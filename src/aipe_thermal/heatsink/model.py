# SPDX-License-Identifier: CC-BY-4.0
"""Plate-fin heat-sink models for natural and forced convection (SI units) -- ``aipe_thermal.heatsink.model``.

Importable, corrected and documented refactor of ``HS_main.py`` (Dr. Fulong Li,
licensed CC BY 4.0, https://creativecommons.org/licenses/by/4.0/). The original
single-script formulation is retained verbatim in the ``legacy_*`` functions at
the end of this module so that every correction can be regression-tested; see
``CHANGELOG.md`` for the file:line of each original defect.

All quantities are SI: metres, kilograms, seconds, kelvin for differences,
``*_C`` suffixes for absolute temperatures in degrees Celsius, pascal, watt.

Geometry (see ``docs/heat_sink_geometry.png``): a rectangular base plate ``W x L``
of thickness ``t_bp`` carries ``n`` straight plate fins of thickness ``t`` and
height ``H`` running along ``L``. ``L`` is the flow direction for forced
convection and the vertical direction for natural convection (vertical base,
vertical fins). ``b`` is the requested minimum gap; the integer fin count is
``n = floor((W - t) / (t + b)) + 1`` (as in the original) and the fins are
spread evenly across ``W`` so the real channel gap is ``b_eff = (W - n t)/(n - 1)``.

Correlations (validity ranges are stated in each docstring and echoed in the
``validity`` list of every evaluation):

* natural convection between vertical isothermal plates: Elenbaas number and
  the Bar-Cohen & Rohsenow (1984) composite ``Nu_b = [576/El^2 + 2.873/El^0.5]^-0.5``;
* forced convection in the fin channels: Teertstra, Yovanovich & Culham (1999)
  composite for developing laminar flow between parallel plates with the
  modified channel Reynolds number ``Re_b* = Re_b b / L``; Gnielinski with the
  hydraulic diameter when the channel flow is transitional or turbulent;
* pressure drop: contraction/expansion (Kays & London asymptotes as used by
  Culham & Muzychka 2001) plus the Muzychka & Yovanovich apparent friction
  factor for developing laminar flow in rectangular ducts; Petukhov friction
  factor for turbulent flow;
* fin efficiency: straight rectangular fin, adiabatic tip (or Harper-Brown
  corrected length), perimeter ``2 (L + t)``;
* base conduction and spreading: Lee, Song, Au & Moran (1995) closed form
  (as popularised by Simons, Electronics Cooling 2004), per device footprint;
* air properties from Sutherland-type fits evaluated at the film temperature.

Nothing here is a CFD or FEM result; the models are analytical screening
correlations with the stated validity ranges.
"""

from __future__ import annotations

import inspect
import json
import math
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Sequence

import numpy as np
from scipy.optimize import NonlinearConstraint, brentq, differential_evolution

from aipe_thermal.units import validate_ambient_C

__version__ = "1.0.0"
__license__ = "CC BY 4.0"
__attribution__ = (
    "Derived from HS_main.py, 'Heat Sink Optimization' by Dr. Fulong Li, CC BY 4.0. "
    "Changes: see CHANGELOG.md (model corrections, SI refactor, Pareto method, tests; packaged as aipe_thermal.heatsink in 0.2.0)."
)

G = 9.80665  # m/s^2
P_ATM = 101325.0  # Pa
R_AIR = 287.05  # J/(kg K)
T0_K = 273.15
SIGMA_SB = 5.670374419e-8  # W/(m^2 K^4)


# =============================================================================
# Fluid properties
# =============================================================================

@dataclass(frozen=True)
class FluidProperties:
    """Thermophysical properties of a fluid at one temperature (SI)."""

    temperature_K: float
    density: float  # kg/m^3
    viscosity: float  # Pa s
    conductivity: float  # W/(m K)
    specific_heat: float  # J/(kg K)
    expansion_coefficient: float  # 1/K
    name: str = "air"

    @property
    def prandtl(self) -> float:
        return self.viscosity * self.specific_heat / self.conductivity

    @property
    def kinematic_viscosity(self) -> float:
        return self.viscosity / self.density

    @property
    def thermal_diffusivity(self) -> float:
        return self.conductivity / (self.density * self.specific_heat)


def air_properties(temperature_K: float, pressure_Pa: float = P_ATM) -> FluidProperties:
    """Dry-air properties at ``temperature_K`` (ideal gas, Sutherland fits).

    * density: ideal gas ``p / (R T)``;
    * viscosity: Sutherland, ``mu0 = 1.716e-5 Pa s`` at 273.15 K, ``S = 110.4 K``;
    * conductivity: Sutherland form with ``k0 = 0.0241 W/(m K)``, ``S = 194 K``;
    * specific heat: quadratic fit ``1002.5 + 2.75e-4 (T - 200)^2``;
    * expansion coefficient: ``1/T`` (ideal gas).

    Agreement with Incropera & DeWitt Table A.4 is within about 0.5 % for
    250 K <= T <= 450 K, the range relevant to electronics cooling.
    """
    T = float(temperature_K)
    if not (150.0 < T < 1500.0):
        raise ValueError(f"air_properties: temperature {T} K outside 150-1500 K")
    density = pressure_Pa / (R_AIR * T)
    viscosity = 1.716e-5 * (T / T0_K) ** 1.5 * (T0_K + 110.4) / (T + 110.4)
    conductivity = 0.0241 * (T / T0_K) ** 1.5 * (T0_K + 194.0) / (T + 194.0)
    specific_heat = 1002.5 + 2.75e-4 * (T - 200.0) ** 2
    return FluidProperties(T, density, viscosity, conductivity, specific_heat, 1.0 / T, "air")


def film_temperature_K(t_ambient_C: float, delta_T_K: float) -> float:
    """Film temperature ``T_amb + dT/2`` in kelvin."""
    return t_ambient_C + T0_K + 0.5 * delta_T_K


# =============================================================================
# Materials
# =============================================================================

@dataclass(frozen=True)
class Material:
    """Isotropic solid with thermal conductivity (W/m K) and density (kg/m^3)."""

    name: str
    conductivity: float
    density: float
    note: str = ""


MATERIALS: dict[str, Material] = {
    "Al-6063-T5": Material("Al-6063-T5", 209.0, 2700.0, "extrusion alloy; ASM typical 209 W/(m K)"),
    "Al-6061-T6": Material("Al-6061-T6", 167.0, 2700.0, "ASM typical 167 W/(m K)"),
    "Al-pure": Material("Al-pure", 237.0, 2702.0, "Incropera & DeWitt Table A.1 at 300 K"),
    "Cu-C11000": Material("Cu-C11000", 391.0, 8940.0, "ETP copper; ASM typical 391 W/(m K)"),
    "Cu-pure": Material("Cu-pure", 401.0, 8933.0, "Incropera & DeWitt Table A.1 at 300 K"),
}
MATERIAL_ALIASES = {"aluminium": "Al-6063-T5", "aluminum": "Al-6063-T5", "al": "Al-6063-T5",
                    "copper": "Cu-C11000", "cu": "Cu-C11000"}


def get_material(name: str | Material) -> Material:
    """Look up a material by key or alias (case-insensitive)."""
    if isinstance(name, Material):
        return name
    key = str(name).strip()
    if key in MATERIALS:
        return MATERIALS[key]
    alias = MATERIAL_ALIASES.get(key.lower())
    if alias:
        return MATERIALS[alias]
    for k in MATERIALS:
        if k.lower() == key.lower():
            return MATERIALS[k]
    raise KeyError(f"unknown material {name!r}; known: {sorted(MATERIALS)} or aliases {sorted(MATERIAL_ALIASES)}")


# =============================================================================
# Geometry
# =============================================================================

@dataclass(frozen=True)
class PlateFinHeatSink:
    """Plate-fin heat sink geometry (metres) and material.

    ``fin_gap`` is the requested minimum gap; ``gap_effective`` is the gap after
    spreading the integer number of fins evenly across ``width``.
    """

    length: float  # L, along the fins (flow / vertical direction)
    width: float  # W, across the fins
    fin_height: float  # H
    fin_thickness: float  # t
    fin_gap: float  # b (requested)
    base_thickness: float  # t_bp
    material: Material = MATERIALS["Al-6063-T5"]

    def __post_init__(self) -> None:
        for name in ("length", "width", "fin_height", "fin_thickness", "fin_gap", "base_thickness"):
            v = getattr(self, name)
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
                raise ValueError(f"PlateFinHeatSink.{name} must be a finite positive number, got {v!r}")
        if self.width < 2 * self.fin_thickness + self.fin_gap:
            raise ValueError("width must accommodate at least two fins and one gap: W >= 2 t + b")
        object.__setattr__(self, "material", get_material(self.material))

    # -- derived geometry -------------------------------------------------
    @property
    def n_fins(self) -> int:
        """Integer fin count ``floor((W - t)/(t + b)) + 1`` (>= 2), as in the original script."""
        return max(2, int(math.floor((self.width - self.fin_thickness) / (self.fin_thickness + self.fin_gap))) + 1)

    @property
    def n_channels(self) -> int:
        return self.n_fins - 1

    @property
    def gap_effective(self) -> float:
        """Actual channel gap when ``n_fins`` fins are spread evenly over ``width``."""
        return (self.width - self.n_fins * self.fin_thickness) / (self.n_fins - 1)

    @property
    def base_area(self) -> float:
        return self.width * self.length

    @property
    def exposed_base_area(self) -> float:
        """Base surface between the fins that is wetted by the fluid."""
        return (self.width - self.n_fins * self.fin_thickness) * self.length

    @property
    def fin_perimeter(self) -> float:
        """Cross-section perimeter ``2 (L + t)`` used in the fin parameter."""
        return 2.0 * (self.length + self.fin_thickness)

    @property
    def fin_cross_section(self) -> float:
        return self.length * self.fin_thickness

    @property
    def fin_area_each(self) -> float:
        """Convective area of one fin ``P H`` (both faces plus the two edges)."""
        return self.fin_perimeter * self.fin_height

    @property
    def channel_area(self) -> float:
        """Total free-flow cross-section of all channels ``(n-1) b_eff H``."""
        return self.n_channels * self.gap_effective * self.fin_height

    @property
    def frontal_area(self) -> float:
        """Frontal area of the finned region ``W H`` seen by a ducted flow."""
        return self.width * self.fin_height

    @property
    def free_area_ratio(self) -> float:
        """sigma = channel area / frontal area = 1 - n t / W."""
        return self.channel_area / self.frontal_area

    @property
    def hydraulic_diameter(self) -> float:
        """Hydraulic diameter of one rectangular channel ``2 b H / (b + H)``."""
        b, H = self.gap_effective, self.fin_height
        return 2.0 * b * H / (b + H)

    @property
    def channel_aspect_ratio(self) -> float:
        """Short side over long side of one channel (0 < eps <= 1)."""
        b, H = self.gap_effective, self.fin_height
        return min(b, H) / max(b, H)

    @property
    def mass(self) -> float:
        """Mass of base plate plus fins (kg)."""
        return self.material.density * (self.base_area * self.base_thickness
                                        + self.n_fins * self.fin_cross_section * self.fin_height)

    @property
    def envelope_volume(self) -> float:
        """Bounding-box volume ``W L (t_bp + H)`` (m^3)."""
        return self.base_area * (self.base_thickness + self.fin_height)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["material"] = self.material.name
        d.update(n_fins=self.n_fins, gap_effective=self.gap_effective, mass_kg=self.mass,
                 envelope_volume_m3=self.envelope_volume, hydraulic_diameter=self.hydraulic_diameter,
                 base_area=self.base_area)
        return d


# =============================================================================
# Fin model
# =============================================================================

def fin_parameter(h: float, conductivity: float, fin_thickness: float, fin_length: float) -> float:
    """Fin parameter ``m = sqrt(h P / (k A_c))`` with ``P = 2 (L + t)`` and ``A_c = L t``.

    The original used ``m = sqrt(2 h / (k t))`` which is the ``t << L`` limit
    (perimeter ``2 L``); the full perimeter is used here.
    """
    if h < 0 or conductivity <= 0 or fin_thickness <= 0 or fin_length <= 0:
        raise ValueError("fin_parameter needs h >= 0 and positive k, t, L")
    perimeter = 2.0 * (fin_length + fin_thickness)
    return math.sqrt(h * perimeter / (conductivity * fin_length * fin_thickness))


def fin_efficiency(h: float, conductivity: float, fin_height: float, fin_thickness: float,
                   fin_length: float, tip: str = "adiabatic") -> float:
    """Efficiency ``tanh(m H_c) / (m H_c)`` of a straight rectangular fin.

    ``tip="adiabatic"`` uses ``H_c = H``; ``tip="convective"`` uses the
    Harper-Brown corrected length ``H_c = H + t/2``. Returns 1 for ``m H -> 0``.
    """
    m = fin_parameter(h, conductivity, fin_thickness, fin_length)
    hc = fin_height + (0.5 * fin_thickness if tip == "convective" else 0.0)
    if tip not in ("adiabatic", "convective"):
        raise ValueError("tip must be 'adiabatic' or 'convective'")
    mh = m * hc
    if mh < 1e-9:
        return 1.0
    return math.tanh(mh) / mh


def single_fin_resistance(h: float, conductivity: float, fin_height: float, fin_thickness: float,
                          fin_length: float, tip: str = "adiabatic") -> float:
    """Base-to-fluid resistance of one fin ``1 / (h P H_c eta)`` (K/W).

    For an adiabatic tip this equals the textbook ``1 / (sqrt(h P k A_c) tanh(m H))``.
    """
    eta = fin_efficiency(h, conductivity, fin_height, fin_thickness, fin_length, tip)
    hc = fin_height + (0.5 * fin_thickness if tip == "convective" else 0.0)
    area = 2.0 * (fin_length + fin_thickness) * hc
    return 1.0 / (h * area * eta)


# =============================================================================
# Natural convection (vertical parallel plates)
# =============================================================================

def elenbaas_number(gap: float, length: float, delta_T_K: float, props: FluidProperties) -> float:
    """Elenbaas number ``El = rho^2 g beta c_p b^4 dT / (mu k L) = Ra_b (b / L)``.

    ``gap`` is the plate spacing, ``length`` the vertical plate length and
    ``delta_T_K`` the plate-to-ambient temperature difference.
    """
    if gap <= 0 or length <= 0 or delta_T_K <= 0:
        raise ValueError("elenbaas_number needs positive gap, length and delta_T")
    return (props.density ** 2 * G * props.expansion_coefficient * props.specific_heat
            * gap ** 4 * delta_T_K) / (props.viscosity * props.conductivity * length)


def bar_cohen_rohsenow_nusselt(elenbaas: float) -> float:
    """Composite Nusselt number ``Nu_b = h b / k = [576 / El^2 + 2.873 / El^0.5]^-0.5``.

    Bar-Cohen & Rohsenow (1984) for symmetric isothermal vertical plates. It
    blends the fully developed limit ``El/24`` (small El) and the isolated
    plate limit ``0.59 El^0.25`` (large El) and is quoted to about +/-10 %
    across the whole range; the physical assumptions are laminar flow,
    isothermal plates, ``b << L``, open channel ends and no radiation.
    """
    if elenbaas <= 0:
        raise ValueError("Elenbaas number must be positive")
    return (576.0 / elenbaas ** 2 + 2.873 / math.sqrt(elenbaas)) ** -0.5


def natural_convection_h(gap: float, length: float, delta_T_K: float, props: FluidProperties) -> tuple[float, float]:
    """Return ``(h, El)`` for vertical plates ``gap`` apart with length ``length``."""
    el = elenbaas_number(gap, length, delta_T_K, props)
    return bar_cohen_rohsenow_nusselt(el) * props.conductivity / gap, el


def bar_cohen_optimum_spacing(length: float, delta_T_K: float, props: FluidProperties) -> float:
    """Optimum plate spacing ``b_opt = 2.714 L / Ra_L^0.25`` (maximum heat transfer per base area).

    Equivalent to ``El(b_opt) = 2.714^4 = 54.3``.
    """
    ra_l = G * props.expansion_coefficient * delta_T_K * length ** 3 / (props.kinematic_viscosity * props.thermal_diffusivity)
    return 2.714 * length / ra_l ** 0.25


def radiation_h(emissivity: float, t_surface_C: float, t_ambient_C: float) -> float:
    """Linearised radiation coefficient ``eps sigma (Ts^2 + Ta^2)(Ts + Ta)`` to large surroundings."""
    if not (0.0 <= emissivity <= 1.0):
        raise ValueError("emissivity must lie in [0, 1]")
    ts, ta = t_surface_C + T0_K, t_ambient_C + T0_K
    return emissivity * SIGMA_SB * (ts ** 2 + ta ** 2) * (ts + ta)


# =============================================================================
# Forced convection and pressure drop
# =============================================================================

def tyc_channel_nusselt(re_b_star: float, prandtl: float) -> float:
    """Teertstra-Yovanovich-Culham (1999) composite Nusselt number for a plate-fin channel.

    ``Nu_b = h b / k = [ (Re_b* Pr / 2)^-3 + (0.664 Re_b*^0.5 Pr^(1/3) (1 + 3.65 / Re_b*^0.5)^0.5)^-3 ]^(-1/3)``
    with ``Re_b* = Re_b b / L`` and ``Re_b = rho V_ch b / mu``. ``h`` is referenced to
    the inlet (ambient) temperature so ``R = 1/(h A)`` already includes the air
    temperature rise (fully developed limit ``Re_b* Pr / 2`` is the energy balance).
    Developed for laminar developing flow between isothermal parallel plates
    (``H >> b``); the authors validated it against numerical data over roughly
    ``0.1 <= Re_b* <= 100`` and against measurements on commercial sinks.
    """
    if re_b_star <= 0 or prandtl <= 0:
        raise ValueError("Re_b* and Pr must be positive")
    developed = re_b_star * prandtl / 2.0
    boundary_layer = 0.664 * math.sqrt(re_b_star) * prandtl ** (1.0 / 3.0) * math.sqrt(1.0 + 3.65 / math.sqrt(re_b_star))
    return (developed ** -3 + boundary_layer ** -3) ** (-1.0 / 3.0)


def laminar_duct_nusselt(aspect_ratio: float, boundary: str = "T") -> float:
    """Fully developed laminar Nusselt number (``D_h`` based) in a rectangular duct.

    Shah & London polynomial fits in the aspect ratio ``eps = short/long`` (0..1):
    ``boundary="T"`` (isothermal walls): ``7.541 (1 - 2.610 e + 4.970 e^2 - 5.119 e^3 + 2.702 e^4 - 0.548 e^5)``
    ``boundary="H"`` (uniform heat flux): ``8.235 (1 - 2.0421 e + 3.0853 e^2 - 2.4765 e^3 + 1.0578 e^4 - 0.1861 e^5)``
    giving 7.541/8.235 for parallel plates and 2.976/3.608 for a square duct.
    """
    e = aspect_ratio
    if not (0.0 <= e <= 1.0):
        raise ValueError("aspect_ratio must lie in [0, 1]")
    if boundary == "T":
        return 7.541 * (1 - 2.610 * e + 4.970 * e ** 2 - 5.119 * e ** 3 + 2.702 * e ** 4 - 0.548 * e ** 5)
    if boundary == "H":
        return 8.235 * (1 - 2.0421 * e + 3.0853 * e ** 2 - 2.4765 * e ** 3 + 1.0578 * e ** 4 - 0.1861 * e ** 5)
    raise ValueError("boundary must be 'T' or 'H'")


def laminar_duct_friction_fanning_re(aspect_ratio: float) -> float:
    """Fully developed laminar ``f_Fanning Re_Dh`` for a rectangular duct (Shah & London).

    ``24 (1 - 1.3553 e + 1.9467 e^2 - 1.7012 e^3 + 0.9564 e^4 - 0.2537 e^5)``:
    24 for parallel plates (Darcy 96), 14.23 for a square duct (Darcy 56.9).
    """
    e = aspect_ratio
    if not (0.0 <= e <= 1.0):
        raise ValueError("aspect_ratio must lie in [0, 1]")
    return 24.0 * (1 - 1.3553 * e + 1.9467 * e ** 2 - 1.7012 * e ** 3 + 0.9564 * e ** 4 - 0.2537 * e ** 5)


def petukhov_friction_factor(reynolds: float) -> float:
    """Darcy friction factor ``(0.790 ln Re - 1.64)^-2`` for smooth ducts, ``3000 <= Re <= 5e6``."""
    if reynolds <= 1000:
        raise ValueError("Petukhov friction factor needs Re > 1000")
    return (0.790 * math.log(reynolds) - 1.64) ** -2


def gnielinski_nusselt(reynolds: float, prandtl: float) -> float:
    """Gnielinski correlation ``Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 (f/8)^0.5 (Pr^(2/3) - 1))``.

    Valid for ``3000 <= Re <= 5e6`` and ``0.5 <= Pr <= 2000`` with the Petukhov
    Darcy friction factor; ``h`` is referenced to the bulk mean fluid temperature.
    """
    if reynolds <= 1000:
        raise ValueError("Gnielinski needs Re > 1000 (laminar flow uses laminar_duct_nusselt)")
    f = petukhov_friction_factor(reynolds)
    return (f / 8.0) * (reynolds - 1000.0) * prandtl / (1.0 + 12.7 * math.sqrt(f / 8.0) * (prandtl ** (2.0 / 3.0) - 1.0))


def dittus_boelter_nusselt(reynolds: float, prandtl: float, heating: bool = True) -> float:
    """Dittus-Boelter ``Nu = 0.023 Re^0.8 Pr^n`` (``n = 0.4`` heating, 0.3 cooling).

    Valid for fully developed turbulent flow, ``Re >= 1e4``, ``0.6 <= Pr <= 160``, ``L/D >= 10``.
    """
    if reynolds <= 0 or prandtl <= 0:
        raise ValueError("Re and Pr must be positive")
    return 0.023 * reynolds ** 0.8 * prandtl ** (0.4 if heating else 0.3)


def developing_laminar_fanning_friction(reynolds_dh: float, length: float, hydraulic_diameter: float,
                                        aspect_ratio: float) -> float:
    """Apparent Fanning friction factor for developing laminar flow in a rectangular duct.

    Muzychka & Yovanovich: ``f_app Re = sqrt((3.44 / L*^0.5)^2 + (f Re)_fd^2)`` with
    ``L* = L / (D_h Re)``; reduces to the fully developed value for long ducts.
    """
    if reynolds_dh <= 0:
        raise ValueError("Re must be positive")
    l_star = length / (hydraulic_diameter * reynolds_dh)
    f_re_fd = laminar_duct_friction_fanning_re(aspect_ratio)
    return math.sqrt((3.44 / math.sqrt(l_star)) ** 2 + f_re_fd ** 2) / reynolds_dh


def channel_pressure_drop(sink: PlateFinHeatSink, volumetric_flow: float, props: FluidProperties,
                          include_entrance_exit: bool = True) -> dict[str, float]:
    """Pressure drop across the fin channels for a fully ducted flow (Pa).

    ``dp = (K_c + K_e + 4 f_app L / D_h) rho V_ch^2 / 2`` with ``K_c = 0.42 (1 - sigma^2)``,
    ``K_e = (1 - sigma^2)^2`` (Kays & London laminar asymptotes as used by Culham &
    Muzychka 2001) and the developing-flow apparent friction factor; the Petukhov
    Darcy factor is used when ``Re_Dh >= 2300``. All flow is assumed to pass
    through the channels (no bypass).
    """
    if volumetric_flow <= 0:
        return {"velocity": 0.0, "reynolds_dh": 0.0, "friction": 0.0, "contraction": 0.0, "expansion": 0.0,
                "total": 0.0, "fanning_friction_factor": 0.0, "regime": "no-flow"}
    v = volumetric_flow / sink.channel_area
    dh = sink.hydraulic_diameter
    re = props.density * v * dh / props.viscosity
    dyn = 0.5 * props.density * v ** 2
    if re < 2300.0:
        f_fanning = developing_laminar_fanning_friction(re, sink.length, dh, sink.channel_aspect_ratio)
        regime = "laminar"
    else:
        f_fanning = petukhov_friction_factor(re) / 4.0
        regime = "turbulent"
    friction = 4.0 * f_fanning * sink.length / dh * dyn
    sigma = sink.free_area_ratio
    kc = 0.42 * (1.0 - sigma ** 2) if include_entrance_exit else 0.0
    ke = (1.0 - sigma ** 2) ** 2 if include_entrance_exit else 0.0
    return {"velocity": v, "reynolds_dh": re, "friction": friction, "contraction": kc * dyn,
            "expansion": ke * dyn, "total": friction + (kc + ke) * dyn, "fanning_friction_factor": f_fanning,
            "regime": regime}


def fan_curve_pressure(fan_curve: Sequence[tuple[float, float]], volumetric_flow: float) -> float:
    """Linear interpolation of a fan curve given as ``[(Q m^3/s, dp Pa), ...]`` sorted by Q."""
    pts = sorted((float(q), float(p)) for q, p in fan_curve)
    if len(pts) < 2:
        raise ValueError("fan curve needs at least two (flow, pressure) points")
    q = volumetric_flow
    if q <= pts[0][0]:
        return pts[0][1]
    if q >= pts[-1][0]:
        return pts[-1][1]
    for (q0, p0), (q1, p1) in zip(pts[:-1], pts[1:]):
        if q0 <= q <= q1:
            return p0 + (p1 - p0) * (q - q0) / (q1 - q0)
    return pts[-1][1]


def fan_operating_point(sink: PlateFinHeatSink, fan_curve: Sequence[tuple[float, float]], props: FluidProperties,
                        n_fans_parallel: int = 1) -> dict[str, float]:
    """Flow where the fan (or ``n`` identical fans in parallel) matches the sink pressure drop.

    Solves ``dp_fan(Q / n) = dp_sink(Q)`` with Brent's method on ``[0, n Q_max]``.
    """
    pts = sorted((float(q), float(p)) for q, p in fan_curve)
    q_max = pts[-1][0] * n_fans_parallel
    if q_max <= 0:
        raise ValueError("fan curve must extend to a positive flow")

    def g(q: float) -> float:
        return fan_curve_pressure(pts, q / n_fans_parallel) - channel_pressure_drop(sink, q, props)["total"]

    if g(q_max) >= 0:  # fan still has pressure at free delivery point -> take the curve end
        q_op = q_max
    else:
        q_op = brentq(g, 1e-9 * q_max, q_max, xtol=1e-12, rtol=1e-10, maxiter=200)
    dp = channel_pressure_drop(sink, q_op, props)["total"]
    return {"volumetric_flow": q_op, "pressure_drop": dp, "fan_pressure": fan_curve_pressure(pts, q_op / n_fans_parallel)}


def fan_power(volumetric_flow: float, pressure_drop: float, efficiency: float = 0.25) -> dict[str, float]:
    """Air power ``Q dp`` and electrical input ``Q dp / eta`` (small axial fans: eta ~ 0.15-0.35)."""
    if not (0.0 < efficiency <= 1.0):
        raise ValueError("fan efficiency must be in (0, 1]")
    air = volumetric_flow * pressure_drop
    return {"air_power_W": air, "electrical_power_W": air / efficiency, "efficiency": efficiency}


# =============================================================================
# Spreading resistance (Lee, Song, Au & Moran 1995)
# =============================================================================

def spreading_resistance(source_area: float, plate_area: float, plate_thickness: float,
                         conductivity: float, h_effective: float) -> dict[str, float]:
    """Closed-form spreading plus one-dimensional base resistance for a centred heat source.

    Lee, Song, Au & Moran (1995), in the area form of Simons (2004):
    ``r1 = sqrt(A_s/pi)``, ``r2 = sqrt(A_p/pi)``, ``eps = r1/r2``, ``tau = t/r2``, ``Bi = h r2/k``,
    ``lambda = pi + 1/(eps sqrt(pi))``,
    ``Phi = (tanh(lambda tau) + lambda/Bi) / (1 + (lambda/Bi) tanh(lambda tau))``,
    ``Psi_max = eps tau/sqrt(pi) + (1 - eps) Phi/sqrt(pi)``,
    ``Psi_avg = eps tau/sqrt(pi) + 0.5 (1 - eps)^1.5 Phi``,
    ``R = Psi / (k r1 sqrt(pi))``.

    ``Psi`` includes the one-dimensional term ``t/(k A_p)``; the returned
    ``spreading_*`` values have it removed. ``h_effective`` is the uniform
    heat-transfer coefficient on the back face (for a finned base use
    ``1/(R_fins A_p)``). Limits: ``A_s -> A_p`` gives pure 1-D conduction;
    ``A_p, t -> inf`` gives the isoflux half-space maximum ``1/(pi k r1)``.
    Accuracy quoted by the authors: within ~10 % of the exact series solution.
    """
    if min(source_area, plate_area, plate_thickness, conductivity, h_effective) <= 0:
        raise ValueError("spreading_resistance needs positive areas, thickness, conductivity and h")
    a_s = min(source_area, plate_area)
    r1 = math.sqrt(a_s / math.pi)
    r2 = math.sqrt(plate_area / math.pi)
    eps = r1 / r2
    tau = plate_thickness / r2
    bi = h_effective * r2 / conductivity
    lam = math.pi + 1.0 / (eps * math.sqrt(math.pi))
    th = math.tanh(lam * tau)
    ratio = lam / bi
    phi = (th + ratio) / (1.0 + ratio * th)
    psi_max = eps * tau / math.sqrt(math.pi) + (1.0 - eps) * phi / math.sqrt(math.pi)
    psi_avg = eps * tau / math.sqrt(math.pi) + 0.5 * (1.0 - eps) ** 1.5 * phi
    denom = conductivity * r1 * math.sqrt(math.pi)
    one_d = plate_thickness / (conductivity * plate_area)
    total_max = psi_max / denom
    total_avg = psi_avg / denom
    return {"one_d": one_d, "spreading_max": max(total_max - one_d, 0.0), "spreading_avg": max(total_avg - one_d, 0.0),
            "total_max": total_max, "total_avg": total_avg, "epsilon": eps, "tau": tau, "biot": bi, "phi": phi}


# =============================================================================
# Heat sources and composite resistance
# =============================================================================

@dataclass(frozen=True)
class HeatSource:
    """Device footprint on the base plate: contact area (m^2) and dissipated power (W)."""

    area: float
    power: float
    name: str = ""

    def __post_init__(self) -> None:
        if self.area <= 0 or self.power < 0:
            raise ValueError("HeatSource needs area > 0 and power >= 0")


def heat_sink_resistance(sink: PlateFinHeatSink, h: float, sources: Sequence[HeatSource] | None = None,
                         tip: str = "adiabatic", h_radiation: float = 0.0,
                         caloric_resistance: float = 0.0) -> dict[str, Any]:
    """Compose the base-to-fluid resistance of the finned surface with base conduction and spreading.

    * ``r_convective = 1/((h + h_rad) A_eff)`` with ``A_eff = A_base,exposed + n eta P H`` --
      the fin-root (base surface) to fluid resistance; ``caloric_resistance`` is
      added for correlations referenced to the bulk mean rather than inlet temperature;
    * for each source the base plate is partitioned into a cell of area
      ``A_base P_i / sum P`` (power-weighted so the shared sink node is consistent)
      and the Lee et al. closed form gives ``r_base_i = one-D + spreading`` with the
      back-face coefficient ``h_eff = 1/(r_convective A_base)``;
    * without sources the base is uniformly heated: ``r_base = t_bp/(k A_base)``.

    ``r_total`` is ``(T_hottest_source - T_fluid) / sum P`` =
    ``r_convective + max_i[P_i r_base_i] / sum P``.
    """
    eta = fin_efficiency(h, sink.material.conductivity, sink.fin_height, sink.fin_thickness, sink.length, tip)
    hc = sink.fin_height + (0.5 * sink.fin_thickness if tip == "convective" else 0.0)
    fin_area = sink.n_fins * sink.fin_perimeter * hc
    a_eff = sink.exposed_base_area + eta * fin_area
    h_total = h + h_radiation
    r_conv = 1.0 / (h_total * a_eff) + caloric_resistance
    h_eff = 1.0 / (r_conv * sink.base_area)
    k = sink.material.conductivity
    result: dict[str, Any] = {"h": h, "h_radiation": h_radiation, "fin_efficiency": eta, "effective_area": a_eff,
                              "wetted_area": sink.exposed_base_area + fin_area, "r_convective": r_conv,
                              "h_effective_base": h_eff, "sources": []}
    if not sources:
        r_base = sink.base_thickness / (k * sink.base_area)
        result.update(r_base_1d=r_base, r_spreading_max=0.0, r_total=r_conv + r_base, hottest_source=None)
        return result
    p_tot = sum(s.power for s in sources)
    if p_tot <= 0:
        raise ValueError("sources must dissipate positive total power")
    worst_rise = 0.0
    worst_name = None
    for s in sources:
        cell = sink.base_area * s.power / p_tot if s.power > 0 else sink.base_area / len(sources)
        if s.area > cell:
            spread = {"one_d": sink.base_thickness / (k * cell), "spreading_max": 0.0, "spreading_avg": 0.0,
                      "total_max": sink.base_thickness / (k * cell), "note": "footprint exceeds its base cell"}
        else:
            spread = spreading_resistance(s.area, cell, sink.base_thickness, k, h_eff)
        rise = s.power * spread["total_max"]
        result["sources"].append({"name": s.name, "power": s.power, "area": s.area, "cell_area": cell,
                                  "r_base_1d": spread["one_d"], "r_spreading_max": spread["spreading_max"],
                                  "r_base_total": spread["total_max"], "base_rise_K": rise})
        if rise > worst_rise:
            worst_rise, worst_name = rise, s.name
    result.update(r_base_1d=max(x["r_base_1d"] for x in result["sources"]),
                  r_spreading_max=max(x["r_spreading_max"] for x in result["sources"]),
                  r_total=r_conv + worst_rise / p_tot, hottest_source=worst_name)
    return result


# =============================================================================
# Operating point and evaluation
# =============================================================================

@dataclass
class OperatingPoint:
    """Thermal boundary conditions for one evaluation.

    * ``cooling``: ``"natural"`` (vertical fins, buoyancy) or ``"forced"`` (ducted flow along the fins);
    * exactly one of ``delta_T_K`` (fixed base-to-ambient rise) or ``heat_W`` (solve for the rise);
    * forced convection needs ``volumetric_flow`` (m^3/s through the channels) or ``fan_curve``
      ``[(Q, dp), ...]`` with ``n_fans``; ``fan_efficiency`` converts air power to electrical power;
    * ``sources`` are device footprints for spreading (optional; uniform base heating otherwise);
    * ``emissivity`` > 0 adds a linearised radiation term on the wetted area (default off).
    """

    t_ambient_C: float = 40.0
    cooling: str = "natural"
    heat_W: float | None = None
    delta_T_K: float | None = None
    volumetric_flow: float | None = None
    fan_curve: Sequence[tuple[float, float]] | None = None
    n_fans: int = 1
    fan_efficiency: float = 0.25
    sources: Sequence[HeatSource] = field(default_factory=tuple)
    emissivity: float = 0.0
    fin_tip: str = "adiabatic"
    pressure_Pa: float = P_ATM

    def __post_init__(self) -> None:
        if self.cooling not in ("natural", "forced"):
            raise ValueError("cooling must be 'natural' or 'forced'")
        self.t_ambient_C = validate_ambient_C(self.t_ambient_C, "t_ambient_C")
        if (self.heat_W is None) == (self.delta_T_K is None):
            raise ValueError("specify exactly one of heat_W or delta_T_K")
        if self.heat_W is not None and self.heat_W <= 0:
            raise ValueError("heat_W must be positive")
        if self.delta_T_K is not None and self.delta_T_K <= 0:
            raise ValueError("delta_T_K must be positive")
        if self.cooling == "forced" and self.volumetric_flow is None and not self.fan_curve:
            raise ValueError("forced convection needs volumetric_flow or fan_curve")
        if self.volumetric_flow is not None and self.volumetric_flow <= 0:
            raise ValueError("volumetric_flow must be positive")
        if self.sources and self.heat_W is not None:
            p_src = sum(s.power for s in self.sources)
            if p_src > 0 and abs(p_src - self.heat_W) > 1e-6 * max(1.0, self.heat_W):
                raise ValueError(f"sum of source powers {p_src} W must equal heat_W {self.heat_W} W")

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["sources"] = [asdict(s) for s in self.sources]
        d["fan_curve"] = [list(p) for p in self.fan_curve] if self.fan_curve else None
        return d


@dataclass
class HeatSinkPerformance:
    """Result of evaluating a heat sink at an operating point (all SI)."""

    sink: dict[str, Any]
    operating_point: dict[str, Any]
    delta_T_base_K: float
    t_base_C: float
    t_film_K: float
    heat_W: float
    r_total: float
    r_convective: float
    r_base_1d: float
    r_spreading_max: float
    h: float
    h_radiation: float
    fin_efficiency: float
    effective_area_m2: float
    mass_kg: float
    envelope_volume_m3: float
    elenbaas: float | None
    reynolds_dh: float | None
    reynolds_b_star: float | None
    regime: str
    volumetric_flow: float | None
    channel_velocity: float | None
    pressure_drop_Pa: float | None
    pressure_drop_breakdown: dict[str, float] | None
    fan_power_W: float | None
    fan_electrical_power_W: float | None
    sources: list[dict[str, Any]]
    hottest_source: str | None
    solver: dict[str, Any]
    validity: list[str]
    warnings: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, **kwargs: Any) -> str:
        kwargs.setdefault("indent", 2)
        return json.dumps(self.to_dict(), default=_json_default, **kwargs)


def _json_default(obj: Any) -> Any:
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if isinstance(obj, Material):
        return obj.name
    raise TypeError(f"not JSON serialisable: {type(obj)!r}")


def evaluate_at_delta_T(sink: PlateFinHeatSink, op: OperatingPoint, delta_T_K: float,
                        solver: dict[str, Any] | None = None) -> HeatSinkPerformance:
    """Evaluate the sink for a given base-to-ambient rise (air properties at the film temperature)."""
    if delta_T_K <= 0:
        raise ValueError("delta_T_K must be positive")
    t_film = film_temperature_K(op.t_ambient_C, delta_T_K)
    props = air_properties(t_film, op.pressure_Pa)
    validity: list[str] = []
    warnings: list[str] = []
    h_rad = radiation_h(op.emissivity, op.t_ambient_C + delta_T_K, op.t_ambient_C) if op.emissivity > 0 else 0.0
    if op.emissivity > 0:
        validity.append("radiation: linearised grey-body exchange on the whole wetted area to surroundings at ambient; "
                        "channel view factors neglected (upper bound for the radiative share)")
    el = re_dh = re_b_star = flow = v_ch = dp_total = fan_air = fan_el = None
    dp: dict[str, float] | None = None
    caloric = 0.0
    if op.cooling == "natural":
        b = sink.gap_effective
        h, el = natural_convection_h(b, sink.length, delta_T_K, props)
        regime = "natural-laminar"
        validity.append("natural convection: Bar-Cohen & Rohsenow composite for vertical isothermal parallel plates "
                        "(laminar, b << L, open channel ends); fins vertical, base vertical; air at film temperature")
        if el < 0.1 or el > 1e5:
            warnings.append(f"Elenbaas number {el:.3g} outside the 0.1-1e5 band where the composite is normally applied")
        if b > 0.25 * sink.length:
            warnings.append("fin gap exceeds L/4: the parallel-plate channel assumption (b << L) is weak")
    else:
        if op.fan_curve:
            fop = fan_operating_point(sink, op.fan_curve, props, op.n_fans)
            flow = fop["volumetric_flow"]
        else:
            flow = float(op.volumetric_flow)  # type: ignore[arg-type]
        v_ch = flow / sink.channel_area
        b = sink.gap_effective
        re_b = props.density * v_ch * b / props.viscosity
        re_b_star = re_b * b / sink.length
        re_dh = props.density * v_ch * sink.hydraulic_diameter / props.viscosity
        if re_dh < 2300.0:
            h = tyc_channel_nusselt(re_b_star, props.prandtl) * props.conductivity / b
            regime = "forced-laminar"
            validity.append("forced convection: Teertstra-Yovanovich-Culham composite (developing laminar flow between "
                            "parallel plates, h referenced to inlet temperature); fully ducted flow, no bypass")
            if not (0.1 <= re_b_star <= 100.0):
                warnings.append(f"Re_b* = {re_b_star:.3g} outside the ~0.1-100 range validated by Teertstra et al.")
            if sink.channel_aspect_ratio > 0.5:
                warnings.append("channel aspect ratio b/H > 0.5: parallel-plate assumption of the TYC model is weak")
        else:
            h = gnielinski_nusselt(re_dh, props.prandtl) * props.conductivity / sink.hydraulic_diameter
            m_dot = props.density * flow
            caloric = 1.0 / (2.0 * m_dot * props.specific_heat)
            regime = "forced-turbulent"
            validity.append("forced convection: Gnielinski on the channel hydraulic diameter (Re_Dh >= 2300) with the "
                            "mean-fluid caloric term 1/(2 m_dot c_p) added so the result is referenced to the inlet")
            if re_dh < 3000.0:
                warnings.append(f"Re_Dh = {re_dh:.0f} is transitional (2300-3000); Gnielinski extrapolated")
        dp = channel_pressure_drop(sink, flow, props)
        dp_total = dp["total"]
        fp = fan_power(flow, dp_total, op.fan_efficiency)
        fan_air, fan_el = fp["air_power_W"], fp["electrical_power_W"]
        validity.append("pressure drop: Kays & London contraction/expansion asymptotes plus Muzychka-Yovanovich developing "
                        "laminar friction (Petukhov when turbulent); fan power = Q dp / eta_fan")
    comp = heat_sink_resistance(sink, h, op.sources, op.fin_tip, h_rad, caloric)
    if comp["fin_efficiency"] < 0.8:
        warnings.append(f"fin efficiency {comp['fin_efficiency']:.2f} < 0.8: isothermal-plate correlations lose accuracy")
    validity.append("fin efficiency: straight rectangular fin, perimeter 2(L + t), tip '" + op.fin_tip + "'")
    validity.append("base: 1-D conduction plus Lee-Song-Au-Moran spreading per footprint (power-weighted base cells); "
                    "uniform base heating when no sources are given")
    heat = op.heat_W if op.heat_W is not None else delta_T_K / comp["r_total"]
    return HeatSinkPerformance(
        sink=sink.to_dict(), operating_point=op.to_dict(), delta_T_base_K=delta_T_K,
        t_base_C=op.t_ambient_C + delta_T_K, t_film_K=t_film, heat_W=heat, r_total=comp["r_total"],
        r_convective=comp["r_convective"], r_base_1d=comp["r_base_1d"], r_spreading_max=comp["r_spreading_max"],
        h=h, h_radiation=h_rad, fin_efficiency=comp["fin_efficiency"], effective_area_m2=comp["effective_area"],
        mass_kg=sink.mass, envelope_volume_m3=sink.envelope_volume, elenbaas=el, reynolds_dh=re_dh,
        reynolds_b_star=re_b_star, regime=regime, volumetric_flow=flow, channel_velocity=v_ch,
        pressure_drop_Pa=dp_total, pressure_drop_breakdown=dp, fan_power_W=fan_air, fan_electrical_power_W=fan_el,
        sources=comp["sources"], hottest_source=comp["hottest_source"],
        solver=solver or {"mode": "fixed-delta-T"}, validity=validity, warnings=warnings)


def evaluate(sink: PlateFinHeatSink, op: OperatingPoint, tol_K: float = 1e-4, max_iter: int = 100,
             relaxation: float = 0.5) -> HeatSinkPerformance:
    """Evaluate at a fixed ``delta_T_K`` or solve ``dT = Q R(dT)`` for a given heat load.

    The heat-load solve is a damped fixed-point iteration
    ``dT <- (1 - w) dT + w Q R(dT)`` (``w = relaxation``) with air properties
    re-evaluated at the film temperature every step; it converges for both the
    fully developed (``R ~ 1/dT``) and isolated-plate (``R ~ dT^-1/4``) natural
    convection regimes. If it does not reach ``tol_K`` in ``max_iter`` steps a
    bracketing Brent solve on the same residual is used. The iteration record is
    returned in ``solver``.
    """
    if op.delta_T_K is not None:
        return evaluate_at_delta_T(sink, op, op.delta_T_K)
    q = float(op.heat_W)  # type: ignore[arg-type]
    if not (0.0 < relaxation <= 1.0):
        raise ValueError("relaxation must lie in (0, 1]")
    d_t = 20.0
    history: list[dict[str, float]] = []
    converged = False
    for i in range(max_iter):
        perf = evaluate_at_delta_T(sink, op, d_t)
        d_new = q * perf.r_total
        history.append({"iteration": i, "delta_T_in": d_t, "delta_T_out": d_new, "r_total": perf.r_total})
        if abs(d_new - d_t) < tol_K:
            d_t = d_new
            converged = True
            break
        d_t = (1.0 - relaxation) * d_t + relaxation * d_new
        d_t = min(max(d_t, 1e-3), 2000.0)
    method = "damped-fixed-point"
    if not converged:
        def residual(x: float) -> float:
            return x - q * evaluate_at_delta_T(sink, op, x).r_total

        d_t = brentq(residual, 1e-3, 2000.0, xtol=tol_K / 10, maxiter=200)
        method = "fixed-point-fallback-brentq"
        converged = True
    solver = {"mode": "heat-load", "method": method, "converged": converged, "iterations": len(history),
              "tolerance_K": tol_K, "relaxation": relaxation, "history": history}
    return evaluate_at_delta_T(sink, op, d_t, solver)


# =============================================================================
# Multi-objective optimisation
# =============================================================================

VARIABLE_ORDER = ("length", "width", "fin_height", "fin_thickness", "fin_gap", "base_thickness", "volumetric_flow")
OBJECTIVES = ("r_th", "mass", "fan_power_W", "volume_m3", "delta_T_K")


@dataclass
class DesignBounds:
    """Lower/upper bounds (m, m^3/s) of the design variables; ``volumetric_flow`` only for forced cooling."""

    length: tuple[float, float] = (0.05, 0.1)
    width: tuple[float, float] = (0.05, 0.1)
    fin_height: tuple[float, float] = (0.02, 0.1)
    fin_thickness: tuple[float, float] = (0.001, 0.005)
    fin_gap: tuple[float, float] = (0.002, 0.01)
    base_thickness: tuple[float, float] = (0.002, 0.01)
    volumetric_flow: tuple[float, float] | None = None

    def as_list(self) -> list[tuple[float, float]]:
        out = [getattr(self, k) for k in VARIABLE_ORDER[:-1]]
        if self.volumetric_flow is not None:
            out.append(self.volumetric_flow)
        for lo, hi in out:
            if not (0 < lo <= hi):
                raise ValueError(f"invalid bounds {(lo, hi)}")
        return out

    @property
    def names(self) -> tuple[str, ...]:
        return VARIABLE_ORDER if self.volumetric_flow is not None else VARIABLE_ORDER[:-1]


@dataclass
class ParetoSettings:
    """Optimiser controls. Every differential-evolution run is seeded from ``seed`` + run index."""

    objectives: tuple[str, ...] = ("r_th", "mass")
    seed: int = 1
    n_weights: int = 9
    n_epsilon: int = 7
    maxiter: int = 60
    popsize: int = 12
    tol: float = 1e-6
    reference: dict[str, float] | None = None  # optional user scale per objective
    penalty: float = 1e3

    def __post_init__(self) -> None:
        if len(self.objectives) < 2:
            raise ValueError("at least two objectives are needed for a Pareto front")
        unknown = [o for o in self.objectives if o not in OBJECTIVES]
        if unknown:
            raise ValueError(f"unknown objectives {unknown}; choose from {OBJECTIVES}")
        if self.n_weights < 2 or self.n_epsilon < 0 or self.maxiter < 1 or self.popsize < 4:
            raise ValueError("n_weights >= 2, n_epsilon >= 0, maxiter >= 1, popsize >= 4 required")


@dataclass
class ParetoPoint:
    variables: dict[str, float]
    objectives: dict[str, float]
    method: str
    parameter: float | None
    performance: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ParetoFront:
    points: list[ParetoPoint]
    objectives: tuple[str, ...]
    utopia: dict[str, float]
    nadir: dict[str, float]
    scale: dict[str, float]
    settings: dict[str, Any]
    bounds: dict[str, Any]
    material: str
    operating_point: dict[str, Any]
    n_evaluations: int
    runs: list[dict[str, Any]]
    attribution: str = __attribution__

    def to_dict(self) -> dict[str, Any]:
        return {"schema": "heatsink-pareto-front/1.0", "attribution": self.attribution, "material": self.material,
                "objectives": list(self.objectives), "utopia": self.utopia, "nadir": self.nadir, "scale": self.scale,
                "settings": self.settings, "bounds": self.bounds, "operating_point": self.operating_point,
                "n_evaluations": self.n_evaluations, "runs": self.runs, "n_points": len(self.points),
                "points": [p.to_dict() for p in self.points]}

    def to_json(self, **kwargs: Any) -> str:
        kwargs.setdefault("indent", 2)
        return json.dumps(self.to_dict(), default=_json_default, **kwargs)

    def knee_point(self) -> ParetoPoint | None:
        """Point closest (normalised Euclidean) to the utopia point."""
        if not self.points:
            return None
        return min(self.points, key=lambda p: sum(((p.objectives[o] - self.utopia[o]) / self.scale[o]) ** 2
                                                  for o in self.objectives))


def pareto_filter(points: Sequence[ParetoPoint], objectives: Sequence[str]) -> list[ParetoPoint]:
    """Return the non-dominated subset (all objectives minimised), de-duplicated and sorted by the 2nd objective."""
    unique: dict[tuple[float, ...], ParetoPoint] = {}
    for p in points:
        key = tuple(round(p.objectives[o], 12) for o in objectives)
        unique.setdefault(key, p)
    cands = list(unique.values())
    front: list[ParetoPoint] = []
    for p in cands:
        f = [p.objectives[o] for o in objectives]
        dominated = False
        for q in cands:
            if q is p:
                continue
            g = [q.objectives[o] for o in objectives]
            if all(gi <= fi for gi, fi in zip(g, f)) and any(gi < fi for gi, fi in zip(g, f)):
                dominated = True
                break
        if not dominated:
            front.append(p)
    return sorted(front, key=lambda p: p.objectives[objectives[1]])


def _de_kwargs(seed: int) -> dict[str, Any]:
    params = inspect.signature(differential_evolution).parameters
    return {"rng": seed} if "rng" in params else {"seed": seed}


class _Problem:
    """Vector -> sink/operating point -> objectives, with evaluation counting."""

    def __init__(self, op: OperatingPoint, bounds: DesignBounds, material: Material, settings: ParetoSettings):
        self.op, self.bounds, self.material, self.settings = op, bounds, material, settings
        self.names = bounds.names
        self.n_eval = 0
        self.cache: dict[tuple[float, ...], dict[str, float] | None] = {}
        self.last_perf: dict[tuple[float, ...], HeatSinkPerformance] = {}

    def decode(self, x: Sequence[float]) -> dict[str, float]:
        return {n: float(v) for n, v in zip(self.names, x)}

    def evaluate(self, x: Sequence[float]) -> dict[str, float] | None:
        key = tuple(round(float(v), 15) for v in x)
        if key in self.cache:
            return self.cache[key]
        self.n_eval += 1
        v = self.decode(x)
        result: dict[str, float] | None
        try:
            sink = PlateFinHeatSink(v["length"], v["width"], v["fin_height"], v["fin_thickness"], v["fin_gap"],
                                    v["base_thickness"], self.material)
            op = self.op
            if "volumetric_flow" in v:
                op = OperatingPoint(**{**op.__dict__, "volumetric_flow": v["volumetric_flow"], "fan_curve": None})
            perf = evaluate(sink, op)
            result = {"r_th": perf.r_total, "mass": perf.mass_kg, "fan_power_W": perf.fan_electrical_power_W or 0.0,
                      "volume_m3": perf.envelope_volume_m3, "delta_T_K": perf.delta_T_base_K}
            self.last_perf[key] = perf
        except (ValueError, ZeroDivisionError, OverflowError):
            result = None
        self.cache[key] = result
        return result

    def objective_vector(self, x: Sequence[float]) -> list[float] | None:
        r = self.evaluate(x)
        if r is None:
            return None
        return [r[o] for o in self.settings.objectives]


def _run_de(problem: _Problem, func: Callable[[np.ndarray], float], seed: int, constraints: Any = None) -> np.ndarray:
    kwargs: dict[str, Any] = dict(bounds=problem.bounds.as_list(), maxiter=problem.settings.maxiter,
                                  popsize=problem.settings.popsize, tol=problem.settings.tol, polish=False,
                                  init="latinhypercube", **_de_kwargs(seed))
    if constraints is not None:
        kwargs["constraints"] = constraints
    res = differential_evolution(func, **kwargs)
    return np.asarray(res.x, dtype=float)


def optimise_pareto(op: OperatingPoint, bounds: DesignBounds | None = None, material: str | Material = "Al-6063-T5",
                    settings: ParetoSettings | None = None) -> ParetoFront:
    """Pareto front of the heat-sink design problem by weight sweep *and* epsilon-constraint.

    1. Anchor runs minimise each objective alone (seeded differential evolution);
       they give the utopia point and a nadir estimate, and the scale
       ``scale_i = nadir_i - utopia_i`` (or ``settings.reference``) used to
       normalise the objectives -- the original script summed K/W and kg directly.
    2. Weight sweep: ``min sum_i w_i (f_i - utopia_i)/scale_i`` for ``n_weights`` weights.
    3. Epsilon-constraint (two objectives): ``min f_1 s.t. f_2 <= eps`` for
       ``n_epsilon`` levels between the anchors; recovers non-convex parts of the
       front that a weighted sum cannot reach.
    4. All solutions are pooled and filtered to the non-dominated set.

    Every run is deterministic for a given ``settings.seed``.
    """
    bounds = bounds or DesignBounds()
    settings = settings or ParetoSettings()
    mat = get_material(material)
    if op.cooling == "forced" and bounds.volumetric_flow is None and op.volumetric_flow is None and not op.fan_curve:
        raise ValueError("forced cooling needs a fixed volumetric_flow, a fan_curve or volumetric_flow bounds")
    if op.cooling == "natural" and bounds.volumetric_flow is not None:
        raise ValueError("volumetric_flow is not a design variable for natural convection")
    problem = _Problem(op, bounds, mat, settings)
    objs = settings.objectives
    big = 1e6
    runs: list[dict[str, Any]] = []
    collected: list[ParetoPoint] = []

    def single(index: int) -> Callable[[np.ndarray], float]:
        def f(x: np.ndarray) -> float:
            v = problem.objective_vector(x)
            return big if v is None else v[index]
        return f

    def record(x: np.ndarray, method: str, parameter: float | None) -> ParetoPoint | None:
        v = problem.evaluate(x)
        if v is None:
            return None
        key = tuple(round(float(t), 15) for t in x)
        perf = problem.last_perf.get(key)
        pt = ParetoPoint(problem.decode(x), {o: v[o] for o in objs}, method, parameter,
                         perf.to_dict() if perf else {})
        collected.append(pt)
        return pt

    run_id = 0
    anchors: list[ParetoPoint] = []
    for i, o in enumerate(objs):
        x = _run_de(problem, single(i), settings.seed + run_id)
        run_id += 1
        pt = record(x, f"anchor-min-{o}", None)
        if pt is None:
            raise RuntimeError(f"anchor optimisation for {o} found no feasible design within the bounds")
        anchors.append(pt)
        runs.append({"method": "anchor", "objective": o, "seed": settings.seed + run_id - 1, "x": pt.variables})
    utopia = {o: min(a.objectives[o] for a in anchors) for o in objs}
    nadir = {o: max(a.objectives[o] for a in anchors) for o in objs}
    scale = {}
    for o in objs:
        ref = (settings.reference or {}).get(o)
        s = ref if ref else nadir[o] - utopia[o]
        scale[o] = s if s > 1e-12 else max(abs(utopia[o]), 1e-12)

    def weighted(w: Sequence[float]) -> Callable[[np.ndarray], float]:
        def f(x: np.ndarray) -> float:
            v = problem.objective_vector(x)
            if v is None:
                return big
            return sum(wi * (vi - utopia[o]) / scale[o] for wi, vi, o in zip(w, v, objs))
        return f

    if len(objs) == 2:
        weights = [(w, 1.0 - w) for w in np.linspace(0.0, 1.0, settings.n_weights)[1:-1]]
    else:
        rng = np.random.default_rng(settings.seed)
        weights = [tuple(rng.dirichlet(np.ones(len(objs)))) for _ in range(settings.n_weights)]
    for w in weights:
        x = _run_de(problem, weighted(w), settings.seed + run_id)
        run_id += 1
        pt = record(x, "weight-sweep", float(w[0]))
        runs.append({"method": "weight-sweep", "weights": [float(t) for t in w], "seed": settings.seed + run_id - 1,
                     "x": pt.variables if pt else None})

    if len(objs) == 2 and settings.n_epsilon > 0:
        o1, o2 = objs
        lo, hi = utopia[o2], nadir[o2]
        if hi - lo > 1e-12:
            levels = np.linspace(lo, hi, settings.n_epsilon + 2)[1:-1]

            def constrained_value(x: np.ndarray) -> float:
                v = problem.evaluate(x)
                return big if v is None else v[o2]

            for eps in levels:
                con = NonlinearConstraint(constrained_value, -np.inf, float(eps))
                x = _run_de(problem, single(0), settings.seed + run_id, constraints=con)
                run_id += 1
                pt = record(x, "epsilon-constraint", float(eps))
                feasible = bool(pt and pt.objectives[o2] <= eps * (1 + 1e-9))
                runs.append({"method": "epsilon-constraint", "epsilon": float(eps), "constrained_objective": o2,
                             "seed": settings.seed + run_id - 1, "feasible": feasible,
                             "x": pt.variables if pt else None})
                if pt and not feasible:
                    collected.remove(pt)
    front = pareto_filter(collected, objs)
    return ParetoFront(points=front, objectives=objs, utopia=utopia, nadir=nadir, scale=scale,
                       settings=asdict(settings), bounds=asdict(bounds), material=mat.name,
                       operating_point=op.to_dict(), n_evaluations=problem.n_eval, runs=runs)


# =============================================================================
# Legacy model (verbatim formulas of the original HS_main.py, regression use only)
# =============================================================================

LEGACY_CONSTANTS = {"RHO_AIR": 1.095, "BETA": 1 / 300, "G": 9.81, "C_P_AIR": 1007, "MU_AIR": 1.95e-5,
                    "K_F_AIR": 0.0274, "K_ALUMINUM": 196, "RHO_ALUMINUM": 2700, "K_COPPER": 380,
                    "RHO_COPPER": 8960 * 0.97}


def legacy_thermal_resistance(L: float, W: float, H: float, t: float, b: float, t_bp: float, delta_temp: float,
                              k: float = LEGACY_CONSTANTS["K_COPPER"]) -> float:
    """Original ``thermal_resistance`` (HS_main.py lines 19-53 at commit f9ae29a): convective term only.

    Fixed air properties, ``m = sqrt(2 h / (k t))`` (perimeter 2L), no base
    conduction or spreading (``t_bp`` unused), copper by default.
    """
    c = LEGACY_CONSTANTS
    el = (c["RHO_AIR"] ** 2 * c["BETA"] * c["G"] * c["C_P_AIR"] * b ** 4 * delta_temp) / (c["MU_AIR"] * c["K_F_AIR"] * L)
    h = (c["K_F_AIR"] / b) * ((576 / el ** 2) + (2.873 / math.sqrt(el))) ** (-1 / 2)
    m = math.sqrt(2 * h / (k * t))
    eta = math.tanh(m * H) / (m * H)
    n = int((W - t) / (t + b)) + 1
    total_area = (W * L - n * L * t) + n * (2 * H * L) * eta
    return 1 / (h * total_area)


def legacy_heat_sink_mass(W: float, L: float, H: float, t: float, t_bp: float, b: float,
                          rho: float = LEGACY_CONSTANTS["RHO_COPPER"]) -> float:
    """Original ``heat_sink_mass`` (HS_main.py lines 56-58): note the different argument order."""
    n = int((W - t) / (t + b)) + 1
    return rho * (W * L * t_bp + n * L * H * t)


def legacy_fitness(params: Sequence[float]) -> float:
    """Original weighted objective (HS_main.py lines 61-69): ``0.5 Rha + 0.5 m`` with identity 'normalisation'."""
    L, W, H, t, b, t_bp = params
    return 0.5 * legacy_thermal_resistance(L, W, H, t, b, t_bp, 40) + 0.5 * legacy_heat_sink_mass(W, L, H, t, t_bp, b)


__all__ = [
    "FluidProperties", "air_properties", "film_temperature_K", "Material", "MATERIALS", "get_material",
    "PlateFinHeatSink", "fin_parameter", "fin_efficiency", "single_fin_resistance", "elenbaas_number",
    "bar_cohen_rohsenow_nusselt", "natural_convection_h", "bar_cohen_optimum_spacing", "radiation_h",
    "tyc_channel_nusselt", "laminar_duct_nusselt", "laminar_duct_friction_fanning_re", "petukhov_friction_factor",
    "gnielinski_nusselt", "dittus_boelter_nusselt", "developing_laminar_fanning_friction", "channel_pressure_drop",
    "fan_curve_pressure", "fan_operating_point", "fan_power", "spreading_resistance", "HeatSource",
    "heat_sink_resistance", "OperatingPoint", "HeatSinkPerformance", "evaluate_at_delta_T", "evaluate",
    "DesignBounds", "ParetoSettings", "ParetoPoint", "ParetoFront", "pareto_filter", "optimise_pareto",
    "legacy_thermal_resistance", "legacy_heat_sink_mass", "legacy_fitness", "LEGACY_CONSTANTS",
]
