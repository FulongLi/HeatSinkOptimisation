"""Analytic, limit and regression tests for aipe_thermal.heatsink (run: ``pytest`` from the repo root).

The regression block documents the original HS_main.py objective values next to
the corrected ones for the same geometry; see CHANGELOG.md for the defect list.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

import aipe_thermal.heatsink as hm

# Geometry used for the regression comparison: inside the original bounds, copper
# (the original script's hard-wired default material).
GEOM = dict(L=0.1, W=0.1, H=0.05, t=0.002, b=0.006, t_bp=0.005)


# --------------------------------------------------------------------------- air

def test_air_properties_match_tabulated_values_at_300K():
    p = hm.air_properties(300.0)
    assert p.density == pytest.approx(1.1614, rel=0.015)  # ideal gas vs table 1.1614
    assert p.viscosity == pytest.approx(1.846e-5, rel=0.005)
    assert p.conductivity == pytest.approx(0.0263, rel=0.005)
    assert p.specific_heat == pytest.approx(1007.0, rel=0.003)
    assert p.prandtl == pytest.approx(0.707, rel=0.01)
    assert p.expansion_coefficient == pytest.approx(1 / 300.0)


def test_air_properties_vary_with_film_temperature():
    cold, hot = hm.air_properties(290.0), hm.air_properties(360.0)
    assert hot.viscosity > cold.viscosity and hot.conductivity > cold.conductivity and hot.density < cold.density
    with pytest.raises(ValueError):
        hm.air_properties(50.0)


# --------------------------------------------------------------------------- geometry

def test_fin_count_and_effective_gap():
    s = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="Cu-C11000")
    assert s.n_fins == int((GEOM["W"] - GEOM["t"]) / (GEOM["t"] + GEOM["b"])) + 1 == 13
    assert s.gap_effective == pytest.approx((0.1 - 13 * 0.002) / 12)
    assert s.gap_effective >= s.fin_gap
    assert s.n_fins * s.fin_thickness + s.n_channels * s.gap_effective == pytest.approx(s.width)
    assert s.free_area_ratio == pytest.approx(1 - 13 * 0.002 / 0.1)


def test_mass_and_material_selection():
    cu = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="copper")
    al = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="aluminium")
    expected_volume = 0.1 * 0.1 * 0.005 + 13 * 0.1 * 0.002 * 0.05
    assert cu.mass == pytest.approx(8940.0 * expected_volume)
    assert al.mass / cu.mass == pytest.approx(2700.0 / 8940.0)
    assert cu.material.name == "Cu-C11000" and al.material.name == "Al-6063-T5"
    with pytest.raises(KeyError):
        hm.get_material("unobtainium")
    with pytest.raises(ValueError):
        hm.PlateFinHeatSink(0.1, 0.004, 0.05, 0.002, 0.006, 0.005)  # width cannot hold two fins


# --------------------------------------------------------------------------- fins

def test_fin_efficiency_limits_and_perimeter():
    assert hm.fin_efficiency(0.0, 200.0, 0.05, 0.002, 0.1) == 1.0
    assert hm.fin_efficiency(1e-6, 200.0, 0.05, 0.002, 0.1) == pytest.approx(1.0, abs=1e-6)
    h, k, H, t, L = 5000.0, 10.0, 0.5, 0.001, 0.1
    m = math.sqrt(h * 2 * (L + t) / (k * L * t))
    assert hm.fin_efficiency(h, k, H, t, L) == pytest.approx(1 / (m * H), rel=1e-9)  # tanh -> 1
    # full perimeter 2(L + t) gives a larger fin parameter than the original 2L approximation
    assert hm.fin_parameter(20.0, 200.0, 0.002, 0.1) > math.sqrt(2 * 20.0 / (200.0 * 0.002))
    assert hm.fin_parameter(20.0, 200.0, 0.002, 0.1) == pytest.approx(math.sqrt(20.0 * 2 * 0.102 / (200.0 * 0.1 * 0.002)))
    assert hm.fin_efficiency(20.0, 200.0, 0.04, 0.002, 0.1, tip="convective") < hm.fin_efficiency(20.0, 200.0, 0.04, 0.002, 0.1)


def test_single_fin_resistance_matches_textbook_adiabatic_tip():
    h, k, H, t, L = 20.0, 200.0, 0.04, 0.002, 0.1
    P, Ac = 2 * (L + t), L * t
    m = math.sqrt(h * P / (k * Ac))
    expected = 1.0 / (math.sqrt(h * P * k * Ac) * math.tanh(m * H))
    assert hm.single_fin_resistance(h, k, H, t, L) == pytest.approx(expected, rel=1e-12)


# --------------------------------------------------------------------------- natural convection

def test_bar_cohen_rohsenow_asymptotes():
    assert hm.bar_cohen_rohsenow_nusselt(1e-3) == pytest.approx(1e-3 / 24, rel=1e-6)
    assert hm.bar_cohen_rohsenow_nusselt(1e9) == pytest.approx(0.59 * 1e9 ** 0.25, rel=2e-3)
    with pytest.raises(ValueError):
        hm.bar_cohen_rohsenow_nusselt(0.0)


def test_optimum_spacing_corresponds_to_elenbaas_54():
    props = hm.air_properties(318.15)
    b_opt = hm.bar_cohen_optimum_spacing(0.1, 40.0, props)
    assert hm.elenbaas_number(b_opt, 0.1, 40.0, props) == pytest.approx(2.714 ** 4, rel=1e-9)
    assert 0.004 < b_opt < 0.01


def test_elenbaas_definition():
    p = hm.air_properties(320.0)
    el = hm.elenbaas_number(0.006, 0.1, 40.0, p)
    ra_b = hm.G * p.expansion_coefficient * 40.0 * 0.006 ** 3 / (p.kinematic_viscosity * p.thermal_diffusivity)
    assert el == pytest.approx(ra_b * 0.006 / 0.1, rel=1e-9)


# --------------------------------------------------------------------------- forced convection

def test_tyc_nusselt_limits():
    assert hm.tyc_channel_nusselt(1e-4, 0.7) == pytest.approx(1e-4 * 0.7 / 2, rel=1e-6)
    assert hm.tyc_channel_nusselt(1e8, 0.7) == pytest.approx(0.664 * 1e4 * 0.7 ** (1 / 3), rel=1e-3)


def test_laminar_duct_values_and_turbulent_correlations():
    assert hm.laminar_duct_nusselt(0.0, "T") == pytest.approx(7.541)
    assert hm.laminar_duct_nusselt(0.0, "H") == pytest.approx(8.235)
    assert hm.laminar_duct_nusselt(1.0, "T") == pytest.approx(2.976, abs=0.01)
    assert hm.laminar_duct_nusselt(1.0, "H") == pytest.approx(3.608, abs=0.01)
    assert hm.laminar_duct_friction_fanning_re(0.0) == pytest.approx(24.0)
    assert hm.laminar_duct_friction_fanning_re(1.0) == pytest.approx(14.23, abs=0.01)
    f = hm.petukhov_friction_factor(1e4)
    assert f == pytest.approx((0.79 * math.log(1e4) - 1.64) ** -2)
    nu = hm.gnielinski_nusselt(1e4, 0.7)
    assert nu == pytest.approx((f / 8) * 9000 * 0.7 / (1 + 12.7 * math.sqrt(f / 8) * (0.7 ** (2 / 3) - 1)))
    assert hm.dittus_boelter_nusselt(1e4, 0.7) == pytest.approx(0.023 * 1e4 ** 0.8 * 0.7 ** 0.4)
    assert abs(nu - hm.dittus_boelter_nusselt(1e4, 0.7)) / nu < 0.2  # same order for air at Re = 1e4
    with pytest.raises(ValueError):
        hm.gnielinski_nusselt(500.0, 0.7)


def test_laminar_channel_pressure_drop_reduces_to_poiseuille():
    # near-parallel plates (b/H = 1e-3) and a long channel (L* ~ 100): dp = 12 mu V L / b^2
    sink = hm.PlateFinHeatSink(5.0, 0.1, 1.0, 0.001, 0.001, 0.005, "Al-6063-T5")
    props = hm.air_properties(300.0)
    v = 0.2
    dp = hm.channel_pressure_drop(sink, v * sink.channel_area, props, include_entrance_exit=False)
    poiseuille = 12 * props.viscosity * v * sink.length / sink.gap_effective ** 2
    assert dp["regime"] == "laminar" and dp["reynolds_dh"] < 100
    assert dp["friction"] == pytest.approx(poiseuille, rel=2e-3)
    assert dp["contraction"] == 0.0 and dp["expansion"] == 0.0
    full = hm.channel_pressure_drop(sink, v * sink.channel_area, props)
    sigma, dyn = sink.free_area_ratio, 0.5 * props.density * v ** 2
    assert full["contraction"] == pytest.approx(0.42 * (1 - sigma ** 2) * dyn)
    assert full["expansion"] == pytest.approx((1 - sigma ** 2) ** 2 * dyn)
    assert full["total"] == pytest.approx(full["friction"] + full["contraction"] + full["expansion"])


def test_developing_friction_exceeds_fully_developed():
    f_short = hm.developing_laminar_fanning_friction(500.0, 0.01, 0.005, 0.1)
    f_long = hm.developing_laminar_fanning_friction(500.0, 1000.0, 0.005, 0.1)  # L* = 400: entrance term negligible
    assert f_short > f_long
    assert f_long == pytest.approx(hm.laminar_duct_friction_fanning_re(0.1) / 500.0, rel=1e-3)


def test_fan_operating_point_matches_curve_and_sink():
    sink = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="Al-6063-T5")
    props = hm.air_properties(313.15)
    curve = [(0.0, 60.0), (0.01, 40.0), (0.02, 15.0), (0.03, 0.0)]
    op = hm.fan_operating_point(sink, curve, props, n_fans_parallel=2)
    q = op["volumetric_flow"]
    assert 0 < q < 0.06
    assert hm.fan_curve_pressure(curve, q / 2) == pytest.approx(hm.channel_pressure_drop(sink, q, props)["total"], rel=1e-6)
    one = hm.fan_operating_point(sink, curve, props, n_fans_parallel=1)["volumetric_flow"]
    assert one < q
    power = hm.fan_power(q, op["pressure_drop"], 0.25)
    assert power["electrical_power_W"] == pytest.approx(4 * q * op["pressure_drop"])


# --------------------------------------------------------------------------- spreading

def test_spreading_resistance_limits():
    k = 200.0
    # source covers the whole plate: pure one-dimensional conduction
    a = 0.01
    full = hm.spreading_resistance(a, a, 0.005, k, 50.0)
    assert full["spreading_max"] == 0.0 and full["total_max"] == pytest.approx(0.005 / (k * a))
    # tiny source on a thick, large plate with a strong back-side sink: isoflux half-space maximum 1/(pi k r1)
    r1 = 0.01
    a_s = math.pi * r1 ** 2
    half_space = hm.spreading_resistance(a_s, a_s / 1e-8, 10 * r1, k, 1e6)
    assert half_space["total_max"] == pytest.approx(1 / (math.pi * k * r1), rel=1e-3)
    # more spreading for a smaller source, and the average is below the maximum
    small = hm.spreading_resistance(1e-4, 1e-2, 0.005, k, 50.0)
    large = hm.spreading_resistance(4e-4, 1e-2, 0.005, k, 50.0)
    assert small["spreading_max"] > large["spreading_max"] > 0
    assert small["spreading_avg"] < small["spreading_max"]
    with pytest.raises(ValueError):
        hm.spreading_resistance(1e-4, 1e-2, 0.0, k, 50.0)


def test_heat_sink_resistance_with_sources_adds_spreading_consistently():
    sink = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="Cu-C11000")
    plain = hm.heat_sink_resistance(sink, 30.0)
    assert plain["r_total"] == pytest.approx(plain["r_convective"] + sink.base_thickness / (sink.material.conductivity * sink.base_area))
    sources = [hm.HeatSource(3.2e-4, 10.0, f"S{i}") for i in range(4)]
    with_src = hm.heat_sink_resistance(sink, 30.0, sources)
    assert with_src["r_convective"] == pytest.approx(plain["r_convective"])
    assert with_src["r_total"] > plain["r_total"]
    # equal sources: r_total = r_conv + (1-D + spreading of one cell) / N
    cell = with_src["sources"][0]
    assert with_src["r_total"] == pytest.approx(with_src["r_convective"] + cell["r_base_total"] / 4)
    assert cell["cell_area"] == pytest.approx(sink.base_area / 4)


# --------------------------------------------------------------------------- operating point / solve

def test_heat_load_solve_is_self_consistent_natural_and_forced():
    sink = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="Al-6063-T5")
    nat = hm.evaluate(sink, hm.OperatingPoint(t_ambient_C=40.0, cooling="natural", heat_W=30.0))
    assert nat.solver["converged"] and nat.solver["mode"] == "heat-load"
    assert nat.delta_T_base_K == pytest.approx(30.0 * nat.r_total, abs=1e-3)
    back = hm.evaluate(sink, hm.OperatingPoint(t_ambient_C=40.0, cooling="natural", delta_T_K=nat.delta_T_base_K))
    assert back.heat_W == pytest.approx(30.0, rel=1e-4)
    assert nat.t_film_K == pytest.approx(40.0 + 273.15 + nat.delta_T_base_K / 2)
    forced = hm.evaluate(sink, hm.OperatingPoint(t_ambient_C=40.0, cooling="forced", heat_W=60.0, volumetric_flow=0.01))
    assert forced.regime == "forced-laminar" and forced.pressure_drop_Pa > 0 and forced.fan_electrical_power_W > 0
    assert forced.delta_T_base_K == pytest.approx(60.0 * forced.r_total, abs=1e-3)
    assert forced.r_total < nat.r_total


def test_more_flow_lowers_resistance_and_raises_pressure_drop():
    sink = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="Al-6063-T5")
    low = hm.evaluate(sink, hm.OperatingPoint(cooling="forced", delta_T_K=20.0, volumetric_flow=0.005))
    high = hm.evaluate(sink, hm.OperatingPoint(cooling="forced", delta_T_K=20.0, volumetric_flow=0.02))
    assert high.r_total < low.r_total and high.pressure_drop_Pa > low.pressure_drop_Pa


def test_operating_point_validation():
    with pytest.raises(ValueError):
        hm.OperatingPoint(cooling="natural")  # neither heat nor delta T
    with pytest.raises(ValueError):
        hm.OperatingPoint(cooling="natural", heat_W=10.0, delta_T_K=10.0)
    with pytest.raises(ValueError):
        hm.OperatingPoint(cooling="forced", heat_W=10.0)  # no flow
    with pytest.raises(ValueError):
        hm.OperatingPoint(cooling="natural", heat_W=10.0, sources=[hm.HeatSource(1e-4, 5.0)])  # power mismatch


def test_evaluation_is_json_serialisable_with_validity_notes():
    sink = hm.PlateFinHeatSink(**_sink_kwargs(GEOM))
    perf = hm.evaluate(sink, hm.OperatingPoint(cooling="natural", delta_T_K=40.0, emissivity=0.8))
    data = json.loads(perf.to_json())
    assert data["sink"]["material"] == "Al-6063-T5"
    assert any("Bar-Cohen" in v for v in data["validity"])
    assert data["h_radiation"] > 0


# --------------------------------------------------------------------------- regression: original vs corrected

def test_regression_original_vs_corrected_objective_values():
    """Documented values for L=W=0.1, H=0.05, t=2 mm, b=6 mm, t_bp=5 mm, copper, dT = 40 K.

    Original HS_main.py (commit f9ae29a): R = 1.490181 K/W, m = 1.564416 kg (rho = 8960*0.97),
    weighted objective 0.5 R + 0.5 m = 1.527299 (K/W added to kg).
    Corrected model (ambient 25 C, film 318 K, b_eff = 6.17 mm, perimeter 2(L+t), base 1-D
    conduction, Cu-C11000 rho = 8940): R = 1.426012 K/W (convective 1.424733 + base 0.001279),
    eta = 0.989131, El = 41.40, h = 5.066 W/m^2 K, m = 1.609200 kg.
    """
    g = GEOM
    assert hm.legacy_thermal_resistance(g["L"], g["W"], g["H"], g["t"], g["b"], g["t_bp"], 40) == pytest.approx(1.490181, abs=1e-6)
    assert hm.legacy_heat_sink_mass(g["W"], g["L"], g["H"], g["t"], g["t_bp"], g["b"]) == pytest.approx(1.564416, abs=1e-6)
    assert hm.legacy_fitness([g["L"], g["W"], g["H"], g["t"], g["b"], g["t_bp"]]) == pytest.approx(1.527299, abs=1e-6)
    sink = hm.PlateFinHeatSink(**_sink_kwargs(g), material="Cu-C11000")
    perf = hm.evaluate(sink, hm.OperatingPoint(t_ambient_C=25.0, cooling="natural", delta_T_K=40.0))
    assert perf.r_total == pytest.approx(1.426012, abs=1e-6)
    assert perf.r_convective == pytest.approx(1.424733, abs=1e-6)
    assert perf.r_base_1d == pytest.approx(1.278772e-3, rel=1e-5)
    assert perf.fin_efficiency == pytest.approx(0.989131, abs=1e-6)
    assert perf.elenbaas == pytest.approx(41.4041, abs=1e-3)
    assert perf.mass_kg == pytest.approx(1.609200, abs=1e-6)
    # the original ignored t_bp entirely: its resistance does not change with base thickness
    thick = hm.legacy_thermal_resistance(g["L"], g["W"], g["H"], g["t"], g["b"], 0.05, 40)
    assert thick == hm.legacy_thermal_resistance(g["L"], g["W"], g["H"], g["t"], g["b"], g["t_bp"], 40)
    corrected_thick = hm.evaluate(hm.PlateFinHeatSink(**{**_sink_kwargs(g), "base_thickness": 0.05}, material="Cu-C11000"),
                                  hm.OperatingPoint(t_ambient_C=25.0, cooling="natural", delta_T_K=40.0))
    assert corrected_thick.r_total > perf.r_total and corrected_thick.mass_kg > perf.mass_kg


def test_regression_with_device_footprints_and_heat_load():
    """Four 16 x 20 mm footprints at 10 W each on the same copper sink, 25 C ambient: dT = 51.84 K,
    R = 1.2959 K/W including 0.0757 K/W spreading per footprint cell (values documented 2026-10)."""
    sink = hm.PlateFinHeatSink(**_sink_kwargs(GEOM), material="Cu-C11000")
    src = [hm.HeatSource(16e-3 * 20e-3, 10.0, f"S{i}") for i in range(4)]
    perf = hm.evaluate(sink, hm.OperatingPoint(t_ambient_C=25.0, cooling="natural", heat_W=40.0, sources=src))
    assert perf.delta_T_base_K == pytest.approx(51.8361, abs=2e-3)
    assert perf.r_total == pytest.approx(1.295904, abs=1e-4)
    assert perf.r_spreading_max == pytest.approx(0.075681, abs=1e-5)
    assert perf.hottest_source in {"S0", "S1", "S2", "S3"}


# --------------------------------------------------------------------------- Pareto

def _tiny_settings(seed: int = 3) -> hm.ParetoSettings:
    return hm.ParetoSettings(n_weights=4, n_epsilon=2, maxiter=10, popsize=6, seed=seed)


def test_pareto_front_is_non_dominated_sorted_and_reproducible():
    op = hm.OperatingPoint(t_ambient_C=25.0, cooling="natural", delta_T_K=40.0)
    a = hm.optimise_pareto(op, hm.DesignBounds(), "Cu-C11000", _tiny_settings())
    b = hm.optimise_pareto(op, hm.DesignBounds(), "Cu-C11000", _tiny_settings())
    assert a.points and [p.objectives for p in a.points] == [p.objectives for p in b.points]
    masses = [p.objectives["mass"] for p in a.points]
    rths = [p.objectives["r_th"] for p in a.points]
    assert masses == sorted(masses) and rths == sorted(rths, reverse=True)
    assert set(a.utopia) == {"r_th", "mass"} and all(a.scale[o] > 0 for o in a.scale)
    methods = {p.method for p in a.points}
    assert methods & {"weight-sweep", "epsilon-constraint", "anchor-min-mass", "anchor-min-r_th"}
    assert any(r["method"] == "epsilon-constraint" for r in a.runs) and any(r["method"] == "weight-sweep" for r in a.runs)
    knee = a.knee_point()
    assert knee in a.points
    data = json.loads(a.to_json())
    assert data["schema"] == "heatsink-pareto-front/1.0" and data["n_points"] == len(a.points)
    assert "Fulong Li" in data["attribution"]
    different = hm.optimise_pareto(op, hm.DesignBounds(), "Cu-C11000", _tiny_settings(seed=4))
    assert different.runs[0]["seed"] != a.runs[0]["seed"]


def test_pareto_filter_removes_dominated_points():
    mk = lambda r, m: hm.ParetoPoint({}, {"r_th": r, "mass": m}, "t", None, {})  # noqa: E731
    front = hm.pareto_filter([mk(1.0, 1.0), mk(0.5, 2.0), mk(1.2, 1.1), mk(0.5, 2.0)], ("r_th", "mass"))
    assert [(p.objectives["r_th"], p.objectives["mass"]) for p in front] == [(1.0, 1.0), (0.5, 2.0)]


def test_forced_pareto_with_flow_variable_and_three_objectives():
    op = hm.OperatingPoint(t_ambient_C=40.0, cooling="forced", heat_W=50.0, volumetric_flow=0.01)
    bounds = hm.DesignBounds(volumetric_flow=(0.002, 0.02))
    settings = hm.ParetoSettings(objectives=("r_th", "mass", "fan_power_W"), n_weights=3, n_epsilon=0, maxiter=5, popsize=5)
    front = hm.optimise_pareto(op, bounds, "Al-6063-T5", settings)
    assert front.points and all("volumetric_flow" in p.variables for p in front.points)
    assert all(p.objectives["fan_power_W"] > 0 for p in front.points)
    with pytest.raises(ValueError):
        hm.optimise_pareto(hm.OperatingPoint(cooling="natural", delta_T_K=40.0), bounds, "Al-6063-T5", _tiny_settings())


def test_settings_validation():
    with pytest.raises(ValueError):
        hm.ParetoSettings(objectives=("r_th",))
    with pytest.raises(ValueError):
        hm.ParetoSettings(objectives=("r_th", "cost"))


# --------------------------------------------------------------------------- script entry point

def test_heatsink_frontend_runs_headless_and_writes_json(tmp_path):
    from aipe_thermal.heatsink import frontend

    out = tmp_path / "front.json"
    rc = frontend.main(["--quick", "--json", str(out), "--material", "Al-6063-T5", "--no-plot"])
    assert rc == 0
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["n_points"] >= 1 and data["material"] == "Al-6063-T5"
    assert frontend.main(["--cooling", "forced", "--quick", "--no-plot"]) == 2  # forced convection needs --flow or --flow-bounds


def test_hs_main_wrapper_still_runs_the_documented_command(tmp_path):
    out = tmp_path / "front.json"
    script = Path(__file__).resolve().parents[1] / "HS_main.py"
    proc = subprocess.run([sys.executable, str(script), "--quick", "--json", str(out), "--no-plot"],
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stderr
    assert "Pareto front" in proc.stdout
    assert json.loads(out.read_text(encoding="utf-8"))["n_points"] >= 1


def _sink_kwargs(g: dict[str, float]) -> dict[str, float]:
    return dict(length=g["L"], width=g["W"], fin_height=g["H"], fin_thickness=g["t"], fin_gap=g["b"], base_thickness=g["t_bp"])
