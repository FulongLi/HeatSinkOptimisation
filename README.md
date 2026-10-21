# Heat Sink Optimization

![license - CC BY 4.0](https://img.shields.io/badge/license-CC--BY-green)
![type - hardware](https://img.shields.io/badge/type-hardware-blue)
![category - power electronics](https://img.shields.io/badge/category-power%20electronics-lightgrey)
![status - refactored](https://img.shields.io/badge/status-refactored%201.0-blue)

A Python tool for plate-fin heat sink design: analytical natural- and
forced-convection models with base conduction and spreading, a heat-load driven
temperature solve, and a seeded two-method Pareto optimisation (weight sweep and
epsilon-constraint) of thermal resistance against mass (optionally fan power and
volume). Version 1.0.0 (October 2026) refactors the original research script into
the importable module `heatsink_model.py`; every correction is listed with the
original file:line in [`CHANGELOG.md`](CHANGELOG.md).

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Mathematical Models](#mathematical-models)
- [Optimization Parameters](#optimization-parameters)
- [Corrections in 1.0.0](#corrections-in-100)
- [Author and Attribution](#author-and-attribution)
- [Citation](#citation)
- [License](#license)
- [Project Status](#project-status)

## Overview

The models cover the chain **fin surface -> fin root -> base plate -> device
footprint** for a rectangular plate-fin heat sink in SI units. Given a heat
load and ambient temperature the module solves the base temperature rise with
air properties at the film temperature, or evaluates the sink at a fixed rise.
The optimiser returns a non-dominated front (not repeated single-objective
runs) with a JSON record of every design, objective, correlation validity note
and warning.

![Heat Sink Geometry](heat_sink_geometry.png)

## Features

- **Natural convection**: Elenbaas number and the Bar-Cohen & Rohsenow composite
  correlation for vertical isothermal parallel plates; optimum-spacing helper;
  optional linearised radiation term.
- **Forced convection**: Teertstra-Yovanovich-Culham developing-laminar channel
  model (Gnielinski when the channel flow is turbulent), pressure drop
  (entrance, exit, developing friction), fan operating point from a fan curve
  or a fixed flow, fan power.
- **Conduction**: fin efficiency with the full perimeter `2(L + t)`, 1-D base
  conduction and Lee-Song-Au-Moran spreading resistance from device footprints.
- **Materials**: aluminium and copper alloys with sourced properties.
- **Heat-load solve**: damped fixed-point iteration of `dT = Q R(dT)` with
  properties at the film temperature (Brent fallback), iteration record included.
- **Pareto optimisation**: objectives scaled by anchor-run reference values,
  seeded `scipy.optimize.differential_evolution`, weight sweep **and**
  epsilon-constraint, non-dominated filter, knee point, JSON output.
- **Tests**: pytest suite with analytic checks and a regression test documenting
  the original versus corrected objective values.

## Requirements

- Python 3.9 or higher
- NumPy >= 1.19, SciPy >= 1.5 (differential evolution with constraints)
- Matplotlib >= 3.3 only for `--plot`
- pytest for the test suite

## Getting Started

```bash
git clone https://github.com/FulongLi/HeatSinkOptimisation.git
cd HeatSinkOptimisation
pip install -r requirements.txt
pytest            # runs tests/test_heatsink_model.py
```

## Usage

### Command line

```bash
python HS_main.py                                   # natural convection, dT = 40 K, Al-6063-T5
python HS_main.py --material copper --heat-w 40 --ambient-c 40 --json front.json
python HS_main.py --cooling forced --flow 0.01 --objectives r_th mass fan_power_W --plot
python HS_main.py --quick                           # small settings for a smoke test
```

The script prints the front (one row per non-dominated design), the knee
point and the original un-normalised objective at the knee geometry for
comparison, and writes JSON with `--json`.

### As a module

```python
import heatsink_model as hm

sink = hm.PlateFinHeatSink(length=0.1, width=0.1, fin_height=0.05, fin_thickness=0.002,
                           fin_gap=0.006, base_thickness=0.005, material="Al-6063-T5")
devices = [hm.HeatSource(area=16e-3 * 20e-3, power=10.0, name=f"S{i}") for i in range(4)]
op = hm.OperatingPoint(t_ambient_C=40.0, cooling="natural", heat_W=40.0, sources=devices)
perf = hm.evaluate(sink, op)          # solves the base temperature rise
print(perf.delta_T_base_K, perf.r_total, perf.r_spreading_max, perf.validity)

forced = hm.OperatingPoint(t_ambient_C=40.0, cooling="forced", heat_W=40.0, sources=devices,
                           fan_curve=[(0.0, 50.0), (0.005, 35.0), (0.009, 15.0), (0.011, 0.0)], n_fans=2)
print(hm.evaluate(sink, forced).to_json())

front = hm.optimise_pareto(op, hm.DesignBounds(), "Al-6063-T5", hm.ParetoSettings(seed=1))
print(front.knee_point().objectives, len(front.points))
```

## Project Structure

```
HeatSinkOptimisation/
├── heatsink_model.py                   # models, solver, Pareto optimiser, JSON output
├── HS_main.py                          # thin CLI (argparse, optional plotting)
├── tests/test_heatsink_model.py        # pytest suite
├── CHANGELOG.md                        # corrections with original file:line
├── heat_sink_optimization_summary.html # original mathematical documentation (historic)
├── heat_sink_geometry.png              # reference geometry image
├── requirements.txt, pyproject.toml    # dependencies, pytest configuration
└── README.md
```

## Mathematical Models

Symbols: `L` fin length (flow/vertical direction), `W` base width, `H` fin
height, `t` fin thickness, `b` requested gap, `b_eff = (W - n t)/(n - 1)`
effective gap, `n = floor((W - t)/(t + b)) + 1` fins, `t_bp` base thickness.

1. **Air properties** at the film temperature `T_amb + dT/2`: ideal-gas
   density, Sutherland viscosity and conductivity, quadratic `c_p` fit, `beta = 1/T`.
2. **Natural convection** (vertical isothermal plates, laminar, `b << L`):
   `El = rho^2 g beta c_p b^4 dT/(mu k L)`, `Nu_b = [576/El^2 + 2.873/El^0.5]^-0.5`,
   `h = Nu_b k/b`. Optimum spacing `b_opt = 2.714 L/Ra_L^0.25` (El = 54.3).
3. **Forced convection** (ducted flow, no bypass): `Re_b* = Re_b b/L`,
   `Nu_b = [(Re_b* Pr/2)^-3 + (0.664 Re_b*^0.5 Pr^(1/3) (1 + 3.65/Re_b*^0.5)^0.5)^-3]^(-1/3)`
   referenced to the inlet temperature; Gnielinski on the hydraulic diameter when
   `Re_Dh >= 2300` (plus the caloric term `1/(2 m_dot c_p)`).
4. **Pressure drop**: `dp = (K_c + K_e + 4 f_app L/D_h) rho V^2/2`,
   `K_c = 0.42(1 - sigma^2)`, `K_e = (1 - sigma^2)^2`,
   `f_app Re = sqrt((3.44/L*^0.5)^2 + (f Re)_fd^2)`; fan operating point where the
   fan curve meets `dp(Q)`; fan power `Q dp/eta_fan`.
5. **Fin efficiency**: `m = sqrt(h P/(k A_c))`, `P = 2(L + t)`, `A_c = L t`,
   `eta = tanh(m H_c)/(m H_c)` (adiabatic tip, or `H_c = H + t/2`).
6. **Resistance chain**: `R_conv = 1/(h (A_base,exposed + n eta P H))`; base
   `t_bp/(k A)` plus Lee-Song-Au-Moran spreading per footprint cell
   (power-weighted partition of the base, back-face `h_eff = 1/(R_conv A_base)`);
   `R_total = R_conv + max_i P_i (R_1D,i + R_sp,i)/sum P`.
7. **Mass** `rho (W L t_bp + n L H t)`, envelope volume `W L (t_bp + H)`.
8. **Optimisation**: anchors (each objective alone) give utopia/nadir and the
   scale; weight sweep `min sum w_i (f_i - u_i)/s_i`; epsilon-constraint
   `min f_1 s.t. f_2 <= eps`; pooled solutions filtered to the non-dominated set.

All correlations are analytical screening models with the validity ranges
stated in the docstrings and echoed in each result; nothing here is a CFD or
FEM result.

## Optimization Parameters

| Parameter | Symbol | Default range | Description |
|-----------|--------|---------------|-------------|
| Length | L | 0.05 - 0.1 m | Fin length (flow / vertical direction) |
| Width | W | 0.05 - 0.1 m | Base width across the fins |
| Fin Height | H | 0.02 - 0.1 m | Height of the fins |
| Fin Thickness | t | 0.001 - 0.005 m | Thickness of each fin |
| Fin Spacing | b | 0.002 - 0.01 m | Requested gap between fins |
| Base Plate Thickness | t_bp | 0.002 - 0.01 m | Base thickness (now enters conduction and spreading) |
| Volumetric flow | Q | user bounds | Optional design variable for forced convection |

Objectives: `r_th` (K/W, hottest footprint), `mass` (kg), optionally
`fan_power_W`, `volume_m3`, `delta_T_K`.

## Corrections in 1.0.0

See [`CHANGELOG.md`](CHANGELOG.md) for the full table. In short: the identity
"normalisation" of the weighted objective, the missing base conduction and
spreading, the `2L` fin perimeter, the hard-wired copper, the fixed 40 K rise,
the unseeded optimiser and the mislabelled "Pareto front" of repeated runs were
corrected; forced convection, pressure drop, fan power, materials, JSON output
and tests were added.

## Author and Attribution

Original work: **Dr. Fulong Li** -- "Heat Sink Optimization", a Python-based
optimization tool for plate fin heat sink design, created for heat sink
optimization research in power electronics applications (CC BY 4.0).

Version 1.0.0 derivative work (October 2026): refactor into `heatsink_model.py`,
model corrections, forced-convection and spreading models, Pareto method, tests
and documentation for the AIPE 10 kW DAB workflow. Changes are indicated in
`CHANGELOG.md` as required by the licence.

## Citation

```
Heat Sink Optimization - A Python-based optimization tool for plate fin heat sink design
Created by Dr. Fulong Li (original work, CC BY 4.0); version 1.0.0 refactor 2026, see CHANGELOG.md
```

## License

This project is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). You are free to:
- Share -- copy and redistribute the material in any medium or format
- Adapt -- remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- Attribution -- You must give appropriate credit, provide a link to the license, and indicate if changes were made.

For more details, see the [LICENSE](LICENSE.md) file or visit [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/).

## Project Status

The original script is archived upstream; the 1.0.0 refactor is maintained as
part of the AIPE ecosystem (consumed by `AIPE-Design-Agent/pea/thermal`).
