# Changelog

This repository is licensed under CC BY 4.0. The original work is
`HS_main.py` by Dr. Fulong Li (archived, commit f9ae29a). The 1.0.0 refactor
below is a derivative work: the corrections, the module split and the tests were
made in October 2026 for the AIPE 10 kW DAB workflow. Attribution to the original
author is retained in `README.md`, `heatsink_model.__attribution__` and here, as
required by the licence (https://creativecommons.org/licenses/by/4.0/).

## 1.0.0 -- 2026-10-21

### Structure

- New importable module `heatsink_model.py` (SI units, docstrings with correlation
  validity ranges, JSON output). `HS_main.py` is now a thin CLI with an
  `if __name__ == "__main__"` guard; plotting is opt-in (`--plot`).
- `tests/test_heatsink_model.py` (pytest): analytic checks (fin efficiency limits,
  single-fin resistance, spreading-resistance limits, Poiseuille pressure drop,
  Bar-Cohen/TYC asymptotes, fan operating point) and a regression test that
  documents the original versus corrected objective values.
- The original formulas are kept verbatim as `legacy_thermal_resistance`,
  `legacy_heat_sink_mass` and `legacy_fitness` for regression use only.

### Corrections (file:line refer to the original `HS_main.py` at commit f9ae29a)

| # | Original defect | Location | Correction |
|---|-----------------|----------|------------|
| 1 | "Normalised" weighted objective is the identity (`Rha_norm = Rha`, `mh_norm = mh`, since commit 45004a0), so K/W and kg are summed directly | `HS_main.py:65-69` | Objectives scaled by `(f - utopia)/(nadir - utopia)` from seeded anchor runs (or user reference values) in `optimise_pareto`; weight sweep **and** epsilon-constraint; final non-dominated filter |
| 2 | `thermal_resistance()` contains only the fin-to-air convective term; base-plate conduction and spreading are absent although `t_bp` is an argument | `HS_main.py:39-53` | `heat_sink_resistance()` adds 1-D base conduction `t_bp/(k A)` and Lee-Song-Au-Moran spreading per device footprint (`spreading_resistance`, `HeatSource`) |
| 3 | Fin parameter `m = sqrt(2 h/(k t))` uses perimeter `2L` (fin thickness ignored) | `HS_main.py:28` | `fin_parameter` uses `P = 2 (L + t)`, `A_c = L t`; optional Harper-Brown tip correction |
| 4 | Material fixed to copper through default arguments (`k=K_COPPER`, `rho=RHO_COPPER`) while the README advertises Al/Cu; unexplained density fudge `8960*0.97` | `HS_main.py:16, 27, 47, 56` | `Material` table (`MATERIALS`: Al-6063-T5, Al-6061-T6, Al-pure, Cu-C11000, Cu-pure) selected per sink; sourced property values |
| 5 | Temperature rise fixed to 40 K in three places; no heat load | `HS_main.py:63, 85, 106` | `OperatingPoint(heat_W=...)` solves `dT = Q R(dT)` by damped fixed-point iteration (Brent fallback); `delta_T_K` mode kept |
| 6 | No forced-convection or pressure-drop model | -- | Teertstra-Yovanovich-Culham channel Nusselt (Gnielinski when turbulent), Kays & London entrance/exit losses, Muzychka-Yovanovich developing friction, fan curve operating point, fan power |
| 7 | `differential_evolution` unseeded | `HS_main.py:84` | Every run seeded (`ParetoSettings.seed` + run index); reproducibility tested |
| 8 | 100 repeated single-objective runs of the same scalar objective labelled "Pareto front" | `HS_main.py:82-101` | True Pareto method: anchors, weight sweep, epsilon-constraint, non-dominated filter (`pareto_filter`) |
| 9 | No `__main__` guard, no tests, no machine-readable output | whole file | Guarded CLI, pytest suite, `ParetoFront.to_json()` / `HeatSinkPerformance.to_json()` |
| 10 | Air properties fixed at one (unstated, ~320 K) temperature | `HS_main.py:7-12` | `air_properties(T)` Sutherland/ideal-gas fits evaluated at the film temperature every iteration |
| 11 | Fin count truncation leaves an unused strip of base; requested gap `b` used in the correlation regardless | `HS_main.py:33, 57` | Fins spread evenly: `gap_effective = (W - n t)/(n - 1)` used consistently for convection, pressure drop and areas |
| 12 | Elenbaas/Bar-Cohen correlation applied without stating validity | `HS_main.py:19-24` | Validity notes in every evaluation (`validity`, `warnings`): vertical isothermal plates, laminar, `b << L`, El band, fin efficiency > 0.8 |

### Regression values (documented in `tests/test_heatsink_model.py`)

Geometry L = W = 0.1 m, H = 0.05 m, t = 2 mm, b = 6 mm, t_bp = 5 mm, copper, dT = 40 K:

| Quantity | Original | Corrected (ambient 25 C) |
|----------|----------|--------------------------|
| Thermal resistance | 1.490181 K/W | 1.426012 K/W (1.424733 convective + 0.001279 base) |
| Mass | 1.564416 kg (rho 8691) | 1.609200 kg (rho 8940) |
| Weighted objective | 1.527299 (K/W + kg) | dimensionless, scaled per front |

The corrected resistance is lower mainly because the 13 fins are spread evenly
(gap 6.17 mm instead of 6 mm, higher Elenbaas number); the perimeter fix and
the base term act in the opposite direction. With four 16 x 20 mm footprints at
10 W each the corrected model reports 0.0757 K/W of spreading per footprint cell
and a solved rise of 51.84 K for 40 W.

### Not changed

- The six design variables and their default bounds (Table 7.1 of the thesis
  referenced by the original) are unchanged (`DesignBounds`).
- The Bar-Cohen & Rohsenow composite correlation is unchanged; only its inputs
  (properties, gap) and its documentation were corrected.
- No CFD or FEM results are produced; all models are analytical correlations.
