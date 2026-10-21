# Changelog

Two numbering schemes live in this repository, on purpose:

* **Distribution version** (`aipe-thermal`, `aipe_thermal.__version__`): starts at **0.2.0** with the packaged
  release below. "0.1" was the unpackaged state consumed by the first executed AIPE DAB workflow run
  (2026-10-21): the loose `heatsink_model.py` here plus the `pea.thermal` modules inside AIPE-Electrical-Design
  (distribution `pea` 0.1.0, which located this repository through a sibling-directory lookup). 0.2.0 is the first
  version that installs, has a manifest and CI, and bundles both halves; it is a Prototype (correlations tested
  against analytic limits and recorded runs, no CFD/FEM/measurement), so it stays below 1.0.
* **Heat-sink model version** (`aipe_thermal.heatsink.__version__`, formerly `heatsink_model.__version__`): the
  correlation model keeps its own version **1.0.0** from the October 2026 refactor, because the AIPE evidence
  records of the executed workflow reference it (`tool.version` / `model.version` = 1.0.0) and the model physics did
  not change in 0.2.0.

Licensing: the heat-sink correlation model is CC BY 4.0 (original work `HS_main.py` by Dr. Fulong Li, archived commit
f9ae29a); the modules migrated from AIPE-Electrical-Design are MIT. See `LICENSES.md` for the per-file table.
Attribution to the original author is retained in `README.md`, `aipe_thermal.heatsink.__attribution__`, `LICENSES.md`
and here, as required by CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).

## 0.2.0 -- 2026-10-21 (distribution `aipe-thermal`)

### Structure

- New installable package `aipe_thermal` (src layout, `pyproject.toml`, dependencies numpy >= 1.24 and scipy >= 1.10
  only, Python >= 3.10). Console script `aipe-thermal` and `python -m aipe_thermal`.
- `heatsink_model.py` moved (with history) to `src/aipe_thermal/heatsink/model.py`; `aipe_thermal.heatsink`
  re-exports every public name, so `import aipe_thermal.heatsink as hm` replaces the module that
  `pea.thermal._sibling.load_heatsink_model()` returned. `HS_main.py`'s logic moved to
  `aipe_thermal.heatsink.frontend` (also `aipe-thermal heatsink-front ...`).
- Migrated from AIPE-Electrical-Design `pea/thermal` (commit 614ecd5, MIT): `network.py`, `electrothermal.py`,
  `cooling.py`, `heatsink_optimisation.py` -> `screening.py`, `dab_case.py` -> `cases/dab_10kw.py`,
  `__main__.py` -> `cli.py`. Imports now point at the package; `_sibling.py` (sibling-checkout locator) is gone.
- `aipe.yaml` capability manifest (`aipe.thermal-design`, validated with the AIPE-Registry validator), GitHub Actions
  (pytest on Python 3.10 / 3.11 / 3.12, documented commands, manifest validation), `LICENSES.md`, `LICENSE-MIT.md`,
  `docs/models.md`, `docs/validation.md`, `examples/dab_10kw_cooling.py` and CLI inputs under `examples/dab-10kw/`.
- The historic documentation moved to `docs/legacy/heat_sink_optimization_summary.html` and
  `docs/heat_sink_geometry.png`; the tracked `.idea/` IDE files were removed (already in `.gitignore`).

### Behaviour

- Input validation (`aipe_thermal.units`): absolute temperatures are degrees Celsius with the `_C` suffix (or
  `_K` keys in dictionaries / JSON, exactly one of the two); ambient and coolant temperatures above 300 degC and any
  absolute temperature above 600 degC are rejected with the hint that the value looks like kelvin; volumetric flows
  are m^3/s and values above 1 m^3/s are rejected with the hint to convert from L/min or CFM; geometry must be
  positive and finite. Applied to `OperatingPoint`, `SharedHeatSink`, `DeviceThermalPath`, `CoolingDevice`,
  `PlateFinGeometry`, `ForcedAirOption`, `LiquidColdPlateOption`, `StageATask`, `iterate()` and the CLI.
- `HEATSINK_OPT_DIR` and the `../HeatSinkOptimisation` directory are no longer consulted anywhere. Setting the variable
  only emits a `DeprecationWarning` (`aipe_thermal._compat.check_legacy_environment`). There is no repository-root
  `heatsink_model.py` any more: path-based loaders (the legacy `pea.thermal._sibling`) must import
  `aipe_thermal.heatsink` instead.
- Result `schema` identifiers renamed: `pea-thermal-cooling-comparison/0.1` -> `aipe-thermal-cooling-comparison/0.1`,
  `pea-thermal-stage-a/0.1` -> `aipe-thermal-stage-a/0.1` (`aipe-dab-thermal-case/0.1` and
  `heatsink-pareto-front/1.0` unchanged); `method` strings name `aipe_thermal.heatsink`. All numerical outputs are
  unchanged (see Regression below).
- CLI: new `heatsink-front` subcommand (the former `HS_main.py` options), `--version`, kelvin inputs; exit codes 0 /
  1 / 3 as before.
- `aipe_thermal.cases.dab_10kw.json_default` is public (`_json_default` kept as an alias for existing callers).

### Compatibility

- `HS_main.py` at the repository root is a thin wrapper around `aipe_thermal.heatsink.frontend`; every documented
  `python HS_main.py ...` command still runs (requires `pip install -e .`). It may be removed in 0.3.
- `requirements.txt` now just installs the package (`-e .[dev]`) plus the optional matplotlib.

### Regression against the recorded AIPE runs (all reproduced to floating-point precision)

| Quantity | Recorded (before packaging) | 0.2.0 |
| --- | --- | --- |
| Heat-sink model, L = W = 0.1 m, H = 0.05 m, t = 2 mm, b = 6 mm, t_bp = 5 mm, Cu, dT = 40 K | R = 1.426012 K/W, m = 1.609200 kg | identical |
| Cooling comparison, secondary bridge, database losses (wf-dab-10kw-first) | liquid 63.040 / forced 64.291 / natural 79.230 degC | identical (63.040292612666 / 64.290973018107 / 79.229518630824) |
| Stage A forced-air screening, secondary bridge (seed 1) | 8 candidates, 6 feasible, knee A-001: 0.4405 K/W, 0.2783 kg, 0.294 W | identical (all 8 candidates to 1e-9 relative) |
| DAB case electrothermal closure (dab-thermal-case.json) | primary 70.57 / 71.08 / 86.40 degC, secondary 66.91 / 69.54 / 87.39 degC, 4 iterations each | identical |

No migrated number differs: the modules were copied verbatim apart from imports, validation and labels, and the
seeded differential evolution reproduces bit-for-bit with the same scipy generator (1.15.3 in the recording and in the
verification environment).

### Tests

84 tests (`python -m pytest`, about 45 s): the 27 heat-sink model tests (the `HS_main` test now runs the wrapper as a
subprocess and the front end directly, +1), the 31 migrated `pea.thermal` tests (the sibling-locator test is replaced
by a bundled-model / deprecation-warning test; four CLI tests added), and new tests for units/validation (6),
recorded-run regression (11) and the manifest (4). See `docs/validation.md`.

## 1.0.0 -- 2026-10-21 (heat-sink correlation model, `heatsink_model.py`)

The original work is `HS_main.py` by Dr. Fulong Li (archived, commit f9ae29a). The 1.0.0 refactor below is a
derivative work: the corrections, the module split and the tests were made in October 2026 for the AIPE 10 kW DAB
workflow.

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
