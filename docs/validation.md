# Validation of `aipe_thermal` 0.2.0

"Validated" here means: analytic limit checks of every correlation, migrated unit tests, and regression against
numbers recorded by earlier AIPE runs. No CFD, FEM or measured value backs any heat-sink, cold-plate or device
temperature; Stage B (FEM/CFD) is flagged, not performed.

## Test suite (`python -m pytest`, 84 tests, about 45 s)

| File | Tests | Origin | What it checks |
| --- | --- | --- | --- |
| `tests/test_heatsink_model.py` | 28 | 27 from the 1.0.0 refactor (+1: wrapper and front end tested separately) | air-property table values, fin count / effective gap, mass, fin efficiency limits and textbook single-fin resistance, Bar-Cohen asymptotes and optimum spacing, Elenbaas definition, TYC limits, Shah-London values, Gnielinski / Dittus-Boelter, Poiseuille limit of the pressure drop, developing friction, fan operating point, spreading limits (1-D and half-space), source bookkeeping, heat-load solve self-consistency, JSON output, original-vs-corrected regression (1.426012 K/W, 1.609200 kg), footprint regression (51.84 K, 0.0757 K/W spreading), Pareto reproducibility / non-dominance / knee, `HS_main.py` wrapper |
| `tests/test_network.py` | 8 | migrated from `pea/thermal` | steady-state coupling through the shared node, power override, allowable power, required resistance, electrothermal solver callable, Foster limits and superposition, exact Foster -> Cauer with simulation within 0.5 %, transient Zth, dict round trip |
| `tests/test_electrothermal.py` | 7 | migrated | two-point R_DS(on) models, switching-energy scaling, datasheet-default losses, convergence to the closed-form fixed point, runaway detection, shared-sink iteration, relaxation settings |
| `tests/test_cooling.py` | 8 | migrated (sibling-locator test replaced by a bundled-model / deprecation-warning test) | bundled model and `HEATSINK_OPT_DIR` deprecation warning (no root `heatsink_model.py`), coolant table, option validation, three-option comparison with review items, electrothermal closure lowers losses, cold-plate physics details, option dict round trip, DAB loss estimates and TI cross-check |
| `tests/test_screening.py` | 4 | migrated | non-dominated flagged candidates with junction checks, seed reproducibility, forced screening with the flow variable, task validation and round trip |
| `tests/test_cli.py` | 8 | 4 migrated + 4 new | network / electrothermal / cooling-compare / heatsink-pareto / dab-case commands, error exit codes, `python -m aipe_thermal` with kelvin inputs and `--version`, console-script registration, unit-mix rejection, `heatsink-front` |
| `tests/test_units.py` | 6 | new | conversions, Celsius-vs-kelvin plausibility, mapping reader (`_C` or `_K`, not both), flow hints (L/min, CFM), conventions applied by every dataclass |
| `tests/test_regression_recorded.py` | 11 | new | see the table below |
| `tests/test_manifest.py` | 4 | new | `aipe.yaml` identity, capabilities, evidence/documentation files, README version and licence statements |

Run a subset without the 20 s seeded Stage A reproduction with `python -m pytest -m "not slow"`.

## Regression against the recorded AIPE runs

Sources (workspace `AIPE-Full-Ecosystem-Integration`, 2026-10-21): the AIPE-Workflows thermal task output
`runs/wf-dab-10kw-first/tasks/thermal-cooling-design/` (cooling comparison and Stage A screening with the
semiconductor-database devices: 17.446 W per secondary device at Tj = 125 degC, Rth(j-c) 0.45758 K/W, four per sink,
40 degC ambient) and the full `pea.thermal.dab_case` run `runs/dab-l1-nominal-ic/thermal/dab-thermal-case.json`
(TI fact-sheet devices with electrothermal closure). The same inputs are shipped as CLI examples in
`examples/dab-10kw/`.

| Item required by the task | Recorded value (pea.thermal + root heatsink_model.py) | aipe_thermal 0.2.0 | Test |
| --- | --- | --- | --- |
| Natural convection (Al 200 x 200 x 88 mm, 69.78 W) | Tj,max 79.2295186308237 degC; h 4.0487 W/m^2K; El 46.2007; eta 0.96025; dT 29.5019 K; R_conv 0.378418 K/W | identical to 1e-10 relative | `test_natural_convection_and_heat_spreading_match_recorded_values` |
| Forced air (Al 100 x 100 x 46 mm, two fans) | Tj,max 64.29097301810704 degC; Q 0.0128179 m^3/s; dp 27.9553 Pa; h 34.394; Re_Dh 1604.09; Re_b* 27.630; fan 1.4333 W | identical | `test_forced_air_and_fan_curve_intersection_match_recorded_values` |
| Fan-curve intersection | dp_fan(Q/2) = dp_sink(Q) = 27.955339 Pa on the 60 mm curve segment (0.005, 35)-(0.009, 15) Pa | identical (Brent solve re-run at the film temperature) | same test |
| Heat spreading (Lee-Song-Au-Moran per footprint cell) | natural 0.173518 K/W, forced 0.125885 K/W, liquid 0.159348 K/W; model regression 0.075681 K/W | identical | natural / forced / liquid tests, `test_regression_with_device_footprints_and_heat_load` |
| Shared-device network | T_sink 52.1667468168599 degC, Tj 64.29097301810704 degC for r_sa 0.174346, r_spread 0.137368; TI cross-check T_sink 57.035979 / Tj 77.479154 degC | identical | `test_shared_sink_network_reproduces_recorded_junction_temperatures`, `test_ti_reference_network_cross_check_matches_dab_thermal_case` |
| Electrothermal convergence (both bridges, three options) | converged in 4 iterations each; primary 70.574 / 71.076 / 86.396 degC (61.02 / 61.08 / 63.09 W), secondary 66.911 / 69.541 / 87.392 degC (82.96 / 83.54 / 87.44 W) | identical to 1e-9 | `test_electrothermal_closure_reproduces_dab_thermal_case_for_both_bridges` |
| Cooling comparison headline (secondary bridge, database losses) | liquid 63.0 / forced 64.3 / natural 79.2 degC | 63.040292612666 / 64.290973018107 / 79.229518630824 degC | `test_cooling_comparison_secondary_bridge_reproduces_wf_dab_10kw_first` |
| Pareto candidate generation (Stage A forced, seed 1, 5 weights, 3 epsilon, maxiter 25, popsize 10) | 8 candidates, 6 feasible, flagged A-001 / A-007 / A-004; knee A-001: R_th 0.4405345 K/W, m 0.2782541 kg, fan 0.293925 W, Tj 80.470 degC; utopia / nadir as recorded | identical for all 8 candidates (objectives, variables, flags, distances) at 1e-6 relative; bit-for-bit in this environment | `test_stage_a_seeded_screening_reproduces_wf_dab_10kw_first_candidates` (slow), `test_stage_a_knee_candidate_geometry_evaluates_to_recorded_values` (optimiser-free) |
| Heat-sink model regression (L = W = 0.1 m, H = 0.05 m, t = 2 mm, b = 6 mm, t_bp = 5 mm, Cu, dT = 40 K) | R = 1.426012 K/W (1.424733 + 0.001279), m = 1.609200 kg; original script 1.490181 K/W, 1.564416 kg | identical | `test_regression_original_vs_corrected_objective_values` |

Where a migrated number could have differed and why it did not: the modules were copied verbatim apart from imports,
input validation, result-schema labels and the removal of the sibling-checkout loader; the heat-sink correlation
model file was moved, not edited (its only functional addition is the ambient-temperature plausibility check). The
Stage A result depends on scipy's seeded differential evolution; it reproduced bit-for-bit because the recording and
the verification used the same environment (Python 3.10.11, numpy 2.2.6, scipy 1.15.3). Before the migration the frozen
legacy code (`pea.thermal` with `HEATSINK_OPT_DIR` pointing at the root module) was run against the same recorded
files and also reproduced every value (differences 0.0), which established the baseline.

Label changes that are not numerical differences: `schema` `pea-thermal-cooling-comparison/0.1` ->
`aipe-thermal-cooling-comparison/0.1`, `pea-thermal-stage-a/0.1` -> `aipe-thermal-stage-a/0.1`; the `method` strings
name `aipe_thermal.heatsink` instead of "HeatSinkOptimisation heatsink_model".

## Continuous integration

`.github/workflows/ci.yml`: `pip install -e ".[dev]"` and `python -m pytest` on Ubuntu with Python 3.10, 3.11 and
3.12, followed by every documented command (`aipe-thermal network | electrothermal | cooling-compare`,
`python -m aipe_thermal heatsink-front --quick`, `python HS_main.py --quick`, `examples/dab_10kw_cooling.py --quick`),
and a `manifest` job that validates `aipe.yaml` together with the AIPE-Core manifest using the AIPE-Registry
validator (`scripts/validate.py`). Dependencies are pinned by lower bound only (numpy >= 1.24, scipy >= 1.10); a scipy
release that changes the differential-evolution random stream would surface as a failure of the slow Stage A
reproduction test, which is intended.

## Known limitations

* All correlations are analytical screening models with the validity ranges listed in `docs/models.md`; results
  outside those ranges carry warnings but are not refused.
* The DAB case uses stated placeholders (reference-design turn-off loss, same-family R_DS(on) temperature factor,
  assumed cooling geometries, generic fan curve, coolant inlet, fan and pump efficiencies).
* No measured or CFD/FEM reference exists yet for any configuration; Stage B verification is the recommended next step
  for the flagged candidates.
