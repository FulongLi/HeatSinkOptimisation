# Licences of the files in this repository

The repository combines two licences. Nothing here changed licence during the
0.2.0 consolidation: the heat-sink correlation model keeps the CC BY 4.0 of the
original research work, and the modules migrated from AIPE-Electrical-Design
(`pea.thermal`, MIT there) keep MIT. The distribution metadata declares the
SPDX expression `CC-BY-4.0 AND MIT`; every source file carries a
`SPDX-License-Identifier` header.

| Path | Licence | Origin and copyright |
| --- | --- | --- |
| `LICENSE.md` | CC BY 4.0 (licence text) | Creative Commons Attribution 4.0 International |
| `LICENSE-MIT.md` | MIT (licence text) | Copyright (c) 2026 Fulong Li (as in AIPE-Electrical-Design) |
| `src/aipe_thermal/heatsink/model.py` | CC BY 4.0 | Dr. Fulong Li, "Heat Sink Optimization" (`HS_main.py`, archived commit f9ae29a); derivative refactor 1.0.0 (2026-10) with corrections listed in `CHANGELOG.md`; moved from the repository root `heatsink_model.py` in 0.2.0 (history preserved by `git mv`) |
| `src/aipe_thermal/heatsink/frontend.py` | CC BY 4.0 | Command-line front end derived from the original `HS_main.py`; moved from the repository root in 0.2.0 |
| `src/aipe_thermal/heatsink/__init__.py` | CC BY 4.0 | Re-export of the model (trivial glue, kept under the model's licence) |
| `HS_main.py` (repository root) | CC BY 4.0 | Thin wrapper for the documented `python HS_main.py` commands |
| `docs/heat_sink_geometry.png`, `docs/legacy/heat_sink_optimization_summary.html` | CC BY 4.0 | Original project documentation by Dr. Fulong Li (historic) |
| `src/aipe_thermal/network.py` | MIT | Migrated from AIPE-Electrical-Design `pea/thermal/network.py` (commit 614ecd5) |
| `src/aipe_thermal/electrothermal.py` | MIT | Migrated from `pea/thermal/electrothermal.py` (commit 614ecd5) |
| `src/aipe_thermal/cooling.py` | MIT | Migrated from `pea/thermal/cooling.py` (commit 614ecd5); sibling-checkout loader replaced by the bundled model |
| `src/aipe_thermal/screening.py` | MIT | Migrated from `pea/thermal/heatsink_optimisation.py` (commit 614ecd5) |
| `src/aipe_thermal/cases/dab_10kw.py` | MIT | Migrated from `pea/thermal/dab_case.py` (commit 614ecd5) |
| `src/aipe_thermal/cli.py` | MIT | Migrated from `pea/thermal/__main__.py` (commit 614ecd5), extended (`heatsink-front`, kelvin inputs, `--version`) |
| `src/aipe_thermal/__init__.py`, `__main__.py`, `_version.py`, `_compat.py`, `units.py`, `cases/__init__.py` | MIT | New in 0.2.0 |
| `tests/test_heatsink_model.py` | CC BY 4.0 | Tests of the correlation model (from the 1.0.0 refactor) |
| `tests/test_network.py`, `test_electrothermal.py`, `test_cooling.py`, `test_screening.py`, `test_cli.py` | MIT | Migrated from AIPE-Electrical-Design `tests/test_thermal_*.py` |
| `tests/test_units.py`, `tests/test_regression_recorded.py`, `tests/test_manifest.py` | MIT | New in 0.2.0 |
| `examples/` | MIT | New in 0.2.0 (inputs reproduce the AIPE DAB workflow run) |
| `README.md`, `CHANGELOG.md`, `docs/models.md`, `docs/validation.md`, `aipe.yaml`, CI configuration | CC BY 4.0 | Project documentation (the heat-sink documentation is a derivative of the original CC BY 4.0 material; the migrated-module documentation is contributed under the same terms) |

## Attribution notice (CC BY 4.0 material)

Original work: Dr. Fulong Li, "Heat Sink Optimization", a Python-based
optimization tool for plate-fin heat-sink design, licensed under the Creative
Commons Attribution 4.0 International licence
(https://creativecommons.org/licenses/by/4.0/). Changes relative to the
original script are listed in `CHANGELOG.md` (1.0.0 corrections, 0.2.0
packaging), as the licence requires. The attribution string is also available
at run time as `aipe_thermal.heatsink.__attribution__` and is embedded in every
Pareto-front JSON record.

## What the combination means for users

* Using the whole package: comply with both licences -- keep the MIT copyright
  and permission notice, and give attribution with a licence link and change
  notice for the heat-sink model.
* Using only `aipe_thermal.network` / `electrothermal` (standard library, MIT):
  MIT alone applies.
* `aipe_thermal.cooling`, `screening` and the DAB case import the CC BY 4.0
  model, so distributions that include them carry the attribution obligation.
