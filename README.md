# AIPE — Thermal Design

`aipe_thermal` evaluates the cooling of power semiconductors with analytical models in SI units: a corrected plate-fin
heat-sink correlation model, shared-sink thermal networks with Foster/Cauer transients, an electrothermal
loss-temperature loop, a natural / forced-air / liquid cold-plate comparison and Stage A heat-sink screening. It is
the thermal tool of the AIPE ecosystem, used by AIPE-Workflows for the 10 kW DAB converter, replacing the former
`pea.thermal` modules of AIPE-Electrical-Design and the loose `heatsink_model.py` of this repository. Every result is
a correlation estimate with recorded validity notes, never a CFD, FEM or measured value.

## Features

- Plate-fin heat-sink model: Bar-Cohen-Rohsenow natural convection, Teertstra-Yovanovich-Culham / Gnielinski channel
  flow, pressure drop with fan-curve operating point, fin efficiency, Lee-Song-Au-Moran spreading, heat-load solve.
- Seeded Pareto optimisation (weight sweep and epsilon-constraint) of resistance, mass, fan power and volume.
- Shared-sink networks with junction margins, exact Foster-to-Cauer conversion and transient junction temperatures.
- Electrothermal iteration with pluggable loss models and runaway detection.
- Cooling-concept comparison (quantitative outputs, qualitative review items) and Stage A screening with junction
  checks and FEM-verification flags.
- Explicit units (degrees Celsius vs kelvin, flows in m^3/s) validated on every input; JSON command line.

## Installation

```bash
git clone https://github.com/FulongLi/AIPE-Thermal-Design.git
cd AIPE-Thermal-Design
pip install -e ".[dev]"      # Python >= 3.10; numpy >= 1.24, scipy >= 1.10, pytest, pyyaml
python -m pytest             # 84 tests, about 45 s
```

## Usage

```bash
aipe-thermal cooling-compare --input examples/dab-10kw/cooling-compare-secondary.json --no-electrothermal
aipe-thermal heatsink-pareto --input examples/dab-10kw/stage-a-forced-secondary.json --output stage-a.json
aipe-thermal network --input examples/dab-10kw/network-ti-reference.json
aipe-thermal electrothermal --input examples/dab-10kw/electrothermal-primary-ti-reference.json
aipe-thermal dab-case --quick --output out/             # the AIPE 10 kW DAB case
python examples/dab_10kw_cooling.py --quick             # API example: comparison + Stage A screening
python HS_main.py --quick --json front.json --no-plot   # former research script, now a thin wrapper
```

```python
import aipe_thermal.heatsink as hm
from aipe_thermal import cooling as cl
from aipe_thermal.cases import dab_10kw as dc

sink = hm.PlateFinHeatSink(0.1, 0.1, 0.05, 0.002, 0.006, 0.005, "Cu-C11000")
perf = hm.evaluate(sink, hm.OperatingPoint(t_ambient_C=25.0, cooling="natural", delta_T_K=40.0))
print(round(perf.r_total, 6), round(perf.mass_kg, 6))        # 1.426012 K/W, 1.6092 kg

devices = [cl.CoolingDevice(f"S{i}", 17.446227, 0.45758, 0.1, 3.2e-4, tj_max_C=175, tj_design_C=125) for i in range(4)]
table = cl.compare_cooling_options(devices, dc.cooling_options(), 40.0, electrothermal=False)["summary_table"]
print([(row["kind"], round(row["tj_max_C"], 1)) for row in table])   # liquid 63.0, forced 64.3, natural 79.2 degC
```

## Project Information

- Version: 0.2.0 (heat-sink correlation model 1.0.0, see CHANGELOG.md)
- Status: Prototype
- Developer: Fulong Li
- Organisation: Spirit Connect / AIPE Labs
- License: CC BY 4.0 for the heat-sink correlation model (`src/aipe_thermal/heatsink/`, `HS_main.py`; original
  research by Dr. Fulong Li, attribution required) and MIT for the migrated network, electrothermal, cooling,
  screening, case and CLI modules; per-file table in [LICENSES.md](LICENSES.md), texts in [LICENSE.md](LICENSE.md)
  and [LICENSE-MIT.md](LICENSE-MIT.md).

## Documentation

- [docs/models.md](docs/models.md) — correlations, validity ranges, network / electrothermal / cold-plate / screening models
- [docs/validation.md](docs/validation.md) — test suite, regression against the recorded AIPE DAB runs, CI
- [CHANGELOG.md](CHANGELOG.md) — 0.2.0 packaging; 1.0.0 corrections of the original script
- [LICENSES.md](LICENSES.md), [aipe.yaml](aipe.yaml) (capability manifest)
- [docs/legacy/heat_sink_optimization_summary.html](docs/legacy/heat_sink_optimization_summary.html) — original mathematical summary (historic)
