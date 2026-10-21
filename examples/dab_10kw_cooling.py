#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""AIPE 10 kW DAB: cooling comparison and Stage A heat-sink screening with the ``aipe_thermal`` API alone.

Two device sets are available:

* default -- the secondary-bridge devices of the AIPE-Workflows run ``wf-dab-10kw-first`` (semiconductor-database
  loss 17.446 W per device at Tj = 125 degC, Rth(j-c) 0.45758 K/W, four devices on one sink, 40 degC ambient), read
  from ``examples/dab-10kw/cooling-compare-secondary.json`` and ``stage-a-forced-secondary.json``. Recorded results:
  liquid 63.04 / forced 64.29 / natural 79.23 degC and Stage A knee A-001 (R_th 0.4405 K/W, 0.2783 kg).
* ``--ti-case`` -- the TI fact-sheet devices of :mod:`aipe_thermal.cases.dab_10kw` for both bridges with the
  electrothermal closure (R_DS(on)(Tj) + turn-off placeholder), i.e. the ``aipe-thermal dab-case`` inputs.

    python examples/dab_10kw_cooling.py                      # comparison + recorded Stage A task (about 20 s)
    python examples/dab_10kw_cooling.py --quick              # smaller optimiser settings (a few seconds)
    python examples/dab_10kw_cooling.py --skip-stage-a       # comparison only
    python examples/dab_10kw_cooling.py --ti-case --output out/

Everything printed is an analytical correlation estimate with stated validity; no CFD, FEM or measured value.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from aipe_thermal import cooling as cl
from aipe_thermal import screening
from aipe_thermal.cases import dab_10kw as dc
from aipe_thermal.cli import json_default

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "dab-10kw"


def _print_table(title: str, comparison: dict) -> None:
    print(f"\n{title}")
    print(f"{'option':62s} {'Tj,max degC':>12s} {'T_sink degC':>12s} {'R_sink K/W':>11s} {'R_hot K/W':>10s} {'mass kg':>8s} {'aux W':>8s}")
    for row in comparison["summary_table"]:
        print(f"{row['option']:62s} {row['tj_max_C']:12.2f} {row['t_sink_C']:12.2f} {row['r_sink_to_fluid_K_per_W']:11.4f} "
              f"{row['r_hottest_total_K_per_W']:10.4f} {row['mass_kg']:8.3f} {row['auxiliary_power_W']:8.3f}")
    print("ranking by Tj:", " > ".join(comparison["ranking_by_tj"]))
    for result in comparison["results"]:
        if result["warnings"]:
            print(f"  warnings ({result['option_name']}): " + "; ".join(result["warnings"]))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ti-case", action="store_true", help="TI fact-sheet devices for both bridges with electrothermal closure")
    p.add_argument("--quick", action="store_true", help="small Stage A optimiser settings")
    p.add_argument("--skip-stage-a", action="store_true", help="cooling comparison only")
    p.add_argument("--output", type=Path, default=None, help="directory for the JSON results")
    args = p.parse_args(argv)
    results: dict[str, dict] = {}

    if args.ti_case:
        est = dc.loss_estimates()
        print("Loss estimates at Tj = 125 degC (W per device): primary %.3f, secondary %.3f; converter total %.1f W"
              % (est["primary"]["total_W"], est["secondary"]["total_W"], est["converter_total_W"]))
        for bridge in ("primary", "secondary"):
            comp = cl.compare_cooling_options(dc.build_devices(bridge), dc.cooling_options(), dc.TI_AMBIENT_C, electrothermal=True)
            results[f"cooling_comparison_{bridge}"] = comp
            _print_table(f"{bridge} bridge, electrothermal closure, {dc.TI_AMBIENT_C:.0f} degC ambient", comp)
            for res in comp["results"]:
                et_summary = res["electrothermal"]
                print(f"  {res['option_name']}: converged in {et_summary['iterations']} iterations, total loss {et_summary['total_loss_W']:.2f} W")
        devices = dc.build_devices("secondary", with_loss_models=False)
        task_dict = {"devices": [d.to_dict() for d in devices], "t_ambient_C": dc.TI_AMBIENT_C, "cooling": "forced",
                     "materials": ["Al-6063-T5"], "flow_bounds": [0.003, 0.03], "objectives": ["r_th", "mass", "fan_power_W"],
                     "bounds": {"length": [0.06, 0.20], "width": [0.06, 0.20], "fin_height": [0.02, 0.08], "fin_thickness": [0.001, 0.004],
                                "fin_gap": [0.002, 0.008], "base_thickness": [0.003, 0.010]},
                     "seed": 21, "n_weights": 6, "n_epsilon": 0, "maxiter": 25, "popsize": 10, "tj_limit_C": dc.TJ_DESIGN_C,
                     "label": "secondary-bridge-forced-air-variable-flow (TI devices)"}
    else:
        spec = json.loads((INPUTS / "cooling-compare-secondary.json").read_text(encoding="utf-8"))
        devices = [cl.CoolingDevice.from_dict(d) for d in spec["devices"]]
        options = [cl.option_from_dict(o) for o in spec["options"]]
        comp = cl.compare_cooling_options(devices, options, float(spec["t_ambient_C"]), electrothermal=False)
        results["cooling_comparison_secondary"] = comp
        _print_table(f"secondary bridge, database losses (fixed), {spec['t_ambient_C']:.0f} degC ambient", comp)
        task_dict = json.loads((INPUTS / "stage-a-forced-secondary.json").read_text(encoding="utf-8"))

    if not args.skip_stage_a:
        if args.quick:
            task_dict = {**task_dict, "n_weights": 3, "n_epsilon": 1, "maxiter": 6, "popsize": 6}
        task = screening.StageATask.from_dict(task_dict)
        stage_a = screening.run_stage_a(task)
        results["stage_a"] = stage_a
        print(f"\nStage A {task.cooling} screening ({task.label}): {stage_a['n_feasible']}/{stage_a['n_candidates']} candidates feasible, "
              f"flagged for FEM verification: {stage_a['flagged_candidates']}")
        for row in stage_a["pareto_set"]:
            if row["fem_verification_recommended"]:
                v, o = row["variables"], row["objectives"]
                print(f"  {row['candidate_id']} ({row['flag_reason']}): R_th {o['r_th']:.4f} K/W, mass {o['mass']:.4f} kg"
                      + (f", fan {o['fan_power_W']:.3f} W" if "fan_power_W" in o else "")
                      + f"; L {v['length'] * 1e3:.1f} W {v['width'] * 1e3:.1f} H {v['fin_height'] * 1e3:.1f} t {v['fin_thickness'] * 1e3:.2f} "
                        f"b {v['fin_gap'] * 1e3:.2f} t_bp {v['base_thickness'] * 1e3:.2f} mm"
                      + (f", Q {v['volumetric_flow'] * 1e3:.2f} L/s" if "volumetric_flow" in v else "")
                      + f"; Tj,max {row['junction_check']['tj_max_C']:.1f} degC")
        print("next stage:", stage_a["next_stage"])

    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        for name, data in results.items():
            path = args.output / f"{name}.json"
            path.write_text(json.dumps(data, indent=2, default=json_default) + "\n", encoding="utf-8")
            print("written", path)
    print("\nAll values are analytical correlation estimates (validity notes are recorded per evaluation); not CFD, FEM or measured.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
