# SPDX-License-Identifier: MIT
"""Command line for aipe_thermal (JSON in, JSON out). ``python -m aipe_thermal`` is the same interface.

    aipe-thermal dab-case [--quick] [--output DIR]
    aipe-thermal cooling-compare --input spec.json [--output out.json] [--no-electrothermal]
    aipe-thermal heatsink-pareto --input task.json [--output out.json]
    aipe-thermal network --input spec.json [--output out.json]
    aipe-thermal electrothermal --input spec.json [--output out.json]
    aipe-thermal heatsink-front [--material M] [--cooling natural|forced] [--heat-w W] [--json out.json] ...

Input layouts (all SI; every absolute temperature may be given as ``<name>_C`` (degrees Celsius) or ``<name>_K``
(kelvin), flows are m^3/s, lengths m, areas m^2):

cooling-compare: {"t_ambient_C": 40, "devices": [{"name", "power_W", "r_jc", "r_cs", "footprint_area_m2", "tj_max_C", "tj_design_C"}],
                  "options": [{"kind": "natural"|"forced"|"liquid", ...}]}   (see aipe_thermal.cooling.option_from_dict)
heatsink-pareto: StageATask fields (see aipe_thermal.screening.StageATask.from_dict)
network:         {"sink": {"r_sa", "t_ambient_C"}, "devices": [{"name", "power_W", "r_jc", "r_cs", "r_spread", "tj_max_C", "tj_design_C",
                  "zth_jc": {"type": "foster", "stages": [{"resistance", "time_constant"}]}}], "transient": {"device", "power_steps": [[t, P]], "times": [...]}}
electrothermal:  {"sink": {...}, "tj_initial_C": 25, "devices": [{"name", "i_rms_A", "f_sw_Hz", "r_jc", "r_cs", "r_spread", "tj_max_C", "tj_design_C",
                  "rds_on": {"r_ref_ohm", "t_ref_C", "r_hot_ohm", "t_hot_C", "kind"}, "fixed_switching_W",
                  "e_off": {"e_ref_J", "t_ref_C", "v_ref_V", "i_ref_A", "t_coeff_per_K"}, "v_off_V", "i_off_A"}], "settings": {...}}
heatsink-front:  the former HS_main.py options (single-material Pareto front; see aipe_thermal.heatsink.frontend)

Exit codes: 0 success, 1 invalid input / missing file, 3 electrothermal iteration did not converge.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from aipe_thermal import electrothermal as et
from aipe_thermal import network as net
from aipe_thermal._version import __version__
from aipe_thermal.heatsink import frontend
from aipe_thermal.units import MAX_AMBIENT_C, temperature_C_from_mapping


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _dump(data: dict[str, Any], output: Path | None) -> None:
    text = json.dumps(data, indent=2, default=json_default, allow_nan=False)
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text + "\n", encoding="utf-8")
        print(f"written {output}")
    else:
        print(text)


def json_default(obj: Any) -> Any:
    """``json.dumps`` fallback for dataclasses with ``to_dict`` and numpy scalars/arrays."""
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if hasattr(obj, "item"):
        return obj.item()
    if hasattr(obj, "tolist"):
        return obj.tolist()
    return str(obj)


def cmd_dab_case(args: argparse.Namespace) -> int:
    from aipe_thermal.cases.dab_10kw import run_dab_case

    result = run_dab_case(quick=args.quick, output_dir=args.output)
    print(json.dumps({"loss_estimates": result["loss_estimates"],
                      "cooling_summary": {b: result["cooling_comparison"][b]["summary_table"] for b in result["cooling_comparison"]},
                      "stage_a_flagged": {k: v["flagged_candidates"] for k, v in result["heat_sink_stage_a"].items()},
                      "written_to": result.get("written_to")}, indent=2, default=json_default))
    return 0


def cmd_cooling_compare(args: argparse.Namespace) -> int:
    from aipe_thermal.cooling import CoolingDevice, compare_cooling_options, option_from_dict

    spec = _load(args.input)
    devices = [CoolingDevice.from_dict(d) for d in spec["devices"]]
    options = [option_from_dict(o) for o in spec["options"]]
    t_ambient_C = temperature_C_from_mapping(spec, "t_ambient", default=40.0, max_C=MAX_AMBIENT_C)
    result = compare_cooling_options(devices, options, t_ambient_C, electrothermal=not args.no_electrothermal)
    _dump(result, args.output)
    return 0


def cmd_heatsink_pareto(args: argparse.Namespace) -> int:
    from aipe_thermal.screening import StageATask, run_stage_a

    task = StageATask.from_dict(_load(args.input))
    _dump(run_stage_a(task), args.output)
    return 0


def cmd_network(args: argparse.Namespace) -> int:
    spec = _load(args.input)
    sink = net.SharedHeatSink.from_dict(spec["sink"])
    devices = [net.DeviceThermalPath.from_dict(d) for d in spec["devices"]]
    result: dict[str, Any] = {"steady_state": net.solve_steady_state(devices, sink)}
    if spec.get("transient"):
        tr = spec["transient"]
        device = next(d for d in devices if d.name == tr["device"])
        others = sum(d.power_W for d in devices if d.name != device.name)
        result["transient"] = net.transient_junction_temperature(device, sink, tr["power_steps"], tr["times"], others)
    _dump(result, args.output)
    return 0


def cmd_electrothermal(args: argparse.Namespace) -> int:
    spec = _load(args.input)
    sink = net.SharedHeatSink.from_dict(spec["sink"])
    paths, models = [], []
    for d in spec["devices"]:
        rds = et.RdsOnModel(float(d["rds_on"]["r_ref_ohm"]), float(d["rds_on"]["t_ref_C"]), float(d["rds_on"]["r_hot_ohm"]),
                            float(d["rds_on"]["t_hot_C"]), str(d["rds_on"].get("kind", "linear")))
        e_off = et.SwitchingEnergyModel(**{k: v for k, v in d["e_off"].items()}) if d.get("e_off") else None
        models.append(et.DeviceLossModel.datasheet_default(d["name"], float(d["i_rms_A"]), rds, float(d["f_sw_Hz"]), e_off,
                                                           float(d.get("v_off_V", 0.0)), float(d.get("i_off_A", 0.0)),
                                                           fixed_switching_W=float(d.get("fixed_switching_W", 0.0)),
                                                           extra_W=float(d.get("extra_W", 0.0))))
        paths.append(net.DeviceThermalPath(d["name"], 0.0, float(d["r_jc"]), float(d.get("r_cs", 0.0)), float(d.get("r_spread", 0.0)),
                                           temperature_C_from_mapping(d, "tj_max", required=False),
                                           temperature_C_from_mapping(d, "tj_design", required=False)))
    settings = et.ElectrothermalSettings(**spec.get("settings", {}))
    tj_initial_C = temperature_C_from_mapping(spec, "tj_initial", default=25.0)
    result = et.iterate(models, net.make_network_solver(paths, sink), tj_initial_C, settings)
    final = net.solve_steady_state(paths, sink, {n: v["total_W"] for n, v in result.losses_W.items()})
    _dump({"electrothermal": result.to_dict(), "network_at_converged_losses": final}, args.output)
    return 0 if result.converged else 3


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="aipe-thermal", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"aipe-thermal {__version__}")
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("dab-case", help="run the AIPE 10 kW DAB thermal case study")
    d.add_argument("--quick", action="store_true")
    d.add_argument("--output", type=Path, default=None, help="directory for dab-thermal-case.json")
    d.set_defaults(func=cmd_dab_case)
    for name, func, help_text in (("cooling-compare", cmd_cooling_compare, "compare natural / forced / liquid options"),
                                  ("heatsink-pareto", cmd_heatsink_pareto, "Stage A heat-sink Pareto screening"),
                                  ("network", cmd_network, "steady-state (and transient) shared-sink network"),
                                  ("electrothermal", cmd_electrothermal, "loss(Tj) <-> Tj iteration with datasheet defaults")):
        s = sub.add_parser(name, help=help_text)
        s.add_argument("--input", type=Path, required=True)
        s.add_argument("--output", type=Path, default=None)
        if name == "cooling-compare":
            s.add_argument("--no-electrothermal", action="store_true")
        s.set_defaults(func=func)
    f = sub.add_parser("heatsink-front", help="single-material Pareto front (former HS_main.py options)",
                       description=frontend.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    frontend.add_arguments(f)
    f.set_defaults(func=frontend.run)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ValueError, KeyError, OSError, TypeError, ImportError, StopIteration) as exc:
        parser.exit(1, f"error: {exc}\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
