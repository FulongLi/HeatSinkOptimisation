# SPDX-License-Identifier: CC-BY-4.0
"""Single-material Pareto front from the command line -- the former ``HS_main.py`` interface.

Original research script by Dr. Fulong Li (CC BY 4.0). The physics, the
optimiser and the JSON output live in :mod:`aipe_thermal.heatsink.model`; this
module only parses arguments, runs ``optimise_pareto`` and prints or plots the
result. It is reachable as ``aipe-thermal heatsink-front ...`` and through the
thin repository-root wrapper ``python HS_main.py ...``. See ``CHANGELOG.md``
for every change relative to the original script.

Examples::

    aipe-thermal heatsink-front                              # natural convection, dT = 40 K, Al-6063-T5, seeded
    aipe-thermal heatsink-front --material copper --heat-w 40 --ambient-c 40 --json front.json
    aipe-thermal heatsink-front --cooling forced --flow 0.01 --objectives r_th mass fan_power_W --plot
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import aipe_thermal.heatsink as hm


def add_arguments(p: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Attach the front-end options to ``p`` (shared by ``HS_main.py`` and the ``heatsink-front`` subcommand)."""
    p.add_argument("--material", default="Al-6063-T5", help=f"one of {sorted(hm.MATERIALS)} or aluminium/copper")
    p.add_argument("--cooling", choices=["natural", "forced"], default="natural")
    load = p.add_mutually_exclusive_group()
    load.add_argument("--delta-t", type=float, default=None, help="fixed base-to-ambient rise in K (default 40 K)")
    load.add_argument("--heat-w", type=float, default=None, help="heat load in W; the rise is solved")
    p.add_argument("--ambient-c", type=float, default=40.0, help="ambient temperature in degrees Celsius")
    p.add_argument("--flow", type=float, default=None, help="forced convection: channel volumetric flow in m^3/s")
    p.add_argument("--flow-bounds", type=float, nargs=2, default=None, metavar=("MIN", "MAX"),
                   help="forced convection: make the flow a design variable within these bounds (m^3/s)")
    p.add_argument("--objectives", nargs="+", default=["r_th", "mass"], help=f"subset of {hm.OBJECTIVES}")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--n-weights", type=int, default=9)
    p.add_argument("--n-epsilon", type=int, default=7)
    p.add_argument("--maxiter", type=int, default=60)
    p.add_argument("--popsize", type=int, default=12)
    p.add_argument("--quick", action="store_true", help="small optimiser settings for smoke tests")
    p.add_argument("--json", type=Path, default=None, help="write the Pareto front as JSON to this path")
    p.add_argument("--plot", action="store_true", help="show the Pareto front and the knee-point geometry (matplotlib)")
    p.add_argument("--no-plot", action="store_true", help="explicitly disable plotting (default)")
    return p


def build_parser() -> argparse.ArgumentParser:
    return add_arguments(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))


def run(args: argparse.Namespace) -> int:
    """Execute the front end for parsed ``args``; returns the process exit code."""
    delta_t = args.delta_t
    heat = args.heat_w
    if delta_t is None and heat is None:
        delta_t = 40.0
    op_kwargs = dict(t_ambient_C=args.ambient_c, cooling=args.cooling, heat_W=heat, delta_T_K=delta_t)
    bounds = hm.DesignBounds()
    if args.cooling == "forced":
        if args.flow_bounds:
            bounds = hm.DesignBounds(volumetric_flow=(args.flow_bounds[0], args.flow_bounds[1]))
            op_kwargs["volumetric_flow"] = args.flow or args.flow_bounds[1]
        elif args.flow:
            op_kwargs["volumetric_flow"] = args.flow
        else:
            print("forced cooling needs --flow or --flow-bounds", file=sys.stderr)
            return 2
    op = hm.OperatingPoint(**op_kwargs)
    if args.quick:
        settings = hm.ParetoSettings(objectives=tuple(args.objectives), seed=args.seed, n_weights=3, n_epsilon=2,
                                     maxiter=8, popsize=6)
    else:
        settings = hm.ParetoSettings(objectives=tuple(args.objectives), seed=args.seed, n_weights=args.n_weights,
                                     n_epsilon=args.n_epsilon, maxiter=args.maxiter, popsize=args.popsize)
    front = hm.optimise_pareto(op, bounds, args.material, settings)
    _print_front(front)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(front.to_json(), encoding="utf-8")
        print(f"\nJSON written to {args.json}")
    if args.plot and not args.no_plot:
        _plot(front)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


def _print_front(front: hm.ParetoFront) -> None:
    objs = front.objectives
    print(f"Pareto front ({len(front.points)} points, {front.n_evaluations} evaluations, material {front.material})")
    print("objectives scaled by", {k: f"{v:.4g}" for k, v in front.scale.items()}, "utopia", {k: f"{v:.4g}" for k, v in front.utopia.items()})
    header = " | ".join(f"{o:>12s}" for o in objs) + " | method              | L      W      H      t      b      t_bp"
    print(header)
    for p in front.points:
        v = p.variables
        row = " | ".join(f"{p.objectives[o]:12.5g}" for o in objs)
        geom = f"{v['length']:.4f} {v['width']:.4f} {v['fin_height']:.4f} {v['fin_thickness']:.4f} {v['fin_gap']:.4f} {v['base_thickness']:.4f}"
        extra = f" flow={v['volumetric_flow']:.4g}" if "volumetric_flow" in v else ""
        print(f"{row} | {p.method:19s} | {geom}{extra}")
    knee = front.knee_point()
    if knee:
        print("\nKnee point (closest to utopia):")
        print(json.dumps({"variables": knee.variables, "objectives": knee.objectives}, indent=2))
        legacy = hm.legacy_fitness([knee.variables[k] for k in ("length", "width", "fin_height", "fin_thickness", "fin_gap", "base_thickness")])
        print(f"original un-normalised objective 0.5*R + 0.5*m at the knee geometry (copper, dT 40 K): {legacy:.4f}")


def _plot(front: hm.ParetoFront) -> None:
    import matplotlib.pyplot as plt  # optional dependency (pip install aipe-thermal[plot])

    objs = front.objectives
    fig = plt.figure(figsize=(14, 7))
    ax1 = fig.add_subplot(121)
    ax1.plot([p.objectives[objs[1]] for p in front.points], [p.objectives[objs[0]] for p in front.points], "o-")
    ax1.set_xlabel(objs[1])
    ax1.set_ylabel(objs[0])
    ax1.set_title("Pareto front (weight sweep + epsilon-constraint, non-dominated)")
    ax1.grid(True)
    knee = front.knee_point()
    if knee:
        v = knee.variables
        sink = hm.PlateFinHeatSink(v["length"], v["width"], v["fin_height"], v["fin_thickness"], v["fin_gap"], v["base_thickness"], front.material)
        ax2 = fig.add_subplot(122, projection="3d")
        ax2.bar3d(0, 0, 0, sink.width, sink.length, sink.base_thickness, color="gray", alpha=0.6)
        for i in range(sink.n_fins):
            x = i * (sink.fin_thickness + sink.gap_effective)
            ax2.bar3d(x, 0, sink.base_thickness, sink.fin_thickness, sink.length, sink.fin_height, color="tab:blue", alpha=0.7)
        ax2.set_xlabel("Width (m)")
        ax2.set_ylabel("Length (m)")
        ax2.set_zlabel("Height (m)")
        ax2.set_title("Knee-point geometry")
    plt.tight_layout()
    plt.show()


__all__ = ["add_arguments", "build_parser", "run", "main"]


if __name__ == "__main__":
    raise SystemExit(main())
