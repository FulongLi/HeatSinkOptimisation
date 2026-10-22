# SPDX-License-Identifier: MIT
"""Regression comparisons against the numbers recorded by earlier AIPE runs, before the ``aipe_thermal`` packaging.

Recorded sources (workspace ``AIPE-Full-Ecosystem-Integration``, 2026-10-21):

* ``runs/wf-dab-10kw-first/tasks/thermal-cooling-design/`` -- AIPE-Workflows thermal task on the selected secondary
  devices (C3M0025065K class: database loss 17.446 W per device at Tj = 125 degC, Rth(j-c) 0.45758 K/W, four devices
  per sink, 40 degC ambient): cooling comparison 63.0 / 64.3 / 79.2 degC (liquid / forced / natural) and the Stage A
  forced-air screening whose knee candidate A-001 has R_th = 0.4405 K/W and m = 0.2783 kg.
* ``runs/dab-l1-nominal-ic/thermal/dab-thermal-case.json`` -- ``pea.thermal.dab_case`` run with the TI fact-sheet
  devices: electrothermal cooling comparison for both bridges, converged in 4 iterations for every option.
* The heat-sink model regression (R = 1.426012 K/W, m = 1.609200 kg for L = W = 0.1 m, H = 0.05 m, t = 2 mm, b = 6 mm,
  t_bp = 5 mm copper, dT = 40 K) stays in ``tests/test_heatsink_model.py``.

The recorded values were produced by the same correlations (``pea.thermal`` + the root ``heatsink_model.py``) and the
migrated package reproduces them to floating-point precision. The Stage A screening is a seeded differential
evolution: it reproduced bit-for-bit in the recording environment (scipy 1.15.3) and is compared at 1e-6 relative so
that a different but equivalent scipy random stream still counts as agreement only if it lands on the same designs.
"""

import json
from pathlib import Path

import pytest

from aipe_thermal import cli
from aipe_thermal import cooling as cl
from aipe_thermal import heatsink as hm
from aipe_thermal import network as net
from aipe_thermal import screening
from aipe_thermal.cases import dab_10kw as dc

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "dab-10kw"

# --- wf-dab-10kw-first (database losses, secondary bridge) -------------------------------------------------------------
WF_T_AMBIENT_C = 40.0
WF_DEVICE = {"power_W": 17.446227058701883, "r_jc": 0.45758000000000004, "r_cs": 0.1, "footprint_area_m2": 0.00032,
             "tj_max_C": 175.0, "tj_design_C": 125.0}
WF_TOTAL_POWER_W = 4 * WF_DEVICE["power_W"]  # 69.78490823480753
# option -> (tj_max_C, t_sink_C, r_sink_to_fluid, r_hottest_total, mass_kg, auxiliary_power_W, pressure_drop_Pa)
WF_SUMMARY = {
    "liquid-cold-plate-Al-100x100-16-minichannels-water-2lpm": (63.040292612666306, 50.399048769504844, 0.07736699676294606,
                                                                0.11911780841361214, 0.27012853953439103, 0.004829297666022693,
                                                                43.46367899420424),
    "forced-air-Al-100x100x46-2-fans": (64.29097301810704, 52.1667468168599, 0.17434638985154288, 0.20868847008746383, 0.5784,
                                        1.4333109797810175, 27.955339236723404),
    "natural-convection-Al-200x200x88": (79.2295186308237, 66.40783427157584, 0.37841755387455045, 0.42275403226392266,
                                         2.5920000000000005, 0.0, None),
}
WF_FAN_CURVE = [(0.0, 50.0), (0.005, 35.0), (0.009, 15.0), (0.011, 0.0)]  # per fan, two fans in parallel
# Stage A forced-air screening (seed 1, 5 weights, 3 epsilon levels, maxiter 25, popsize 10, Al-6063-T5)
WF_STAGE_A_UTOPIA = {"r_th": 0.16961535507482098, "mass": 0.0874456161947565, "fan_power_W": 0.0032306302303881716}
WF_STAGE_A_NADIR = {"r_th": 1.202539081782396, "mass": 1.327960409184956, "fan_power_W": 1.2113278765935378}
WF_STAGE_A_CANDIDATES = {  # candidate -> (r_th K/W, mass kg, fan_power_W), recorded to 9 significant digits
    "A-000": (1.20253908, 0.0874456162, 1.14043613), "A-001": (0.44053454, 0.278254076, 0.293924775),
    "A-002": (0.352796249, 0.29763139, 1.04873689), "A-003": (1.1175474, 0.363467815, 0.00323063023),
    "A-004": (0.430009833, 0.539023753, 0.0206105698), "A-005": (0.405959954, 0.566941343, 0.0249506446),
    "A-006": (0.267566968, 0.900922839, 0.114248933), "A-007": (0.169615355, 1.32796041, 1.21132788),
}
WF_STAGE_A_FLAGGED = ["A-001", "A-007", "A-004"]
WF_KNEE_VARIABLES = {"length": 0.08099197532160482, "width": 0.06462391382662835, "fin_height": 0.04517684721720107,
                     "fin_thickness": 0.002025704984332154, "fin_gap": 0.004084806218350192, "base_thickness": 0.004112578190442571,
                     "volumetric_flow": 0.006398004819063193}
WF_KNEE_OBJECTIVES = {"r_th": 0.4405345400652951, "mass": 0.2782540764293634, "fan_power_W": 0.29392477468682326}

# --- dab-thermal-case.json (TI fact-sheet devices, electrothermal closure) ----------------------------------------------
# bridge -> option kind -> (tj_max_C, t_sink_C, total_power_W at convergence, conduction_W per device at convergence)
CASE_SUMMARY = {
    "primary": {"liquid": (70.57442242886376, 49.721146688056415, 61.01893310063836, 8.50473327515959),
                "forced": (71.07554138924051, 50.65240488214567, 61.08454107726919, 8.521135269317298),
                "natural": (86.39623406681312, 64.67179907602286, 63.09021664570295, 9.022554161425738)},
    "secondary": {"liquid": (66.91122888934885, 51.41794770668751, 82.96238691649486, 13.990596729123716),
                  "forced": (69.540587668888, 54.55894091094716, 83.53675383679378, 14.134188459198445),
                  "natural": (87.3917814819305, 70.83747487280148, 87.43610615151621, 15.109026537879052)},
}


def _wf_devices() -> list[cl.CoolingDevice]:
    return [cl.CoolingDevice.from_dict({"name": f"secondary-{i + 1}", **WF_DEVICE}) for i in range(4)]


def _wf_sources() -> list[hm.HeatSource]:
    return [hm.HeatSource(WF_DEVICE["footprint_area_m2"], WF_DEVICE["power_W"], f"secondary-{i + 1}") for i in range(4)]


# ------------------------------------------------------------------------------------------------ cooling comparison


def test_cooling_comparison_secondary_bridge_reproduces_wf_dab_10kw_first():
    comp = cl.compare_cooling_options(_wf_devices(), dc.cooling_options(), WF_T_AMBIENT_C, electrothermal=False)
    assert comp["ranking_by_tj"] == list(WF_SUMMARY)
    for row in comp["summary_table"]:
        tj, t_sink, r_sink, r_hot, mass, aux, dp = WF_SUMMARY[row["option"]]
        assert row["tj_max_C"] == pytest.approx(tj, abs=1e-9)
        assert row["t_sink_C"] == pytest.approx(t_sink, abs=1e-9)
        assert row["r_sink_to_fluid_K_per_W"] == pytest.approx(r_sink, rel=1e-12)
        assert row["r_hottest_total_K_per_W"] == pytest.approx(r_hot, rel=1e-12)
        assert row["mass_kg"] == pytest.approx(mass, rel=1e-12)
        assert row["auxiliary_power_W"] == pytest.approx(aux, rel=1e-9, abs=1e-15)
        assert row["pressure_drop_Pa"] == (None if dp is None else pytest.approx(dp, rel=1e-12))
        assert row["total_power_W"] == pytest.approx(WF_TOTAL_POWER_W) and row["within_design_limit"]
    # the headline numbers quoted in the task specification: 63.0 / 64.3 / 79.2 degC for liquid / forced / natural
    by_kind = {r["kind"]: r["tj_max_C"] for r in comp["summary_table"]}
    assert (round(by_kind["liquid"], 1), round(by_kind["forced"], 1), round(by_kind["natural"], 1)) == (63.0, 64.3, 79.2)
    assert comp["schema"] == "aipe-thermal-cooling-comparison/0.1" and "aipe_thermal.heatsink" in comp["method"]


def test_cooling_comparison_example_input_runs_through_the_cli(tmp_path):
    out = tmp_path / "cmp.json"
    assert cli.main(["cooling-compare", "--input", str(EXAMPLES / "cooling-compare-secondary.json"), "--output", str(out),
                     "--no-electrothermal"]) == 0
    data = json.loads(out.read_text(encoding="utf-8"))
    for row in data["summary_table"]:
        assert row["tj_max_C"] == pytest.approx(WF_SUMMARY[row["option"]][0], abs=1e-9)


def test_natural_convection_and_heat_spreading_match_recorded_values():
    natural = dc.cooling_options()[0]
    res = cl.evaluate_option(natural, _wf_devices(), WF_T_AMBIENT_C, electrothermal=False)
    d = res.details
    assert res.tj_max_C == pytest.approx(79.2295186308237, abs=1e-9)
    assert d["heat_sink"]["n_fins"] == 20 and d["heat_sink"]["gap_effective"] == pytest.approx(0.008421052631578947, rel=1e-12)
    assert d["h"] == pytest.approx(4.048658839144143, rel=1e-10)
    assert d["elenbaas"] == pytest.approx(46.20066690818446, rel=1e-10)
    assert d["fin_efficiency"] == pytest.approx(0.9602504759312321, rel=1e-10)
    assert d["effective_area_m2"] == pytest.approx(0.6527059076419485, rel=1e-10)
    assert d["delta_T_base_K"] == pytest.approx(29.50187013405834, abs=1e-9)
    assert d["regime"] == "natural-laminar" and d["solver"]["converged"] and d["solver"]["iterations"] == 11
    cell = d["sources"][0]
    assert cell["cell_area"] == pytest.approx(0.01) and cell["r_base_1d"] == pytest.approx(0.0038277511961722485, rel=1e-10)
    assert cell["r_spreading_max"] == pytest.approx(0.17351816236131662, rel=1e-10)  # Lee-Song-Au-Moran spreading per footprint
    assert cell["r_base_total"] == pytest.approx(0.17734591355748885, rel=1e-10)
    assert cell["base_rise_K"] == pytest.approx(3.094017075856867, rel=1e-10)
    # the same numbers straight from the heat-sink model
    sink = hm.PlateFinHeatSink(0.20, 0.20, 0.08, 0.002, 0.008, 0.008, "Al-6063-T5")
    perf = hm.evaluate(sink, hm.OperatingPoint(t_ambient_C=WF_T_AMBIENT_C, cooling="natural", heat_W=WF_TOTAL_POWER_W, sources=_wf_sources()))
    assert perf.r_convective == pytest.approx(0.37841755387455045, rel=1e-12)
    assert perf.r_total == pytest.approx(0.42275403226392266, rel=1e-12)
    assert perf.r_spreading_max == pytest.approx(0.17351816236131662, rel=1e-10)
    assert perf.mass_kg == pytest.approx(2.592, rel=1e-12) and perf.envelope_volume_m3 == pytest.approx(3.52e-3, rel=1e-12)


def test_forced_air_and_fan_curve_intersection_match_recorded_values():
    forced = dc.cooling_options()[1]
    assert forced.fan_curve == WF_FAN_CURVE and forced.n_fans == 2
    res = cl.evaluate_option(forced, _wf_devices(), WF_T_AMBIENT_C, electrothermal=False)
    d = res.details
    assert res.tj_max_C == pytest.approx(64.29097301810704, abs=1e-9)
    assert d["heat_sink"]["n_fins"] == 22 and d["heat_sink"]["gap_effective"] == pytest.approx(0.0031904761904761906, rel=1e-12)
    assert d["heat_sink"]["hydraulic_diameter"] == pytest.approx(0.005909592061742007, rel=1e-12)
    assert res.volumetric_flow == pytest.approx(0.012817864305310907, rel=1e-10)  # fan-curve operating point, two fans
    assert res.pressure_drop_Pa == pytest.approx(27.955339236723404, rel=1e-10)
    assert d["h"] == pytest.approx(34.39405736957006, rel=1e-10)
    assert d["reynolds_dh"] == pytest.approx(1604.0911930369675, rel=1e-10)
    assert d["reynolds_b_star"] == pytest.approx(27.63010706056221, rel=1e-10)
    assert d["channel_velocity"] == pytest.approx(4.7827851885488455, rel=1e-10)
    assert d["fin_efficiency"] == pytest.approx(0.8960171810250009, rel=1e-10)
    assert d["delta_T_base_K"] == pytest.approx(14.563305605899286, abs=1e-9)
    assert d["regime"] == "forced-laminar" and d["solver"]["iterations"] == 17
    br = d["pressure_drop_breakdown"]
    assert br["friction"] == pytest.approx(21.21239183362834, rel=1e-10)
    assert br["contraction"] == pytest.approx(2.9163195441251424, rel=1e-10)
    assert br["expansion"] == pytest.approx(3.8266278589699185, rel=1e-10)
    assert br["fanning_friction_factor"] == pytest.approx(0.024873204765692564, rel=1e-10)
    assert d["fan_air_power_W"] == pytest.approx(0.3583277449452544, rel=1e-10)
    assert res.auxiliary_power_W == pytest.approx(1.4333109797810175, rel=1e-10)  # Q dp / 0.25
    assert d["sources"][0]["r_spreading_max"] == pytest.approx(0.1258850673551671, rel=1e-10)
    assert d["sources"][0]["r_base_total"] == pytest.approx(0.13736832094368384, rel=1e-10)
    # fan-curve intersection: the fan pressure at Q/2 per fan equals the channel pressure drop at the film temperature
    sink = hm.PlateFinHeatSink(0.10, 0.10, 0.04, 0.0015, 0.003, 0.006, "Al-6063-T5")
    props = hm.air_properties(hm.film_temperature_K(WF_T_AMBIENT_C, d["delta_T_base_K"]))
    op = hm.fan_operating_point(sink, WF_FAN_CURVE, props, n_fans_parallel=2)
    assert op["volumetric_flow"] == pytest.approx(res.volumetric_flow, rel=1e-9)
    assert hm.fan_curve_pressure(WF_FAN_CURVE, op["volumetric_flow"] / 2) == pytest.approx(op["pressure_drop"], rel=1e-9)
    assert op["pressure_drop"] == pytest.approx(27.955339236723404, rel=1e-9)
    # hand check of the linear fan-curve segment between (0.005, 35) and (0.009, 15) Pa
    q_per_fan = res.volumetric_flow / 2
    assert 35.0 - 20.0 * (q_per_fan - 0.005) / 0.004 == pytest.approx(res.pressure_drop_Pa, rel=1e-9)


def test_liquid_cold_plate_matches_recorded_values():
    liquid = dc.cooling_options()[2]
    res = cl.evaluate_option(liquid, _wf_devices(), WF_T_AMBIENT_C, electrothermal=False)
    d = res.details
    assert res.tj_max_C == pytest.approx(63.040292612666306, abs=1e-9) and res.fluid_inlet_C == 45.0
    assert d["coolant"]["temperature_C"] == pytest.approx(45.25300270137734, rel=1e-10)
    assert d["coolant"]["density"] == pytest.approx(989.7308139329559, rel=1e-10)
    assert d["mass_flow_kg_s"] == pytest.approx(0.032991027131098534, rel=1e-10)
    assert d["coolant_rise_K"] == pytest.approx(0.5060549640943711, rel=1e-10)
    assert d["reynolds_dh"] == pytest.approx(860.2908773716506, rel=1e-10) and d["regime"] == "laminar"
    assert d["prandtl"] == pytest.approx(3.9266665211100573, rel=1e-10)
    assert d["nusselt"] == pytest.approx(3.8963387318400007, rel=1e-10)
    assert d["h"] == pytest.approx(662.735064292824, rel=1e-10)
    assert d["wall_fin_efficiency"] == pytest.approx(0.9788806066110491, rel=1e-10)
    assert d["r_convective"] == pytest.approx(0.07374117724913634, rel=1e-10)
    assert d["r_caloric"] == pytest.approx(0.003625819513809717, rel=1e-10)
    assert d["darcy_friction_factor"] == pytest.approx(0.11449014626052975, rel=1e-10)
    assert res.pressure_drop_Pa == pytest.approx(43.46367899420424, rel=1e-10)
    assert d["pump_hydraulic_power_W"] == pytest.approx(0.0014487892998068079, rel=1e-10)
    assert res.auxiliary_power_W == pytest.approx(0.004829297666022693, rel=1e-10)
    assert d["sources"][0]["r_spreading_max"] == pytest.approx(0.15934774421031983, rel=1e-10)
    assert res.mass_kg == pytest.approx(0.27012853953439103, rel=1e-10) and res.volume_m3 == pytest.approx(1.2e-4, rel=1e-12)


# ------------------------------------------------------------------------------------------------ shared-sink network


def test_shared_sink_network_reproduces_recorded_junction_temperatures():
    r_sa, r_base_total = 0.17434638985154288, 0.13736832094368384  # forced-air option, recorded
    sink = net.SharedHeatSink(r_sa, WF_T_AMBIENT_C)
    paths = [net.DeviceThermalPath(f"secondary-{i + 1}", WF_DEVICE["power_W"], WF_DEVICE["r_jc"], WF_DEVICE["r_cs"], r_base_total,
                                   WF_DEVICE["tj_max_C"], WF_DEVICE["tj_design_C"]) for i in range(4)]
    result = net.solve_steady_state(paths, sink)
    assert result["t_sink_C"] == pytest.approx(52.1667468168599, abs=1e-9)
    assert result["t_j_max_C"] == pytest.approx(64.29097301810704, abs=1e-9)
    assert result["devices"][0]["margin_to_tj_design_K"] == pytest.approx(125.0 - 64.29097301810704, abs=1e-9)
    assert result["all_within_design_limit"] and result["all_within_absolute_limit"]
    # allowable power and required resistance are consistent with the same network
    p_allow = net.allowable_power(paths[0], sink, 3 * WF_DEVICE["power_W"], 125.0)
    assert net.solve_steady_state(paths, sink, {"secondary-1": p_allow})["t_j_max_C"] == pytest.approx(125.0)
    r_req = net.required_sink_resistance(paths, WF_T_AMBIENT_C, 125.0)
    assert net.solve_steady_state(paths, net.SharedHeatSink(r_req, WF_T_AMBIENT_C))["t_j_max_C"] == pytest.approx(125.0)


def test_ti_reference_network_cross_check_matches_dab_thermal_case(tmp_path):
    est = dc.loss_estimates()
    assert est["primary"]["total_W"] == pytest.approx(17.035978999999998, rel=1e-12)
    assert est["secondary"]["total_W"] == pytest.approx(23.91275, rel=1e-12)
    assert est["converter_total_W"] == pytest.approx(163.794916, rel=1e-12)
    ti = dc.ti_reference_network(est)
    assert ti["primary"]["t_sink_C"] == pytest.approx(57.035979, rel=1e-12)
    assert ti["primary"]["t_j_C"] == pytest.approx(77.47915379999999, rel=1e-12)
    assert ti["primary"]["allowable_power_per_device_W_at_125C"] == pytest.approx(44.98941794117647, rel=1e-12)
    assert ti["secondary"]["t_sink_C"] == pytest.approx(63.91275, rel=1e-12)
    assert ti["secondary"]["t_j_C"] == pytest.approx(77.782145, rel=1e-12)
    assert ti["secondary"]["allowable_power_per_device_W_at_125C"] == pytest.approx(67.63298611111111, rel=1e-12)
    # the same network through the CLI example (steady state plus a Foster transient that ends at the static value)
    out = tmp_path / "net.json"
    assert cli.main(["network", "--input", str(EXAMPLES / "network-ti-reference.json"), "--output", str(out)]) == 0
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["steady_state"]["t_sink_C"] == pytest.approx(57.035979, rel=1e-9)
    assert data["steady_state"]["t_j_max_C"] == pytest.approx(77.47915379999999, rel=1e-9)
    assert data["transient"][-1]["t_j_C"] == pytest.approx(data["steady_state"]["devices"][0]["t_j_C"], rel=1e-6)


# ------------------------------------------------------------------------------------------------ electrothermal loop


def test_electrothermal_closure_reproduces_dab_thermal_case_for_both_bridges():
    for bridge, expected in CASE_SUMMARY.items():
        comp = cl.compare_cooling_options(dc.build_devices(bridge), dc.cooling_options(), dc.TI_AMBIENT_C, electrothermal=True)
        assert comp["ranking_by_tj"] == list(WF_SUMMARY)  # liquid, forced, natural
        for row in comp["summary_table"]:
            tj, t_sink, total, _ = expected[row["kind"]]
            assert row["tj_max_C"] == pytest.approx(tj, abs=1e-9)
            assert row["t_sink_C"] == pytest.approx(t_sink, abs=1e-9)
            assert row["total_power_W"] == pytest.approx(total, rel=1e-12)
        for res in comp["results"]:
            _, _, total, conduction = expected[res["kind"]]
            et_summary = res["electrothermal"]
            assert et_summary["converged"] and et_summary["reason"] == "converged" and et_summary["iterations"] == 4
            assert et_summary["total_loss_W"] == pytest.approx(total, rel=1e-12)
            first = next(iter(et_summary["losses_W"].values()))
            assert first["conduction_W"] == pytest.approx(conduction, rel=1e-12) and first["switching_W"] == 6.75
            assert len(et_summary["history"]) == 4
    # the recorded junction temperatures of the natural option (both bridges) stay below the 125 degC design target
    assert CASE_SUMMARY["primary"]["natural"][0] < 125.0 and CASE_SUMMARY["secondary"]["natural"][0] < 125.0


def test_electrothermal_example_input_runs_through_the_cli(tmp_path):
    out = tmp_path / "et.json"
    assert cli.main(["electrothermal", "--input", str(EXAMPLES / "electrothermal-primary-ti-reference.json"), "--output", str(out)]) == 0
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["electrothermal"]["converged"]
    # two devices on the 0.5 K/W TI sink with R_DS(on)(Tj): the converged Tj is consistent with the network at the converged losses
    p = data["electrothermal"]["losses_W"]["primary-1"]["total_W"]
    tj = data["network_at_converged_losses"]["devices"][0]["t_j_C"]
    assert tj == pytest.approx(40.0 + 2 * p * 0.5 + p * 1.2, rel=1e-9)
    assert p < 17.035979  # below the fixed 125 degC estimate because the junction settles far below 125 degC


# ------------------------------------------------------------------------------------------------ Stage A screening


def _wf_stage_a_task() -> screening.StageATask:
    return screening.StageATask.from_dict(json.loads((EXAMPLES / "stage-a-forced-secondary.json").read_text(encoding="utf-8")))


def test_stage_a_knee_candidate_geometry_evaluates_to_recorded_values():
    v = WF_KNEE_VARIABLES
    sink = hm.PlateFinHeatSink(v["length"], v["width"], v["fin_height"], v["fin_thickness"], v["fin_gap"], v["base_thickness"], "Al-6063-T5")
    assert sink.n_fins == 11 and sink.gap_effective == pytest.approx(0.004234115899897466, rel=1e-12)
    op = hm.OperatingPoint(t_ambient_C=WF_T_AMBIENT_C, cooling="forced", heat_W=WF_TOTAL_POWER_W, sources=_wf_sources(),
                           volumetric_flow=v["volumetric_flow"], fan_efficiency=0.25)
    perf = hm.evaluate(sink, op)
    assert perf.r_total == pytest.approx(WF_KNEE_OBJECTIVES["r_th"], rel=1e-10)
    assert perf.mass_kg == pytest.approx(WF_KNEE_OBJECTIVES["mass"], rel=1e-12)
    assert perf.fan_electrical_power_W == pytest.approx(WF_KNEE_OBJECTIVES["fan_power_W"], rel=1e-10)
    assert perf.delta_T_base_K == pytest.approx(30.74266271247645, abs=1e-9)
    assert perf.h == pytest.approx(31.134716692341737, rel=1e-10)
    assert perf.fin_efficiency == pytest.approx(0.9086620360356972, rel=1e-10)
    assert perf.pressure_drop_Pa == pytest.approx(11.485016931022734, rel=1e-10)
    assert perf.r_convective == pytest.approx(0.40965656127906913, rel=1e-10)
    assert perf.r_spreading_max == pytest.approx(0.10847382702310296, rel=1e-10)
    assert perf.envelope_volume_m3 == pytest.approx(0.0002579817611765843, rel=1e-12)
    assert perf.regime == "forced-laminar" and not perf.warnings
    check = screening._junction_check(_wf_stage_a_task(), perf.to_dict())
    assert check["tj_max_C"] == pytest.approx(80.47032973611077, abs=1e-9)
    assert check["t_sink_C"] == pytest.approx(68.58784553664665, abs=1e-9)
    assert check["feasible"] and check["limit_C"] == 125.0 and check["margin_K"] == pytest.approx(44.52967026388923, abs=1e-9)
    # the headline numbers quoted in the task specification
    assert (round(perf.r_total, 4), round(perf.mass_kg, 4)) == (0.4405, 0.2783)


@pytest.mark.slow
def test_stage_a_seeded_screening_reproduces_wf_dab_10kw_first_candidates():
    task = _wf_stage_a_task()
    assert task.seed == 1 and task.n_weights == 5 and task.n_epsilon == 3 and task.maxiter == 25 and task.popsize == 10
    result = screening.run_stage_a(task)
    assert result["schema"] == "aipe-thermal-stage-a/0.1" and result["stage"] == "A"
    assert result["n_candidates"] == 8 and result["n_feasible"] == 6
    assert result["flagged_candidates"] == WF_STAGE_A_FLAGGED
    for key in ("r_th", "mass", "fan_power_W"):
        assert result["utopia"][key] == pytest.approx(WF_STAGE_A_UTOPIA[key], rel=1e-6)
        assert result["nadir"][key] == pytest.approx(WF_STAGE_A_NADIR[key], rel=1e-6)
    rows = {r["candidate_id"]: r for r in result["pareto_set"]}
    assert list(rows) == list(WF_STAGE_A_CANDIDATES)
    for cid, (r_th, mass, fan) in WF_STAGE_A_CANDIDATES.items():
        o = rows[cid]["objectives"]
        assert (o["r_th"], o["mass"], o["fan_power_W"]) == pytest.approx((r_th, mass, fan), rel=1e-6)
    knee = rows["A-001"]
    assert "knee" in knee["flag_reason"] and knee["fem_verification_recommended"]
    for key, value in WF_KNEE_OBJECTIVES.items():
        assert knee["objectives"][key] == pytest.approx(value, rel=1e-6)
    for key, value in WF_KNEE_VARIABLES.items():
        assert knee["variables"][key] == pytest.approx(value, rel=1e-6)
    assert knee["junction_check"]["tj_max_C"] == pytest.approx(80.47032973611077, rel=1e-6)
    assert knee["distance_to_utopia"] == pytest.approx(0.3877502108940394, rel=1e-6)
    assert rows["A-007"]["flag_reason"] == "lowest thermal resistance" and rows["A-004"]["flag_reason"] == "next closest to utopia"
    assert not rows["A-000"]["junction_check"]["feasible"] and not rows["A-003"]["junction_check"]["feasible"]
