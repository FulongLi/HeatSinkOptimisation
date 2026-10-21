"""Tests for the ``aipe-thermal`` / ``python -m aipe_thermal`` command line."""

import json
import subprocess
import sys

import pytest

from aipe_thermal import cli


def _write(tmp_path, name, data):
    p = tmp_path / name
    p.write_text(json.dumps(data), encoding="utf-8")
    return p


def test_cli_network_steady_and_transient(tmp_path):
    spec = {"sink": {"r_sa": 0.5, "t_ambient_C": 40.0},
            "devices": [{"name": "S1", "power_W": 17.0, "r_jc": 1.1, "r_cs": 0.1, "tj_design_C": 125,
                         "zth_jc": {"type": "foster", "stages": [{"resistance": 0.6, "time_constant": 0.01}, {"resistance": 0.5, "time_constant": 1.0}]}},
                        {"name": "S2", "power_W": 17.0, "r_jc": 1.1, "r_cs": 0.1, "tj_design_C": 125}],
            "transient": {"device": "S1", "power_steps": [[0.0, 17.0]], "times": [0.0, 0.01, 100.0]}}
    out = tmp_path / "net.json"
    assert cli.main(["network", "--input", str(_write(tmp_path, "net-in.json", spec)), "--output", str(out)]) == 0
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["steady_state"]["t_sink_C"] == pytest.approx(57.0)
    assert data["transient"][-1]["t_j_C"] == pytest.approx(data["steady_state"]["devices"][0]["t_j_C"], rel=1e-6)


def test_cli_electrothermal_with_datasheet_defaults(tmp_path):
    spec = {"sink": {"r_sa": 0.5, "t_ambient_C": 40.0}, "tj_initial_C": 50.0,
            "devices": [{"name": "S1", "i_rms_A": 9.67, "f_sw_Hz": 100e3, "r_jc": 1.1, "r_cs": 0.1, "tj_design_C": 125,
                         "rds_on": {"r_ref_ohm": 0.075, "t_ref_C": 25, "r_hot_ohm": 0.1275, "t_hot_C": 175}, "fixed_switching_W": 6.75},
                        {"name": "S2", "i_rms_A": 9.67, "f_sw_Hz": 100e3, "r_jc": 1.1, "r_cs": 0.1, "tj_design_C": 125,
                         "rds_on": {"r_ref_ohm": 0.075, "t_ref_C": 25, "r_hot_ohm": 0.1275, "t_hot_C": 175},
                         "e_off": {"e_ref_J": 60e-6, "t_ref_C": 25, "v_ref_V": 800, "i_ref_A": 20, "t_coeff_per_K": 0.001}, "v_off_V": 800, "i_off_A": 14.26}],
            "settings": {"tol_K": 0.001}}
    out = tmp_path / "et.json"
    assert cli.main(["electrothermal", "--input", str(_write(tmp_path, "et-in.json", spec)), "--output", str(out)]) == 0
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["electrothermal"]["converged"] and data["electrothermal"]["history"]
    assert data["network_at_converged_losses"]["t_j_max_C"] > 40.0


def test_cli_errors_exit_nonzero(tmp_path):
    bad = _write(tmp_path, "bad.json", {"sink": {"r_sa": 0.5}})  # missing t_ambient_C
    with pytest.raises(SystemExit) as exc:
        cli.main(["network", "--input", str(bad)])
    assert exc.value.code == 1


def test_cli_cooling_compare_pareto_and_dab_case(tmp_path):
    from aipe_thermal import cooling as cl
    from aipe_thermal.cases import dab_10kw as dc

    devices = [d.to_dict() for d in dc.build_devices("primary", with_loss_models=False)]
    for d in devices:
        d.pop("loss_model")
    spec = {"t_ambient_C": 40.0, "devices": devices, "options": [cl.option_to_dict(o) for o in dc.cooling_options()]}
    out = tmp_path / "cmp.json"
    assert cli.main(["cooling-compare", "--input", str(_write(tmp_path, "cmp-in.json", spec)), "--output", str(out), "--no-electrothermal"]) == 0
    cmp_data = json.loads(out.read_text(encoding="utf-8"))
    assert len(cmp_data["summary_table"]) == 3 and cmp_data["ranking_by_tj"]

    task = {"devices": devices, "t_ambient_C": 40.0, "cooling": "natural", "materials": ["Al-6063-T5"], "tj_limit_C": 125.0,
            "bounds": {"length": [0.1, 0.3], "width": [0.1, 0.3], "fin_height": [0.03, 0.12], "fin_thickness": [0.0015, 0.005],
                       "fin_gap": [0.004, 0.015], "base_thickness": [0.004, 0.012]},
            "n_weights": 3, "n_epsilon": 1, "maxiter": 5, "popsize": 6}
    out2 = tmp_path / "pareto.json"
    assert cli.main(["heatsink-pareto", "--input", str(_write(tmp_path, "task.json", task)), "--output", str(out2)]) == 0
    pareto = json.loads(out2.read_text(encoding="utf-8"))
    assert pareto["stage"] == "A" and pareto["flagged_candidates"]

    assert cli.main(["dab-case", "--quick", "--output", str(tmp_path / "dab")]) == 0
    case = json.loads((tmp_path / "dab" / "dab-thermal-case.json").read_text(encoding="utf-8"))
    assert case["quick_mode"] is True and set(case["cooling_comparison"]) == {"primary", "secondary"}
    assert set(case["heat_sink_stage_a"]) == {"natural", "forced"}
    assert case["loss_estimates"]["primary"]["total_W"] == pytest.approx(dc.loss_estimates()["primary"]["total_W"])
    assert any("proxy" in a for a in case["assumptions"])


def test_python_dash_m_entry_point_accepts_kelvin_inputs_and_reports_version(tmp_path):
    spec = {"sink": {"r_sa": 0.5, "t_ambient_K": 313.15},
            "devices": [{"name": "S1", "power_W": 17.0, "r_jc": 1.1, "r_cs": 0.1, "tj_design_K": 398.15}]}
    inp = _write(tmp_path, "net-k.json", spec)
    out = tmp_path / "net-k.out.json"
    proc = subprocess.run([sys.executable, "-m", "aipe_thermal", "network", "--input", str(inp), "--output", str(out)],
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["steady_state"]["t_ambient_C"] == pytest.approx(40.0)  # kelvin input converted to degrees Celsius
    row = data["steady_state"]["devices"][0]
    assert row["margin_to_tj_design_K"] == pytest.approx(125.0 - row["t_j_C"])
    version = subprocess.run([sys.executable, "-m", "aipe_thermal", "--version"], capture_output=True, text=True, timeout=120)
    assert version.returncode == 0 and "aipe-thermal 0.2.0" in version.stdout + version.stderr


def test_console_script_is_registered():
    from importlib.metadata import entry_points

    eps = entry_points(group="console_scripts")
    assert any(ep.name == "aipe-thermal" and ep.value == "aipe_thermal.cli:main" for ep in eps)


def test_cli_rejects_mixed_temperature_units_and_kelvin_passed_as_celsius(tmp_path):
    devices = [{"name": "S1", "power_W": 1.0, "r_jc": 1.0}]
    mixed = _write(tmp_path, "mixed.json", {"sink": {"r_sa": 0.5, "t_ambient_C": 40.0, "t_ambient_K": 313.15}, "devices": devices})
    with pytest.raises(SystemExit) as exc:
        cli.main(["network", "--input", str(mixed)])
    assert exc.value.code == 1
    kelvin_as_celsius = _write(tmp_path, "kelvin.json", {"sink": {"r_sa": 0.5, "t_ambient_C": 313.15}, "devices": devices})
    with pytest.raises(SystemExit) as exc:
        cli.main(["network", "--input", str(kelvin_as_celsius)])
    assert exc.value.code == 1


def test_cli_heatsink_front_subcommand(tmp_path):
    out = tmp_path / "front.json"
    assert cli.main(["heatsink-front", "--quick", "--json", str(out), "--no-plot", "--seed", "3"]) == 0
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["n_points"] >= 1 and "Fulong Li" in data["attribution"]
