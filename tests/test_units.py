# SPDX-License-Identifier: MIT
"""Unit conventions: degrees Celsius versus kelvin made explicit, flows in m^3/s, positive geometry."""

import pytest

from aipe_thermal import cooling as cl
from aipe_thermal import electrothermal as et
from aipe_thermal import heatsink as hm
from aipe_thermal import network as net
from aipe_thermal import screening, units


def test_conversions():
    assert units.kelvin_from_celsius(40.0) == pytest.approx(313.15)
    assert units.celsius_from_kelvin(313.15) == pytest.approx(40.0)
    assert units.m3_s_from_lpm(2.0) == pytest.approx(2.0 / 60000.0)
    assert units.m3_s_from_cfm(1.0) == pytest.approx(4.719e-4, rel=1e-3)


def test_temperature_validation_distinguishes_celsius_from_kelvin():
    assert units.validate_temperature_C(125.0, "tj_design_C") == 125.0
    assert units.validate_ambient_C(40.0, "t_ambient_C") == 40.0
    with pytest.raises(ValueError, match="kelvin"):
        units.validate_ambient_C(313.15, "t_ambient_C")  # 40 degC written in kelvin
    with pytest.raises(ValueError, match="kelvin"):
        units.validate_temperature_C(773.15, "tj_max_C")
    with pytest.raises(ValueError, match="absolute zero"):
        units.validate_temperature_C(-300.0, "t_C")
    with pytest.raises(ValueError, match="finite"):
        units.validate_temperature_C(float("nan"), "t_C")
    with pytest.raises(ValueError, match="real number"):
        units.validate_temperature_C("warm", "t_C")
    assert units.validate_temperature_K(313.15, "t_K") == 313.15
    with pytest.raises(ValueError, match="positive"):
        units.validate_temperature_K(-1.0, "t_K")


def test_mapping_reader_accepts_exactly_one_unit():
    assert units.temperature_C_from_mapping({"t_ambient_C": 40.0}, "t_ambient") == 40.0
    assert units.temperature_C_from_mapping({"t_ambient_K": 313.15}, "t_ambient") == pytest.approx(40.0)
    assert units.temperature_C_from_mapping({}, "t_ambient", default=25.0) == 25.0
    assert units.temperature_C_from_mapping({"tj_max_C": None}, "tj_max", required=False) is None
    with pytest.raises(KeyError):
        units.temperature_C_from_mapping({}, "t_ambient")
    with pytest.raises(ValueError, match="not both"):
        units.temperature_C_from_mapping({"t_ambient_C": 40.0, "t_ambient_K": 313.15}, "t_ambient")
    with pytest.raises(ValueError, match="kelvin"):
        units.temperature_C_from_mapping({"t_ambient_C": 313.15}, "t_ambient", max_C=units.MAX_AMBIENT_C)


def test_flow_validation_hints_at_lpm_and_cfm():
    assert units.validate_flow_m3_s(2.0 / 60000.0, "q") == pytest.approx(3.3333e-5, rel=1e-4)
    with pytest.raises(ValueError, match="L/min"):
        units.validate_flow_m3_s(2.0, "volumetric_flow")  # 2 L/min typed as 2 m^3/s
    with pytest.raises(ValueError, match="CFM"):
        units.validate_flow_m3_s(30.0, "volumetric_flow")  # 30 CFM typed as 30 m^3/s
    with pytest.raises(ValueError, match="positive"):
        units.validate_flow_m3_s(0.0, "volumetric_flow")
    assert units.require_positive(1e-3, "x") == 1e-3
    assert units.require_non_negative(0.0, "x") == 0.0
    with pytest.raises(ValueError):
        units.require_non_negative(-1e-9, "x")


def test_network_and_heatsink_dataclasses_apply_the_conventions():
    with pytest.raises(ValueError, match="kelvin"):
        net.SharedHeatSink(0.5, 313.15)
    with pytest.raises(ValueError, match="kelvin"):
        net.DeviceThermalPath("S1", 10.0, 1.0, tj_max_C=848.15)
    assert net.SharedHeatSink.from_dict({"r_sa": 0.5, "t_ambient_K": 313.15}).t_ambient_C == pytest.approx(40.0)
    dev = net.DeviceThermalPath.from_dict({"name": "S1", "power_W": 1.0, "r_jc": 1.0, "tj_design_K": 398.15})
    assert dev.tj_design_C == pytest.approx(125.0) and dev.tj_max_C is None
    with pytest.raises(ValueError, match="kelvin"):
        net.required_sink_resistance([dev], 313.15, 125.0)
    with pytest.raises(ValueError, match="kelvin"):
        hm.OperatingPoint(t_ambient_C=313.15, cooling="natural", delta_T_K=40.0)
    with pytest.raises(ValueError, match="kelvin"):
        et.iterate([et.DeviceLossModel("D", conduction=lambda tj: 1.0, switching=lambda tj: 0.0)], lambda p: {"D": 50.0}, 698.15)


def test_cooling_and_screening_dataclasses_apply_the_conventions():
    with pytest.raises(ValueError, match="kelvin"):
        cl.LiquidColdPlateOption(0.1, 0.1, 0.004, 16, 0.003, 0.005, 0.0025, 318.15, 3.3e-5)
    with pytest.raises(ValueError, match="L/min"):
        cl.LiquidColdPlateOption(0.1, 0.1, 0.004, 16, 0.003, 0.005, 0.0025, 45.0, 2.0)  # 2 L/min given as 2 m^3/s
    liquid = cl.LiquidColdPlateOption.from_dict({"length": 0.1, "width": 0.1, "base_thickness": 0.004, "n_channels": 16,
                                                  "channel_width": 0.003, "channel_height": 0.005, "wall_thickness": 0.0025,
                                                  "coolant_inlet_K": 318.15, "volumetric_flow": 3.3e-5})
    assert liquid.coolant_inlet_C == pytest.approx(45.0)
    with pytest.raises(ValueError, match="positive"):
        cl.PlateFinGeometry(0.1, 0.1, -0.04, 0.0015, 0.003, 0.006)
    geometry = cl.PlateFinGeometry(0.1, 0.1, 0.04, 0.0015, 0.003, 0.006)
    with pytest.raises(ValueError, match="CFM"):
        cl.ForcedAirOption(geometry, volumetric_flow=25.0)  # CFM given as m^3/s
    with pytest.raises(ValueError, match="non-negative"):
        cl.ForcedAirOption(geometry, fan_curve=[(0.0, 50.0), (0.01, -5.0)])
    assert cl.ForcedAirOption(geometry, fan_curve=[(0.0, 50.0), (0.01, 0.0)]).fan_curve[-1] == (0.01, 0.0)
    assert geometry.to_sink().n_fins == hm.PlateFinHeatSink(0.1, 0.1, 0.04, 0.0015, 0.003, 0.006).n_fins
    dev = cl.CoolingDevice.from_dict({"name": "S1", "power_W": 10.0, "r_jc": 0.5, "footprint_area_m2": 3.2e-4, "tj_max_K": 448.15})
    assert dev.tj_max_C == pytest.approx(175.0)
    with pytest.raises(ValueError, match="footprint"):
        cl.CoolingDevice("S1", 10.0, 0.5, 0.1, 0.0)
    task = screening.StageATask.from_dict({"devices": [{"name": "S1", "power_W": 10.0, "r_jc": 0.5, "footprint_area_m2": 3.2e-4}],
                                           "t_ambient_K": 313.15, "tj_limit_K": 398.15})
    assert task.t_ambient_C == pytest.approx(40.0) and task.tj_limit_C == pytest.approx(125.0)
    with pytest.raises(ValueError, match="flow_bounds"):
        screening.StageATask(devices=task.devices, cooling="forced", flow_bounds=(0.02, 0.01))
    with pytest.raises(ValueError, match="L/min"):
        screening.StageATask(devices=task.devices, cooling="forced", volumetric_flow=3.0)
