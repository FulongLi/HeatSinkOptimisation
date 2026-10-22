# SPDX-License-Identifier: MIT
"""Unit conventions and input validation for ``aipe_thermal`` (SI).

Conventions used throughout the package:

* absolute temperatures carry an explicit unit suffix: ``*_C`` for degrees
  Celsius, ``*_K`` for kelvin; temperature *differences* are kelvin
  (``delta_T_K``, ``margin_K``, ``rise_K``);
* lengths in metres, areas in square metres, volumetric flows in cubic metres
  per second (``1 L/min = 1/60000 m^3/s``, ``1 CFM = 4.719e-4 m^3/s``),
  pressures in pascal, powers in watt, thermal resistances in K/W.

Dictionary inputs (CLI JSON, ``from_dict`` constructors) may give an absolute
temperature either as ``<name>_C`` or as ``<name>_K``; exactly one of the two
is accepted and kelvin values are converted. Plain float arguments are always
degrees Celsius for absolute temperatures; values that are implausible for
electronics cooling (ambient or coolant above ``MAX_AMBIENT_C``, any absolute
temperature above ``MAX_PLAUSIBLE_C``) are rejected with the hint that they
look like kelvin, and flows above ``MAX_PLAUSIBLE_FLOW_M3_S`` with the hint
that they look like L/min or CFM.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

T0_K = 273.15
ABSOLUTE_ZERO_C = -T0_K
MAX_AMBIENT_C = 300.0  # ambient / coolant inlet above this is almost certainly a kelvin value (313.15 "degC" = 40 degC)
MAX_PLAUSIBLE_C = 600.0  # above this any "absolute temperature in degC" is almost certainly a kelvin value
MAX_PLAUSIBLE_FLOW_M3_S = 1.0  # 60 000 L/min / 2100 CFM: above this a flow was almost certainly given in L/min or CFM
M3_S_PER_LPM = 1.0 / 60000.0
M3_S_PER_CFM = 0.3048 ** 3 / 60.0


def kelvin_from_celsius(t_C: float) -> float:
    """``T[K] = T[degC] + 273.15``."""
    return float(t_C) + T0_K


def celsius_from_kelvin(t_K: float) -> float:
    """``T[degC] = T[K] - 273.15``."""
    return float(t_K) - T0_K


def m3_s_from_lpm(q_lpm: float) -> float:
    """Litres per minute to cubic metres per second."""
    return float(q_lpm) * M3_S_PER_LPM


def m3_s_from_cfm(q_cfm: float) -> float:
    """Cubic feet per minute to cubic metres per second."""
    return float(q_cfm) * M3_S_PER_CFM


def require_finite(value: Any, name: str) -> float:
    """Return ``value`` as a finite float or raise ``ValueError`` naming the quantity."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a real number, got {value!r}")
    try:
        v = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a real number, got {value!r}") from exc
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return v


def require_positive(value: Any, name: str) -> float:
    v = require_finite(value, name)
    if v <= 0.0:
        raise ValueError(f"{name} must be positive, got {v!r}")
    return v


def require_non_negative(value: Any, name: str) -> float:
    v = require_finite(value, name)
    if v < 0.0:
        raise ValueError(f"{name} must be non-negative, got {v!r}")
    return v


def _kelvin_key(name: str) -> str:
    return name[:-2] + "_K" if name.endswith("_C") else name + "_K"


def validate_temperature_C(value: Any, name: str, max_C: float = MAX_PLAUSIBLE_C) -> float:
    """Validate an absolute temperature given in degrees Celsius.

    Rejects non-finite values, values at or below absolute zero and values
    above ``max_C`` (default 600 degC; use ``MAX_AMBIENT_C`` for ambient and
    coolant temperatures), the latter with the hint that the number looks like
    a kelvin value.
    """
    v = require_finite(value, name)
    if v <= ABSOLUTE_ZERO_C:
        raise ValueError(f"{name} = {v} degC is at or below absolute zero")
    if v > max_C:
        raise ValueError(
            f"{name} = {v} is not a plausible temperature in degrees Celsius (> {max_C} degC). Absolute "
            f"temperatures are degrees Celsius; if the value is in kelvin pass '{_kelvin_key(name)}' or convert "
            f"with celsius_from_kelvin()."
        )
    return v


def validate_ambient_C(value: Any, name: str) -> float:
    """Validate an ambient / coolant-inlet temperature in degrees Celsius (at most ``MAX_AMBIENT_C``)."""
    return validate_temperature_C(value, name, MAX_AMBIENT_C)


def validate_temperature_K(value: Any, name: str, max_K: float = MAX_PLAUSIBLE_C + T0_K) -> float:
    """Validate an absolute temperature given in kelvin (positive, at most ``max_K``)."""
    v = require_positive(value, name)
    if v > max_K:
        raise ValueError(f"{name} = {v} K is not a plausible temperature (> {max_K} K)")
    return v


def validate_flow_m3_s(value: Any, name: str, max_m3_s: float = MAX_PLAUSIBLE_FLOW_M3_S) -> float:
    """Validate a volumetric flow given in m^3/s (positive, finite, plausible)."""
    v = require_positive(value, name)
    if v > max_m3_s:
        raise ValueError(
            f"{name} = {v} m^3/s is not a plausible volumetric flow (> {max_m3_s} m^3/s). Flows are m^3/s; convert "
            f"with m3_s_from_lpm() (1 L/min = {M3_S_PER_LPM:.4e} m^3/s) or m3_s_from_cfm() (1 CFM = {M3_S_PER_CFM:.4e} m^3/s)."
        )
    return v


def temperature_C_from_mapping(data: Mapping[str, Any], base: str, *, default: float | None = None,
                               required: bool = True, max_C: float = MAX_PLAUSIBLE_C) -> float | None:
    """Read ``<base>_C`` or ``<base>_K`` (exactly one) from ``data`` and return degrees Celsius.

    ``default`` (degrees Celsius) is used when neither key is present; without a
    default a missing required temperature raises ``KeyError`` and an optional
    one returns ``None``. ``max_C`` is the plausibility limit applied to the
    Celsius value (``MAX_AMBIENT_C`` for ambient and coolant temperatures).
    """
    key_c, key_k = f"{base}_C", f"{base}_K"
    has_c = data.get(key_c) is not None
    has_k = data.get(key_k) is not None
    if has_c and has_k:
        raise ValueError(f"give either {key_c} (degrees Celsius) or {key_k} (kelvin), not both")
    if has_c:
        return validate_temperature_C(data[key_c], key_c, max_C)
    if has_k:
        return validate_temperature_C(celsius_from_kelvin(validate_temperature_K(data[key_k], key_k)), key_c, max_C)
    if default is not None:
        return validate_temperature_C(default, key_c, max_C)
    if required:
        raise KeyError(f"missing temperature: give {key_c} (degrees Celsius) or {key_k} (kelvin)")
    return None


__all__ = ["T0_K", "ABSOLUTE_ZERO_C", "MAX_AMBIENT_C", "MAX_PLAUSIBLE_C", "MAX_PLAUSIBLE_FLOW_M3_S", "M3_S_PER_LPM",
           "M3_S_PER_CFM", "kelvin_from_celsius", "celsius_from_kelvin", "m3_s_from_lpm", "m3_s_from_cfm",
           "require_finite", "require_positive", "require_non_negative", "validate_temperature_C", "validate_ambient_C",
           "validate_temperature_K", "validate_flow_m3_s", "temperature_C_from_mapping"]
