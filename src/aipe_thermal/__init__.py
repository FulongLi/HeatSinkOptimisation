# SPDX-License-Identifier: MIT
"""Thermal design evaluation for power semiconductors and their cooling (SI units).

Modules:

- ``aipe_thermal.heatsink`` -- plate-fin heat-sink correlation model: natural
  (Elenbaas / Bar-Cohen-Rohsenow) and forced (Teertstra-Yovanovich-Culham,
  Gnielinski) convection, pressure drop and fan operating point, fin efficiency,
  Lee-Song-Au-Moran spreading, heat-load solve and seeded Pareto optimiser.
  CC BY 4.0, original research by Fulong Li (``heatsink_model`` 1.0.0).
- ``aipe_thermal.network`` -- steady-state junction -> case -> interface ->
  shared heat sink -> ambient networks, Foster/Cauer transient impedances (MIT).
- ``aipe_thermal.electrothermal`` -- loss(Tj) <-> Tj fixed-point iteration with
  pluggable loss callables and datasheet-based defaults (MIT).
- ``aipe_thermal.cooling`` -- natural convection, forced air and liquid cold
  plate comparison for devices sharing one sink (MIT).
- ``aipe_thermal.screening`` -- Stage A analytical Pareto screening of plate-fin
  geometries with junction checks and FEM-verification flags (MIT).
- ``aipe_thermal.cases.dab_10kw`` -- the AIPE 10 kW DAB case study (MIT).
- ``aipe_thermal.units`` -- unit conventions and input validation (MIT).
- ``aipe-thermal`` / ``python -m aipe_thermal`` -- JSON command-line interface.

Results are analytical correlation estimates with stated validity; none of them
is a CFD, FEM or measured value. See ``LICENSES.md`` for the per-file licences.
"""

from __future__ import annotations

from aipe_thermal._version import __version__
from aipe_thermal import electrothermal, network, units  # noqa: E402  (standard library only)
from aipe_thermal import heatsink  # noqa: E402  (numpy + scipy; CC BY 4.0 correlation model)
from aipe_thermal import cooling, screening  # noqa: E402
from aipe_thermal._compat import check_legacy_environment  # noqa: E402

check_legacy_environment()

__all__ = ["__version__", "cooling", "electrothermal", "heatsink", "network", "screening", "units",
           "check_legacy_environment"]
