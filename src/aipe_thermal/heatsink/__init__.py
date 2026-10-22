# SPDX-License-Identifier: CC-BY-4.0
"""Plate-fin heat-sink correlation model (CC BY 4.0; original research by Fulong Li).

This package is the authoritative home of the former stand-alone
``heatsink_model.py`` (model version ``__version__`` = 1.0.0, see
``CHANGELOG.md`` for the corrections relative to the original ``HS_main.py``).
Everything public of :mod:`aipe_thermal.heatsink.model` is re-exported here, so
``import aipe_thermal.heatsink as hm`` is a drop-in replacement for the module
that ``pea.thermal._sibling.load_heatsink_model()`` used to return.
"""

from aipe_thermal.heatsink.model import *  # noqa: F401,F403
from aipe_thermal.heatsink.model import (  # noqa: F401
    G, MATERIAL_ALIASES, OBJECTIVES, P_ATM, R_AIR, SIGMA_SB, T0_K, VARIABLE_ORDER, __all__, __attribution__,
    __license__, __version__,
)
