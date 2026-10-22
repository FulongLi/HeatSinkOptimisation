# SPDX-License-Identifier: MIT
"""Deprecation shim for the pre-0.2 layout (``pea.thermal`` plus a sibling ``HeatSinkOptimisation`` checkout).

The heat-sink correlation model is bundled as ``aipe_thermal.heatsink``: no
sibling checkout, directory name or environment variable is consulted. The
legacy variable ``HEATSINK_OPT_DIR`` is only *inspected* so that users who still
set it get a ``DeprecationWarning``; its value is never used.
"""

from __future__ import annotations

import os
import warnings
from typing import Mapping

LEGACY_ENV_VAR = "HEATSINK_OPT_DIR"


def check_legacy_environment(environ: Mapping[str, str] | None = None) -> str | None:
    """Warn (``DeprecationWarning``) when the legacy ``HEATSINK_OPT_DIR`` variable is set; return the message or ``None``."""
    env = os.environ if environ is None else environ
    value = env.get(LEGACY_ENV_VAR)
    if not value:
        return None
    message = (f"{LEGACY_ENV_VAR}={value!r} is ignored: since aipe-thermal 0.2 the heat-sink correlation model is "
               "bundled as aipe_thermal.heatsink and no sibling checkout is searched. Unset the variable.")
    warnings.warn(message, DeprecationWarning, stacklevel=3)
    return message


__all__ = ["LEGACY_ENV_VAR", "check_legacy_environment"]
