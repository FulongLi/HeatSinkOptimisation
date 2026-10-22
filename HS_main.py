# SPDX-License-Identifier: CC-BY-4.0
"""Thin wrapper kept for the documented ``python HS_main.py ...`` commands (Dr. Fulong Li, CC BY 4.0).

The implementation lives in :mod:`aipe_thermal.heatsink.frontend` (same options);
``aipe-thermal heatsink-front ...`` is the packaged equivalent. Requires
``pip install -e .`` (or an installed ``aipe-thermal`` distribution).
"""

from aipe_thermal.heatsink.frontend import build_parser, main  # noqa: F401

if __name__ == "__main__":
    raise SystemExit(main())
