# SPDX-License-Identifier: MIT
"""``python -m aipe_thermal`` -- same interface as the ``aipe-thermal`` console script."""

import sys

from aipe_thermal.cli import main

if __name__ == "__main__":
    sys.exit(main())
