"""Compatibility entry point for the complete analysis."""

from __future__ import annotations

import sys

from abca4_avi.cli import main


if __name__ == "__main__":
    main(["analysis", *sys.argv[1:]])
