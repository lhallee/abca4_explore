"""Compatibility entry point for AlphaGenome AVI annotation."""

from __future__ import annotations

import sys

from abca4_avi.cli import main


if __name__ == "__main__":
    main(["annotate", *sys.argv[1:]])
