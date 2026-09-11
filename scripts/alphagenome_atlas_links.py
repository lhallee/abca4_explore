"""Compatibility entry point for AlphaGenome Atlas links."""

from __future__ import annotations

import sys

from abca4_avi.cli import main


if __name__ == "__main__":
    main(["links", *sys.argv[1:]])
