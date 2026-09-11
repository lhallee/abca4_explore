"""Compatibility entry point for ClinVar commands."""

from __future__ import annotations

import sys

from abca4_avi.cli import main


if __name__ == "__main__":
    main(["clinvar", *sys.argv[1:]])
