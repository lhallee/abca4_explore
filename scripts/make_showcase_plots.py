"""Build leakage-aware ABCA4 AVI presentation figures and assessment."""

from __future__ import annotations

import argparse
import json

from pathlib import Path

from abca4_avi.showcase import build_showcase


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    arguments = parser.parse_args()
    result = build_showcase(arguments.root)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
