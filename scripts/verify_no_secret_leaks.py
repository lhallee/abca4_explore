"""Verify that loaded credentials do not occur in project artifacts."""

from __future__ import annotations

import argparse

from pathlib import Path

from abca4_avi.config import load_credentials
from abca4_avi.security import find_credential_leaks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    arguments = parser.parse_args()
    root = arguments.root.resolve()
    matches = find_credential_leaks(root, load_credentials(root))
    if matches:
        print(f"Credential leak check failed in {len(matches)} file(s).")
        for path in matches:
            print(path.relative_to(root))
        raise SystemExit(1)
    print("Credential leak check passed: 0 matching files.")


if __name__ == "__main__":
    main()
