"""Credential leak checks that never expose credential values."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from .config import Credentials


EXCLUDED_NAMES = frozenset(
    {
        ".git",
        ".pytest_cache",
        ".ruff_cache",
        ".secrets.env",
        ".venv",
        "__pycache__",
    }
)


def _secret_bytes(credentials: Credentials) -> tuple[bytes, ...]:
    values = (
        credentials.ncbi_api_key,
        credentials.alphagenome_api_key,
        credentials.ncbi_email,
    )
    return tuple(
        value.encode("utf-8")
        for value in values
        if value is not None and len(value) >= 8
    )


def _candidate_files(root: Path) -> Iterable[Path]:
    for path in root.rglob("*"):
        relative = path.relative_to(root)
        if any(
            part in EXCLUDED_NAMES or part.startswith(".uv-cache")
            for part in relative.parts
        ):
            continue
        if path.is_file():
            yield path


def find_credential_leaks(root: Path, credentials: Credentials) -> list[Path]:
    """Return files containing loaded credentials without exposing the values."""
    secrets = _secret_bytes(credentials)
    if not secrets:
        return []

    matches: list[Path] = []
    for path in _candidate_files(root.resolve()):
        try:
            content = path.read_bytes()
        except OSError:
            continue
        if any(secret in content for secret in secrets):
            matches.append(path)
    return matches
