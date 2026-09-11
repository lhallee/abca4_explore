"""Project configuration and secret loading."""

from __future__ import annotations

import os

from dataclasses import dataclass, field
from pathlib import Path

from dotenv import dotenv_values


@dataclass(frozen=True, kw_only=True)
class Credentials:
    """API credentials held only in process memory."""

    ncbi_api_key: str | None = field(default=None, repr=False)
    alphagenome_api_key: str | None = field(default=None, repr=False)
    ncbi_email: str | None = field(default=None, repr=False)


@dataclass(frozen=True, kw_only=True)
class ProjectPaths:
    """Canonical project output paths."""

    root: Path
    raw: Path
    cache: Path
    processed: Path
    tables: Path
    plots: Path
    logs: Path
    reference: Path

    @classmethod
    def from_root(cls, root: Path) -> ProjectPaths:
        root = root.resolve()
        return cls(
            root=root,
            raw=root / "data" / "raw",
            cache=root / "data" / "cache",
            processed=root / "data" / "processed",
            tables=root / "results" / "tables",
            plots=root / "plots",
            logs=root / "results" / "logs",
            reference=root / "data" / "reference",
        )

    def create(self) -> None:
        """Create output directories."""
        for path in (
            self.raw,
            self.cache,
            self.processed,
            self.tables,
            self.plots,
            self.logs,
            self.reference,
        ):
            path.mkdir(parents=True, exist_ok=True)


def _clean_secret(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def redact_text(text: str, credentials: Credentials) -> str:
    """Remove known credential values from text before it is persisted."""
    redacted = text
    for secret in (
        credentials.ncbi_api_key,
        credentials.alphagenome_api_key,
        credentials.ncbi_email,
    ):
        if secret:
            redacted = redacted.replace(secret, "[REDACTED]")
    return redacted


def load_credentials(root: Path) -> Credentials:
    """Load canonical variables and project aliases without printing values."""
    secret_path = root.resolve() / ".secrets.env"
    secret_values = dotenv_values(secret_path) if secret_path.exists() else {}

    ncbi_api_key = _clean_secret(
        os.environ.get("NCBI_API_KEY")
        or secret_values.get("NCBI_API_KEY")
        or secret_values.get("NCBI")
    )
    alphagenome_api_key = _clean_secret(
        os.environ.get("ALPHAGENOME_API_KEY")
        or secret_values.get("ALPHAGENOME_API_KEY")
        or secret_values.get("ALPHA_GENOME")
    )
    ncbi_email = _clean_secret(
        os.environ.get("NCBI_EMAIL") or secret_values.get("NCBI_EMAIL")
    )

    return Credentials(
        ncbi_api_key=ncbi_api_key,
        alphagenome_api_key=alphagenome_api_key,
        ncbi_email=ncbi_email,
    )
