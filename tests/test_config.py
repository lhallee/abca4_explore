from __future__ import annotations

from abca4_avi.config import Credentials, redact_text


def test_secret_redaction_removes_all_known_values() -> None:
    credentials = Credentials(
        ncbi_api_key="ncbi-secret",
        alphagenome_api_key="alpha-secret",
        ncbi_email="private@example.org",
    )

    output = redact_text("ncbi-secret alpha-secret private@example.org", credentials)

    assert output == "[REDACTED] [REDACTED] [REDACTED]"
    assert "secret" not in output
