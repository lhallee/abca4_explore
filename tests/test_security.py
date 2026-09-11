from __future__ import annotations

from pathlib import Path

from abca4_avi.config import Credentials
from abca4_avi.security import find_credential_leaks


def test_secret_file_is_excluded_but_artifact_leak_is_found(tmp_path: Path) -> None:
    credential = "test-secret-value"
    credentials = Credentials(alphagenome_api_key=credential)
    (tmp_path / ".secrets.env").write_text(
        f"ALPHA_GENOME={credential}\n", encoding="utf-8"
    )

    assert find_credential_leaks(tmp_path, credentials) == []

    artifact = tmp_path / "result.txt"
    artifact.write_text(credential, encoding="utf-8")

    assert find_credential_leaks(tmp_path, credentials) == [artifact]
