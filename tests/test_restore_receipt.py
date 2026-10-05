"""Regression tests for multi-call restoration receipts."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_restore_groups_merge_verified_receipts(tmp_path):
    target = tmp_path / "data"
    target.mkdir()
    entries = []
    for file_id, group, payload in (
        ("core_sample", "core", b"core raster fixture"),
        ("external_grid", "external", b"external raster fixture"),
    ):
        dest = f"{file_id}.bin"
        (target / dest).write_bytes(payload)
        entries.append(
            {
                "id": file_id,
                "dest": dest,
                "group": group,
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        )
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"files": entries}), encoding="utf-8")

    for group in ("core", "external"):
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "restore_data.py"),
                "--manifest",
                str(manifest),
                "--target-dir",
                str(target),
                "--group",
                group,
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    receipt = json.loads((target / "restore_receipt.json").read_text(encoding="utf-8"))
    assert receipt["schema_version"] == 2
    assert receipt["verified_files_count"] == 2
    assert {item["id"] for item in receipt["files"]} == {"core_sample", "external_grid"}
    assert receipt["last_requested_group"] == "external"
    assert Path(receipt["storage_root"]) == target


def test_build_pipeline_refuses_default_holdout_rescore():
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_pipeline.py")],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "holdout is exhausted" in result.stderr
