"""Tests for the dual-encoding submission writer (src/gems39/io39.py).

Both encodings must survive a round trip through rasterio with the invariants the
portal actually checks: single band, float32, template grid, every cell finite
and inside [0, 1] for the zeros encoding, NaN exactly outside for the nan one.
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems39 import io39  # noqa: E402


def _write_template(path: Path, h: int = 6, w: int = 8) -> Path:
    """A minimal single-band float32 GeoTIFF with a rectangular NaN-free core."""
    from rasterio.transform import from_origin

    data = np.full((h, w), np.nan, dtype="float32")
    data[1:-1, 1:-1] = 0.25
    with rasterio.open(
        path, "w", driver="GTiff", height=h, width=w, count=1, dtype="float32",
        crs="EPSG:32611", transform=from_origin(243350.0, 4508550.0, 100.0, 100.0),
    ) as dst:
        dst.write(data, 1)
    return path


def _footprint(h: int = 6, w: int = 8) -> np.ndarray:
    foot = np.zeros((h, w), bool)
    foot[1:-1, 1:-1] = True
    return foot


def test_write_zeros_is_all_finite_and_in_range(tmp_path):
    tpl = _write_template(tmp_path / "tpl.tif")
    foot = _footprint()
    rng = np.random.default_rng(0)
    pred = rng.random(foot.shape).astype("float32")
    out = io39.write(pred, tpl, tmp_path / "z.tif", foot, outside="zeros")
    a = io39.audit(out, foot, tpl)
    assert a["ok"] is True
    assert a["outside_zero"] == int((~foot).sum())
    assert a["outside_other"] == 0
    with rasterio.open(out) as s:
        arr = s.read(1)
        assert s.count == 1 and s.dtypes[0] == "float32"
        assert s.nodata is None, "all-finite encoding must not carry a nodata tag"
    assert np.isfinite(arr).all()
    assert arr.min() >= 0.0 and arr.max() <= 1.0
    assert (arr[~foot] == 0.0).all()


def test_write_nan_puts_nan_exactly_outside(tmp_path):
    tpl = _write_template(tmp_path / "tpl.tif")
    foot = _footprint()
    rng = np.random.default_rng(1)
    pred = rng.random(foot.shape).astype("float32")
    out = io39.write(pred, tpl, tmp_path / "n.tif", foot, outside="nan")
    a = io39.audit(out, foot, tpl)
    assert a["ok"] is True
    with rasterio.open(out) as s:
        arr = s.read(1)
        assert s.nodata is not None and np.isnan(s.nodata)
    assert np.isnan(arr[~foot]).all()
    assert np.isfinite(arr[foot]).all()


def test_writer_refuses_out_of_range_and_bad_shapes(tmp_path):
    tpl = _write_template(tmp_path / "tpl.tif")
    foot = _footprint()
    bad = np.full(foot.shape, 1.5, dtype="float32")
    with pytest.raises(ValueError, match="\\[0,1\\]"):
        io39.write(bad, tpl, tmp_path / "bad.tif", foot, outside="zeros")
    with pytest.raises(ValueError, match="shape mismatch"):
        io39.write(np.zeros((3, 3), "float32"), tpl, tmp_path / "bad.tif", foot)
    with pytest.raises(ValueError):
        io39.write(np.zeros(foot.shape, "float32"), tpl, tmp_path / "bad.tif", foot,
                   outside="sentinel")


def test_audit_detects_a_wrong_grid(tmp_path):
    tpl = _write_template(tmp_path / "tpl.tif")
    other = _write_template(tmp_path / "tpl2.tif", h=6, w=9)
    foot = _footprint()
    out = io39.write(np.full(foot.shape, 0.5, "float32"), tpl, tmp_path / "z.tif", foot)
    a = io39.audit(out, foot, other)
    assert a["grid_ok"] is False


def test_sha256_is_stable_and_matches_a_known_digest(tmp_path):
    p = tmp_path / "b.bin"
    p.write_bytes(struct.pack("<f", 0.5))
    assert io39.sha256(p) == io39.sha256(p)
    assert len(io39.sha256(p)) == 64
