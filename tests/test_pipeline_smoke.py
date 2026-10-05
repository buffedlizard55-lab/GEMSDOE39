"""Smoke tests: metric identity, SPRT boundaries, validator on the shipped tif."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems39 import metric, grid  # noqa: E402
from gems39.sprt_select import sprt_pairwise  # noqa: E402


def test_metric_perfect():
    pred = np.zeros((5, 5), np.float32)
    pred[2, 2] = 1.0
    truth = np.zeros((5, 5), bool)
    truth[2, 2] = True
    r = metric.dti_components(pred, truth)
    assert r["tp"] == 1.0
    assert r["fn"] == 0.0
    assert abs(r["dti"] - 1.0) < 1e-6


def test_metric_miss():
    pred = np.zeros((5, 5), np.float32)
    truth = np.zeros((5, 5), bool)
    truth[2, 2] = True
    r = metric.dti_components(pred, truth)
    assert r["tp"] == 0.0
    assert r["dti"] == 0.0


def test_sprt_boundaries():
    r = sprt_pairwise([True] * 8, p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    # 8 wins under p0=0.5,p1=0.7 is LLR = 8*log(0.7/0.5) = 8*0.3365 = 2.69
    assert abs(r["llr"] - 8 * np.log(0.7 / 0.5)) < 1e-6
    assert r["decision"] == "continue"  # upper ~ 2.89
    r2 = sprt_pairwise([True] * 9, p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    assert r2["decision"] == "accept_H1"


def test_submission_format():
    sub = ROOT / "docs" / "downloads"
    tifs = list(sub.glob("gemsdoe39-*-zeros.tif"))
    assert tifs, "no submission built"
    p = tifs[-1]
    import rasterio
    with rasterio.open(p) as s:
        assert s.count == 1
        assert s.dtypes[0] == "float32"
        assert str(s.crs).upper().startswith("EPSG:32611")
        assert s.height == 3730 and s.width == 3292
        a = s.read(1)
    sample = ROOT / "data" / "sample_submission.tif"
    if sample.exists():
        with rasterio.open(sample) as r:
            foot = np.isfinite(r.read(1))
            assert tuple(s.transform)[:6] == tuple(r.transform)[:6]
        assert np.isfinite(a[foot]).all()
        assert (a[foot] >= 0).all() and (a[foot] <= 1).all()
        assert (a[~foot] == 0).all()


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_") and callable(v):
            v()
            print("PASS", k)
