"""Unit/integration tests for DTI, SPRT, blocked folds, and submission format."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems39 import emission, grid, holdout, metric  # noqa: E402
from gems39.sprt_select import PairwiseSPRT, sprt_pairwise  # noqa: E402


def test_metric_perfect_and_miss():
    pred = np.zeros((5, 5), np.float32)
    pred[2, 2] = 1.0
    truth = np.zeros((5, 5), bool)
    truth[2, 2] = True
    result = metric.dti_components(pred, truth)
    assert result["tp"] == 1.0
    assert result["fn"] == 0.0
    assert abs(result["dti"] - 1.0) < 1e-6
    missed = metric.dti_components(np.zeros_like(pred), truth)
    assert missed["tp"] == 0.0
    assert missed["dti"] == 0.0


def test_binary_and_soft_metric_agree_for_binary_predictions():
    truth = np.zeros((9, 11), bool)
    truth[3:6, 4] = True
    truth[6, 5:8] = True
    pred = np.zeros_like(truth)
    pred[3:6, 5] = True
    pred[6, 5:8] = True
    soft = metric.dti_components(pred.astype(np.float32), truth)
    binary = metric.dti_binary(pred, truth)
    for key in ("tp", "fp", "fn", "dti"):
        assert abs(soft[key] - binary[key]) < 1e-9, (key, soft[key], binary[key])


def test_sprt_anytime_boundaries_and_generator_accounting():
    result = sprt_pairwise((True for _ in range(8)), p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    assert result["decision"] == "continue"
    assert result["n"] == 8 and result["wins"] == 8
    assert abs(result["llr"] - 8 * np.log(0.7 / 0.5)) < 1e-12
    assert abs(result["upper"] - np.log(20.0)) < 1e-12
    assert abs(result["lower"] - np.log(0.10)) < 1e-12
    assert sprt_pairwise([True] * 9)["decision"] == "accept_H1"
    assert sprt_pairwise([False] * 5)["decision"] == "accept_H0"


def test_incremental_sprt_refuses_post_boundary_peeking():
    test = PairwiseSPRT()
    for _ in range(9):
        state = test.update(True)
    assert state["decision"] == "accept_H1"
    try:
        test.update(True)
    except RuntimeError:
        pass
    else:
        raise AssertionError("SPRT must refuse observations after crossing a boundary")


def test_exact_known_fault_mask_zero_buffer():
    foot = np.ones((9, 9), bool)
    catalogue = np.zeros_like(foot)
    catalogue[4, 4] = True
    blocked = emission._catalogue_block(catalogue, foot, 0)
    assert blocked.sum() == 1
    assert blocked[4, 4]
    assert not blocked[4, 5]


def test_read_bands_selects_metadata_names_and_zeros_outside(tmp_path):
    profile = {
        "driver": "GTiff",
        "width": 4,
        "height": 3,
        "count": 2,
        "dtype": "float32",
        "crs": "EPSG:32611",
        "transform": from_origin(0, 300, 100, 100),
        "nodata": np.float32(-3.4028234663852886e38),
    }
    path = tmp_path / "features.tif"
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(np.arange(12, dtype=np.float32).reshape(3, 4), 1)
        dst.write(np.arange(12, dtype=np.float32).reshape(3, 4) + 100, 2)
        dst.set_band_description(1, "geod_test - test band")
        dst.set_band_description(2, "other_test - another band")
        dst.update_tags(1, band_name="geod_test")
        dst.update_tags(2, band_name="other_test")
    footprint = np.ones((3, 4), bool)
    footprint[0, 0] = False
    bands = grid.read_bands(path, footprint, names=["GEOD_TEST"])
    assert set(bands) == {"geod_test"}
    assert np.isfinite(bands["geod_test"]).all()
    assert bands["geod_test"][0, 0] == 0.0
    try:
        grid.read_bands(path, footprint, names=["not_a_band"])
    except ValueError:
        pass
    else:
        raise AssertionError("Unknown metadata band must be rejected, not guessed")


def _write_tiny_holdout(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "width": 40,
        "height": 40,
        "count": 1,
        "crs": "EPSG:32611",
        "transform": from_origin(0, 4000, 100, 100),
        "compress": "deflate",
    }
    sample = np.zeros((40, 40), np.float32)
    with rasterio.open(root / "sample_submission.tif", "w", dtype="float32", nodata=np.nan, **profile) as dst:
        dst.write(sample, 1)
    labels = np.zeros((40, 40), np.int8)
    labels[7:9, 7] = 1
    labels[7:9, 32] = 1
    labels[32, 7:9] = 1
    labels[32, 32:34] = 1
    with rasterio.open(root / "labels.tif", "w", dtype="int8", nodata=-1, **profile) as dst:
        dst.write(labels, 1)


def test_spatial_holdout_hides_whole_components_without_label_leakage(tmp_path):
    _write_tiny_holdout(tmp_path)
    ctx = holdout.load_holdout(tmp_path, hide_fraction=0.5, seed=7, collar_px=2)
    assert [fold.name for fold in ctx.folds] == ["NW", "NE", "SW", "SE"]
    for fold in ctx.folds:
        assert fold.n_truth > 0
        assert not np.any(fold.truth & fold.known)
        assert np.all(fold.truth[~fold.domain] == 0)
        pred = np.zeros(ctx.footprint.shape, bool)
        pred[fold.bbox] = fold.truth
        result = holdout.evaluate_fold(pred, fold)
        assert abs(result["dti"] - 1.0) < 1e-6


def test_submission_writer_nan_outside_and_audit(tmp_path):
    profile = {
        "driver": "GTiff",
        "width": 5,
        "height": 4,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:32611",
        "transform": from_origin(0, 400, 100, 100),
        "nodata": np.nan,
    }
    sample_path = tmp_path / "sample.tif"
    sample = np.full((4, 5), np.nan, np.float32)
    sample[1:3, 1:4] = 0
    with rasterio.open(sample_path, "w", **profile) as dst:
        dst.write(sample, 1)
    foot = np.isfinite(sample)
    pred = np.zeros((4, 5), np.float32)
    pred[2, 3] = 1
    out_path = tmp_path / "submission.tif"
    grid.write_submission(pred, sample_path, out_path, foot)
    audit = grid.audit_submission(out_path, foot, sample_path)
    assert audit["ok"]
    with rasterio.open(out_path) as ds:
        written = ds.read(1)
        assert np.isnan(written[~foot]).all()
        assert ds.nodata is not None and np.isnan(ds.nodata)
    try:
        grid.write_submission(pred, sample_path, tmp_path / "bad.tif", foot, outside="zero")
    except ValueError:
        pass
    else:
        raise AssertionError("Zero-filled outside output must be rejected")
    bad_range = pred.copy()
    bad_range[2, 3] = 1.01
    try:
        grid.write_submission(bad_range, sample_path, tmp_path / "bad_range.tif", foot)
    except ValueError:
        pass
    else:
        raise AssertionError("In-footprint values above 1 must be rejected")


def test_current_candidate_submission_format_if_built():
    candidates = sorted((ROOT / "docs" / "downloads").glob("gemsdoe39-h39x01-*-nan.tif"))
    if not candidates:
        return
    sample_path = ROOT / "data" / "sample_submission.tif"
    if not sample_path.exists():
        return
    path = candidates[-1]
    with rasterio.open(path) as ds, rasterio.open(sample_path) as ref:
        values = ds.read(1)
        footprint = np.isfinite(ref.read(1))
        assert ds.count == 1 and ds.dtypes == ("float32",)
        assert ds.crs == ref.crs
        assert ds.shape == ref.shape
        assert tuple(ds.transform)[:6] == tuple(ref.transform)[:6]
        assert np.isfinite(values[footprint]).all()
        assert ((values[footprint] >= 0) & (values[footprint] <= 1)).all()
        assert np.isnan(values[~footprint]).all()
        assert ds.nodata is not None and np.isnan(ds.nodata)


if __name__ == "__main__":
    for name, value in list(globals().items()):
        if name.startswith("test_") and callable(value):
            value()
            print("PASS", name)
