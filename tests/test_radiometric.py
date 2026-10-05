from __future__ import annotations

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from gems39.radiometric import (
    RATIO_BANDS,
    radiometric_conductivity_score,
    read_conductivity,
    read_ratio_bands,
)


def test_parallel_cross_gradient_is_bounded_and_outside_is_zero():
    n = 128
    yy, xx = np.mgrid[:n, :n]
    edge_x = np.tanh((xx - 64) / 4.0).astype(np.float32)
    edge_y = np.tanh((yy - 64) / 4.0).astype(np.float32)
    footprint = np.ones((n, n), bool)
    footprint[:3, :5] = False
    ratios = {
        "thk": edge_x,
        "uk": np.zeros_like(edge_x),
        "uth": np.ones_like(edge_x),
    }

    parallel = radiometric_conductivity_score(ratios, 2.0 * edge_x, footprint)
    perpendicular = radiometric_conductivity_score(ratios, 2.0 * edge_y, footprint)

    assert parallel.shape == footprint.shape
    assert np.isfinite(parallel).all()
    assert parallel.min() >= 0.0 and parallel.max() <= 1.0
    assert parallel[:, 56:72].max() > 0.8
    assert perpendicular.max() < 1e-5
    assert np.all(parallel[~footprint] == 0.0)
    assert np.all(perpendicular[~footprint] == 0.0)


def test_score_rejects_missing_or_nonfinite_inputs():
    footprint = np.ones((12, 12), bool)
    field = np.ones((12, 12), np.float32)
    with pytest.raises(ValueError, match="exactly the ratio bands"):
        radiometric_conductivity_score({"thk": field}, field, footprint)

    ratios = {name: field.copy() for name in RATIO_BANDS}
    ratios["uk"][3, 4] = np.nan
    with pytest.raises(ValueError, match="Invalid ratio field"):
        radiometric_conductivity_score(ratios, field, footprint)


def _write_one_band(path, data, *, dtype, transform, description, nodata=None):
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=data.shape[1],
        height=data.shape[0],
        count=1,
        dtype=dtype,
        crs="EPSG:32611",
        transform=transform,
        nodata=nodata,
    ) as dst:
        dst.write(data, 1)
        dst.set_band_description(1, description)


def test_ratio_reader_uses_descriptions_not_numeric_order(tmp_path):
    h, w = 5, 6
    transform = from_origin(100, 500, 100, 100)
    ref = tmp_path / "reference.tif"
    _write_one_band(
        ref,
        np.ones((h, w), np.float32),
        dtype="float32",
        transform=transform,
        description="template",
    )
    ratio_path = tmp_path / "ratios.tif"
    values = {"thk": 64, "uk": 128, "uth": 255}
    with rasterio.open(
        ratio_path,
        "w",
        driver="GTiff",
        width=w,
        height=h,
        count=3,
        dtype="uint8",
        crs="EPSG:32611",
        transform=transform,
    ) as dst:
        for band, name in enumerate(("uth", "thk", "uk"), 1):
            dst.write(np.full((h, w), values[name], np.uint8), band)
            dst.set_band_description(band, name.upper())

    footprint = np.ones((h, w), bool)
    footprint[0, 0] = False
    loaded = read_ratio_bands(ratio_path, footprint, reference_path=ref)

    assert set(loaded) == set(RATIO_BANDS)
    assert loaded["thk"][2, 2] == pytest.approx(64 / 255.0)
    assert loaded["uk"][2, 2] == pytest.approx(128 / 255.0)
    assert loaded["uth"][2, 2] == pytest.approx(1.0)
    assert all(layer[0, 0] == 0.0 for layer in loaded.values())


def test_conductivity_reader_keeps_raw_values_and_fills_nodata(tmp_path):
    h, w = 6, 7
    transform = from_origin(100, 600, 100, 100)
    ref = tmp_path / "reference.tif"
    _write_one_band(
        ref,
        np.ones((h, w), np.float32),
        dtype="float32",
        transform=transform,
        description="template",
    )
    cond_path = tmp_path / "features.tif"
    data = np.arange(h * w, dtype=np.float32).reshape(h, w) + 100.0
    data[2, 3] = -9999.0
    _write_one_band(
        cond_path,
        data,
        dtype="float32",
        transform=transform,
        description="cond_surf - Surface conductivity",
        nodata=-9999.0,
    )
    footprint = np.ones((h, w), bool)
    footprint[0, 0] = False

    loaded = read_conductivity(cond_path, footprint, reference_path=ref)

    assert loaded[1, 1] == pytest.approx(data[1, 1])
    assert loaded[2, 3] != -9999.0
    assert np.isfinite(loaded[footprint]).all()
    assert loaded[0, 0] == 0.0


def test_ratio_reader_rejects_grid_mismatch(tmp_path):
    h, w = 4, 5
    ref = tmp_path / "reference.tif"
    ratio_path = tmp_path / "ratios.tif"
    _write_one_band(
        ref,
        np.ones((h, w), np.float32),
        dtype="float32",
        transform=from_origin(0, 1000, 100, 100),
        description="template",
    )
    with rasterio.open(
        ratio_path,
        "w",
        driver="GTiff",
        width=w,
        height=h,
        count=3,
        dtype="uint8",
        crs="EPSG:32611",
        transform=from_origin(100, 1000, 100, 100),
    ) as dst:
        for band, name in enumerate(("ThK", "UK", "UTh"), 1):
            dst.write(np.ones((h, w), np.uint8), band)
            dst.set_band_description(band, name)

    with pytest.raises(ValueError, match="grid mismatch"):
        read_ratio_bands(
            ratio_path,
            np.ones((h, w), bool),
            reference_path=ref,
        )
