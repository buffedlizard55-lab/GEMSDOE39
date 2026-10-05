"""Tests for the session-2 science module and the shared emitter.

Each test names the code path it executes so the check is traceable:

* ``emit_sequence``/``emit_prefix`` must reproduce ``emit_fast`` exactly --
  the whole budget sweep depends on that equivalence.
* ``rel_density_profile`` must be a density per unit area, not a raw count, and
  must zero the bins inside the catalogue exclusion.
* ``profile_weights`` must never authorise a dot inside the exclusion radius.
* ``calibrated_score`` must be monotone in per-dot credit and must saturate in
  the dot count (the property that lets it explain rho(n, score) = -0.86).
* ``expected_credit`` must be 1.0 for dots sitting on truth and 0.0 for dots
  farther than the 300 m kernel radius.
* ``hide_components`` must hide whole components and must not leak the hidden
  pixels into the visible catalogue that the candidate is allowed to see.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gems39 import emission as em  # noqa: E402
from gems39 import h40  # noqa: E402


def _toy(seed=0, shape=(400, 400)):
    """A tiny deterministic context so the tests do not need the 419 MB raster."""
    rng = np.random.default_rng(seed)
    H, W = shape
    foot = np.zeros(shape, bool)
    foot[20:H - 20, 20:W - 20] = True
    cat = np.zeros(shape, bool)
    yy = np.arange(60, H - 60, 37)
    for y in yy:
        cat[y, 40:W - 40] = True
    off = np.zeros(shape, bool)
    fy, fx = np.nonzero(foot & ~cat)
    idx = rng.integers(0, fy.size, 4000)
    off[fy[idx], fx[idx]] = True
    from scipy.ndimage import distance_transform_edt

    return dict(foot=foot, cat=cat, sgmc=off | cat, off=off,
                dcat=distance_transform_edt(~cat), blocks=h40.make_blocks(foot, 2, 2),
                shape=shape)


# --------------------------------------------------------------------- emitter
def test_emit_sequence_prefix_equals_emit_fast():
    rng = np.random.default_rng(7)
    shape = (300, 300)
    foot = np.zeros(shape, bool)
    foot[10:290, 10:290] = True
    field = np.where(foot, rng.random(shape), 0.0).astype(np.float32)
    cat = np.zeros(shape, bool)
    cat[150, 10:290] = True
    for budget in (200, 900, 2500):
        fast = em.emit_fast(field, foot, budget, min_dist=3.0, catalogue=cat,
                            cat_buffer_px=2)
        seq = em.emit_sequence(field, foot, max_n=budget, min_dist=3.0,
                               catalogue=cat, cat_buffer_px=2)
        pref = em.emit_prefix(seq, budget, shape)
        assert fast.sum() == pref.sum() == budget, (budget, int(fast.sum()), int(pref.sum()))
        assert np.array_equal(fast, pref)


def test_emit_sequence_respects_min_dist_and_catalogue_buffer():
    shape = (200, 200)
    foot = np.ones(shape, bool)
    field = np.ones(shape, np.float32)
    cat = np.zeros(shape, bool)
    cat[100, 50:150] = True
    seq = em.emit_sequence(field, foot, max_n=400, min_dist=3.0, catalogue=cat,
                           cat_buffer_px=2)
    ys, xs = seq[:, 0], seq[:, 1]
    from scipy.spatial import cKDTree
    from scipy.ndimage import distance_transform_edt

    d, _ = cKDTree(seq.astype(float)).query(seq.astype(float), k=2)
    assert d[:, 1].min() >= 3.0 - 1e-9
    assert distance_transform_edt(~cat)[ys, xs].min() > 2.0


# -------------------------------------------------------------------- profile
def test_rel_density_profile_is_area_normalised_and_respects_exclusion():
    ctx = _toy()
    prof = h40.rel_density_profile(ctx, exclude_px=2.0)
    for b in prof["bins"]:
        if b["area_px"]:
            assert b["density"] == pytest.approx(b["off_px"] / b["area_px"])
        if b["lo_m"] < 200.0:
            assert b["usable"] is False
    # the weights derived from it must be exactly zero inside the exclusion
    w = h40.profile_weights(prof, ctx["shape"], ctx["dcat"], cap=1.25, floor=0.15)
    assert (w[ctx["dcat"] <= 2.0] == 0).all()
    assert (w[ctx["dcat"] > 2.0][ctx["foot"][ctx["dcat"] > 2.0]] > 0).any()
    assert w.max() <= 1.25 + 1e-6
    assert w[w > 0].min() >= 0.15 - 1e-6


def test_rel_density_profile_mask_restricts_the_measurement():
    ctx = _toy()
    full = h40.rel_density_profile(ctx)
    half = h40.rel_density_profile(ctx, mask=(ctx["blocks"] != 0))
    tot_full = sum(b["area_px"] for b in full["bins"])
    tot_half = sum(b["area_px"] for b in half["bins"])
    assert tot_half < tot_full


# ------------------------------------------------------------------ instrument
def test_expected_credit_bounds():
    ctx = _toy()
    truth = ctx["off"]
    on = truth.copy()
    w_on, n = h40.expected_credit(on, truth)
    assert w_on == pytest.approx(1.0)
    assert n == int(truth.sum())
    far = np.zeros(ctx["shape"], bool)
    from scipy.ndimage import distance_transform_edt

    far[distance_transform_edt(~truth) > 3.0] = True
    w_far, _ = h40.expected_credit(far, truth)
    assert w_far == 0.0
    empty, n0 = h40.expected_credit(np.zeros(ctx["shape"], bool), truth)
    assert empty == 0.0 and n0 == 0


def test_calibrated_score_monotone_in_credit_and_saturating_in_budget():
    K, a = 5365.0, 2.18
    lo = h40.calibrated_score(40000, 0.04, K, a)
    hi = h40.calibrated_score(40000, 0.06, K, a)
    assert hi > lo
    # saturating: at a fixed per-dot credit the score must eventually fall as the
    # budget grows (this is what reproduces rho(n, live score) = -0.86)
    s = [float(h40.calibrated_score(n, 0.056, K, a))
         for n in (20000, 40000, 80000, 200000, 500000)]
    assert s.index(max(s)) < len(s) - 1
    assert s[-1] < max(s)
    # TP is capped at K, so the score can never exceed K / (0.2K + 0.8K) = 1.0
    assert float(h40.calibrated_score(10_000_000, 1.0, K, a)) <= 1.0


def test_fit_calibrated_recovers_a_saturating_generative_process():
    K_true, a_true = 6000.0, 2.5
    n = np.array([20000, 30000, 45000, 60000, 90000, 130000, 180000, 250000], float)
    w = np.array([0.050, 0.048, 0.047, 0.046, 0.044, 0.043, 0.042, 0.041])
    s = np.asarray(h40.calibrated_score(n, w, K_true, a_true), float)
    fit = h40.fit_calibrated(n, w, s)
    assert fit["rho_in_sample"] > 0.99
    assert abs(fit["K"] - K_true) / K_true < 0.25
    assert fit["rmse"] < 0.005


# -------------------------------------------------------------------- holdout
def test_hide_components_hides_whole_components_without_leakage():
    ctx = _toy()
    comps = h40.label_components(ctx["cat"])
    truth, visible, domain = h40.hide_components(ctx, 0, hide_frac=0.30, seed=1,
                                                 components=comps)
    assert truth.any()
    assert not (truth & visible).any()
    # every hidden pixel belongs to a *globally* labelled component that is
    # wholly hidden -- no half-trace left visible across a block boundary
    hidden_ids = set(np.unique(comps[truth])) - {0}
    assert hidden_ids
    for i in hidden_ids:
        assert not (visible & (comps == i)).any()
        # and the whole component is in the hidden set, not just its inner part
        assert ((comps == i) & ~truth).sum() == 0
    # the visible catalogue is never part of the evaluation domain
    assert not (domain & visible).any()
    assert (truth & domain).sum() > 0


def test_hide_components_covers_every_block_exactly_once():
    ctx = _toy()
    comps = h40.label_components(ctx["cat"])
    seen = np.zeros(ctx["shape"], bool)
    for b in range(4):
        t, _, _ = h40.hide_components(ctx, b, hide_frac=0.25, seed=b, components=comps)
        assert not (t & seen).any(), "a component was hidden in two blocks"
        seen |= t
    assert seen.sum() > 0


def test_rank_lift_is_0p5_for_a_constant_field():
    ctx = _toy()
    truth, visible, domain = h40.hide_components(ctx, 1, hide_frac=0.30, seed=3)
    const = np.ones(ctx["shape"], np.float32)
    assert h40.rank_lift(const, truth, domain) == pytest.approx(0.5, abs=0.02)


def test_rank_lift_orders_a_perfect_field_above_an_antiperfect_one():
    ctx = _toy()
    comps = h40.label_components(ctx["cat"])
    truth, visible, domain = h40.hide_components(ctx, 2, hide_frac=0.30, seed=5,
                                                 components=comps)
    good = np.where(truth, 1.0, 0.0).astype(np.float32)
    bad = np.where(truth, 0.0, 1.0).astype(np.float32)
    lg = h40.rank_lift(good, truth, domain)
    lb = h40.rank_lift(bad, truth, domain)
    assert lg > 0.9 and lb < 0.1


# ------------------------------------------------------- real-raster regression
_REAL = ROOT / "data" / "training_features.tif"


@pytest.mark.skipif(not _REAL.is_file(), reason="competition raster not restored")
def test_read_bands_returns_finite_clipped_values():
    """Regression: the nearest-fill used EDT of the *valid* mask, which returns
    the index of the nearest invalid pixel and turned every band into NaN.
    """
    with rasterio.open(ROOT / "data" / "labels.tif") as src:
        foot = src.read(1) >= 0
    B = h40._read_bands(str(_REAL), ["det_elev", "tmi", "cond_surf"], foot)
    for nm, v in B.items():
        assert v.shape == foot.shape
        assert np.isfinite(v[foot]).all(), f"{nm} has non-finite cells inside the footprint"
        assert v[~foot].max() == 0.0
        assert v[foot].std() > 0.0, f"{nm} is degenerate"


@pytest.mark.skipif(not _REAL.is_file(), reason="competition raster not restored")
def test_structural_field_is_not_degenerate():
    with rasterio.open(ROOT / "data" / "labels.tif") as src:
        foot = src.read(1) >= 0
    f = h40.structural_field(ROOT / "data", foot)
    assert np.isfinite(f).all()
    assert f[foot].std() > 1e-3
    assert (f[foot] > 0).mean() > 0.05
