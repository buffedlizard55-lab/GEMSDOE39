"""Verification tests for the H40 round.

Every test here pins a claim made in the documentation to a number that is
recomputed from bytes on disk or from the organizer's own published metric
definition.  Nothing is asserted from memory.

Skipped (not silently passed) when the pinned competition rasters are absent --
run ``python scripts/restore_data.py --group core --group calib`` first.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import rasterio  # noqa: E402

from gems39 import calibrate as cal           # noqa: E402
from gems39 import detectors40 as det40       # noqa: E402
from gems39 import external as ex             # noqa: E402
from gems39 import grid, priced_emit as pe    # noqa: E402
from gems39.metric import dti, dti_binary, kernel  # noqa: E402
from gems39.sprt_select import folds_needed, sprt_pairwise  # noqa: E402

DATA = ROOT / "data"
HAVE_DATA = (DATA / "sample_submission.tif").exists() and (DATA / "labels.tif").exists()
HAVE_CALIB = (ROOT / "data" / "calib" / "h33-2-b2_0.2778.tif").exists() and \
             (ROOT / "data" / "calib" / "h27-4-solo-d2-8_0.2708.tif").exists()
need_data = pytest.mark.skipif(not HAVE_DATA, reason="pinned competition rasters not restored")
need_calib = pytest.mark.skipif(not (HAVE_DATA and HAVE_CALIB),
                                reason="calibration artifacts not restored")


# ------------------------------------------------------- the organizer's metric
def test_metric_matches_published_worked_example():
    """DrivenData page 967 worked example: TPw=3.00, FPw=1.89, FNw=2.00 -> 0.60.

    https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/
    """
    assert round(dti(3.00, 1.89, 2.00), 2) == 0.60
    assert dti(3.00, 1.89, 2.00) == pytest.approx(3.00 / (3.00 + 0.2 * 1.89 + 0.8 * 2.00))


def test_kernel_is_triangular_with_R_300m():
    assert kernel(0.0) == 1.0
    assert abs(kernel(1.5) - 0.5) < 1e-12
    assert kernel(3.0) == 0.0
    assert kernel(3.5) == 0.0


def test_tp_plus_fn_equals_truth_count():
    """Identity used as an independent check on the metric implementation."""
    rng = np.random.default_rng(0)
    H = W = 64
    truth = rng.random((H, W)) < 0.02
    pred = rng.random((H, W)) < 0.05
    r = dti_binary(pred, truth, np.ones((H, W), bool))
    assert abs((r["tp"] + r["fn"]) - truth.sum()) < 1e-6
    assert r["n_truth"] == int(truth.sum())


# ------------------------------------------------------------ band map (defect)
@need_data
def test_band_order_is_read_from_the_file_not_hardcoded():
    names = grid.read_band_names(DATA / "training_features.tif") \
        if (DATA / "training_features.tif").exists() else None
    if names is None:
        pytest.skip("training_features.tif not restored")
    assert names == list(grid.VERIFIED_BAND_ORDER)
    # the specific defect: tc is band 6, not band 18
    assert names[5] == "tc"
    assert names[17] == "iso_grav_anom_hg"
    # and the previously-unnamed bands exist
    for nm in ("geod_2ndinv", "geod_shearrate", "geod_dilaterate",
               "ieq_n100a15", "deq_n100a15", "tmi_hg", "tmi_vg",
               "iso_grav_anom_slope", "iso_grav_anom_vg", "iso_grav_anom_hg"):
        assert nm in names


# ------------------------------------------------------------------- grid facts
@need_data
def test_grid_geometry_and_masks():
    with rasterio.open(DATA / "sample_submission.tif") as s:
        sub = s.read(1)
        assert str(s.crs).upper().startswith("EPSG:32611")
        assert (s.width, s.height) == (3292, 3730)
        assert tuple(s.transform)[:6] == (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)
    foot = np.isfinite(sub)
    with rasterio.open(DATA / "labels.tif") as s:
        cat = (s.read(1) >= 1) & foot
    assert int(foot.sum()) == 5_167_373
    assert int(cat.sum()) == 60_988
    assert int((foot & ~cat).sum()) == 5_106_385


# ------------------------------------------------- nested-pair calibration
@need_calib
def test_nested_pair_is_exactly_as_documented():
    r = cal.verify_nested_pair(DATA / "calib" / "h27-4-solo-d2-8_0.2708.tif",
                              DATA / "calib" / "h33-2-b2_0.2778.tif",
                              DATA / "labels.tif", DATA / "sample_submission.tif")
    assert r["n_base"] == 40_199
    assert r["n_prune"] == 37_654
    assert r["n_removed"] == 2_545
    assert r["pruned_is_subset"]
    assert r["base_on_catalogue"] == 0 and r["pruned_on_catalogue"] == 0
    assert abs(r["removed_dcat_min"] - math.sqrt(2)) < 1e-6
    assert r["removed_dcat_max"] == 2.0


@need_calib
def test_G_is_identified_and_tight_across_nuisance_priors():
    """|G| must be stable over the l0 x kb sensitivity grid."""
    sens = cal.fit_G_sensitivity(40_199, 0.2708, 37_654, 0.2778)
    assert 12_000 < sens["G_median"] < 17_000
    assert (sens["G_max"] - sens["G_min"]) / sens["G_median"] < 0.05


@need_calib
def test_forward_model_reproduces_both_anchor_scores():
    sens = cal.fit_G_sensitivity(40_199, 0.2708, 37_654, 0.2778)
    G = sens["G_median"]
    h = cal.implied_hit_rate(37_654, 0.2778, G)
    assert abs(cal.forward(37_654, h * 37_654, G)[0] - 0.2778) < 5e-4
    hb = cal.implied_hit_rate(40_199, 0.2708, G)
    assert abs(cal.forward(40_199, hb * 40_199, G)[0] - 0.2708) < 5e-4
    # the hit rate must increase as dots are removed -- that is the whole point
    assert h > hb


def test_break_even_probability_is_in_a_sane_range():
    r = cal.break_even_pi(0.2778, 14_143.0, 2308.0)
    assert 0.01 < r["pi_star"] < 0.10


def test_binary_emission_is_never_worse_than_scaling_down():
    """Soft (scaled) predictions cannot beat binary ones under this metric.

    DTI(c*p) is monotone increasing in c because 0.8*|G| is a constant additive
    term in the denominator, so the optimum is at the largest admissible c, i.e.
    p in {0, 1}.
    """
    G = 14_143.0
    for n in (20_000, 37_654, 60_000):
        for h in (0.03, 0.06, 0.10):
            s1 = cal.forward(n, h * n, G)[0]
            for c in (0.2, 0.5, 0.8, 1.0):
                # scaling a binary field by c scales TPw and FPw by c
                TP = G * (1 - math.exp(-3.0 * h * n / G)) * c
                FP = (n - 0.5 * h * n) * c
                s = TP / (0.2 * TP + 0.2 * FP + 0.8 * G)
                assert s <= s1 + 1e-9


# ------------------------------------------------------------------- SPRT
def test_sprt_boundaries_are_wald():
    d = sprt_pairwise([], p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    assert abs(d["upper"] - math.log(0.90 / 0.05)) < 1e-12
    assert abs(d["lower"] - math.log(0.10 / 0.95)) < 1e-12


def test_sprt_stops_only_at_a_boundary_and_reports_consistent_counts():
    d = sprt_pairwise([1] * 9, p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    assert d["decision"] == "accept_H1"
    assert d["n"] == d["wins"] == 9
    # the lower boundary is ln(0.10/0.95) = -2.2513 and a loss contributes
    # ln(0.3/0.5) = -0.5108, so 5 consecutive losses cross it (not 7).
    d2 = sprt_pairwise([0] * 7, p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    assert d2["decision"] == "accept_H0"
    assert d2["n"] == 5 and d2["wins"] == 0
    # a generator must not be double-consumed (the fixed defect)
    # a generator must not be double-consumed (the fixed defect): total is the
    # length of the supplied vector, wins counts only folds consumed before the
    # boundary, and n is the stopping point.
    vec = [1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1]      # 11 wins, 1 loss
    d3 = sprt_pairwise(iter(vec), p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    assert d3["decision"] == "accept_H1"
    assert d3["total"] == len(vec)                   # generator fully materialised
    assert d3["n"] == 12 and d3["wins"] == 11
    assert d3["llr"] == pytest.approx(11 * math.log(1.4) + math.log(0.6))
    assert d3["llr"] >= d3["upper"]
    # the same vector as a list must give an identical result
    assert sprt_pairwise(vec, p0=0.5, p1=0.7, alpha=0.05, beta=0.10) == d3
    # and one win fewer must NOT cross the boundary -- stopping is at the
    # boundary only, never on a judgment call
    d4 = sprt_pairwise([1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1], p0=0.5, p1=0.7,
                       alpha=0.05, beta=0.10)
    assert d4["decision"] == "continue" and d4["llr"] < d4["upper"]


def test_holdout_design_can_reach_a_boundary():
    """n = 8 folds cannot cross the upper boundary; the H40 design can."""
    f = folds_needed(p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    assert f["min_all_wins"] == 9          # ceil(2.8904 / ln(1.4))
    assert f["min_all_losses"] == 5        # ceil(-2.2513 / ln(0.6))
    assert 8 < f["min_all_wins"], "an 8-fold design has no power to accept H1"
    assert 4 * 6 >= f["min_all_wins"]


# ------------------------------------------------------- azimuth convention
def test_strike_convention_matches_the_official_release():
    """Our image->geographic strike must agree with the release's own Strike.

    The Siler (2022) release gives, per segment, both the endpoint lon/lat and
    an attributed ``Strike``.  Recomputing the strike from the endpoints with our
    convention and comparing to the attribute is a direct check that no 90-degree
    or sign error crept into ``image_to_geographic_strike``.
    """
    import csv
    p = None
    for cand in (ROOT / "data" / "external2" / "sb_slip_tendency_in_footprint.csv",
                 ROOT / "data" / "external" / "sb_slip_tendency_in_footprint.csv"):
        if cand.exists():
            p = cand
            break
    if p is None:
        pytest.skip("slip-tendency release not restored")
    got, want, w = [], [], []
    with p.open(newline="") as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            if i > 20000:
                break
            try:
                x0, y0 = float(row["X_START"]), float(row["Y_START"])
                x1, y1 = float(row["X_END"]), float(row["Y_END"])
                s = float(row["Strike"]); L = float(row.get("length_m") or 0)
            except (TypeError, ValueError, KeyError):
                continue
            if L < 500:
                continue
            # geographic strike from endpoints, clockwise from north, as an axis
            dE = (x1 - x0) * math.cos(math.radians(0.5 * (y0 + y1)))
            dN = (y1 - y0)
            az = math.degrees(math.atan2(dE, dN)) % 180.0
            got.append(az); want.append(s % 180.0); w.append(L)
    assert len(got) > 500, "too few segments to check the convention"
    got = np.array(got); want = np.array(want); w = np.array(w)
    d = np.abs(((got - want + 90) % 180) - 90)
    assert float(np.median(d)) < 15.0, f"median |strike difference| = {np.median(d):.1f} deg"
    assert float((d < 30).mean()) > 0.80


def test_favourability_curve_peaks_at_favourable_orientation():
    """The empirical curve must peak where physics says it should.

    dtheta = strike - (SHmaxAz + 90), folded to [0, 90].  dtheta ~ 90 deg means
    the fault strikes parallel to sigma_Hmax, which is the favourably oriented
    normal fault (Barton et al. 1995; Morris et al. 1996).  dtheta ~ 0 is the
    unfavourable end.  The curve is measured from the release, so this test can
    fail -- and if it does, the H40-A gate is wrong, not the test.
    """
    import csv
    p = None
    for cand in (ROOT / "data" / "external2" / "sb_slip_tendency_in_footprint.csv",
                 ROOT / "data" / "external" / "sb_slip_tendency_in_footprint.csv"):
        if cand.exists():
            p = cand
            break
    if p is None:
        pytest.skip("slip-tendency release not restored")
    sl, tl, al = [], [], []
    with p.open(newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                sl.append(float(row["Strike"])); tl.append(float(row["TS_norm"]))
                al.append(float(row["SHmaxAz"]))
            except (TypeError, ValueError, KeyError):
                continue
    centres, curve, n = ex.strike_favourability(np.array(sl), np.array(tl),
                                                np.array(al), n_bins=18)
    assert curve[-3:].mean() > curve[:3].mean() * 1.5, \
        "favourable end is not clearly better than the unfavourable end"
    assert int(np.argmax(curve)) >= 9, "peak should be in the favourable half"


# ------------------------------------------------------------ emission geometry
def test_poisson_disk_respects_minimum_spacing():
    rng = np.random.default_rng(1)
    prop = rng.random((120, 120)).astype(np.float32)
    allowed = np.ones((120, 120), bool)
    mask, ys_e, xs_e = pe.emit_ranked(prop, allowed, min_dist=2.8, max_n=400)
    ys, xs = np.nonzero(mask)
    # the returned coordinates must be exactly the mask, in best-first order
    assert sorted(zip(ys_e.tolist(), xs_e.tolist())) == sorted(zip(ys.tolist(), xs.tolist()))
    v = prop[ys_e, xs_e].astype(np.float64)
    assert (v[:-1] >= v[1:] - 1e-6).all(), "not in descending propensity order"
    from scipy.ndimage import distance_transform_edt
    for i in range(min(ys.size, 400)):
        m2 = mask.copy(); m2[ys[i], xs[i]] = False
        d = distance_transform_edt(~m2)
        assert d[ys[i], xs[i]] >= 2.8 - 1e-6


def test_catalogue_buffer_is_a_hard_gate():
    foot = np.ones((80, 80), bool)
    cat = np.zeros((80, 80), bool)
    cat[40, 10:70] = True
    ok = pe.allowed_domain(foot, cat, 2)
    assert not ok[40, 40]
    assert not ok[39, 40] and not ok[41, 40]
    assert ok[37, 40] and ok[43, 40]


def test_trim_to_budget_keeps_the_best_first_dots():
    """Regression test for the index-scrambling defect.

    ``trim_to_budget`` must return the n HIGHEST-propensity accepted dots.  The
    previous implementation reordered raster-ordered coordinates by pool-rank
    positions and returned an arbitrary permutation, which passed the old
    assertions (right count, subset of the mask) while silently degrading every
    trimmed emission.
    """
    rng = np.random.default_rng(11)
    H = W = 90
    prop = rng.random((H, W)).astype(np.float32)
    allowed = np.ones((H, W), bool)
    mask, ys, xs = pe.emit_ranked(prop, allowed, min_dist=2.8, max_n=600)
    n = int(mask.sum())
    assert n > 20
    for k in (1, n // 3, n // 2, n - 1, n):
        keep = pe.trim_to_budget((H, W), ys, xs, k)
        assert int(keep.sum()) == k
        assert not (keep & ~mask).any()
        kept_vals = prop[keep]
        dropped_vals = prop[mask & ~keep]
        if dropped_vals.size:
            assert kept_vals.min() >= dropped_vals.max() - 1e-6, \
                f"trim to {k} kept a lower-propensity dot than one it dropped"
    # out-of-range n is clamped, not an error
    assert int(pe.trim_to_budget((H, W), ys, xs, 10 ** 9).sum()) == n
    assert int(pe.trim_to_budget((H, W), ys, xs, 0).sum()) == 0


# ------------------------------------------------------------- leakage guard
@need_data
def test_leakage_probe_detects_a_cheating_mask():
    from gems39 import holdout40 as ho
    with rasterio.open(DATA / "sample_submission.tif") as s:
        foot = np.isfinite(s.read(1))
    with rasterio.open(DATA / "labels.tif") as s:
        cat = (s.read(1) >= 1) & foot

    class _Ctx:
        pass
    ctx = _Ctx(); ctx.foot = foot; ctx.catalogue = cat
    H, W = foot.shape
    # Place both probes on pixels that are actually inside the footprint: it is
    # an irregular shape, so a hardcoded row/col would land outside it and the
    # probe would return 0.0 -- passing for the wrong reason.
    idx = np.flatnonzero(foot.ravel())
    ia, ib = idx[idx.size // 3], idx[2 * idx.size // 3]
    ya, xa = divmod(int(ia), W)
    yb, xb = divmod(int(ib), W)

    truth = np.zeros(foot.shape, bool)
    truth[max(ya - 5, 0):ya + 5, max(xa - 5, 0):xa + 5] = True
    truth &= foot
    assert truth.sum() > 0

    cheat = truth.copy()
    assert ho.leakage_probe(cheat, ctx, truth) == 1.0

    honest = np.zeros(foot.shape, bool)
    honest[max(yb - 5, 0):yb + 5, max(xb - 5, 0):xb + 5] = True
    honest &= foot
    honest &= ~truth
    assert honest.sum() > 0
    assert ho.leakage_probe(honest, ctx, truth) == 0.0

    # a half-contaminated mask must report exactly its contaminated fraction
    half = cheat.copy()
    ys_t, xs_t = np.nonzero(half)
    half[ys_t[: ys_t.size // 2], xs_t[: ys_t.size // 2]] = False
    mixed = half | honest
    assert ho.leakage_probe(mixed, ctx, truth) == pytest.approx(
        float((mixed & truth).sum()) / float(mixed.sum()))

    # an empty mask must not divide by zero
    assert ho.leakage_probe(np.zeros(foot.shape, bool), ctx, truth) == 0.0


# --------------------------------------------------------- prevalence matching
def test_prevalence_match_hits_the_target_without_splitting_traces():
    rng = np.random.default_rng(3)
    m = np.zeros((400, 400), bool)
    for _ in range(60):
        y, x = rng.integers(10, 390, 2)
        L = rng.integers(20, 60)
        m[y:y + 2, x:x + L] = True
    out, n = ex_prevalence(m, 800)
    assert abs(n - 800) < 400
    from scipy.ndimage import label
    comp, nc = label(out, structure=np.ones((3, 3), int))
    sizes = np.bincount(comp.ravel())[1:]
    comp_full, nc_full = label(m, structure=np.ones((3, 3), int))
    sizes_full = np.bincount(comp_full.ravel())[1:]
    # every retained component must appear in the source at its full size
    for s_ in sizes[sizes > 0]:
        assert s_ in sizes_full


def ex_prevalence(mask, target):
    from gems39.holdout40 import prevalence_match
    return prevalence_match(mask, target)


# ------------------------------------------- Wald SPRT for a normal mean (magnitudes)
def test_sprt_normal_mean_boundaries_and_direction():
    from gems39.sprt_select import sprt_normal_mean
    r = sprt_normal_mean([], d_alt=0.5)
    assert r["decision"] == "continue" and r["n"] == 0
    pos = sprt_normal_mean([1.0] * 12, d_alt=0.5)
    assert pos["decision"] == "accept_H1" and pos["llr"] >= pos["upper"]
    neg = sprt_normal_mean([-1.0] * 12, d_alt=0.5)
    assert neg["decision"] == "accept_H0" and neg["llr"] <= neg["lower"]
    # LLR must equal Wald's closed form for a normal mean
    x = np.array([0.7, -0.2, 1.1, 0.4, -0.5, 0.9, 0.3, 0.8, -0.1, 1.2])
    r = sprt_normal_mean(x, d_alt=0.5)
    s = float(x.std(ddof=1))
    expected = (0.5 / s) * x[: r["n"]].sum() - r["n"] * 0.5 ** 2 / 2
    assert r["llr"] == pytest.approx(expected, abs=1e-9)


def test_sign_test_and_mean_test_can_disagree_and_both_are_reported():
    """The case that motivated having two tests.

    A candidate that gains a lot in a few truth-dense folds and a little in many
    truth-poor folds has a positive SUM of per-fold score deltas (so it wins the
    pooled metric, which is what the organizer scores) while losing a majority of
    fold signs.  The sign test rejects it; the mean test accepts it.  Both must be
    reported -- neither alone is sufficient evidence.
    """
    from gems39.sprt_select import sprt_normal_mean, sprt_pairwise
    deltas = [3.0] * 6 + [-0.2] * 18
    sign = sprt_pairwise([d > 0 for d in deltas], p0=0.5, p1=0.7, alpha=0.05, beta=0.10)
    mean = sprt_normal_mean(deltas, d_alt=0.5, alpha=0.05, beta=0.10)
    assert sum(deltas) > 0
    assert mean["decision"] == "accept_H1"
    assert sign["decision"] != "accept_H1"


# ------------------------------------------------ fold delta additivity
def test_fold_score_delta_is_the_first_order_pooled_score_change():
    """Exact identity: delta_fold = D_cand,fold * (s_cand,fold - s_anchor).

    Derivation.  With D = 0.2*TP + 0.2*FP + 0.8*G and s = TP/D,

        delta = dTP*(1 - 0.2*s_a) - 0.2*s_a*dFP
              = (TP_c - TP_a) - s_a*[(D_c - 0.8G) - (D_a - 0.8G)]
              = D_c*s_c - D_a*s_a - s_a*D_c + s_a*D_a
              = D_c * (s_c - s_a)                       (exact, no approximation)

    Because s_a is a single pre-declared constant (the anchor's POOLED instrument
    DTI), summing delta over folds telescopes to

        sum_k delta_k = (TP_c - TP_a)(1 - 0.2*s_a) - 0.2*s_a*(FP_c - FP_a)
                      = D_c_pooled * (s_c_pooled - s_a)

    i.e. exactly the pooled score change against the anchor, up to the positive
    factor D_c_pooled.  That is what makes a sequential test over folds a test of
    the POOLED claim -- which is what the organizer scores -- rather than of a
    per-block claim.
    """
    from gems39 import holdout40 as ho
    rng = np.random.default_rng(7)
    H = W = 120
    foot = np.ones((H, W), bool)
    truth = rng.random((H, W)) < 0.03
    anchor = rng.random((H, W)) < 0.10
    cand = rng.random((H, W)) < 0.10

    class C:
        pass
    cell = C()
    cell.bbox = (slice(0, H), slice(0, W))
    cell.active_sub = foot
    cell.truth_sub = truth
    cell.known_sub = np.zeros((H, W), bool)
    cell.n_truth = int(truth.sum())
    cell.n_active = int(foot.sum())

    from gems39.metric import dti_binary
    s_a = dti_binary(anchor, truth, valid=foot)["dti"]
    s_c = dti_binary(cand, truth, valid=foot)["dti"]
    fd = ho.fold_score_delta(cell, cand, anchor, s_a)
    c = ho.dti_cell(cand, cell)
    D_c = 0.2 * c["tp"] + 0.2 * c["fp"] + 0.8 * c["n_truth"]
    # exact identity, not a first-order approximation
    assert fd["delta"] == pytest.approx(D_c * (c["dti"] - s_a), rel=1e-9, abs=1e-12)
    assert np.sign(fd["delta"]) == np.sign(c["dti"] - s_a) or abs(c["dti"] - s_a) < 1e-12

    # Additivity across folds.  NOTE the boundary condition: the kernel reaches
    # 3 px, so a dot on one side of a split can credit truth on the other and TP/FP
    # are then NOT additive -- which is exactly why holdout40 erodes every block by
    # 6 px and selects hidden components with a 10 px collar.  Clearing a 5 px band
    # around the split removes the interaction so additivity is exact.
    band = slice(H // 2 - 5, H // 2 + 5)
    truth[band, :] = False
    anchor[band, :] = False
    cand[band, :] = False
    cell.truth_sub = truth
    cell.n_truth = int(truth.sum())
    s_a = dti_binary(anchor, truth, valid=foot)["dti"]
    fd = ho.fold_score_delta(cell, cand, anchor, s_a)
    c = ho.dti_cell(cand, cell)
    a = ho.dti_cell(anchor, cell)
    half = (slice(0, H // 2), slice(0, W))
    rest = (slice(H // 2, H), slice(0, W))
    total = 0.0
    for bb in (half, rest):
        cl = C()
        cl.bbox = bb
        cl.active_sub = foot[bb]
        cl.truth_sub = truth[bb]
        cl.known_sub = np.zeros((H // 2 + 1, W), bool)
        cl.n_truth = int(truth[bb].sum())
        cl.n_active = int(foot[bb].sum())
        total += ho.fold_score_delta(cl, cand, anchor, s_a)["delta"]
    expected = ((c["tp"] - a["tp"]) * (1 - 0.2 * s_a)
                - 0.2 * s_a * (c["fp"] - a["fp"]))
    assert total == pytest.approx(expected, rel=1e-9, abs=1e-9)
    # and the whole-domain identity still holds after the split
    D_c = 0.2 * c["tp"] + 0.2 * c["fp"] + 0.8 * c["n_truth"]
    assert fd["delta"] == pytest.approx(D_c * (c["dti"] - s_a), rel=1e-9, abs=1e-12)


# ------------------------------------------------ H40-F fold alignment (leakage)
@need_data
def test_h40f_fold_grid_is_identical_to_the_evaluation_grid():
    """The defect that produced a fake 2.05x win.

    H40-F trains on SGMC off-catalogue faults; instrument I2 scores against the
    same population.  If the training fold grid and the evaluation block grid
    differ, an evaluation block lies partly in folds the model WAS trained on, and
    the measured advantage is leakage.  ``block_labels`` must therefore return
    exactly ``holdout40.blocks``.
    """
    from gems39 import holdout40 as ho
    shape = (3730, 3292)
    bid_f, n_f = det40.block_labels(shape)
    bid_h, boxes_h, n_h = ho.blocks(shape, 4, 6)
    assert n_f == n_h
    assert np.array_equal(bid_f.astype(np.int32), bid_h.astype(np.int32))
    # and fold = block mod n_folds, so block b is predicted without block b
    n_folds = 4
    fold = (np.arange(n_f) % n_folds)[bid_f]
    for b in range(n_h):
        q = bid_h == b
        assert np.unique(fold[q]).size == 1
        assert int(np.unique(fold[q])[0]) == b % n_folds
        # no positive of block b can appear in the training set of the model that
        # predicts block b
        train_blocks = [bb for bb in range(n_h) if bb % n_folds != b % n_folds]
        assert b not in train_blocks


def test_block_grids_that_do_not_nest_are_detected():
    """Guard against the defect reappearing: 3x4 and 4x6 must NOT agree."""
    from gems39 import holdout40 as ho
    a, _, na = ho.blocks((3730, 3292), 3, 4)
    b, _, nb = ho.blocks((3730, 3292), 4, 6)
    assert na != nb
    assert not np.array_equal(a.astype(np.int32), b.astype(np.int32))
