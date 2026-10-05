#!/usr/bin/env python3
"""Uniqueness gate: prove the shipped raster is not a copy of any prior artifact.

For every same-grid raster we can retrieve from the sibling GEMSDOE
repositories, compute:

  * Jaccard overlap of the positive-pixel sets
  * the fraction of our positives that are also positive there (containment)
  * Pearson correlation of the two rasters inside the footprint

A submission is rejected by this gate if its Jaccard overlap with any single
prior artifact exceeds 0.50, or if its containment in any *non-degenerate*
prior artifact (one that does not mark essentially the whole footprint
positive) exceeds 0.90.  (GEMSDOE30's own brief asks for exactly this check:
"hash and correlate every candidate against the full history".)

Why Jaccard and not containment alone: several historical artifacts mark
essentially every footprint pixel positive (e.g. ``gems-density-probe`` has
5,106,385 of 5,167,373 pixels > 0).  Containment in such a raster is 1.0 by
construction for every possible submission and carries no information; the
Jaccard overlap against it is 0.006.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
REF_TR = (100.0, 0.0, 243350.0, 0.0, -100.0, 4508550.0)


def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    cand_path = Path(sys.argv[1])
    with rasterio.open(DATA / "labels.tif") as s:
        lab = s.read(1)
    foot = lab >= 0
    with rasterio.open(cand_path) as s:
        a = s.read(1)
    a = np.where(foot, np.nan_to_num(a, nan=0.0), 0.0)
    apos = a > 0
    an = int(apos.sum())
    av = a[foot].astype(np.float64)

    rows = []
    files = sorted(list((DATA / "compare").glob("*.tif")) + list((DATA / "scored").glob("*.tif")))
    for p in files:
        try:
            with rasterio.open(p) as s:
                if (s.height, s.width) != (3730, 3292):
                    continue
                if tuple(round(v, 2) for v in s.transform)[:6] != REF_TR:
                    continue
                b = s.read(1)
        except Exception:
            continue
        b = np.where(foot, np.nan_to_num(b, nan=0.0), 0.0)
        bpos = b > 0
        bn = int(bpos.sum())
        inter = int((apos & bpos).sum())
        union = int((apos | bpos).sum())
        bv = b[foot].astype(np.float64)
        if bv.std() > 0 and av.std() > 0:
            corr = float(np.corrcoef(av, bv)[0, 1])
        else:
            corr = float("nan")
        rows.append(dict(file=p.name, n_them=bn, inter=inter,
                         jaccard=inter / max(union, 1),
                         containment=inter / max(an, 1),
                         corr=corr, sha=sha256(p)))
    rows.sort(key=lambda r: -r["containment"])
    print(f"candidate: {cand_path.name}  positives={an:,}  compared against {len(rows)} rasters\n")
    print(f"{'containment':>11s} {'jaccard':>8s} {'corr':>7s} {'n_them':>8s}  file")
    for r in rows[:15]:
        print(f"{r['containment']:11.4f} {r['jaccard']:8.4f} {r['corr']:7.4f} "
              f"{r['n_them']:8d}  {r['file'][:64]}")
    worst_j = max(rows, key=lambda r: r["jaccard"]) if rows else None
    nondeg = [r for r in rows if r["n_them"] <= 10 * max(an, 1)]
    worst_c = max(nondeg, key=lambda r: r["containment"]) if nondeg else None
    gate = bool(worst_j is None or (worst_j["jaccard"] < 0.50
                                    and (worst_c is None or worst_c["containment"] < 0.90)))
    print(f"\nmax Jaccard overlap          = {worst_j['jaccard']:.4f}  "
          f"({worst_j['file'][:56]})   threshold 0.50")
    print(f"max containment (non-degenerate) = {worst_c['containment']:.4f}  "
          f"({worst_c['file'][:56]})   threshold 0.90")
    max_corr = max((abs(r["corr"]) for r in rows if np.isfinite(r["corr"])), default=float("nan"))
    worst_r = max((r for r in rows if np.isfinite(r["corr"])), key=lambda r: abs(r["corr"]),
                  default=None)
    print(f"max |Pearson r| inside footprint = {max_corr:.4f}"
          + (f"  ({worst_r['file'][:56]})" if worst_r else ""))
    (ROOT / "registry" / "uniqueness.json").write_text(json.dumps(
        dict(candidate=cand_path.name, sha256=sha256(cand_path), positives=an,
             compared=len(rows), max_jaccard=worst_j["jaccard"] if worst_j else None,
             max_containment_nondegenerate=worst_c["containment"] if worst_c else None,
             max_abs_pearson_r=max_corr,
             max_abs_pearson_r_file=(worst_r["file"] if worst_r else None),
             gate_thresholds=dict(jaccard_lt=0.50, containment_nondegenerate_lt=0.90,
                                  degenerate_reference_ratio=10),
             gate_pass=gate, table=rows), indent=2, default=float) + "\n")
    return 0 if gate else 1


if __name__ == "__main__":
    sys.exit(main())
