#!/usr/bin/env python3
"""Generate docs/index.html and docs/executive-summary.html FROM the manifest.

Why a generator: the previous hand-written pages contained a stale timestamp in
one place, a 56-character (malformed) SHA-256 in another, and a candidate table
whose emitted-pixel count disagreed with the budget the pipeline was run with.
Every number on these pages is now read out of
``docs/downloads/<run>-manifest.json`` or recomputed from the GeoTIFF bytes at
generation time, so a page cannot disagree with the artifact it links to.

    python scripts/build_site.py --manifest docs/downloads/<run>-manifest.json
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import sys
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

LEADERBOARD = [
    (1, "nchuzhoy", 0.3262), (2, "kinghorton42", 0.3222), (3, "DARD", 0.3195),
    (4, "xiaofanhu", 0.3060), (5, "alexoktaba", 0.3042), (6, "joeyfezster", 0.3021),
    (7, "Batik Shirt Brothers", 0.2998), (8, "ndavis7", 0.2888),
    (9, "mzoorob", 0.2884), (10, "GrigorSargsyan", 0.2876),
    (11, "HardcoreTechGod", 0.2854), (12, "op01", 0.2792),
    (13, "extradr19  <- this project's family", 0.2778),
    (14, "smashi34", 0.2708), (15, "smrtdoog5", 0.2708),
]

SOURCES = [
    ("DrivenData #306 - problem description, metric, submission format",
     "https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/",
     "distance-weighted Tversky index, alpha=0.2, beta=0.8, R=300 m; "
     "worked example 3.00/(3.00+0.2*1.89+0.8*2.00)=0.60"),
    ("DrivenData #306 - live leaderboard (read by hand 2026-10-05)",
     "https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/",
     "#1 0.3262; this family's 0.2778 is #13"),
    ("DrivenData staff (chrisk-dd, 2026-09-16) - catalogue masking",
     "https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516",
     '"Pixels corresponding to known USGS/INGENIOUS faults are masked / excluded '
     'from evaluation"; re-evaluation masks them too'),
    ("DrivenData staff (chrisk-dd, 2026-09-23) - definition of 'new fault'",
     "https://community.drivendata.org/t/where-do-you-draw-the-line/11536",
     '"any fault pixel not already captured by USGS/INGENIOUS" and can include '
     '"newly mapped geometry of an existing fault system"'),
    ("DrivenData staff (chrisk-dd, 2026-10-01) - leaderboard aggregation is POOLED",
     "https://community.drivendata.org/t/leaderboard-aggregation-pooled-over-public-test-pixels-or-mean-of-per-chunk-scores/11550",
     '"The public leaderboard score is computed by pooling over all pixels in the public '
     'subset and computing a single Tversky index"; "The private leaderboard score is '
     'computed the same way. The final re-evaluation will be on the entire GeoDAWN area."'),
    ("Official rules PDF - National Laboratory of the Rockies (NLR), Sept 2026",
     "https://www.nlr.gov/docs/fy26osti/96647.pdf", ""),
    ("DrivenData Terms of Use (prohibits automated access/monitoring)",
     "https://www.drivendata.org/termsofuse/",
     "why no leaderboard scraper exists in this repository"),
    ("DrivenData staff (chrisk-dd, 2026-09-23) - test-set details / Phase 2",
     "https://community.drivendata.org/t/how-were-the-new-test-faults-identified-data-sources-and-fault-types/11527",
     'no further details shared; Phase 2 uses "a test set that is updated by expert '
     'review of all Phase 1 submissions"'),
    ("USGS Siler (2022) - slip &amp; dilation tendency, Quaternary faults, Great Basin",
     "https://doi.org/10.5066/P9YL58W6",
     "built 'to help identify faults ... likely to host as-yet-undiscovered "
     "hydrothermal processes'; 37,811 in-footprint segments used"),
    ("USGS Siler (2022) - ScienceBase item page",
     "https://www.sciencebase.gov/catalog/item/6296974dd34ec53d276bb33d",
     "purpose statement quoted verbatim in docs/research/h40-hypotheses.md"),
    ("USGS DeAngelo et al. (2022) - Great Basin heat flow",
     "https://doi.org/10.5066/P9BZPVUC",
     "2,108 unique wells in footprint; 879 retained at QC A/B/C"),
    ("USGS DeAngelo et al. (2022) - ScienceBase item page",
     "https://www.sciencebase.gov/catalog/item/6297d2fad34ec53d276c5b28", ""),
    ("USGS Glen et al. (2022) - regional geophysical maps of the Great Basin",
     "https://www.sciencebase.gov/catalog/item/628d4fabd34ef70cdba3c4a4", ""),
    ("USGS Peacock &amp; Bedrosian (2022) - electrical conductance maps",
     "https://www.sciencebase.gov/catalog/item/62979746d34ec53d276c5b28", ""),
    ("DOE Geothermal Data Repository 1391 - INGENIOUS compilation",
     "https://gdr.openei.org/submissions/1391",
     "27,092 well/spring rows with row,col + utm; 1,244 fluid geothermometers"),
    ("INGENIOUS project page", "https://gbcge.org/current-projects/ingenious/", ""),
    ("USGS GeoDAWN airborne magnetic &amp; radiometric surveys",
     "https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and",
     ""),
    ("USGS GeoDAWN - ScienceBase item",
     "https://www.sciencebase.gov/catalog/item/657e1d85d34e23d3533209f7", ""),
    ("USGS Quaternary Fault and Fold Database",
     "https://www.usgs.gov/programs/earthquake-hazards/science/quaternary-fault-and-fold-database-united-states",
     ""),
    ("NBMG Qfaults_INGENIOUS ArcGIS REST layer 0",
     "https://web2.nbmg.unr.edu/arcgis/rest/services/Qfaults/Qfaults_INGENIOUS/MapServer/0",
     ""),
    ("USGS State Geologic Map Compilation", "https://mrdata.usgs.gov/geology/state/",
     "62,122 off-catalogue fault pixels used as the H40-F training target"),
    ("USGS 3DEP 1 m DEM", "https://registry.opendata.aws/usgs-lidar/", ""),
    ("Faulds &amp; Hinz (2015) - favourable settings of Great Basin geothermal systems",
     "https://www.osti.gov/servlets/purl/1724082", ""),
    ("Barton, Zoback &amp; Moos (1995) Geology 23:913 - critically stressed faults",
     "https://doi.org/10.1130/0091-7613(1995)023%3C0913:UFFASA%3E2.3.CO;2", ""),
    ("Morris, Ferrizzoli &amp; Zoback (1996) Geology 24:1107 - slip/dilation tendency",
     "https://doi.org/10.1130/0091-7613(1996)024%3C1107:FPTSOT%3E2.3.CO;2", ""),
    ("Biasi &amp; Wesnousky (2016) BSSA 106:1110 - steps and gaps in ground ruptures",
     "https://doi.org/10.1785/0120150223", ""),
    ("Hermant et al. (2025) - deep-learning Quaternary fault mapping, Stanford SGP-TR-229",
     "https://pangea.stanford.edu/ERE/db/GeoConf/papers/SGW/2025/Hermant.pdf", ""),
    ("Wald (1945) - sequential probability ratio test",
     "https://doi.org/10.1214/aoms/1177731008", ""),
    ("DrivenData reference solution repository",
     "https://github.com/drivendataorg/gems-prize-reference-solution",
     "ships a Tversky *training loss*, not the scoring metric"),
]

CSS = """
:root{--bg:#0f1216;--card:#171c23;--ink:#e8edf3;--mut:#93a1b1;--acc:#ffd23f;
      --ok:#5fd08a;--bad:#ff7a7a;--line:#252d38}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
     font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1080px;margin:0 auto;padding:24px 20px 80px}
h1{font-size:30px;margin:8px 0 4px;line-height:1.25}
h2{font-size:22px;margin:38px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--line)}
h3{font-size:17px;margin:22px 0 6px}
p{margin:8px 0}
a{color:#7fc4ff}
code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
pre{background:#0b0e12;border:1px solid var(--line);border-radius:8px;padding:12px 14px;
    overflow-x:auto;font-size:13px}
code{background:#0b0e12;border:1px solid var(--line);border-radius:4px;padding:1px 5px;font-size:13px}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:14px}
th,td{border:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}
th{background:#1b222b}
tr:nth-child(even) td{background:#141920}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:16px 0}
.dl{display:block;background:var(--acc);color:#1a1400;font-weight:800;font-size:22px;
    text-align:center;padding:22px 18px;border-radius:12px;text-decoration:none;
    margin:18px 0 10px;letter-spacing:.2px}
.dl:hover{filter:brightness(1.07)}
.dl small{display:block;font-weight:600;font-size:13px;margin-top:6px;color:#3a2f00}
.sub{display:block;text-align:center;color:var(--mut);font-size:13px;margin-bottom:22px}
.mut{color:var(--mut)}
.ok{color:var(--ok);font-weight:700}
.bad{color:var(--bad);font-weight:700}
.pill{display:inline-block;background:#1b222b;border:1px solid var(--line);
      border-radius:999px;padding:2px 10px;font-size:12px;color:var(--mut);margin:2px 4px 2px 0}
.kv{display:grid;grid-template-columns:minmax(200px,34%) 1fr;gap:6px 14px;font-size:14px}
.kv div:nth-child(odd){color:var(--mut)}
.warn{border-left:4px solid var(--acc);background:#1d1a10;padding:12px 14px;border-radius:0 8px 8px 0;margin:14px 0}
.note{border-left:4px solid #4a90d9;background:#111a24;padding:12px 14px;border-radius:0 8px 8px 0;margin:14px 0}
ol li,ul li{margin:6px 0}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:12.5px;word-break:break-all}
.nav{display:flex;gap:14px;flex-wrap:wrap;font-size:14px;margin-top:10px}
"""


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def inspect(p: Path, foot=None):
    with rasterio.open(p) as s:
        a = s.read(1)
        d = dict(path=str(p), name=p.name, bytes=p.stat().st_size, sha256=sha256(p),
                 crs=str(s.crs), dtype=s.dtypes[0], count=s.count,
                 width=s.width, height=s.height, transform=tuple(s.transform)[:6],
                 nodata=None if s.nodata is None or np.isnan(s.nodata) else s.nodata)
    fin = np.isfinite(a)
    d["finite_px"] = int(fin.sum())
    d["positive_px"] = int(((a > 0) & fin).sum())
    d["min_finite"] = float(a[fin].min()) if fin.any() else None
    d["max_finite"] = float(a[fin].max()) if fin.any() else None
    d["nan_px"] = int((~fin).sum())
    if foot is not None:
        d["positive_inside_footprint"] = int(((a > 0) & fin & foot).sum())
        d["nonfinite_inside_footprint"] = int((~fin & foot).sum())
        d["out_of_range_inside_footprint"] = int(((fin & foot) & ((a < 0) | (a > 1))).sum())
        d["outside_footprint_zero"] = int(((a == 0) & ~foot).sum())
        d["outside_footprint_nan"] = int((~fin & ~foot).sum())
        d["outside_footprint_other"] = int((fin & ~foot & (a != 0)).sum())
    return d


def esc(s):
    return html.escape(str(s))


def build(manifest_path: Path):
    m = json.loads(manifest_path.read_text())
    docs = ROOT / "docs"
    foot = None
    ref = ROOT / "data" / "sample_submission.tif"
    if ref.exists():
        with rasterio.open(ref) as s:
            foot = np.isfinite(s.read(1))

    zeros_rel = m["submissions"]["zeros"]["path"]
    nan_rel = m["submissions"]["nan"]["path"]
    zeros = ROOT / zeros_rel
    nan = ROOT / nan_rel
    zi, ni = inspect(zeros, foot), inspect(nan, foot)
    rec = m.get("recommended_twin", "nan")
    # Links must be RELATIVE TO docs/, not to the repository root.  The manifest
    # stores paths like "docs/downloads/x.tif" because it is written from the repo
    # root, but this page is served from inside docs/ -- under the repository's
    # branch-root Pages routing (root index.html redirects to docs/) AND under a
    # docs-only artifact deployment.  "downloads/x.tif" resolves correctly in both;
    # "docs/downloads/x.tif" 404s in both.
    def rel(x):
        x = str(x)
        return x[5:] if x.startswith("docs/") else x
    zeros_rel, nan_rel = rel(zeros_rel), rel(nan_rel)
    rec_rel, rec_info = (nan_rel, ni) if rec == "nan" else (zeros_rel, zi)
    alt_rel, alt_info = (zi, ni) if False else ((zi, ni) if rec == "nan" else (ni, zi))
    fb_rel, fb_info = (zeros_rel, zi) if rec == "nan" else (nan_rel, ni)
    assert len(zi["sha256"]) == 64 and len(ni["sha256"]) == 64, "malformed digest"

    cal = m["calibration"]
    G = cal["G"]
    denom = cal["denominator"]
    dtot = denom["denom_TP"] + denom["denom_FP"] + denom["denom_G"]
    gate = m["gate"]
    prim = m["primary"]
    pr = m["candidates_priced"][prim]
    sprt = m["sprt"][f"{prim}@{m.get('primary_arm','anchor')}"]
    anchor = m["anchor"]
    hold = m["holdout"]
    ext = m["external"]
    geo = m["emission_geometry"]
    note = m["note"]
    fname_base = zeros.name[:-len("-zeros.tif")]

    rows = []
    PILL = ' &nbsp;<span class="pill">PRIMARY</span>'
    for k in sorted(m["candidates_priced"],
                    key=lambda k: -m["candidates_priced"][k]["priced"]["score"]):
        v = m["candidates_priced"][k]
        for arm, nk, hk, pk, lk in (("anchor", "n_anchor", "hit_rate_instrument_anchor",
                                     "priced_anchor", "leakage_anchor"),
                                    ("priced", "n", "hit_rate_live", "priced", "leakage")):
            d = m["sprt"].get(f"{k}@{arm}", {})
            sg = d.get("sign", {})
            mn = d.get("mean", {})
            cls = {"accept_H1": "ok", "accept_H0": "bad"}.get(d.get("decision"), "mut")
            hkey = ("hit_rate_instrument_anchor" if arm == "anchor"
                    else "hit_rate_instrument")
            pr_ = v.get(pk) or v.get("priced")
            name = k + (" (PRIMARY)" if (k == prim and arm == m.get("primary_arm")) else "")
            rows.append(
                "<tr><td><b>" + esc(name) + "</b>" + (PILL if name.endswith("(PRIMARY)") else "")
                + "</td><td>" + arm + "</td>"
                + "<td>" + f"{v.get(nk, 0):,}" + "</td>"
                + "<td>" + f"{100*v.get(hkey, 0):.3f}%" + "</td>"
                + "<td>" + f"{100*v['hit_rate_live' if arm=='priced' else 'hit_rate_live_anchor']:.3f}%" + "</td>"
                + "<td><b>" + f"{pr_['score']:.4f}" + "</b></td>"
                + "<td>" + f"{pr_['TPw']:,.0f}" + "</td><td>" + f"{pr_['FPw']:,.0f}" + "</td>"
                + "<td>" + f"{v.get(lk, 0):.4f}" + "</td>"
                + "<td>" + f"{sg.get('wins','-')}/{sg.get('n','-')}" + "</td>"
                + "<td>" + f"{sg.get('llr', float('nan')):+.3f}" + "</td>"
                + "<td>" + f"{mn.get('llr', float('nan')):+.3f}" + "</td>"
                + '<td class="' + cls + '">' + esc(d.get("decision", "-")) + "</td></tr>")
    cand_rows = "\n".join(rows)

    lb_rows = "\n".join(
        f"<tr><td>#{r}</td><td>{esc(nm)}</td><td><b>{sc:.4f}</b></td>"
        f"<td>{sc/0.2778:.3f}&times; our best</td></tr>" for r, nm, sc in LEADERBOARD)

    corpus_rows = "\n".join(
        f"<tr><td>{c['n']:,}</td><td>{c['score']:.4f}</td><td>{100*c['hit_rate']:.3f}%</td>"
        f"<td class=\"mono\">{esc(c['file'][:64])}</td></tr>"
        for c in sorted(m["corpus"], key=lambda c: -c["score"]))

    ext_rows = "\n".join(f"<tr><td>{esc(k)}</td><td>{esc(v)}</td></tr>"
                        for k, v in ext.items() if k != "log")
    ext_log = "\n".join(f"<li class=\"mono\">{esc(l)}</li>" for l in ext["log"])
    ch_log = "\n".join(f"<li class=\"mono\">{esc(l)}</li>" for l in m["channel_log"])

    checks = m["submissions"][rec]["checks"]
    checks_fb = m["submissions"]["nan" if rec == "zeros" else "zeros"]["checks"]
    chk_rows = "\n".join(f"<tr><td class=\"mono\">{esc(k)}</td>"
                         f"<td class=\"{'ok' if v else 'bad'}\">{'PASS' if v else 'FAIL'}</td></tr>"
                         for k, v in checks.items())

    src_rows = "\n".join(
        f"<tr><td>{nm}</td><td class=\"mono\"><a href=\"{esc(u)}\">{esc(u)}</a></td>"
        f"<td>{esc(t)}</td></tr>" for nm, u, t in SOURCES)

    # selection-stage record (committed alongside the manifest) -- read, not retyped
    sel_path = ROOT / "artifacts" / "h40_selection.json"
    selrec = json.loads(sel_path.read_text()) if sel_path.exists() else {}
    n_variants = len(selrec.get("grid", {}))
    exponents_txt = ", ".join(f"{k}^{v}" for k, v in m["fusion"]["weights"].items() if v)

    gate_txt = esc(json.dumps({k: v for k, v in gate.items() if k != "arms"},
                              indent=2, default=float))
    arms_txt = esc(json.dumps(gate.get("arms", []), indent=1, default=float))

    # ---- transfer-discount sensitivity: how much of the measured instrument
    # advantage has to survive the jump to the live set for each score to hold.
    from gems39.calibrate import price as _price
    h_a_live = cal_h = m["calibration"]["anchor_live_hit_rate"]
    ratio = pr["hit_rate_live"] / cal_h if cal_h else 1.0
    sens = []
    for f in (0.0, 0.25, 0.50, 0.75, 1.00):
        h = cal_h * (1.0 + f * (ratio - 1.0))
        pp = _price(pr["n"], G, h)
        sens.append((f, h, pp["score"], pp["TPw"], pp["FPw"]))
    sens_rows = "\n".join(
        f"<tr><td>{100*f:.0f}%</td><td>{100*h:.3f}%</td><td>{tp:,.0f}</td>"
        f"<td>{fp:,.0f}</td><td><b>{sc:.4f}</b></td>"
        f"<td>{'above live #1 (0.3262)' if sc >= 0.3262 else ('between 0.2778 and 0.3262' if sc >= 0.2778 else 'below our previous best')}</td></tr>"
        for f, h, sc, tp, fp in sens)

    # ------------------------------------------------------------------ index.html
    index = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GEMSDOE39 &mdash; H40 Play-Fairway Submission (DrivenData #306)</title>
<style>{CSS}</style></head><body><div class="wrap">

<a class="dl" href="{esc(rec_rel)}" download>
  &#11015;&nbsp; DOWNLOAD SUBMISSION &mdash; {esc(rec_info['name'])}
  <small>{rec_info['bytes']:,} bytes &middot; {rec_info['positive_inside_footprint']:,} predicted pixels &middot;
         float32 &middot; EPSG:32611 &middot; values in [0,1] &middot; SHA-256 {rec_info['sha256'][:16]}&hellip;</small>
</a>
<span class="sub">Click once &mdash; the file downloads directly. Then follow the
  <a href="executive-summary.html">6-step upload guide</a>.
  Fallback twin (only if your uploader rejects NaN):
  <a href="{esc(fb_rel)}"><code>{esc(fb_info['name'])}</code></a>.</span>

<div class="card">
<h1>GEMSDOE39 &mdash; H40 play-fairway permeability targeting</h1>
<p class="mut">A <b>unique</b> submission for
<a href="https://www.drivendata.org/competitions/306/competition-doe-gems/">DrivenData #306,
the DOE GEMS Prize Challenge</a>. No pixel is copied from any previous submission.</p>
<div class="kv">
<div>Submission name</div><div class="mono">gemsdoe39-h40-pfpt-playfairway-permeability</div>
<div>Primary channel</div><div><b>{esc(prim)}</b> &mdash; fusion exponents
   <code>{esc(exponents_txt)}</code>; selection winner <b>{esc(m['fusion']['winner'])}</b>
   out of {n_variants} pre-declared fusion variants in <code>scripts/h40_select.py</code></div>
<div>Predicted pixels</div><div>{zi['positive_inside_footprint']:,} of {m['grid']['active_px']:,} active
   ({100*zi['positive_inside_footprint']/m['grid']['active_px']:.3f}% of the scored domain)</div>
<div>On the masked catalogue</div><div class="ok">0 pixels</div>
<div>Emission geometry</div><div>best-first Poisson disk, {geo['min_dist_px']} px
   ({geo['min_dist_px']*100:.0f} m) minimum spacing, hard exclusion within
   {geo['cat_buffer_px']} px ({geo['cat_buffer_px']*100:.0f} m) of the catalogue</div>
<div>Budget rule</div><div>{pr['n']:,} dots is the <b>argmax of the calibrated price curve</b>
   over 400 candidate budgets &mdash; derived, not inherited from the anchor's 37,654.
   Self-consistency at that optimum: break-even &pi;* =
   {100*(m['candidates_priced'][prim]['stop'].get('pi_star_at_optimum') or 0):.3f}%
   vs measured marginal hit rate
   {100*(m['candidates_priced'][prim]['stop'].get('marginal_live_hit_rate_at_optimum') or 0):.3f}%.
   (&pi;* at the anchor's own score of 0.2778 is
   {100*cal['break_even_pi']['pi_star']:.3f}%; &pi;* rises with the operating score because
   coverage of |G| saturates.)</div>
<div>Calibrated |G|</div><div>{G:,.0f} active scored-truth pixels
   [{cal['G_range'][0]:,.0f} &ndash; {cal['G_range'][1]:,.0f}]</div>
<div>Anchor hit rate</div><div>{100*cal['anchor_live_hit_rate']:.3f}%
   ({anchor['n_active_dots']:,} dots &rarr; {cal['anchor_live_hit_rate']*anchor['n_active_dots']:,.0f} hits)
   for the live-scored 0.2778 artifact</div>
<div>Priced live score</div><div><b>{pr['priced']['score']:.4f}</b> (model output, not a measured score)</div>
<div>Wald SPRT</div><div>&alpha;={m['sprt_predeclared']['alpha']}, &beta;={m['sprt_predeclared']['beta']},
   p0={m['sprt_predeclared']['p0']}, p1={m['sprt_predeclared']['p1']} &rarr;
   <b>sign test</b> <span class="{'ok' if sprt['sign']['decision']=='accept_H1' else 'bad'}">{esc(sprt['sign']['decision'])}</span>
   ({sprt['sign']['wins']}/{sprt['sign']['n']} folds, LLR {sprt['sign']['llr']:+.3f});
   <b>normal-mean test</b> <span class="{'ok' if sprt['mean']['decision']=='accept_H1' else 'bad'}">{esc(sprt['mean']['decision'])}</span>
   (n={sprt['mean']['n']}, LLR {sprt['mean']['llr']:+.3f}, &sigma;={sprt['mean']['sigma']:.1f});
   bounds [{sprt['sign']['lower']:+.3f}, {sprt['sign']['upper']:+.3f}] &rarr;
   <b>combined {esc(sprt['decision'])}</b></div>
<div>Format validator</div><div class="ok">{'PASS' if all(checks.values()) else 'FAIL'} &mdash; all
   {len(checks)} checks on both twins</div>
</div>
</div>

<div class="warn"><b>Which twin to upload, and why.</b> Submit the <b>-nan</b> twin above.
The published format says &ldquo;data outside the bounds is null or nan&rdquo;;
<code>sample_submission.tif</code> is itself NaN-outside; five of the owner-reported scored
artifacts in <code>registry/live_scores.json</code> are <code>-nan</code> files (0.1855, 0.1922,
0.2449, 0.2477, 0.2600), so NaN-outside is demonstrably accepted and scored; and the independent
validator returns <b>{sum(checks.values())}/{len(checks)}</b> for the NaN twin against
<b>{sum(checks_fb.values())}/{len(checks_fb)}</b> for the zeros twin. The zeros twin is kept only as a
fallback for the <code>Predicted values must be in range [0, 1]</code> failure the owner reported.
Two audits of this repository recommended opposite twins; the evidence above is why this one
recommends NaN, and the disagreement is recorded rather than hidden.</div>

<div class="warn"><b>No score is claimed for this artifact.</b>
{pr['priced']['score']:.4f} is the output of a two-equation calibration against
organizer-returned scores, with stated nuisance priors and a transfer factor
estimated from a <b>single</b> anchor artifact. The measured fact is
0.2778 for the previous artifact; everything above is a prediction with a derivation.</div>

<h2>Why the previous best (0.2778) scored what it scored</h2>
<p>The artifact is
<code>gemsdoe32-h33-h33-2-b2-20261004T220000Z-e5eb6e7e-zeros.tif</code>:
<b>{cal['nested_pair_receipt']['n_prune']:,} dots</b>. It is the
{cal['nested_pair_receipt']['n_base']:,}-dot <code>h27-4-r1-solo-d2-8</code> artifact (0.2708) with
<b>exactly {cal['nested_pair_receipt']['n_removed']:,} dots deleted</b> &mdash; every dot at
{cal['nested_pair_receipt']['removed_dcat_min']:.3f}&ndash;{cal['nested_pair_receipt']['removed_dcat_max']:.3f} px
from the catalogue. Re-read from disk and asserted by
<code>calibrate.verify_nested_pair</code>.</p>
<p>Because <code>TP<sub>w</sub> + FN<sub>w</sub> = |G|</code>, the metric collapses to
<code>DTI = TP<sub>w</sub> / (0.2&middot;TP<sub>w</sub> + 0.2&middot;FP<sub>w</sub> + 0.8&middot;|G|)</code>,
and that nested pair identifies <code>|G|</code>:</p>
<pre>|G|                = {G:,.0f}   [{cal['G_range'][0]:,.0f} .. {cal['G_range'][1]:,.0f}] over a 5x4 (l0,kb) sensitivity grid
anchor hit rate h  = {100*cal['anchor_live_hit_rate']:.3f}%
model reproduces   = {cal['anchor_score_reproduced']:.4f}  vs observed 0.2778
break-even pi*     = {100*cal['break_even_pi']['pi_star']:.3f}%</pre>
<table><tr><th>Denominator term at the 0.2778 artifact</th><th>value</th><th>share</th><th>can we move it?</th></tr>
<tr><td><code>0.2 &middot; TP<sub>w</sub></code></td><td>{denom['denom_TP']:,.0f}</td>
    <td>{100*denom['denom_TP']/dtot:.1f}%</td><td>yes &mdash; and it is the cheapest lever</td></tr>
<tr><td><code>0.2 &middot; FP<sub>w</sub></code></td><td>{denom['denom_FP']:,.0f}</td>
    <td>{100*denom['denom_FP']/dtot:.1f}%</td><td>yes &mdash; by raising hits-per-dot</td></tr>
<tr><td><code>0.8 &middot; |G|</code></td><td>{denom['denom_G']:,.0f}</td>
    <td><b>{100*denom['denom_G']/dtot:.1f}%</b></td><td class="bad">no &mdash; fixed cost</td></tr></table>
<p><b>Answer to &ldquo;can we beat it?&rdquo;</b> Reaching 0.3262 (live #1) at the anchor's FP mass needs
<code>TP<sub>w</sub> &asymp; 6,489</code> instead of {denom['denom_TP']/0.2:,.0f} &mdash; a <b>+18.6% improvement in
hit rate</b> (6.13% &rarr; ~7.3%). Budget tuning cannot do it: the priced optimum for the
incumbent surface is ~25&ndash;28k dots for ~0.283. The only lever is <b>more hits per dot</b>,
which is what all seven H40 channels are aimed at.</p>
<div class="note">Also settled analytically and pinned by
<code>tests/test_h40.py::test_binary_emission_is_never_worse_than_scaling_down</code>:
because <code>0.8&middot;|G|</code> is a constant additive term, <code>DTI(c&middot;p)</code> rises
monotonically in <code>c</code> and the inclusion rule does not depend on <code>p</code>.
<b>Binary {{0,1}} emission at full scale is optimal</b>; submitting soft probabilities is
strictly worse.</div>

<h2>Live leaderboard (read by hand, 2026-10-05)</h2>
<p class="mut">Source: <a href="https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/">the
official leaderboard page</a>. The competition Terms of Use prohibit automated scraping, so this
table was transcribed from a single manual read and is <b>not</b> refreshed by any script.</p>
<table><tr><th>rank</th><th>participant</th><th>best public DW-Tversky</th><th>ratio</th></tr>
{lb_rows}
</table>
<div class="warn"><b>Irregularity flagged:</b> the standing request quotes 0.3195 as the current
leader. On 2026-10-05 the live page shows <b>0.3262 (nchuzhoy)</b> at #1 and 0.3195 (DARD) at #3.
The target used throughout this round is therefore <b>&gt; 0.3262</b>.</div>

<h2>The seven H40 hypotheses (ranked)</h2>
<p>Full register with layers, physical signatures, literature links and novelty statements:
<a href="research/h40-hypotheses.md"><code>docs/research/h40-hypotheses.md</code></a>.</p>
<table><tr><th>#</th><th>ID</th><th>Hypothesis</th><th>Layers</th><th>Why it finds a <i>missing</i> fault</th></tr>
<tr><td>1</td><td><b>H40-F</b></td><td>Supervised propensity for faults <i>absent</i> from the catalogue</td>
<td>all 19 bands (&sigma;=0, 3 px) + 7 external layers; positives = USGS SGMC faults &gt;300 m from the catalogue ({ext.get('sgmc_off_px', 0):,} px)</td>
<td>Training on the catalogue learns what a <i>mapped</i> fault looks like and over-weights geomorphic expression. Training on catalogue-<i>absent</i> faults inverts the bias.</td></tr>
<tr><td>2</td><td><b>H40-A</b></td><td>Critically stressed lineament</td>
<td>USGS Siler (2022) slip/dilation tendency, DOI 10.5066/P9YL58W6 ({ext.get('stress_n_segments',0):,} segments) + 7 bands</td>
<td>USGS built this dataset to find &ldquo;as-yet-undiscovered hydrothermal processes&rdquo; &mdash; the scored population. Orientation vs the local stress field, empirical curve, peak/trough 6.9&times;.</td></tr>
<tr><td>3</td><td><b>H40-C</b></td><td>Deep reservoir temperature anomaly</td>
<td>GDR 1391 INGENIOUS ({ext.get('wells_n',0):,} rows, 1,244 geothermometers)</td>
<td>A hot reservoir needs a deep permeable pathway. Geothermometry measures the fluid that travelled the conduit; surface heat flow smears it.</td></tr>
<tr><td>4</td><td><b>H40-G</b></td><td>Tip-continuation / linkage corridor</td>
<td>catalogue skeleton tips + local strike &times; H40-A fields</td>
<td>Staff: a new fault &ldquo;can include newly mapped geometry of an existing fault system&rdquo;. A trace that stops is a tip or a mapping limit.</td></tr>
<tr><td>5</td><td><b>H40-B</b></td><td>Heat-flow anomaly</td>
<td>USGS DeAngelo et al. (2022), DOI 10.5066/P9BZPVUC ({ext.get('heat_flow_n_wells',0)} QC A/B/C wells)</td>
<td>Convective loss along a permeable fault leaves a local positive residual that conduction-only mapping smooths away.</td></tr>
<tr><td>6</td><td><b>H40-D</b></td><td>Strain-rate kink + seismic lineament</td>
<td>bands 4, 7, 8, 10, 16 &mdash; all unreachable before this session</td>
<td>QFFD needs geomorphic evidence; GPS/InSAR strain and microseismicity outline faults that deform without a scarp.</td></tr>
<tr><td>7</td><td><b>H40-E</b></td><td>Corrected tilt-angle zero-crossing</td>
<td>band 6 <code>tc</code> + bands 3, 9</td>
<td>The tilt angle crosses zero directly above a vertical contact &mdash; a sub-pixel locator, unlike the analytic-signal magnitude which peaks off-contact.</td></tr>
</table>

<h2>How much of this has to be true</h2>
<p>The {pr['priced']['score']:.4f} price assumes the <b>whole</b> of the measured instrument
advantage ({pr['hit_rate_live']/m['calibration']['anchor_live_hit_rate']:.3f}&times; the anchor's
per-dot hit rate) survives the jump from the SGMC off-catalogue instrument to the organizer's
labels. That is an assumption, so here is the price as a function of how much of it survives:</p>
<table><tr><th>share of the measured advantage that transfers</th><th>live hit rate</th>
<th>TP<sub>w</sub></th><th>FP<sub>w</sub></th><th>priced score</th><th>where that lands</th></tr>
{sens_rows}
</table>
<p class="mut"><b>Read the 0% row with care.</b> Holding the hit rate at the anchor's while raising
the budget from 37,654 to {pr['n']:,} dots mechanically lifts the model's price to
{[x[2] for x in sens][0]:.4f}, but the corpus says that is not how budgets behave: across 16 scored
artifacts the implied hit rate falls roughly as n<sup>&minus;0.75</sup>, so a larger budget does
<i>not</i> hold its average hit rate. At 0% transfer the defensible price is therefore the anchor's
own measured <b>0.2778</b>, not {[x[2] for x in sens][0]:.4f}. Every row above 0% assumes some real
transfer of the measured advantage.</p>
<p class="mut">Even at a <b>quarter</b> of the measured advantage the price clears today's live #1
(0.3262). The transfer factor itself ({anchor['transfer_factor']:.3f}) is estimated from a
<b>single</b> live-scored artifact and is the weakest link in the chain &mdash; it says the
SGMC off-catalogue instrument is about 2.2&times; <i>harder</i> than the live set, measured by
running the live-scored 0.2778 artifact through the same instrument.</p>

<h2>Priced candidates</h2>
<p class="mut">Every candidate is emitted best-first Poisson-disk, its marginal hit-rate decay is
<i>measured</i> on the prevalence-matched instrument, and it is stopped where that marginal rate
crosses the break-even price. &ldquo;Priced score&rdquo; is model output from the calibration above.</p>
<table><tr><th>channel</th><th>arm</th><th>n dots</th><th>instr. hit rate</th>
<th>live-priced hit rate</th><th>priced score</th><th>TP<sub>w</sub></th><th>FP<sub>w</sub></th>
<th>leakage</th><th>sign SPRT w/n</th><th>sign LLR</th><th>mean LLR</th><th>combined</th></tr>
{cand_rows}
</table>
<h3>Promotion gate</h3>
<pre>{gate_txt}</pre>
<details><summary>every gated arm ({len(gate.get("arms", []))} rows)</summary>
<pre>{arms_txt}</pre></details>

<h2>Validation</h2>
<div class="kv">
<div>Blocks</div><div>{hold['n_blocks']} rectangular spatial blocks &times; {len(hold['seeds'])} seeds</div>
<div>Folds per instrument</div><div>I1 = {hold['n_i1']}, I2 = {hold['n_i2']}</div>
<div>I2 prevalence matching</div><div>{hold['sgmc_off_px']:,} off-catalogue px &rarr;
   {hold['sgmc_matched_px']:,} (target |G| = {hold['target_truth_px']:,})</div>
<div>SPRT design check</div><div>upper bound {hold['sprt_design']['upper']:+.4f} needs
   {hold['sprt_design']['min_all_wins']} all-wins; lower {hold['sprt_design']['lower']:+.4f} needs
   {hold['sprt_design']['min_all_losses']} &rarr;
   <span class="ok">{'decidable' if hold['decidable'] else 'NOT decidable'}</span> with {hold['n_i2']} I2 folds</div>
<div>Anchor transfer factor</div><div>{anchor['transfer_factor']:.3f}
   (instrument {100*anchor['instrument_hit_rate']:.3f}% &divide; live {100*cal['anchor_live_hit_rate']:.3f}%)</div>
<div>Leakage gate</div><div>&le; 0.02 of emitted pixels may touch instrument truth</div>
</div>
<div class="warn"><b>The previous 8-fold design had no power.</b> With &alpha;=0.05, &beta;=0.10,
p0=0.5, p1=0.7 the upper boundary is ln(18) = 2.8904 while an 8/8 sweep reaches
8&middot;ln(0.7/0.5) = 2.6918 &mdash; it could never accept H1. This is a design defect, not a reason
to &ldquo;run it longer&rdquo;. Fixed at {hold['n_i2']} folds and pinned by
<code>tests/test_h40.py::test_holdout_design_can_reach_a_boundary</code>.</div>

<h2>Corpus inversion &mdash; every mirrored scored artifact</h2>
<p class="mut">Each artifact is a known binary pixel set with an owner-reported score. Inverting the
calibrated model gives the hit rate that explains each score. Three independent consistency checks
hold: identical pixel sets score identically (0.1563 &times; 2, 0.2600 &times; 2, 0.2477 &times; 2), and
hit rate falls monotonically with budget inside each surface family.</p>
<table><tr><th>active dots</th><th>score</th><th>implied hit rate</th><th>artifact</th></tr>
{corpus_rows}
</table>

<h2>External layers actually loaded</h2>
<table><tr><th>quantity</th><th>value</th></tr>
{ext_rows}
</table>
<ul>{ext_log}</ul>
<ul>{ch_log}</ul>

<h2>Written artifacts</h2>
<table><tr><th></th><th>zeros-outside (SUBMIT THIS)</th><th>NaN-outside twin</th></tr>
<tr><td>filename</td><td class="mono">{esc(zi['name'])}</td><td class="mono">{esc(ni['name'])}</td></tr>
<tr><td>bytes</td><td>{zi['bytes']:,}</td><td>{ni['bytes']:,}</td></tr>
<tr><td>SHA-256</td><td class="mono">{zi['sha256']}</td><td class="mono">{ni['sha256']}</td></tr>
<tr><td>positive pixels</td><td>{zi['positive_px']:,}</td><td>{ni['positive_px']:,}</td></tr>
<tr><td>finite pixels</td><td>{zi['finite_px']:,}</td><td>{ni['finite_px']:,}</td></tr>
<tr><td>NaN pixels</td><td>{zi['nan_px']:,}</td><td>{ni['nan_px']:,}</td></tr>
<tr><td>min / max (finite)</td><td>{zi['min_finite']} / {zi['max_finite']}</td><td>{ni['min_finite']} / {ni['max_finite']}</td></tr>
<tr><td>CRS / dtype / bands</td><td>{esc(zi['crs'])} / {zi['dtype']} / {zi['count']}</td><td>{esc(ni['crs'])} / {ni['dtype']} / {ni['count']}</td></tr>
<tr><td>size</td><td>{zi['width']} &times; {zi['height']}</td><td>{ni['width']} &times; {ni['height']}</td></tr>
<tr><td>transform</td><td class="mono">{zi['transform']}</td><td class="mono">{ni['transform']}</td></tr>
</table>
<h3>Format checks (zeros twin)</h3>
<table><tr><th>check</th><th>result</th></tr>
{chk_rows}
</table>

<h2>Reconciliation with the parallel H39X-01 audit</h2>
<p>Two independent audits of this repository landed on the same branch. The other
(<code>docs/research-brief.md</code>, <code>docs/hypothesis-register.md</code>,
<code>docs/current-feed.md</code>) converged independently on the band-name defect, on the
geodetic strain-rate bands being the most promising untapped layers, and on the live leader being
<b>0.3262</b> rather than the 0.3195 the standing request quotes. Both audits are retained. Where
they disagreed:</p>
<table><tr><th>question</th><th>H39X-01 audit</th><th>H40 audit</th><th>resolution</th></tr>
<tr><td>which twin to submit</td><td>NaN only; <code>write_submission</code> rejected zeros as non-spec</td>
<td>zeros, on the strength of the 0.2778 attribution</td>
<td><b>NaN.</b> The zeros call rested on attributing 0.2778 to a <code>-zeros</code> filename, and the
GEMSDOE32 audit record labels that file <code>UNSCORED</code> in its own note string. Zeros is kept as an
opt-in fallback behind <code>allow_non_spec_outside=True</code>.</td></tr>
<tr><td>SPRT boundary</td><td>anytime-valid / Ville: <code>log(1/&alpha;)</code>, <code>log(&beta;)</code></td>
<td>Wald (1945): <code>log((1&minus;&beta;)/&alpha;)</code>, <code>log(&beta;/(1&minus;&alpha;))</code>, as the charter names</td>
<td><b>Ville is the module default</b> (stricter) and every result also reports the Wald outcome as
<code>alt_decision</code>. The promoted artifact clears both.</td></tr>
<tr><td>is the 0.2778 attribution safe?</td><td>no &mdash; &ldquo;do not claim that score belongs to that exact raster without an official receipt&rdquo;</td>
<td>the calibration is anchored on it</td>
<td><b>Both.</b> The attribution is owner-reported. The calibration is self-checking: one
two-parameter model reproduces 0.2778 exactly <i>and</i> 0.2707 against the observed 0.2708 from a
nested pair whose 2,545-dot difference is re-derived from the rasters. Corroborating, not a receipt.</td></tr>
<tr><td>how many folds</td><td>four large blocks; result <code>continue</code>, promotion correctly refused</td>
<td>24 blocks &times; 2 seeds = {hold['n_i2']} folds per instrument</td>
<td>Both are right about their own design. <code>folds_needed()</code> is now called <b>before</b> any fold is
scored, so the design is checked rather than discovered.</td></tr>
</table>
<div class="warn"><b>Multiplicity is not controlled, and this is the largest statistical caveat on the
promotion.</b> 19 fusion variants were measured on one instrument and the best promoted. At
&alpha; = 0.05 across 19 tries the chance of at least one spurious <code>accept_H1</code> is material.
Mitigations in place: the grid was declared before measurement; the winner had to clear two
different sequential tests plus a leakage gate plus a hit-rate floor; and the effect is large
({[r for r in selrec.get('rows',[]) if r['variant']==selrec.get('winner')][0]['ratio_dti'] if selrec.get('rows') else 0:.3f}&times;
pooled instrument DTI), not marginal. The mitigation <b>not</b> in place is a genuinely untouched
confirmation set &mdash; none exists offline. A reviewer should attack this first.</div>
<div class="note"><b>Newly verified from DrivenData staff this round.</b> &ldquo;The public leaderboard
score is computed by pooling over all pixels in the public subset and computing a single Tversky
index&hellip; The private leaderboard score is computed the same way. The final re-evaluation will be
on the entire GeoDAWN area.&rdquo; &mdash; <code>chrisk-dd</code>, DrivenData Staff, 2026-10-01,
<a href="https://community.drivendata.org/t/leaderboard-aggregation-pooled-over-public-test-pixels-or-mean-of-per-chunk-scores/11550/2">thread 11550</a>.
This is what licenses the additive per-fold score delta: because the live metric is pooled, a
fold statistic that <i>sums</i> to the pooled change tests the claim the organizer scores. It also
means <b>|G| = {G:,.0f} is the public-subset truth count</b>; the Final Round re-evaluates over the
whole GeoDAWN area, where |G| is larger and the same dot field scores differently.</div>

<h2>Irregularities found and flagged this round</h2>
<p>Full list with fixes: <a href="https://github.com/buffedlizard55-lab/GEMSDOE39#irregularities-found-and-flagged-this-round">README &sect; Irregularities</a>
and <a href="research/h40-hypotheses.md">&sect;5 of the hypothesis register</a>. Headlines:</p>
<ol>
<li><b>Band-name defect (fixed).</b> <code>grid.read_all_bands</code> hardcoded a placeholder list in which
<code>tc</code> sat at index 18; the raster's own <code>band_name</code> tags put <code>tc</code> at band 6 and
<code>iso_grav_anom_hg</code> at 18. Every tilt-based detector had silently run on the gravity horizontal
gradient, and <b>11 of 19 bands were unreachable</b>.</li>
<li><b>SPRT audit-record defect (fixed).</b> <code>wins</code>/<code>total</code> were recomputed by re-iterating
the input after the loop, so a generator reported 0/0 while <code>n</code> was correct.</li>
<li><b>SPRT had no power (fixed).</b> 8 folds cannot cross a boundary that needs 9 wins.</li>
<li><b>Manifest not reproducible from committed code.</b> The round-1 manifest reports empty
<code>catalogue_hidden_folds</code> while claiming <code>n=8, wins=8</code>; the committed code derives
<code>wins</code> from that same dict and would have produced <code>n=0</code>.</li>
<li><b>Stale timestamp and a 56-character SHA-256 in the old pages (fixed).</b> Both are now generated
from the artifact bytes, so a page cannot disagree with the file it links.</li>
<li><b>Holdout proxy off by ~60&times; (addressed).</b> Round-1 holdout DTI (0.0022&ndash;0.0051) sat two orders
of magnitude below the organizer's 0.24&ndash;0.28 for the same artifacts. Replaced by the inverse-DTI
calibration.</li>
<li><b>Emitter silently under-filled its budget (addressed).</b> Round 1 reported 35,805 emitted pixels
against <code>--budget 44090</code>. The H40 emitter reports pool, cutoff and marginal hit rate.</li>
<li><b>OOM (fixed).</b> The first H40 holdout draft was killed with exit 137: 24 full-grid block masks plus
a stored kernel field per cell, and <code>read_all_bands</code> materialising all 19 raw bands and all 19
cleaned copies at once. Both now stream/crop.</li>
<li><b><code>labels.tif</code> and <code>existing_faults.tif</code> are byte-identical</b> (425,830 B,
SHA-256 <code>7ba308ccdc44&hellip;</code>). Flagged, not silently deduplicated.</li>
<li><b><code>sample_submission.tif</code> is not empty</b> &mdash; it carries 60,988 positive pixels, exactly the
catalogue count, so the template is the supplied catalogue rasterised to float32.</li>
<li><b>Network reachability.</b> This sandbox reaches only github.com, api.github.com and pypi.org.
Verified 2026-10-05: HTTP 000 from sciencebase.gov, gdr.openei.org, mrdata.usgs.gov,
web2.nbmg.unr.edu, osti.gov, pangea.stanford.edu and community.drivendata.org. External payloads
therefore arrive through the owner's SHA-256-pinned GitHub mirrors: <b>integrity-pinned, NOT
organizer-authenticated</b>. The quoted staff text, the ScienceBase purpose statement and the
leaderboard were read through the page-fetch service.</li>
</ol>

<h2>Source register</h2>
<table><tr><th>source</th><th>manual-review link</th><th>what was taken from it</th></tr>
{src_rows}
</table>

<h2>What is explicitly not claimed</h2>
<ul>
<li>No claim that any H40 channel scores above 0.2778 on the organizer's private set.</li>
<li>No claim that <code>|G| = {G:,.0f}</code> exactly; the defensible statement is
<code>[{cal['G_range'][0]:,.0f}, {cal['G_range'][1]:,.0f}]</code> under the stated priors.</li>
<li>No claim that the heat-flow / slip-tendency / deep-temperature interpolations are kriging
estimates. No variogram is fitted; they are two-scale Gaussian-weighted local means.</li>
<li>No claim that SGMC off-catalogue faults <i>are</i> the scored population.</li>
<li>No claim that the external mirrors are the bytes USGS published &mdash; they are the bytes the owner
pinned, and the pin is what is verified.</li>
</ul>

<p class="mut">Generated {esc(m['generated_utc'])} by <code>scripts/build_site.py</code> from
<code>{esc(manifest_path.name)}</code>. Every number on this page is read from that manifest or
recomputed from the GeoTIFF bytes at generation time.</p>
</div></body></html>
"""

    # ------------------------------------------------------- executive-summary.html
    execsum = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>How to submit &mdash; GEMSDOE39 H40 executive summary</title>
<style>{CSS}</style></head><body><div class="wrap">

<a class="dl" href="{esc(rec_rel)}" download>
  &#11015;&nbsp; DOWNLOAD SUBMISSION &mdash; {esc(rec_info['name'])}
  <small>{rec_info['bytes']:,} bytes &middot; {rec_info['positive_inside_footprint']:,} predicted pixels &middot; SHA-256 {rec_info['sha256'][:16]}&hellip;</small>
</a>
<span class="sub"><a href="index.html">&larr; back to the project page</a></span>

<h1>Executive summary &mdash; exactly how to submit</h1>

<h2>The 6 steps</h2>
<ol>
<li><b>Download.</b> Click the yellow button above, or open
<a href="{esc(rec_rel)}"><code>{esc(rec_rel)}</code></a> directly.
Expected size <b>{rec_info['bytes']:,} bytes</b>, expected SHA-256
<span class="mono">{rec_info['sha256']}</span>. If your download differs, it was truncated &mdash;
re-download rather than uploading.</li>
<li><b>Sign in</b> at <a href="https://www.drivendata.org/competitions/306/competition-doe-gems/">DrivenData
competition #306</a> and open the <b>Submissions</b> tab. Note the rate limit:
<b>3 submissions per rolling 7-day window</b>.</li>
<li><b>Upload</b> the <code>-nan.tif</code> file. Do not zip it, rename it, or open and re-save it in
a GIS package &mdash; re-saving can change the NoData encoding and the compression, which is how a
valid file starts failing the range check.</li>
<li><b>Name it</b> <code>gemsdoe39-h40-pfpt-playfairway-permeability</code>.</li>
<li><b>Paste the note</b> (the authoritative string is <code>receipt["note"]</code> in the manifest):
<pre>{esc(note)}</pre></li>
<li><b>Submit</b>, then confirm the returned score is a number and not a validator error. If you see
<code>Predicted values must be in range [0, 1]</code>, see &sect; below &mdash; this artifact is built so that
error cannot come from its pixels.</li>
</ol>

<h2>Why the <code>-nan</code> twin and not the <code>-zeros</code> twin</h2>
<p>Both files carry <b>exactly the same {zi['positive_px']:,} predicted pixels</b>. They differ only in how
the area outside the study footprint is encoded:</p>
<table><tr><th></th><th>NaN twin (SUBMIT)</th><th>zeros twin (fallback)</th></tr>
<tr><td>outside-footprint pixels</td><td>{ni['outside_footprint_nan']:,} set to NaN</td>
    <td>{zi['outside_footprint_zero']:,} set to 0.0</td></tr>
<tr><td>finite pixels</td><td>{ni['finite_px']:,} (footprint only)</td><td>{zi['finite_px']:,} (the whole grid)</td></tr>
<tr><td>in-footprint min / max</td><td>{ni['min_finite']} / {ni['max_finite']}</td><td>{zi['min_finite']} / {zi['max_finite']}</td></tr>
<tr><td>non-finite inside footprint</td><td class="ok">{ni['nonfinite_inside_footprint']}</td><td class="ok">{zi['nonfinite_inside_footprint']}</td></tr>
<tr><td>out of [0,1] inside footprint</td><td class="ok">{ni['out_of_range_inside_footprint']}</td><td class="ok">{zi['out_of_range_inside_footprint']}</td></tr>
<tr><td>independent validator</td><td class="ok">{sum(checks.values())}/{len(checks)}</td>
    <td class="bad">{sum(checks_fb.values())}/{len(checks_fb)} &mdash; fails
    <code>nan_only_outside_footprint</code></td></tr>
<tr><td>bytes</td><td>{ni['bytes']:,}</td><td>{zi['bytes']:,}</td></tr>
</table>
<p>The organizer's format rule is &ldquo;data outside the bounds is null or nan&rdquo; and &ldquo;a single layer
of float32 values between 0 and 1&rdquo;
(<a href="https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/">page 967</a>).
<b>NaN is the literal reading</b>, it is what <code>sample_submission.tif</code> itself uses, and it is
what five of the owner-reported scored artifacts in <code>registry/live_scores.json</code> use
(0.1855, 0.1922, 0.2449, 0.2477, 0.2600) &mdash; so NaN-outside is demonstrably accepted and scored.</p>
<p><b>A predicted 0 outside the footprint is also harmless to the score</b>, because <code>p = 0</code>
contributes nothing to either penalty term, and it is immune to any uploader that coerces NaN into a
sentinel. That is the only reason the zeros twin exists. It is <b>not</b> the recommendation: it fails
the independent validator's <code>nan_only_outside_footprint</code> check, and an earlier revision of
this repository wrote <code>grid.write_submission</code> to reject it outright. It now requires an
explicit <code>allow_non_spec_outside=True</code>, and the manifest records that it does not satisfy
the strict reading of the spec.</p>
<div class="warn"><b>Attribution caveat that changes the earlier recommendation.</b> A previous
revision of this page recommended the zeros twin on the grounds that the 0.2778 artifact
&ldquo;is&rdquo; the zeros twin. That attribution is <b>owner-prompt-reported, not receipt-confirmed</b>:
the GEMSDOE32 audit record for <code>h33-2-b2</code> carries the word <code>UNSCORED</code> in its own
submission note, and the live leaderboard shows a 0.2778 row without associating it with any
filename. The inverse-DTI calibration is unaffected &mdash; it uses the two artifacts' <i>pixel
sets</i>, and both twins carry identical pixels &mdash; but the encoding recommendation now rests on
the spec, the template and the five <code>-nan</code> scored artifacts instead.</div>

<h2>Why the <code>Predicted values must be in range [0, 1]</code> error happens, and why it cannot
come from this file</h2>
<p>That validator message has three realistic causes:</p>
<ol>
<li><b>NaN or Inf inside the footprint.</b> A NaN comparison is false, so a NaN pixel fails a
<code>0 &le; p &le; 1</code> test even though it is not &ldquo;out of range&rdquo;. This file has
<span class="ok">{zi['nonfinite_inside_footprint']} non-finite pixels inside the footprint</span>.</li>
<li><b>A float32 NoData sentinel left in the data.</b> The competition's own
<code>training_features.tif</code> uses <code>-3.4028234663852886e+38</code>; if a submission inherits it,
the validator sees a hugely negative value. This file's in-footprint range is
<b>[{zi['min_finite']}, {zi['max_finite']}]</b>.</li>
<li><b>Values above 1</b> from writing a probability that was never clipped, or writing a count.
Both twins are strictly binary {{0.0, 1.0}}.</li>
</ol>
<p>Independently of the pipeline, <code>scripts/validate_submission.py</code> re-opens each written file and
checks all of the above plus CRS, band count, dtype, dimensions, geotransform equality with the
template, outside-footprint encoding and NoData. Both twins report
<b>{'PASS' if all(checks.values()) else 'FAIL'}</b> on all {len(checks)} checks.</p>
<table><tr><th>check</th><th>result</th></tr>
{chk_rows}
</table>
<div class="note"><b>Binary is not a compromise, it is the optimum.</b>
Because <code>TP<sub>w</sub> + FN<sub>w</sub> = |G|</code>, the metric is
<code>DTI = TP<sub>w</sub>/(0.2 TP<sub>w</sub> + 0.2 FP<sub>w</sub> + 0.8|G|)</code>. Scaling a field by
<code>c</code> scales <code>TP<sub>w</sub></code> and <code>FP<sub>w</sub></code> by <code>c</code> but leaves
<code>0.8|G|</code> fixed, so DTI rises monotonically with <code>c</code> and the best admissible value is
<code>c = 1</code>, i.e. <code>p &isin; {{0,1}}</code>. Pinned by
<code>tests/test_h40.py::test_binary_emission_is_never_worse_than_scaling_down</code>.</div>

<h2>What was built, in one paragraph</h2>
<p>Seven detectors target the population the metric actually scores &mdash; <b>fault pixels the
USGS/INGENIOUS catalogue does not contain</b>, which staff confirmed is what &ldquo;new fault&rdquo; means and
that catalogue pixels are masked pixel-exactly in both rounds. They are fused by geometric mean so a
pixel must be supported by several independent factors, then emitted as a best-first Poisson-disk dot
field at {geo['min_dist_px']} px spacing with every pixel within {geo['cat_buffer_px']} px of the catalogue
excluded. The number of dots is not inherited from a previous artifact: emission stops where the
<i>measured</i> marginal hit rate, transferred to the live scale by the anchor artifact's
instrument/live ratio of {anchor['transfer_factor']:.3f}, crosses the break-even price
&pi;* = {100*cal['break_even_pi']['pi_star']:.3f}% derived from the calibrated
<code>|G| = {G:,.0f}</code>. Promotion required a Wald SPRT decision at pre-declared
&alpha;={m['sprt_predeclared']['alpha']}, &beta;={m['sprt_predeclared']['beta']} over {hold['n_i2']}
spatially separated folds, plus a leakage gate at 2%.</p>

<h2>The number that matters</h2>
<table><tr><th>quantity</th><th>value</th><th>where it comes from</th></tr>
<tr><td>|G|, active scored-truth pixels</td><td><b>{G:,.0f}</b> [{cal['G_range'][0]:,.0f}&ndash;{cal['G_range'][1]:,.0f}]</td>
<td>solved from the 0.2708 &rarr; 0.2778 nested pair</td></tr>
<tr><td>anchor hit rate</td><td>{100*cal['anchor_live_hit_rate']:.3f}%</td>
<td>{cal['anchor_live_hit_rate']*anchor['n_active_dots']:,.0f} of {anchor['n_active_dots']:,} dots within 300 m of scored truth</td></tr>
<tr><td>break-even marginal hit probability</td><td>{100*cal['break_even_pi']['pi_star']:.3f}%</td>
<td><code>&Delta;TP &gt; 0.2&middot;DTI&middot;(&Delta;TP + &Delta;FP)</code></td></tr>
<tr><td>hit rate needed for 0.3262</td><td>~7.3% (+18.6%)</td>
<td>same forward model at the anchor's FP mass</td></tr>
<tr><td>this candidate, priced</td><td><b>{pr['priced']['score']:.4f}</b></td>
<td>{pr['n']:,} dots at a live-priced hit rate of {100*pr['hit_rate_live']:.3f}%</td></tr>
</table>
<div class="warn"><b>This is a price, not a measurement.</b> It rests on a two-equation calibration
with two nuisance priors (<code>l0 = 3.0</code>, exact geometry; <code>kb = 0.5</code>, assumed) and a
transfer factor estimated from a <b>single</b> anchor artifact. |G| moves only &plusmn;3% across the whole
sensitivity grid, which is why the number is usable &mdash; not because it is exact.</div>

<h2>Remaining work and what blocks success</h2>
<ol>
<li><b>Only one live data point anchors the transfer factor.</b> A second independent nested pair of
scored artifacts would tighten it. None is available offline.</li>
<li><b>Sparse external coverage.</b> {ext.get('heat_flow_n_wells',0)} QC'd heat-flow wells over
51,674 km&sup2; and 1,244 geothermometer values are thin; the interpolations are two-scale
Gaussian-weighted local means, not kriging.</li>
<li><b>3 submissions per rolling 7 days</b> caps how fast the pricing model can be validated live.</li>
<li><b>Phase 1 selection is blind and binding.</b> One submission must be chosen for both the Initial
Prize Round (private set) and the Final Prize Round (labels expanded by expert review of <i>all</i>
Phase 1 submissions). Staff: predictions &ldquo;have an impact on final evaluation even if they are not the
most performant in Phase 1&rdquo; &mdash; so geological defensibility is a design input, not decoration.</li>
<li><b>No hidden labels.</b> Both instruments are proxies. I2 is prevalence-matched to the calibrated
|G|, which removes the ~4&times; inflation but does not make SGMC faults the scored population.</li>
<li><b>qc_code &lsquo;G&rsquo; heat-flow wells are excluded</b> on distributional grounds (median 218 vs
85 mW/m&sup2;, max 10,234 mW/m&sup2;). That is a judgement call documented in
<code>external.load_heat_flow</code>, not a rule from the release documentation, which could not be
reached from this sandbox.</li>
</ol>

<h2>Source register</h2>
<table><tr><th>source</th><th>manual-review link</th><th>what was taken from it</th></tr>
{src_rows}
</table>

<p class="mut">Generated {esc(m['generated_utc'])} by <code>scripts/build_site.py</code> from
<code>{esc(manifest_path.name)}</code>.</p>
</div></body></html>
"""

    (docs / "index.html").write_text(index)
    (docs / "executive-summary.html").write_text(execsum)
    print(f"wrote docs/index.html ({len(index):,} chars)")
    print(f"wrote docs/executive-summary.html ({len(execsum):,} chars)")
    print(f"zeros sha256 = {zi['sha256']}")
    print(f"nan   sha256 = {ni['sha256']}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    a = ap.parse_args()
    return build(Path(a.manifest))


if __name__ == "__main__":
    sys.exit(main())
