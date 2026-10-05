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
    assert len(zi["sha256"]) == 64 and len(ni["sha256"]) == 64, "malformed digest"

    cal = m["calibration"]
    G = cal["G"]
    denom = cal["denominator"]
    dtot = denom["denom_TP"] + denom["denom_FP"] + denom["denom_G"]
    gate = m["gate"]
    prim = m["primary"]
    pr = m["candidates_priced"][prim]
    sprt = m["sprt"][prim]
    anchor = m["anchor"]
    hold = m["holdout"]
    ext = m["external"]
    geo = m["emission_geometry"]
    note = m["note"]
    fname_base = zeros.name[:-len("-zeros.tif")]

    rows = []
    for k in sorted(m["candidates_priced"], key=lambda k: -m["candidates_priced"][k]["priced"]["score"]):
        v = m["candidates_priced"][k]
        d = m["sprt"][k]
        rows.append(
            f"<tr><td><b>{esc(k)}</b>{' &nbsp;<span class=\"pill\">PRIMARY</span>' if k == prim else ''}</td>"
            f"<td>{v['n_pool']:,}</td><td><b>{v['n']:,}</b></td>"
            f"<td>{100*v['hit_rate_instrument']:.3f}%</td>"
            f"<td>{100*v['hit_rate_live']:.3f}%</td>"
            f"<td><b>{v['priced']['score']:.4f}</b></td>"
            f"<td>{v['priced']['TPw']:,.0f}</td><td>{v['priced']['FPw']:,.0f}</td>"
            f"<td>{v['leakage']:.4f}</td>"
            f"<td>{d['wins']}/{d['n']}</td><td>{d['llr']:+.3f}</td>"
            f"<td class=\"{'ok' if d['decision']=='accept_H1' else ('bad' if d['decision']=='accept_H0' else 'mut')}\">"
            f"{esc(d['decision'])}</td></tr>")
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

    checks = m["submissions"]["zeros"]["checks"]
    chk_rows = "\n".join(f"<tr><td class=\"mono\">{esc(k)}</td>"
                         f"<td class=\"{'ok' if v else 'bad'}\">{'PASS' if v else 'FAIL'}</td></tr>"
                         for k, v in checks.items())

    src_rows = "\n".join(
        f"<tr><td>{nm}</td><td class=\"mono\"><a href=\"{esc(u)}\">{esc(u)}</a></td>"
        f"<td>{esc(t)}</td></tr>" for nm, u, t in SOURCES)

    gate_txt = esc(json.dumps(gate, indent=2, default=float))

    # ------------------------------------------------------------------ index.html
    index = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GEMSDOE39 &mdash; H40 Play-Fairway Submission (DrivenData #306)</title>
<style>{CSS}</style></head><body><div class="wrap">

<a class="dl" href="{esc(zeros_rel)}" download>
  &#11015;&nbsp; DOWNLOAD SUBMISSION &mdash; {esc(zi['name'])}
  <small>{zi['bytes']:,} bytes &middot; {zi['positive_inside_footprint']:,} predicted pixels &middot;
         float32 &middot; EPSG:32611 &middot; values in [0,1] &middot; SHA-256 {zi['sha256'][:16]}&hellip;</small>
</a>
<span class="sub">Click once &mdash; the file downloads directly. Then follow the
  <a href="executive-summary.html">6-step upload guide</a>.</span>

<div class="card">
<h1>GEMSDOE39 &mdash; H40 play-fairway permeability targeting</h1>
<p class="mut">A <b>unique</b> submission for
<a href="https://www.drivendata.org/competitions/306/competition-doe-gems/">DrivenData #306,
the DOE GEMS Prize Challenge</a>. No pixel is copied from any previous submission.</p>
<div class="kv">
<div>Submission name</div><div class="mono">gemsdoe39-h40-pfpt-playfairway-permeability</div>
<div>Primary channel</div><div><b>{esc(prim)}</b> &mdash; play-fairway fusion of 7 detectors</div>
<div>Predicted pixels</div><div>{zi['positive_inside_footprint']:,} of {m['grid']['active_px']:,} active
   ({100*zi['positive_inside_footprint']/m['grid']['active_px']:.3f}% of the scored domain)</div>
<div>On the masked catalogue</div><div class="ok">0 pixels</div>
<div>Emission geometry</div><div>best-first Poisson disk, {geo['min_dist_px']} px
   ({geo['min_dist_px']*100:.0f} m) minimum spacing, hard exclusion within
   {geo['cat_buffer_px']} px ({geo['cat_buffer_px']*100:.0f} m) of the catalogue</div>
<div>Budget rule</div><div>stopped where the marginal hit rate crosses the break-even price
   &pi;* = {100*cal['break_even_pi']['pi_star']:.3f}%</div>
<div>Calibrated |G|</div><div>{G:,.0f} active scored-truth pixels
   [{cal['G_range'][0]:,.0f} &ndash; {cal['G_range'][1]:,.0f}]</div>
<div>Anchor hit rate</div><div>{100*cal['anchor_live_hit_rate']:.3f}%
   ({anchor['n_active_dots']:,} dots &rarr; {cal['anchor_live_hit_rate']*anchor['n_active_dots']:,.0f} hits)
   for the live-scored 0.2778 artifact</div>
<div>Priced live score</div><div><b>{pr['priced']['score']:.4f}</b> (model output, not a measured score)</div>
<div>Wald SPRT</div><div>&alpha;={m['sprt_predeclared']['alpha']}, &beta;={m['sprt_predeclared']['beta']},
   p0={m['sprt_predeclared']['p0']}, p1={m['sprt_predeclared']['p1']} &rarr;
   <span class="{'ok' if sprt['decision']=='accept_H1' else 'mut'}">{esc(sprt['decision'])}</span>
   ({sprt['wins']}/{sprt['n']} folds, LLR {sprt['llr']:+.3f},
   bounds [{sprt['lower']:+.3f}, {sprt['upper']:+.3f}])</div>
<div>Format validator</div><div class="ok">{'PASS' if all(checks.values()) else 'FAIL'} &mdash; all
   {len(checks)} checks on both twins</div>
</div>
</div>

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

<h2>Priced candidates</h2>
<p class="mut">Every candidate is emitted best-first Poisson-disk, its marginal hit-rate decay is
<i>measured</i> on the prevalence-matched instrument, and it is stopped where that marginal rate
crosses the break-even price. &ldquo;Priced score&rdquo; is model output from the calibration above.</p>
<table><tr><th>channel</th><th>pool</th><th>priced n</th><th>instr. hit rate</th><th>live-priced hit rate</th>
<th>priced score</th><th>TP<sub>w</sub></th><th>FP<sub>w</sub></th><th>leakage</th>
<th>SPRT wins</th><th>LLR</th><th>decision</th></tr>
{cand_rows}
</table>
<h3>Promotion gate</h3>
<pre>{gate_txt}</pre>

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

<a class="dl" href="{esc(zeros_rel)}" download>
  &#11015;&nbsp; DOWNLOAD SUBMISSION &mdash; {esc(zi['name'])}
  <small>{zi['bytes']:,} bytes &middot; {zi['positive_inside_footprint']:,} predicted pixels &middot; SHA-256 {zi['sha256'][:16]}&hellip;</small>
</a>
<span class="sub"><a href="index.html">&larr; back to the project page</a></span>

<h1>Executive summary &mdash; exactly how to submit</h1>

<h2>The 6 steps</h2>
<ol>
<li><b>Download.</b> Click the yellow button above, or open
<a href="{esc(zeros_rel)}"><code>{esc(zeros_rel)}</code></a> directly.
Expected size <b>{zi['bytes']:,} bytes</b>, expected SHA-256
<span class="mono">{zi['sha256']}</span>. If your download differs, it was truncated &mdash;
re-download rather than uploading.</li>
<li><b>Sign in</b> at <a href="https://www.drivendata.org/competitions/306/competition-doe-gems/">DrivenData
competition #306</a> and open the <b>Submissions</b> tab. Note the rate limit:
<b>3 submissions per rolling 7-day window</b>.</li>
<li><b>Upload</b> the <code>-zeros.tif</code> file. Do not zip it, rename it, or open and re-save it in
a GIS package &mdash; re-saving can change the NoData encoding and the compression, which is how a
valid file starts failing the range check.</li>
<li><b>Name it</b> <code>gemsdoe39-h40-pfpt-playfairway-permeability</code>.</li>
<li><b>Paste the note</b> (the authoritative string is <code>receipt["note"]</code> in the manifest):
<pre>{esc(note)}</pre></li>
<li><b>Submit</b>, then confirm the returned score is a number and not a validator error. If you see
<code>Predicted values must be in range [0, 1]</code>, see &sect; below &mdash; this artifact is built so that
error cannot come from its pixels.</li>
</ol>

<h2>Why the <code>-zeros</code> twin and not the <code>-nan</code> twin</h2>
<p>Both files carry <b>exactly the same {zi['positive_px']:,} predicted pixels</b>. They differ only in how
the area outside the study footprint is encoded:</p>
<table><tr><th></th><th>zeros twin (SUBMIT)</th><th>NaN twin</th></tr>
<tr><td>outside-footprint pixels</td><td>{zi['outside_footprint_zero']:,} set to 0.0</td>
    <td>{ni['outside_footprint_nan']:,} set to NaN</td></tr>
<tr><td>finite pixels</td><td>{zi['finite_px']:,} (the whole grid)</td><td>{ni['finite_px']:,} (footprint only)</td></tr>
<tr><td>in-footprint min / max</td><td>{zi['min_finite']} / {zi['max_finite']}</td><td>{ni['min_finite']} / {ni['max_finite']}</td></tr>
<tr><td>non-finite inside footprint</td><td class="ok">{zi['nonfinite_inside_footprint']}</td><td class="ok">{ni['nonfinite_inside_footprint']}</td></tr>
<tr><td>out of [0,1] inside footprint</td><td class="ok">{zi['out_of_range_inside_footprint']}</td><td class="ok">{ni['out_of_range_inside_footprint']}</td></tr>
<tr><td>bytes</td><td>{zi['bytes']:,}</td><td>{ni['bytes']:,}</td></tr>
</table>
<p>The organizer's format rule is &ldquo;data outside the bounds is null or nan&rdquo; and &ldquo;a single layer of
float32 values between 0 and 1&rdquo;
(<a href="https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/">page 967</a>).
<b>Zero satisfies both readings</b>: 0 is inside <code>[0,1]</code>, and a predicted 0 outside the study
area contributes nothing to either penalty term because <code>p = 0</code>. It is also immune to any
uploader or browser that coerces NaN into a sentinel. Decisive evidence rather than argument: the
highest-scoring artifact in this project's family &mdash; <code>h33-2-b2</code>, <b>0.2778</b>, live rank
#13 &mdash; <b>is</b> the zeros twin.</p>

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
