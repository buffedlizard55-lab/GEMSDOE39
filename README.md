# GEMSDOE39 — GEMS fault discovery, holdout-first

> **Session rule:** Read this entire README before changing the project. Start with the constraints and evidence below; do not treat prior repository prose, owner-reported scores, or old artifacts as verified facts without re-checking them.

## Current deliverable (2026-10-05)

**A new, locally generated, format-verified GeoTIFF is available at** [`docs/downloads/gemsdoe39-h39x01-strain-corridor-20261005T023209Z-d2f7e140-nan.tif`](docs/downloads/gemsdoe39-h39x01-strain-corridor-20261005T023209Z-d2f7e140-nan.tif). It is a fresh geodetic-strain/dilatation candidate; no prior prediction pixels were copied. Its audit manifest is next to the TIFF. The [live GEMS download site](https://buffedlizard55-lab.github.io/GEMSDOE39/docs/) (source file: `docs/index.html`) has the prominent download button; the [executive submission guide](docs/executive-summary.html) explains format, upload steps, and why this candidate is not currently approved.

- 44,090 positive pixels; single-band float32; EPSG:32611; 3,730×3,292; 100 m; exact sample geotransform.
- Every in-footprint value is finite and in `[0,1]`; all 7,111,787 cells outside the sample footprint are NaN; nodata is NaN.
- SHA-256: `a612f66fe640cc6d5dd2333a330a454c2d8affef04f4672961d07bd4415a0b83`.
- Independent validator: **14/14 checks passed**.
- Pixel-identity check against 17 locally available prior TIFFs (including the downloaded H33-2-B2 reference): **no exact duplicate**. This does not prove uniqueness against every file hosted by every sibling site.
- **Do not spend a weekly competition submission slot on this candidate.** The four-block spatial proxy comparison produced two wins and two losses. Although fixed SPRT settings were used, the exact H39X-01 detector recipe was recorded in the hypothesis register only after the holdout result; this is exploratory, not preregistered or confirmatory evidence. Pooled proxy DTI was 0.052944 vs 0.024763 for the locked incumbent, but the sequential likelihood-ratio test returned `continue`, not acceptance. Local format-validity does not establish private-test improvement. No competition score is claimed.

This project uses the official file format: data outside the sample bounds are NaN/null. The legacy zeros-outside TIFFs remain only as historical artifacts; they are not the recommended download and the current validator rejects them as not matching the published outside-bounds requirement.

### Public leaderboard snapshot

The official leaderboard page was checked on **2026-10-05**. At that snapshot, the top public score was **0.3262** (`nchuzhoy`), followed by 0.3222 and 0.3195. Thus 0.3195 is not currently the leader in the fetched page. A 0.2778 leaderboard row was visible at rank 13, but the leaderboard does not associate that row with a particular GeoTIFF filename. The GEMSDOE32 site itself labels its H33-2-B2 download “UNSCORED”; exact attribution of 0.2778 to that file is therefore unresolved. See [the dated feed/evidence register](docs/current-feed.md) and the [live leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/).

## Mission and mandatory operating constraints

The mission is to improve the probability of winning the DOE GEMS Prize (#306) by generating scientifically defensible, reproducible fault predictions and an obvious-to-download, valid GeoTIFF. No code or copy can guarantee a score above 0.3262, because the private expert labels are unavailable.

1. **Unique artifact.** Generate each candidate from the data and documented transforms; never copy a previous submission as the deliverable. Use previous artifacts only for educational analysis, baseline comparison, or duplicate checks, and disclose that use.
2. **Format before score.** Validate CRS, shape, transform, dtype, footprint, finite/range constraints, outside-NaN behavior, and SHA-256 after writing. DrivenData requires one float32 layer, EPSG:32611, 100 m, matching bounds, values `[0,1]`, and null/NaN outside bounds.
3. **Holdout before any weekly slot.** Preregister one candidate, its transforms, budget, fold construction, comparator, metric, and stopping rule before examining its holdout result. Do not upload an untested idea. The current candidate is explicitly **not approved** by its gate.
4. **SPRT, not peeking.** Compare one preregistered candidate with one locked incumbent, one block at a time. Use fixed fold order, predeclared `p0=0.50`, `p1=0.70`, `alpha=0.05`, `beta=0.10`, and conservative likelihood-ratio thresholds `upper=ln(1/alpha)`, `lower=ln(beta)`. Stop at a boundary. If fixed folds are exhausted without crossing, record `continue` and stop this experiment; do not keep re-scoring the same holdout with new seeds.
5. **State assumptions.** The SPRT's nominal error control is conditional on independent Bernoulli block outcomes or a valid conditional supermartingale. A spatial collar reduces leakage but does not prove geologic independence. Do not claim unconditional error control where this assumption is doubtful.
6. **Separate multiplicity from sequential stopping.** Test only the single registered candidate in a holdout cycle. Do not pick a winner from many hypotheses that share the same holdout. Future portfolios need a separate multiplicity policy (for example, alpha allocation or a genuinely untouched confirmation set).
7. **No unsupported causal claims.** Public leaderboard scores, owner-reported historical scores, proxy-holdout DTI, model estimates, and private-test scores are different evidence classes. Label them correctly. A proxy DTI is not a leaderboard prediction.
8. **Sources and data provenance.** Prefer official DrivenData, USGS, DOE/NLR, and peer-reviewed/USGS research links. The challenge data page is login-gated in this environment; restored owner mirrors are SHA-256 checked but are not organizer-authenticated. Store that limitation with every data claim.
9. **Respect platform terms.** DrivenData's Terms of Use prohibit robot/spider or other automatic access for monitoring/copying. Provide a live leaderboard link and timestamped snapshots; do not deploy an automated DrivenData scraper. The current page can be rechecked in a future session through an allowed human/agent review.
10. **Core values:** **Maximize P(Win)** by evidence-based choices; **Own the Outcome** by fixing defects and publishing limitations rather than hiding them. Generative AI use must be disclosed in the competition narrative as required by the official rules.

## Protocol deviation — H39X-01 is exploratory, not preregistered

The exact H39X-01 feature recipe was entered into `docs/hypothesis-register.md` only after the four-block result was observed. Fixed SPRT parameters and a blocked fold design do not retroactively preregister the candidate formula. Therefore the scores and `continue` result are exploratory only; this test cannot confirm the hypothesis, and the artifact is not eligible for a weekly slot regardless of a later rerun on the same folds. Do not tune on this exposed holdout. A future confirmatory experiment must freeze the full formula and analysis plan before opening a genuinely untouched holdout. See the dated [hypothesis register](docs/hypothesis-register.md) and [research brief](docs/research-brief.md).

## Restore data and run local checks

Competition data are not committed to Git. The present sandbox restored the core and external rasters from the SHA-pinned owner mirrors. On a new machine:

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python scripts/restore_data.py --group core
python scripts/restore_data.py --group external
python scripts/prepare_data.py
python scripts/validate_submission.py docs/downloads/gemsdoe39-h39x01-strain-corridor-20261005T023209Z-d2f7e140-nan.tif
pytest -q
```

The H39X-01 spatial folds are exhausted and cannot be used to test new variants. `scripts/build_pipeline.py` refuses to score them unless the explicit `--replay-exposed-holdout` audit flag is supplied; that deterministic replay still can never approve a submission. Register future formulas before evaluating a genuinely untouched holdout. `restore_data.py` can also restore the historical artifacts for educational comparison with `--group scored`; they are not inputs to the H39X-01 detector. The final upload still requires the competitor's own DrivenData account; this project does not log in, submit a file, or consume a weekly slot.

## Current hypothesis and holdout result

The highest-ranked current research candidate is **H39X-01: strain-localized, dilatational fault corridor**. It combines robust-scaled geodetic second-invariant and shear-rate amplitudes, multi-scale invariant line response, and paired/sign-transition dilatation support. It uses `geod_2ndinv`, `geod_shearrate`, and `geod_dilaterate`; the locked local incumbent comparator uses detrended elevation, total magnetic intensity, and detrended-elevation slope. Its geodetic-strain input family differs from the checked-in H39-A…E detector families. **The exact detector formula was documented only after holdout scoring**, so it must not be called preregistered or confirmatory. See the post-hoc record in `docs/hypothesis-register.md`.

| Spatial block | H39X-01 DTI proxy | Incumbent DTI proxy | Candidate win? |
|---|---:|---:|---|
| NW | 0.010223 | 0.029566 | No |
| NE | 0.000000 | 0.014530 | No |
| SW | 0.069819 | 0.033305 | Yes |
| SE | 0.060368 | 0.022261 | Yes |
| **Pooled** | **0.052944** | **0.024763** | descriptive proxy only |

The sequential log-likelihood ratio after the four fixed blocks is `−0.348707`, between the lower boundary `−2.302585` and upper boundary `2.995732`. The correct statistical outcome is **continue / inconclusive**, not acceptance. Four blocks provide little power, and geologic dependence weakens the Bernoulli-independence assumption. Because the exact candidate recipe was recorded post-run, this outcome is exploratory and the candidate is not eligible for a leaderboard slot.

See the [ranked hypothesis register and post-hoc formula record](docs/hypothesis-register.md), [audit manifest](docs/downloads/gemsdoe39-h39x01-strain-corridor-20261005T023209Z-d2f7e140-manifest.json), and [dated evidence feed](docs/current-feed.md).

## Research hypotheses not run this cycle

These are ranked research ideas, not measured performance claims. Only H39X-01 was tested; the other candidates were not run against the same holdout.

1. **Strain-localized dilatational corridor** — provided geodetic strain-rate/invariant layers; coherent high-shear band plus a bipolar dilatation transition. May expose active blind structures absent from surface catalogues; no new external source required. Tested this cycle; inconclusive by SPRT.
2. **Radiometric alteration front × conductivity edge** — GeoDAWN K/Th/U ratios and `cond_surf`; co-located radiometric ratio gradients and an electrical-conductivity transition. Could reflect hydrothermal alteration but is non-unique and near-surface confounded. Official USGS GeoDAWN DOI is [10.5066/P93LGLVQ](https://doi.org/10.5066/P93LGLVQ), checked obtainable 2026-10-05; no new external data required for a first test because a derived layer already exists in the mirror.
3. **Earthquake-density edge aligned with strain** — provided `ieq_n100a15`, `deq_n100a15`, `geod_shearrate`; lineament or density transition rather than a raw seismicity hotspot. May expose active faults, but completeness, catalog smoothing, and unrelated events are major risks.
4. **Gravity derivative / basement inflection** — provided `iso_grav_anom_hg`, `iso_grav_anom_vg`, `depth_to_base_surf`; co-located horizontal-gradient ridge and vertical-gradient sign transition at basement-depth inflection. May detect buried basin-bounding faults but overlaps existing gravity/basement evidence.

Full layer-level description, physical signature, novelty, qualitative expected DTI direction, relative implementation cost, and future multiplicity controls are in `docs/hypothesis-register.md`. It explicitly marks the H39X-01 formula as post-hoc rather than a preregistration.

## What the historical 0.2778 result does — and does not — tell us

The challenge metric is a distance-weighted Tversky index with α=0.2 for false positives, β=0.8 for false negatives, and a 300 m triangular distance kernel. It rewards predictions near truth and penalizes broad, unsupported mass; a sparse, well-ranked dot set can therefore outperform an over-thick map. The GEMSDOE32 site describes H33-2-B2 as a 37,654-dot flank-B=2-pruned file based on a dotted-ridge family. This is a plausible mechanism for a high public score, but the same page labels that specific file “UNSCORED.” The public leaderboard displays a 0.2778 row but not a prediction filename. Therefore the file-to-score attribution, and the causal explanation for that particular score, cannot be verified from these pages alone.

The maximum in the current public snapshot is 0.3262, so 0.2778 is not the current target. A score above either number cannot be inferred from the historical series or a proxy holdout. The candidate generated here has no DrivenData score.

## Verified sources and knowledge register

- [Competition overview](https://www.drivendata.org/competitions/306/competition-doe-gems/)
- [Problem description, DTI metric, datasets, and submission format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [About page](https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/)
- [Competition data page](https://www.drivendata.org/competitions/306/competition-doe-gems/data/) (redirects to login)
- [Official leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/) (timestamped snapshot 2026-10-05)
- [Official rules PDF — National Laboratory of the Rockies (NLR)](https://www.nlr.gov/docs/fy26osti/96647.pdf)
- [DrivenData reference solution](https://github.com/drivendataorg/gems-prize-reference-solution)
- [USGS GeoDAWN data release](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) · [DOI / data](https://doi.org/10.5066/P93LGLVQ)
- [USGS LiDAR download FAQ / 3DEP entry points](https://www.usgs.gov/faqs/what-lidar-data-and-where-can-i-download-it)
- [USGS: structural discontinuities and hydrothermal systems in the Great Basin](https://www.usgs.gov/publications/structural-discontinuities-and-their-control-hydrothermal-systems-great-basin-usa)
- [USGS: 3D geologic mapping and geothermal potential in Nevada/Oregon](https://www.usgs.gov/publications/three-dimensional-geologic-mapping-assess-geothermal-potential-examples-nevada-and)
- [DrivenData staff clarification: known-fault pixels are masked](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2)
- [DrivenData staff clarification: public score pools all public-subset pixels](https://community.drivendata.org/t/leaderboard-aggregation-pooled-over-public-test-pixels-or-mean-of-per-chunk-scores/11550/2)
- [DrivenData Terms of Use](https://www.drivendata.org/termsofuse/) (prohibits automated access/monitoring)

The line-by-line source notes, leaderboard snapshot, provenance flags, and irregularities are in [docs/current-feed.md](docs/current-feed.md). The detailed scientific decision register is [docs/research-brief.md](docs/research-brief.md).

## Known limitations / review flags

1. **No private labels:** local CV uses visible catalogue faults as imperfect proxies; it cannot predict the private expert labels or Phase 2 expert review.
2. **SPRT independence:** the four large spatial blocks have an inward 1.5 km collar, but remaining spatial dependence cannot be excluded. The nominal error guarantee is conditional on an independent-block or valid conditional-supermartingale model; current result did not cross either boundary. The detector recipe itself was not registered before evaluation.
3. **No submission authentication:** DrivenData requires an account. This agent generated and validated a file but did not upload it or use a weekly slot.
4. **Input provenance:** the competition data page redirects to login. Files here were restored from owner-controlled GitHub mirrors with SHA-256 pins, not downloaded directly from the organizer. Hashes certify mirror integrity, not organizer authenticity.
5. **Leaderboard feed:** DrivenData Terms prohibit automated access for monitoring/copying; the page is linked and a dated snapshot is stored, not auto-scraped.
6. **Historical score mapping:** user-supplied scores and owner pages do not always agree. In particular, the GEMSDOE32 page calls H33-2-B2 unscored while a 0.2778 leaderboard row exists. Do not claim that score belongs to that exact raster without an official receipt.
7. **Artifact scope:** duplicate testing covered 17 local/retrieved TIFFs, not every public artifact from all GEMSDOE sites. The detector is newly computed and its pixel set differs from the retrieved H33 reference, but the scope limit is explicit.
8. **Model training:** this run is a deterministic, unsupervised geological feature transform, not a trained neural network. Training on the visible catalogue can teach mapped-fault morphology but cannot directly supervise “new fault” truth. The DrivenData reference U-Net needs compute/data and still cannot remove that label-shift problem.
9. **AI disclosure:** official rules allow generative AI but require the competitor to describe its extent/use in the narrative and remain responsible for every claim. Include that disclosure in any prize package.

## Three-pass implementation and review record

1. **Pass 1 — build:** refactor the feature reader, candidate pipeline, spatial holdout/SPRT path, NaN-outside GeoTIFF writer, strict validator, and first public-facing documentation. The H39X-01 formula was not recorded in the hypothesis register before scoring; the result is therefore explicitly exploratory.
2. **Pass 2 — independent review:** inspect protocol leakage, post-hoc wording, range/footprint behavior, uniqueness scope, score attribution, restore/preparation helpers, and download/guide links; correct issues rather than promoting the candidate.
3. **Pass 3 — final re-check:** rerun unit/integration tests and smoke checks for restore receipts, source inspection, baseline writing and validation; re-audit the candidate checksum/manifest, all local download links, and the final diff against the charter. Do not consume a slot; the SPRT result remains `continue`.

## Historical score snapshots supplied by the project owner

These values are preserved for longitudinal context, **not independently verified as exact file/score pairs**. The current organizer leaderboard snapshot takes precedence for live rank.

| Site / project | Score(s) reported in the supplied project history |
|---|---|
| GEMSDOE | 0.1563 |
| 6GEMSDOE | 0.0286 |
| GEMSDOE3 | 0.1193, 0.0830, 0.1152 |
| GEMSDOE2 | 0.1560 |
| GEMSDOE4 | 0.0343 |
| 5GEMSDOE | 0.1563 |
| 7GEMSDOE | 0.1461 |
| 8GEMSDOE | 0.1563 |
| GEMSDOE9 | 0.0107 |
| 11GEMSDOE | 0.0202 |
| 12GEMSDOE | 0.1294 |
| 15GEMSDOE | 0.0782 |
| 14GEMSDOE | 0.0020 |
| 17GEMSDOE | 0.0187 |
| 18GEMSDOE | 0.0297 |
| GEMSDOE19 | 0.1894, 0.1922 |
| GEMSDOE10 | 0.0461, 0.0921, 0.1280, 0.1839 |
| 13GEMSDOE | 0.0904 |
| 16GEMSDOE | 0.1855, 0.0976, 0.0360 |
| GEMSDOE21 | 0.1894 |
| 20GEMSDOE | 0.1890, 0.1859 |
| GEMSDOE22 | 0.1002, 0.0748 |
| GEMSDOE23 | 0.1352 |
| GEMSDOE24 | 0.2477 |
| GEMSDOE25 | 0.2600 |
| GEMSDOE26 | 0.1223 |
| GEMSDOE27 | 0.2449 |
| GEMSDOE28 | 0.2708 |
| GEMSDOE29 | 0.2600; several other candidates listed without scores |
| GEMSDOE30 | 0.2600 |
| GEMSDOE31 | 0.2708 |
| GEMSDOE32 | 0.2778 (file attribution not verified) |
| GEMSDOE33–38 | no score supplied for listed candidates |
| GEMSDOE39 / 40GEMSDOE | no prior score supplied |
