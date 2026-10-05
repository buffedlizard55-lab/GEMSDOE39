# Research brief and decision register

**Evidence snapshot:** 2026-10-05 UTC. This is a scientific working document, not a claim of private-test performance. Read with the [ranked hypothesis register and post-hoc formula record](hypothesis-register.md) and [source feed](current-feed.md).

## Executive decision

**Do not spend a weekly DrivenData submission slot on H39X-01.** The candidate passed the local format audit (14/14) and was not an exact pixelwise duplicate of 17 locally available/retrieved TIFFs. Its four-block visible-label proxy comparison produced two wins and two losses; the Bernoulli likelihood-ratio test returned `continue`, not acceptance. More importantly, the exact H39X-01 detector formula was entered into the hypothesis register only after the holdout result. This result is exploratory and cannot be retroactively called preregistered or confirmatory.

The TIFF is a genuine, newly computed artifact for audit and inspection; it is not an approved leaderboard prediction. No upload, public score, private score, or weekly slot consumption occurred.

## Problem framing and official requirements

The task is to predict fault structures relevant to geothermal systems. The official problem page specifies distance-weighted Tversky scoring (α=0.2, β=0.8, 300 m triangular support), and one Float32 GeoTIFF on the supplied EPSG:32611, 100 m grid. Predictions inside the bounds must be in `[0,1]`; cells outside the competition bounds must be null/NoData (represented here as NaN). See the [official problem and submission instructions](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/).

The public leaderboard snapshot checked 2026-10-05 showed 0.3262 at rank 1, 0.3222 at rank 2, and 0.3195 at rank 3; 0.2778 appeared at rank 13. The leaderboard does not show TIFF filenames. The owner-controlled GEMSDOE32 page calls H33-2-B2 “UNSCORED,” so attribution of the public 0.2778 row to that exact file remains unverified. These public values do not identify any private-set score. See [current-feed.md](current-feed.md).

## Candidate, recipe provenance, and observed result

**H39X-01 — geodetic strain / dilatation corridor.** The candidate uses the supplied mirror's `geod_2ndinv`, `geod_shearrate`, and `geod_dilaterate` bands. It combines a robust-scaled strain amplitude, a multi-scale line response, and dilatation sign-transition/gradient support, then emits a sparse 44,090-pixel prediction using best-first Poisson-disk thinning at 2.7-pixel minimum spacing and exact-known-pixel masking. Full formula and limitations are documented in the [hypothesis register](hypothesis-register.md).

> **Protocol deviation:** the exact H39X-01 detector formula/parameters were not in the hypothesis register before holdout scoring; they were recorded after the run. The fixed four-block layout and the used SPRT settings do not make the detector retrospectively preregistered. Treat every result below as exploratory. Do not tune or re-evaluate variants on these exposed folds. A new confirmation requires a genuinely untouched holdout or external new-fault labels and a complete, timestamped plan written before scoring.

The current artifact is `downloads/gemsdoe39-h39x01-strain-corridor-20261005T023209Z-d2f7e140-nan.tif`; checksum and format details are in its adjacent [manifest](downloads/gemsdoe39-h39x01-strain-corridor-20261005T023209Z-d2f7e140-manifest.json). Its pooled score below is an offline DTI proxy, not a DrivenData score.

| Spatial block | H39X-01 proxy DTI | Locked incumbent proxy DTI | Candidate win? |
|---|---:|---:|---|
| NW | 0.010223 | 0.029566 | No |
| NE | 0.000000 | 0.014530 | No |
| SW | 0.069819 | 0.033305 | Yes |
| SE | 0.060368 | 0.022261 | Yes |
| **Pooled** | **0.052944** | **0.024763** | Descriptive proxy only |

The fixed simple-hypothesis Bernoulli SPRT used `p0=0.50`, `p1=0.70`, `alpha=0.05`, `beta=0.10`, upper boundary `ln(1/alpha)=2.995732`, and lower boundary `ln(beta)=−2.302585`. After four outcomes (loss, loss, win, win), the LLR was `−0.348707`, so the result is `continue`. This is neither acceptance nor evidence of a leaderboard improvement. Four blocks provide little power; block-level geologic independence is uncertain, and the nominal error argument is conditional on an independent Bernoulli or valid conditional-supermartingale assumption.

The comparator was a deterministic repository baseline reconstructed from `det_elev`, `tmi`, and `det_elev_slope`, with a matched 44,090-pixel emission budget and the same exact-known-cell mask. It is a local reference, not the best public leaderboard score.

## Spatial proxy design and limitations

- The sample raster's finite mask defines the competition footprint. The visible label raster is used only to construct proxy truth and known-fault masks; these labels are not the private newly identified-fault set.
- Whole 8-connected label components are assigned by centroid to NW/NE/SW/SE quadrants. Each evaluation domain is eroded inward by 15 pixels (1.5 km). A fixed-seed sample of whole components totaling approximately 20% of each fold's label pixels is hidden in that fold; the candidate's exact known-cell mask removes the remaining visible catalogue cells.
- Candidate and incumbent were emitted at equal budgets. Block DTI is diagnostic; pooled counts are combined before calculating the pooled DTI, consistent with organizer staff guidance that leaderboard aggregation is pooled.
- Spatial collars reduce direct proximity leakage but do not prove geological independence. The hidden public-label components are imperfect proxies and may differ in morphology, completeness, and class balance from the competition's private faults.
- Because the exact detector recipe is post-hoc and only four blocks were evaluated, the outcome is not a clean confirmatory performance estimate. Do not make inferential or score forecasts from the pooled value.

## Ranked hypotheses for prospective follow-up

The current ranking below is a research prioritization, not a claim that the four ideas were all time-stamped before H39X-01 was scored. Only H39X-01 was evaluated on the exposed holdout. The other three have no measured DTI or public/private scores. Detailed input layers, physical signatures, novelty, expected direction, relative cost, and multiplicity plan are in [hypothesis-register.md](hypothesis-register.md).

| Rank | Distinct geological signal | Expected DTI direction (hypothesis only) | Status |
|---:|---|---|---|
| 1 | Geodetic second invariant + shear corridor supported by dilatation transition | Up vs topographic/geophysical incumbent if hidden faults localize strain; magnitude unquantified | Exploratory result: 2/4 wins, SPRT continue |
| 2 | GeoDAWN radiometric-ratio alteration front co-located with conductivity transition | Up if alteration and structure jointly improve fault localization; strong lithology/surface confounding | Not tested |
| 3 | Earthquake-intensity/distance edge aligned with a geodetic shear corridor | Up if density transitions track active structures; completeness/smoothing confounding | Not tested |
| 4 | Gravity horizontal/vertical-gradient ridge at a basement-depth inflection | Up if basin-bounding structures are captured; partly redundant with existing potential-field/basement detectors | Not tested |

No numeric improvement or leaderboard score is forecast for any follow-up. Rank reflects expected information gain and distinctness versus implementation/redundancy risk, not measured DTI.

## Multiplicity and future test policy

The current code compared one H39X-01 candidate against one locked comparator, but that does **not** erase post-selection risk: the exact detector was documented after results, and the number of exploratory feature variants considered during development is not fully logged. The current SPRT should therefore be treated as a descriptive screen only. Applying an alpha correction now cannot restore the holdout's untouched status.

For a future confirmatory cycle, choose one of two routes before labels are opened:

1. **Single selected candidate:** select and tune on development data only, then hash/register one final formula, comparator, emission budget, folds, metric, and sequential rule before opening a genuinely new confirmation set. Use the single-test alpha allocation only once.
2. **Four-candidate screening family:** if all four hypotheses are tested on the same fresh holdout, preallocate family-wise α=0.05 as αᵢ=0.0125 per candidate (Bonferroni/union bound) and use the corresponding per-test sequential upper boundary `ln(1/0.0125)=ln(80)≈4.382027`, conditional on each test's own validity assumptions. Do not promote the best-looking candidate from the family without the corrected rule; confirm any selected winner on an additional untouched set. The correction must be declared before the new holdout, and it does not authorize reuse of the current folds.

For either route, stop at the first declared boundary; if fixed folds are exhausted without crossing, record `continue` and stop. No seed search, formula swaps, manual visual selection, or repeated interim result review on the same holdout.

## Data provenance and external source status

The competition data page is login-gated in this environment. The core 19-band raster, labels, sample template, and optional external rasters were restored from owner-controlled GitHub mirrors at pinned commits and matched SHA-256 hashes in `registry/data_manifest.json`. Those hashes verify consistency with those mirrors, **not organizer authenticity**. H39X-01 requires no additional external dataset. The prospective radiometric hypothesis uses existing owner-derived GeoDAWN layers sourced from the official [USGS GeoDAWN data release](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) and [DOI 10.5066/P93LGLVQ](https://doi.org/10.5066/P93LGLVQ); layer alignment/quantization must be independently validated before use.

DrivenData's [Terms of Use](https://www.drivendata.org/termsofuse/) restrict automated site access for monitoring or copying. This repository stores a dated manual snapshot and a live link, not an automated leaderboard scraper.

## Research sources

- [DrivenData competition overview](https://www.drivendata.org/competitions/306/competition-doe-gems/)
- [Official task, metric, and submission format](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/)
- [Official about / geoscience context](https://www.drivendata.org/competitions/306/competition-doe-gems/page/968/)
- [Official leaderboard](https://www.drivendata.org/competitions/306/competition-doe-gems/leaderboard/) — snapshot 2026-10-05
- [National Laboratory of the Rockies official rules PDF](https://www.nlr.gov/docs/fy26osti/96647.pdf) — includes submission limit and AI disclosure requirement
- [DrivenData reference solution](https://github.com/drivendataorg/gems-prize-reference-solution)
- [USGS GeoDAWN release](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) and [DOI](https://doi.org/10.5066/P93LGLVQ)
- [USGS: structural discontinuities and hydrothermal systems in the Great Basin](https://www.usgs.gov/publications/structural-discontinuities-and-their-control-hydrothermal-systems-great-basin-usa)
- [USGS: 3DEP LiDAR data access](https://www.usgs.gov/faqs/what-lidar-data-and-where-can-i-download-it)
- [DrivenData staff clarification: known faults are masked](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2)
- [DrivenData staff clarification: public score aggregates pooled pixels](https://community.drivendata.org/t/leaderboard-aggregation-pooled-over-public-test-pixels-or-mean-of-per-chunk-scores/11550/2)
- [DrivenData Terms of Use](https://www.drivendata.org/termsofuse/)

## Next scientifically valid actions

1. Do not run new variants on the exposed four-block set and do not submit H39X-01.
2. Obtain a new spatially independent label source or reserve a new geographic confirmation region. If none exists, keep current work exploratory and do not claim statistical confirmation.
3. Freeze the next detector formula and its complete candidate family, comparator, budget, holdout construction/order, DTI implementation, exclusions, stopping boundaries, and multiplicity rule in a time-stamped register before scoring.
4. Independently inspect the external GeoDAWN band metadata, spatial alignment, quantization, and use rights before testing H39X-02.
5. Rerun all local code/data/file checks, maintain the mirror-provenance caveat, and include the required AI-use disclosure in any eventual prize submission.
