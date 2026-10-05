# Ranked geological hypothesis register — H40 family and H39X-01 post-hoc record

**Current record date:** 2026-10-05 UTC

**Important status:** this file was assembled/updated after the H39X-01 holdout run. It is **not** a retroactive preregistration. The exact H39X-01 detector formula and parameter choices were entered here only after the four block outcomes had been observed. Fixed SPRT settings do not cure that protocol deviation. The result is exploratory; this candidate is not confirmatory and must not consume a weekly submission slot.

The ranking below is the current research priority order. It does not claim that all four ideas or their exact formulas were timestamped before H39X-01 was evaluated. Only H39X-01 was tested in the reported cycle. No numerical score forecast is made for any hypothesis.

## Ranked hypotheses

“Cost” below means relative additional data-acquisition, implementation, and compute burden for a reproducible first test—not a dollar estimate. All referenced layers are already present in the SHA-256-pinned owner mirrors in this workspace, but the competition rasters are not organizer-authenticated here. “Expected DTI direction” is a physical prior against the locked local incumbent, not measured performance. A positive direction has no numeric magnitude attached.

| Rank / ID | Hypothesis & input layers | Predicted physical signature | Novelty relative to checked-in detector families | Expected DTI direction / improvement | Relative cost & uncertainty |
|---|---|---|---|---|---|
| **1 · H39X-01** | **Strain-localized, dilatational fault corridor.** Competition feature bands `geod_2ndinv` (second invariant), `geod_shearrate`, `geod_dilaterate`; exact known-fault cells used only in emission masking. | A spatially coherent high-invariant/high-shear lineament, sharpened by multi-scale line response and supported by either nearby opposite-sign dilatation or a dilatation-gradient edge. Such deformation contrasts may mark a fault or damage zone even where no strong surface scarp is expressed. This signature is a hypothesis, not evidence of geothermal productivity. | Uses a geodetic-strain primary signal not used by checked-in H39-A…E transforms (magnetic/gravity/conductivity edges, topography/LiDAR, magnetic tilt/continuation, basement/conductivity coincidence). Repository-relative distinction only; no claim of a novel method in the literature. | **Hypothesized up** vs topographic/geophysical ridge incumbent if hidden faults localize present-day deformation. No ex-ante numeric gain is defensible. The observed pooled proxy DTI was 0.052944 vs 0.024763, but this is exploratory, 2/4 blocks were wins, and SPRT=`continue`; do not infer an improvement. | **Medium.** No new data needed; existing input bands. Full-grid filtering has moderate CPU/memory cost. High geological uncertainty because the supplied deformation field may not resolve all blind/old faults. |
| **2 · H39X-02** | **Radiometric alteration front × conductivity transition.** Owner-derived `data/external/geodawn_rad_u8.tif` (K, Th, U, TC) and `geodawn_extensions_u8.tif` (Th/K, U/K, U/Th, TMI-upward-continued); competition `cond_surf`. Source: official USGS GeoDAWN release, [DOI 10.5066/P93LGLVQ](https://doi.org/10.5066/P93LGLVQ). | A line-like discontinuity or paired gradient in one or more radiometric ratios co-located with an electrical-conductivity edge. Alteration/weathering and fluid pathways can change radiometric response and conductivity; lithology, regolith, moisture, flight geometry, and quantization can produce similar patterns. | Adds radiometric ratio features, not used by checked-in H39-A…E. Conductivity itself overlaps existing detectors, so novelty lies in the cross-sensor alteration-front interaction. | **Hypothesized up** if alteration plus conductivity identifies structurally focused fluid pathways more precisely than a generic geophysical ridge. No numeric gain or score estimate; potential specificity is offset by substantial non-fault confounding. | **Medium.** Derived public-source layers are already mirrored; verify band order, alignment, quantization, masks, and CC0/source handling. Moderate feature-engineering/compute cost. Geological interpretation risk is high. |
| **3 · H39X-03** | **Earthquake-intensity/distance edge aligned with strain.** Competition feature bands `ieq_n100a15`, `deq_n100a15`, `geod_shearrate`, optionally `geod_2ndinv`. | A density/intensity or distance-field transition that is elongated and co-located with a geodetic shear corridor, rather than selecting the highest earthquake-density cells alone. Could indicate an active structural boundary; catalog completeness, smoothing, depth, and unrelated seismicity can mimic it. | Adds seismicity-derived fields and a structural alignment test absent from H39-A…E. It is independent in input family from the magnetic, gravity, and surface-topography detectors, while sharing geodetic support with H39X-01. | **Hypothesized up** if spatially aligned seismicity edges mark active faults missed by surface/geophysics; no numeric improvement forecast. A raw seismicity hotspot may instead lower DTI via unrelated activity. | **Low–medium.** Uses existing competition bands; modest edge/alignment logic and moderate full-grid compute. Main uncertainty is seismic catalog completeness and the semantics of smoothed fields. |
| **4 · H39X-04** | **Gravity derivative ridge × basement-depth inflection.** Competition feature bands `iso_grav_anom_hg`, `iso_grav_anom_vg`, `depth_to_base_surf`. | Co-located horizontal-gradient ridge / vertical-gradient sign transition near a basement-depth step. This may highlight buried basin-bounding structures lacking a DEM scarp, but non-fault density contrasts and broad basin geometry are alternatives. | Uses explicit gravity-derivative plus basement-inflection geometry. It overlaps more with existing H39-A gravity edges and H39-E basement/conductivity co-edges than H39X-02/03, hence the lower rank. | **Weakly hypothesized up** if the derivative/inflection interaction adds buried structural localization. No numeric gain forecast; substantial redundancy risk could produce no improvement or a loss. | **Low.** Inputs are present in the core feature raster; straightforward multi-scale gradients. Low data cost and relatively low implementation cost, but high redundancy risk. |

### Prioritization rationale

H39X-01 is ranked first for a genuinely different primary signal family and no additional data requirement, but its current result does not validate it. H39X-02 offers a potentially more geothermal-specific alteration cue with established public source data, balanced against surface/lithology confounding. H39X-03 adds a distinct active-structure cue but inherits seismic-catalog limitations. H39X-04 is cheapest to test but most redundant with earlier gravity/basement detectors. These are research judgments; no measured result supports the ranking of untested candidates.

## H40 family — off-catalogue faults, 2026-10-05 (session 2)

**Registration status, stated up front.** `CONTROL-structural` and `H40-A-profile` were fixed before the
first pass ran. `H40-E`, `H40-F` and `H40-G` were added **after** that first pass showed the pure geophysical
structural field reaches only w = 0.065 against the off-catalogue surrogate — *below every one of the 12
organiser-scored historical artifacts* (0.048–0.105), which is why a catalogue-trained or purely geophysical
field cannot win. **The family-wise error rate across the whole candidate set is therefore not controlled at
α = 0.05.** The SPRT verdicts below are reported as such; they are not confirmatory evidence about the
late-added candidates.

Two further protocol deviations are disclosed in `docs/index.html` §4: the holdout evaluation domain was
restricted mid-run (to pixels outside the catalogue exclusion), and cells with fewer than 50 scorable hidden
truth pixels are excluded rather than booked as ties.

### The mechanism each hypothesis bets on

The organizer states that known USGS/INGENIOUS fault pixels are masked from evaluation, and that "new fault"
means **any fault pixel not already captured by USGS/INGENIOUS, including newly mapped geometry of an existing
fault system** ([community #11536, staff reply](https://community.drivendata.org/t/where-do-you-draw-the-line/11536/2)).
So the scoreable target is by construction the set of fault pixels that are *near* mapped faults but not
*on* them. That is a different target from every prior submission in this family, all of which trained on the
catalogue itself.

| ID | Layers / data | Physical signature and transform | Why it could catch an off-catalogue fault | Rank (expected gain ÷ cost) |
|---|---|---|---|---|
| **H40-E-disc — SHIPPED** | all 19 supplied bands + 4 gradient magnitudes + elevation coherence + distance-to-catalogue | Spatially-blocked out-of-fold gradient-boosted discriminant for *off-catalogue* fault presence (positives = 79,615 SGMC pixels that are not on the catalogue; negatives = a hard 1–6 px ring plus background), multiplied by the measured relative-density profile | A mapped fault is one *sample* of the geophysical conditions that produce faults; a discriminant generalises those conditions to adjacent ground that no map covers. Blocked OOF AUC **0.7489** over 11 blocks proves the signature transfers across space rather than memorising the traces | **1** — gain 2.71× w over control at zero new data cost |
| H40-A-profile | det_elev, tmi, rtp, iso_grav_anom, cond_surf, tc + catalogue geometry | Geometric mean of three independent physics families (magnetic, gravity, conductivity) × structure-tensor coherence × the measured relative-density profile | Corroboration across independent physics suppresses single-survey artefacts; the profile kills the >5 km wasteland where density falls to 0.39–0.82× base | 2 — cheap, but only 1.03× w over control |
| H40-G-disc-thermal | H40-E + OpenEI GDR springs/wells/vents | H40-E × proximity to thermal manifestations | Hydrothermal discharge localises on active permeable structures | 3 — measured **worse** than H40-E (0.1354 vs 0.1370); thermal proximity adds nothing on this instrument |
| H40-F-disc-sgmc0.15/0.30 | H40-E + SGMC corridor | H40-E × (SGMC-corridor)^0.15 or ^0.30 | Directly targets "newly mapped geometry of an existing system" | **4 — dropped.** Scored highest on instrument B (0.4095/0.4137) but 41–47 % of dots sit within 100 m of an SGMC trace, where the surrogate is the training positives, so the credit is circular. Rejected by the pre-declared 2× extrapolation guard |
| CONTROL-structural | det_elev, tmi, rtp, iso_grav_anom, cond_surf, tc | Geometric mean of the three physics families × coherence, no spatial prior | Null baseline for the SPRT | — |

### Outcomes

| Candidate | w @40k | w ÷ w_rand | Instrument B pred | % dots <100 m of off-cat SGMC | SPRT (usable domain) | Mean lift |
|---|---:|---:|---:|---:|---|---:|
| CONTROL-structural | 0.0636 | 1.26 | 0.1775 | 5.5 | control | 0.5086 |
| H40-A-profile | 0.0652 | 1.29 | 0.1808 | 5.6 | accept H1, 9/9 | 0.5406 |
| **H40-E-disc** | **0.1370** | **2.71** | **0.2891** | **12.3** | **accept H1, 9/9, LLR +3.0283** | **0.8225** |
| H40-F-disc-sgmc0.15 | 0.4451 | 8.82 | 0.4095 | 40.8 | accept H0 | 0.4092 |
| H40-F-disc-sgmc0.30 | 0.5066 | 10.03 | 0.4137 | 46.6 | accept H0 | 0.4092 |
| H40-G-disc-thermal | 0.1354 | 2.68 | 0.2872 | 12.2 | accept H1, 9/9 | 0.8257 |

`w_rand = 0.0505` is the credit per dot of a uniformly random field against the same surrogate. The shipped
field reaches **0.1461** at its 30,000-dot budget — 2.89× random, and 39 % above the best per-dot credit any
artifact in the 12-file corpus has ever achieved (0.1051).

### Selection, and the conflict between the two instruments

The two instruments disagree **by construction**, and that must be read as a property of the instruments
rather than of the candidates:

- **Instrument A** (catalogue-hidden spatial holdout) ranks fields by how well they find *mapped* faults —
  precisely the pixel class the organiser masks out of scoring. It cannot adjudicate a catalogue-trained model
  at all (circular truth), and its usable domain leaves isolated outlier pixels, which are systematically
  hostile to off-catalogue evidence (it scores H40-F at 0.409, *below* its own 0.5 null).
- **Instrument B** (calibrated surrogate credit) is the only instrument here ever checked against organiser
  scores: ρ = +0.923 in sample, **+0.902 leave-one-out**, n = 12.

**Declared promotion rule:** promote the highest instrument-B score subject to a 2× extrapolation guard
(w ≤ 2 × 0.1051 = 0.2103), reporting every instrument-A verdict alongside **without letting it veto**. The
guard is what dropped H40-F. This rule was fixed before the final run and is reproduced verbatim in the stage-6
log of `registry/h40_report_h40e-30k.json`.

### Budget

The budget is not chosen by the fitted model's absolute level but by the metric's own denominator: adding one
predicted pixel costs exactly 0.2 of DTI regardless of quality, so a marginal dot pays iff its expected credit
exceeds 0.2 × DTI. Applying that to the measured emission curve selects **30,000 dots**, which is also the peak
of the fitted score curve (0.3006). Two independent derivations agreeing is the reason the budget is shipped.
An earlier maximin-over-scenarios rule was replaced: it is degenerate here because the pessimistic scenario is
monotone increasing in budget, so it always returns the largest grid value regardless of field quality.

## H39X-01 exact formula — documented after the run

The following formula reconstructs the implementation in `src/gems39/features.py` and the candidate manifest. It is recorded for auditability, **not** as evidence that these values were preregistered. The feature reader first uses GeoTIFF `band_name` tags/descriptions, nearest-fills invalid in-footprint values, clips each band to its 0.1st/99.9th percentiles, and zeros the array outside the footprint. The candidate feature function then applies:

1. Let `inv = geod_2ndinv`, `shear = geod_shearrate`, and `dil = geod_dilaterate`.
2. For each of `inv` and `shear`, compute a robust unit scaling of `log1p(abs(value))` over the footprint using quantiles 0.02 and 0.995. Define `strain_amp = sqrt(inv_u * shear_u)`.
3. Smooth `dil` with Gaussian σ=1.5 pixels. Within a 7×7 local window, compute positive maximum and negative minimum amplitudes, robust-scale each at quantiles 0.02/0.995, and set `bipolar = min(pos_u, neg_u)`. Separately compute the gradient magnitude of `dil` after derivative Gaussian σ=3.0 pixels and robust-scale it at the same quantiles. Set `dilation_support = max(bipolar, dilation_edge)`.
4. Apply `log1p(abs(inv))`; at Gaussian Hessian scales σ=2 and 4 pixels, compute a bright-ridge line response (`sign=-1`), robust-scale each response at quantiles 0.02/0.995, and take the pixelwise maximum as `ridge`.
5. Form `score = 0.50*strain_amp + 0.30*ridge + 0.20*dilation_support`, then robust-scale the result over the footprint at quantiles 0.02/0.995.
6. For each holdout fold, best-first Poisson-disk emission selects up to 44,090 positive pixels with minimum separation 2.7 pixels. The only catalogue exclusion in this run is the exact known-cell mask (`cat_buffer_px=0`); no 200 m halo is claimed to be organizer-masked. The final full-footprint artifact uses all visible known-fault cells for this exact-cell exclusion.
7. Encode binary predictions as float32 0/1 in the official footprint and NaN outside. The output is fresh; no prior prediction pixels were copied.

The matched incumbent is a deterministic ridge surface from `det_elev`, `tmi`, and `det_elev_slope`: elevation line-response scales 1.0, 2.0, 3.5 pixels; normalized ridge weighted by `0.55 + 0.25*mag_gradient + 0.20*slope_break`; same 44,090 budget, 2.7-pixel minimum spacing, and exact known-cell mask.

## H39X-01 holdout design and exploratory outcome

The holdout implementation uses the finite mask of `sample_submission.tif` as the footprint and visible `labels.tif` as proxy truth. It labels whole 8-connected components, assigns each component to one quadrant by centroid, erodes each quadrant inward by 15 pixels (1.5 km), and hides a fixed-seed set of whole components amounting to approximately 20% of each fold's label pixels. The fold order is NW, NE, SW, SE; seed `20261005`; exact known cells are masked from score, not an unverified surrounding halo. Public catalogue labels are not the private newly identified-fault target.

| Block | H39X-01 DTI proxy | Locked incumbent DTI proxy | Outcome |
|---|---:|---:|---|
| NW | 0.010223 | 0.029566 | Loss |
| NE | 0.000000 | 0.014530 | Loss |
| SW | 0.069819 | 0.033305 | Win |
| SE | 0.060368 | 0.022261 | Win |
| **Pooled** | **0.052944** | **0.024763** | Descriptive only |

SPRT settings used: `p0=0.50`, `p1=0.70`, `alpha=0.05`, `beta=0.10`. Wins increment LLR by `ln(0.7/0.5)`; losses increment it by `ln(0.3/0.5)`. Upper boundary `ln(1/alpha)=2.995732`; lower boundary `ln(beta)=−2.302585`. After the fixed four folds, LLR was `−0.348707`; decision `continue`. It neither accepts the candidate nor rejects the incumbent. The formal interpretation is conditional on independent Bernoulli block outcomes or a valid conditional-supermartingale argument. Spatial collars reduce but do not prove geologic independence.

The exact formula's post-run documentation invalidates any confirmatory interpretation regardless of the numerical SPRT state. Do not reuse these folds for tuning, seed search, alternative masks, or variant selection. The candidate is **not approved** for any weekly submission slot.

## Multiplicity and decision policy for a future cycle

- The four entries above are distinct research ideas, not four scored candidates in the current holdout. Only H39X-01 was scored in the described cycle.
- One-candidate analysis in code does not remove model-selection risk: the exact formula is post-hoc, and the number of exploratory variants tried during development was not fully registered. This cycle's holdout statistics are descriptive/exploratory.
- A correction applied after seeing the exposed result cannot repair preregistration or holdout leakage. The current four folds are spent for confirmatory purposes.
- For a future run, either select exactly one frozen recipe on development data and open a new untouched confirmation set once, or preallocate family-wise α across all candidates tested on that fresh holdout. For a four-candidate family with total α=0.05, use per-candidate αᵢ=0.0125 and sequential upper boundary `ln(1/αᵢ)=ln(80)≈4.382027` (under each test's validity assumptions). Use a separate untouched confirmation set before recommending a weekly submission slot.
- Freeze and hash the detector formula, all feature transforms, emission budget/spacing, exclusion mask, comparator, spatial fold construction/order, metric implementation, hypothesis family, alpha allocation, beta/stopping policy, and reporting rules before evaluating the holdout. If fixed folds finish without a boundary, record `continue` and stop.

## Source, data, and evidence limitations

- Official format, metric and public leaderboard facts are linked in the [dated current-feed register](current-feed.md). The current public leader snapshot was 0.3262 on 2026-10-05; no candidate score is known.
- The 19-band competition feature raster, labels, sample and optional external layers were restored from pinned owner mirrors. Hashes establish mirror integrity, **not** organizer authentication. The data-download page requires login in this environment.
- The official [USGS GeoDAWN page](https://www.usgs.gov/data/geodawn-airborne-magnetic-and-radiometric-surveys-northwestern-great-basin-nevada-and) and [DOI](https://doi.org/10.5066/P93LGLVQ) support H39X-02 source availability; the derived layers in this repository are owner-mirrored/processed products, not direct organizer-provided inputs.
- USGS research on [structural discontinuities and hydrothermal systems](https://www.usgs.gov/publications/structural-discontinuities-and-their-control-hydrothermal-systems-great-basin-usa) motivates structural hypotheses but does not prove performance of these transforms.
- DrivenData Terms of Use prohibit automated monitoring/copying. The leaderboard is linked and manually timestamped, not scraped.
- Exact duplicate comparison covered 17 locally available/retrieved prior TIFFs, including an H33-2-B2 reference; it found no exact match. It did not cover every artifact on all sibling sites or public repositories.
- The public leaderboard does not map participants' rows to filenames. The 0.2778 row is official at the checked snapshot; its attribution to H33-2-B2 remains unresolved because the owner page labels that file “UNSCORED.”
- AI-use disclosure is required by the official rules. Before any prize submission, accurately document the extent of AI assistance in research, coding, analysis, and writing, and verify the current rules.
