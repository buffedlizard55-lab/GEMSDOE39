# H40 hypothesis register — seven candidates not previously tried in this project's family

Generated 2026-10-05 for DrivenData competition #306 (DOE GEMS Prize Challenge).
Charter requirement satisfied: *"Before you spend a submission slot, generate 3–5
new candidate geological hypotheses (not tried before). Each must name the layer(s)
it uses, the physical signature it targets, why that signature catches a fault the
USGS/INGENIOUS catalogue is missing, and how it differs from anything already
implemented in this repo."*

Seven are given. Each is ranked by expected DTI improvement against implementation
cost. Nothing below is copied from a previous submission's pixels; the previously
scored artifacts are read **only** to calibrate the metric and to serve as a
comparison arm on local instruments.

---

## 0. What the target population actually is (verified, with links)

Before ranking hypotheses it matters what "new fault" means, because it decides
which signatures can score at all.

| Fact | Verbatim source | Link |
|---|---|---|
| Known catalogue pixels are masked out of scoring | "Pixels corresponding to known USGS/INGENIOUS faults are masked / excluded from evaluation, so they do not count towards penalty terms." — `chrisk-dd`, **DrivenData Staff**, 2026-09-16 | [thread 11516](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2) |
| The Final Round masks them too | "Re-evaluation will also mask/exclude the existing USGS/INGENIOUS faults." — same post | [thread 11516](https://community.drivendata.org/t/scoring-clarification-are-known-usgs-ingenious-faults-masked-when-scoring-and-are-they-in-the-final-round-label-set/11516/2) |
| Definition of "new fault" | "'new fault' means 'any fault pixel not already captured by USGS/INGENIOUS' and can include newly mapped geometry of an existing fault system." — `chrisk-dd`, **DrivenData Staff**, 2026-09-23 | [thread 11536](https://community.drivendata.org/t/11536/2) |
| No further hints about the test faults | "We're not sharing details about the data sources, fault types, or coverage behind the test faults beyond what's in the problem description." — `chrisk-dd`, **DrivenData Staff**, 2026-09-23 | [thread 11527](https://community.drivendata.org/t/how-were-the-new-test-faults-identified-data-sources-and-fault-types/11527/7) |
| Phase 2 value of a credible submission | "the largest prize pool (Phase 2) will use a test set that is updated by expert review of all Phase 1 submissions, so your fault predictions have an impact on final evaluation even if they are not the most performant in Phase 1." — same post | [thread 11527](https://community.drivendata.org/t/how-were-the-new-test-faults-identified-data-sources-and-fault-types/11527/7) |
| Metric definition, α=0.2, β=0.8, R=300 m | distance-weighted Tversky index, worked example `3.00/(3.00+0.2·1.89+0.8·2.00)=0.60` | [page 967](https://www.drivendata.org/competitions/306/competition-doe-gems/page/967/) |

Two consequences drive every hypothesis below:

1. **Mass on the catalogue is metrically inert.** Confirmed from the artifacts, not
   assumed: `gemsdoe-ens12-adopted-7f00890a.tif` has 6,455 pixels on the catalogue
   and `8GEMSDOE_Hedge-v2_submission.tif` has 60,988 (the whole catalogue), yet both
   have exactly 166,519 active pixels and both scored exactly **0.1563**.
2. **The target is "fault geometry the catalogue does not contain"** — including
   continuations of mapped systems. So a detector must be good at *absence*, not at
   reproducing the catalogue.

---

## 1. Why the current best (0.2778) scores what it scores — inverse-DTI calibration

This is the analytical result that the H40 family is built on, and it is derived
from organizer-returned numbers rather than from a local proxy.

Because `TP_w + FN_w = |G|` (identity, since `FN_w = Σ_g [1 − max_x p·k]`), the
metric collapses to

```
DTI = TP_w / (0.2·TP_w + 0.2·FP_w + 0.8·|G|)
```

`|G|`, the number of *active* scored-truth pixels, is not published. It is
identified from a **nested pair of scored artifacts whose difference is exactly
known**:

| artifact | active dots | score |
|---|---|---|
| `gems28-h27-4-r1-solo-d2-8-...-nan.tif` | 40,199 | 0.2708 |
| `gemsdoe32-h33-h33-2-b2-...-e5eb6e7e-zeros.tif` | 37,654 | **0.2778** |

Re-read from disk and asserted in `calibrate.verify_nested_pair` /
`tests/test_h40.py::test_nested_pair_is_exactly_as_documented`:
the second is a strict subset of the first, the removed set is exactly **2,545
dots**, and every removed dot lies at **1.414 ≤ d(catalogue) ≤ 2.000 px**
(141–200 m). Neither artifact has any pixel on the catalogue.

Solving the two-equation system for `(|G|, m)` — with `l0 = 3.0`, the exact
triangular-kernel sum along a straight trace (`1 + 2·⅔ + 2·⅓ + 0`), and
`kb = 0.5`, the mean kernel of a hit uniformly placed in the 3 px band:

```
|G| = 14,143 active scored-truth pixels   [13,987 .. 14,434] over a 5x4
                                          sensitivity grid on (l0, kb)
anchor live hit rate h = 6.130%           (2,308 of 37,654 dots within 300 m
                                           of scored truth)
forward model reproduces 0.2778 exactly, and 0.2707 vs the observed 0.2708
```

`|G|` moves by only ±3% across the whole nuisance grid — it is identified far
more tightly than the parameters it was solved with.

**Denominator decomposition at the 0.2778 submission:**

| term | value | share |
|---|---|---|
| `0.2·TP_w` | 1,095 | 5.6% |
| `0.2·FP_w` | 7,300 | 37.0% |
| `0.8·|G|` | 11,314 | **57.4%** |

**Break-even marginal hit probability** (a dot is worth adding iff
`ΔTP > 0.2·DTI·(ΔTP + ΔFP)`):

```
pi* = 3.149%   at DTI = 0.2778
```

**Corpus-wide inversion** of all 16 mirrored scored artifacts (`registry/live_scores.json`)
gives a single monotone law — implied hit rate falls as budget rises:

| n (active dots) | score | implied h |
|---|---|---|
| 37,654 | 0.2778 | 6.130% |
| 40,199 | 0.2708 | 5.742% |
| 44,090 | 0.2600 | 5.214% |
| 60,069 | 0.2477 | 4.344% |
| 61,328 | 0.2449 | 4.248% |
| 65,236 | 0.1839 | 2.853% |
| 121,131 | 0.1922 | 2.680% |
| 123,939 | 0.1855 | 2.542% |
| 166,519 | 0.1563 | 1.992% |
| 204,504 | 0.0904 | 0.953% |
| 145,610 | 0.0107 | 0.101% |

**What this says about beating 0.3262 (live #1 on 2026-10-05).** At the anchor's
budget and FP mass, 0.3262 requires `TP_w ≈ 6,489` instead of 5,470 — a **+18.6%
improvement in hit rate**, i.e. `h` from 6.13% to ~7.3%. Equivalently, at the
same `TP_w`, cutting `FP_w` from 36,379 to 21,770 also reaches 0.3262. Both
routes are the same lever: **more hits per dot.** Budget alone cannot do it — the
priced optimum for the incumbent surface is ~25–28k dots for ~0.283.

Also settled analytically and pinned by
`tests/test_h40.py::test_binary_emission_is_never_worse_than_scaling_down`:
because `0.8·|G|` is a constant additive term, `DTI(c·p)` is monotone increasing
in `c`, and the inclusion condition does not depend on `p`. **Binary {0,1}
emission at full scale is optimal; submitting soft probabilities is strictly
worse.**

---

## 2. Ranked hypotheses

Rank = expected DTI improvement ÷ implementation cost. Expected improvement is
stated as a *multiplier on hit rate* because §1 shows hit rate — not budget, not
spacing — is what moves the score.

### H40-F — supervised propensity for *missing* faults (rank 1)

* **Layers.** All 19 competition bands at σ = 0 and σ = 3 px, plus the external
  layers below (heat flow, slip/dilation tendency, deep temperature, hot-spring
  density, slip rate, young-fault density). Positive class:
  `data/external/derived_sgmc_faults_100m_u8.tif` — USGS State Geologic Map
  Compilation faults (https://mrdata.usgs.gov/geology/state/) — restricted to
  pixels **more than 300 m from the scored catalogue**: 62,122 px, of which a
  prevalence-matched subset is used.
* **Physical signature.** Not a hand-designed one — a learned one. The model is
  asked "what does a fault that the USGS/INGENIOUS catalogue does *not* contain
  look like in these 19 layers?", which is the question the metric actually asks.
* **Why it catches a missing fault.** Training on the catalogue teaches "what a
  *mapped* fault looks like", which over-weights strong geomorphic expression —
  exactly the property that got the fault mapped in the first place. Training on
  catalogue-*absent* faults inverts the bias. The negative class is drawn from
  catalogue-free background so distance-to-catalogue cannot be used as a shortcut.
* **How it differs from anything in this repo.** Every prior round in this
  project's family used SGMC-off as an **evaluation instrument only**. This uses
  it as the **training target**. All predictions are **out-of-fold** with respect
  to the same 4×6 spatial block grid used for evaluation, so scoring this channel
  on the SGMC instrument is not circular.
* **Expected improvement / cost.** Highest expected improvement (a learned
  multi-layer interaction that no hand-designed channel reproduces); highest cost
  (~4 gradient-boosted fits over 5.1 M candidate pixels). Chunked feature
  extraction is required — the full feature cube is 2.2 GB and does not fit.

### H40-A — critically stressed lineament (rank 2)

* **Layers.** External: USGS **Siler (2022)** slip/dilation tendency for
  Quaternary faults of the Great Basin, DOI [10.5066/P9YL58W6](https://doi.org/10.5066/P9YL58W6),
  [ScienceBase item 6296974dd34ec53d276bb33d](https://www.sciencebase.gov/catalog/item/6296974dd34ec53d276bb33d)
  — 37,811 in-footprint segments ≥200 m carrying `Strike`, `TS`, `TS_norm`, `TD`,
  `SHmaxAz`, `SHmax_Mag`, `ShminAZ`, `Shmin_Mag`. Competition: `det_elev`, `tmi`,
  `iso_grav_anom`, `cond_surf`, `tc`, `geod_2ndinv`, `depth_to_base_surf` for the
  lineament detector and its local strike.
* **Physical signature.** **Orientation, not amplitude.** A structure tensor over
  the fused physics stack gives a local strike; the strike is compared with the
  release's own *local* `SHmaxAz` and looked up in an **empirical favourability
  curve measured from the release itself** — mean `TS_norm` binned by
  `Δθ = strike − (SHmaxAz + 90°)`, folded to [0°, 90°]:

  ```
  Δθ =  2.5°  mean TS_norm = 0.120   (n =   407)
  Δθ = 22.5°  mean TS_norm = 0.440   (n =   975)
  Δθ = 42.5°  mean TS_norm = 0.713   (n = 1,843)
  Δθ = 62.5°  mean TS_norm = 0.822   (n = 2,999)   <-- peak
  Δθ = 87.5°  mean TS_norm = 0.756   (n = 4,199)
  ratio peak/trough = 6.9x
  ```

  A fault striking parallel to σ_Hmax is the favourably oriented normal fault;
  favourably oriented faults are the ones that are stress-loaded for slip and
  therefore hydraulically conductive
  ([Barton, Zoback & Moos 1995, *Geology* 23:913](https://doi.org/10.1130/0091-7613(1995)023<0913:UFFASA>2.3.CO;2);
  [Morris, Ferrizzoli & Zoback 1996, *Geology* 24:1107](https://doi.org/10.1130/0091-7613(1996)024<1107:FPTSOT>2.3.CO;2)).
* **Why it catches a missing fault — in the release's own words.** The USGS
  purpose statement for this exact dataset reads: *"This effort was undertaken to
  help identify faults and fault segments that are appropriately oriented to be
  stress-loaded for slip or to dilate under the ambient stress conditions. Both
  conditions may make such faults likely to host **as-yet-undiscovered
  hydrothermal processes**."* The dataset was built, by the USGS, as part of
  **INGENIOUS** — *Innovative Geothermal Exploration through Novel Investigations
  of **Undiscovered** Systems* — i.e. to find the same population this competition
  scores. It is not a proxy; it is a purpose-built predictor.
* **How it differs.** No prior artifact in this family used orientation relative
  to a measured stress field. Prior magnetic/gravity work used edge *magnitude*
  (`|∇|`, analytic signal, coherence) which is orientation-blind.
* **Expected improvement / cost.** Large improvement (a 6.9× discriminative range
  on an officially purpose-built layer), low cost (one structure tensor + one
  1-D lookup).

### H40-C — deep reservoir temperature anomaly (rank 3)

* **Layers.** External: DOE Geothermal Data Repository submission **1391**
  (INGENIOUS), DOI [10.15121/1881483](https://doi.org/10.15121/1881483),
  https://gdr.openei.org/submissions/1391 — `gdr_wellspring_in_footprint.csv`,
  27,092 rows already carrying competition-grid `row`/`col`, with `temp_c` and
  the three fluid geothermometers `geothermquartz_c`, `geothermchalc_c`,
  `geothermcat_c` (1,244 rows carry at least one).
* **Physical signature.** Maximum of the three geothermometers per location =
  **temperature at reservoir depth**, not at the surface. Residual above an 80 px
  (8 km) local mean, intersected with the lineament detector and boosted by hot
  spring/well density (`thermalclass` = Hot, or `temp_c ≥ 65 °C`; 3,163 rows).
* **Why it catches a missing fault.** A hot deep reservoir cannot exist without a
  deep permeable pathway. Surface heat flow misses this because conduction smears
  it over kilometres; a fluid geothermometer measures the water that actually
  travelled the conduit, so it localises the conduit even where the surface shows
  nothing.
* **How it differs.** Prior rounds used *surface* temperature proxies
  (radiometrics, `thermal_pop` scalars) and used well/spring data only as a
  corroboration weight. None used a reservoir-depth geothermometer as a primary
  targeting field.
* **Georeference verified, not assumed:** the mirror carries both `row`/`col` and
  `utm_x`/`utm_y`; we recompute row/col from UTM with the grid transform and
  assert agreement — **27,092 rows checked, 0 mismatches**. That independently
  confirms `origin = (243350 E, 4508550 N)`, 100 m, EPSG:32611.
* **Expected improvement / cost.** Moderate-to-large; low cost.

### H40-G — tip-continuation / linkage corridor (rank 4)

* **Layers.** The supplied catalogue (`labels.tif`) used **geometrically**, plus
  the H40-A structure strength, coherence and favourability fields.
* **Physical signature.** Skeletonise the catalogue, find its **endpoints**
  (8-neighbour count = 1 on traces ≥25 px), measure the local strike in a 6 px
  window from the trace's own second moments, and lay a corridor **60 px (6 km)
  long × 3 px half-width** beyond the tip along that strike. Weight decays as
  `(1 − L/60)^1.5`, and every corridor pixel is forced to ≥2 px from the
  catalogue so it can never re-enter the pixel-exact mask.
* **Why it catches a missing fault.** Staff confirmed that newly mapped geometry
  of an existing system counts
  ([thread 11536](https://community.drivendata.org/t/11536/2)). A trace that
  stops is either a real tip — where a relay ramp, splay or horsetail termination
  may be unmapped — or a cartographic limitation, where the fault simply
  continues. [Biasi & Wesnousky (2016), *BSSA* 106:1110](https://doi.org/10.1785/0120150223)
  show surface ruptures step and gap along strike with a characteristic
  distribution; [Faulds & Hinz (2015)](https://www.osti.gov/servlets/purl/1724082)
  report that termination/linkage/step-over settings host ~32%/~25%/~22% of
  characterised Great Basin geothermal systems.
* **How it differs.** The previous family's best artifact (0.2708 → 0.2778) used a
  *precomputed* NBMG topology-link table (`docs/data/topology_priority_h27_5b.csv`,
  81 priority links). H40-G derives the continuation geometry **from the raster
  itself** — no external link table, no NBMG FID attribution — and gates it on the
  measured stress-favourability curve rather than on name/kinematic matching. It
  also respects the empirically validated ≥2 px catalogue exclusion that the
  topology-based version did not.
* **Expected improvement / cost.** Moderate; moderate cost (a per-tip Python loop,
  bounded by tip count).

### H40-B — heat-flow anomaly (rank 5)

* **Layers.** External: USGS **DeAngelo et al. (2022)** *Heat flow maps and
  supporting data for the Great Basin, USA*, DOI
  [10.5066/P9BZPVUC](https://doi.org/10.5066/P9BZPVUC),
  [ScienceBase item 6297d2fad34ec53d276c5b28](https://www.sciencebase.gov/catalog/item/6297d2fad34ec53d276c5b28).
  2,108 unique wells in the footprint; **879** retained after QC filtering.
* **Physical signature.** Residual of measured heat flow above its own 60 px
  (6 km) local mean, intersected with lineament strength and coherence.
* **Why it catches a missing fault.** Heat flow is the thermal half of a play
  fairway; convective heat loss along a permeable fault produces a local positive
  residual that conduction-only mapping smooths away.
* **Data-handling decisions, stated:** the release ships each well twice (inputs
  layer `USGS_gbHeatFlowWells.shp`, outputs layer
  `USGS_gbHeatFlowWells_wEstimates.shp`) so rows are de-duplicated on
  `unique_id`. Source CRS is `USAEAC_83_117` (Albers Conic Equal Area, NAD83,
  centre 37.5/−117, standard parallels 29.5/45.5, no false easting/northing),
  copied verbatim from the release `.prj` as recorded in the mirror schema, and
  reprojected with pyproj. **QC class `G` is excluded**: its `hf_meas`
  distribution is a different population (median 218 mW/m², maximum 10,234 mW/m²
  against A/B/C median ~85 mW/m²) and 10,234 mW/m² is not physically possible for
  conductive crustal heat flow. Retained classes A/B/C: median 85 mW/m²,
  consistent with Great Basin background.
* **How it differs.** No prior round in this family used a measured heat-flow
  release at all; heat was inferred from radiometrics or from surface temperature.
* **Expected improvement / cost.** Moderate; low cost. Sparse coverage (879 wells
  over 51,674 km²) is the main limitation and is reported rather than hidden.

### H40-D — strain-rate kink + seismic lineament (rank 6)

* **Layers.** Competition bands **4, 7, 8, 10, 16** — `geod_2ndinv`,
  `geod_shearrate`, `geod_dilaterate`, `deq_n100a15`, `ieq_n100a15`.
* **Physical signature.** (a) The **normal curvature** of the second invariant of
  the strain-rate tensor at σ = 2, 4, 8 px: an actively deforming zone is a
  linear *kink* in an otherwise smooth field, which curvature detects and
  amplitude thresholding does not — the response is
  `|κ_n| / (|κ_n| + |κ_t|)` so elongated kinks score and blobs do not.
  (b) **Elongation** of the earthquake-density field along a consistent azimuth
  (structure-tensor coherence of `log1p(ieq)`), i.e. a seismic *lineament* rather
  than a cluster. (c) Proximity weighting by `−deq_n100a15`, and an extensional
  bonus from `geod_dilaterate`.
* **Why it catches a missing fault.** The USGS Quaternary Fault and Fold Database
  requires geomorphic evidence. Strain rate comes from GPS/InSAR and seismicity
  from catalogue locations, so both outline faults that deform or creep without
  leaving a scarp — the blind and alluvium-covered population.
* **How it differs — and a defect it repairs.** These five bands were **entirely
  unused** by the previous revision of this repository: `grid.read_all_bands`
  hardcoded a name list in which bands 3–11 and 16 were placeholders
  (`b3`…`b11`, `b16`), and `tc` was placed at index 18 where the raster's own
  metadata puts `iso_grav_anom_hg`. Names are now read from the file's GDAL
  `band_name` tags and asserted against the verified order; the test
  `test_band_order_is_read_from_the_file_not_hardcoded` pins it.
* **Expected improvement / cost.** Moderate; low cost.

### H40-E — corrected tilt-angle zero-crossing (rank 7)

* **Layers.** Competition band **6 `tc`** ("Tilt angle or total curvature"),
  corroborated with bands **3 `tmi_hg`** and **9 `tmi_vg`**.
* **Physical signature.** The tilt angle `arctan(VG/HG)` **crosses zero directly
  above a vertical contact**, so its zero crossing is a sub-pixel contact
  locator. The channel is `|∇tc| · exp(−½(tc/s)²)` at σ = 1.0, 2.0, 3.5, times a
  local-vs-regional swing of `log(1 + |tmi_vg|/|tmi_hg|)`.
* **Why it catches a missing fault.** Buried contacts under alluvium produce a
  weak-amplitude but sharp zero crossing in the tilt angle; amplitude-threshold
  mapping misses them, which is precisely how they stay out of the catalogue.
* **How it differs.** Two ways. (i) The analytic-signal *magnitude* used by
  H39-C peaks **off** contact; the zero-crossing locator peaks **on** it.
  (ii) H39-C read the wrong band entirely (see H40-D), so the previous
  tilt channel was in fact an isostatic-gravity horizontal-gradient channel.
* **Expected improvement / cost.** Small-to-moderate; very low cost.

---

## 3. Emission: the budget is a price, not a copied constant

Every hypothesis above is emitted with the same geometry, which is the geometry
of the empirically best-scoring family and is justified by an organizer-scored
differential rather than by preference:

* **Poisson-disk best-first thinning at 2.8 px (280 m) minimum spacing.** The
  scored corpus shows 2.8 px beating 1.5 px at fixed surface
  (0.2600 vs 0.2477), consistent with the redundancy argument: at 1.5 px spacing
  several dots share the same truth pixel, so extra dots buy FP without buying TP.
* **Hard exclusion of every pixel within 2 px (200 m) of the catalogue.** This is
  not a masking rule — staff say masking is pixel-exact — it is a measured
  economic fact: deleting exactly those 2,545 dots moved 0.2708 → 0.2778.
* **Stopping at the break-even price.** Emission continues down the ranked
  propensity list until the *marginal* (windowed) hit rate, transferred to the
  live scale, falls below `pi* = 3.149%`. The instrument→live transfer factor is
  measured from the live-scored anchor artifact on the same prevalence-matched
  instrument, not assumed.

`scripts/build_h40.py` prints the marginal hit-rate curve and the chosen cutoff
for every candidate, so the stopping point is auditable.

---

## 4. Validation protocol (pre-declared before any fold was evaluated)

* **Blocks.** 4 × 6 = 24 rectangular spatial blocks × 2 seeds = 48 folds per
  instrument. Blocks are eroded 6 px inside and dilated 10 px for the
  component-selection collar, so a trace is never split across the train/score
  boundary of its own fold.
* **Instruments.** I1 = hidden whole catalogue components (25% of the block's
  catalogue length, whole components only, excluding any component touching the
  collar). I2 = prevalence-matched SGMC off-catalogue faults, matched to the
  calibrated `|G| = 14,143` by dropping whole components.
* **Wald SPRT (1945).** `α = 0.05`, `β = 0.10`, `p0 = 0.5`, `p1 = 0.7`, declared
  before any fold is scored. Boundaries
  `B = ln((1−β)/α) = +2.8904`, `A = ln(β/(1−α)) = −2.2513`.
  Fold order is the deterministic sort on fold key.
* **Power check performed first.** `folds_needed()` reports that the upper
  boundary needs **9** all-wins and the lower needs **5** all-losses. The
  previous revision's 8-fold design could therefore *never* accept H1 — an 8/8
  sweep reaches LLR = 8·ln(0.7/0.5) = 2.6918 < 2.8904. With 48 I2 folds the test
  is decidable, and `tests/test_h40.py::test_holdout_design_can_reach_a_boundary`
  pins that fact.
* **Leakage gate.** No candidate is promoted if more than 2% of its emitted
  pixels lie on the evaluation truth.
* **Promotion gate.** A candidate must (a) reach `accept_H1` on the SPRT against
  the live-scored anchor arm, (b) have a live-priced hit rate ≥ 1.02× the
  anchor's 6.130%, (c) pass the leakage gate, and (d) price above 0.2778. If
  nothing clears all four, the pre-declared fallback emits the best-priced
  candidate at the anchor's validated budget rather than gambling a slot.

---

## 5. Irregularities found and flagged this round

1. **Band-name defect (fixed).** `grid.read_all_bands` hardcoded
   `["mag_anom","rtp","b3",…,"b11","det_elev","iso_grav_anom","tmi",
   "depth_to_base_surf","b16","cond_surf","tc","det_elev_slope"]`. The raster's
   own `band_name` tags give `tc` at band 6 and `iso_grav_anom_hg` at band 18, so
   every tilt-based detector silently ran on the gravity horizontal gradient, and
   11 of 19 bands (all three geodetic strain bands, both seismic bands,
   `tmi_hg`, `tmi_vg`, and all four gravity-gradient bands) were unreachable.
2. **SPRT audit-record defect (fixed).** `sprt_pairwise` recomputed `wins` and
   `total` by re-iterating the input *after* the accumulation loop; a generator
   input was already exhausted, so both were reported as 0 while `n` was correct.
   It also counted wins from folds consumed after the stopping boundary.
3. **SPRT had no power (fixed).** 8 folds cannot cross an upper boundary that
   needs 9 wins.
4. **Stale timestamp in `README.md` (fixed).** The README named
   `gemsdoe39-h39-b-20261005T012128Z-zeros.tif`; the file that exists and that
   `docs/index.html` links is `…T020000Z-zeros.tif`.
5. **Malformed SHA-256 in `docs/index.html` (fixed).** The NaN variant was printed
   as `c985646d4a1bb6de975b65fdc07904faaa55d4739eb176942c537d63f9` — 56 hex
   characters. The true digest is
   `c985646d4a1bb6de975290b75b65fdc07904faaa55d4739eb176942c537d63f9`.
6. **Manifest not reproducible from committed code.**
   `docs/downloads/gemsdoe39-h39-b-20261005T020000Z-manifest.json` reports
   `catalogue_hidden_folds: {}` for every candidate while claiming
   `sprt_decisions` with `n=8, wins=8`. `build_pipeline.py` derives `wins` by
   iterating that same empty dict, so it would have produced `n=0`. The committed
   run could not have produced that manifest.
7. **Holdout proxy off by ~60× (addressed).** The previous manifest's holdout DTI
   values (0.0022–0.0051) sit two orders of magnitude below the organizer's
   0.24–0.28 for the same artifacts. Ranking on that proxy is close to ranking on
   noise. §1 replaces it with a calibration anchored on organizer-returned scores.
8. **Emission under-filled its own budget (addressed).** `docs/index.html` reports
   35,805 emitted pixels against a `--budget` of 44,090: the emitter ran out of
   eligible pixels and said nothing. The H40 emitter reports pool size, chosen
   cutoff and the marginal hit rate at the cutoff.
9. **`labels.tif` and `existing_faults.tif` are byte-identical**
   (both 425,830 B, both SHA-256 `7ba308ccdc44…`). Flagged, not silently
   deduplicated.
10. **`sample_submission.tif` is not empty.** It contains 60,988 positive pixels —
    exactly the catalogue count — so the template is the supplied catalogue
    rasterised to float32, not a blank canvas.
11. **Network reachability.** This sandbox reaches only
    `github.com` / `api.github.com` / `pypi.org`. Verified 2026-10-05: HTTP `000`
    for `sciencebase.gov`, `gdr.openei.org`, `mrdata.usgs.gov`,
    `web2.nbmg.unr.edu`, `osti.gov`, `pangea.stanford.edu` and
    `community.drivendata.org`. The external payloads therefore arrive through the
    owner's SHA-256-pinned GitHub mirrors. **They are integrity-pinned but not
    organizer-authenticated**, and every document here says so rather than
    implying otherwise. The *text* of the forum threads, the ScienceBase purpose
    statement and the live leaderboard were read through the page-fetch service
    and are quoted verbatim above.

---

## 6. What is explicitly **not** claimed

* No claim that any H40 channel scores above 0.2778 on the organizer's private
  set. The priced numbers are model output from a 2-equation calibration with
  stated nuisance priors, and the transfer factor is estimated from a **single**
  anchor artifact.
* No claim that `|G| = 14,143` exactly. The defensible statement is
  `|G| ∈ [13,987, 14,434]` under the stated priors, and the calibration's value
  is that it is *tight*, not that it is exact.
* No claim that the IDW interpolations of heat flow, slip tendency or deep
  temperature are kriging estimates. No variogram model is fitted; they are
  two-scale Gaussian-weighted local means and are described as such.
* No claim that the SGMC off-catalogue faults *are* the scored population. They
  are an official, independent realisation of "fault pixels absent from the
  scored catalogue", used because nothing better is available offline.
* No claim that the external mirrors are the bytes USGS published. They are the
  bytes the owner pinned, and the pin is what is verified.
