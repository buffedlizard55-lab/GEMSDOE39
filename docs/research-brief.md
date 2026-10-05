# Research brief — GEMSDOE39

**This brief is superseded.** The round-1 (H39) register is retained below for
provenance only; it is not the basis of the current artifact.

The current round is **H40 — play-fairway permeability targeting**. Its register,
calibration, validation protocol, ranked hypotheses and irregularity list live in:

* **[`docs/research/h40-hypotheses.md`](research/h40-hypotheses.md)** — the seven
  ranked H40 hypotheses, the inverse-DTI calibration that identifies `|G|`, the
  break-even emission price, the pre-declared SPRT design, and every flagged
  irregularity with its fix.
* [`docs/index.html`](index.html) — the landing page and download button,
  generated from the run manifest by `scripts/build_site.py`.
* [`docs/executive-summary.html`](executive-summary.html) — the six-step upload
  procedure, the zeros-vs-NaN decision, and the source register.
* [`README.md`](../README.md) — the permanent project charter, status, quick
  start and limitations.

Why the round-1 register is not used: its holdout proxy returned DTI values of
0.0022–0.0051 for artifacts the organizer scored at 0.24–0.28, a ~60× scale
error, and its `catalogue_hidden_folds` records were empty while its SPRT block
claimed `n=8, wins=8`. Ranking on that proxy is close to ranking on noise. See
irregularities 6 and 7 in the H40 register.

---

## Round-1 (H39) register — retained for provenance

| Rank | ID | Hypothesis | Layers |
|---|---|---|---|
| 1 | H39-B | Scarp-curvature step-over / horsetail splay | `det_elev`, `det_elev_slope`, 1 m LiDAR scarp composite |
| 2 | H39-A | Cross-gradient tensor multi-physics edge eigen-coherence | `rtp`, `tmi`, `iso_grav_anom`, `cond_surf` |
| 3 | H39-C | Tilt-derivative analytic-signal edge | `tc`, `tmi`, `rtp` |
| 4 | H39-D | Magnetic upward-continued "worm" ridges | `rtp`, `tmi` |
| 5 | H39-E | Basement-depth × conductivity co-located parallel edge | `depth_to_base_surf`, `cond_surf` |

Two of these were invalidated by defects found in the H40 round:

* **H39-C read the wrong band.** `grid.read_all_bands` placed `tc` at index 18,
  where the raster's own `band_name` metadata puts `iso_grav_anom_hg`. `tc` is
  band 6. H39-C therefore measured the isostatic-gravity horizontal gradient.
* **H39-A/B/D/E could not reach 11 of the 19 bands**, because bands 3–11 and 16
  were placeholders (`b3`…`b11`, `b16`) in the same hardcoded list. The
  geodetic strain-rate bands, both seismic bands, `tmi_hg`, `tmi_vg` and three
  gravity-gradient bands were unreachable.

Both defects are fixed and pinned by
`tests/test_h40.py::test_band_order_is_read_from_the_file_not_hardcoded`.
