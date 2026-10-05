from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.evaluate_h39y01 import _build_report, _fold_masks, _make_cells
from gems39.sprt_select import PairwiseSPRT


def test_locked_grid_order_erosion_and_seed_sequence():
    footprint = np.ones((800, 800), bool)
    catalogue = np.zeros_like(footprint)
    catalogue[50:61, 50:61] = True
    catalogue[10, 50:751] = True
    catalogue[250:261, 250:261] = True

    cells, grid_record, components = _make_cells(footprint, catalogue)

    assert [(cell.row, cell.col) for cell in cells] == [
        (0, 0), (0, 2), (0, 4), (0, 6),
        (2, 0), (2, 2), (2, 4), (2, 6),
        (4, 0), (4, 2), (4, 4), (4, 6),
        (6, 0), (6, 2), (6, 4), (6, 6),
    ]
    assert [cell.seed for cell in cells] == list(range(20261005, 20261021))
    assert grid_record["component_connectivity"] == 8
    assert cells[0].active.sum() == 70 * 70
    assert cells[0].n_truth == 121
    assert cells[5].n_truth == 121

    train_mask, train_catalogue, hidden_map = _fold_masks(
        cells[0], footprint, catalogue, components
    )
    # Component 1 crosses far beyond the holdout tile; the complete component
    # is hidden from model fitting, not just the pixels inside the collar.
    assert hidden_map[10, 700]
    assert not train_catalogue[10, 700]
    assert not train_mask[10, 700]
    assert not train_mask[50, 50]
    assert train_mask[200, 200]


def test_report_builder_preserves_grid_fields_without_duplicate_keyword(tmp_path):
    footprint = np.ones((800, 800), bool)
    catalogue = np.zeros_like(footprint)
    catalogue[50:61, 50:61] = True
    cells, grid_record, components = _make_cells(footprint, catalogue)
    paths = {}
    for name in ("sample", "labels", "features", "ratios"):
        path = tmp_path / name
        path.write_text(name)
        paths[name] = path
    input_info = dict(
        reference=dict(height=800, width=800, crs="EPSG:32611", transform=[100, 0, 0, 0, -100, 0]),
        footprint_pixels=int(footprint.sum()),
        catalogue_pixels=int(catalogue.sum()),
        raw_labels_outside_footprint=0,
        paths=paths,
    )
    sprt = PairwiseSPRT()
    report = _build_report(
        "inconclusive",
        input_info,
        grid_record,
        cells,
        int(components.max()),
        {},
        [],
        sprt,
        "2026-10-05T00:00:00+00:00",
        None,
    )

    assert report["geometry"]["full_catalogue_component_count"] == int(components.max())
    assert report["geometry"]["grid_n"] == 8
    assert report["promotion_gate"]["weekly_slot_recommendation"] is False
