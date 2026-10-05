#!/usr/bin/env python3
"""Create a one-pixel TIFF used only for a low-level writer smoke test.

This file is intentionally NOT a GEMS submission: it has no competition grid,
CRS, bounds, footprint, or NaN-outside mask. The competition-ready candidate is
built by scripts/build_pipeline.py instead.
"""
from pathlib import Path
import struct


# Tiny little-endian TIFF: 1×1 float32 sample with value 0.731.
entries = [
    (256, 3, 1, struct.pack("<H", 1) + b"\0\0"),
    (257, 3, 1, struct.pack("<H", 1) + b"\0\0"),
    (258, 3, 1, struct.pack("<H", 32) + b"\0\0"),
    (259, 3, 1, struct.pack("<H", 1) + b"\0\0"),
    (262, 3, 1, struct.pack("<H", 1) + b"\0\0"),
    (273, 4, 1, struct.pack("<I", 134)),
    (277, 3, 1, struct.pack("<H", 1) + b"\0\0"),
    (278, 4, 1, struct.pack("<I", 1)),
    (279, 4, 1, struct.pack("<I", 4)),
    (339, 3, 1, struct.pack("<H", 3) + b"\0\0"),
]
content = bytearray(b"II" + struct.pack("<HI", 42, 8))
content += struct.pack("<H", len(entries))
for tag, typ, count, value in entries:
    content += struct.pack("<HHI4s", tag, typ, count, value)
content += struct.pack("<I", 0) + struct.pack("<f", 0.731)
output = Path("artifacts/format_smoke_test.tif")
output.parent.mkdir(parents=True, exist_ok=True)
output.write_bytes(content)
print(f"Wrote {output}; one-pixel smoke test only, NOT a submission")
