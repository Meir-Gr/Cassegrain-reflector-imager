"""Verify docs/FINDINGS.pdf pagination: every table's first and last row must
land on the same page (checked via text pairs unique to each table)."""

from pathlib import Path

from pypdf import PdfReader

PDF_FILE = Path(__file__).resolve().parent.parent / "docs" / "FINDINGS.pdf"

# (label, text in first row / header, text in last row)
CHECKS = [
    ("baseline system table", "Primary diameter", "550 nm"),
    ("scales table", "0.567 arcsec", "2.75"),
    ("EE80 table", "EE80 diameter", "5104"),
    ("survival hierarchy table", "Finest edge fringes", "Dies last"),
    ("ring contrast table", "Measured contrast ratio", "broadband residual"),
    ("arago spot table", "core invisible", "merged into hole"),
    ("vibration table", "measured", "0.74"),
    ("size laws table", "Donut outer diameter", "see 5.3"),
    ("vane width table", "Shadow depth", "2.6e-04"),
    ("formula sheet", "Defocus phase", "rotation of"),
    ("files table", "cassegrain_pupil_and_psfs.png",
     "cassegrain_arago_vs_source_metrics.png"),
]

pages = [(p.extract_text() or "") for p in PdfReader(PDF_FILE).pages]


def pages_with(snippet):
    return [i + 1 for i, t in enumerate(pages) if snippet in t]


bad = 0
for label, first, last in CHECKS:
    p1, p2 = pages_with(first), pages_with(last)
    if not p1 or not p2:
        print(f"  ? {label:26s} snippet not found ({first!r}: {p1}, {last!r}: {p2})")
        bad += 1
    elif set(p1) & set(p2):
        print(f"  OK {label:26s} page {sorted(set(p1) & set(p2))[0]}")
    else:
        print(f"  SPLIT {label:23s} rows on pages {p1} vs {p2}")
        bad += 1

print(f"\n{len(CHECKS) - bad}/{len(CHECKS)} tables intact; total pages: {len(pages)}")
