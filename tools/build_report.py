"""Build a PDF from a Markdown report (default docs/FINDINGS.md).

Usage:  python tools/build_report.py [file.md]     ->  file.pdf
        (a bare name such as THEORY.md is looked up in docs/)

Math in $...$ / $$...$$ is typeset with matplotlib mathtext using the STIX
fonts (professional serif math faces, the family Office equations emulate)
and embedded as baseline-aligned images.
"""

import hashlib
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import math_to_image

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (CondPageBreak, Image, KeepTogether,
                                ListFlowable, ListItem, Paragraph,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

rcParams["mathtext.fontset"] = "stix"

ROOT = Path(__file__).resolve().parent.parent
_arg = Path(sys.argv[1] if len(sys.argv) > 1 else "FINDINGS.md")
MD_FILE = (_arg if _arg.is_file() else ROOT / "docs" / _arg).resolve()
PDF_FILE = MD_FILE.with_suffix(".pdf")
MATH_DIR = ROOT / ".cache" / "mathtext"
MATH_DIR.mkdir(parents=True, exist_ok=True)
MATH_DPI = 600.0
MATH_FALLBACKS = []

# Windows fonts with full Unicode coverage (lambda, arrows, primes...).
pdfmetrics.registerFont(TTFont("Uni", "C:/Windows/Fonts/segoeui.ttf"))
pdfmetrics.registerFont(TTFont("Uni-Bold", "C:/Windows/Fonts/segoeuib.ttf"))
pdfmetrics.registerFont(TTFont("Uni-Italic", "C:/Windows/Fonts/segoeuii.ttf"))
pdfmetrics.registerFontFamily("Uni", normal="Uni", bold="Uni-Bold",
                              italic="Uni-Italic", boldItalic="Uni-Bold")

S = {
    "h1": ParagraphStyle("h1", fontName="Uni-Bold", fontSize=17, leading=21,
                         spaceAfter=10, spaceBefore=4,
                         textColor=colors.HexColor("#1a3550")),
    "h2": ParagraphStyle("h2", fontName="Uni-Bold", fontSize=13.5, leading=17,
                         spaceBefore=14, spaceAfter=6,
                         textColor=colors.HexColor("#1a3550")),
    "h3": ParagraphStyle("h3", fontName="Uni-Bold", fontSize=11.5, leading=15,
                         spaceBefore=10, spaceAfter=4,
                         textColor=colors.HexColor("#2c567e")),
    "body": ParagraphStyle("body", fontName="Uni", fontSize=9.5, leading=13,
                           spaceAfter=5, autoLeading="max"),
    "bullet": ParagraphStyle("bullet", fontName="Uni", fontSize=9.5, leading=13,
                             autoLeading="max"),
    "cell": ParagraphStyle("cell", fontName="Uni", fontSize=8, leading=10.5,
                           autoLeading="max"),
    "cellh": ParagraphStyle("cellh", fontName="Uni-Bold", fontSize=8,
                            leading=10.5, autoLeading="max"),
    "caption": ParagraphStyle("caption", fontName="Uni-Italic", fontSize=8.5,
                              leading=11, alignment=1, spaceBefore=2,
                              spaceAfter=8, textColor=colors.HexColor("#555555")),
}

DOC_WIDTH = A4[0] - 3.0 * cm


def render_math(tex, size):
    """Typeset a math snippet; returns (png path, width pt, height pt, depth pt)."""
    key = hashlib.md5(f"{size}|{tex}".encode()).hexdigest()[:16]
    png = MATH_DIR / f"{key}.png"
    meta = MATH_DIR / f"{key}.depth"
    if not (png.exists() and meta.exists()):
        depth_px = math_to_image(f"${tex}$", str(png), dpi=MATH_DPI,
                                 prop=FontProperties(size=size), format="png")
        meta.write_text(str(depth_px))
    depth_px = float(meta.read_text())
    w_px, h_px = ImageReader(str(png)).getSize()
    k = 72.0 / MATH_DPI
    return png, w_px * k, h_px * k, depth_px * k


def math_tag(tex, size):
    """Inline <img> markup for a math snippet, baseline-aligned."""
    try:
        png, w, h, depth = render_math(tex, size)
    except Exception:
        MATH_FALLBACKS.append(tex)
        return f'<font face="Courier" size="8.5">{tex}</font>'
    return (f'<img src="{png.as_posix()}" width="{w:.2f}" height="{h:.2f}" '
            f'valign="{-depth:.2f}"/>')


def inline(text, msize=10.0):
    parts = re.split(r"(\$[^$]+\$)", text)
    out = []
    for p in parts:
        if len(p) > 2 and p.startswith("$") and p.endswith("$"):
            out.append(math_tag(p[1:-1].strip(), msize))
        else:
            p = p.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            p = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", p)
            p = re.sub(r"`(.+?)`",
                       r'<font face="Courier" size="8.5">\1</font>', p)
            out.append(p)
    return "".join(out)


def make_table(rows):
    header = [Paragraph(inline(c, 8.5), S["cellh"]) for c in rows[0]]
    body = [[Paragraph(inline(c, 8.5), S["cell"]) for c in r] for r in rows[1:]]
    ncols = len(rows[0])
    tbl = Table([header] + body, colWidths=[DOC_WIDTH / ncols] * ncols,
                repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dbe6f0")),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#9fb3c4")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]))
    return tbl


def make_image(path, caption):
    img_w, img_h = ImageReader(str(path)).getSize()
    aspect = img_h / img_w
    width = DOC_WIDTH if aspect < 0.55 else 0.72 * DOC_WIDTH
    flow = [Image(str(path), width=width, height=width * aspect, hAlign="CENTER")]
    if caption:
        flow.append(Paragraph(caption, S["caption"]))
    return flow


def build():
    lines = MD_FILE.read_text(encoding="utf-8").splitlines()
    story, para, bullets, table = [], [], [], []
    pending = []          # heading(s) waiting to be glued to the next block

    def emit(*flowables):
        """Append flowables, gluing any pending heading to them."""
        if pending:
            story.append(KeepTogether([*pending, *flowables]))
            pending.clear()
        else:
            story.extend(flowables)

    def flush():
        nonlocal para, bullets, table
        if para:
            emit(Paragraph(inline(" ".join(para)), S["body"]))
            para = []
        if bullets:
            emit(ListFlowable(
                [ListItem(Paragraph(inline(b), S["bullet"]),
                          bulletFontName="Uni") for b in bullets],
                bulletType="bullet", bulletFontSize=8, leftIndent=14))
            story.append(Spacer(1, 5))
            bullets = []
        if table:
            rows = [[c.strip() for c in row.strip().strip("|").split("|")]
                    for row in table]
            rows = [r for r in rows if not all(set(c) <= set("-: ") for c in r)]
            # Table flowables split by default; force them atomic
            emit(KeepTogether([make_table(rows)]))
            story.append(Spacer(1, 8))
            table = []

    for line in lines:
        stripped = line.strip()
        img = re.match(r"!\[(.*?)\]\((.*?)\)", stripped)
        if stripped.startswith("|"):
            table.append(stripped)
            continue
        if table:
            flush()
        if not stripped:
            flush()
        elif stripped.startswith("$$") and stripped.endswith("$$") and len(stripped) > 4:
            flush()
            try:
                png, w, h, _ = render_math(stripped[2:-2].strip(), 12.5)
                emit(KeepTogether([
                    Spacer(1, 4),
                    Image(str(png), width=w, height=h, hAlign="CENTER"),
                    Spacer(1, 6)]))
            except Exception:
                MATH_FALLBACKS.append(stripped)
                emit(Paragraph(inline(stripped.strip("$")), S["body"]))
        elif img:
            flush()
            path = MD_FILE.parent / img.group(2)
            if path.exists():
                # image and caption stay on the same page
                emit(KeepTogether(make_image(path, img.group(1))))
            else:
                print(f"WARNING: image not found: {img.group(2)}")
        elif stripped.startswith("### "):
            flush()
            # require room below, and glue the heading to its first block
            story.append(CondPageBreak(6.0 * cm))
            pending.append(Paragraph(inline(stripped[4:]), S["h3"]))
        elif stripped.startswith("## "):
            flush()
            story.append(CondPageBreak(8.0 * cm))
            pending.append(Paragraph(inline(stripped[3:]), S["h2"]))
        elif stripped.startswith("# "):
            flush()
            pending.append(Paragraph(inline(stripped[2:]), S["h1"]))
        elif stripped.startswith("- "):
            para and flush()
            bullets.append(stripped[2:])
        else:
            if bullets:                      # continuation of a bullet line
                bullets[-1] += " " + stripped
            else:
                para.append(stripped)
    flush()
    story.extend(pending)                    # heading at EOF (defensive)

    doc_title = next((ln[2:].strip() for ln in lines if ln.startswith("# ")),
                     MD_FILE.stem)
    doc = SimpleDocTemplate(str(PDF_FILE), pagesize=A4,
                            leftMargin=1.5 * cm, rightMargin=1.5 * cm,
                            topMargin=1.6 * cm, bottomMargin=1.6 * cm,
                            title=doc_title)
    doc.build(story)
    if MATH_FALLBACKS:
        print(f"WARNING: {len(MATH_FALLBACKS)} math segment(s) fell back to "
              f"plain text: {MATH_FALLBACKS[:4]}")
    print(f"saved {PDF_FILE.name} ({PDF_FILE.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    build()
