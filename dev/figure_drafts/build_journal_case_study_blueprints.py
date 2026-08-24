#!/usr/bin/env python3
"""Render the canonical posteriorHCA case-study figure blueprints.

Only the current SAVI and four-gene IF/TC HCC designs are retained. Output is
PDF-only; temporary PNG previews are written to /tmp for visual quality control.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

DEPS = Path("/tmp/posteriorhca_pptx_deps")
if str(DEPS) not in sys.path:
    sys.path.insert(0, str(DEPS))

from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas as pdfcanvas

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "dev" / "figure_drafts"
MM = 72.0 / 25.4
W = 183 * MM
H = 170 * MM
SCALE = 2.5

FONT_REG = "/usr/share/fonts/dejavu/DejaVuSans.ttf"
FONT_BOLD = "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf"
FONT_ITALIC = "/usr/share/fonts/dejavu/DejaVuSans-Oblique.ttf"
pdfmetrics.registerFont(TTFont("DV", FONT_REG))
pdfmetrics.registerFont(TTFont("DV-Bold", FONT_BOLD))
pdfmetrics.registerFont(TTFont("DV-Italic", FONT_ITALIC))

COL = {
    "black": "#111111", "dark": "#3F3F3F", "mid": "#777777",
    "light": "#D9D9D9", "pale": "#F2F2F2", "white": "#FFFFFF",
    "healthy": "#8C8C8C", "ctrl": "#009E73", "treated": "#7A5195",
    "savi": "#D55E00", "nonrec": "#0072B2", "rec": "#D55E00",
    "gold": "#E69F00", "blue_pale": "#EAF3F8",
    "orange_pale": "#FCEFE8", "green_pale": "#E8F5F0",
    "purple_pale": "#F0EDF7",
}

SAVI_PDF = OUT / "posteriorHCA_Fig2_SAVI_Nature_blueprint.pdf"
HCC_PDF = OUT / "posteriorHCA_Fig3_HCC_Nature_blueprint.pdf"
SAVI_PREVIEW = Path("/tmp/posteriorHCA_Fig2_SAVI_Nature_blueprint.png")
HCC_PREVIEW = Path("/tmp/posteriorHCA_Fig3_HCC_Nature_blueprint.png")
MIN_PRINT_TEXT_PT = 5.0


def _rgb(hex_color: str) -> tuple[int, int, int]:
    value = hex_color.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


class Surface:
    """Top-left coordinate drawing surface backed by PDF and a QC bitmap."""

    def __init__(self, pdf_path: Path, preview_path: Path, height: float = H):
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        self.pdf_path = pdf_path
        self.preview_path = preview_path
        self.height = height
        self.c = pdfcanvas.Canvas(
            str(pdf_path),
            pagesize=(W, self.height),
            pageCompression=1,
            initialFontName="DV",
            initialFontSize=6,
            initialLeading=7.2,
        )
        self.c.setTitle(pdf_path.stem)
        self.img = Image.new("RGB", (round(W * SCALE), round(self.height * SCALE)), "white")
        self.d = ImageDraw.Draw(self.img)

    def _pil_font(self, size: float, bold: bool = False, italic: bool = False):
        path = FONT_BOLD if bold else FONT_ITALIC if italic else FONT_REG
        return ImageFont.truetype(path, max(1, round(size * SCALE)))

    def text(
        self,
        x: float,
        y: float,
        value: str,
        size: float = 6.0,
        color: str = COL["black"],
        bold: bool = False,
        italic: bool = False,
        align: str = "left",
    ) -> None:
        font_name = "DV-Bold" if bold else "DV-Italic" if italic else "DV"
        width = pdfmetrics.stringWidth(value, font_name, size)
        xx = x - width / 2 if align == "center" else x - width if align == "right" else x
        self.c.setFont(font_name, size)
        self.c.setFillColor(HexColor(color))
        self.c.drawString(xx, self.height - y - size * .80, value)

        font = self._pil_font(size, bold, italic)
        box = self.d.textbbox((0, 0), value, font=font)
        xp = x * SCALE
        if align == "center":
            xp -= (box[2] - box[0]) / 2
        elif align == "right":
            xp -= box[2] - box[0]
        self.d.text((xp, y * SCALE), value, font=font, fill=_rgb(color))

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        color: str = COL["black"],
        width: float = .6,
        dash: tuple[float, float] | None = None,
    ) -> None:
        self.c.setStrokeColor(HexColor(color))
        self.c.setLineWidth(width)
        self.c.setDash(*(dash or ()))
        self.c.line(x1, self.height - y1, x2, self.height - y2)
        self.c.setDash()
        if dash:
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy)
            if length == 0:
                return
            ux, uy = dx / length, dy / length
            position = 0.0
            while position < length:
                end = min(length, position + dash[0])
                self.d.line(
                    ((x1 + ux * position) * SCALE, (y1 + uy * position) * SCALE,
                     (x1 + ux * end) * SCALE, (y1 + uy * end) * SCALE),
                    fill=_rgb(color), width=max(1, round(width * SCALE)),
                )
                position += dash[0] + dash[1]
        else:
            self.d.line(
                (x1 * SCALE, y1 * SCALE, x2 * SCALE, y2 * SCALE),
                fill=_rgb(color), width=max(1, round(width * SCALE)),
            )

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        fill: str | None = None,
        stroke: str | None = None,
        width: float = .5,
    ) -> None:
        self.c.setLineWidth(width)
        self.c.setStrokeColor(HexColor(stroke or fill or COL["white"]))
        self.c.setFillColor(HexColor(fill or COL["white"]))
        self.c.rect(x, self.height - y - h, w, h, fill=1 if fill else 0, stroke=1 if stroke else 0)
        self.d.rectangle(
            (x * SCALE, y * SCALE, (x + w) * SCALE, (y + h) * SCALE),
            fill=_rgb(fill) if fill else None,
            outline=_rgb(stroke) if stroke else None,
            width=max(1, round(width * SCALE)),
        )

    def circle(
        self,
        x: float,
        y: float,
        r: float = 2.2,
        fill: str = COL["white"],
        stroke: str | None = COL["black"],
        width: float = .6,
    ) -> None:
        self.c.setLineWidth(width)
        self.c.setFillColor(HexColor(fill))
        self.c.setStrokeColor(HexColor(stroke or fill))
        self.c.circle(x, self.height - y, r, fill=1, stroke=1 if stroke else 0)
        self.d.ellipse(
            ((x - r) * SCALE, (y - r) * SCALE, (x + r) * SCALE, (y + r) * SCALE),
            fill=_rgb(fill), outline=_rgb(stroke) if stroke else None,
            width=max(1, round(width * SCALE)),
        )

    def square(self, x: float, y: float, r: float, fill: str, stroke: str | None = None):
        self.rect(x - r, y - r, 2 * r, 2 * r, fill=fill, stroke=stroke)

    def diamond(self, x: float, y: float, r: float, fill: str, stroke: str | None = None):
        points = [(x, y - r), (x + r, y), (x, y + r), (x - r, y)]
        path = self.c.beginPath()
        path.moveTo(points[0][0], self.height - points[0][1])
        for px, py in points[1:]:
            path.lineTo(px, self.height - py)
        path.close()
        self.c.setFillColor(HexColor(fill))
        self.c.setStrokeColor(HexColor(stroke or fill))
        self.c.drawPath(path, fill=1, stroke=1 if stroke else 0)
        self.d.polygon(
            [(px * SCALE, py * SCALE) for px, py in points],
            fill=_rgb(fill), outline=_rgb(stroke) if stroke else None,
        )

    def finish(self) -> None:
        self.c.showPage()
        self.c.save()
        self.img.save(self.preview_path, dpi=(300, 300))


def panel_label(s: Surface, label: str, y: float) -> None:
    s.text(5, y + 1, label, size=7, bold=True)


def axis_x(
    s: Surface,
    x0: float,
    x1: float,
    y: float,
    lo: float,
    hi: float,
    ticks: list[float],
    label: str,
    fmt=lambda z: f"{z:g}",
) -> None:
    s.line(x0, y, x1, y, width=.55)
    for tick in ticks:
        xx = x0 + (tick - lo) / (hi - lo) * (x1 - x0)
        s.line(xx, y, xx, y + 3, width=.45)
        s.text(xx, y + 5, fmt(tick), size=5.2, align="center")
    s.text((x0 + x1) / 2, y + 15, label, size=6.0, align="center")


def legend_mark(s: Surface, x: float, y: float, kind: str, color: str, label: str) -> float:
    if kind == "circle":
        s.circle(x, y, 2.1, color, None)
    elif kind == "square":
        s.square(x, y, 2.1, color)
    else:
        s.diamond(x, y, 2.4, color)
    s.text(x + 5, y - 3, label, size=5.5)
    return x + 8 + pdfmetrics.stringWidth(label, "DV", 5.5)


class PublicationSurface(Surface):
    """Surface with the 5-pt production floor used by the journal blueprints."""

    def text(
        self,
        x: float,
        y: float,
        value: str,
        size: float = 6.0,
        color: str = COL["black"],
        bold: bool = False,
        italic: bool = False,
        align: str = "left",
    ) -> None:
        super().text(
            x,
            y,
            value,
            size=max(MIN_PRINT_TEXT_PT, size),
            color=color,
            bold=bold,
            italic=italic,
            align=align,
        )


def read_reactable_table(path: Path, required_columns: set[str]) -> list[dict[str, str]]:
    """Read one reactable JSON payload embedded as a single HTML line."""
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped.startswith('{"columns":'):
            continue
        try:
            payload = json.loads(stripped)
        except json.JSONDecodeError:
            continue
        columns = {
            str(col["name"][0]): str(col["label"][0])
            for col in payload.get("columns", [])
        }
        labels = set(columns.values())
        if not required_columns.issubset(labels):
            continue
        return [
            {columns[key]: str(value) for key, value in row.items() if key in columns}
            for row in payload.get("data", [])
        ]
    raise RuntimeError(f"No reactable table with {sorted(required_columns)} in {path}")


def dashed_frame(s: Surface, x: float, y: float, w: float, h: float) -> None:
    dash = (3, 2)
    for x1, y1, x2, y2 in (
        (x, y, x + w, y),
        (x + w, y, x + w, y + h),
        (x + w, y + h, x, y + h),
        (x, y + h, x, y),
    ):
        s.line(x1, y1, x2, y2, color=COL["mid"], width=0.55, dash=dash)


def planned_tag(s: Surface, x: float, y: float, text: str = "ANALYSIS NEEDED") -> None:
    width = pdfmetrics.stringWidth(text, "DV-Bold", 5.0) + 9
    s.rect(x - width, y, width, 11, fill=COL["orange_pale"], stroke=COL["gold"], width=0.45)
    s.text(x - width / 2, y + 2.1, text, size=5.0, bold=True, align="center")


def compact_axis_x(
    s: Surface,
    x0: float,
    x1: float,
    y: float,
    lo: float,
    hi: float,
    ticks: tuple[float, ...],
    label: str,
) -> None:
    """Five-point x axis for shallow panels without crossing the panel frame."""
    s.line(x0, y, x1, y, width=.55)
    for tick in ticks:
        xx = x0 + (tick - lo) / (hi - lo) * (x1 - x0)
        s.line(xx, y, xx, y + 2.5, width=.4)
        s.text(xx, y + 4, f"{tick:g}", size=5.0, align="center")
    s.text((x0 + x1) / 2, y + 11.5, label, size=5.0, align="center")


def draw_savi_design_support_compact(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    """Compact, source-faithful collection schedule and expression-model support."""
    female_fill = True
    female_hist = "#666666"
    male_hist = "#C2C2C2"
    state_color = {"CTRL": COL["ctrl"], "U": COL["savi"], "T": COL["treated"]}

    # Left: exact-age collections plus a separate non-numeric Adult category.
    lx0, lx1 = x + 37, x + 217
    rows_y = {"STING1 · v2": y + 25, "STING2 · v3": y + 38, "STING3 · v3": y + 51}
    to_age = lambda age: lx0 + age / 40 * (lx1 - lx0)
    s.text(x + 5, y + 3, "Observed collections", size=5.2, bold=True)
    cursor = x + 85
    for state, label, mark in (("CTRL", "CTRL", "circle"), ("U", "untreated", "diamond"), ("T", "JAK inhibitor", "square")):
        color = state_color[state]
        if mark == "circle":
            s.circle(cursor, y + 7, 1.6, color, None)
        elif mark == "diamond":
            s.diamond(cursor, y + 7, 1.8, color, None)
        else:
            s.square(cursor, y + 7, 1.6, color)
        s.text(cursor + 4, y + 4, label, size=5.0)
        cursor += 10 + pdfmetrics.stringWidth(label, "DV", 5.0)
    s.circle(x + 251, y + 7, 1.5, COL["dark"], None)
    s.text(x + 255, y + 4, "F", size=5.0)
    s.circle(x + 266, y + 7, 1.5, COL["white"], COL["dark"], .45)
    s.text(x + 270, y + 4, "M", size=5.0)
    s.text(x + 5, y + 12, "P1 family: P2 uncle · C8 paternal aunt · C9 mother", size=5.0)
    s.text(
        x + 5,
        y + 18,
        "Tx: P1 ruxolitinib 52 mo / tofacitinib 2 mo · P2 ruxolitinib 1 mo · P4/P6 ruxolitinib 13 mo",
        size=5.0,
        color=COL["mid"],
    )
    for label, yy in rows_y.items():
        s.text(x + 4, yy - 3, label, size=5.0)
        s.line(lx0, yy, lx1, yy, color=COL["light"], width=.3)

    def collection_mark(age: float, row: str, lane: float, state: str, sex: str) -> tuple[float, float]:
        xx, yy = to_age(age), rows_y[row] + lane
        color = state_color[state]
        fill = color if sex == "F" else COL["white"]
        if state == "CTRL":
            s.circle(xx, yy, 1.55, fill, color, .5)
        elif state == "U":
            s.diamond(xx, yy, 1.8, fill, color)
        else:
            s.square(xx, yy, 1.55, fill, color)
        return xx, yy

    # Controls; all exact labels come from the source/GEO metadata.
    for _, age, row, lane, sex in (
        ("C7", 5, "STING1 · v2", -3.2, "F"),
        ("C8", 36, "STING1 · v2", -3.2, "F"),
        ("C9", 35, "STING1 · v2", 3.2, "F"),
        ("C10", 10, "STING2 · v3", -3.2, "M"),
        ("C11", 16, "STING2 · v3", -3.2, "M"),
    ):
        collection_mark(age, row, lane, "CTRL", sex)

    # Patient trajectories. Labels at treated points give drug and duration.
    trajectories = (
        ("P1", "STING1 · v2", "F", ((4, 3.2, "U", "P1"), (8, 3.2, "T", "R52"), (8.5, -3.2, "T", "T2"))),
        ("P2", "STING1 · v2", "M", ((31, 3.2, "U", "P2"), (34, 3.2, "T", "R1"))),
        ("P4", "STING2 · v3", "M", ((7.5, 3.2, "U", "P4"), (8.5, 3.2, "T", "R13"))),
        ("P6", "STING2 · v3", "M", ((8, -3.2, "U", "P6"), (9.5, -3.2, "T", "R13"))),
    )
    for _, row, sex, visits in trajectories:
        coords: list[tuple[float, float]] = []
        for age, lane, state, _ in visits:
            coords.append(collection_mark(age, row, lane, state, sex))
        for start, end in zip(coords[:-1], coords[1:]):
            s.line(start[0], start[1], end[0], end[1], color=COL["mid"], width=.4)
    p5x, p5y = collection_mark(10.5, "STING3 · v3", 0, "U", "M")
    s.text(p5x, p5y + 4, "P5", size=5.0, color=COL["dark"], align="center")

    # C1/C2 are categorical Adult observations; no finite age interval is invented.
    adult_x = x + 229
    s.text(adult_x + 33, rows_y["STING1 · v2"] - 3, "C7/8/9 · P1/2", size=5.0, color=COL["mid"], align="center")
    s.text(adult_x + 33, rows_y["STING2 · v3"] - 3, "C10/11 · P4/6", size=5.0, color=COL["mid"], align="center")
    s.rect(adult_x, rows_y["STING3 · v3"] - 6, 66, 12, fill=COL["white"], stroke=COL["ctrl"], width=.45)
    s.text(adult_x + 33, rows_y["STING3 · v3"] - 2.8, "C1/C2 · Adult", size=5.0, color=COL["dark"], bold=True, align="center")
    s.text(adult_x + 33, rows_y["STING3 · v3"] + 3, "exact age unavailable", size=5.0, color=COL["mid"], align="center")
    compact_axis_x(s, lx0, lx1, y + h - 13, 0, 40, (0, 10, 20, 30, 40), "exact age at collection (years)")
    # Right: exact current expression-reference support; greys avoid cohort-colour reuse.
    rx0, rx1 = x + 321, x + w - 5
    base, hist_top = y + h - 14, y + 27
    female_bins = (9, 4, 47, 92, 73, 50, 31, 8, 3)
    male_bins = (11, 4, 41, 62, 56, 39, 24, 12, 3)
    max_total = max(a + b for a, b in zip(female_bins, male_bins))
    s.text(rx0, y + 3, "HCA expression support", size=5.0, bold=True)
    s.text(rx0, y + 10, "569 donors · 973 samples · 13 datasets", size=5.0)
    s.text(rx0, y + 16, "46 samples in coarse harmonised 10x 3′ group", size=5.0, color=COL["mid"])
    s.text(rx0, y + 22, "query status: observed · integrated over · outside support", size=4.8, color=COL["mid"])
    s.text(rx1, y + 10, "unique healthy-HCA donors", size=4.8, color=COL["mid"], align="right")
    bw = (rx1 - rx0) / len(female_bins)
    for i, (nf, nm) in enumerate(zip(female_bins, male_bins)):
        hf = (base - hist_top) * nf / max_total
        hm = (base - hist_top) * nm / max_total
        bx = rx0 + i * bw
        s.rect(bx + .4, base - hm, bw - .8, hm, fill=male_hist)
        s.rect(bx + .4, base - hm - hf, bw - .8, hf, fill=female_hist)
    for age in (4, 5, 7.5, 8, 8.5, 9.5, 10, 10.5, 16, 31, 34, 35, 36):
        xx = rx0 + age / 90 * (rx1 - rx0)
        s.line(xx, hist_top - 2, xx, hist_top + 3, color=COL["savi"], width=.55)
    s.line(rx0, base, rx1, base, color=COL["mid"], width=.35)
    for age in (0, 30, 60, 90):
        xx = rx0 + age / 90 * (rx1 - rx0)
        s.line(xx, base, xx, base + 2, color=COL["mid"], width=.3)
        s.text(xx, base + 3, str(age), size=5.0, align="center")
    s.text((rx0 + rx1) / 2, y + h - 4, "unique donor age · orange ticks, exact SAVI ages", size=5.0, align="center")
    s.text(x + w - 4, y + 3, "overlap audit pending", size=4.3, color=COL["dark"], bold=True, align="right")


def draw_savi_composition_compact(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
    rows: list[dict[str, str]],
) -> None:
    """Blank two-contrast frame for the final patient-level composition model."""
    dashed_frame(s, x, y, w, h)
    planned_tag(s, x + w - 3, y + 3, "PATIENT MODEL NEEDED")
    names = (
        ("Naive CD4 T", "↑"),
        ("Eff./memory CD4", "↓"),
        ("Naive CD8 T", "↑"),
        ("Eff./memory CD8", "↓"),
        ("NK", "↓"),
        ("MAIT", "↓"),
        ("gamma-delta T", "↓"),
    )
    s.text(x + 7, y + 18, "source-reported direction", size=5.0, bold=True)
    legend_mark(s, x + 11, y + 29, "circle", COL["mid"], "untreated − study control")
    legend_mark(s, x + 106, y + 29, "diamond", COL["healthy"], "untreated − HCA")
    x0, x1 = x + 85, x + w - 10
    top, bottom = y + 41, y + h - 40
    lo, hi = -15.0, 15.0
    zero = x0 + (0 - lo) / (hi - lo) * (x1 - x0)
    s.line(zero, top, zero, bottom, color=COL["mid"], width=.4, dash=(2, 2))
    dy = (bottom - top) / len(names)
    for i, (display, direction) in enumerate(names):
        yy = top + (i + .5) * dy
        s.text(x + 7, yy - 2.4, display, size=5.0)
        s.text(x + 75, yy - 2.4, direction, size=5.0, bold=True, align="center")
        s.line(x0, yy, x1, yy, color=COL["light"], width=.25)
    compact_axis_x(s, x0, x1, y + h - 35, lo, hi, (-15, 0, 15), "composition difference (percentage points)")
    s.text(x + 7, y + h - 15, "blank: eight-category patient model + HCA composition-support audit", size=4.4, color=COL["mid"])
    s.text(x + 7, y + h - 7, "HCA = unseen-study mean · broad NK (source strongest: CD56-bright)", size=4.0)


def draw_savi_signature_compact(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
    rows: list[dict[str, str]],
) -> None:
    """Blank aligned slot for covered-subset and constituent posterior results."""
    dashed_frame(s, x, y, w, h)
    planned_tag(s, x + w - 3, y + 3, "ALIGNED SCORE + POSTERIORS NEEDED")
    genes = (
        "ADRB2", "PTGS2", "AQP9", "OTUD1", "PTX3", "TGIF1", "TSC22D2", "ZBTB43",
        "ZFP36L2", "JARID2", "KLF4", "TREM1", "ANKRD28", "ANXA1", "PLAUR", "IL1B",
    )
    s.text(x + 7, y + 18, "Published 21-gene IFN-independent STING-activation signature", size=5.1, bold=True)
    s.text(
        x + 7,
        y + 25,
        "16-gene covered subset: median + 95% subject-coupled stability envelope",
        size=4.2,
        color=COL["dark"],
    )
    plot_left = x + 58
    gap = 6
    col_w = (w - 72 - 3 * gap) / 4
    columns = (
        ("untreated −", "study control"),
        ("untreated −", "unseen-study HCA"),
        ("paired untreated −", "treated · 4 patients"),
        ("study control −", "unseen-study HCA"),
    )
    lo, hi = -4.0, 7.0
    top, bottom = y + 55, y + h - 25
    for k, (line1, line2) in enumerate(columns):
        cx0 = plot_left + k * (col_w + gap)
        cx1 = cx0 + col_w
        zero = cx0 + (0 - lo) / (hi - lo) * (cx1 - cx0)
        s.rect(cx0, y + 31, col_w, 7, fill=COL["pale"], stroke=COL["light"], width=.3)
        s.line(zero, y + 32, zero, y + 37, color=COL["mid"], width=.3, dash=(1.5, 1.5))
        s.text((cx0 + cx1) / 2, y + 41, line1, size=3.8, bold=True, align="center")
        s.text((cx0 + cx1) / 2, y + 47, line2, size=3.8, align="center")
        s.line(zero, top, zero, bottom, color=COL["mid"], width=.35, dash=(2, 2))
    dy = (bottom - top) / len(genes)
    for i, gene in enumerate(genes):
        yy = top + (i + .5) * dy
        s.text(x + 7, yy - 2.4, gene, size=5.0, italic=True)
        for k in range(len(columns)):
            cx0 = plot_left + k * (col_w + gap)
            cx1 = cx0 + col_w
            s.line(cx0, yy, cx1, yy, color=COL["light"], width=.25)
    axis_y = y + h - 24
    for k in range(len(columns)):
        cx0 = plot_left + k * (col_w + gap)
        cx1 = cx0 + col_w
        s.line(cx0, axis_y, cx1, axis_y, color=COL["mid"], width=.3)
        for tick in (-4, 0, 4):
            px = cx0 + (tick - lo) / (hi - lo) * (cx1 - cx0)
            s.text(px, axis_y + 2, str(tick), size=5.0, align="center")
    s.text(x + w / 2, y + h - 13, "log2 fold change in offset-centred marginal expected abundance", size=5.0, align="center")
    s.text(
        x + 7,
        y + h - 5,
        "blank: cell-aligned target models + joint HCA/target draws required · 5 source genes not covered",
        size=5.0,
        color=COL["mid"],
    )


def draw_savi_density_blank_compact(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    """Final-size vector slots for the requested comprehensive/personalised densities."""
    planned_tag(s, x + w - 3, y + 1, "VECTOR · CELL-ALIGNED RERUN NEEDED")
    labels = ("external + internal direction", "study-control/HCA shift", "largest subject influence")
    gap = 13
    tile_w = (w - 30 - 2 * gap) / 3
    for i, label in enumerate(labels):
        xx = x + 15 + i * (tile_w + gap)
        s.text(xx + tile_w / 2, y + 15, label, size=4.8, bold=True, align="center")
        for row, row_label in enumerate(("comprehensive cohort query", "metadata-conditioned healthy predictive query")):
            py = y + 23 + row * 31
            s.rect(xx, py, tile_w, 25, fill=COL["white"], stroke=COL["light"], width=.35)
            s.text(xx + 4, py + 3, row_label, size=5.0, color=COL["mid"])
            if row == 0:
                s.line(xx + 8, py + 18, xx + tile_w - 5, py + 18, color=COL["mid"], width=.35)
                s.text(xx + tile_w - 4, py + 18, "offset-zero log2 marginal expected count", size=4.1, color=COL["mid"], align="right")
            else:
                for j, patient in enumerate(("P1", "P2", "P4", "P5", "P6")):
                    ridge_y = py + 8 + j * 3.2
                    s.text(xx + 3, ridge_y - 1.7, patient, size=4.0, color=COL["mid"])
                    s.line(xx + 13, ridge_y, xx + tile_w - 5, ridge_y, color=COL["light"], width=.25)
                s.text(xx + tile_w - 4, py + 21, "offset-centred log2 count", size=4.1, color=COL["mid"], align="right")
    s.text(
        x + w / 2,
        y + h - 5,
        "gene names fixed by ordered non-duplicating rules after c/f · blank axes are planned analyses",
        size=5.0,
        align="center",
    )


def draw_savi_predictive_percentile_heatmap_compact(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    """Final-size slot for SAVI healthy-prediction percentiles."""
    dashed_frame(s, x, y, w, h)
    planned_tag(s, x + w - 3, y + 3, "AGE + CELL-ALIGNMENT RERUN NEEDED")
    genes = (
        "ADRB2", "PTGS2", "AQP9", "OTUD1", "PTX3", "TGIF1", "TSC22D2", "ZBTB43",
        "ZFP36L2", "JARID2", "KLF4", "TREM1", "ANKRD28", "ANXA1", "PLAUR", "IL1B",
    )
    sample_specs = (
        ("C10_17. Disease-associated monocytes", "10"),
        ("C11_17. Disease-associated monocytes", "11"),
        ("C1_17. Disease-associated monocytes", "1"),
        ("C2_17. Disease-associated monocytes", "2"),
        ("C7_17. Disease-associated monocytes", "7"),
        ("C8-aunt-P1-STING_17. Disease-associated monocytes", "8"),
        ("C9-mother-P1-STING_17. Disease-associated monocytes", "9"),
        ("P1-STING-ht_17. Disease-associated monocytes", "1"),
        ("P2-STING-ht_17. Disease-associated monocytes", "2"),
        ("P4-STING-ht_17. Disease-associated monocytes", "4"),
        ("P5-STING_17. Disease-associated monocytes", "5"),
        ("P6-STING-ht_17. Disease-associated monocytes", "6"),
        ("P1-STING-ht-T_17. Disease-associated monocytes", "1R"),
        ("P1-STING-ht-T2_17. Disease-associated monocytes", "1T"),
        ("P2-STING-ht-T_17. Disease-associated monocytes", "2R"),
        ("P4-STING-ht-T_17. Disease-associated monocytes", "4R"),
        ("P6-STING-ht-T_17. Disease-associated monocytes", "6R"),
    )
    samples = tuple(item[0] for item in sample_specs)
    gx, gy = x + 43, y + 38
    heat_w, heat_h = w - 82, h - 60
    cw, ch = heat_w / len(samples), heat_h / len(genes)
    for start, count, color, label in (
        (0, 7, COL["ctrl"], "CTRL"),
        (7, 5, COL["savi"], "untreated"),
        (12, 5, COL["treated"], "treated"),
    ):
        bx = gx + start * cw
        s.rect(bx, y + 23, count * cw - .5, 3, fill=color)
        s.text(bx + count * cw / 2, y + 15, label, size=5.0, align="center")
    for c, (_, short) in enumerate(sample_specs):
        s.text(gx + (c + .5) * cw, y + 29, short, size=3.4, align="center")
    for r, gene in enumerate(genes):
        yy = gy + r * ch
        s.text(x + 5, yy, gene, size=3.8, italic=True)
        for c, _sample in enumerate(samples):
            s.rect(
                gx + c * cw,
                yy,
                cw - .35,
                ch - .22,
                fill=COL["white"],
                stroke=COL["light"],
                width=.25,
            )
    cx0, cx1 = gx + heat_w + 5, x + w - 5
    s.text((cx0 + cx1) / 2, y + 29, "7 CTRLs", size=3.8, bold=True, align="center")
    s.line((cx0 + cx1) / 2, gy, (cx0 + cx1) / 2, gy + heat_h, color=COL["mid"], width=.35, dash=(2, 2))
    for r, _gene in enumerate(genes):
        yy = gy + (r + .5) * ch
        s.line(cx0, yy, cx1, yy, color=COL["light"], width=.25)
    axis_y = gy + heat_h + 3
    s.line(cx0, axis_y, cx1, axis_y, color=COL["mid"], width=.3)
    for value in (0, .5, 1):
        xx = cx0 + value * (cx1 - cx0)
        s.line(xx, axis_y, xx, axis_y + 2, color=COL["mid"], width=.3)
        s.text(xx, axis_y + 2, f"{value:g}", size=4.5, align="center")
    s.text(
        gx + heat_w / 2,
        gy + heat_h / 2 - 2,
        "blank: corrected ages + cell-aligned predictions required",
        size=4.0,
        color=COL["mid"],
        align="center",
    )
    s.text(
        x + 5,
        y + h - 13,
        "IDs shown · R ruxolitinib · T tofacitinib · right: 7-control min–max + median",
        size=3.4,
        color=COL["mid"],
    )
    s.text(
        x + 5,
        y + h - 6,
        "tile = fraction of same-gene matched-healthy simulations below observed count",
        size=3.5,
    )


def draw_savi_influence_compact_blank(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    """Two-estimand subject-influence slot; full 16-gene grid goes to Extended Data."""
    dashed_frame(s, x, y, w, h)
    planned_tag(s, x + w - 3, y + 3, "TARGET-ONLY DELETIONS NEEDED")
    labels = ("16-gene subset",) + tuple(f"selected gene {i}" for i in range(1, 6))
    plot_left = x + 60
    gap = 12
    facet_w = (w - 75 - gap) / 2
    top, bottom = y + 47, y + h - 31
    s.text(x + 7, y + 19, "rows", size=4.6, bold=True)
    facet_specs = (
        ("untreated / unseen-study HCA", "log2(untreated / unseen HCA)"),
        ("untreated / study control", "log2(untreated / control)"),
    )
    for k, (title, _axis_label) in enumerate(facet_specs):
        x0 = plot_left + k * (facet_w + gap)
        x1 = x0 + facet_w
        if "unseen" in title:
            line1, line2 = "untreated /", "unseen-study HCA"
        else:
            line1, line2 = "untreated /", "study control"
        s.text((x0 + x1) / 2, y + 14, line1, size=4.6, bold=True, align="center")
        s.text((x0 + x1) / 2, y + 20, line2, size=4.6, bold=True, align="center")
        zero = x0 + (4 / 11) * facet_w
        s.line(zero, top, zero, bottom, color=COL["mid"], width=.4, dash=(2, 2))
    s.text(x + 7, y + 29, "full: diamond + interval · deletion: open point", size=4.1, color=COL["mid"])
    s.text(x + 7, y + 36, "thin: min–max · purple: direction category changed", size=4.1, color=COL["mid"])
    dy = (bottom - top) / len(labels)
    for i, label in enumerate(labels):
        yy = top + (i + .5) * dy
        s.text(x + 7, yy - 2.4, label, size=5.0, italic=True)
        for k in range(2):
            x0 = plot_left + k * (facet_w + gap)
            s.line(x0, yy, x0 + facet_w, yy, color=COL["light"], width=.25)
    for k, (_title, axis_label) in enumerate(facet_specs):
        x0 = plot_left + k * (facet_w + gap)
        compact_axis_x(s, x0, x0 + facet_w, y + h - 30, -4, 7, (-4, 0, 4), axis_label)
    s.text(
        x + 7,
        y + h - 12,
        "12 individual deletions in both effects · categories: positive / negative / unresolved",
        size=4.0,
    )
    s.text(x + 7, y + h - 5, "family deletions + all 16 genes: Extended Data", size=4.0, color=COL["mid"])


def draw_savi() -> None:
    s = PublicationSurface(SAVI_PDF, SAVI_PREVIEW, height=170 * MM)
    report = ROOT / "dev" / "SAVI" / "SAVI_case_study_universal_baseline.html"
    gene_rows = read_reactable_table(
        report,
        {"gene", "Category", "baseline_mean_log_mu", "cohort_log_mu", "cohort_se"},
    )
    composition_rows = read_reactable_table(
        report,
        {"paper_compartment", "Category", "cohort_mean", "baseline_mean"},
    )

    # a: observed source design and expression-reference support.
    ya, ha = 5, 78
    panel_label(s, "a", ya)
    draw_savi_design_support_compact(s, 25, ya + 5, 488, ha - 8)

    # b/c: compact composition check beside the signature-first scientific centre.
    ybc, hbc = 91, 144
    panel_label(s, "b", ybc)
    draw_savi_composition_compact(s, 25, ybc + 5, 185, hbc - 8, composition_rows)
    s.text(216, ybc + 1, "c", size=7, bold=True)
    draw_savi_signature_compact(s, 235, ybc + 5, 278, hbc - 8, gene_rows)

    # d: final-size vector density slots; current raster crops are not reused.
    yd, hd = 243, 92
    panel_label(s, "d", yd)
    draw_savi_density_blank_compact(s, 5, yd, 500, hd)

    # e/f: metadata-conditioned target-domain diagnostic and two-estimand influence.
    yef, hef = 343, 133
    panel_label(s, "e", yef)
    draw_savi_predictive_percentile_heatmap_compact(s, 25, yef + 5, 235, hef - 8)
    s.text(266, yef + 1, "f", size=7, bold=True)
    draw_savi_influence_compact_blank(s, 285, yef + 5, 228, hef - 8)

    s.finish()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def draw_hcc_case_structure(s: Surface, x: float, y: float, w: float, h: float) -> None:
    """Actual four-gene IF/TC design; no unused compartment taxonomy."""
    # The figure starts at the measurements that actually enter the analysis.
    sx = x + 8
    boxes = (
        (sx, 91, "20 GeoMx ROI", "measurements"),
        (sx + 113, 105, "aggregate within", "patient × IF/TC"),
        (sx + 240, 103, "5 spatial profiles", "3 unique patients"),
    )
    for bx, bw, line1, line2 in boxes:
        s.rect(bx, y + 12, bw, 31, fill=COL["pale"], stroke=COL["light"], width=.45)
        s.text(bx + bw / 2, y + 20, line1, size=5.1, bold=True, align="center")
        s.text(bx + bw / 2, y + 29, line2, size=4.8, color=COL["mid"], align="center")
    for ax in (sx + 99, sx + 226):
        s.line(ax, y + 27.5, ax + 10, y + 27.5, color=COL["mid"], width=.7)
        s.line(ax + 7, y + 24.5, ax + 10, y + 27.5, color=COL["mid"], width=.7)
        s.line(ax + 7, y + 30.5, ax + 10, y + 27.5, color=COL["mid"], width=.7)

    # Patient/compartment table with the number of ROI measurements contributing.
    px = x + 371
    s.text(px + 31, y + 4, "IF", size=5.0, bold=True, align="center")
    s.text(px + 74, y + 4, "TC", size=5.0, bold=True, align="center")
    patients = [
        ("A24  41F", "non-REC", 2, 3),
        ("A6  57M", "non-REC", 7, None),
        ("A29  57M", "REC", 4, 4),
    ]
    for i, (patient, outcome, nif, ntc) in enumerate(patients):
        yy = y + 17 + i * 14
        color = COL["rec"] if outcome == "REC" else COL["nonrec"]
        s.text(px, yy - 3, patient, size=5.0, bold=True)
        s.circle(px + 31, yy, 3.8, color, None)
        s.text(px + 31, yy - 3, str(nif), size=5.0, color=COL["white"], bold=True, align="center")
        if ntc is not None:
            s.circle(px + 74, yy, 3.8, color, None)
            s.text(px + 74, yy - 3, str(ntc), size=5.0, color=COL["white"], bold=True, align="center")
        else:
            s.text(px + 74, yy - 3, "—", size=5.0, color=COL["mid"], align="center")
    legend_mark(s, x + 373, y + h - 10, "circle", COL["nonrec"], "non-REC")
    legend_mark(s, x + 435, y + h - 10, "circle", COL["rec"], "REC")
    s.text(x + 8, y + h - 9, "SPON2 · VIM · ZFP36 · ZFP36L2", size=5.0, italic=True)
    s.text(x + 205, y + h - 9, "patient = biological unit · symbols show ROI count", size=4.5)


def draw_hcc_positive_control(s: Surface, x: float, y: float, w: float, h: float, ground: list[dict[str, str]]) -> None:
    positive = [r for r in ground if r["comparison_kind"] == "recurrence_direction"]
    values = {(r["gene"], r["compartment"]): float(r["observed_metric"]) for r in positive}
    genes = ["SPON2", "VIM", "ZFP36", "ZFP36L2"]
    x0, x1, lo, hi = x + 78, x + w - 12, 0, 4.5
    xx = lambda v: x0 + v / hi * (x1 - x0)
    xleg = legend_mark(s, x + 87, y + 10, "circle", COL["nonrec"], "IF")
    legend_mark(s, xleg + 10, y + 10, "square", COL["rec"], "TC")
    s.text(x + w - 8, y + 7, "extraction: 8/8 source directions", size=4.8, bold=True, align="right")
    dy = (h - 42) / len(genes)
    for i, gene in enumerate(genes):
        yy = y + 25 + (i + .5) * dy
        s.text(x + 8, yy - 3, gene, size=5.0, italic=True)
        s.line(x0, yy, x1, yy, color=COL["light"], width=0.3)
        vif, vtc = values[(gene, "IF")], values[(gene, "TC")]
        s.line(xx(0), yy - 3, xx(vif), yy - 3, color=COL["nonrec"], width=0.7)
        s.line(xx(0), yy + 3, xx(vtc), yy + 3, color=COL["rec"], width=0.7)
        s.circle(xx(vif), yy - 3, 2.1, COL["nonrec"], None)
        s.square(xx(vtc), yy + 3, 2.0, COL["rec"])
    axis_x(s, x0, x1, y + h - 16, lo, hi, [0, 2, 4], "log2 fold change (non-REC / REC)")


def draw_hcc_comprehensive_density_blueprint(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    """Four-gene cohort densities against the comprehensive HCA reference."""
    dashed_frame(s, x, y, w, h)
    planned_tag(s, x + w - 3, y + 3, "DENSITIES TO REGENERATE")
    genes = ("SPON2", "VIM", "ZFP36", "ZFP36L2")
    x0, x1 = x + 55, x + w - 9
    top, bottom = y + 28, y + h - 22
    for i, gene in enumerate(genes):
        yy = top + (i + .5) * (bottom - top) / len(genes)
        s.text(x + 7, yy - 2.5, gene, size=5.0, italic=True)
        s.line(x0, yy, x1, yy, color=COL["light"], width=.3)
    s.text(x + 8, y + 13, "grey density: comprehensive healthy-liver NK reference", size=4.7)
    s.text(x + 8, y + 20, "marks: IF/TC × non-REC/REC cohort means", size=4.7)
    compact_axis_x(s, x0, x1, y + h - 21, 2, 10, (2, 6, 10), "shared offset-centred log expected abundance")
    s.text(x + 7, y + h - 5, "blank: regenerate native-vector posterior densities on the frozen common scale", size=4.3, color=COL["mid"])


def draw_hcc_personalised_density_blueprint(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    """All four genes and five observed patient-by-compartment profiles."""
    dashed_frame(s, x, y, w, h)
    planned_tag(s, x + w - 3, y + 3, "20 DENSITIES TO REGENERATE")
    genes = ("SPON2", "VIM", "ZFP36", "ZFP36L2")
    profiles = ("A24 IF", "A6 IF", "A29 IF", "A24 TC", "A29 TC")
    gx, gy = x + 50, y + 31
    grid_w, grid_h = w - 59, h - 55
    cw, ch = grid_w / len(profiles), grid_h / len(genes)
    for c, profile in enumerate(profiles):
        color = COL["rec"] if profile.startswith("A29") else COL["nonrec"]
        s.text(gx + (c + .5) * cw, y + 19, profile, size=4.4, bold=True, color=color, align="center")
    for r, gene in enumerate(genes):
        yy = gy + r * ch
        s.text(x + 7, yy + ch / 2 - 2.5, gene, size=5.0, italic=True)
        for c in range(len(profiles)):
            xx = gx + c * cw
            s.rect(xx, yy, cw - .8, ch - .7, fill=COL["white"], stroke=COL["light"], width=.3)
            s.line(xx + 4, yy + ch / 2, xx + cw - 5, yy + ch / 2, color=COL["light"], width=.3)
    s.text(x + 7, y + h - 13, "cell: metadata-conditioned healthy predictive density + observed ROI-aggregated mark", size=4.3)
    s.text(x + 7, y + h - 6, "descriptive cross-platform query; no GeoMx calibration or predictive P value", size=4.3, color=COL["mid"])


def draw_hcc_comprehensive_effect_blueprint(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
    cohort_rows: list[dict[str, str]],
) -> None:
    """Current four-cohort point estimates; joint posterior intervals remain blank."""
    dashed_frame(s, x, y, w, h)
    planned_tag(s, x + w - 3, y + 3, "JOINT INTERVALS NEEDED")
    genes = ("SPON2", "VIM", "ZFP36", "ZFP36L2")
    cohorts = (
        ("IF_non_REC", "circle", COL["nonrec"], -3.0),
        ("IF_REC", "circle", COL["rec"], -1.0),
        ("TC_non_REC", "square", COL["nonrec"], 1.0),
        ("TC_REC", "square", COL["rec"], 3.0),
    )
    index = {(r["gene"], r["cohort"]): r for r in cohort_rows}
    x0, x1, lo, hi = x + 81, x + w - 14, -3.0, 1.5
    xx = lambda v: x0 + (v - lo) / (hi - lo) * (x1 - x0)
    top, bottom = y + 30, y + h - 27
    s.line(xx(0), top, xx(0), bottom, color=COL["mid"], width=.5, dash=(2, 2))
    legend_x = x + 13
    for label, mark, color in (
        ("IF non-REC", "circle", COL["nonrec"]),
        ("IF REC", "circle", COL["rec"]),
        ("TC non-REC", "square", COL["nonrec"]),
        ("TC REC", "square", COL["rec"]),
    ):
        legend_x = legend_mark(s, legend_x, y + 15, mark, color, label) + 8
    dy = (bottom - top) / len(genes)
    for i, gene in enumerate(genes):
        centre_y = top + (i + .5) * dy
        s.text(x + 8, centre_y - 2.5, gene, size=5.0, italic=True)
        s.line(x0, centre_y, x1, centre_y, color=COL["light"], width=.3)
        for cohort, mark, color, offset in cohorts:
            value = float(index[(gene, cohort)]["delta_log_mu_vs_baseline"])
            yy = centre_y + offset
            if mark == "circle":
                s.circle(xx(value), yy, 2.0, color, None)
            else:
                s.square(xx(value), yy, 1.9, color)
    compact_axis_x(s, x0, x1, y + h - 26, lo, hi, (-3, 0, 1), "HCC cohort − comprehensive HCA log expected abundance")
    s.text(x + 8, y + h - 6, "points are audited estimates; final lines require joint HCC/HCA posterior draws", size=4.4, color=COL["mid"])


def draw_hcc_personalised_score(
    s: Surface,
    x: float,
    y: float,
    w: float,
    h: float,
    gene_rows: list[dict[str, str]],
) -> None:
    """Mean absolute predictive deviation with arithmetic leave-one-gene-out range."""
    sample_order = ("A24_IF", "A6_IF", "A29_IF", "A24_TC", "A29_TC")
    by_sample: dict[str, list[dict[str, str]]] = {sample: [] for sample in sample_order}
    for row in gene_rows:
        if row["sample_id"] in by_sample and row["gene"] in {"SPON2", "VIM", "ZFP36", "ZFP36L2"}:
            by_sample[row["sample_id"]].append(row)
    values: dict[str, tuple[float, float, float, str]] = {}
    for sample in sample_order:
        rows = by_sample[sample]
        z = [abs(float(row["predictive_z"])) for row in rows]
        full = sum(z) / len(z)
        loo = [(sum(z) - value) / (len(z) - 1) for value in z]
        values[sample] = (full, min(loo), max(loo), rows[0]["recurrence"])

    x0, x1, lo, hi = x + 67, x + w - 13, 0.0, 2.5
    xx = lambda v: x0 + (v - lo) / (hi - lo) * (x1 - x0)
    top, bottom = y + 38, y + h - 29
    s.text(x + 8, y + 8, "four-gene metadata-conditioned predictive deviation", size=5.0, bold=True)
    s.text(x + w - 8, y + 18, "thin line: leave-one-gene-out range", size=4.5, color=COL["mid"], align="right")
    dy = (bottom - top) / len(sample_order)
    for i, sample in enumerate(sample_order):
        yy = top + (i + .5) * dy
        full, low, high, recurrence = values[sample]
        color = COL["rec"] if recurrence == "REC" else COL["nonrec"]
        s.text(x + 8, yy - 2.5, sample.replace("_", " "), size=5.0, bold=True)
        s.line(x0, yy, x1, yy, color=COL["light"], width=.3)
        s.line(xx(low), yy, xx(high), yy, color=color, width=1.0)
        s.line(xx(low), yy - 2.2, xx(low), yy + 2.2, color=color, width=.45)
        s.line(xx(high), yy - 2.2, xx(high), yy + 2.2, color=color, width=.45)
        s.circle(xx(full), yy, 2.4, color, None)
    compact_axis_x(s, x0, x1, y + h - 28, lo, hi, (0, 1, 2), "mean absolute predictive z")
    s.text(x + 8, y + h - 6, "descriptive: 3 patients, 1 REC · range is not a confidence interval", size=4.4, color=COL["mid"])


def draw_hcc() -> None:
    s = PublicationSurface(HCC_PDF, HCC_PREVIEW, height=170 * MM)
    data_dir = ROOT / "dev" / "Jie_HCC" / "data" / "processed"
    ground = read_csv(data_dir / "Jie_HCC_IF_TC_liver_baseline_ground_truth_comparison.csv")
    cohort_rows = read_csv(data_dir / "Jie_HCC_four_cohort_comprehensive_liver_baseline_mu_tests.csv")
    personalised_rows = read_csv(data_dir / "Jie_HCC_personalised_baseline_gene_patient_tests.csv")

    # a: only the ROI measurements, profiles, outcomes and genes actually used.
    ya, ha = 5, 76
    panel_label(s, "a", ya)
    draw_hcc_case_structure(s, 25, ya + 3, 488, ha - 4)

    # b: source-direction extraction check; c: comprehensive HCA density query.
    ybc, hbc = 88, 105
    panel_label(s, "b", ybc)
    draw_hcc_positive_control(s, 25, ybc + 5, 235, hbc - 8, ground)
    s.text(266, ybc + 1, "c", size=7, bold=True)
    draw_hcc_comprehensive_density_blueprint(s, 285, ybc + 5, 228, hbc - 8)

    # d: formal comprehensive-baseline cohort contrasts; intervals still needed.
    yd, hd = 201, 115
    panel_label(s, "d", yd)
    draw_hcc_comprehensive_effect_blueprint(s, 25, yd + 5, 488, hd - 8, cohort_rows)

    # e/f: all personalised densities and the leave-one-gene-out score summary.
    yef, hef = 324, 152
    panel_label(s, "e", yef)
    draw_hcc_personalised_density_blueprint(s, 25, yef + 5, 235, hef - 8)
    s.text(266, yef + 1, "f", size=7, bold=True)
    draw_hcc_personalised_score(s, 285, yef + 5, 228, hef - 8, personalised_rows)

    s.finish()


if __name__ == "__main__":
    draw_savi()
    draw_hcc()
    print(SAVI_PDF)
    print(HCC_PDF)
