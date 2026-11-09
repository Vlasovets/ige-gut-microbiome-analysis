#!/usr/bin/env python
"""
Module 06r - assemble the manuscript Figure 4 from its panel files, so the
composite is rebuilt by script whenever a panel changes (the previous composite
was assembled by hand and fell behind the panel fixes of the co-author review).

Layout (as the manuscript figure):
  left   A  alpha diversity                06o -> fig4A_alpha.png
         B  family composition             06n -> fig4B_family_composition.png
  right  C  LinDA volcanoes (2x2)          06k --set score --compact
            ANCOM-BC2 volcanoes (2x2)      06k --set score --compact --method ANCOMBC2
All panels are 300 dpi with the same font sizes; the taller column is scaled
down (never up) so both columns end at the same height.

Run the four panel scripts first. Runs in scont_plots.
Output: results/figures/review/fig4_composite.{png,pdf}
Usage:  python revision/analysis/06r_fig4_assemble.py
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

REV = Path(__file__).resolve().parent.parent
FIG = REV / "results" / "figures" / "review"
FONTS = Path.home() / ".local/share/fonts/roboto"
DPI = 300
NAVY = "#2a3f5f"
PANELS = {"A": "fig4A_alpha.png", "B": "fig4B_family_composition.png",
          "LinDA": "fig4_volcanoes_contrasts_score_compact.png",
          "ANCOM-BC2": "fig4_volcanoes_contrasts_score_ancombc2_compact.png"}
MARGIN, GAP, COL_GAP, LETTER_W, HEADER_H = 40, 70, 110, 90, 90

im = {k: Image.open(FIG / f).convert("RGB") for k, f in PANELS.items()}
letter_font = ImageFont.truetype(str(FONTS / "Roboto-700-normal.ttf"), 72)
header_font = ImageFont.truetype(str(FONTS / "Roboto-500-normal.ttf"), 50)


def scaled(img, factor):
    if abs(factor - 1) < 1e-3:
        return img
    return img.resize((round(img.width * factor), round(img.height * factor)), Image.LANCZOS)


# left column: A and B at a common width (they share the stratum x-axis)
wl = max(im["A"].width, im["B"].width)
left = [scaled(im[k], wl / im[k].width) for k in ("A", "B")]
# right column: the two volcano blocks at a common width, each under a header
wr = max(im["LinDA"].width, im["ANCOM-BC2"].width)
right = [scaled(im[k], wr / im[k].width) for k in ("LinDA", "ANCOM-BC2")]

hl = sum(i.height for i in left) + GAP
hr = sum(i.height for i in right) + 2 * HEADER_H + GAP
if hr > hl:   # shrink the taller column's images to match
    f = (hl - 2 * HEADER_H - GAP) / sum(i.height for i in right)
    right = [scaled(i, f) for i in right]
else:
    f = (hr - GAP) / sum(i.height for i in left)
    left = [scaled(i, f) for i in left]
h = max(sum(i.height for i in left) + GAP, sum(i.height for i in right) + 2 * HEADER_H + GAP)
wl, wr = max(i.width for i in left), max(i.width for i in right)

W = MARGIN + LETTER_W + wl + COL_GAP + LETTER_W + wr + MARGIN
H = MARGIN + h + MARGIN
canvas = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(canvas)

x0, y = MARGIN + LETTER_W, MARGIN
for letter, img in zip("AB", left):
    canvas.paste(img, (x0, y))
    draw.text((MARGIN, y), letter, font=letter_font, fill="black")
    y += img.height + GAP

x1, y = x0 + wl + COL_GAP + LETTER_W, MARGIN
draw.text((x1 - LETTER_W, y), "C", font=letter_font, fill="black")
for name, img in zip(("LinDA", "ANCOM-BC2"), right):
    draw.text((x1 + wr / 2, y + HEADER_H / 2), name, font=header_font, fill=NAVY, anchor="mm")
    canvas.paste(img, (x1, y + HEADER_H))
    y += HEADER_H + img.height + GAP

for ext in ("png", "pdf"):
    canvas.save(FIG / f"fig4_composite.{ext}", dpi=(DPI, DPI))
print(f"wrote {FIG}/fig4_composite.{{png,pdf}}  {W}x{H} px = {W / DPI:.1f} x {H / DPI:.1f} in;"
      f" right column scale {f:.2f}" if hr > hl else f" left column scale {f:.2f}")
