"""flow_RNA_RNA_inter: interaction filtering and sequence support.
Adds the final split into mRNA-associated and ncRNA-ncRNA records (EC).
The release-dependent counts are read from data/paper_numbers_<date>.json
(05_paper_numbers.py); the two RNAInter counts (raw, curated) are fixed.
Palette matches the draw.io schemas: cream input, orange intermediate,
green terminal."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

CREAM="#FFF2CC"; CREAM_E="#D6B656"
ORANGE="#FFE6CC"; ORANGE_E="#D79B00"
GREEN="#D5E8D4"; GREEN_E="#82B366"
BLUE="#DAE8FC"; BLUE_E="#6C8EBF"

import glob, json
# interaction and entity counts of the current release (paper_figures/05_paper_numbers.py)
N = json.load(open(sorted(glob.glob("data/paper_numbers_*.json"))[-1]))
I, E = N["interactions"], N["entities"]
print("release:", N["release"], I, E)

fig, ax = plt.subplots(figsize=(9.2, 7.4))
ax.set_xlim(0, 10); ax.set_ylim(0, 8.6); ax.axis("off")

def box(x, y, w, h, title, value, fc, ec, fs=13, vfs=15):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor=fc, edgecolor=ec, linewidth=1.6))
    ax.text(x, y+h*0.16, title, ha="center", va="center", fontsize=fs)
    ax.text(x, y-h*0.22, value, ha="center", va="center",
            fontsize=vfs, fontweight="bold")

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=17, linewidth=1.5, color="#444444"))

W, H = 5.6, 1.15
box(5, 7.85, W, H, "Raw RNA–RNA interactions in RNAInter v4.0", "9,484,609", CREAM, CREAM_E)
arrow(5, 7.27, 5, 6.63)
ax.text(5.12, 6.95, "filtering of ambiguous, non-RNA\nand unmappable records",
        ha="left", va="center", fontsize=9.5, color="#555555")

box(5, 6.05, W, H, "Curated interactions", "9,429,404", CREAM, CREAM_E)
arrow(5, 5.47, 5, 4.83)
ax.text(5.12, 5.15, "both partners resolved to\n$\\geq$ 1 transcript sequence",
        ha="left", va="center", fontsize=9.5, color="#555555")

box(5, 4.25, W, H, "Sequence-supported RNA–RNA interactions", f"{I['total']:,}", ORANGE, ORANGE_E)

# split
ax.plot([5, 5], [3.67, 3.25], color="#444444", linewidth=1.5)
ax.plot([2.6, 7.4], [3.25, 3.25], color="#444444", linewidth=1.5)
arrow(2.6, 3.25, 2.6, 2.72)
arrow(7.4, 3.25, 7.4, 2.72)

box(2.6, 2.15, 4.3, 1.15, "mRNA-associated", f"{I['mrna']:,}", GREEN, GREEN_E, fs=13, vfs=15)
ax.text(2.6, 1.42, f"{I['mrna_pct']:.1f} %", ha="center", va="center", fontsize=11, color="#555555")
box(7.4, 2.15, 4.3, 1.15, "ncRNA–ncRNA", f"{I['ncrna']:,}", BLUE, BLUE_E, fs=13, vfs=15)
ax.text(7.4, 1.42, f"{I['ncrna_pct']:.1f} %", ha="center", va="center", fontsize=11, color="#555555")

ax.text(5, 0.72, f"{E['total']:,} unique sequence-resolved RNA entities\n"
                 f"({E['mrna']:,} mRNA · {E['ncrna']:,} ncRNA)",
        ha="center", va="center", fontsize=10.5, color="#333333")

plt.tight_layout()
import os, sys
OUT = sys.argv[1] if len(sys.argv) > 1 else "flow_RNA_RNA_inter.png"
os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print("written:", OUT)
