#!/usr/bin/env python
"""
Module 06f - prove which ordering the published Figure 1A used, from the
published PNG itself.

The figure was built with dash-bio's Clustergram, and that code is not in the
repository. Reading the 104 x 275 presence/absence matrix back out of the
published heatmap and matching it to the recomputed data (fig1a_order.py)
shows:

  IgEs         = the notebook's complete-Jaccard dendrogram, leaf 0 at the BOTTOM
  individuals  = the notebook's complete-Euclidean LAC dendrogram with some
                 branches flipped (same tree, same heights, different
                 left-right drawing)

Writes revision/lib/fig1A_lac_flips.csv: the linkage rows whose children were
swapped. That file holds no participant IDs; 06g turns it back into the exact
order. Also checks the published LAC strip against the recomputed LAC scores.

Input:  the published figure PNG (manuscript repo), KORA_data
Output: revision/lib/fig1A_lac_flips.csv, results/figures/review/fig1A_order_qc.txt
Env:    scont_plots
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from PIL import Image
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig1a_order import apply_flips, leaf_ids, recompute  # noqa: E402

REV = Path(__file__).resolve().parent.parent
PUBLISHED = Path("/home/itg/oleg.vlasovets/slr_example/The-role-of-IgE-patterns-and-their-"
                 "link-to-the-gut-microbiome-in-allergic-sensitization/figures/png/"
                 "ige_complete_LA_euclidean.png")
OUT_FLIPS = REV / "lib" / "fig1A_lac_flips.csv"
OUT_QC = REV / "results" / "figures" / "review" / "fig1A_order_qc.txt"
PURPLE, YELLOW = np.array([68, 1, 84]), np.array([253, 231, 37])   # viridis ends
SEARCH = (900, 850, 8700, 4250)      # px window around panel A's heatmap
NR, NC = 104, 275

log = []
def say(*a):
    s = " ".join(map(str, a)); print(s); log.append(s)


def read_heatmap(img):
    a = np.asarray(img.crop(SEARCH)).astype(int)
    hm = (np.abs(a - PURPLE).sum(2) < 60) | (np.abs(a - YELLOW).sum(2) < 80)
    r = np.where(hm.mean(1) > 0.5)[0]; c = np.where(hm.mean(0) > 0.5)[0]
    x0, y0 = SEARCH[0] + c.min(), SEARCH[1] + r.min()
    x1, y1 = SEARCH[0] + c.max() + 1, SEARCH[1] + r.max() + 1
    a = np.asarray(img.crop((x0, y0, x1, y1))).astype(int)
    ys = ((np.arange(NR) + 0.5) * a.shape[0] / NR).astype(int)
    xs = ((np.arange(NC) + 0.5) * a.shape[1] / NC).astype(int)
    M = np.zeros((NR, NC), int)
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            w = a[y - 4:y + 5, x - 4:x + 5].reshape(-1, 3)   # 9x9 px, dodges the dotted cluster lines
            M[i, j] = int((np.abs(w - YELLOW).sum(1) < 120).mean() >
                          (np.abs(w - PURPLE).sum(1) < 120).mean())
    say(f"published heatmap at x {x0}-{x1}, y {y0}-{y1} px; {M.sum()} IgE-positive cells")
    return M, (x0, x1, y1)


def read_lac_strip(img, x0, x1, y_below):
    """Map each LAC-strip cell back onto the viridis scale (rank-level check only)."""
    from matplotlib import colormaps
    lut = (colormaps["viridis"](np.linspace(0, 1, 256))[:, :3] * 255)
    a = np.asarray(img.crop((x0, y_below, x1, y_below + 1200))).astype(np.int32)
    xs = ((np.arange(NC) + 0.5) * a.shape[1] / NC).astype(int)
    # distance to the nearest viridis colour, at cell-centre columns only
    d = np.min(((a[:, xs, None, :] - lut[None, None]) ** 2).sum(-1), axis=2)
    rows = np.where((d < 400).mean(1) > 0.9)[0]
    # the 3-row strip has thin anti-aliased gaps between rows: take the full span
    y0, y1 = rows.min(), rows.max() + 1
    ys = y0 + ((np.arange(3) + 0.5) * (y1 - y0) / 3).astype(int)
    V = np.zeros((3, NC))
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            px = a[y - 3:y + 4, x - 3:x + 4].reshape(-1, 3).mean(0)
            V[i, j] = np.argmin(((lut - px) ** 2).sum(1)) / 255
    return V


def main():
    Image.MAX_IMAGE_PIXELS = None          # the published PNG is 12917 x 8770
    img = Image.open(PUBLISHED).convert("RGB")
    r = recompute()
    vis, lac, link = r["vis"], r["lac"], r["lac_link"]
    say(f"recomputed: {vis.shape[0]} individuals x {vis.shape[1]} IgEs, "
        f"{int(vis.values.sum())} positive cells")

    M, (x0, x1, y1) = read_heatmap(img)

    # 1. IgE rows: leaf order, drawn bottom-up (plotly puts row 0 at the bottom).
    rows_top_down = list(vis.columns)[::-1]
    D = vis[rows_top_down].values
    mism = (M.T[:, None, :] != D[None, :, :]).sum(2)       # image column x individual
    say(f"IgE rows = Jaccard leaf order, bottom-up: every published column has an "
        f"exact match = {bool((mism.min(1) == 0).all())}")
    assert (mism.min(1) == 0).all()

    # 2. Individuals: rotate the LAC tree so its leaves follow the published
    #    columns. Individuals with identical profiles are interchangeable in the
    #    picture, so use the mean published position of each profile.
    prof = vis[rows_top_down].astype(str).agg("".join, axis=1)
    img_prof = pd.Series(["".join(map(str, c)) for c in M.T])
    key = {s: img_prof.index[img_prof == prof[s]].to_numpy().mean() for s in vis.index}
    ids = lac.index.to_numpy()
    n = len(ids)
    flips = []

    def place(node):
        if node < n:
            return [ids[node]]
        i = node - n
        a, b = place(int(link[i, 0])), place(int(link[i, 1]))
        if np.mean([key[s] for s in a]) > np.mean([key[s] for s in b]):
            flips.append(i); return b + a
        return a + b
    place(2 * n - 2)

    order = leaf_ids(apply_flips(link, flips), ids)
    diff = int((vis.loc[order, rows_top_down].values.T != M).sum())
    say(f"LAC tree with {len(flips)} of {n - 1} branches flipped -> "
        f"{diff} of {M.size} cells differ from the published heatmap")
    assert diff == 0
    default = leaf_ids(link, ids)
    say(f"(unflipped scipy order: "
        f"{int((vis.loc[default, rows_top_down].values.T != M).sum())} cells differ)")
    runs = "".join(pd.Series(r["cluster"].reindex(order).values).loc[
        lambda s: s != s.shift()].values)
    say(f"cluster sequence along published columns: {runs}")
    assert runs == "ABCDEFGHI"

    # 3. LAC strip: rank agreement with the recomputed scores, per LAC.
    V = read_lac_strip(img, x0, x1, y1 + 200)
    for k in range(3):
        rho = spearmanr(V[k], lac.loc[order].iloc[:, k]).correlation
        say(f"published LAC strip row {k + 1} vs recomputed LAC{k + 1}: Spearman {rho:.3f}")

    pd.DataFrame({"linkage_row": sorted(flips)}).to_csv(OUT_FLIPS, index=False)
    OUT_QC.parent.mkdir(parents=True, exist_ok=True)
    OUT_QC.write_text("\n".join(log) + "\n")
    print(f"wrote {OUT_FLIPS}\nwrote {OUT_QC}")


if __name__ == "__main__":
    main()
