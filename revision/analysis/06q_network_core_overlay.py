#!/usr/bin/env python
"""
Module 06q - genus network figure with the prevalence core overlaid (R1).

Styled after the submitted paper's network figure (figures/ige_networks.png):
one row per partition, three panels per row

  non-IgE  |  IgE-sensitized  |  Differential network

all drawn on ONE shared layout, so every panel and every partition is directly
comparable.

  layout        union of edges over all five partitions. The main connected
                component is laid out with Kamada-Kawai, then spread with
                Fruchterman-Reingold so it fills the panel; small components
                and unconnected genera sit on an outer ring (as in the paper,
                where unconnected taxa surround the connected core)
  node fill     phylum, as in the paper
  core overlay  prevalence-core genera (>= 80% of the cohort) are drawn
                saturated with a dark outline; non-core genera pale
  node size     cohort mean relative abundance (sqrt scale)
  edges         SPIEC-EASI partial correlations; green positive, red negative,
                width ~ |association|
  differential  edges with nominal permutation p < 0.05, coloured by how the
                association changes from non-IgE to IgE (paper's categories);
                thick black-cased = significant after Benjamini-Hochberg (q < 0.05) over
                the edges present in at least one group. NB: the column is named
                q_adaptBH, but module 03 fills it with p.adjust(method = "BH").
                Nominal edges are drawn for orientation only: they do not
                survive multiplicity correction.

DA is NOT encoded on the nodes: no genus is differentially abundant with a
calibrated method (LinDA, ALDEx2; full cohort, adjusted, 5% prevalence). The
per-node DA statistics are in results/networks_full/network_node_annotation.csv
(module 03b).

Runs in scont_plots (cnsplots), Roboto font.
Output: results/figures/review/fig_network_core_<partitions joined by _>.{png,svg,pdf}
Usage:  python revision/analysis/06q_network_core_overlay.py [--partitions any_ige ...]
"""

import argparse
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib import font_manager as fm
from matplotlib.colors import to_rgb
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D

REV = Path(__file__).resolve().parent.parent
RES = REV / "results"
NETF = RES / "networks_full"
OUT = RES / "figures" / "review"
ALL_PARTS = ["any_ige", "multiple", "medium_score", "low_score", "cluster_I"]
ALPHA_Q = 0.05          # FDR threshold (BH over testable edges; see docstring)
ALPHA_NOMINAL = 0.05    # differential panel: edges drawn for orientation
SEED = 42

NAVY = "#2a3f5f"
PHYLUM = {"Firmicutes": "#2a9d8f", "Bacteroidota": "#d1495b",
          "Actinobacteriota": "#56b4e9", "Proteobacteria": "#7a2418",
          "Verrucomicrobiota": "#8e6cc0"}
OTHER_PHYLUM = "#9aa0a6"
POS, NEG = "#1a7f1a", "#e8312f"
# change classes, "non-IgE -> IgE"; Okabe-Ito, matching the paper's hues
CHANGE = {"+→0": ("#009e73", "+ → 0  (positive in non-IgE only)"),
          "0→+": ("#e6c619", "0 → +  (positive in IgE only)"),
          "−→0": ("#0072b2", "− → 0  (negative in non-IgE only)"),
          "0→−": ("#d55e00", "0 → −  (negative in IgE only)"),
          "+→+": ("#cc79a7", "+ → +  (strength changed)"),
          "−→−": ("#56b4e9", "− → −  (strength changed)"),
          "+→−": ("#000000", "+ → −  (sign reversed)"),
          "−→+": ("#999999", "− → +  (sign reversed)")}
PART_LABEL = {"any_ige": "any IgE sensitization", "multiple": "multiple sensitizations",
              "medium_score": "medium sensitization score", "low_score": "low sensitization score",
              "cluster_I": "cluster I"}

ap = argparse.ArgumentParser()
ap.add_argument("--partitions", nargs="+", default=["any_ige"], choices=ALL_PARTS)
PARTS = ap.parse_args().partitions

# ── data ──────────────────────────────────────────────────────────────────
edges = {p: pd.read_csv(NETF / p / f"edges_{p}.csv") for p in ALL_PARTS}
for p, e in edges.items():
    assert (e.design == "full").all(), p
ann_all = pd.read_csv(NETF / "network_node_annotation.csv")
ann = ann_all[ann_all.partition == "any_ige"].set_index("genus")   # node attributes are cohort-level
assert len(ann) == 50, "annotation incomplete; rerun 03b"
ann["core80"] = ann.core80.astype(bool)
genera = list(ann.index)


def sign(x):
    return "+" if x > 0 else ("−" if x < 0 else "0")


# ── shared layout ─────────────────────────────────────────────────────────
G = nx.Graph()
G.add_nodes_from(genera)
for e in edges.values():
    u = e[e.present_ige | e.present_non_ige]
    G.add_edges_from(zip(u.taxon_i, u.taxon_j))
comps = sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)))
main = G.subgraph(comps[0])
p_main = nx.kamada_kawai_layout(main)
m = np.array([p_main[g] for g in main])
m = m - m.mean(0)
m = m / np.linalg.norm(m, axis=1).max()
# The densely connected centre packs into a knot under any force layout, so:
# (1) radial fisheye r -> r^0.6 opens the centre, keeping angular order;
r = np.linalg.norm(m, axis=1, keepdims=True)
m = m / np.maximum(r, 1e-9) * r ** 0.6
# (2) scale by the 85th-percentile radius (a few long leaves must not set the
#     scale) and pull the extreme leaves in;
r = np.linalg.norm(m, axis=1, keepdims=True)
m = m / np.quantile(r, 0.85) * 0.72
r = np.linalg.norm(m, axis=1, keepdims=True)
m = m / np.maximum(r, 1e-9) * np.minimum(r, 0.95)
# (3) deterministic overlap removal in final units: push apart pairs < D_MIN.
D_MIN = 0.155
for _ in range(600):
    d = m[:, None, :] - m[None, :, :]
    dist = np.linalg.norm(d, axis=2) + np.eye(len(m))
    push = np.clip(D_MIN - dist, 0, None) / 2
    np.fill_diagonal(push, 0)
    if not push.any():
        break
    m += (d / dist[..., None] * push[..., None]).sum(1)
p_main = dict(zip(main, m))
R_RING = np.linalg.norm(m, axis=1).max() + 0.17
pos = dict(zip(p_main, m))
# remaining genera on an outer ring; members of one small component stay adjacent,
# then by phylum so colours cluster
rest = [sorted(c, key=lambda g: (str(ann.at[g, "phylum"]), g)) for c in comps[1:]]
rest = sorted(rest, key=lambda c: (str(ann.at[c[0], "phylum"]), -len(c)))
ring = [g for c in rest for g in c]
ang = np.pi / 2 + 2 * np.pi * (np.arange(len(ring)) + 0.5) / len(ring)
for g, t in zip(ring, ang):
    pos[g] = np.array([R_RING * np.cos(t), R_RING * np.sin(t)])
xy = np.array([pos[g] for g in genera])

# ── style ─────────────────────────────────────────────────────────────────
cns.setup_matplotlib()
for f in Path.home().joinpath(".local/share/fonts/roboto").glob("*.ttf"):
    fm.fontManager.addfont(str(f))
plt.rcParams.update({"font.family": ["Roboto", "DejaVu Sans"], "savefig.transparent": False,
                     "savefig.facecolor": "white", "figure.facecolor": "white",
                     "axes.facecolor": "white", "font.size": 9, "svg.fonttype": "none"})


def pale(c, w=0.35):
    """Blend a colour towards white (opaque, so edges do not show through)."""
    r = np.array(to_rgb(c))
    return tuple(w * r + (1 - w) * np.ones(3))


base = [PHYLUM.get(ann.at[g, "phylum"], OTHER_PHYLUM) for g in genera]
fill = [c if ann.at[g, "core80"] else pale(c) for g, c in zip(genera, base)]
edgec = [NAVY if ann.at[g, "core80"] else "white" for g in genera]
lw_node = [1.1 if ann.at[g, "core80"] else 0.5 for g in genera]
size = 40 + 380 * np.sqrt(ann.loc[genera, "mean_rel_abundance"] / ann.mean_rel_abundance.max())


def edge_w(a):
    return 0.4 + abs(a) * 14


def draw_nodes(ax, bold):
    ax.scatter(xy[:, 0], xy[:, 1], s=size, c=fill, edgecolors=edgec,
               linewidths=lw_node, zorder=3)


LABEL_KW = dict(fontsize=6.5, color="#333333", zorder=5,
                path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
R_OUTER = 0.80          # nodes beyond this radius get radial, side-aware labels


def label_text(g, wrap=False):
    """Long names on the outer ring wrap onto two lines so they stay in the panel."""
    t = g.replace("_", " ")
    if ((wrap and len(t) > 16) or len(t) > 24) and " " in t:
        sp = [i for i, c in enumerate(t) if c == " "]
        cut = min(sp, key=lambda i: abs(i - len(t) / 2))
        t = t[:cut] + "\n" + t[cut + 1:]
    return t


def relax_outer(items, line_h=0.062, iters=200):
    """items: [genus, x, y, n_lines] on one side. Spread labels vertically so
    neighbours do not overlap; deterministic, keeps the order by y."""
    items = sorted(items, key=lambda t: t[2])
    ys = np.array([t[2] for t in items], float)
    hs = np.array([t[3] * line_h for t in items])
    for _ in range(iters):
        moved = False
        for a in range(len(ys) - 1):
            need = (hs[a] + hs[a + 1]) / 2 + 0.008
            gap = ys[a + 1] - ys[a]
            if gap < need:
                d = (need - gap) / 2
                ys[a] -= d; ys[a + 1] += d
                moved = True
        if not moved:
            break
    return {t[0]: (t[1], y) for t, y in zip(items, ys)}


def place_labels(ax):
    """Outer genera: label placed radially on the outer side (long names
    wrapped), then spaced vertically per side. Inner genera: de-overlapped with
    adjustText, outer labels as obstacles. Returns {genus: (x, y, ha, x0, y0)},
    (x0, y0) = the undisplaced position, used to decide on a leader line.
    Computed once and reused in every panel: the layout is shared."""
    out, outer_txt = {}, []
    radius = np.linalg.norm(xy, axis=1)
    sides = {"left": [], "right": []}
    for g, (x, y), r in zip(genera, xy, radius):
        if r > R_OUTER:
            lx, ly = np.array([x, y]) * (1 + 0.055 / r)
            ha = "left" if x >= 0 else "right"
            sides[ha].append((g, lx, ly, label_text(g, True).count("\n") + 1))
            out[g] = (lx, ly, ha, lx, ly)
    for ha, items in sides.items():
        for g, (lx, ly) in relax_outer(items).items():
            x0, y0 = out[g][3], out[g][4]
            out[g] = (lx, ly, ha, x0, y0)
            outer_txt.append(ax.text(lx, ly, label_text(g, True), ha=ha, va="center",
                                     linespacing=0.95, **LABEL_KW))
    inner = [g for g, r in zip(genera, radius) if r <= R_OUTER]
    texts = [ax.text(pos[g][0], pos[g][1] + 0.05, label_text(g), ha="center", va="center",
                     linespacing=0.95, **LABEL_KW) for g in inner]
    adjust_text(texts, x=xy[:, 0], y=xy[:, 1], objects=outer_txt, ax=ax,
                expand=(1.05, 1.2), force_text=(0.2, 0.35), force_static=(0.3, 0.5),
                force_pull=(0.02, 0.02), max_move=(20, 20), ensure_inside_axes=True)
    for g, t in zip(inner, texts):
        x, y = t.get_position()
        out[g] = (x, y, "center", pos[g][0], pos[g][1] + 0.05)
    for t in texts + outer_txt:
        t.remove()
    return out


def draw_labels(ax, lab_pos, bold):
    for g in genera:
        x, y, ha, x0, y0 = lab_pos[g]
        nx_, ny_ = pos[g]
        moved = np.hypot(x - x0, y - y0) > 0.04
        if (ha == "center" and np.hypot(x - nx_, y - ny_) > 0.075) or (ha != "center" and moved):
            ax.plot([nx_, x], [ny_, y], color="#9a9a9a", lw=0.45, zorder=4)
        ax.text(x, y, label_text(g, ha != "center"), ha=ha, va="center", linespacing=0.95,
                fontweight="bold" if g in bold else "normal", **LABEL_KW)


def frame(ax, title):
    ax.set_title(title, fontsize=11, fontweight="bold", color=NAVY, pad=4)
    h = R_RING + 0.08
    ax.set_xlim(-h - 0.5, h + 0.5); ax.set_ylim(-h, h)
    ax.set_aspect("equal"); ax.set_axis_off()


def draw_group(ax, e, present, asso):
    for r in e[e[present]].itertuples():
        (x1, y1), (x2, y2) = pos[r.taxon_i], pos[r.taxon_j]
        a = getattr(r, asso)
        ax.plot([x1, x2], [y1, y2], color=POS if a > 0 else NEG, lw=edge_w(a),
                alpha=0.85, zorder=2, solid_capstyle="round")


def draw_diff(ax, e):
    t = e[e.testable & (e.p_perm < ALPHA_NOMINAL)].copy()
    t["cls"] = [f"{sign(c)}→{sign(i)}" for c, i in zip(t.asso_non_ige, t.asso_ige)]
    t["fdr"] = t.q_adaptBH < ALPHA_Q
    for r in t.sort_values("fdr").itertuples():
        (x1, y1), (x2, y2) = pos[r.taxon_i], pos[r.taxon_j]
        col = CHANGE[r.cls][0]
        if r.fdr:
            ax.plot([x1, x2], [y1, y2], color="black", lw=6.0, zorder=2, solid_capstyle="round")
            ax.plot([x1, x2], [y1, y2], color=col, lw=4.0, zorder=2, solid_capstyle="round")
        else:
            ax.plot([x1, x2], [y1, y2], color=col, lw=1.6, alpha=0.9, zorder=2,
                    solid_capstyle="round")
    return set(t.cls)


# ── figure ────────────────────────────────────────────────────────────────
nrow = len(PARTS)
ROW_H, LEG_H = 4.4, 1.05
fig, axes = plt.subplots(nrow, 3, figsize=(16, ROW_H * nrow + LEG_H), squeeze=False)
classes_used = set()
rows = []
for k, p in enumerate(PARTS):
    e = edges[p]
    sig = e[e.testable & (e.q_adaptBH < ALPHA_Q)]
    n_ige, n_ctrl = int(e.n_ige.iloc[0]), int(e.n_ctrl.iloc[0])
    titles = [f"non-IgE (n = {n_ctrl}, {int(e.present_non_ige.sum())} edges)",
              f"IgE-sensitized (n = {n_ige}, {int(e.present_ige.sum())} edges)",
              f"Differential network ({len(sig)} edge{'s' * (len(sig) != 1)} at q < {ALPHA_Q})"]
    for j_, t in enumerate(titles):
        frame(axes[k, j_], t)
    if nrow > 1:
        axes[k, 0].text(-0.04, 1.0, "ABCDE"[k], transform=axes[k, 0].transAxes,
                        fontsize=20, color=NAVY, va="bottom")
    axes[k, 1].text(0.5, 1.075, PART_LABEL[p], transform=axes[k, 1].transAxes,
                    ha="center", va="bottom", fontsize=9, color="#555555", style="italic")
    rows.append((e, set(sig.taxon_i) | set(sig.taxon_j)))

ROW_FRAC = LEG_H / (ROW_H * nrow + LEG_H)
fig.tight_layout(rect=(0, ROW_FRAC, 1, 0.985), w_pad=0.3, h_pad=1.2)
# labels: identical layout in every panel, so de-overlap once on the final axes
draw_nodes(axes[0, 0], set())
lab_pos = place_labels(axes[0, 0])
for c in list(axes[0, 0].collections):
    c.remove()

for k, (e, bold) in enumerate(rows):
    draw_group(axes[k, 0], e, "present_non_ige", "asso_non_ige")
    draw_group(axes[k, 1], e, "present_ige", "asso_ige")
    classes_used.update(draw_diff(axes[k, 2], e))
    for ax in axes[k]:
        draw_nodes(ax, bold)
        draw_labels(ax, lab_pos, bold)

# grouped legend, as in the paper
dot = dict(marker="o", ls="", ms=9)
leg_node = [Line2D([], [], mfc=pale(PHYLUM["Firmicutes"]), mec="white", label="genus", **dot),
            Line2D([], [], mfc=PHYLUM["Firmicutes"], mec=NAVY, mew=1.3,
                   label="core genus (prevalence ≥ 80%)", **dot)]
leg_edge = [Line2D([], [], color=POS, lw=2.5, label="positive association"),
            Line2D([], [], color=NEG, lw=2.5, label="negative association")]
phyla_used = [ph for ph in PHYLUM if (ann.phylum == ph).any()]
leg_phy = [Line2D([], [], marker="s", ls="", ms=9, mfc=PHYLUM[ph], mec=PHYLUM[ph], label=ph)
           for ph in phyla_used]
leg_chg = [Line2D([], [], color=CHANGE[c][0], lw=2.5, label=CHANGE[c][1])
           for c in CHANGE if c in classes_used]
leg_chg.append(Line2D([], [], color="#bbbbbb", lw=3.5,
                      path_effects=[pe.Stroke(linewidth=5.5, foreground="black"), pe.Normal()],
                      label=f"significant, q < {ALPHA_Q} (Benjamini–Hochberg)"))
specs = [("Node:", leg_node, 1, 0.03), ("Edge:", leg_edge, 1, 0.215),
         ("Phylum:", leg_phy, 2, 0.36),
         (f"Differential edge, non-IgE → IgE (nominal p < {ALPHA_NOMINAL}):", leg_chg, 2, 0.585)]
leg_top = (LEG_H - 0.02) / (ROW_H * nrow + LEG_H)
for title, h, ncol, x in specs:
    lg = fig.legend(handles=h, title=title, loc="upper left", bbox_to_anchor=(x, leg_top),
                    ncol=ncol, frameon=False, fontsize=8.5, title_fontsize=10,
                    alignment="left", columnspacing=1.4, handlelength=2.0, handletextpad=0.6)
    lg.get_title().set_color(NAVY)


OUT.mkdir(parents=True, exist_ok=True)
stem = OUT / f"fig_network_core_{'_'.join(PARTS)}"
for ext in ("png", "svg", "pdf"):
    fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight")
print(f"wrote {stem}.{{png,svg,pdf}}  (core {int(ann.core80.sum())}/50; partitions {', '.join(PARTS)})")
