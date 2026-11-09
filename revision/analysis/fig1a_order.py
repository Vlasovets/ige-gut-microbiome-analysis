"""
Recompute the orderings behind published Figure 1A, step for step as in
notebooks/01_IgE_analysis.ipynb. Shared by 06f (proof against the published PNG)
and 06g (rebuilt figure).

  cell 17   sensitized = the 275 individuals with any non-zero IgE
  cell 18   complete-linkage Jaccard trees for individuals and for IgEs; the
            reordered binary matrix is the notebook's ige_clustered_275.csv
  cell 37   PCA of the IgE co-occurrence matrix of that reordered matrix
  cell 43   varimax of the first 3 loadings -> latent allergy components (LACs)
  cell 21   complete-linkage Euclidean tree on the LACs orders the individuals

Neither ige_clustered_275.csv nor latent_allergy_components.csv survives on
disk; both are rebuilt here. Check that they are right: every individual lands
at its saved KORA_data/cluster_df_latent.csv row_ix_org (asserted).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from sklearn.decomposition import PCA

KORA = Path("/home/itg/oleg.vlasovets/slr_example/KORA_data")


def varimax(Phi, gamma=1.0, q=20, tol=1e-6):
    """Verbatim from notebook cell 36."""
    p, k = Phi.shape
    R = np.eye(k)
    d = 0
    for _ in range(q):
        d_old = d
        L = Phi @ R
        u, s, vh = np.linalg.svd(Phi.T @ (L ** 3 - (gamma / p) * L @ np.diag(np.diag(L.T @ L))))
        R = u @ vh
        d = s.sum()
        if d_old != 0 and d / d_old < 1 + tol:
            break
    return Phi @ R


def recompute(kora: Path = KORA):
    ige = pd.read_csv(kora / "ige_kora_subset_all.csv", index_col="u3_16s_id",
                      low_memory=False)
    test = ige.drop(ige[ige.eq(0).all(axis=1)].index)

    row_order = sch.dendrogram(sch.linkage(test, method="complete", metric="jaccard"),
                               no_plot=True)["leaves"]
    ige_link = sch.linkage(test.T, method="complete", metric="jaccard")
    col_order = sch.dendrogram(ige_link, no_plot=True)["leaves"]
    vis = (test.iloc[row_order, col_order] > 0).astype(int)

    cl = pd.read_csv(kora / "cluster_df_latent.csv", index_col=0)
    pos = pd.Series(np.arange(len(vis)), index=vis.index)
    assert (pos.reindex(cl.index).values == cl.row_ix_org.values).all(), \
        "recomputed ige_clustered_275 does not match cluster_df_latent.row_ix_org"

    co = vis.T @ vis
    pca = PCA(n_components=len(co)).fit(co)
    lac = pd.DataFrame(vis.values @ varimax(pca.components_.T[:, :3]),
                       index=vis.index, columns=["LAC1", "LAC2", "LAC3"])
    # PCA/varimax signs are arbitrary and sklearn changed its svd_flip
    # convention in 1.5; orient each LAC so its scores are mostly positive,
    # as published (scores 0-3).
    lac = lac * np.sign(lac.sum())
    lac_link = sch.linkage(lac.values, method="complete", metric="euclidean")

    return dict(vis=vis, lac=lac, lac_link=lac_link, cluster=cl.cluster,
                ige_link=ige_link, ige_names=list(test.columns))


def leaf_ids(link, ids):
    return np.asarray(ids)[sch.dendrogram(link, no_plot=True)["leaves"]]


def apply_flips(link, flips):
    """Swap the two children of the given linkage rows. This changes only the
    left-right drawing order of the tree, never its topology or heights."""
    out = link.copy()
    for i in flips:
        out[i, [0, 1]] = out[i, [1, 0]]
    return out
