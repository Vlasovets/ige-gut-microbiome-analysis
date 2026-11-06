#!/usr/bin/env python
"""
Reconstruct clustering_validation.png (and .svg) from available data files.

Replicates the logic from notebook cell 094eeb10 of 01_IgE_analysis.ipynb.

Inputs (from KORA_data/):
  ige_kora_subset_all.csv  — IgE sensitization values (508 samples × 104 allergens)

Outputs (thesis/figures/abrogate/):
  clustering_validation.png
  clustering_validation.svg
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from numpy import eye, asarray, dot, sum, diag
from numpy.linalg import svd
from pathlib import Path

KORA_DATA = Path("/home/itg/oleg.vlasovets/slr_example/KORA_data")
OUT_DIR   = Path("/home/itg/oleg.vlasovets/thesis/figures/abrogate")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── varimax rotation ────────────────────────────────────────────────────────
def varimax(Phi, gamma=1.0, q=20, tol=1e-6):
    p, k = Phi.shape
    R = eye(k)
    d = 0
    for _ in range(q):
        d_old = d
        Lambda = dot(Phi, R)
        u, s, vh = svd(dot(Phi.T, asarray(Lambda)**3
                           - (gamma / p) * dot(Lambda, diag(diag(dot(Lambda.T, Lambda))))))
        R = dot(u, vh)
        d = sum(s)
        if d_old != 0 and d / d_old < 1 + tol:
            break
    return dot(Phi, R)


# ── clustering validation ────────────────────────────────────────────────────
def validate_latent_clustering(data, linkage_matrix, method='complete', metric='euclidean'):
    if isinstance(data, pd.DataFrame):
        data_array = data.values
    else:
        data_array = data

    distance_matrix = pdist(data_array, metric=metric)
    distance_matrix_square = squareform(distance_matrix)

    silhouette_scores, calinski_scores, davies_bouldin_scores = [], [], []
    inertias, valid_clusters, cluster_assignments = [], [], {}

    print("Evaluating cluster quality for k=2 to k=15 ...")
    for n_clusters in range(2, 16):
        try:
            cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            actual_n = len(np.unique(cluster_labels))
            print(f"  k={n_clusters}: {actual_n} actual clusters", end="")
            if actual_n < 2 or actual_n > n_clusters:
                print(" | skipped")
                continue
            cluster_assignments[n_clusters] = cluster_labels
            sil = silhouette_score(distance_matrix_square, cluster_labels, metric='precomputed')
            ch  = calinski_harabasz_score(data_array, cluster_labels)
            db  = davies_bouldin_score(data_array, cluster_labels)
            wcss = 0.0
            for cid in np.unique(cluster_labels):
                pts = data_array[cluster_labels == cid]
                if len(pts) > 1:
                    wcss += np.sum((pts - pts.mean(axis=0)) ** 2)
            silhouette_scores.append(sil)
            calinski_scores.append(ch)
            davies_bouldin_scores.append(db)
            inertias.append(float(wcss))
            valid_clusters.append(n_clusters)
            print(f" | Sil={sil:.4f}  CH={ch:.2f}  DB={db:.4f}  WCSS={wcss:.2f}")
        except Exception as e:
            print(f" | error: {e}")

    return (valid_clusters, silhouette_scores, calinski_scores,
            davies_bouldin_scores, inertias, cluster_assignments)


# ── load & prepare data ──────────────────────────────────────────────────────
print("Loading IgE data ...")
ige_data = pd.read_csv(KORA_DATA / "ige_kora_subset_all.csv",
                       sep=",", index_col="u3_16s_id", low_memory=False)

control_samples = ige_data[ige_data.eq(0).all(axis=1)]
test_samples    = ige_data.drop(control_samples.index)
print(f"  Total: {len(ige_data)}  |  test (IgE+): {len(test_samples)}  |  control: {len(control_samples)}")

# Binarize: IgE > 0 → 1 (presence/absence)
vis_ige = (test_samples > 0).astype(int)

# ── IgE hierarchical clustering (for column/row ordering) ───────────────────
print("Clustering IgE presence/absence ...")
row_linkage = sch.linkage(vis_ige.values, method='complete', metric='jaccard')
col_linkage = sch.linkage(vis_ige.values.T, method='complete', metric='jaccard')
row_order = sch.dendrogram(row_linkage, no_plot=True)['leaves']
col_order = sch.dendrogram(col_linkage, no_plot=True)['leaves']
vis_ige_clustered = vis_ige.iloc[row_order, col_order]

# ── PCA + varimax on co-occurrence matrix ────────────────────────────────────
print("Running PCA + varimax ...")
N = len(vis_ige_clustered)
co_occurrence = vis_ige_clustered.T.dot(vis_ige_clustered)
n_iges = len(co_occurrence)

pca = PCA(n_components=n_iges)
pca.fit(co_occurrence)
loadings_df = pd.DataFrame(
    pca.components_.T,
    index=co_occurrence.columns,
    columns=[f'PC{i+1}' for i in range(n_iges)]
)

# Varimax rotation of top 3 PCs
top_3_pcs = loadings_df.iloc[:, :3].values
sparse_loadings = varimax(top_3_pcs)

rotated_proj = vis_ige_clustered.values.dot(sparse_loadings)
latent_allergy_components = pd.DataFrame(
    rotated_proj,
    index=vis_ige_clustered.index,
    columns=['LA1', 'LA2', 'LA3']
)
print(f"  latent_allergy_components shape: {latent_allergy_components.shape}")

# ── linkage on latent components ─────────────────────────────────────────────
print("Computing linkage on latent components ...")
lac_arr = np.asarray(latent_allergy_components, dtype=np.float64)
lac_arr = np.nan_to_num(lac_arr, nan=0.0, posinf=0.0, neginf=0.0)
print(f"  NaN: {np.isnan(lac_arr).any()}  Inf: {np.isinf(lac_arr).any()}")
row_linkage_test = sch.linkage(lac_arr, method='complete', metric='euclidean')

# ── clustering validation ────────────────────────────────────────────────────
(cluster_range, silhouette_scores, calinski_scores,
 davies_bouldin_scores, inertias, cluster_assignments) = validate_latent_clustering(
    lac_arr, row_linkage_test, 'complete', 'euclidean')

if not cluster_range:
    print("No valid clustering results. Aborting.")
    raise SystemExit(1)

cluster_range      = np.array(cluster_range)
silhouette_scores  = np.array(silhouette_scores)
calinski_scores    = np.array(calinski_scores)
davies_bouldin_scores = np.array(davies_bouldin_scores)
inertias           = np.array(inertias)

# ── figure ───────────────────────────────────────────────────────────────────
print("Plotting ...")
plt.rcParams.update({'font.size': 11, 'font.weight': 'normal'})

fig = plt.figure(figsize=(20, 15))
gs  = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)

# 1. Dendrogram (top row, spanning all columns)
ax1 = fig.add_subplot(gs[0, :])
dendrogram(row_linkage_test, truncate_mode='level', p=15, ax=ax1)
ax1.set_title('Hierarchical Clustering Dendrogram — Latent Allergy Components', fontsize=14)
ax1.set_xlabel('Sample index')
ax1.set_ylabel('Euclidean distance')

nine_mask = cluster_range == 9
nine_idx  = int(np.where(nine_mask)[0][0]) if np.any(nine_mask) else None

# 2. Silhouette
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(cluster_range, silhouette_scores, 'ro-', linewidth=2, markersize=8)
ax2.set_xlabel('Number of clusters')
ax2.set_ylabel('Silhouette score')
ax2.set_title('Silhouette analysis')
ax2.grid(True, alpha=0.3)
best_idx = int(np.argmax(silhouette_scores))
ax2.plot(cluster_range[best_idx], silhouette_scores[best_idx],
         'g*', markersize=15, label=f'Best (k={cluster_range[best_idx]})')
if nine_idx is not None:
    ax2.plot(9, silhouette_scores[nine_idx], 'bo', markersize=12, label='k=9')
ax2.legend()

# 3. Elbow
ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(cluster_range, inertias, 'bo-', linewidth=2, markersize=8)
ax3.set_xlabel('Number of clusters')
ax3.set_ylabel('Within-cluster sum of squares')
ax3.set_title('Elbow method')
ax3.grid(True, alpha=0.3)
if nine_idx is not None:
    ax3.plot(9, inertias[nine_idx], 'ro', markersize=12, label='k=9')
    ax3.legend()
ax3.text(0.68, 0.82, 'Inspect visually\nfor elbow point',
         transform=ax3.transAxes, fontsize=10,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))

# 4. Calinski-Harabasz
ax4 = fig.add_subplot(gs[1, 2])
ax4.plot(cluster_range, calinski_scores, 'go-', linewidth=2, markersize=8)
ax4.set_xlabel('Number of clusters')
ax4.set_ylabel('Calinski-Harabasz index')
ax4.set_title('Calinski-Harabasz index\n(higher = better)')
ax4.grid(True, alpha=0.3)
best_idx = int(np.argmax(calinski_scores))
ax4.plot(cluster_range[best_idx], calinski_scores[best_idx],
         'r*', markersize=15, label=f'Best (k={cluster_range[best_idx]})')
if nine_idx is not None:
    ax4.plot(9, calinski_scores[nine_idx], 'bo', markersize=12, label='k=9')
ax4.legend()

# 5. Davies-Bouldin
ax5 = fig.add_subplot(gs[2, 0])
ax5.plot(cluster_range, davies_bouldin_scores, 'mo-', linewidth=2, markersize=8)
ax5.set_xlabel('Number of clusters')
ax5.set_ylabel('Davies-Bouldin index')
ax5.set_title('Davies-Bouldin index\n(lower = better)')
ax5.grid(True, alpha=0.3)
best_idx = int(np.argmin(davies_bouldin_scores))
ax5.plot(cluster_range[best_idx], davies_bouldin_scores[best_idx],
         'r*', markersize=15, label=f'Best (k={cluster_range[best_idx]})')
if nine_idx is not None:
    ax5.plot(9, davies_bouldin_scores[nine_idx], 'bo', markersize=12, label='k=9')
ax5.legend()

# 6. Combined normalised metrics
ax6 = fig.add_subplot(gs[2, 1])

def safe_normalize(arr):
    arr = np.array(arr, dtype=float)
    if arr.max() == arr.min():
        return np.full_like(arr, 0.5)
    return (arr - arr.min()) / (arr.max() - arr.min())

ax6.plot(cluster_range, safe_normalize(silhouette_scores), 'r-',  lw=2, label='Silhouette',             alpha=0.8)
ax6.plot(cluster_range, safe_normalize(calinski_scores),   'g-',  lw=2, label='Calinski-Harabasz',      alpha=0.8)
ax6.plot(cluster_range, 1 - safe_normalize(davies_bouldin_scores), 'm-', lw=2, label='Davies-Bouldin (inv)', alpha=0.8)
ax6.plot(cluster_range, 1 - safe_normalize(inertias),      'b-',  lw=2, label='Elbow (inv)',            alpha=0.8)
if 9 in cluster_range:
    ax6.axvline(x=9, color='orange', linestyle='--', lw=3, alpha=0.8, label='k=9')
ax6.set_xlabel('Number of clusters')
ax6.set_ylabel('Normalised score (0–1)')
ax6.set_title('Combined normalised metrics')
ax6.legend(fontsize=9)
ax6.grid(True, alpha=0.3)

# 7. Summary table
ax7 = fig.add_subplot(gs[2, 2])
ax7.axis('off')
summary_data = [
    [int(k),
     f"{silhouette_scores[i]:.4f}",
     f"{calinski_scores[i]:.2f}",
     f"{davies_bouldin_scores[i]:.4f}",
     f"{inertias[i]:.1f}"]
    for i, k in enumerate(cluster_range)
]
tbl = ax7.table(cellText=summary_data,
                colLabels=['k', 'Silhouette', 'Calinski-H', 'Davies-B', 'WCSS'],
                cellLoc='center', loc='center',
                bbox=[0.05, 0.05, 0.90, 0.90])
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 1.5)
if nine_idx is not None:
    for j in range(5):
        tbl[(nine_idx + 1, j)].set_facecolor('#FFE4B5')
ax7.set_title('Clustering metrics summary', fontsize=12, fontweight='bold')

fig.suptitle(
    'Clustering validation: latent allergy components (Euclidean distance)',
    fontsize=16, fontweight='bold'
)

for fmt in ('png', 'svg'):
    out_path = OUT_DIR / f"clustering_validation.{fmt}"
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved {out_path}")

plt.close(fig)
print("Done.")
