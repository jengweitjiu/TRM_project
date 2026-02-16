import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

sc.settings.verbosity = 2

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "raw_data"
OUT = BASE / "processed_data"
FIG = BASE / "figures_results"

FILES = [
    (RAW / "GSM6190732_BCH01-NL1.h5ad", "NL"),
    (RAW / "GSM6190734_BCH01-L1.h5ad",  "L"),
]

# ---------- helpers ----------
def find_gene(adata, candidates):
    """Return first gene present in adata.var_names (case-sensitive), else None."""
    for g in candidates:
        if g in adata.var_names:
            return g
    return None

def add_basic_qc(adata):
    # Common mt gene patterns: human MT-*, mouse mt-*
    mt_flag = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
    adata.var["mt"] = mt_flag
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

def ensure_raw_counts(adata):
    # If counts layer exists, use it; else keep X
    # We'll just proceed with X as-is; many GEO h5ad already has log1p normalized.
    pass

def preprocess_for_umap(adata):
    # If data already normalized/logged, running normalize_total again can distort.
    # We'll do a conservative path:
    # 1) select HVGs
    # 2) scale
    # 3) PCA/Neighbors/UMAP
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)

def score_sets(adata):
    # T cell / CD4 / CD8 / TRM modules (minimal, robust)
    tcell = ["TRAC", "CD3D", "CD3E"]
    cd4   = ["CD4", "IL7R"]
    cd8   = ["CD8A", "CD8B"]
    trm   = ["ITGAE", "CD69", "CXCR6", "ITGA1", "RGS1"]
    circ  = ["CCR7", "S1PR1", "SELL", "LTB"]

    # only keep genes that exist
    def filt(gs): return [g for g in gs if g in adata.var_names]
    sc.tl.score_genes(adata, filt(tcell), score_name="score_Tcell")
    sc.tl.score_genes(adata, filt(cd4),   score_name="score_CD4")
    sc.tl.score_genes(adata, filt(cd8),   score_name="score_CD8")
    sc.tl.score_genes(adata, filt(trm),   score_name="score_TRMcore")
    sc.tl.score_genes(adata, filt(circ),  score_name="score_Circ")

def gate_tcells(adata):
    # Soft gate: T cell score > threshold OR TRAC present >0
    # Threshold is dataset-dependent; start modest.
    trac = find_gene(adata, ["TRAC"])
    if trac is not None:
        keep = (adata.obs["score_Tcell"] > 0.2) | (adata[:, trac].X.toarray().flatten() > 0 if hasattr(adata[:, trac].X, "toarray") else adata[:, trac].X > 0)
    else:
        keep = (adata.obs["score_Tcell"] > 0.2)
    return adata[keep].copy()

def label_lineage(adata):
    # CD4-like vs CD8-like by module scores
    adata.obs["lineage"] = np.where(adata.obs["score_CD8"] > adata.obs["score_CD4"], "CD8-like", "CD4-like")

def gate_trm(adata):
    # Pragmatic TRM: TRMcore high + circulation low
    keep = (adata.obs["score_TRMcore"] > 0.2) & (adata.obs["score_Circ"] < 0.2)
    return adata[keep].copy()

# ---------- main ----------
adatas = []
for fn, cond in FILES:
    print("\n=== Loading:", fn, "===")
    a = sc.read_h5ad(str(fn))
    a.obs["condition"] = cond
    add_basic_qc(a)
    score_sets(a)
    a = gate_tcells(a)
    score_sets(a)  # re-score after subsetting
    label_lineage(a)
    adatas.append(a)

print("\n=== Concatenating NL + L T cells ===")
adata = ad.concat(adatas, join="outer", label="batch", keys=[x[1] for x in FILES], fill_value=0)
# fill NaNs in scores if any
for col in ["score_Tcell","score_CD4","score_CD8","score_TRMcore","score_Circ"]:
    if col in adata.obs:
        adata.obs[col] = adata.obs[col].astype(float)

print("Cells (T cells):", adata.n_obs, " Genes:", adata.n_vars)

# UMAP on T cells
preprocess_for_umap(adata)

# TRM subset
print("\n=== TRM gating ===")
trm = gate_trm(adata)
print("Cells (TRM):", trm.n_obs)

if trm.n_obs > 200:
    preprocess_for_umap(trm)

# Save outputs
adata.write(str(OUT / "BCH01_Tcells_NL_vs_L.h5ad"))
trm.write(str(OUT / "BCH01_TRM_NL_vs_L.h5ad"))

# Make quick plots
sc.settings.figdir = str(FIG)
sc.pl.umap(adata, color=["condition", "lineage", "score_TRMcore", "score_Circ"], wspace=0.4, save="_Tcells_overview.png")
if trm.n_obs > 0:
    cols = ["condition", "lineage", "score_TRMcore", "score_Circ"]
    sc.pl.umap(trm, color=cols, wspace=0.4, save="_TRM_overview.png")

print(f"\nDONE. Saved to {OUT}")
print(f"UMAP images saved under {FIG}")
# 重新分群 TRM
sc.pp.highly_variable_genes(trm, n_top_genes=3000, flavor="seurat_v3", subset=True)
sc.pp.scale(trm, max_value=10)
sc.tl.pca(trm)
sc.pp.neighbors(trm)
sc.tl.leiden(trm, resolution=0.5)
sc.tl.umap(trm)

sc.pl.umap(trm, color=["condition","lineage","leiden"], save="_TRM_clusters.png")
sc.tl.rank_genes_groups(trm, groupby="condition", method="wilcoxon")
sc.pl.rank_genes_groups(trm, n_genes=20, save="_TRM_L_vs_NL.png")
pd.crosstab(trm.obs["condition"], trm.obs["lineage"])

