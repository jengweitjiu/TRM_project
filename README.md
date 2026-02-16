# TRM_project

Single-cell RNA-seq analysis of tissue-resident memory T cells (TRM) in psoriatic skin, combining Scanpy preprocessing, pySCENIC gene regulatory network inference, and basin stability modeling to identify transcription factors that maintain the IL-22-high TRM attractor state.

## Data Source

- **GEO Accession**: [GSE204765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204765)
- **Samples**: BCH01 non-lesional (NL) vs lesional (L) skin biopsies
- **Format**: 10x scRNA-seq (h5ad)

## Project Structure

```
TRM_project/
├── raw_data/                  GEO h5ad files (not tracked by git)
├── processed_data/            T cell / TRM subsets, expression matrix, metadata
├── scripts/
│   ├── trm_scanpy.py          Scanpy pipeline: QC → T cell gating → TRM gating → UMAP
│   ├── setup_pyscenic.sh      Conda environment setup for pySCENIC
│   ├── IL22TRMPOU2F1.ipynb    Main analysis: basin stability & TF perturbation
│   ├── locking_tf_rank.ipynb  Delta AUC ranking of locking TFs
│   └── depth_analysis_simulation.ipynb  Depth-based analysis & power simulations
├── GRN_SCENIC/                Regulons, GRN, AUC matrix, motif/ranking databases
├── figures_results/           Figures (pdf/png/svg), tables, and result CSVs
└── manuscript/                Cell Systems submission materials
```

## Analysis Pipeline

1. **Preprocessing** (`trm_scanpy.py`)
   - Load NL and L h5ad from GEO
   - QC and module scoring (T cell, CD4, CD8, TRM core, circulating)
   - Gate T cells (TRAC+ or score > 0.2), then gate TRM (TRM-core high, circulating low)
   - Concatenate, HVG selection, PCA, UMAP, Leiden clustering
   - Differential expression: L vs NL (Wilcoxon)

2. **GRN Inference** (pySCENIC)
   - Regulon detection with SCENIC motif databases (hg38, v9)
   - AUC scoring per cell (`trm_auc.csv`)

3. **Basin Stability Modeling** (`IL22TRMPOU2F1.ipynb`)
   - Define IL-22-high attractor state from TRM metadata
   - Fit logistic regression classifier in regulon AUC space
   - Compute basin stability curves under noise perturbation
   - Single-TF leave-one-out: identify TFs most critical for stability
   - Pairwise TF removal: identify cooperative TF pairs (e.g., POU2F1 + CHD1/RBBP5)
   - Permutation test for synergy significance

4. **Depth-Based Reanalysis** (`depth_analysis_simulation.ipynb`)
   - Cross-validated attractor depth estimation
   - Power simulations for cooperative TF detection

## Key Findings

- **POU2F1** is a central locking TF for the IL-22-high TRM state
- Cooperative TF pairs (e.g., CHD1 + RBBP5) show synergistic effects on basin stability
- TRM in lesional skin exhibit distinct regulon activity compared to non-lesional skin

## Requirements

```
scanpy
anndata
pyscenic
scikit-learn
pandas
numpy
matplotlib
seaborn
```

Setup pySCENIC environment:
```bash
bash scripts/setup_pyscenic.sh
```

## Large Files (not tracked)

The following are excluded via `.gitignore` and must be obtained separately:

| File | Source |
|------|--------|
| `raw_data/GSM6190732_BCH01-NL1.h5ad` | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6190732) |
| `raw_data/GSM6190734_BCH01-L1.h5ad` | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6190734) |
| `GRN_SCENIC/hg38__refseq-r80__*_rankings.feather` | [SCENIC databases](https://resources.aertslab.org/cistarget/) |
| `GRN_SCENIC/motifs-v9-nr.hgnc-*.tbl` | [SCENIC databases](https://resources.aertslab.org/cistarget/) |
