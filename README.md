# BIOT 670I — Project 2a: Single Cell RNA-Seq Analysis

## Overview

This project reimplements the Scanpy legacy workflow on the 3k PBMC dataset as a single master pipeline script. The goal is to reproduce a standard scRNA-seq analysis, from raw count data through quality control, normalization, dimensionality reduction, clustering, and cell type annotation, in a reproducible and self-contained pipeline.

The dataset used is the 10x Genomics PBMC 3k dataset

## Project Structure

```text
BIFS-618---Project-2a-Single-Cell-RNA-Seq/
│
├── README.md
│
├── 00_raw_data/
│   ├── README.md
│   ├── qc_metrics.csv
│   └── figures/
│       ├── README.md
│       ├── highest_expr_genes.png
│       ├── scatter.png
│       └── violin.png
│
├── 01_clustering/
│   ├── README.md
│   └── figures/
│       ├── README.md
│       ├── hvg_plot.png
│       ├── pca_variance_ratio.png
│       └── umap_qc_overlay.png
│
├── 02_annotation/
│   ├── README.md
│   └── figures/
│       ├── README.md
│       ├── dotplot_markers.png
│       ├── umap_cell_type.png
│       └── umap_marker_genes.pdf
│
├── Draft_Code/
│   ├── README.md
│   ├── Draft Code.txt
│   ├── DraftCode_DialogBoxes.txt
│   ├── The_Gruesome_Details.txt
│   ├── draftcode_v1.txt
│   ├── draftcodedb.txt
│   ├── pipeline.sh
│   ├── scanpy_master.py
│   ├── scanpy_master_logs_seed_v2.py
│   └── scanpy_master_logs_v2.py
│
└── Pipeline/
    ├── README.md
    ├── scanpy_master_logs_seed_v3.py
    └── run_scanpy_pipeline.sh
```

## Pipeline Summary

| Step                     | Method                            | Input → Output                        |
| ------------------------ | --------------------------------- | ------------------------------------- |
| Quality Control          | Scanpy QC metrics                 | Raw matrix → filtered matrix          |
| Normalization            | Total count normalization + log1p | Filtered → normalized                 |
| HVG Selection            | Seurat-flavor dispersion          | 32,738 genes → 2,000 genes            |
| Dimensionality Reduction | PCA                               | 2,000 genes → 50 principal components |
| Neighbor Graph           | k-nearest neighbors               | 50 PCs → connectivity graph           |
| UMAP                     | Fuzzy set optimization            | Graph → 2D embedding                  |
| Clustering               | Leiden algorithm                  | Graph → cluster labels                |
| Cell Type Annotation     | Marker gene scoring               | Clusters → cell type labels           |

***

## Reproducibility

All stochastic steps in the pipeline are controlled using a shared seed value:

```python
SEED = 594
```

This seed is applied to:

*   NumPy
*   Python `random` module
*   Scrublet
*   PCA
*   Neighbor graph construction
*   UMAP
*   Leiden clustering

Each execution produces a `run_log.txt` file documenting execution metadata, package versions, seed value, and matrix dimensions at every stage of the pipeline, ensuring reproducibility and auditability.

***

## Group

**Group 2a — BIOT 670I**
