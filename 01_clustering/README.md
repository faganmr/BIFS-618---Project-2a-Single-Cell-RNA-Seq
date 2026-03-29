---
# 01_clustering

This folder contains outputs from the normalization, highly variable gene (HVG) selection, dimensionality reduction, and Leiden clustering stages of the pipeline. These steps transform the filtered expression matrix into a structured set of cell clusters that reflect underlying transcriptional similarity.

---

## Contents

| File | Description |
|------|-------------|
| `figures/hvg_plot.png` | Scatter plot of mean expression vs. normalized dispersion, highlighting the top 2,000 HVGs selected for downstream analysis |
| `figures/pca_variance_ratio.png` | Elbow plot of variance explained by each principal component, used to assess the informativeness of the PCA decomposition |
| `figures/umap_qc_overlay.png` | Four-panel UMAP plot showing Leiden cluster assignments alongside `log1p_total_counts`, `pct_counts_mt`, and `log1p_n_genes_by_counts` overlays |

---

## Workflow Overview

```
Raw counts (filtered)
        │
        ▼
Normalization (library-size + log1p)
        │
        ▼
HVG Selection (top 2,000 genes, Seurat flavor)
        │
        ▼
PCA (dimensionality reduction)
        │
        ▼
k-Nearest Neighbor Graph
        │
        ▼
UMAP (visualization)
        │
        ▼
Leiden Clustering (community detection)
```

---

## Step 1 — Normalization

Raw counts are preserved before any transformation:

```python
adata.layers["counts"] = adata.X.copy()
```

Each cell is then normalized to a total count of 10,000 molecules, followed by log transformation:

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
```

Storing the normalized matrix in `adata.raw` ensures that downstream marker gene scoring and annotation reference a consistent, full-gene expression state even after the dataset is later subset to HVGs.

---

## Step 2 — Highly Variable Gene (HVG) Selection

The top 2,000 most variable genes are identified using the Seurat dispersion-based method:

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", batch_key="sample")
adata = adata[:, adata.var["highly_variable"]].copy()
```

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `n_top_genes` | 2,000 | Captures major immune lineage markers while reducing noise |
| `flavor` | `"seurat"` | Bins genes and z-scores dispersion for normalized comparison across expression levels |
| `batch_key` | `"sample"` | Computes variability within samples before comparing across genes; retained for future multi-sample compatibility |

---

## Step 3 — Principal Component Analysis (PCA)

PCA is run on the HVG-subset matrix to reduce dimensionality while preserving dominant axes of variation:

```python
sc.tl.pca(adata, random_state=594)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```

Variance explained by the first 5 principal components from this run:

| PC | Variance Ratio |
|----|---------------|
| 1 | 0.1007 |
| 2 | 0.0359 |
| 3 | 0.0241 |
| 4 | 0.0112 |
| 5 | 0.0070 |

---

## Step 4 — Neighbor Graph Construction

A k-nearest neighbor (kNN) graph is constructed in PCA space, where each cell is a node and edges connect cells with similar transcriptional profiles:

```python
sc.pp.neighbors(adata, random_state=594)
```

| Metric | Value |
|--------|-------|
| Neighbor graph shape | (2700, 2700) |
| Total neighbor edges | 62,096 |

---

## Step 5 — UMAP

UMAP provides a two-dimensional visualization of the structure encoded in the neighbor graph. It does not determine clusters itself but reflects the relationships already defined by the kNN graph:

```python
sc.tl.umap(adata, random_state=594)
```

---

## Step 6 — Leiden Clustering

Community detection is performed on the neighbor graph using the Leiden algorithm with the igraph backend:

```python
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, random_state=594)
sc.pl.umap(adata, color="leiden")
```

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `flavor` | `"igraph"` | Preferred backend per current Scanpy documentation |
| `n_iterations` | `2` | Sufficient for stable community detection with minimal overhead |
| `random_state` | `594` | Global seed for reproducibility |

This run produced **9 Leiden clusters** (labeled 0–8), with cluster assignments stored in `adata.obs["leiden"]` until cell type annotation is applied in `02_annotation`.

---

## Matrix Dimensions

All values confirmed from `run_log.txt` (timestamp: 2026-03-25T19:41:44):

| Stage | n_obs (cells) | n_vars (genes) |
|-------|--------------|----------------|
| Raw load | 2,700 | 32,738 |
| After filtering | 2,700 | 13,714 |
| After normalization | 2,700 | 13,714 |
| After HVG subset | 2,700 | 2,000 |
| After PCA | 2,700 | 2,000 |
| After UMAP | 2,700 | 2,000 |
| After clustering | 2,700 | 2,000 |

---

## Reproducibility

All stochastic steps in this module use `SEED = 594`:

| Step | Seed applied |
|------|-------------|
| PCA | `random_state=594` |
| Neighbor graph | `random_state=594` |
| UMAP | `random_state=594` |
| Leiden clustering | `random_state=594` |

Additional reproducibility metrics from `run_log.txt`:

| Metric | Value |
|--------|-------|
| HVG count | 2,000 |
| Number of clusters | 9 |
| Predicted doublets | 38 |
| NumPy seed snapshot | 2857075577 |

---

## References

- McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform manifold approximation and projection for dimension reduction. *arXiv*. https://arxiv.org/abs/1802.03426
- Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: Guaranteeing well connected communities. *Scientific Reports, 9*, Article 5233. https://doi.org/10.1038/s41598-019-41695-z
- Wolf, F. A., Angerer, P., & Theis, F. J. (2018). Scanpy: Large scale single cell gene expression data analysis. *Genome Biology, 19*, Article 15. https://doi.org/10.1186/s13059-017-1382-0
```
