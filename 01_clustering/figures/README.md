# figures

This folder contains all visualizations generated during the normalization, HVG selection, dimensionality reduction, and Leiden clustering stages of the pipeline.

---

## Contents

| File | Description |
|------|-------------|
| `hvg_plot.png` | Scatter plot of mean expression vs. normalized dispersion, highlighting the top 2,000 highly variable genes selected for downstream analysis |
| `pca_variance_ratio.png` | Log-scale elbow plot of variance explained by each of the top 50 principal components, used to assess how many PCs carry meaningful biological signal |
| `umap_qc_overlay.png` | Four-panel UMAP displaying Leiden cluster assignments (top left), `log1p_total_counts` (top right), `pct_counts_mt` (bottom left), and `log1p_n_genes_by_counts` (bottom right) |

---

## How These Figures Are Generated

```python
# HVG plot
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", batch_key="sample")
sc.pl.highly_variable_genes(adata)

# PCA variance ratio
sc.tl.pca(adata, random_state=594)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# UMAP + QC overlay
sc.tl.umap(adata, random_state=594)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, random_state=594)
sc.pl.umap(adata, color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"], wspace=0.5, ncols=2)
```

Figures are saved automatically to this directory because the pipeline sets `sc.settings.autosave = True` and `sc.settings.autoshow = False` at runtime.

---

## Interpreting the Figures

**HVG Plot** — Genes highlighted in the upper portion of the dispersion plot are the 2,000 selected HVGs. The remaining genes shown in gray are excluded from PCA and clustering. A clear separation between highly variable and lowly variable genes indicates a well-structured dataset.

**PCA Variance Ratio** — The steep drop after the first few principal components is expected. PCs 1 and 2 explain the largest share of variance (10.07% and 3.59% respectively), with the curve flattening beyond PC 5. This confirms that the top PCs capture the dominant axes of biological variation in the dataset.

**UMAP + QC Overlay** — The four-panel layout serves as a sanity check that cluster structure reflects true biological variation rather than technical artifacts. Ideally:
- Total counts and gene counts should be distributed broadly across clusters rather than concentrated in one region
- Mitochondrial percentage should be uniformly low across all clusters; elevated mitochondrial signal concentrated in a specific cluster may indicate residual low-quality cells

This run produced **9 Leiden clusters** (labeled 0–8), which are visually separated across the UMAP embedding.

---

## References

- McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform manifold approximation and projection for dimension reduction. *arXiv*. https://arxiv.org/abs/1802.03426
- Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: Guaranteeing well connected communities. *Scientific Reports, 9*, Article 5233. https://doi.org/10.1038/s41598-019-41695-z
- Wolf, F. A., Angerer, P., & Theis, F. J. (2018). Scanpy: Large scale single cell gene expression data analysis. *Genome Biology, 19*, Article 15. https://doi.org/10.1186/s13059-017-1382-0
