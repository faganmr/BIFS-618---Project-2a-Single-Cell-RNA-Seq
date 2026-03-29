---
# figures

This folder contains all visualizations generated during the annotation stage of the pipeline. These figures are used to validate and display the cell type labels assigned to each Leiden cluster based on curated PBMC marker gene expression.

---

## Contents

| File | Description |
|------|-------------|
| `dotplot_markers.png` | Dot plot of curated PBMC marker gene expression across Leiden clusters, showing both mean expression level and fraction of cells expressing each gene |
| `umap_marker_genes.png` | UMAP panels displaying the expression of individual marker genes across all cells, arranged in a 3-column grid |
| `umap_cell_type.png` | Final UMAP plot colored by assigned cell type labels, representing the completed output of the annotation workflow |

---

## How These Figures Are Generated

All three figures are produced during the annotation stage of the pipeline using the following calls:

```python
sc.pl.dotplot(
    adata,
    safe_markers,
    groupby="leiden",
    standard_scale="var",
    use_raw=True,
    save="_markers",
)

sc.pl.umap(
    adata,
    color=marker_list,
    ncols=3,
    size=10,
    use_raw=True,
    save="_marker_genes",
)

sc.pl.umap(adata, color="cell_type", save="_cell_type")
```

All figures use expression values from `adata.raw` — the normalized, full-gene dataset stored before HVG subsetting — to ensure that marker genes not included in the HVG subset remain accessible for visualization.

---

## Interpreting the Figures

**Dot Plot** — Each row corresponds to a Leiden cluster and each column corresponds to a marker gene. Dot size reflects the percentage of cells in that cluster expressing the gene, and dot color reflects mean expression level standardized per gene. Strong, cluster-specific dot patterns indicate clean separation between cell types.

**Marker Gene UMAP Panels** — Each panel shows the expression of a single marker gene overlaid on the UMAP embedding. Marker genes for a given cell type should be concentrated in the same spatial region of the UMAP, confirming that cluster boundaries align with biological identity.

**Annotated Cell Type UMAP** — Displays the final result of the annotation workflow. Each cell is colored by its assigned cell type label. The following populations are identified in the PBMC3k dataset:

| Cell Type | Key Marker Genes |
|-----------|-----------------|
| CD4 T cells | CD3D, CD3E, IL7R, CCR7 |
| CD8 T cells | CD3D, CD8A, CD8B, GZMK |
| NK cells | NKG7, GNLY, KLRD1 |
| B cells | MS4A1, CD79A, CD79B |
| Classical Monocytes | LYZ, S100A8, S100A9 |
| Nonclassical Monocytes | FCGR3A, MS4A7 |
| Dendritic cells | FCER1A, CST3 |
| Megakaryocytes | PPBP, PF4 |

---

## References

- Scanpy Developers. (2024). *Scanpy documentation*. https://scanpy.readthedocs.io/en/stable/
- Wolf, F. A., Angerer, P., & Theis, F. J. (2018). Scanpy: Large scale single cell gene expression data analysis. *Genome Biology, 19*, Article 15. https://doi.org/10.1186/s13059-017-1382-0
```
