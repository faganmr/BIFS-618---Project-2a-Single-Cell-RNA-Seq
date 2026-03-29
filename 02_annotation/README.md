# 02_annotation

This folder contains outputs from the cell type annotation stage of the pipeline. Annotation assigns biological labels to each Leiden cluster based on the expression of known PBMC marker genes. This step does not create new clusters — it maps numerical cluster identities to biologically meaningful cell type names.

---

## Contents

| File | Description |
|------|-------------|
| `figures/umap_cell_type.png` | Final UMAP plot colored by assigned cell type labels |
| `figures/dotplot_markers.png` | Dot plot showing marker gene expression level and detection rate across Leiden clusters |
| `figures/umap_marker_genes.png` | UMAP panels displaying expression of individual marker genes across cells |

---

## Annotation Strategy

Annotation in this pipeline is performed at the **cluster level** using a manually curated marker gene dictionary tailored to the expected immune cell populations in the PBMC3k dataset. The pipeline does not use supervised classification; instead, it computes mean marker gene expression per cluster and assigns the label with the highest average expression score.

This approach uses expression values from `adata.raw` — the normalized, full-gene dataset stored before HVG subsetting — to ensure that marker genes not included in the HVG subset are still accessible for scoring.

---

## Step 1 — Marker Gene Ranking

Differentially expressed genes are ranked per cluster using the Wilcoxon rank-sum test:

```python
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

This step provides an unsupervised view of which genes most distinguish each cluster before formal labels are applied.

---

## Step 2 — PBMC Marker Dictionary

The following marker panels are defined directly in the pipeline script and used for cluster scoring:

| Cell Type | Marker Genes |
|-----------|-------------|
| CD4 T cells | CD3D, CD3E, IL7R, CCR7 |
| CD8 T cells | CD3D, CD8A, CD8B, GZMK |
| NK cells | NKG7, GNLY, KLRD1 |
| B cells | MS4A1, CD79A, CD79B |
| Classical Monocytes | LYZ, S100A8, S100A9 |
| Nonclassical Monocytes | FCGR3A, MS4A7 |
| Dendritic cells | FCER1A, CST3 |
| Megakaryocytes | PPBP, PF4 |

Only marker genes confirmed to be present in `adata.raw.var_names` are used for scoring. This safe marker filtering step prevents runtime errors if any expected genes are absent from the dataset:

```python
safe_markers = {}
if adata.raw is not None:
    for k, genes in pbmc_markers.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if present:
            safe_markers[k] = present
```

---

## Step 3 — Cluster Scoring and Label Assignment

Each Leiden cluster is scored against every marker panel by computing the mean expression of that panel's genes across all cells in the cluster, using raw normalized values:

```python
def mean_expr_raw(mask, genes):
    x = adata.raw[mask, genes].X
    try:
        return float(x.mean())
    except Exception:
        return float(x.toarray().mean())
```

The cell type label with the highest mean expression score is assigned to each cluster. Clusters that do not exceed a positive score for any panel receive the label `"Unknown"` as a fallback.

Label assignments are then mapped back to individual cells:

```python
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_labels).astype("category")
```

---

## Step 4 — Annotation Visualization

Two visualizations are generated to support and validate the annotation:

**Dot plot** — shows both the mean expression level and the fraction of cells expressing each marker gene across clusters:

```python
sc.pl.dotplot(adata, safe_markers, groupby="leiden", standard_scale="var", use_raw=True, save="_markers")
```

**Marker gene UMAP panels** — overlays individual marker gene expression on the UMAP to confirm that marker-positive cells co-localize with their expected clusters:

```python
sc.pl.umap(adata, color=marker_list, ncols=3, size=10, use_raw=True, save="_marker_genes")
```

**Final annotated UMAP** — displays cell type labels directly on the UMAP:

```python
sc.pl.umap(adata, color="cell_type", save="_cell_type")
```

---

## Cell Types Identified

The pipeline identified the following cell populations in the PBMC3k dataset:

- B cells
- CD4 T cells
- CD8 T cells
- Classical Monocytes
- Dendritic cells
- Megakaryocytes
- NK cells
- Nonclassical Monocytes

---

## Notes on Generalizability

The current annotation approach is tailored to PBMC cell populations. Applying this pipeline to other datasets would require replacing the `pbmc_markers` dictionary with panels appropriate to the tissue or condition under study. Automated annotation tools such as CellTypist were explored during development and represent a promising direction for future pipeline extensions.

---

## References

- Luecken, M. D., & Theis, F. J. (2019). Current best practices in single cell RNA seq analysis: A tutorial. *Molecular Systems Biology, 15*(6), e8746. https://doi.org/10.15252/msb.20188746
- Scanpy Developers. (2024). *Scanpy documentation*. https://scanpy.readthedocs.io/en/stable/
- Wolf, F. A., Angerer, P., & Theis, F. J. (2018). Scanpy: Large scale single cell gene expression data analysis. *Genome Biology, 19*, Article 15. https://doi.org/10.1186/s13059-017-1382-0
