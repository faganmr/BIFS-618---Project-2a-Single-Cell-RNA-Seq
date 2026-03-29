---
# figures

This folder contains all quality control (QC) visualizations generated during the initial data loading and preprocessing stage of the pipeline, before any filtering or normalization is applied.

---

## Contents

| File | Description |
|------|-------------|
| `highest_expr_genes.png` | Bar plot of the top 20 most highly expressed genes across all cells before filtering, used to identify dominant transcripts and rule out technical artifacts |
| `violin.png` | Three-panel violin plot displaying the distribution of `n_genes_by_counts`, `total_counts`, and `pct_counts_mt` across all cells, used to assess overall cell quality and guide filtering thresholds |
| `scatter.png` | Scatter plot of total UMI counts vs. number of detected genes per cell, colored by mitochondrial percentage, used to identify low-quality cells and outliers prior to filtering |

---

## How These Figures Are Generated

All three figures are produced during the QC stage of the pipeline using the following calls:

```python
sc.pl.highest_expr_genes(adata, n_top=20)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(
    adata,
    x="total_counts",
    y="n_genes_by_counts",
    color="pct_counts_mt",
)
```

Figures are saved automatically to this directory because the pipeline sets `sc.settings.autosave = True` and `sc.settings.autoshow = False` at runtime.

---

## Interpreting the Figures

**Highest Expressed Genes** — Genes like `MALAT1`, `B2M`, and ribosomal genes (`RPL`, `RPS` prefixes) dominating the top of this plot are expected in PBMC data and indicate normal biological signal rather than technical artifacts.

**Violin Plots** — The bulk of cells should show moderate, normally distributed gene and count values. Cells in the low tail of `n_genes_by_counts` or the high tail of `pct_counts_mt` are candidates for removal in the filtering step.

**Scatter Plot** — Cells should follow a roughly linear relationship between total counts and gene counts. Cells that deviate significantly from this trend, particularly those with high mitochondrial percentage, are likely low quality.

---

## References

- Wolf, F. A., Angerer, P., & Theis, F. J. (2018). Scanpy: Large scale single cell gene expression data analysis. *Genome Biology, 19*, Article 15. https://doi.org/10.1186/s13059-017-1382-0
