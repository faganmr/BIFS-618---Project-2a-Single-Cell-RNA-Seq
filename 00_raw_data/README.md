---
# 00_raw_data

This folder contains the raw input data, quality control (QC) metrics, and preprocessing outputs generated at the start of the pipeline before any filtering or normalization is applied.

---

## Contents

| File | Description |
|------|-------------|
| `qc_metrics.csv` | Per-cell QC metrics table exported after `sc.pp.calculate_qc_metrics()`, including `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo`, and `pct_counts_hb` |
| `figures/highest_expr_genes.png` | Bar plot of the top 20 most highly expressed genes before filtering |
| `figures/violin.png` | Violin plots of gene counts, total counts, and mitochondrial percentage |
| `figures/scatter.png` | Scatter plot of total counts vs. gene counts, colored by mitochondrial percentage |

---

## Dataset

The pipeline loads the **PBMC3k dataset** directly via Scanpy using:

```python
adata = sc.datasets.pbmc3k()
```

This dataset contains approximately 2,700 peripheral blood mononuclear cells (PBMCs) sequenced using 10x Genomics Chromium v1 chemistry and is a standard benchmark dataset for scRNA-seq analysis pipelines (10x Genomics, 2017).

---

## Data Sources

### Primary Access Method (Used in This Project)

The dataset is loaded programmatically via Scanpy:

https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html

### Original Experimental Source

The underlying experimental data was generated and released by **10x Genomics**:

https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0

---

## QC Metrics Computed

The following metrics are calculated per cell using `sc.pp.calculate_qc_metrics()`:

- **`n_genes_by_counts`** — Number of genes detected per cell
- **`total_counts`** — Total UMI counts per cell
- **`pct_counts_mt`** — Percentage of counts from mitochondrial genes (genes prefixed `MT-`)
- **`pct_counts_ribo`** — Percentage of counts from ribosomal genes (prefixed `RPS` or `RPL`)
- **`pct_counts_hb`** — Percentage of counts from hemoglobin genes

---

## Filtering Thresholds Applied

After QC visualization, the following thresholds are applied to remove low-quality observations:

```python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| Minimum genes per cell | 200 | Removes empty droplets and severely damaged cells |
| Minimum cells per gene | 3 | Removes extremely rare transcripts that add noise without biological signal |

---

## Doublet Detection

Doublet detection is performed using **Scrublet** immediately after filtering:

```python
sc.pp.scrublet(adata, batch_key="sample", random_state=594)
```

Scrublet generates simulated artificial doublets and assigns each observed cell a `doublet_score` and a `predicted_doublet` flag. This step is wrapped in a `try/except` block to ensure pipeline stability across environments where Scrublet may be unavailable.

---

## Matrix Dimensions (Logged)

The run log records matrix dimensions at this stage. Expected values for PBMC3k:

| Stage | n_obs (cells) | n_vars (genes) |
|-------|--------------|----------------|
| Raw load | 2,700 | 32,738 |
| After filtering | 2,700 | 13,714 |

---

## References

- 10x Genomics. (2017). *3k peripheral blood mononuclear cells (PBMCs) from a healthy donor (v1 chemistry)*. https://www.10xgenomics.com/resources/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-0-0
- Wolf, F. A., Angerer, P., & Theis, F. J. (2018). Scanpy: Large scale single cell gene expression data analysis. *Genome Biology, 19*, Article 15. https://doi.org/10.1186/s13059-017-1382-0
