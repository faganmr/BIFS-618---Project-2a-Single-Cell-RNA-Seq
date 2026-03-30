# Code — Master Pipeline

This folder contains the two scripts that make up the pipeline: a Bash launcher that handles environment setup, and a Python script that runs the full scRNA-seq analysis.

---

## Files

| File | Description |
|---|---|
| `run_scanpy_pipeline.sh` | Bash script — creates/activates the conda environment and launches the pipeline |
| `scanpy_master_logs_seed_v3.py` | Python script — full scRNA-seq analysis pipeline |

---

## Usage

```bash
bash run_scanpy_pipeline.sh <env_name> [script_name]
```

**Example:**
```bash
bash run_scanpy_pipeline.sh scanpy scanpy_master_logs_seed_v3.py
```

The `env_name` argument is required. The `script_name` argument is optional and defaults to `scanpy_master_logs_seed_v3.py`.

---

## Bash Script — `run_scanpy_pipeline.sh`

The bash script handles environment setup before handing off to Python. It performs the following steps in order:

1. **Validates inputs** — exits with a usage message if no environment name is provided
2. **Checks for conda** — exits if conda is not found on the system
3. **Creates or reuses the conda environment** — if the named environment does not exist, it creates one with Python 3.10.19; if it already exists, it skips creation
4. **Activates the environment** — sources the conda profile and activates the named environment
5. **Validates the Python version** — confirms the active Python matches the expected version (3.10.19) and exits if there is a mismatch
6. **Installs base requirements** — runs `pip install numpy` to satisfy the import at the top of the Python script before the self-installer runs
7. **Launches the pipeline** — runs the Python script with the active environment's interpreter

---

## Python Pipeline — `scanpy_master_logs_seed_v3.py`

### Environment Check

On startup the script verifies that a conda environment is active, then checks all required packages against expected versions. If any packages are missing or mismatched it presents a GUI dialog asking whether to install or update them automatically. After installation the script restarts itself so the newly installed packages are loaded cleanly.

Required packages and expected versions:

| Package | Version |
|---|---|
| scanpy | 1.11.5 |
| anndata | 0.11.4 |
| numpy | 2.2.6 |
| pandas | 2.3.3 |
| scipy | 1.15.2 |
| matplotlib | 3.10.8 |
| seaborn | 0.13.2 |
| scikit-learn | 1.7.2 |
| umap-learn | 0.5.11 |
| python-igraph | 0.10.8 |
| leidenalg | 0.10.2 |
| scikit-misc | 0.5.2 |

### Output Directory

A file dialog opens at startup prompting the user to select an output directory. All figures, tables, and logs are written there.

### Analysis Steps

#### 1. Load Data
The PBMC 3k dataset is loaded via `sc.datasets.pbmc3k()`. Raw counts are preserved in `adata.layers["counts"]` before any transformations.

#### 2. Quality Control
Mitochondrial, ribosomal, and hemoglobin genes are flagged. Per-cell QC metrics are computed with `sc.pp.calculate_qc_metrics()` and saved to `tables/qc_metrics.csv`. QC plots are generated and a summary dialog allows the user to review before continuing.

#### 3. Filtering
Low-quality cells (fewer than 200 detected genes) and low-prevalence genes (detected in fewer than 3 cells) are removed.

#### 4. Doublet Detection
Scrublet is run via `sc.pp.scrublet()` to flag likely doublets. This step is wrapped in a try/except and the pipeline continues if it fails.

#### 5. Normalization
Total counts per cell are normalized to 10,000 and log-transformed with `sc.pp.log1p()`. The normalized full-gene object is stored in `adata.raw` for use in downstream marker gene plots.

#### 6. Highly Variable Gene Selection
The top 2,000 highly variable genes are identified using the Seurat-flavor dispersion method (`sc.pp.highly_variable_genes()`). The matrix is then subset to these genes for PCA and clustering.

#### 7. PCA
PCA is run with `sc.tl.pca()`. A variance ratio plot is saved to assess how many principal components capture meaningful variance.

#### 8. Neighbor Graph
A k-nearest-neighbor graph is constructed in PC space with `sc.pp.neighbors()`. This graph is shared by both UMAP and Leiden clustering downstream.

#### 9. UMAP
UMAP coordinates are computed with `sc.tl.umap()` and plotted coloured by sample.

#### 10. Leiden Clustering
Leiden clustering is run with `sc.tl.leiden()` using the igraph flavor. Clusters are plotted on the UMAP.

#### 11. Marker Genes
Wilcoxon rank-sum marker genes are computed per cluster with `sc.tl.rank_genes_groups()`. A curated set of PBMC marker genes is used to generate a dot plot and per-marker UMAP overlays.

#### 12. Cell Type Annotation
Each Leiden cluster is assigned a cell type label based on mean marker gene expression from the curated panel. Labels are stored in `adata.obs["cell_type"]` and plotted on the UMAP.

### Outputs

| File | Description |
|---|---|
| `run_log.txt` | Environment info, seed, package versions, matrix dimensions, reproducibility metrics |
| `tables/qc_metrics.csv` | Per-cell QC values |
| `figures/` | All plots saved automatically by Scanpy |
| `pbmc3k_processed.h5ad` | Final processed AnnData object |

### Seed

All stochastic steps use `SEED = 594`:

```python
np.random.seed(SEED)
random.seed(SEED)
sc.tl.pca(adata, random_state=SEED)
sc.pp.neighbors(adata, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
sc.tl.leiden(adata, random_state=SEED)
sc.pp.scrublet(adata, random_state=SEED)
```
