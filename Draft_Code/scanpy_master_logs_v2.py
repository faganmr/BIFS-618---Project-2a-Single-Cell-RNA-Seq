#!/usr/bin/env python3
"""
Group 2 Master Scanpy scRNA seq Pipeline
Default dataset is PBMC3k from Scanpy
"""

from __future__ import annotations

import os
import sys
import subprocess
from datetime import datetime

import anndata as ad
import scanpy as sc
import scrublet as scr
import numpy as np

import tkinter as tk
from tkinter import filedialog, messagebox


# ============================
# PROMPTS
# ============================

# Open file explorer for directory selection
root = tk.Tk()
root.withdraw()
output_dir = filedialog.askdirectory(title="Select Output Directory")
if not output_dir:
    sys.exit("No directory selected")
root.destroy()  # Clean up tkinter root

sc.settings.figdir = f"{output_dir}/figures"
sc.settings.autosave = True
sc.settings.autoshow = False

sc.settings.set_figure_params(dpi=50, facecolor="white")


# ============================
# ENVIRONMENT AND MODULE CHECKS
# ============================

def ensure_conda_env() -> None:
    # Make sure conda is active so the pipeline runs in a controlled environment
    current_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    if not current_env:
        print("Error: no conda environment detected.")
        print("Run: conda activate <your_env_name>")
        sys.exit(1)


def check_required_modules() -> None:
    # Confirm required Python modules are installed before running analysis
    required = [
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "sklearn",
        "umap",
        "igraph",
        "leidenalg",
        "scrublet",
        "skmisc",  
    ]

    missing = []
    for m in required:
        try:
            __import__(m)
        except ImportError:
            missing.append(m)

    if missing:
        print("Missing modules detected:")
        for m in missing:
            print("  -", m)

        print("")
        print("Suggested install commands:")

        # Map import names to install names
        conda_pkgs = []
        pip_pkgs = []

        for m in missing:
            if m in {"igraph", "leidenalg"}:
                conda_pkgs.extend(["python-igraph", "leidenalg"])
            elif m == "scrublet":
                pip_pkgs.append("scrublet")
            elif m == "umap":
                pip_pkgs.append("umap-learn")
            elif m == "sklearn":
                conda_pkgs.append("scikit-learn")
            elif m == "skmisc":
                pip_pkgs.append("scikit-misc")  # correct package name
            else:
                pip_pkgs.append(m)

        if conda_pkgs:
            print("")
            print("Conda:")
            print("conda install -c conda-forge " + " ".join(sorted(set(conda_pkgs))) + " -y")

        if pip_pkgs:
            print("")
            print("Pip:")
            print("pip install " + " ".join(sorted(set(pip_pkgs))))

        sys.exit(1)


def write_run_log(output_dir: str) -> None:
    # Save a run log so results are traceable during presentation
    path = os.path.join(output_dir, "run_log.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write("Group 2 Scanpy Pipeline Run Log\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write(f"Python executable: {sys.executable}\n")
        f.write(f"Python version: {sys.version}\n")
        f.write(f"Conda env: {os.environ.get('CONDA_DEFAULT_ENV', '')}\n")

    # Add an addendum to the run_log to note module versions.
    required = [
        "scanpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "sklearn",
        "umap",
        "igraph",
        "leidenalg",
        "scrublet",
        "skmisc",
    ]

    dim_file = os.path.join(output_dir, "run_log.txt")

    with open(dim_file, "a", encoding="utf-8") as f:
        f.write("\n\n=== MODULE VERSIONS ===\n")
        for r in required:
            try:
                mod = __import__(r)
                version = getattr(mod, "__version__", "unknown")
                f.write(f"{r} version: {version}\n")
            except ImportError:
                f.write(f"{r} version: not installed\n")    



# ============================
# MAIN PIPELINE
# ============================

def main() -> None:
    ensure_conda_env()
    check_required_modules()

    # Create output folders
    fig_dir = os.path.join(output_dir, "figures")
    table_dir = os.path.join(output_dir, "tables")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(table_dir, exist_ok=True)

    write_run_log(output_dir)

    # Configure Scanpy plotting behavior
    sc.settings.figdir = fig_dir
    sc.settings.autosave = True
    sc.settings.autoshow = False
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    print("Saving results to:")
    print(output_dir)

    # ============================
    # LOAD DATA
    # ============================

    # Load PBMC3k dataset
    adata = sc.datasets.pbmc3k()

    # Add a sample label for consistency with batch parameters
    adata.obs["sample"] = "pbmc3k"
    adata.obs["sample"] = adata.obs["sample"].astype("category")

    # Preserve raw counts before any transformations
    adata.layers["counts"] = adata.X.copy()

    # ============================
    # QC METRICS
    # ============================

    # Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    # Identify ribosomal genes
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

    # Identify hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains(r"^HB(?!P)", regex=True)

    # Compute QC metrics per cell
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    # Save a focused QC table
    qc_cols = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
    ]
    keep_cols = [c for c in qc_cols if c in adata.obs.columns]
    adata.obs[keep_cols].to_csv(os.path.join(table_dir, "qc_metrics.csv"), index=True)

    # QC Summary Dialog
    msg = f"""
    QC Summary:

    Cells: {adata.n_obs}
    Genes: {adata.n_vars}

    Mitochondrial %:

    Mean: {adata.obs.pct_counts_mt.mean():.2f}
    Max:  {adata.obs.pct_counts_mt.max():.2f}

    Hemoglobin %:
    Mean: {adata.obs.pct_counts_hb.mean():.2f}
    Max:  {adata.obs.pct_counts_hb.max():.2f}

    Ribosomal %:
    Mean: {adata.obs.pct_counts_ribo.mean():.2f}
    Max:  {adata.obs.pct_counts_ribo.max():.2f}

    Continue with analysis?
    """

    if not messagebox.askyesno("QC Check", msg):
        sys.exit("Analysis stopped by user")

   
   # Starting a log of adata dimensions for tracking
    dim_file = os.path.join(output_dir, "run_log.txt")

    with open(dim_file, "a", encoding="utf-8") as f:
        f.write(f"\n\n=== Matrix Dimensions ===\n")
        f.write(f"n_obs: {adata.n_obs}\n")
        f.write(f"n_vars: {adata.n_vars}\n")
        f.write(f"\n")

    # ============================
    # QC PLOTS
    # ============================

    # Plot highest expressed genes
    sc.pl.highest_expr_genes(adata, n_top=20)

    # Plot QC violin plots before filtering
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
    )

    # Plot QC scatter before filtering
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
    )

    # ============================
    # FILTERING
    # ============================

    # Filter low quality cells
    sc.pp.filter_cells(adata, min_genes=200)

    # Filter low prevalence genes
    sc.pp.filter_genes(adata, min_cells=3)

    #Tracking after Filetering
    # Starting a log of adata dimensions for tracking
    with open(dim_file, "a", encoding="utf-8") as f:
        f.write(f"Filtering:\n")
        f.write(f"n_obs: {adata.n_obs}\n")
        f.write(f"n_vars: {adata.n_vars}\n")
        f.write(f"\n")

    # ============================
    # DOUBLET DETECTION (optional - runs automatically)
    # ============================

    try:
        sc.pp.scrublet(adata, batch_key="sample", random_state=42)
    except Exception as e:
        print("Scrublet failed, continuing without doublet calls.")
        print(str(e))

    with open(dim_file, "a", encoding="utf-8") as f:
        f.write(f"Doublet Detection:\n")
        f.write(f"n_obs: {adata.n_obs}\n")
        f.write(f"n_vars: {adata.n_vars}\n")
        f.write(f"\n")
    # ============================
    # NORMALIZATION
    # ============================
    
    # Save original counts
    adata.layers["counts"] = adata.X.copy()

    # Normalize total counts per cell to 10,000
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform normalized data
    sc.pp.log1p(adata)

    # Save normalized full dataset
    adata.raw = adata

    # Normalization Summary Prompt
    msg = f"""
    Normalization Summary:
    Cells: {adata.n_obs}
    Genes: {adata.n_vars}
    Continue with analysis?
    """

    if not messagebox.askyesno("Normalization", msg):
        sys.exit("Analysis stopped by user")

    print("Normalization Complete")

    with open(dim_file, "a", encoding="utf-8") as f:
        f.write(f"Normalization:\n")
        f.write(f"n_obs: {adata.n_obs}\n")
        f.write(f"n_vars: {adata.n_vars}\n")
        f.write(f"\n")
    # ============================
    # HIGHLY VARIABLE GENES
    # ============================

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")

    # Plot HVG selection
    sc.pl.highly_variable_genes(adata)

    # Subset to HVGs for PCA and clustering
    if "highly_variable" in adata.var.columns:
        adata = adata[:, adata.var["highly_variable"]].copy()

        

    # ============================
    # PCA
    # ============================

    # Run PCA
    sc.tl.pca(adata)

    # Plot PCA variance ratio
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

    with open(dim_file, "a", encoding="utf-8") as f:
        f.write(f"PCA:\n")
        f.write(f"n_obs: {adata.n_obs}\n")
        f.write(f"n_vars: {adata.n_vars}\n")
        f.write(f"\n")

    # ============================
    # NEIGHBORS AND UMAP
    # ============================

    # Compute neighbor graph
    sc.pp.neighbors(adata)

    # Compute UMAP
    sc.tl.umap(adata)

    # Plot UMAP colored by sample
    sc.pl.umap(adata, color="sample", size=2)

    with open(dim_file, "a", encoding="utf-8") as f:
        f.write(f"UMAP:\n")
        f.write(f"n_obs: {adata.n_obs}\n")
        f.write(f"n_vars: {adata.n_vars}\n")

        f.write(f"\n")

    # ============================
    # CLUSTERING
    # ============================

    # Run Leiden clustering using igraph flavor
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    # Plot clusters on UMAP
    sc.pl.umap(adata, color="leiden")

    with open(dim_file, "a", encoding="utf-8") as f:
        f.write(f"Clustering:\n")
        f.write(f"n_obs: {adata.n_obs}\n")
        f.write(f"n_vars: {adata.n_vars}\n")

    # ============================
    # MARKER GENES
    # ============================

    # Rank marker genes per cluster
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # Plot ranked marker genes
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    # Plot PBMC marker panels using the full gene space
    pbmc_markers = {
    "CD4_T_cells": ["CD3D", "CD3E", "IL7R", "CCR7"],
    "CD8_T_cells": ["CD3D", "CD8A", "CD8B", "GZMK"],
    "NK_cells": ["NKG7", "GNLY", "KLRD1"],
    "B_cells": ["MS4A1", "CD79A", "CD79B"],
    "Classical_Monocytes": ["LYZ", "S100A8", "S100A9"],
    "Nonclassical_Monocytes": ["FCGR3A", "MS4A7"],
    "Dendritic_cells": ["FCER1A", "CST3"],
    "Megakaryocytes": ["PPBP", "PF4"],
}

    safe_markers = {}
    if adata.raw is not None:
        for k, genes in pbmc_markers.items():
            present = [g for g in genes if g in adata.raw.var_names]
            if present:
                safe_markers[k] = present

    # Plot dotplot for marker genes across Leiden clusters
    if safe_markers:
        sc.pl.dotplot(
            adata,
            safe_markers,
            groupby="leiden",
            standard_scale="var",
            use_raw=True,
            save="_markers",
        )

        # Plot marker genes on UMAP for quick visual confirmation
        marker_list = []
        for gene_list in safe_markers.values():
            marker_list.extend(gene_list)

        sc.pl.umap(
            adata,
            color=marker_list,
            ncols=3,
            size=10,
            use_raw=True,
            save="_marker_genes",
        )

    # ============================
    # CELL TYPE ANNOTATION AND UMAP
    # ============================

    # Assign a cell type label to each Leiden cluster using marker gene expression
    def mean_expr_raw(mask, genes):
        x = adata.raw[mask, genes].X
        try:
            return float(x.mean())
        except Exception:
            return float(x.toarray().mean())

    cluster_labels = {}
    clusters = sorted(adata.obs["leiden"].unique().tolist())

    for cl in clusters:
        mask = adata.obs["leiden"] == cl
        best_label = "Unknown"
        best_score = -1.0

        for label, genes in safe_markers.items():
            score = mean_expr_raw(mask, genes)
            if score > best_score:
                best_score = score
                best_label = label

        cluster_labels[cl] = best_label

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_labels).astype("category")

    # Plot UMAP colored by annotated cell types
    sc.pl.umap(adata, color="cell_type", save="_cell_type")

    #=============================
    #Reproducibility Check 
    #=============================



    repro_file = os.path.join(output_dir, "run_log.txt")

    with open(repro_file, "a", encoding="utf-8") as f:

        f.write("\n=== REPRODUCIBILITY LOG ===\n")

        # Basic dimensions
        f.write(f"Cells (n_obs): {adata.n_obs}\n")
        f.write(f"Genes (n_vars): {adata.n_vars}\n")

        # HVG count
        if "highly_variable" in adata.var.columns:
            f.write(f"HVG count: {adata.var['highly_variable'].sum()}\n")

        # Cluster info
        if "leiden" in adata.obs.columns:
            f.write(f"Number of clusters: {adata.obs['leiden'].nunique()}\n")

        # PCA variance
        if "pca" in adata.uns:
            var_ratio = adata.uns["pca"]["variance_ratio"][:5]
            f.write(f"PCA variance ratio (first 5 PCs): {var_ratio.tolist()}\n")

        # Neighbor graph info
        if "connectivities" in adata.obsp:
            f.write(f"Neighbor graph shape: {adata.obsp['connectivities'].shape}\n")
            f.write(f"Total neighbor edges: {adata.obsp['connectivities'].nnz}\n")

        # Doublet info
        if "predicted_doublet" in adata.obs.columns:
            f.write(f"Predicted doublets: {adata.obs['predicted_doublet'].sum()}\n")

        # Random seed snapshot (numpy)
        f.write(f"NumPy seed snapshot: {np.random.get_state()[1][0]}\n")

        f.write("=== END LOG ===\n\n")

        print("Reproducibility log written to:", repro_file)

        # ============================
        # FINAL OVERLAYS
        # ============================

        # Overlay doublet calls if scrublet was run
        overlay_cols = ["leiden"]
        for c in ["predicted_doublet", "doublet_score"]:
            if c in adata.obs.columns:
                overlay_cols.append(c)

        sc.pl.umap(adata, color=overlay_cols, wspace=0.5, size=3)

        # Overlay QC metrics for a final sanity check
        qc_overlay = []
        for c in ["log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"]:
            if c in adata.obs.columns:
                qc_overlay.append(c)

        if qc_overlay:
            sc.pl.umap(adata, color=["leiden"] + qc_overlay, wspace=0.5, ncols=2)

        # ============================
        # SAVE OUTPUTS
        # ============================

        # Save processed object
        out_h5ad = os.path.join(output_dir, "pbmc3k_processed.h5ad")
        adata.write(out_h5ad)

        print("Pipeline finished successfully.")
        print("Saved:")
        print(out_h5ad)

    # Open figures directory on Linux if possible
    if sys.platform.startswith("linux"):
        try:
            subprocess.run(["xdg-open", fig_dir], check=False)
        except Exception:
            pass


if __name__ == "__main__":
    main()
