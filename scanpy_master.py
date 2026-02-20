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


# ============================
# PROMPTS
# ============================

def prompt_yes_no(question: str, default_yes: bool = False) -> bool:
    # Ask a yes or no question in the terminal
    suffix = " [Y n]: " if default_yes else " [y N]: "
    ans = input(question + suffix).strip().lower()
    if ans == "":
        return default_yes
    return ans in {"y", "yes"}


def prompt_output_dir() -> str:
    # Ask the user where results should be saved
    out = input("Enter output folder path to save results: ").strip()
    if not out:
        print("No output folder provided. Exiting.")
        sys.exit(1)
    return os.path.abspath(os.path.expanduser(out))


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
            else:
                # Most packages can be installed by their standard name via pip
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


# ============================
# MAIN PIPELINE
# ============================

def main() -> None:
    ensure_conda_env()
    check_required_modules()

    output_dir = prompt_output_dir()

    run_now = prompt_yes_no("Run the Scanpy pipeline now?", default_yes=True)
    if not run_now:
        print("Pipeline canceled by user.")
        sys.exit(0)

    run_scrublet = prompt_yes_no("Run doublet detection using scrublet?", default_yes=False)

    import scanpy as sc

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

    # ============================
    # DOUBLET DETECTION
    # ============================

    # Run scrublet if the user opts in
    if run_scrublet:
        try:
            sc.pp.scrublet(adata, batch_key="sample")
        except Exception as e:
            print("Scrublet failed, continuing without doublet calls.")
            print(str(e))

    # ============================
    # NORMALIZATION
    # ============================

    # Normalize total counts per cell
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform normalized data
    sc.pp.log1p(adata)

    # Store the full gene space for marker plots after HVG subsetting
    adata.raw = adata

    # ============================
    # HIGHLY VARIABLE GENES
    # ============================

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

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

    # ============================
    # NEIGHBORS AND UMAP
    # ============================

    # Compute neighbor graph
    sc.pp.neighbors(adata)

    # Compute UMAP
    sc.tl.umap(adata)

    # Plot UMAP colored by sample
    sc.pl.umap(adata, color="sample", size=2)

    # ============================
    # CLUSTERING
    # ============================

    # Run Leiden clustering using igraph flavor
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    # Plot clusters on UMAP
    sc.pl.umap(adata, color="leiden")

    # ============================
    # MARKER GENES
    # ============================

    # Rank marker genes per cluster
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # Plot ranked marker genes
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    # Plot PBMC marker panels using the full gene space
    pbmc_markers = {
        "T_cells": ["CD3D", "CD3E"],
        "CD4_T": ["IL7R", "CCR7"],
        "CD8_NK": ["NKG7", "GNLY"],
        "B_cells": ["MS4A1", "CD79A"],
        "Monocytes": ["LYZ", "S100A8", "S100A9"],
        "Dendritic": ["FCER1A", "CST3"],
        "Platelets": ["PPBP", "PF4"],
    }

    safe_markers = {}
    if adata.raw is not None:
        for k, genes in pbmc_markers.items():
            present = [g for g in genes if g in adata.raw.var_names]
            if present:
                safe_markers[k] = present

    if safe_markers:
        sc.pl.dotplot(adata, safe_markers, groupby="leiden", standard_scale="var", use_raw=True)

        marker_list = []
        for gene_list in safe_markers.values():
            marker_list.extend(gene_list)

        sc.pl.umap(adata, color=marker_list, ncols=3, size=10, use_raw=True)

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
