#!/bin/bash

echo "Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sc

echo ""
echo "Python being used:"
which python
python --version
echo ""

echo "Checking Python modules..."

python - <<EOF
import sys

modules = [
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "sklearn",
    "umap",
    "igraph",
    "scrublet"
]

missing = []

for m in modules:
    try:
        __import__(m)
        print(f"OK: {m}")
    except ImportError:
        print(f"MISSING: {m}")
        missing.append(m)

if missing:
    print("\nMissing modules detected:", ", ".join(missing))
    sys.exit(1)
else:
    print("\nAll required modules are installed.")
EOF

STATUS=$?

if [ $STATUS -ne 0 ]; then
    echo ""
    echo "Environment check failed. Pipeline will NOT run."
    exit 1
fi

# ==============================
# Confirmation Box
# ==============================

echo ""
echo "========================================"
echo " Environment check PASSED"
echo " Conda environment: sc"
echo " Python pipeline ready to run"
echo "========================================"
echo ""

read -p "Do you want to run the Scanpy pipeline now? (y/n): " CONFIRM

case "$CONFIRM" in
    y|Y|yes|YES)
        echo ""
        echo "Starting Scanpy pipeline..."
        ;;
    *)
        echo ""
        echo "Pipeline execution canceled by user."
        exit 0
        ;;
esac

echo ""
python /home/tyler/Desktop/10XGenome_Scanpy.py
