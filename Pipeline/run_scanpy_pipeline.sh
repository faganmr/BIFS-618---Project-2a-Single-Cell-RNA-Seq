#!/bin/bash

ENV_NAME="$1"
SCRIPT_NAME="${2:-scanpy_master_logs_seed_v3.py}"
PYTHON_VERSION="3.10.19"

if [ -z "$ENV_NAME" ]; then
    echo "ERROR: Please provide a Conda environment name."
    echo "Usage:"
    echo "bash run_scanpy_pipeline.sh <env_name> [script_name]"
    exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: Conda was not found."
    exit 1
fi

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    echo "Environment already exists: $ENV_NAME"
else
    echo "Environment not found. Creating $ENV_NAME with Python $PYTHON_VERSION..."
    conda create -y -n "$ENV_NAME" python="$PYTHON_VERSION"

    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to create environment $ENV_NAME"
        exit 1
    fi
fi

conda activate "$ENV_NAME"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate environment $ENV_NAME"
    exit 1
fi

CURRENT_PYTHON="$(python -c 'import sys; print(sys.version.split()[0])')"

echo "Active Python version: $CURRENT_PYTHON"

if [ "$CURRENT_PYTHON" != "$PYTHON_VERSION" ]; then
    echo "ERROR: Python version mismatch."
    echo "Expected: $PYTHON_VERSION"
    echo "Found: $CURRENT_PYTHON"
    exit 1
fi

if [ ! -f "$SCRIPT_NAME" ]; then
    echo "ERROR: Script not found: $SCRIPT_NAME"
    exit 1
fi

echo "Running pipeline script: $SCRIPT_NAME"
python "$SCRIPT_NAME"
