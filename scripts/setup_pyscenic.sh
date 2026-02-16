#!/bin/bash

echo "=========================================="
echo "PyScenic Setup Script"
echo "=========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install Miniconda or Anaconda first."
    exit 1
fi

# Environment name
ENV_NAME="pyscenic_env"

echo ""
echo "Step 1: Removing existing environment (if any)..."
conda env remove -n $ENV_NAME -y 2>/dev/null

echo ""
echo "Step 2: Creating new conda environment with Python 3.9..."
conda create -n $ENV_NAME python=3.9 -y

echo ""
echo "Step 3: Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

echo ""
echo "Step 4: Installing compatible NumPy..."
pip install numpy==1.23.5

echo ""
echo "Step 5: Installing PyScenic..."
pip install pyscenic

echo ""
echo "Step 6: Verifying installation..."
echo "Python version:"
python --version
echo ""
echo "NumPy version:"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"
echo ""
echo "PyScenic version:"
python -c "import pyscenic; print(f'PyScenic: {pyscenic.__version__}')"

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "To use PyScenic, activate the environment with:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Test PyScenic with:"
echo "  pyscenic --help"
echo "=========================================="
