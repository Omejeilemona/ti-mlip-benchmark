#!/bin/bash
# Download the Ti DFT dataset (Röcken & Zavadlav, npj Comput. Mater. 10, 69 (2024); Apache-2.0)
set -e
mkdir -p data && cd data
git clone --depth 1 https://github.com/tummfm/Fused-EXP-DFT-MLP.git src
mkdir -p dft && unzip -oq "src/Dataset/Data DFT and EXP/DFT_data/InitAndBulk_256atoms_curatedData.zip" -d dft
echo "DFT data in data/dft/261022_AllInitAndBulk_256atoms_with_types_curatedData/"
