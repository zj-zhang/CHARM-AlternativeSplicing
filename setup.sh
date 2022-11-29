#!/usr/bin/env bash -l
set -euo pipefail

# set-up environment for reproducing
conda init bash
conda env create -n charm_altsp -f conda_rnaseq.yml
conda activate charm_altsp
pip install gdown
pip install jemm==0.0.1
pip install SciencePlots
pip install gseapy
# fix legacy folder names
ln -sf notebooks notebook
# download data; ~12GB
gdown 1Pe7H8TqjLBDqCKPC245kvMo_4fKGknFM
# unzip; ~52GB
tar -xvzf CHARM-AlternativeSplicing.data-V7.1.tar.gz


