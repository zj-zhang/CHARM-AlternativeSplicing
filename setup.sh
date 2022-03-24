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
# download data; ~11GB
gdown 1ThJRxusyaLrH3Wof6AYPJQJL8LZJQIuX
tar -xvzf CHARM-AlternativeSplicing.data-V7.tar.gz


