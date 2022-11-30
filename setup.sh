#!/usr/bin/env bash -i
# set-up environment for reproducing
conda init bash
conda create -n charm_altsp python jupyter numpy statsmodels h5py rpy2 pandas=1.1.5 matplotlib=3.2.1 seaborn=0.11.1
conda activate charm_altsp
pip install gdown
pip install jemm==0.0.2 tqdm joblib
pip install scikit-learn==1.1.3
pip install SciencePlots==1.0.9
pip install gseapy==1.0.0
# fix legacy folder names
#ln -sf notebooks notebook
# download data; ~12GB
if [ ! -f "CHARM-AlternativeSplicing.data-V7.1.tar.gz" ]; then 
    gdown 1Pe7H8TqjLBDqCKPC245kvMo_4fKGknFM ;
else 
    echo "Found previous downloaded tarball";
fi
# unzip; ~52GB
if [ ! -d "data-V7" ]; then 
    tar -xvzf CHARM-AlternativeSplicing.data-V7.1.tar.gz ;
else
    echo "Found previous unzipped datafolder";
fi


