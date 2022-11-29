# CHARM-AlternativeSplicing
## Developing Alternative Splicing Diagnostic Biomarkers using COVID-19 Health Action Response for Marines

A collection of jupyter notebooks, for Blood Transcriptome Analysis of Alternative Splicing in CHARM Covid-19 datasets

## Installation
First clone this github repository to your local machine.

Run the provided `setup.sh` bash script to automate the installation process. 
This will create a new conda environment named `charm_altsp` with all packages specified in `conda_rnaseq.yml` file. 

The script will also automatically download the data files, which requires ~50GB of local storage space. `gdown` may fail to download the file due to network issues, in which case you may either upgrade gdown by `pip install --upgrade gdown`, or follow the "You may still be able to access the file from the browser" link and download from your browser (not recommended as the file is large).

The inclusion of all intermediate files is to faciliate re-run and reproducibility of our results without starting from scratch (e.g., re-aligning all FASTQ reads). Therefore, one can run any individual Jupyter notebook in a matter a few minutes (depending on your local disk I/O speed).


## Reproduce the Whole-blood RNA-seq Biomarker Discovery

To reproduce CHARM DAS/DEG predictive model, run `notebooks/05-Navy_predictive-model.V2.ipynb`. This also generates all comparisons between CHARM DEG/DAS and six public available gene signatures.

If you'd like to extract the trained model based on CHARM DAS/DEG, the trained classifiers are stored in a pickled Python dictionary in file `data-V7/das_classifier2/predictive_model_dict.pkl`. To read in, use the following Python code:
```python
(clfs, model_eval_df, train_stats, train_psi_dfs, train_dfs, duke_psi_dfs, duke_preds, duke_y) = pickle.load(
    open('%s/predictive_model_dict.pkl' % CLFDIR, 'rb'))
print('reloaded previously trained model')
```
The usage of the above code snippet can also be found in the notebook `notebooks/05-Navy_predictive-model.V2.ipynb`.


## Reproduce the Fluidigm Biomarker Validation Results

To reproduce the Fluidigm diagnostic AS biomarker validation performance using our preset AS biomarkers, run `notebooks/06-Fluidigm.V2.ipynb`.

To reproduce the forward selection of Fluidigm AS biomarkers, run `notebooks/06-Fluidigm.FZZgen.ipynb`.

This notebook will require the data files in `data-V7/fluidigm/`.


## Reproduce the Statistical Testing Results

Our statistical inference is based on a Python Package `jemm`, which should have been installed if you followed the `setup.sh` script.

By default, the pre-generated regression tables are stored in `data-V7/joint_[SE|A5|A3|RI]` folders for each AS type. This should faciliate most downstream analyses unless you want to re-generate these files from scratch.

There are two ways to reproduce the Statistical Testing Results:
- You can run each notebook from `01-Navy_joint-SE.ipynb` to `04-Navy_joint-RI.ipynb`. The alternative splicing type is as indicated in the filename. You will need to set the flag `FORCE_RERUN=True`
- Use `run_jem_tests.py` script. This script can be called by `sbatch_jem.sh` to submit slurm jobs to compute nodes.

## Reproduce All Other Results

Follow the notebooks. The code should be straightforward (though messy, as I don't have enough time to clean up everything) and runs fast with all the intermediate files available.

## Contact

Still have questions? Leave a GitHub Issue and I will get back to you.