"""cli for running LMM tests

FZZ, 10/20/2021
"""


from jemm.junction import JunctionCountTable
from jemm.transcript import TranscriptMeasureTable
from jemm.data_table import DataTableIndexed
from jemm.model import JemmLinearRegression, JemmLinearMixedModel
from jemm.plots import facet_boxplot
from jemm.covariate import Contrasts, Covariate
from jemm.utils import setup_logger
from jemm.genomic_annotation import ExonSet

import os
import pickle
import joblib


import argparse


PCS_TO_INCL = ['PC%i'%i for i in [0, 9]]
LOGGER = setup_logger()

USE_RE = True
MIN_RE_VAR = 0.001

USE_ANNOT = True


def get_covs(DATA_VER):
    # build constrasts
    contrasts = Contrasts(name="final", levels=[
        'Control',
        'Pre',
        'First',
        'Mid',
        'Post',
    ])
    covs = Covariate(fp="./%s/charm_master.clean_w_10pc.csv" % DATA_VER, sep=",",
                     index_col=0,
                     contrasts=contrasts,
                     main_effects=['final', 'pid','Sex', 'plateNum'] + PCS_TO_INCL if USE_RE is True else ['final', 'Sex', 'plateNum'] + PCS_TO_INCL,
                     factor_conversion={
                         'final': {
                             'First': 'final@First',
                             'Mid': 'final@Mid',
                             'Post': 'final@Post',
                             'Pre': 'final@Pre',
                         },
                         'plateNum': {
                             'P2R': 'plate@P2R',
                             'P3': 'plate@P3',
                             'P4': 'plate@P4',
                             'P5': 'plate@P5',
                             'P6': 'plate@P6',
                             'P7': 'plate@P7',
                             'P8': 'plate@P8',
                             'P9': 'plate@P9',
                             'P10': 'plate@P10',
                             'P12': 'plate@P12',
                             'P13': 'plate@P13',
                             'P14': 'plate@P14',
                             'P15': 'plate@P15',
                             'P16': 'plate@P16',
                             'P17': 'plate@P17',
                             'P18': 'plate@P18',
                             'P19': 'plate@P19',
                             'P20': 'plate@P20',
                             'P21': 'plate@P21',
                             'P22': 'plate@P22',
                             'P23': 'plate@P23',
                             'P24': 'plate@P24',
                             'P25': 'plate@P25'
                         },
                         'Sex': {'M': 'Sex@M'}
                     },
                     verbose=True
                 )
    return covs


def load_data(use_jct=False, use_txr=True):
    if use_jct:
        #JCT = pickle.load(open("./%s/compiled/jct_%s.pkl" % (DATA_VER, EVENT_TYPE), "rb"))
        JCT = DataTableIndexed("./%s/compiled/jct_%s.txt" % (DATA_VER, EVENT_TYPE))
    else:
        JCT = None
    if use_txr:
        #TXR = pickle.load(open("./%s/compiled/txr_%s.pkl" % (DATA_VER, EVENT_TYPE), "rb"))
        TXR = DataTableIndexed("./%s/compiled/txr_%s.txt" % (DATA_VER, EVENT_TYPE))
    else:
        TXR = None
    return JCT, TXR


def read_exonset():
    exonset1 = ExonSet.from_rmats("./data/rmats/fromGTF.%s.txt"%EVENT_TYPE, event_type=EVENT_TYPE)
    exonset2 = ExonSet.from_suppa("/mnt/ceph/users/zzhang/SUPPA/index/hg38/suppa_gencodev34_%s_strict.ioe"%EVENT_TYPE,
                                 cache_dir='./data/',
                                 event_type=EVENT_TYPE)
    exonset = exonset1 + exonset2
    if USE_ANNOT:
        exonset.read_coding_region_from_gtf_by_splice_site('/mnt/home/zzhang/ceph/genome_annotation/gencode.v34.annotation.gtf')
        exonset.read_pfam_domain_from_bed_by_splice_site('./data/pfam/hg38_exon_pfam_annot.bed.gz',
                                                         clan_annot_fp='./data/pfam/Pfam-A.clans.tsv.gz')
        exonset.read_rbp_clip([
            './data/encode_eclip/hg38_gencodev34.all_splice_sites.encodeIDR.bed',
            './data/encode_eclip/hg38_gencodev34.all_splice_sites.clipdb.bed'
            ])
    return exonset


def run_ctrl(JCT, TXR, exonset=None):
    contrasts = Contrasts(name="final", levels=['Control'])
    cov_ctrl = Covariate(fp="./%s/charm_master.clean_w_10pc.csv" % DATA_VER, sep=",", index_col=0,
                         contrasts=contrasts,
                         main_effects=['Sex', 'tp', 'plateNum', 'pid'] + PCS_TO_INCL if USE_RE is True \
                             else ['Sex', 'tp', 'plateNum'] + PCS_TO_INCL,
                         interaction_effects=['Sex|tp'],
                         verbose=True
                         )
    LOGGER.info("CTRL formula: %s" % cov_ctrl.formula)
    jem_ctrl = JemmLinearMixedModel(
        junction_measure=JCT,
        transcript_measure=TXR,
        covariates=cov_ctrl,
        diff_intercept_by_measure=True,
        group_varname='pid',
        min_groupvar=MIN_RE_VAR,
        optimizer='bfgs'
    )

    _ = jem_ctrl.run_tests(test_type='Wald',
                      data="txr",
                      force_diff_intercept=False,
                      force_rerun=True,
                      pval_adjust_method="fdr",
                      nthreads=12
                      )
    jem_ctrl.save_regression_table(ctrl_fp,
                           exonset=exonset
                           )
    black_list_exons = [e for e in jem_ctrl.stats_tests if 
                    jem_ctrl.stats_tests[e].loc['tp', 'qvals']<FDR_THRESH or
                    jem_ctrl.stats_tests[e].loc['Sex@M|tp', 'qvals']<FDR_THRESH
                   ]
    return black_list_exons


def run_disease(JCT, TXR, black_list_exons=None, exonset=None):
    global DATA_VER
    covs = get_covs(DATA_VER=DATA_VER)
    LOGGER.info("Disease model: %s" % covs.formula)
    # init jemm
    Jemm = JemmLinearMixedModel if USE_RE else JemmLinearRegression
    jem = Jemm(junction_measure=JCT, 
               transcript_measure=TXR, 
               covariates=covs,
               diff_intercept_by_measure=True,
               group_varname='pid',
               min_groupvar=MIN_RE_VAR, 
               optimizer='bfgs',
    )
    black_list_exons = black_list_exons or []
    LOGGER.debug("Before filter: %i"%len(jem.event_index))
    LOGGER.debug("Black list=%i"%len(black_list_exons))
    jem.event_index = [e for e in jem.event_index if e not in black_list_exons]
    LOGGER.debug("After filter: %i"%len(jem.event_index))

    # data stats
    try:
        LOGGER.debug('unique pids = %i' % len(jem.covariates.covariate.pid.unique()))
        LOGGER.debug('total sids = %i' % len(jem.covariates.covariate.index.unique()))
    except:
        pass


    # run tests if cannot find previous results
    _ = jem.run_tests(test_type='Wald',
                      data="txr",
                      force_diff_intercept=False,
                      force_rerun=True,
                      pval_adjust_method="fdr",
                      nthreads=24
                     )
    # if new run, save the results for future re-use
    jem.munge_covariates([
            'final@Pre',
            'final@First',
            'final@Mid',
            'final@Post',
        ],
        meta_name = "final4cond"
        )
    outfp = os.path.join(OUTDIR, "joint.%s.reg_table.tsv"%EVENT_TYPE)
    reg_table = jem.save_regression_table(outfp,
                              exonset=exonset,
                              annotations=['coding_region', 'pfam_domain', 'rbp_clip'] if USE_ANNOT else None,
                              order_by_covariate='final4cond',
                              order_by='logP'
                             )

    return jem


def main():
    LOGGER.info("load data")
    jct, txr = load_data(use_jct=(EVENT_TYPE=="SE"), use_txr=True)
    LOGGER.info("load annot")
    exonset = read_exonset()
    #LOGGER.info("run ctrl")
    #black_list_exons = run_ctrl(JCT=jct, TXR=txr, exonset=exonset)
    LOGGER.info("run tests")
    run_disease(JCT=jct, TXR=txr, exonset=exonset)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', type=str, required=True)
    ap.add_argument('--event-type', type=str, choices=['SE', 'A3SS', 'A5SS', 'RI'], required=True)
    try:
        args = ap.parse_args()

        DATA_VER = args.data
        EVENT_TYPE = args.event_type
        OUTDIR = './%s/joint_%s/' % (DATA_VER, EVENT_TYPE)
        os.makedirs(OUTDIR, exist_ok=True)

        main()
    except Exception as e:
        print(e)
        pass
