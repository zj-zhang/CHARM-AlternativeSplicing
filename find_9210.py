
import os
from tqdm import tqdm
import pickle
import pandas as pd
import numpy as np
import pickle
from jemm import rmats_helper, suppa_helper
from jemm import kallisto_helper
from jemm.data_table import DataTableIndexed
from jemm.junction import JunctionCountTable
from jemm.transcript import TranscriptMeasureTable

from jemm.covariate import Covariate, Contrasts
from jemm.model import JemmLinearMixedModel, JemmLinearRegression


DATA_VER = 'data-V7'
EVENT_TYPE = 'SE'
FDR_THRESH = 0.05
PCS_TO_INCL = []


# Model method to use
USE_RE = True
MIN_RE_VAR = 0.001  # sqrt(0.001) = 0.03, 0.03*1.96=0.06~5% change
#JCT = JunctionCountTable.from_plaintext("./%s/compiled/jct_%s.JemmPlainText.txt" % (DATA_VER, EVENT_TYPE) )
#TXR = TranscriptMeasureTable.from_plaintext("./%s/compiled/txr_%s.JemmPlainText.txt" % (DATA_VER, EVENT_TYPE))
JCT = pickle.load(open("./%s/compiled/jct_%s.reduced.pkl" % (DATA_VER, EVENT_TYPE), "rb"))
TXR = pickle.load(open("./%s/compiled/txr_%s.reduced.pkl" % (DATA_VER, EVENT_TYPE), "rb"))
#JCT = DataTableIndexed("./%s/compiled/jct_%s.txt" % (DATA_VER, EVENT_TYPE) )
#TXR = DataTableIndexed("./%s/compiled/txr_%s.txt" % (DATA_VER, EVENT_TYPE))


black_list_exons = pickle.load(open("./%s/joint_%s/black_list_exons.pkl"%(DATA_VER, EVENT_TYPE), "rb"))

contrasts = Contrasts(name="final", levels=[
    'Control',
    'Pre',
    'First',
    'Mid',
    'Post',
    'Asymptomatic',
    'Mild',
    'Moderate',
    'Exposed',
    'False Negative',
    'Immune',
    'Reinfection'
])

covs = Covariate(fp="./%s/charm_master.clean_w_10pc.csv" % DATA_VER, sep=",", 
                 index_col=0,
                 contrasts=contrasts,
                 main_effects=['final', 'pid','Sex'] + PCS_TO_INCL if USE_RE is True else ['final', 'Sex'] + PCS_TO_INCL,
                 factor_conversion={
                     'final': {
                         'Asymptomatic': 'final@Asymptomatic',
                         'Exposed': 'final@Exposed',
                         'False Negative': 'final@False Negative',
                         'First': 'final@First',
                         'Immune': 'final@Immune',
                         'Mid': 'final@Mid',
                         'Mild': 'final@Mild',
                         'Moderate': 'final@Moderate',
                         'Post': 'final@Post',
                         'Pre': 'final@Pre',
                         'Reinfection': 'final@Reinfection'
                     },
                     'Sex': {'M': 'Sex@M'}
                 },  
                 verbose=True
             )


jmm = JemmLinearRegression(junction_measure=JCT,
           transcript_measure=TXR,
           covariates=covs,
           diff_intercept_by_measure=True)
print("Before filter: %i"%len(jmm.event_index))
print("Black list=%i"%len(black_list_exons))
jmm.event_index = [e for e in jmm.event_index if e not in black_list_exons]
print("After filter: %i"%len(jmm.event_index))

_ = jmm.run_tests(test_type='Wald',
                  data="both",
                  force_diff_intercept=True,
                  force_rerun=True,
                  pval_adjust_method="fdr",
                  nthreads=1
                 )

tested_eids = [e for e in jmm.event_index if e in jmm.stats_tests]


from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import scipy.stats as ss
import numpy as np

txr_dfs = {'SE': jmm.transcript_measure.applymap(lambda x: x.psi)}

pca = PCA(n_components=10, random_state=777)
a = txr_dfs['SE'].loc[tested_eids]
pca.fit(scale(a.fillna(a.mean(axis=0)).fillna(0.5), axis=1))


meta = jmm.covariates.meta
pcs = []
for i in range(len(pca.components_)):
    meta["PC%i"%i] = pca.components_[i]
    sp = ss.spearmanr(pca.components_[i], meta['final'])
    sp2 = ss.spearmanr(pca.components_[i], meta['Sex'])
    sp3 = ss.spearmanr(pca.components_[i], meta['Age'])
    sp4 = ss.spearmanr(pca.components_[i], meta['tp'])

    print('PC%i, Final Spearman=%.3f, pval=%.3f' % (i, sp.correlation, sp.pvalue))
    print('      Sex   Spearman=%.3f, pval=%.3f' % (sp2.correlation, sp2.pvalue))
    print('      Age   Spearman=%.3f, pval=%.3f' % (sp3.correlation, sp3.pvalue))
    print('      TP    Spearman=%.3f, pval=%.3f' % (sp4.correlation, sp4.pvalue))
    if sp.pvalue > 0.01 and sp2.pvalue > 0.01 and sp3.pvalue > 0.01:
        pcs.append(i)
print('PCs w/o strong corr.: %s' % pcs)




