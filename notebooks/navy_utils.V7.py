from jemm.covariate import Contrasts, Covariate
import sys
DATA_VER = sys.argv[1]

# Plot style
from matplotlib import rcParams
from matplotlib import font_manager
#font_manager.fontManager.addfont('./data/Helvetica.ttf')
#rcParams['font.family'] = 'Helvetica'
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False
rcParams['xtick.top'] = False
rcParams['ytick.right'] = False
rcParams['figure.figsize'] = (4,4)
rcParams['figure.dpi'] = 600
rcParams['axes.titlesize'] = 14
rcParams['axes.labelsize'] = 14
rcParams['lines.linewidth'] = 1.5
rcParams['lines.markersize'] = 8
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
rcParams['legend.fontsize'] = 14

try:
    PCS_TO_INCL = sys.argv[2]
    PCS_TO_INCL = ['PC%i'%int(i) for i in PCS_TO_INCL.split(',')]
except:
    PCS_TO_INCL = []
try:
    USE_RE = sys.argv[3]
    USE_RE = True if USE_RE=='True' else False
except:
    USE_RE = False

print('navy_utils.py {data} {pcs} {re}'.format(data=DATA_VER, pcs=str(PCS_TO_INCL), re=USE_RE))

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

joint_data_files = {
    'SE': {
        'jct': './%s/compiled/jct_SE.txt' % DATA_VER,
        'txr': './%s/compiled/txr_SE.txt' % DATA_VER
    },
    'A5SS': {
        'jct': './%s/compiled/jct_A5SS.reduced.pkl' % DATA_VER,
        'txr': './%s/compiled/txr_A5SS.pkl' % DATA_VER
    },
    'A3SS': {
        'jct': './%s/compiled/jct_A3SS.reduced.pkl'% DATA_VER,
        'txr': './%s/compiled/txr_A3SS.pkl'% DATA_VER
    },
    'RI': {
        'jct': './%s/compiled/jct_RI.reduced.pkl'% DATA_VER,
        'txr': './%s/compiled/txr_RI.pkl'% DATA_VER
    }
}

joint_reg_tables = {
    'SE': './%s/joint_SE/joint.SE.reg_table.tsv'% DATA_VER,
    'A5SS': './%s/joint_A5SS/joint.A5SS.reg_table.tsv'% DATA_VER,
    'A3SS': './%s/joint_A3SS/joint.A3SS.reg_table.tsv'% DATA_VER,
    'RI': './%s/joint_RI/joint.RI.reg_table.tsv'% DATA_VER
}

m_reg_tables = {
    'SE': './%s/gender_SE/male.SE.reg_table.tsv'% DATA_VER,
    'A5SS': './%s/gender_A5SS/male.A5SS.reg_table.tsv'% DATA_VER,
    'A3SS': './%s/gender_A3SS/male.A3SS.reg_table.tsv'% DATA_VER,
    'RI': './%s/gender_RI/male.RI.reg_table.tsv'% DATA_VER
}

f_reg_tables = {
    'SE': './%s/gender_SE/female.SE.reg_table.tsv'% DATA_VER,
    'A5SS': './%s/gender_A5SS/female.A5SS.reg_table.tsv'% DATA_VER,
    'A3SS': './%s/gender_A3SS/female.A3SS.reg_table.tsv'% DATA_VER,
    'RI': './%s/gender_RI/female.RI.reg_table.tsv'% DATA_VER
}
