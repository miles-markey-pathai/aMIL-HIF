import numpy as np
import pandas as pd
from scipy import stats 
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests 
from typing import Tuple

def compare_hif_distributions(hifs: pd.DataFrame,
                              group_label: pd.Series,
                              multipletests_method: str = 'fdr_bh',
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Performs Kolmogorov-Smirnov and Mann-Whitney U test for each HIF."""
    
    ks_dict = {}
    mwh_dict = {}
    groups = np.unique(group_label)
    
    for hif in hifs.columns:
        ks = stats.ks_2samp(hifs.loc[group_label == groups[0], hif].dropna(),
                            hifs.loc[group_label == groups[1], hif].dropna())
        ks_dict[hif] = [ks[0], ks[1]]
        u_stat, pval = stats.mannwhitneyu(hifs.loc[group_label == groups[0], hif].dropna(),
                                          hifs.loc[group_label == groups[1], hif].dropna())
        mwh_dict[hif] = [u_stat, pval]
        
    df_ks = pd.DataFrame.from_dict(ks_dict, orient='index', columns=['statistic', 'pvalue'])
    df_ks['pvalue_corrected'] = multipletests(df_ks['pvalue'], method=multipletests_method)[1]
    df_ks = df_ks.sort_values(by='pvalue')
    
    df_mwh = pd.DataFrame.from_dict(mwh_dict, orient='index', columns=['statistic', 'pvalue'])
    df_mwh['pvalue_corrected'] = multipletests(df_mwh['pvalue'], method=multipletests_method)[1]
    df_mwh = df_mwh.sort_values(by='pvalue')

    return df_ks, df_mwh
