import re
import pandas as pd

def get_property(feature):
    """Extract property from feature."""
    return feature.split(']')[0].split('NUCLEUS_')[-1]

def get_property_type(feature):
    """Extract property type from feature."""
    
    prop = get_property(feature)
    
    if 'CHANNEL' in prop:
        return 'Color'
    
    keywords = ['AREA','LENGTH','PERIMETER']
    for k in keywords:
        if k in feature:
            return 'Size'
    return 'Shape'


def get_hif_type(feature):
    return feature.split('[')[0]
    
def get_cell(feature):
    """Extract cell from feature."""
    
    # nuclear hif 
    if any(substring in feature for substring in ['IQR', 'MEAN', 'MEDIAN', 'STD']):  
        return feature.split('[')[1].split('_')[0]
    # tissue/cell hif 
    elif ('SELECTED BY' in feature):
        match = re.search(r'^[a-z ]{1,}\[\[\[(?P<cell>[- \w]+) [CELLS]+\] SELECTED BY', 
                  feature, re.IGNORECASE).group('cell')
        return match 
    elif ('DENSITY RATIO' in feature) :
        return None
    elif 'CELL' in feature:
        if ('DENSITY' not in feature) or \
           ((('DENSITY') in feature) and ('OVER' in feature)):
            match = re.search(r'^[a-z ]{1,}\[{1,2}(?P<cell>[- \w]+) [CELLS]+\] OVER', 
                              feature, re.IGNORECASE).group('cell')
        else: 
            match = re.search(r'^[a-z ]{1,}\[{1,2}(?P<cell>[- \w]+) [CELLS]+\] IN', 
                              feature, re.IGNORECASE).group('cell')
        return match
    else:
        return None

    
    
def get_tissue(feature):
    if ('COUNT PROP' in feature) or ('DENSITY' in feature):
        match = re.search(r'\] IN \[(?P<tissue>.*)\]', 
                          feature, re.IGNORECASE).group('tissue')
        return match
    elif 'AREA PROP' in feature:
        match = re.search(r'^AREA PROP \[\[(?P<tissue>.*)\] OVER', 
                          feature, re.IGNORECASE).group('tissue')
        return match

def standardize_hif_names(hif_df):
    for feature in list(hif_df):
        newFeat = feature
        #newFeat = feature[:-11]
        if ('PD-L1-HIGH') in newFeat or ('PD-L1-LOW') in newFeat:
            newFeatList = newFeat.split('[')
            newFeatList[-1] = 'HIGH ATTENTION REGIONS]_HE'
            newFeat = '['.join(newFeatList)
        else:
            newFeatList = newFeat.split('_')
            newFeatList[-1] = 'HE'
            newFeat = '_'.join(newFeatList)
        hif_df.rename(columns={feature:newFeat}, inplace=True)
    return hif_df
def standardize_area_hif_names(hif_df):
    for feature in list(hif_df):
        newFeat = feature[:-11]
        if ('EXCITATORY REGIONS') in newFeat or ('INHIBITORY REGIONS') in newFeat:
            newFeatList = newFeat.split('[')
            newFeatList[3] = 'HIGH ATTENTION REGIONS]] IN '
            newFeat = '['.join(newFeatList)
        hif_df.rename(columns={feature:newFeat}, inplace=True)
    return hif_df

def get_hif_df(mil_df, tissue_low, tissue_high, cell_low, cell_high, drop_slides, drop_columns):
    mil_df.set_index(["slideID"], inplace=True)
    mil_df = mil_df[['mean_score', 'model_prediction', 'tgfbcaf_pos', '__set_TGFb-CAF']]
    mil_df.rename(columns={'pdl1_pos':'Class'}, inplace=True)
    low_df = mil_df[mil_df['model_prediction']==0]
    high_df = mil_df[mil_df['model_prediction']==1]
    tissue_high.set_index(["H & E_ID"], inplace=True)
    drop_abs_prefix = 'MM2'
    drop_abs_list = list()
    for hif in list(tissue_high):    
        if drop_abs_prefix in hif.split("[")[0].strip():
            drop_abs_list.append(hif) 
    tissue_high.drop(columns=drop_abs_list, inplace=True)
    tissue_high.drop(columns=["CASE_ID"], inplace=True)
    tissue_high = standardize_area_hif_names(tissue_high)
    tissue_low.set_index(["H & E_ID"], inplace=True)
    drop_abs_prefix = 'MM2'
    drop_abs_list = list()
    for hif in list(tissue_low):    
        if drop_abs_prefix in hif.split("[")[0].strip():
            drop_abs_list.append(hif) 
    tissue_low.drop(columns=drop_abs_list, inplace=True)
    tissue_low.drop(columns=["CASE_ID"], inplace=True)
    tissue_low = standardize_area_hif_names(tissue_low)
    cell_high.set_index(["H & E_ID"], inplace=True)
    drop_abs_prefix = ['TOTAL']
    drop_abs_list = list()
    for hif in list(cell_high):
        if hif.split("[")[0].strip() in drop_abs_prefix:
            drop_abs_list.append(hif)
        if 'NORMAL' in hif: 
            drop_abs_list.append(hif)
        if 'NECROSIS' in hif: 
            drop_abs_list.append(hif)
        if 'PD-L1-LOW' in hif: 
            drop_abs_list.append(hif)
    cell_high.drop(columns=drop_abs_list, inplace=True)
    cell_high.drop(columns=["CASE_ID"], inplace=True)
    cell_high = standardize_hif_names(cell_high)
    cell_low.set_index(["H & E_ID"], inplace=True)
    drop_abs_prefix = ['TOTAL']
    drop_abs_list = list()
    for hif in list(cell_low):
        if hif.split("[")[0].strip() in drop_abs_prefix:
            drop_abs_list.append(hif)
        if 'NORMAL' in hif: 
            drop_abs_list.append(hif)
        if 'NECROSIS' in hif: 
            drop_abs_list.append(hif)
        if 'PD-L1-HIGH' in hif: 
            drop_abs_list.append(hif)
    cell_low.drop(columns=drop_abs_list, inplace=True)
    cell_low.drop(columns=["CASE_ID"], inplace=True)
    cell_low = standardize_hif_names(cell_low)
    hif_df_e = tissue_high.join(cell_high)
    hif_df_i = tissue_low.join(cell_low)
    hif_df_e = hif_df_e.loc[high_df.index]
    hif_df_i = hif_df_i.loc[low_df.index]
    hif_df = pd.concat([hif_df_e, hif_df_i])
    hif_df.drop(index=drop_slides, inplace=True)
    hif_df.drop(columns=drop_columns, inplace=True)
    hif_mil_df = hif_df.join(mil_df)
    return hif_mil_df
    