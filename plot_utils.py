import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D 
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import pandas as pd
import numpy as np 
import seaborn as sns 
from statannotations.Annotator import Annotator
import textwrap
from typing import List, Tuple

import utils 

# pathAI palette 
pa_colors = pd.Series({'spacecadet': '#2D2849',
 'darkgreen': '#195866',
 'mediumgreen': '#2AA095',
 'violet': '#793EF8',
 'turquoise': '#4CEAD3',
 'raspberry': '#C5304B',
 'darkgray': '#333132',
 'lightgray': '#EEEEEE',
 'lightviolet': '#EBE2FE',
 'lightgreen': '#E1EFED',
 '_darkblue': '#173D52',
 '_lightblue': '#16B2E1',
 '_purple': '#852693',
 '_red': '#E11C29'})


def boxplots_hifs_by_group(hifs: pd.DataFrame,
                           df_pvals: pd.DataFrame, 
                           regex: str, 
                           hue: pd.Series,
                           hue_colors: dict, 
                           groupx: List[str],
                           fname: str = None,
                           groupx_name: str = 'cell',
                           groupx_labels: List[str] = None,
                           ylabel: str = None,
                           figsize=(16, 7),
) -> None:
    """Plots boxplot of a given feature by cell type.
    
    Args:
        hifs: A dataframe whose rows represent samples and columns represent features. 
        df_pvals: A dataframe whose rows represent features and 
            columns include pvalues for annotation.
        regex: A string (regular expression) that subsets features for plotting.
        hue: A list of group labels for each row in `hifs` dataframe. 
        hue_colors: A dictionary of colors for each group label. 
        groupx: A list of strings indicating cell types. 
        fname: A string for file name/path. 
            Defaults to None, in which case figure is not saved. 
        groupx_name: A string indicating cell or tissue. 
            Use `cell` for cell or `tissue` for tissue. 
        groupx_labels: A list of strings for setting x-axis label. 
            Defaults to None, in which case `groupx` is used for x-axis label.
        ylabel: A string for setting y-axis label. 
            Defaults to None, in which case `regex` is used for y-axis label. 
        figsize: A tuple of floats specifying width, height in inches.
            Defaults to (16, 7).

    Returns:
        None
    """
    
    if groupx_labels is None: 
        groupx_labels = groupx
    if ylabel is None:
        ylabel = regex 
        
    hifs_filtered = hifs.filter(regex=regex)
    #print(hifs_filtered)
    names_hifs = list(hifs_filtered)
    hifs_filtered = hifs_filtered.merge(hue, on='H & E_ID')
    hifs_filtered_melt = pd.melt(hifs_filtered, id_vars=hue.name, 
                                 value_vars=names_hifs)
    if groupx_name == 'cell':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_cell, hifs_filtered_melt.variable))
    elif groupx_name == 'tissue':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_tissue, hifs_filtered_melt.variable))
    groups = list(hue_colors.keys())
    legend_elements = [Line2D([0], [0], color=hue_colors[groups[0]], marker='o', 
                              alpha=0.5, lw=15, label=groups[0]),
                       Line2D([0], [0], color=hue_colors[groups[1]], marker='o', 
                              alpha=0.5, lw=15, label=groups[1])]
    pairs = [((cell, groups[0]), (cell, groups[1])) for cell in groupx] 
    pvalues_corrected = []
    for groupx_val in groupx:
        hif_name = hifs_filtered.columns[hifs_filtered.columns.str.contains(f'\[{groupx_val}')][0]
        pval_tmp = df_pvals.loc[df_pvals.index == hif_name, 'pvalue_corrected'].squeeze()
        
        if pval_tmp < 1e-20:
            pvalues_corrected.append('****') 
        elif pval_tmp < 1e-10:
            pvalues_corrected.append('***')
        elif pval_tmp < 1e-5:
            pvalues_corrected.append('**')
        elif pval_tmp < 0.05:
            pvalues_corrected.append('*')
        elif pval_tmp >= 0.05:
            pvalues_corrected.append('ns')
            
    fig, ax = plt.subplots(figsize=figsize)
    ax = sns.boxplot(x=groupx_name, y='value', hue=hue.name, 
                hue_order=list(hue_colors.keys()),
                data=hifs_filtered_melt, palette=hue_colors,
                order=groupx)
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    sns.swarmplot(x=groupx_name, y='value', hue=hue.name, dodge=True,
                  hue_order=list(hue_colors.keys()),
                  order=groupx, data=hifs_filtered_melt, ax=ax, 
                  palette=hue_colors, size=2).set(ylabel="", xlabel="")

    annotator=Annotator(ax, pairs, data=hifs_filtered_melt,
                        x=groupx_name, y='value', hue=hue.name,
                        order=groupx)
    annotator.configure(test=None, text_format='simple', fontsize=15) #, pvalue_format='simple')
    annotator.set_custom_annotations(pvalues_corrected)
    annotator.annotate()

    ax.legend(handles=legend_elements, loc='best', fontsize=15, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #ax.set_xticklabels(groupx_labels)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlabel('')
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500)

    return None

def boxplots_hifs_by_group3(hifs1: pd.DataFrame,
                           hifs2: pd.DataFrame,
                           df_pvals1: pd.DataFrame, 
                           df_pvals2: pd.DataFrame, 
                           regex: str, 
                           hue: pd.Series,
                           hue_colors: dict, 
                           groupx: List[str],
                           fname: str = None,
                           groupx_name: str = 'cell',
                           groupx_labels: List[str] = None,
                           ylabel: str = None,
                           figsize=(16, 7),
) -> None:
    """Plots boxplot of a given feature by cell type.
    
    Args:
        hifs: A dataframe whose rows represent samples and columns represent features. 
        df_pvals: A dataframe whose rows represent features and 
            columns include pvalues for annotation.
        regex: A string (regular expression) that subsets features for plotting.
        hue: A list of group labels for each row in `hifs` dataframe. 
        hue_colors: A dictionary of colors for each group label. 
        groupx: A list of strings indicating cell types. 
        fname: A string for file name/path. 
            Defaults to None, in which case figure is not saved. 
        groupx_name: A string indicating cell or tissue. 
            Use `cell` for cell or `tissue` for tissue. 
        groupx_labels: A list of strings for setting x-axis label. 
            Defaults to None, in which case `groupx` is used for x-axis label.
        ylabel: A string for setting y-axis label. 
            Defaults to None, in which case `regex` is used for y-axis label. 
        figsize: A tuple of floats specifying width, height in inches.
            Defaults to (16, 7).

    Returns:
        None
    """
    
    if groupx_labels is None: 
        groupx_labels = groupx
    if ylabel is None:
        ylabel = regex 
        
    hifs1_filtered = hifs1.filter(regex=regex)
    hifs2_filtered = hifs2.filter(regex=regex)
    #print(hifs_filtered)
    names_hifs = list(hifs1_filtered)
    hifs1_filtered = hifs1_filtered.merge(hue, on='H & E_ID')
    hifs1_filtered_melt = pd.melt(hifs1_filtered, id_vars=hue.name, 
                                 value_vars=names_hifs)
    hifs2_filtered = hifs2_filtered.merge(hue, on='H & E_ID')
    hifs2_filtered_melt = pd.melt(hifs2_filtered, id_vars=hue.name, 
                                 value_vars=names_hifs)
    if groupx_name == 'cell':
        hifs1_filtered_melt[groupx_name] = list(map(utils.get_cell, hifs1_filtered_melt.variable))
    elif groupx_name == 'tissue':
        hifs1_filtered_melt[groupx_name] = list(map(utils.get_tissue, hifs1_filtered_melt.variable))
    groups = list(hue_colors.keys())
    legend_elements = [Line2D([0], [0], color=hue_colors[groups[0]], marker='o', 
                              alpha=0.5, lw=15, label=groups[0]),
                       Line2D([0], [0], color=hue_colors[groups[1]], marker='o', 
                              alpha=0.5, lw=15, label=groups[1])]
    pairs = [((cell, groups[0]), (cell, groups[1])) for cell in groupx] 
    pvalues_corrected = []
    for groupx_val in groupx:
        hif_name = hifs1_filtered.columns[hifs1_filtered.columns.str.contains(f'\[{groupx_val}')][0]
        pval_tmp = df_pvals1.loc[df_pvals1.index == hif_name, 'pvalue_corrected'].squeeze()
        
        if pval_tmp < 1e-20:
            pvalues_corrected.append('****') 
        elif pval_tmp < 1e-10:
            pvalues_corrected.append('***')
        elif pval_tmp < 1e-5:
            pvalues_corrected.append('**')
        elif pval_tmp < 0.05:
            pvalues_corrected.append('*')
        elif pval_tmp >= 0.05:
            pvalues_corrected.append('ns')
            
    fig, ax = plt.subplots(figsize=figsize)
    ax = sns.boxplot(x=groupx_name, y='value', hue=hue.name, 
                hue_order=list(hue_colors.keys()),
                data=hifs1_filtered_melt, palette=hue_colors,
                order=groupx)
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    sns.swarmplot(x=groupx_name, y='value', hue=hue.name, dodge=True,
                  hue_order=list(hue_colors.keys()),
                  order=groupx, data=hifs1_filtered_melt, ax=ax, 
                  palette=hue_colors, size=2).set(ylabel="", xlabel="")

    annotator=Annotator(ax, pairs, data=hifs1_filtered_melt,
                        x=groupx_name, y='value', hue=hue.name,
                        order=groupx)
    annotator.configure(test=None, text_format='simple', fontsize=15) #, pvalue_format='simple')
    annotator.set_custom_annotations(pvalues_corrected)
    annotator.annotate()

    ax.legend(handles=legend_elements, loc='best', fontsize=15, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #ax.set_xticklabels(groupx_labels)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlabel('')
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500)

    return None

def boxplots_hifs_by_group_two_region(hifs: pd.DataFrame,
                           df_pvals: pd.DataFrame, 
                           regex: str, 
                           hue: pd.Series,
                           hue_colors: dict, 
                           groupx: List[str],
                           fname: str = None,
                           groupx_name: str = 'cell',
                           groupx_labels: List[str] = None,
                           ylabel: str = None,
                           figsize=(16, 7),
) -> None:
    """Plots boxplot of a given feature by cell type.
    
    Args:
        hifs: A dataframe whose rows represent samples and columns represent features. 
        df_pvals: A dataframe whose rows represent features and 
            columns include pvalues for annotation.
        regex: A string (regular expression) that subsets features for plotting.
        hue: A list of group labels for each row in `hifs` dataframe. 
        hue_colors: A dictionary of colors for each group label. 
        groupx: A list of strings indicating cell types. 
        fname: A string for file name/path. 
            Defaults to None, in which case figure is not saved. 
        groupx_name: A string indicating cell or tissue. 
            Use `cell` for cell or `tissue` for tissue. 
        groupx_labels: A list of strings for setting x-axis label. 
            Defaults to None, in which case `groupx` is used for x-axis label.
        ylabel: A string for setting y-axis label. 
            Defaults to None, in which case `regex` is used for y-axis label. 
        figsize: A tuple of floats specifying width, height in inches.
            Defaults to (16, 7).

    Returns:
        None
    """
    
    if groupx_labels is None: 
        groupx_labels = groupx
    if ylabel is None:
        ylabel = regex 
        
    hifs_filtered = hifs.filter(regex=regex)
    #print(hifs_filtered)
    names_hifs = list(hifs_filtered)
    hifs_filtered = hifs_filtered.merge(hue, on='H & E_ID')
    hifs_filtered_melt = pd.melt(hifs_filtered, id_vars=hue.name, 
                                 value_vars=names_hifs)
    if groupx_name == 'cell':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_cell, hifs_filtered_melt.variable))
    elif groupx_name == 'tissue':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_tissue, hifs_filtered_melt.variable))
    groups = list(hue_colors.keys())
    legend_elements = [Line2D([0], [0], color=hue_colors[groups[0]], marker='o', 
                              alpha=0.5, lw=15, label=groups[0]),
                       Line2D([0], [0], color=hue_colors[groups[1]], marker='o', 
                              alpha=0.5, lw=15, label=groups[1])]
    pairs = [((cell, groups[0]), (cell, groups[1])) for cell in groupx] 
    pvalues_corrected = []
    for groupx_val in groupx:
        hif_name = hifs_filtered.columns[hifs_filtered.columns.str.contains(f'\[{groupx_val}')][0]
        pval_tmp = df_pvals.loc[df_pvals.index == hif_name, 'pvalue_corrected'].squeeze()
        
        if pval_tmp < 1e-20:
            pvalues_corrected.append('****') 
        elif pval_tmp < 1e-10:
            pvalues_corrected.append('***')
        elif pval_tmp < 1e-5:
            pvalues_corrected.append('**')
        elif pval_tmp < 0.05:
            pvalues_corrected.append('*')
        elif pval_tmp >= 0.05:
            pvalues_corrected.append('ns')
            
    fig, ax = plt.subplots(figsize=figsize)
    ax = sns.boxplot(x=groupx_name, y='value', hue=hue.name, 
                hue_order=list(hue_colors.keys()),
                data=hifs_filtered_melt, palette=hue_colors,
                order=groupx)
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    sns.swarmplot(x=groupx_name, y='value', hue=hue.name, dodge=True,
                  hue_order=list(hue_colors.keys()),
                  order=groupx, data=hifs_filtered_melt, ax=ax, 
                  palette=hue_colors, size=2).set(ylabel="", xlabel="")

    annotator=Annotator(ax, pairs, data=hifs_filtered_melt,
                        x=groupx_name, y='value', hue=hue.name,
                        order=groupx)
    annotator.configure(test=None, text_format='simple', fontsize=15) #, pvalue_format='simple')
    annotator.set_custom_annotations(pvalues_corrected)
    annotator.annotate()

    ax.legend(handles=legend_elements, loc='best', fontsize=15, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #ax.set_xticklabels(groupx_labels)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlabel('')
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500)

    return None

def boxplots_hifs_by_group2(hifs: pd.DataFrame,
                           df_pvals: pd.DataFrame, 
                           regex: str, 
                           hue: pd.Series,
                           hue_colors: dict, 
                           groupx: List[str],
                           fname: str = None,
                           groupx_name: str = 'cell',
                           groupx_labels: List[str] = None,
                           ylabel: str = None,
                           figsize=(16, 7),
) -> None:
    """Plots boxplot of a given feature by cell type.
    
    Args:
        hifs: A dataframe whose rows represent samples and columns represent features. 
        df_pvals: A dataframe whose rows represent features and 
            columns include pvalues for annotation.
        regex: A string (regular expression) that subsets features for plotting.
        hue: A list of group labels for each row in `hifs` dataframe. 
        hue_colors: A dictionary of colors for each group label. 
        groupx: A list of strings indicating cell types. 
        fname: A string for file name/path. 
            Defaults to None, in which case figure is not saved. 
        groupx_name: A string indicating cell or tissue. 
            Use `cell` for cell or `tissue` for tissue. 
        groupx_labels: A list of strings for setting x-axis label. 
            Defaults to None, in which case `groupx` is used for x-axis label.
        ylabel: A string for setting y-axis label. 
            Defaults to None, in which case `regex` is used for y-axis label. 
        figsize: A tuple of floats specifying width, height in inches.
            Defaults to (16, 7).

    Returns:
        None
    """
    
    if groupx_labels is None: 
        groupx_labels = groupx
    if ylabel is None:
        ylabel = regex 
        
    hifs_filtered = hifs.filter(regex=regex)
    names_hifs = list(hifs_filtered)
    hifs_filtered = hifs_filtered.merge(hue, on='H & E_ID')
    hifs_filtered_melt = pd.melt(hifs_filtered, id_vars=hue.name, 
                                 value_vars=names_hifs)
    if groupx_name == 'cell':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_cell, hifs_filtered_melt.variable))
    elif groupx_name == 'tissue':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_tissue, hifs_filtered_melt.variable))
    groups = list(hue_colors.keys())
    legend_elements = [Line2D([0], [0], color=hue_colors[groups[0]], marker='o', 
                              alpha=0.5, lw=15, label=groups[0]),
                       Line2D([0], [0], color=hue_colors[groups[1]], marker='o', 
                              alpha=0.5, lw=15, label=groups[1])]
    pairs = [((cell, groups[0]), (cell, groups[1])) for cell in groupx] 
    pvalues_corrected = []
    for groupx_val in groupx:
        hif_name = hifs_filtered.columns[hifs_filtered.columns.str.contains(f'\[{groupx_val}')][0]
        pval_tmp = df_pvals.loc[df_pvals.index == hif_name, 'pvalue_corrected'].squeeze()
        
        if pval_tmp < 1e-20:
            pvalues_corrected.append('****') 
        elif pval_tmp < 1e-10:
            pvalues_corrected.append('***')
        elif pval_tmp < 1e-5:
            pvalues_corrected.append('**')
        elif pval_tmp < 0.05:
            pvalues_corrected.append('*')
        elif pval_tmp >= 0.05:
            pvalues_corrected.append('ns')
            
    fig, ax = plt.subplots(figsize=figsize)
    ax = sns.boxplot(x=groupx_name, y='value', hue=hue.name, 
                hue_order=list(hue_colors.keys()),
                data=hifs_filtered_melt, palette=hue_colors,
                order=groupx)
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    sns.swarmplot(x=groupx_name, y='value', hue=hue.name, dodge=True,
                  hue_order=list(hue_colors.keys()),
                  order=groupx, data=hifs_filtered_melt, ax=ax, 
                  palette=hue_colors, size=2).set(ylabel="", xlabel="")

    annotator=Annotator(ax, pairs, data=hifs_filtered_melt,
                        x=groupx_name, y='value', hue=hue.name,
                        order=groupx)
    annotator.configure(test=None, text_format='simple', fontsize=15) #, pvalue_format='simple')
    annotator.set_custom_annotations(pvalues_corrected)
    annotator.annotate()

    ax.legend(handles=legend_elements, loc='best', fontsize=15, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xticklabels(groupx_labels)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlabel('')
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500)

    return None


def boxplots_pvals_by_groupx_and_groupy(df_pvals: pd.DataFrame, 
                                        groupx_vals: List[str],
                                        groupy_vals: List[str],
                                        groupx_col: str = 'tissue',
                                        groupy_col: str = 'cell',
                                        title: str = '',
                                        groupx_labels: List[str] = None,
                                        groupy_labels: List[str] = None,
                                        boxplot: bool = True,
                                        axvline: float = 0.05,
                                        colors: List[str] = pa_colors,
                                        fname: str = None,
                                        palette: dict = {},
                                        figsize: Tuple = (15, 7),
                                        xlim: List[float] = None,
                                        **kwargs
) -> None:
    """Plots boxplots of pvalues by group x and group y.
    
    Args:
        df_pvals: A dataframe whose rows represent features and 
            columns include pvalues for annotation.
        groupx_vals: A list of strings specifying categories on x-axis.
        groupy_vals: A list of strings specifying categories on y-axis.
        groupx_col: A string specifying column name for groups on x-axis. 
        groupy_col: A string specifying column name for groups on y-axis. 
        title: A string specifying figure title.
            Defaults to ''.
        groupx_labels: A list of strings for setting x-axis label. 
            Defaults to None, in which case `groupx_vals` is used for x-axis label.
        groupy_labels: A list of strings for setting x-axis label. 
            Defaults to None, in which case `groupy_vals` is used for x-axis label.
        boxplot: A bool specifying whether boxplots are plotted.
            Defaults to True.
        axvline: A float specifying x position in data coordinates of the vertical line.
            Defaults to 0.01.
        colors: A list of strings specifying colors. 
        fname: A string for file name/path. 
            Defaults to None, in which case figure is not saved. 

    Returns:
        None
    """
    
    fig, axes = plt.subplots(1, len(groupx_vals), figsize=figsize, sharey=True, sharex=True)
    PROPS = {
        'boxprops':{'edgecolor':'black'},
        'medianprops':{'color':'black'},
        'whiskerprops':{'color':'black'},
        'capprops':{'color':'black'},
        #'flierprops': {"marker": "x"}
    }
    if groupx_labels is None: 
        groupx_labels = groupx_vals 
    if groupy_labels is None: 
        groupy_labels = groupy_vals 

    for idx, groupx_val, color in zip(np.arange(6), groupx_vals, colors):

        df_groupx = df_pvals.loc[df_pvals[groupx_col] == groupx_val, :]
        dict_pvals = {groupy_val: -np.log10(df_groupx.loc[df_groupx[groupy_col] == groupy_val, 'pvalue_corrected']) for groupy_val in groupy_vals}
        sortedTypesDescending = sorted(dict_pvals.items(), key=lambda x: x[1].size)
        sortedTypesDescending.reverse()
        df_groupy = pd.concat([key[1].to_frame(name=key[0]).reset_index(drop=True) for key in sortedTypesDescending], axis=1, join='outer')
        
        if len(palette) == 0:
            if boxplot: 
                sns.boxplot(data=df_groupy, orient='h', ax=axes[idx], linewidth=1, order=groupy_vals, color=color, **PROPS)
            sns.swarmplot(data=df_groupy, orient='h', ax=axes[idx], dodge=True, order=groupy_vals, color=color,  size=4.5, alpha=0.8)
        else:
            if boxplot: 
                sns.boxplot(data=df_groupy, orient='h', ax=axes[idx], linewidth=1, order=groupy_vals, palette=palette, **PROPS)
            sns.swarmplot(data=df_groupy, orient='h', ax=axes[idx], dodge=True, order=groupy_vals, palette=palette,  size=4.5)
        
        for patch in axes[idx].patches:
            r, g, b, a = patch.get_facecolor()
            patch.set_facecolor((r, g, b, .5))

        axes[idx].annotate(text=groupx_labels[idx], xy=(0.5, 1), xycoords='axes fraction', fontsize=15, ha='center', va='bottom')
        axes[idx].axvline(x = -math.log10(axvline), color='tab:red', linestyle='-', linewidth=0.5)
        axes[idx].grid(axis='y')
        axes[idx].tick_params(axis='both', which='major', labelsize=15)
        
        axes[idx].set_yticklabels(groupy_labels)
        if xlim:
            axes[idx].set_xlim(xlim)
    fig.supxlabel('-log10(FDR adjusted p-value)', fontsize=15)
    fig.suptitle(title, fontsize=15)
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500, **kwargs)
        
    return None


def boxplots_pvals_by_group(df_pvals: pd.DataFrame, 
                            group_vals: List[str],
                            group_col: str = 'cell',
                            title: str = '',
                            group_labels: List[str] = None,
                            boxplot: bool = True,
                            color =pa_colors['spacecadet'],
                            fname: str = None
) -> None:
    """Plots boxplots of pvalues for each group."""
    
    fig, ax = plt.subplots(1, 1, figsize=(15,7), sharey=True, sharex=True)
    PROPS = {
        'boxprops':{'edgecolor':'black'},
        'medianprops':{'color':'black'},
        'whiskerprops':{'color':'black'},
        'capprops':{'color':'black'},
    }

    dict_pvals = {group_val: -np.log10(df_pvals.loc[df_pvals[group_col] == group_val, 'pvalue_corrected']) for group_val in group_vals}
    sorted_group = sorted(dict_pvals.items(), key=lambda x: x[1].size)
    sorted_group.reverse()
    df_group = pd.concat([key[1].to_frame(name=key[0]).reset_index(drop=True) for key in sorted_group], axis=1, join='outer')
    
    if boxplot:
        sns.boxplot(data=df_group, orient='h', ax=ax, linewidth=1, order=group_vals, color=color, **PROPS)
    sns.swarmplot(data=df_group, orient='h', ax=ax, dodge=True, color=color, order=group_vals, size=4.5, alpha=0.8)
    
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))

    ax.axvline(x = -math.log10(0.01), color='tab:red', linestyle='-', linewidth=0.5)
    ax.grid(axis='y')
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel('-log10(FDR adjusted p-value)', fontsize=15)
    ax.set_title(title, fontsize=15)
    ax.set_yticklabels(group_labels)
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500)
    
    return None


def plot_stacked_barplots(hifs: pd.DataFrame, 
                          group_labels: pd.Series,  
                          feats: List[str], 
                          title: str ='', 
                          width: float =0.85, 
                          colors = pa_colors, 
                          xticks = False, 
                          figsize = (20, 5), 
                          fname: str = None, 
                          legend: bool = True,
                          legend_labels: List[str] = None,
                          ymax: float = 100,
                          xlab: str = '', 
                          ylab: str = '',
                          **kwargs):
    """Plots stacked barplots."""
    
    (fig, ax) = plt.subplots(1, 1, figsize=figsize)
    tt = pd.DataFrame(columns = feats)
    for group in np.unique(group_labels):
        id_slides = hifs.index[group_labels == group]
        tt = pd.concat([tt, hifs.loc[id_slides, feats].sort_values(by=feats, ascending=False)])
    ind = np.arange(tt.shape[0])
    bottom = np.zeros(len(ind)) 

    for i, col in enumerate(feats):
        ax.bar(ind, tt.loc[:, col], bottom=bottom, width=width, color=colors[i], align='center')
        bottom = bottom + tt.loc[:, col]
    if legend_labels is None:
        labels = [re.search(r'^.*\[\[(?P<stroma>.*?)\]', feat, re.IGNORECASE).group('stroma') for feat in feats]
    else: 
        labels = legend_labels
    if legend:
        ax.legend(labels=labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), frameon=False,
                 fontsize=15, ncol=len(np.unique(feats)))
    ax.set_ylabel(ylab, fontsize=15)
    ax.set_xlabel(xlab, fontsize=15)
    ax.set_ylim(0, ymax)
    ax.set_xlim(0.8, ind[-1] + 0.3)
    ax.tick_params(axis='both', which='major', labelsize=15)
    
    if fname != None:
        plt.savefig(fname, dpi=1200, facecolor='white', bbox_inches='tight', **kwargs)

def violinplots_hifs_by_group(hifs: pd.DataFrame,
                           df_pvals: pd.DataFrame, 
                           regex: str, 
                           hue: pd.Series,
                           hue_colors: dict, 
                           groupx: List[str],
                           fname: str = None,
                           groupx_name: str = 'cell',
                           groupx_labels: List[str] = None,
                           ylabel: str = None,
                           figsize: Tuple[int]=(16, 7),
                           ylim: List[float]=[-0.17, 1.17],
                           yticks=np.arange(0, 1.2, step=0.2),
                           **kwargs
) -> None:
    """Plots violin plot of a given feature by cell type.
    
    Args:
        hifs: A dataframe whose rows represent samples and columns represent features. 
        df_pvals: A dataframe whose rows represent features and 
            columns include pvalues for annotation.
        regex: A string (regular expression) that subsets features for plotting.
        hue: A list of group labels for each row in `hifs` dataframe. 
        hue_colors: A dictionary of colors for each group label. 
        groupx: A list of strings indicating cell types. 
        fname: A string for file name/path. 
            Defaults to None, in which case figure is not saved. 
        groupx_name: A string indicating cell or tissue. 
            Use `cell` for cell or `tissue` for tissue. 
        groupx_labels: A list of strings for setting x-axis label. 
            Defaults to None, in which case `groupx` is used for x-axis label.
        ylabel: A string for setting y-axis label. 
            Defaults to None, in which case `regex` is used for y-axis label. 
        figsize: A tuple of floats specifying width, height in inches.
            Defaults to (16, 7).

    Returns:
        None
    """
    
    if groupx_labels is None: 
        groupx_labels = groupx
    if ylabel is None:
        ylabel = regex 
        
    hifs_filtered = hifs.filter(regex=regex)
    names_hifs = list(hifs_filtered)
    hifs_filtered = hifs_filtered.merge(hue, on='H & E_ID')
    hifs_filtered_melt = pd.melt(hifs_filtered, id_vars=hue.name, 
                                 value_vars=names_hifs)
    if groupx_name == 'cell':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_cell, hifs_filtered_melt.variable))
    elif groupx_name == 'tissue':
        hifs_filtered_melt[groupx_name] = list(map(utils.get_tissue, hifs_filtered_melt.variable))
    groups = list(hue_colors.keys())
    legend_elements = [Line2D([0], [0], color=hue_colors[groups[0]], marker='o', 
                              alpha=0.5, lw=15, label=groups[0]),
                       Line2D([0], [0], color=hue_colors[groups[1]], marker='o', 
                              alpha=0.5, lw=15, label=groups[1])]
    pairs = [((cell, groups[0]), (cell, groups[1])) for cell in groupx] 
    pvalues_corrected = []
    for groupx_val in groupx:
        hif_name = hifs_filtered.columns[hifs_filtered.columns.str.contains(f'\[{groupx_val}')][0]
        pval_tmp = df_pvals.loc[df_pvals.index == hif_name, 'pvalue_corrected'].squeeze()
        
        if pval_tmp < 1e-20:
            pvalues_corrected.append('****') 
        elif pval_tmp < 1e-10:
            pvalues_corrected.append('***')
        elif pval_tmp < 1e-5:
            pvalues_corrected.append('**')
        elif pval_tmp < 0.05:
            pvalues_corrected.append('*')
        elif pval_tmp >= 0.05:
            pvalues_corrected.append('ns')
            
    fig, ax = plt.subplots(figsize=figsize)
    ax = sns.violinplot(x=groupx_name, y='value', hue=hue.name, 
                hue_order=list(hue_colors.keys()),
                data=hifs_filtered_melt, palette=hue_colors,
                order=groupx)
    plt.setp(ax.collections, alpha=.4)

    sns.swarmplot(x=groupx_name, y='value', hue=hue.name, dodge=True,
                  hue_order=list(hue_colors.keys()),
                  order=groupx, data=hifs_filtered_melt, ax=ax, 
                  palette=hue_colors, size=2).set(ylabel="", xlabel="")

    annotator=Annotator(ax, pairs, data=hifs_filtered_melt,
                        x=groupx_name, y='value', hue=hue.name,
                        order=groupx)
    annotator.configure(test=None, text_format='simple', fontsize=15) #, pvalue_format='simple')
    annotator.set_custom_annotations(pvalues_corrected)
    annotator.annotate()

    ax.legend(handles=legend_elements, loc='best', fontsize=15, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xticklabels(groupx_labels)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_ylim(ylim)
    ax.set_xlabel('')
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500, **kwargs)

    return None


#def plot_stacked_barplots(hifs: pd.DataFrame, 
#                          group_labels: pd.Series,  
#                          feats: List[str], 
#                          title: str ='', 
#                          width: float =0.25, 
#                          colors = pa_colors, 
#                          xticks = False, 
#                          figsize = (20, 5), 
#                          fname: str = None, 
#                          legend: bool = True,
#                          legend_labels: List[str] = None,
#                          ymax: float = 1,
#                          xlab: str = '', 
#                          ylab: str = ''):
#    """Plots stacked barplots."""
#    
#    (fig, ax) = plt.subplots(1, 1, figsize=figsize)
#    tt = pd.DataFrame(columns = feats)
#    #print(group_labels)
#    for group in np.unique(group_labels):
#        id_slides = hifs.index[group_labels == group]
#        tt = pd.concat([tt, hifs.loc[id_slides, feats].sort_values(by=feats, ascending=False)])
#    
#    ind = np.arange(tt.shape[0]) + 1 #/ 2
#    bottom = np.zeros(len(ind))
#    legend_labels.append('Other tissue')
#    for i, col in enumerate(feats):
#        ax.bar(ind, tt.loc[:, col], bottom=bottom, width=0.85, color=colors[i], align='center')
#        bottom = bottom + tt.loc[:, col]
#    if legend_labels is None:
#        labels = [re.search(r'^.*\[\[(?P<stroma>.*?)\]', feat, re.IGNORECASE).group('stroma') for feat in feats]
#    else: 
#        labels = legend_labels
#    if legend:
#        ax.legend(labels=labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), frameon=False,
#                 fontsize=15, ncol=len(np.unique(feats)))
#    ax.set_ylabel(ylab, fontsize=15)
#    ax.set_xlabel(xlab, fontsize=15)
#    ax.set_title(title, fontsize=15)
#    ax.set_ylim(0, ymax)
#    ax.set_xlim(0.8, ind[-1] + 0.3)
#    ax.tick_params(axis='both', which='major', labelsize=15)
#    
#    if fname != None:
#        plt.savefig(fname, dpi=1200, facecolor='white', bbox_inches='tight')
#

def plot_overlapping_histograms(hifs: pd.DataFrame,
                                group_label: pd.Series, 
                                list_hifs: List[str],
                                colors: dict,
                                alpha: float = 0.5, 
                                figsize=(12, 10),
                                **kwargs
) -> None:
    """Plot overlapping histograms."""
    
    groups = np.unique(group_label)
    
    n_rows = math.ceil(len(list_hifs)/4)
    fig, axes = plt.subplots(n_rows, 4, figsize=figsize, constrained_layout = True)
    
    # when the subsets have unequal numbers of observations, comparing their distributions in terms of counts may not be ideal.
    for idx, ax in enumerate(axes.flat):
        if idx >= len(list_hifs):
            fig.delaxes(ax)
            continue 
        hif = list_hifs[idx]
        ax.hist(hifs.loc[group_label == groups[0], hif], 
                label=groups[0], alpha=alpha, color=colors[groups[0]], 
                density=True, **kwargs)
        ax.hist(hifs.loc[group_label == groups[1], hif],
                label=groups[1], alpha=alpha, color=colors[groups[1]], 
                density=True, **kwargs)
        ax.legend(loc='upper right')
        ax.set_title(textwrap.fill(hif.split('_HE')[0], 50))
        if idx % 4 == 0:
            ax.set_ylabel('Density')
    
    return None


def boxplots_pvals_by_cell(df_pvals: pd.DataFrame, 
                           cell_types: List[str],
                           title: str = '',
                           cell_types_labels: List[str] = None,
) -> None:
    
    fig, ax = plt.subplots(1, 1, figsize=(15,7), sharey=True, sharex=True)
    PROPS = {
        'boxprops':{'facecolor':'none', 'edgecolor':'black'},
        'medianprops':{'color':'black'},
        'whiskerprops':{'color':'black'},
        'capprops':{'color':'black'},
        #'flierprops': {"marker": "x"}
    }

    dict_pvals = {celltype: -np.log10(df_pvals.loc[df_pvals['cell'] == celltype, 'pvalue_corrected']) for celltype in cell_types}
    sortedCellTypesDescending = sorted(dict_pvals.items(), key=lambda x: x[1].size)
    sortedCellTypesDescending.reverse()
    dfCellType = pd.concat([key[1].to_frame(name=key[0]).reset_index(drop=True) for key in sortedCellTypesDescending], axis=1, join='outer')

    sns.boxplot(data=dfCellType, orient='h', ax=ax, linewidth=1, order=cell_types, **PROPS)
    sns.swarmplot(data=dfCellType, orient='h', ax=ax, dodge=True, color=pa_colors[3], order=cell_types, size=4.5, alpha=0.8)

    ax.axvline(x = -math.log10(0.01), color='tab:red', linestyle='-', linewidth=0.5)
    ax.grid(axis='y')
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel('-log10(FDR corrected p-value)', fontsize=15)
    ax.set_title(title, fontsize=15)
    ax.set_yticklabels(cell_types_labels)
    
    return None


def boxplots_pvals_by_tissue_and_cell(df_pvals: pd.DataFrame, 
                                      tissue_types: List[str],
                                      cell_types: List[str],
                                      title: str = '',
                                      cell_types_labels: List[str] = None,
) -> None:

    fig, ax = plt.subplots(1, len(tissue_types), figsize=(15,7), sharey=True, sharex=True)
    PROPS = {
        'boxprops':{'facecolor':'none', 'edgecolor':'black'},
        'medianprops':{'color':'black'},
        'whiskerprops':{'color':'black'},
        'capprops':{'color':'black'},
        #'flierprops': {"marker": "x"}
    }
    if cell_types_labels is None: 
        cell_types_labels = cell_types

    for idx, tissue_type in zip(np.arange(6), tissue_types):

        df_tissue = df_pvals.loc[df_pvals['tissue'] == tissue_type, :]
        dict_pvals = {celltype: -np.log10(df_tissue.loc[df_tissue['cell'] == celltype, 'pvalue_corrected']) for celltype in cell_types}
        sortedCellTypesDescending = sorted(dict_pvals.items(), key=lambda x: x[1].size)
        sortedCellTypesDescending.reverse()
        dfCellType = pd.concat([key[1].to_frame(name=key[0]).reset_index(drop=True) for key in sortedCellTypesDescending], axis=1, join='outer')

        sns.boxplot(data=dfCellType, orient='h', ax=ax[idx], linewidth=1, order=cell_types, **PROPS)
        sns.swarmplot(data=dfCellType, orient='h', ax=ax[idx], dodge=True, order=cell_types, color=pa_colors[3],  size=4.5, alpha=0.8)

        ax[idx].annotate(text=tissue_type, xy=(0.5, 1), xycoords='axes fraction', fontsize=15, ha='center', va='bottom')
        ax[idx].axvline(x = -math.log10(0.01), color='tab:red', linestyle='-', linewidth=0.5)
        ax[idx].grid(axis='y')
        ax[idx].tick_params(axis='both', which='major', labelsize=15)
        
        ax[idx].set_yticklabels(cell_types_labels)
    fig.supxlabel('-log10(FDR corrected p-value)', fontsize=15)
    fig.suptitle(title, fontsize=15)
        
    return None 


def boxplots_hifs_by_cell(hifs: pd.DataFrame,
                          df_pvals: pd.DataFrame, 
                          group_labels: pd.Series,
                          cell_types: List[str],
                          region: str,
                          color_dict: dict, 
                          fname: str = None,
                          cell_types_labels: List[str] = None,
) -> None:
    
    if cell_types_labels is None: 
        cell_types_labels = cell_types
        
    hifs_filtered = hifs.filter(regex=f'^COUNT PROP.*IN \[{region}\].*')
    names_hifs = list(hifs_filtered)
    hifs_filtered = hifs_filtered.merge(group_labels, on='slideID')
    hifs_filtered_melt = pd.melt(hifs_filtered, id_vars=group_labels.name, 
                                 value_vars=names_hifs)
    hifs_filtered_melt['cell'] = list(map(utils.get_cell, hifs_filtered_melt.variable))
    
    groups = list(color_dict.keys())

    legend_elements = [Line2D([0], [0], color=color_dict[groups[0]], marker='o', 
                              alpha=0.5, lw=15, label=groups[0]),
                       Line2D([0], [0], color=color_dict[groups[1]], marker='o', 
                              alpha=0.5, lw=15, label=groups[1])]
    pairs = [((cell, groups[0]), (cell, groups[1])) for cell in cell_types] 
    pvalues_corrected = []
    for cell in cell_types:
        hif_name = f'COUNT PROP [[{cell} CELLS] OVER [ALL CELLS]] IN [{region}]_HE'
        pval_tmp = df_pvals.loc[df_pvals.index == hif_name, 'pvalue_corrected'].squeeze()
        if pval_tmp < 1e-20:
            pvalues_corrected.append('***') 
        elif pval_tmp < 1e-10:
            pvalues_corrected.append('**')
        elif pval_tmp < 1e-2:
            pvalues_corrected.append('*')
        elif pval_tmp > 0.01:
            pvalues_corrected.append('ns')
            
    fig, ax = plt.subplots(figsize=(16, 7))
    ax = sns.boxplot(x='cell', y='value', hue=group_labels.name, 
                data=hifs_filtered_melt, palette=list(color_dict.values()),
                order=cell_types)
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    sns.swarmplot(x='cell', y='value', hue=group_labels.name, dodge=True,
                  order=cell_types, data=hifs_filtered_melt, ax=ax, 
                  palette=list(color_dict.values()), size=2).set(ylabel="", xlabel="")

    annotator=Annotator(ax, pairs, data=hifs_filtered_melt,
                        x='cell', y='value', hue=group_labels.name)
    annotator.configure(test=None, text_format='simple') #, pvalue_format='simple')
    annotator.set_custom_annotations(pvalues_corrected)
    annotator.annotate()

    ax.legend(handles=legend_elements, loc='best', fontsize=15, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xticklabels(cell_types_labels)
    ax.set_ylabel(f'Proportion in {region.lower()}', fontsize=15)
    ax.set_xlabel('')
    
    if fname is not None:
        fig.savefig(fname, bbox_inches='tight',
                    dpi=500)

    return None 


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html
    
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
