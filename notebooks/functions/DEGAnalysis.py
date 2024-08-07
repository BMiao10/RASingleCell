############ LIBRARIES ############
# single cell packages
import pandas as pd
import anndata
import scanpy as sc

# plotting
import matplotlib.pyplot as plt
import seaborn as sns

# stats and math
from statannot import add_stat_annotation
import itertools
from itertools import product
import statsmodels.stats as sm
import scipy
import numpy as np

############ UNPAIRED ANALYSIS ############
def getCellSubsetProportions(adata, cluster_col = "singleR_blueprint", sample_col="Sample"):
    
    '''
    Counts number of each cell type in each sample and provides
    
    Parameters
    ----------
        adata (anndata.AnnData): adata
        cluster_col (str): col in adata.obs containing clusters to count cells in
        sample_col (str): sample col in adata.obs to count cells for
        
    Output
    ----------
        Returns dataframe containing cell counts
    
    '''
    
    # get annotated adata 
    if cluster_col not in adata.obs.columns: print("cell annotations not found")
    
    # segment metadata
    sample_counts = adata.obs[[sample_col,cluster_col]]
    
    # add cell counts column
    sample_counts['cell_count'] = 1
    
    # get cell counts for each label
    sample_counts = sample_counts.groupby([sample_col,cluster_col], as_index=False).count()
    sample_counts = sample_counts.replace(np.nan, 0)
    
    # get total cell counts per sample
    total_counts = adata.obs[[sample_col]]
    total_counts["count"] = 1
    total_counts = total_counts.groupby([sample_col], as_index=False).count()
    total_counts = dict(zip(total_counts[sample_col],total_counts["count"]))
    
    # map total cell counts to per sample count data frame
    sample_counts["total"] = [total_counts[cell] for cell in sample_counts[sample_col]]
    sample_counts = sample_counts[sample_counts["total"] != 0]
    sample_counts["proportion"] = 100*(sample_counts["cell_count"] / sample_counts["total"])
  
    return sample_counts

def boxplotByCellType(cell_count_df, hue=None, palette=None, height=10, width=16, x_type='singleR_blueprint', 
                      y_type='cell_count', offset=True, order=None, hue_order=None,
                      x_label='', y_label='Proportion of total cells'):
    '''
    Box plots comparing cell counts between groups
    Can run getCellSubsetProportions() to get cell_count_df
    
    Parameters
    ----------
        cell_count_df (pd.DataFrame): dataframe containing cell counts (y_type) for each group (x_type)
        x_type (str): column name in cell_count_df containing cell type labels
        y_type (str): column name in cell_count_df containing cell counts
        x_label, y_label (str): axis object labels
        offset (bool): offset x tick label if label length > 14 (default = True)
        order (list): list containing order of x ticks
        hue : column name in cell_count_df to group samples by for comparison
        (deprecated) palette
        height, width : passed to sns.boxplot
    
    Output
    ----------
        Box plots comparing cell counts between groups
    
    '''
    # create 
    fig, ax = plt.subplots(figsize=(width, height))
    
    """# update palette
    if palette is None: 
        if hue == "sex":  palette = {"M": "#9cf6fd", "F": "#fed7fe"}
        elif hue == "sampleType":  
            print("here") 
            palette = {"Control": "#9cf6fd", "RA": "#fed7fe"}"""
    
    # plot
    sns.stripplot(data=cell_count_df, x=x_type, y=y_type, hue=hue, jitter=True, hue_order=hue_order,
                  linewidth=0.5,edgecolor='gray', size=5, palette="dark:.25", dodge=True) #split=True,

    # Get the ax object to use later.
    ax = sns.boxplot(data=cell_count_df, x=x_type, y=y_type, hue=hue, fliersize=0, width=0.8,
                     order=order, linewidth=1.2, saturation=0.7, hue_order=hue_order)
    
    # Get the handles and labels
    handles, labels = ax.get_legend_handles_labels()

    # When creating the legend, only use the first two elements
    # to effectively remove the last two.
    n = len(cell_count_df[hue].unique())
    l = plt.legend(handles[n:], labels[n:], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
    # set params
    ax.set(xlabel=x_label, ylabel=y_label)
    ax.tick_params(labelsize=15)
    
    # offset tick labels
    
    if offset: 
        if isinstance(order, list): tick_labels = order
        else: tick_labels = [label.get_text() for label in ax.get_xticklabels()]
        
        ax.set_xticklabels(["\n" + label if i%2==0 else label for i, label in  enumerate(tick_labels)])
    
    return ax

def unpairedTestByCellType(cell_type_df, subset_col="PBMC_subset", sample_col="Sample", 
                      groupby="dx", groups = ["control", "RA"], 
                      value_col="proportion", stats_test=scipy.stats.kruskal, corr="s", **kwargs):
    
    """
    Args:
        cell_type_df(pd.DataFrame): dataframe containing cell counts/proportions, cell subset labels, and metadata. If paired, must have index with (Sample, pair)
        subset_col (str): column containing cell subsets
        groupby (str): column containing groups to compare between
        groupA, groupB (str): regex to select groups to compare between
        sample_col (str): column containing samples IDs 
        pair_index (str): column containing pair index for each sample 
        stats_test (test): what scipy test to apply
        value_col (str): column with values to compare between groups
        **kwargs: passed to stats_test method
        
    Returns:
        pd.DataFrame()
    """
    pvals = {}
    
    for s in cell_type_df[subset_col].unique():
        # subset to cell type
        curr_df = cell_type_df[cell_type_df[subset_col] == s]
        values = []
        
        for n in groups:
            # subset to group within cell type
            curr_df_sub = curr_df[curr_df[groupby].str.contains(n)]
            
            # get values
            values.append(list(curr_df_sub[value_col]))
            
        if len(values) > 2:
            pvals[s] = stats_test(*values, **kwargs).pvalue
        else:
            pvals[s] = stats_test(values[0], values[1], **kwargs).pvalue

    if corr is not None:
        padj = sm.multitest.multipletests(list(pvals.values()), method=corr)[1] # get adjusted pvals
        pvals = dict(zip(pvals.keys(), padj))
        
    return pvals

def _boxPairs(groupA, groupB):
    """
    Creates pairs of samples to add statistical annotations for 
    
    Args:
        groupA, groupB (list): list of items to create pairs of
        
    Returns:
        list<tuples>
        
    Example:
        boxPairs(groupA = ["Classical monocytes"], groupB = ["Non-remission", "Remission"])
        >> [(('Classical monocytes', 'Non-remission'), ('Classical monocytes', 'Remission'))]
    """
    groups = []
    box_pairs = []

    # create initial groups
    for a in groupA:
        groups.append([x for x in product([a], groupB)])
    
    # create pairs of groups, ignoring duplicates
    for curr in groups:
        for c in curr:
            curr.remove(c) 
            box_pairs = box_pairs + [x for x in product([c], curr)]

    return box_pairs

def compositionalAnalysis(adata, hue_groups, hue="cohort", subset_col="PBMC", 
                          sample_col="Sample", height=6,width=14, single_pair=False,
                              stats_test=scipy.stats.mannwhitneyu, stat_annot_test='Mann-Whitney',
                          corr="fdr_bh", density=False, return_df=True, filename=None,min_cells=1500,
                          **kwargs):
    '''
    Calculates relative density of group in umap embeddings
    
    Args:
        adata (anndata.AnnData): contains raw or normalized counts
        hue_groups (list): hue types
        hue (str): adata.obs column with groups to compare statistically
        subset_col (str): adata.obs column with cell types
        stats_test (test): scipy stats test to apply
        min_cells (int): min cells per subset to keep
        stat_annot_test (str): if using multiple hues, pass this to stat_annot
        single_pair (bool): only allow tests to be compared between 2 hue groups
        **kwargs
    
    Returns:
        anndata.AnnData
    
    '''

    # compositional analysis
    plt.rcParams["figure.figsize"]=(6,6)
    if density:
        sc.tl.embedding_density(adata, basis='umap', groupby=hue)
        sc.pl.embedding_density(adata, basis='umap', key="umap_density_%s"%hue, group=hue_groups)

    # plot cell counts
    cell_type_df = getCellSubsetProportions(adata, cluster_col=subset_col)
    cell_type_df[hue] = cell_type_df[sample_col].map(dict(zip(adata.obs[sample_col], adata.obs[hue])))
    cell_type_df = cell_type_df[cell_type_df[hue] != 'nan']
    cell_type_df = cell_type_df[~cell_type_df[hue].isna()]
    
    # remove subsets that with less than x# of cells
    cells_to_keep = cell_type_df.groupby(subset_col).sum()
    cells_to_keep = list(cells_to_keep[cells_to_keep["cell_count"]>=min_cells].index)
    cell_type_df = cell_type_df[cell_type_df[subset_col].isin(cells_to_keep)]
    
    # boxplot with annotations
    cell_type_df = cell_type_df[cell_type_df[hue].isin(hue_groups)]
    cell_type_df[subset_col] = pd.Categorical(cell_type_df[subset_col], ordered=True, categories=cell_type_df[subset_col].unique())
    ax = boxplotByCellType(cell_type_df, x_type=subset_col, hue=hue, y_type="proportion", 
                           height=height,width=width, hue_order=hue_groups)
    
    #cell_type_df["stats_pair"] = cell_type_df[sample_col].map(dict(zip(adata.obs[sample_col], adata.obs["stats_pair"])))
    cell_type_df["sample_id"] = [s.split("_")[0][-2:] for s in cell_type_df[sample_col]]
    
    # get pvals
    pvals = unpairedTestByCellType(cell_type_df, subset_col=subset_col, groupby=hue, 
                                   groups = cell_type_df[hue].unique(), 
                                   stats_test=stats_test, corr=corr, **kwargs)
    
    if single_pair: # len(hue_groups) > 2:
        box_pairs = _boxPairs(cell_type_df[subset_col].unique(), [hue_groups[0], hue_groups[-1]])
        
        add_stat_annotation(ax, data=cell_type_df, x=subset_col, y="proportion", hue=hue,  pvalues=list(pvals.values()), 
                            box_pairs=box_pairs, hue_order=hue_groups, perform_stat_test=False,  text_format='star', verbose=0) 
    else:
        box_pairs = _boxPairs(cell_type_df[subset_col].unique(), hue_groups)
        
        add_stat_annotation(ax, data=cell_type_df, x=subset_col, y="proportion",
                            hue = hue, box_pairs=box_pairs, comparisons_correction=corr,
                        test=stat_annot_test, text_format='star', loc='inside', verbose=2)
    
        
    if filename is not None: plt.savefig(filename)
    
    if return_df: return cell_type_df
    elif return_df=="ax": return ax
    elif return_df=="fig": return fig
    
    return pvals

############ _idED ANALYSIS ############
def pairedTestByCellType(cell_type_df, subset_col="PBMC_subset", groupby="dx", 
                 groupA="RA", groupB="control", sample_col="Sample", 
                 pair_index="stats_pair", compare_col="proportion", 
                 stats_test=scipy.stats.wilcoxon, corr="fdr_bh", **kwargs):
    """
    Args:
        cell_type_df(pd.DataFrame): dataframe containing cell counts/proportions, cell subset labels, and metadata. If paired, must have index with (Sample, pair)
        subset_col (str): column containing cell subsets
        groupby (str): column containing groups to compare between
        groupA, groupB (str): regex to select groups to compare between
        sample_col (str): column containing samples IDs 
        pair_index (str): column containing pair index for each sample 
        stats_test (test): scipy stats test to apply
        compare_col (str): column with values to compare between groups
        **kwargs: passed to stats_test method
        
    Returns:
        pd.DataFrame()
    """
    
    # get "control" samples
    A_df = cell_type_df[cell_type_df[groupby].str.contains(groupA)]
    A_df.columns = [groupA+"_"+s for s in A_df.columns]

    # get "case" samples
    B_df = cell_type_df[cell_type_df[groupby].str.contains(groupB)]
    B_df.columns = [groupB+"_"+s for s in B_df.columns]
    
    # perform statistical testing
    pvals = {}
    for s in cell_type_df[subset_col].unique():
        curr_A_df = A_df[A_df[groupA+"_"+subset_col] == s]
        curr_B_df = B_df[B_df[groupB+"_"+subset_col] == s]
            
        curr_A_df = curr_A_df.set_index(groupA+"_"+pair_index)
        curr_B_df = curr_B_df.set_index(groupB+"_"+pair_index)
            
        curr_df = curr_A_df.merge(curr_B_df, right_index=True, left_index=True, how="inner")
        
        pvals[s] =  stats_test(curr_df[groupA+"_"+compare_col],
                               curr_df[groupB+"_"+compare_col],  **kwargs).pvalue 
    
    if corr is not None:
        padj = sm.multitest.multipletests(list(pvals.values()), method=corr)[1] # get adjusted pvals
        pvals = dict(zip(pvals.keys(), padj))
        
    return pvals




