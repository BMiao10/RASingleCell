"""
### Adapted from: M.D. Luecken, F.J. Theis, "Current best practices in single-cell RNA-seq analysis: a tutorial", Molecular Systems Biology 15(6) (2019): e8746

https://github.com/theislab/single-cell-tutorial
"""

############ LIBRARIES ############
# single cell packages
import pandas as pd
import anndata
import scanpy as sc
import harmonypy as hm

# plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mpl_colors
from numpy import linspace
import seaborn as sns

# stats & math
import numpy as np

# settings
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 2 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white')

############ SEURAT I/O ############
def loadFromSeurat(mtx_path, genes_path, barcodes_path, metadata_path, output=None):
    
    """
    Converts R seurat/singlecellexperiment data into gene x barcodes anndata object
    
    Args:
        mtx_path (str): path to .mtx file for count matrix
        genes_path (str): path to .csv file for gene/row labels
        barcodes_path (str): path to .csv with barcodes/column labels
        metadata_path (str): path to .csv file with metadata
        output (str): Path to write anndata object to or None if anndata should be returned (default=None)
    
    Returns:
        anndata.Anndata
    """
    
    #loading in seurat data matrix into anndata
    _adata = sc.read(mtx_path, cache=False)
    _adata = _adata.transpose()

    #loading in adata index/column headers
    _genes = pd.read_csv(genes_path, header=0)
    _barcodes = pd.read_csv(barcodes_path, header=0)
    
    _genes, barcodes = _genes[['x']], barcodes[['x']]
    _genes.columns, barcodes.columns = ['genes'], ['barcodes']
    _genes = _genes.set_index('genes',drop=True)
    _barcodes = _barcodes.set_index('barcodes',drop=True)

    #loading metadata to be added to adata.obs
    _metadata = pd.read_csv(metadata_path, header=0)
    _metadata.columns = ['barcodes' if 'Unnamed' in col else col for col in _metadata.columns]
    _metadata = _metadata.set_index('barcodes', drop=True)
    _barcodes = _barcodes.join(_metadata)

    #setting adata var/obs
    _adata.var = genes
    _adata.obs = barcodes
    
    if type(output) == str: _adata.write(output)

    return _adata


############ PREPROCESSING ############
def labelSampleType(adata):
    """
    Adds "sampleType", "sex", "race" columns to anndata with samples labeled as RA, Control, or AS
    
    Args:
        adata (Dataframe): gene x barcodes anndata or pandas dataframe with "Sample" label metadata
    
    Returns:
        Dataframe
        
    """
    
    #Define a variable that stores sample type labels
    AS_pts = ['RA11', 'RA12', 'RA13', 'RA14']
    
    if isinstance(adata, pd.DataFrame):
        adata['sampleType'] = ['AS' if any(substring in sample for substring in AS_pts) 
                                   else 'RA' if 'RA' in sample 
                                   else 'Control' for sample in adata['Sample']]
    else: 
        adata.obs['sampleType'] = ['AS' if any(substring in sample for substring in AS_pts) 
                                   else 'RA' if 'RA' in sample 
                                   else 'Control' for sample in adata.obs['Sample']]
    
    # manually annotate sex & ethnicity
    adata.obs['sex'] = [sample.split('_')[-1] for sample in adata.obs['Sample']]
    adata.obs['race'] = [sample.split('_')[2] for sample in adata.obs['Sample']]

    return adata

def cleanRAData(adata, plot=False, remove_batches=[]):
    """
    Light preprocessing **specific to this RA dataset**
    
    Args:
        adata (anndata.Anndata): gene x barcode anndata containing raw read counts 
        plot (bool): whether show metric plots (default=False)
        remove_batches (list): batches to remove (default=[])
        singlets_only (bool): bool to remove non-singlet cells (default=True)
      
    Returns:
        anndata.Anndata
    """

    if type(adata) == str: adata = sc.read(adata, cache=False)
    
    # annotate batches
    adata.obs['batch'] = [item.split('-')[0].split('run')[-1] for item in adata.obs['Lane']] 

    #remove any '.' from obs names
    adata.obs.columns = [item.replace('.', '_') for item in adata.obs.columns]
    
    # remove nan samples
    adata = adata[adata.obs["Sample"]!="nan",:]

    # remove batches from adata
    for batch in remove_batches:
        adata = adata[~(adata.obs['batch'] == batch)]
        
    # add sample labels
    adata = labelSampleType(adata)
        
    return adata

############ SCANPY METRICS ############   
def plotQCMetrics(adata, groupby="batch", out_path=None):
    """
    Plotting standard QC plots
    
    Args:
        adata (anndata.AnnData): object containing single cell value counts
        groupby (str): column to group by when showing metrics
        out_path (str): folder to save QC metrics to 
    
    Returns:
        anndata.AnnData
    
    """
    
    # add QC metrics for mitochondrial and ribosomal genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', "ribo", "hb"], percent_top=None, log1p=False, inplace=True)
    
    # value counts for figure sizing
    unique_size = adata.obs[groupby].nunique()
    
    # plot QC metrics
    show = True if out_path is None else False
    sc.set_figure_params(scanpy=False,figsize=(unique_size/3, 4))
    ax = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
             jitter=0.1, groupby = groupby, scale="width", ha='center', 
                      rotation_mode='anchor', rotation=90,cut=0, show=show,
                     **{"wspace":0})
    
    if out_path is not None:
        ax[0].figure.savefig(out_path+"/qc_overall_violinplots.pdf", bbox_inches="tight")
    
    sc.set_figure_params(scanpy=False,figsize=(unique_size/3, 4))
    ax = sc.pl.violin(adata, ['pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'],
             jitter=0.1, groupby = groupby, scale="width", rotation= 90,
                      ha='center',  rotation_mode='anchor', cut=0, show=show, **{"wspace":0})
    
    if out_path is not None:
        ax[0].figure.savefig(out_path+"/qc_mt-ribo-hb-counts_violinplots.pdf", bbox_inches="tight")
        
    #Data quality summary plots
    fig, axs = plt.subplots(ncols =3, figsize = (15,4), sharey= False, sharex=False);
    sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_mt', show=show, ax = axs[0])
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color=groupby, show=show, ax = axs[1],);
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color=groupby, show=show, ax = axs[2]);
    
    if out_path is not None:
        plt.savefig(out_path+"/genes_by_total_counts_scatterplot.pdf", bbox_inches="tight")
    else:
        plt.show()
        
    # histplots
    fig, axs = plt.subplots(ncols =3, figsize = (20,4), sharey= False);
    sns.histplot(adata.obs['total_counts'], kde = False, bins=60, ax = axs[0]);
    sns.histplot(adata.obs['n_genes_by_counts'], kde = False, bins=60, ax = axs[1]);
    sns.histplot(adata.obs['pct_counts_mt'], kde = False, bins=80, ax = axs[2]);
    
    if out_path is not None:
        plt.savefig(out_path+"/qc_histplots.pdf", bbox_inches="tight")
    else:
        plt.show()
    
    return adata

    
def set_filters(adata, min_cells=3, min_genes=100, max_genes=None, min_count=None, max_count=None, 
               mtpct=25, rbpct=None, remove_mt=True, remove_rb=True, remove_hb=True, remove_genes=["MALAT1"]):
    
    """
    Filters cells by number of gene/read counts and remove mitochondrial/ribosomal genes
    
    Args:
        adata (anndata.Anndata): anndata containing qc metrics
        min_cells (int): min num cells to keep gene (default=3)
        min_genes (int): min num genes to keep cell (default=100)
        min_count (int): min read count to keep gene 
        max_count (int): max read count to keep gene 
        mtpct (int): max percent of mitochondrial genes to keep cell
        rbpct (int): max percent of ribosomal genes to keep cell
        remove_mt (bool): removing mitochondiral genes
        remove_rb (bool): removing ribosomal genes 
        remove_hb (bool): removing hemoglobin genes
        remove_genes (list<str>): custom list of genes to remove
    
    Returns:
        anndata.Anndata
    
    """
    
    # filter out genes found in less than min_cells
    print('Total number of genes: {:d}'.format(adata.n_vars))
    if isinstance(min_cells, int): 
        sc.pp.filter_genes(adata, min_cells = min_cells)
        print('Genes after min cells filter: {:d}'.format(adata.n_vars))

    # filter out cells with < min_genes
    print()
    if isinstance(min_genes, int): 
        sc.pp.filter_cells(adata, min_genes= min_genes)
        print('Cells after min genes filter: {:d}'.format(adata.n_obs))
    
    # filter out cells with > max_genes
    print()
    if isinstance(max_genes, int): 
        sc.pp.filter_cells(adata, max_genes= max_genes)
        print('Cells after max genes filter: {:d}'.format(adata.n_obs))

    # filter out cells with < min_count reads
    print('\nTotal number of cells: {:d}'.format(adata.n_obs))
    if isinstance(min_count, int): 
        sc.pp.filter_cells(adata, min_counts = min_count)
        print('Cells after min count filter: {:d}'.format(adata.n_obs))

    # filter out cells with > max_count reads
    print()
    if isinstance(max_count, int): 
        sc.pp.filter_cells(adata, max_counts = max_count)
        print('Cells after max count filter: {:d}'.format(adata.n_obs))
    
    # filter out cells with mt % > mtpct
    if mtpct is not None: 
        adata = adata[adata.obs.pct_counts_mt < mtpct, :]
        print('\nCells after mito filter: {:d}'.format(adata.n_obs))

    # filter out cells with ribo % < rbpct
    if rbpct is not None: 
        adata = adata[adata.obs[ 'pct_counts_ribo'] > rbpct, :]
        #adata = adata[adata.obs[ 'pct_counts_ribo'] < 60, :]
        print('Cells after ribo filter: {:d}'.format(adata.n_obs))
        
    # filter out mitochondrial genes
    if remove_mt: adata = adata[:,~adata.var_names.str.startswith('MT-')]
        
    # filter out ribosomal genes
    if remove_rb: 
        adata = adata[:,~adata.var_names.str.startswith('RPL')]
        adata = adata[:,~adata.var_names.str.startswith('RPS')]
    
    # filter out hemoglobin genes
    if remove_hb: adata = adata[:,~adata.var_names.str.contains('^HB[^(P)]')]
    
    # filter out unwanted gene
    if len(remove_genes) > 0: 
        adata = adata[:,~adata.var_names.isin(remove_genes)]
    
    return adata

def printMetrics(adata, metric_list):
    
    """
    Print out count value metrics from anndata
    
    Args: 
        adata : anndata object
        metricList : list of metrics to print value counts for
    
    Args:
        prints out metrics to console
    """

    # print the total size of the data set
    print('cells x genes = ' + str(adata.shape) + '\n')

    # print some metrics about the dataset
    for _metric in metric_list:
        if _metric in adata.obs.columns:
                print(adata.obs[_metric].value_counts())
                print('')
        else: print(_metric + ' not found')


############ DIMENSIONALITY REDUCTION ############   
def batchCorrect(adata, output_path=None, max_iter_h=3, epsilon=-np.inf,
                 PCS="X_pca", vars_use = ['batch', 'Lane'], **kwargs):
    """
    Batch correction using harmonypy package
    
    Args:
        adata (anndata.Anndata): AnnData with PCA calculated and stored in adata.obsm
        output_path (str, pathlike): File to save corrected PCs to, use None to return
        max_iter (int): Number of iterations to perform. If None, runs default harmonypy
        PCS (str): Key in obsm containing PCs to correct
        vars_use (list<str>): List of strings to use in harmonypy analysis
        
    Returns:
        pandas.Dataframe
        
    """
        
    # extract PCA components
    metadata = adata.obs
    data_mat = adata.obsm[PCS]

    # Run Harmony
    if max_iter_h == None: ho = hm.run_harmony(data_mat, metadata, vars_use) 
    else: ho = hm.run_harmony(data_mat, metadata, vars_use, 
                              max_iter_harmony = max_iter_h, epsilon_harmony=epsilon, 
                              **kwargs) 

    # extract corrected pcs
    res = pd.DataFrame(ho.Z_corr)

    # update PCA embeddings
    adata.obsm['X_pca'] = res.T.values

    # plot pca
    sc.pl.pca(adata, color=vars_use)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)
    
    # return batch corrected or save to file
    if output_path is None: 
        return adata.obsm['X_pca']
    else: 
        cols = ["PC"+str(i+1) for i in range(adata.obsm['X_pca'].shape[1])]
        corr_PCA = pd.DataFrame(adata.obsm['X_pca'], columns = cols, index=adata.obs_names)
        corr_PCA.to_csv(output_path)



def calculateUMAP(adata, flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5, hvg_batch='batch',
                 n_neighbors=10, n_pcs=60, return_full=False, umap=False):
    """
    Processes anndata for UMAP calculation with HVG segmentation, scaling, PCA, bbknn, and umap
    
    Args:
        adata (anndata.AnnData): 
        flavor (str): 
        min_mean (float): 
        max_mean (int): 
        min_disp (float): 
        hvg_batch (str): 
        neighbors (int): 
        bbknn_batch (str): 
        n_pcs (int): 
        bbknn_pcs (int): 
        return_full (bool): 
    
    """
    
    if return_full: 
        adata_all = adata.copy()
        adata_all.raw = adata_all
    
    # select HVG
    sc.pp.highly_variable_genes(adata, flavor=flavor, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, batch_key=hvg_batch)
    sc.pl.highly_variable_genes(adata)

    adata = adata[:, adata.var.highly_variable]
    print("Highly variable genes: %s"%len(adata.var_names))
    
    # regress out unwanted variables
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

    # scale data
    sc.pp.scale(adata, zero_center = True, max_value=10)

    ### dimensional reduction and plotting
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs = n_pcs)
    
    # return the full data
    if return_full: 
        adata = adata_all.raw.toarray()

    # calculate neighbors and UMAP embeddings
    if umap:
        sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
        sc.tl.umap(adata, min_dist=0.1, spread=1, )
        sc.pl.umap(adata, color=['sampleType', 'batch', 'Lane'], wspace=0.35)
        
    return adata
            
def cluster_small_multiples(adata, clust_key, size=10, frameon=False, legend_loc=None, **kwargs):
    """
    Plot subset of umap 
    Adapted from: https://github.com/theislab/scanpy/issues/955
    
    Args:
        adata (anndata.AnnData): anndata object 
        clust_key (str): adata.obs column containing categories to plot
        size (int): 
        frameon (bool): frame around subplots
        legend_loc (str, None): location of legend
        **kwargs: passed to sc.pl.umap

    """
    # load embeddings
    sc.pl.umap(adata, color=clust_key, show=False,return_fig=True)
    
    # split colors
    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        adata.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        adata.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(adata, groups=adata.obs[clust].cat.categories[1:].values, 
               color=adata.obs[clust_key].cat.categories.tolist(), 
               size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)
    
    
    
    
    


