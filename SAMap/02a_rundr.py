# conda create -n sam-amp_v3 python=3.7
# conda activate sam-amp_v3 
# pip install sam-algorithm
# pip install hnswlib
#  pip install igraph
#  pip install leidenalg
#cd /athena/abc/scratch/paz2005/projects/czi_samap_myeloid

#conda activate sam-amp_v3 

#######.  Running first steps outside of the main SAM function
from samalg import SAM
import pandas as pd

# https://github.com/atarashansky/SAMap/blob/main/samap/mapping.py
def prepare_SAMap_loadings(sam, npcs=300):
    """ Prepares SAM object to contain the proper PC loadings associated with its manifold.
    Deposits the loadings in `sam.adata.varm['PCs_SAMap']`.
    
    Parameters
    ----------    
    sam - SAM object
    
    npcs - int, optional, default 300
        The number of PCs to calculate loadings for.
    
    """
    ra = sam.adata.uns["run_args"]
    preprocessing = ra.get("preprocessing", "StandardScaler")
    weight_PCs = ra.get("weight_PCs", False)
    A, _ = sam.calculate_nnm(
        n_genes=sam.adata.shape[1],
        preprocessing=preprocessing,
        npcs=npcs,
        weight_PCs=weight_PCs,
        sparse_pca=False,
        update_manifold=False,
        weight_mode='dispersion'
    )
    sam.adata.varm["PCs_SAMap"] = A


data="/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/sce_zf_myeloidOnly_counts_forSAMap_filtered.h5ad"
sam = SAM()
sam.load_data(data)
sam.preprocess_data(sum_norm="cell_median",norm="log",thresh_low=0.0,thresh_high=0.96,min_expression=1)
sam.run(preprocessing="StandardScaler",npcs=100,weight_PCs=False,k=20,n_genes=3000,weight_mode='rms')
sam.leiden_clustering(res=3)
prepare_SAMap_loadings(sam) 
sam.save_anndata(data.split('.h5ad')[0]+'_pr.h5ad')

