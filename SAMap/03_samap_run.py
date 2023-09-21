#  srun --mem=250G --ntasks=1 --cpus-per-task=1 --time=96:0:0 --pty bash -l 
# conda activate sam-amp_v3
# pip install samap

from samap.mapping import SAMAP
from samap.utils import save_samap
from samap.utils import load_samap
from samap.analysis import (get_mapping_scores, GenePairFinder, transfer_annotations,sankey_plot, chord_plot, CellTypeTriangles,ParalogSubstitutions, FunctionalEnrichment,convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd


fn1 = '/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/sce_mm_myeloidOnly_counts_forSAMap_filtered_pr.h5ad'
fn2 = '/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/sce_zf_myeloidOnly_counts_forSAMap_filtered_pr.h5ad'

#filenames = {'mm':fn1,'dr':fn2}

sam1=SAM()
sam1.load_data(fn1)

sam2=SAM()
sam2.load_data(fn2)

sam1.adata
sam2.adata

n1 = pd.read_csv("/athena/abc/scratch/paz2005/samap_directory/mm.names.tsv", sep='\t')
n1 = [tuple(row) for row in n1.values]


n2 = pd.read_csv("/athena/abc/scratch/paz2005/samap_directory/dr.names.tsv", sep='\t')
n2 = [tuple(row) for row in n2.values]

names = {'mm':n1,'dr':n2}

sams = {'mm':sam1,'dr':sam2}

# run SAMap
sm = SAMAP(sams, f_maps = '/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/maps/', names=names, save_processed=False)

save_samap(sm, '/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/samap_out/tmp.h5')

sm.run(pairwise=True)
save_samap(sm, '/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/samap_out/samap_after_pairwise.h5')

sm.samap.adata.write("/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/samap_out/samap_adata.h5")

p = sm.scatter()
import matplotlib.pyplot as plt
plt.savefig('/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/samap_out/samap_after_pairwise_scatter.png')




####
## cell type mapping scores
sm = load_samap('/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/samap_out/samap_after_pairwise.h5')

keys = {'mm':'clst.MonoM0Only_k200','dr':'M0_k50_cl9.15k50'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
D.head()
D.to_csv('/athena/abc/scratch/paz2005/projects/czi_samap_myeloid/samap_out/D_mappingTable_cellLabels.csv')