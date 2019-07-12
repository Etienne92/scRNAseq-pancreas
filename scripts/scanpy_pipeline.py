import scanpy as sc 
import bbknn
import numpy as np
from sklearn import mixture
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Takes a raw anndata (not filtered/normalized etc) and applies the scanpy pipeline.')
parser.add_argument("-i",type=str,help="Input file")
parser.add_argument("-o",type=str,help="Output file")
parser.add_argument("-datasets",type=str,default="all",help="Datasets to use")

args=parser.parse_args()



def detect_doublets(adata,marker_genes=["GCG","INS","SST","PPY","COL3A1","CFTR","PRSS2","GHRL"],inplace=True):
    """Filter out cells which express more than 1 marker gene. The threshold between expressed and not expressed is defined by fitting a gaussian mixture model."""
    counts=np.zeros((1,adata.shape[0]))
    for gene in marker_genes:
        gm = mixture.GaussianMixture(n_components=2, covariance_type='full',reg_covar=0.3)
        expressions = (adata[:,gene].X).reshape(-1,1)
        gm.fit(expressions)
        predictions = gm.predict(expressions)
        if gm.predict([[0]]):
            predictions = 1 - predictions
        counts= counts + predictions
    if inplace:
        adata._inplace_subset_obs((counts <=1)[0])
    else: 
        #In that case, the doublets won't be removed, but the "doublet score" will be added to the anndata. This is useful for testing that this filter correctly identifies the doublets.
        adata.obs["doublets"] = counts[0] 

adata = sc.read(args.i)
if args.datasets !="all":
    datasets = args.datasets.split(",")
    adata._inplace_subset_obs(adata.obs["dataset"].isin(datasets))
sc.pp.filter_cells(adata,min_counts = 80000)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
sc.pp.filter_cells(adata,min_genes=3700)
sc.pp.filter_genes(adata,min_cells=3)

sc.pp.normalize_total(adata,target_sum = 10000)

#Filter out cells which have a high percentage of mitochondrial genes
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

adata._inplace_subset_obs(adata.obs['percent_mito'] < 0.25)

        
sc.pp.log1p(adata)
detect_doublets(adata,inplace=True)

sc.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean = 10000000,min_disp=1,max_disp=1000000,flavor="cell_ranger")
adata.var["highly_variable"] = (adata.var["dispersions"] > 1) & (adata.var["means"]>0.0125) #use dispersions and not normalized dispersions

sc.tl.pca(adata,n_comps = 20)
bbknn.bbknn(adata,batch_key="dataset",n_pcs=20)
sc.tl.louvain(adata,resolution=0.1,key_added="louvain0.1")
sc.tl.louvain(adata,resolution=0.3,key_added="louvain0.3")
sc.tl.louvain(adata,resolution=0.5,key_added="louvain0.5")
sc.tl.louvain(adata,resolution=0.7,key_added="louvain0.7")
sc.tl.louvain(adata,resolution=1.0,key_added="louvain1.0")
sc.tl.louvain(adata,resolution=2.0,key_added="louvain2.0")

sc.tl.umap(adata)
sc.tl.tsne(adata,n_pcs=20)

adata.write(args.o)

