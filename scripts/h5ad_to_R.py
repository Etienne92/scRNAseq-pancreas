import scanpy as sc 
import pandas as pd
import numpy as np
import argparse



parser = argparse.ArgumentParser(description='Converts a raw anndata (containing the read counts, without normalization) and returns csv files to be used by DESeq2 and edgeR.')
parser.add_argument("-celltype",type=str,default= "beta cell",help="Only returns the cells which have this inferred cell type.")
parser.add_argument("-datasets",type=str,help="Datasets to use")
parser.add_argument("-i",type=str,help="Input file")
parser.add_argument('-pseudobulk', type=str, default = "True", help="If true, do pseudobulk (aggregate the expressions by individual).")
parser.add_argument('-normalize',default=False,dest='normalize',action='store_true')
args=parser.parse_args()

datasets = args.datasets.split(",")
datasets = [x.strip() for x in datasets]
do_pseudobulk = (args.pseudobulk == "True")

adata = sc.read(args.i)
adata._inplace_subset_obs(adata.obs["dataset"].isin(datasets))
adata._inplace_subset_obs(adata.obs["inferred_cell_type"]==args.celltype)

sc.pp.filter_cells(adata,min_genes=3700)
min_cells = adata.obs.index.shape[0] / 4
sc.pp.filter_genes(adata,min_cells=min_cells)

adata2 = adata.copy()
sc.pp.normalize_total(adata2,target_sum = 10000)
sc.pp.log1p(adata2) 
sc.pp.highly_variable_genes(adata2,min_mean=0.1,max_mean = 100000,min_disp=0.2,max_disp=1000000,flavor="cell_ranger")
selected_genes = adata2.var[(adata2.var["means"]>0.1) & (adata2.var["dispersions"]>0.2)].index
adata = adata[:,selected_genes]
if args.normalize:
    adata = adata2[:,selected_genes]
if do_pseudobulk:
    df = pd.DataFrame(np.log1p(adata.X.todense()))
    df.index = adata.obs.index
    df.columns = adata.var.index
    df["individual"] = adata.obs["individual"].astype(str)
    #Filter out individuals who do not have enough cells of this cell type
    

    
    if args.normalize:
        means = df.groupby("individual").mean()
    else:
        means = np.expm1(df.groupby("individual").mean())
        counts = df.groupby("individual").count()
        means = means[counts.iloc[:,0]>5]
        #Multiply the mean value by the number of cells of each individual, to get something closer to the number of counts
        for individual in means.index:
            means.loc[individual,:] = means.loc[individual,:] *  counts.loc[individual,counts.columns[0]]

    attributes = {"dataset":{},"sex":{},"disease":{},"ethnic_group":{},"age":{},"body_mass_index":{}}
    for cell in df.index:
            individual = df.loc[cell,"individual"]
            if individual in means.index and (not individual in attributes["dataset"]):
                for attribute in attributes:
                        attributes[attribute][individual] = adata.obs.loc[cell,attribute]

    coldata=pd.DataFrame(attributes)
    coldata.index.name = "individual"
    cells = pd.Series(coldata.index)
    if args.normalize:
        df = means
    else:
        df = means.astype(int)
else: #No pseudobulk
    if args.normalize:
        df = pd.DataFrame(adata.X.todense())
    else:
        df = pd.DataFrame(adata.X.todense().astype(int)) #DESeq2 expects counts (integers)
    df.columns = list(adata.var.index)
    df.index = list(adata.obs.index)
    coldata = pd.DataFrame({"disease":adata.obs["disease"],"dataset":adata.obs["dataset"],"sex":adata.obs["sex"],"individual":adata.obs["individual"],"ethnic_group":adata.obs["ethnic_group"],
    "age":adata.obs["age"],"body_mass_index":adata.obs["body_mass_index"]})
    cells = pd.Series(adata.obs.index)

genes = pd.Series(adata.var.index)


df = df.transpose()
df.index.name="gene"

df.to_csv("counts.csv",index_label=False,index=False)
coldata.to_csv("coldata.csv")
genes.to_csv("genes.csv",index_label=False,index=False,header=False)
cells.to_csv("cells.csv",index_label=False,index=False,header=False)
