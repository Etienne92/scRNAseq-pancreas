import scanpy as sc 
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import argparse

parser = argparse.ArgumentParser(description='Uses Scanpy to normalize and filter cells/genes, and then uses statsmodels to identify differentially expressed genes.')
parser.add_argument("-i",type=str,help="Input file (counts matrix)")
parser.add_argument("-celltype",type=str,default= "beta cell",help="Only considers the cells which have this inferred cell type.")
parser.add_argument("-datasets",type=str,help="Datasets to use.")
parser.add_argument("-design",type=str,help="Design for the linear model, for example dataset+sex+disease")
parser.add_argument('-mean', dest='mean', action='store_true')
parser.add_argument('-nomean', dest='mean', action='store_false')
args=parser.parse_args()

datasets = args.datasets.split(",")
datasets = [x.strip() for x in datasets]

adata = sc.read(args.i)
adata = adata[adata.obs["dataset"].isin(datasets)]
adata = adata[adata.obs["inferred_cell_type"] == args.celltype]

sc.pp.filter_cells(adata,min_counts = 80000)
sc.pp.filter_cells(adata,min_genes=3700)
min_cells = adata.obs.index.shape[0] / 4
sc.pp.filter_genes(adata,min_cells=min_cells)
sc.pp.normalize_total(adata,target_sum = 10000,fraction=0.3)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,min_mean=0.1,max_mean = 10000000,min_disp=0.2,max_disp=100000,flavor="cell_ranger")
selected_genes = adata.var[(adata.var["means"]>0.1) & (adata.var["dispersions"]>0.2)].index
adata = adata[:,selected_genes]


df = pd.DataFrame(adata.X.todense())
df.index = adata.obs.index
df.columns = [x.replace("-","").replace(".","") for x in adata.var.index]
df["individual"] = adata.obs["individual"]
df["dataset"] = adata.obs["dataset"]
df["disease"] = adata.obs["disease"]
df["sex"] = adata.obs["sex"]
df["ethnic_group"] = adata.obs["ethnic_group"]
df["body_mass_index"] = adata.obs["body_mass_index"]
df["age"] = adata.obs["age"]

if args.mean:
        #Filter out individuals who do not have enough cells of this cell type
        counts_donor = df.groupby("individual").count()
        df["cells_donor"] = [counts_donor.loc[df.loc[x,"individual"],"disease"] for x in df.index]
        df = df[df["cells_donor"]>5]

def linear_model(df,design):
        pvalues=[]
        for gene in df.columns:
                if (df[gene].dtype == "float32" or df[gene].dtype == "float64"):  #make sure the column is a gene, not metadata
                        results = smf.ols(gene +design, data=df).fit()
                        pvalues.append((gene,results.pvalues["disease[T.normal]"]))
        pvalues.sort(key = lambda x : float(x[1]))
        return(pvalues)

def adjust_pvalues(pvalues,FDR=0.05):
        """Correct p values with the Benjamini-Hochberg correction method."""
        values = [x[1] for x in pvalues]
        adjusted_values = multipletests(values,alpha=FDR,method="fdr_bh")
        return [(pvalues[i][0],adjusted_values[1][i]) for i in range(len(pvalues))]

def means_individual(df):
        means = df.groupby("individual").mean()
        attributes = {"dataset":{},"sex":{},"disease":{}}
        for cell in df.index:
                individual = df.loc[cell,"individual"]
                for attribute in attributes:
                        attributes[attribute][individual] = df.loc[cell,attribute]

        df_attributes=pd.DataFrame(attributes)
        means = means.join(df_attributes)
        return means

if args.mean:
    means = means_individual(df)
    pvalues= linear_model(means,args.design) 
else:
    pvalues= linear_model(df,args.design)
apvalues = adjust_pvalues(pvalues)

genes = [x[0] for x in apvalues]
adjpvalues = [x[1] for x in apvalues]
nonadjpvalues = [x[1] for x in pvalues]
df_results=pd.DataFrame({"gene":genes,"adj p-value":adjpvalues,"pvalues":nonadjpvalues})

df_results.to_csv("python_pvalues.csv")
