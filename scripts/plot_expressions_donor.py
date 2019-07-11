import scanpy as sc 
import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import statsmodels.api as sm
import statsmodels.formula.api as smf
import argparse

parser = argparse.ArgumentParser(description='Creates boxplots of the expressions of the most differentially expressed genes.')
parser.add_argument("-anndata",type=str,help="Raw anndata object")
parser.add_argument("-celltype",type=str,default= "beta cell",help="Only considers the cells which have this inferred cell type.")
parser.add_argument("-results",type=str,help = "Results file of the most differentially expressed genes")
parser.add_argument("-datasets",type=str,default="all",help="Datasets to use")
args=parser.parse_args()

adata = sc.read(args.i)
cell_type = args.celltype
datasets = args.datasets.split(",")


adata = adata[adata.obs["dataset"] in datasets]
adata = adata[adata.obs["inferred_cell_type"] == args.celltype]
sc.pp.filter_cells(adata,min_counts = 80000)
sc.pp.filter_cells(adata,min_genes=3700)
#min_cells = adata.obs.index.shape[0] / 4
sc.pp.filter_genes(adata,min_cells=5)
sc.pp.normalize_total(adata,target_sum = 10000,fraction=0.3)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,min_mean=0.1,max_mean = 10000000,min_disp=0.3,max_disp=100000,flavor="cell_ranger")
selected_genes = adata.var[(adata.var["means"]>0.01) & (adata.var["dispersions"]>0.03)].index
adata = adata[:,selected_genes]


df = pd.DataFrame(adata.X.todense())
df.index = adata.obs.index
df.columns = adata.var.index
df["individual"] = adata.obs["individual"]
df["dataset"] = adata.obs["dataset"]
df["disease"] = adata.obs["disease"]
df["sex"] = adata.obs["sex"]
dataset_mapping = {"E-MTAB-5061":"Donors from Segerstolpe et al","E-GEOD-81608":"Donors from Xin et al","E-ENAD-27":"Donors from Lowler et al"}


#df["individual"] = df["individual"].map(lambda x : x.split("_")[-1])
#df["individual"] = df2["individual"].map(lambda x : x.split("_")[-1])

def boxplot(df,gene,params=None):
    gene_name = ""
    if "gene" in adata.var.columns:
        transcript_id = gene[:-1] + "." + gene[-1]
        gene_name = " (" + adata.var.loc[transcript_id,"gene"] + ")"
    datasets = df["dataset"].unique()
    nb_datasets = len(datasets)
    width_ratios = [len(df.loc[df["dataset"]==dataset,:]["individual"].unique()) for dataset in datasets]

    fig, axs = plt.subplots(1, nb_datasets, sharey=True,gridspec_kw={'width_ratios': width_ratios})
    for i,dataset in enumerate(datasets):
        subset = df.loc[df["dataset"]==dataset,:]
        subset["individual"] = subset["individual"].astype(str)
        subset = subset.sort_values(by = "individual")
        seaborn.boxplot(y=gene,x="individual",data=subset,hue="disease",showmeans=True,ax=axs[i]).set(xlabel=dataset_mapping[dataset])

    axs[0].tick_params(top=False, bottom=False, left=True, right=False, labelleft=True, labelbottom=False)
    axs[0].spines['right'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    axs[0].yaxis.set_label_text("Log expression of " + gene + gene_name)
    for i in [1,2]:
        for pos in ["top","bottom","left","right"]:
            axs[i].spines[pos].set_visible(False)
            axs[i].tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
            #axs[i].set_yticklabels([])
            axs[i].yaxis.set_label_text("")
    handles, _ = axs[2].get_legend_handles_labels()
    axs[2].legend(handles, ["Type II Diabetes","Healthy"],loc = 9)
    axs[0].get_legend().set_visible(False)
    axs[1].get_legend().set_visible(False)
    if params is not None:
        for i,dataset in enumerate(datasets):
                axs[i].hlines(params["means"][dataset],-1,width_ratios[i],color="red")
                axs[i].hlines(params["healthy"][dataset],-1,width_ratios[i],color="green")
                axs[i].hlines(params["T2D"][dataset],-1,width_ratios[i],color="yellow")
    fig.savefig("expression_donors_"+gene+".png",dpi=1000)

def means_individual(df):
        means = df.groupby("individual").mean()
        means = means
        attributes = {"dataset":{},"sex":{},"disease":{}}
        for cell in df.index:
                individual = df.loc[cell,"individual"]
                for attribute in attributes:
                        attributes[attribute][individual] = df.loc[cell,attribute]

        df_attributes=pd.DataFrame(attributes)
        means = means.join(df_attributes)
        return means


def correct_covariates(df,gene):
        means = means_individual(df)
        df_modified = df.copy()
        results = smf.ols(gene + '~ disease + sex + dataset', data=means).fit()

        for param in results.params.index:
                if param[:7]=="dataset":
                        dataset = param[param.find("[")+3:param.find("]")]
                        for cell in df_modified.index:
                                if df_modified.loc[cell,"dataset"] == dataset:
                                        df_modified.loc[cell,gene] = df_modified.loc[cell,gene] - results.params[param]
                elif param[:3] == "sex":
                        sex = param[param.find("[")+3:param.find("]")]
                        for cell in df_modified.index:
                                if df_modified.loc[cell,"sex"] == sex:
                                        df_modified.loc[cell,gene] = df_modified.loc[cell,gene] - results.params[param]
        return df_modified


boxplot(df,gene)

#df_modified = correct_covariates(df,gene)
#boxplot(df_modified,gene)