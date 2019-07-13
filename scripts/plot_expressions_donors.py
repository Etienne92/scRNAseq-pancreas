import scanpy as sc 
import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import statsmodels.formula.api as smf
import os
import argparse

parser = argparse.ArgumentParser(description='Creates boxplots of the expressions of the most differentially expressed genes.')
parser.add_argument("-anndata",type=str,help="Raw anndata object")
parser.add_argument("-celltype",type=str,default= "beta cell",help="Only considers the cells which have this inferred cell type.")
parser.add_argument("-results",type=str,help = "Results file of the most differentially expressed genes")
parser.add_argument("-datasets",type=str,default="all",help="Datasets to use")
parser.add_argument("-correct_covariates",type=str,default="False",help="If True, remove the unwanted variations before plotting.")
args=parser.parse_args()

adata = sc.read(args.anndata)
cell_type = args.celltype
datasets = args.datasets.split(",")
datasets.sort(reverse=True)
use_correct_covariates = (args.correct_covariates == "True")


adata._inplace_subset_obs(adata.obs["dataset"].isin(datasets))
adata._inplace_subset_obs(adata.obs["inferred_cell_type"] == args.celltype)
sc.pp.filter_cells(adata,min_counts = 80000)
sc.pp.filter_cells(adata,min_genes=3700)
sc.pp.filter_genes(adata,min_cells=3)
sc.pp.normalize_total(adata,target_sum = 10000)
sc.pp.log1p(adata)

df = pd.DataFrame(adata.X.todense())
df.index = adata.obs.index
df.columns = adata.var.index
df["individual"] = adata.obs["individual"]
df["dataset"] = adata.obs["dataset"]
df["disease"] = adata.obs["disease"]
df["sex"] = adata.obs["sex"]
dataset_mapping = {"E-MTAB-5061":"Donors from Segerstolpe et al","E-GEOD-81608":"Donors from Xin et al","E-ENAD-27":"Donors from Lawlor et al"}

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
        subset = df.loc[df["dataset"]==dataset,:].copy()
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
    fig.savefig("boxplots/"+gene+".png",dpi=1000)

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


def correct_covariates(df,gene,covariates=["dataset","sex"]):
        """Fit a linear model and remove the variation due to the covariates."""
        means = means_individual(df)
        df_modified = df.copy()
        model = gene  + '~ disease +' + "+".join(covariates)
        results = smf.ols(model, data=means).fit()

        for param in results.params.index:
                for covariate in covariates:
                        if param[:len(covariate)]== covariate:
                                factor_value = param[param.find("[")+3:param.find("]")]
                                for cell in df_modified.index:
                                        if df_modified.loc[cell,covariate] == factor_value:
                                                df_modified.loc[cell,gene] = df_modified.loc[cell,gene] - results.params[param]
        return df_modified


#Generate the list of genes which have an adjusted p-value < 0.05 for at least one method
results = pd.read_csv(args.results,index_col=0)
genes = {}
for method in results.columns:
        for i in range(15):
                entry = results.loc[i,method]
                gene = entry[:entry.find("(")-1]
                #adj_pvalue = float(entry[entry.find("(")+1:entry.find(")")])
                #if adj_pvalue>0.05:
                #        break
                if gene in adata.var.index:
                        genes[gene] = 1
                elif (gene[:-1] + "." + gene[-1]) in adata.var.index: #The . in the gene names are removed for the differential expression because they interfere with the R syntax.
                        genes[gene[:-1] + "." + gene[-1]] = 1
if not os.path.exists('boxplots'):
    os.makedirs('boxplots')

if use_correct_covariates:
        for gene in genes:
                df_modified = correct_covariates(df,gene)
                boxplot(df_modified,gene)
else:
        for gene in genes:
                boxplot(df,gene)