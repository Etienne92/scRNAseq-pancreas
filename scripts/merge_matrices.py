import scanpy as sc 
import pandas as pd
import scipy
import anndata
import numpy as np
import sys
import os
from ID2name import ID2name_converter

import argparse



parser = argparse.ArgumentParser(description='Merge several expression matrices and add metadata.')
parser.add_argument("-datasets",type=str,help="Datasets to use")
parser.add_argument("-folder",type=str,help="Folder containing the datasets")
args=parser.parse_args()

datasets = args.datasets.split(",")

def add_metadata(adata,sdrfpath):
    """Parse the sdrf file and add the metadata to the anndata"""
    with open(sdrfpath) as sdrf:
        header=next(sdrf).split("\t")
        enaIndex=-1
        infoIndex={}
        infoData={}
        useful_info = ["body mass index","bmi","age","sex","inferred cell type","ancestry category","ethnic group","disease","individual","inferred cell type"]
        for i in range(len(header)):
            if header[i][:15]=="Characteristics" or header[i][:12]=="Factor Value":
                s=header[i]
                info = s[s.find("[")+1:s.find("]")] #select the substring between the brackets
                if info in useful_info:
                    infoIndex[info]=i
                    infoData[info]={}
            if header[i]=="Comment[ENA_RUN]" or header[i] == "Comment [ENA_RUN]" or header[i] == "Source Name":
                enaIndex=i
        for line in sdrf:
            linesplit=line.split("\t")
            ena=linesplit[enaIndex]
            if ena in adata.obs.index:
                for info in infoIndex:
                    infoData[info][ena]=linesplit[infoIndex[info]].rstrip()

        #Rename some metadata names (remove spaces and make consistent across datasets)
        if "bmi" in infoData:
            infoData["body_mass_index"] = infoData.pop("bmi")
            infoIndex["body_mass_index"] = infoIndex.pop("bmi")
        if "body mass index" in infoData:
            infoData["body_mass_index"] = infoData.pop("body mass index")
            infoIndex["body_mass_index"] = infoIndex.pop("body mass index")
        if "ancestry category" in infoData:
            infoData["ethnic_group"] = infoData.pop("ancestry category")
            infoIndex["ethnic_group"] = infoIndex.pop("ancestry category")
        if "ethnic group" in infoData:
            infoData["ethnic_group"] = infoData.pop("ethnic group")
            infoIndex["ethnic_group"] = infoIndex.pop("ethnic group")
        if "inferred cell type" in infoData:
            infoData["inferred_cell_type"] = infoData.pop("inferred cell type")
            infoIndex["inferred_cell_type"] = infoIndex.pop("inferred cell type")

        #Quantitative information is converted to numeric values, so that cellxgene can use it as continuous metadata, instead of categorical.
        quantitative_info=["age","body_mass_index"]
        for info in infoIndex:
            if not info in quantitative_info:
                adata.obs[info]=pd.Series(infoData[info])
            else:
                series = pd.Series(infoData[info])
                adata.obs[info]=pd.to_numeric(series.replace("not available",np.nan))


cell_types = {"pancreatic A cell": "alpha cell","pancreatic D cell":"delta cell","pancreatic ductal cell":"ductal cell","pancreatic epsilon cell":"epsilon cell","pancreatic PP cell": "PP cell",
 "pancreatic stellate cell":"stellate cell","type B pancreatic cell": "beta cell"}

adatas = []
for dataset in datasets:
    folder = args.folder + "/" + dataset
    adata = sc.read_10x_mtx(folder)
    sc.pp.filter_cells(adata,min_counts = 0) #just for getting the counts for each cell
    add_metadata(adata,folder+"/sdrf.txt")
    adata = adata[adata.obs["inferred_cell_type"]!="not applicable"] #Filter out samples which are not single cells

    #Change the names of the donors to reflect their disease: H1,H2,H3... for healthy and D1,D2,D3 for patients with T2D
    individual_mapping={} 
    counts_healthy=1
    counts_diabetic=1
    for i in range(len(adata.obs["individual"])):
        if not adata.obs["individual"][i] in individual_mapping:
            if adata.obs["disease"][i] == "normal":
                individual_mapping[adata.obs["individual"][i]] = dataset + "_H"+str(counts_healthy)
                counts_healthy+=1
            else:
                individual_mapping[adata.obs["individual"][i]] = dataset + "_D"+str(counts_diabetic)
                counts_diabetic+=1
    adata.obs["individual"].replace(individual_mapping,inplace=True)
    adata.obs["inferred_cell_type"].replace(cell_types,inplace=True)  #unify cell type names across datasets
    adata.obs["disease"].replace({"type II diabetes mellitus":"T2D"},inplace=True)
    adata.var.drop("gene_ids",axis=1,inplace=True)

    adatas.append(adata)

merged = adatas[0].concatenate(adatas[1:],batch_key = "dataset",batch_categories = datasets)


#Use gene names instead of gene IDs. 
converter=ID2name_converter(os.path.dirname(sys.argv[0]) +"/../data/Homo_sapiens.GRCh38.95.gtf.gz")
merged.var.index = merged.var.index.map(lambda x: converter.id2name(x))
merged.var_names_make_unique()

merged.write("merged_counts_matrix.h5ad")
