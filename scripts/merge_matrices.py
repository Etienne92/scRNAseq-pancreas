import scanpy as sc 
import pandas as pd
import scipy
import anndata
import numpy as np
import sys
from ID2name import ID2name_converter



def add_metadata(adata,sdrfpath):
    """Parse the sdrf file and add the metadata to the anndata"""
    with open(sdrfpath) as sdrf:
        header=next(sdrf).split("\t")
        enaIndex=-1
        infoIndex={}
        infoData={}
        for i in range(len(header)):
            if header[i][:15]=="Characteristics" or header[i][:12]=="Factor Value":
                s=header[i]
                info = s[s.find("[")+1:s.find("]")] #select the substring between the brackets
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

        #Make metadata names consistent across datasets
        if "bmi" in infoData:
            infoData["body mass index"] = infoData.pop("bmi")
            infoIndex["body mass index"] = infoIndex.pop("bmi")
        if "ancestry category" in infoData:
            infoData["ethnic group"] = infoData.pop("ancestry category")
            infoIndex["ethnic group"] = infoIndex.pop("ancestry category")

        #Quantitative information is converted to numeric values, so that cellxgene can use it as continuous metadata, instead of categorical.
        quantitative_info=["age","body mass index"]
        for info in infoIndex:
            if not info in quantitative_info:
                adata.obs[info]=pd.Series(infoData[info])
            else:
                series = pd.Series(infoData[info])
                adata.obs[info]=pd.to_numeric(series.replace("not available",np.nan))


cell_types = {"pancreatic A cell": "alpha cell","pancreatic D cell":"delta cell","pancreatic ductal cell":"ductal cell","pancreatic epsilon cell":"epsilon cell","pancreatic PP cell": "PP cell",
 "pancreatic stellate cell":"stellate cell","type B pancreatic cell": "beta cell"}


datasets = sys.argv[1:]

dfs=[]
attributes = {
    "dataset":[],
    "individual":[],
    "inferred cell type":[],
    "disease":[],
    "body mass index":[],
    "age":[],
    "sex":[],
    "ethnic group":[],
    "n_counts":[]
}
for dataset in datasets:
    adata = sc.read_10x_mtx(dataset)
    sc.pp.filter_cells(adata,min_counts = 0) #just for getting the counts for each cell
    add_metadata(adata,dataset+"/sdrf.txt")
    adata = adata[adata.obs["inferred cell type"]!="not applicable"] #Filter out samples which are not single cells
    matrix = adata.X.todense()
    df = pd.DataFrame(matrix)
    df.index = adata.obs.index
    df.columns = adata.var.index
    dfs.append(df)
    attributes["dataset"] = attributes["dataset"] + [dataset] * df.shape[0]
    for attribute in attributes:
        if attribute != "dataset" and attribute!="individual":
            if attribute in adata.obs:
                attributes[attribute] = attributes[attribute] + list(adata.obs[attribute])
            else:
                attributes[attribute] = attributes[attribute] + ["not available"] * len(adata.obs)

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
    attributes["individual"] = attributes["individual"] + [individual_mapping.get(individual,individual) for individual in adata.obs["individual"]]

attributes["inferred cell type"] = [cell_types.get(celltype,celltype) for celltype in attributes["inferred cell type"]] #unify cell type names across datasets
attributes["disease"] = [x if x=="normal" else "T2D" for x in attributes["disease"]]
merged_df = pd.concat(dfs,sort=True)
merged_df.fillna(0, inplace=True)

sparse_merged_matrix = scipy.sparse.csr.csr_matrix(merged_df)

adata = anndata.AnnData(sparse_merged_matrix)
adata.obs.index = merged_df.index
adata.var.index = merged_df.columns
for attribute in attributes:
    corrected_name = attribute.replace(" ","_")
    adata.obs[corrected_name] = attributes[attribute]

#Use gene names instead of gene IDs. 
converter=ID2name_converter("../data/Homo_sapiens.GRCh38.95.gtf.gz")
adata.var.index = adata.var.index.map(lambda x: converter.id2name(x))
adata.var_names_make_unique()

adata.write("merged_counts_matrix.h5ad")
