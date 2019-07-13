# scRNAseq-pancreas

This is a Nextflow worfklow which reproduces the main results of my project during my internship at the EBI. It is used to cluster and visualize cells from several single-cell RNA-seq human pancreas datasets, and performs differential expression to identify genes linked to T2D. An important point of this differential expression is that it avoids pseudoreplication, which is very common with scRNAseq.

## Inputs
The pipeline takes as inputs several datasets in one directory. For each dataset, we need the result of the aggregation workflow and the sdrf file. The datasets I used can be downloaded [here](https://www.dropbox.com/s/eaaljveie79k2ec/datasets.zip?dl=0). They were obtained with my modified pipeline, but we would probably get similar results with the standard SCXA workflow (I only made the quality control filters stricter and I added a filter on the percentage of reads that can be pseudoaligned).

## Outputs
The pipeline outputs are:
- a processed anndata object, which contains clusters, UMAP coordinates and metadata. This is meant to be visualized with a tool like cellxgene.
- The ranked list of differentially expressed genes with their adjusted p-values, for each method.
- Boxplots of the most differentially expressed genes.
They are located in a directory called results.

## Execution
The workflow can be run from this repository by executing:
```
nextflow run Etienne92/scRNAseq-pancreas --datasets <datasets_dir>
```
where <datasets_dir> is the path to the directory which contains the datasets (typically a path to [this directory](https://www.dropbox.com/s/eaaljveie79k2ec/datasets.zip?dl=0)).
