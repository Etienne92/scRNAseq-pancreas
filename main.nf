#!/usr/bin/env nextflow

scripts="${workflow.projectDir}/scripts"

DATASETS = Channel.fromPath("datasets/*", checkIfExists: true,type:'dir').collect()

process merge_matrices{
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    input:
        file datasets from DATASETS
        
    output:
        file "merged_counts_matrix.h5ad" into MERGED_COUNTS_MATRIX

    """
        python ${scripts}/merge_matrices.py ${datasets} 
    """
}

MERGED_COUNTS_MATRIX.into{
    COUNTS_FOR_SCANPY
    COUNTS_FOR_PYTHON
    COUNTS_FOR_R
    COUNTS_FOR_R_NORM
}

process scanpy_pipeline {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_FOR_SCANPY
    output:
        file "processed_anndata.h5ad" into PROCESSED_ANNDATA

    """
        python ${scripts}/scanpy_pipeline.py -i ${counts} -o processed_anndata.h5ad -datasets ${params.datasets}
    """
}

process python_de {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_FOR_PYTHON
        
    output:
        file "python_pvalues.csv" into PYTHON_PVALUES

    """
        python ${scripts}/python_de.py -i ${counts} -celltype ${params.celltype} -datasets ${params.datasets} -design ${params.design} -${params.mean}
    """
}

process h5ad_to_R {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    input:
        file counts from COUNTS_FOR_R
        
    output:
        file "counts.csv" into COUNTS_R
        file "coldata.csv" into COLDATA_R
        file "genes.csv" into GENES_R
        file "cells.csv" into CELLS_R

    """
        python ${scripts}/h5ad_to_R.py -i ${counts} -celltype ${params.celltype} -datasets ${params.datasets} -${params.mean}
    """
}
COUNTS_R.into{
    COUNTS_DESEQ2
    COUNTS_EDGER
}
COLDATA_R.into{
    COLDATA_DESEQ2
    COLDATA_EDGER
}
GENES_R.into{
    GENES_DESEQ2
    GENES_EDGER
}
CELLS_R.into{
    CELLS_DESEQ2
    CELLS_EDGER
}

process DESeq2 {
    conda "${workflow.projectDir}/envs/DESeq2.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_DESEQ2
        file coldata from COLDATA_DESEQ2
        file genes from GENES_DESEQ2
        file cells from CELLS_DESEQ2
        
    output:
        file "DESeq2_results.csv" into RESULTS_DESEQ2

    """
        Rscript ${scripts}/script_DESeq2.R ${counts} ${coldata} ${genes} ${cells} ${params.design}
    """
}

process edgeR {
    conda "${workflow.projectDir}/envs/edgeR.yml"
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_EDGER
        file coldata from COLDATA_EDGER
        file genes from GENES_EDGER
        file cells from CELLS_EDGER
        
    output:
        file "edgeR_results.csv" into RESULTS_EDGER

    """
        Rscript ${scripts}/script_edgeR.R ${counts} ${coldata} ${genes} ${cells} ${params.design}
    """
}

process h5ad_to_R_norm {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    input:
        file counts from COUNTS_FOR_R_NORM
        
    output:
        file "counts.csv" into COUNTS_R_NORM
        file "coldata.csv" into COLDATA_R_NORM
        file "genes.csv" into GENES_R_NORM
        file "cells.csv" into CELLS_R_NORM

    """
        python ${scripts}/h5ad_to_R.py -i ${counts} -celltype ${params.celltype} -datasets ${params.datasets} -nomean -normalize
    """
}

process mixed_model {
    conda "${workflow.projectDir}/envs/lme4.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_R_NORM
        file coldata from COLDATA_R_NORM
        file genes from GENES_R_NORM
        file cells from CELLS_R_NORM
        
    output:
        file "mixed_results.csv" into RESULTS_MIXED

    """
        Rscript ${scripts}/mixed_model.R ${counts} ${coldata} ${genes} ${cells} ${params.design}
    """
}

process merge_results {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "results", mode: 'copy', overwrite: true
    input:
        file python_results from PYTHON_PVALUES
        file deseq2_results from RESULTS_DESEQ2
        file edger_results from RESULTS_EDGER
        file mixed_results from RESULTS_MIXED
        
    output:
        file "results_complete.csv" into RESULTS_COMPLETE

    """
        python ${scripts}/merge_results.py -python ${python_results} -deseq2 ${deseq2_results} -edger ${edger_results} -mixed ${mixed_results}
    """
}