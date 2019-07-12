#!/usr/bin/env nextflow

scripts="${workflow.projectDir}/scripts"

process merge_matrices{
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    output:
        file "merged_counts_matrix.h5ad" into MERGED_COUNTS_MATRIX

    """
        python ${scripts}/merge_matrices.py -datasets ${params.pipeline.datasets} -folder ${params.datasets} 
    """
}

MERGED_COUNTS_MATRIX.into{
    COUNTS_FOR_SCANPY
    COUNTS_FOR_OLS
    COUNTS_FOR_R
    COUNTS_FOR_R_NORM
    COUNTS_FOR_BOXPLOT
}

process scanpy_pipeline {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${PWD}/results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_FOR_SCANPY
    output:
        file "processed_anndata.h5ad" into PROCESSED_ANNDATA

    """
        python ${scripts}/scanpy_pipeline.py -i ${counts} -o processed_anndata.h5ad -datasets ${params.pipeline.datasets}
    """
}

process ordinary_least_squares {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${PWD}/results/partial_results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_FOR_OLS
        
    output:
        file "results_ols.csv" into RESULTS_OLS

    """
        python ${scripts}/ols.py -i ${counts} -celltype ${params.de.celltype} -datasets ${params.de.datasets} -design ${params.de.design} -pseudobulk ${params.de.pseudobulk}
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
        python ${scripts}/h5ad_to_R.py -i ${counts} -celltype ${params.de.celltype} -datasets ${params.de.datasets} -pseudobulk ${params.de.pseudobulk}
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
    conda "${workflow.projectDir}/envs/Renv.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${PWD}/results/partial_results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_DESEQ2
        file coldata from COLDATA_DESEQ2
        file genes from GENES_DESEQ2
        file cells from CELLS_DESEQ2
        
    output:
        file "results_DESeq2.csv" into RESULTS_DESEQ2

    """
        Rscript ${scripts}/script_DESeq2.R ${counts} ${coldata} ${genes} ${cells} ${params.de.design}
    """
}

process edgeR {
    conda "${workflow.projectDir}/envs/Renv.yml"
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${PWD}/results/partial_results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_EDGER
        file coldata from COLDATA_EDGER
        file genes from GENES_EDGER
        file cells from CELLS_EDGER
        
    output:
        file "results_edgeR.csv" into RESULTS_EDGER

    """
        Rscript ${scripts}/script_edgeR.R ${counts} ${coldata} ${genes} ${cells} ${params.de.design}
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
        python ${scripts}/h5ad_to_R.py -i ${counts} -celltype ${params.de.celltype} -datasets ${params.de.datasets} -pseudobulk False -normalize
    """
}

process mixed_model {
    conda "${workflow.projectDir}/envs/Renv.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${PWD}/results/partial_results", mode: 'copy', overwrite: true
    input:
        file counts from COUNTS_R_NORM
        file coldata from COLDATA_R_NORM
        file genes from GENES_R_NORM
        file cells from CELLS_R_NORM
        
    output:
        file "results_mixed.csv" into RESULTS_MIXED

    """
        Rscript ${scripts}/script_mixed.R ${counts} ${coldata} ${genes} ${cells} ${params.de.design}
    """
}

process merge_results {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${PWD}/results", mode: 'copy', overwrite: true
    input:
        file results_ols from RESULTS_OLS
        file results_DESeq2 from RESULTS_DESEQ2
        file results_edgeR from RESULTS_EDGER
        file results_mixed from RESULTS_MIXED
        
    output:
        file "results_complete.csv" into RESULTS_COMPLETE

    """
        python ${scripts}/merge_results.py -ols ${results_ols} -deseq2 ${results_DESeq2} -edger ${results_edgeR} -mixed ${results_mixed}
    """
}

process boxplots {
    conda "${workflow.projectDir}/envs/scanpy.yml"
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 2
    publishDir "${PWD}/results", mode: 'copy', overwrite: true
    input:
        file results from RESULTS_COMPLETE
        file anndata from COUNTS_FOR_BOXPLOT

        
    output:
        file "boxplots" into BOXPLOTS

    """
        python ${scripts}/plot_expressions_donors.py -anndata ${anndata} -results ${results} -celltype ${params.de.celltype} -datasets ${params.de.datasets} -correct_covariates ${params.plot.correct_covariates}
    """
}
