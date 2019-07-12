library(DESeq2)
library("BiocParallel")
register(MulticoreParam(8))
args = commandArgs(trailingOnly=TRUE)

counts = read.csv(args[1],header=TRUE,sep=",")
matrix = data.matrix(counts)
cells = scan(args[4], sep=',', what = "", quiet = TRUE)
genes = scan(args[3], sep=',', what = "", quiet = TRUE)
colnames(matrix) <- cells
rownames(matrix) <- genes

coldata = read.csv(args[2])
dds <- DESeqDataSetFromMatrix(countData=matrix,colData=coldata,design = as.formula(args[5]))
dds <- DESeq(dds,parallel=TRUE) #minReplicatesForReplace=Inf
res <- results(dds,contrast = c("disease","normal","T2D"),cooksCutoff=FALSE) #cooksCutoff=FALSE, independentFiltering=FALSE
print(dds)
resOrdered <- res[order(res$pvalue),]
print(resOrdered)
print(summary(res))
write.csv(as.data.frame(resOrdered), file="results_DESeq2.csv")
