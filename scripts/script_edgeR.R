#!/usr/bin/Rscript

library("edgeR")
args = commandArgs(trailingOnly=TRUE)

counts = read.csv(args[1],header=TRUE,sep=",")
matrix = data.matrix(counts)
coldata = read.csv(args[2])
genes = scan(args[3], sep=',', what = "", quiet = TRUE)
cells = scan(args[4], sep=',', what = "", quiet = TRUE)

colnames(matrix) <- cells
rownames(matrix) <- genes
#diseases = scan(args[3], sep=',', what = "", quiet = TRUE)
#group = factor(diseases)
#disease = factor(disease)

#et <- exactTest(y)
#topTags(et)
variables = strsplit(substring(args[5],2),split="\\+")
factors = vector("list",length(variables[[1]]))
formula = "~"
count=1
disease_index=1
for (x in variables[[1]]){
    factors[[count]] = factor(coldata[x][,1])
    formula = paste(formula,"factors[[",count,"]]+",sep="")
    if (x=="disease"){
        disease_index=count
    }
    count = count +1
}
formula = substring(formula,1,nchar(formula)-1)

#f_disease = factor(coldata["disease"][,1])
#f_sex = factor(coldata["sex"][,1])
#f_dataset = factor(coldata["dataset"][,1])
#design <- model.matrix(~f_sex + f_dataset + f_disease)
print(formula)
print(factors)
print(factors[[1]])
design <- model.matrix(as.formula(formula))

y = DGEList(counts=matrix,group=factors[[disease_index]],genes=genes)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design) #robust True ?
qlf <- glmQLFTest(fit)
topTags(qlf)
results <- topTags(qlf, n=Inf, adjust.method="BH", sort.by="PValue")
write.csv(as.data.frame(results), file="results_edgeR.csv")