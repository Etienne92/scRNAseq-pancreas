library(lme4)

args = commandArgs(trailingOnly=TRUE)

counts = read.csv(args[1],header=TRUE,sep=",")
matrix = data.matrix(counts)
cells = scan(args[4], sep=',', what = "", quiet = TRUE)
genes = scan(args[3], sep=',', what = "", quiet = TRUE)
colnames(matrix) <- cells
rownames(matrix) <- genes

coldata = read.csv(args[2])
f_disease = factor(coldata["disease"][,1])
f_sex = factor(coldata["sex"][,1])
f_dataset = factor(coldata["dataset"][,1])
f_individual = factor(coldata["individual"][,1])


l=list()
for (i in 1:length(genes)){
    expression = matrix[genes[i],]
    mres1 <- lmer(expression ~ f_sex + f_dataset + (1 | f_individual))
    mres2 <- lmer(expression ~ f_sex + f_dataset + (1 | f_disease/f_individual))
    anova_res <- anova(mres1, mres2)
    pval <- anova_res$`Pr(>Chisq)`[2]
    l[genes[i]] = pval
}
or = order(unlist(l))
results = l[or]

adj_pvals = p.adjust(results,method="BH")

df = data.frame("adj p-val"=adj_pvals,"p-val"=unlist(results))
colnames(df)[1] <- "gene"

write.csv(df, file="mixed_results.csv")