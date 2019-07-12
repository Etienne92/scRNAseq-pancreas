import pandas as pd 
import argparse

parser = argparse.ArgumentParser(description='Merge the results of all methods into one single csv file.')
parser.add_argument("-ols",type=str)
parser.add_argument("-edger",type=str)
parser.add_argument("-deseq2",type=str)
parser.add_argument("-mixed",type=str)
args=parser.parse_args()


results={}

ols_results = pd.read_csv(args.ols)
results["OLS"]= [f"{ols_results.loc[i,'gene']} ({float(ols_results.loc[i,'adj p-value']):.3})" for i in range(ols_results.shape[0])][:1000]

edger_results = pd.read_csv(args.edger)
results["edgeR"]= [f"{edger_results.loc[i,'genes']} ({float(edger_results.loc[i,'FDR']):.3})" for i in range(edger_results.shape[0])][:1000]

deseq2_results = pd.read_csv(args.deseq2)
deseq2_results.columns = ["gene","baseMean","logFold2change","lfcSE","stat","pvalue","padj"]
results["DESeq2"]= [f"{deseq2_results.loc[i,'gene']} ({float(deseq2_results.loc[i,'padj']):.3})" for i in range(deseq2_results.shape[0])][:1000]

mixed_results = pd.read_csv(args.mixed)
mixed_results.columns = ["gene","adj p-val","p-val"]
results["mixed model"] =[f"{mixed_results.loc[i,'gene']} ({float(mixed_results.loc[i,'adj p-val']):.3})" for i in range(mixed_results.shape[0])][:1000]

df = pd.DataFrame(results)
df.to_csv("results_complete.csv")

