library(glue)
library(Seurat)
library(stringr)

# load data
load("/data/gusev/USERS/scarver/Light_Hvratin/data/GSE102827_MATRIX.RDat")

genes.keep = readRDS("/data/gusev/USERS/scarver/Light_Hvratin/data/genes_names_QC.rds")
cells.keep = readRDS("/data/gusev/USERS/scarver/Light_Hvratin/data/cell_names_QC.rds")

# subset metadata and expression data by QC'd genes and cells
filtered_exp = exp[genes.keep,cells.keep ]

# parse args
args = commandArgs(trailingOnly=TRUE)

name = args[1]
method = args[2] 
data_param = args[3] 
num_cells = args[4] 
seed = args[5] 
dim = args[6] 
vargenes = args[7] 

# read in embeddings
embeddings = read.table(glue("/data/gusev/USERS/scarver/Light_Hvratin/cNMF_final/cNMF_output/{name}/{name}.usages.k_{dim}.dt_2_0.consensus.txt"),sep="\t",row.names=1,header=TRUE)


print(glue("embedding dim: {dim(embeddings)}"))
# obtain time metadata and convert to time continuum
meta$time = as.numeric(substr((meta$stim),1,1))

# subset metadata time to match the embeddings
final_meta = meta[match(rownames(embeddings), rownames(meta)),]

# compute adjusted r squared
r_sq = summary(lm(scale(final_meta$time) ~ scale(embeddings[,1:dim(embeddings)[2]])) )$adj.r.sq

# compute correlation per PC
corr = rep(NA, dim(embeddings)[2])
for (i in 1:dim(embeddings)[2]){
   corr[i] = summary(lm(scale(final_meta$time) ~ scale(embeddings[,i])) )$adj.r.sq    
}

df = data.frame(method= method, data_param = data_param, num_cells = num_cells, seed = seed, dim=dim, vargenes=vargenes, adj.rsq = r_sq, max.rsq =max(corr))

# append run to output dataframe
write.table(df, file = "/data/gusev/USERS/scarver/Light_Hvratin/output_final_format.txt", append = TRUE, quote = FALSE, sep = " ", row.names = FALSE,col.names = FALSE)
