setwd("/data/gusev/SCRNA/HRVATIN")

args = commandArgs(trailingOnly=TRUE)

method = args[1] # c("scaled_NMF", "PCA", "Cond_PCA")
data_param = args[2] # c("all_cell_types","int_ex")
num_cells = args[3] # c("all", # cells to downsample)
seed = as.integer(args[4])
dim = as.integer(args[5]) # rank of matrix in DR method

library(Seurat)
library(glue)
print(method)
if (method == "NMF"){
    library(NMF)
}


if (method == "PCA" | method == "Cond_PCA"){
    #load("QC_GSE102827_MATRIX.RData") # log normalized data
    load("QC_GSE102827_MATRIX_20k_vargenes.RData")
    print("loaded data for pca or cond pca")
} else {
    load("counts_norm_sub_GSE102827_MATRIX_20k_vargenes.RData") # loads meta and exp, normalized data
    print("loaded count normalized data")
}

print(glue("method: {method}"))
print(glue("data_param: {data_param}"))
print(glue("num_cells: {num_cells}"))
print(glue("seed: {seed}"))
print(glue("dim: {dim}"))

exp_seu <- CreateSeuratObject(counts = exp, meta.data = meta)

print("start subsetting")
print(glue("dim data: {dim(exp_seu)}"))
if (data_param == "all_cell_types"){
    print("subsetting to all cell types exluding NA cells")
    cells.use = rownames(exp_seu@meta.data[!is.na(exp_seu@meta.data$maintype),])
    sub <- subset(exp_seu, cells = cells.use)

} else {
    print("subsetting to Excitatory and Interneurons")
    sub = subset(x = exp_seu, subset = maintype == "Excitatory" | maintype == 'Interneurons')
}

# subset data if necessary
if (num_cells == "all"){
    print("no downsampling")
    print(glue("downsampled dims: {dim(sub)}"))
} else {
    print("performing downsampling") 
    set.seed(seed)
    num_cells = as.numeric(num_cells)
    cells.use = sample(colnames(sub), size = num_cells)
    sub <- subset(sub, cells = cells.use)
    print("finished downsampling")
    print(glue("downsampled dims: {dim(sub)}"))
}

# properly scale data

if (method == "PCA" | method == "UMAP") {
    print("not conditional")
    sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures =20000)
    print('done variable features')

    var_feat = sub@assays$RNA@var.features
    sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000) 

    sub <- ScaleData(sub, features=var_feat, scale.max=9999999999)
    print("finished scaling data")

} else if (method == "Cond_PCA"){
    print("conditioning out cell type")
    print(glue("cell types: {unique(sub@meta.data$maintype)}"))
    sub@meta.data$maintype = as.factor(sub@meta.data$maintype)
    
    sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures =20000)
    print('done variable features')

    var_feat = sub@assays$RNA@var.features
    sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000) 
    
    sub <- ScaleData(sub, features=var_feat, vars.to.regress = "maintype" ,scale.max=9999999999)
    
    print("finished scaling data")
    print(dim(sub))
} else {
    print("no scaling, method is standard NMF")

}

print(glue("dimensions of scaled data: {dim(sub@assays$RNA@scale.data)}"))


if (method == "NMF"){
    print("performing standard NMF")
    unscaled_data = sub@assays$RNA@counts
    # remove rows or columns with all zeros
    unscaled_data = unscaled_data[rowSums(unscaled_data != 0) > 0, colSums(unscaled_data != 0) > 0]
    unscaled_data = t(round(as.matrix(unscaled_data),4) )
    print(glue("transposed data dim: {dim(unscaled_data)}"))
    print("performing NMF")
    out = nmf(unscaled_data, rank=dim)
    w_cells <- out@fit@W
    h_genes <- out@fit@H

    output <- list()
    output[["embeddings"]] = w_cells
    output[["loadings"]] = h_genes
    embeddings = output[["embeddings"]]
    #saveRDS(output, glue("/data/gusev/SCRNA/HRVATIN/downsampling/QC/vargenes_20k/NMF/{method}_cells_{num_cells}_seed_{seed}_{data_param}_rank_{dim}.rds"))
} else if (method == "UMAP") {
    print("performing UMAP")
    sub <- RunPCA(sub, npcs = dim)
    sub <- RunUMAP(sub, dims = 1:10)
    #sub <- RunUMAP(sub, dims = 1:10, umap.method ="umap-learn")

    embeddings = sub[["umap"]]@cell.embeddings
    loadings = NA

    print(glue("dim emb {dim(embeddings)}"))
    print(glue("load dim NA"))

    output <- list()
    output[["embeddings"]] = embeddings
    output[["loadings"]] = loadings
    
    #saveRDS(output, file = glue("/data/gusev/SCRNA/HRVATIN/downsampling/QC/vargenes_20k/UMAP/{method}_cells_{num_cells}_seed_{seed}_{data_param}.rds"))

} else if (method == "cNMF") {
    print("obtaining data for cNMF")
    unscaled_data = sub@assays$RNA@counts
    # remove rows or columns with all zeros
    unscaled_data = unscaled_data[rowSums(unscaled_data != 0) > 0, colSums(unscaled_data != 0) > 0]
    unscaled_data = t(round(as.matrix(unscaled_data),4) )
    print(glue("transposed data dim: {dim(unscaled_data)}"))
    print("outputting data for cNMF")
    write.table(unscaled_data, file=glue("./{method}_{num_cells}_{seed}_{data_param}_{dim}.tab"), append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
    
} else {
    print("performing PCA")
    # CHNGE VARIABLE NAME
    sub <- RunPCA(sub, npcs = dim)

    embeddings = Embeddings(sub, reduction = "pca")
    loadings = Loadings(sub, reduction = "pca")

    print(glue("dim emb {dim(embeddings)}"))
    print(glue("dim load {dim(loadings)}"))

    output <- list()
    output[["embeddings"]] = embeddings
    output[["loadings"]] = loadings
    
    #saveRDS(output, file = glue("/data/gusev/SCRNA/HRVATIN/downsampling/QC/vargenes_20k/PCA/{method}_cells_{num_cells}_seed_{seed}_{data_param}.rds"))
}

if (method != "cNMF"){
    # filter metadata to match embeddings
    filtered_meta = meta[rownames(embeddings) ,]
    # extract numbers from simulation vector to obtain time (in hours) THEN convert to numeric
    filtered_meta$time = as.numeric(substr((filtered_meta$stim),1,1))

    # compute adjusted r squared
    r_sq = summary(lm(scale(filtered_meta$time) ~ scale(embeddings[,1:dim(embeddings)[2]])) )$adj.r.sq

    # compute correlation per PC
    corr = rep(NA, dim(embeddings)[2])
    for (i in 1:dim(embeddings)[2]){
       corr[i] = summary(lm(scale(filtered_meta$time) ~ scale(embeddings[,i])) )$adj.r.sq    

    }

    df = data.frame(method= method, data_param = data_param, num_cells = num_cells, seed = seed, adj.rsq = r_sq, max.rsq =max(corr)  )

    # append run to output dataframe
    write.table(df, file = "/PHShome/sv433/Gusev_Lab/Light_Hvratin/output.txt", append = TRUE, quote = FALSE, sep = " ",
                row.names = FALSE,
                col.names = FALSE)
}

