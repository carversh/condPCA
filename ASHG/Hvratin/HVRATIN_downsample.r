library(Seurat)
library(glue)

args = commandArgs(trailingOnly=TRUE)

method = args[1] # c("scaled_NMF", "PCA", "Cond_PCA")
data_param = args[2] # c("all_cell_types","int_ex")
num_cells = args[3] # c("all", # cells to downsample)
seed = as.integer(args[4])
dim = as.integer(args[5]) # rank of matrix in DR method
vargenes = as.integer(args[6])

print(glue("method: {method}"))
print(glue("data_param: {data_param}"))
print(glue("num_cells: {num_cells}"))
print(glue("seed: {seed}"))
print(glue("dim: {dim}"))
print(glue("vargenes: {vargenes}"))

# load NMF library if necessary
if (method == "NMF" | method == "scaled_NMF"){
    library(NMF)
}

# load data
load("~/Gusev_Lab/Light_Hvratin/data/GSE102827_MATRIX.RDat")

genes.keep = readRDS("~/Gusev_Lab/Light_Hvratin/data/genes_names_QC.rds")
cells.keep = readRDS("~/Gusev_Lab/Light_Hvratin/data/cell_names_QC.rds")

# subset metadata and expression data by QC'd genes and cells
filtered_exp = exp[genes.keep,cells.keep ]
print(glue("QC'd data: {dim(filtered_exp)}"))

filtered_meta = meta[cells.keep,]
print(glue("dim meta: {dim(filtered_meta)}"))


exp_seu <- CreateSeuratObject(counts = filtered_exp, meta.data = filtered_meta)

# SUBSETTING TO INT/EX IF NECESSARY
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

print(dim(sub))

# FIND VARIABLE FEATURES
sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = vargenes)
print('done variable features')
var_feat = sub@assays$RNA@var.features
print(glue("length of variable features: {length(var_feat)}"))

# PERFORM DOWNSAMPLING IF NECESSARY 
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

print(dim(sub))

if (method == "PCA"){
    sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
    sub <- ScaleData(object = sub, features=var_feat, scale.max=9999999999999, assay='RNA')
    sub <- RunPCA(object = sub, npcs = dim, verbose = FALSE, assay = 'RNA')
    embeddings = Embeddings(sub, reduction = "pca")
    final_meta = sub@meta.data[rownames(embeddings),]

    print(glue("dim emb {dim(embeddings)}"))  
    print(glue('dim emb {dim(Loadings(sub, reduction = "pca"))}'))
    
} else if (method == "Cond_PCA") {
    sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
    sub@meta.data$maintype = as.factor(sub@meta.data$maintype)
    print("conditioning")
    sub <- ScaleData(object = sub, features=var_feat, vars.to.regress = "maintype", scale.max=9999999999999, assay='RNA')
    sub <- RunPCA(object = sub, npcs = dim, verbose = FALSE, assay = 'RNA')
    embeddings = Embeddings(sub, reduction = "pca")
    final_meta = sub@meta.data[rownames(embeddings),]
    
    print(glue("dim emb {dim(embeddings)}"))  
    print(glue('dim emb {dim(Loadings(sub, reduction = "pca"))}'))
    
} else if (method == "UMAP") {
    sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
    sub <- ScaleData(object = sub, features=var_feat, scale.max=9999999999999, assay='RNA')
    sub <- RunPCA(object = sub, npcs = dim, verbose = FALSE, assay = 'RNA')
    sub <- RunUMAP(sub, dims = 1:10)
    embeddings = sub[["umap"]]@cell.embeddings
    final_meta = sub@meta.data[rownames(embeddings),]
    
    print(glue("dim emb {dim(embeddings)}"))
       
} else if (method == "scaled_NMF"){
    sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
    sub <- ScaleData(object = sub, features=var_feat, scale.max=9999999999999, assay='RNA')
    scaled_data = sub@assays$RNA@scale.data
    # replace negative values with 0 
    scaled_data[scaled_data<0] <- 0
    # remove rows or columns with all zeros
    scaled_data = scaled_data[rowSums(scaled_data != 0) > 0, colSums(scaled_data != 0) > 0]
    scaled_data = t(round(as.matrix(scaled_data),4) )
    print("performing scaled NMF")
    out = nmf(scaled_data, rank=dim)
    w_cells <- out@fit@W
    h_genes <- out@fit@H

    output <- list()
    output[["embeddings"]] = w_cells
    output[["loadings"]] = h_genes
    embeddings = output[["embeddings"]]
    final_meta = sub@meta.data[rownames(embeddings),]
    
    print(glue("dim emb {dim(embeddings)}"))
    print(glue('dim load {dim(output$loadings)}'))
    
} else if (method == "NMF"){
    print("performing standard NMF")  
    unscaled_data = sub@assays$RNA@counts[var_feat,]
    print(glue("dim input NMF: {dim(unscaled_data)}"))
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
    final_meta = sub@meta.data[rownames(unscaled_data),]
    print(dim(final_meta))
    
    print(glue("dim emb {dim(embeddings)}"))
    print(glue('dim load {dim(output[["loadings"]])}'))
    
} else { # method should be cNMF
    print(method)
    print("obtaining data for cNMF")
    unscaled_data = sub@assays$RNA@counts[var_feat,]
    # remove rows or columns with all zeros
    unscaled_data = unscaled_data[rowSums(unscaled_data != 0) > 0, colSums(unscaled_data != 0) > 0]
    unscaled_data = t(round(as.matrix(unscaled_data),4) )
    print(glue("transposed data dim: {dim(unscaled_data)}"))
    print("outputting data for cNMF")
    write.table(unscaled_data, file=glue("./{method}_{num_cells}_{seed}_{data_param}_{dim}_{vargenes}.tab"), append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
}

if (method != "cNMF"){
    print(glue("embedding dim: {dim(embeddings)}"))
    # obtain time metadata
    final_meta$time = as.numeric(substr((final_meta$stim),1,1))

    # compute adjusted r squared
    r_sq = summary(lm(scale(final_meta$time) ~ scale(embeddings[,1:dim(embeddings)[2]])) )$adj.r.sq

    # compute correlation per PC
    corr = rep(NA, dim(embeddings)[2])
    for (i in 1:dim(embeddings)[2]){
       corr[i] = summary(lm(scale(final_meta$time) ~ scale(embeddings[,i])) )$adj.r.sq    

    }

    df = data.frame(method= method, data_param = data_param, num_cells = num_cells, seed = seed, dim=dim, vargenes=length(var_feat), adj.rsq = r_sq, max.rsq =max(corr))

    # append run to output dataframe
    write.table(df, file = "/PHShome/sv433/Gusev_Lab/Light_Hvratin/output_final_format.txt", append = TRUE, quote = FALSE, sep = " ", row.names = FALSE,col.names = FALSE)
}

df = data.frame(method=method,data_param=data_param,num_cells=num_cells,seed=seed, dim=dim,vargenes=vargenes )
write.table(df, file = "/PHShome/sv433/Gusev_Lab/Light_Hvratin/finished_downsampling_scripts.txt", append = TRUE, quote = FALSE, sep = " ", row.names = FALSE,col.names = FALSE)
