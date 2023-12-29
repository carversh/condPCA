library(Seurat)
library(glue)

args = commandArgs(trailingOnly=TRUE)

cells = args[1]
genes = args[2]
seed = 12345
counts_matrix <- readRDS(glue("/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_{seed}.rds")) 

metadata = readRDS(glue("/data/gusev/USERS/scarver/create_datasets/data/METADATA_cells{cells}_genes_{genes}_seed_{seed}.rds"))
metadata$cell_type = as.factor(metadata$cell_type)

dim=50


set.seed(seed)


seurat_object <- CreateSeuratObject(counts = counts_matrix )
seurat_object = AddMetaData(object = seurat_object, metadata = metadata, col.name = "cell_type")
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = length(rownames(seurat_object)))
var_feat = seurat_object@assays$RNA@var.features
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000,vars.to.regress = "cell_type") 
seurat_object <- ScaleData(seurat_object, features = var_feat , scale.max=9999999999, vars.to.regress="cell_type")
seurat_object <- RunPCA(seurat_object, npcs=dim)

library(pryr)
options(scipen = 999)
current_memory_usage_mb <- pryr::mem_used()
current_memory_usage_bytes = current_memory_usage_mb * 1e6
cat("memory_used_bytes:", current_memory_usage_bytes, "\n")
