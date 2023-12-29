library(Seurat)
library(glue)

# read in cell type proportions
ct_prop = read.table("/data/gusev/USERS/scarver/Gusev_Lab/Simulation/Gamma-Poisson/ct_proportions.txt",header=TRUE, row.names=1)

# read in overdispersion parameters
output = readRDS("/data/gusev/USERS/scarver/Gusev_Lab/Simulation/Gamma-Poisson/fitted_parameters_Gamma_poiss.RDS")
overdispersion = data.frame(output$overdispersion)
Mus = data.frame(output$Mu)

args = commandArgs(trailingOnly=TRUE)

seed = 12345
set.seed(seed)

cells= as.integer(args[1])
N = as.integer(1.05 * cells)
#2526 total genes
genes=as.integer(args[2]) # num genes
M = as.integer(1.5 * genes)
cell_counts = data.frame(t(as.integer(ct_prop$prop * N ) ) )
colnames(cell_counts) = rownames(ct_prop)

# increasing the number of genes from ~2500 to M
ct = as.list(colnames(cell_counts) )
sample_from <- function(CT, Mus_or_overdispersion, missing_genes) {
    sample_output = sample(Mus_or_overdispersion[, CT], size = missing_genes, replace = TRUE)
    return(matrix(sample_output))
}

# add or subtract genes if the number of simulated genes is greater than the input
if (M > nrow(overdispersion)){
    # missing genes 
    missing_genes = M - nrow(overdispersion)
    sample_overdisp_matrices = lapply(ct, sample_from, Mus_or_overdispersion = overdispersion, missing_genes = missing_genes)
    sample_overdisp <- do.call(cbind, sample_overdisp_matrices)
    sample_Mu_matrices = lapply(ct, sample_from, Mus_or_overdispersion = Mus, missing_genes = missing_genes)
    sample_Mu <- do.call(cbind, sample_Mu_matrices)
    # add rownames to sample_Mu and sample_overdisp
    
    colnames(sample_overdisp) = colnames(sample_Mu) = ct

}
Mus = rbind(Mus,sample_Mu)
overdispersion = rbind(overdispersion,sample_overdisp)

rownames(Mus) = rownames(overdispersion) = paste0("gene", 1:nrow(overdispersion))

sample_matrix <- function(ct){
    cell_counts[,ct]
    n_samples = cell_counts[,ct]
    means = Mus[,ct]
    sizes = 1/overdispersion[,ct]

    samples_matrix <- matrix(rnbinom(n_samples * length(means), size = rep(sizes, each = n_samples), mu = rep(means, each = n_samples)), ncol = length(means))
    return(samples_matrix)
}

matrices_list <- lapply(ct, sample_matrix)
counts_matrix <- do.call(rbind, matrices_list)
rownames(counts_matrix) <- paste0("cell_", 1:dim(counts_matrix)[1])
colnames(counts_matrix) <- rownames(Mus)

# subset to columns and rows that are not all zero
counts_matrix = counts_matrix[rowSums(counts_matrix)  != 0,colSums(counts_matrix)  != 0]

# Find columns with zero variance
zero_variance_columns <- apply(counts_matrix, 2, function(col) var(col) <= 0.02)

# Remove columns with zero variance
counts_matrix <- counts_matrix[, !zero_variance_columns]

counts_matrix <- as.matrix(t(counts_matrix) )

counts_matrix <- counts_matrix[1:genes,1:cells]

print(dim(counts_matrix))
# Convert the dataframe to a named vector
freq_vector <- unlist(cell_counts)
# Create a vector repeating each element by its frequency
result_vector <- rep(names(freq_vector), times = freq_vector)
metadata_ct = data.frame(result_vector)
num_cells_sim = rowSums(cell_counts)
rownames(metadata_ct) = paste0("cell_", 1:num_cells_sim)
metadata_ct = metadata_ct[colnames(counts_matrix),]
metadata_ct = data.frame(metadata_ct)
rownames(metadata_ct) = colnames(counts_matrix)
colnames(metadata_ct) = "cell_type"

saveRDS(metadata_ct, glue("/data/gusev/USERS/scarver/create_datasets/data/METADATA_cells{cells}_genes_{genes}_seed_{seed}.rds") )

saveRDS(counts_matrix, glue("/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_{seed}.rds") )
