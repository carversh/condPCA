import os
import psutil
import time
import argparse

seed_value = 12345  # You can use any integer value as the seed
parser = argparse.ArgumentParser()
parser.add_argument('arg1', help='Number of cells', type=str)
parser.add_argument('arg2', help='Number of genes', type=str)

args = parser.parse_args()

cells = args.arg1
genes = args.arg2
print(cells)
print(genes)

from condPCA_comp import *
test = condPCA(count_matrix_path=f'/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_{seed_value}.rds', metadata_path=f'/data/gusev/USERS/scarver/create_datasets/data/METADATA_cells{cells}_genes_{genes}_seed_12345.rds', object_columns=['celltype'], save_image_outputs=False) # object_columns are columns that must be factors
#start time
test.Normalize()
test.Standardize()
test.CondPCA_fit()

# retrieve used memory information
current_memory_usage_bytes = psutil.virtual_memory().used
print("memory_used_bytes:", current_memory_usage_bytes)

print(test.CondPCA_gene_loadings.shape)
print(test.CondPCA_cell_embeddings.shape)
