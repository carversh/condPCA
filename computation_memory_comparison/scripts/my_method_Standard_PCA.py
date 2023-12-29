import pandas as pd
import numpy as np
import scanpy
from sklearn.decomposition import PCA
import pyreadr
import anndata
import argparse
import os
import psutil
import time


seed_value = 12345  # You can use any integer value as the seed
parser = argparse.ArgumentParser()
parser.add_argument('arg1', help='Number of cells', type=str)
parser.add_argument('arg2', help='Number of genes', type=str)

args = parser.parse_args()

cells = args.arg1
genes = args.arg2
count_matrix_path = f'/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_{seed_value}.rds'

try:
    count_matrix = anndata.AnnData(pyreadr.read_r(count_matrix_path)[None].T)
except:
    print("reading txt file instead of rds")
    count_matrix_path = f'/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_{seed_value}.txt'
    count_matrix = anndata.AnnData(pd.read_csv(count_matrix_path, sep='\t', index_col=0).T)

scanpy.pp.normalize_total(count_matrix, target_sum = 10000)
# log transform
scanpy.pp.log1p(count_matrix)

# compute means of genes or covariates
mean_vector = np.mean(count_matrix.X, axis=0)

# compute standard deviation of genes or covariates
std_vector = np.std(count_matrix.X, axis=0)

if sum(std_vector == 0) > 0 :
    indices_std_0 = np.where(std_vector == 0)[0]
    print("WARNING: some std values are 0")
    print(indices_std_0)
    arr_without_column = np.delete(count_matrix.X, indices_std_0, axis=1)
    # recompute means and std
    mean_vector = np.mean(arr_without_column, axis=0)
    std_vector = np.std(arr_without_column, axis=0)
    stand_mat = (arr_without_column - mean_vector) / std_vector
else:
    # standardize by gene or covariates
    stand_mat = (count_matrix.X - mean_vector) / std_vector

# instantiate PCA with hyperparameters
pca = PCA(n_components=50, random_state=seed_value)

# projections (of input data onto eigenvectors)
pca.fit(stand_mat)

# retrieve used memory information
current_memory_usage_bytes = psutil.virtual_memory().used
print("memory_used_bytes:", current_memory_usage_bytes)

# retrieve eigenvectors/gene loadings
gene_loadings = pca.components_

