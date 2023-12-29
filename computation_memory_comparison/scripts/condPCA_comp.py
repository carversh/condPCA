#!/usr/bin/python

import pandas as pd
import numpy as np
import scanpy
import scipy
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pyreadr
import anndata

class condPCA(object):
    def __init__(self, count_matrix_path, metadata_path, object_columns, vars_to_regress=True, n_PCs=50, random_seed=12345, vargenes_IterPCA=500, vargenes_Stand_Cond="all", save_image_outputs = False, BIC=True):
        """
        Parameters
        ----------
        count_matrix:
            Count matrix that must be log-normalized and standardized

        metadata:
            metadata containing cell type labels named "celltype"

        object_columns:
            columns that will be one hot encoded/columns that are factors

        vars_to_regress:
            list of variables to regress out

        """
        
        #self.count_matrix = scanpy.read(count_matrix_path) # cells x genes, pd.read_csv(count_matrix_path, sep='\t', header=0, index_col=0)
       # self.count_matrix = anndata.AnnData(pyreadr.read_r(count_matrix_path)[None].T)
        try:
            self.count_matrix = anndata.AnnData(pyreadr.read_r(count_matrix_path)[None].T)
        except:
            print("reading txt file instead of rds")
            # count_matrix_path = f'/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_{seed_value}.txt'
            count_matrix_path = count_matrix_path.replace(".rds", ".txt")
            self.count_matrix = anndata.AnnData(pd.read_csv(count_matrix_path, sep='\t', index_col=0).T)        

        # create a separate count matrix that won't be log normalized for iterative PCA
        #self.iter_count_matrix = scanpy.read(count_matrix_path)
       
        #self.metadata = pd.read_csv(metadata_path, sep='\t', header=0, index_col=0)
        self.metadata = pyreadr.read_r(metadata_path)[None]
        self.metadata = self.metadata.rename(columns={'cell_type': 'celltype'})
        
        if vars_to_regress:
            
            self.vars_to_regress = self.metadata.columns
        
        else: # if vars_to_regress is a list, convert to pandas core Index object
            
            self.vars_to_regress = pd.Index(vars_to_regress)

        # one hot encode necessary metadata variables
        self.object_columns = object_columns # obtain columns that must be one hot encoded
        
        self.metadata[self.object_columns] = self.metadata[self.object_columns].astype(object) # convert these columns to objects

        self.random_seed = random_seed # set random seed
        
        self.n_PCs = n_PCs # number of PCs to extract
        
        self.vargenes_IterPCA = vargenes_IterPCA # number of variable genes for iterative pca
        
        if vargenes_Stand_Cond == "all":
            self.vargenes_Stand_Cond = self.count_matrix.shape[1]
        else:
            self.vargenes_Stand_Cond = vargenes_Stand_Cond # number of variable genes for standard and conditional pca

        self.save_image_outputs = save_image_outputs # save image outputs
        
        self.BIC = BIC # compute BIC

    def Normalize(self):
        """
        Normalize and take log1p of count data (for conditional and standard, not iterative)
        """
        
        # update scanpy object to normalize all rows, so every cell sums to 10k
        scanpy.pp.normalize_total(self.count_matrix, target_sum = 10000)
       
        # log transform
        scanpy.pp.log1p(self.count_matrix)
        
        # plot mean variance relationship if specified by user
        if self.save_image_outputs:
            self._plot_mean_variance_relationship(self.count_matrix.X, label="All Cells")

        # find highly variable genes
        scanpy.pp.highly_variable_genes(self.count_matrix, n_top_genes=self.vargenes_Stand_Cond)

    def Standardize(self):
        """
        Standardize count data AND metadata (for conditional and standard, not iterative)
        """
        
        # Standardize count data

        # if matrix is sparse, make dense
        if scipy.sparse.issparse(self.count_matrix.X):
            
            self.count_matrix.X = self.count_matrix.X.todense()
        
        # only subset the matrix to the most variable genes
        self.standardized_count_data = self._standardize(self.count_matrix.X[:, self.count_matrix.var['highly_variable']])

        # Process metadata/covariates for standardization:

        # subset to only variables that you want to regress out
        self.metadata = self.metadata[self.vars_to_regress]
       
        # WARNING IN FOLLOWING LINE BECAUSE CONVERTING OBJECT THAT LOOKS NUMERIC TO BE ONE HOT ENCODED, this is batch
        self.IterPCA_metadata = pd.get_dummies(self.metadata, drop_first=False)
        
        # Convert factor covariates to dummy variables dropping one column
        self.metadata = pd.get_dummies(self.metadata, drop_first=True)
        
        self.standardized_metadata = self._standardize(self.metadata)
    
    def _standardize(self, mat): # simple function performing standardization
        
        # compute means of genes or covariates
        mean_vector = np.mean(mat, axis=0)
       
       # compute standard deviation of genes or covariates
        std_vector = np.std(mat, axis=0)
        
        
        if sum(std_vector == 0) > 0 :
            indices_std_0 = np.where(std_vector == 0)[0]
            print("WARNING: some std values are 0")
            print(indices_std_0)
            arr_without_column = np.delete(mat, indices_std_0, axis=1)
            # recompute means and std
            mean_vector = np.mean(arr_without_column, axis=0)
            std_vector = np.std(arr_without_column, axis=0)
            stand_mat = (arr_without_column - mean_vector) / std_vector
        else:
            # standardize by gene or covariates
            stand_mat = (mat - mean_vector) / std_vector
        return stand_mat
    
    def _regress_covariates(self, standardized_metadata, standardized_count_data): # function regressing set of covariates
        
        # append ones to standardized meta for intercept
        standardized_metadata = np.c_[np.ones((standardized_metadata.shape[0], 1)), standardized_metadata]
        
        # compute inverse of np.matmul(A^T, A) where A is the standardized metadata or covariates
        inv_cov = np.linalg.pinv(np.matmul(standardized_metadata.T, standardized_metadata) )
        
        # compute betas per gene
        betas = inv_cov @ standardized_metadata.T @ standardized_count_data
        
        # compute prediction
        prediction = np.matmul(standardized_metadata, betas) # compute prediction
        
        # compute residual
        residual = standardized_count_data - prediction
        
        standardized_residual = self._standardize(residual)
        
        return standardized_residual
    
    def _fit_pca(self, mat, standardPCA, condPCA, iterPCA, iterPCA_genenames, iterPCA_cellnames, iterPCA_CT): # fitting PCA
        
        # instantiate PCA with hyperparameters
        pca = PCA(n_components=self.n_PCs, random_state=self.random_seed)
        
        # projections (of input data onto eigenvectors)
        pca.fit(mat)
       
        # retrieve eigenvectors/gene loadings
        gene_loadings = pca.components_

        # retrive cell embeddings
        cell_embeddings = pca.transform(mat)

        # retrieve eigenvalues
        eigenvalues = pca.explained_variance_
        
        # if iterative PCA
        if iterPCA:
            
            # convert gene loadings to dataframe
            gene_loadings = pd.DataFrame(gene_loadings.T, index = list(iterPCA_genenames ), columns = [f'PC_{i}' for i in range(1, (gene_loadings.T.shape[1]+1))])
                
            # convert cell embeddings to dataframe
            cell_embeddings = pd.DataFrame(cell_embeddings, index = list(iterPCA_cellnames), columns = [f'PC_{i}' for i in range(1, (cell_embeddings.shape[1]+1))])

        # convert eigenvalues to dataframe
        eigenvalues = pd.DataFrame(eigenvalues, index = [f'PC_{i}' for i in range(1, (eigenvalues.shape[0]+1))], columns=["Eigenvalues"])

        # if Standard or Conditional PCA, construct dataframes based on gene and cell list from original count matrix
        if not iterPCA:

            # convert gene loadings to dataframe
            gene_loadings = pd.DataFrame(gene_loadings.T, index = list(self.count_matrix.var_names[self.count_matrix.var['highly_variable']] ), columns = [f'PC_{i}' for i in range(1, (gene_loadings.T.shape[1]+1))])
                
            # convert cell embeddings to dataframe
            cell_embeddings = pd.DataFrame(cell_embeddings, index = list(self.count_matrix.obs_names), columns = [f'PC_{i}' for i in range(1, (cell_embeddings.shape[1]+1))])
      
        if self.BIC:
            if standardPCA:
                # compute BIC
                min_BIC_index = self._compute_BIC(eigenvalues, self.standardized_residual, "Standard PCA")

            if condPCA:
                # compute BIC
                min_BIC_index = self._compute_BIC(eigenvalues, self.standardized_residual, "CondPCA")

            if iterPCA:
                # compute BIC
                min_BIC_index = self._compute_BIC(eigenvalues, self.standardized_residual, "IterPCA", iterPCA_CT)

            BIC_cutoff = "PC_" + (min_BIC_index + 1).astype(str)
            return cell_embeddings, gene_loadings, eigenvalues, BIC_cutoff

        return cell_embeddings, gene_loadings, eigenvalues, "Not Calculated"
    
    def _fit_model(self, standardized_metadata, standardized_count_data, standardPCA=False, condPCA = False, iterPCA=False, iterPCA_genenames=False, iterPCA_cellnames=False, iterPCA_CT=False): # regress out covariates and then input into PCA

        # regress out covariates (including celltype) and retrieve standardized residual
        self.standardized_residual = self._regress_covariates(standardized_metadata = standardized_metadata, standardized_count_data= standardized_count_data) # REMOVE SELF
        
        #np.save('/Users/shayecarver/condPCA/final_method/standardized_residual.npy', standardized_residual) # THIS CAN BE DELETED, IT WAS FOR DEBUGGING
        
        # if not iterative PCA, able to add gene names and cell names here, but must subset if IterPCA
        if not iterPCA:
            # return standardized residual as a dataframe with gene and cell names:
            self.standardized_residual = pd.DataFrame(self.standardized_residual, index = list(self.count_matrix.obs_names), columns = list(self.count_matrix.var_names[self.count_matrix.var['highly_variable']]))# REMOVE SELF

        if iterPCA:
            # return standardized residual as a dataframe with gene and cell names of the given subset:
            self.standardized_residual = pd.DataFrame(self.standardized_residual, index = list(iterPCA_cellnames), columns = list(iterPCA_genenames))# REMOVE SELF

        
        # perform PCA on residualized matrix
        return self._fit_pca(self.standardized_residual, standardPCA=standardPCA, condPCA=condPCA, iterPCA=iterPCA, iterPCA_genenames=iterPCA_genenames, iterPCA_cellnames=iterPCA_cellnames, iterPCA_CT=iterPCA_CT)# REMOVE SELF

    def _mapping_IterPCA_subset_dataframes_to_PCA(self, metadata, CT_exp_column): # function that subsets count matrix to a particular cell type and then performs PCA on that subset
        # remove "celltype_" from the string CT_exp_column
        CT_label = CT_exp_column.replace('celltype_', '')

        # extract indices of the cells that belong to the particular cell type of interest (indicated by CT_column, which is a column name)
        indices_given_ct = self.dataframe_CT_indices[CT_exp_column]

        # Check if the sum of indices is less than or equal to 200
        if sum(indices_given_ct) <= 200:
            # Return empty output
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), None
 
        # subset the count data to the cells belonging to the cell type
        metadata_subset_to_CT = metadata[indices_given_ct]

        # Re-process from log-normalized data to standadization of the matrix (find new set of variable genes for the subset)
        
        # make a tmp copy and subset the matrix to cells in the particular cell type identify highly variable genes from log normalized count matrix
        # Subset the Scanpy object to the specified row/cell indices
        tmp_copy_counts = self.iter_count_matrix[indices_given_ct].copy()

        # update scanpy object to normalize all rows, so every cell sums to 10k
        scanpy.pp.normalize_total(tmp_copy_counts, target_sum = 10000)
       
        # log transform
        scanpy.pp.log1p(tmp_copy_counts)
        
        # plot mean variance relationship if specified by user
        if self.save_image_outputs:

            # if white space, replace with underscore
            CT_label = CT_label.replace(' ', '_')

            self._plot_mean_variance_relationship(tmp_copy_counts.X, label=CT_label)
        
        # find highly variable genes after subsetting count matrix to the cells in a particular cell type
        scanpy.pp.highly_variable_genes(tmp_copy_counts, n_top_genes=self.vargenes_IterPCA)

        # subset the scanpy object to the most variable genes
        tmp_copy_counts_subset = tmp_copy_counts[:, tmp_copy_counts.var['highly_variable']]
        
        # Re-standardize count databecause it has just been subset
        count_data_subset_to_CT = self._standardize(tmp_copy_counts_subset.X)
        
        # Re-standardize metadata because it has just been subset
        metadata_subset_to_CT = self._standardize(metadata_subset_to_CT)

        # extract the cell names or barcodes of the cells belonging to the cell type of interest
        cellnames = tmp_copy_counts_subset.obs_names

        # extract the gene names of the genes belonging to the most variable genes within that subset
        genenames = tmp_copy_counts_subset.var_names

        # fit the given model by regressing out covariates and performing PCA
        return self._fit_model(standardized_metadata=metadata_subset_to_CT,standardized_count_data = count_data_subset_to_CT, iterPCA=True, iterPCA_genenames=genenames, iterPCA_cellnames = cellnames, iterPCA_CT=CT_label)

    
    def CondPCA_fit(self):
       
        # fit linear model (regress out covariates) and fit PCA -- covariates contain cell type
        self.CondPCA_cell_embeddings, self.CondPCA_gene_loadings, self.CondPCA_eigenvalues, self.CondPCA_BIC_cutoff = self._fit_model(standardized_metadata=self.standardized_metadata,standardized_count_data= self.standardized_count_data, condPCA = True)


    def StandardPCA_fit(self):
        
        # remove celltype from covariate space
        standardized_metadata_minus_celltype = self.standardized_metadata.drop(columns = self.standardized_metadata.filter(like="celltype", axis=1).columns )
        
        # fit linear model (regress out covariates) and fit PCA -- covariates do not contain cell type
        self.StandardPCA_cell_embeddings, self.StandardPCA_gene_loadings, self.StandardPCA_eigenvalues, self.StandardPCA_BIC_cutoff = self._fit_model(standardized_metadata=standardized_metadata_minus_celltype,standardized_count_data= self.standardized_count_data,standardPCA=True)

    def Iter_PCA_fit(self):

        # remove celltype from standardized covariate space
        metadata_minus_celltype = self.standardized_metadata.drop(columns = self.standardized_metadata.filter(like="celltype", axis=1).columns )
        
        # get dataframe with boolean indices for each cell type
        self.dataframe_CT_indices = self.IterPCA_metadata.filter(like="celltype", axis=1).astype(bool)
        
        # get the name of the columns that indicate a cell type
        celltype_colnames = self.dataframe_CT_indices.columns

        # create empty dictionaries to store results of iterative PCA per cell type
        self.IterPCA_cell_embeddings = {}
        self.IterPCA_gene_loadings = {}
        self.IterPCA_eigenvalues = {}
        self.IterPCA_BIC_cutoff = {}

        # iterate through each cell type and perform iterative PCA, storing results in dictionaries
        for celltype_column in celltype_colnames:
            # obtain cell type name, replace spaces with underscores
            tmp_CT = celltype_column.replace("celltype_", "").replace(" ", "_")
            tmp_result = self._mapping_IterPCA_subset_dataframes_to_PCA(metadata_minus_celltype, celltype_column)
            # append results to appropriate dictionary
            self.IterPCA_cell_embeddings[tmp_CT] = tmp_result[0]
            self.IterPCA_gene_loadings[tmp_CT] = tmp_result[1]
            self.IterPCA_eigenvalues[tmp_CT] = tmp_result[2]
            self.IterPCA_BIC_cutoff[tmp_CT] = tmp_result[3]

    def _plot_mean_variance_relationship(self, log_normed_data, label):
        # compute the mean of every gene/column of the matrix
        mean_vector = np.mean(log_normed_data, axis=0)
        # compute the variance of every gene/column of the matrix
        var_vector = np.var(log_normed_data, axis=0)
        plt.scatter(mean_vector, var_vector)
        plt.xlabel("Mean")
        plt.ylabel("Variance")
        plt.title(f'Mean-Variance Relationship ({label}, {log_normed_data.shape[0]} Cells)')
        # add dashed red line with slope 1 and y intercept 0
        plt.plot([0, 5], [0, 5], 'r--')
        # save plot to current directory
        plt.savefig(f'Mean_Variance_Relationship_{label}.png')
        # close plot
        plt.close()

    def _compute_BIC(self, eigenvalues, standardized_residual, label, ct=False):
        
        # extarct the standardized residuals as a numpy array
        X = standardized_residual.values

        # Compute the sample covariance matrix
        cov_matrix = (X.T @ X) / (X.shape[0] - 1)

        # Compute the trace of the sample covariance matrix (this is equal to the sum of the eigenvalues)
        trace = np.trace(cov_matrix)

        # compute BIC

        p = X.shape[1]
        n = X.shape[0]

        # Initialize an array to store BIC values
        BIC_values = np.zeros(self.n_PCs)

        # Perform calculations for each value of j
        for j in range(0, self.n_PCs ):

            ell_bar = (trace - np.sum(eigenvalues.iloc[0:j, 0]) ) / (p - j)

            dj = (j + 1) * (p + 1 - j / 2)

            term_1 = n * np.log(np.prod(eigenvalues.iloc[0:j, 0]))

            term_2 = n * (p - j) * np.log(ell_bar)

            term_3 = np.log(n) * dj

            term_4 = n * (np.log((n - 1) / n )**p) + n * p * (1 + np.log(2 * np.pi))

            BIC_j = term_1 + term_2 + term_3 + term_4

            # Store BIC value in the array
            BIC_values[j] = BIC_j

        # Find the index corresponding to the minimum BIC value
        min_BIC_index = np.argmin(BIC_values)
        #min_BIC_value = BIC_values[min_BIC_index]

        if self.save_image_outputs:
            # Plot BIC values
            plt.plot(range(1, self.n_PCs + 1), BIC_values, marker='o', linestyle='-', color='b')
            plt.axvline(x=min_BIC_index + 1, color='r', linestyle='--', label=f'Min BIC at $j = {min_BIC_index + 1}$')
            plt.xlabel('Number of Principal Components (j)')
            plt.ylabel('BIC Value')
            
            plt.legend()
            plt.grid(True)
            if not ct:
                plt.title(f'{label} BIC Values vs. Number of Principal Components')
                plt.savefig(f'BIC_plot_{label}.png')
            if ct:
                plt.title(f'{label} BIC Values vs. Number of Principal Components ({ct})')
                plt.savefig(f'BIC_plot_{label}_{ct}.png')
            plt.close()

        return min_BIC_index
        
def main():
    pass  # Empty main function

if __name__ == "__main__":
    main()
