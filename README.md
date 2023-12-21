# condPCA Package v1.0.0

condPCA is an analysis pipeline for inferring cell states within and across cell types by leveraging cell type labels. 

It takes a count matrix (N cells X G genes) as input and produces ...

![Inline Image](https://github.com/carversh/condPCA/blob/main/images/pipeline_visual.png?raw=true)

# Installation

```
pip install CondPCA
```

# Running CondPCA
test

# Step by Step Guide (Quickstart)

### Step 1 - initialize a class with your input which consists of a (1) QC'd counts matrix and a (2) covariates matrix containing a column titled "celltype" with respective cell type labels of each cell

Example command:

```
scExp = condPCA(count_matrix_path="matrix.txt", metadata_path="metadata.txt", object_columns=['Batch', 'Sex','celltype'], save_image_outputs=False, BIC=True)

```

Input Data
  - count_matrix_path - points to the tab delimited file that is cells by genes dimensional, where rownames correspond to the barcoded cells and columnnames correspond to the genes, or a Scanpy object.
  - metadata_path - points to the tab delimited file that is cells by covariates/features dimensional.There must be a single column with a column name "celltype" that contains the cell type labels corresponding to each barcoded cell in the count matrix.

Parameters
  - object_columns - A list that specifies column names whose column should be treated as a factor or object and not a numeric column. The method will one hot encode these columns. Celltype must be in this list. Commonly, batch is another covariate that is included in this list.
  - n_PCs - number of PCs to output per method. Default is 200.
  - random_seed - Defailt is 998999.
  - vargenes_Stand_Cond - The number of variable genes to analyze during Standard and Conditional PCA. This can be set to an integer or "all" if all genes should be considered for analyses.
  - vargenes_IterPCA - The number of variable genes to analyze during Iterative PCA. This can be set to an integer or "all" if all genes should be considered for analyses.
  -  
  - save_image_outputs - True/False. Specifies whether to save the graphical outputs during the pipeline processing (i.e. Mean-Variance Log-Normalization Graphs, BIC Cutoff Graphs,...)
  - BIC - True/False. Specified whether to compute the Bayesian Information Criterion after each method so that the set of states has a statistical cutoff.
### Step 2 - log-normalize the count data

Example command:

```
scExp.Normalize()
```

### Step 3 - standardize the count data

Example command:

```
scExp.Standardize()
```

### Step 4 - perform Standard PCA

Example command:

```
scExp.StandardPCA_fit()
```

Outputs
  - scExp.StandardPCA_cell_embeddings - cell embeddings outputted by Standard PCA
  - scExp.StandardPCA_gene_loadings - gene loadings or eigenvectors outputted by Standard PCA
  - scExp.StandardPCA_eigenvalues - eigenvalues outputted by Standard PCA
  - scExp.StandardPCA_BIC_cutoff - PC cutoff that specifies the maximum state that is significant. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable

### Step 4 - perform Conditional PCA

Example command:

```
scExp.CondPCA_fit()
```

Outputs
  - scExp.CondPCA_cell_embeddings - cell embeddings outputted by Conditional PCA
  - scExp.CondPCA_gene_loadings - gene loadings or eigenvectors outputted by Conditional PCA
  - scExp.CondPCA_eigenvalues - eigenvalues outputted by Conditional PCA
  - scExp.CondPCA_BIC_cutoff - PC cutoff that specifies the maximum state that is significant. For significant states, subset the cell embeddings and gene loadings from PC1 to the PC specified in this variable
