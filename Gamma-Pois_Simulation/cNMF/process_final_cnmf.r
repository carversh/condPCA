#!/usr/bin/env Rscript
library(glue)
library(Seurat)
library(stringr)



# parse args
args = commandArgs(trailingOnly=TRUE)


flag_sim_100_cts = args[9]
name = args[8]


# determine which dataset to read in
if (is.na(flag_sim_100_cts)){
    print("performing simulation on standard 7 Morabito cell types")
    
    # read in cell type proportions
    ct_prop = read.table("../ct_proportions.txt",header=TRUE, row.names=1)

    # read in overdispersion parameters
    output = readRDS("/PHShome/sv433/scratch/MORABITO/RNA/fitted_parameters_Gamma_poiss.RDS")
    overdispersion = output$overdispersion
    Mus = output$Mu
    
} else {
    print("performing simulation on 100 cell types, synthetic")
    
    # read in cell type proportions
    ct_prop = read.table("/PHShome/sv433/scratch/MORABITO/RNA/ct_proportions_100_cts.txt",header=TRUE, row.names=1)
    
    # read in overdispersion parameters
    output = readRDS("/PHShome/sv433/scratch/MORABITO/RNA/fitted_parameters_Gamma_poiss_100_cts.RDS")
    overdispersion = output$overdispersion
    Mus = output$Mu
}

# parameters inputted 
seed = as.integer(args[1])
set.seed(seed)

method = args[2]

state_type = args[3] # "within_one_ct", "across_cts"

if (state_type == "within_one_ct"){
    ct_in_state = sample(rownames(ct_prop),1)
    
} else if (state_type == "cts_6") {
    ct_in_state = sample(rownames(ct_prop),6)
    
} else { #"across_cts"
    ct_in_state = unique(rownames(ct_prop) )      
}

# dimensionality of DR output
dim = as.integer(args[4])

# compute cells per cell type according to total_cells hyperparameter and proportions in real data
total_cells = as.integer(args[5])
ct_counts = as.list( ceiling(ct_prop$prop * total_cells) )
names(ct_counts) = rownames(ct_prop)

perc_genes = as.numeric(args[6])
perc_cells = as.numeric(args[7])

print(paste("seed:", seed))
print(paste("Method:", method))
print(paste("state_type:", state_type))
print(paste("dim:", dim))
print(paste("total_cells:", total_cells))
print(paste("perc_genes:", perc_genes))
print(paste("perc_cells:", perc_cells))
print(paste("ct_in_state:", ct_in_state))
print(paste("flag_sim_100_cts:", flag_sim_100_cts))

state_info = list("State_1" = c(perc_genes, perc_cells))
#state_info_ct = list("State_1" = c("MG","ODC"))
state_info_ct = list("State_1" = ct_in_state)

# total cells
N = do.call(sum, ct_counts)
# total genes
M = length(rownames(Mus) )
# gene names
genes = rownames(Mus)
# total number of states
num_state = length(names(state_info))
# create cell type label vector
ct_labels = c()
for (ct in names(ct_counts)){
    ct_labels = append(ct_labels, rep(ct, ct_counts[[ct]]) )
}

total_cells_can_occupy_state = sum(ct_labels %in% ct_in_state)
total_cells_occupy_state = floor(total_cells_can_occupy_state * perc_cells)
total_genes_in_state = floor(M * perc_genes)

# if simulating states, 
  # if state has not been mapped to speciic cell types based on input, peform this mapping
if (is.list(state_info_ct)){
    
    # if cell type-state assignment has NOT been inputted, make assignment based on number of cell types in state
#     if (is.numeric(state_info_ct[[1]]) ){
#         # iterate over states
#         for (state in names(state_info)){
#             # sample the number of cell types in the state from total cell types
#             ct_in_state = sample(unique(ct_labels), size = state_info_ct[[state]], replace = FALSE)
#             # rewrite this into state_info_ct
#             state_info_ct[[state]] = ct_in_state   
#         } 
#     } 
  # CREATE CELL-STATE INDICATOR MATRIX (sample cells in state)
    cell_state_ind = matrix(0,N,num_state)
    colnames(cell_state_ind) = names(state_info)
    
    
  # CREATE GENE-STATE INDICATOR MATRIX (sample genes in state)
    gene_state_ind = matrix(FALSE,M, num_state)
    rownames(gene_state_ind) = genes
    colnames(gene_state_ind) = names(state_info)
    
    # iterate over states to populate the CELL-STATE INDICATOR MATRIX and the GENE-STATE INDICATOR MATRIX
    for (state in names(state_info)){
        # FOR CELL-STATE INDICATOR MATRIX
        # compute total number of cells that can possibly occupy state 
        total_cells_state = length(which(ct_labels %in% state_info_ct[[state]]) )
        # sample total cells that can occupy a state based on % cells in state inputted
        cell_indices_in_state = sample(which(ct_labels %in% state_info_ct[[state]]), size = floor(total_cells_state * state_info[[state]][2]), replace = FALSE)
        # add 1 to indicator matrix for cell in state
        cell_state_ind[cell_indices_in_state, state] = 1
    #####################################################################      
        # FOR GENE-STATE INDICATOR MATRIX
        # sample genes in a given state based on % genes in state inputted
        gene_indices_in_state = sample(genes, size = floor(M * state_info[[state]][1]))
        # add 1 to indicator matrix for gene in state
        gene_state_ind[gene_indices_in_state, state] = TRUE  
    }
}


# CREATE FOLD CHANGE PER STATE
log_FC = replicate(n=num_state, {runif(n=M, 0.5, 3) })
# create vector or -1 and +1 
pos_neg = replicate(n=num_state, {rbinom(M, 1, 0.5)})
pos_neg = replace(pos_neg,pos_neg==0,-1)
FC = 2^(log_FC * pos_neg)

rownames(FC) = genes
colnames(FC) = names(state_info)
# if gene is not in state, assign FC of 0 --> this will not add additional expression to baseline
FC[!gene_state_ind] = 0


compute_y_j_i <- function(i){ # i is the cell index
    ct = ct_labels[i]
    b_j_i = Mus[,ct]
    # add cell type continuum (for now it's 1)
    W_j = matrix(1,M,1)
    # W_j = runif(M, min = 0, max = 1) # for cell type continuum

    if (!is.list(state_info)){
        #print("no states being simulated")
        y_j_i = (W_j * b_j_i)


    } else {
        #print(glue("{length(names(state_info))} state(s) being simulated") )
        
        # sample state continuums from standard uniform
        state_cont = runif(n=num_state, min = 0, max = 1)
        # W_j_z should be 1 by state
        W_j_z = state_cont * cell_state_ind[i,,drop=FALSE]
        S_j_z_i = b_j_i * FC

        # duplicate state continuum to be M x State dimensional
        state_cont_M_dimensional = matrix(rep(t(W_j_z),M),ncol=ncol(W_j_z),byrow=TRUE)
        # set genes not in the continuum to 0
        state_cont_M_dimensional[!gene_state_ind] = 0 

        sum_continuums = (W_j + matrix(apply(state_cont_M_dimensional,1,sum) ) )

        y_j_i = ((W_j * b_j_i) + (S_j_z_i %*% t(W_j_z)) ) / sum_continuums


        # continuums don't make sense because some states shouldn't be a part of the continuum so this weighting should change
        # sum of continuums is also wrong
        continuums = sweep(cbind(W_j,state_cont_M_dimensional), 1, sum_continuums, "/")
        rownames(continuums) = genes
        colnames(continuums) = c("Cell Type", names(state_info))


    }
    return(list("y_j_i"=y_j_i, "continuums"=continuums) ) 
    #return(list("y_j_i"=y_j_i, "state_cont"=state_cont) ) 
} 

# obtain y_j_i output over all cells
out = lapply(1:N, FUN=compute_y_j_i) 


# parse output into cell and state continuums
continuums = lapply(out, function(l) l[[2]])


# obtain mean continuums per cell
max_continuums = do.call("rbind", lapply(continuums, function(l) t(data.frame(apply(l,2,max)) )) )
rownames(max_continuums) = c()
###########################
# state_cont = max_continuums[,"State_1"]


# read in embeddings
embeddings = read.table(glue("./cNMF_output/{name}/{name}.usages.k_{dim}.dt_0_5.consensus.txt"),sep="\t",row.names=1,header=TRUE)


print(glue("embedding dim: {dim(embeddings)}"))

# # compute adjusted r squared
# r_sq = summary(lm(scale(state_cont) ~ scale(embeddings[,1:dim(embeddings)[2]])) )$adj.r.sq

# # compute correlation per PC
# corr = rep(NA, dim(embeddings)[2])
# for (i in 1:dim(embeddings)[2]){
#    corr[i] = summary(lm(scale(state_cont) ~ scale(embeddings[,i])) )$adj.r.sq    
# }
# print(state_info[[state]])
# print(state_info[[state]][2])

# # df = data.frame(method = method, state_type = state_type, ct_in_state = paste(ct_in_state,collapse='_'), perc_genes= perc_genes, perc_cells = perc_cells, seed = seed, adj.rsq = r_sq, max.rsq =max(corr) , total_cells_can_occupy_state = total_cells_can_occupy_state, total_cells_occupy_state = floor(total_cells_state * state_info[[state]][2]), total_genes = M, total_genes_in_state = total_genes_in_state )
                                         
# df = data.frame(seed=seed, method=method, state_type=state_type, dim=dim, total_cells=total_cells, perc_genes=perc_genes, perc_cells=perc_cells, flag_sim_100_cts=flag_sim_100_cts,ct_in_state = "all", adj.rsq = r_sq, max.rsq =max(corr), max.rsq_sub = max(sub_corr_sq) , total_cells_can_occupy_state = total_cells_can_occupy_state, total_cells_occupy_state = floor(total_cells_state * state_info[[state]][2]), total_genes = M, total_genes_in_state = total_genes_in_state)
    
# # append run to output dataframe
# write.table(df, file = "../output.txt", append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
                                         
                                         
###############
                                         

# compute adjusted r squared
r_sq = summary(lm(scale(max_continuums[,"State_1"]) ~ scale(embeddings[,1:dim(embeddings)[2]])) )$adj.r.sq

# compute squared correlation per PC
corr = rep(NA, dim(embeddings)[2])
for (i in 1:dim(embeddings)[2]){
   corr[i] = summary(lm(scale(max_continuums[,"State_1"]) ~ scale(embeddings[,i])) )$adj.r.sq    
}
# compute squared correlation when subsetting to the cell type that is undergoing the state
print(ct_in_state)
if (state_type == "within_one_ct"){

    sub_corr_sq = rep(NA, dim(embeddings)[2])

    for (i in 1:dim(embeddings)[2]){
        sub_corr_sq[i] = cor(embeddings[ct_labels == ct_in_state,i], max_continuums[,"State_1"][ct_labels == ct_in_state])^2
    }
} else {
    sub_corr_sq = corr

}
print(dim(embeddings))                          
# determine which dataset to read in
if (is.na(flag_sim_100_cts)){
    print("saving output standard 7 Morabito cell types")
    df = data.frame(seed=seed, method=method, state_type=state_type, dim=dim, total_cells=total_cells, perc_genes=perc_genes, perc_cells=perc_cells, flag_sim_100_cts=flag_sim_100_cts,ct_in_state = paste(ct_in_state,collapse='_'), adj.rsq = r_sq, max.rsq =max(corr) , max.rsq_sub = max(sub_corr_sq), total_cells_can_occupy_state = total_cells_can_occupy_state, total_cells_occupy_state = floor(total_cells_state * state_info[[state]][2]), total_genes = M, total_genes_in_state = total_genes_in_state)

    # append run to output dataframe
    write.table(df, file = "../output.txt", append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

    output = list(embeddings=embeddings, max_continuums = max_continuums,ct_labels=ct_labels)
    # saveRDS(output, glue("{method}_{state_type}_allcts_{perc_genes}_{perc_cells}_seed{seed}_{total_cells}_dim{dim}.rds"))    
} else {
    print("saving output on 100 cell types, synthetic")
    if (state_type == "within_one_ct"){
        df = data.frame(seed=seed, method=method, state_type=state_type, dim=dim, total_cells=total_cells, perc_genes=perc_genes, perc_cells=perc_cells, flag_sim_100_cts=flag_sim_100_cts,ct_in_state = "all", adj.rsq = r_sq, max.rsq =max(corr), max.rsq_sub = max(sub_corr_sq) , total_cells_can_occupy_state = total_cells_can_occupy_state, total_cells_occupy_state = floor(total_cells_state * state_info[[state]][2]), total_genes = M, total_genes_in_state = total_genes_in_state)

        print("100 cts within one ct")
            # append run to output dataframe
        write.table(df, file = "../output.txt", append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
        #output = list(embeddings=embeddings, max_continuums = max_continuums,ct_labels=ct_labels,state_info_ct=state_info_ct)
        #saveRDS(output, glue("{method}_{state_type}_allcts_{perc_genes}_{perc_cells}_seed{seed}_{total_cells}_{flag_sim_100_cts}_dim{dim}.rds"))
    } else { #"across_cts"
        print("100 cts across cts")
        df = data.frame(seed=seed, method=method, state_type=state_type, dim=dim, total_cells=total_cells, perc_genes=perc_genes, perc_cells=perc_cells, flag_sim_100_cts=flag_sim_100_cts,ct_in_state = "all", adj.rsq = r_sq, max.rsq =max(corr), max.rsq_sub = max(sub_corr_sq) , total_cells_can_occupy_state = total_cells_can_occupy_state, total_cells_occupy_state = floor(total_cells_state * state_info[[state]][2]), total_genes = M, total_genes_in_state = total_genes_in_state)

        print("100 cts within one ct")
            # append run to output dataframe
        write.table(df, file = "../output.txt", append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
        #output = list(embeddings=embeddings, max_continuums = max_continuums, ct_labels=ct_labels,state_info_ct=state_info_ct)
        #saveRDS(output, glue("{method}_{state_type}_allcts_{perc_genes}_{perc_cells}_seed{seed}_{total_cells}_{flag_sim_100_cts}_dim{dim}.rds"))
    }    
}
