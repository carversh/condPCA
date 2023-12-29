library(glue)
args = commandArgs(trailingOnly=TRUE)

cells = args[1]
genes = args[2]
print(cells)
print(genes)
counts = readRDS(glue("/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_12345.rds") )
write.table(counts, file =glue( "/data/gusev/USERS/scarver/create_datasets/data/cells{cells}_genes_{genes}_seed_12345.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

