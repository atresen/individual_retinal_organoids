#!/usr/bin/env Rscript
# R version 4.3.0
# Running under: macOS 13.3.1 with M1 Chip

suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra)
  library(viridis)
  library(devtools)
  library(monocle3)
  library(data.table)
  library(gridExtra)
  library(gghighlight)
  library(Matrix)
  library(RColorBrewer)
  source("monocle3_help_functions.R")
  source("chiSq_test_functions.R")
  
})

arg = commandArgs(trailingOnly=TRUE)

# test to ensure the file arg is not empty
if (length(args)==0) {
  stop("No input file supplied", call.=FALSE)
}

# Read in the unfiltered CDS (w.r.t. hashes)
# At this point count matrix has been made with all cells with greater than 100 UMIs
# if starting from GEO data...
cds = readRDS(paste0("data/R_objects/", arg[1], "_cds.RDS"))
head(colData(cds))



# Filter CDS based on a upper and lower treshold  -------------------------

sum(colData(cds)$n.umi > 100)

sum(colData(cds)$n.umi > 20000)

# filter for cells with more than 3 standard deviations above the mean
upper_cutoff = mean(colData(cds)$n.umi) + 3* sqrt(var(colData(cds)$n.umi))

ggplot(colData(cds) %>%
         as.data.frame()) +
  geom_histogram(aes(x = log10(n.umi)),
                 fill = NA,
                 color = "black") +
  geom_vline(xintercept = log10(upper_cutoff),
             color = "red") +
  theme_classic()

cds = cds[,colData(cds)$n.umi <= upper_cutoff]
dim(cds) #XXXXX cells in total


# Import and clean BBI hash output  ----------------------------------------------------------

# read cell names
cell_list = fread(paste0("data/", arg[1], ".hashumis_cells.txt"),
                  header = FALSE, 
                  data.table = F)[,1]
head(cell_list)
length(cell_list)

# read hash names
hash_list = fread(paste0("data/", arg[1], ".hashumis_hashes.txt"),
                  header = FALSE, 
                  data.table = F)[,1]
head(hash_list)
length(hash_list)

# read hash matrix and add cell and hash names (check that the lengths of each match)
hash_mtx = readMM(paste0("data/", arg[1], ".hashumis.mtx"))
hash_mtx[1:3, 1:3]
dim(hash_mtx)
rownames(hash_mtx) = hash_list
colnames(hash_mtx) = cell_list

# Flip it so that the hashes are columns and the 
hash_mtx = t(hash_mtx)
# Comes in as a sparse matrix
hash_mtx[1:3, 1:3]


# Read in UMIs.per.cell to get background cells
umis_per_cell = 
  read.table(paste0("data/umis_per_cell_barcode.txt"),
             col.names = c("Cell","n.umi"),
             sep = "\t")
#head(umis_per_cell)

# Apply Background Correction to hashTable --------------------------------
background_cell_hashes  =
  umis_per_cell$Cell[umis_per_cell$n.umi < 5] %>%
  as.character()

test_cells = 
  colnames(cds) %>%
  as.character()

corrected_hash_table = 
  assign_hash_labels(hash_matrix = hash_mtx,
                     test_cell_hashes = test_cells,
                     background_cell_hashes = background_cell_hashes)

corrected_hash_table

new_coldata =
  colData(cds) %>%
  as.data.frame() %>%
  left_join(corrected_hash_table,
            by = c("cell"="Cell"))

rownames(new_coldata) = new_coldata$cell
cds_with_hashes =
  new_cell_data_set(expression_data = counts(cds),
                    cell_metadata = new_coldata,
                    gene_metadata = 
                      rowData(cds) %>% 
                      as.data.frame(row.names = rownames(cds)))

fwrite(corrected_hash_table, paste0("data/", arg[1], "_corrected_hash_table.csv", sep = ","))

saveRDS(object = cds_with_hashes,
        file = paste0("data/R_objects/", arg[1], "_cds_with_hash_data.RDS"))
