require(data.table)
library(tidyr)
library(dplyr)
setwd("~/GitHub/scTerra")

express_mat_51 <- read.delim("inputs/SRR10018151_GE_matrix_filtered.txt", row.names = 1)
#1802 cell barcodes
redi_cos_51 <- read.delim("inputs/SRR10018151_ind_ss_redi_cos.txt")
#17210 snvs each with corresponding cell barcodes
#many SNVs found in different cells
barcodelist <- colnames(express_mat_51) 
#list of barcodes to use for filtering

redi_separate_51 <- separate_rows(redi_cos_51, Cells, convert = TRUE)
#separate out/expand the Cells column
#note this breaks the N_cells column now
filt_redit_51<- redi_separate_51 %>% filter(Cells %in% barcodelist)
#keep only cells found in the GE barcode list