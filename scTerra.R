require(data.table)

setwd("~/GitHub/scTerra")

express_mat_51 <- read.delim("inputs/SRR10018151_GE_matrix_filtered.txt", row.names = 1)
#1802 cell barcodes
redi_cos_51 <- read.delim("inputs/SRR10018151_ind_ss_redi_cos.txt")
#17210 snvs each with corresponding cell barcodes
#many SNVs found in different cells
cells_51 <- as.data.frame(redi_cos_51$Cells)
#only barcodes (comma seperated)
barcodelist <- colnames(express_mat_51) 
#list of barcodes