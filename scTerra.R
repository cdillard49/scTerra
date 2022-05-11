require(data.table)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(tidyr)
library(dplyr)
library(DESeq2)
setwd("~/GitHub/scTerra")

end_reads_50 <- read.delim("inputs/end_reads_barcodes50")
express_mat_51 <- read.delim("inputs/SRR10018151_GE_matrix_filtered.txt", row.names = 1)
#1802 cell barcodes
redi_cos_51 <- read.delim("inputs/SRR10018151_ind_ss_redi_cos.txt")
#17210 snvs each with corresponding cell barcodes
#many SNVs found in different cells
barcodelist <- colnames(express_mat_51) 
#list of barcodes to use for filtering
redi_separate_51 <- separate_rows(redi_cos_51, Cells, convert = TRUE)
#separate out/expand the Cells column
#note this breaks the N_cells column for now
filt_redit_51<- redi_separate_51 %>% filter(Cells %in% barcodelist)
#keep only cells found in the GE barcode list



#DESEQ2 GE using random cells (WIP: will use cells with specific scTerra SNVs)
letters <- LETTERS[1:2]
fac <- sample(letters, 1802, replace = TRUE)
metadata_51 <- data.frame(
  sample = barcodelist,
  condition = fac,
  row.names = "sample")
metadata_51$condition <- as.factor(metadata_51$condition)
#metadata creation, change the condition column to match SNV condition

pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)