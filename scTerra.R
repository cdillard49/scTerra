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

dds <- DESeqDataSetFromMatrix(countData = express_mat_51, colData = metadata_51, 
                              design = ~ condition)
#create the deseq object
dds$condition <- relevel(dds$condition, ref = "A")
#set condition
#here we set as A, but in the future will set condition appropriately
dds <- DESeq(dds)
#initialize the differential expression
res <- results(dds, independentFiltering = FALSE)
#initialize the results, with no filtering
res <- res[order(res$padj),]
#results by order of adjusted p-value
summary(results(dds, alpha=0.05))
#view summary of results that match alpha <0.05
normalized_counts <- counts(dds, normalized=TRUE)
#normalize counts variable
head(normalized_counts)
#view the normalized counts

resLFC<- lfcShrink(dds, coef="condition_B_vs_A", type = "normal", lfcThreshold = 1, alpha = 0.05)
plotMA(resLFC, ylim=c(-2,2), cex=.4)
plotMA(res, ylim=c(-2,2),cex=.4)
#visualize with MA plot
#logfold change shrinkage vs un-shrunk