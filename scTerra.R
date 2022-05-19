require(data.table)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(tidyr)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
setwd("~/GitHub/scTerra")
#end_reads_50 <- read.delim("inputs/end_reads_barcodes50")

snv_input <- read.delim("example.txt/tsv") ##Not Ran
#file needs to have column with SNV and row with SNV location/coordinates
snv_list <- snv_input$SNV ##Not Ran

express_mat_51 <- read.delim("inputs/SRR10018151_GE_matrix_filtered.txt", row.names = 1)
#1802 cell barcodes
redi_cos_51 <- read.delim("inputs/SRR10018151_ind_ss_redi_cos.txt")
#17210 snvs each with corresponding cell barcodes
#many snvss found in different cells
barcodelist <- colnames(express_mat_51) 
#list of barcodes to use for filtering
----
#redi_cos_51 <- subset(redi_cos_51, subset=(N_cells >= 10))
#remove low cell count (10)
----
redi_separate_51 <- separate_rows(redi_cos_51, Cells, convert = TRUE)
#separate out/expand the Cells column
#note this breaks the N_cells column for now
filt_sub_redi_51 <- redi_separate_51 %>% filter(Cells %in% barcodelist)
#keep only cells found in the GE barcode list


#DESEQ2 GE using cells with specific SNVs
#line 33: change the SNV to the desired SNV of choice
snv <- filt_sub_redi_51[filt_sub_redi_51$SNV %in% c("1:100291620_T>G"), ]
snv_cells_A <- snv$Cells
#A is chosen SNV
snv_cells_B <- barcodelist[!barcodelist %in% snv_cells_A]
#B is all other SNVs
#filtering out only those snvs selected and setting the condition for dseq as A
snv_cells_A <- as.data.frame(snv_cells_A)
snv_cells_B <- as.data.frame(snv_cells_B)
snv_cells_A$condition <- "A"
snv_cells_B$condition <- "B"
names(snv_cells_A)[1] <- "barcode"
names(snv_cells_B)[1] <- "barcode"
dseq_meta <- rbind(snv_cells_A, snv_cells_B)
#binding the two conditional dataframes together 
dseq_meta$condition <- as.factor(dseq_meta$condition)

coldata <- data.frame(barcode = barcodelist, row.names = "barcode")
dseq_meta_index <- match(rownames(coldata), dseq_meta$barcode)
#create the GE matrix barcode dataframe for reordering 
#match the GE barcode to the metadata barcode and get an index for reordering
dseq_meta_reorder <- as.data.frame(dseq_meta[dseq_meta_index, ])
dseq_meta_reorder <- data.frame(dseq_meta_reorder, row.names = "barcode")
#reordering and dataframe manipulation to get into deseq2 format

#metadata creation

dds <- DESeqDataSetFromMatrix(countData = express_mat_51, colData = dseq_meta_reorder, 
                              design = ~ condition)
#create the deseq object
dds$condition <- relevel(dds$condition, ref = "B")
#set reference condition
#here we set B reference/control
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

resLFC<- lfcShrink(dds, coef="condition_A_vs_B", type = "normal", lfcThreshold = 1)
plotMA(resLFC, ylim=c(-2,2), cex=.4)
plotMA(res, ylim=c(-2,2),cex=.4)
EnhancedVolcano(resLFC, lab=rownames(resLFC), x="log2FoldChange", y= "pvalue", 
                title = "Normal vs. SNV" ,xlab = bquote(~Log[2]~ "fold change"), 
                colAlpha = .5, pointSize = 1.0, labSize = 3.0, col = c("black", "blue", "green", "red"), 
                legendPosition = "right" )

#visualize with MA plot
#logfold change shrinkage vs un-shrunk
#volcanoPlots