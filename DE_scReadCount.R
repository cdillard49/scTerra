if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
devtools::install_github('kevinblighe/EnhancedVolcano')

require(data.table)
library(tidyr)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
setwd("~/GitHub/scTerra")
sample_name <- "SRR10018151"

#use to generate variable for naming later, only use sample name that matches
# the one from line 19
snv_input <- read.delim("inputs/example.txt")
#file needs to have column with "SNV" and row with SNV location/coordinates
snv_list <- snv_input$SNV 
express_mat_51 <- read.delim("inputs/SRR10018151_GE_matrix_filtered.txt", row.names = 1)
screadcounts <- na.omit(read.delim("inputs/SRR10018151.tsv"))
barcodelist <- colnames(express_mat_51) 
#list of barcodes to use for filtering

##redi_cos_51 <- subset(redi_cos_51, subset=(N_cells >= 10))
##remove low cell count (10)


#separate out/expand the Cells column
#note this breaks the N_cells column for now
filt_screadcounts <- screadcounts %>% filter(ReadGroup %in% barcodelist)
vaf_cells <- filt_screadcounts[filt_screadcounts$VAF == 0, ]
rm(screadcounts)
#keep only cells found in the GE barcode list

####


#DESEQ2 GE using cells with specific SNVs
for (input in snv_list) {
  for (string in input) {
    match <- as.data.frame(str_split(input, ":|_|>", simplify = TRUE))
  }
  snv <- filt_screadcounts[filt_screadcounts$CHROM %in% match$V1, ]
  snv <- snv[snv$POS %in% match$V2, ]
  snv <- snv[snv$REF %in% match$V3, ]
  snv <- snv[snv$ALT %in% match$V4, ]
  snv_vaf <- snv[snv$VAF > 0, ]

  snv_cells_A <- snv_vaf$ReadGroup
  #A is chosen SNV
  snv_cells_B <- vaf_cells$ReadGroup
  #B is all other SNVs
  #filtering out only those snvs selected and setting the condition for deseq as A
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
  output_gene <- as.data.frame(res[order(res$log2FoldChange, decreasing=TRUE),])
  #output genes ordered by upregulated genes, set decreasing=FALSE for downregulated
  
  file_name_genes <- paste0("outputs/", sample_name, "_DE_genes_", input, ".tsv")
  file_name_genes <- gsub(":|>", "_", file_name_genes)
  
  file_name_plotMA_shrink <- paste0("outputs/", sample_name, "_DE_MA_Shrink_", input, ".pdf")
  file_name_plotMA_shrink <- gsub(":|>", "_",file_name_plotMA_shrink)
  
  file_name_plotMA <- paste0("outputs/", sample_name, "_DE_MA_", input, ".pdf")
  file_name_plotMA <- gsub(":|>", "_", file_name_plotMA)
  
  file_name_volcano_shrink <- paste0("outputs/", sample_name, "_DE_Volcano_", input, ".pdf")
  file_name_volcano_shrink <- gsub(":|>", "_", file_name_volcano_shrink)
  #file name format for all plots and gene table, removing ":" and ">"
  
  write.table(x = output_gene, file = file_name_genes, quote = FALSE, sep = "\t", col.names = NA)
  #write gene table to file
  
  resLFC<- lfcShrink(dds, coef="condition_A_vs_B", type = "normal", lfcThreshold = 1)
  #apply lfcShrink to results
  
  pdf(file = file_name_volcano_shrink)
  p <- EnhancedVolcano(resLFC, lab=rownames(resLFC), x="log2FoldChange", y= "pvalue", 
                       pCutoff = 1, title = "Normal vs. SNV" ,xlab = bquote(~Log[2]~ "fold change"), 
                       colAlpha = .5, pointSize = 1.0, labSize = 3.0, col = c("black", "blue", "green", "red"), 
                       legendPosition = "right",
                       legendLabSize = 10.0,
                       legendIconSize = 3.0)
  print(p)
  dev.off()
  
  pdf(file = file_name_plotMA_shrink)
  plotMA(resLFC, ylim=c(-2,2), cex=.4)
  dev.off()
  
  pdf(file = file_name_plotMA)
  plotMA(res, ylim=c(-2,2),cex=.4)
  dev.off()
  
  
  print(paste0("SNV ", input, " finished."))
  #visualize with MA plot
  #logfold change shrinkage vs un-shrunk
  #volcanoPlots
  #write plots to pdf format
  
}

####
print("DONE")
