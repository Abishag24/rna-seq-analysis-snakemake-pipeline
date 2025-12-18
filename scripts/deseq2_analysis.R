###############################################################################
# Script: deseq2_analysis.R
# Purpose: Differential gene expression analysis using DESeq2
# Organism: Saccharomyces cerevisiae
#
# Inputs:
#   - counts/gene_counts.txt       (featureCounts output)
#   - counts/Column_Data_new.csv   (sample metadata)
#
# Outputs:
#   - results/All_DEGs.csv
#   - results/Upregulated_genes.csv
#   - results/Downregulated_genes.csv
#   - QC and visualization plots (PCA, MA, Volcano, Heatmaps)
#   - GO enrichment results
#
# Execution:
#   This script is executed via Snakemake.
#   Not intended to be run manually.
#    
#  Author: Abishag Jacquline
###############################################################################

#load the library 
library(DESeq2)
library(matrixStats)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(enrichplot)

counts <- read.table("counts/gene_counts.txt", 
                     header=TRUE, sep="\t", comment.char = "#", check.names = FALSE)

# Paths are fixed for standalone reproducibility


#remove annotation columns and make count_matrix ready for DESeq2
count_matrix <- counts[, 7:ncol(counts)]
rownames(count_matrix) <- counts$Geneid
Count_data <- as.matrix(count_matrix)

#load column data (what each column represents)
Col_data <- read.csv(file = "counts/Column_Data_new.csv", header = T, sep = ",",row.names = 1)

# Clean sample names (remove path and BAM suffix)

colnames(Count_data) <- sub("^aligned/", "", colnames(Count_data))
colnames(Count_data) <- sub("\\.Aligned.sortedByCoord.out\\.bam$", "", colnames(Count_data))

# Ensure sample names match between count matrix and metadata

if (all(rownames(Col_data) == colnames(Count_data))) {
  message("✅ Sample names match between Count_data and Col_data")
} else {
  stop("❌ Sample names do NOT match between Count_data and Col_data")
}

#count number of NA values in matrix
sum(is.na(Count_data)) #how many NAs are there


#remove or replace NA values, if there (optional)
# Filtering applied at raw count level (zero counts) and DESeq2 level (low-expression genes)
#replace NA values in matrix with zero
Count_data[is.na(Count_data)] <- 0


# Filter genes with zero counts across all samples
Count_data <- Count_data[rowSums(Count_data) > 0, ]

##################################################################################################

#creating DESeq object
dds <- DESeqDataSetFromMatrix(countData = Count_data,
                              colData = Col_data,
                              design = ~condition) # we're testing for the different conditions




dds


#set the factor level (normal first n then condition)
dds$condition <- relevel(dds$condition, ref = "WT") #check from the col_data files, write the name based on that


# additional filtering (optional)
filtered_count <- rowSums(counts(dds)) >=10 #filter those rows with at least 10 reads in total
dds <- dds[filtered_count,]



# Run DESeq
dds <- DESeq(dds)


# Remove rows with NA padj
res <- results(dds)
res_clean <- res[complete.cases(res$padj, res$log2FoldChange), ]
summary(res_clean)
#contrast argument
resultsNames(dds)

#for multiple group comparison (when you have more than 2 condition) 
#results(dds, contrast = ("normal",  "condition1", "condition2"))

###keeping only significant results, padj<0.05 and log2FoldChange >1
resSigUp <- subset(res_clean, padj < 0.05 & log2FoldChange >1) #Cutoffs can be adjusted depending on study design and stringency requirements


###keep only sig results, padj<0.05 and log2FoldChange < -1
resSigDown <- subset(res_clean, padj < 0.05 & log2FoldChange < -1) 


write.csv(resSigUp, "results/Upregulated_genes.csv")
write.csv(resSigDown, "results/Downregulated_genes.csv")
### Combined results in a single file Upregulated genes and Downregulated genes
resSig_combined <- subset(
  res_clean,  
  padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)
)

write.csv(resSig_combined, "results/All_DEGs.csv")
##################################################
#Visualization

#1. Volcano plot

png("results/volcano.png", width = 1200, height = 1000, res = 150)


EnhancedVolcano(res_clean,
                lab = rownames(res_clean),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c( 'YBL006C','YIR005W','YNR057C', 'YIR035C','YHR092C',       #up
                               'YHL046W-A','YHL047C','YOR382W','YKR035W-A','YHL040C'),    #down
                title = NULL,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,15),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                labSize = 3,
                border = 'full' )



dev.off()


#2. MA PLOT

png("results/MA_plot.png",
    width = 1200, height = 1000, res = 150)
plotMA(res_clean)
dev.off()


# normalized counts
norm_counts <- counts(dds, normalized = TRUE) #normalized data
head(norm_counts)



##### Heatmap (normalized data)#######

# Upregulated genes
up_genes <- rownames(resSigUp)


# Downregulated genes
down_genes <- rownames(resSigDown)

# Variance-stabilizing transformation (for visualization)
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)


# Significant genes (padj < 0.05, |log2FC| > 1)
sig_genes_heatmap <- rownames(
  res_clean[res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1, ]
)
sig_genes_enrich <- rownames(resSig_combined)

# Extract vst values for significant genes
mat_sig <- mat[sig_genes_heatmap, ]

length(sig_genes_heatmap)
length(sig_genes_enrich)

# ---- Heatmap of all significant DEGs ----

set.seed(123)
pheatmap(mat_sig,
         scale = "row",                      # row-wise normalization
         show_rownames = FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = Col_data)


# ---- Heatmap of top 100 most variable DEGs ----

set.seed(123)
if (nrow(mat_sig) > 5) {
  top_genes <- head(order(rowVars(mat_sig), decreasing=TRUE), 50)
  pheatmap(mat_sig[top_genes, ],
           scale = "row",
           show_rownames = TRUE,
           annotation_col = Col_data)
  
}

#################################################################################

#PCA plot (VST-normalized DESeq2 data)

# PCA on variance-stabilized data
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))


png("results/PCA.png",
    width = 1200, height = 1000, res = 150)

set.seed(123)
ggplot(pca_data, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA plot")

dev.off()


##################################################################################

# read csv file from g-profiler 

gene_df <- read.csv("gProfiler_scerevisiae_12-9-2025.csv")
# keep as character vector
gene_names <- as.character(gene_df$name)

#clean
gene_names <- gene_names[!is.na(gene_names)]
gene_names <- gene_names[gene_names != ""]

## gProfiler results loaded for external validation only
# Not used directly in clusterProfiler enrichment


# Take significant genes from DESeq2

# Map ORFs to ENTREZ IDs
entrez_ids <- mapIds(org.Sc.sgd.db, 
                     keys = sig_genes_enrich, 
                     column = "ENTREZID", 
                     keytype = "ORF", 
                     multiVals = "first")

# Remove NAs

entrez_ids <- na.omit(entrez_ids)



# GO ENRICHMENT

# for Biological Process
go_bp <- enrichGO(gene        = entrez_ids,
                  OrgDb         = org.Sc.sgd.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = FALSE)
#Because readable=TRUE tries to convert ENTREZ → SYMBOL, but 
#yeast annotation package does not have a SYMBOL column, so the conversion fails


# Create the plot object , display and save
if(nrow(as.data.frame(go_bp)) > 0){ 
  p <-dotplot(go_bp, showCategory = 20) +
    ggtitle("GO BP Enrichment") +
    theme(axis.text.y = element_text(size = 10))
  
  print(p)
  
  ggsave("results/GO_BP_dotplot.png",
         plot = p, width = 8, height = 7, dpi = 150)
}
#-----
# for Molecular Function

# Now run GO MF enrichment
go_mf <- enrichGO(gene        = entrez_ids,
                  OrgDb       = org.Sc.sgd.db,
                  keyType     = "ENTREZID",
                  ont         = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.2,
                  readable    = FALSE)

# create the plot object , print nd save

if(nrow(as.data.frame(go_mf)) > 0){ 
  p <- dotplot(go_mf, showCategory = 20) +
    ggtitle("GO MF Enrichment") +
    theme(axis.text.y = element_text(size = 10))
  
  print(p)
  
  ggsave("results/GO_MF_dotplot.png",
         plot = p, width = 8, height = 7, dpi = 150)
}
#-------------



# for Cellular Component


go_cc <- enrichGO(gene        = entrez_ids,
                  OrgDb         = org.Sc.sgd.db,
                  keyType       = "ENTREZID",
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = FALSE)

# Cellular Component (CC) terms in yeast are:
# Broad (e.g. nucleus, cytoplasm, ribosome) ,Highly shared across many genes
# and are spread across multiple locations
# Because of this: Enrichment statistics cannot find a strong signal

# NOTE:
# CC terms showed no significant enrichment after multiple testing correction
# (adjusted p-values > 0.9). Cutoffs were relaxed to visualize CC distribution
# without claiming statistical significance.



# to save the plot

if(nrow(as.data.frame(go_cc)) > 0) {
  p <-dotplot(go_cc, showCategory = 20) +
    ggtitle("GO CC Enrichment") +
    labs(subtitle = "No terms significant after BH correction")+
    theme(axis.text.y = element_text(size = 10))
  
  print(p)
  
  ggsave("results/GO_CC_dotplot.png",
         plot = p, width = 8, height = 7, dpi = 150)
}
#_---------------------

# KEGG

# For Saccharomyces cerevisiae, KEGG uses systematic ORF IDs
# (e.g., YLR070C) rather than Entrez IDs.

#  “For yeast, KEGG pathways are indexed using systematic ORF IDs from SGD. Unlike human data, Entrez IDs are 
#not the primary identifier for S. cerevisiae, so ORF IDs are the correct input.”

# KEGG enrichment was attempted but skipped
# Reason: KEGG REST API does not reliably support ORF-based mapping for yeast
# GO enrichment (BP, MF, CC) used instead




# Export enrichment results

write.csv(as.data.frame(go_bp), "results/GO_BP.csv")
write.csv(as.data.frame(go_mf), "results/GO_MF.csv")
write.csv(as.data.frame(go_cc), "results/GO_CC.csv")



# Save session info for reproducibility
writeLines(capture.output(sessionInfo()),
           "results/sessionInfo.txt")








