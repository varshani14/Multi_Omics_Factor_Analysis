
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
install.packages(c(
  "pheatmap",
  "ggplot2",
  "plotly",
  "corrplot"
))

# Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "EDASeq"
))

library(pheatmap) # for drawing heatmaps of gene expression and clustering samples/genes
library(ggplot2) # for creating publication-quality plots and visualizations
library(plotly) # for interactive plots such as interactive PCA or 3D visualizations
library(corrplot) # for visualizing correlation matrices between samples
library(DESeq2) # for differential gene expression analysis of RNA-seq count data
library(EDASeq) # for exploratory data analysis and normalization of sequencing data

# load the saved RNA-seq dataset (counts and sample metadata) into the R environment
load("colorectal_rnaseq.RData")

# basic commands to explore the dataset and inspect objects in the R environment
str(counts)
head(counts)
tail(counts)
colnames(counts)
row.names(counts)
dim(counts)
summary(counts[,1:3])
counts[1:5, 1:5]
ls()
rm(list=ls())

##### ------------ Exploring the normalization approaches ####------------------

### To compute the CPM values for each sample (excluding the width column)
cpm <- apply(subset(counts, select = c(-width)), 2, 
             function(x) x/sum(as.numeric(x)) * 10^6)

colSums(cpm)

# create a vector of gene lengths 
geneLengths <- as.vector(subset(counts, select = c(width)))

### Computing RPKM
rpkm <- apply(X = subset(counts, select = c(-width)),
              MARGIN = 2, 
              FUN = function(x) {
                10^9 * x / geneLengths / sum(as.numeric(x))
              })

colSums(rpkm)

#find gene length normalized values 
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000))

#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

colSums(tpm)

##### ------------ Exploratory analysis of the read count table ####------------

### Clustering

#compute the variance of each gene across samples
V <- apply(cpm, 1, var)
#sort the results by variance in decreasing order 
#and select the top 100 genes 
selectedGenes <- names(V[order(V, decreasing = T)][1:50])

pheatmap(tpm[selectedGenes,], scale = 'row', show_rownames = FALSE)

pheatmap(tpm[selectedGenes,], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = coldata)

### PCA

#transpose the matrix 
M <- t(cpm[selectedGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1)

#compute PCA 

pcaResults <- prcomp(M, scale. = TRUE)

pca_df <- as.data.frame(pcaResults$x)

# add sample group information
pca_df$group <- coldata$group
pca_df$sample <- rownames(pca_df)

# 3D PCA plot
plot_ly(pca_df,
        x = ~PC1,
        y = ~PC2,
        z = ~PC3,
        color = ~group,
        colors = c("CTRL" = "blue", "CASE" = "red"),
        text = ~sample,
        type = "scatter3d",
        mode = "markers")

plot_ly(pca_df,
        x = ~PC1,
        y = ~PC2,
        color = ~group,
        colors = c("CTRL" = "blue", "CASE" = "red"),
        text = ~sample,
        type = "scatter",
        mode = "markers")

var_exp <- (pcaResults$sdev^2 / sum(pcaResults$sdev^2)) * 100

plot_ly(pca_df,
        x = ~PC1,
        y = ~PC2,
        color = ~group,
        colors = c("CTRL" = "blue", "CASE" = "red"),
        text = ~sample,
        type = "scatter",
        mode = "markers") %>%
  layout(xaxis = list(title = paste0("PC1 (", round(var_exp[1],1), "%)")),
         yaxis = list(title = paste0("PC2 (", round(var_exp[2],1), "%)")))

summary(pcaResults)

### Correlation plots
correlationMatrix <- cor(tpm)
corrplot(correlationMatrix, order = 'hclust', 
         addrect = 2, addCoef.col = 'white', 
         number.cex = 0.5) 

# split the clusters into two based on the clustering similarity 
pheatmap(correlationMatrix,  
         annotation_col = coldata, 
         cutree_cols = 2)

##### ---------------- Differential expression analysis ####--------------------

countData <- as.matrix(subset(counts, select = c(-width)))
designFormula <- "~ group"
coldata$group <- factor(coldata$group)

#create a DESeq dataset object from the count matrix and the colData 
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = coldata,
  design = ~ group
)
#print dds object to see the contents
print(dds)

#For each gene, we count the total number of reads for that gene in all samples 
#and remove those that don't have at least 1 read. 
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

dds <- DESeq(dds)

#compute the contrast for the 'group' variable where 'CTRL' 
#samples are used as the control group. 
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]

#shows a summary of the results
print(DEresults)
DEresults_df <- as.data.frame(DEresults)
#library(writexl)
#write_xlsx(DEresults_df, "DESeq2_results.xlsx")

##### ------------------------ Diagnostic plots ####----------------------------

#MA plot
DESeq2::plotMA(object = dds, ylim = c(-5, 5))

#Volcano plot
res <- results(dds)
res_df <- as.data.frame(res)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.5) +
  xlab("log2 Fold Change") +
  ylab("-log10(p-value)") +
  theme_minimal()

res_df$significant <- res_df$padj < 0.05

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  xlab("log2 Fold Change") +
  ylab("-log10(p-value)") +
  theme_minimal()


res_df$regulation <- "NS"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Down" = "blue",
                                "NS" = "grey",
                                "Up" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(x = "log2 Fold Change",
       y = "-log10(p-value)",
       color = "Regulation")

# P-value distribution
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + 
  geom_histogram(bins = 100)


# extract normalized counts
countsNormalized <- DESeq2::counts(dds, normalized = TRUE)

# select top 500 most variable genes
selectedGenes <- names(sort(apply(countsNormalized, 1, var),
                            decreasing = TRUE)[1:500])

# transpose so samples are rows
M <- t(countsNormalized[selectedGenes,])

# log transform
M <- log2(M + 1)

# compute PCA
pcaResults <- prcomp(M, scale. = TRUE)

# convert to dataframe
pca_df <- as.data.frame(pcaResults$x)

# add metadata
pca_df$group <- coldata$group
pca_df$sample <- rownames(pca_df)

plot_ly(pca_df,
        x = ~PC1,
        y = ~PC2,
        color = ~group,
        text = ~sample,
        type = "scatter",
        mode = "markers")

plot_ly(pca_df,
        x = ~PC1,
        y = ~PC2,
        z = ~PC3,
        color = ~group,
        text = ~sample,
        type = "scatter3d",
        mode = "markers")

#Relative Log Expression (RLE) plot
par(mfrow = c(1, 2))

plotRLE(countData,
        outline = FALSE,
        ylim = c(-4, 4),
        col = as.numeric(coldata$group),
        main = "Raw Counts",
        cex.axis = 0.7,
        las = 2)

plotRLE(DESeq2::counts(dds, normalized = TRUE),
        outline = FALSE,
        ylim = c(-4, 4),
        col = as.numeric(coldata$group),
        main = "Normalized Counts",
        cex.axis = 0.7,
        las = 2)
