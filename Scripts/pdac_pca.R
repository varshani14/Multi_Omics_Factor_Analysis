install.packages("ggplot2")
install.packages("ggfortify")
install.packages("dplyr")
library(ggplot2)
library(ggfortify)
library(dplyr)
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")
list.files("data")
dir.create("plots", showWarnings = FALSE)
tumor_purity <- read.csv("data/tumor_purity.csv", row.names = 1)
head(tumor_purity)
sample_ids <- rownames(tumor_purity)
table(substr(sample_ids, 1, 3))
pca_data <- tumor_purity
pca_data <- na.omit(pca_data)
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
sample_type <- ifelse(substr(rownames(pca_data), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:length(variance_explained), Variance = variance_explained)

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - tumor_purity",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_tumor_purity.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - tumor_purity",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("plots/scatter_tumor_purity.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
xcell <- read.csv("data/xcell.csv", row.names = 1)
head(xcell)
dim(xcell)
sample_type <- ifelse(substr(rownames(xcell), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- xcell
pca_data <- na.omit(pca_data)
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - xcell",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_xcell.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - xcell",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_xcell.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
cibersort <- read.csv("data/cibersort.csv", row.names = 1)
head(cibersort)
dim(cibersort)
sample_type <- ifelse(substr(rownames(cibersort), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- cibersort
pca_data <- na.omit(pca_data)
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - cibersort",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_cibersort.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - cibersort",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_cibersort.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
miRNA <- read.csv("data/miRNA.csv", row.names = 1)
head(miRNA)
dim(miRNA)
sample_type <- ifelse(substr(rownames(miRNA), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- miRNA
pca_data <- na.omit(pca_data)
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - miRNA",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_miRNA.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - miRNA",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_miRNA.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
circular_RNA <- read.csv("data/circular_RNA.csv", row.names = 1)
head(circular_RNA)
dim(circular_RNA)
sample_type <- ifelse(substr(rownames(circular_RNA), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- circular_RNA
pca_data <- na.omit(pca_data)
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
str(circular_RNA[, 1:5])
str(circular_RNA[, 1:15])
dim(circular_RNA)circular_RNA <- read.csv("data/circular_RNA.csv", row.names = 1, skip = 4)
dim(circular_RNA)
str(circular_RNA[, 1:5])
sample_type <- ifelse(substr(rownames(circular_RNA), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- circular_RNA
pca_data <- na.omit(pca_data)
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - circular_RNA",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_circular_RNA.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - circular_RNA",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_circular_RNA.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
sample_type <- ifelse(substr(rownames(pca_data), 1, 3) == "C3L", "Tumor", "Normal")
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - circular_RNA",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_circular_RNA.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
proteomics <- read.csv("data/proteomics.csv", row.names = 1)
head(proteomics[, 1:5])
dim(proteomics)
sample_type <- ifelse(substr(rownames(proteomics), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- proteomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
proteomics <- read.csv("data/proteomics.csv", row.names = 1, skip = 1)
dim(proteomics)
sample_type <- ifelse(substr(rownames(proteomics), 1, 3) == "C3L", "Tumor", "Normal")
pca_data <- proteomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - proteomics",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_proteomics.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - proteomics",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_proteomics.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
pca_scores_clean <- pca_scores[pca_scores$PC1 > -500, ]
scatter_plot <- ggplot(pca_scores_clean, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - proteomics",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_proteomics.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
CNV <- read.csv("data/CNV.csv", row.names = 1)
head(CNV[, 1:5])
dim(CNV)
CNV <- read.csv("data/CNV.csv", row.names = 1, skip = 1)
str(CNV[, 1:5])
sample_type <- ifelse(substr(rownames(CNV), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- CNV
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - CNV",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_CNV.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - CNV",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_CNV.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
library(data.table)
transcriptomics <- fread("data/transcriptomics.csv")
transcriptomics <- as.data.frame(transcriptomics)
rownames(transcriptomics) <- transcriptomics[, 1]
transcriptomics <- transcriptomics[, -1]
sample_type <- ifelse(substr(rownames(transcriptomics), 1, 3) == "C3L", "Tumor", "Normal")
pca_data <- transcriptomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
transcriptomics <- fread("data/transcriptomics.csv", skip = 1)
str(transcriptomics[, 1:5])
transcriptomics <- as.data.frame(transcriptomics)
rownames(transcriptomics) <- transcriptomics$Database_ID
transcriptomics$Database_ID <- NULL
dim(transcriptomics)
sample_type <- ifelse(substr(rownames(transcriptomics), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- transcriptomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
nrow(pca_result$x)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - transcriptomics",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_transcriptomics.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - transcriptomics",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_transcriptomics.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
pca_scores_clean <- pca_scores[pca_scores$PC1 > -500, ]
scatter_plot <- ggplot(pca_scores_clean, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - transcriptomics",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_transcriptomics.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(data.table)
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")
phospho <- fread("data/phosphoproteomics.csv", skip = 1)
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]
dim(phospho)
str(phospho[, 1:5])
phospho <- fread("data/phosphoproteomics.csv", skip = 3)
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]
str(phospho[, 1:5])
sample_type <- ifelse(substr(rownames(phospho), 1, 3) == "C3L", "Tumor", "Normal")
table(sample_type)
pca_data <- phospho
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
dim(pca_data)
pca_data <- phospho
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.7]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
dim(pca_data)
str(pca_data[, 1:5])
pca_data <- pca_data[, sapply(pca_data, function(x) sd(x) > 0)]
sum(is.na(colnames(pca_data)))
pca_data_matrix <- as.matrix(pca_data)
keep_cols <- which(apply(pca_data_matrix, 2, sd) > 0)
pca_data_matrix <- pca_data_matrix[, keep_cols]
dim(pca_data_matrix)
pca_result <- prcomp(pca_data_matrix, center = TRUE, scale. = TRUE)
nrow(pca_result$x)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = 1:10, Variance = variance_explained[1:10])

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "#E67E22") +
  geom_line(color = "#16A085", size = 1) +
  geom_point(color = "#16A085", size = 3) +
  labs(title = "Scree Plot - phosphoproteomics",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("plots/scree_phosphoproteomics.png", scree_plot, width = 7, height = 5, dpi = 300)
print(scree_plot)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

scatter_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - phosphoproteomics",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_phosphoproteomics.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
summary(pca_scores$PC1)
pca_scores_clean <- pca_scores[pca_scores$PC1 > -100, ]
scatter_plot <- ggplot(pca_scores_clean, aes(x = PC1, y = PC2, color = SampleType)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "PCA Scatter Plot - phosphoproteomics",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Sample Type") +
  theme_minimal()

ggsave("plots/scatter_phosphoproteomics.png", scatter_plot, width = 7, height = 5, dpi = 300)
print(scatter_plot)
clinical <- read.csv("data/clinical.csv", row.names = 1)
head(rownames(clinical))
clinical <- read.csv("data_csv/clinical.csv", row.names = 1)
head(rownames(clinical))
pca_scores_stage <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- substr(rownames(pca_scores_stage), 1, 9)
pca_scores_stage$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]
library(data.table)

transcriptomics <- fread("data/transcriptomics.csv", skip = 1)
transcriptomics <- as.data.frame(transcriptomics)
rownames(transcriptomics) <- transcriptomics$Database_ID
transcriptomics$Database_ID <- NULL

sample_type <- ifelse(substr(rownames(transcriptomics), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- transcriptomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type
transcriptomics <- as.data.frame(transcriptomics)
rownames(transcriptomics) <- transcriptomics$Database_ID
transcriptomics$Database_ID <- NULL

sample_type <- ifelse(substr(rownames(transcriptomics), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- transcriptomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type
dim(pca_data)
install.packages("irlba")
library(irlba)
pca_data_matrix <- as.matrix(pca_data)
pca_data_scaled <- scale(pca_data_matrix, center = TRUE, scale = TRUE)

pca_result_irlba <- irlba::prcomp_irlba(pca_data_scaled, n = 10, center = FALSE, scale. = FALSE)

variance_explained <- (pca_result_irlba$sdev^2) / sum(apply(pca_data_scaled, 2, var)) * 100

pca_scores <- as.data.frame(pca_result_irlba$x)
colnames(pca_scores) <- paste0("PC", 1:10)
pca_scores$SampleType <- sample_type
head(pca_scores)
sum(is.na(pca_data_matrix))
rowSums(is.na(pca_data_matrix))
pca_data_matrix <- pca_data_matrix[-1, ]
sample_type <- sample_type[-1]
pca_data_scaled <- scale(pca_data_matrix, center = TRUE, scale = TRUE)
pca_result_irlba <- irlba::prcomp_irlba(pca_data_scaled, n = 10, center = FALSE, scale. = FALSE)

variance_explained <- (pca_result_irlba$sdev^2) / sum(apply(pca_data_scaled, 2, var)) * 100

pca_scores <- as.data.frame(pca_result_irlba$x)
colnames(pca_scores) <- paste0("PC", 1:10)
pca_scores$SampleType <- sample_type
head(pca_scores)
