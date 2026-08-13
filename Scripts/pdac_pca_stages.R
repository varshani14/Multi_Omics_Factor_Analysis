library(ggplot2)
library(dplyr)
library(data.table)

setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")

clinical <- read.csv("data_csv/clinical.csv", row.names = 1)
dir.create("plots_pca_stages", showWarnings = FALSE)

library(data.table)

transcriptomics <- fread("data_csv/transcriptomics.csv", skip = 1)
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
head(pca_scores[, 1:3])
tail(rownames(transcriptomics), 5)
sum(rownames(transcriptomics) == "Patient_ID")
transcriptomics <- transcriptomics[rownames(transcriptomics) != "Patient_ID", ]

sample_type <- ifelse(substr(rownames(transcriptomics), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- transcriptomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- rownames(pca_scores_tumor)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

table(pca_scores_stage$Stage)
stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)
table(pca_scores_stage$Stage_Label)
stage_plot <- ggplot(pca_scores_stage, aes(x = PC1, y = PC2, color = Stage_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Scatter Plot - Transcriptomics - Coloured by Pathological Stage",
       subtitle = paste0("CPTAC PDAC Patients (n = ", nrow(pca_scores_stage), ")"),
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal()

ggsave("plots_pca_stages/pca_transcriptomics_stage.png", stage_plot, width = 8, height = 5, dpi = 300)
print(stage_plot)
nrow(pca_scores_stage[pca_scores_stage$Stage == "Stage I", ])
library(data.table)

CNV <- fread("data_csv/CNV.csv", skip = 1)
CNV <- as.data.frame(CNV)
rownames(CNV) <- CNV[, 1]
CNV <- CNV[, -1]

sample_type <- ifelse(substr(rownames(CNV), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- CNV
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")
getwd()
library(data.table)

CNV <- fread("data_csv/CNV.csv", skip = 1)
CNV <- as.data.frame(CNV)
rownames(CNV) <- CNV[, 1]
CNV <- CNV[, -1]

sample_type <- ifelse(substr(rownames(CNV), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- CNV
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
CNV <- CNV[rownames(CNV) != "Patient_ID", ]

sample_type <- ifelse(substr(rownames(CNV), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- CNV
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- rownames(pca_scores_tumor)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)

table(pca_scores_stage$Stage_Label)
clinical <- read.csv("data_csv/clinical.csv", row.names = 1)
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- rownames(pca_scores_tumor)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)

table(pca_scores_stage$Stage_Label)
stage_plot <- ggplot(pca_scores_stage, aes(x = PC1, y = PC2, color = Stage_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Scatter Plot - CNV - Coloured by Pathological Stage",
       subtitle = paste0("CPTAC PDAC Patients (n = ", nrow(pca_scores_stage), ")"),
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal()

ggsave("plots_pca_stages/pca_CNV_stage.png", stage_plot, width = 8, height = 5, dpi = 300)
print(stage_plot)
library(ggplot2)
stage_plot <- ggplot(pca_scores_stage, aes(x = PC1, y = PC2, color = Stage_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Scatter Plot - CNV - Coloured by Pathological Stage",
       subtitle = paste0("CPTAC PDAC Patients (n = ", nrow(pca_scores_stage), ")"),
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal()

ggsave("plots_pca_stages/pca_CNV_stage.png", stage_plot, width = 8, height = 5, dpi = 300)
print(stage_plot)
proteomics <- read.csv("data_csv/proteomics.csv", row.names = 1, skip = 1)

sample_type <- ifelse(substr(rownames(proteomics), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- proteomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
proteomics <- proteomics[rownames(proteomics) != "Patient_ID", ]

sample_type <- ifelse(substr(rownames(proteomics), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- proteomics
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.3]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- rownames(pca_scores_tumor)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)

table(pca_scores_stage$Stage_Label)
head(rownames(pca_scores_stage))
nrow(proteomics)
library(ggplot2)

stage_plot <- ggplot(pca_scores_stage, aes(x = PC1, y = PC2, color = Stage_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Scatter Plot - Proteomics - Coloured by Pathological Stage",
       subtitle = paste0("CPTAC PDAC Patients (n = ", nrow(pca_scores_stage), ")"),
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal()

ggsave("plots_pca_stages/pca_proteomics_stage.png", stage_plot, width = 8, height = 5, dpi = 300)
print(stage_plot)
library(data.table)

phospho <- fread("data_csv/phosphoproteomics.csv", skip = 3)
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]

sample_type <- ifelse(substr(rownames(phospho), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- phospho
pca_data <- pca_data[, colMeans(is.na(pca_data)) < 0.7]
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
str(pca_data[, 1:5])
sd_values <- apply(pca_data, 2, sd)
sum(is.na(sd_values))
sd_values <- apply(pca_data, 2, sd)
sum(is.na(sd_values))
sum(sd_values == 0, na.rm = TRUE)
nrow(pca_data)
pca_data <- as.data.frame(sapply(pca_data, as.numeric))
rownames(pca_data) <- rownames(phospho)

sd_values <- apply(pca_data, 2, sd)
sum(is.na(sd_values))
bad_cols <- which(is.na(sd_values))
length(bad_cols)

head(phospho[, bad_cols[1]], 10)
colnames(pca_data)[bad_cols[1]]
pca_data <- as.data.frame(sapply(phospho, as.numeric))
rownames(pca_data) <- rownames(phospho)

pca_data[is.na(pca_data)] <- 0

sd_values <- apply(pca_data, 2, sd)
sum(is.na(sd_values))
pca_matrix <- as.matrix(sapply(phospho, as.numeric))
rownames(pca_matrix) <- rownames(phospho)

pca_matrix[is.na(pca_matrix)] <- 0

sd_values <- apply(pca_matrix, 2, sd)
sum(is.na(sd_values))
library(data.table)

phospho <- fread("data_csv/phosphoproteomics.csv", skip = 3)
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]

sample_type <- ifelse(substr(rownames(phospho), 1, 3) == "C3L", "Tumor", "Normal")

pca_matrix <- as.matrix(sapply(phospho, as.numeric))
rownames(pca_matrix) <- rownames(phospho)

pca_matrix[is.na(pca_matrix)] <- 0
pca_matrix <- pca_matrix[, apply(pca_matrix, 2, sd) > 0]
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")
library(data.table)

phospho <- fread("data_csv/phosphoproteomics.csv", skip = 3)
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]

sample_type <- ifelse(substr(rownames(phospho), 1, 3) == "C3L", "Tumor", "Normal")

pca_matrix <- as.matrix(sapply(phospho, as.numeric))
rownames(pca_matrix) <- rownames(phospho)

pca_matrix[is.na(pca_matrix)] <- 0
pca_matrix <- pca_matrix[, apply(pca_matrix, 2, sd) > 0]
pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
sd_check <- apply(pca_matrix, 2, sd)
sum(sd_check == 0, na.rm = TRUE)
sum(is.na(sd_check))
ncol(pca_matrix)
pca_matrix <- pca_matrix[, !is.na(sd_check) & sd_check > 0]

ncol(pca_matrix)
pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
pca_matrix <- pca_matrix[rownames(pca_matrix) != "Patient_ID", ]
sample_type <- ifelse(substr(rownames(pca_matrix), 1, 3) == "C3L", "Tumor", "Normal")
pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- rownames(pca_scores_tumor)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)


table(pca_scores_stage$Stage_Label)
clinical <- read.csv("data_csv/clinical.csv", row.names = 1)
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- rownames(pca_scores_tumor)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)

table(pca_scores_stage$Stage_Label)
head(rownames(pca_scores_stage))
nrow(phospho)
library(ggplot2)

stage_plot <- ggplot(pca_scores_stage, aes(x = PC1, y = PC2, color = Stage_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Scatter Plot - Phosphoproteomics - Coloured by Pathological Stage",
       subtitle = paste0("CPTAC PDAC Patients (n = ", nrow(pca_scores_stage), ")"),
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal()

ggsave("plots_pca_stages/pca_phosphoproteomics_stage.png", stage_plot, width = 8, height = 5, dpi = 300)
print(stage_plot)
miRNA <- read.csv("data_csv/miRNA.csv", row.names = 1)

sample_type <- ifelse(substr(rownames(miRNA), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- miRNA
pca_data <- na.omit(pca_data)
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- rownames(pca_scores_tumor)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)

table(pca_scores_stage$Stage_Label)
identical(sort(rownames(pca_scores_stage)), sort(patient_id[patient_id %in% rownames(clinical)]))
length(unique(patient_id))
proteomics_patients <- rownames(proteomics)[substr(rownames(proteomics), 1, 3) == "C3L"]
miRNA_patients <- rownames(miRNA)[substr(rownames(miRNA), 1, 3) == "C3L"]

length(intersect(proteomics_patients, miRNA_patients))
length(proteomics_patients)
length(miRNA_patients)
miRNA_patients <- rownames(miRNA)[substr(rownames(miRNA), 1, 3) == "C3L"]
length(miRNA_patients)

length(intersect(miRNA_patients, rownames(clinical)))
library(ggplot2)

stage_plot <- ggplot(pca_scores_stage, aes(x = PC1, y = PC2, color = Stage_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Scatter Plot - miRNA - Coloured by Pathological Stage",
       subtitle = paste0("CPTAC PDAC Patients (n = ", nrow(pca_scores_stage), ")"),
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal()

ggsave("plots_pca_stages/pca_miRNA_stage.png", stage_plot, width = 8, height = 5, dpi = 300)
print(stage_plot)
circular_RNA <- read.csv("data_csv/circular_RNA.csv", row.names = 1, skip = 4)

sample_type <- ifelse(substr(rownames(circular_RNA), 1, 3) == "C3L", "Tumor", "Normal")

pca_data <- circular_RNA
pca_data <- na.omit(pca_data)
pca_data <- pca_data[, apply(pca_data, 2, sd) > 0]
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
sample_type <- ifelse(substr(rownames(pca_data), 1, 3) == "C3L", "Tumor", "Normal")

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleType <- sample_type

head(pca_scores[, 1:3])
pca_scores_tumor <- pca_scores[pca_scores$SampleType == "Tumor", ]
patient_id <- substr(rownames(pca_scores_tumor), 1, 9)
pca_scores_tumor$Stage <- clinical$tumor_stage_pathological[match(patient_id, rownames(clinical))]

pca_scores_stage <- pca_scores_tumor[!is.na(pca_scores_tumor$Stage), ]

stage_counts <- table(pca_scores_stage$Stage)

pca_scores_stage$Stage_Label <- factor(
  pca_scores_stage$Stage,
  levels = names(stage_counts),
  labels = paste0(names(stage_counts), " (n=", stage_counts, ")")
)

table(pca_scores_stage$Stage_Label)
library(ggplot2)

stage_plot <- ggplot(pca_scores_stage, aes(x = PC1, y = PC2, color = Stage_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Scatter Plot - circular_RNA - Coloured by Pathological Stage",
       subtitle = paste0("CPTAC PDAC Patients (n = ", nrow(pca_scores_stage), ")"),
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal()

ggsave("plots_pca_stages/pca_circular_RNA_stage.png", stage_plot, width = 8, height = 5, dpi = 300)
print(stage_plot)
