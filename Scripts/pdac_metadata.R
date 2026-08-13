library(ggplot2)
library(dplyr)
library(data.table)
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")
clinical <- read.csv("data/clinical.csv", row.names = 1)
dim(clinical)
colnames(clinical)
table(clinical$sex)
install.packages("patchwork")
library(patchwork)
table(clinical$sex)
table(clinical$tumor_stage_pathological)
table(clinical$histologic_grade)
table(clinical$tumor_site)
table(clinical$vital_status_at_date_of_last_contact)
summary(clinical$age)
summary(clinical$age)
p1 <- ggplot(clinical, aes(x = sex, fill = sex)) +
  geom_bar(width = 0.6) +
  scale_fill_manual(values = c("Female" = "#E67E22", "Male" = "#16A085")) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 4) +
  labs(title = "Sex", x = NULL, y = "Patients") +
  theme_minimal() + theme(legend.position = "none")

p2 <- ggplot(clinical, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "#E67E22", color = "white") +
  labs(title = "Age Distribution", x = "Age (years)", y = "Patients") +
  theme_minimal()

p3 <- ggplot(clinical, aes(x = tumor_stage_pathological, fill = tumor_stage_pathological)) +
  geom_bar(width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 3.5) +
  labs(title = "Tumor Stage", x = NULL, y = "Patients") +
  theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p4 <- ggplot(clinical, aes(x = histologic_grade, fill = histologic_grade)) +
  geom_bar(width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 3.5) +
  labs(title = "Histologic Grade", x = NULL, y = "Patients") +
  theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p5 <- ggplot(clinical, aes(x = tumor_site, fill = tumor_site)) +
  geom_bar(width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 4) +
  labs(title = "Tumor Site", x = NULL, y = "Patients") +
  theme_minimal() + theme(legend.position = "none")

p6 <- ggplot(clinical, aes(x = vital_status_at_date_of_last_contact, fill = vital_status_at_date_of_last_contact)) +
  geom_bar(width = 0.6) +
  scale_fill_manual(values = c("Deceased" = "#E67E22", "Living" = "#16A085")) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 4) +
  labs(title = "Vital Status", x = NULL, y = "Patients") +
  theme_minimal() + theme(legend.position = "none")
class(p1)
clinical_overview <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(title = "Clinical Overview - CPTAC PDAC",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

ggsave("plots/clinical_overview.png", clinical_overview, width = 14, height = 8, dpi = 300)
print(clinical_overview)
followup <- read.csv("data_csv/followup.csv")
dim(followup)
colnames(followup)
table(followup$is_this_patient_lost_to_follow.up)
medical_history <- read.csv("data_csv/medical_history.csv")
dim(medical_history)
colnames(medical_history)
ancestry <- read.csv("data_csv/ancestry_prediction.csv")
dim(ancestry)
colnames(ancestry)
table(ancestry$consensus_pred_ancestry)
hla <- read.csv("data_csv/hla_typing.csv")
dim(hla)
colnames(hla)
ancestry_data <- as.data.frame(table(ancestry$consensus_pred_ancestry))
colnames(ancestry_data) <- c("Ancestry", "Count")

bar_ancestry <- ggplot(ancestry_data, aes(x = reorder(Ancestry, -Count), y = Count, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  labs(title = "Patient Ancestry Distribution - CPTAC PDAC",
       subtitle = "Consensus predicted ancestry",
       x = "Ancestry Group",
       y = "Number of Patients") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave("plots_pca_metadata/ancestry_distribution.png", bar_ancestry, width = 8, height = 6, dpi = 300)
print(bar_ancestry)
hla <- read.csv("data_csv/hla_typing.csv")
table(hla$A1)
table(hla$B1)
table(hla$C1)
hla_A <- as.data.frame(table(hla$A1)) %>% top_n(10, Freq) %>% arrange(desc(Freq))
colnames(hla_A) <- c("HLA_Type", "Count")

hla_B <- as.data.frame(table(hla$B1)) %>% top_n(10, Freq) %>% arrange(desc(Freq))
colnames(hla_B) <- c("HLA_Type", "Count")

hla_C <- as.data.frame(table(hla$C1)) %>% top_n(10, Freq) %>% arrange(desc(Freq))
colnames(hla_C) <- c("HLA_Type", "Count")

hp1 <- ggplot(hla_A, aes(x = reorder(HLA_Type, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#E67E22", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  labs(title = "Top HLA-A Types", x = NULL, y = "Patients") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

hp2 <- ggplot(hla_B, aes(x = reorder(HLA_Type, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#16A085", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  labs(title = "Top HLA-B Types", x = NULL, y = "Patients") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

hp3 <- ggplot(hla_C, aes(x = reorder(HLA_Type, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#8E44AD", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  labs(title = "Top HLA-C Types", x = NULL, y = "Patients") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

hla_overview <- hp1 / hp2 / hp3 +
  plot_annotation(title = "Top 10 HLA Types - CPTAC PDAC Patients",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

ggsave("plots_pca_metadata/hla_top_types.png", hla_overview, width = 12, height = 10, dpi = 300)
print(hla_overview)
readLines("data_csv/somatic_mutation.csv", n = 10)
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")
clinical <- read.csv("data_csv/clinical.csv", row.names = 1)
table(clinical$tumor_stage_pathological)
numeric_cols <- sapply(clinical, is.numeric)
sum(numeric_cols)
clinical_numeric <- clinical[, numeric_cols]
colnames(clinical_numeric)
colSums(is.na(clinical_numeric))
clinical_pca_data <- clinical_numeric[, colSums(is.na(clinical_numeric)) < 30]
ncol(clinical_pca_data)
for(i in 1:ncol(clinical_pca_data)) {
  clinical_pca_data[is.na(clinical_pca_data[,i]), i] <- mean(clinical_pca_data[,i], na.rm = TRUE)
}
sum(is.na(clinical_pca_data))
pca_clinical <- prcomp(clinical_pca_data, center = TRUE, scale. = TRUE)
variance_explained <- (pca_clinical$sdev^2) / sum(pca_clinical$sdev^2) * 100
round(variance_explained[1:5], 2)
pca_scores <- as.data.frame(pca_clinical$x)
pca_scores$Stage <- clinical$tumor_stage_pathological[match(rownames(pca_scores), rownames(clinical))]
pca_scores_filtered <- pca_scores[pca_scores$Stage != "Staging is not applicable or unknown", ]
table(pca_scores_filtered$Stage)
scatter_clinical_stages <- ggplot(pca_scores_filtered, aes(x = PC1, y = PC2, color = Stage)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Stage I" = "#16A085",
                                "Stage II" = "#E67E22",
                                "Stage III" = "#8E44AD",
                                "Stage IV" = "#C0392B")) +
  labs(title = "PCA Scatter Plot - Clinical Data by Cancer Stage",
       subtitle = "CPTAC PDAC Patients (n = 135)",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11))

ggsave("plots_metadata/clinical_pca_stages.png", scatter_clinical_stages, width = 9, height = 6, dpi = 300)
print(scatter_clinical_stages)
library(ggplot2)
scatter_clinical_stages <- ggplot(pca_scores_filtered, aes(x = PC1, y = PC2, color = Stage)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Stage I" = "#16A085",
                                "Stage II" = "#E67E22",
                                "Stage III" = "#8E44AD",
                                "Stage IV" = "#C0392B")) +
  labs(title = "PCA Scatter Plot - Clinical Data by Cancer Stage",
       subtitle = "CPTAC PDAC Patients (n = 135)",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
       color = "Cancer Stage") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11))

ggsave("plots_metadata/clinical_pca_stages.png", scatter_clinical_stages, width = 9, height = 6, dpi = 300)
print(scatter_clinical_stages)
