"MOFA2" %in% rownames(installed.packages())
"reticulate" %in% rownames(installed.packages())
"basilisk" %in% rownames(installed.packages())
library(MOFA2)
library(ggplot2)
library(dplyr)
library(data.table)
library(reticulate)
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")

exists("tumor_purity")
exists("xcell")
exists("cibersort")
exists("miRNA")
exists("circular_RNA")
exists("proteomics")
exists("CNV")
exists("transcriptomics_num")
exists("phospho_num")
exists("clinical")
tumor_purity <- read.csv("data_csv/tumor_purity.csv", row.names = 1)
xcell <- read.csv("data_csv/xcell.csv", row.names = 1)
cibersort <- read.csv("data_csv/cibersort.csv", row.names = 1)
miRNA <- read.csv("data_csv/miRNA.csv", row.names = 1)
circular_RNA <- read.csv("data_csv/circular_RNA.csv", row.names = 1, skip = 4)
proteomics <- read.csv("data_csv/proteomics.csv", row.names = 1, skip = 1)
CNV <- read.csv("data_csv/CNV.csv", row.names = 1, skip = 1)
clinical <- read.csv("data_csv/clinical.csv", row.names = 1)
transcriptomics <- fread("data_csv/transcriptomics.csv")
transcriptomics <- as.data.frame(transcriptomics)
rownames(transcriptomics) <- transcriptomics[, 1]
transcriptomics <- transcriptomics[, -1]
transcriptomics_num <- as.data.frame(sapply(transcriptomics, as.numeric))
rownames(transcriptomics_num) <- rownames(transcriptomics)

phospho <- fread("data_csv/phosphoproteomics.csv")
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]
phospho_num <- as.data.frame(sapply(phospho, as.numeric))
rownames(phospho_num) <- rownames(phospho)
dim(tumor_purity)
dim(xcell)
dim(cibersort)
dim(miRNA)
dim(circular_RNA)
dim(proteomics)
dim(CNV)
dim(transcriptomics_num)
dim(phospho_num)
dim(clinical)
class(unlist(phospho_num))
common_samples <- Reduce(intersect, list(
  rownames(tumor_purity), rownames(xcell), rownames(cibersort),
  rownames(miRNA), rownames(circular_RNA), rownames(proteomics),
  rownames(CNV), rownames(transcriptomics_num), rownames(phospho_num)
))
length(common_samples)
table(substr(common_samples, 1, 3))
tumor_purity_f <- tumor_purity[common_samples, ]
xcell_f <- xcell[common_samples, ]
cibersort_f <- cibersort[common_samples, ]
miRNA_f <- miRNA[common_samples, ]
circular_RNA_f <- circular_RNA[common_samples, ]
proteomics_f <- proteomics[common_samples, ]
CNV_f <- CNV[common_samples, ]
transcriptomics_f <- transcriptomics_num[common_samples, ]
phospho_f <- phospho_num[common_samples, ]

dim(tumor_purity_f)
dim(phospho_f)
tumor_purity_log <- log2(as.matrix(tumor_purity_f) + 1)
xcell_log <- log2(as.matrix(xcell_f) + 1)
cibersort_log <- log2(as.matrix(cibersort_f) + 1)
miRNA_log <- log2(as.matrix(miRNA_f) + 1)
circular_RNA_log <- log2(as.matrix(circular_RNA_f) + 1)
transcriptomics_log <- log2(as.matrix(transcriptomics_f) + 1)

proteomics_mat <- as.matrix(proteomics_f)
CNV_mat <- as.matrix(CNV_f)
phospho_mat <- as.matrix(phospho_f)
tumor_purity_t <- t(tumor_purity_log)
xcell_t <- t(xcell_log)
cibersort_t <- t(cibersort_log)
miRNA_t <- t(miRNA_log)
circular_RNA_t <- t(circular_RNA_log)
transcriptomics_t <- t(transcriptomics_log)
proteomics_t <- t(proteomics_mat)
CNV_t <- t(CNV_mat)
phospho_t <- t(phospho_mat)

dim(tumor_purity_t)
dim(phospho_t)
data_list <- list(
  tumor_purity = tumor_purity_t,
  xcell = xcell_t,
  cibersort = cibersort_t,
  miRNA = miRNA_t,
  circular_RNA = circular_RNA_t,
  proteomics = proteomics_t,
  CNV = CNV_t,
  transcriptomics = transcriptomics_t,
  phosphoproteomics = phospho_t
)

length(data_list)
names(data_list)
sapply(data_list, dim)
MOFAobject <- create_mofa(data_list)
MOFAobject
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)
model_opts$num_factors <- 15
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42
MOFAobject <- prepare_mofa(MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAobject_trained <- run_mofa(MOFAobject, outfile = "MOFA_PDAC_model.hdf5", use_basilisk = TRUE)
reticulate::py_config()
MOFAobject_trained <- run_mofa(MOFAobject, outfile = "MOFA_PDAC_model.hdf5", use_basilisk = TRUE)
reticulate::py_install("mofapy2", pip = TRUE)
MOFAobject_trained <- run_mofa(MOFAobject, outfile = "MOFA_PDAC_model.hdf5", use_basilisk = FALSE)
clean_view <- function(mat) {
  mat[is.infinite(mat)] <- NA
  all_na_rows <- apply(mat, 1, function(x) all(is.na(x)))
  mat <- mat[!all_na_rows, ]
  row_var <- apply(mat, 1, function(x) var(x, na.rm = TRUE))
  mat <- mat[!is.na(row_var) & row_var > 0, ]
  return(mat)
}

data_list_clean <- lapply(data_list, clean_view)
sapply(data_list_clean, dim)
MOFAobject <- create_mofa(data_list_clean)

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)

model_opts$num_factors <- 15
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42

MOFAobject <- prepare_mofa(MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAobject_trained <- run_mofa(MOFAobject, outfile = "MOFA_PDAC_model.hdf5", use_basilisk = FALSE)
MOFAobject_trained
sample_metadata <- data.frame(
  sample = samples_names(MOFAobject_trained)[[1]],
  tumor_type = ifelse(substr(samples_names(MOFAobject_trained)[[1]], 1, 3) == "C3L", "Tumor", "Normal"),
  stage = clinical$tumor_stage_pathological[match(samples_names(MOFAobject_trained)[[1]], rownames(clinical))],
  age = clinical$age[match(samples_names(MOFAobject_trained)[[1]], rownames(clinical))],
  sex = clinical$sex[match(samples_names(MOFAobject_trained)[[1]], rownames(clinical))],
  vital_status = clinical$vital_status_at_date_of_last_contact[match(samples_names(MOFAobject_trained)[[1]], rownames(clinical))]
)

samples_metadata(MOFAobject_trained) <- sample_metadata
head(sample_metadata)
table(sample_metadata$tumor_type)
p_var <- plot_variance_explained(MOFAobject_trained, x = "view", y = "factor")
ggsave("plots_mofa/variance_explained_per_factor.png", p_var, width = 10, height = 8, dpi = 300)
print(p_var)
p_var_total <- plot_variance_explained(MOFAobject_trained, plot_total = TRUE)[[2]]
ggsave("plots_mofa/total_variance_per_view.png", p_var_total, width = 10, height = 6, dpi = 300)
print(p_var_total)
p_scatter_tumor <- plot_factors(MOFAobject_trained, 
                                factors = c(1, 2), 
                                color_by = "tumor_type",
                                dot_size = 3) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "MOFA Factor 1 vs Factor 2 - by Tumor Type") +
  theme_minimal()

ggsave("plots_mofa/factor_scatter_tumor_type.png", p_scatter_tumor, width = 8, height = 6, dpi = 300)
print(p_scatter_tumor)

p_scatter_stage <- plot_factors(MOFAobject_trained, 
                                factors = c(1, 2), 
                                color_by = "stage",
                                dot_size = 3) +
  scale_color_manual(values = c("Stage I" = "#16A085",
                                "Stage II" = "#E67E22",
                                "Stage III" = "#8E44AD",
                                "Stage IV" = "#C0392B",
                                "Staging is not applicable or unknown" = "gray")) +
  labs(title = "MOFA Factor 1 vs Factor 2 - by Cancer Stage") +
  theme_minimal()

ggsave("plots_mofa/factor_scatter_stage.png", p_scatter_stage, width = 8, height = 6, dpi = 300)
print(p_scatter_stage)
p_violin <- plot_factor(MOFAobject_trained, 
                        factors = 1:5,
                        color_by = "tumor_type",
                        add_violin = TRUE,
                        dodge = TRUE) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  scale_fill_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085"))

ggsave("plots_mofa/factor_violin_tumor_type.png", p_violin, width = 12, height = 6, dpi = 300)
print(p_violin)
p_weights_F1 <- plot_top_weights(MOFAobject_trained, 
                                 view = "transcriptomics", 
                                 factor = 1, 
                                 nfeatures = 15)
ggsave("plots_mofa/top_weights_factor1_transcriptomics.png", p_weights_F1, width = 8, height = 6, dpi = 300)
print(p_weights_F1)

p_weights_F2_prot <- plot_top_weights(MOFAobject_trained, 
                                      view = "proteomics", 
                                      factor = 2, 
                                      nfeatures = 15)
ggsave("plots_mofa/top_weights_factor2_proteomics.png", p_weights_F2_prot, width = 8, height = 6, dpi = 300)
print(p_weights_F2_prot)
p_corr <- correlate_factors_with_covariates(MOFAobject_trained,
                                            covariates = c("tumor_type", "stage", "age", "sex", "vital_status"),
                                            plot = "log_pval")
ggsave("plots_mofa/factor_metadata_correlation.png", p_corr, width = 8, height = 8, dpi = 300)
print(p_corr)
p_corr <- correlate_factors_with_covariates(MOFAobject_trained,
                                            covariates = c("tumor_type", "stage", "age", "sex", "vital_status"),
                                            plot = "log_pval") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white", color = "white"))

ggsave("plots_mofa/factor_metadata_correlation.png", p_corr, width = 8, height = 8, dpi = 300, bg = "white")
library(ggplot2)

p_corr <- correlate_factors_with_covariates(MOFAobject_trained,
                                            covariates = c("tumor_type", "stage", "age", "sex", "vital_status"),
                                            plot = "log_pval")

ggsave("plots_mofa/factor_metadata_correlation.png", p_corr, width = 8, height = 8, dpi = 300, bg = "white")
p_heatmap_F1 <- plot_data_heatmap(MOFAobject_trained,
                                  view = "transcriptomics",
                                  factor = 1,
                                  features = 25,
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  show_colnames = FALSE)
ggsave("plots_mofa/heatmap_factor1_transcriptomics.png", p_heatmap_F1, width = 10, height = 7, dpi = 300, bg = "white")

p_heatmap_F2 <- plot_data_heatmap(MOFAobject_trained,
                                  view = "proteomics",
                                  factor = 2,
                                  features = 25,
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  show_colnames = FALSE)
ggsave("plots_mofa/heatmap_factor2_proteomics.png", p_heatmap_F2, width = 10, height = 7, dpi = 300, bg = "white")
p_beeswarm <- plot_factor(MOFAobject_trained,
                          factors = 1:5,
                          color_by = "tumor_type",
                          dot_size = 2,
                          dodge = TRUE,
                          add_violin = TRUE) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  scale_fill_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085"))

ggsave("plots_mofa/factor_beeswarm.png", p_beeswarm, width = 12, height = 6, dpi = 300, bg = "white")
print(p_beeswarm)
save(MOFAobject_trained, sample_metadata, file = "MOFA_PDAC_results.RData")
data_list_v2 <- list(
  circular_RNA = circular_RNA_t,
  CNV = CNV_t,
  miRNA = miRNA_t,
  phosphoproteomics = phospho_t,
  proteomics = proteomics_t,
  transcriptomics = transcriptomics_t
)

length(data_list_v2)
names(data_list_v2)
sapply(data_list_v2, dim)
circular_RNA <- read.csv("data_csv/circular_RNA.csv", row.names = 1, skip = 4)
CNV <- read.csv("data_csv/CNV.csv", row.names = 1, skip = 1)
miRNA <- read.csv("data_csv/miRNA.csv", row.names = 1)
proteomics <- read.csv("data_csv/proteomics.csv", row.names = 1, skip = 1)

transcriptomics <- fread("data_csv/transcriptomics.csv", skip = 1)
transcriptomics <- as.data.frame(transcriptomics)
rownames(transcriptomics) <- transcriptomics[, 1]
transcriptomics <- transcriptomics[, -1]
transcriptomics_num <- as.data.frame(sapply(transcriptomics, as.numeric))
rownames(transcriptomics_num) <- rownames(transcriptomics)

phospho <- fread("data_csv/phosphoproteomics.csv", skip = 1)
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]
phospho_num <- as.data.frame(sapply(phospho, as.numeric))
rownames(phospho_num) <- rownames(phospho)

clinical <- read.csv("data_csv/clinical.csv", row.names = 1)

common_samples <- Reduce(intersect, list(
  rownames(circular_RNA), rownames(CNV), rownames(miRNA),
  rownames(proteomics), rownames(transcriptomics_num), rownames(phospho_num)
))
length(common_samples)
circular_RNA_f <- circular_RNA[common_samples, ]
CNV_f <- CNV[common_samples, ]
miRNA_f <- miRNA[common_samples, ]
proteomics_f <- proteomics[common_samples, ]
transcriptomics_f <- transcriptomics_num[common_samples, ]
phospho_f <- phospho_num[common_samples, ]

circular_RNA_log <- log2(as.matrix(circular_RNA_f) + 1)
miRNA_log <- log2(as.matrix(miRNA_f) + 1)
transcriptomics_log <- log2(as.matrix(transcriptomics_f) + 1)

CNV_mat <- as.matrix(CNV_f)
proteomics_mat <- as.matrix(proteomics_f)
phospho_mat <- as.matrix(phospho_f)

circular_RNA_t <- t(circular_RNA_log)
CNV_t <- t(CNV_mat)
miRNA_t <- t(miRNA_log)
proteomics_t <- t(proteomics_mat)
transcriptomics_t <- t(transcriptomics_log)
phospho_t <- t(phospho_mat)
data_list_v2 <- list(
  circular_RNA = circular_RNA_t,
  CNV = CNV_t,
  miRNA = miRNA_t,
  proteomics = proteomics_t,
  transcriptomics = transcriptomics_t,
  phosphoproteomics = phospho_t
)

length(data_list_v2)
names(data_list_v2)
sapply(data_list_v2, dim)
clean_view <- function(mat) {
  mat[is.infinite(mat)] <- NA
  all_na_rows <- apply(mat, 1, function(x) all(is.na(x)))
  mat <- mat[!all_na_rows, ]
  row_var <- apply(mat, 1, function(x) var(x, na.rm = TRUE))
  mat <- mat[!is.na(row_var) & row_var > 0, ]
  return(mat)
}

data_list_v2_clean <- lapply(data_list_v2, clean_view)
sapply(data_list_v2_clean, dim)

MOFAobject_v2 <- create_mofa(data_list_v2_clean)
MOFAobject_v2
library(MOFA2)
MOFAobject_v2 <- create_mofa(data_list_v2_clean)
MOFAobject_v2
data_opts <- get_default_data_options(MOFAobject_v2)
model_opts <- get_default_model_options(MOFAobject_v2)
train_opts <- get_default_training_options(MOFAobject_v2)

model_opts$num_factors <- 15
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42

MOFAobject_v2 <- prepare_mofa(MOFAobject_v2, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAobject_v2_trained <- run_mofa(MOFAobject_v2, outfile = "MOFA_PDAC_v2_model.hdf5", use_basilisk = TRUE)
library(reticulate)
reticulate::py_config()
reticulate::py_install("mofapy2", pip = TRUE)
reticulate::py_require(c("mofapy2"))
MOFAobject_v2_trained <- run_mofa(MOFAobject_v2, outfile = "MOFA_PDAC_v2_model.hdf5", use_basilisk = FALSE)
sample_metadata_v2 <- data.frame(
  sample = samples_names(MOFAobject_v2_trained)[[1]],
  tumor_type = ifelse(substr(samples_names(MOFAobject_v2_trained)[[1]], 1, 3) == "C3L", "Tumor", "Normal"),
  stage = clinical$tumor_stage_pathological[match(samples_names(MOFAobject_v2_trained)[[1]], rownames(clinical))],
  age = clinical$age[match(samples_names(MOFAobject_v2_trained)[[1]], rownames(clinical))],
  sex = clinical$sex[match(samples_names(MOFAobject_v2_trained)[[1]], rownames(clinical))],
  vital_status = clinical$vital_status_at_date_of_last_contact[match(samples_names(MOFAobject_v2_trained)[[1]], rownames(clinical))]
)

samples_metadata(MOFAobject_v2_trained) <- sample_metadata_v2
head(sample_metadata_v2)
table(sample_metadata_v2$tumor_type)
p_overview <- plot_data_overview(MOFAobject_v2_trained)
ggsave("plots_mofa_refined/data_overview.png", p_overview, width = 10, height = 6, dpi = 300)
print(p_overview)
p_var <- plot_variance_explained(MOFAobject_v2_trained, x = "view", y = "factor")
ggsave("plots_mofa_refined/variance_explained_per_factor.png", p_var, width = 10, height = 8, dpi = 300)
print(p_var)
p_var_total <- plot_variance_explained(MOFAobject_v2_trained, plot_total = TRUE)[[2]]
ggsave("plots_mofa_refined/total_variance_per_view.png", p_var_total, width = 10, height = 6, dpi = 300)
print(p_var_total)
p_scatter_tumor_v2 <- plot_factors(MOFAobject_v2_trained, 
                                   factors = c(1, 2), 
                                   color_by = "tumor_type",
                                   dot_size = 3) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "MOFA Factor 1 vs Factor 2 - by Tumor Type (Refined 6-view model)") +
  theme_minimal()

ggsave("plots_mofa_refined/factor_scatter_tumor_type.png", p_scatter_tumor_v2, width = 8, height = 6, dpi = 300)
print(p_scatter_tumor_v2)
class(sample_metadata_v2$tumor_type)
unique(sample_metadata_v2$tumor_type)
head(samples_metadata(MOFAobject_v2_trained))
unique(samples_metadata(MOFAobject_v2_trained)$tumor_type)
factors_v2 <- get_factors(MOFAobject_v2_trained, factors = "all")[[1]]
factors_df <- as.data.frame(factors_v2)
factors_df$sample <- rownames(factors_df)
factors_df <- merge(factors_df, sample_metadata_v2, by = "sample")

p_scatter_tumor_v2 <- ggplot(factors_df, aes(x = Factor1, y = Factor2, color = tumor_type)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "MOFA Factor 1 vs Factor 2 - by Tumor Type (Refined 6-view model)",
       x = "Factor 1", y = "Factor 2", color = "Sample Type") +
  theme_minimal()

ggsave("plots_mofa_refined/factor_scatter_tumor_type.png", p_scatter_tumor_v2, width = 8, height = 6, dpi = 300)
print(p_scatter_tumor_v2)
p_scatter_stage_v2 <- ggplot(factors_df, aes(x = Factor1, y = Factor2, color = stage)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Stage I" = "#16A085",
                                "Stage II" = "#E67E22",
                                "Stage III" = "#8E44AD",
                                "Stage IV" = "#C0392B",
                                "Staging is not applicable or unknown" = "gray")) +
  labs(title = "MOFA Factor 1 vs Factor 2 - by Cancer Stage (Refined 6-view model)",
       x = "Factor 1", y = "Factor 2", color = "Stage") +
  theme_minimal()

ggsave("plots_mofa_refined/factor_scatter_stage.png", p_scatter_stage_v2, width = 8, height = 6, dpi = 300)
print(p_scatter_stage_v2)
p_violin_v2 <- plot_factor(MOFAobject_v2_trained, 
                           factors = 1:5,
                           color_by = "tumor_type",
                           add_violin = TRUE,
                           dodge = TRUE) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  scale_fill_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085"))

ggsave("plots_mofa_refined/factor_violin_tumor_type.png", p_violin_v2, width = 12, height = 6, dpi = 300)
print(p_violin_v2)
p_beeswarm_v2 <- plot_factor(MOFAobject_v2_trained,
                             factors = 1:5,
                             color_by = "tumor_type",
                             dot_size = 2,
                             dodge = TRUE,
                             add_violin = TRUE) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  scale_fill_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085"))

ggsave("plots_mofa_refined/factor_beeswarm.png", p_beeswarm_v2, width = 12, height = 6, dpi = 300, bg = "white")
print(p_beeswarm_v2)
p_corr_v2 <- correlate_factors_with_covariates(MOFAobject_v2_trained,
                                               covariates = c("tumor_type", "stage", "age", "sex", "vital_status"),
                                               plot = "log_pval") +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white", color = "white"))

ggsave("plots_mofa_refined/factor_metadata_correlation.png", p_corr_v2, width = 8, height = 8, dpi = 300, bg = "white")
print(p_corr_v2)
p_weights_F1_v2 <- plot_top_weights(MOFAobject_v2_trained, 
                                    view = "phosphoproteomics", 
                                    factor = 1, 
                                    nfeatures = 15) +
  labs(title = "Top 15 Weights - Factor 1 (Phosphoproteomics)")

ggsave("plots_mofa_refined/top_weights_factor1_phosphoproteomics.png", p_weights_F1_v2, width = 8, height = 6, dpi = 300)
print(p_weights_F1_v2)
p_weights_F2_v2 <- plot_top_weights(MOFAobject_v2_trained, 
                                    view = "proteomics", 
                                    factor = 2, 
                                    nfeatures = 15) +
  labs(title = "Top 15 Weights - Factor 2 (Proteomics)")

ggsave("plots_mofa_refined/top_weights_factor2_proteomics.png", p_weights_F2_v2, width = 8, height = 6, dpi = 300)
print(p_weights_F2_v2)
p_heatmap_F1_v2 <- plot_data_heatmap(MOFAobject_v2_trained,
                                     view = "phosphoproteomics",
                                     factor = 1,
                                     features = 25,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     show_colnames = FALSE)

png("plots_mofa_refined/heatmap_factor1_phosphoproteomics.png", width = 10, height = 7, units = "in", res = 300)
print(p_heatmap_F1_v2)
dev.off()

print(p_heatmap_F1_v2)
p_heatmap_F2_v2 <- plot_data_heatmap(MOFAobject_v2_trained,
                                     view = "proteomics",
                                     factor = 2,
                                     features = 25,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     show_colnames = FALSE)

png("plots_mofa_refined/heatmap_factor2_proteomics.png", width = 10, height = 7, units = "in", res = 300)
print(p_heatmap_F2_v2)
dev.off()

print(p_heatmap_F2_v2)
png("plots_mofa_refined/factor_metadata_correlation.png", width = 8, height = 8, units = "in", res = 300)
correlate_factors_with_covariates(MOFAobject_v2_trained,
                                  covariates = c("tumor_type", "stage", "age", "sex", "vital_status"),
                                  plot = "log_pval")
dev.off()
p_corr_v2_final <- correlate_factors_with_covariates(MOFAobject_v2_trained,
                                                     covariates = c("tumor_type", "stage", "age", "sex", "vital_status"),
                                                     plot = "log_pval")

png("plots_mofa_refined/factor_metadata_correlation.png", width = 8, height = 8, units = "in", res = 300)
print(p_corr_v2_final)
dev.off()
p_factor_corr <- plot_factor_cor(MOFAobject_v2_trained)
png("plots_mofa_refined/factor_correlation_matrix.png", width = 8, height = 8, units = "in", res = 300)
plot_factor_cor(MOFAobject_v2_trained)
dev.off()

plot_factor_cor(MOFAobject_v2_trained)
graphics.off()

p_corr_v2_final <- correlate_factors_with_covariates(MOFAobject_v2_trained,
                                                     covariates = c("tumor_type", "stage", "age", "sex", "vital_status"),
                                                     plot = "log_pval")

png("plots_mofa_refined/factor_metadata_correlation.png", width = 8, height = 8, units = "in", res = 300)
print(p_corr_v2_final)
dev.off()
set.seed(42)
MOFAobject_v2_trained <- run_umap(MOFAobject_v2_trained)
p_umap_v2 <- plot_dimred(MOFAobject_v2_trained,
                         method = "UMAP",
                         color_by = "tumor_type",
                         dot_size = 3) +
  labs(title = "UMAP of All MOFA Factors - by Tumor Type (Refined 6-view model)")

ggsave("plots_mofa_refined/umap_factors_tumor_type.png", p_umap_v2, width = 8, height = 6, dpi = 300)
print(p_umap_v2)
p_umap_v2 <- plot_dimred(MOFAobject_v2_trained,
                         method = "UMAP",
                         color_by = "tumor_type",
                         dot_size = 3) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "UMAP of All MOFA Factors - by Tumor Type (Refined 6-view model)")

ggsave("plots_mofa_refined/umap_factors_tumor_type.png", p_umap_v2, width = 8, height = 6, dpi = 300)
print(p_umap_v2)
umap_coords <- MOFAobject_v2_trained@dim_red$UMAP
umap_df <- merge(umap_coords, sample_metadata_v2, by = "sample")

p_umap_v2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = tumor_type)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Tumor" = "#E67E22", "Normal" = "#16A085")) +
  labs(title = "UMAP of All MOFA Factors - by Tumor Type (Refined 6-view model)",
       color = "Sample Type") +
  theme_minimal()

ggsave("plots_mofa_refined/umap_factors_tumor_type.png", p_umap_v2, width = 8, height = 6, dpi = 300)
print(p_umap_v2)
weights_F1 <- get_weights(MOFAobject_v2_trained, views = "phosphoproteomics", factors = 1)[[1]]
weights_df <- data.frame(feature = rownames(weights_F1), weight = weights_F1[,1])

feature_values <- data_list_v2_clean$phosphoproteomics
factor1_scores <- get_factors(MOFAobject_v2_trained, factors = 1)[[1]][,1]

pvals <- apply(feature_values, 1, function(x) {
  if (sum(!is.na(x)) > 3) {
    cor.test(x, factor1_scores[names(x)], method = "spearman")$p.value
  } else {
    NA
  }
})

weights_df$pval <- pvals[weights_df$feature]
weights_df$padj <- p.adjust(weights_df$pval, method = "BH")
weights_df$neg_log_padj <- -log10(weights_df$padj)

head(weights_df)
weights_df$neg_log_padj[is.infinite(weights_df$neg_log_padj)] <- max(weights_df$neg_log_padj[is.finite(weights_df$neg_log_padj)], na.rm = TRUE) + 1

weights_df$significant <- ifelse(weights_df$padj < 0.05 & abs(weights_df$weight) > 0.3, "Significant", "Not significant")

p_volcano_v2 <- ggplot(weights_df, aes(x = weight, y = neg_log_padj, color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "#E67E22", "Not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano Plot - Factor 1 (Phosphoproteomics)",
       x = "MOFA Weight (Factor 1)",
       y = "-log10(adjusted p-value)",
       color = "") +
  theme_minimal()

ggsave("plots_mofa_refined/volcano_factor1_phosphoproteomics.png", p_volcano_v2, width = 8, height = 6, dpi = 300)
print(p_volcano_v2)
weights_df$direction <- case_when(
  weights_df$padj < 0.05 & weights_df$weight > 0.3 ~ "Positive (sig)",
  weights_df$padj < 0.05 & weights_df$weight < -0.3 ~ "Negative (sig)",
  TRUE ~ "Not significant"
)

p_volcano_v2 <- ggplot(weights_df, aes(x = weight, y = neg_log_padj, color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Positive (sig)" = "#E67E22", 
                                "Negative (sig)" = "#16A085", 
                                "Not significant" = "grey80")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano Plot - Factor 1 (Phosphoproteomics)",
       x = "MOFA Weight (Factor 1)",
       y = "-log10(adjusted p-value)",
       color = "Direction") +
  theme_minimal()

ggsave("plots_mofa_refined/volcano_factor1_phosphoproteomics.png", p_volcano_v2, width = 8, height = 6, dpi = 300)
print(p_volcano_v2)
