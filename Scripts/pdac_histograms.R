library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)
setwd("C:/Users/varsh/Documents/R/cptac_pdac_pca")
circular_RNA <- read.csv("data_csv/circular_RNA.csv", row.names = 1, skip = 4)

vals_circ <- unlist(circular_RNA)
vals_circ <- vals_circ[!is.na(vals_circ)]
vals_circ <- log1p(vals_circ[vals_circ >= 0])

h1 <- ggplot(data.frame(value = vals_circ), aes(x = value)) +
  geom_histogram(bins = 50, fill = "#E67E22", color = "white") +
  labs(title = "circular_RNA", x = "log1p(Value)", y = "Frequency") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11))

print(h1)
CNV <- read.csv("data_csv/CNV.csv", row.names = 1, skip = 1)

vals_cnv <- unlist(CNV)
vals_cnv <- vals_cnv[!is.na(vals_cnv)]

h2 <- ggplot(data.frame(value = vals_cnv), aes(x = value)) +
  geom_histogram(bins = 50, fill = "#E67E22", color = "white") +
  labs(title = "CNV", x = "Value", y = "Frequency") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11))

print(h2)
miRNA <- fread("data_csv/miRNA.csv")
miRNA <- as.data.frame(miRNA)
rownames(miRNA) <- miRNA[, 1]
miRNA <- miRNA[, -1]

vals_mirna <- unlist(miRNA)
vals_mirna <- vals_mirna[!is.na(vals_mirna)]
vals_mirna <- log1p(vals_mirna[vals_mirna >= 0])

h3 <- ggplot(data.frame(value = vals_mirna), aes(x = value)) +
  geom_histogram(bins = 50, fill = "#E67E22", color = "white") +
  labs(title = "miRNA", x = "log1p(Value)", y = "Frequency") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11))

print(h3)
proteomics <- fread("data_csv/proteomics.csv", skip = 1)
proteomics <- as.data.frame(proteomics)
rownames(proteomics) <- proteomics[, 1]
proteomics <- proteomics[, -1]

vals_prot <- unlist(proteomics)
vals_prot <- vals_prot[!is.na(vals_prot)]

h4 <- ggplot(data.frame(value = vals_prot), aes(x = value)) +
  geom_histogram(bins = 50, fill = "#E67E22", color = "white") +
  labs(title = "proteomics", x = "Value", y = "Frequency") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11))

print(h4)
phospho <- fread("data_csv/phosphoproteomics.csv", skip = 1)
phospho <- as.data.frame(phospho)
rownames(phospho) <- phospho[, 1]
phospho <- phospho[, -1]

vals_phospho <- unlist(phospho)
phospho_num <- as.data.frame(sapply(phospho, as.numeric))

vals_phospho <- unlist(phospho_num)
vals_phospho <- vals_phospho[!is.na(vals_phospho)]

h5 <- ggplot(data.frame(value = vals_phospho), aes(x = value)) +
  geom_histogram(bins = 50, fill = "#E67E22", color = "white") +
  labs(title = "phosphoproteomics", x = "Value", y = "Frequency") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11))

print(h5)
transcriptomics <- fread("data_csv/transcriptomics.csv", skip = 1)
transcriptomics <- as.data.frame(transcriptomics)
rownames(transcriptomics) <- transcriptomics[, 1]
transcriptomics <- transcriptomics[, -1]
transcriptomics_num <- as.data.frame(sapply(transcriptomics, as.numeric))

vals_trans <- unlist(transcriptomics_num)
vals_trans <- vals_trans[!is.na(vals_trans)]
vals_trans <- log1p(vals_trans[vals_trans >= 0])

h6 <- ggplot(data.frame(value = vals_trans), aes(x = value)) +
  geom_histogram(bins = 50, fill = "#E67E22", color = "white") +
  labs(title = "transcriptomics", x = "log1p(Value)", y = "Frequency") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11))

print(h6)
histograms_all <- (h1 | h2 | h3) / (h4 | h5 | h6) +
  plot_annotation(
    title = "Distribution of the 6 CPTAC PDAC Omics Datasets",
    subtitle = "Log1p-transformed for zero-inflated layers (circular_RNA, miRNA, transcriptomics); raw scale for pre-normalized layers (CNV, proteomics, phosphoproteomics)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5))
  )

ggsave("plots_histograms/histograms_all_omics.png", histograms_all, width = 16, height = 10, dpi = 300)
print(histograms_all)
