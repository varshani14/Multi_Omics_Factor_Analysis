# Multi-Omics Factor Analysis (MOFA) of CPTAC PDAC

## Overview

This project explores pancreatic ductal adenocarcinoma (PDAC) using CPTAC multi-omics data.

The main aim is to move from individual-omics exploration to multi-omics integration using Multi-Omics Factor Analysis (MOFA). The workflow examines the individual data layers, integrates them through MOFA, investigates clinical associations, and performs downstream biological interpretation of an important MOFA factor.

The repository is organized so that the complete analysis can be followed from data acquisition and preprocessing through PCA, MOFA, clinical association analysis, and downstream biological interpretation.

---

## Analysis Workflow

```text
CPTAC PDAC data
       │
       ▼
Data preparation & metadata exploration
       │
       ├── Feature distributions
       ├── Individual-omics PCA
       │      ├── Scree plots
       │      └── Tumor vs Normal plots
       └── Stage-based PCA
       │
       ▼
MOFA Model A
Original 9-view model
       │
       ▼
MOFA Model B
Refined 6-view model
       │
       ├── Factor visualization
       ├── Variance explained
       ├── Factor correlations
       ├── Clinical/metadata associations
       ├── Feature weights
       ├── Heatmaps
       └── UMAP
       │
       ├───────────────┐
       ▼               ▼
Clinical           Downstream
Association        Biological
Analysis            Interpretation
                       │
                       ├── CNV annotation
                       ├── GO enrichment
                       ├── KEGG enrichment
                       └── GSEA
```

---

## Data

The project uses multiple molecular and clinical data layers derived from the CPTAC PDAC dataset.

### Original MOFA Model A — 9 views

- Tumor purity
- xCell
- CIBERSORT
- miRNA
- Circular RNA
- Proteomics
- CNV
- Transcriptomics
- Phosphoproteomics

### Refined MOFA Model B — 6 views

- Circular RNA
- CNV
- miRNA
- Phosphoproteomics
- Proteomics
- Transcriptomics

Model A is retained as the original 9-view analysis, while Model B represents the refined 6-view analysis used for the subsequent interpretation.

---

## Analysis

### 1. Data Download and Preparation

The data acquisition workflow is documented in:

`Scripts/cptac_pdac_data_download.ipynb`

The downloaded CSV datasets are stored in:

`data_csv/`

The R scripts then prepare, transform, and align the different data layers using common sample identifiers.

### 2. Metadata and Clinical Exploration

Clinical and metadata characteristics are explored before multi-omics integration.

This includes:

- Clinical variables
- Ancestry distribution
- HLA typing
- Clinical PCA
- Pathological-stage visualization

The corresponding outputs are organized under:

`plots_metadata/`

and

`clinical_association_analysis/`

### 3. Multi-Omics Feature Distributions

Feature distributions are examined across the different data views before PCA and MOFA integration.

The combined visualization is stored in:

`plots_histograms/`

### 4. Individual-Omics PCA

PCA is performed separately for the individual data views to examine:

- Variance captured by principal components
- Internal structure of each data layer
- Tumor vs Normal separation

The outputs are organized as:

```text
plots_pca/
├── scree/
└── scatter/
```

### 5. Stage-Based PCA

Additional PCA visualizations examine molecular variation across pathological cancer stages.

The outputs are stored in:

`plots_pca_stages/`

---

## MOFA Analysis

### MOFA Model A — Original 9-View Model

Model A represents the original multi-omics integration using nine data views.

Its outputs include:

- Variance explained
- Total variance per view
- Tumor vs Normal factor visualization
- Cancer-stage visualization
- Factor distributions
- Top feature weights
- Factor/metadata correlations
- Heatmaps
- Factor beeswarm visualization

The corresponding plots are stored in:

`mofa_model_A/`

### MOFA Model B — Refined 6-View Model

Model B is the refined six-view model and forms the main basis for the subsequent interpretation.

Its outputs include:

- Data overview
- Variance explained
- Tumor vs Normal factor visualization
- Cancer-stage visualization
- Factor violin and beeswarm plots
- Factor correlation matrix
- Factor/metadata correlation
- Top feature weights
- Feature heatmaps
- UMAP visualization
- Factor 1 phosphoproteomics analysis
- Factor 2 proteomics analysis

These outputs are stored in:

`mofa_model_B/`

---

## Clinical Association Analysis

The clinical analysis investigates associations between MOFA factor results and clinical characteristics, including:

- Histological grade
- Pathological stage
- Survival status
- Recurrence
- Perineural invasion

The results are separated into:

```text
clinical_association_analysis/
├── plots/
└── tables/
```

---

## Downstream Biological Interpretation

The downstream analysis focuses on the biological interpretation of **MOFA Factor 5**.

### CNV Annotation

CNV-associated features are annotated with gene symbols and gene names where available.

### Gene Ontology Enrichment

GO enrichment is used to investigate biological processes associated with the selected gene set.

### KEGG Pathway Enrichment

KEGG analysis is used to examine pathway-level associations.

### Gene Set Enrichment Analysis

GSEA is performed using a ranked gene list to investigate coordinated enrichment of biological gene sets.

The resulting outputs are organized as:

```text
downstream_analysis/
├── plots/
└── tables/
```

---

## Repository Structure

```text
Multi_Omics_Factor_Analysis/
│
├── Scripts/
│   ├── cptac_pdac_data_download.ipynb
│   ├── pdac_metadata.R
│   ├── pdac_histograms.R
│   ├── pdac_pca.R
│   ├── pdac_pca_stages.R
│   ├── pdac_mofa.R
│   ├── pdac_clinical_association.R
│   └── pdac_downstream_analysis.R
│
├── data_csv/
│
├── plots_metadata/
├── plots_histograms/
│
├── plots_pca/
│   ├── scree/
│   └── scatter/
│
├── plots_pca_stages/
│
├── mofa_model_A/
├── mofa_model_B/
│
├── clinical_association_analysis/
│   ├── plots/
│   └── tables/
│
└── downstream_analysis/
    ├── plots/
    └── tables/
```

---

## Reproducibility

The analysis is implemented primarily in R, with a Python notebook used for the initial data-download workflow.

The repository keeps the analysis scripts, input data tables, visual outputs, MOFA results, clinical association results, and downstream biological interpretation results together in a structured layout.

The scripts are separated by analytical stage so that the workflow can be followed from data preparation through multi-omics integration and biological interpretation.

---

