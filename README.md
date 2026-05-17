# Multi-Omics Factor Analysis (MOFA) — CPTAC PDAC

Project analyzing pancreatic ductal adenocarcinoma (PDAC) using CPTAC multi-omics data. The goal is to perform PCA analysis on individual omics layers to assess Tumor vs Normal segregation, followed by integrative MOFA analysis.

## Workflow

1. **Data download** (`CPTAC_pdac_PCA.ipynb`) — Python notebook downloading all 15 PDAC datasets from the CPTAC repository.
2. **PCA analysis** (`CPTAC_pdac_PCA.R`) — Principal Component Analysis on 9 numeric omics datasets: tumor_purity, xcell, cibersort, miRNA, circular_RNA, proteomics, CNV, transcriptomics, phosphoproteomics. Output: scree and scatter plots colored by Tumor (orange) and Normal (teal).
3. **Metadata analysis** (`CPTAC_pdac_Metadata.R`) — Summary visualizations of clinical, ancestry, and HLA typing data.

