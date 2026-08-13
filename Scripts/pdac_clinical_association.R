library(tidyverse)

clinical <- read_csv(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/data_csv/clinical.csv",
  show_col_types = FALSE
)

dim(clinical)
names(clinical)
table(clinical$tumor_stage_pathological, useNA = "ifany")

table(clinical$histologic_grade, useNA = "ifany")

table(clinical$`Recurrence status (1, yes; 0, no)`, useNA = "ifany")

table(clinical$`Survival status (1, dead; 0, alive)`, useNA = "ifany")
table(clinical$tumor_necrosis, useNA = "ifany")

table(clinical$margin_status, useNA = "ifany")

table(clinical$perineural_invasion, useNA = "ifany")

table(clinical$tumor_site, useNA = "ifany")

table(clinical$tumor_focality, useNA = "ifany")
table(clinical$pathologic_staging_primary_tumor_pt, useNA = "ifany")

table(clinical$pathologic_staging_regional_lymph_nodes_pn, useNA = "ifany")

table(clinical$pathologic_staging_distant_metastasis_pm, useNA = "ifany")
summary(clinical$age)

summary(clinical$tumor_size_cm)

summary(clinical$bmi)

summary(clinical$`Overall survival, days`)

summary(clinical$`Recurrence-free survival, days`)
clinical$pT_clean <- str_extract(
  toupper(clinical$pathologic_staging_primary_tumor_pt),
  "T[1-4]"
)

table(clinical$pT_clean, useNA = "ifany")
clinical$pN_clean <- str_extract(
  toupper(clinical$pathologic_staging_regional_lymph_nodes_pn),
  "N[0-2]"
)

table(clinical$pN_clean, useNA = "ifany")
clinical_mofa <- clinical %>%
  select(
    Patien,
    tumor_stage_pathological,
    histologic_grade,
    pT_clean,
    pN_clean,
    tumor_necrosis,
    margin_status,
    perineural_invasion,
    tumor_site,
    tumor_size_cm,
    age,
    `Recurrence status (1, yes; 0, no)`,
    `Survival status (1, dead; 0, alive)`,
    `Overall survival, days`
  )
dir.create(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis",
  showWarnings = FALSE
)

write_csv(
  clinical_mofa,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/clinical_mofa_variables.csv"
)
file.exists(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/clinical_mofa_variables.csv"
)
file.exists(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/MOFA_PDAC_v2_model.hdf5"
)
library(MOFA2)

MOFA_B <- load_model(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/MOFA_PDAC_v2_model.hdf5"
)

MOFA_B
mofa_samples <- samples_names(MOFA_B)[[1]]

clinical_mofa_matched <- clinical_mofa[
  match(mofa_samples, clinical_mofa$Patien),
]

clinical_mofa_matched$sample <- mofa_samples

dim(clinical_mofa_matched)
sum(is.na(clinical_mofa_matched$Patien))
mofa_factors <- get_factors(
  MOFA_B,
  factors = "all",
  as.data.frame = TRUE
)

head(mofa_factors)
mofa_factors_wide <- mofa_factors %>%
  select(sample, factor, value) %>%
  pivot_wider(
    names_from = factor,
    values_from = value
  )

dim(mofa_factors_wide)
id="d4n8vs"
mofa_clinical <- clinical_mofa_matched %>%
  left_join(
    mofa_factors_wide,
    by = "sample"
  )

dim(mofa_clinical)
stage_results <- lapply(
  paste0("Factor", 1:15),
  function(f) {
    model <- lm(
      mofa_clinical[[f]] ~ mofa_clinical$tumor_stage_pathological
    )
    
    p <- anova(model)[1, "Pr(>F)"]
    
    data.frame(
      factor = f,
      p_value = p
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(q_value)

stage_results
recurrence_results <- lapply(
  paste0("Factor", 1:15),
  function(f) {
    model <- lm(
      mofa_clinical[[f]] ~ mofa_clinical$`Recurrence status (1, yes; 0, no)`
    )
    
    p <- anova(model)[1, "Pr(>F)"]
    
    data.frame(
      factor = f,
      p_value = p
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(q_value)

recurrence_results
write_csv(
  stage_results,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/stage_MOFA_association.csv"
)
write_csv(
  recurrence_results,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/recurrence_MOFA_association.csv"
)
grade_results <- lapply(
  paste0("Factor", 1:15),
  function(f) {
    model <- lm(
      mofa_clinical[[f]] ~ mofa_clinical$histologic_grade
    )
    
    p <- anova(model)[1, "Pr(>F)"]
    
    data.frame(
      factor = f,
      p_value = p
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(q_value)

grade_results
survival_results <- lapply(
  paste0("Factor", 1:15),
  function(f) {
    model <- lm(
      mofa_clinical[[f]] ~ mofa_clinical$`Survival status (1, dead; 0, alive)`
    )
    
    p <- anova(model)[1, "Pr(>F)"]
    
    data.frame(
      factor = f,
      p_value = p
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(q_value)

survival_results
perineural_results <- lapply(
  paste0("Factor", 1:15),
  function(f) {
    model <- lm(
      mofa_clinical[[f]] ~ mofa_clinical$perineural_invasion
    )
    
    p <- anova(model)[1, "Pr(>F)"]
    
    data.frame(
      factor = f,
      p_value = p
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(q_value)

perineural_results
margin_results <- lapply(
  paste0("Factor", 1:15),
  function(f) {
    model <- lm(
      mofa_clinical[[f]] ~ mofa_clinical$margin_status
    )
    
    p <- anova(model)[1, "Pr(>F)"]
    
    data.frame(
      factor = f,
      p_value = p
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(q_value)

margin_results
necrosis_results <- lapply(
  paste0("Factor", 1:15),
  function(f) {
    model <- lm(
      mofa_clinical[[f]] ~ mofa_clinical$tumor_necrosis
    )
    
    p <- anova(model)[1, "Pr(>F)"]
    
    data.frame(
      factor = f,
      p_value = p
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(q_value)

necrosis_results
write_csv(
  grade_results,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/grade_MOFA_association.csv"
)
write_csv(
  survival_results,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/survival_MOFA_association.csv"
)
write_csv(
  perineural_results,
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/perineural_MOFA_association.csv"
)
p_stage <- ggplot(
  mofa_clinical,
  aes(
    x = tumor_stage_pathological,
    y = Factor5
  )
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(
    title = "MOFA Factor 5 across pathological stages",
    x = "Pathological stage",
    y = "MOFA Factor 5"
  ) +
  theme_classic()

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Pathological_Stage.png",
  p_stage,
  width = 8,
  height = 6,
  dpi = 300
)

p_stage
stage_plot_data <- mofa_clinical %>%
  filter(
    tumor_stage_pathological %in% c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  ) %>%
  mutate(
    Stage = factor(
      tumor_stage_pathological,
      levels = c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      )
    )
  )

stage_counts <- stage_plot_data %>%
  count(Stage) %>%
  mutate(
    label = paste0(Stage, "\n(n=", n, ")")
  )

p_stage <- ggplot(
  stage_plot_data,
  aes(
    x = Stage,
    y = Factor5,
    fill = Stage
  )
) +
  geom_violin(
    alpha = 0.65,
    trim = FALSE,
    color = "grey25"
  ) +
  geom_boxplot(
    width = 0.16,
    alpha = 0.9,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.65,
    color = "black"
  ) +
  scale_x_discrete(
    labels = stage_counts$label
  ) +
  scale_fill_brewer(
    palette = "Set2"
  ) +
  labs(
    title = "MOFA Factor 5 across PDAC Pathological Stages",
    subtitle = "Exploratory association of a latent multi-omics factor with tumour stage",
    x = "Pathological stage",
    y = "MOFA Factor 5"
  ) +
  annotate(
    "text",
    x = 2.5,
    y = max(stage_plot_data$Factor5, na.rm = TRUE) * 0.95,
    label = "ANOVA: p = 0.0091 | FDR q = 0.136",
    size = 4
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Pathological_Stage.png",
  p_stage,
  width = 9,
  height = 7,
  dpi = 300
)

p_stage
stage_plot_data <- mofa_clinical %>%
  filter(
    tumor_stage_pathological %in% c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  ) %>%
  mutate(
    Stage = factor(
      tumor_stage_pathological,
      levels = c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      )
    )
  )

p_stage <- ggplot(
  stage_plot_data,
  aes(
    x = Stage,
    y = Factor5,
    fill = Stage
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.70,
    color = "grey25",
    linewidth = 0.7
  ) +stage_plot_data <- mofa_clinical %>%
  filter(
    tumor_stage_pathological %in% c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  ) %>%
  mutate(
    Stage = factor(
      tumor_stage_pathological,
      levels = c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      )
    )
  )

p_stage <- ggplot(
  stage_plot_data,
  aes(
    x = Stage,
    y = Factor5,
    fill = Stage
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.70,
    color = "grey25",
    linewidth = 0.7
  ) +
  geom_boxplot(
    width = 0.16,
    fill = "white",
    color = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.055,
    size = 2.1,
    alpha = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Stage I" = "#66C2A5",
      "Stage II" = "#FC8D62",
      "Stage III" = "#8DA0CB",
      "Stage IV" = "#E78AC3"
    )
  ) +
  labs(
    title = "MOFA Factor 5 across PDAC Pathological Stages",
    subtitle = "Distribution of latent multi-omics factor scores across tumour stages",
    x = "Pathological stage",
    y = "MOFA Factor 5",
    fill = "Pathological stage"
  ) +
  annotate(
    "label",
    x = 2.5,
    y = 1.05,
    label = "ANOVA p = 0.0091\nFDR q = 0.136",
    size = 4,
    fontface = "bold",
    fill = "white",
    color = "black"
  ) +
  coord_cartesian(
    ylim = c(
      min(stage_plot_data$Factor5, na.rm = TRUE) - 0.05,
      1.15
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Pathological_Stage.png",
  p_stage,
  width = 10,
  height = 7,
  dpi = 300
)

p_stage

  geom_boxplot(
    width = 0.16,
    fill = "white",
    color = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.055,
    size = 2.1,
    alpha = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Stage I" = "#66C2A5",
      "Stage II" = "#FC8D62",
      "Stage III" = "#8DA0CB",
      "Stage IV" = "#E78AC3"
    )
  ) +
  labs(
    title = "MOFA Factor 5 across PDAC Pathological Stages",
    subtitle = "Distribution of latent multi-omics factor scores across tumour stages",
    x = "Pathological stage",
    y = "MOFA Factor 5",
    fill = "Pathological stage"
  ) +
  annotate(
    "label",
    x = 2.5,
    y = 1.05,
    label = "ANOVA p = 0.0091\nFDR q = 0.136",
    size = 4,
    fontface = "bold",
    fill = "white",
    color = "black"
  ) +
  coord_cartesian(
    ylim = c(
      min(stage_plot_data$Factor5, na.rm = TRUE) - 0.05,
      1.15
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Pathological_Stage.png",
  p_stage,
  width = 10,
  height = 7,
  dpi = 300
)

p_stage
grade_plot_data <- mofa_clinical %>%
  filter(
    !is.na(histologic_grade),
    histologic_grade %in% c(
      "G1 Well differentiated",
      "G2 Moderately differentiated",
      "G3 Poorly differentiated",
      "G4 Undifferentiated"
    )
  ) %>%
  mutate(
    Grade = factor(
      histologic_grade,
      levels = c(
        "G1 Well differentiated",
        "G2 Moderately differentiated",
        "G3 Poorly differentiated",
        "G4 Undifferentiated"
      )
    )
  )
grade_plot_data <- grade_plot_data %>%
  mutate(
    Grade_short = factor(
      Grade,
      levels = levels(Grade),
      labels = c(
        "G1\nWell differentiated",
        "G2\nModerately differentiated",
        "G3\nPoorly differentiated",
        "G4\nUndifferentiated"
      )
    )
  )
grade_counts <- grade_plot_data %>%
  count(Grade_short)

grade_counts
p_grade <- ggplot(
  grade_plot_data,
  aes(
    x = Grade_short,
    y = Factor5,
    fill = Grade_short
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_violin(
    data = grade_plot_data %>%
      filter(Grade != "G4 Undifferentiated"),
    trim = FALSE,
    alpha = 0.70,
    color = "grey25",
    linewidth = 0.7
  ) +
  geom_boxplot(
    data = grade_plot_data %>%
      filter(Grade != "G4 Undifferentiated"),
    width = 0.16,
    fill = "white",
    color = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.055,
    size = 2.1,
    alpha = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "G1\nWell differentiated" = "#66C2A5",
      "G2\nModerately differentiated" = "#FC8D62",
      "G3\nPoorly differentiated" = "#8DA0CB",
      "G4\nUndifferentiated" = "#E78AC3"
    )
  ) +
  labs(
    title = "MOFA Factor 5 across PDAC Histological Grades",
    subtitle = "Distribution of latent multi-omics factor scores across tumour grades",
    x = "Histological grade",
    y = "MOFA Factor 5",
    fill = "Histological grade"
  ) +
  annotate(
    "label",
    x = 2.5,
    y = max(grade_plot_data$Factor5, na.rm = TRUE) * 0.95,
    label = "ANOVA p = 0.0091\nFDR q = 0.136",
    size = 4,
    fontface = "bold",
    fill = "white",
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Histological_Grade.png",
  p_grade,
  width = 10,
  height = 7,
  dpi = 300
)

p_grade
p_stage <- ggplot(
  stage_plot_data,
  aes(
    x = Stage,
    y = Factor5,
    fill = Stage
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.70,
    color = "grey25",
    linewidth = 0.7
  ) +
  geom_boxplot(
    width = 0.16,
    fill = "white",
    color = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.055,
    size = 2.1,
    alpha = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Stage I" = "#66C2A5",
      "Stage II" = "#FC8D62",
      "Stage III" = "#8DA0CB",
      "Stage IV" = "#E78AC3"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Stage I\n(n = 13)",
      "Stage II\n(n = 55)",
      "Stage III\n(n = 27)",
      "Stage IV\n(n = 6)"
    )
  ) +
  labs(
    title = "MOFA Factor 5 across PDAC Pathological Stages",
    subtitle = "Distribution of latent multi-omics factor scores across tumour stages",
    x = "Pathological stage",
    y = "MOFA Factor 5",
    fill = "Pathological stage"
  ) +
  annotate(
    "label",
    x = 2.5,
    y = 1.05,
    label = "ANOVA p = 0.0091\nFDR q = 0.136",
    size = 4,
    fontface = "bold",
    fill = "white",
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Pathological_Stage.png",
  p_stage,
  width = 10,
  height = 7,
  dpi = 300
)

p_stage
p_grade <- ggplot(
  grade_plot_data,
  aes(
    x = Grade_short,
    y = Factor5,
    fill = Grade_short
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_violin(
    data = grade_plot_data %>%
      filter(Grade != "G4 Undifferentiated"),
    trim = FALSE,
    alpha = 0.70,
    color = "grey25",
    linewidth = 0.7
  ) +
  geom_boxplot(
    data = grade_plot_data %>%
      filter(Grade != "G4 Undifferentiated"),
    width = 0.16,
    fill = "white",
    color = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.055,
    size = 2.1,
    alpha = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "G1\nWell differentiated" = "#66C2A5",
      "G2\nModerately differentiated" = "#FC8D62",
      "G3\nPoorly differentiated" = "#8DA0CB",
      "G4\nUndifferentiated" = "#E78AC3"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "G1\nWell differentiated\n(n = 6)",
      "G2\nModerately differentiated\n(n = 72)",
      "G3\nPoorly differentiated\n(n = 26)",
      "G4\nUndifferentiated\n(n = 1)"
    )
  ) +
  labs(
    title = "MOFA Factor 5 across PDAC Histological Grades",
    subtitle = "Distribution of latent multi-omics factor scores across tumour grades",
    x = "Histological grade",
    y = "MOFA Factor 5",
    fill = "Histological grade"
  ) +
  annotate(
    "label",
    x = 2.5,
    y = max(grade_plot_data$Factor5, na.rm = TRUE) * 0.95,
    label = "ANOVA p = 0.0091\nFDR q = 0.136",
    size = 4,
    fontface = "bold",
    fill = "white",
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Histological_Grade.png",
  p_grade,
  width = 10,
  height = 7,
  dpi = 300
)

p_grade
survival_plot_data <- mofa_factors_wide %>%
  select(sample, Factor5) %>%
  left_join(
    clinical_mofa_matched %>%
      select(
        sample,
        `Survival status (1, dead; 0, alive)`
      ),
    by = "sample"
  ) %>%
  filter(
    !is.na(`Survival status (1, dead; 0, alive)`)
  ) %>%
  mutate(
    Survival_Status = factor(
      `Survival status (1, dead; 0, alive)`,
      levels = c(0, 1),
      labels = c("Alive", "Dead")
    )
  )

survival_counts <- survival_plot_data %>%
  count(Survival_Status)

survival_counts
survival_stats <- survival_results %>%
  filter(factor == "Factor5")

p_survival <- ggplot(
  survival_plot_data,
  aes(
    x = Survival_Status,
    y = Factor5,
    fill = Survival_Status
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.70,
    color = "grey25",
    linewidth = 0.7
  ) +
  geom_boxplot(
    width = 0.16,
    fill = "white",
    color = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.055,
    size = 2.1,
    alpha = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Alive" = "#66C2A5",
      "Dead" = "#FC8D62"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Alive\n(n = 28)",
      "Dead\n(n = 75)"
    )
  ) +
  labs(
    title = "MOFA Factor 5 across PDAC Survival Status",
    subtitle = "Distribution of latent multi-omics factor scores by survival status",
    x = "Survival status",
    y = "MOFA Factor 5",
    fill = "Survival status"
  ) +
  annotate(
    "label",
    x = 1.5,
    y = max(survival_plot_data$Factor5, na.rm = TRUE) * 0.95,
    label = paste0(
      "p = ",
      formatC(
        survival_stats$p_value,
        format = "f",
        digits = 4
      ),
      "\nFDR q = ",
      formatC(
        survival_stats$q_value,
        format = "f",
        digits = 3
      )
    ),
    size = 4,
    fontface = "bold",
    fill = "white",
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  "C:/Users/varsh/Documents/R/cptac_pdac_pca/clinical_association_analysis/Factor5_Survival_Status.png",
  p_survival,
  width = 10,
  height = 7,
  dpi = 300
)

p_survival
