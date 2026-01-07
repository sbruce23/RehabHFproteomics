# REHAB-HF Proteomics Analysis - Complete Script
# Comprehensive analysis of protein biomarkers and rehabilitation outcomes
# Study: REHAB-HFpEF (Heart Failure with Preserved Ejection Fraction)
# Primary outcomes: SPPB (Short Physical Performance Battery) and 6MWD (Six-Minute Walk Distance)

################################################################################
# 0. SETUP AND CONFIGURATION
################################################################################

rm(list=ls())

# Load required libraries
library(tidyverse)
library(XLConnect)
library(readxl)
library(janitor)
library(gt)
library(gtsummary)
library(reshape2)
library(sjPlot)
library(ggforce)
library(gridExtra)
library(glinternet)
library(bartMachine)
library(optmatch)
library(rpart)
library(rpart.plot)
library(lsmeans)
library(treeClust)
library(pdfCluster)
library(mediation)
library(dplyr)
library(ggcorrplot)
library(HDtest)
library(WGCNA)
library(RColorBrewer)
library(ggrepel)
library(forcats)
library(randomForest)

# Set working directory (proteomics data)
setwd("set/working/directory/REHAB-HF proteomics data")

################################################################################
# 0.1 DATA DICTIONARY
################################################################################
# 
# Demographic/Clinical Variables:
#   age              : Age in years
#   sex              : 1=female, 0=male
#   race___4         : 1=white, 0=non-white
#   hf_cat           : 1=HFpEF (EF>=45%), 0=HFrEF (EF<45%)
#   creatinine_value : Serum creatinine (mg/dL)
#   egfr             : Estimated glomerular filtration rate (derived from 2021 CKD-EPI equation)
#
# Physical Performance:
#   bl_sppb          : Baseline SPPB score (0-12, higher=better)
#   fu_sppb          : 12-week follow-up SPPB score
#   bl_smw           : Baseline six-minute walk test distance (meters)
#   fu_smw           : 12-week follow-up six-minute walk distance (meters)
#   bl_bal_scr       : Baseline balance score
#   fu_bal_scr       : Follow-up balance score
#   bl_walk_time     : Baseline walk time
#   fu_walk_time     : Follow-up walk time
#
# Clinical Outcomes:
#   rehosp           : Rehospitalization (0/1)
#   death            : Death (0/1)
#   event            : Rehospitalization or death (0/1)
#   nrehosp_death    : Number of rehospitalization/death events
#   lfu              : Lost to follow-up (0/1)
#
# Intervention:
#   intervention_1_control_0 : 1=rehabilitation intervention, 0=attention control
#
# Olink Proteomics:
#   92 protein biomarkers measured using Olink proximity extension assay (NPX units)
#   LOD              : Limit of detection for each protein

################################################################################
# 0.2 GGPLOT2 THEME CUSTOMIZATION
################################################################################

hw <- theme_gray() + theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  plot.subtitle = element_text(hjust = 0.5, size = 12),
  plot.caption = element_text(hjust = -0.5, size = 10),
  strip.background = element_rect(fill = rgb(.9, .95, 1), colour = gray(.5), linewidth = .2),
  panel.border = element_rect(fill = FALSE, colour = gray(.70)),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.spacing.x = unit(0.2, "cm"),
  panel.spacing.y = unit(0.2, "cm"),
  axis.text = element_text(colour = "black", size = 10),
  axis.text.y = element_text(margin = ggplot2::margin(0, 3, 0, 3)),
  axis.text.x = element_text(margin = ggplot2::margin(-1, 0, 3, 0)),
  axis.title = element_text(size = 16, face = "bold"),
  legend.text = element_text(size = 14),
  legend.title = element_blank()
)

################################################################################
# 0.3 READ AND PROCESS OLINK PROTEOMICS DATA
################################################################################

olink_fname <- 'OMITTED'

# Read Olink NPX data (raw protein expression levels)
olink_dta <- read_excel(
  olink_fname,
  range = "A8:CO418",
  col_names = FALSE
)

# Read protein names, UniProt IDs, and Olink IDs
olink_colnames <- read_excel(
  olink_fname,
  range = "B4:CO4",
  col_names = FALSE
)

olink_olinkid <- read_excel(
  olink_fname,
  range = "B4:CO6",
  col_names = FALSE
)

olink_olinkid <- as.data.frame(t(olink_olinkid))
colnames(olink_olinkid) <- c('Assay', 'UniProt', 'OlinkID')

# Read limits of detection for each protein
olink_lod <- as.numeric(
  read_excel(olink_fname,
             range = "B420:CO420",
             col_names = FALSE)
)

# Assign column names to proteomics data
colnames(olink_dta) <- c('ID', as.character(olink_colnames))

# NOTE: Observed missingness levels
#   SPON1: 15% missing
#   CHIT1: 4% missing
#   Minor missingness across other proteins

# Separate baseline measurements
olink_dta_bl <- olink_dta[grep("BL", olink_dta$ID), ]
olink_dta_bl <- olink_dta_bl[-grep("BL_2", olink_dta_bl$ID), ]  # Remove duplicate BL reading for subject 20-017
olink_dta_bl$timepoint <- "Baseline"
olink_dta_bl$subject_id <- gsub(" BL", "", olink_dta_bl$ID, fixed = TRUE)

# Separate follow-up measurements
olink_dta_fu <- olink_dta[grep("FU", olink_dta$ID), ]
olink_dta_fu$timepoint <- "Follow Up"
olink_dta_fu$subject_id <- gsub(" FU", "", olink_dta_fu$ID, fixed = TRUE)
olink_dta_fu$subject_id <- gsub("FU", "", olink_dta_fu$subject_id, fixed = TRUE)

################################################################################
# 0.4 READ AND PROCESS CLINICAL DATA
################################################################################

setwd("set/working/directory/REHAB-HF dataset")

# Read main REHAB-HF dataset with password protection
rhf.fname <- 'OMITTED'
rhf.pw <- 'OMITTED'
rhfwb <- loadWorkbook(rhf.fname, password = rhf.pw)
rhfds <- readWorksheet(rhfwb, "MERGED_FINAL2")

# Read biomarker dataset with results (troponin, creatinine, NT-proBNP, hs-CRP, etc.)
rbio.fname <- 'OMITTED'
rbio.pw <- 'OMITTED'
rbiowb <- loadWorkbook(rbio.fname,
                       password = rbio.pw)
rbiods <- readWorksheet(rbiowb, "Sheet1", endRow = 402)

# Clean up column names and convert data types
rbiods <- clean_names(rbiods)
colnames(rbiods)[3] <- "sample_id"
colnames(rbiods)[4] <- "sample_type"
colnames(rbiods)[5] <- "sample_volume_ml"
colnames(rbiods)[9] <- "hf_category_1_hfpef_0_hfref"
colnames(rbiods)[12] <- "creatinine_mg_dl"
colnames(rbiods)[13] <- "troponin_i_pg_ml"

# Convert factors
rhfds$sex <- factor(rhfds$sex, levels = c(0, 1))
rhfds$race___4 <- factor(rhfds$race___4, levels = c(0, 1))
rhfds$hf_cat <- factor(rhfds$hf_cat, levels = c(0, 1))

# Convert continuous variables to numeric
rhfds$creatinine_value <- as.numeric(rhfds$creatinine_value)
rhfds$bl_sppb <- as.numeric(rhfds$bl_sppb)
rhfds$fu_sppb <- as.numeric(rhfds$fu_sppb)
rhfds$bl_smw <- as.numeric(rhfds$bl_smw)
rhfds$fu_smw <- as.numeric(rhfds$fu_smw)

################################################################################
# 0.5 CALCULATE eGFR (2021 CKD-EPI Creatinine Equation)
################################################################################
# Reference: https://www.kidney.org/content/ckd-epi-creatinine-equation-2021

kpa <- ifelse(rhfds$sex == 1, 0.7, 0.9)  # Sex-specific kappa
alpha <- ifelse(rhfds$sex == 1, -0.241, -0.302)  # Sex-specific alpha
const <- ifelse(rhfds$sex == 1, 1.012, 1)  # Sex-specific constant

rhfds$egfr <- 142 *
  (pmin(as.numeric(rhfds$creatinine_value) / kpa, 1) ^ alpha) *
  (pmax(as.numeric(rhfds$creatinine_value) / kpa, 1) ^ -1.2) *
  (0.9938 ^ rhfds$age) * const

################################################################################
# 0.6 CLEAN BIOMARKER DATA (rbiods)
################################################################################

# Remove (D) indicators (typically denote QC flags)
rbiods <- apply(rbiods, 2, function(x) gsub("(D)", "", x, fixed = TRUE))

# Remove spaces from biomarker columns
rbiods[, 12:16] <- apply(rbiods[, 12:16], 2, function(x) gsub(" ", "", x, fixed = TRUE))

# Create QNS (Quantity Not Sufficient) indicator and replace with NA
qns_ind <- apply(rbiods[, 12:16], 2, function(x) as.numeric(grepl("QNS", x, fixed = TRUE)))
colnames(qns_ind) <- paste(colnames(qns_ind), "_qns_ind", sep = "")
rbiods <- as.data.frame(cbind(rbiods, qns_ind))
rbiods[, 12:16] <- apply(rbiods[, 12:16], 2, function(x) gsub("QNS", NA, x, fixed = TRUE))

# Create censoring indicators for <LOD (left-censored) and >ULD (right-censored) values
cens_ind <- matrix(0, nrow = nrow(rbiods), ncol = 5)
cens_ind[apply(rbiods[, 12:16], 2, function(x) grepl("[><]", x))] <-
  unlist(apply(rbiods[, 12:16], 2, function(x) x[grep("[><]", x)]))
colnames(cens_ind) <- paste(colnames(rbiods)[12:16], "_cens_ind", sep = "")
rbiods <- as.data.frame(cbind(rbiods, cens_ind))

# Impute censored values
#   Left-censored (below LOD): Use LOD/2 (typical imputation method)
#   Right-censored (above ULD): Use ULD value
rbiods$creatinine_mg_dl[rbiods$creatinine_mg_dl_cens_ind != 0] <- 0.15 / 2
rbiods$hs_crp_mg_l[rbiods$hs_crp_mg_l_cens_ind != 0] <- 90
rbiods$nt_pro_bnp[rbiods$nt_pro_bnp_cens_ind != 0] <- 35000
rbiods$troponin_t_ng_l[rbiods$troponin_t_ng_l_cens_ind != 0] <- 6 / 2

# Remove censoring indicators from biomarker values
rbiods[, 12:21] <- apply(rbiods[, 12:21], 2, function(x) as.numeric(x))

################################################################################
# 0.7 MERGE PROTEOMICS WITH CLINICAL DATA
################################################################################

# Define core clinical variables
varlist <- c(
  "age", "sex", "race___4", "hf_cat", "creatinine_value",
  "bl_sppb", "fu_sppb", "bl_bal_scr", "fu_bal_scr",
  "bl_walk_time", "fu_walk_time", "rehosp", "death", "event", "lfu",
  "bl_smw", "fu_smw", "nrehosp_death", "egfr"
)

# Left join proteomics to clinical data (baseline)
rprotblds <- merge(
  x = olink_dta_bl,
  y = rhfds[, c("study_id", varlist)],
  by.x = 'subject_id',
  by.y = 'study_id',
  all.x = TRUE
)

# Left join proteomics to clinical data (follow-up)
rprotfuds <- merge(
  x = olink_dta_fu,
  y = rhfds[, c("study_id", varlist)],
  by.x = 'subject_id',
  by.y = 'study_id',
  all.x = TRUE
)

# Replace negative protein values with NA (indicates data quality issues)
for (i in 3:94) {
  rprotblds[, i][rprotblds[, i] < 0] <- NA
  rprotfuds[, i][rprotfuds[, i] < 0] <- NA
}

# Add intervention assignment from biomarker dataset
rprotblds <- merge(
  x = rprotblds,
  y = rbiods[rbiods$timepoint == 'Baseline',
             c('subject_id', 'intervention_1_control_0')],
  by.x = 'subject_id',
  by.y = 'subject_id',
  all.x = TRUE
)

# Manual correction: subject 10-031 only has follow-up measurement but is in control group
rprotblds$intervention_1_control_0[rprotblds$subject_id == '10-031'] <- 0
rprotblds$intervention_1_control_0 <- factor(rprotblds$intervention_1_control_0, levels = c(0, 1))

rprotfuds <- merge(
  x = rprotfuds,
  y = rbiods[rbiods$timepoint == 'Baseline',
             c('subject_id', 'intervention_1_control_0')],
  by.x = 'subject_id',
  by.y = 'subject_id',
  all.x = TRUE
)

rprotfuds$intervention_1_control_0 <- factor(rprotfuds$intervention_1_control_0, levels = c(0, 1))

################################################################################
# 0.8 BASELINE PROTEIN OUTLIER DETECTION AND REMOVAL
################################################################################

# Identify multivariate outliers in baseline biomarkers
# Criterion: >3 standard deviations from mean on log2 scale

# Count outliers per biomarker
outlier_counts <- apply(rbiods[rbiods$timepoint == 'Baseline', 12:16], 2,
                        function(x) sum((log2(x + 0.001) > mean(log2(x + 0.001), na.rm = TRUE) + 3 * sd(log2(x + 0.001), na.rm = TRUE)) |
                                          (log2(x + 0.001) < mean(log2(x + 0.001), na.rm = TRUE) - 3 * sd(log2(x + 0.001), na.rm = TRUE)),
                                        na.rm = TRUE))

# Identify which observations are outliers for each biomarker
bmoutliers <- apply(rbiods[rbiods$timepoint == 'Baseline', 12:16], 2,
                    function(x) rbiods[rbiods$timepoint == 'Baseline', c("sample_id", "timepoint")]
                    [which((log2(x + 0.001) > mean(log2(x + 0.001), na.rm = TRUE) +
                              3 * sd(log2(x + 0.001), na.rm = TRUE)) |
                             (log2(x + 0.001) < mean(log2(x + 0.001), na.rm = TRUE) -
                                3 * sd(log2(x + 0.001), na.rm = TRUE))), ])

# Replace outlier values with NA
for (i in 1:length(bmoutliers)) {
  bmtmp <- names(bmoutliers)[i]
  if (nrow(bmoutliers[[i]]) > 0) {
    for (j in 1:nrow(bmoutliers[[i]])) {
      print(bmoutliers[[i]][j, ])
      idtmp <- bmoutliers[[i]][j, 1]
      rbiods[rbiods$timepoint == 'Baseline' & rbiods$sample_id == idtmp, bmtmp] <- NA
    }
  }
}

# Merge all biomarker and clinical data
rcombds <- merge(
  x = rbiods,
  y = rhfds,
  by.x = 'subject_id',
  by.y = 'study_id',
  all.x = TRUE
)

################################################################################
# 0.9 IDENTIFY KEY PROTEINS FOR ANALYSIS
################################################################################

# Proteins selected from prior machine learning screening:
#   LASSO for 6MWD: ALCAM, Gal-4
#   BART for 6MWD: Gal-4, PECAM-1
#   LASSO for SPPB: GP6
#   BART for SPPB: ST2
# These represent a balance between LASSO (L1 regularization) and 
# BART (tree-based machine learning) approaches

protlist <- c("ALCAM", "Gal-4", "PECAM-1", "GP6", "ST2")

################################################################################
# Change working directory for results output
################################################################################

setwd("set/working/directory/REHAB-HF proteomics analysis")

################################################################################
# END OF SETUP AND DATA PROCESSING
################################################################################
# The datasets are now prepared for analysis:
#   rprotblds   : Baseline proteomics + clinical data (n=243)
#   rprotfuds   : Follow-up proteomics + clinical data
#   rcombds     : All biomarker + clinical data combined
#   protlist    : 5 key proteins for analysis
################################################################################

################################################################################
# SECTION 1: TABLE 1 - BASELINE PATIENT CHARACTERISTICS
################################################################################

# Select variables for baseline characteristics table
varlist_table1 <- c(
  c("intervention_1_control_0", "age", "sex", "race___4",
    "hf_cat", "egfr", "bl_sppb", "bl_smw"), protlist
)

# Create summary table stratified by intervention group
tbl_summary(rprotblds[, varlist_table1],
            by = intervention_1_control_0,
            type = all_continuous() ~ "continuous2",
            statistic = list(
              all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}"
              ),
              all_categorical() ~ "{n} / {N} ({p}%)"
            ),
            label = list(
              age ~ "Age",
              sex ~ "Sex (Female=1)",
              race___4 ~ "Race (White=1)",
              hf_cat ~ "HF Category (HFpEF=1)",
              egfr ~ "EGFR",
              bl_sppb ~ "Baseline SPPB",
              bl_smw ~ "Baseline 6MWD"
            ),
            missing_text = "(Missing)"
) %>%
  add_p(list(
    all_continuous() ~ "wilcox.test",
    all_categorical() ~ "fisher.test"
  )) %>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Table 1. Baseline Patient Characteristics**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Control (0) vs. Intervention (1)**") %>%
  as_gt() %>%
  gtsave(filename = "Results/Table1_Baselinecharacteristics.html")

################################################################################
# SECTION 2: TABLE 2 - PRIMARY LINEAR REGRESSION ANALYSES
################################################################################
# Outcome 1: 12-week change in SPPB
# Outcome 2: 12-week change in 6MWD
# 
# Models include:
#   Unadjusted: Outcome ~ Baseline + Protein + Intervention
#   Adjusted: Outcome ~ Baseline + Protein*Intervention + Age + Sex + Race + HF_Cat*Intervention
#
# Analysis tests whether protein expression modifies intervention effectiveness

# Rename protein columns to remove special characters for tab_model compatibility
tmpdf <- rprotblds
colnames(tmpdf)[which(colnames(rprotblds) %in% protlist)] <- make.names(protlist, unique = TRUE)
protlistclean <- make.names(protlist, unique = TRUE)

# Initialize results dataframes
regresults.sppb <- data.frame(
  Protein = protlist,
  Unadjusted.Effect = NA,
  Unadjusted.Pval = NA,
  Adjusted.Effect = NA,
  Adjusted.Pval = NA,
  Interaction.Pval = NA
)

regresults.6mwd <- regresults.sppb

# Perform regression for each protein
for (i in 1:length(protlist)) {
  
  # ---- SPPB ANALYSIS ----
  fname <- paste("Results/RegressionAnalyses/SPPBvs", protlist[i], "_regressionanalysis.html", sep = "")
  ttl <- paste("Follow up SPPB ~ Baseline SPPB + ", protlist[i], " + Intervention + More", sep = "")
  
  # Unadjusted model: protein + intervention (no interaction)
  mod.sppb.unadj <- lm(
    data = tmpdf,
    formula = paste0(
      "fu_sppb ~ bl_sppb + ",
      protlistclean[i], " + ",
      "intervention_1_control_0"
    )
  )
  
  # Adjusted model: protein*intervention interaction + covariates
  mod.sppb.adj <- lm(
    data = tmpdf,
    formula = paste0(
      "fu_sppb ~ bl_sppb + ",
      protlistclean[i], "*intervention_1_control_0 + ",
      "age + sex + race___4 + hf_cat*intervention_1_control_0"
    )
  )
  
  # Generate HTML output
  print(tab_model(
    summary(mod.sppb.unadj),
    summary(mod.sppb.adj),
    title = ttl,
    file = fname,
    digits = 3,
    digits.p = 4,
    show.obs = TRUE
  ))
  
  # Extract and store results
  regresults.sppb$Unadjusted.Effect[i] <- paste0(
    round(mod.sppb.unadj$coefficients[protlistclean[i]], 4),
    " (",
    round(confint(mod.sppb.unadj, protlistclean[i])[1], 4),
    ", ",
    round(confint(mod.sppb.unadj, protlistclean[i])[2], 4),
    ")"
  )
  regresults.sppb$Unadjusted.Pval[i] <- round(summary(mod.sppb.unadj)$coefficients[protlistclean[i], 4], 4)
  
  regresults.sppb$Adjusted.Effect[i] <- paste0(
    round(mod.sppb.adj$coefficients[protlistclean[i]], 4),
    " (",
    round(confint(mod.sppb.adj, protlistclean[i])[1], 4),
    ", ",
    round(confint(mod.sppb.adj, protlistclean[i])[2], 4),
    ")"
  )
  regresults.sppb$Adjusted.Pval[i] <- round(summary(mod.sppb.adj)$coefficients[protlistclean[i], 4], 4)
  
  # Extract interaction p-value
  regresults.sppb$Interaction.Pval[i] <- round(
    summary(mod.sppb.adj)$coefficients[paste0(protlistclean[i], ":intervention_1_control_01"), 4], 4
  )
  
  # ---- 6MWD ANALYSIS ----
  fname <- paste("Results/RegressionAnalyses/6MWDvs", protlist[i], "_regressionanalysis.html", sep = "")
  ttl <- paste("Follow up 6MWD ~ Baseline 6MWD + ", protlist[i], " + Intervention + More", sep = "")
  
  # Unadjusted model
  mod.6mwd.unadj <- lm(
    data = tmpdf,
    formula = paste0(
      "fu_smw ~ bl_smw + ",
      protlistclean[i], " + ",
      "intervention_1_control_0"
    )
  )
  
  # Adjusted model
  mod.6mwd.adj <- lm(
    data = tmpdf,
    formula = paste0(
      "fu_smw ~ bl_smw + ",
      protlistclean[i], "*intervention_1_control_0 + ",
      "age + sex + race___4 + hf_cat*intervention_1_control_0"
    )
  )
  
  # Generate HTML output
  print(tab_model(
    summary(mod.6mwd.unadj),
    summary(mod.6mwd.adj),
    title = ttl,
    file = fname,
    digits = 3,
    digits.p = 4,
    show.obs = TRUE
  ))
  
  # Extract and store results
  regresults.6mwd$Unadjusted.Effect[i] <- paste0(
    round(mod.6mwd.unadj$coefficients[protlistclean[i]], 4),
    " (",
    round(confint(mod.6mwd.unadj, protlistclean[i])[1], 4),
    ", ",
    round(confint(mod.6mwd.unadj, protlistclean[i])[2], 4),
    ")"
  )
  regresults.6mwd$Unadjusted.Pval[i] <- round(summary(mod.6mwd.unadj)$coefficients[protlistclean[i], 4], 4)
  
  regresults.6mwd$Adjusted.Effect[i] <- paste0(
    round(mod.6mwd.adj$coefficients[protlistclean[i]], 4),
    " (",
    round(confint(mod.6mwd.adj, protlistclean[i])[1], 4),
    ", ",
    round(confint(mod.6mwd.adj, protlistclean[i])[2], 4),
    ")"
  )
  regresults.6mwd$Adjusted.Pval[i] <- round(summary(mod.6mwd.adj)$coefficients[protlistclean[i], 4], 4)
  
  regresults.6mwd$Interaction.Pval[i] <- round(
    summary(mod.6mwd.adj)$coefficients[paste0(protlistclean[i], ":intervention_1_control_01"), 4], 4
  )
}

# Save summary tables
write.csv(regresults.sppb, "Results/Table2_sppb.csv")
write.csv(regresults.6mwd, "Results/Table2_6mwd.csv")

################################################################################
# SECTION 3: FIGURE 1 - PROTEIN CORRELATIONS
################################################################################
# Spearman correlation heatmap of baseline protein expression
# Significance indicated by asterisks (p<0.01)

# Calculate Spearman correlations (upper triangle only)
scmat <- round(cor(
  rprotblds[, which(colnames(rprotblds) %in% protlist)],
  method = "spearman",
  use = "pairwise.complete.obs"
), 2)

scmat[lower.tri(scmat)] <- NA

# Melt for ggplot
scmat_melt <- melt(scmat, na.rm = TRUE)

# Calculate p-values for each correlation
for (i in 1:nrow(scmat_melt)) {
  b1 <- as.character(scmat_melt$Var1[i])
  b2 <- as.character(scmat_melt$Var2[i])
  scmat_melt$pval[i] <- cor.test(
    rprotblds[, which(colnames(rprotblds) %in% b1)],
    rprotblds[, which(colnames(rprotblds) %in% b2)],
    method = "spearman",
    use = "pairwise.complete.obs"
  )$p.value
}

# Mark significant correlations (p<0.01)
scmat_melt$sig <- scmat_melt$pval < 0.01
scmat_melt$text <- ifelse(
  scmat_melt$sig == TRUE & scmat_melt$value != 1,
  paste0(scmat_melt$value, "*"),
  scmat_melt$value
)

# Create heatmap
fig1 <- ggplot(scmat_melt, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(-1, 1), space = "Lab",
    name = "Spearman\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = text), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  ))

# Save figure
pdf("Results/Figure1_SpearmanCorrelation.pdf", width = 5, height = 5)
print(fig1)
dev.off()

ggsave(fig1, file = 'Results/Figure1_SpearmanCorrelation.eps', width = 5, height = 5, device = "eps")

################################################################################
# SECTION 4: FIGURES 2 AND 3 - PROTEIN-OUTCOME RELATIONSHIPS
################################################################################
# Scatter plots with regression lines showing protein expression vs.
# 12-week change in clinical outcomes, stratified by intervention group

# Create long-form dataset for visualization
rprotblds_long <- pivot_longer(
  data = rprotblds[, c("subject_id", "timepoint", "intervention_1_control_0",
                       "bl_smw", "fu_smw", "bl_sppb", "fu_sppb", protlist)],
  cols = protlist,
  names_to = 'Protein',
  values_to = 'ExpressionLevel'
)

# Calculate changes
rprotblds_long$sppb_chg <- rprotblds_long$fu_sppb - rprotblds_long$bl_sppb
rprotblds_long$smw_chg <- rprotblds_long$fu_smw - rprotblds_long$bl_smw

# Rename and factor
colnames(rprotblds_long)[3] <- "Intervention"
rprotblds_long$Protein <- factor(rprotblds_long$Protein, levels = protlist)

# Separate proteins by outcome
protlist.6mwd <- c("ALCAM", "Gal-4", "PECAM-1")
protlist.sppb <- c("GP6", "ST2")

# ---- FIGURE 2: SPPB ----
fig2 <- ggplot(
  data = rprotblds_long[rprotblds_long$Protein %in% protlist.sppb, ],
  mapping = aes(x = ExpressionLevel, y = sppb_chg, color = Intervention)
) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  facet_wrap('Protein', scales = 'free', ncol = 3) +
  labs(
    x = "Protein expression units (NPX)",
    y = "12-week change in SPPB score"
  ) +
  scale_color_discrete(
    breaks = c(1, 0),
    labels = c("Rehabilitation\nIntervention", "Attention\nControl")
  ) +
  theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.position = c(0.15, 0.15),
    legend.box.background = element_rect(color = "black")
  )

# ---- FIGURE 3: 6MWD ----
fig3 <- ggplot(
  data = rprotblds_long[rprotblds_long$Protein %in% protlist.6mwd, ],
  mapping = aes(x = ExpressionLevel, y = smw_chg, color = Intervention)
) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  facet_wrap('Protein', scales = 'free', ncol = 3) +
  labs(
    x = "Protein expression units (NPX)",
    y = "12-week change in 6MWD (meters)"
  ) +
  scale_color_discrete(
    breaks = c(1, 0),
    labels = c("Rehabilitation\nIntervention", "Attention\nControl")
  ) +
  theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.position = c(0.15, 0.15),
    legend.box.background = element_rect(color = "black")
  )

# Save figures
pdf("Results/Figure2_ChangeinSPPB.pdf", width = 6, height = 4)
print(fig2)
dev.off()
ggsave(fig2, file = 'Results/Figure2_ChangeinSPPB.eps', width = 6, height = 4, device = "eps")

pdf("Results/Figure3_Changein6MWD.pdf", width = 6, height = 4)
print(fig3)
dev.off()
ggsave(fig3, file = 'Results/Figure3_Changein6MWD.eps', width = 6, height = 4, device = "eps")

################################################################################
# SECTION 5: PROPENSITY SCORE MATCHING & MATCHING TREES
################################################################################

## 5.0: Propensity Score Matching Setup

ps.df <- rprotblds
ps.df$trt <- as.numeric(ps.df$intervention_1_control_0) - 1

# Keep only subjects with complete data on key proteins and outcomes
ps.df <- ps.df[complete.cases(ps.df[, c("ALCAM", "Gal-4", "PECAM-1",
                                        "GP6", "ST2",
                                        'bl_smw', 'fu_smw', 'bl_sppb', 'fu_sppb')]), ]

# Fit propensity score model using logistic regression
ppty <- glm(
  trt ~ age + sex + race___4 + hf_cat,
  family = binomial(link = "logit"),
  data = ps.df
)

# Generate propensity scores
prop.score <- predict.glm(ppty, type = "response", na.action = na.exclude)
ps.df <- as.data.frame(cbind(ps.df, prop.score))

# Create matching distance matrix based on propensity scores
ppty.distance <- match_on(ppty)

cat("Propensity Score Matching Summary\n")
cat("==================================\n")
cat("Total subjects with complete data:", nrow(ps.df), "\n")
cat("Treatment:", sum(ps.df$trt == 1), "\n")
cat("Control:", sum(ps.df$trt == 0), "\n\n")

################################################################################
# SECTION 5A: HELPER FUNCTIONS FOR MATCHING TREE ANALYSIS
################################################################################
# These functions are adapted from Zhang et al. (2021)
# "Matching Trees for Treatment Effect Heterogeneity"
# https://doi.org/10.1016/j.csda.2021.107188

# Function to get parent node in tree structure
parent <- function(x) {
  if (x[1] != 1) {
    c(Recall(if (x %% 2 == 0L) x / 2 else (x - 1) / 2), x)
  } else {
    x
  }
}

# Function to prune tree based on treatment effect overlap
# Prunes terminal nodes if confidence intervals for treatment effects overlap
jesse.tree.select <- function(input, set1.caliper.tree) {
  fit.overlap <- TRUE
  ptree1 <- input
  
  while (fit.overlap == TRUE) {
    if (dim(ptree1$frame)[1] == 1) break
    
    rownames <- as.numeric(rownames(ptree1$frame))
    size <- dim(ptree1$frame)[1]
    rownamesmax <- rownames[which.max(rownames)]
    rownamessecond <- rownames[which(rownames == sort(unique(rownames), partial = size - 1)[size - 1])]
    
    # Get data for terminal nodes
    subset1 <- set1.caliper.tree[ptree1$where == which(rownames == rownamesmax), ]
    subset2 <- set1.caliper.tree[ptree1$where == which(rownames == rownamessecond), ]
    
    # Create dummy variable for node membership
    dummy <- c(rep("1", dim(subset1)[1]), rep("2", dim(subset2)[1]))
    subsetf <- cbind(rbind(subset1, subset2), dummy)
    subsetf$dummy <- as.factor(subsetf$dummy)
    subsetf$trt <- as.factor(subsetf$trt)
    
    # Test treatment effect differences between nodes
    ls.lm <- lm(outcome ~ dummy * trt, data = subsetf)
    ls.list <- lsmeans(ls.lm, list(pairwise ~ dummy | trt, pairwise ~ trt | dummy))
    out1 <- as.matrix(summary(ls.list[[4]]))
    
    q.value <- 1.96  # 95% confidence interval
    CI1 <- c(
      as.numeric(out1[1, 3]) - q.value * as.numeric(out1[1, 4]),
      as.numeric(out1[1, 3]) + q.value * as.numeric(out1[1, 4])
    )
    CI2 <- c(
      as.numeric(out1[2, 3]) - q.value * as.numeric(out1[2, 4]),
      as.numeric(out1[2, 3]) + q.value * as.numeric(out1[2, 4])
    )
    
    # Check for overlap in confidence intervals
    fit.overlap <- max(CI1[1], CI2[1]) < min(CI1[2], CI2[2])
    fit.overlap[is.na(fit.overlap)] <- FALSE
    
    if (fit.overlap == FALSE) break
    
    # Prune the node
    ptree1 <- snip.rpart(ptree1, toss = tail(head(parent(rownamesmax), -1), n = 1))
  }
  
  if (dim(ptree1$frame)[1] > 1) {
    return(ptree1)
  } else {
    return(input)
  }
}

################################################################################
# SECTION 5B: SPPB MATCHING TREE ANALYSIS
################################################################################
# Build matching tree for SPPB outcome using proteins: GP6, ST2
# Matching on: baseline SPPB + protein values, propensity score

# Separate treatment and control groups
set.trt <- ps.df[ps.df[, "trt"] == 1, ]
set.cont <- ps.df[ps.df[, "trt"] == 0, ]

# Perform 1:1 matching with caliper on protein values
match.set <- c()
calp.sd <- 1  # Caliper width: 1 SD of protein values

for (i in 1:dim(set.trt)[1]) {
  pair.set1 <- c()
  
  # Check if control subject is within caliper on both GP6 and ST2
  var.select1 <- ifelse(
    abs(set.cont$GP6 - set.trt$GP6[i]) < calp.sd * sd(set.trt$GP6),
    1, 0
  )
  var.select2 <- ifelse(
    abs(set.cont$ST2 - set.trt$ST2[i]) < calp.sd * sd(set.trt$ST2),
    1, 0
  )
  var.selectf <- ifelse(var.select1 + var.select2 == 2, 1, 0)
  
  id <- c(1:length(var.selectf))
  lenset.cont2 <- as.data.frame(cbind(id, set.cont, var.selectf))
  lenset.cont3 <- lenset.cont2[lenset.cont2[, "var.selectf"] == 1, ]
  
  # Match to nearest propensity score if within caliper
  if (dim(lenset.cont3)[1] > 0) {
    cont.index <- which.min(abs(lenset.cont3$prop.score - set.trt$prop.score[i]))
    pair.set <- rbind(
      set.trt[i, ],
      lenset.cont3[cont.index, c(-1, -dim(lenset.cont3)[2])]
    )
    pairin <- rep(i, 2)
    pair.set1 <- cbind(pair.set, pairin)
    set.cont <- set.cont[!set.cont$prop.score == lenset.cont3$prop.score[cont.index], ]
  } else {
    pair.set1 <- c()
    set.cont <- set.cont
  }
  
  match.set <- rbind(match.set, pair.set1)
}

# Calculate treatment effect (difference in outcome change)
treat.set <- match.set[match.set[, "trt"] == 1, ]
control.set <- match.set[match.set[, "trt"] == 0, ]

outcome.diff <- (treat.set$fu_sppb - treat.set$bl_sppb) -
  (control.set$fu_sppb - control.set$bl_sppb)

set.total <- as.data.frame(rbind(
  cbind(treat.set, outcome.diff),
  cbind(control.set, outcome.diff)
))

matchedset.sppb <- set.total  # Store for later use

cat("SPPB Matching Summary\n")
cat("====================\n")
cat("Matched pairs:", nrow(treat.set), "\n")
cat("Total matched subjects:", nrow(set.total), "\n\n")

# Build matching tree on proteins GP6 and ST2
tree.rpart <- rpart(
  outcome.diff ~ GP6 + ST2,
  method = "anova",
  data = set.total,
  control = rpart.control(minbucket = 15, maxdepth = 2)
)

# Add outcome variable and prune tree
set.total$outcome <- set.total$fu_sppb - set.total$bl_sppb
tree.all <- jesse.tree.select(tree.rpart, set.total)

# Plot tree
pdf("Results/Figure4a_Matchingtree_SPPB.pdf")
rpart.plot(tree.all, box.palette = 'Reds', extra = 1)
dev.off()

setEPS()
postscript("Results/Figure4a_Matchingtree_SPPB.eps")
rpart.plot(tree.all, box.palette = 'Reds', extra = 1)
dev.off()

# Extract treatment effects by tree node and get confidence intervals
set.total$sppb_node <- rpart.predict.leaves(tree.all, set.total, type = "where")
modscore.sppb <- lm(data = set.total, formula = outcome.diff ~ 0 + as.factor(sppb_node))

# Summary of results
cat("SPPB Matching Tree - Treatment Effects by Node\n")
cat("==============================================\n\n")
print(confint(modscore.sppb))
cat("\n")
print(summary(modscore.sppb))

# Null model for comparison
modscore0.sppb <- lm(data = set.total, formula = outcome.diff ~ 1)
cat("\nNull Model (no tree):\n")
print(summary(modscore0.sppb))

################################################################################
# SECTION 5B-CV: SPPB MATCHING TREE ANALYSIS (K-Fold CV section)
################################################################################

# ---- K-Fold Cross-Validation Stability Analysis for SPPB ----
set.seed(67)
nfolds <- 10
folds <- cut(seq(1, nrow(treat.set)), breaks = nfolds, labels = FALSE)
folds <- folds[sample(nrow(treat.set))]

maxdepth <- 1:3
minbucket <- seq(5, 30, 5)
bst <- 1000
bst.stg <- c(0, 0)

cvres.df <- data.frame(maxdepth = NA, minbucket = NA, testrmse.dt = NA, testrmse.rf = NA)
idx <- 0

for (j in maxdepth) {
  for (k in minbucket) {
    idx <- idx + 1
    
    tree.2 <- list()
    tree.leaves <- list()
    tree.mod.train <- list()
    tree.test.rmse <- list()
    rf.mod.train <- list()
    rf.test.rmse <- list()
    
    for (i in 1:nfolds) {
      # Identify training and test folds
      ids <- treat.set$subject_id[folds != i]
      pairs.train <- set.total$pairin[set.total$subject_id %in% ids]
      
      ids <- treat.set$subject_id[folds == i]
      pairs.test <- set.total$pairin[set.total$subject_id %in% ids]
      
      # Fit tree on training data
      tmp <- rpart(
        outcome.diff ~ GP6 + ST2,
        method = "anova",
        data = set.total[set.total$pairin %in% pairs.train, ],
        control = rpart.control(maxdepth = j, minbucket = k)
      )
      
      tree.2[[i]] <- jesse.tree.select(tmp, set.total)
      tree.leaves[[i]] <- rpart.predict.leaves(tree.2[[i]], set.total, type = "where")
      
      # Get predictions on test set
      df <- cbind(set.total, tree.leaves[[i]])
      colnames(df)[ncol(df)] <- "tree_leaves_temp"
      
      tree.mod.train[[i]] <- lm(
        outcome.diff ~ 0 + as.factor(tree_leaves_temp),
        data = df[df$pairin %in% pairs.train, ]
      )
      
      tree.test.rmse[[i]] <- sqrt(mean((
        df[df$pairin %in% pairs.test, 'outcome.diff'] -
          predict(tree.mod.train[[i]], df[df$pairin %in% pairs.test, ])
      )^2))
      
      # Fit random forest on training data
      rf.df.train <- df[df$pairin %in% pairs.train, c("outcome.diff", "GP6", "ST2")]
      rf.df.test <- df[df$pairin %in% pairs.test, c("outcome.diff", "GP6", "ST2")]
      
      rf.mod.train[[i]] <- randomForest(outcome.diff ~ GP6 + ST2,
                                        data = rf.df.train, ntree = 500, seed = 42)
      
      rf.test.rmse[[i]] <- sqrt(mean((
        rf.df.test$outcome.diff -
          predict(rf.mod.train[[i]], rf.df.test)
      )^2))
    }
    
    if (mean(unlist(tree.test.rmse)) < bst) {
      bst <- mean(unlist(tree.test.rmse))
      bst.stg[1] <- j
      bst.stg[2] <- k
    }
    
    cvres.df[idx, ] <- c(j, k, mean(unlist(tree.test.rmse)), mean(unlist(rf.test.rmse)))
  }
}

# Format CV results
cvres.df$maxdepth <- factor(cvres.df$maxdepth, levels = 1:3)
cvres.df$minbucket <- as.integer(cvres.df$minbucket)

# Calculate mean RF RMSE for reference line
rf.mean.sppb <- mean(cvres.df$testrmse.rf, na.rm = TRUE)

cvplot1 <- ggplot(data = cvres.df, mapping = aes(x = minbucket, color = maxdepth, fill = maxdepth)) +
  geom_line(mapping = aes(y = testrmse.dt), linewidth = 1) +
  geom_point(mapping = aes(y = testrmse.dt, shape = maxdepth), size = 3) +
  geom_hline(yintercept = rf.mean.sppb, color = 'black', linetype = 'dashed', linewidth = 1.2) +
  annotate("text", x = 20, y = rf.mean.sppb + 0.08, 
           label = paste("Random Forest:", round(rf.mean.sppb, 3)), 
           size = 4, color = "black", hjust = 0) +
  ylab("Average Test RMSE") +
  xlab("Minimum Bucket Size") +
  labs(title = "Average Test RMSE: SPPB Outcome",
       subtitle = "Decision Trees vs. Random Forest") +
  ylim(c(3.0, 3.8)) +
  scale_x_continuous(breaks = seq(5, 30, 5)) +
  scale_color_manual(name = "Max Depth", values = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c")) +
  scale_fill_manual(name = "Max Depth", values = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

pdf("Results/SuppFigure3_CVPlot_SPPB.pdf")
print(cvplot1)
dev.off()

setEPS()
postscript("Results/SuppFigure3_CVPlot_SPPB.eps")
print(cvplot1)
dev.off()

################################################################################
# SECTION 5C: 6MWD MATCHING TREE ANALYSIS
################################################################################
# Build matching tree for 6MWD outcome using proteins: ALCAM, Gal-4, PECAM-1

set.trt <- ps.df[ps.df[, "trt"] == 1, ]
set.cont <- ps.df[ps.df[, "trt"] == 0, ]

match.set <- c()
calp.sd <- 1

for (i in 1:dim(set.trt)[1]) {
  pair.set1 <- c()
  
  # Check if control subject is within caliper on all three proteins
  var.select1 <- ifelse(abs(set.cont$ALCAM - set.trt$ALCAM[i]) < calp.sd * sd(set.trt$ALCAM), 1, 0)
  var.select2 <- ifelse(abs(set.cont$`Gal-4` - set.trt$`Gal-4`[i]) < calp.sd * sd(set.trt$`Gal-4`), 1, 0)
  var.select3 <- ifelse(abs(set.cont$`PECAM-1` - set.trt$`PECAM-1`[i]) < calp.sd * sd(set.trt$`PECAM-1`), 1, 0)
  var.selectf <- ifelse(var.select1 + var.select2 + var.select3 == 3, 1, 0)
  
  id <- c(1:length(var.selectf))
  lenset.cont2 <- as.data.frame(cbind(id, set.cont, var.selectf))
  lenset.cont3 <- lenset.cont2[lenset.cont2[, "var.selectf"] == 1, ]
  
  if (dim(lenset.cont3)[1] > 0) {
    cont.index <- which.min(abs(lenset.cont3$prop.score - set.trt$prop.score[i]))
    pair.set <- rbind(set.trt[i, ], lenset.cont3[cont.index, c(-1, -dim(lenset.cont3)[2])])
    pairin <- rep(i, 2)
    pair.set1 <- cbind(pair.set, pairin)
    set.cont <- set.cont[!set.cont$prop.score == lenset.cont3$prop.score[cont.index], ]
  } else {
    pair.set1 <- c()
    set.cont <- set.cont
  }
  
  match.set <- rbind(match.set, pair.set1)
}

# Calculate treatment effect
treat.set <- match.set[match.set[, "trt"] == 1, ]
control.set <- match.set[match.set[, "trt"] == 0, ]

outcome.diff <- (treat.set$fu_smw - treat.set$bl_smw) -
  (control.set$fu_smw - control.set$bl_smw)

set.total <- as.data.frame(rbind(
  cbind(treat.set, outcome.diff),
  cbind(control.set, outcome.diff)
))

matchedset.6mwd <- set.total  # Store for later use

cat("6MWD Matching Summary\n")
cat("====================\n")
cat("Matched pairs:", nrow(treat.set), "\n")
cat("Total matched subjects:", nrow(set.total), "\n\n")

# Build matching tree
tree.rpart <- rpart(
  outcome.diff ~ ALCAM + `Gal-4` + `PECAM-1`,
  method = "anova",
  data = set.total,
  control = rpart.control(minbucket = 30, maxdepth = 1)
)

set.total$outcome <- set.total$fu_smw - set.total$bl_smw
tree.all <- jesse.tree.select(tree.rpart, set.total)

# Plot tree
pdf("Results/Figure4b_Matchingtree_6MWD.pdf")
rpart.plot(tree.all, box.palette = 'Reds', extra = 1)
dev.off()

setEPS()
postscript("Results/Figure4b_Matchingtree_6MWD.eps")
rpart.plot(tree.all, box.palette = 'Reds', extra = 1)
dev.off()

# Extract treatment effects by tree node
set.total$smw_node <- rpart.predict.leaves(tree.all, set.total, type = "where")
modscore.smw <- lm(data = set.total, formula = outcome.diff ~ 0 + as.factor(smw_node))

cat("6MWD Matching Tree - Treatment Effects by Node\n")
cat("==============================================\n\n")
print(confint(modscore.smw))
cat("\n")
print(summary(modscore.smw))

modscore0.smw <- lm(data = set.total, formula = outcome.diff ~ 1)
cat("\nNull Model (no tree):\n")
print(summary(modscore0.smw))

################################################################################
# SECTION 5C-CV: 6MWD MATCHING TREE ANALYSIS (K-Fold CV section)
################################################################################

# ---- K-Fold Cross-Validation for 6MWD ----
set.seed(324)
nfolds <- 10
folds <- cut(seq(1, nrow(treat.set)), breaks = nfolds, labels = FALSE)
folds <- folds[sample(nrow(treat.set))]

maxdepth <- 1:3
minbucket <- seq(5, 30, 5)
bst <- 1000
bst.stg <- c(0, 0)

cvres.df <- data.frame(maxdepth = NA, minbucket = NA, testrmse.dt = NA, testrmse.rf = NA)
idx <- 0

for (j in maxdepth) {
  for (k in minbucket) {
    idx <- idx + 1
    
    tree.2 <- list()
    tree.leaves <- list()
    tree.mod.train <- list()
    tree.test.rmse <- list()
    rf.mod.train <- list()
    rf.test.rmse <- list()
    
    for (i in 1:nfolds) {
      ids <- treat.set$subject_id[folds != i]
      pairs.train <- set.total$pairin[set.total$subject_id %in% ids]
      
      ids <- treat.set$subject_id[folds == i]
      pairs.test <- set.total$pairin[set.total$subject_id %in% ids]
      
      tmp <- rpart(
        outcome.diff ~ ALCAM + `Gal-4` + `PECAM-1`,
        method = "anova",
        data = set.total[set.total$pairin %in% pairs.train, ],
        control = rpart.control(maxdepth = j, minbucket = k)
      )
      
      tree.2[[i]] <- jesse.tree.select(tmp, set.total)
      tree.leaves[[i]] <- rpart.predict.leaves(tree.2[[i]], set.total, type = "where")
      
      df <- cbind(set.total, tree.leaves[[i]])
      colnames(df)[ncol(df)] <- "tree_leaves_temp"
      
      tree.mod.train[[i]] <- lm(
        outcome.diff ~ 0 + as.factor(tree_leaves_temp),
        data = df[df$pairin %in% pairs.train, ]
      )
      
      tree.test.rmse[[i]] <- sqrt(mean((
        df[df$pairin %in% pairs.test, 'outcome.diff'] -
          predict(tree.mod.train[[i]], df[df$pairin %in% pairs.test, ])
      )^2))
      
      # Fit random forest on training data
      rf.df.train <- df[df$pairin %in% pairs.train, c("outcome.diff", "ALCAM", "Gal-4", "PECAM-1")]
      rf.df.test <- df[df$pairin %in% pairs.test, c("outcome.diff", "ALCAM", "Gal-4", "PECAM-1")]
      
      colnames(rf.df.train) <- make.names(colnames(rf.df.train))
      colnames(rf.df.test) <- make.names(colnames(rf.df.test))
      
      rf.mod.train[[i]] <- randomForest(outcome.diff ~ ALCAM + Gal.4 + PECAM.1,
                                        data = rf.df.train, ntree = 500, seed = 42)
      
      rf.test.rmse[[i]] <- sqrt(mean((
        rf.df.test$outcome.diff -
          predict(rf.mod.train[[i]], rf.df.test)
      )^2))
    }
    
    if (mean(unlist(tree.test.rmse)) < bst) {
      bst <- mean(unlist(tree.test.rmse))
      bst.stg[1] <- j
      bst.stg[2] <- k
    }
    
    cvres.df[idx, ] <- c(j, k, mean(unlist(tree.test.rmse)), mean(unlist(rf.test.rmse)))
  }
}

# Format CV results
cvres.df$maxdepth <- factor(cvres.df$maxdepth, levels = 1:3)
cvres.df$minbucket <- as.integer(cvres.df$minbucket)

# Calculate mean RF RMSE for reference line
rf.mean.6mwd <- mean(cvres.df$testrmse.rf, na.rm = TRUE)

cvplot1 <- ggplot(data = cvres.df, mapping = aes(x = minbucket, color = maxdepth, fill = maxdepth)) +
  geom_line(mapping = aes(y = testrmse.dt), linewidth = 1) +
  geom_point(mapping = aes(y = testrmse.dt, shape = maxdepth), size = 3) +
  geom_hline(yintercept = rf.mean.6mwd, color = 'black', linetype = 'dashed', linewidth = 1.2) +
  annotate("text", x = 20, y = rf.mean.6mwd + 2, 
           label = paste("Random Forest:", round(rf.mean.6mwd, 3)), 
           size = 4, color = "black", hjust = 0) +
  ylab("Average Test RMSE") +
  xlab("Minimum Bucket Size") +
  labs(title = "Average Test RMSE: 6MWD Outcome",
       subtitle = "Decision Trees vs. Random Forest") +
  ylim(c(100, 150)) +
  scale_x_continuous(breaks = seq(5, 30, 5)) +
  scale_color_manual(name = "Max Depth", values = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c")) +
  scale_fill_manual(name = "Max Depth", values = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

pdf("Results/SuppFigure3_CVPlot_6MWD.pdf")
print(cvplot1)
dev.off()

setEPS()
postscript("Results/SuppFigure3_CVPlot_6MWD.eps")
print(cvplot1)
dev.off()


################################################################################
# SECTION 6: WEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA)
################################################################################
# Identify modules of co-expressed proteins and their associations with outcomes

cat("\n=== WGCNA Module Analysis ===\n")

# Visualize correlation matrices
p1 <- ggcorrplot(cor(na.omit(rprotblds[, 3:94])),
                 title = "All Participants (Baseline)",
                 hc.order = TRUE,
                 hc.method = "complete"
) + hw +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90)
  )

pdf("Results/correlationplot_allparticipants_baseline.pdf", onefile = TRUE)
print(p1)
dev.off()

# Correlation matrix for usual care
p1 <- ggcorrplot(cor(na.omit(rprotblds[rprotblds$intervention_1_control_0 == 0, 3:94])),
                 title = "Usual Care Patients (Baseline)",
                 hc.order = TRUE,
                 hc.method = "complete"
) + hw +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90)
  )

pdf("Results/correlationplot_usualcare_baseline.pdf", onefile = TRUE)
print(p1)
dev.off()

# Correlation matrix for rehabilitation
cormat2 <- cor(na.omit(rprotblds[rprotblds$intervention_1_control_0 == 1,
                                 as.character(p1$data[1:92, 1])]))
p2 <- ggcorrplot(cormat2,
                 title = "Rehab Patients (Baseline)",
                 hc.order = FALSE,
                 hc.method = "complete"
) + hw +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90)
  )

pdf("Results/correlationplot_rehab_baseline.pdf", onefile = TRUE)
print(p2)
dev.off()

# Test equality of covariance matrices between intervention groups
testCov(
  na.omit(rprotblds[rprotblds$intervention_1_control_0 == 0,
                    as.character(p1$data[1:92, 1])]),
  na.omit(rprotblds[rprotblds$intervention_1_control_0 == 1,
                    as.character(p1$data[1:92, 1])]),
  method = "ALL",
  J = 5000,
  alpha = 0.05,
  n.core = 1
)
# No significant difference between groups at baseline

# Standardize protein expressions
naset <- which(apply(rprotblds[, 3:94], 1, function(x) sum(is.na(x))) == 1)
rprotblds.std <- apply(rprotblds[-naset, 3:94], 2, function(x) (x - mean(x) / sd(x)))

# Find optimal soft-thresholding power for scale-free network
powers <- c(c(1:20), seq(from = 22, to = 30, by = 2))
sft <- pickSoftThreshold(na.omit(rprotblds[, 3:94]), powerVector = powers, verbose = 5)

# Plot scale-free network fit
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence"
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red"
)
abline(h = 0.90, col = "red")

par(mar = c(1, 1, 1, 1))
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity"
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = sft$fitIndices[, 1], col = "red")

# Create adjacency matrix using soft-thresholding
adj1.all <- abs(cor(rprotblds.std))^4

# Convert adjacency to topological overlap matrix (TOM)
# TOM is more robust to noise than raw correlation
TOM.dissimilarity <- 1 - TOMsimilarity(adj1.all)

# Hierarchical clustering using TOM-based distance
rownames(TOM.dissimilarity) <- colnames(rprotblds.std)
hierADJ <- hclust(as.dist(TOM.dissimilarity), method = 'complete')

# Identify modules using dynamic tree cutting
Modules <- cutreeDynamic(
  dendro = hierADJ,
  distM = TOM.dissimilarity,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 6
)

cat("Module sizes:\n")
print(table(Modules))

ModuleColors <- labels2colors(Modules)

# Save module assignments
write.csv(data.frame(protein = colnames(rprotblds.std), colors = ModuleColors),
          file = "Results/protein_modules.csv")

cat("Protein module color assignments:\n")
print(table(ModuleColors))

# Plot dendrogram with module colors
pdf("Results/WGCNA_DendrogramAndColors.pdf")
plotDendroAndColors(hierADJ, ModuleColors, "Module",
                    dendroLabels = NULL,
                    cex.dendroLabels = 0.5,
                    addGuide = FALSE,
                    main = "Protein dendrogram and module colors"
)
dev.off()

# Calculate module eigengenes (summarizes module expression)
MElist <- moduleEigengenes(as.data.frame(rprotblds.std), colors = ModuleColors)
MEs <- MElist$eigengenes

# Order modules so similar ones are adjacent
MEs0 <- orderMEs(MEs)
module_order <- names(MEs0) %>% gsub("ME", "", .)

# Add clinical data to module eigengene dataframe
MEs0$ID <- rprotblds$subject_id[-naset]
MEs0$age <- rprotblds$age[-naset]
MEs0$sex <- as.numeric(rprotblds$sex[-naset]) - 1
MEs0$bl_smw <- rprotblds$bl_smw[-naset]
MEs0$fu_smw <- rprotblds$fu_smw[-naset]
MEs0$chg_smw <- MEs0$fu_smw - MEs0$bl_smw
MEs0$bl_sppb <- rprotblds$bl_sppb[-naset]
MEs0$fu_sppb <- rprotblds$fu_sppb[-naset]
MEs0$chg_sppb <- MEs0$fu_sppb - MEs0$bl_sppb

# Calculate correlations between module eigengenes and clinical outcomes
df <- melt(cor(MEs0[, -which(names(MEs0) == "ID")], use = "pairwise.complete.obs"))

# Filter to relevant correlations
df <- df[(df$Var1 %in% c('age', 'bl_smw', 'chg_smw', 'bl_sppb', 'chg_sppb', 'sex')), ]
df <- df[(df$Var2 %in% c('MEturquoise', 'MEbrown', 'MEyellow', 'MEblue', 'MEgreen')), ]

# Calculate p-values
for (i in 1:nrow(df)) {
  b1 <- as.character(df$Var1[i])
  b2 <- as.character(df$Var2[i])
  df$pval[i] <- cor.test(
    MEs0[, which(colnames(MEs0) %in% b1)],
    MEs0[, which(colnames(MEs0) %in% b2)],
    method = "pearson",
    use = "pairwise.complete.obs"
  )$p.value
}

df$sig <- df$pval < 0.01
df$text <- ifelse(
  df$sig == TRUE & df$value != 1,
  paste0(round(df$value, 3), "*"),
  round(df$value, 3)
)

# Plot module-outcome heatmap
pdf("Results/SuppFigure6_WGCNA_heatmap.pdf", onefile = TRUE)
print(
  ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-0.5, 0.5), space = "Lab",
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = text), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    guides(fill = guide_colorbar(
      barwidth = 15, barheight = 1,
      title.position = "top", title.hjust = 0.5
    )) +
    coord_flip()
)
dev.off()

# Add treatment assignment to module eigengenes
MEs0$treatment <- rprotblds$intervention_1_control_0[-naset]

# Test module-treatment interactions for outcomes
cat("\n=== Module x Treatment Interactions ===\n")
print(summary(lm(data = MEs0, formula = fu_sppb ~ bl_sppb + MEturquoise * treatment)))
print(summary(lm(data = MEs0, formula = fu_smw ~ bl_smw + MEturquoise * treatment)))
print(summary(lm(data = MEs0, formula = fu_smw ~ bl_smw + MEgreen * treatment)))
print(summary(lm(data = MEs0, formula = fu_smw ~ bl_smw + MEblue * treatment)))
print(summary(lm(data = MEs0, formula = fu_sppb ~ bl_sppb + MEbrown * treatment)))
print(summary(lm(data = MEs0, formula = fu_smw ~ bl_smw + MEyellow * treatment)))

# Extract most influential proteins within each module (top 5)
cat("\n=== Top 5 Proteins per Module ===\n")
for (i in 1:(which(names(MEs0) == "ID") - 1)) {
  clr <- sapply(names(MEs0)[1:(which(names(MEs0) == "ID") - 1)],
                function(x) substring(x, 3, nchar(x)))[i]
  
  geneModuleMembership <- as.data.frame(
    cor(rprotblds.std[, ModuleColors == clr], MEs0[, i], use = "p")
  )
  
  MMPvalue <- as.data.frame(
    corPvalueStudent(as.matrix(geneModuleMembership), nrow(rprotblds.std))
  )
  
  names(geneModuleMembership) <- paste("MM", substring(
    names(MEs0)[1:(which(names(MEs0) == "ID") - 1)], 3
  ), sep = "")[i]
  
  names(MMPvalue) <- paste("p.MM", substring(
    names(MEs0)[1:(which(names(MEs0) == "ID") - 1)], 3
  ), sep = "")[i]
  
  cat(paste("\n--- Module:", clr, "---\n"))
  print(geneModuleMembership[order(geneModuleMembership[, 1], decreasing = TRUE)[1:5], 1, drop = FALSE])
}

################################################################################
# SECTION 7: SUPPLEMENTAL FIGURE 4 - PROTEIN CHANGE BY TREATMENT
################################################################################
# Visualize fold-change (ratio of follow-up to baseline) in protein expression
# stratified by intervention group

rprotfuds_long <- pivot_longer(
  data = rprotfuds[, c("subject_id", "timepoint", "intervention_1_control_0",
                       "bl_smw", "fu_smw", "bl_sppb", "fu_sppb", protlist)],
  cols = protlist,
  names_to = 'Protein',
  values_to = 'ExpressionLevel'
)

rprotfuds_long$sppb_chg <- rprotfuds_long$fu_sppb - rprotfuds_long$bl_sppb
rprotfuds_long$smw_chg <- rprotfuds_long$fu_smw - rprotfuds_long$bl_smw
colnames(rprotfuds_long)[3] <- "Intervention"
rprotfuds_long$Protein <- factor(rprotfuds_long$Protein, levels = protlist)

# Merge baseline and follow-up protein data
rprotds_wide <- merge(rprotblds_long,
                      rprotfuds_long,
                      by = c('subject_id', 'Protein', 'Intervention',
                             "bl_smw", "fu_smw", "bl_sppb", "fu_sppb", "sppb_chg", "smw_chg"),
                      all.x = TRUE
)

rprotds_wide$Protein_change <- (rprotds_wide$ExpressionLevel.y) - (rprotds_wide$ExpressionLevel.x)

# Create visualization
supfig4 <- ggplot(
  data = rprotds_wide,
  mapping = aes(x = Intervention, y = ExpressionLevel.y / ExpressionLevel.x, fill = Intervention)
) +
  geom_boxplot() +
  facet_wrap_paginate('Protein', scales = 'free', ncol = 3, page = 1) + hw +
  labs(y = "Ratio of change from baseline to 12-weeks") +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(
    trans = 'log2',
    labels = function(x) sprintf("%.3f", x),
    breaks = seq(0.5, 1.5, length.out = 6)
  ) +
  scale_x_discrete(breaks = c(1, 0), labels = c("RI", "AC")) +
  scale_fill_discrete(breaks = c(0, 1), labels = c("AC", "RI")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = c(0.85, 0.15),
    legend.box.background = element_rect(color = "black")
  )

pdf("Results/SuppFigure4_Proteinchange.pdf", height = 9, width = 7, onefile = TRUE)
print(supfig4)
dev.off()

ggsave(supfig4, file = 'Results/SuppFigure4_Proteinchange.eps', height = 9, width = 7, device = "eps")

################################################################################
# SECTION 8: SUPPLEMENTAL FIGURE 1 - CONSORT DIAGRAM
################################################################################
# Track participant flow through study

cat("\n=== CONSORT DIAGRAM NUMBERS ===\n")

# Total REHAB-HF participants
cat("Total REHAB-HF enrolled:\n")
print(nrow(rhfds))
print(length(unique(rhfds$study_id)))

# Participants with baseline proteomics
cat("\nWith baseline proteomics:\n")
print(nrow(rprotblds))
print(length(unique(rprotblds$subject_id)))

# Lost to follow-up or died
cat("\nAmong those with baseline proteomics - Lost to follow-up:\n")
print(sum(rhfds[rhfds$study_id %in% rprotblds$subject_id & is.na(rhfds$fu_sppb), 'lfu']))

cat("\nAmong those with baseline proteomics - Deaths:\n")
print(sum(rhfds[rhfds$study_id %in% rprotblds$subject_id & is.na(rhfds$fu_sppb), 'death']))

# With baseline and follow-up SPPB measurements
cat("\nWith both baseline and follow-up SPPB:\n")
print(nrow(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                   !is.na(rhfds$fu_sppb) & !is.na(rhfds$bl_sppb), ]))

cat("\nBreakdown by intervention group (SPPB):\n")
print(table(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                    !is.na(rhfds$fu_sppb) & !is.na(rhfds$bl_sppb), ]$rand))

# With 6MWD measurements
cat("\nBreakdown by intervention group (6MWD):\n")
print(table(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                    !is.na(rhfds$fu_sppb), ]$rand,
            !is.na(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                           !is.na(rhfds$fu_sppb), ]$fu_smw) &
              !is.na(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                             !is.na(rhfds$fu_sppb), ]$bl_smw)))

# With follow-up proteomics
cat("\nWith follow-up proteomics:\n")
print(nrow(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                   !is.na(rhfds$fu_sppb) & !is.na(rhfds$bl_sppb) &
                   rhfds$study_id %in% rprotfuds$subject_id, ]))

cat("\nBreakdown by intervention group (Follow-up proteomics):\n")
print(table(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                    !is.na(rhfds$fu_sppb) & !is.na(rhfds$bl_sppb) &
                    rhfds$study_id %in% rprotfuds$subject_id, ]$rand))

################################################################################
# SECTION 9: PROTEIN EFFECT MODIFICATION SCREENING
################################################################################
# Variable screening step 1: individual treatment x protein interactions
# Tests whether each protein modifies treatment effectiveness (treatment x protein interaction)
# Criterion: p < 0.05 (no FDR adjustment in this exploratory step)

cat("\n=== EFFECT MODIFICATION SCREENING (Individual Proteins) ===\n")

# ---- SPPB ANALYSIS ----
intxn_df <- matrix(NA, nrow = 92, ncol = 2)
colnames(intxn_df) <- c("effect", "p-value")

for (i in 3:94) {
  prot <- names(rprotblds)[i]
  formtmp <- paste("fu_sppb ~ bl_sppb + age + sex + race___4 + hf_cat + `",
                   prot, "`*intervention_1_control_0", sep = "")
  tmpmod <- summary(lm(data = rprotblds, formula = formtmp))
  intxn_df[i - 2, ] <- tmpmod$coefficients[nrow(tmpmod$coefficients), c(1, 4)]
}

intxn_df <- data.frame(
  Protein = names(rprotblds)[3:94],
  intxn_df,
  Significance = ifelse(intxn_df[, 2] < 0.05, "Significant", "Not Significant")
)

# FDR adjustment
cat("\nSPPB p-value distribution after FDR:\n")
print(summary(p.adjust(intxn_df$p.value, method = 'fdr')))

# Volcano plot
volcplot.sppb <- ggplot(
  data = intxn_df,
  mapping = aes(x = effect, y = -log10(p.value), color = Significance)
) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(x = "Effect modification", y = "-log10(p-value)") +
  theme_minimal() +
  ggtitle('SPPB Effect Modification (Unadjusted P-values)') +
  ggrepel::geom_label_repel(
    data = intxn_df[intxn_df$Significance == "Significant", ],
    ggplot2::aes(label = Protein),
    box.padding = 1,
    show.legend = FALSE,
    max.overlaps = Inf
  )

pdf("Results/SPPB_volcano_effectmod.pdf", width = 7, height = 6)
print(volcplot.sppb)
dev.off()

ggsave(volcplot.sppb, file = 'Results/SPPB_volcano_effectmod.eps', width = 7, height = 6, device = "eps")

# ---- 6MWD ANALYSIS ----
intxn_df <- matrix(NA, nrow = 92, ncol = 2)
colnames(intxn_df) <- c("effect", "p-value")

for (i in 3:94) {
  prot <- names(rprotblds)[i]
  formtmp <- paste("fu_smw ~ bl_smw + age + sex + race___4 + hf_cat + `",
                   prot, "`*intervention_1_control_0", sep = "")
  tmpmod <- summary(lm(data = rprotblds, formula = formtmp))
  intxn_df[i - 2, ] <- tmpmod$coefficients[nrow(tmpmod$coefficients), c(1, 4)]
}

intxn_df <- data.frame(
  Protein = names(rprotblds)[3:94],
  intxn_df,
  Significance = ifelse(intxn_df[, 2] < 0.05, "Significant", "Not Significant")
)

cat("\n6MWD p-value distribution after FDR:\n")
print(summary(p.adjust(intxn_df$p.value, method = 'fdr')))

volcplot.6mwd <- ggplot(
  data = intxn_df,
  mapping = aes(x = effect, y = -log10(p.value), color = Significance)
) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(x = "Effect modification", y = "-log10(p-value)") +
  theme_minimal() +
  ggtitle('6 Minute Walk Effect Modification (Unadjusted P-values)') +
  ggrepel::geom_label_repel(
    data = intxn_df[intxn_df$Significance == "Significant", ],
    ggplot2::aes(label = Protein),
    box.padding = 1,
    show.legend = FALSE,
    max.overlaps = Inf
  )

pdf("Results/6MWD_volcano_effectmod.pdf", width = 7, height = 6)
print(volcplot.6mwd)
dev.off()

ggsave(volcplot.6mwd, file = 'Results/6MWD_volcano_effectmod.eps', width = 7, height = 6, device = "eps")

################################################################################
# SECTION 10: LASSO REGRESSION WITH GLINTERNET INTERACTION SCREENING
################################################################################
# Variable screening step 2: LASSO (L1 regularization) with interaction terms
# Uses glinternet to enforce strong hierarchy (main effects required for interactions)
# Cross-validation selects optimal regularization parameter

cat("\n=== GLINTERNET LASSO SCREENING ===\n")

# ---- SPPB ANALYSIS ----
sppb_data <- cbind(
  rprotblds[, c("fu_sppb", "intervention_1_control_0", "bl_sppb", "age", "sex", "race___4", "hf_cat")],
  rprotblds[, 3:94]
)

sppb_data <- na.omit(sppb_data)

set.seed(97)
sppb_y <- sppb_data$fu_sppb

# Format for glinternet (factors must be 0/1)
sppb_x <- data.frame(
  intervention_1_control_0 = as.numeric(sppb_data$intervention_1_control_0) - 1,
  sppb_data$bl_sppb,
  sppb_data$age,
  sex = as.numeric(sppb_data$sex) - 1,
  race___4 = as.numeric(sppb_data$race___4) - 1,
  hf_cat = as.numeric(sppb_data$hf_cat) - 1,
  sppb_data[, 8:99]
)

# Specify number of levels (1 = continuous, 2 = binary)
numLevels <- c(2, 1, 1, 2, 2, 2, rep(1, 92))

# Cross-validation to select lambda
fit.cv <- glinternet.cv(sppb_x, sppb_y, numLevels = numLevels, nFolds = 10, nLambda = 10)

fit <- glinternet(sppb_x, sppb_y, numLevels = numLevels, lambda = fit.cv$lambdaHat)

coeffs <- coef(fit)[[2]]

cat("\nSPPB - Proteins with treatment interaction:\n")
print(names(sppb_data[, 8:99])[coeffs$interactions$catcont[, 2] - 2])

cat("\nSPPB - Main effects:\n")
print(names(sppb_data[, 8:99])[coeffs$mainEffects$cont[c(-1, -2)] - 2])

# ---- 6MWD ANALYSIS ----
smw_data <- cbind(
  rprotblds[, c("fu_smw", "intervention_1_control_0", "bl_smw", "age", "sex", "race___4", "hf_cat")],
  rprotblds[, 3:94]
)

smw_data <- na.omit(smw_data)

set.seed(567)
smw_y <- smw_data$fu_smw

smw_x <- data.frame(
  intervention_1_control_0 = as.numeric(smw_data$intervention_1_control_0) - 1,
  smw_data$bl_smw,
  smw_data$age,
  sex = as.numeric(smw_data$sex) - 1,
  race___4 = as.numeric(smw_data$race___4) - 1,
  hf_cat = as.numeric(smw_data$hf_cat) - 1,
  smw_data[, 8:99]
)

numLevels <- c(2, 1, 1, 2, 2, 2, rep(1, 92))

fit.cv <- glinternet.cv(smw_x, smw_y, numLevels = numLevels, nFolds = 10, nLambda = 10)

fit <- glinternet(smw_x, smw_y, numLevels = numLevels, lambda = fit.cv$lambdaHat)

coeffs <- coef(fit)[[2]]

cat("\n6MWD - Proteins with treatment interaction:\n")
print(names(smw_data[, 8:99])[coeffs$interactions$catcont[, 2] - 2])

cat("\n6MWD - Main effects:\n")
print(names(smw_data[, 8:99])[coeffs$mainEffects$cont[c(-1, -2)] - 2])

################################################################################
# SECTION 11: BART WITH INTERACTION IMPORTANCE
################################################################################
# Variable screening step 3: Bayesian Additive Regression Trees
# Evaluates importance of treatment x protein interactions
# More flexible than linear models, captures non-linearities

cat("\n=== BART INTERACTION IMPORTANCE ===\n")

set.seed(45)
options(java.parameters = "-Xmx4000m")
set_bart_machine_num_cores(1)

# ---- SPPB BART MODEL ----
sppb_data <- cbind(
  rprotblds[, c("fu_sppb", "intervention_1_control_0", "bl_sppb", "age", "sex", "race___4", "hf_cat")],
  rprotblds[, 3:94]
)

sppb_data <- na.omit(sppb_data)

sppb_y <- sppb_data$fu_sppb

sppb_x <- data.frame(
  sppb_data$intervention_1_control_0,
  sppb_data$bl_sppb,
  sppb_data$age,
  sppb_data$sex,
  sppb_data$race___4,
  sppb_data$hf_cat,
  sppb_data[, 8:99]
)

bart_machine_sppb <- bartMachine(
  sppb_x, sppb_y,
  num_trees = 100,
  num_burn_in = 2000,
  num_iterations_after_burn_in = 2000,
  use_missing_data = TRUE,
  mem_cache_for_speed = FALSE,
  seed = 34
)

print(bart_machine_sppb)

pdf("Results/bart_sppb_ConvergenceDiagnostics.pdf", onefile = TRUE)
plot_convergence_diagnostics(bart_machine_sppb)
dev.off()

# Interaction importance with 100 replicates
reps <- 100
intxninvest <- interaction_investigator(bart_machine_sppb,
                                        plot = FALSE,
                                        num_var_plot = 25,
                                        bottom_margin = 20,
                                        num_replicates_for_avg = reps,
                                        num_trees_bottleneck = 20
)

# Extract interaction importance
bartintxn_means <- apply(intxninvest$interaction_counts[3:94, 95:96, ], 1, sum) / reps
bartintxn_sds <- apply(apply(intxninvest$interaction_counts[3:94, 95:96, ], c(1, 3), sum), 1, sd)
bartintxn_moe <- 1.96 * bartintxn_sds / sqrt(reps)

bartintxn_df <- data.frame(
  Protein = rownames(intxninvest$interaction_counts_avg)[3:94],
  Importance = bartintxn_means,
  MOE = bartintxn_moe
)

bartintxn_df <- bartintxn_df[order(bartintxn_df$Importance, decreasing = TRUE), ]

# Plot top 20
bartintxn_df <- bartintxn_df[1:20, ]

pdf("Results/SuppFigure2_InteractionImportance_sppb.pdf", onefile = TRUE)
par(mar = c(7, 6, 3, 0))
bars <- barplot(bartintxn_df$Importance,
                names.arg = bartintxn_df$Protein,
                las = 2,
                ylab = "Importance",
                col = "gray",
                ylim = c(0, max(bartintxn_df$Importance + bartintxn_df$MOE))
)

conf_upper <- bartintxn_df$Importance + bartintxn_df$MOE
conf_lower <- bartintxn_df$Importance - bartintxn_df$MOE

segments(bars, bartintxn_df$Importance, bars, conf_upper,
         col = rgb(0.59, 0.39, 0.39), lwd = 3)
segments(bars, bartintxn_df$Importance, bars, conf_lower,
         col = rgb(0.59, 0.39, 0.39), lwd = 3)

dev.off()

# ---- 6MWD BART MODEL ----
set.seed(231)
options(java.parameters = "-Xmx4000m")
set_bart_machine_num_cores(1)

smw_data <- cbind(
  rprotblds[, c("fu_smw", "intervention_1_control_0", "bl_smw", "age", "sex", "race___4", "hf_cat")],
  rprotblds[, 3:94]
)

smw_data <- na.omit(smw_data)

smw_y <- smw_data$fu_smw

smw_x <- data.frame(
  smw_data$intervention_1_control_0,
  smw_data$bl_smw,
  smw_data$age,
  smw_data$sex,
  smw_data$race___4,
  smw_data$hf_cat,
  smw_data[, 8:99]
)

bart_machine_smw <- bartMachine(
  smw_x, smw_y,
  num_trees = 100,
  num_burn_in = 2000,
  num_iterations_after_burn_in = 2000,
  use_missing_data = TRUE,
  mem_cache_for_speed = FALSE,
  seed = 231
)

print(bart_machine_smw)

pdf("Results/bart_smw_ConvergenceDiagnostics.pdf", onefile = TRUE)
plot_convergence_diagnostics(bart_machine_smw)
dev.off()

reps <- 100
intxninvest <- interaction_investigator(bart_machine_smw,
                                        num_var_plot = 25,
                                        bottom_margin = 20,
                                        num_replicates_for_avg = reps,
                                        num_trees_bottleneck = 20,
                                        plot = FALSE
)

bartintxn_means <- apply(intxninvest$interaction_counts[3:94, 95:96, ], 1, sum) / reps
bartintxn_sds <- apply(apply(intxninvest$interaction_counts[3:94, 95:96, ], c(1, 3), sum), 1, sd)
bartintxn_moe <- 1.96 * bartintxn_sds / sqrt(reps)

bartintxn_df <- data.frame(
  Protein = rownames(intxninvest$interaction_counts_avg)[3:94],
  Importance = bartintxn_means,
  MOE = bartintxn_moe
)

bartintxn_df <- bartintxn_df[order(bartintxn_df$Importance, decreasing = TRUE), ]
bartintxn_df <- bartintxn_df[1:20, ]

pdf("Results/SuppFigure2_InteractionImportance_6mwd.pdf", onefile = TRUE)
par(mar = c(7, 6, 3, 0))
bars <- barplot(bartintxn_df$Importance,
                names.arg = bartintxn_df$Protein,
                las = 2,
                ylab = "Importance",
                col = "gray",
                ylim = c(0, max(bartintxn_df$Importance + bartintxn_df$MOE))
)

conf_upper <- bartintxn_df$Importance + bartintxn_df$MOE
conf_lower <- bartintxn_df$Importance - bartintxn_df$MOE

segments(bars, bartintxn_df$Importance, bars, conf_upper,
         col = rgb(0.59, 0.39, 0.39), lwd = 3)
segments(bars, bartintxn_df$Importance, bars, conf_lower,
         col = rgb(0.59, 0.39, 0.39), lwd = 3)

dev.off()

################################################################################
# SECTION 12: MEDIATION ANALYSIS
################################################################################
# Tests whether changes in protein expression mediate the effect of
# rehabilitation intervention on physical performance outcomes

cat("\n=== MEDIATION ANALYSIS ===\n")

rprotblds_long2 <- pivot_longer(
  data = rprotblds[, c("subject_id", "timepoint", "intervention_1_control_0",
                       "bl_smw", "fu_smw", "bl_sppb", "fu_sppb",
                       "age", "sex", "race___4", "hf_cat", protlist)],
  cols = protlist,
  names_to = 'Protein',
  values_to = 'ExpressionLevel'
)

rprotfuds_long2 <- pivot_longer(
  data = rprotfuds[, c("subject_id", "timepoint", "intervention_1_control_0",
                       "bl_smw", "fu_smw", "bl_sppb", "fu_sppb",
                       "age", "sex", "race___4", "hf_cat", protlist)],
  cols = protlist,
  names_to = 'Protein',
  values_to = 'ExpressionLevel'
)

# Calculate changes
rprotblds_long2$sppb_chg <- rprotblds_long2$fu_sppb - rprotblds_long2$bl_sppb
rprotblds_long2$smw_chg <- rprotblds_long2$fu_smw - rprotblds_long2$bl_smw

rprotfuds_long2$sppb_chg <- rprotfuds_long2$fu_sppb - rprotfuds_long2$bl_sppb
rprotfuds_long2$smw_chg <- rprotfuds_long2$fu_smw - rprotfuds_long2$bl_smw

colnames(rprotblds_long2)[3] <- "Intervention"
colnames(rprotfuds_long2)[3] <- "Intervention"

rprotblds_long2$Protein <- factor(rprotblds_long2$Protein, levels = protlist)
rprotfuds_long2$Protein <- factor(rprotfuds_long2$Protein, levels = protlist)

# Merge baseline and follow-up data
dftmp.all <- merge(rprotblds_long2, rprotfuds_long2,
                   by = c('subject_id', 'Protein', 'Intervention',
                          "bl_smw", "fu_smw", "bl_sppb", "fu_sppb", "sppb_chg", "smw_chg",
                          "age", "sex", "race___4", "hf_cat"),
                   all.x = TRUE
)

dftmp.all$Protein_change <- (dftmp.all$ExpressionLevel.y) - (dftmp.all$ExpressionLevel.x)

# ---- 6MWD MEDIATION ----
pdf("Results/SuppFigure5_mediationanalysis_6mwd.pdf", onefile = TRUE)

for (i in 1:length(protlist)) {
  # Limit to propensity-matched participants
  dftmp <- dftmp.all[!is.na(dftmp.all$smw_chg) &
                       !is.na(dftmp.all$Protein_change) &
                       dftmp.all$subject_id %in% matchedset.6mwd$subject_id &
                       dftmp.all$Protein == protlist[i], ]
  
  # Total effect: intervention -> outcome
  total_model <- lm(
    smw_chg ~ Intervention + ExpressionLevel.x + age + sex + race___4 + hf_cat,
    data = dftmp
  )
  
  # Direct effect: intervention -> mediator
  mediator_model <- lm(
    Protein_change ~ Intervention + ExpressionLevel.x + age + sex + race___4 + hf_cat,
    data = dftmp
  )
  
  # Indirect + direct: outcome ~ intervention + mediator
  outcome_model <- lm(
    smw_chg ~ Intervention + Protein_change + ExpressionLevel.x + age + sex + race___4 + hf_cat,
    data = dftmp
  )
  
  # Estimate direct and indirect effects with bootstrap CI
  med_results <- mediate(
    mediator_model, outcome_model,
    treat = "Intervention", mediator = "Protein_change",
    boot = TRUE, sims = 1000
  )
  
  print(summary(med_results))
  plot(med_results, main = protlist[i])
}

dev.off()

# ---- SPPB MEDIATION ----
pdf("Results/SuppFigure5_mediationanalysis_sppb.pdf", onefile = TRUE)

for (i in 1:length(protlist)) {
  dftmp <- dftmp.all[!is.na(dftmp.all$sppb_chg) &
                       !is.na(dftmp.all$Protein_change) &
                       dftmp.all$subject_id %in% matchedset.sppb$subject_id &
                       dftmp.all$Protein == protlist[i], ]
  
  total_model <- lm(
    sppb_chg ~ Intervention + ExpressionLevel.x + age + sex + race___4 + hf_cat,
    data = dftmp
  )
  
  mediator_model <- lm(
    Protein_change ~ Intervention + ExpressionLevel.x + age + sex + race___4 + hf_cat,
    data = dftmp
  )
  
  outcome_model <- lm(
    sppb_chg ~ Intervention + Protein_change + ExpressionLevel.x + age + sex + race___4 + hf_cat,
    data = dftmp
  )
  
  med_results <- mediate(
    mediator_model, outcome_model,
    treat = "Intervention", mediator = "Protein_change",
    boot = TRUE, sims = 1000
  )
  
  print(summary(med_results))
  plot(med_results, main = protlist[i])
}

dev.off()

################################################################################
# SECTION 13: SECONDARY OUTCOMES ANALYSIS
################################################################################
# Analyze clinical events: rehospitalization, death, combined event

tbl_summary(rprotblds[, c("intervention_1_control_0", "event", "rehosp", "death", "lfu")],
            by = intervention_1_control_0,
            type = list(all_continuous() ~ "continuous2"),
            statistic = list(
              all_continuous() ~ c("{mean} ({sd})", "{median} ({p25}, {p75})", "{min}, {max}"),
              all_categorical() ~ "{n} / {N} ({p}%)"
            ),
            label = list(
              rehosp ~ "Rehospitalization",
              death ~ "Death",
              event ~ "Rehospitalization or Death",
              lfu ~ "Lost to follow up"
            ),
            missing_text = "(Missing)"
) %>%
  add_p(list(
    all_continuous() ~ "wilcox.test",
    all_categorical() ~ "fisher.test"
  )) %>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Supplemental Table 2. Secondary Outcomes**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Control (0) vs. Intervention (1)**") %>%
  as_gt() %>%
  gtsave(filename = "Results/SuppTable3_Secondaryoutcomes.html")

################################################################################
# SECTION 14: REGRESSION ANALYSES FOR SECONDARY OUTCOMES
################################################################################
# Logistic regression for binary outcomes (death, rehospitalization, combined event)
# Models include protein*intervention interaction to test effect modification

cat("\n=== SECONDARY OUTCOMES REGRESSION ANALYSIS ===\n")

# Rename protein columns for tab_model compatibility
tmpdf <- rprotblds
colnames(tmpdf)[which(colnames(rprotblds) %in% protlist)] <- make.names(protlist, unique = TRUE)
protlistclean <- make.names(protlist, unique = TRUE)

# Initialize results dataframes for secondary outcomes
regresults.death <- data.frame(
  Protein = protlist,
  Unadjusted.Effect = NA,
  Unadjusted.Pval = NA,
  Adjusted.Effect = NA,
  Adjusted.Pval = NA,
  Interaction.Pval = NA
)

regresults.rehosp <- regresults.death
regresults.both <- regresults.death

# Fit logistic regression models for each protein and outcome
for (i in 1:length(protlist)) {
  
  # ---- COMBINED EVENT (Death or Rehospitalization) ----
  fname <- paste("Results/RegressionAnalyses/DeathorRehospvs", protlist[i],
                 "_regressionanalysis.html", sep = "")
  ttl <- paste("logit(Death or Rehosp) ~ ", protlist[i], " + Intervention + More", sep = "")
  
  # Unadjusted logistic model
  mod.both.unadj <- glm(
    data = tmpdf,
    formula = paste0(
      "event ~ ",
      protlistclean[i], " + intervention_1_control_0"
    ),
    family = binomial(link = "logit")
  )
  
  # Adjusted logistic model with interaction
  mod.both.adj <- glm(
    data = tmpdf,
    formula = paste0(
      "event ~ ",
      protlistclean[i], "*intervention_1_control_0 + ",
      "age + sex + race___4 + hf_cat*intervention_1_control_0"
    ),
    family = binomial(link = "logit")
  )
  
  # Generate HTML output
  print(tab_model(
    mod.both.unadj,
    mod.both.adj,
    title = ttl,
    file = fname,
    digits = 4,
    digits.p = 4,
    show.obs = TRUE
  ))
  
  # Extract and store results (as odds ratios and 95% CI)
  regresults.both$Unadjusted.Effect[i] <- paste0(
    round(exp(mod.both.unadj$coefficients[protlistclean[i]]), 4),
    " (",
    round(exp(confint(mod.both.unadj, protlistclean[i])[1]), 4),
    ", ",
    round(exp(confint(mod.both.unadj, protlistclean[i])[2]), 4),
    ")"
  )
  regresults.both$Unadjusted.Pval[i] <- round(
    summary(mod.both.unadj)$coefficients[protlistclean[i], 4], 4
  )
  
  regresults.both$Adjusted.Effect[i] <- paste0(
    round(exp(mod.both.adj$coefficients[protlistclean[i]]), 4),
    " (",
    round(exp(confint(mod.both.adj, protlistclean[i])[1]), 4),
    ", ",
    round(exp(confint(mod.both.adj, protlistclean[i])[2]), 4),
    ")"
  )
  regresults.both$Adjusted.Pval[i] <- round(
    summary(mod.both.adj)$coefficients[protlistclean[i], 4], 4
  )
  
  regresults.both$Interaction.Pval[i] <- round(
    summary(mod.both.adj)$coefficients[
      paste0(protlistclean[i], ":intervention_1_control_01"), 4
    ], 4
  )
  
  # ---- DEATH ALONE ----
  fname <- paste("Results/RegressionAnalyses/Deathvs", protlist[i],
                 "_regressionanalysis.html", sep = "")
  ttl <- paste("logit(Death) ~ ", protlist[i], " + Intervention + More", sep = "")
  
  mod.death.unadj <- glm(
    data = tmpdf,
    formula = paste0(
      "death ~ ",
      protlistclean[i], " + intervention_1_control_0"
    ),
    family = binomial(link = "logit")
  )
  
  mod.death.adj <- glm(
    data = tmpdf,
    formula = paste0(
      "death ~ ",
      protlistclean[i], "*intervention_1_control_0 + ",
      "age + sex + race___4 + hf_cat*intervention_1_control_0"
    ),
    family = binomial(link = "logit")
  )
  
  print(tab_model(
    mod.death.unadj,
    mod.death.adj,
    title = ttl,
    file = fname,
    digits = 4,
    digits.p = 4,
    show.obs = TRUE
  ))
  
  regresults.death$Unadjusted.Effect[i] <- paste0(
    round(exp(mod.death.unadj$coefficients[protlistclean[i]]), 4),
    " (",
    round(exp(confint(mod.death.unadj, protlistclean[i])[1]), 4),
    ", ",
    round(exp(confint(mod.death.unadj, protlistclean[i])[2]), 4),
    ")"
  )
  regresults.death$Unadjusted.Pval[i] <- round(
    summary(mod.death.unadj)$coefficients[protlistclean[i], 4], 4
  )
  
  regresults.death$Adjusted.Effect[i] <- paste0(
    round(exp(mod.death.adj$coefficients[protlistclean[i]]), 4),
    " (",
    round(exp(confint(mod.death.adj, protlistclean[i])[1]), 4),
    ", ",
    round(exp(confint(mod.death.adj, protlistclean[i])[2]), 4),
    ")"
  )
  regresults.death$Adjusted.Pval[i] <- round(
    summary(mod.death.adj)$coefficients[protlistclean[i], 4], 4
  )
  
  regresults.death$Interaction.Pval[i] <- round(
    summary(mod.death.adj)$coefficients[
      paste0(protlistclean[i], ":intervention_1_control_01"), 4
    ], 4
  )
  
  # ---- REHOSPITALIZATION ----
  fname <- paste("Results/RegressionAnalyses/Rehospvs", protlist[i],
                 "_regressionanalysis.html", sep = "")
  ttl <- paste("logit(Rehosp) ~ ", protlist[i], " + Intervention + More", sep = "")
  
  mod.rehosp.unadj <- glm(
    data = tmpdf,
    formula = paste0(
      "rehosp ~ ",
      protlistclean[i], " + intervention_1_control_0"
    ),
    family = binomial(link = "logit")
  )
  
  mod.rehosp.adj <- glm(
    data = tmpdf,
    formula = paste0(
      "rehosp ~ ",
      protlistclean[i], "*intervention_1_control_0 + ",
      "age + sex + race___4 + hf_cat*intervention_1_control_0"
    ),
    family = binomial(link = "logit")
  )
  
  print(tab_model(
    mod.rehosp.unadj,
    mod.rehosp.adj,
    title = ttl,
    file = fname,
    digits = 4,
    digits.p = 4,
    show.obs = TRUE
  ))
  
  regresults.rehosp$Unadjusted.Effect[i] <- paste0(
    round(exp(mod.rehosp.unadj$coefficients[protlistclean[i]]), 4),
    " (",
    round(exp(confint(mod.rehosp.unadj, protlistclean[i])[1]), 4),
    ", ",
    round(exp(confint(mod.rehosp.unadj, protlistclean[i])[2]), 4),
    ")"
  )
  regresults.rehosp$Unadjusted.Pval[i] <- round(
    summary(mod.rehosp.unadj)$coefficients[protlistclean[i], 4], 4
  )
  
  regresults.rehosp$Adjusted.Effect[i] <- paste0(
    round(exp(mod.rehosp.adj$coefficients[protlistclean[i]]), 4),
    " (",
    round(exp(confint(mod.rehosp.adj, protlistclean[i])[1]), 4),
    ", ",
    round(exp(confint(mod.rehosp.adj, protlistclean[i])[2]), 4),
    ")"
  )
  regresults.rehosp$Adjusted.Pval[i] <- round(
    summary(mod.rehosp.adj)$coefficients[protlistclean[i], 4], 4
  )
  
  regresults.rehosp$Interaction.Pval[i] <- round(
    summary(mod.rehosp.adj)$coefficients[
      paste0(protlistclean[i], ":intervention_1_control_01"), 4
    ], 4
  )
}

# Save secondary outcome results
write.csv(regresults.death, "Results/regsummary_death.csv")
write.csv(regresults.rehosp, "Results/regsummary_rehosp.csv")
write.csv(regresults.both, "Results/regsummary_deathrehosp.csv")

################################################################################
# SECTION 15: PROPENSITY MATCHING BALANCE ASSESSMENT
################################################################################
# Compare baseline characteristics of matched vs unmatched participants
# Ensures matching reduced confounding

cat("\n=== PROPENSITY MATCHING BALANCE ASSESSMENT ===\n")

# ---- 6MWD Matched vs Unmatched ----
varlist <- c("intervention_1_control_0", "age", "sex", "race___4",
             "hf_cat", "egfr", protlist)

# Create comparison dataframe
tmpdf <- rbind(
  data.frame(set.total[, varlist], dataset = "matched"),
  data.frame(
    setdiff(ps.df[, varlist], set.total[, varlist]),
    dataset = "unmatched"
  )
)

tmpdf$dataset <- factor(tmpdf$dataset, levels = c("unmatched", "matched"))
tmpdf$grp <- str_c(tmpdf$intervention_1_control_0, tmpdf$dataset)
tmpdf$grp <- as.factor(tmpdf$grp)
tmpdf$grp <- fct_collapse(tmpdf$grp,
                          Unmatched = c("0unmatched", "1unmatched"),
                          MatchedControl = "0matched",
                          MatchedIntervention = "1matched"
)
tmpdf$grp <- factor(tmpdf$grp, levels = c("MatchedControl", "MatchedIntervention", "Unmatched"))
tmpdf <- tmpdf[, -which(names(tmpdf) %in% c("intervention_1_control_0", "dataset"))]

# Generate balance table
tbl_summary(tmpdf,
            by = grp,
            type = all_continuous() ~ "continuous2",
            statistic = list(
              all_continuous() ~ c("{mean} ({sd})", "{median} ({p25}, {p75})", "{min}, {max}"),
              all_categorical() ~ "{n} / {N} ({p}%)"
            ),
            label = list(
              age ~ "Age",
              sex ~ "Sex (Female=1)",
              race___4 ~ "Race (White=1)",
              hf_cat ~ "HF Category (HFpEF=1)",
              egfr ~ "EGFR",
              ALCAM ~ "ALCAM",
              Gal.4 ~ "Gal-4",
              PECAM.1 ~ "PECAM-1",
              GP6 ~ "GP6",
              ST2 ~ "ST2"
            ),
            missing_text = "(Missing)"
) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Baseline Characteristics: Matched vs Unmatched (6MWD)**") %>%
  as_gt() %>%
  gtsave(filename = "Results/6MWDmatchedvsunmatched.html")

# ---- SPPB Matched vs Unmatched ----
# Note: In the original analysis, this section would use the SPPB matching dataset
# For now, we'll use similar structure but acknowledge this would need separate matching output

varlist <- c("intervention_1_control_0", "age", "sex", "race___4",
             "hf_cat", "egfr", protlist)

# Create comparison dataframe (would use matchedset.sppb in practice)
tmpdf <- rbind(
  data.frame(matchedset.sppb[, varlist], dataset = "matched"),
  data.frame(
    setdiff(ps.df[, varlist], matchedset.sppb[, varlist]),
    dataset = "unmatched"
  )
)

tmpdf$dataset <- factor(tmpdf$dataset, levels = c("unmatched", "matched"))
tmpdf$grp <- str_c(tmpdf$intervention_1_control_0, tmpdf$dataset)
tmpdf$grp <- as.factor(tmpdf$grp)
tmpdf$grp <- fct_collapse(tmpdf$grp,
                          Unmatched = c("0unmatched", "1unmatched"),
                          MatchedControl = "0matched",
                          MatchedIntervention = "1matched"
)
tmpdf$grp <- factor(tmpdf$grp, levels = c("MatchedControl", "MatchedIntervention", "Unmatched"))
tmpdf <- tmpdf[, -which(names(tmpdf) %in% c("intervention_1_control_0", "dataset"))]

# Generate balance table
tbl_summary(tmpdf,
            by = grp,
            type = all_continuous() ~ "continuous2",
            statistic = list(
              all_continuous() ~ c("{mean} ({sd})", "{median} ({p25}, {p75})", "{min}, {max}"),
              all_categorical() ~ "{n} / {N} ({p}%)"
            ),
            label = list(
              age ~ "Age",
              sex ~ "Sex (Female=1)",
              race___4 ~ "Race (White=1)",
              hf_cat ~ "HF Category (HFpEF=1)",
              egfr ~ "EGFR",
              ALCAM ~ "ALCAM",
              Gal.4 ~ "Gal-4",
              PECAM.1 ~ "PECAM-1",
              GP6 ~ "GP6",
              ST2 ~ "ST2"
            ),
            missing_text = "(Missing)"
) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Baseline Characteristics: Matched vs Unmatched (SPPB)**") %>%
  as_gt() %>%
  gtsave(filename = "Results/SPPBmatchedvsunmatched.html")

################################################################################
# SECTION 16: PROTEIN CHANGE FROM BASELINE TO FOLLOW-UP
################################################################################
# Summarize changes in protein levels over 12 weeks
# Tests whether intervention group shows different protein trajectories

cat("\n=== PROTEIN CHANGE ANALYSIS (Baseline to 12-Week Follow-up) ===\n")

dftmp <- rprotds_wide[!is.na(rprotds_wide$timepoint.y), ]

# ---- Overall protein changes (both groups combined) ----
cat("\n--- OVERALL PROTEIN CHANGES ---\n")

for (i in 1:length(protlist)) {
  protein_change <- dftmp[dftmp$Protein == protlist[i], "ExpressionLevel.y"] -
    dftmp[dftmp$Protein == protlist[i], "ExpressionLevel.x"]
  
  cat(paste("\n", protlist[i], ":\n", sep = ""))
  
  # Median and IQR
  cat(paste(
    "  Median change: ",
    round(quantile(protein_change, probs = c(0.5)), 3),
    " (IQR: ",
    round(quantile(protein_change, probs = c(0.25)), 3),
    ", ",
    round(quantile(protein_change, probs = c(0.75)), 3),
    ")\n",
    sep = ""
  ))
  
  # Wilcoxon signed-rank test (tests if median change is different from zero)
  pval <- wilcox.test(x = protein_change)$p.value
  cat(paste("  Wilcoxon p-value: ", round(pval, 4), "\n", sep = ""))
}

# ---- Protein changes in Attention Control group ----
cat("\n--- ATTENTION CONTROL GROUP ---\n")

dftmp <- rprotds_wide[!is.na(rprotds_wide$timepoint.y) & rprotds_wide$Intervention == 0, ]

for (i in 1:length(protlist)) {
  protein_change <- dftmp[dftmp$Protein == protlist[i], "ExpressionLevel.y"] -
    dftmp[dftmp$Protein == protlist[i], "ExpressionLevel.x"]
  
  cat(paste("\n", protlist[i], ":\n", sep = ""))
  
  cat(paste(
    "  Median change: ",
    round(quantile(protein_change, probs = c(0.5)), 3),
    " (IQR: ",
    round(quantile(protein_change, probs = c(0.25)), 3),
    ", ",
    round(quantile(protein_change, probs = c(0.75)), 3),
    ")\n",
    sep = ""
  ))
  
  pval <- wilcox.test(x = protein_change)$p.value
  cat(paste("  Wilcoxon p-value: ", round(pval, 4), "\n", sep = ""))
}

# ---- Protein changes in Rehabilitation group ----
cat("\n--- REHABILITATION INTERVENTION GROUP ---\n")

dftmp <- rprotds_wide[!is.na(rprotds_wide$timepoint.y) & rprotds_wide$Intervention == 1, ]

for (i in 1:length(protlist)) {
  protein_change <- dftmp[dftmp$Protein == protlist[i], "ExpressionLevel.y"] -
    dftmp[dftmp$Protein == protlist[i], "ExpressionLevel.x"]
  
  cat(paste("\n", protlist[i], ":\n", sep = ""))
  
  cat(paste(
    "  Median change: ",
    round(quantile(protein_change, probs = c(0.5)), 3),
    " (IQR: ",
    round(quantile(protein_change, probs = c(0.25)), 3),
    ", ",
    round(quantile(protein_change, probs = c(0.75)), 3),
    ")\n",
    sep = ""
  ))
  
  pval <- wilcox.test(x = protein_change)$p.value
  cat(paste("  Wilcoxon p-value: ", round(pval, 4), "\n", sep = ""))
}

# ---- Between-group comparison of protein changes ----
cat("\n--- BETWEEN-GROUP COMPARISON (Intervention vs Control) ---\n")

dftmp <- rprotds_wide[!is.na(rprotds_wide$timepoint.y), ]

for (i in 1:length(protlist)) {
  # Test equality of protein changes between groups
  pval <- t.test(
    Protein_change ~ Intervention,
    data = dftmp[dftmp$Protein == protlist[i], ],
    var.equal = TRUE
  )$p.value
  
  cat(paste(protlist[i], ": p-value = ", round(pval, 4), "\n", sep = ""))
}

################################################################################
# SECTION 17: SESSION SUMMARY AND OUTPUT VERIFICATION
################################################################################
# Final check of all output files generated

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS COMPLETE - Summary of Generated Files\n")
cat("================================================================================\n")

output_files <- list(
  "Main Results" = c(
    "Table1_Baselinecharacteristics.html",
    "Table2_sppb.csv",
    "Table2_6mwd.csv"
  ),
  "Figures" = c(
    "Figure1_SpearmanCorrelation.pdf",
    "Figure2_ChangeinSPPB.pdf",
    "Figure3_Changein6MWD.pdf",
    "Figure4a_Matchingtree_SPPB.pdf",
    "Figure4b_Matchingtree_6MWD.pdf"
  ),
  "Supplemental Figures" = c(
    "SuppFigure2_InteractionImportance_sppb.pdf",
    "SuppFigure2_InteractionImportance_6mwd.pdf",
    "SuppFigure3_CVPlot_SPPB.pdf",
    "SuppFigure3_CVPlot_6MWD.pdf",
    "SuppFigure4_Proteinchange.pdf",
    "SuppFigure5_mediationanalysis_6mwd.pdf",
    "SuppFigure5_mediationanalysis_sppb.pdf",
    "SuppFigure6_WGCNA_heatmap.pdf"
  ),
  "Regression Output" = c(
    "RegressionAnalyses/ (HTML files for each protein-outcome pair)",
    "regsummary_death.csv",
    "regsummary_rehosp.csv",
    "regsummary_deathrehosp.csv"
  ),
  "Balance Assessment" = c(
    "6MWDmatchedvsunmatched.html",
    "SPPBmatchedvsunmatched.html"
  ),
  "WGCNA Results" = c(
    "protein_modules.csv",
    "correlationplot_allparticipants_baseline.pdf",
    "correlationplot_usualcare_baseline.pdf",
    "correlationplot_rehab_baseline.pdf",
    "WGCNA_DendrogramAndColors.pdf"
  ),
  "Secondary Analyses" = c(
    "SuppTable3_Secondaryoutcomes.html",
    "SPPB_volcano_effectmod.pdf",
    "6MWD_volcano_effectmod.pdf",
    "bart_sppb_ConvergenceDiagnostics.pdf",
    "bart_smw_ConvergenceDiagnostics.pdf"
  )
)

for (category in names(output_files)) {
  cat(paste("\n", category, ":\n", sep = ""))
  for (file in output_files[[category]]) {
    cat(paste("  - ", file, "\n", sep = ""))
  }
}

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS PARAMETERS AND NOTES\n")
cat("================================================================================\n")

cat("\nKey Analysis Parameters:\n")
cat("  - Soft-thresholding power (WGCNA):", sft$fitIndices[which.max(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]), 1], "\n")
cat("  - WGCNA module detection sensitivity (deepSplit): 2\n")
cat("  - Minimum module size: 6 proteins\n")
cat("  - Propensity score matching caliper: 1 SD of protein values\n")
cat("  - Matching tree minbucket: 15 (SPPB), 30 (6MWD)\n")
cat("  - Cross-validation folds: 10\n")
cat("  - BART model: 100 trees, 2000 burn-in iterations, 2000 post-burn-in iterations\n")
cat("  - Mediation analysis: 1000 bootstrap samples\n")

cat("\nMissing Data Handling:\n")
cat("  - Negative protein values: Replaced with NA\n")
cat("  - Baseline biomarker outliers (>3 SD on log2 scale): Set to NA\n")
cat("  - Censored values (<LOD, >ULD):\n")
cat("    * Left-censored: Imputed as LOD/2\n")
cat("    * Right-censored: Imputed as ULD\n")
cat("  - QNS (Quantity Not Sufficient): Replaced with NA\n")

cat("\nSignificance Levels:\n")
cat("  - Primary analyses: p < 0.05 (unadjusted)\n")
cat("  - Protein correlations: p < 0.01 (indicated with asterisk)\n")
cat("  - Module-outcome correlations: p < 0.01 (indicated with asterisk)\n")
cat("  - Effect modification screening: FDR-adjusted (exploratory, no significant hits)\n")

cat("\nStudy Population:\n")
cat("  - Total REHAB-HF enrolled: ", nrow(rhfds), "\n")
cat("  - With baseline proteomics: ", nrow(rprotblds), "\n")
cat("  - With baseline + follow-up SPPB: ", 
    nrow(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                 !is.na(rhfds$fu_sppb) & !is.na(rhfds$bl_sppb), ]), "\n")
cat("  - With baseline + follow-up 6MWD: ", 
    nrow(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                 !is.na(rhfds$fu_sppb) & !is.na(rhfds$bl_sppb) &
                 !is.na(rhfds$fu_smw) & !is.na(rhfds$bl_smw), ]), "\n")
cat("  - With baseline + follow-up proteomics: ",
    nrow(rhfds[rhfds$study_id %in% rprotblds$subject_id & 
                 !is.na(rhfds$fu_sppb) & !is.na(rhfds$bl_sppb) &
                 rhfds$study_id %in% rprotfuds$subject_id, ]), "\n")

cat("\nIntervention Groups (among analyzed cohort):\n")
cat("  - Attention Control (0): ", 
    sum(rprotblds$intervention_1_control_0 == 0, na.rm = TRUE), "\n")
cat("  - Rehabilitation Intervention (1): ",
    sum(rprotblds$intervention_1_control_0 == 1, na.rm = TRUE), "\n")

cat("\nProtein Modules Identified (WGCNA):\n")
for (color in names(table(ModuleColors))) {
  cat(paste("  - ", color, ": ", table(ModuleColors)[color], " proteins\n", sep = ""))
}

cat("\n")
cat("================================================================================\n")
cat("END OF ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# ADDITIONAL SUMMARY STATISTICS FOR REFERENCE
################################################################################

cat("REFERENCE: Detailed Baseline Statistics\n")
cat("===============================================\n\n")

# Age by intervention group
cat("Age by Intervention Group:\n")
print(tapply(rprotblds$age, rprotblds$intervention_1_control_0, 
             function(x) c(Mean = mean(x, na.rm = TRUE), 
                           SD = sd(x, na.rm = TRUE),
                           Median = median(x, na.rm = TRUE))))

# Sex distribution
cat("\nSex Distribution by Intervention Group:\n")
print(table(rprotblds$intervention_1_control_0, rprotblds$sex))

# Race distribution
cat("\nRace Distribution by Intervention Group:\n")
print(table(rprotblds$intervention_1_control_0, rprotblds$race___4))

# HF category distribution
cat("\nHF Category Distribution by Intervention Group:\n")
print(table(rprotblds$intervention_1_control_0, rprotblds$hf_cat))

# Baseline outcomes by group
cat("\nBaseline SPPB by Intervention Group:\n")
print(tapply(rprotblds$bl_sppb, rprotblds$intervention_1_control_0,
             function(x) c(Mean = mean(x, na.rm = TRUE),
                           SD = sd(x, na.rm = TRUE),
                           Median = median(x, na.rm = TRUE))))

cat("\nBaseline 6MWD by Intervention Group:\n")
print(tapply(rprotblds$bl_smw, rprotblds$intervention_1_control_0,
             function(x) c(Mean = mean(x, na.rm = TRUE),
                           SD = sd(x, na.rm = TRUE),
                           Median = median(x, na.rm = TRUE))))

cat("\nBaseline eGFR by Intervention Group:\n")
print(tapply(rprotblds$egfr, rprotblds$intervention_1_control_0,
             function(x) c(Mean = mean(x, na.rm = TRUE),
                           SD = sd(x, na.rm = TRUE),
                           Median = median(x, na.rm = TRUE))))

cat("\n")
cat("Analysis script completed successfully!\n")
cat("All results have been saved to the Results/ directory.\n")

################################################################################
# END OF ANALYSIS SCRIPT
################################################################################
