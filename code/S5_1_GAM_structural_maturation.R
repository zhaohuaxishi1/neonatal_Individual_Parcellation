# run_GAM_structural_maturation.R
# Generalized Additive Model (GAM) analyses of structural maturation
# Author: jlzhao
# Based on: "Neonatal functional topographic refinement and cortical structure maturation"
#
# This script quantifies age-dependent changes in cortical structural features
# (myelination, thickness, curvature, and sulcal depth)
# during the early postnatal period in 301 term-born neonates (Subset 2).
#
# GAM model: feature ~ s(scan_age, bs="cs", k=3) + sex + mFD + TimeInterval
# FDR-corrected voxelwise significance maps are saved as .mat files.

library(R.matlab)
library(mgcv)
library(nlme)
library(visreg)

cat("Starting GAM analyses for cortical structure maturation...\n")

# --------------------------- Configuration ---------------------------
data_dir <- "data/structure_basis"
output_dir <- file.path("output", "gam_structurematuration")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_subjects <- 301
n_vertices <- 5124

# Structural metrics:
# 1=myelin, 2=thickness, 3=curvature, 4=sulcal depth, others for derived measures
metric_labels <- c("SulcalDepth","Curvature", "Myelination", "Thickness")

# --------------------------- Load Behavioral Data ---------------------------
age_data <- readMat(file.path(data_dir, "CombinedAgeScoreStructure.mat"))
scan_age <- as.numeric(age_data$CombinedAgeScore[[1]])
birth_age <- as.numeric(age_data$CombinedAgeScore[[2]])
sex <- as.factor(age_data$CombinedAgeScore[[3]])
mFD <- as.numeric(age_data$CombinedAgeScore[[4]])
TimeInterval <- as.numeric(age_data$CombinedAgeScore[[5]])

Behavior <- data.frame(scan_age, birth_age, sex, mFD, TimeInterval)

# --------------------------- Run GAM for Each Metric ---------------------------
metric_file <- file.path(data_dir, "sub_metric_all_term.mat")
metric_data <- readMat(metric_file)
all_metrics <- metric_data$sub.metric.all

for (metric_idx in seq_along(all_metrics)) {
  metric_name <- metric_labels[metric_idx]
  cat(sprintf("\nRunning GAM for %s (Metric %d of %d)\n",
              metric_name, metric_idx, length(all_metrics)))

  Measure_I <- all_metrics[[metric_idx]][[1]]
  valid_cols <- which(colSums(Measure_I) != 0)
  Measure_nonzero <- Measure_I[, valid_cols]
  n_features <- ncol(Measure_nonzero)

  # Initialize vectors
  p_vec <- z_vec <- f_vec <- numeric(n_features)

  # --------------------------- GAM Loop ---------------------------
  for (v in seq_len(n_features)) {
    if (v %% 500 == 0) cat("  Vertex", v, "of", n_features, "\n")

    y <- as.numeric(Measure_nonzero[, v])
    gam_model <- gam(y ~ s(scan_age, bs = "cs", k = 3) + sex + TimeInterval + mFD,
                     data = Behavior, method = "REML", na.action = na.omit)
    gam_summary <- summary(gam_model)

    p_vec[v] <- gam_summary$s.table["s(scan_age)", "p-value"]
    f_vec[v] <- gam_summary$s.table["s(scan_age)", "F"]
    z_vec[v] <- qnorm(p_vec[v] / 2, lower.tail = FALSE)

    # Determine direction from linear model
    lm_model <- lm(y ~ scan_age + sex + TimeInterval + mFD, data = Behavior)
    if (summary(lm_model)$coefficients[2, 3] < 0) {
      z_vec[v] <- -z_vec[v]
      f_vec[v] <- -f_vec[v]
    }
  }

  # --------------------------- Multiple Comparison Correction ---------------------------
  p_fdr <- p.adjust(p_vec, method = "fdr")
  p_bonf <- p.adjust(p_vec, method = "bonferroni")

  z_fdr <- ifelse(p_fdr < 0.05, z_vec, 0)
  f_fdr <- ifelse(p_fdr < 0.05, f_vec, 0)
  z_bonf <- ifelse(p_bonf < 0.05, z_vec, 0)
  f_bonf <- ifelse(p_bonf < 0.05, f_vec, 0)

  # Fill full-length vectors (with zeros for non-brain vertices)
  p_all <- z_all <- p_fdr_all <- z_fdr_all <- p_bonf_all <- z_bonf_all <-
    f_all <- f_fdr_all <- f_bonf_all <- rep(0, n_vertices)

  p_all[valid_cols] <- p_vec
  z_all[valid_cols] <- z_vec
  p_fdr_all[valid_cols] <- p_fdr
  z_fdr_all[valid_cols] <- z_fdr
  p_bonf_all[valid_cols] <- p_bonf
  z_bonf_all[valid_cols] <- z_bonf
  f_all[valid_cols] <- f_vec
  f_fdr_all[valid_cols] <- f_fdr
  f_bonf_all[valid_cols] <- f_bonf

  # --------------------------- Save Results ---------------------------
  out_file <- file.path(output_dir, paste0("AgeEffect_Structure_", metric_name, ".mat"))
  writeMat(out_file,
           Gam_P_Vector_All = p_all,
           Gam_Z_Vector_All = z_all,
           Gam_P_FDR_Vector_All = p_fdr_all,
           Gam_Z_FDR_Sig_Vector_All = z_fdr_all,
           Gam_P_Bonf_Vector_All = p_bonf_all,
           Gam_Z_Bonf_Sig_Vector_All = z_bonf_all,
           Gam_F_Vector_All = f_all,
           Gam_F_FDR_Sig_Vector_All = f_fdr_all,
           Gam_F_Bonf_Sig_Vector_All = f_bonf_all)

  cat(sprintf("✅ Saved results for %s → %s\n", metric_name, out_file))
}

cat("\n✅ Finished GAM analyses for all 8 structural metrics.\n")
