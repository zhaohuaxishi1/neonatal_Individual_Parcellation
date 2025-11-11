# run_GAM_voxelwise_networks.R
# Perform voxel-wise GAM analyses for each probabilistic network
# Author: jlzhao

library(R.matlab)
library(mgcv)
library(nlme)
library(visreg)

# -------------------------------------------------------------
# Configuration
# -------------------------------------------------------------
cat("Starting voxel-wise GAM analysis ...\n")

# Relative paths (consistent with GitHub repo structure)
data_dir <- "data"
atlas_dir <- file.path(data_dir, "network_loading")
output_dir <- "output/gam_voxelwise"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Network labels (to name outputs)
network_labels <- c(
  "FroLim", "aDMN", "DorFro", "FootMot", "HandMot",
  "SupPar", "rDMN", "lDMN", "pDMN", "LatVis", "PriVis"
)

# Parameters
n_subjects <- 315
n_voxels <- 6905

# -------------------------------------------------------------
# Load behavioral variables
# -------------------------------------------------------------
behavior_data <- readMat(file.path(data_dir, "AgeScore.mat"))

scan_age     <- as.numeric(behavior_data$AgeScore[[1]])
birth_age    <- as.numeric(behavior_data$AgeScore[[2]])
sex          <- as.factor(behavior_data$AgeScore[[3]])
time_interval <- as.numeric(behavior_data$AgeScore[[4]])
mean_fd      <- as.numeric(behavior_data$AgeScore[[5]])

Behavior <- data.frame(
  scan_age = scan_age,
  birth_age = birth_age,
  sex = sex,
  TimeInterval = time_interval,
  mFD = mean_fd
)

# -------------------------------------------------------------
# Load network loading matrix
# -------------------------------------------------------------
atlas_file <- file.path(atlas_dir, "V_vector_age.mat")
atlas_data <- readMat(atlas_file)
V_vector <- atlas_data$V.vector

# -------------------------------------------------------------
# GAM analysis loop for each network
# -------------------------------------------------------------
for (net_idx in 1:length(network_labels)) {
  network_name <- network_labels[net_idx]
  cat(sprintf("Running voxel-wise GAM for network: %s\n", network_name))

  # Extract network-specific voxel loadings
  net_data <- V_vector[, ((net_idx - 1) * n_voxels + 1):(net_idx * n_voxels)]

  # Remove zero columns (non-informative voxels)
  valid_cols <- which(colSums(net_data) != 0)
  net_data_nonzero <- net_data[, valid_cols]
  n_features <- ncol(net_data_nonzero)

  # Initialize result vectors
  p_vec <- z_vec <- f_vec <- numeric(n_features)

  for (v in 1:n_features) {
    if (v %% 500 == 0) cat("  Processing voxel", v, "of", n_features, "\n")

    y <- as.numeric(net_data_nonzero[, v])
    gam_model <- gam(
      y ~ s(scan_age, bs = "cs", k = 3) + sex + TimeInterval + mFD,
      data = Behavior, method = "REML", na.action = na.omit
    )

    gam_summary <- summary(gam_model)
    p_vec[v] <- gam_summary$s.table["s(scan_age)", "p-value"]
    f_vec[v] <- gam_summary$s.table["s(scan_age)", "F"]
    z_vec[v] <- qnorm(p_vec[v] / 2, lower.tail = FALSE)

    # Determine direction (sign) of effect
    lm_model <- lm(y ~ scan_age + sex + TimeInterval + mFD, data = Behavior)
    if (summary(lm_model)$coefficients[2, 3] < 0) {
      z_vec[v] <- -z_vec[v]
    }
  }

  # ---------------------------------------------------------
  # Multiple comparison correction
  # ---------------------------------------------------------
  p_fdr <- p.adjust(p_vec, method = "fdr")
  p_bonf <- p.adjust(p_vec, method = "bonferroni")

  z_fdr <- ifelse(p_fdr < 0.05, z_vec, 0)
  f_fdr <- ifelse(p_fdr < 0.05, f_vec, 0)
  z_bonf <- ifelse(p_bonf < 0.05, z_vec, 0)
  f_bonf <- ifelse(p_bonf < 0.05, f_vec, 0)

  # ---------------------------------------------------------
  # Write full-length vectors (with zeros for removed voxels)
  # ---------------------------------------------------------
  p_all <- z_all <- p_fdr_all <- z_fdr_all <- p_bonf_all <- z_bonf_all <- f_all <- f_fdr_all <- f_bonf_all <- rep(0, n_voxels)

  p_all[valid_cols] <- p_vec
  z_all[valid_cols] <- z_vec
  p_fdr_all[valid_cols] <- p_fdr
  z_fdr_all[valid_cols] <- z_fdr
  p_bonf_all[valid_cols] <- p_bonf
  z_bonf_all[valid_cols] <- z_bonf
  f_all[valid_cols] <- f_vec
  f_fdr_all[valid_cols] <- f_fdr
  f_bonf_all[valid_cols] <- f_bonf

  # ---------------------------------------------------------
  # Save results
  # ---------------------------------------------------------
  out_file <- file.path(output_dir, paste0("AgeEffect_Voxelwise_Network_", network_name, ".mat"))
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
}

cat("✅ Finished voxel-wise GAM analyses. Results saved to:", output_dir, "\n")
