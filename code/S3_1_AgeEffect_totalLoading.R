# scripts/run_GAM_total_representation.R
# GAM analysis of age-related changes in total network representation
# Author: jlzhao

library(R.matlab)
library(mgcv)
library(visreg)
library(ggplot2)

# ------------------------ Config ------------------------
network_count <- 11
result_dir <- file.path("output", "gam_totalloading")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
network_labels <- c('FroLim', 'aDMN', 'DorFro', 'FootMot', 'HandMot',
                    'SupPar', 'rDMN', 'lDMN', 'pDMN', 'LatVis', 'PriVis')

# Load data
representation <- readMat("data/TotalNetworkRepresentation.mat")$TotalNetworkRepresentation
meta <- readMat("data/AgeScore.mat")$AgeScore
scan_age <- as.numeric(meta[[1]])
birth_age <- as.numeric(meta[[2]])
sex <- as.factor(meta[[3]])
time_interval <- as.numeric(meta[[4]])
mFD <- as.numeric(meta[[5]])

# ------------------------ Init ------------------------
p_values <- numeric(network_count)
F_values <- numeric(network_count)

# Color palette for lines and points
color_line <- "#7F7F7F"
color_point <- c(
  '#E4A74C','#9A80CB','#EF7FCC','#108B7A','#84BC7D',
  '#FFEA84','#AA50B8','#AA50B8','#E33C33','#007CC8','#52AAC8'
)

# ------------------------ GAM Loop ------------------------
for (i in 1:network_count) {
  network_name <- network_labels[i]
  y <- representation[, i]
  df <- data.frame(
    scan_age = scan_age,
    sex = sex,
    time_interval = time_interval,
    mFD = mFD,
    y = y
  )

  # Fit GAM
  gam_model <- gam(y ~ s(scan_age, bs = "cs", k = 3) + sex + time_interval + mFD,
                   data = df, method = "REML")
  s <- summary(gam_model)
  p_values[i] <- s$s.table["s(scan_age)", "p-value"]
  F_values[i] <- s$s.table["s(scan_age)", "F"]

  # Visualization
  vis_data <- visreg(gam_model, "scan_age", gg = TRUE, type = "conditional",
                     scale = "response", overlay = TRUE, rug = FALSE, plot = FALSE)

  smooth_df <- data.frame(
    x = vis_data$fit$scan_age,
    y = vis_data$fit$visregFit,
    lower = vis_data$fit$visregLwr,
    upper = vis_data$fit$visregUpr
  )

  point_df <- data.frame(
    x = vis_data$res$scan_age,
    y = vis_data$res$visregRes
  )

  fig <- ggplot() +
    geom_point(data = point_df, aes(x, y), color = color_point[i], size = 2, alpha = 0.8) +
    geom_line(data = smooth_df, aes(x, y), color = color_line, size = 1.5) +
    geom_ribbon(data = smooth_df, aes(x = x, ymin = lower, ymax = upper), fill = color_line, alpha = 0.2) +
    labs(x = "Scan age (weeks)", y = "Total representation") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    scale_x_continuous(breaks = seq(37, 45, 2)) +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.25)))

  ggsave(
    filename = file.path(result_dir, paste0("GAM_Network_", i, '_', network_name, ".tiff")),
    plot = fig, width = 17, height = 15, dpi = 75, units = "cm"
  )
}

# ------------------------ FDR Correction & Save ------------------------
p_values_fdr <- p.adjust(p_values, method = "fdr")
F_values_fdr <- F_values
F_values_fdr[p_values_fdr >= 0.05] <- 0

write.table(p_values,        file = file.path(result_dir, "network_p_values.txt"),       row.names = FALSE, col.names = FALSE)
write.table(F_values,        file = file.path(result_dir, "network_F_values.txt"),       row.names = FALSE, col.names = FALSE)
write.table(p_values_fdr,    file = file.path(result_dir, "network_p_values_FDR.txt"),   row.names = FALSE, col.names = FALSE)
write.table(F_values_fdr,    file = file.path(result_dir, "network_F_values_FDR.txt"),   row.names = FALSE, col.names = FALSE)

cat("Finished GAM modeling and saved results to", result_dir, "\n")
