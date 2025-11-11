# plot_network_loading_vs_variability.R
# Visualize correlation between network loading and individual variability using hexbin
# Author: jlzhao

library(R.matlab)
library(ggplot2)
library(hexbin)

# ------------------- Configuration -------------------
network_index <- 1                  # Select network ID (1–11)
vertex_count <- 6905                # Number of cortical voxels/vertices
data_dir <- file.path("output", "variability_loading_corr")
figure_dir <- file.path("output", "variability_loading_corr")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------- Load Data -------------------
loadings_mat <- readMat(file.path(data_dir, "Network_Loading.mat"))
variability_mat <- readMat(file.path(data_dir, "Individual_Variability.mat"))

network_loading <- loadings_mat$Network.Loading[, network_index]
individual_variability <- variability_mat$Individual.Variability


# ------------------- Correlation ------------------
cor_result <- cor.test(network_loading, individual_variability, method = "pearson")
cat(sprintf("Network %d | Pearson r = %.4f | p = %.4g\n",
            network_index, cor_result$estimate, cor_result$p.value))

# ------------------- Prepare Data -------------------
df <- data.frame(
  loading = as.numeric(network_loading),
  variability = as.numeric(individual_variability)
)

hex_info <- hexbin(df$loading, df$variability)
df_hex <- data.frame(hcell2xy(hex_info), count = hex_info@count)

# ------------------- Plot -------------------
# ------------------- Color Settings -------------------
hex_palette <- c("#CCCCCC", "#4C4C4C", "#4C4C4C", "#4C4C4C")

plot_limits <- list(
  x = c(0.30, 1.00),
  y = c(0.01, 0.08)
)

p <- ggplot() +
  geom_hex(data = df_hex, aes(x, y, fill = count), stat = "identity") +
  scale_fill_gradientn(colours = hex_palette) +
  geom_smooth(data = df, aes(x = loading, y = variability),
              method = lm, color = "white", linetype = "dashed") +
  theme_classic() +
  labs(x = "Network Loading", y = "Individual Variability") +
  theme(legend.position = 'none') +
  theme(
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 16),
    aspect.ratio = 1
  ) +
  scale_x_continuous(limits = plot_limits$x, breaks = c(0.3, 0.6, 1.0)) +
  scale_y_continuous(limits = plot_limits$y, breaks = c(0.01, 0.04, 0.08))
p

# ------------------- Save Figure -------------------
output_name <- sprintf("variability_loading_corr.tiff", network_index)
ggsave(
  filename = file.path(figure_dir, output_name),
  plot = p, width = 17, height = 15, dpi = 150, units = "cm"
)
