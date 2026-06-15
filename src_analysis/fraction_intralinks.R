# Compute mixing parameters for KEGG pathways in FunCoup network

source(file.path("src_analysis", "paths.R"))

# Libraries
library(fastmatch)
library(ggplot2)
library(cowplot)
library(scales)

# --- 1. Load data ---
KEGG <- read.delim(benchmark_file("data", "KEGG_pathways"), header = FALSE)
net <- read.delim(benchmark_file("data", "fc3.tsv"), header = TRUE)
# Filter edges by weight threshold
net <- subset(net, V3 >= 0.75)[, 1:3]

paths <- unique(KEGG$V2)

# Compute degree for all nodes
deg_table <- table(c(net$V1, net$V2))
deg_df <- as.data.frame(deg_table, stringsAsFactors = FALSE)
names(deg_df) <- c("gene", "degree")

# --- 2. Function to compute mixing per pathway ---
compute_mixing <- function(path_id) {
  genes_in_path <- KEGG$V1[KEGG$V2 == path_id]
  # Subgraph edges incident to pathway
  sub_edges <- subset(net, V1 %in% genes_in_path | V2 %in% genes_in_path)
  # Edges fully inside pathway
  inside_edges <- subset(sub_edges, V1 %in% genes_in_path & V2 %in% genes_in_path)
  
  # Initialize vectors
  mixing <- numeric(); degrees <- numeric(); weighted <- numeric(); abs_diff <- numeric()
  
  for (g in genes_in_path) {
    # Edges incident to g
    incident <- subset(sub_edges, V1 == g | V2 == g)
    # Inside edges for g
    inside <- subset(inside_edges, V1 == g | V2 == g)
    total_n <- nrow(incident)
    inside_n <- nrow(inside)
    # Mixing parameter: fraction of edges leaving pathway
    mix_param <- if (total_n == 0) 0 else (total_n - inside_n) / total_n
    # Node degree from global degree table
    deg <- deg_df$degree[match(g, deg_df$gene)]
    if (is.na(deg)) deg <- 0
    # Weighted contribution
    weighted_val <- deg * mix_param
    # Absolute difference
    abs_val <- abs(inside_n - total_n)
    
    mixing <- c(mixing, mix_param)
    degrees <- c(degrees, deg)
    weighted <- c(weighted, weighted_val)
    abs_diff <- c(abs_diff, abs_val)
  }
  
  # Weighted mixing (normalized by total degree)
  total_degree <- sum(degrees)
  weight_norm <- if (total_degree == 0) rep(0, length(degrees)) else degrees / total_degree
  weighted_mixing <- mixing * weight_norm
  
  data.frame(
    gene = genes_in_path,
    mixing = mixing,
    degree = degrees,
    weighted = weighted,
    weighted_norm = weighted_mixing,
    abs_diff = abs_diff,
    pathway = path_id,
    stringsAsFactors = FALSE
  )
}

# Compute for all pathways
mixing_list <- lapply(paths, compute_mixing)
names(mixing_list) <- paths

# --- 3. Summaries per pathway ---
# Mean unweighted mixing
mean_mix <- sapply(mixing_list, function(df) mean(df$mixing, na.rm = TRUE))
# Total links inside
total_links <- sapply(mixing_list, function(df) sum(df$degree))
# Total outward mixing (sum of weighted)
total_out <- sapply(mixing_list, function(df) sum(df$weighted, na.rm = TRUE))
# Fraction of intralinks (1 - outward mixing / total links)
frac_intra <- 1 - (total_out / total_links)

# Load functional performance (FPR) from supplementary
fps <- read.delim(benchmark_file("Supplementary_data", "Supplementary_Data_3"), header = TRUE)
fps_sub <- fps[match(paths, fps$pathway), ]

# Combine into data frame
results <- data.frame(
  pathway = paths,
  frac_intra = frac_intra,
  FPR = fps_sub$FPR,
  total_links = total_links,
  stringsAsFactors = FALSE
)

# Spearman correlation
corr_spear <- cor(results$frac_intra, results$FPR, method = "spearman")
message("Spearman correlation: ", round(corr_spear, 2))

# --- 4. Scatter plot ---
results$size <- results$total_links / 1e5
p <- ggplot(results, aes(x = frac_intra, y = FPR, size = size)) +
  geom_point() +
  scale_size(name = "Total links", breaks = c(0.15, 0.6, 1.2, 3.5)) +
  scale_x_continuous(
    name = "Fraction of intralinks",
    trans = "log",
    breaks = trans_breaks("log", function(x) exp(x)),
    labels = trans_format("log", math_format(e^.x))
  ) +
  ylab("FPR, BinoX enrichment") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Add panel label
library(cowplot)
p <- ggdraw(p) + draw_plot_label("B", size = 14)

# Save plot
tiff(output_file("fraction_intralinks.tiff"), width = 5.5, height = 4.5, units = "in", res = 400)
print(p)
dev.off()
