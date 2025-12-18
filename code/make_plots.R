# make_plots.R

suppressPackageStartupMessages({
  library(data.table)    
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(uwot)          
  library(Rtsne)
  library(matrixStats)   
  library(cowplot)       
  library(patchwork)
})

project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
processed_dir <- file.path(data_dir, "analyzed")
plot_dir <- file.path(project_dir, "plots")

if (!dir.exists(plot_dir)) dir.create(plot_dir)

input_data <- file.path(processed_dir, "GSE140686_Final_Mvals.rds")
input_meta <- file.path(processed_dir, "metadata_complete.csv")

run_umap <- function(meth_df, selected_ids, center=FALSE, top_var_n=10000, 
                     run_pca=TRUE, n_pcs=50, n_neighbors=15, 
                     min_dist=0.1, n_components=2, metric="euclidean", 
                     seed=0) {
  
  valid_ids <- intersect(selected_ids, rownames(meth_df))
  X_subset <- meth_df[valid_ids, , drop = FALSE]
  
  vars <- colVars(as.matrix(X_subset))
  
  # Get top variable probes
  top_indices <- order(vars, decreasing = TRUE)[1:min(top_var_n, length(vars))]
  X_analysis <- X_subset[,top_indices]
    
  # Center Data 
  if (center) {
    X_analysis <- scale(X_analysis, center = TRUE, scale = FALSE)
  }
  
  # PCA
  if (run_pca) {
    pca_res <- prcomp(X_analysis, center = TRUE, scale. = FALSE, rank. = n_pcs)
    X_input <- pca_res$x
  } else {
    X_input <- X_analysis
  }
  
  # UMAP
  umap_res <- uwot::umap(
    X_input,
    n_neighbors = n_neighbors,
    n_components = n_components,
    min_dist = min_dist,
    metric = metric,
    ret_model = FALSE,
    seed = seed
  )
  
  res_df <- as.data.frame(umap_res)
  colnames(res_df) <- paste0("UMAP", 1:n_components)
  rownames(res_df) <- rownames(X_analysis)
  
  return(res_df)
}

run_tsne <- function(meth_df, selected_ids, top_var_n=10000, 
                     n_pcs=50, perplexity=30, n_components=2, seed=0) {
  
  valid_ids <- intersect(selected_ids, rownames(meth_df))
  X_subset <- meth_df[valid_ids, , drop = FALSE]
  
  vars <- colVars(as.matrix(X_subset))
  top_indices <- order(vars, decreasing = TRUE)[1:min(top_var_n, length(vars))]
  X_analysis <- X_subset[,top_indices]
  
  set.seed(seed)
  tsne_res <- Rtsne(
    X_analysis, 
    dims = n_components, 
    perplexity = perplexity,
    initial_dims = n_pcs,
    pca = TRUE,             
    check_duplicates = FALSE,
    verbose = FALSE
  )
  
  res_df <- as.data.frame(tsne_res$Y)
  colnames(res_df) <- paste0("tSNE", 1:n_components)
  rownames(res_df) <- rownames(X_analysis)
  
  return(res_df)
}

plot_umap_gg <- function(df, x_col="UMAP1", y_col="UMAP2", color_col="Color_dx",
                        label_col = "Dx", show_legend=FALSE,
                        title_text = NULL, facet_col=NULL, 
                        xlim=NULL, ylim=NULL, rect_coords=NULL,
                        point_size=1.0) {
  
  if (!is.null(label_col)) {
    unique_map <- unique(df[, c(label_col, color_col)])
    
    if (any(duplicated(unique_map[[label_col]]))) {
      
      dupes <- unique_map[[label_col]][duplicated(unique_map[[label_col]])]
      stop(paste0("Ambiguous Legend Mapping! The label '", dupes[1], 
                  "' is assigned to multiple different colors in column '", color_col, "'. ",
                  "Please ensure each label corresponds to exactly one color."))
    }
    
    my_palette <- setNames(unique_map[[color_col]], unique_map[[label_col]])
    
    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[label_col]])) +
      scale_color_manual(values = my_palette, 
                         guide = if (show_legend) "legend" else "none")
      
  } else {
    p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_col]])) +
      scale_color_identity(guide = if (show_legend) "legend" else "none")
  }

  p <- p +
    geom_point(size = point_size, alpha = 1) +
    theme_classic() +
    coord_fixed(ratio = 1, xlim = xlim, ylim = ylim) +
    labs(title = title_text) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)), nrow = 1, scales = "free_y")
  }
  
  if (!is.null(rect_coords)) {
    p <- p + geom_rect(
      xmin = rect_coords$xmin, xmax = rect_coords$xmax,
      ymin = rect_coords$ymin, ymax = rect_coords$ymax,
      fill = NA, color = "red", linetype = "dashed", linewidth = 0.8,
      inherit.aes = FALSE
    )
  }
  
  return(p)
}

message("=== Generating Plots ===")

# Load Data
# ------------------------------------------------------------------------------
if (!file.exists(input_data)) stop(paste("Data not found at:", input_data))

message("   [LOADING] Methylation data (M-values)...")
meth_df <- readRDS(input_data)

message(sprintf("   Loaded dimensions: %d Samples x %d Probes", nrow(meth_df), ncol(meth_df)))

message("   [LOADING] Metadata...")
meta <- fread(input_meta, data.table = FALSE)

if (!"ID" %in% colnames(meta)) stop("Metadata must contain an 'ID' column.")
rownames(meta) <- meta$ID

common_ids <- intersect(rownames(meth_df), rownames(meta))
message(sprintf("   Matched %d samples between data and metadata.", length(common_ids)))

if (length(common_ids) == 0) stop("No matching IDs found between data and metadata.")

meth_df <- meth_df[common_ids, ]
meta <- meta[common_ids, ]

# Save Legend
# ------------------------------------------------------------------------------
message("   Plotting Legend...")

# legend_key <- meta %>%
#   select(Dx, Color_dx) %>%
#   distinct() %>%
#   arrange(Dx)

# dummy_plot <- ggplot(legend_key, aes(x = 1, y = 1, color = Color_dx)) +
#   geom_point(alpha = 0) + 
#   scale_color_identity(
#     breaks = legend_key$Color_dx,     
#     labels = legend_key$Dx,            
#     guide = guide_legend(
#       ncol = 2,                        
#       override.aes = list(
#         alpha = 1,                     
#         size = 5,                      
#         shape = 15                     
#       )
#     )
#   ) +
#   theme_void() +
#   theme(
#     legend.position = "right",
#     legend.text = element_text(size = 10),
#     legend.key = element_blank()      
#   )

# legend_only <- get_legend(dummy_plot)
# n_rows <- ceiling(nrow(legend_key) / 2)
# calc_height <- max(4, n_rows * 0.25)
# print(dummy_plot)
# ggsave(
#   file.path(plot_dir, "dx_legend.pdf"), 
#   plot = ggdraw(legend_only), 
#   width = 6,      
#   height = calc_height
# )

# Random State Comparison
# ------------------------------------------------------------------------------
message("   [PLOTTING] UMAP Random State Comparison...")
emb_list <- list()
plot_list <- list()

for (s in c(0, 1)) {
  tmp <- run_umap(meth_df, common_ids, seed=s)
  tmp$random_state <- paste0("random state = ", s)
  tmp <- tmp %>% rownames_to_column("ID")
  tmp <- left_join(tmp, meta, by="ID")
  
  emb_list[[length(emb_list)+1]] <- tmp
  
  p <- plot_umap_gg(tmp, title_text = paste0("random state = ", s))
  plot_list[[length(plot_list)+1]] <- p
}

rng_plots <- wrap_plots(plot_list, nrow = 1)
ggsave(file.path(plot_dir, "umap_GPL13534_M_random_state_comparison_R.pdf"), rng_plots, width=14, height=7)

combined <- bind_rows(emb_list)

dat0 <- combined %>% filter(random_state == "random state = 0")
p_zoom0 <- plot_umap_gg(dat0, xlim=c(0.5, 7), ylim=c(-3.0, 3.0))
ggsave(file.path(plot_dir, "umap_GPL13534_M_rng0_zoom1.pdf"), p_zoom0, width=7, height=7)

p_box0 <- plot_umap_gg(dat0, rect_coords=list(xmin=0.0, xmax=7.5, ymin=-3.5, ymax=3.5))
ggsave(file.path(plot_dir, "umap_GPL13534_M_rng0_box1.pdf"), p_box0, width=7, height=7)

dat1 <- combined %>% filter(random_state == "random state = 1")
p_zoom1 <- plot_umap_gg(dat1, xlim=c(1, 7.2), ylim=c(-3.5, 2.5))
ggsave(file.path(plot_dir, "umap_GPL13534_M_rng1_zoom1.pdf"), p_zoom1, width=7, height=7)

p_box1 <- plot_umap_gg(dat1, rect_coords=list(xmin=0.5, xmax=7.7, ymin=-4.0, ymax=3.0))
ggsave(file.path(plot_dir, "umap_GPL13534_M_rng1_box1.pdf"), p_box1, width=7, height=7)

# Neighbors Comparison
# ------------------------------------------------------------------------------
# message("   [PLOTTING] UMAP Neighbors Comparison...")
# plot_list <- list()

# for (nn in c(5, 15, 50)) {
#   tmp <- run_umap(meth_df, common_ids, n_neighbors=nn, seed=0)

#   tmp$n_neighbors <- paste0("n_neighbors = ", nn)

#   tmp <- tmp %>% rownames_to_column("ID")

#   tmp <- left_join(tmp, meta, by="ID")

#   p <- plot_umap_gg(tmp, title_text = paste0("n_neighbors = ", nn))

#   plot_list[[length(plot_list)+1]] <- p
# }

# nn_plots <- wrap_plots(plot_list, nrow = 1)
# ggsave(file.path(plot_dir, "umap_GPL13534_M_n_neighbors_comparison.pdf"), nn_plots, width=21, height=7, dpi=300)

# # 5. Perplexity Comparison (tSNE)
# # ------------------------------------------------------------------------------
# message("   [PLOTTING] tSNE Perplexity Comparison...")
# plot_list <- list()

# for (perp in c(5, 15, 50)) {
#   tmp <- run_tsne(meth_df, common_ids, perplexity=perp, seed=0)

#   tmp$perplexity <- paste0("perplexity = ", perp)

#   tmp <- tmp %>% rownames_to_column("ID")

#   tmp <- left_join(tmp, meta, by="ID")

#   p <- plot_umap_gg(tmp, x_col="tSNE1", y_col="tSNE2", title_text = paste0("perplexity = ", perp))

#   plot_list[[length(plot_list)+1]] <- p
# }

# perp_plots <- wrap_plots(plot_list, nrow = 1)
# ggsave(file.path(plot_dir, "tsne_GPL13534_M_perplexity_comparison.pdf"), perp_plots, width=21, height=7, dpi=300)
# # Diagnosis Subsets
# # ------------------------------------------------------------------------------
# message("   [PLOTTING] Diagnosis Subsets...")

# plot_list <- list()

# for (c in unique(meta$WHO_differentiation)) {
#   sub_meta <- meta %>% filter(WHO_differentiation == c | Dx %in% c("CTRL1", "CTRL3"))
#   sub_ids <- sub_meta$ID

#   if (length(sub_ids) > 15){
#     tmp <- run_umap(meth_df, sub_ids, n_pcs=25, seed=0)

#     tmp$WHO_differentiation <- c

#     tmp <- tmp %>% rownames_to_column("ID")

#     tmp <- left_join(tmp, meta, by="ID")

#     p <- plot_umap_gg(tmp, title_text = c, show_legend=TRUE)
#     new_name <- gsub(" ", "_", c) %>% gsub("/", "_", .)
#     ggsave(file.path(plot_dir, sprintf("umap_GPL13534_M_%s_noBlood.pdf", new_name)), p, width=7, height=7)

#     plot_list[[length(plot_list)+1]] <- p

#   }
# }

message("=== Plots Generated in /plots ===")