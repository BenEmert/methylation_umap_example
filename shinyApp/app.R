library(shiny)
library(bslib)
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(uwot)        
library(Rtsne)       
library(matrixStats)
library(shinycssloaders)
library(shinyWidgets)
getwd()
# ==============================================================================
# Global Setup & Data Loading
# ==============================================================================

data_dir <- "../data/analyzed"
file_betas <- file.path(data_dir, "GSE140686_Final_Betas_corrected.rds")
file_mvals <- file.path(data_dir, "GSE140686_Final_Mvals.rds")
file_meta  <- file.path(data_dir, "metadata_complete.csv")

betas_mat <- NULL
mvals_mat <- NULL
meta_df   <- NULL
global_vars_beta <- NULL
global_vars_mval <- NULL

load_data <- function() {
  if (file.exists(file_betas) && file.exists(file_mvals) && file.exists(file_meta)) {
    message("Loading data...")
    
    betas_mat <<- readRDS(file_betas)
    mvals_mat <<- readRDS(file_mvals)
    meta_df   <<- fread(file_meta, data.table = FALSE)
    
    # Ensure IDs match
    common_ids <- intersect(rownames(betas_mat), meta_df$ID)
    betas_mat <<- betas_mat[common_ids, ]
    mvals_mat <<- mvals_mat[common_ids, ]
    meta_df   <<- meta_df %>% filter(ID %in% common_ids)
    rownames(meta_df) <<- meta_df$ID
    
    message("Calculating global variance...")
    global_vars_beta <<- colVars(betas_mat)
    names(global_vars_beta) <<- colnames(betas_mat)
    global_vars_beta <<- names(sort(global_vars_beta, decreasing = TRUE))
    
    global_vars_mval <<- colVars(mvals_mat)
    names(global_vars_mval) <<- colnames(mvals_mat)
    global_vars_mval <<- names(sort(global_vars_mval, decreasing = TRUE))
    
    message("Data loaded successfully.")
    return(TRUE)
  } else {
    warning("Data files not found. Please ensure pipeline has run.")
    return(FALSE)
  }
}

data_loaded <- load_data()

unique_diagnoses <- if(data_loaded) sort(unique(meta_df$Diagnosis)) else c()

# ==============================================================================
# UI Definition
# ==============================================================================

ui <- page_sidebar(
  theme = bs_theme(version = 5, bootswatch = "zephyr"),
  title = "Methylation UMAP/t-SNE Explorer",
  
  sidebar = sidebar(
    width = 400,
    h4("Configuration"),
    
    # Target Plot Selection
    radioButtons("active_plot", "Target Plot:", 
                 choices = c("Plot A (Left)", "Plot B (Right)"), 
                 inline = TRUE),
    hr(),
    
    # Algorithm Selection
    selectInput("algo_type", "Algorithm", choices = c("UMAP", "t-SNE")),
    
    # Input Parameters
    selectInput("preset", "Presets", 
                choices = c("Custom", "Fast Preview", "High Detail", "Standard Analysis"), 
                selected = "Custom"),
    
    radioButtons("norm_method", "Normalization", 
                 choices = c("M-values", "Beta values"), inline = TRUE),
    
    sliderInput("top_n", "Number of Variable Features (CpGs)", 
                min = 0, max = 100000, value = 5000, step = 500),
    
    # Dynamic UI for UMAP vs t-SNE parameters
    uiOutput("algo_params"),
    
    accordion(
      accordion_panel(
        "Advanced Settings",
        sliderInput("pca_n", "PCA Components (Pre-processing)", 
                    min = 10, max = 100, value = 50),
        selectInput("metric", "Distance Metric", 
                    choices = c("euclidean", "cosine", "correlation", "manhattan")),
        radioButtons("var_mode", "Variance Calculation", 
                     choices = c("Subset variance", "Global variance")),
        actionButton("randomize", "Random Seed"),
        verbatimTextOutput("seed_display", placeholder = TRUE)
      ),
      open = FALSE
    ),
    
    hr(),
    selectizeInput("filter_dx", "Filter by Diagnosis", 
                   choices = c("All", unique_diagnoses), 
                   selected = "All", multiple = TRUE),
    textOutput("sample_count"),
    
    hr(),
    radioButtons("color_by", "Color Mapping", 
                 choices = c("Diagnosis (Dx)", "WHO Category"), inline = TRUE),
    
    hr(),
    actionButton("run", "RUN ANALYSIS", class = "btn-primary btn-lg w-100"),
    div(style = "margin-top: 10px;", textOutput("status_text")),
    textOutput("compute_time")
  ),
  
  # Main Layout
  layout_columns(
    col_widths = c(6, 6),
    card(
      card_header("Plot A"),
      withSpinner(plotlyOutput("plot_a", height = "600px"))
    ),
    card(
      card_header("Plot B"),
      withSpinner(plotlyOutput("plot_b", height = "600px"))
    )
  ),
  
  # Legend Section
  card(
    card_header("Legend"),
    uiOutput("dynamic_legend")
  )
)

# ==============================================================================
# Server Logic
# ==============================================================================

server <- function(input, output, session) {
  
  create_ggplot_interactive <- function(df, x_col, y_col, color_col) {
    if (is.null(df)) return(NULL)
  
    p <- ggplot(df, aes(
        x = .data[[x_col]], 
        y = .data[[y_col]], 
        color = .data[[color_col]],
        text = paste0(
            "<b>ID:</b> ", ID, "<br>",
            "<b>Dx:</b> ", Diagnosis, "<br>",
            "<b>WHO:</b> ", WHO_differentiation
        )
        )) +
        geom_point(size = 2, alpha = 0.8) + 
        scale_color_identity() +
        theme_classic() +
        theme(
            axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none" # Hide ggplot legend (we use the HTML one)
            )
    

    ggplotly(p, tooltip = "text") %>%
        config(displayModeBar = FALSE) # Optional: Hides the floating toolbar
    }
  
  # Reactive State
  vals <- reactiveValues(
    plot_a_data = NULL,
    plot_b_data = NULL,
    seed = 0
  )
  
  # Random Seed Logic
  observeEvent(input$randomize, {
    vals$seed <- sample.int(10000, 1)
  })
  
  output$seed_display <- renderText({
    paste("Seed:", vals$seed)
  })
  
  # Dynamic Parameters based on Algorithm
  output$algo_params <- renderUI({
    if (input$algo_type == "UMAP") {
      tagList(
        sliderInput("n_neighbors", "Neighbors (n_neighbors)", 
                    min = 0, max = 100, value = 15),
        sliderInput("min_dist", "Minimum Distance", 
                    min = 0.0, max = 1.0, value = 0.1, step = 0.05)
      )
    } else {
      # t-SNE Parameters
      tagList(
        sliderInput("perplexity", "Perplexity", 
                    min = 2, max = 100, value = 30),
        numericInput("tsne_iter", "Max Iterations", value = 1000, min = 250, max = 5000)
      )
    }
  })
  
    observeEvent(input$top_n, {
    # If user drags to 0, snap back to 500
    if (input$top_n < 500) {
      updateSliderInput(session, "top_n", value = 500)
      showNotification("Minimum probes set to 500.", type = "message", duration = 2)
    }
  })

  observeEvent(input$n_neighbors, {
    req(input$algo_type == "UMAP") 
    # If user drags to < 2, snap back to 2
    if (input$n_neighbors < 2) {
      updateSliderInput(session, "n_neighbors", value = 2)
      showNotification("Minimum neighbors set to 2.", type = "message", duration = 2)
    }
  })

  # Presets Logic
  observeEvent(input$preset, {
    req(input$preset != "Custom")
    
    if (input$preset == "Fast Preview") {
      updateSliderInput(session, "top_n", value = 1000)
      if (input$algo_type == "UMAP") {
        updateSliderInput(session, "n_neighbors", value = 15)
        updateSliderInput(session, "min_dist", value = 0.5)
      } else {
        updateSliderInput(session, "perplexity", value = 30)
      }
    } else if (input$preset == "High Detail") {
      updateSliderInput(session, "top_n", value = 10000)
      if (input$algo_type == "UMAP") {
        updateSliderInput(session, "n_neighbors", value = 30)
        updateSliderInput(session, "min_dist", value = 0.1)
      } else {
        updateSliderInput(session, "perplexity", value = 50)
      }
    } else if (input$preset == "Standard Analysis") {
      updateSliderInput(session, "top_n", value = 5000)
      if (input$algo_type == "UMAP") {
        updateSliderInput(session, "n_neighbors", value = 15)
        updateSliderInput(session, "min_dist", value = 0.1)
      } else {
        updateSliderInput(session, "perplexity", value = 30)
      }
    }
  })
  
  # Sample Filtering
  get_selected_ids <- reactive({
    req(data_loaded)
    if ("All" %in% input$filter_dx) {
      return(rownames(meta_df))
    } else {
      return(meta_df$ID[meta_df$Diagnosis %in% input$filter_dx])
    }
  })
  
  output$sample_count <- renderText({
    n <- length(get_selected_ids())
    paste(n, "samples selected")
  })
  
  # MAIN COMPUTATION
  observeEvent(input$run, {
    req(data_loaded)
    
    output$status_text <- renderText("Computing...")
    start_time <- Sys.time()
    
    sel_ids <- get_selected_ids()
    if (length(sel_ids) < 3) {
      showNotification("Please select at least 3 samples.", type = "error")
      output$status_text <- renderText("Error: Too few samples")
      return()
    }
    
    # 1. Select Data Matrix
    if (input$norm_method == "M-values") {
      full_mat <- mvals_mat
      global_order <- global_vars_mval
    } else {
      full_mat <- betas_mat
      global_order <- global_vars_beta
    }
    
    # 2. Subset Samples
    X <- full_mat[sel_ids, , drop = FALSE]
    
    # 3. Select Features (Variance)
    top_n <- min(input$top_n, ncol(X))
    
    if (input$var_mode == "Global variance") {
      # Use pre-calculated order
      features <- global_order[1:top_n]
    } else {
      # Calculate variance on subset
      vars <- colVars(X)
      # Get indices of top variance
      top_idx <- order(vars, decreasing = TRUE)[1:top_n]
      features <- colnames(X)[top_idx]
    }
    
    X_sub <- X[, features]
    
    target_plot_name <- ifelse(input$active_plot == "Plot A (Left)", "Plot A", "Plot B")
    
    algo_detail <- if(input$algo_type == "UMAP") {
      paste0(", n_neighbors: ", n_neigh_safe, ", min_dist: ", input$min_dist)
    } else {
      paste0(", perp: ", input$perplexity, ", iter: ", input$tsne_iter)
    }

    title_str <- paste0(
      target_plot_name,
      ", N probes: ", top_n,
      ", ", input$norm_method,
      ", n_PCs: ", input$pca_n,
      algo_detail,
      ", seed: ", vals$seed
    )
    
    # Run Algorithm
    set.seed(vals$seed)
    
    res_coords <- tryCatch({
      if (input$algo_type == "UMAP") {
        # UMAP using 'uwot'
        umap_res <- uwot::umap(
          X_sub,
          n_neighbors = min(input$n_neighbors, nrow(X_sub) - 1),
          min_dist = input$min_dist,
          n_components = 2,
          metric = input$metric,
          pca = min(input$pca_n, min(dim(X_sub))),
          ret_model = FALSE,
          n_threads = 1 
        )
        colnames(umap_res) <- c("UMAP1", "UMAP2")
        as.data.frame(umap_res)
        
      } else {
        # t-SNE using 'Rtsne'
        tsne_res <- Rtsne(
          X_sub,
          dims = 2,
          perplexity = min(input$perplexity, floor((nrow(X_sub) - 1) / 3)),
          theta = 0.5,
          pca = TRUE,
          initial_dims = min(input$pca_n, min(dim(X_sub))),
          max_iter = input$tsne_iter,
          check_duplicates = FALSE
        )
        res_df <- as.data.frame(tsne_res$Y)
        colnames(res_df) <- c("tSNE1", "tSNE2")
        res_df
      }
    }, error = function(e) {
      showNotification(paste("Algorithm Error:", e$message), type = "error")
      return(NULL)
    })
    
    if (is.null(res_coords)) {
      output$status_text <- renderText("Computation Failed")
      return()
    }
    
    # Prepare Plot Data
    res_coords$ID <- rownames(X_sub)
    plot_df <- left_join(res_coords, meta_df, by = "ID")
    
    # Store result in the active slot
    if (input$active_plot == "Plot A (Left)") {
      vals$plot_a_data <- plot_df
      vals$plot_a_title <- title_str
    } else {
      vals$plot_b_data <- plot_df
      vals$plot_b_title <- title_str
    }
    
    end_time <- Sys.time()
    elapsed <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
    output$compute_time <- renderText(paste("Time:", elapsed, "s"))
    output$status_text <- renderText("✅ Complete")
  })
  
  get_coord_cols <- function(df) {
    if ("UMAP1" %in% names(df)) return(c("UMAP1", "UMAP2"))
    if ("tSNE1" %in% names(df)) return(c("tSNE1", "tSNE2"))
    return(c("Dim1", "Dim2")) # Fallback
  }

  output$plot_a <- renderPlotly({
    df <- vals$plot_a_data
    if (is.null(df)) return(plotly_empty() %>% layout(title = "Waiting for run..."))
    
    coords <- get_coord_cols(df)
    c_col <- if (input$color_by == "Diagnosis (Dx)") "Color_dx" else "Color_WHO"
    
    create_ggplot_interactive(df, coords[1], coords[2], c_col)
  })
  
  output$plot_b <- renderPlotly({
    df <- vals$plot_b_data
    if (is.null(df)) return(plotly_empty() %>% layout(title = "Waiting for run..."))
    
    coords <- get_coord_cols(df)
    c_col <- if (input$color_by == "Diagnosis (Dx)") "Color_dx" else "Color_WHO"
    
    create_ggplot_interactive(df, coords[1], coords[2], c_col)
  })
  
  # STATIC HTML LEGEND
  output$dynamic_legend <- renderUI({
    # Combine data from both plots to show comprehensive legend
    df_list <- list(vals$plot_a_data, vals$plot_b_data)
    combined_df <- bind_rows(df_list)
    
    if (nrow(combined_df) == 0) return(div("Run analysis to see legend."))
    
    # Determine which columns to show
    # We always show Diagnosis vs Color_dx and WHO vs Color_WHO
    
    legend_html <- tagList(
      div(style = "display: flex; flex-wrap: wrap; gap: 20px;",
          # Legend Block 1: Diagnosis
          div(
            h6("Diagnosis"),
            div(style = "display: flex; flex-wrap: wrap; gap: 5px; max-height: 200px; overflow-y: auto;",
                lapply(unique(combined_df$Dx), function(dx) {
                  # Find color for this Dx
                  col <- combined_df$Color_dx[combined_df$Dx == dx][1]
                  full_name <- combined_df$Diagnosis[combined_df$Dx == dx][1]
                  div(style = "display: flex; align-items: center; margin-right: 10px;",
                      span(style = paste0("width: 12px; height: 12px; background-color:", col, 
                                          "; display: inline-block; margin-right: 5px; border: 1px solid #ccc;")),
                      span(style = "font-size: 0.8em;", dx)
                  )
                })
            )
          ),
          # Legend Block 2: WHO (Optional, nice to have)
          div(
            h6("WHO Category"),
            div(style = "display: flex; flex-wrap: wrap; gap: 5px;",
                lapply(unique(combined_df$WHO_differentiation), function(who) {
                  # Find color is trickier if color is not strictly mapped to WHO in data
                  # But usually WHO categories don't have direct color columns in your pipeline
                  # If you don't have Color_WHO in your data, remove this block.
                  # Assuming Color_WHO DOES NOT exist based on your previous scripts, 
                  # I will comment this out to prevent errors, unless you added it.
                  NULL 
                })
            )
          )
      )
    )
    legend_html
  })
}

# Run the application 
shinyApp(ui = ui, server = server)