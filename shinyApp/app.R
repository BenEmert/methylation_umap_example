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
# ==============================================================================
# Global Setup & Data Loading
# ==============================================================================

data_dir <- "../data/analyzed"

file_betas <- file.path(data_dir, "GSE140686_Final_Betas_corrected.rds")
file_mvals <- file.path(data_dir, "GSE140686_Final_Mvals.rds")
file_meta  <- file.path(data_dir, "metadata_complete.csv")
# mvals_mat <<- readRDS('data/analyzed/GSE140686_Final_Mvals.rds')
betas_mat <- NULL
mvals_mat <- NULL
meta_df   <- NULL
global_vars_beta <- NULL
global_vars_mval <- NULL

load_data <- function() {
  if (file.exists(file_betas) && file.exists(file_mvals) && file.exists(file_meta)) {
    message("Loading Metadata...")
    meta_df <<- fread(file_meta, data.table = FALSE)
    rownames(meta_df) <<- meta_df$ID
    
    message("Loading Beta Matrix...")
    betas_mat <<- readRDS(file_betas)
    
    message("Loading M-value Matrix...")
    mvals_mat <<- readRDS(file_mvals)
    
    common_ids <- intersect(rownames(betas_mat), meta_df$ID)
    common_ids <- intersect(common_ids, rownames(mvals_mat))
    
    betas_mat <<- betas_mat[common_ids, ]
    mvals_mat <<- mvals_mat[common_ids, ]
    meta_df   <<- meta_df %>% filter(ID %in% common_ids)
    
    message("Calculating Global Variances...")
    
    beta_v <- colVars(betas_mat)
    names(beta_v) <- colnames(betas_mat)
    global_vars_beta <<- names(sort(beta_v, decreasing = TRUE))
    
    mval_v <- colVars(mvals_mat)
    names(mval_v) <- colnames(mvals_mat)
    global_vars_mval <<- names(sort(mval_v, decreasing = TRUE))
    
    message("Data loaded successfully.")
    return(TRUE)
  } else {
    warning("Data files not found.")
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

  tags$head(
    tags$style(HTML("
      /* Target the selected item in the picker dropdown */
      .bootstrap-select .dropdown-menu li.selected a.dropdown-item {
        background-color: #0d6efd !important; /* Blue background */
        color: white !important;              /* White text */
        font-weight: bold !important;
      }
      /* Ensure the checkmark is white to match the text */
      .bootstrap-select .dropdown-menu li.selected a.dropdown-item span.check-mark {
        color: white !important;
        opacity: 1 !important;
      }
      /* Hover state for better UX */
      .bootstrap-select .dropdown-menu li a.dropdown-item:hover {
        background-color: #e9ecef;
      }
    "))
  ),

  sidebar = sidebar(
    width = 400,
    h4("Configuration"),
    
    radioButtons("active_plot", "Target Plot:", 
                 choices = c("Plot A (Left)", "Plot B (Right)"), 
                 inline = TRUE),
    hr(),
        
    div(class = "d-grid gap-2",
      selectInput("algo_type", NULL, choices = c("UMAP", "t-SNE"), width = "100%"),
      selectInput("preset", NULL, choices = c("Custom", "Fast Preview", "High Detail", "Standard Analysis"), selected = "Custom", width = "100%")
    ),
    
    # --- Collapsible Settings ---
    accordion(
      open = FALSE, 
      accordion_panel(
        "Data Parameters",
        radioButtons("norm_method", "Normalization", choices = c("M-values", "Beta values")),
        sliderInput("top_n", "Var. CpGs", min = 500, max = 50000, value = 5000, step = 500),
        radioButtons("var_mode", "Variance Mode", choices = c("Subset", "Global"))
      ),

    accordion_panel(
        "Algorithm Parameters",
        uiOutput("algo_params"),
        selectInput("metric", "Distance Metric", choices = c("euclidean", "cosine", "correlation", "manhattan")),
        sliderInput("pca_n", "PCA Components", min = 10, max = 100, value = 50),
        
        # --- Editable Seed with Sync Button ---
        div(style = "display: flex; gap: 5px; align-items: flex-end;",
            div(style = "flex-grow: 1;", numericInput("seed_input", "Random Seed", value = 42, min = 1)),
            div(actionButton("randomize", NULL, icon = icon("dice"), class = "btn-secondary mb-3"))
        )
      ),
    ),

    hr(),
    
    pickerInput(
      inputId = "filter_dx",
      label = "Filter Diagnoses", 
      choices = unique_diagnoses,
      selected = unique_diagnoses,
      options = list(
        `actions-box` = TRUE, 
        `live-search` = TRUE,
        `show-tick` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} selected",
        `none-selected-text` = "No samples selected"
      ), 
      multiple = TRUE
    ),

    div(style="font-size: 0.8rem; color: #777; margin-bottom: 10px;", textOutput("sample_count")),
    
    hr(),
    radioButtons("color_by", "Color Mapping", 
                 choices = c("Diagnosis (Dx)", "Methylation Class", "WHO Category")),
    
    hr(),
    actionButton("run", "RUN ANALYSIS", class = "btn-primary btn-lg w-100"),
    div(style = "margin-top: 10px;", textOutput("status_text")),
    textOutput("compute_time"),

    hr(),
    div(class = "d-grid gap-2",
        downloadButton("dl_a", "Download Data (Plot A)", class = "btn-outline-primary btn-sm"),
        downloadButton("dl_b", "Download Data (Plot B)", class = "btn-outline-primary btn-sm")
    )
  ),
  
  # Plot Layout
  layout_columns(
    col_widths = c(6, 6),
    card(
      card_header(uiOutput("header_a", container = span)),
      withSpinner(plotlyOutput("plot_a", height = "600px"))
    ),
    card(
      card_header(uiOutput("header_b", container = span)),
      withSpinner(plotlyOutput("plot_b", height = "600px"))
    )
  ),
  
  # Legend 
  accordion(
    open = FALSE, 
    accordion_panel(
      "View Legend",
      div(
        style = "max-height: 30vh; overflow-y: auto; padding-right: 10px;", 
        uiOutput("dynamic_legend")
        )
      )
    )
)
# ==============================================================================
# Server Logic
# ==============================================================================

server <- function(input, output, session) {
  
  output$selected_dx_tags <- renderUI({
    req(input$filter_dx)
    if (length(input$filter_dx) == length(unique_diagnoses)) {
      return(span(class = "badge bg-secondary", "All Diagnoses Selected"))
    }
    if (length(input$filter_dx) > 10) {
      return(span(class = "badge bg-info", paste(length(input$filter_dx), "Diagnoses Selected")))
    }
    
    tagList(
      lapply(input$filter_dx, function(x) {
        span(class = "badge bg-light text-dark border", style="margin-right: 3px; margin-bottom: 3px;", x)
      })
    )
  })

  create_ggplot_interactive <- function(df, x_col, y_col, color_col) {
    if (is.null(df)) return(NULL)
    if (!color_col %in% colnames(df)) return(NULL)
  
    p <- ggplot(df, aes(
        x = .data[[x_col]], y = .data[[y_col]], color = .data[[color_col]],
        text = paste0(
            "<b>ID:</b> ", ID, "<br>",
            "<b>Dx:</b> ", Diagnosis, "<br>",
            "<b>Methyl:</b> ", Methylation_Class, "<br>",
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
        config(displayModeBar = TRUE, 
             modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d"))
    }
  vals <- reactiveValues(
    plot_a_data = NULL, plot_a_params = NULL, plot_a_title = NULL,
    plot_b_data = NULL, plot_b_params = NULL, plot_b_title = NULL
  )
  
  observeEvent(input$randomize, {
    updateNumericInput(session, "seed_input", value = sample.int(10000, 1))
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
        sliderInput("tsne_iter", "Max Iterations", value = 1000, 
                    min = 100, max = 10000, step = 100)
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
    if (input$n_neighbors < 2) {
      updateSliderInput(session, "n_neighbors", value = 2)
      showNotification("Minimum neighbors set to 2.", type = "message", duration = 2)
    }
  })

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
  
  # If a user changes a parameter, check if it deviates from the preset
  
  check_preset_deviation <- function(param, value) {
    req(input$preset != "Custom")
    
    # Define expected values for each preset
    presets <- list(
      "Fast Preview"      = list(top_n = 1000,  n_neighbors = 15, min_dist = 0.5, perplexity = 30),
      "High Detail"       = list(top_n = 10000, n_neighbors = 30, min_dist = 0.1, perplexity = 50),
      "Standard Analysis" = list(top_n = 5000,  n_neighbors = 15, min_dist = 0.1, perplexity = 30)
    )
    
    expected <- presets[[input$preset]][[param]]
    
    # If the parameter exists in the preset definition, check for deviation
    if (!is.null(expected)) {
      if (abs(value - expected) > 0.001) {
      updateSelectInput(session, "preset", selected = "Custom")
      }
    }
  }

  observeEvent(input$top_n,       { check_preset_deviation("top_n", input$top_n) })
  observeEvent(input$n_neighbors, { check_preset_deviation("n_neighbors", input$n_neighbors) })
  observeEvent(input$min_dist,    { check_preset_deviation("min_dist", input$min_dist) })
  observeEvent(input$perplexity,  { check_preset_deviation("perplexity", input$perplexity) })

  # Sample Filtering
  get_selected_ids <- reactive({
    req(data_loaded)
    if (is.null(input$filter_dx)) return(character(0))
    return(meta_df$ID[meta_df$Diagnosis %in% input$filter_dx])
  })
  
  output$sample_count <- renderText({
    n <- length(get_selected_ids())
    paste(n, "samples selected")
  })

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
    
    withProgress(message = 'Initializing...', value = 0, {
      
      # Prepare Matrix
      incProgress(0.1, detail = "Preparing Matrix...")
      if (input$norm_method == "M-values") {
        X_full <- mvals_mat
        global_order <- global_vars_mval
      } else {
        X_full <- betas_mat
        global_order <- global_vars_beta
      }
      
      # Subset Samples
      incProgress(0.1, detail = "Subsetting...")
      X_sel <- X_full[sel_ids, , drop = FALSE]
      
      # Feature Selection
      incProgress(0.1, detail = "Selecting Features...")
      
      top_n <- input$top_n

      if (input$var_mode == "Global") {
        features <- global_order[1:top_n]
      } else {
        vars <- colVars(X_sel, na.rm=TRUE)
        top_idx <- order(vars, decreasing = TRUE)[1:top_n]
        features <- colnames(X_sel)[top_idx]
      }
      X_final <- X_sel[, features, drop=FALSE]
      
      # --- PARAMETER VALIDATION & WARNINGS ---
      
      # Handle potentially NULL inputs
      in_neighbors <- if(is.null(input$n_neighbors)) 15 else input$n_neighbors
      in_min_dist  <- if(is.null(input$min_dist)) 0.1 else input$min_dist
      in_perp      <- if(is.null(input$perplexity)) 30 else input$perplexity
      in_iter      <- if(is.null(input$tsne_iter)) 1000 else input$tsne_iter
      
      # PCA Check
      max_possible_pca <- min(nrow(X_final), ncol(X_final)) - 1
      
      if (input$pca_n > max_possible_pca) {
        safe_pca_n <- max_possible_pca
        showNotification(
          paste("Warning: PCA components reduced to", safe_pca_n, "due to sample size limits."),
          type = "warning", duration = 6
        )
      } else {
        safe_pca_n <- input$pca_n
      }
      
      # Robust Neighbors Check (UMAP only)
      n_neigh_safe <- in_neighbors
      if (input$algo_type == "UMAP") {
        max_neighbors <- nrow(X_final) - 1
        if (n_neigh_safe > max_neighbors) {
          n_neigh_safe <- max_neighbors
          if (n_neigh_safe < 2) n_neigh_safe <- 2
          showNotification(
            paste("Warning: Neighbors reduced to", n_neigh_safe, "due to sample size limits."),
            type = "warning", duration = 6
          )
        }
      }
      
      # Robust Perplexity Check (t-SNE only)
      safe_perp <- in_perp
      if (input$algo_type == "t-SNE") {
        # Limit is (N - 1) / 3
        max_perp <- floor((nrow(X_final) - 1) / 3)
        if (max_perp < 1) max_perp <- 1
        
        if (safe_perp > max_perp) {
          safe_perp <- max_perp
          showNotification(
            paste("Warning: Perplexity reduced to", safe_perp, "due to sample size limits."),
            type = "warning", duration = 6
          )
        }
      }
      
      # Store current params to attach to the plot later
      current_params <- list(
        algo = input$algo_type,
        norm = input$norm_method,
        var_mode = input$var_mode,
        top_n = top_n,
        distance = input$metric,
        pca = safe_pca_n,
        seed = input$seed_input,
        neigh = if(input$algo_type == "UMAP"){ n_neigh_safe } else { NULL },
        min_dist = if(input$algo_type == "UMAP"){ in_min_dist } else { NULL },
        perp = if(input$algo_type == "t-SNE"){ safe_perp } else { NULL },
        iter = if(input$algo_type == "t-SNE"){ in_iter } else { NULL }
      )
      # --- RUN ALGORITHM ---
      
      target_plot_name <- ifelse(input$active_plot == "Plot A (Left)", "Plot A", "Plot B")
      
      # Format: UMAP | N=5000 | M-values | nn=15 md=0.1 | seed=42
      algo_str <- if(input$algo_type == "UMAP") {
        paste0("nn=", n_neigh_safe, " md=", in_min_dist)
      } else {
        paste0("perp=", safe_perp, " iter=", in_iter)
      }
      
      title_str <- paste(
        input$algo_type,
        paste0("N=", top_n),
        input$norm_method,
        paste0("D=", input$metric),
        algo_str,
        paste0("PCA=", safe_pca_n),
        paste0("seed=", input$seed_input),
        sep = " | "
      )
      
      incProgress(0.3, detail = paste0("Running ", input$algo_type, "..."))
      
      res_coords <- tryCatch({
        
        if (input$algo_type == "UMAP") {
          umap_res <- uwot::umap(
            X_final,
            n_neighbors = n_neigh_safe,
            min_dist = in_min_dist,
            n_components = 2,
            metric = input$metric,
            pca = safe_pca_n,
            ret_model = FALSE,
            n_threads = 1,
            seed = input$seed_input 
          )
          colnames(umap_res) <- c("UMAP1", "UMAP2")
          as.data.frame(umap_res)
          
        } else {
          set.seed(input$seed_input)
          
          tsne_res <- Rtsne(
            X_final,
            dims = 2,
            perplexity = safe_perp,
            theta = 0.5,
            pca = TRUE,
            initial_dims = safe_pca_n,  
            max_iter = in_iter,
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
      
      if (!is.null(res_coords)) {
        res_coords$ID <- rownames(X_final)
        plot_df <- left_join(res_coords, meta_df, by = "ID")
        
        if (input$active_plot == "Plot A (Left)") {
          vals$plot_a_data <- plot_df
          vals$plot_a_title <- title_str
          vals$plot_a_params <- current_params
        } else {
          vals$plot_b_data <- plot_df
          vals$plot_b_title <- title_str
          vals$plot_b_params <- current_params
        }
        
        output$status_text <- renderText("Complete")
        end_time <- Sys.time()
        elapsed <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
        output$compute_time <- renderText(paste("Time:", elapsed, "s"))
      } else {
        output$status_text <- renderText("Failed")
      }
    })
  })
  
  get_coord_cols <- function(df) {
    if ("UMAP1" %in% names(df)) return(c("UMAP1", "UMAP2"))
    if ("tSNE1" %in% names(df)) return(c("tSNE1", "tSNE2"))
    return(c("Dim1", "Dim2")) 
  }

  render_plot <- function(df) {
    if (is.null(df)) return(plotly_empty() %>% layout(title = "Waiting for run..."))
    
    coords <- get_coord_cols(df)
    c_col <- switch(input$color_by,
                    "Diagnosis (Dx)" = "Color_dx",
                    "Methylation Class" = "Color_Methyl",
                    "WHO Category" = "Color_WHO")
    
    # Pass title to layout
    create_ggplot_interactive(df, coords[1], coords[2], c_col)
  }
  
  output$plot_a <- renderPlotly({ render_plot(vals$plot_a_data) })
  output$plot_b <- renderPlotly({ render_plot(vals$plot_b_data) })

  # --- Dynamic Card Headers ---
  output$header_a <- renderUI({
    if (is.null(vals$plot_a_title)) {
      return(strong("Plot A (Left)"))
    } else {
      tagList(
        strong("Plot A"), 
        span(" | ", style = "color: #ccc;"),
        span(vals$plot_a_title, style = "font-size: 0.85rem; font-weight: normal;")
      )
    }
  })

  output$header_b <- renderUI({
    if (is.null(vals$plot_b_title)) {
      return(strong("Plot B (Right)"))
    } else {
      tagList(
        strong("Plot B"), 
        span(" | ", style = "color: #ccc;"),
        span(vals$plot_b_title, style = "font-size: 0.85rem; font-weight: normal;")
      )
    }
  })
  
  output$dynamic_legend <- renderUI({
    df_list <- list(vals$plot_a_data, vals$plot_b_data)
    combined_df <- bind_rows(df_list)
    
    if (nrow(combined_df) == 0) return(div("Run analysis to see legend."))
    
    create_legend_block <- function(title, label_col, color_col) {
      if (!all(c(label_col, color_col) %in% colnames(combined_df))) return(NULL)
      
      # Extract unique pairs of Label + Color
      pairs <- combined_df %>% 
        select(all_of(c(label_col, color_col))) %>% 
        distinct() %>% 
        arrange(.data[[label_col]]) %>%
        na.omit()
      
      div(
        h6(title),
        div(style = "display: flex; flex-wrap: wrap; gap: 8px; max-height: 200px; overflow-y: auto; margin-bottom: 15px;",
            lapply(1:nrow(pairs), function(i) {
              lbl <- pairs[[label_col]][i]
              col <- pairs[[color_col]][i]
              div(style = "display: flex; align-items: center; margin-right: 10px;",
                  span(style = paste0("width: 12px; height: 12px; background-color:", col, 
                                      "; display: inline-block; margin-right: 5px; border: 1px solid #ccc; border-radius: 2px;")),
                  span(style = "font-size: 0.8em;", lbl)
              )
            })
        )
      )
    }
    
    tagList(
      div(style = "display: flex; flex-direction: column;",
          create_legend_block("Diagnosis", "Diagnosis", "Color_dx"),
          create_legend_block("Methylation Class", "Methylation_Class", "Color_Methyl"),
          create_legend_block("WHO Differentiation", "WHO_differentiation", "Color_WHO")
      )
    )
  })

  generate_download <- function(target) {
    downloadHandler(
      filename = function() {
        paste0("methylation_embedding_", tolower(target), "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        # Select data based on target
        if (target == "A") {
          df <- vals$plot_a_data
          p  <- vals$plot_a_params
        } else {
          df <- vals$plot_b_data
          p  <- vals$plot_b_params
        }
        
        req(df)
        
        # Create Metadata Header Lines

        header_lines <- c(
          paste0("# Algorithm: ", p$algo),
          paste0("# Normalization: ", p$norm),
          paste0("# Top Variable Features: ", p$top_n),
          paste0("# Distance metric: ", p$distance),
          paste0("# Variance mode: ", p$var_mode),
          paste0("# PCA Components: ", p$pca),
          paste0("# Random Seed: ", p$seed)
        )
        
        if (p$algo == "UMAP") {
          header_lines <- c(header_lines, paste0("# Neighbors: ", p$neigh), paste0("# Min Dist: ", p$min_dist))
        } else {
          header_lines <- c(header_lines, paste0("# Perplexity: ", p$perp), paste0("# Iterations: ", p$iter))
        }
        
        # Write header then append CSV data
        writeLines(header_lines, file)
        suppressWarnings(write.table(df, file, append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE))
      }
    )
  }
  
  output$dl_a <- generate_download("A")
  output$dl_b <- generate_download("B")
}

shinyApp(ui = ui, server = server)