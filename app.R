# app.R — FPCA visualiser
# -----------------------------------
# Expected fpca_bundle (.rds) fields:
# scores (n x npc matrix), functions (funData), values (npc),
# fit (optional), meta (optional), file_order (length n),
# map_col (key column name in meta; default "File")

library(shiny)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(plotly)
library(DT)
library(scales)
library(MFPCA)
library(funData)

# ---- helpers ----

# Shared color palette for consistency between plotly scatter and ggplot reconstruction
# Using plotly's default color sequence
get_color_palette <- function(n) {
  # Plotly's default qualitative color sequence (D3 category10)
  plotly_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                     "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  if (n <= length(plotly_colors)) {
    return(plotly_colors[seq_len(n)])
  }
  # If more colors needed, use colorRampPalette
  colorRampPalette(plotly_colors)(n)
}

get_mu_vec <- function(fpca_bundle) {
  stopifnot(!is.null(fpca_bundle$mu), methods::is(fpca_bundle$mu, "funData"))
  as.numeric(fpca_bundle$mu@X)
}

make_pc_table <- function(values) {
  varprop <- values / sum(values)
  tibble(
    PC = paste0("PC", seq_along(values)),
    eigenvalue = values,
    var_prop = varprop,
    var_prop_pct = 100 * varprop,
    cum_var_prop = cumsum(varprop),
    cum_var_prop_pct = 100 * cumsum(varprop)
  )
}

fun_pc_long <- function(functions_fd, pcs_to_plot) {
  # functions_fd@X: npc x T
  phi_mat <- functions_fd@X
  time_vals <- functions_fd@argvals[[1]]
  
  pcs_to_plot <- sort(unique(pcs_to_plot))
  pcs_to_plot <- pcs_to_plot[pcs_to_plot >= 1 & pcs_to_plot <= nrow(phi_mat)]
  
  plot_df <- tibble(Time = time_vals)
  for (pc in pcs_to_plot) {
    plot_df[[paste0("PC", pc)]] <- as.numeric(phi_mat[pc, ])
  }
  
  plot_long <- plot_df |>
    pivot_longer(cols = -Time, names_to = "Curve", values_to = "Value")
  
  plot_long
}

reconstruct_curve <- function(mu_vec, phi_mat, scores_row, K = NULL, pcs = NULL) {
  # If K provided: use first K PCs
  # If pcs provided: use specific PCs
  if (!is.null(K)) {
    use <- seq_len(min(K, length(scores_row)))
  } else if (!is.null(pcs)) {
    use <- pcs
  } else {
    use <- seq_along(scores_row)
  }
  use <- use[use >= 1 & use <= length(scores_row)]
  
  mu_vec + as.numeric(colSums(phi_mat[use, , drop = FALSE] * scores_row[use]))
}

# ---- UI ----

ui <- fluidPage(
  # JavaScript for copy to clipboard and custom click handling
 tags$head(tags$script(HTML("
    Shiny.addCustomMessageHandler('copyToClipboard', function(text) {
      navigator.clipboard.writeText(text).then(function() {
        alert('Copied to clipboard!');
      }).catch(function(err) {
        // Fallback for older browsers
        var textarea = document.createElement('textarea');
        textarea.value = text;
        document.body.appendChild(textarea);
        textarea.select();
        document.execCommand('copy');
        document.body.removeChild(textarea);
        alert('Copied to clipboard!');
      });
    });
  "))),
  titlePanel("FPCA Visualiser"),

  sidebarLayout(
    sidebarPanel(
      fileInput("rds_file", "Load fpca_bundle (.rds)", accept = c(".rds")),
      actionButton("load_example", "Load Example Dataset", width = "100%"),
      tags$hr(),

      # Only show these controls when on the "Scores + Reconstruction" tab
      conditionalPanel(
        condition = "input.main_tabs == 'Scores + Reconstruction'",

        tags$hr(),
        h4("Scores scatter"),
        uiOutput("xy_select_ui"),
        uiOutput("aes_select_ui"),

        tags$hr(),
        h4("Interactive Reconstruction"),

        # Slider-based reconstruction controls
        checkboxInput("show_slider_recon", "Show slider-based reconstruction", value = TRUE),
        p("You can also move the red dot in the figure, but note that the real-time reconstruction is slow online."),
        uiOutput("slider_x_ui"),
        uiOutput("slider_y_ui"),

        tags$hr(),
        h4("Plot Settings"),
        fluidRow(
          column(6, numericInput("scatter_height", "Scatter height (px)", value = 480, min = 200, max = 1000, step = 20)),
          column(6, numericInput("recon_height", "Recon height (px)", value = 320, min = 150, max = 800, step = 20))
        ),
        uiOutput("recon_ylim_ui")
      )
    ),

    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Instruction",
                 tags$div(
                   style = "max-width: 900px; padding: 20px;",
                   
                   h2("Purpose of the tool"),
                   p("This tool provides interactive visualisation of FPCA results, enabling examination of PC score distributions and direct visualisation of reconstructed curves of selected observations. It is particularly useful for interpreting results and identifying outliers."),
                   p("Currently, it supports one-dimensional FPCA computed via the PACE algorithm."),
                   
                   tags$hr(),
                   
                   h2("How to save an FPCA bundle?"),
                   p("This tool requires uploading an FPCA result bundle as an ", tags$code(".rds"), " file. The steps below describe how to prepare and save one."),
                   p("You can also load the example dataset to proceed."),

                   tags$hr(),

                   h3("1. Run FPCA analysis"),
                   tags$pre(
                     style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto; font-size: 12px;",
                     'raw.data <- read.csv("data/combine.csv") %>% dplyr::select("File", "Timepoint", "time", "Speaker", "gender", "Tone", "Syllable", "Intonation", "f0_semi")

curves <- raw.data %>%
  arrange(File, Timepoint)

stopifnot(!anyNA(curves$File), !anyNA(curves$Timepoint), !anyNA(curves$f0_semi))

# Functional data object, only the essential matrix
curvesFun <- long2irregFunData(curves, id = "File", time = "Timepoint", value = "f0_semi") %>% 
  as.funData()

# explicit file order
file_order <- curves %>%
  distinct(File) %>%
  pull(File)

# Compute FPCA
fpca_pace <- PACE(curvesFun)'
                   ),
                   
                   h3("2. Prepare file_order"),
                   p(tags$code("file_order"), " defines the mapping between rows of ", tags$code("fpca_pace$scores"),
                     " and your original curves. It must be in the same order as the data used to build the FPCA input."),
                   tags$pre(
                     style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto; font-size: 12px;",
                     'file_order <- curves %>%
  dplyr::distinct(File) %>%
  dplyr::pull(File)'
                   ),
                   tags$div(
                     style = "background-color: #fff3cd; padding: 10px; border-radius: 5px; margin: 10px 0;",
                     tags$strong("Important: "), "If you filtered or reordered curves before creating the FPCA object, do the same here."
                   ),
                   
                   h3("3. Prepare meta_df (optional)"),
                   p(tags$code("meta_df"), " should:"),
                   tags$ul(
                     tags$li("Have one row per curve."),
                     tags$li("Contain no time-varying columns (no ", tags$code("time"), ", ", tags$code("f0"), ", etc.)."),
                     tags$li("Convert categorical variables to ", tags$code("factor"), " (recommended).")
                   ),
                   tags$pre(
                     style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto; font-size: 12px;",
                     'meta <- curves %>%
  dplyr::select(-any_of(c("Timepoint", "time", "f0_semi"))) %>%
  dplyr::distinct(File, .keep_all = TRUE) %>%
  mutate(
    Speaker    = as.factor(Speaker),
    gender     = as.factor(gender),
    Tone       = as.factor(Tone),
    Syllable   = as.factor(Syllable),
    Intonation = as.factor(Intonation)
  )'
                   ),
                   

                   h3("4. Save the bundle"),
                   
                   h4("The save function"),
                   p("Run the following save function in R. The goal of ", tags$code("save_fpca_pace_bundle()"), " is to:"),
                   tags$ul(
                     tags$li("Collect all ", tags$strong("essential FPCA outputs"), " (scores, eigenfunctions, eigenvalues, mean curve)."),
                     tags$li("Preserve a ", tags$strong("stable mapping"), " between FPCA scores and original curve IDs."),
                     tags$li("Optionally attach ", tags$strong("metadata"), " (speaker, condition, tone, etc.) in a safe and mergeable way."),
                     tags$li("Save everything into ", tags$strong("one portable ", tags$code(".rds"), " object"), ".")
                   ),
                   
                   tags$pre(
                     style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto; font-size: 12px;",
                     'save_fpca_pace_bundle <- function(
  fpca_pace,
  file_order,
  meta_df = NULL,
  map_col = "File",
  out_file = "fpca_pace_bundle.rds"
) {
  # ---- basic checks ----
  stopifnot(
    is.matrix(fpca_pace$scores),
    is.character(file_order),
    length(file_order) == nrow(fpca_pace$scores),
    is.character(map_col),
    length(map_col) == 1
  )

  # ---- FPCA object checks ----
  stopifnot(
    !is.null(fpca_pace$functions),
    methods::is(fpca_pace$functions, "funData"),
    !is.null(fpca_pace$values),
    length(fpca_pace$values) == ncol(fpca_pace$scores),
    !is.null(fpca_pace$mu),
    methods::is(fpca_pace$mu, "funData")
  )

  # ---- optional metadata checks ----
  if (!is.null(meta_df)) {
    stopifnot(is.data.frame(meta_df))
    stopifnot(map_col %in% colnames(meta_df))

    key <- meta_df[[map_col]]

    # enforce exactly one row per key
    stopifnot(nrow(meta_df) == length(unique(key)))

    # ensure metadata covers all keys
    missing_keys <- setdiff(file_order, key)
    if (length(missing_keys) > 0) {
      stop(
        "meta_df is missing key(s) in column \'", map_col, "\': ",
        paste(missing_keys, collapse = ", ")
      )
    }
  }

  # ---- bundle ----
  fpca_bundle <- list(
    # core FPCA outputs
    scores    = fpca_pace$scores,
    functions = fpca_pace$functions,
    values    = fpca_pace$values,
    mu        = fpca_pace$mu,

    # additional outputs
    fit       = fpca_pace$fit,
    estVar    = fpca_pace$estVar,
    npc       = fpca_pace$npc,
    sigma2    = fpca_pace$sigma2,

    # bookkeeping
    file_order = file_order,
    meta       = meta_df,
    map_col    = map_col,
    created_at = Sys.time()
  )

  saveRDS(fpca_bundle, out_file)
  invisible(fpca_bundle)
}'
                   ),

                   p("The saved ", tags$code(".rds"), " file is a named list with (at minimum):"),
                   tags$ul(
                     tags$li(tags$code("scores"), " \u2014 Numeric matrix ", tags$code("[n_observations \u00d7 npc]"), " of FPCA scores."),
                     tags$li(tags$code("functions"), " \u2014 ", tags$code("funData"), " object containing FPCA eigenfunctions (PC dimensions)."),
                     tags$li(tags$code("values"), " \u2014 Numeric vector of eigenvalues (length = ", tags$code("npc"), ")."),
                     tags$li(tags$code("mu"), " \u2014 ", tags$code("funData"), " object representing the ", tags$strong("mean curve"), ". This is required for meaningful reconstruction: ",
                             tags$em("\u0177(t) = \u03bc(t) + \u03a3 s_k \u03c6_k(t)")),
                     tags$li(tags$code("fit, estVar, npc, sigma2"), " \u2014 Additional FPCA outputs saved as-is."),
                     tags$li(tags$code("file_order"), " \u2014 Character vector (length = ", tags$code("n_observations"), ") mapping each row of ", tags$code("scores"), " to a curve ID."),
                     tags$li(tags$code("meta"), " (optional) \u2014 A data frame with one row per curve, containing metadata."),
                     tags$li(tags$code("map_col"), " \u2014 Name of the ID column in ", tags$code("meta"), " used to match ", tags$code("file_order"), " (default: ", tags$code('"File"'), ")."),
                     tags$li(tags$code("created_at"), " \u2014 Timestamp of bundle creation.")
                   ),

                   h4("Running the save function"),
                   tags$pre(
                     style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto; font-size: 12px;",
'save_fpca_pace_bundle(
  fpca_pace  = fpca_pace,
  file_order = file_order,
  meta_df    = meta,
  map_col    = "File",
  out_file   = "data/fpca_pace_bundle.rds"
)'
                   )
                 )
        ),
        tabPanel("Summary",
                 verbatimTextOutput("bundle_info"),
                 plotlyOutput("var_bar", height = 260),
                 DTOutput("pc_table")
        ),
        tabPanel("Eigenfunctions",
                 plotOutput("eigen_plot", height = 320),
                 plotOutput("pc_effect_plot", height = 320)
        ),
        tabPanel("Scores + Reconstruction",
                 uiOutput("scatter_plot_ui"),
                 tags$hr(),
                 fluidRow(
                   column(4,
                     h4("Point Selection"),
                     helpText("Click points on the scatter plot to select them. Click again to deselect."),
                    uiOutput("selection_info"),
                    fluidRow(
                      column(6, actionButton("clear_selection", "Clear All", width = "100%")),
                      column(6, actionButton("copy_selection", "Copy Names", width = "100%"))
                    ),
                    tags$hr(),
                     radioButtons(
                       "recon_mode", "Reconstruction mode",
                       choices = c("Observation (use K PCs)" = "obs",
                                   "Axis-only (use selected x/y PCs)" = "axis"),
                       selected = "axis"
                     ),
                     conditionalPanel(
                       condition = "input.recon_mode == 'obs'",
                       numericInput("K", "K (number of PCs)", value = 5, min = 1, max = 10, step = 1)
                     ),
                   ),
                   column(8,
                     uiOutput("recon_plot_ui")
                   )
                 )
        )
      )
    )
  )
)

# ---- Server ----

server <- function(input, output, session) {

  fpca_bundle <- reactiveVal(NULL)
  selected_rowids <- reactiveVal(integer(0))  # multi-selection
  live_crosshair <- reactiveVal(NULL)  # For real-time drag updates
  drag_end_time <- reactiveVal(0)  # Timestamp when drag ended (to ignore late events)
  plot_zoom <- reactiveVal(NULL)  # Preserve zoom state across re-renders

  observeEvent(input$rds_file, {
    req(input$rds_file)
    obj <- readRDS(input$rds_file$datapath)
    fpca_bundle(obj)
    selected_rowids(integer(0))  # clear selection on new file
    plot_zoom(NULL)  # reset zoom on new file
  })

  observeEvent(input$load_example, {
    path <- "data/example_bundle.rds"
    if (!file.exists(path)) {
      showNotification("Example file not found: data/example_bundle.rds", type = "error")
      return()
    }
    obj <- readRDS(path)
    fpca_bundle(obj)
    selected_rowids(integer(0))
    plot_zoom(NULL)
    showNotification("Example dataset loaded.", type = "message", duration = 4)
  })

  # Clear selection button
  observeEvent(input$clear_selection, {
    selected_rowids(integer(0))
  })

  # Copy selection to clipboard
  observeEvent(input$copy_selection, {
    ids <- selected_rowids()
    if (length(ids) == 0) return()

    df <- pcscores()
    key_name <- map_col()
    selected_keys <- df[[key_name]][df$.rowid %in% ids]
    text_to_copy <- paste(selected_keys, collapse = "\n")

    session$sendCustomMessage("copyToClipboard", text_to_copy)
  })

  # Derived reactive pieces
  bundle_ok <- reactive({
    b <- fpca_bundle()
    if (is.null(b)) return(FALSE)
    
    needed <- c("scores", "functions", "values", "file_order", "mu")
    
    has_needed <- all(needed %in% names(b))
    if (!has_needed) return(FALSE)
    
    if (!is.matrix(b$scores)) return(FALSE)
    if (length(b$file_order) != nrow(b$scores)) return(FALSE)
    if (length(b$values) != ncol(b$scores)) return(FALSE)
    if (!methods::is(b$functions, "funData")) return(FALSE)
    if (!methods::is(b$mu, "funData")) return(FALSE)
    
    TRUE
  })
  
  npc <- reactive({
    req(bundle_ok())
    ncol(fpca_bundle()$scores)
  })
  
  map_col <- reactive({
    b <- fpca_bundle()
    if (is.null(b)) return("File")
    if (!is.null(b$map_col) && is.character(b$map_col) && length(b$map_col) == 1) {
      b$map_col
    } else {
      "File"
    }
  })
  
  pcscores <- reactive({
    req(bundle_ok())
    b <- fpca_bundle()
    n <- nrow(b$scores)
    p <- ncol(b$scores)
    
    key_name <- map_col()
    
    df <- as_tibble(b$scores)
    names(df) <- paste0("PC", seq_len(p))
    df <- df %>%
      mutate(.rowid = seq_len(n)) %>%
      mutate(!!key_name := b$file_order)
    
    if (!is.null(b$meta) && is.data.frame(b$meta) && (key_name %in% names(b$meta))) {
      # ensure distinct keys in meta
      meta2 <- b$meta %>% distinct(.data[[key_name]], .keep_all = TRUE)
      df <- df %>% left_join(meta2, by = key_name)
    }
    
    df
  })
  
  # UI: PC selectors
  output$pc_select_ui <- renderUI({
    req(bundle_ok())
    p <- npc()
    selectInput(
      "pcs_to_plot", "Eigenfunctions: PCs to plot",
      choices = as.list(setNames(seq_len(p), paste0("PC", seq_len(p)))),
      selected = as.list(setNames(seq_len(min(3, p)), paste0("PC", seq_len(min(3, p))))),
      multiple = TRUE
    )
  })
  
  output$xy_select_ui <- renderUI({
    req(bundle_ok())
    p <- npc()
    tagList(
      selectInput("x_pc", "X axis", choices = paste0("PC", 1:p), selected = "PC1"),
      selectInput("y_pc", "Y axis", choices = paste0("PC", 1:p), selected = if (p >= 2) "PC2" else "PC1")
    )
  })
  
  output$aes_select_ui <- renderUI({
    req(bundle_ok())
    b <- fpca_bundle()
    key_name <- map_col()

    meta_cols <- character(0)
    if (!is.null(b$meta) && is.data.frame(b$meta)) {
      meta_cols <- setdiff(names(b$meta), key_name)
    }

    choices <- c("None" = "", meta_cols)

    tagList(
      selectInput("color_by", "Color by (meta column)", choices = choices, selected = ""),
      selectInput("shape_by", "Shape by (meta column)", choices = choices, selected = "")
    )
  })

  # Dynamic slider for X axis PC value
  output$slider_x_ui <- renderUI({
    req(bundle_ok(), input$x_pc)
    df <- pcscores()
    xcol <- input$x_pc
    if (!(xcol %in% names(df))) return(NULL)

    xvals <- df[[xcol]]
    xmin <- floor(min(xvals, na.rm = TRUE) * 10) / 10
    xmax <- ceiling(max(xvals, na.rm = TRUE) * 10) / 10
    xmid <- round(mean(xvals, na.rm = TRUE), 2)

    sliderInput(
      "slider_x", paste0(xcol, " value"),
      min = xmin, max = xmax, value = xmid, step = (xmax - xmin) / 100
    )
  })

  # Dynamic slider for Y axis PC value
  output$slider_y_ui <- renderUI({
    req(bundle_ok(), input$y_pc)
    df <- pcscores()
    ycol <- input$y_pc
    if (!(ycol %in% names(df))) return(NULL)

    yvals <- df[[ycol]]
    ymin <- floor(min(yvals, na.rm = TRUE) * 10) / 10
    ymax <- ceiling(max(yvals, na.rm = TRUE) * 10) / 10
    ymid <- round(mean(yvals, na.rm = TRUE), 2)

    sliderInput(
      "slider_y", paste0(ycol, " value"),
      min = ymin, max = ymax, value = ymid, step = (ymax - ymin) / 100
    )
  })

  # Display selection info with inline tag/chip layout for deselection
  output$selection_info <- renderUI({
    ids <- selected_rowids()
    if (length(ids) == 0) {
      return(helpText("No points selected."))
    }
    df <- pcscores()
    key_name <- map_col()
    color_by <- input$color_by
    use_color_by <- is.character(color_by) &&
      length(color_by) == 1 &&
      !is.na(color_by) &&
      nzchar(color_by) &&
      (color_by %in% names(df))

    # Get color palette for consistent coloring
    if (isTRUE(use_color_by)) {
      all_color_levels <- unique(as.character(df[[color_by]]))
      all_colors <- get_color_palette(length(all_color_levels))
      names(all_colors) <- all_color_levels
    }

    # Create inline tag/chip items
    selected_items <- lapply(ids, function(rid) {
      idx <- which(df$.rowid == rid)
      if (length(idx) == 0) return(NULL)

      key_val <- df[[key_name]][idx]

      # Get color for this item
      if (isTRUE(use_color_by)) {
        color_val <- as.character(df[[color_by]][idx])
        bg_color <- all_colors[color_val]
      } else {
        bg_color <- "#1f77b4"
      }

      # Determine text color based on background brightness
      # Simple heuristic: use white text for darker backgrounds
      tags$span(
        style = paste0(
          "display: inline-flex; align-items: center; ",
          "margin: 2px 4px 2px 0; padding: 4px 8px; ",
          "background-color: ", bg_color, "; ",
          "color: white; ",
          "border-radius: 16px; ",
          "font-size: 12px; ",
          "white-space: nowrap;"
        ),
        tags$span(as.character(key_val), style = "margin-right: 6px;"),
        tags$span(
          "\u00d7",
          style = "cursor: pointer; font-weight: bold; font-size: 14px; line-height: 1; opacity: 0.8;",
          onclick = sprintf("Shiny.setInputValue('deselect_id', %d, {priority: 'event'});", rid),
          onmouseover = "this.style.opacity='1';",
          onmouseout = "this.style.opacity='0.8';"
        )
      )
    })

    tagList(
      tags$div(
        style = "margin-bottom: 8px;",
        tags$strong(paste0("Selected ", key_name, " (", length(ids), "):"))
      ),
      tags$div(
        style = "display: flex; flex-wrap: wrap; align-items: center; max-height: 200px; overflow-y: auto;",
        selected_items
      )
    )
  })

  # Observer for deselection via JS onclick
  observeEvent(input$deselect_id, {
    rid <- input$deselect_id
    if (!is.null(rid)) {
      current_sel <- selected_rowids()
      selected_rowids(setdiff(current_sel, rid))
    }
  })

  # Compute extreme reconstruction range (for default axis limits)
  extreme_recon_range <- reactive({
    req(bundle_ok())
    b <- fpca_bundle()
    mu_vec <- get_mu_vec(b)
    phi_mat <- b$functions@X
    scores_mat <- b$scores
    p <- ncol(scores_mat)

    # Get min and max scores for each PC
    score_mins <- apply(scores_mat, 2, min, na.rm = TRUE)
    score_maxs <- apply(scores_mat, 2, max, na.rm = TRUE)

    # Reconstruct extreme curves using all combinations of min/max scores
    # For efficiency, compute contribution of each PC at its extremes
    extreme_vals <- mu_vec
    for (k in seq_len(p)) {
      phi_k <- as.numeric(phi_mat[k, ])
      contrib_min <- score_mins[k] * phi_k
      contrib_max <- score_maxs[k] * phi_k
      # Take the more extreme contribution at each time point
      extreme_vals <- c(
        extreme_vals,
        mu_vec + contrib_min,
        mu_vec + contrib_max
      )
    }

    # Also reconstruct using actual extreme observations
    for (i in seq_len(nrow(scores_mat))) {
      recon_i <- reconstruct_curve(mu_vec, phi_mat, as.numeric(scores_mat[i, ]))
      extreme_vals <- c(extreme_vals, recon_i)
    }

    list(
      y_min = min(extreme_vals, na.rm = TRUE),
      y_max = max(extreme_vals, na.rm = TRUE)
    )
  })

  # Dynamic UI for scatter plot with adjustable height
  output$scatter_plot_ui <- renderUI({
    height <- input$scatter_height
    if (is.null(height)) height <- 480
    plotlyOutput("score_scatter", height = paste0(height, "px"))
  })

  # Dynamic UI for reconstruction plot with adjustable height
  output$recon_plot_ui <- renderUI({
    height <- input$recon_height
    if (is.null(height)) height <- 320
    plotOutput("recon_plot", height = paste0(height, "px"))
  })

  # Dynamic y-axis range for reconstruction plot
  output$recon_ylim_ui <- renderUI({
    req(bundle_ok())
    extremes <- extreme_recon_range()

    # Round to reasonable precision
    default_ymin <- floor(extremes$y_min * 10) / 10
    default_ymax <- ceiling(extremes$y_max * 10) / 10

    tagList(
      tags$label("Reconstruction Y-axis range:", style = "font-weight: normal; margin-bottom: 5px;"),
      fluidRow(
        column(6, numericInput("recon_ylim_low", "Min", value = default_ymin, step = 0.1)),
        column(6, numericInput("recon_ylim_high", "Max", value = default_ymax, step = 0.1))
      )
    )
  })


  # Debounced slider values for scatter plot (avoid too many redraws)
  slider_x_debounced <- reactive({
    input$slider_x
  }) %>% debounce(100)

  slider_y_debounced <- reactive({
    input$slider_y
  }) %>% debounce(100)
  
  # Summary text
  output$bundle_info <- renderPrint({
    req(bundle_ok())
    b <- fpca_bundle()
    p <- npc()
    cat("Loaded bundle keys:\n")
    cat(paste(sort(names(b)), collapse = ", "), "\n\n")
    cat("n observations:", nrow(b$scores), "\n")
    cat("npc:", p, "\n")
    cat("map_col:", map_col(), "\n")
    if (!is.null(b$meta)) {
      cat("meta present: TRUE\n")
      col_types <- sapply(b$meta, function(x) {
        if (is.factor(x))    "factor"
        else if (is.integer(x)) "integer"
        else if (is.numeric(x)) "numeric"
        else if (is.logical(x)) "logical"
        else                    "character"
      })
      cat("meta variables:\n")
      for (nm in names(col_types)) {
        extra <- if (col_types[nm] == "factor") {
          lvls <- levels(b$meta[[nm]])
          paste0("  [", length(lvls), " levels: ", paste(head(lvls, 5), collapse = ", "),
                 if (length(lvls) > 5) ", ..." else "", "]")
        } else ""
        cat(sprintf("  %-20s %s%s\n", nm, col_types[nm], extra))
      }
    } else {
      cat("meta present: FALSE\n")
    }
    if (!is.null(b$created_at)) cat("created_at:", as.character(b$created_at), "\n")
  })
  
  # PC table + variance bar
  output$pc_table <- renderDT({
    req(bundle_ok())
    tbl <- make_pc_table(fpca_bundle()$values) %>%
      mutate(
        var_prop_pct = round(var_prop_pct, 2),
        cum_var_prop_pct = round(cum_var_prop_pct, 2)
      ) %>%
      select(PC, eigenvalue, var_prop_pct, cum_var_prop_pct)
    datatable(tbl, rownames = FALSE, options = list(pageLength = 10))
  })
  
  output$var_bar <- renderPlotly({
    req(bundle_ok())
    tbl <- make_pc_table(fpca_bundle()$values)
    plot_ly(
      data = tbl,
      x = ~PC,
      y = ~var_prop_pct,
      type = "bar"
    ) %>%
      layout(
        yaxis = list(title = "Variance explained (%)"),
        xaxis = list(title = "")
      )
  })
  
  # Eigenfunctions plot
  output$eigen_plot <- renderPlot({
    req(bundle_ok())
    pcs <- suppressWarnings(as.integer(gsub("^PC", "", input$pcs_to_plot)))
    pcs <- pcs[!is.na(pcs)]
    if (length(pcs) == 0) pcs <- 1:min(3, npc())
    
    b <- fpca_bundle()
    plot_long <- fun_pc_long(b$functions, pcs)
    
    ggplot(plot_long, aes(x = Time, y = Value, color = Curve, linetype = Curve)) +
      geom_line(linewidth = 1) +
      labs(title = "Eigenfunctions", x = "Time", y = "") +
      theme_minimal(base_size = 11)
  })
  
  # PC effect plot (mean ± SD)
  output$pc_effect_plot <- renderPlot({
    req(bundle_ok())
    b <- fpca_bundle()
    phi_mat <- b$functions@X
    time_vals <- b$functions@argvals[[1]]
    eigenvals <- b$values
    p <- npc()
    
    nPC <- min(3, p)
    sd_scores <- sqrt(eigenvals)
    varprop <- eigenvals / sum(eigenvals)
    var_labels <- paste0("PC", 1:nPC, " (", round(varprop[1:nPC] * 100), "%)")
    
    mu_vec <- get_mu_vec(b)
    
    fractions <- seq(-1, 1, by = 0.5)
    
    pc_curves <- expand_grid(
      PC = 1:nPC,
      fraction = fractions
    ) %>%
      mutate(label = factor(PC, levels = 1:nPC, labels = var_labels)) %>%
      rowwise() %>%
      mutate(y = list(mu_vec + fraction * sd_scores[PC] * as.numeric(phi_mat[PC, ]))) %>%
      ungroup() %>%
      unnest_longer(y) %>%
      group_by(PC, fraction, label) %>%
      mutate(time = time_vals) %>%
      ungroup()
    
    # mark mean
    pc_curves <- pc_curves %>%
      mutate(is_mean = (fraction == 0))
    
    ggplot(pc_curves, aes(x = time, y = y, group = fraction, color = fraction)) +
      geom_line(data = pc_curves %>% filter(!is_mean), linewidth = 0.8) +
      geom_line(data = pc_curves %>% filter(is_mean), linewidth = 1.2) +
      facet_wrap(~ label, nrow = 1) +
      labs(title = "PC effects (mean ± k·SD)", x = "Time", y = "") +
      theme_minimal(base_size = 11)
  })
  
  # Scores table
  output$scores_table <- renderDT({
    req(bundle_ok())
    datatable(pcscores(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Scatter plot (interactive) + click handling
  output$score_scatter <- renderPlotly({
    req(bundle_ok())
    # Reference height to trigger re-render on resize
    input$scatter_height
    df <- pcscores()
    req(input$x_pc, input$y_pc)

    xcol <- input$x_pc
    ycol <- input$y_pc
    key_name <- map_col()

    color_by <- input$color_by
    shape_by <- input$shape_by

    # Mark selected points
    sel_ids <- selected_rowids()
    df$is_selected <- df$.rowid %in% sel_ids

    # base tooltip
    tt <- paste0(
      key_name, ": ", df[[key_name]],
      "<br>", xcol, ": ", signif(df[[xcol]], 4),
      "<br>", ycol, ": ", signif(df[[ycol]], 4)
    )

    # Convert color_by and shape_by to factors (categorical)
    if (!is.null(color_by) && nzchar(color_by) && (color_by %in% names(df))) {
      df[[color_by]] <- as.factor(df[[color_by]])
      tt <- paste0(tt, "<br>", color_by, ": ", df[[color_by]])
    }
    if (!is.null(shape_by) && nzchar(shape_by) && (shape_by %in% names(df))) {
      df[[shape_by]] <- as.factor(df[[shape_by]])
      tt <- paste0(tt, "<br>", shape_by, ": ", df[[shape_by]])
    }
    tt <- paste0(tt, "<br>Selected: ", ifelse(df$is_selected, "Yes", "No"))

    # Determine marker sizes and outline for selection indication
    marker_size <- ifelse(df$is_selected, 12, 8)
    marker_line <- lapply(seq_len(nrow(df)), function(i) {
      if (df$is_selected[i]) list(color = "black", width = 2) else list(color = "rgba(0,0,0,0)", width = 0)
    })

    # plotly scatter with a stable key
    # Build color and symbol mappings
    use_color <- !is.null(color_by) && nzchar(color_by) && (color_by %in% names(df))
    use_shape <- !is.null(shape_by) && nzchar(shape_by) && (shape_by %in% names(df))

    # Map shape_by to plotly symbols
    if (use_shape) {
      shape_levels <- levels(df[[shape_by]])
      plotly_symbols <- c("circle", "square", "diamond", "cross", "x", "triangle-up", "triangle-down", "star")
      symbol_map <- setNames(plotly_symbols[seq_along(shape_levels)], shape_levels)
      df$.symbol <- symbol_map[as.character(df[[shape_by]])]
    } else {
      df$.symbol <- "circle"  # default symbol
    }

    # Create consistent color mapping for color_by
    if (use_color) {
      color_levels <- levels(df[[color_by]])
      color_palette <- get_color_palette(length(color_levels))
      names(color_palette) <- color_levels
      df$.color <- color_palette[as.character(df[[color_by]])]
    } else {
      df$.color <- "#1f77b4"  # default blue
    }

    # Base scatter plot
    if (use_color && use_shape) {
      p <- plot_ly(
        data = df,
        x = df[[xcol]],
        y = df[[ycol]],
        type = "scatter",
        mode = "markers",
        color = df[[color_by]],
        colors = color_palette,
        symbol = df[[shape_by]],
        symbols = plotly_symbols[seq_along(shape_levels)],
        text = tt,
        hoverinfo = "text",
        key = df$.rowid,
        customdata = df[[key_name]],
        marker = list(size = 8),
        source = "scatter_main"
      )
    } else if (use_color) {
      p <- plot_ly(
        data = df,
        x = df[[xcol]],
        y = df[[ycol]],
        type = "scatter",
        mode = "markers",
        color = df[[color_by]],
        colors = color_palette,
        text = tt,
        hoverinfo = "text",
        key = df$.rowid,
        customdata = df[[key_name]],
        marker = list(size = 8),
        source = "scatter_main"
      )
    } else if (use_shape) {
      p <- plot_ly(
        data = df,
        x = df[[xcol]],
        y = df[[ycol]],
        type = "scatter",
        mode = "markers",
        symbol = df[[shape_by]],
        symbols = plotly_symbols[seq_along(shape_levels)],
        text = tt,
        hoverinfo = "text",
        key = df$.rowid,
        customdata = df[[key_name]],
        marker = list(size = 8),
        source = "scatter_main"
      )
    } else {
      p <- plot_ly(
        data = df,
        x = df[[xcol]],
        y = df[[ycol]],
        type = "scatter",
        mode = "markers",
        text = tt,
        hoverinfo = "text",
        key = df$.rowid,
        customdata = df[[key_name]],
        marker = list(size = 8, color = "#1f77b4"),
        source = "scatter_main"
      )
    }

    # Add highlight trace for selected points (with same source for click events)
    # Use same symbol as original with transparent fill and black outline
    df_selected <- df[df$is_selected, , drop = FALSE]
    if (nrow(df_selected) > 0) {
      selected_tt <- tt[df$is_selected]
      p <- p %>% add_trace(
        x = df_selected[[xcol]],
        y = df_selected[[ycol]],
        type = "scatter",
        mode = "markers",
        marker = list(
          size = 14,
          symbol = df_selected$.symbol,
          color = "rgba(0,0,0,0)",
          line = list(color = "black", width = 3)
        ),
        text = selected_tt,
        hoverinfo = "text",
        key = df_selected$.rowid,
        name = "Selected",
        showlegend = FALSE,
        inherit = FALSE
      )
    }

    # Add crosshairs - use slider values (dot moves via JS during drag)
    show_slider <- isTRUE(input$show_slider_recon)
    slider_x_val <- slider_x_debounced()
    slider_y_val <- slider_y_debounced()

    xmin <- min(df[[xcol]], na.rm = TRUE)
    xmax <- max(df[[xcol]], na.rm = TRUE)
    ymin <- min(df[[ycol]], na.rm = TRUE)
    ymax <- max(df[[ycol]], na.rm = TRUE)

    shapes <- list()
    annotations <- list()

    if (show_slider && !is.null(slider_x_val) && !is.null(slider_y_val)) {
      # Vertical line as trace
      p <- p %>% add_segments(
        x = slider_x_val, xend = slider_x_val,
        y = ymin - 0.1 * abs(ymax - ymin), yend = ymax + 0.1 * abs(ymax - ymin),
        line = list(color = "red", width = 1.5, dash = "dash"),
        showlegend = FALSE,
        hoverinfo = "skip",
        inherit = FALSE
      )
      # Horizontal line as trace
      p <- p %>% add_segments(
        x = xmin - 0.1 * abs(xmax - xmin), xend = xmax + 0.1 * abs(xmax - xmin),
        y = slider_y_val, yend = slider_y_val,
        line = list(color = "red", width = 1.5, dash = "dash"),
        showlegend = FALSE,
        hoverinfo = "skip",
        inherit = FALSE
      )

      # Crosshair point - NOT editable, we handle dragging ourselves
      shapes[[1]] <- list(
        type = "circle",
        xref = "x",
        yref = "y",
        xsizemode = "pixel",
        ysizemode = "pixel",
        xanchor = slider_x_val,
        yanchor = slider_y_val,
        x0 = -10,
        x1 = 10,
        y0 = -10,
        y1 = 10,
        fillcolor = "rgba(255, 0, 0, 0.8)",
        line = list(color = "darkred", width = 2),
        editable = FALSE
      )
    }

    # Apply saved zoom state if available
    zoom_state <- plot_zoom()
    x_range <- if (!is.null(zoom_state)) zoom_state$x else NULL
    y_range <- if (!is.null(zoom_state)) zoom_state$y else NULL

    p %>%
      layout(
        xaxis = list(title = xcol, range = x_range),
        yaxis = list(title = ycol, range = y_range),
        dragmode = "zoom",
        shapes = shapes,
        showlegend = TRUE
      ) %>%
      config(edits = list(shapePosition = FALSE)) %>%
      htmlwidgets::onRender("
        function(el, x) {
          var gd = el;

          // Capture zoom state changes
          el.on('plotly_relayout', function(eventData) {
            if (eventData) {
              var xRange = eventData['xaxis.range[0]'] !== undefined ?
                [eventData['xaxis.range[0]'], eventData['xaxis.range[1]']] : null;
              var yRange = eventData['yaxis.range[0]'] !== undefined ?
                [eventData['yaxis.range[0]'], eventData['yaxis.range[1]']] : null;

              // Also check for autorange reset
              if (eventData['xaxis.autorange'] || eventData['yaxis.autorange']) {
                Shiny.setInputValue('plot_zoom_change', {
                  x: null,
                  y: null,
                  timestamp: new Date().getTime()
                }, {priority: 'event'});
              } else if (xRange || yRange) {
                Shiny.setInputValue('plot_zoom_change', {
                  x: xRange,
                  y: yRange,
                  timestamp: new Date().getTime()
                }, {priority: 'event'});
              }
            }
          });

          // Handle point clicks for selection
          el.on('plotly_click', function(eventData) {
            if (eventData && eventData.points && eventData.points.length > 0) {
              var pt = eventData.points[0];
              Shiny.setInputValue('scatter_click', {
                x: pt.x,
                y: pt.y,
                timestamp: new Date().getTime()
              }, {priority: 'event'});
            }
          });

          // Custom crosshair dragging
          var isDragging = false;
          var lastUpdate = 0;
          var throttleMs = 30;
          var lastDragPos = null;
          var savedDragmode = null;

          // Helper to check if click is near crosshair
          function isNearCrosshair(evt) {
            var shapes = gd.layout.shapes;
            if (!shapes || shapes.length === 0) return false;

            var xaxis = gd._fullLayout.xaxis;
            var yaxis = gd._fullLayout.yaxis;
            if (!xaxis || !yaxis) return false;

            var rect = el.getBoundingClientRect();
            var clickX = evt.clientX - rect.left;
            var clickY = evt.clientY - rect.top;

            var shape = shapes[0];
            var crossX = shape.xanchor;
            var crossY = shape.yanchor;

            var crossPixelX = xaxis.d2p(crossX) + xaxis._offset;
            var crossPixelY = yaxis.d2p(crossY) + yaxis._offset;
            var pixelDist = Math.sqrt(Math.pow(clickX - crossPixelX, 2) + Math.pow(clickY - crossPixelY, 2));

            return pixelDist < 20;
          }

          // Mousedown handler - capture phase to intercept before plotly
          el.addEventListener('mousedown', function(evt) {
            if (!isNearCrosshair(evt)) return;

            var shapes = gd.layout.shapes;
            var shape = shapes[0];

            isDragging = true;
            lastDragPos = {x: shape.xanchor, y: shape.yanchor};

            // Disable plotly's drag mode temporarily
            savedDragmode = gd.layout.dragmode;
            Plotly.relayout(gd, {dragmode: false});

            evt.preventDefault();
            evt.stopImmediatePropagation();
          }, true);

          document.addEventListener('mousemove', function(evt) {
            if (!isDragging) return;
            evt.preventDefault();

            var now = Date.now();
            if (now - lastUpdate < throttleMs) return;
            lastUpdate = now;

            var xaxis = gd._fullLayout.xaxis;
            var yaxis = gd._fullLayout.yaxis;
            if (!xaxis || !yaxis) return;

            var rect = el.getBoundingClientRect();
            var clickX = evt.clientX - rect.left;
            var clickY = evt.clientY - rect.top;

            var dataX = xaxis.p2d(clickX - xaxis._offset);
            var dataY = yaxis.p2d(clickY - yaxis._offset);

            // Clamp to axis range
            dataX = Math.max(xaxis.range[0], Math.min(xaxis.range[1], dataX));
            dataY = Math.max(yaxis.range[0], Math.min(yaxis.range[1], dataY));

            lastDragPos = {x: dataX, y: dataY};

            // Update shape position visually
            Plotly.relayout(gd, {
              'shapes[0].xanchor': dataX,
              'shapes[0].yanchor': dataY
            });

            // Update reconstruction plot
            Shiny.setInputValue('crosshair_dragging', {
              x: dataX,
              y: dataY
            }, {priority: 'event'});
          });

          document.addEventListener('mouseup', function(evt) {
            if (!isDragging) return;
            isDragging = false;

            // Restore dragmode
            if (savedDragmode !== null) {
              Plotly.relayout(gd, {dragmode: savedDragmode});
              savedDragmode = null;
            }

            // Send final position
            if (lastDragPos) {
              Shiny.setInputValue('crosshair_drag_end', {
                x: lastDragPos.x,
                y: lastDragPos.y,
                timestamp: new Date().getTime()
              }, {priority: 'event'});
            }
          });
        }
      ")
  })

  # Handle click to toggle point selection
  # Uses custom JS handler with timestamp to ensure every click is registered
  observeEvent(input$scatter_click, {
    click_data <- input$scatter_click
    if (is.null(click_data)) return()

    clicked_x <- click_data$x
    clicked_y <- click_data$y

    # Skip if no valid coordinates
    if (is.null(clicked_x) || is.null(clicked_y)) return()

    df <- pcscores()
    xcol <- input$x_pc
    ycol <- input$y_pc

    if (!(xcol %in% names(df)) || !(ycol %in% names(df))) return()

    # Find the point closest to the clicked location
    x_range <- max(df[[xcol]], na.rm = TRUE) - min(df[[xcol]], na.rm = TRUE)
    y_range <- max(df[[ycol]], na.rm = TRUE) - min(df[[ycol]], na.rm = TRUE)

    # Normalize distances to handle different axis scales
    x_dist <- (df[[xcol]] - clicked_x) / max(x_range, 1e-10)
    y_dist <- (df[[ycol]] - clicked_y) / max(y_range, 1e-10)
    distances <- sqrt(x_dist^2 + y_dist^2)
    min_idx <- which.min(distances)

    # Only accept if the distance is reasonably small (within 5% of plot range)
    if (length(min_idx) == 0 || distances[min_idx] >= 0.05) return()

    clicked_id <- df$.rowid[min_idx]

    # Toggle selection
    current_sel <- selected_rowids()

    if (clicked_id %in% current_sel) {
      # Deselect
      selected_rowids(setdiff(current_sel, clicked_id))
    } else {
      # Select (add to current selection)
      selected_rowids(c(current_sel, clicked_id))
    }
  })

  # Handle zoom state changes - preserve across re-renders
  observeEvent(input$plot_zoom_change, {
    zoom_data <- input$plot_zoom_change
    if (is.null(zoom_data)) return()
    plot_zoom(list(x = zoom_data$x, y = zoom_data$y))
  })

  # Handle real-time crosshair dragging - update live position only (fast)
  observeEvent(input$crosshair_dragging, {
    drag_data <- input$crosshair_dragging
    if (is.null(drag_data)) return()

    # Ignore late events that arrive within 150ms of drag_end
    # These are queued events from before mouseup
    time_since_drag_end <- as.numeric(Sys.time()) - drag_end_time()
    if (time_since_drag_end < 0.15) return()

    live_crosshair(list(x = drag_data$x, y = drag_data$y))
  })

  # Handle drag end - update sliders, keep live_crosshair until sliders catch up
  observeEvent(input$crosshair_drag_end, {
    drag_data <- input$crosshair_drag_end
    if (is.null(drag_data)) return()

    # Record drag end time to ignore late dragging events
    drag_end_time(as.numeric(Sys.time()))

    new_x <- drag_data$x
    new_y <- drag_data$y

    if (!is.null(new_x) && !is.null(new_y)) {
      # Keep live_crosshair at final position to prevent flash
      # It will be cleared by slider observer once sliders catch up
      live_crosshair(list(x = new_x, y = new_y))
      updateSliderInput(session, "slider_x", value = new_x)
      updateSliderInput(session, "slider_y", value = new_y)
    }
  })


  # Clear live_crosshair when sliders catch up after drag, or user manually adjusts

  observeEvent(c(input$slider_x, input$slider_y), {
    live_pos <- live_crosshair()
    if (is.null(live_pos) || is.null(input$slider_x) || is.null(input$slider_y)) return()

    values_match <- abs(input$slider_x - live_pos$x) <= 0.01 &&
                    abs(input$slider_y - live_pos$y) <= 0.01

    # Check if we recently ended a drag (within 500ms)
    time_since_drag_end <- as.numeric(Sys.time()) - drag_end_time()
    recently_ended_drag <- time_since_drag_end < 0.5 && time_since_drag_end > 0

    if (values_match && recently_ended_drag) {
      # Sliders caught up after drag end - clear live_crosshair
      live_crosshair(NULL)
    } else if (!values_match) {
      # User manually adjusted sliders - clear live_crosshair
      live_crosshair(NULL)
    }
  }, ignoreInit = TRUE)
  
  output$recon_plot <- renderPlot({
    req(bundle_ok())
    req(input$recon_mode, input$x_pc, input$y_pc)
    # Reference height to trigger re-render on resize
    input$recon_height
    df <- pcscores()
    b <- fpca_bundle()

    # data for reconstruction
    phi_mat <- b$functions@X
    time_vals <- b$functions@argvals[[1]]
    mu_vec <- get_mu_vec(b)
    p <- ncol(b$scores)
    key_name <- map_col()

    # Get color_by setting from scatter plot
    color_by <- input$color_by
    use_color_by <- is.character(color_by) &&
      length(color_by) == 1 &&
      !is.na(color_by) &&
      nzchar(color_by) &&
      (color_by %in% names(df))

    # parse selected x/y PCs
    pc_index <- function(pcname) as.integer(gsub("^PC", "", pcname))
    ax <- pc_index(input$x_pc)
    ay <- pc_index(input$y_pc)
    ax <- if (!is.na(ax) && ax >= 1 && ax <= p) ax else 1
    ay <- if (!is.na(ay) && ay >= 1 && ay <= p) ay else min(2, p)

    # Start building plot data
    plot_data <- tibble(time = time_vals, y = mu_vec, label = "Mean", type = "mean", color_group = "Mean")

    # Crosshair-based reconstruction
    # Use live_crosshair if dragging, otherwise use slider values
    show_slider <- isTRUE(input$show_slider_recon)
    live_pos <- live_crosshair()

    if (!is.null(live_pos)) {
      # Use live drag position (fast path)
      crosshair_x <- live_pos$x
      crosshair_y <- live_pos$y
    } else {
      # Use slider values
      crosshair_x <- input$slider_x
      crosshair_y <- input$slider_y
    }

    if (show_slider && !is.null(crosshair_x) && !is.null(crosshair_y)) {
      # Create a synthetic scores row with crosshair values for selected PCs
      slider_scores <- rep(0, p)
      slider_scores[ax] <- crosshair_x
      slider_scores[ay] <- crosshair_y

      slider_recon <- reconstruct_curve(mu_vec, phi_mat, slider_scores, pcs = c(ax, ay))
      slider_label <- paste0("Crosshair (", input$x_pc, "=", round(crosshair_x, 2),
                             ", ", input$y_pc, "=", round(crosshair_y, 2), ")")
      plot_data <- bind_rows(
        plot_data,
        tibble(time = time_vals, y = slider_recon, label = slider_label, type = "slider", color_group = "Crosshair")
      )
    }

    # Selected points reconstruction
    sel_ids <- selected_rowids()
    if (length(sel_ids) > 0) {
      for (rid in sel_ids) {
        i <- which(df$.rowid == rid)
        if (length(i) != 1) next

        scores_row <- as.numeric(b$scores[i, ])
        key_val <- df[[key_name]][i]

        # Get color_by value for this point
        if (isTRUE(use_color_by)) {
          color_val <- as.character(df[[color_by]][i])
        } else {
          color_val <- as.character(key_val)
        }

        # decide reconstruction mode
        if (isTRUE(input$recon_mode == "obs")) {
          K <- min(input$K, p)
          y <- reconstruct_curve(mu_vec, phi_mat, scores_row, K = K)
        } else {
          y <- reconstruct_curve(mu_vec, phi_mat, scores_row, pcs = c(ax, ay))
        }

        plot_data <- bind_rows(
          plot_data,
          tibble(time = time_vals, y = y, label = as.character(key_val), type = "selected", color_group = color_val)
        )
      }
    }

    # Check if we have anything to plot besides mean
    has_content <- nrow(plot_data) > length(time_vals)
    validate(need(has_content, "Use the sliders or click points on the scatter plot to reconstruct curves."))

    # Build title
    if (isTRUE(input$recon_mode == "obs")) {
      mode_text <- paste0("K = ", min(input$K, p), " PCs")
    } else {
      mode_text <- paste0(input$x_pc, " + ", input$y_pc, " only")
    }

    # Define colors: mean=gray dashed, slider=red, selected points by color_group
    plot_data$type <- factor(plot_data$type, levels = c("mean", "slider", "selected"))

    # Create endpoint data for line labels (file names)
    endpoint_data <- plot_data %>%
      filter(type == "selected") %>%
      group_by(label) %>%
      slice_max(order_by = time, n = 1, with_ties = FALSE) %>%
      ungroup()

    # Create a color mapping based on color_group using same palette as scatter plot
    unique_color_groups <- unique(plot_data$color_group[plot_data$type == "selected"])

    # Get all levels of color_by to ensure consistent color mapping
    if (isTRUE(use_color_by)) {
      all_color_levels <- unique(as.character(df[[color_by]]))
      all_colors <- get_color_palette(length(all_color_levels))
      names(all_colors) <- all_color_levels
      # Extract only the colors we need
      group_colors <- all_colors[unique_color_groups]
    } else if (length(unique_color_groups) > 0) {
      group_colors <- get_color_palette(length(unique_color_groups))
      names(group_colors) <- unique_color_groups
    } else {
      group_colors <- character(0)
    }

    # Build full color map
    color_map <- c("Mean" = "gray40", "Crosshair" = "red")
    color_map <- c(color_map, group_colors)

    # Create linetype map
    linetype_map <- c("Mean" = "dashed", "Crosshair" = "solid")
    for (grp in unique_color_groups) {
      linetype_map[grp] <- "solid"
    }

    # Get y-axis limits
    ylim_low <- input$recon_ylim_low
    ylim_high <- input$recon_ylim_high

    # Determine legend title
    legend_title <- if (isTRUE(use_color_by)) color_by else key_name

    plt <- ggplot(plot_data, aes(x = time, y = y, color = color_group, linetype = color_group, group = label)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = color_map, name = legend_title) +
      scale_linetype_manual(values = linetype_map, name = legend_title) +
      labs(
        title = paste0("Reconstructed Curves (", mode_text, ")"),
        x = "Time",
        y = ""
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        plot.margin = margin(5.5, 40, 5.5, 5.5)  # Extra right margin for labels
      )

    # Add direct line labels for file names (simple fixed position, no repel)
    if (nrow(endpoint_data) > 0) {
      plt <- plt +
        geom_text(
          data = endpoint_data,
          aes(x = time, y = y, label = label),
          hjust = 0,
          nudge_x = 0.5,
          size = 3,
          show.legend = FALSE
        )
    }

    # Apply y-axis limits if set, use clip="off" to allow labels outside panel
    ylim_low  <- input$recon_ylim_low
    ylim_high <- input$recon_ylim_high

    if (is.finite(ylim_low) && is.finite(ylim_high)) {
      plt <- plt + coord_cartesian(ylim = c(ylim_low, ylim_high), clip = "off")
    } else {
      plt <- plt + coord_cartesian(clip = "off")
    }

    plt
  })
  
  # Update K slider max when npc changes
  observeEvent(npc(), {
    updateNumericInput(session, "K", max = npc(), value = min(5, npc()))
  })
  
}

shinyApp(ui, server)
