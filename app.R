# ============================================================
# CHIPMUNK - Shoe Search with Principal Component Analysis
# Combined version of app.R and app(1).R
# All interface text, messages, and comments are in English.
# ============================================================

# 1. PACKAGES AND SETTINGS -----------------------------------
required_packages <- c("shiny", "bslib", "png", "magick")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages first: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

library(shiny)
library(bslib)
library(png)
library(magick)
options(shiny.maxRequestSize = 30 * 1024^2)

IMAGE_SIDE <- 64L
PIXEL_COUNT <- IMAGE_SIDE * IMAGE_SIDE
NEIGHBORS <- 5L
MAX_PCS <- 100L
CACHE_FILE <- "pca_cache_unified_v2.rds"
CACHE_VERSION <- "column-major-png-grayscale-v2"

# 2. CATALOGUE DATA ------------------------------------------
if (!file.exists("shoes_meta.csv")) {
  stop("Missing shoes_meta.csv. Put it in the same folder as app.R.")
}

shoes_meta <- read.csv("shoes_meta.csv", stringsAsFactors = FALSE)
if (!all(c("id", "category") %in% names(shoes_meta))) {
  stop("shoes_meta.csv must contain columns named 'id' and 'category'.")
}

shoes_meta$id <- suppressWarnings(as.integer(as.character(shoes_meta$id)))
shoes_meta$category <- trimws(as.character(shoes_meta$category))
if (anyNA(shoes_meta$id) || anyDuplicated(shoes_meta$id) ||
    anyNA(shoes_meta$category) || any(!nzchar(shoes_meta$category))) {
  stop("Catalogue IDs must be unique integers and categories must not be empty.")
}
if (nrow(shoes_meta) <= NEIGHBORS) {
  stop("The catalogue must have at least six images for five-neighbor search.")
}

catalogue_paths <- file.path("shoes", sprintf("%04d.png", shoes_meta$id))
missing_images <- catalogue_paths[!file.exists(catalogue_paths)]
if (length(missing_images) > 0) {
  stop("Catalogue images are missing. First missing file: ", missing_images[1])
}

category_levels <- sort(unique(shoes_meta$category))
category_colors <- setNames(
  grDevices::hcl.colors(length(category_levels), palette = "Dark 3"),
  category_levels
)

# 3. IMAGE PREPROCESSING -------------------------------------
# IMPORTANT: Both catalogue and external images become column-major
# 64 x 64 grayscale vectors, using png::readPNG as the final decoder.
# Do not mix this representation with a row-major magick vector.
read_png_vector <- function(path) {
  img <- png::readPNG(path)
  dims <- dim(img)
  
  if (length(dims) == 3L) {
    if (dims[3] >= 3L) {
      # Convert color PNGs to grayscale; catalogue files are normally gray.
      gray <- 0.2126 * img[, , 1] +
        0.7152 * img[, , 2] +
        0.0722 * img[, , 3]
      if (dims[3] >= 4L) {
        alpha <- img[, , 4]
        gray <- gray * alpha + (1 - alpha)  # White background
      }
      img <- gray
    } else {
      img <- img[, , 1]
    }
  }
  
  if (!identical(dim(img), c(IMAGE_SIDE, IMAGE_SIDE))) {
    stop("Expected a 64 x 64 PNG: ", path)
  }
  as.numeric(round(pmax(0, pmin(1, img)) * 255))
}

read_external_vector <- function(path) {
  # Convert real photographs in a consistent way: orient, flatten alpha,
  # pad to a square without stretching, resize and convert to grayscale.
  img <- magick::image_read(path)[1]
  dims <- magick::image_info(img)
  
  # An existing 64x64 PNG should remain pixel-exact. This also works for
  # Shiny temporary uploads whose paths have no .png file extension.
  if (dims$width[1] == IMAGE_SIDE && dims$height[1] == IMAGE_SIDE) {
    exact_png <- tryCatch(png::readPNG(path), error = function(e) NULL)
    if (!is.null(exact_png)) return(read_png_vector(path))
  }
  
  img <- magick::image_orient(img)
  img <- magick::image_background(img, "white", flatten = TRUE)
  dims <- magick::image_info(img)
  side <- max(dims$width[1], dims$height[1])
  img <- magick::image_extent(
    img, geometry = paste0(side, "x", side),
    gravity = "center", color = "white"
  )
  img <- magick::image_resize(img, "64x64!")
  img <- magick::image_convert(img, colorspace = "gray")
  
  temporary_png <- tempfile(fileext = ".png")
  on.exit(unlink(temporary_png), add = TRUE)
  magick::image_write(img, path = temporary_png, format = "png", depth = 8)
  read_png_vector(temporary_png)
}

save_vector_png <- function(vector) {
  stopifnot(length(vector) == PIXEL_COUNT)
  values <- pmax(0, pmin(255, as.numeric(vector))) / 255
  image_matrix <- matrix(values, nrow = IMAGE_SIDE, ncol = IMAGE_SIDE)
  destination <- tempfile(fileext = ".png")
  png::writePNG(image_matrix, destination)
  destination
}

# Match a self-check filename to its optional ground-truth label.
read_selfcheck_label <- function(filename) {
  path <- file.path("selfcheck", "labels.csv")
  if (!file.exists(path)) return(NA_character_)
  
  labels <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                     error = function(e) NULL)
  if (is.null(labels) || nrow(labels) == 0L) return(NA_character_)
  
  headers <- tolower(names(labels))
  file_col <- match(TRUE, headers %in% c(
    "file", "filename", "file_name", "image", "photo", "name", "id"
  ))
  cat_col <- match(TRUE, headers %in% c("category", "label", "class"))
  if (is.na(file_col) || is.na(cat_col)) return(NA_character_)
  
  normalize_filename <- function(x) {
    x <- tolower(tools::file_path_sans_ext(basename(as.character(x))))
    sub("^selfcheck_0*", "", x)
  }
  row_index <- match(normalize_filename(filename),
                     normalize_filename(labels[[file_col]]))
  if (is.na(row_index)) return(NA_character_)
  
  value <- trimws(as.character(labels[[cat_col]][row_index]))
  if (is.na(value) || !nzchar(value)) return(NA_character_)
  known <- match(tolower(value), tolower(category_levels))
  if (is.na(known)) return(value)
  category_levels[known]
}

# 4. PCA WITH A VERSIONED, VALIDATED CACHE --------------------
# Changing the processing algorithm or source files invalidates the cache.
info <- file.info(catalogue_paths)
cache_manifest <- list(
  version = CACHE_VERSION,
  ids = shoes_meta$id,
  categories = shoes_meta$category,
  sizes = as.numeric(info$size),
  modified = as.numeric(info$mtime)
)

cached <- if (file.exists(CACHE_FILE)) {
  tryCatch(readRDS(CACHE_FILE), error = function(e) NULL)
} else {
  NULL
}

cache_ok <- is.list(cached) && identical(cached$manifest, cache_manifest) &&
  is.matrix(cached$X) && nrow(cached$X) == nrow(shoes_meta) &&
  ncol(cached$X) == PIXEL_COUNT &&
  is.matrix(cached$loadings) && is.matrix(cached$scores) &&
  length(cached$x_bar) == PIXEL_COUNT &&
  nrow(cached$scores) == nrow(shoes_meta) &&
  ncol(cached$scores) == ncol(cached$loadings) &&
  length(cached$var_explained) == ncol(cached$loadings)

# Spot-check actual pixels even when the cache metadata appears unchanged.
# This catches a cache built with a different pixel order or decoder.
if (cache_ok) {
  probes <- unique(c(1L, ceiling(nrow(shoes_meta) / 2), nrow(shoes_meta)))
  cache_ok <- all(vapply(probes, function(i) {
    identical(as.numeric(cached$X[i, ]), read_png_vector(catalogue_paths[i]))
  }, logical(1)))
}

if (cache_ok) {
  message("Loading validated PCA cache: ", CACHE_FILE)
  pca_data <- cached
} else {
  message("Building PCA catalogue from ", nrow(shoes_meta), " images...")
  X <- t(vapply(catalogue_paths, read_png_vector, numeric(PIXEL_COUNT)))
  x_bar <- colMeans(X)
  centered <- sweep(X, 2, x_bar, "-")
  components <- min(MAX_PCS, nrow(X) - 1L, ncol(X))
  
  # PCA without variance scaling: brightness/shape variation is preserved.
  decomposition <- svd(centered, nu = 0L, nv = components)
  loadings <- decomposition$v[, seq_len(components), drop = FALSE]
  scores <- centered %*% loadings
  explained <- 100 * decomposition$d[seq_len(components)]^2 /
    sum(decomposition$d^2)
  
  pca_data <- list(
    manifest = cache_manifest,
    X = X,
    x_bar = x_bar,
    loadings = loadings,
    scores = scores,
    var_explained = explained
  )
  saveRDS(pca_data, CACHE_FILE)
  message("PCA cache saved to ", CACHE_FILE)
}

K_MAX <- ncol(pca_data$loadings)
stopifnot(K_MAX >= 2L)

# 5. MATCHING, RECONSTRUCTION, AND ACCURACY -------------------
predict_from_categories <- function(categories_in_distance_order) {
  counts <- table(categories_in_distance_order)
  tied <- names(counts)[counts == max(counts)]
  # Because neighbors are already sorted, the nearest tied class wins.
  categories_in_distance_order[
    which(categories_in_distance_order %in% tied)[1]
  ]
}

find_similar_shoes <- function(query_vector, K, exclude_id = NULL) {
  stopifnot(length(query_vector) == PIXEL_COUNT,
            K >= 1L, K <= K_MAX)
  indices <- seq_len(K)
  V <- pca_data$loadings[, indices, drop = FALSE]
  query_scores <- as.numeric((query_vector - pca_data$x_bar) %*% V)
  catalogue_scores <- pca_data$scores[, indices, drop = FALSE]
  distance_squared <- rowSums(sweep(catalogue_scores, 2, query_scores, "-")^2)
  
  if (!is.null(exclude_id)) {
    excluded_row <- match(as.integer(exclude_id), shoes_meta$id)
    if (!is.na(excluded_row)) distance_squared[excluded_row] <- Inf
  }
  
  nearest <- order(distance_squared)[seq_len(NEIGHBORS)]
  results <- data.frame(
    Rank = seq_len(NEIGHBORS),
    ID = shoes_meta$id[nearest],
    Category = shoes_meta$category[nearest],
    Distance = round(sqrt(distance_squared[nearest]), 2),
    row.names = NULL
  )
  list(
    best_id = results$ID[1],
    top5 = results,
    prediction = predict_from_categories(results$Category),
    query_scores = query_scores
  )
}

reconstruct_shoe <- function(vector, K) {
  V <- pca_data$loadings[, seq_len(K), drop = FALSE]
  q <- as.numeric((vector - pca_data$x_bar) %*% V)
  as.numeric(pca_data$x_bar + V %*% q)
}

# Plot in the same orientation as png::writePNG: top-left pixel stays top-left.
plot_shoe <- function(vector, title = "", signed = FALSE) {
  m <- matrix(as.numeric(vector), IMAGE_SIDE, IMAGE_SIDE)
  z <- t(m[IMAGE_SIDE:1L, , drop = FALSE])
  
  if (signed) {
    bound <- max(abs(vector))
    if (!is.finite(bound) || bound == 0) bound <- 1
    palette <- grDevices::colorRampPalette(
      c("#245A9C", "#FFFFFF", "#C64342")
    )(201)
    limits <- c(-bound, bound)
  } else {
    palette <- gray(seq(0, 1, length.out = 256))
    limits <- c(0, 255)
  }
  
  image(
    x = seq_len(IMAGE_SIDE), y = seq_len(IMAGE_SIDE), z = z,
    col = palette, zlim = limits,
    axes = FALSE, xlab = "", ylab = "", asp = 1,
    main = title, useRaster = TRUE
  )
}

# Optional, expensive computation: leave one catalogue image out, then
# classify it using its five closest remaining neighbours at each K.
accuracy_cache <- new.env(parent = emptyenv())
compute_loo_accuracy <- function(progress_callback = function(k) NULL) {
  if (exists("results", envir = accuracy_cache, inherits = FALSE)) {
    return(get("results", envir = accuracy_cache, inherits = FALSE))
  }
  
  S <- pca_data$scores
  n <- nrow(S)
  truth_factor <- factor(shoes_meta$category, levels = category_levels)
  truth <- as.integer(truth_factor)
  category_count <- length(category_levels)
  D2 <- matrix(0, n, n)
  overall <- numeric(K_MAX)
  by_category <- matrix(
    NA_real_, nrow = K_MAX, ncol = category_count,
    dimnames = list(NULL, category_levels)
  )
  
  for (k in seq_len(K_MAX)) {
    D2 <- D2 + outer(S[, k], S[, k], "-")^2
    diag(D2) <- Inf
    
    neighbor_rows <- t(vapply(seq_len(n), function(i) {
      row_distances <- D2[i, ]
      threshold <- sort.int(row_distances, partial = NEIGHBORS)[NEIGHBORS]
      candidates <- which(row_distances <= threshold)
      candidates[order(row_distances[candidates], candidates)][seq_len(NEIGHBORS)]
    }, integer(NEIGHBORS)))
    
    predicted <- vapply(seq_len(n), function(i) {
      neighbor_classes <- truth[neighbor_rows[i, ]]
      counts <- tabulate(neighbor_classes, nbins = category_count)
      tied <- which(counts == max(counts))
      neighbor_classes[which(neighbor_classes %in% tied)[1]]
    }, integer(1))
    
    correct <- predicted == truth
    overall[k] <- mean(correct)
    by_category[k, ] <- vapply(seq_len(category_count), function(j) {
      mean(correct[truth == j])
    }, numeric(1))
    progress_callback(k)
  }
  
  result <- list(K = seq_len(K_MAX), overall = overall, by_cat = by_category)
  assign("results", result, envir = accuracy_cache)
  result
}

# 6. PERSONAL SHOE DATA --------------------------------------
parse_our_filename <- function(filename) {
  parts <- strsplit(
    tools::file_path_sans_ext(basename(filename)), "_", fixed = TRUE
  )[[1]]
  member <- if (length(parts) >= 1L) parts[1] else "Unknown"
  category <- if (length(parts) >= 2L) parts[2] else NA_character_
  matching <- match(tolower(category), tolower(category_levels))
  category <- if (is.na(matching)) NA_character_ else category_levels[matching]
  list(name = member, category = category)
}

empty_our <- function() {
  list(
    df = data.frame(
      file = character(), name = character(), category = character(),
      stringsAsFactors = FALSE
    ),
    X = matrix(numeric(0), nrow = 0L, ncol = PIXEL_COUNT)
  )
}

load_folder_shoes <- function() {
  result <- empty_our()
  if (!dir.exists("ourshoes")) return(result)
  filenames <- sort(list.files(
    "ourshoes", pattern = "\\.(png|jpg|jpeg)$", ignore.case = TRUE
  ))
  
  for (filename in filenames) {
    vector <- tryCatch(
      read_external_vector(file.path("ourshoes", filename)),
      error = function(e) {
        warning("Skipping unreadable personal shoe photo: ", filename)
        NULL
      }
    )
    if (is.null(vector)) next
    details <- parse_our_filename(filename)
    result$df <- rbind(result$df, data.frame(
      file = filename, name = details$name,
      category = details$category, stringsAsFactors = FALSE
    ))
    result$X <- rbind(result$X, vector)
  }
  result
}

initial_our <- load_folder_shoes()
selfcheck_files <- if (dir.exists("selfcheck")) {
  sort(list.files("selfcheck", pattern = "\\.(png|jpg|jpeg)$",
                  ignore.case = TRUE))
} else {
  character(0)
}

# 7. USER INTERFACE ------------------------------------------
# A scrollable content panel for each top-level navigation tab.
# The navbar remains visible; only the selected tab content scrolls.
scroll_tab <- function(...) {
  div(class = "chipmunk-scroll-pane", ...)
}

shoe_image_card <- function(title, output_id) {
  card(
    fill = FALSE,
    card_header(title),
    card_body(fill = FALSE, fillable = FALSE,
              imageOutput(output_id, width = "100%", height = "auto"))
  )
}

ui <- page_navbar(
  title = "Chipmunk | Shoe PCA",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  # Keep the navbar and active tab fitted to the browser viewport.
  fillable = TRUE,
  fillable_mobile = TRUE,
  padding = 0,
  header = tags$head(tags$style(HTML(" 
    /* App window: top navigation is always visible. */
    html { height: 100%; }
    body.bslib-page-navbar {
      height: 100vh;
      height: 100dvh;
      max-height: 100dvh;
      min-height: 0;
      overflow: hidden !important;
      display: flex;
      flex-direction: column;
    }
    body.bslib-page-navbar > .navbar {
      flex: 0 0 auto;
      z-index: 1030;
    }
    /* bslib's content container fills the space below the navbar. */
    body.bslib-page-navbar > .container,
    body.bslib-page-navbar > .container-fluid {
      min-height: 0;
      flex: 1 1 auto;
      overflow: hidden;
    }
    body.bslib-page-navbar .tab-content {
      min-height: 0;
      flex: 1 1 auto;
      overflow: hidden;
    }
    body.bslib-page-navbar .tab-content > .tab-pane {
      min-height: 0;
      max-height: 100%;
      overflow: hidden;
    }
    body.bslib-page-navbar .tab-content > .tab-pane.active {
      display: flex;
      flex-direction: column;
      height: 100%;
    }
    /* Independent scrollbar on the right of the selected tab. */
    .chipmunk-scroll-pane {
      flex: 1 1 auto;
      width: 100%;
      height: 100%;
      min-height: 0;
      overflow-y: auto;
      overflow-x: hidden;
      overscroll-behavior-y: contain;
      scrollbar-gutter: stable;
      scrollbar-width: thin;
      scrollbar-color: #93a7b4 #edf2f5;
      padding: 18px 22px 28px;
    }
    .chipmunk-scroll-pane::-webkit-scrollbar { width: 10px; }
    .chipmunk-scroll-pane::-webkit-scrollbar-track { background: #edf2f5; }
    .chipmunk-scroll-pane::-webkit-scrollbar-thumb {
      background: #93a7b4;
      border-radius: 8px;
      border: 2px solid #edf2f5;
    }
    @media (max-width: 767px) {
      .chipmunk-scroll-pane { padding: 12px 12px 24px; }
    }
    .bslib-card { margin-bottom: 12px; }
    .shiny-image-output img { object-fit: contain; }
    .shoe-caption { color: #576574; font-size: 0.93rem; }
  "))),
  
  nav_panel(
    "Shoe Search",
    scroll_tab(
      layout_sidebar(
        fill = FALSE, fillable = FALSE,
        sidebar = sidebar(
          title = "Search Controls", width = 310,
          radioButtons("src", "Photo Source",
                       choices = c("Catalogue", "Self-check", "Upload"),
                       selected = "Catalogue", inline = TRUE),
          conditionalPanel(
            "input.src == 'Catalogue'",
            numericInput("cat_id", "Catalogue ID",
                         value = shoes_meta$id[1],
                         min = min(shoes_meta$id), max = max(shoes_meta$id))
          ),
          conditionalPanel(
            "input.src == 'Self-check'",
            selectInput("sc_file", "Self-check Photo", choices = selfcheck_files),
            p(class = "shoe-caption", "If labels.csv is available, the true category is displayed.")
          ),
          conditionalPanel(
            "input.src == 'Upload'",
            fileInput("photo", "Upload Photo (JPG/PNG)",
                      accept = c(".jpg", ".jpeg", ".png"))
          ),
          hr(),
          sliderInput("K", "Number of Principal Components (K)",
                      min = 1, max = K_MAX, value = min(20L, K_MAX), step = 1),
          sliderInput("image_size", "Image Size", min = 96, max = 240,
                      value = 160, step = 8, post = " px"),
          uiOutput("image_size_css")
        ),
        card(
          fill = FALSE,
          card_header("Search Results and Prediction"),
          card_body(
            fill = FALSE, fillable = FALSE,
            uiOutput("prediction"),
            hr(),
            layout_columns(
              col_widths = c(5, 7), fill = FALSE, fillable = FALSE,
              shoe_image_card("Best Match", "best_match"),
              card(fill = FALSE, card_header("Top 5 Similar Shoes"),
                   tableOutput("top5"))
            ),
            hr(),
            h4("PCA Image Comparison"),
            layout_columns(
              col_widths = c(4, 4, 4), fill = FALSE, fillable = FALSE,
              shoe_image_card("Selected Shoe", "original_shoe"),
              shoe_image_card("Average Shoe", "average_shoe"),
              shoe_image_card("PCA Reconstruction", "reconstructed_shoe")
            ),
            h4("Top 5 Nearest Neighbors"),
            plotOutput("search_top5_gallery", height = "230px")
          )
        )
      )
      
    )
  ),
  
  nav_panel(
    "Explain PCA",
    scroll_tab(
      p("Explore how the photo selected on the Shoe Search tab is represented and matched."),
      layout_columns(
        col_widths = c(6, 6), fill = FALSE, fillable = FALSE,
        card(card_header("Step 1 | Average shoe"),
             p("The mean image of all catalogue photos (x-bar)."),
             plotOutput("pca_step1", height = "250px")),
        card(card_header("Step 2 | Photo minus average"),
             p("Blue means below the average pixel value; red means above it."),
             plotOutput("pca_step2", height = "250px"))
      ),
      card(card_header("Step 3 | PCA score map (PC1 and PC2)"),
           p("Catalogue points are colored by category; your selected photo is highlighted."),
           plotOutput("pca_step3", height = "410px")),
      card(card_header("Step 4 | First 10 eigenshoes"),
           p("Principal component loadings, with their individual explained variance."),
           plotOutput("pca_step4", height = "440px")),
      card(card_header("Step 5 | Query scores"),
           p("The first K PCA coordinates used to calculate nearest-neighbor distances."),
           tableOutput("pca_step5_table")),
      card(card_header("Step 6 | Reconstruct using different K values"),
           p("Compare the same photo reconstructed using different numbers of components."),
           plotOutput("pca_step6", height = "300px")),
      card(card_header("Classification accuracy versus K"),
           p("Leave-one-out evaluation of the catalogue: each photo is classified using the five nearest other photos. Computing all K values can be slow."),
           actionButton("run_acc", "Compute leave-one-out accuracy", class = "btn-primary"),
           plotOutput("acc_plot", height = "400px"),
           textOutput("acc_text"),
           p(class = "shoe-caption",
             "Choosing K based on this same catalogue gives an optimistic estimate of performance on new photos. Test personal shoes separately."))
      
    )
  ),
  
  nav_panel(
    "Our Shoes",
    scroll_tab(
      p("Test images the PCA model has not seen. Use filenames such as Alex_Sandals_1.jpg (Name_Category_Number)."),
      p("Uploaded images are temporary. To keep photos across restarts, save them in the ourshoes folder."),
      fileInput("our_upload", "Add Our Photos", multiple = TRUE,
                accept = c(".jpg", ".jpeg", ".png")),
      actionButton("our_reset", "Remove Uploaded Photos"),
      uiOutput("our_upload_msg"),
      uiOutput("our_empty_msg"),
      card(
        card_header("Analyze a Selected Shoe"),
        selectInput("our_pick", "Choose a Shoe", choices = character(0), width = "100%"),
        uiOutput("our_prediction"),
        plotOutput("our_best_plot", height = "300px"),
        h5("Five Closest Catalogue Photos"),
        tableOutput("our_top5_table"),
        plotOutput("our_top5_plot", height = "220px")
      ),
      card(card_header("All Our Photos"), plotOutput("our_grid", height = "auto")),
      card(card_header("Overall Results at the Selected K"),
           textOutput("our_summary"), tableOutput("our_all_table"))
      
    )
  ),
  
  nav_panel(
    "Data 1",
    scroll_tab(
      p("Descriptive exploration of a Project 1 CSV dataset. Upload the original data; no results are invented."),
      fileInput("data1_file", "Upload Project 1 CSV", accept = ".csv"),
      layout_columns(
        col_widths = c(4, 4, 4), fill = FALSE, fillable = FALSE,
        selectInput("data1_group", "Group By", choices = c("No grouping" = "")),
        selectInput("data1_filter", "Filter Column", choices = c("No filter" = "")),
        selectInput("data1_value", "Filter Value", choices = c("All values" = ""))
      ),
      selectInput("data1_measure", "Numeric Measure", choices = character(0)),
      h4("Descriptive Results"),
      tableOutput("data1_summary"),
      textOutput("data1_description"),
      h4("Data Preview"),
      tableOutput("data1_preview"),
      p(class = "shoe-caption",
        "Complete the Project 1 research question, interpretation, and next-analysis plan using your actual dataset.")
      
    )
  ),
  
  nav_panel(
    "AI",
    scroll_tab(
      card(
        card_header("Ethics Statement"),
        p("We understand the rules on data and AI use for this project, including following the license and privacy rules of every dataset we used, understanding all of our code, and knowing that presenting AI output as our own work may be academic misconduct. We certify that we have followed all of them. Any violation is our responsibility, and we agree to receive an F in this course if we have not."),
        p(strong("Signed: TODO - all team members must review and sign."))
      ),
      card(
        card_header("What AI Did and What the Team Did"),
        p("Draft development notes. Verify and edit these statements before submitting."),
        tableOutput("ai_table")
      ),
      card(
        card_header("Individual Reflections"),
        p("Each member must write or confirm their own reflection."),
        tableOutput("ai_reflection")
      )
      
    )
  )
)

# 8. REACTIVE SERVER -----------------------------------------
server <- function(input, output, session) {
  
  # Search controls: resize output images without changing the 64 x 64 PCA data.
  output$image_size_css <- renderUI({
    size <- input$image_size
    best_size <- round(size * 4 / 3)
    tags$style(HTML(sprintf(" 
      #best_match { width: min(100%%, %dpx) !important; }
      #original_shoe, #average_shoe, #reconstructed_shoe {
        width: min(100%%, %dpx) !important;
      }
      #best_match, #original_shoe, #average_shoe, #reconstructed_shoe {
        height: auto !important; aspect-ratio: 1 / 1;
      }
      #best_match img, #original_shoe img, #average_shoe img,
      #reconstructed_shoe img { width: 100%% !important; height: 100%% !important; }
    ", best_size, size)))
  })
  
  # One query representation is shared by both the Search and Explain PCA tabs.
  query_data <- reactive({
    req(input$src)
    
    if (identical(input$src, "Catalogue")) {
      req(input$cat_id)
      index <- match(as.integer(input$cat_id), shoes_meta$id)
      validate(need(!is.na(index), "Catalogue ID not found."))
      return(list(
        vector = as.numeric(pca_data$X[index, ]),
        actual = shoes_meta$category[index],
        exclude_id = shoes_meta$id[index]
      ))
    }
    
    if (identical(input$src, "Self-check")) {
      req(input$sc_file)
      path <- file.path("selfcheck", basename(input$sc_file))
      validate(need(file.exists(path), "Self-check photo not found."))
      vector <- tryCatch(read_external_vector(path), error = function(e) NULL)
      validate(need(!is.null(vector), "Could not read the self-check photo."))
      return(list(
        vector = vector,
        actual = read_selfcheck_label(input$sc_file),
        exclude_id = NULL
      ))
    }
    
    req(input$photo)
    vector <- tryCatch(
      read_external_vector(input$photo$datapath), error = function(e) NULL
    )
    validate(need(!is.null(vector), "Could not read the uploaded image."))
    list(vector = vector, actual = NA_character_, exclude_id = NULL)
  })
  
  shoe_search <- reactive({
    query <- query_data()
    find_similar_shoes(query$vector, input$K, query$exclude_id)
  })
  
  output$prediction <- renderUI({
    result <- shoe_search()
    actual <- query_data()$actual
    known <- !is.na(actual) && nzchar(actual)
    verdict <- if (!known) "Unknown (no ground-truth category)" else if (
      tolower(result$prediction) == tolower(actual)
    ) "Correct" else "Incorrect"
    
    tagList(
      p(strong("Predicted Category: "), result$prediction),
      p(strong("Actual Category: "), if (known) actual else "Not available"),
      p(strong("Result: "), verdict)
    )
  })
  
  output$top5 <- renderTable(shoe_search()$top5,
                             striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$best_match <- renderImage({
    id <- shoe_search()$best_id
    list(src = file.path("shoes", sprintf("%04d.png", id)),
         contentType = "image/png", alt = paste("Catalogue shoe", id))
  }, deleteFile = FALSE)
  
  output$original_shoe <- renderImage({
    list(src = save_vector_png(query_data()$vector),
         contentType = "image/png", alt = "Selected shoe")
  }, deleteFile = TRUE)
  
  output$average_shoe <- renderImage({
    list(src = save_vector_png(pca_data$x_bar),
         contentType = "image/png", alt = "Average shoe")
  }, deleteFile = TRUE)
  
  output$reconstructed_shoe <- renderImage({
    list(src = save_vector_png(reconstruct_shoe(query_data()$vector, input$K)),
         contentType = "image/png",
         alt = paste("PCA reconstruction using", input$K, "components"))
  }, deleteFile = TRUE)
  
  output$search_top5_gallery <- renderPlot({
    result <- shoe_search()$top5
    old <- par(mfrow = c(1, NEIGHBORS), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (j in seq_len(NEIGHBORS)) {
      row <- match(result$ID[j], shoes_meta$id)
      plot_shoe(pca_data$X[row, ],
                paste0("#", j, " ID ", result$ID[j], "\n", result$Category[j]))
    }
  })
  
  # EXPLAIN PCA ----------------------------------------------
  output$pca_step1 <- renderPlot({
    plot_shoe(pca_data$x_bar, "Average shoe")
  })
  
  output$pca_step2 <- renderPlot({
    difference <- query_data()$vector - pca_data$x_bar
    plot_shoe(difference, "Selected photo - average", signed = TRUE)
  })
  
  output$pca_step3 <- renderPlot({
    query <- query_data()$vector
    two_PC <- pca_data$loadings[, 1:2, drop = FALSE]
    position <- as.numeric((query - pca_data$x_bar) %*% two_PC)
    coordinates <- pca_data$scores[, 1:2, drop = FALSE]
    
    plot(
      coordinates[, 1], coordinates[, 2],
      col = grDevices::adjustcolor(
        category_colors[shoes_meta$category], alpha.f = 0.5
      ),
      pch = 16, cex = 0.8,
      xlim = range(c(coordinates[, 1], position[1])),
      ylim = range(c(coordinates[, 2], position[2])),
      xlab = "PC1 Score", ylab = "PC2 Score", main = "PCA Score Map"
    )
    points(position[1], position[2], pch = 23, bg = "black",
           col = "#FFD43B", cex = 2.5, lwd = 2)
    legend(
      "topright", legend = c(category_levels, "Selected photo"),
      col = c(unname(category_colors), "#FFD43B"),
      pt.bg = c(rep(NA_character_, length(category_levels)), "black"),
      pch = c(rep(16, length(category_levels)), 23),
      cex = 0.8, bg = "white"
    )
  })
  
  output$pca_step4 <- renderPlot({
    n_show <- min(10L, K_MAX)
    old <- par(mfrow = c(2, 5), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (k in seq_len(n_show)) {
      plot_shoe(
        pca_data$loadings[, k],
        sprintf("PC %d (%.1f%%)", k, pca_data$var_explained[k]),
        signed = TRUE
      )
    }
  })
  
  output$pca_step5_table <- renderTable({
    score_vector <- shoe_search()$query_scores
    data.frame(
      PC = paste0("PC", seq_along(score_vector)),
      Score = round(score_vector, 2)
    )
  }, striped = TRUE, spacing = "s")
  
  output$pca_step6 <- renderPlot({
    ks <- sort(unique(c(1L, 5L, 20L, min(50L, K_MAX), input$K)))
    ks <- ks[ks <= K_MAX]
    old <- par(mfrow = c(1, length(ks)), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (k in ks) {
      plot_shoe(reconstruct_shoe(query_data()$vector, k),
                paste0("K = ", k, if (k == input$K) "\n(selected)" else ""))
    }
  })
  
  accuracy_result <- eventReactive(input$run_acc, {
    withProgress(message = "Computing leave-one-out accuracy", value = 0, {
      compute_loo_accuracy(function(k) {
        incProgress(1 / K_MAX, detail = paste("K =", k))
      })
    })
  })
  
  output$acc_plot <- renderPlot({
    validate(need(input$run_acc > 0,
                  "Click the button to calculate accuracy across K values."))
    result <- accuracy_result()
    best_k <- which.max(result$overall)
    plot(
      result$K, 100 * result$overall, type = "l", lwd = 3,
      ylim = c(0, 100),
      xlab = "Number of PCs (K)", ylab = "Accuracy (%)",
      main = "Catalogue leave-one-out accuracy"
    )
    for (j in seq_along(category_levels)) {
      lines(result$K, 100 * result$by_cat[, j],
            col = category_colors[j], lwd = 1.5, lty = 2)
    }
    abline(v = input$K, col = "red", lty = 3)
    points(input$K, 100 * result$overall[input$K], col = "red", pch = 19)
    points(best_k, 100 * result$overall[best_k], pch = 8, cex = 2)
    legend(
      "bottomright", bg = "white", cex = 0.8,
      legend = c("Overall", category_levels,
                 paste("Selected K =", input$K), paste("Best K =", best_k)),
      col = c("black", unname(category_colors), "red", "black"),
      lty = c(1, rep(2, length(category_levels)), 3, NA),
      lwd = c(3, rep(1.5, length(category_levels)), 1, NA),
      pch = c(NA, rep(NA, length(category_levels)), NA, 8)
    )
  })
  
  output$acc_text <- renderText({
    req(input$run_acc > 0)
    result <- accuracy_result()
    best_k <- which.max(result$overall)
    sprintf(
      "Selected K = %d: %.1f%% accuracy (%d of %d). Highest catalogue leave-one-out accuracy: %.1f%% at K = %d.",
      input$K, 100 * result$overall[input$K],
      as.integer(round(result$overall[input$K] * nrow(shoes_meta))),
      nrow(shoes_meta), 100 * result$overall[best_k], best_k
    )
  })
  
  # OUR SHOES ------------------------------------------------
  ours <- reactiveValues(
    df = initial_our$df,
    X = initial_our$X,
    message = "",
    preferred = if (nrow(initial_our$df)) initial_our$df$file[1] else ""
  )
  
  observeEvent(input$our_upload, {
    uploads <- input$our_upload
    req(nrow(uploads) > 0)
    new <- empty_our()
    messages <- character(0)
    
    for (j in seq_len(nrow(uploads))) {
      vector <- tryCatch(
        read_external_vector(uploads$datapath[j]),
        error = function(e) NULL
      )
      if (is.null(vector)) {
        messages <- c(messages,
                      paste("Skipped unreadable image:", uploads$name[j]))
        next
      }
      
      details <- parse_our_filename(uploads$name[j])
      new$df <- rbind(new$df, data.frame(
        file = uploads$name[j], name = details$name,
        category = details$category, stringsAsFactors = FALSE
      ))
      new$X <- rbind(new$X, vector)
      if (is.na(details$category)) {
        messages <- c(messages, paste(
          uploads$name[j], "has no valid category in its filename; ground truth is unknown."
        ))
      }
    }
    
    if (nrow(new$df) > 0L) {
      # Replace images with the same filename rather than duplicating them.
      keep <- !(ours$df$file %in% new$df$file)
      ours$df <- rbind(ours$df[keep, , drop = FALSE], new$df)
      ours$X <- rbind(ours$X[keep, , drop = FALSE], new$X)
      ours$preferred <- new$df$file[1]
    }
    ours$message <- paste(
      c(sprintf("Added %d photo(s).", nrow(new$df)), messages),
      collapse = "\n"
    )
  })
  
  observeEvent(input$our_reset, {
    ours$df <- initial_our$df
    ours$X <- initial_our$X
    ours$preferred <- if (nrow(initial_our$df)) initial_our$df$file[1] else ""
    ours$message <- "Removed temporary uploads. Folder photos are preserved."
  })
  
  output$our_upload_msg <- renderUI({
    if (!nzchar(ours$message)) return(NULL)
    p(style = "white-space: pre-line;", ours$message)
  })
  
  output$our_empty_msg <- renderUI({
    if (nrow(ours$df) > 0L) return(NULL)
    p("No personal shoes yet. Upload photos or add JPG/PNG files to the ourshoes folder.")
  })
  
  observe({
    data <- ours$df
    
    # A new installation may have no photos in the ourshoes folder.
    # Do not construct dropdown labels from an empty data frame:
    # paste0() with literal separators can return one label even when
    # data$file has zero elements, making setNames() fail.
    if (nrow(data) == 0L) {
      updateSelectInput(
        session, "our_pick",
        choices = character(0),
        selected = character(0)
      )
      return(invisible(NULL))
    }
    
    labels <- paste0(
      data$name, " - ",
      ifelse(is.na(data$category), "Unknown", data$category),
      " (", data$file, ")"
    )
    choices <- stats::setNames(data$file, labels)
    chosen <- if (ours$preferred %in% data$file) {
      ours$preferred
    } else {
      data$file[1]
    }
    updateSelectInput(session, "our_pick", choices = choices, selected = chosen)
  })
  
  our_index <- reactive({
    req(input$our_pick)
    index <- match(input$our_pick, ours$df$file)
    validate(need(!is.na(index), "Choose a personal shoe first."))
    index
  })
  
  our_result <- reactive({
    find_similar_shoes(as.numeric(ours$X[our_index(), ]), input$K)
  })
  
  output$our_prediction <- renderUI({
    index <- our_index()
    result <- our_result()
    actual <- ours$df$category[index]
    verdict <- if (is.na(actual)) "Unknown" else if (
      tolower(actual) == tolower(result$prediction)
    ) "Correct" else "Incorrect"
    tagList(
      p(strong("Predicted: "), result$prediction),
      p(strong("Actual: "), if (is.na(actual)) "Not available" else actual),
      p(strong("Result: "), verdict)
    )
  })
  
  output$our_best_plot <- renderPlot({
    index <- our_index()
    result <- our_result()
    best_row <- match(result$best_id, shoes_meta$id)
    old <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    plot_shoe(ours$X[index, ], "Selected personal shoe")
    plot_shoe(pca_data$X[best_row, ],
              paste("Closest twin | ID", result$best_id))
  })
  
  output$our_top5_table <- renderTable(our_result()$top5,
                                       striped = TRUE, bordered = TRUE)
  
  output$our_top5_plot <- renderPlot({
    results <- our_result()$top5
    old <- par(mfrow = c(1, NEIGHBORS), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (j in seq_len(NEIGHBORS)) {
      row <- match(results$ID[j], shoes_meta$id)
      plot_shoe(pca_data$X[row, ],
                paste0("#", j, " ID ", results$ID[j],
                       "\n", results$Category[j]))
    }
  })
  
  output$our_grid <- renderPlot({
    validate(need(nrow(ours$df) > 0L, "No personal shoes to display."))
    n <- nrow(ours$df)
    columns <- min(n, 6L)
    rows <- ceiling(n / columns)
    old <- par(mfrow = c(rows, columns), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (i in seq_len(n)) {
      plot_shoe(ours$X[i, ], paste0(ours$df$name[i], "\n",
                                    ifelse(is.na(ours$df$category[i]),
                                           "Unknown", ours$df$category[i])))
    }
  }, height = function() {
    n <- max(nrow(ours$df), 1L)
    200 * ceiling(n / min(n, 6L))
  })
  
  our_all_results <- reactive({
    req(nrow(ours$df) > 0L)
    predictions <- vapply(seq_len(nrow(ours$df)), function(i) {
      find_similar_shoes(as.numeric(ours$X[i, ]), input$K)$prediction
    }, character(1))
    actual <- ours$df$category
    verdict <- ifelse(
      is.na(actual), "Unknown",
      ifelse(tolower(predictions) == tolower(actual), "Correct", "Incorrect")
    )
    data.frame(
      Name = ours$df$name, File = ours$df$file,
      Actual = ifelse(is.na(actual), "Unknown", actual),
      Predicted = predictions, Result = verdict,
      row.names = NULL
    )
  })
  
  output$our_summary <- renderText({
    results <- our_all_results()
    known <- results$Result != "Unknown"
    sprintf(
      "K = %d: %d of %d personal shoes with known categories were classified correctly (%d unknown).",
      input$K, sum(results$Result == "Correct"), sum(known), sum(!known)
    )
  })
  
  output$our_all_table <- renderTable(our_all_results(),
                                      striped = TRUE, bordered = TRUE)
  
  # DATA 1 ---------------------------------------------------
  data1_raw <- reactive({
    req(input$data1_file)
    result <- tryCatch(
      read.csv(input$data1_file$datapath, stringsAsFactors = FALSE,
               check.names = FALSE),
      error = function(e) NULL
    )
    validate(need(!is.null(result) && ncol(result) > 0L,
                  "Please upload a readable CSV dataset."))
    result
  })
  
  observeEvent(data1_raw(), {
    data <- data1_raw()
    cat_columns <- names(data)[vapply(data, function(x) {
      is.character(x) || is.factor(x) || is.logical(x)
    }, logical(1))]
    number_columns <- names(data)[vapply(data, is.numeric, logical(1))]
    updateSelectInput(session, "data1_group",
                      choices = c("No grouping" = "", setNames(names(data), names(data))),
                      selected = "")
    updateSelectInput(session, "data1_filter",
                      choices = c("No filter" = "", setNames(cat_columns, cat_columns)),
                      selected = "")
    updateSelectInput(session, "data1_measure",
                      choices = c("Row count only" = "", setNames(number_columns, number_columns)),
                      selected = "")
  })
  
  observeEvent(list(data1_raw(), input$data1_filter), {
    data <- data1_raw()
    column <- input$data1_filter
    values <- if (!is.null(column) && nzchar(column) && column %in% names(data)) {
      sort(unique(as.character(data[[column]][!is.na(data[[column]])])))
    } else {
      character(0)
    }
    updateSelectInput(session, "data1_value",
                      choices = c("All values" = "", setNames(values, values)),
                      selected = "")
  })
  
  data1_filtered <- reactive({
    data <- data1_raw()
    column <- input$data1_filter
    value <- input$data1_value
    if (!is.null(column) && nzchar(column) && column %in% names(data) &&
        !is.null(value) && nzchar(value)) {
      data <- data[!is.na(data[[column]]) &
                     as.character(data[[column]]) == value, , drop = FALSE]
    }
    data
  })
  
  output$data1_preview <- renderTable({
    head(data1_filtered(), 12)
  }, striped = TRUE, bordered = TRUE)
  
  output$data1_summary <- renderTable({
    data <- data1_filtered()
    group_name <- input$data1_group
    measure <- input$data1_measure
    if (is.null(group_name) || !(group_name %in% names(data))) group_name <- ""
    if (is.null(measure) || !(measure %in% names(data)) ||
        !is.numeric(data[[measure]])) measure <- ""
    
    group_values <- if (nzchar(group_name)) {
      ifelse(is.na(data[[group_name]]), "(Missing)",
             as.character(data[[group_name]]))
    } else {
      rep("All rows", nrow(data))
    }
    levels <- unique(group_values)
    if (!length(levels)) return(data.frame(Message = "No rows match the filter."))
    
    do.call(rbind, lapply(levels, function(level) {
      keep <- group_values == level
      values <- if (nzchar(measure)) data[[measure]][keep] else numeric(0)
      observed <- values[is.finite(values)]
      data.frame(
        Group = level,
        Count = sum(keep),
        Mean = if (nzchar(measure) && length(observed)) round(mean(observed), 3) else NA_real_,
        Median = if (nzchar(measure) && length(observed)) round(median(observed), 3) else NA_real_,
        Minimum = if (nzchar(measure) && length(observed)) min(observed) else NA_real_,
        Maximum = if (nzchar(measure) && length(observed)) max(observed) else NA_real_,
        row.names = NULL
      )
    }))
  }, striped = TRUE, bordered = TRUE)
  
  output$data1_description <- renderText({
    data <- data1_filtered()
    sprintf("The current filter contains %d rows and %d columns. Statistics are descriptive only; they do not establish causation.",
            nrow(data), ncol(data))
  })
  
  # AI DOCUMENTATION -----------------------------------------
  output$ai_table <- renderTable({
    data.frame(
      Tab = c("Shoe Search", "Explain PCA", "Our Shoes", "Data 1"),
      `What AI Did (Draft)` = c(
        "Combined catalogue search, caching, image conversion, and result display.",
        "Combined PCA explanation plots, reconstruction, and leave-one-out evaluation.",
        "Combined multi-upload handling, previews, and prediction summaries.",
        "Implemented optional CSV descriptive summaries."
      ),
      `What the Team Did` = rep(
        "TODO: each team member must verify and describe their actual work.", 4
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)
  
  output$ai_reflection <- renderTable({
    data.frame(
      Member = c("박선우", "송현정", "한우영", "신유진", "Flo"),
      Reflection = rep(
        "TODO: write and verify your own reflection on AI assistance, mistakes, and learning.", 5
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)
}

# 9. START APPLICATION ---------------------------------------
shinyApp(ui = ui, server = server)
