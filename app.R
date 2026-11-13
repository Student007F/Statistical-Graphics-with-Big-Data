# ============================================================
# CHIPMUNK - Shoe Search with Principal Component Analysis
# Shoe PCA, K-pop lyrics analysis, and personal shoe testing.
# ============================================================

# 1. PACKAGES AND SETTINGS -----------------------------------
required_packages <- c("shiny", "bslib", "png", "magick", "dplyr", "ggplot2")
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

# png, magick, and ggplot2 are accessed with explicit package namespaces.
library(shiny)
library(bslib)
library(dplyr)
options(shiny.maxRequestSize = 30 * 1024^2)

IMAGE_SIDE <- 64L
PIXEL_COUNT <- IMAGE_SIDE * IMAGE_SIDE
NEIGHBORS <- 5L
MAX_PCS <- 100L
CACHE_VERSION <- "column-major-png-grayscale-v2"
# Uploaded photos are downscaled to this size (longest side) for cropping.
CROP_WORKING_SIDE <- 1000L
MIN_CROP_SIDE <- 8L

# 2. PROJECT FOLDERS -----------------------------------------
# Catalogue and personal shoe images live inside data/ alongside CSV files.
# Self-check images stay in selfcheck/; generated files stay in generated/.
DATA_DIR <- "data"
CATALOGUE_DIR <- file.path(DATA_DIR, "shoes")
SELFCHECK_DIR <- "selfcheck"
PERSONAL_DIR <- file.path(DATA_DIR, "ourshoes")
CACHE_DIR <- file.path("generated", "cache")
TEMP_DIR <- file.path("generated", "temp")
CATALOGUE_CSV <- file.path(DATA_DIR, "shoes_meta.csv")
SELFCHECK_LABEL_CSV <- file.path(DATA_DIR, "selfcheck_labels.csv")
CACHE_FILE <- file.path(CACHE_DIR, "pca_cache_unified_v2.rds")

# The project already uses data/shoes/ and data/ourshoes/.
# Only folders containing generated files need to be created at startup.
for (folder in c(CACHE_DIR, TEMP_DIR)) {
  if (!dir.exists(folder) &&
      !dir.create(folder, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create generated folder: ", folder)
  }
}

# 3. CATALOGUE DATA ------------------------------------------
if (!file.exists(CATALOGUE_CSV)) {
  stop("Missing ", CATALOGUE_CSV,
       ". Place shoes_meta.csv in the data folder next to app.R.")
}

shoes_meta <- read.csv(CATALOGUE_CSV, stringsAsFactors = FALSE)
if (!all(c("id", "category") %in% names(shoes_meta))) {
  stop("data/shoes_meta.csv must contain columns named 'id' and 'category'.")
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

catalogue_paths <- file.path(CATALOGUE_DIR, sprintf("%04d.png", shoes_meta$id))
missing_images <- catalogue_paths[!file.exists(catalogue_paths)]
if (length(missing_images) > 0) {
  stop("Catalogue images are missing. First missing file: ", missing_images[1])
}

category_levels <- sort(unique(shoes_meta$category))
category_colors <- setNames(
  grDevices::hcl.colors(length(category_levels), palette = "Dark 3"),
  category_levels
)

# 4. IMAGE PREPROCESSING -------------------------------------
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

# Apply camera orientation and replace transparency with a white background.
prepare_external_image <- function(img) {
  img <- magick::image_orient(img)
  magick::image_background(img, "white", flatten = TRUE)
}

# Pad a prepared magick image to a square without stretching, resize it to
# 64 x 64, convert to grayscale, and decode it with the same PNG reader
# used for the catalogue.
image_to_vector <- function(img) {
  img <- magick::image_repage(img)
  dims <- magick::image_info(img)
  side <- max(dims$width[1], dims$height[1])
  img <- magick::image_extent(
    img, geometry = paste0(side, "x", side),
    gravity = "center", color = "white"
  )
  img <- magick::image_resize(img, "64x64!")
  img <- magick::image_convert(img, colorspace = "gray")
  
  temporary_png <- tempfile(pattern = "preprocess_", tmpdir = TEMP_DIR, fileext = ".png")
  on.exit(unlink(temporary_png), add = TRUE)
  magick::image_write(img, path = temporary_png, format = "png", depth = 8)
  read_png_vector(temporary_png)
}

read_external_vector <- function(path) {
  # Convert real photographs in a consistent way: orient, flatten alpha,
  # pad to a square without stretching, resize and convert to grayscale.
  img <- magick::image_read(path)[1]
  dims <- magick::image_info(img)
  
  # An existing 64x64 PNG should remain pixel-exact. This also works for
  # Shiny temporary uploads whose paths have no .png file extension.
  if (tolower(dims$format[1]) == "png" &&
      dims$width[1] == IMAGE_SIDE && dims$height[1] == IMAGE_SIDE) {
    exact_vector <- tryCatch(read_png_vector(path), error = function(e) NULL)
    if (!is.null(exact_vector)) return(exact_vector)
  }
  
  image_to_vector(prepare_external_image(img))
}

# Convert a Shiny brush (plot coordinates, origin at bottom-left) into a
# crop region in image pixels (origin at top-left).
brush_to_crop <- function(brush, width, height) {
  x0 <- max(0, floor(brush$xmin))
  x1 <- min(width, ceiling(brush$xmax))
  top <- max(0, floor(height - brush$ymax))
  bottom <- min(height, ceiling(height - brush$ymin))
  crop_width <- x1 - x0
  crop_height <- bottom - top
  if (crop_width < MIN_CROP_SIDE || crop_height < MIN_CROP_SIDE) return(NULL)
  list(
    x = x0, y = top, width = crop_width, height = crop_height,
    geometry = sprintf("%dx%d+%d+%d", crop_width, crop_height, x0, top)
  )
}

save_vector_png <- function(vector) {
  stopifnot(length(vector) == PIXEL_COUNT)
  values <- pmax(0, pmin(255, as.numeric(vector))) / 255
  image_matrix <- matrix(values, nrow = IMAGE_SIDE, ncol = IMAGE_SIDE)
  destination <- tempfile(pattern = "preview_", tmpdir = TEMP_DIR, fileext = ".png")
  png::writePNG(image_matrix, destination)
  destination
}

# Preload optional self-check ground-truth labels once at startup.
selfcheck_labels_df <- if (file.exists(SELFCHECK_LABEL_CSV)) {
  tryCatch(read.csv(SELFCHECK_LABEL_CSV, stringsAsFactors = FALSE),
           error = function(e) NULL)
} else {
  NULL
}

# Match a self-check filename to its optional ground-truth label.
read_selfcheck_label <- function(filename) {
  labels <- selfcheck_labels_df
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

# 5. PCA WITH A VERSIONED, VALIDATED CACHE --------------------
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

# 6. MATCHING, RECONSTRUCTION, AND ACCURACY -------------------
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

# Draw the five nearest matches for both Shoe Search and Our Shoes.
plot_neighbor_gallery <- function(top5) {
  old <- par(mfrow = c(1, NEIGHBORS), mar = c(1, 1, 3, 1))
  on.exit(par(old))
  for (j in seq_len(NEIGHBORS)) {
    row <- match(top5$ID[j], shoes_meta$id)
    plot_shoe(pca_data$X[row, ],
              paste0("#", j, " ID ", top5$ID[j], "\n", top5$Category[j]))
  }
}

# Optional, expensive computation: leave one catalogue image out, then
# classify it using its five closest remaining neighbours at each K
# within the precomputed PCA coordinate space.
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

# 7. PERSONAL SHOE DATA --------------------------------------
parse_our_filename <- function(filename) {
  stem <- tools::file_path_sans_ext(basename(filename))
  # Match category names literally against catalogue categories. This handles
  # underscores and regex-special characters in categories without escaping.
  # Check the longest category first to avoid matching a shorter suffix.
  candidates <- unique(c(stem, sub("_[0-9]+$", "", stem)))
  categories <- category_levels[order(nchar(category_levels), decreasing = TRUE)]
  for (category in categories) {
    suffix <- paste0("_", tolower(category))
    for (candidate in candidates) {
      if (endsWith(tolower(candidate), suffix)) {
        member <- substr(candidate, 1L, nchar(candidate) - nchar(suffix))
        if (nzchar(member)) return(list(name = member, category = category))
      }
    }
  }
  # Keep older, unrecognized photos visible, but never invent a category.
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  member <- if (length(parts) && nzchar(parts[1])) parts[1] else "Unknown"
  list(name = member, category = NA_character_)
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
  if (!dir.exists(PERSONAL_DIR)) return(empty_our())
  filenames <- sort(list.files(
    PERSONAL_DIR, pattern = "\\.(png|jpg|jpeg)$", ignore.case = TRUE
  ))
  if (length(filenames) == 0L) return(empty_our())
  
  rows <- vector("list", length(filenames))
  vectors <- vector("list", length(filenames))
  keep <- 0L
  
  for (filename in filenames) {
    vector <- tryCatch(
      read_external_vector(file.path(PERSONAL_DIR, filename)),
      error = function(e) {
        warning("Skipping unreadable personal shoe photo: ", filename)
        NULL
      }
    )
    if (is.null(vector)) next
    details <- parse_our_filename(filename)
    keep <- keep + 1L
    rows[[keep]] <- data.frame(
      file = filename, name = details$name,
      category = details$category, stringsAsFactors = FALSE
    )
    vectors[[keep]] <- vector
  }
  if (keep == 0L) return(empty_our())
  list(
    df = do.call(rbind, rows[seq_len(keep)]),
    X = do.call(rbind, vectors[seq_len(keep)])
  )
}

selfcheck_files <- if (dir.exists(SELFCHECK_DIR)) {
  sort(list.files(SELFCHECK_DIR, pattern = "\\.(png|jpg|jpeg)$",
                  ignore.case = TRUE))
} else {
  character(0)
}

# 8. K-POP LYRICS DATA FOR DATA 1 -----------------------------
# The source CSV is optional during startup. Missing or invalid data is
# reported inside Data 1, rather than blocking the entire shoe application.
LYRICS_FILE <- file.path(DATA_DIR, "kpop_lyrics_counts.csv")
LYRICS_REQUIRED <- c("year", "is_hive", "pct_main", "pct_sub_1", "pct_sub_2")
LYRICS_RULE_COLUMNS <- c(
  main = "pct_main", sub_1 = "pct_sub_1", sub_2 = "pct_sub_2"
)
LYRICS_RULE_CHOICES <- c(
  "Main: all English words" = "main",
  "Sub 1: exclude English fillers (ad-libs)" = "sub_1",
  "Sub 2: repeated words count once per song" = "sub_2"
)

load_lyrics_csv <- function(path) {
  if (!file.exists(path)) {
    return(list(data = NULL,
                problem = paste0("Missing ", path, ". Add the CSV and restart the app.")))
  }
  
  data <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(data)) {
    return(list(data = NULL, problem = paste0("Cannot read ", path, ".")))
  }
  absent <- setdiff(LYRICS_REQUIRED, names(data))
  if (length(absent)) {
    return(list(data = NULL,
                problem = paste("K-pop CSV is missing columns:",
                                paste(absent, collapse = ", "))))
  }
  
  data$year <- suppressWarnings(as.integer(as.character(data$year)))
  for (column in unname(LYRICS_RULE_COLUMNS)) {
    # Accept both 32.5 and "32.5%", always treating the values as 0-100.
    values <- gsub("%", "", as.character(data[[column]]), fixed = TRUE)
    data[[column]] <- suppressWarnings(as.numeric(values))
    invalid <- !is.na(data[[column]]) &
      (!is.finite(data[[column]]) | data[[column]] < 0 | data[[column]] > 100)
    data[[column]][invalid] <- NA_real_
  }
  flags <- tolower(trimws(as.character(data$is_hive)))
  data$is_hive <- ifelse(flags %in% c("true", "1", "yes"), "true",
                         ifelse(flags %in% c("false", "0", "no"), "false", NA_character_))
  data <- data[!is.na(data$year) & !is.na(data$pct_main), , drop = FALSE]
  if (nrow(data) == 0L) {
    return(list(data = NULL, problem = "No rows with valid year and pct_main values."))
  }
  list(data = data, problem = NULL)
}

lyrics_source <- load_lyrics_csv(LYRICS_FILE)
lyrics_data <- lyrics_source$data
lyrics_first_year <- if (is.null(lyrics_data)) 2014L else min(lyrics_data$year)
lyrics_last_year <- if (is.null(lyrics_data)) 2025L else max(lyrics_data$year)

# Convert selected counting rules to a consistent long-form table.
# Songs with missing data for one rule are omitted from that rule only.
lyrics_make_long <- function(data, rules, split_groups = FALSE) {
  stopifnot(all(rules %in% names(LYRICS_RULE_COLUMNS)))
  if (split_groups) {
    data <- data[!is.na(data$is_hive), , drop = FALSE]
    data$group <- ifelse(data$is_hive == "true", "HYBE", "Non-HYBE")
  } else {
    data$group <- rep("All songs", nrow(data))
  }
  
  long <- dplyr::bind_rows(lapply(rules, function(rule) {
    piece <- data[!is.na(data[[LYRICS_RULE_COLUMNS[[rule]]]]), , drop = FALSE]
    piece$rule <- rep(rule, nrow(piece))
    piece$pct <- piece[[LYRICS_RULE_COLUMNS[[rule]]]]
    piece
  }))
  if (length(rules) > 1L && split_groups) {
    long$series <- paste(long$rule, long$group, sep = " / ")
  } else if (length(rules) > 1L) {
    long$series <- long$rule
  } else {
    long$series <- long$group
  }
  long$series <- factor(long$series, levels = unique(long$series))
  long
}

# Shared statistics for Findings A/B and the optional HYBE comparison.
summarize_lyrics <- function(data, threshold = NULL) {
  if (is.null(threshold)) {
    return(data %>%
             dplyr::group_by(series) %>%
             dplyr::summarise(n = dplyr::n(), mean_pct = mean(pct),
                              median_pct = median(pct), .groups = "drop"))
  }
  data %>%
    dplyr::group_by(series) %>%
    dplyr::summarise(n = dplyr::n(), mean_pct = mean(pct),
                     median_pct = median(pct),
                     share_thr = mean(pct >= threshold) * 100,
                     .groups = "drop")
}

# Compact KPI cards use the same bslib/Bootstrap visual language as the other tabs.
lyrics_stat_cards <- function(summary_table, threshold = NULL) {
  if (!nrow(summary_table)) return(p("No data matches the current selection."))
  tagList(lapply(seq_len(nrow(summary_table)), function(i) {
    values <- list(
      c("Number of songs", format(summary_table$n[i], big.mark = ",")),
      c("Mean English share", sprintf("%.1f%%", summary_table$mean_pct[i])),
      c("Median English share", sprintf("%.1f%%", summary_table$median_pct[i]))
    )
    if (!is.null(threshold)) {
      values[[4]] <- c(
        sprintf("Songs with at least %d%% English", threshold),
        sprintf("%.1f%%", summary_table$share_thr[i])
      )
    }
    div(
      class = "d1-stat-series",
      if (nrow(summary_table) > 1L) h5(as.character(summary_table$series[i])),
      div(class = "d1-stat-grid", lapply(values, function(value) {
        div(class = "d1-stat-tile",
            div(class = "d1-stat-title", value[1]),
            div(class = "d1-stat-value", value[2]))
      }))
    )
  }))
}

# Show data-dependent endpoint comparisons without calling them trends.
# Each series uses its first and last available year, not extrapolated values.
lyrics_year_change_notes <- function(yearly, column, measure_label) {
  # Keep the series order used in the chart legend.
  series_names <- if (is.factor(yearly$series)) {
    levels(droplevels(yearly$series))
  } else {
    unique(as.character(yearly$series))
  }

  cards <- lapply(series_names, function(group_name) {
    rows <- yearly[as.character(yearly$series) == group_name, , drop = FALSE]
    rows <- rows[order(rows$year), , drop = FALSE]

    if (nrow(rows) < 2L) {
      return(div(
        class = "d1-change-card",
        div(class = "d1-change-series", group_name),
        div(class = "d1-change-note",
            "Only one year is available for the current selection.")
      ))
    }

    first <- rows[[column]][1]
    last <- rows[[column]][nrow(rows)]
    change <- last - first
    direction <- if (abs(change) < 0.05) "flat" else if (change > 0) "up" else "down"
    arrow <- c(up = "▲", down = "▼", flat = "■")[[direction]]

    div(
      class = "d1-change-card",
      div(class = "d1-change-series", group_name),
      div(class = "d1-change-measure", measure_label),
      div(
        class = "d1-change-years",
        div(class = "d1-change-point",
            div(class = "d1-change-year", rows$year[1]),
            div(class = "d1-change-value", sprintf("%.1f%%", first))),
        div(class = "d1-change-arrow", "→"),
        div(class = "d1-change-point",
            div(class = "d1-change-year", rows$year[nrow(rows)]),
            div(class = "d1-change-value", sprintf("%.1f%%", last)))
      ),
      div(
        class = paste("d1-change-delta", direction),
        sprintf("%s %+.1f percentage points", arrow, change)
      )
    )
  })
  div(class = "d1-change-grid", cards)
}

# Shared theme for Data 1 plots: visually consistent with Flatly cards.
lyrics_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", color = "#2C3E50"),
      axis.title = ggplot2::element_text(color = "#34495E"),
      legend.position = "bottom"
    )
}

# Avoid showing non-working charts when the required K-pop CSV is missing.
lyrics_unavailable_card <- function() {
  card(
    card_header("K-pop lyrics dataset not available"),
    p(lyrics_source$problem),
    p("Expected path: ", tags$code(LYRICS_FILE)),
    p("Expected columns: ", tags$code(paste(LYRICS_REQUIRED, collapse = ", "))),
    p("The Shoe Search, Explain PCA, and Our Shoes sections still work.")
  )
}

# 9. USER INTERFACE ------------------------------------------
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
    .d1-main-content { min-width: 0; width: 100%; }
    .d1-main-content .card { margin-bottom: 16px; }
    .d1-stat-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(155px, 1fr));
      gap: 12px;
      margin: 10px 0 16px;
    }
    .d1-stat-series { margin-bottom: 8px; }
    .d1-stat-series h5 { color: #2c3e50; font-weight: 650; margin-top: 14px; }
    .d1-stat-tile {
      border: 1px solid #dce6ed;
      border-radius: 10px;
      background: #f7fafc;
      padding: 14px 16px;
    }
    .d1-stat-title { color: #617383; font-size: 0.87rem; }
    .d1-stat-value { color: #2c3e50; font-size: 1.5rem; font-weight: 700; }
    .d1-method-details {
      margin-top: 14px;
      padding: 14px;
      border: 1px solid #dce6ed;
      border-radius: 10px;
    }
    .d1-method-details summary { cursor: pointer; font-weight: 650; }
    .d1-method-details p { margin-top: 12px; }
    .d1-method-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(155px, 1fr));
      gap: 12px;
      margin-top: 14px;
    }
    .d1-method-grid > div {
      padding: 12px;
      background: #f7fafc;
      border-radius: 8px;
    }
    .d1-example-table { margin-top: 12px; }
    .d1-main-content .shiny-plot-output { max-width: 100%; }
    @media (max-width: 767px) {
      .d1-stat-grid { grid-template-columns: 1fr 1fr; }
      .d1-stat-value { font-size: 1.25rem; }
    }
    /* Explain PCA Step 4: colour legend for eigenshoe loadings. */
    .pca-color-legend { display: grid; gap: 6px; margin: 4px 0 8px; font-size: 0.93rem; }
    .pca-swatch {
      display: inline-block; width: 14px; height: 14px; margin-right: 6px;
      border: 1px solid #b8c4cc; border-radius: 3px; vertical-align: -2px;
    }
    /* Data 1: value tooltip for points on the yearly line chart. */
    .d1-plot-wrap { position: relative; }
    .d1-plot-tooltip {
      position: absolute; z-index: 10; pointer-events: none;
      transform: translate(-50%, calc(-100% - 12px));
      background: rgba(44, 62, 80, 0.94); color: #fff;
      padding: 6px 10px; border-radius: 6px;
      font-size: 0.85rem; white-space: nowrap;
      box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2);
    }
    /* Data 1: first-year vs last-year comparison cards. */
    .d1-change-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 12px;
      margin: 6px 0 12px;
    }
    .d1-change-card {
      border: 1px solid #dce6ed;
      border-radius: 10px;
      background: #f7fafc;
      padding: 14px 16px;
    }
    .d1-change-series { color: #2c3e50; font-weight: 700; font-size: 1.02rem; }
    .d1-change-measure { color: #617383; font-size: 0.85rem; margin-bottom: 10px; }
    .d1-change-years { display: flex; align-items: center; gap: 14px; }
    .d1-change-year { color: #617383; font-size: 0.85rem; }
    .d1-change-value { color: #2c3e50; font-size: 1.35rem; font-weight: 700; }
    .d1-change-arrow { color: #93a7b4; font-size: 1.3rem; }
    .d1-change-delta {
      display: inline-block;
      margin-top: 10px;
      padding: 3px 10px;
      border-radius: 999px;
      font-size: 0.87rem;
      font-weight: 650;
    }
    .d1-change-delta.up { background: #e3f6f1; color: #0f7a63; }
    .d1-change-delta.down { background: #fbe9e7; color: #b03a2e; }
    .d1-change-delta.flat { background: #edf2f5; color: #576574; }
    .d1-change-note { color: #576574; font-size: 0.93rem; margin-top: 6px; }
    /* Shoe Search: manual crop controls for uploaded photos. */
    .crop-buttons { display: flex; flex-wrap: wrap; gap: 8px; margin: 8px 0; }
    .crop-status { color: #576574; font-size: 0.93rem; }
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
                         min = min(shoes_meta$id), max = max(shoes_meta$id),
                         step = 1)
          ),
          conditionalPanel(
            "input.src == 'Self-check'",
            selectInput("sc_file", "Self-check Photo", choices = selfcheck_files),
            p(class = "shoe-caption", "If data/selfcheck_labels.csv is available, the true category is displayed.")
          ),
          conditionalPanel(
            "input.src == 'Upload'",
            fileInput("photo", "Upload Photo (JPG/PNG)",
                      accept = c(".jpg", ".jpeg", ".png")),
            p(class = "shoe-caption",
              "After uploading, crop the photo around the shoe in the panel on the right for better matches.")
          ),
          hr(),
          sliderInput("K", "Number of Principal Components (K)",
                      min = 1, max = K_MAX, value = min(20L, K_MAX), step = 1),
          sliderInput("image_size", "Image Size", min = 96, max = 240,
                      value = 160, step = 8, post = " px"),
          uiOutput("image_size_css")
        ),
        conditionalPanel(
          "input.src == 'Upload'",
          card(
            fill = FALSE,
            card_header("Crop the Uploaded Photo"),
            card_body(
              fill = FALSE, fillable = FALSE,
              p(class = "shoe-caption",
                "Drag a box tightly around the shoe, then click 'Apply crop'. Cropping removes most of the background, so the PCA focuses on the shoe itself. A plain, light background inside the box gives the best results."),
              plotOutput("upload_crop_plot", height = "360px",
                         brush = brushOpts("upload_brush", resetOnNew = TRUE,
                                           fill = "#18BC9C", stroke = "#18BC9C",
                                           opacity = 0.2)),
              div(class = "crop-buttons",
                  actionButton("apply_crop", "Apply crop", class = "btn-primary"),
                  actionButton("reset_crop", "Use full photo",
                               class = "btn-outline-secondary")),
              div(class = "crop-status", textOutput("crop_status"))
            )
          )
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
             p("Blue indicates darker than the average shoe (lower pixel value); red indicates brighter than the average shoe (higher pixel value)."),
             plotOutput("pca_step2", height = "250px"))
      ),
      card(card_header("Step 3 | PCA score map (PC1 and PC2)"),
           p("Catalogue points are colored by category; your selected photo is highlighted."),
           plotOutput("pca_step3", height = "410px")),
      card(card_header("Step 4 | First 10 eigenshoes"),
           p("Principal component loadings drawn as 64 x 64 images, with their individual explained variance."),
           div(
             class = "pca-color-legend",
             div(span(class = "pca-swatch", style = "background:#C64342;"),
                 strong("Red (positive weight): "),
                 "a shoe with a high score on this PC is brighter than the average shoe in these pixels."),
             div(span(class = "pca-swatch", style = "background:#245A9C;"),
                 strong("Blue (negative weight): "),
                 "a shoe with a high score on this PC is darker than the average shoe in these pixels."),
             div(span(class = "pca-swatch", style = "background:#FFFFFF;"),
                 strong("White (near zero): "),
                 "these pixels barely affect this PC.")
           ),
           p(class = "shoe-caption",
             "For a negative score the effect is reversed. The sign of a PC is arbitrary: what matters is that red and blue regions change in opposite directions, so each eigenshoe describes a contrast between them. Stronger colors mean a larger influence."),
           plotOutput("pca_step4", height = "440px")),
      card(card_header("Step 5 | Query scores"),
           p("The first K PCA coordinates used to calculate nearest-neighbor distances."),
           tableOutput("pca_step5_table")),
      card(card_header("Step 6 | Reconstruct using different K values"),
           p("Compare the same photo reconstructed using different numbers of components."),
           plotOutput("pca_step6", height = "300px")),
      card(card_header("Classification accuracy versus K"),
           p("Leave-one-out k-NN evaluation within the catalogue PCA space: each photo is classified using the five nearest other photos. Computing all K values can be slow."),
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
      p("Can the app recognize real shoes it has never seen? Explore all team member photos (3 per member) and test individual predictions."),
      uiOutput("our_empty_msg"),
      card(
        card_header("All Our Photos (3 per Member)"),
        plotOutput("our_grid", height = "auto")
      ),
      card(
        card_header("Analyze a Selected Shoe"),
        selectInput("our_pick", "Choose a Shoe", choices = character(0), width = "100%"),
        uiOutput("our_prediction"),
        plotOutput("our_best_plot", height = "300px"),
        h5("Five Closest Catalogue Photos"),
        tableOutput("our_top5_table"),
        plotOutput("our_top5_plot", height = "220px")
      ),
      card(
        card_header("Overall Results at the Selected K"),
        textOutput("our_summary"),
        tableOutput("our_all_table")
      )
    )
  ),
  
  nav_panel(
    "Data 1",
    scroll_tab(
      layout_sidebar(
        fill = FALSE, fillable = FALSE,
        sidebar = sidebar(
          title = "Data 1 Controls", width = 320,
          radioButtons(
            "d1_view", "Select Analysis",
            choices = c("Finding A", "Finding B"),
            selected = "Finding A"
          ),
          hr(),
          conditionalPanel(
            "input.d1_view == 'Finding A'",
            sliderInput(
              "year_range", "Select Year Range",
              min = lyrics_first_year, max = lyrics_last_year,
              value = c(lyrics_first_year, lyrics_last_year),
              step = 1, sep = ""
            ),
            checkboxGroupInput(
              "a_rules", "Counting Rules (one or more)",
              choices = LYRICS_RULE_CHOICES, selected = "main"
            ),
            checkboxInput("a_hive", "Split by HYBE / Non-HYBE", FALSE),
            checkboxInput("a_show_b", "Show detailed HYBE comparison", FALSE),
            conditionalPanel(
              "input.a_show_b",
              selectInput(
                "hive_group", "Group for detailed comparison",
                choices = c("Both groups" = "all", "HYBE" = "true",
                            "Non-HYBE" = "false"),
                selected = "all"
              )
            ),
            p(class = "shoe-caption",
              "A song with an unknown HYBE affiliation is excluded when comparing groups.")
          ),
          conditionalPanel(
            "input.d1_view == 'Finding B'",
            sliderInput(
              "c_thr", "Minimum English Share",
              min = 10, max = 90, value = 50, step = 5, post = "%"
            ),
            checkboxInput("c_hive", "Split by HYBE / Non-HYBE", FALSE),
            p(class = "shoe-caption",
              "Finding B always uses the main counting rule and all available years.")
          )
        ),
        div(
          class = "d1-main-content",
          conditionalPanel(
            "input.d1_view == 'Finding A'",
            if (is.null(lyrics_data)) lyrics_unavailable_card() else tagList(
              card(
                card_header("Finding A | English in K-pop lyrics over time"),
                p("How has the proportion of English words in K-pop lyrics changed across years?"),
                tags$details(
                  class = "d1-method-details",
                  tags$summary("How are English percentages calculated?"),
                  p("For each song, English share = English words / (Korean words + English words) × 100."),
                  div(class = "d1-method-grid",
                      div(strong("Main"), p("All English words, including fillers and repetitions.")),
                      div(strong("Sub 1"), p("English ad-libs such as 'oh' and 'yeah' are excluded.")),
                      div(strong("Sub 2"), p("Each distinct word is counted once per song."))
                  ),
                  p(strong("Illustrative example (not a dataset observation)")),
                  tags$code("baby baby baby 사랑해 yeah yeah 사랑해"),
                  tags$table(
                    class = "table table-sm table-striped d1-example-table",
                    tags$thead(tags$tr(tags$th("Rule"), tags$th("English"),
                                       tags$th("Korean"), tags$th("English share"))),
                    tags$tbody(
                      tags$tr(tags$td("Main"), tags$td("5"), tags$td("2"), tags$td("71.4%")),
                      tags$tr(tags$td("Sub 1"), tags$td("3"), tags$td("2"), tags$td("60.0%")),
                      tags$tr(tags$td("Sub 2"), tags$td("2"), tags$td("1"), tags$td("66.7%"))
                    )
                  )
                )
              ),
              card(
                card_header("Summary Statistics"),
                uiOutput("a_stats")
              ),
              card(
                card_header("Average English Share by Year"),
                div(
                  class = "d1-plot-wrap",
                  plotOutput("year_plot", height = "410px",
                             click = "year_plot_click",
                             hover = hoverOpts("year_plot_hover", delay = 80,
                                               delayType = "debounce")),
                  uiOutput("year_plot_tooltip")
                ),
                p(class = "shoe-caption", "Each point shows the unweighted mean of song-level English percentages for that year and selected group/rule. Hover over or tap a point to see its value.")
              ),
              card(
                card_header("Observed Results for Selected Filters"),
                uiOutput("a_observation"),
                p(class = "shoe-caption",
                  "These comparisons use only the first and last available year in each selected series; they do not establish a continuous trend or causation.")
              ),
              card(
                card_header("Why It Matters and Next Analysis"),
                p(strong("Why it matters: "), "Comparing counting rules can show how fillers and repeated words affect the measured English share. This dataset describes lyric composition, not whether English lyrics cause greater popularity."),
                p(strong("How this helps our next analysis: "), "Compare the selected series before deciding whether their patterns differ. Follow up by checking the role of repetition and ad-libs, and use a separate dataset with chart performance metrics to investigate any association with popularity.")
              ),
              conditionalPanel(
                "input.a_show_b",
                card(
                  card_header("HYBE vs Non-HYBE | Detailed Comparison"),
                  p("Compare average English usage across the two company-affiliation groups within the selected year range."),
                  uiOutput("hive_stats"),
                  plotOutput("hive_plot", height = "380px"),
                  p(class = "shoe-caption",
                    "The selected year range and group in the sidebar determine which songs appear here (using the main counting rule). These are descriptive differences, not evidence of a causal effect.")
                )
              )
            )
          ),
          conditionalPanel(
            "input.d1_view == 'Finding B'",
            if (is.null(lyrics_data)) lyrics_unavailable_card() else tagList(
              card(
                card_header("Finding B | How many songs are English-centered?"),
                p("What proportion of K-pop songs reach or exceed a given minimum English lyric threshold across years?")
              ),
              card(card_header("Summary Statistics"), uiOutput("c_stats")),
              card(
                card_header("English-Centered Songs by Year"),
                plotOutput("c_plot", height = "420px"),
                p(class = "shoe-caption",
                  "For each year, the percentage is the number of songs meeting the threshold divided by the number of songs with valid lyrics data in that year.")
              ),
              card(
                card_header("Observed Results for Selected Threshold"),
                uiOutput("c_observation"),
                p(class = "shoe-caption",
                  "The change compares the first and last available year for each displayed group; it does not establish a continuous trend or causation.")
              ),
              card(
                card_header("Why It Matters and Next Analysis"),
                p(strong("Why it matters: "), "Average English share alone does not reveal how many sampled songs exceed a chosen threshold. The threshold-based percentage describes that proportion within the dataset, not all K-pop releases."),
                p(strong("How this helps our next analysis: "), "Use the selected threshold and year-by-year results to identify releases worth investigating. If the comparison reveals differences between years or affiliation groups, follow up with suitable data on their target markets, release formats, or artists.")
              )
            )
          )
        )
      )
    )
  ),
  
  nav_panel(
    "AI",
    scroll_tab(
      card(
        card_header("Ethics Statement"),
        p("We understand the rules on data and AI use for this project, including following the license and privacy rules of every dataset we used, understanding all of our code, and knowing that presenting AI output as our own work may be academic misconduct. We certify that we have followed all of them. Any violation is our responsibility, and we agree to receive an F in this course if we have not."),
        p(strong("Signed: 박선우, 송현정, 한우영, 신유진, Flo"))
      ),
      card(
        card_header("What AI Did and What the Team Did"),
        p("Overview of how AI tools and team members contributed to each section."),
        tableOutput("ai_table")
      ),
      card(
        card_header("Individual Reflections"),
        p("Personal reflections from each team member on AI assistance, challenges, and learning."),
        tableOutput("ai_reflection")
      )
    )
  )
)

# 10. REACTIVE SERVER -----------------------------------------
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
  
  # UPLOAD CROPPING --------------------------------------------
  # The uploaded photo is oriented, flattened onto white and downscaled once.
  # The same working image is used for display and for cropping, so the
  # brush coordinates always match the pixels that are cropped.
  upload_image <- reactive({
    req(input$photo)
    img <- tryCatch({
      working <- prepare_external_image(magick::image_read(input$photo$datapath)[1])
      magick::image_resize(
        working, paste0(CROP_WORKING_SIDE, "x", CROP_WORKING_SIDE, ">")
      )
    }, error = function(e) NULL)
    validate(need(!is.null(img), "Could not read the uploaded image."))
    img
  })
  
  upload_crop <- reactiveVal(NULL)
  
  # A new photo always starts uncropped.
  observeEvent(input$photo, upload_crop(NULL))
  
  observeEvent(input$apply_crop, {
    brush <- input$upload_brush
    if (is.null(brush)) {
      showNotification("Drag a box around the shoe first.", type = "warning")
      return()
    }
    dims <- magick::image_info(upload_image())
    region <- brush_to_crop(brush, dims$width[1], dims$height[1])
    if (is.null(region)) {
      showNotification("The selected area is too small. Drag a larger box.",
                       type = "warning")
      return()
    }
    upload_crop(region)
  })
  
  observeEvent(input$reset_crop, upload_crop(NULL))
  
  output$upload_crop_plot <- renderPlot({
    validate(need(input$photo, "Upload a photo in the sidebar to start cropping."))
    img <- upload_image()
    dims <- magick::image_info(img)
    width <- dims$width[1]
    height <- dims$height[1]
    
    old <- par(mar = c(0, 0, 0, 0))
    on.exit(par(old))
    plot(NA, xlim = c(0, width), ylim = c(0, height), asp = 1,
         xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    rasterImage(as.raster(img), 0, 0, width, height)
    
    region <- upload_crop()
    if (!is.null(region)) {
      rect(region$x, height - region$y - region$height,
           region$x + region$width, height - region$y,
           border = "#E74C3C", lwd = 3)
    }
  })
  
  output$crop_status <- renderText({
    req(input$photo)
    region <- upload_crop()
    if (is.null(region)) {
      "Currently using the full photo."
    } else {
      sprintf("Currently using the cropped area (red box, %d x %d px).",
              region$width, region$height)
    }
  })
  
  # One query representation is shared by both the Search and Explain PCA tabs.
  query_data <- reactive({
    req(input$src)
    
    if (identical(input$src, "Catalogue")) {
      req(input$cat_id)
      requested_id <- suppressWarnings(as.numeric(input$cat_id))
      validate(need(length(requested_id) == 1L &&
                      is.finite(requested_id) &&
                      requested_id == trunc(requested_id),
                    "Catalogue ID must be a whole number."))
      index <- match(requested_id, shoes_meta$id)
      validate(need(!is.na(index), "Catalogue ID not found."))
      return(list(
        vector = as.numeric(pca_data$X[index, ]),
        actual = shoes_meta$category[index],
        exclude_id = shoes_meta$id[index]
      ))
    }
    
    if (identical(input$src, "Self-check")) {
      validate(need(length(selfcheck_files) > 0L,
                    "No self-check photos found in selfcheck/."))
      req(input$sc_file)
      path <- file.path(SELFCHECK_DIR, basename(input$sc_file))
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
    region <- upload_crop()
    vector <- if (is.null(region)) {
      tryCatch(read_external_vector(input$photo$datapath),
               error = function(e) NULL)
    } else {
      tryCatch(
        image_to_vector(magick::image_crop(upload_image(), region$geometry)),
        error = function(e) NULL
      )
    }
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
    list(src = file.path(CATALOGUE_DIR, sprintf("%04d.png", id)),
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
    plot_neighbor_gallery(shoe_search()$top5)
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
  # Load saved photos at the start of each browser session, so restarting
  # the app or opening a new session picks up everything on disk.
  ours <- load_folder_shoes()
  
  output$our_empty_msg <- renderUI({
    if (nrow(ours$df) > 0L) return(NULL)
    p("No personal shoes found in data/ourshoes/.")
  })
  
  observe({
    data <- ours$df
    
    # An empty folder must produce zero choices, not one label.
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
    updateSelectInput(session, "our_pick", choices = choices, selected = data$file[1])
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
    plot_neighbor_gallery(our_result()$top5)
  })
  
  # Display in a 3-column grid (3 shoes per member row).
  output$our_grid <- renderPlot({
    validate(need(nrow(ours$df) > 0L, "No personal shoes to display."))
    ord <- order(tolower(ours$df$name), ours$df$file)
    n <- length(ord)
    columns <- min(n, 3L)
    rows <- ceiling(n / columns)
    old <- par(mfrow = c(rows, columns), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (i in ord) {
      plot_shoe(ours$X[i, ], paste0(ours$df$name[i], "\n",
                                    ifelse(is.na(ours$df$category[i]),
                                           "Unknown", ours$df$category[i])))
    }
  }, height = function() {
    n <- max(nrow(ours$df), 1L)
    210 * ceiling(n / min(n, 3L))
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
  
  # DATA 1: K-POP LYRICS FINDINGS A AND B ----------------------
  # The K-pop CSV is loaded once on startup. If unavailable, Data 1
  # displays a clear instruction, while all other tabs keep working.
  lyrics_required <- reactive({
    req(!is.null(lyrics_data))
    lyrics_data
  })
  
  a_long <- reactive({
    data <- lyrics_required()
    validate(need(length(input$a_rules) > 0L,
                  "Select at least one counting rule in the sidebar."))
    req(input$year_range)
    data <- data[data$year >= input$year_range[1] &
                   data$year <= input$year_range[2], , drop = FALSE]
    long <- lyrics_make_long(data, input$a_rules, isTRUE(input$a_hive))
    validate(need(nrow(long) > 0L, "No valid songs match the selected options."))
    long
  })
  
  output$a_stats <- renderUI({
    long <- a_long()
    lyrics_stat_cards(summarize_lyrics(long))
  })
  
  a_yearly <- reactive({
    a_long() %>%
      dplyr::group_by(year, series) %>%
      dplyr::summarise(mean_pct = mean(pct), .groups = "drop")
  })
  
  output$a_observation <- renderUI({
    lyrics_year_change_notes(a_yearly(), "mean_pct", "Average English share")
  })
  
  output$year_plot <- renderPlot({
    data <- a_yearly()
    ggplot2::ggplot(data, ggplot2::aes(x = year, y = mean_pct,
                                       colour = series, group = series)) +
      ggplot2::geom_line(linewidth = 1.05) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_x_continuous(breaks = sort(unique(data$year))) +
      ggplot2::scale_y_continuous(limits = c(0, 100)) +
      ggplot2::labs(title = "English Usage in K-pop Lyrics",
                    x = "Year", y = "Average English share (%)", colour = NULL) +
      lyrics_plot_theme() +
      ggplot2::theme(legend.position = if (nlevels(data$series) > 1L)
        "bottom" else "none")
  })
  
  # Tooltip for the yearly line chart: hover on desktop, tap on touch screens.
  output$year_plot_tooltip <- renderUI({
    event <- input$year_plot_hover
    if (is.null(event)) event <- input$year_plot_click
    if (is.null(event) || is.null(event$coords_css)) return(NULL)
    
    yearly <- a_yearly()
    point <- nearPoints(yearly, event, xvar = "year", yvar = "mean_pct",
                        threshold = 15, maxpoints = 1)
    if (nrow(point) == 0L) return(NULL)
    
    div(
      class = "d1-plot-tooltip",
      style = sprintf("left: %.0fpx; top: %.0fpx;",
                      event$coords_css$x, event$coords_css$y),
      strong(point$year),
      if (nlevels(yearly$series) > 1L) span(" · ", as.character(point$series)),
      br(),
      sprintf("Average English share: %.1f%%", point$mean_pct)
    )
  })
  
  hive_filtered <- reactive({
    data <- lyrics_required()
    req(input$year_range, input$hive_group)
    data <- data[data$year >= input$year_range[1] &
                   data$year <= input$year_range[2] &
                   !is.na(data$is_hive), , drop = FALSE]
    if (input$hive_group != "all") {
      data <- data[data$is_hive == input$hive_group, , drop = FALSE]
    }
    validate(need(nrow(data) > 0L, "No songs have a known affiliation for this group and year range."))
    data$group <- ifelse(data$is_hive == "true", "HYBE", "Non-HYBE")
    data
  })
  
  output$hive_stats <- renderUI({
    data <- hive_filtered()
    lyrics_stat_cards(summarize_lyrics(data.frame(
      series = data$group, pct = data$pct_main
    )))
  })
  
  output$hive_plot <- renderPlot({
    data <- hive_filtered() %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(mean_pct = mean(pct_main), .groups = "drop")
    ggplot2::ggplot(data, ggplot2::aes(x = group, y = mean_pct, fill = group)) +
      ggplot2::geom_col(width = 0.55) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f%%", mean_pct)),
        vjust = -0.5, size = 4
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 105)) +
      ggplot2::scale_fill_manual(values = c("HYBE" = "#18BC9C",
                                            "Non-HYBE" = "#3498DB")) +
      ggplot2::labs(title = "English Share by Company Affiliation",
                    x = NULL, y = "Mean English share (%)") +
      lyrics_plot_theme() + ggplot2::theme(legend.position = "none")
  })
  
  c_long <- reactive({
    long <- lyrics_make_long(lyrics_required(), "main", isTRUE(input$c_hive))
    validate(need(nrow(long) > 0L, "No songs match the current selection."))
    long
  })
  
  output$c_stats <- renderUI({
    req(input$c_thr)
    long <- c_long()
    lyrics_stat_cards(summarize_lyrics(long, input$c_thr),
                      threshold = input$c_thr)
  })
  
  c_yearly <- reactive({
    req(input$c_thr)
    c_long() %>%
      dplyr::group_by(year, series) %>%
      dplyr::summarise(share = mean(pct >= input$c_thr) * 100,
                       .groups = "drop")
  })
  
  output$c_observation <- renderUI({
    req(input$c_thr)
    lyrics_year_change_notes(
      c_yearly(), "share",
      sprintf("Songs with at least %d%% English", input$c_thr)
    )
  })
  
  output$c_plot <- renderPlot({
    req(input$c_thr)
    data <- c_yearly()
    data$year_label <- factor(data$year, levels = sort(unique(data$year)))
    ggplot2::ggplot(data, ggplot2::aes(x = year_label, y = share,
                                       fill = series)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                        width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f%%", share)),
        position = ggplot2::position_dodge(width = 0.8),
        vjust = -0.4, size = 3
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 105)) +
      ggplot2::labs(title = sprintf("Songs with at Least %d%% English", input$c_thr),
                    x = "Year", y = "Share of songs (%)", fill = NULL) +
      lyrics_plot_theme() +
      ggplot2::theme(legend.position = if (nlevels(data$series) > 1L)
        "bottom" else "none")
  })
  
  # AI DOCUMENTATION -----------------------------------------
  output$ai_table <- renderTable({
    data.frame(
      Tab = c(
        "Shoe Search",
        "Explain PCA",
        "Our Shoes",
        "Data 1"
      ),
      `What AI Did` = c(
        "Combined catalogue search, caching, image conversion, and result display.",
        "Combined PCA explanation plots, reconstruction, and leave-one-out evaluation.",
        "Implemented personal shoe loading, gallery display, and prediction summaries.",
        "Integrated K-pop Finding A/B charts, responsive filters, and descriptive result cards."
      ),
      `What the Team Did` = c(
        "'SHOE SEARCH' ",
        "'EXPLAIN PCA' ",
        "'OUR SHOES' ",
        "'DATA 1' "
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)
  
  output$ai_reflection <- renderTable({
    data.frame(
      Member = c(
        "박선우",
        "송현정",
        "한우영",
        "신유진",
        "Florian Urthaler"
      ),
      Reflection = c(
        "박선우",
        "송현정",
        "한우영",
        "신유진",
        "During this project, I mainly focused on combining our team's code, improving usability, and organizing the application. I worked on creating a consistent interface, improving navigation, cleaning up the folder structure, and testing the different functions.
        \n I used AI to help merge code, improve the R Shiny layout, remove unnecessary code, and fix errors. However, AI-generated solutions did not always work correctly and needed testing and adjustments. I learned that AI is a useful tool for speeding up development, but understanding and verifying the code remains important."
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)
}

# 11. START APPLICATION ---------------------------------------
shinyApp(ui = ui, server = server)
