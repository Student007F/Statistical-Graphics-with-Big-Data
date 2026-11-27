# ===========================================================
# CHIPMUNK - Project 3 (R Shiny)
# ===========================================================

library(shiny)
library(bslib)
library(png)
library(dplyr)
library(ggplot2)
library(plotly)

IMAGE_SIDE <- 64L
PIXEL_COUNT <- IMAGE_SIDE * IMAGE_SIDE
NEIGHBORS <- 5L
MAX_PCS <- 100L

DATA_DIR <- "data"
CATALOGUE_DIR <- file.path(DATA_DIR, "shoes")
CATALOGUE_CSV <- file.path(DATA_DIR, "shoes_meta.csv")
CACHE_DIR <- file.path("generated", "cache")
CACHE_FILE <- file.path(CACHE_DIR, "pca_cache.rds")
LYRICS_FILE <- file.path(DATA_DIR, "kpop_lyrics_counts.csv")

dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. SHOE CATALOGUE & PCA ------------------------------------
shoes_meta <- read.csv(CATALOGUE_CSV, stringsAsFactors = FALSE)
shoes_meta$id <- as.integer(shoes_meta$id)
shoes_meta$category <- trimws(as.character(shoes_meta$category))
catalogue_paths <- file.path(CATALOGUE_DIR, sprintf("%04d.png", shoes_meta$id))

category_levels <- sort(unique(shoes_meta$category))
category_colors <- setNames(
  grDevices::hcl.colors(length(category_levels), palette = "Dark 3"),
  category_levels
)

read_png_vector <- function(path) {
  img <- png::readPNG(path)
  if (length(dim(img)) == 3L) {
    img <- if (dim(img)[3] >= 3L) {
      0.2126 * img[, , 1] + 0.7152 * img[, , 2] + 0.0722 * img[, , 3]
    } else {
      img[, , 1]
    }
  }
  as.numeric(round(pmax(0, pmin(1, img)) * 255))
}

if (file.exists(CACHE_FILE)) {
  pca_data <- readRDS(CACHE_FILE)
} else {
  X <- t(vapply(catalogue_paths, read_png_vector, numeric(PIXEL_COUNT)))
  x_bar <- colMeans(X)
  centered <- sweep(X, 2, x_bar, "-")
  n_pcs <- min(MAX_PCS, nrow(X) - 1L, ncol(X))
  decomp <- svd(centered, nu = 0L, nv = n_pcs)
  loadings <- decomp$v[, seq_len(n_pcs), drop = FALSE]
  pca_data <- list(
    X = X,
    x_bar = x_bar,
    loadings = loadings,
    scores = centered %*% loadings,
    var_explained = 100 * decomp$d[seq_len(n_pcs)]^2 / sum(decomp$d^2)
  )
  saveRDS(pca_data, CACHE_FILE)
}

K_MAX <- ncol(pca_data$loadings)

DISC_PCS <- min(20L, K_MAX)

# ANOVA F-statistic per PC across all categories
pc_f_stats <- apply(pca_data$scores[, seq_len(DISC_PCS), drop = FALSE], 2, function(s) {
  summary(aov(s ~ shoes_meta$category))[[1]][["F value"]][1]
})

# Two-sample |t|-statistic per PC for each category vs. all other shoes
pc_t_stats <- sapply(category_levels, function(cat_name) {
  is_cat <- shoes_meta$category == cat_name
  apply(pca_data$scores[, seq_len(DISC_PCS), drop = FALSE], 2, function(s) {
    abs(t.test(s[is_cat], s[!is_cat])$statistic)
  })
})

LOO_FILE <- file.path(CACHE_DIR, "loo_accuracy.rds")

compute_loo_accuracy <- function(progress_cb = function(k) NULL) {
  if (file.exists(LOO_FILE)) return(readRDS(LOO_FILE))
  
  S <- pca_data$scores
  n <- nrow(S)
  truth <- as.integer(factor(shoes_meta$category, levels = category_levels))
  n_cat <- length(category_levels)
  D2 <- matrix(0, n, n)
  overall <- numeric(K_MAX)
  by_cat <- matrix(NA_real_, K_MAX, n_cat, dimnames = list(NULL, category_levels))
  
  for (k in seq_len(K_MAX)) {
    D2 <- D2 + outer(S[, k], S[, k], "-")^2
    diag(D2) <- Inf
    
    predicted <- vapply(seq_len(n), function(i) {
      d <- D2[i, ]
      nbrs <- order(d)[seq_len(NEIGHBORS)]
      classes <- truth[nbrs]
      counts <- tabulate(classes, nbins = n_cat)
      tied <- which(counts == max(counts))
      classes[which(classes %in% tied)[1]]
    }, integer(1))
    
    correct <- predicted == truth
    overall[k] <- mean(correct)
    by_cat[k, ] <- vapply(seq_len(n_cat), function(j) mean(correct[truth == j]), numeric(1))
    progress_cb(k)
  }
  
  res <- list(K = seq_len(K_MAX), overall = overall, by_cat = by_cat)
  saveRDS(res, LOO_FILE)
  res
}

predict_category <- function(categories) {
  counts <- table(categories)
  tied <- names(counts)[counts == max(counts)]
  categories[which(categories %in% tied)[1]]
}

find_similar_shoes <- function(query_vector, K, exclude_id = NULL) {
  idx <- seq_len(K)
  V <- pca_data$loadings[, idx, drop = FALSE]
  q_scores <- as.numeric((query_vector - pca_data$x_bar) %*% V)
  cat_scores <- pca_data$scores[, idx, drop = FALSE]
  d2 <- rowSums(sweep(cat_scores, 2, q_scores, "-")^2)
  
  if (!is.null(exclude_id)) {
    ex <- match(as.integer(exclude_id), shoes_meta$id)
    if (!is.na(ex)) d2[ex] <- Inf
  }
  
  nearest <- order(d2)[seq_len(NEIGHBORS)]
  top5 <- data.frame(
    Rank = seq_len(NEIGHBORS),
    ID = shoes_meta$id[nearest],
    Category = shoes_meta$category[nearest],
    Distance = round(sqrt(d2[nearest]), 2)
  )
  list(
    best_id = top5$ID[1],
    top5 = top5,
    prediction = predict_category(top5$Category),
    query_scores = q_scores
  )
}

reconstruct_shoe <- function(vector, K) {
  V <- pca_data$loadings[, seq_len(K), drop = FALSE]
  q <- as.numeric((vector - pca_data$x_bar) %*% V)
  as.numeric(pca_data$x_bar + V %*% q)
}

plot_shoe <- function(vector, title = "", signed = FALSE, bottom = "") {
  m <- matrix(as.numeric(vector), IMAGE_SIDE, IMAGE_SIDE)
  z <- t(m[IMAGE_SIDE:1L, , drop = FALSE])
  if (signed) {
    bound <- max(abs(vector))
    if (!is.finite(bound) || bound == 0) bound <- 1
    pal <- grDevices::colorRampPalette(c("#245A9C", "#FFFFFF", "#C64342"))(201)
    lims <- c(-bound, bound)
  } else {
    pal <- gray(seq(0, 1, length.out = 256))
    lims <- c(0, 255)
  }
  old <- par(mar = c(if (nzchar(bottom)) 2 else 0.5, 0.5, if (nzchar(title)) 2.5 else 0.5, 0.5))
  on.exit(par(old))
  image(
    seq_len(IMAGE_SIDE), seq_len(IMAGE_SIDE), z,
    col = pal, zlim = lims, axes = FALSE,
    xlab = "", ylab = "", asp = 1, main = title, useRaster = TRUE
  )
  if (nzchar(bottom)) mtext(bottom, side = 1, line = 0.6, cex = 0.9)
}

# 2. K-POP LYRICS DATA (DATA 1) ------------------------------
LYRICS_RULE_COLS <- c(main = "pct_main", sub_1 = "pct_sub_1", sub_2 = "pct_sub_2")
LYRICS_RULE_CHOICES <- c(
  "Main: all English words" = "main",
  "Sub 1: exclude English fillers (ad-libs)" = "sub_1",
  "Sub 2: repeated words count once per song" = "sub_2"
)

lyrics_data <- read.csv(LYRICS_FILE, stringsAsFactors = FALSE, check.names = FALSE)
lyrics_data$year <- suppressWarnings(as.integer(lyrics_data$year))
for (col in unname(LYRICS_RULE_COLS)) {
  vals <- suppressWarnings(as.numeric(gsub("%", "", as.character(lyrics_data[[col]]), fixed = TRUE)))
  vals[!is.finite(vals) | vals < 0 | vals > 100] <- NA_real_
  lyrics_data[[col]] <- vals
}
flags <- tolower(trimws(as.character(lyrics_data$is_hive)))
lyrics_data$is_hive <- ifelse(
  flags %in% c("true", "1", "yes"), "true",
  ifelse(flags %in% c("false", "0", "no"), "false", NA_character_)
)
lyrics_data <- lyrics_data[!is.na(lyrics_data$year) & !is.na(lyrics_data$pct_main), ]

lyrics_first_year <- min(lyrics_data$year)
lyrics_last_year <- max(lyrics_data$year)

lyrics_make_long <- function(data, rules, split_groups = FALSE) {
  if (split_groups) {
    data <- data[!is.na(data$is_hive), , drop = FALSE]
    data$group <- ifelse(data$is_hive == "true", "HYBE", "Non-HYBE")
  } else {
    data$group <- "All songs"
  }
  long <- dplyr::bind_rows(lapply(rules, function(r) {
    piece <- data[!is.na(data[[LYRICS_RULE_COLS[[r]]]]), , drop = FALSE]
    piece$rule <- r
    piece$pct <- piece[[LYRICS_RULE_COLS[[r]]]]
    piece
  }))
  long$series <- if (length(rules) > 1L && split_groups) {
    paste(long$rule, long$group, sep = " / ")
  } else if (length(rules) > 1L) {
    long$rule
  } else {
    long$group
  }
  long$series <- factor(long$series, levels = unique(long$series))
  long
}

lyrics_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# 3. USER INTERFACE ------------------------------------------
ui <- page_navbar(
  title = "Chipmunk",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  header = tags$head(tags$style(HTML("
    .shoe-caption { color: #576574; font-size: 0.93rem; }
    .pca-color-legend { display: grid; gap: 6px; margin: 4px 0 8px; font-size: 0.93rem; }
    .pca-swatch {
      display: inline-block; width: 14px; height: 14px; margin-right: 6px;
      border: 1px solid #b8c4cc; border-radius: 3px; vertical-align: -2px;
    }
    .d1-plot-wrap { position: relative; }
    .d1-plot-tooltip {
      position: absolute; z-index: 10; pointer-events: none;
      transform: translate(-50%, calc(-100% - 12px));
      background: rgba(44, 62, 80, 0.94); color: #fff;
      padding: 6px 10px; border-radius: 6px;
      font-size: 0.85rem; white-space: nowrap;
      box-shadow: 0 2px 6px rgba(0, 0, 0, 0.2);
    }
    .d1-method-details {
      margin-top: 14px; padding: 14px;
      border: 1px solid #dce6ed; border-radius: 10px;
    }
    .d1-method-details summary { cursor: pointer; font-weight: 650; }
    .d1-method-details p { margin-top: 12px; }
    .d1-method-grid {
      display: grid; grid-template-columns: repeat(auto-fit, minmax(155px, 1fr));
      gap: 12px; margin-top: 14px;
    }
    .d1-method-grid > div { padding: 12px; background: #f7fafc; border-radius: 8px; }
    .d1-example-table { margin-top: 12px; }
    .pca-flow {
      display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
      gap: 10px; margin: 10px 0;
    }
    .pca-flow > div {
      padding: 10px 12px; background: #f7fafc; border-left: 4px solid #18BC9C;
      border-radius: 8px; font-size: 0.92rem;
    }
    .pca-flow .step-tag { display: block; color: #18BC9C; font-weight: 700; font-size: 0.8rem; }
  "))),
  
  # TAB 1: PCA
  nav_panel(
    "PCA",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        title = "PCA Controls", width = 300,
        numericInput("cat_id", "Catalogue ID",
                     value = shoes_meta$id[1],
                     min = min(shoes_meta$id), max = max(shoes_meta$id), step = 1),
        sliderInput("K", "Number of Principal Components (K)",
                    min = 1, max = K_MAX, value = min(20L, K_MAX), step = 1),
        hr(),
        p(strong("Try a tied vote")),
        p(class = "shoe-caption",
          "At K = 20, these shoes get a 2 : 2 : 1 vote. Click one to load it (K is set to 20) and see how the tie is broken under the Top 5 table."),
        div(
          actionButton("tie_26", "ID 26", class = "btn-outline-primary btn-sm"),
          actionButton("tie_43", "ID 43", class = "btn-outline-primary btn-sm"),
          actionButton("tie_575", "ID 575", class = "btn-outline-primary btn-sm"),
          actionButton("tie_2", "ID 2", class = "btn-outline-primary btn-sm")
        ),
        p(class = "shoe-caption",
          "ID 575 and ID 2: the nearest shoe's category loses because it has only 1 vote (ID 2 ends up wrong).")
      ),
      navset_tab(
        id = "pca_tabs",
        nav_panel(
          "Shoe Search",
          card(
            card_header("Search Results and Prediction"),
            uiOutput("prediction"),
            hr(),
            layout_columns(
              col_widths = c(5, 7),
              card(card_header("Best Match"), plotOutput("best_match", height = "260px")),
              card(card_header("Top 5 Similar Shoes"), tableOutput("top5"), uiOutput("vote_breakdown"))
            ),
            tags$details(
              class = "d1-method-details",
              tags$summary("How are the distance and the predicted category decided?"),
              p(strong("1. Distance: "),
                "Every shoe is turned into K numbers (its PC scores). The distance between the selected shoe (q) and a catalogue shoe (s) is the straight-line (Euclidean) distance between their K scores:"),
              p(HTML("Distance = &radic;[ (q<sub>1</sub> &minus; s<sub>1</sub>)&sup2; + (q<sub>2</sub> &minus; s<sub>2</sub>)&sup2; + &hellip; + (q<sub>K</sub> &minus; s<sub>K</sub>)&sup2; ]")),
              p("A smaller distance means a more similar shoe. The selected shoe itself is left out, so it cannot match itself. Because only the first K scores are used, the distances and even the ranking can change when you move the K slider."),
              p(strong("2. Majority vote: "),
                "The 5 nearest shoes each vote for their own category, and the category with the most votes is the prediction."),
              p(strong("3. Tie-break: "),
                "If two categories get the same number of votes, the one whose shoe appears first in the ranking (i.e., the closest shoe) wins."),
              tags$table(
                class = "table table-sm table-striped d1-example-table",
                tags$thead(tags$tr(tags$th("Top 5 categories (rank 1 → 5)"), tags$th("Votes"), tags$th("Prediction"))),
                tags$tbody(
                  tags$tr(tags$td("Shoes, Shoes, Slippers, Boots, Boots"),
                          tags$td("Shoes 2 · Boots 2 · Slippers 1"),
                          tags$td("Shoes (tied with Boots, but Shoes is at rank 1)")),
                  tags$tr(tags$td("Slippers, Boots, Shoes, Shoes, Boots"),
                          tags$td("Boots 2 · Shoes 2 · Slippers 1"),
                          tags$td("Boots (first tied shoe is at rank 2; Slippers is closest but has only 1 vote)"))
                )
              ),
              p(class = "shoe-caption", "Illustrative examples. For real tied votes, click ID 26 or ID 43 under 'Try a tied vote' in the sidebar.")
            ),
            hr(),
            h4("PCA Image Comparison"),
            layout_columns(
              col_widths = c(4, 4, 4),
              card(card_header("Selected Shoe"), plotOutput("original_shoe", height = "220px")),
              card(card_header("Average Shoe"), plotOutput("average_shoe", height = "220px")),
              card(card_header("PCA Reconstruction"), plotOutput("reconstructed_shoe", height = "220px"))
            ),
                        tags$details(
              class = "d1-method-details",
              tags$summary("What do these three images show?"),
              p(strong("How the app sees a shoe: "),
                "Selected − Average = how this shoe differs from a typical shoe → that difference is compressed into K numbers (PC scores) → Average + those K numbers = Reconstruction."),
              tags$ol(
                tags$li(strong("Selected − Average: "),
                        "Subtracting the Average Shoe from the Selected Shoe leaves only what makes this shoe different from a typical shoe, for example a higher ankle or a darker sole."),
                tags$li(strong("Compressed into K numbers: "),
                        "From the whole catalogue, PCA has learned the most common patterns of difference between shoes (the principal components). Each PC score says how strongly one of these patterns appears in this shoe, so K numbers replace 4,096 pixels."),
                tags$li(strong("Average + those K numbers: "),
                        "Starting from the Average Shoe and adding each pattern in the amount given by its score rebuilds the shoe. The closer the Reconstruction looks to the original, the better the K numbers describe the shoe.")
              ),
              div(
                class = "d1-method-grid",
                div(strong("Selected Shoe"),
                    p("The chosen catalogue photo, stored as a 64 × 64 grayscale image: 4,096 numbers from 0 (black) to 255 (white). This shoe itself is excluded from the search results.")),
                div(strong("Average Shoe"),
                    p("The pixel-by-pixel average of all catalogue photos. Sharp areas are where most shoes look alike; blurry areas are where they differ. Every shoe is described by how it differs from this image.")),
                div(strong("PCA Reconstruction"),
                    p(HTML("The selected shoe rebuilt from only K numbers: Average + score<sub>1</sub> × PC<sub>1</sub> + … + score<sub>K</sub> × PC<sub>K</sub>. Fine details are lost when K is small.")),
                    p(strong(textOutput("recon_caption", inline = TRUE))),
                    p("Shoe Search measures distances using these same K numbers, so details missing here are also ignored by the search."))
              )
            ),
            h4("Top 5 Nearest Neighbors"),
            plotOutput("top5_gallery", height = "260px"),
            p(class = "shoe-caption",
              "Distance = Euclidean distance to the selected shoe using its first K PC scores (smaller = more similar).")
          )
        ),
        nav_panel(
          "Explain PCA",
          card(
            card_header("The big picture: how Shoe Search works"),
            p(strong("Key idea: "),
              "Instead of comparing 4,096 pixels, the app compares a few numbers (K PC scores) that capture the main ways shoes differ from each other. The steps below follow the photo selected in Shoe Search through this process."),
            div(
              class = "pca-flow",
              div(span(class = "step-tag", "START"),
                  "A photo is a 64 × 64 grayscale image = 4,096 numbers (0 = black, 255 = white)."),
              div(span(class = "step-tag", "STEPS 1–2"),
                  "Subtract the average shoe, leaving only how this shoe differs from a typical one."),
              div(span(class = "step-tag", "STEP 4"),
                  "PCA learns the most common patterns of difference across the catalogue (eigenshoes)."),
              div(span(class = "step-tag", "STEPS 5–6"),
                  "Describe the shoe with K scores, one per pattern; rebuilding the photo from them shows what is kept."),
              div(span(class = "step-tag", "SEARCH · STEP 3"),
                  "Compare K scores with every catalogue shoe; the 5 nearest vote on the category."),
              div(span(class = "step-tag", "STEP 7"),
                  "Choose K: enough PCs to tell categories apart, without keeping unnecessary detail.")
            )
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header("Step 1 | Average shoe"),
              p("The mean image of all catalogue photos (x-bar)."),
              plotOutput("pca_step1", height = "250px")
            ),
            card(
              card_header("Step 2 | Photo minus average"),
              p("Blue indicates darker than the average shoe (lower pixel value); red indicates brighter than the average shoe (higher pixel value)."),
              plotOutput("pca_step2", height = "250px")
            )
          ),
          card(
            card_header("Step 3 | 3D PCA score map & discriminative PCs"),
            selectInput(
              "score_map_mode",
              "Select 3D Score Map Axes",
              choices = c(
                "Standard (PC1, PC2, PC3)" = "standard",
                "Best overall separation (Top 3 ANOVA F-stat)" = "all",
                setNames(category_levels, paste("Best for", category_levels, "(Top 3 |t|-stat)"))
              ),
              selected = "standard",
              width = "420px"
            ),
            p(class = "shoe-caption", textOutput("pc_disc_text")),
            tags$details(
              class = "d1-method-details",
              tags$summary("How to read this 3D map"),
              tags$ul(
                tags$li(strong("Each dot is one catalogue shoe, "),
                        "placed by its scores on the 3 chosen PCs and colored by category. Dots that are close together look similar on these 3 patterns."),
                tags$li(strong("Black diamond: "), "the selected shoe. ",
                        strong("Red open circles: "), "its Top 5 matches from Shoe Search."),
                tags$li(strong("Why the Top 5 may not look closest here: "),
                        "the search measures distance with all K PCs, but this map shows only 3 of them. Two shoes can overlap on these 3 axes and still differ on the others, or the reverse."),
                tags$li(
                  strong("Axis options:"),
                  tags$ul(
                    tags$li(strong("Standard: "),
                            "PC1–PC3, the PCs that capture the most overall pixel variation. Variation is not the same as category separation (for example, PC1 may mostly reflect overall brightness)."),
                    tags$li(strong("Best overall separation: "),
                            "the 3 PCs whose category averages differ most compared with the spread within each category (ANOVA F-statistic)."),
                    tags$li(strong("Best for [category]: "),
                            "the 3 PCs that best separate one category from all other shoes (|t|-statistic).")
                  )
                ),
                tags$li(strong("Interaction: "),
                        "drag to rotate, scroll to zoom, hover over a dot to see its ID. Click a category in the legend to hide it; double-click to show only that category.")
              )
            ),
            plotlyOutput("pca_step3", height = "500px"),
            hr(),
            h5("Category Separation Strength (PC1–PC20)"),
            p(class = "shoe-caption",
              "Highlighted bars are the 3 PCs currently plotted in the 3D map above."),
            plotOutput("pc_disc_plot", height = "320px")
          ),
          card(
            card_header("Step 4 | First 10 eigenshoes"),
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
            plotOutput("pca_step4", height = "440px")
          ),
          card(
            card_header("Step 5 | Query scores"),
            p("The first K PCA coordinates used to calculate nearest-neighbor distances. These K numbers are all that Shoe Search knows about the selected shoe."),
            tags$details(
              class = "d1-method-details",
              tags$summary("How to read a score"),
              tags$ul(
                tags$li(strong("Score: "),
                        "how strongly the eigenshoe pattern of that PC (Step 4) appears in this shoe. It measures how well the difference image in Step 2 lines up with that pattern."),
                tags$li(strong("Sign: "),
                        "read it with the colors in Step 4 (the Meaning column applies that rule to each score). A score near 0 means the pattern barely matters for this shoe."),
                tags$li(strong("Catalogue percentile: "),
                        "raw scores are not comparable across PCs (PC1 scores spread much wider than later PCs), so each score is compared with the catalogue. 95% means this shoe's score is higher than or equal to that of 95% of catalogue shoes on that PC.")
              )
            ),
            tableOutput("pca_step5_table")
          ),
          card(
            card_header("Step 6 | Reconstruct using different K values"),
            p("Compare the same photo reconstructed using different numbers of components."),
            plotOutput("pca_step6", height = "300px")
          ),
          card(
            card_header("Step 7 | Choosing K: picture quality vs. classification accuracy"),
            p("There are two ways to judge K. ",
              strong("Picture quality"), " (left): how much of the pixel variation the first K PCs keep. ",
              strong("Classification accuracy"), " (right): how often the 5 nearest neighbors guess a shoe's category correctly. ",
              "A larger K always keeps more variation, but later PCs mostly add fine details, so accuracy can stop improving or even drop."),
            p("Each catalogue shoe is classified by its 5 nearest neighbors across K = 1 to 100."),
            actionButton("run_acc", "Compute / Load Accuracy", class = "btn-primary"),
            layout_columns(
              col_widths = c(6, 6),
              plotOutput("var_plot", height = "400px"),
              plotOutput("acc_plot", height = "400px")
            ),
            textOutput("acc_text"),
            uiOutput("k_compare")
          )
        )
      )
    )
  ),
  
  # TAB 2: DATA 1
  nav_panel(
    "Data 1",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        title = "Data 1 Controls", width = 300,
        radioButtons("d1_view", "Select Analysis",
                     choices = c("Analysis 1", "Analysis 2"),
                     selected = "Analysis 1"),
        hr(),
        conditionalPanel(
          "input.d1_view == 'Analysis 1'",
          sliderInput("year_range", "Select Year Range",
                      min = lyrics_first_year, max = lyrics_last_year,
                      value = c(lyrics_first_year, lyrics_last_year),
                      step = 1, sep = ""),
          checkboxGroupInput("a_rules", "Counting Rules (one or more)",
                             choices = LYRICS_RULE_CHOICES, selected = "main"),
          checkboxInput("a_hive", "Split by HYBE / Non-HYBE", FALSE),
          checkboxInput("a_show_b", "Show detailed HYBE comparison", FALSE),
          conditionalPanel(
            "input.a_show_b",
            selectInput("hive_group", "Group for detailed comparison",
                        choices = c("Both groups" = "all", "HYBE" = "true",
                                    "Non-HYBE" = "false"),
                        selected = "all")
          ),
          p(class = "shoe-caption",
            "A song with an unknown HYBE affiliation is excluded when comparing groups.")
        ),
        conditionalPanel(
          "input.d1_view == 'Analysis 2'",
          sliderInput("c_thr", "Minimum English Share",
                      min = 10, max = 90, value = 50, step = 5, post = "%"),
          checkboxInput("c_hive", "Split by HYBE / Non-HYBE", FALSE),
          p(class = "shoe-caption",
            "Analysis 2 always uses the main counting rule and all available years.")
        )
      ),
      conditionalPanel(
        "input.d1_view == 'Analysis 1'",
        card(
          card_header("Analysis 1 | English in K-pop lyrics over time"),
          p("How has the proportion of English words in charted K-pop lyrics changed across years?"),
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
        card(card_header("Summary Statistics"), tableOutput("a_stats")),
        card(
          card_header("Average English Share by Year"),
          div(
            class = "d1-plot-wrap",
            plotOutput("year_plot", height = "400px",
                       click = "year_plot_click",
                       hover = hoverOpts("year_plot_hover", delay = 80,
                                         delayType = "debounce")),
            uiOutput("year_plot_tooltip")
          ),
          p(class = "shoe-caption",
            "Each point shows the unweighted mean of song-level English percentages for that year and selected group/rule. Hover over or tap a point to see its value.")
        ),
        card(
          card_header("Why It Matters and Next Analysis"),
          p(
            strong("Why it matters: "),
            "In these charted songs, the average English share rose from 12.2% in 2014 to 35.6% in 2025 (Main rule, all songs). The same rise appears under Sub 1 and Sub 2, so it is not just an effect of ad-libs or repeated words. This dataset describes lyric composition, not whether English lyrics cause greater popularity."
          ),
          p(strong("How This Helps Our Next Analysis: ")),
          p(
            strong("1. Pin Down When the Change Happens"), br(),
            "The yearly average stays between 11% and 20% through 2022, then rises to 25.3% (2023), 34.0% (2024) and 35.6% (2025). We will use 2022-2025 as the window in which to compare other data, such as search interest or tourist arrivals."
          ),
          p(
            strong("2. Test the HYBE-First Pattern"), br(),
            "In 2021, 4 of 7 HYBE songs (57.1%) were at least 50% English versus 8.7% of other songs; by 2025 the figures were 60.0% and 39.3%. HYBE has only 6-15 songs per year since 2018, so we will repeat this with each song counted once."
          )
        ),
        conditionalPanel(
          "input.a_show_b",
          card(
            card_header("HYBE vs Non-HYBE | Detailed Comparison"),
            p("Compare average English usage across the two company-affiliation groups within the selected year range."),
            tableOutput("hive_stats"),
            plotOutput("hive_plot", height = "360px"),
            p(class = "shoe-caption",
              "The selected year range and group in the control panel determine which songs appear here (using the main counting rule). These are descriptive differences, not evidence of a causal effect.")
          )
        )
      ),
      conditionalPanel(
        "input.d1_view == 'Analysis 2'",
        card(
          card_header("Analysis 2 | How many songs are English-centered?"),
          p("What proportion of charted K-pop songs reach or exceed a given minimum English lyric threshold across years?")
        ),
        card(card_header("Summary Statistics"), tableOutput("c_stats")),
        card(
          card_header("English-Centered Songs by Year"),
          plotOutput("c_plot", height = "400px"),
          p(class = "shoe-caption",
            "For each year, the percentage is the number of songs meeting the threshold divided by the number of songs with valid lyrics data in that year.")
        ),
        card(
          card_header("Why It Matters and Next Analysis"),
          p(
            strong("Why it matters: "),
            "At the 50% threshold, the share of English-centered songs rose from 12.0% in 2022 to 41.4% in 2025 (Main rule, all songs), so the rise in the average comes with a larger group of songs that are mostly English. Average English share alone does not reveal this, because it does not show how many sampled songs exceed a chosen threshold. The threshold-based percentage describes that proportion within the dataset, not all K-pop releases."
          ),
          p(
            strong("How this helps our next analysis: "),
            "The share of songs with at least 50% English rose from 12.0% in 2022 to 41.4% in 2025. In a later project, we plan to bring in other data, such as search interest or tourist arrivals, and check whether they change during the same 2022-2025 period. We would also add artist information, such as group or solo, to describe which artists make these English-centered songs."
          )
        )
      )
    )
  ),
  
  # TAB 3: DATA 2
  nav_panel(
    "Data 2",
    card(
      card_header("Data 2 | Exploratory Data Analysis"),
      p("Add EDA for a new dataset here (at least one of Data 1 or Data 2 must be high-dimensional: p > n).")
    )
  ),
  
  # TAB 4: AI
  nav_panel(
    "AI",
    card(
      card_header("Ethics Statement"),
      p("We understand the rules on data and AI use for this project, including following the license and privacy rules of every dataset we used, understanding all of our code, and knowing that presenting AI output as our own work may be academic misconduct. We certify that we have followed all of them. Any violation is our responsibility, and we agree to receive an F in this course if we have not."),
      p(strong("Signed: 박선우, 송현정, 한우영, 신유진, Florian Urthaler"))
    ),
    card(card_header("What AI Did and What We Did"), tableOutput("ai_table")),
    card(card_header("Reflections"), tableOutput("ai_reflection")),
    card(card_header("Contributions"), tableOutput("ai_contribution"))
  )
)

# 4. SERVER --------------------------------------------------
server <- function(input, output, session) {
  
  # PCA TAB
  query_data <- reactive({
    req(input$cat_id)
    idx <- match(as.integer(input$cat_id), shoes_meta$id)
    validate(need(!is.na(idx), "Catalogue ID not found."))
    list(
      vector = as.numeric(pca_data$X[idx, ]),
      actual = shoes_meta$category[idx],
      id = shoes_meta$id[idx]
    )
  })
  
  shoe_search <- reactive({
    q <- query_data()
    find_similar_shoes(q$vector, input$K, exclude_id = q$id)
  })
  
  output$prediction <- renderUI({
    res <- shoe_search()
    actual <- query_data()$actual
    verdict <- if (tolower(res$prediction) == tolower(actual)) "Correct" else "Incorrect"
    tagList(
      p(strong("Predicted Category: "), res$prediction),
      p(strong("Actual Category: "), actual),
      p(strong("Result: "), verdict)
    )
  })
  
  output$top5 <- renderTable(shoe_search()$top5, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$vote_breakdown <- renderUI({
    res <- shoe_search()
    cats <- res$top5$Category
    counts <- table(factor(cats, levels = unique(cats)))   # rank order of first appearance
    max_n <- max(counts)
    tied <- names(counts)[counts == max_n]
    vote_text <- paste(sprintf("%s %d", names(counts), as.integer(counts)), collapse = " · ")
    
    explain <- if (length(tied) > 1L) {
      sprintf("Tie between %s (%d votes each). %s is chosen because its closest shoe is at rank #%d, ahead of the other tied categories.",
              paste(tied, collapse = " and "), max_n, res$prediction, match(res$prediction, cats))
    } else {
      sprintf("%s has the most votes (%d of %d).", res$prediction, max_n, NEIGHBORS)
    }
    
    tagList(
      p(strong("Votes: "), vote_text),
      p(class = "shoe-caption", explain)
    )
  })
  
  load_tie_example <- function(id) {
    updateSliderInput(session, "K", value = 20)
    updateNumericInput(session, "cat_id", value = id)
    nav_select("pca_tabs", "Shoe Search")
  }
  observeEvent(input$tie_26, load_tie_example(26))
  observeEvent(input$tie_43, load_tie_example(43))
  observeEvent(input$tie_575, load_tie_example(575))
  observeEvent(input$tie_2, load_tie_example(2))
  
  output$best_match <- renderPlot({
    id <- shoe_search()$best_id
    row <- match(id, shoes_meta$id)
    plot_shoe(pca_data$X[row, ], paste("Catalogue ID", id))
  })
  
  output$original_shoe <- renderPlot({
    plot_shoe(query_data()$vector)
  })
  
  output$average_shoe <- renderPlot({
    plot_shoe(pca_data$x_bar)
  })
  
  output$reconstructed_shoe <- renderPlot({
    plot_shoe(reconstruct_shoe(query_data()$vector, input$K))
  })
    output$recon_caption <- renderText({
    k <- input$K
    sprintf("K = %d keeps %.1f%% of the pixel variation across catalogue photos.",
            k, cumsum(pca_data$var_explained)[k])
  })
  
  output$top5_gallery <- renderPlot({
    top5 <- shoe_search()$top5
    old <- par(mfrow = c(1, NEIGHBORS), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (j in seq_len(NEIGHBORS)) {
      row <- match(top5$ID[j], shoes_meta$id)
      plot_shoe(pca_data$X[row, ], paste0("#", j, " ID ", top5$ID[j], "\n", top5$Category[j]),
                bottom = sprintf("Distance = %.2f", top5$Distance[j]))
    }
  })
  
  output$pca_step1 <- renderPlot({
    plot_shoe(pca_data$x_bar, "Average shoe")
  })
  
  output$pca_step2 <- renderPlot({
    plot_shoe(query_data()$vector - pca_data$x_bar, "Selected photo - average", signed = TRUE)
  })
  
  selected_pcs <- reactive({
    mode <- input$score_map_mode
    if (is.null(mode) || mode == "standard") return(1:3)
    if (mode == "all") return(order(pc_f_stats, decreasing = TRUE)[1:3])
    order(pc_t_stats[, mode], decreasing = TRUE)[1:3]
  })
  
  output$pc_disc_text <- renderText({
    pcs <- selected_pcs()
    mode <- input$score_map_mode
    basis <- if (mode == "standard") {
      "highest overall pixel variance"
    } else if (mode == "all") {
      "highest ANOVA F-statistic across all categories"
    } else {
      paste0("highest |t|-statistic separating ", mode, " from other shoes")
    }
    sprintf("Plotting PC%d, PC%d, and PC%d (%s).", pcs[1], pcs[2], pcs[3], basis)
  })
  
  output$pca_step3 <- renderPlotly({
    pcs <- selected_pcs()
    mode <- input$score_map_mode
    q <- query_data()$vector
    pos <- as.numeric((q - pca_data$x_bar) %*% pca_data$loadings[, pcs, drop = FALSE])
    coords <- pca_data$scores[, pcs, drop = FALSE]
    axis_labels <- sprintf("PC%d (%.1f%%)", pcs, pca_data$var_explained[pcs])
    
    p <- plot_ly()
    
    # Draw each category separately so we can dim non-target categories
    for (cat in category_levels) {
      idx <- which(shoes_meta$category == cat)
      is_focus <- mode %in% c("standard", "all") || mode == cat
      p <- p %>%
        add_trace(
          x = coords[idx, 1], y = coords[idx, 2], z = coords[idx, 3],
          type = "scatter3d", mode = "markers", name = cat,
          marker = list(
            size = if (is_focus) 3.5 else 2.2,
            opacity = if (is_focus) 0.75 else 0.45,
            color = unname(category_colors[cat])
          ),
          text = paste("ID:", shoes_meta$id[idx], "<br>Category:", cat),
          hoverinfo = "text"
        )
    }
    
    # Highlight Top 5 nearest neighbors from Shoe Search
    top5_rows <- match(shoe_search()$top5$ID, shoes_meta$id)
    p <- p %>%
      add_trace(
        x = coords[top5_rows, 1], y = coords[top5_rows, 2], z = coords[top5_rows, 3],
        type = "scatter3d", mode = "markers", name = "Top 5 matches",
        marker = list(size = 6, color = "#E74C3C", symbol = "circle-open",
                      line = list(width = 3)),
        text = paste0("Match #", seq_len(NEIGHBORS), " | ID: ", shoe_search()$top5$ID,
                      "<br>Category: ", shoe_search()$top5$Category),
        hoverinfo = "text"
      ) %>%
      # Highlight the selected query shoe
      add_trace(
        x = pos[1], y = pos[2], z = pos[3],
        type = "scatter3d", mode = "markers", name = "Selected photo",
        marker = list(size = 8, color = "black", symbol = "diamond",
                      line = list(color = "#FFD43B", width = 2)),
        text = paste("Selected ID:", query_data()$id),
        hoverinfo = "text"
      ) %>%
      layout(
        margin = list(l = 0, r = 0, t = 10, b = 0),
        legend = list(itemsizing = "constant"),
        scene = list(
          aspectmode = "cube",
          camera = list(eye = list(x = 1.45, y = 1.45, z = 1.15)),
          xaxis = list(title = axis_labels[1]),
          yaxis = list(title = axis_labels[2]),
          zaxis = list(title = axis_labels[3])
        )
      ) %>%
      config(displaylogo = FALSE)
    
    p
  })
  
  output$pc_disc_plot <- renderPlot({
    pcs <- selected_pcs()
    mode <- input$score_map_mode
    
    if (mode %in% c("standard", "all")) {
      vals <- pc_f_stats
      y_lab <- "ANOVA F-statistic (All Categories)"
    } else {
      vals <- pc_t_stats[, mode]
      y_lab <- paste0("|t|-statistic (", mode, " vs. Others)")
    }
    
    df <- data.frame(
      PC = factor(paste0("PC", seq_len(DISC_PCS)), levels = paste0("PC", seq_len(DISC_PCS))),
      Stat = as.numeric(vals),
      Active = seq_len(DISC_PCS) %in% pcs
    )
    
    ggplot2::ggplot(df, ggplot2::aes(x = PC, y = Stat, fill = Active)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#18BC9C", "FALSE" = "#BDC3C7"),
        labels = c("TRUE" = "Plotted in 3D map", "FALSE" = "Other PCs"),
        name = NULL
      ) +
      ggplot2::labs(x = NULL, y = y_lab) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom"
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
    s <- shoe_search()$query_scores
    pct <- vapply(seq_along(s), function(k) 100 * mean(pca_data$scores[, k] <= s[k]), numeric(1))
    data.frame(
      PC = paste0("PC", seq_along(s)),
      Score = round(s, 2),
      `Catalogue percentile` = sprintf("%.0f%%", pct),
      Meaning = ifelse(
        seq_along(s) > 10L, "(eigenshoe not shown in Step 4)",
        ifelse(s >= 0, "Red areas brighter, blue areas darker",
               "Red areas darker, blue areas brighter")
      ),
      check.names = FALSE
    )
  }, striped = TRUE, spacing = "s")
  
  output$pca_step6 <- renderPlot({
    ks <- sort(unique(c(1L, 5L, 20L, min(50L, K_MAX), input$K)))
    ks <- ks[ks <= K_MAX]
    old <- par(mfrow = c(1, length(ks)), mar = c(1, 1, 3, 1))
    on.exit(par(old))
    for (k in ks) {
      plot_shoe(
        reconstruct_shoe(query_data()$vector, k),
        paste0("K = ", k, if (k == input$K) "\n(selected)" else "")
      )
    }
  })
  
  accuracy_result <- eventReactive(input$run_acc, {
    if (file.exists(LOO_FILE)) return(readRDS(LOO_FILE))
    if (is.null(input$run_acc) || input$run_acc == 0) return(NULL)
    withProgress(message = "Computing leave-one-out accuracy", value = 0, {
      compute_loo_accuracy(function(k) incProgress(1 / K_MAX, detail = paste("K =", k)))
    })
  }, ignoreNULL = FALSE)
  
  output$acc_plot <- renderPlot({
    res <- accuracy_result()
    validate(need(!is.null(res), "Click the button to calculate accuracy across K."))
    best_k <- which.max(res$overall)
    plot(
      res$K, 100 * res$overall, type = "l", lwd = 3, ylim = c(0, 100),
      xlab = "Number of PCs (K)", ylab = "Accuracy (%)",
      main = "Catalogue Leave-One-Out Accuracy"
    )
    for (j in seq_along(category_levels)) {
      lines(res$K, 100 * res$by_cat[, j], col = category_colors[j], lwd = 1.5, lty = 2)
    }
    abline(v = input$K, col = "red", lty = 3)
    points(input$K, 100 * res$overall[input$K], col = "red", pch = 19)
    points(best_k, 100 * res$overall[best_k], pch = 8, cex = 2)
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
    res <- accuracy_result()
    req(res)
    best_k <- which.max(res$overall)
    sprintf(
      "Selected K = %d: %.1f%% accuracy. Best K = %d: %.1f%% accuracy.",
      input$K, 100 * res$overall[input$K], best_k, 100 * res$overall[best_k]
    )
  })
  
  output$var_plot <- renderPlot({
    cum <- cumsum(pca_data$var_explained)
    k <- input$K
    plot(
      seq_len(K_MAX), pca_data$var_explained, type = "h", col = "#BDC3C7", lwd = 2,
      ylim = c(0, 100), xlab = "Number of PCs (K)", ylab = "Variance explained (%)",
      main = "Pixel Variation Kept by K PCs"
    )
    lines(seq_len(K_MAX), cum, lwd = 3, col = "#2C3E50")
    abline(v = k, col = "red", lty = 3)
    points(k, cum[k], col = "red", pch = 19)
    legend(
      "right", bg = "white", cex = 0.8,
      legend = c("Cumulative", "Each PC", paste("Selected K =", k)),
      col = c("#2C3E50", "#BDC3C7", "red"),
      lty = c(1, 1, 3), lwd = c(3, 2, 1)
    )
  })
  
  output$k_compare <- renderUI({
    k <- input$K
    cum <- cumsum(pca_data$var_explained)
    base <- sprintf("K = %d keeps %.1f%% of the pixel variation.", k, cum[k])
    res <- accuracy_result()
    if (is.null(res)) {
      return(p(class = "shoe-caption", base,
               " Click the button above to compare this with classification accuracy."))
    }
    best_k <- which.max(res$overall)
    p(class = "shoe-caption", sprintf(
      "%s Accuracy is highest at K = %d (%.1f%%), where %.1f%% of the variation is kept. Beyond that point, extra PCs still improve the reconstruction but do not improve classification.",
      base, best_k, 100 * res$overall[best_k], cum[best_k]
    ))
  })
  
  # DATA 1 TAB
  a_long <- reactive({
    validate(need(length(input$a_rules) > 0L, "Select at least one counting rule."))
    req(input$year_range)
    df <- lyrics_data[lyrics_data$year >= input$year_range[1] &
                        lyrics_data$year <= input$year_range[2], , drop = FALSE]
    long <- lyrics_make_long(df, input$a_rules, isTRUE(input$a_hive))
    validate(need(nrow(long) > 0L, "No valid songs match the selection."))
    long
  })
  
  a_yearly <- reactive({
    a_long() %>%
      dplyr::group_by(year, series) %>%
      dplyr::summarise(mean_pct = mean(pct), .groups = "drop")
  })
  
  output$a_stats <- renderTable({
    a_long() %>%
      dplyr::group_by(Series = series) %>%
      dplyr::summarise(
        Songs = dplyr::n(),
        `Mean (%)` = round(mean(pct), 1),
        `Median (%)` = round(median(pct), 1),
        .groups = "drop"
      )
  }, striped = TRUE, bordered = TRUE)
  
  output$year_plot <- renderPlot({
    data <- a_yearly()
    ggplot2::ggplot(data, ggplot2::aes(x = year, y = mean_pct, colour = series, group = series)) +
      ggplot2::geom_line(linewidth = 1.05) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_x_continuous(breaks = sort(unique(data$year))) +
      ggplot2::scale_y_continuous(limits = c(0, 100)) +
      ggplot2::labs(title = "English Usage in K-pop Lyrics",
                    x = "Year", y = "Average English share (%)", colour = NULL) +
      lyrics_plot_theme() +
      ggplot2::theme(legend.position = if (nlevels(data$series) > 1L) "bottom" else "none")
  })
  
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
    req(input$year_range, input$hive_group)
    df <- lyrics_data[lyrics_data$year >= input$year_range[1] &
                        lyrics_data$year <= input$year_range[2] &
                        !is.na(lyrics_data$is_hive), , drop = FALSE]
    if (input$hive_group != "all") df <- df[df$is_hive == input$hive_group, , drop = FALSE]
    validate(need(nrow(df) > 0L, "No songs match this group and year range."))
    df$group <- ifelse(df$is_hive == "true", "HYBE", "Non-HYBE")
    df
  })
  
  output$hive_stats <- renderTable({
    hive_filtered() %>%
      dplyr::group_by(Group = group) %>%
      dplyr::summarise(
        Songs = dplyr::n(),
        `Mean (%)` = round(mean(pct_main), 1),
        `Median (%)` = round(median(pct_main), 1),
        .groups = "drop"
      )
  }, striped = TRUE, bordered = TRUE)
  
  output$hive_plot <- renderPlot({
    data <- hive_filtered() %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(mean_pct = mean(pct_main), .groups = "drop")
    ggplot2::ggplot(data, ggplot2::aes(x = group, y = mean_pct, fill = group)) +
      ggplot2::geom_col(width = 0.55) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", mean_pct)), vjust = -0.5, size = 4) +
      ggplot2::scale_y_continuous(limits = c(0, 105)) +
      ggplot2::scale_fill_manual(values = c("HYBE" = "#18BC9C", "Non-HYBE" = "#3498DB")) +
      ggplot2::labs(title = "English Share by Company Affiliation",
                    x = NULL, y = "Mean English share (%)") +
      lyrics_plot_theme() +
      ggplot2::theme(legend.position = "none")
  })
  
  c_long <- reactive({
    lyrics_make_long(lyrics_data, "main", isTRUE(input$c_hive))
  })
  
  output$c_stats <- renderTable({
    req(input$c_thr)
    c_long() %>%
      dplyr::group_by(Series = series) %>%
      dplyr::summarise(
        Songs = dplyr::n(),
        `Mean (%)` = round(mean(pct), 1),
        `Median (%)` = round(median(pct), 1),
        `Above Threshold (%)` = round(mean(pct >= input$c_thr) * 100, 1),
        .groups = "drop"
      )
  }, striped = TRUE, bordered = TRUE)
  
  output$c_plot <- renderPlot({
    req(input$c_thr)
    data <- c_long() %>%
      dplyr::group_by(year, series) %>%
      dplyr::summarise(share = mean(pct >= input$c_thr) * 100, .groups = "drop")
    data$year_label <- factor(data$year, levels = sort(unique(data$year)))
    ggplot2::ggplot(data, ggplot2::aes(x = year_label, y = share, fill = series)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", share)),
                         position = ggplot2::position_dodge(width = 0.8), vjust = -0.4, size = 3) +
      ggplot2::scale_y_continuous(limits = c(0, 105)) +
      ggplot2::labs(title = sprintf("Songs with at Least %d%% English", input$c_thr),
                    x = "Year", y = "Share of songs (%)", fill = NULL) +
      lyrics_plot_theme() +
      ggplot2::theme(legend.position = if (nlevels(data$series) > 1L) "bottom" else "none")
  })
  
  # AI TAB
  output$ai_table <- renderTable({
    data.frame(
      Tab = c("PCA", "Data 1", "Data 2", "AI"),
      `R Shiny (AI vs Team)` = rep("TODO: summarize AI vs team work for R Shiny", 4),
      `Streamlit (AI vs Team)` = rep("TODO: summarize AI vs team work for Streamlit", 4),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)
  
  output$ai_reflection <- renderTable({
    data.frame(
      Member = c("박선우", "송현정", "한우영", "신유진", "Florian Urthaler"),
      Reflection = rep("TODO: 1-2 paragraphs (what AI helped with, where it failed, what you learned)", 5)
    )
  }, striped = TRUE, bordered = TRUE)
  
  output$ai_contribution <- renderTable({
    data.frame(
      Member = c("박선우", "송현정", "한우영", "신유진", "Florian Urthaler"),
      Contribution = rep("TODO: exact contribution of each member", 5)
    )
  }, striped = TRUE, bordered = TRUE)
}

shinyApp(ui = ui, server = server)
