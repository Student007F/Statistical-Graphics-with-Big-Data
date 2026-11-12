# Project 2 template: the tabs and the inputs only. Everything else, and the look of the app, is yours.
library(shiny)
library(dplyr)
library(ggplot2)
library(bslib)

# Example: read one photo as 4,096 numbers (64 x 64 grayscale, row by row, 0 = black to 255 = white)
library(magick)

# 1. 메타데이터 읽기 (id와 category가 들어있는 파일)
meta <- read.csv("shoes_meta.csv")

# 2. 2,000장의 신발 이미지를 한 번에 읽어서 행렬(Matrix)로 만드는 함수
load_all_shoes <- function() {
  img_list <- lapply(1:2000, function(i) {
    path <- file.path("shoes", sprintf("%04d.png", i))
    img <- image_read(path)
    # 64x64 흑백 이미지로 변환 후 숫자 벡터(0~255)로 바꿈
    as.numeric(t(matrix(as.integer(image_data(img, channels = "gray")), 64, 64)))
  })
  # ★ 핵심: 리스트들을 2000행 x 4096열의 행렬(Matrix)로 묶어서 반환합니다!
  shoe_mat <- do.call(rbind, img_list)
  return(shoe_mat)
}

#앱이 켜질 때 딱 한 번만 실행되도록 여기에 두는 것입니다.
cat("Loading shoe images... please wait.\n")
shoe_matrix <- load_all_shoes()
cat("Images loaded.\n")

# 3. 평균 신발(x_bar) 계산 및 PCA(주성분 분석) 실행
x_bar <- colMeans(shoe_matrix)
pca_result <- prcomp(shoe_matrix, center = TRUE, scale. = FALSE)

# [수정 1] 어떤 사진(JPG/PNG, 가로·세로 무엇이든)이든 카탈로그와 같은 형식의 벡터로 바꾸는 함수
#   - 휴대폰 사진의 회전 정보 적용 → 투명 배경은 흰색으로 → 흰색을 채워 정사각형으로
#     (신발이 안 찌그러지고 안 잘림) → 64x64 → 흑백 → 0~255 숫자 4096개
load_any <- function(path) {
  img <- image_read(path)[1]                              # 여러 프레임이면 첫 장만
  img <- image_orient(img)                                # 휴대폰 사진의 회전(EXIF) 적용
  img <- image_background(img, "white", flatten = TRUE)   # 투명 배경 -> 흰색
  info <- image_info(img)
  s <- max(info$width, info$height)                       # 긴 변 기준 정사각형
  img <- image_extent(img, paste0(s, "x", s), gravity = "center", color = "white")  # 남는 곳은 흰색으로 채움 (자르지 않음)
  img <- image_resize(img, "64x64!")                      # 이미 정사각형이라 비율이 안 깨짐
  img <- image_convert(img, colorspace = "gray")          # 흑백으로
  as.numeric(t(matrix(as.integer(image_data(img, channels = "gray")), 64, 64)))
}

# [수정 5] 예측 동점 규칙: 최다 득표 카테고리가 여러 개면, 그중 가장 가까운 사진의 카테고리를 선택
#   top5는 거리순으로 정렬되어 있으므로, 위에서부터 처음 만나는 동점 카테고리가 가장 가까운 사진의 것
predict_cat <- function(top5) {
  tb   <- table(top5$Category)
  tied <- names(tb)[tb == max(tb)]                # 최다 득표 카테고리들 (동점이면 2개)
  top5$Category[top5$Category %in% tied][1]
}

# [수정 9] K별 예측 정확도 (leave-one-out): 카탈로그 2,000장 각각을 "자기 자신만 빼고" 검색해서
#   top 5 다수결(동점이면 가장 가까운 사진)로 카테고리를 예측하고, 실제 카테고리와 비교
#   -> 정답(카테고리)을 아는 카탈로그 사진에만 가능. 한 번 계산하면 acc_cache에 저장해서 다시 계산하지 않음.
acc_cache <- new.env()

compute_loo_acc <- function(progress = function(k) NULL, max_k = 100) {
  if (!is.null(acc_cache$res)) return(acc_cache$res)
  S     <- pca_result$x
  n     <- nrow(S)
  cat_f <- factor(meta$category)
  truth <- as.integer(cat_f)
  n_cat <- nlevels(cat_f)
  D2      <- matrix(0, n, n)                       # 제곱 거리 (PC를 하나씩 더하며 누적)
  overall <- numeric(max_k)
  by_cat  <- matrix(NA_real_, max_k, n_cat, dimnames = list(NULL, levels(cat_f)))
  
  for (k in 1:max_k) {
    D2 <- D2 + outer(S[, k], S[, k], "-")^2        # PC k까지 쓴 제곱 거리
    Dk <- D2
    diag(Dk) <- Inf                                # 자기 자신 제외
    # 각 사진마다 가장 가까운 5장의 위치 (거리순)
    idx <- t(apply(Dk, 1, function(d) {
      thr  <- sort.int(d, partial = 5)[5]
      cand <- which(d <= thr)
      cand[order(d[cand])][1:5]
    }))
    N <- matrix(truth[idx], nrow = n)              # 이웃 5장의 카테고리 번호
    pred <- vapply(seq_len(n), function(i) {       # 다수결, 동점이면 가장 가까운 이웃의 카테고리
      nb   <- N[i, ]
      cnt  <- tabulate(nb, nbins = n_cat)
      tied <- which(cnt == max(cnt))
      as.integer(nb[nb %in% tied][1])
    }, integer(1))
    ok <- (pred == truth)
    overall[k] <- mean(오케이)
    by_cat[k, ] <- tapply(ok, factor(truth, levels = 1:n_cat), mean)
    progress(k)
  }
  acc_cache$res <- list(K = 1:max_k, overall = overall, by_cat = by_cat)
  acc_cache$res
}

# [수정 11] Our Shoes ------------------------------------------------------------
# 임의의 사진 벡터(길이 4096)에 대해 가장 가까운 카탈로그 5장을 찾는 함수 (Shoe Search와 같은 계산)
#   (이름을 search로 하면 R 기본 함수 search()를 가리므로 find_top5로 함)
find_top5 <- function(vec, K) {
  q  <- as.vector((vec - x_bar) %*% pca_result$rotation[, 1:K])
  cs <- pca_result$x[, 1:K, drop = FALSE]
  d  <- sqrt(rowSums((cs - matrix(q, nrow(cs), K, byrow = TRUE))^2))
  top <- order(d)[1:5]
  data.frame(Rank = 1:5, ID = top, Category = meta$category[top], Distance = round(d[top], 2))
}

# 0~255 벡터를 64x64 흑백 그림으로 그리는 함수 (정사각형 비율, 밝기 범위 고정)
draw_gray <- function(vec, main = "", sub = "") {
  m <- matrix(vec, 64, 64, byrow = TRUE)
  image(t(apply(m, 2, rev)), col = gray(0:100/100), zlim = c(0, 255), axes = FALSE, asp = 1,
        main = main, sub = sub)
}

# ourshoes/ 폴더의 사진 읽기. 파일 이름 규칙: 이름_카테고리_번호.jpg (예: Whae_Sandals_1.jpg)
cat_levels <- sort(unique(meta$category))
our_files  <- list.files("ourshoes", pattern = "\\.(jpg|jpeg|png)$", ignore.case = TRUE)
if (length(our_files) > 0) {
  parts <- strsplit(tools::file_path_sans_ext(our_files), "_")
  our <- data.frame(
    file     = our_files,
    name     = vapply(parts, function(p) p[1], character(1)),
    category = vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, character(1)),
    stringsAsFactors = FALSE)
  our$category <- cat_levels[match(tolower(our$category), tolower(cat_levels))]   # Boots/Sandals/Shoes/Slippers로 통일 (오타면 NA)
  our_mat <- t(sapply(file.path("ourshoes", our$file), load_any))
  rownames(our_mat) <- NULL
} else {
  our <- data.frame(file = character(), name = character(), category = character(), stringsAsFactors = FALSE)
  our_mat <- NULL
}

# img <- image_read("shoes/0001.png")
# x <- as.numeric(t(matrix(as.integer(image_data(img, channels = "gray")), 64, 64)))

# ============================================================
# Finding A / B 용 데이터 읽기 + 공통 보조 함수
#   CSV(kpop_lyrics_counts.csv)는 app.R 과 같은 폴더, 또는 그 안의 data 폴더에 둡니다.
# ============================================================
LYRICS_FILE <- if (file.exists("kpop_lyrics_counts.csv")) "kpop_lyrics_counts.csv" else "data/kpop_lyrics_counts.csv"
if (!file.exists(LYRICS_FILE)) stop("kpop_lyrics_counts.csv 가 없습니다. app.R 과 같은 폴더(또는 data 폴더)에 두세요.")

# Load K-pop lyrics data
lyrics <- read.csv(LYRICS_FILE)

# Remove rows with missing English percentages
lyrics <- lyrics %>%
  mutate(
    is_hive = tolower(trimws(as.character(is_hive)))
  ) %>%
  filter(!is.na(pct_main))

# 카운팅 규칙 3가지와 CSV 열 이름의 연결 (Finding A, B에서 사용)
RULE_COLS <- c(main = "pct_main", sub_1 = "pct_sub_1", sub_2 = "pct_sub_2")
if (!all(RULE_COLS %in% names(lyrics))) {
  stop("CSV에 없는 열: ", paste(setdiff(RULE_COLS, names(lyrics)), collapse = ", "))
}
RULE_CHOICES <- c("main: all English words" = "main",
                  "sub_1: without English fillers (ad-libs)" = "sub_1",
                  "sub_2: a repeated word counts once per song" = "sub_2")

# 곡 표를 "규칙 x 집단" 기준으로 만들어 주는 함수
#   d    : 연도 등으로 이미 걸러 둔 lyrics
#   rules: 선택한 규칙 이름 (예: c("main", "sub_2"))
#   split: TRUE 면 HYBE / Non-HYBE 로 나누고, FALSE 면 전체 곡을 하나로 봄
#   반환 : 규칙마다 곡을 한 번씩 복사한 표. pct = 그 규칙의 영어 비중, series = 그래프에서 선·막대를 나누는 이름
make_long <- function(d, rules, split) {
  if (split) {
    d$group <- ifelse(d$is_hive == "true", "HYBE", "Non-HYBE")
  } else {
    d$group <- rep("All songs", nrow(d))
  }
  long <- bind_rows(lapply(rules, function(r) {
    x <- d
    x$rule <- rep(r, nrow(x))
    x$pct  <- x[[RULE_COLS[[r]]]]
    x
  }))
  if (length(rules) > 1 && split) {
    long$series <- paste0(long$rule, " / ", long$group)
  } else if (length(rules) > 1) {
    long$series <- long$rule
  } else {
    long$series <- long$group
  }
  ord <- order(match(long$rule, names(RULE_COLS)), long$group)
  long$series <- factor(long$series, levels = unique(long$series[ord]))
  long
}

# 상단 통계 상자(대시보드). s = series별 요약표, thr 가 있으면 "기준 이상 곡 비율" 상자를 하나 더 보여 줌
stat_box <- function(title, value, w) {
  column(w, wellPanel(h4(title), h3(value)))
}
stat_rows <- function(s, thr = NULL) {
  w <- if (is.null(thr)) 4 else 3
  rows <- lapply(seq_len(nrow(s)), function(i) {
    boxes <- list(
      stat_box("Number of Songs", s$n[i], w),
      stat_box("Mean English Percentage", round(s$mean_pct[i], 2), w),
      stat_box("Median English Percentage", round(s$median_pct[i], 2), w)
    )
    if (!is.null(thr)) {
      boxes[[4]] <- stat_box(paste0("Songs with at least ", thr, "% English"),
                             paste0(round(s$share_thr[i], 1), "%"), w)
    }
    tagList(
      if (nrow(s) > 1) h5(as.character(s$series[i])),
      do.call(fluidRow, boxes)
    )
  })
  tagList(rows)
}

ui <- navbarPage("Team_chipmunk",
                 tabPanel("Shoe Search",
                          radioButtons("src", "Photo", c("Catalogue", "Self-check", "Upload"), inline = TRUE),
                          numericInput("cat_id", "Catalogue id", 1, 1, 2000),
                          selectInput("sc_file", "Self-check photo", sprintf("selfcheck_%02d.png", 1:4)),
                          fileInput("photo", "Upload", accept = c(".jpg", ".jpeg", ".png")),
                          sliderInput("K", "Number of PCs (K)", 1, 100, 20),
                          hr(),
                          h3("Prediction Result"),
                          textOutput("pred_text"),       # 예측 결과 및 정답 여부
                          h4("Top 5 Closest Shoes"),
                          tableOutput("top5_table"),     # Top 5 목록 (ID, 카테고리, 거리)
                          h4("Best Match Details"),
                          plotOutput("best_match_plot")   #[완료] TODO: best match, top 5, prediction, right or wrong
                 ),
                 tabPanel("Explain PCA",
                          p("For the photo chosen in Shoe Search, show how PCA finds the match."),
                          
                          h4("Step 1. Average shoe"), 
                          p("x̅, the mean of the 2,000 photos."),
                          plotOutput("pca_step1"),
                          
                          h4("Step 2. Photo − average"), 
                          p("x − x̅, where the photo is brighter or darker than the average shoe."),
                          plotOutput("pca_step2"),
                          
                          h4("Step 3. Score map"), 
                          p("All catalogue photos on PC1 and PC2, colored by category, with your photo marked."),
                          plotOutput("pca_step3"),
                          
                          h4("Step 4. Eigenshoes"), 
                          p("Show PC1, PC2, … as pictures (their loadings drawn as 64 × 64 images), so we can see what shape each PC captures, and how much variance each explains."),
                          plotOutput("pca_step4", height = "420px"),   # [수정 3] 격자(2x5)가 보이도록 높이 키움
                          
                          h4("Step 5. Scores"), 
                          p("The photo’s first K PC scores, the numbers used to find the match."),
                          tableOutput("pca_step5_table"),
                          
                          h4("Step 6. Rebuild"), 
                          p("The photo rebuilt from the average shoe and the first K eigenshoes, for several K."),
                          plotOutput("pca_step6"),
                          
                          h4("More"), 
                          p("Anything else you want to show or explain related to it."),
                          
                          # [수정 9] K별 예측 정확도
                          h5("Prediction accuracy for each K"),
                          p("How well does the app predict the shoe type for different numbers of PCs (K)? We test it on the 2,000 catalogue photos,",
                            "because only they have a known true category. Each photo is searched against the other 1,999 photos",
                            "(itself left out), the most common category among its top 5 is the prediction, and we compare it with the true category.",
                            "This is called leave-one-out testing."),
                          actionButton("run_acc", "Compute accuracy for K = 1 to 100 (about 10-30 seconds, only needed once)"),
                          plotOutput("acc_plot", height = "420px"),
                          textOutput("acc_text")
                 ),
                 tabPanel("Our Shoes",
                          p("Can your app recognize real shoes it has never seen?"),
                          p("Each member photographed 3 different shoes of their own (shoe only, plain background). The app has never seen these photos."),
                          p("Choose a shoe below to see its best match, its top 5 and whether the prediction is right. The number of PCs (K) is the one set in the Shoe Search tab."),
                          uiOutput("our_shoes_ui")     # [수정 11] 사진 목록, 선택 위젯, 결과가 여기에 표시됨
                 ),
                 tabPanel("Data 1",
                          sidebarLayout(
                            sidebarPanel(
                              width = 3,
                              tags$style(HTML("#d1_nav { flex-direction: column; }")),   # 왼쪽 메뉴를 세로로
                              
                              # ---- 왼쪽 위: Finding 탭 (내용은 오른쪽, 이 탭은 이름만 가진 메뉴) ----
                              tabsetPanel(
                                id = "d1_nav",
                                type = "pills",
                                tabPanel("Finding A"),
                                tabPanel("Finding B")
                              ),
                              hr(),
                              
                              # ---- 왼쪽 아래: 선택한 탭의 조정 위젯 ----
                              conditionalPanel(
                                condition = "input.d1_nav == 'Finding A'",
                                sliderInput(
                                  "year_range",
                                  "Select year range",
                                  min = min(lyrics$year),
                                  max = max(lyrics$year),
                                  value = c(min(lyrics$year), max(lyrics$year)),
                                  step = 1,
                                  sep = ""
                                ),
                                checkboxGroupInput(
                                  "a_rules",
                                  "Counting rule (choose one or more)",
                                  choices = RULE_CHOICES,
                                  selected = "main"
                                ),
                                checkboxInput("a_hive", "Compare HYBE vs Non-HYBE", value = FALSE),
                                checkboxInput("a_show_b", "Show the HYBE vs Non-HYBE comparison (at the bottom)", value = FALSE),
                                helpText("Two or more rules draw one line per rule. HYBE has only a few songs in some years, so read its line with care.")
                              ),
                              conditionalPanel(
                                condition = "input.d1_nav == 'Finding B'",
                                sliderInput(
                                  "c_thr",
                                  "English-centered = at least this share of English words",
                                  min = 10, max = 90, value = 50, step = 5, post = "%"
                                ),
                                checkboxInput("c_hive", "Compare HYBE vs Non-HYBE", value = FALSE),
                                helpText("Counting rule: main (all English words). The boxes cover all years, 2014-2025.")
                              )
                            ),
                            
                            mainPanel(
                              width = 9,
                              
                              # ================= Finding A =================
                              conditionalPanel(
                                condition = "input.d1_nav == 'Finding A'",
                                
                                h3("Has English Usage in K-pop Lyrics Increased Over Time?"),
                                
                                p("Explore how the proportion of English words in K-pop lyrics has changed from 2014 to 2025."),
                                
                                # Plain-language explanation of the three counting rules (click the title to fold it)
                                tags$details(
                                  open = "open",
                                  style = "margin-bottom: 20px;",
                                  tags$summary(tags$b("How do the counting rules work?")),
                                  p(style = "margin-top: 10px;",
                                    "For every song we count the Korean words (K) and the English words (E). The English share of a song is E divided by (K + E), times 100. The three rules differ only in which words are counted."),
                                  fluidRow(
                                    column(4, wellPanel(
                                      h4("main"),
                                      p(strong("Every word counts.")),
                                      p("Fillers such as 'oh' and 'yeah' are counted, and so is every repeat."),
                                      p(em("Use: the basic measure."))
                                    )),
                                    column(4, wellPanel(
                                      h4("sub_1"),
                                      p(strong("English fillers are left out.")),
                                      p("Ad-libs such as oh, yeah, la and lalala are not counted."),
                                      p(em("Answers: 'Do ad-libs make English look bigger than it is?'"))
                                    )),
                                    column(4, wellPanel(
                                      h4("sub_2"),
                                      p(strong("A repeated word counts once per song.")),
                                      p("Inside one song, a word that appears many times is counted a single time."),
                                      p(em("Answers: 'Is it just the chorus repeating the same English words?'"))
                                    ))
                                  ),
                                  h5("Example (an invented line, only to show how the counting works)"),
                                  p(tags$code("baby baby baby 사랑해 yeah yeah 사랑해")),
                                  tags$table(
                                    class = "table table-sm",
                                    tags$thead(tags$tr(tags$th("Rule"), tags$th("English words (E)"), tags$th("Korean words (K)"), tags$th("English share"))),
                                    tags$tbody(
                                      tags$tr(tags$td("main"),  tags$td("5  (baby x3, yeah x2)"), tags$td("2  (사랑해 x2)"), tags$td("5 / 7 = 71.4%")),
                                      tags$tr(tags$td("sub_1"), tags$td("3  (baby x3; yeah is a filler)"), tags$td("2  (사랑해 x2)"), tags$td("3 / 5 = 60.0%")),
                                      tags$tr(tags$td("sub_2"), tags$td("2  (baby, yeah; each once)"), tags$td("1  (사랑해, once)"), tags$td("2 / 3 = 66.7%"))
                                    )
                                  ),
                                  helpText("sub_1 and sub_2 count fewer words than main, so their English share is usually lower. What we compare is how it changes over the years.")
                                ),
                                
                                # Summary statistics (one row of boxes per selected line)
                                uiOutput("a_stats"),
                                
                                # Graph
                                h4("Key Finding"),
                                
                                h3("English Lyrics Do Not Determine Song Popularity in the Korean Market"),
                                
                                p("Our results suggest that the proportion of English lyrics alone cannot explain a song's popularity on the Melon chart. Since Melon is a Korean music streaming platform primarily used by Korean listeners, our findings mainly reflect music popularity and listening preferences in the Korean market."),
                                
                                p("Therefore, these results should not be directly generalized to international K-pop audiences."),
                                
                                hr(),
                                
                                h4("Finding A: English Lyrics Usage by Year"),
                                
                                plotOutput("year_plot", height = "400px"),
                                
                                # Interpretation
                                h4("Why Does It Matter?"),
                                p("This analysis helps us understand how the use of English in K-pop lyrics has changed over time."),
                                
                                h4("Next Analysis"),
                                p("We can further investigate whether English usage differs between different groups of artists."),
                                
                                # ---- HYBE vs Non-HYBE comparison (the old Finding B tab): shown only when its box is ticked ----
                                conditionalPanel(
                                  condition = "input.a_show_b",
                                  hr(),
                                  
                                  
                                  h3("Does English Usage Differ Between HYBE and Non-HYBE Artists?"),
                                  
                                  # Why did we choose HYBE?
                                  h4("Why Did We Choose HYBE?"),
                                  
                                  p("HYBE is one of the major entertainment companies in the K-pop industry, with artists who have gained significant international recognition."),
                                  
                                  p("As K-pop continues to expand globally, we became curious about whether English usage in song lyrics differs between artists affiliated with HYBE and those from other entertainment companies."),
                                  
                                  wellPanel(
                                    h4("Our Research Question"),
                                    p("Does English usage differ between HYBE and Non-HYBE artists?")
                                  ),
                                  
                                  p("To answer this question, we divided the songs into two groups: HYBE and Non-HYBE, and compared the proportion of English words in their lyrics."),
                                  
                                  # Input: Select artist group
                                  h4("Explore the Data"),
                                  
                                  selectInput(
                                    "hive_group",
                                    "Select Group",
                                    choices = c(
                                      "All Groups" = "all",
                                      "HYBE" = "true",
                                      "Non-HYBE" = "false"
                                    ),
                                    selected = "all"
                                  ),
                                  
                                  # Summary statistics
                                  fluidRow(
                                    column(
                                      4,
                                      wellPanel(
                                        h4("Number of Songs"),
                                        h3(textOutput("hive_count"))
                                      )
                                    ),
                                    column(
                                      4,
                                      wellPanel(
                                        h4("Mean English Percentage"),
                                        h3(textOutput("hive_mean"))
                                      )
                                    ),
                                    column(
                                      4,
                                      wellPanel(
                                        h4("Median English Percentage"),
                                        h3(textOutput("hive_median"))
                                      )
                                    )
                                  ),
                                  
                                  # Graph
                                  h4("English Usage by Group"),
                                  plotOutput("hive_plot", height = "400px"),
                                  
                                  # Interpretation
                                  h4("Why Does It Matter?"),
                                  
                                  p("This comparison helps us explore whether English usage is associated with an artist's entertainment company affiliation. It may also provide insights into different approaches to the global K-pop market."),
                                  
                                  # Next Analysis
                                  h4("Next Analysis"),
                                  
                                  p("We can further investigate whether English usage is related to chart rankings and whether these patterns are consistent across different years.")
                                )
                              ),
                              
                              # ================= Finding B =================
                              conditionalPanel(
                                condition = "input.d1_nav == 'Finding B'",
                                
                                h3("How Many Songs Are English-Centered?"),
                                
                                p("Explore the share of songs whose lyrics are mostly English, year by year. You decide how much English a song needs to count as English-centered."),
                                
                                # Summary statistics (one row of boxes per selected group)
                                uiOutput("c_stats"),
                                
                                h4("Share of English-Centered Songs by Year"),
                                plotOutput("c_plot", height = "400px"),
                                
                                h4("Why Does It Matter?"),
                                
                                p("A rising average could mean that every song uses a little more English, or that more songs are written mostly in English. Counting the songs above a cut-off separates the two. With the default 50% cut-off, English-centered songs grew from about 4% of the chart in 2014 to about 41% in 2025."),
                                
                                p("The cut-off changes the story. With a stricter 70% cut-off, English-centered songs stay at or below 4% through 2022 and reach about 28% in 2025, so the growth is recent."),
                                
                                h4("Next Analysis"),
                                
                                p("We can look at which artists and songs are behind the increase, and compare it with our other datasets over 2022-2025, the years when the change is steepest.")
                              )
                            )
                          )
                 ),
                 tabPanel("AI",
                          h3("Ethics Statement"),
                          tags$blockquote("“We understand the rules on data and AI use for this project, including following the license and privacy rules of every dataset we used, understanding all of our code, and knowing that presenting AI output as our own work may be academic misconduct. We certify that we have followed all of them. Any violation is our responsibility, and we agree to receive an F in this course if we have not.”"),
                          
                          hr(),
                          h3("What AI did and what you did"),
                          p("A short summary for each tab of the project."),
                          tableOutput("ai_summary_table"),
                          
                          hr(),
                          h3("Individual Reflections"),
                          p("Each member's reflection on AI assistance, failures, and learnings (including the required keyword)."),
                          tableOutput("ai_reflection_table")
                 ),
                 
                 theme = bs_theme(
                   version = 5,
                   bootswatch = "minty",
                   primary = "#8064D8",
                   secondary = "#B9A5F5",
                   bg = "#F8F7FC",
                   fg = "#29243A"
                 ),
                 
                 header = tags$head(
                   tags$style(HTML("
      body {
        background-color: #F8F7FC;
        font-family: 'Arial', sans-serif;
      }

      .navbar {
        background-color: #FFFFFF !important;
        border-bottom: 1px solid #EAE5F5;
        padding: 12px;
      }

      .navbar-brand {
        color: #7153C6 !important;
        font-weight: bold;
        font-size: 22px;
      }

      .navbar-nav .nav-link {
        color: #625B70 !important;
        font-weight: 600;
        border-radius: 10px;
        margin: 3px;
      }

      .navbar-nav .nav-link.active {
        background-color: #8064D8 !important;
        color: white !important;
      }

      h3 {
        color: #6047A8;
        font-weight: bold;
        margin-top: 20px;
      }

      h4 {
        color: #51466D;
        font-weight: bold;
      }

      .well {
        background-color: white;
        border: 1px solid #EAE5F5;
        border-radius: 15px;
        box-shadow: 0 3px 12px rgba(80, 60, 130, 0.06);
        padding: 20px;
      }

      .btn {
        background-color: #8064D8;
        border-color: #8064D8;
        border-radius: 10px;
      }

      .form-control {
        border-radius: 10px;
        border-color: #DDD5EF;
      }
    "))
                 )
)

server <- function(input, output, session) {
  #tab1
  # 1. 사용자가 선택한 사진의 픽셀 벡터를 가져오는 반응형 함수 (Reactive)
  selected_shoe_vec <- reactive({
    src_type <- input$src
    
    if (src_type == "Catalogue") {
      # 카탈로그에서 고른 경우: 해당 ID의 행 데이터 가져오기 (1부터 시작하므로 그대로 인덱스 사용)
      idx <- input$cat_id
      return(shoe_matrix[idx, ])
      
    } else if (src_type == "Self-check") {
      # 셀프체크 사진인 경우
      return(load_any(file.path("selfcheck", input$sc_file)))   # [수정 1] 공통 변환 함수 사용
      
    } else {
      # 파일을 직접 업로드한 경우
      req(input$photo) # 업로드 파일이 있을 때만 실행
      # [수정 1] 정사각형으로 맞추고(흰색으로 채움) 64x64 흑백으로 변환하는 공통 함수 사용
      return(load_any(input$photo$datapath))
    }
  })
  
  # 2. PCA를 이용해 가장 가까운 Top 5 찾기 및 거리 계산
  search_results <- reactive({
    query_vec <- selected_shoe_vec()
    K <- input$K
    
    # 쿼리 이미지의 평균 중심 맞추기
    query_centered <- query_vec - x_bar
    
    # 선택된 K개의 주성분(Eigenshoes)만 사용하여 쿼리의 PC 스코어 계산
    query_scores <- query_centered %*% pca_result$rotation[, 1:K]
    
    # 카탈로그 전체의 PC 스코어 (이미 계산된 pca_result$x 활용)
    catalog_scores <- pca_result$x[, 1:K]
    
    # 유클리드 거리(Euclidean Distance) 계산
    distances <- sqrt(rowSums((catalog_scores - matrix(query_scores, nrow = 2000, ncol = K, byrow = TRUE))^2))
    
    # 만약 카탈로그 ID로 검색한 경우, 자기 자신과의 거리를 무한대로 만들어 검색 결과에서 제외
    if (input$src == "Catalogue") {
      distances[input$cat_id] <- Inf
    }
    
    # 거리가 가장 가까운 상위 5개 인덱스 추출
    top5_idx <- order(distances)[1:5]
    
    # 결과 데이터프레임 만들기
    res_df <- data.frame(
      Rank = 1:5,
      ID = top5_idx,
      Category = meta$category[top5_idx],
      Distance = round(distances[top5_idx], 2)
    )
    
    return(res_df)
  })
  
  # 3. 화면 출력 1: Top 5 표 (Table)
  output$top5_table <- renderTable({
    search_results()
  })
  
  # 4. 화면 출력 2: 예측 결과 (Text)
  output$pred_text <- renderText({
    top5 <- search_results()
    # Top 5 중 가장 많이 나온 카테고리를 예측값으로 선정
    # [수정 5] 동점이면 가장 가까운 사진의 카테고리를 선택
    pred_cat <- predict_cat(top5)
    
    # 만약 카탈로그 검색이라면 실제 정답 카테고리와 비교
    if (input$src == "Catalogue") {
      true_cat <- meta$category[input$cat_id]
      is_correct <- ifelse(pred_cat == true_cat, "Right", "Wrong")
      paste0("Predicted category: ", pred_cat, " | Actual category: ", true_cat, " -> ", is_correct)
    } else {
      paste0("Predicted category: ", pred_cat)
    }
  })
  
  # 5. 화면 출력 3: 베스트 매치 시각화 (Plot - 원본, 평균, 복원 이미지 비교)
  output$best_match_plot <- renderPlot({
    top5 <- search_results()
    best_id <- top5$ID[1] # 가장 가까운 1등 신발 ID
    K <- input$K
    
    # 1등 신발의 원본 벡터를 64x64 이미지로 변환
    best_mat <- matrix(shoe_matrix[best_id, ], 64, 64, byrow = TRUE)
    
    # K개의 PC로 복원된 이미지 계산
    mean_mat <- matrix(x_bar, 64, 64, byrow = TRUE)
    query_centered <- selected_shoe_vec() - x_bar
    scores <- query_centered %*% pca_result$rotation[, 1:K]
    rebuilt_vec <- x_bar + scores %*% t(pca_result$rotation[, 1:K])
    rebuilt_mat <- matrix(rebuilt_vec, 64, 64, byrow = TRUE)
    
    # 화면에 3개의 이미지를 나란히 그리기 (원본, 평균, 복원)
    par(mfrow = c(1, 3), mar = c(2, 2, 2, 2))
    image(t(apply(best_mat, 2, rev)), col = gray(0:100/100), axes = FALSE, asp = 1, main = paste0("Best Match (ID: ", best_id, ")"))
    image(t(apply(mean_mat, 2, rev)), col = gray(0:100/100), axes = FALSE, asp = 1, main = "Average Shoe")
    image(t(apply(rebuilt_mat, 2, rev)), col = gray(0:100/100), axes = FALSE, asp = 1, main = paste0("Rebuilt (K=", K, ")"))
  })
  
  # AI 탭 테이블 (TODO 제거 및 예시 채우기)
  output$ai_table <- renderTable(data.frame(
    Tab = c("Shoe Search", "Explain PCA", "Our Shoes", "Data 1"),
    #수정하기!! TODO
    `What AI did` = c("Implemented the distance calculation and PCA matching algorithm", "Helped with the step-by-step PCA visualizations", "Connected our own shoe data", "Wrote data analysis queries"),
    `What we did` = c("UI design and final verification", "Layout arrangement", "Photographed and uploaded our own shoes", "Drew business insights"),
    check.names = FALSE))
  
  # --- Part 3: AI 탭 출력 로직 ---
  
  # 1. 탭별 AI vs 인간 역할 분담 요약 표
  output$ai_summary_table <- renderTable({
    data.frame(
      Tab = c("Shoe Search", "Explain PCA", "Our Shoes", "Finding A", "Finding B"),
      `What AI did` = c(
        "Wrote the skeleton code for the distance calculation and PCA matching algorithm",
        "Helped with the code for the 6-step PCA visualization (average, deviation, score map, eigenshoes, etc.)",
        "Designed the upload and preprocessing logic for personal shoe images",
        "TODO: write what AI really did for Finding A",
        "TODO: write what AI really did for Finding B"
      ),
      `What we did` = c(
        "UI layout design and final result verification",
        "Refined graph layouts and text explanations",
        "Photographed our own shoes and managed the directory structure",
        "TODO: write what we did for Finding A",
        "TODO: write what we did for Finding B"
      ),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)
  
  # 2. 팀원별 소감 표 
  output$ai_reflection_table <- renderTable({
    data.frame(
      Member = c("박선우", "송현정", "한우영", "신유진", "Flo"),
      #수정하기 TODO, 각자 역할 쓰기
      Reflection = c(
        "AI was a great help in writing the skeleton code for complex matrix operations and distance calculations. However, I studied the formulas and code line by line to understand the concept of dimensionality reduction properly.",
        "AI provided an excellent layout structure while we visualized the 6 steps of PCA. I edited details such as grayscale image handling myself and made the code my own.",
        "AI was useful for debugging errors that came up during data preprocessing and text handling. It shortened overall development time, so I could focus more on interpreting and verifying the data.",
        "I worked out how to connect the UI with the reactive server logic together with AI. It was a reliable helper for turning a complex algorithm into code.",
        "Together with my teammates, I reviewed and verified the code AI suggested. I built my development skills by fixing the exception handling and detailed settings that AI missed."
      ),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)
  
  # --- Tab 2: Explain PCA 시각화 출력들 ---
  
  # Step 1. 평균 신발 이미지
  output$pca_step1 <- renderPlot({
    mean_mat <- matrix(x_bar, 64, 64, byrow = TRUE)
    image(t(apply(mean_mat, 2, rev)), col = gray(0:100/100), axes = FALSE, asp = 1, main = "Average Shoe (x̄)")
  })
  
  # Step 2. 사진 - 평균 (편차 이미지)
  output$pca_step2 <- renderPlot({
    query_vec <- selected_shoe_vec()
    diff_vec <- query_vec - x_bar
    diff_mat <- matrix(diff_vec, 64, 64, byrow = TRUE)
    # 편차는 음수/양수가 있으므로 색상 팔레트 조정 가능 (여기선 흑백으로 표현)
    image(t(apply(diff_mat, 2, rev)), col = gray(0:100/100), axes = FALSE, asp = 1, main = "Photo − Average (x - x̄)")
  })
  
  # Step 3. Score map (PC1 vs PC2 산점도 및 현재 사진 마킹)
  output$pca_step3 <- renderPlot({
    # 카탈로그 전체의 PC1, PC2 스코어 가져오기
    pc1 <- pca_result$x[, 1]
    pc2 <- pca_result$x[, 2]
    
    # 현재 선택한 사진의 PC 스코어 계산
    query_vec <- selected_shoe_vec()
    query_centered <- query_vec - x_bar
    query_pc1 <- sum(query_centered * pca_result$rotation[, 1])
    query_pc2 <- sum(query_centered * pca_result$rotation[, 2])
    
    # [수정 2] 카테고리별 색상 (반투명, 검정·노랑은 쓰지 않음 -> 내 사진 표시와 겹치지 않음)
    cat_cols <- c(Boots = "#1f77b4", Sandals = "#2ca02c", Shoes = "#ff7f0e", Slippers = "#9467bd")
    cat_cols <- cat_cols[sort(unique(meta$category))]
    
    plot(pc1, pc2, col = adjustcolor(cat_cols[meta$category], alpha.f = 0.45), pch = 16, cex = 0.8,
         main = "Score Map (PC1 vs PC2)", xlab = "PC1 Score", ylab = "PC2 Score")
    
    # [수정 2] 내 사진: 카테고리에 쓰지 않은 색(검정 바탕 + 노란 테두리)의 큰 마름모 + 글자
    points(query_pc1, query_pc2, pch = 23, bg = "black", col = "yellow", cex = 3, lwd = 3)
    text(query_pc1, query_pc2, labels = "Your photo", pos = 3, offset = 1.2, font = 2, cex = 1.1)
    
    legend("topright", legend = c(names(cat_cols), "Your photo"),
           col = c(cat_cols, "yellow"), pt.bg = c(rep(NA, length(cat_cols)), "black"),
           pch = c(rep(16, length(cat_cols)), 23), pt.cex = c(rep(1, length(cat_cols)), 1.6),
           cex = 0.8, bg = "white")
  })
  
  # Step 4. Eigenshoes (PC1, PC2의 고유벡터 이미지)
  output$pca_step4 <- renderPlot({
    # [수정 3] 처음 10개 PC를 2x5 격자로, 제목에 설명 분산(%) 표시, 부호를 보여주는 파랑-흰색-빨강 팔레트
    n_show <- 10
    var_exp <- pca_result$sdev^2 / sum(pca_result$sdev^2)       # 각 PC의 설명 분산 비율
    pal <- colorRampPalette(c("blue", "white", "red"))(101)
    
    par(mfrow = c(2, 5), mar = c(1, 1, 3, 1))
    for (k in 1:n_show) {
      v <- pca_result$rotation[, k]
      m <- matrix(v, 64, 64, byrow = TRUE)
      lim <- max(abs(v))                                        # 0이 흰색이 되도록 좌우 대칭 범위
      image(t(apply(m, 2, rev)), col = pal, zlim = c(-lim, lim), axes = FALSE, asp = 1,
            main = sprintf("PC%d (%.1f%%)", k, 100 * var_exp[k]))
    }
  })
  
  # Step 5. Scores (선택한 사진의 첫 K개 PC 스코어 표)
  output$pca_step5_table <- renderTable({
    query_vec <- selected_shoe_vec()
    query_centered <- query_vec - x_bar
    K <- input$K          # [수정 4] 슬라이더 K개 전부 표시 (이전: 최대 10개)
    scores <- query_centered %*% pca_result$rotation[, 1:K]
    
    data.frame(
      PC = paste0("PC", 1:K),
      Score = round(as.vector(scores), 2)
    )
  }, striped = TRUE, hover = TRUE)
  
  # Step 6. Rebuild (다양한 K에 따른 복원 이미지 비교)
  output$pca_step6 <- renderPlot({
    query_vec <- selected_shoe_vec()
    query_centered <- query_vec - x_bar
    
    # [수정 6] 비교할 K 값들: 고정값(1, 5, 20, 50) + 슬라이더로 고른 K (중복이면 한 번만, 작은 순서로 정렬)
    k_values <- sort(unique(c(1, 5, 20, min(50, ncol(pca_result$rotation)), input$K)))
    
    par(mfrow = c(1, length(k_values)), mar = c(2, 2, 3, 2))
    for (k in k_values) {
      scores <- query_centered %*% pca_result$rotation[, 1:k]
      rebuilt_vec <- x_bar + scores %*% t(pca_result$rotation[, 1:k])
      rebuilt_mat <- matrix(rebuilt_vec, 64, 64, byrow = TRUE)
      
      is_sel <- (k == input$K)                   # 내가 고른 K인지
      image(t(apply(rebuilt_mat, 2, rev)), col = gray(0:100/100), axes = FALSE, asp = 1,
            main = paste0("Rebuilt (K=", k, ")", if (is_sel) "\n(your K)" else ""),
            col.main = if (is_sel) "red" else "black")   # 내가 고른 K는 제목을 빨간색으로 강조
    }
  })
  
  # [수정 9] K별 예측 정확도 (버튼을 누르면 계산, 결과는 저장되어 다시 계산하지 않음)
  acc_result <- eventReactive(input$run_acc, {
    withProgress(message = "Computing accuracy for K = 1 to 100", value = 0, {
      compute_loo_acc(function(k) incProgress(1/100, detail = paste("K =", k)))
    })
  })
  
  output$acc_plot <- renderPlot({
    validate(need(!is.null(input$run_acc) && input$run_acc > 0, "Click the button above to compute the accuracy."))
    r <- acc_result()
    cat_cols <- c(Boots = "#1f77b4", Sandals = "#2ca02c", Shoes = "#ff7f0e", Slippers = "#9467bd")
    cats <- colnames(r$by_cat)
    
    par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
    plot(r$K, 100 * r$overall, type = "l", lwd = 3, ylim = c(0, 100),
         xlab = "Number of PCs (K)", ylab = "Prediction accuracy (%)",
         main = "Leave-one-out accuracy on the 2,000 catalogue photos")
    for (j in seq_along(cats))                                   # categories: thin dashed lines
      lines(r$K, 100 * r$by_cat[, j], col = cat_cols[cats[j]], lwd = 1.5, lty = 2)
    abline(v = input$K, col = "red", lty = 3)                    # selected K
    points(input$K, 100 * r$overall[input$K], pch = 19, col = "red", cex = 1.6)
    
    # [수정 10] 정확도가 가장 높은 K를 별(★)과 글자로 표시 (같은 값이 여러 개면 가장 작은 K)
    best     <- which.max(r$overall)
    best_acc <- 100 * r$overall[best]
    points(best, best_acc, pch = 8, col = "black", cex = 2.2, lwd = 2)
    text(best, best_acc, labels = sprintf("Best K = %d (%.1f%%)", best, best_acc),
         pos = if (best_acc > 90) 1 else 3, offset = 1, font = 2, cex = 1)
    
    legend("bottomright", bg = "white", cex = 0.85,
           legend = c("Overall", cats, paste0("Selected K = ", input$K), paste0("Best K = ", best)),
           col = c("black", cat_cols[cats], "red", "black"),
           lwd = c(3, rep(1.5, length(cats)), NA, NA),
           lty = c(1, rep(2, length(cats)), NA, NA),
           pch = c(NA, rep(NA, length(cats)), 19, 8))
  })
  
  output$acc_text <- renderText({
    req(input$run_acc)
    r <- acc_result()
    best <- which.max(r$overall)
    sprintf("At the selected K = %d, the app predicts the correct category for %.1f%% of the catalogue photos (%d of 2,000). The highest accuracy is %.1f%% at K = %d.",
            input$K, 100 * r$overall[input$K], round(r$overall[input$K] * 2000),
            100 * r$overall[best], best)
  })
  
  # [수정 11] Our Shoes 탭 ---------------------------------------------------------
  output$our_shoes_ui <- renderUI({
    if (nrow(our) == 0) {
      return(p(strong("No photos found. Put the photos in a folder named 'ourshoes' next to app.R, named Name_Category_Number.jpg (e.g., Whae_Sandals_1.jpg), then restart the app.")))
    }
    tagList(
      h4("All our shoes"),
      p("Each photo is shown as the app sees it (squared to 64 × 64, grayscale), with the member's name and the shoe's category."),
      plotOutput("our_grid", height = "auto"),
      hr(),
      selectInput("our_pick", "Choose a shoe",
                  choices = setNames(our$file, paste0(our$name, " - ", our$category, " (", our$file, ")")),
                  width = "500px"),
      textOutput("our_pred"),
      plotOutput("our_best_plot", height = "320px"),
      h4("Top 5 Closest Shoes"),
      tableOutput("our_top5_table"),
      plotOutput("our_top5_plot", height = "220px"),
      hr(),
      h4("Result for all our shoes at the selected K"),
      textOutput("our_summary"),
      tableOutput("our_all_table")
    )
  })
  
  # 전체 사진 목록 (6열 격자, 행 수에 맞춰 높이가 자동으로 정해짐)
  output$our_grid <- renderPlot({
    n  <- nrow(our)
    nc <- min(n, 6)
    nr <- ceiling(n / nc)
    par(mfrow = c(nr, nc), mar = c(1, 1, 3, 1))
    for (i in seq_len(n)) draw_gray(our_mat[i, ], main = paste0(our$name[i], "\n", our$category[i]))
  }, height = function() 200 * ceiling(max(nrow(our), 1) / min(max(nrow(our), 1), 6)))
  
  our_idx <- reactive({ req(input$our_pick); match(input$our_pick, our$file) })
  our_res <- reactive({ find_top5(our_mat[our_idx(), ], input$K) })
  
  output$our_pred <- renderText({
    pred   <- predict_cat(our_res())
    actual <- our$category[our_idx()]
    verdict <- if (is.na(actual)) "Actual category unknown (check the file name)"
    else if (pred == actual) "Right" else "Wrong"
    paste0("Predicted category: ", pred, " | Actual category: ", actual, " -> ", verdict)
  })
  
  output$our_best_plot <- renderPlot({
    i <- our_idx(); r <- our_res()
    par(mfrow = c(1, 2), mar = c(2, 1, 3, 1))
    draw_gray(our_mat[i, ], main = paste0(our$name[i], " - ", our$category[i]))
    draw_gray(shoe_matrix[r$ID[1], ], main = "Closest tiwn",       # 제목은 과제 문구 그대로 (오타 아님, 수정 금지)
              sub = paste0("ID ", r$ID[1], " - ", r$Category[1]))
  })
  
  output$our_top5_table <- renderTable({ our_res() })
  
  output$our_top5_plot <- renderPlot({
    r <- our_res()
    par(mfrow = c(1, 5), mar = c(1, 1, 3, 1))
    for (j in 1:5) draw_gray(shoe_matrix[r$ID[j], ], main = paste0("#", j, "  ID ", r$ID[j], "\n", r$Category[j]))
  })
  
  our_all <- reactive({
    K <- input$K
    pred <- vapply(seq_len(nrow(our)), function(i) predict_cat(find_top5(our_mat[i, ], K)), character(1))
    data.frame(Name = our$name, File = our$file, Actual = our$category, Predicted = pred,
               Result = ifelse(is.na(our$category), "Unknown", ifelse(pred == our$category, "Right", "Wrong")),
               stringsAsFactors = FALSE)
  })
  
  output$our_summary <- renderText({
    d <- our_all()
    known <- d$Result != "Unknown"
    sprintf("K = %d: %d of %d shoes predicted correctly.", input$K, sum(d$Result == "Right"), sum(known))
  })
  
  output$our_all_table <- renderTable({ our_all() })
  
  # ==================================
  # Finding A: English Usage Over Time
  # ==================================
  
  # 선택한 연도 범위 + 규칙(여러 개 가능) + HYBE 선택을 반영한 표
  a_long <- reactive({
    validate(need(length(input$a_rules) > 0, "Please choose at least one counting rule."))
    rules <- names(RULE_COLS)[names(RULE_COLS) %in% input$a_rules]
    d <- lyrics %>%
      filter(
        year >= input$year_range[1],
        year <= input$year_range[2]
      )
    make_long(d, rules, isTRUE(input$a_hive))
  })
  
  output$a_stats <- renderUI({
    long <- a_long()
    validate(need(nrow(long) > 0, "No songs match the selected options."))
    s <- long %>%
      group_by(series) %>%
      summarise(
        n = n(),
        mean_pct = mean(pct),
        median_pct = median(pct),
        .groups = "drop"
      )
    stat_rows(s)
  })
  
  output$year_plot <- renderPlot({
    long <- a_long()
    validate(need(nrow(long) > 0, "No songs match the selected options."))
    plot_data <- long %>%
      group_by(year, series) %>%
      summarise(
        mean_pct = mean(pct),
        .groups = "drop"
      )
    
    ggplot(plot_data, aes(x = year, y = mean_pct, colour = series, group = series)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      scale_x_continuous(breaks = sort(unique(plot_data$year))) +
      labs(
        title = "Average English Percentage in K-pop Lyrics",
        x = "Year",
        y = "English Percentage (%)",
        colour = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = if (length(unique(plot_data$series)) > 1) "bottom" else "none")
  })
  
  # ==================================
  # HYBE vs Non-HYBE comparison (shown inside Finding A)
  # ==================================
  
  hive_data <- reactive({
    if (input$hive_group == "all") {
      return(lyrics)
    }
    
    lyrics %>%
      filter(is_hive == input$hive_group)
  })
  
  output$hive_count <- renderText({
    nrow(hive_data())
  })
  
  output$hive_mean <- renderText({
    round(mean(hive_data()$pct_main), 2)
  })
  
  output$hive_median <- renderText({
    round(median(hive_data()$pct_main), 2)
  })
  
  output$hive_plot <- renderPlot({
    plot_data <- lyrics %>%
      group_by(is_hive) %>%
      summarise(
        mean_pct = mean(pct_main),
        .groups = "drop"
      ) %>%
      mutate(
        group = ifelse(is_hive == "true", "HYBE", "Non-HYBE")
      )
    
    ggplot(plot_data, aes(x = group, y = mean_pct, fill = group)) +
      geom_col(width = 0.6) +
      geom_text(
        aes(label = paste0(round(mean_pct, 2), "%")),
        vjust = -0.5,
        size = 5
      ) +
      labs(
        title = "Average English Percentage by Group",
        x = "Group",
        y = "Mean English Percentage (%)"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })
  
  # ==================================
  # Finding B: English-centered songs by year (bar chart)
  # ==================================
  
  c_long <- reactive({
    make_long(lyrics, "main", isTRUE(input$c_hive))
  })
  
  output$c_stats <- renderUI({
    long <- c_long()
    validate(need(nrow(long) > 0, "No songs match the selected options."))
    thr <- input$c_thr
    s <- long %>%
      group_by(series) %>%
      summarise(
        n = n(),
        mean_pct = mean(pct),
        median_pct = median(pct),
        share_thr = mean(pct >= thr) * 100,
        .groups = "drop"
      )
    stat_rows(s, thr = thr)
  })
  
  output$c_plot <- renderPlot({
    long <- c_long()
    validate(need(nrow(long) > 0, "No songs match the selected options."))
    thr <- input$c_thr
    plot_data <- long %>%
      group_by(year, series) %>%
      summarise(
        share = mean(pct >= thr) * 100,
        .groups = "drop"
      )
    plot_data$year_f <- factor(plot_data$year)
    
    ggplot(plot_data, aes(x = year_f, y = share, fill = series)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(
        aes(label = paste0(round(share, 1), "%")),
        position = position_dodge(width = 0.8),
        vjust = -0.5,
        size = 3.3
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(
        title = paste0("Songs with at least ", thr, "% English Words, by Year"),
        x = "Year",
        y = "Share of songs (%)",
        fill = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = if (length(unique(plot_data$series)) > 1) "bottom" else "none")
  })
}


shinyApp(ui, server)