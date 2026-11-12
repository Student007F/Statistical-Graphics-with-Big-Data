# Project 2 template: the tabs and the inputs only. Everything else, and the look of the app, is yours.
library(shiny)

# Example: read one photo as 4,096 numbers (64 x 64 grayscale, row by row, 0 = black to 255 = white)
library(magick)

# [수정 12] 업로드 최대 용량 30MB (기본값 5MB는 휴대폰 사진에 작을 수 있음)
options(shiny.maxRequestSize = 30 * 1024^2)

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
    overall[k] <- mean(ok)
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
    # [수정 12] 앱에서 우리 사진 업로드 (여러 장 선택 가능, 파일 이름 규칙: 이름_카테고리_번호.jpg)
    fileInput("our_upload", "Add our photos (select one or several; name them Name_Category_Number, e.g., Whae_Sandals_1.jpg)",
              multiple = TRUE, accept = c(".jpg", ".jpeg", ".png")),
    p("Photos added here stay only until the page is refreshed. To keep them permanently, put the files in the 'ourshoes' folder next to app.R."),
    uiOutput("our_upload_msg"),
    actionButton("our_reset", "Remove uploaded photos"),
    # [수정 13] 선택한(방금 올린) 사진의 분석을 맨 위에, 전체 사진 목록과 전체 결과는 그 아래에 표시
    uiOutput("our_empty_msg"),
    h4("Analysis of the selected shoe"),
    p("The photo you just added is selected automatically. Choose another shoe from the list to analyse it."),
    selectInput("our_pick", "Choose a shoe", choices = character(0), width = "500px"),
    textOutput("our_pred"),
    plotOutput("our_best_plot", height = "320px"),
    h4("Top 5 Closest Shoes"),
    tableOutput("our_top5_table"),
    plotOutput("our_top5_plot", height = "220px"),
    hr(),
    h4("All our shoes"),
    p("Each photo is shown as the app sees it (squared to 64 x 64, grayscale), with the member's name and the shoe's category."),
    plotOutput("our_grid", height = "auto"),
    hr(),
    h4("Result for all our shoes at the selected K"),
    textOutput("our_summary"),
    tableOutput("our_all_table")
  ),
  tabPanel("Data 1",
    p("One Project 1 dataset, descriptive results only (filter, sort, count, mean, median, rate, group summary)."),
    p("It should include:"),
    tags$ul(
      tags$li("the question you are trying to answer;"),
      tags$li("the inputs that show it (e.g., the filter, grouping or time range you choose);"),
      tags$li("the result in numbers;"),
      tags$li("why it matters;"),
      tags$li("how this finding helps your next analysis of this dataset."))
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
      Tab = c("Shoe Search", "Explain PCA", "Our Shoes", "Data 1"),
      `What AI did` = c(
        "Wrote the skeleton code for the distance calculation and PCA matching algorithm",
        "Helped with the code for the 6-step PCA visualization (average, deviation, score map, eigenshoes, etc.)",
        "Designed the upload and preprocessing logic for personal shoe images",
        "Suggested ideas for descriptive statistics (filter, sort, group summary)"
      ),
      `What we did` = c(
        "UI layout design and final result verification",
        "Refined graph layouts and text explanations",
        "Photographed our own shoes and managed the directory structure",
        "Chose the analysis topic and drew business insights"
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

  # [수정 11, 12] Our Shoes 탭 ---------------------------------------------------------
  # ours: 화면에 보이는 우리 신발 사진들. 처음에는 ourshoes/ 폴더의 사진이고, 업로드하면 추가됨 (새로고침하면 사라짐)
  ours <- reactiveValues(df = our, mat = our_mat, msg = character(), pick = NULL)

  # 업로드한 사진 추가: 파일 이름(이름_카테고리_번호)에서 이름과 카테고리를 읽음
  observeEvent(input$our_upload, {
    up <- input$our_upload
    new_df <- NULL; new_mat <- NULL; msgs <- character()
    for (i in seq_len(nrow(up))) {
      parts  <- strsplit(tools::file_path_sans_ext(up$name[i]), "_")[[1]]
      person <- parts[1]
      cat_i  <- if (length(parts) >= 2) cat_levels[match(tolower(parts[2]), tolower(cat_levels))] else NA_character_
      v <- tryCatch(load_any(up$datapath[i]), error = function(e) NULL)
      if (is.null(v)) {
        msgs <- c(msgs, paste0(up$name[i], ": could not be read as an image, so it was skipped."))
        next
      }
      if (is.na(cat_i))
        msgs <- c(msgs, paste0(up$name[i], ": no category found in the file name (expected Name_Category_Number, ",
                               "with Category one of ", paste(cat_levels, collapse = ", "),
                               "), so its true category is unknown."))
      new_df  <- rbind(new_df, data.frame(file = up$name[i], name = person, category = cat_i, stringsAsFactors = FALSE))
      new_mat <- rbind(new_mat, v)
    }
    if (!is.null(new_df)) {
      keep    <- !(ours$df$file %in% new_df$file)          # 같은 이름의 파일을 다시 올리면 교체
      old_mat <- if (is.null(ours$mat)) NULL else ours$mat[keep, , drop = FALSE]
      ours$mat  <- unname(rbind(old_mat, new_mat))
      ours$df   <- rbind(ours$df[keep, , drop = FALSE], new_df)
      ours$pick <- new_df$file[1]                          # 방금 올린 첫 사진을 자동 선택
      msgs <- c(sprintf("Added %d photo(s).", nrow(new_df)), msgs)
    }
    ours$msg <- msgs
  })

  # 업로드한 사진 지우기 (ourshoes/ 폴더의 사진만 남김)
  observeEvent(input$our_reset, {
    ours$df <- our; ours$mat <- our_mat; ours$msg <- character(); ours$pick <- NULL
  })

  output$our_upload_msg <- renderUI({ tagList(lapply(ours$msg, function(m) p(m))) })

  # 사진이 하나도 없을 때 안내 문구
  output$our_empty_msg <- renderUI({
    if (nrow(ours$df) == 0)
      p(strong("No photos yet. Add photos above, or put them in a folder named 'ourshoes' next to app.R (named Name_Category_Number.jpg) and restart the app."))
  })

  # 드롭다운 목록: 사진이 추가/삭제될 때마다 갱신하고, 방금 올린 첫 사진을 선택
  observe({
    d <- ours$df
    ch <- if (nrow(d) == 0) character(0) else setNames(d$file, paste0(d$name, " - ", d$category, " (", d$file, ")"))
    updateSelectInput(session, "our_pick", choices = ch, selected = isolate(ours$pick))
  })

  # 전체 사진 목록 (6열 격자, 행 수에 맞춰 높이가 자동으로 정해짐)
  output$our_grid <- renderPlot({
    req(nrow(ours$df) > 0)
    n  <- nrow(ours$df)
    nc <- min(n, 6)
    nr <- ceiling(n / nc)
    par(mfrow = c(nr, nc), mar = c(1, 1, 3, 1))
    for (i in seq_len(n)) draw_gray(ours$mat[i, ], main = paste0(ours$df$name[i], "\n", ours$df$category[i]))
  }, height = function() {
    n <- max(nrow(ours$df), 1)
    200 * ceiling(n / min(n, 6))
  })

  our_idx <- reactive({
    req(input$our_pick)
    i <- match(input$our_pick, ours$df$file)
    req(!is.na(i))
    i
  })
  our_res <- reactive({ find_top5(ours$mat[our_idx(), ], input$K) })

  output$our_pred <- renderText({
    pred   <- predict_cat(our_res())
    actual <- ours$df$category[our_idx()]
    verdict <- if (is.na(actual)) "Actual category unknown (check the file name)"
               else if (pred == actual) "Right" else "Wrong"
    paste0("Predicted category: ", pred, " | Actual category: ", actual, " -> ", verdict)
  })

  output$our_best_plot <- renderPlot({
    i <- our_idx(); r <- our_res()
    par(mfrow = c(1, 2), mar = c(2, 1, 3, 1))
    draw_gray(ours$mat[i, ], main = paste0(ours$df$name[i], " - ", ours$df$category[i]))
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
    req(nrow(ours$df) > 0)
    K <- input$K
    pred <- vapply(seq_len(nrow(ours$df)), function(i) predict_cat(find_top5(ours$mat[i, ], K)), character(1))
    data.frame(Name = ours$df$name, File = ours$df$file, Actual = ours$df$category, Predicted = pred,
               Result = ifelse(is.na(ours$df$category), "Unknown", ifelse(pred == ours$df$category, "Right", "Wrong")),
               stringsAsFactors = FALSE)
  })

  output$our_summary <- renderText({
    d <- our_all()
    known <- d$Result != "Unknown"
    sprintf("K = %d: %d of %d shoes predicted correctly.", input$K, sum(d$Result == "Right"), sum(known))
  })

  output$our_all_table <- renderTable({ our_all() })
}


shinyApp(ui, server)
