# ============================================================
# Bootstrap uncertainty for seasonal NOK (Mouillot et al. 2005)
#
# 方法：对每个标准层的原始观测值有放回重采样
#       每次重采样后重新计算层均值，传入nicheoverlap核心函数
#       重复N_BOOT次，取SD和95% CI
#
# 输入：AOA Group A 和 Prochlorococcus 原始数据
# 输出：四个季节的 NOK point estimate + bootstrap SD + 95% CI
# ============================================================

library(readxl)

# ── 文件路径（修改为你本地路径）──────────────────────────────
PATH_AOA <- "C:/Paper and funding/2024 Autotrophy niche/Submission versions/Text/submission records/20251016 Revised for NC comments/Evolution data compilation/HOT_AmoA_gene_abundances_1692386884939.xlsx"
PATH_HOT <- "C:/Paper and funding/2024 Autotrophy niche/Submission versions/Text/submission records/20251016 Revised for NC comments/Evolution data compilation/HOT_2013-2019 processed.xlsx"

N_BOOT <- 1000   # bootstrap次数，可改为2000提高精度
set.seed(42)

# ── 共同标准层（去掉200m：Pro在200m无观测）───────────────────
STD_DEPTHS <- c(5, 25, 45, 75, 100, 125, 150, 175)

# ════════════════════════════════════════════════════════════
# 辅助函数
# ════════════════════════════════════════════════════════════
get_season <- function(month) {
  if (month %in% c(12, 1, 2)) return("Winter")
  if (month %in% c(3, 4, 5))  return("Spring")
  if (month %in% c(6, 7, 8))  return("Summer")
  return("Fall")
}

snap_to_std <- function(press,
                        std = c(5, 25, 45, 75, 100, 125, 150, 175, 200),
                        tol = 15) {
  press <- as.numeric(press)
  if (is.na(press)) return(NA)
  dists <- abs(std - press)
  idx   <- which.min(dists)
  if (dists[idx] <= tol) std[idx] else NA
}

# ════════════════════════════════════════════════════════════
# nicheoverlap 核心函数
# 完全复现 Mouillot et al. 2005 的 KDE 逻辑：
#   1. 相对丰度 × 1000 取整 → 展开成样本向量
#   2. Silverman bandwidth → gaussian KDE（1024点）
#   3. spline 积分校正 → 两曲线最小值积分 = NOK
# ════════════════════════════════════════════════════════════
nicheoverlap_core <- function(trab) {
  tr  <- trab[, 1]
  abN <- as.matrix(trab[, -1, drop = FALSE])
  N   <- ncol(abN)
  if (N < 2) return(NA)
  
  nichedens <- function(tr1, ab1) {
    ab1    <- as.numeric(ab1)
    abrel  <- ab1 / sum(ab1, na.rm = TRUE)
    abrelW <- round(abrel, 3) * 1000
    trrep  <- c()
    for (i in seq_along(abrelW)) {
      if (!is.na(abrelW[i]) && abrelW[i] > 0)
        trrep <- c(trrep, rep(tr1[i], as.integer(abrelW[i])))
    }
    if (length(trrep) < 2) return(NULL)
    bw  <- 1.06 * sd(trrep) * (length(trrep)^(-0.2))
    res <- density(trrep, bw = bw, kernel = "gaussian",
                   from = min(tr1), to = max(tr1), n = 1024)
    cbind(res$x, res$y)
  }
  
  dens_list <- vector("list", N)
  for (i in 1:N) {
    d <- nichedens(tr, abN[, i])
    if (is.null(d)) return(NA)
    dens_list[[i]] <- d
  }
  
  x     <- dens_list[[1]][, 1]
  densN <- sapply(dens_list, function(d) d[, 2])
  
  # 积分校正：使每条曲线在 [min(tr), max(tr)] 上积分 = 1
  densNcorr <- densN
  for (i in 1:N) {
    f    <- splinefun(x, densN[, i])
    intg <- tryCatch(
      integrate(f, min(x), max(x), subdivisions = 1000)$value,
      error = function(e) NA
    )
    if (is.na(intg) || intg <= 0) return(NA)
    densNcorr[, i] <- densN[, i] / intg
  }
  
  # NOK = integral of min(f1, f2)
  minN  <- apply(densNcorr, 1, min)
  fmin  <- splinefun(x, minN)
  overlap <- tryCatch(
    integrate(fmin, min(x), max(x), subdivisions = 1000)$value,
    error = function(e) NA
  )
  round(overlap, 4)
}

# ════════════════════════════════════════════════════════════
# 1. 加载 AOA Group A 原始数据
# ════════════════════════════════════════════════════════════
cat("加载 AOA 数据...\n")
df_aoa <- read_excel(PATH_AOA, sheet = "data")
df_aoa$month    <- as.integer(format(as.POSIXct(df_aoa$time), "%m"))
df_aoa$season   <- sapply(df_aoa$month, get_season)
df_aoa$depth_std <- sapply(as.numeric(df_aoa$depth), snap_to_std)
df_aoa$AOA_A    <- ifelse(as.integer(df_aoa$AmoA_GrpA_flag) == 0,
                          as.numeric(df_aoa$AmoA_GrpA_copies), NA)
df_aoa <- df_aoa[!is.na(df_aoa$depth_std) & df_aoa$depth_std <= 200, ]

# ════════════════════════════════════════════════════════════
# 2. 加载 Prochlorococcus 原始数据
# ════════════════════════════════════════════════════════════
cat("加载 Pro 数据...\n")
df_raw <- read_excel(PATH_HOT, sheet = "Sheet1", col_names = FALSE)
colnames(df_raw) <- c("date","time","press","Nitrate","LLN",
                      "HetBact","Pro","Syn","Euk")
df_raw <- df_raw[3:nrow(df_raw), ]   # 跳过列名行和单位行

df_raw$press <- as.numeric(df_raw$press)
df_raw$Pro   <- as.numeric(df_raw$Pro)
df_raw$date  <- as.numeric(df_raw$date)
df_raw$date_str <- formatC(as.integer(df_raw$date), width = 6, flag = "0")
df_raw$month    <- as.integer(substr(df_raw$date_str, 1, 2))
df_raw$season   <- sapply(df_raw$month, get_season)
df_raw$depth    <- df_raw$press
df_raw$Pro[!is.na(df_raw$Pro) & df_raw$Pro == -9] <- NA   # 缺失值标记

df_hot <- df_raw[!is.na(df_raw$depth) &
                   df_raw$depth > 0 &
                   df_raw$depth <= 200, ]

# ════════════════════════════════════════════════════════════
# 3. 逐季节：point estimate + bootstrap
# ════════════════════════════════════════════════════════════
SEASONS <- c("Winter", "Spring", "Summer", "Fall")
results <- data.frame(
  Season = character(),
  NOK_pct = numeric(),
  SD_pct  = numeric(),
  CI_lo_pct = numeric(),
  CI_hi_pct = numeric(),
  n_boot_valid = integer(),
  stringsAsFactors = FALSE
)

for (s in SEASONS) {
  cat("\n══ 季节:", s, "══\n")
  
  # 各层原始观测值
  sub_aoa <- df_aoa[df_aoa$season == s &
                      df_aoa$depth_std %in% STD_DEPTHS, ]
  sub_pro <- df_hot[df_hot$season == s, ]
  
  aoa_obs <- lapply(STD_DEPTHS, function(d) {
    v <- sub_aoa$AOA_A[sub_aoa$depth_std == d]
    v[!is.na(v)]
  })
  pro_obs <- lapply(STD_DEPTHS, function(d) {
    v <- sub_pro$Pro[sub_pro$depth >= (d - 10) & sub_pro$depth <= (d + 10)]
    v[!is.na(v)]
  })
  
  # Point estimate（季节均值）
  aoa_mean <- sapply(aoa_obs, function(x) if (length(x) > 0) mean(x) else NA)
  pro_mean <- sapply(pro_obs, function(x) if (length(x) > 0) mean(x) else NA)
  valid    <- !is.na(aoa_mean) & !is.na(pro_mean)
  
  trab_pt <- cbind(STD_DEPTHS[valid], aoa_mean[valid], pro_mean[valid])
  nok_pt  <- nicheoverlap_core(trab_pt)
  cat("  Point estimate NOK =", round(nok_pt * 100, 1), "%\n")
  
  # Bootstrap
  cat("  Running bootstrap (n =", N_BOOT, ")...\n")
  nok_boot <- numeric(N_BOOT)
  
  for (b in 1:N_BOOT) {
    aoa_b <- sapply(aoa_obs, function(x) {
      if (length(x) == 0) NA
      else mean(sample(x, length(x), replace = TRUE))
    })
    pro_b <- sapply(pro_obs, function(x) {
      if (length(x) == 0) NA
      else mean(sample(x, length(x), replace = TRUE))
    })
    valid_b <- !is.na(aoa_b) & !is.na(pro_b)
    if (sum(valid_b) < 3) { nok_boot[b] <- NA; next }
    
    trab_b      <- cbind(STD_DEPTHS[valid_b], aoa_b[valid_b], pro_b[valid_b])
    nok_boot[b] <- tryCatch(nicheoverlap_core(trab_b), error = function(e) NA)
    
    if (b %% 200 == 0) cat("   ", b, "/", N_BOOT, "\n")
  }
  
  nok_valid <- nok_boot[!is.na(nok_boot)]
  nok_sd    <- sd(nok_valid)
  ci_lo     <- quantile(nok_valid, 0.025)
  ci_hi     <- quantile(nok_valid, 0.975)
  
  cat("  Bootstrap SD =", round(nok_sd * 100, 1), "%\n")
  cat("  95% CI = [", round(ci_lo * 100, 1), "%,",
      round(ci_hi * 100, 1), "%]\n")
  cat("  有效bootstrap次数:", length(nok_valid), "/", N_BOOT, "\n")
  
  results <- rbind(results, data.frame(
    Season       = s,
    NOK_pct      = round(nok_pt  * 100, 1),
    SD_pct       = round(nok_sd  * 100, 1),
    CI_lo_pct    = round(ci_lo   * 100, 1),
    CI_hi_pct    = round(ci_hi   * 100, 1),
    n_boot_valid = length(nok_valid)
  ))
}

# ════════════════════════════════════════════════════════════
# 4. 输出结果
# ════════════════════════════════════════════════════════════
cat("\n\n══════════════════════════════════════════════════\n")
cat("季节性 NOK 结果汇总 (Group A vs Prochlorococcus)\n")
cat("══════════════════════════════════════════════════\n")
print(results)

write.csv(results, "NOK_seasonal_bootstrap_results.csv", row.names = FALSE)
cat("\n结果已保存至 NOK_seasonal_bootstrap_results.csv\n")