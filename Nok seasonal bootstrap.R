# ============================================================
# Bootstrap uncertainty for seasonal NOK 
# Niche overlap calculation is based on Mouillot et al. 2005
#
# Method：Resample the original observations for each standard layer with replacement.
#         Recalculate the layer mean after each resampling iteration and pass it to the nicheoverlap core function.
#         Repeat for N_BOOT iterations to calculate the Standard Deviation (SD) and 95% Confidence Interval (CI).
#
# Input：AOA Group A 和 Prochlorococcus Original observed data
# Output：NOK point estimate + bootstrap SD + 95% CI in diferent seasons
# ============================================================

library(readxl)

# ── File path（change to your own path）──────────────────────────────
PATH_AOA <- "your path/HOT_AmoA_gene_abundances_1692386884939.xlsx"
PATH_HOT <- "your path/HOT_2013-2019 processed.xlsx"

N_BOOT <- 1000   # Increase bootstrap replicates to 2000 to improve statistical precision
set.seed(42)

# ── Snap to standard depths (delete 200m, no obs data of Pro in this layer) ───────────────────
STD_DEPTHS <- c(5, 25, 45, 75, 100, 125, 150, 175)

# ════════════════════════════════════════════════════════════
# Helper Functions
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
# nicheoverlap Core Function
# reproduces the Kernel Density Estimation (KDE) logic described by Mouillot et al. (2005):
#   1. Data Expansion: Transforms relative abundance into a sample vector by multiplying by 1000 and rounding to the nearest integer.
#   2. Kernel Density Estimation: Applies a Gaussian KDE using Silverman’s Rule of Thumb for bandwidth selection (evaluated at 1024 points).
#   3. Integration & Overlap (NOK): Performs Spline Interpolation for integration correction. 
#      The final niche overlap value NOK is calculated as the area under the intersection of the two density curves.
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
  
  # Integral Normalization: Ensures that each probability density curve integrates to 1 over the interval [min(tr), max(tr)].
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
# 1. Upload AOA Group A
# ════════════════════════════════════════════════════════════
cat("uploading AOA ...\n")
df_aoa <- read_excel(PATH_AOA, sheet = "data")
df_aoa$month    <- as.integer(format(as.POSIXct(df_aoa$time), "%m"))
df_aoa$season   <- sapply(df_aoa$month, get_season)
df_aoa$depth_std <- sapply(as.numeric(df_aoa$depth), snap_to_std)
df_aoa$AOA_A    <- ifelse(as.integer(df_aoa$AmoA_GrpA_flag) == 0,
                          as.numeric(df_aoa$AmoA_GrpA_copies), NA)
df_aoa <- df_aoa[!is.na(df_aoa$depth_std) & df_aoa$depth_std <= 200, ]

# ════════════════════════════════════════════════════════════
# 2. Upload Prochlorococcus
# ════════════════════════════════════════════════════════════
cat("uploading Pro ...\n")
df_raw <- read_excel(PATH_HOT, sheet = "Sheet1", col_names = FALSE)
colnames(df_raw) <- c("date","time","press","Nitrate","LLN",
                      "HetBact","Pro","Syn","Euk")
df_raw <- df_raw[3:nrow(df_raw), ]   # Skip the header row (column names) and the units row

df_raw$press <- as.numeric(df_raw$press)
df_raw$Pro   <- as.numeric(df_raw$Pro)
df_raw$date  <- as.numeric(df_raw$date)
df_raw$date_str <- formatC(as.integer(df_raw$date), width = 6, flag = "0")
df_raw$month    <- as.integer(substr(df_raw$date_str, 1, 2))
df_raw$season   <- sapply(df_raw$month, get_season)
df_raw$depth    <- df_raw$press
df_raw$Pro[!is.na(df_raw$Pro) & df_raw$Pro == -9] <- NA   # mark the missing value

df_hot <- df_raw[!is.na(df_raw$depth) &
                   df_raw$depth > 0 &
                   df_raw$depth <= 200, ]

# ════════════════════════════════════════════════════════════
# 3. Each season：point estimate + bootstrap
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
  cat("\n══ Season:", s, "══\n")
  
  # Original obs value in each layer
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
  
  # Point estimate（seasonal mean）
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
  cat("  effective times of bootstrap:", length(nok_valid), "/", N_BOOT, "\n")
  
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
# 4. Output
# ════════════════════════════════════════════════════════════
cat("\n\n══════════════════════════════════════════════════\n")
cat("Summary of seasonal NOK (Group A vs Prochlorococcus)\n")
cat("══════════════════════════════════════════════════\n")
print(results)

write.csv(results, "NOK_seasonal_bootstrap_results.csv", row.names = FALSE)
cat("\nresults are saved in NOK_seasonal_bootstrap_results.csv\n")
