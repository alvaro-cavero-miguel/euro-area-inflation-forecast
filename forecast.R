# =====================================================================
# Econ 460 - Forecast report
# Euro area HICP inflation (annual rate, %) - forecast of October 2026
# Author: Álvaro Cavero Miguel
#
# Steps:
#   1. Data (Eurostat, prc_hicp_minr): annual rate + price index
#   2. Exploration: ACF / PACF and the change in the seasonal pattern in 2011
#   3. Models and the model selection rule (Box-Pierce + BIC)
#   4. Out-of-sample comparison (RMSPE, recursive window, 2015 - present)
#   5. Final model: selection table, parameters, 12-month forecast, 50% / 90% intervals
#   6. Diagnostics of the final model (Box-Jenkins check)
# =====================================================================

library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------------
# 0. Settings
# ---------------------------------------------------------------------
DOWNLOAD       <- FALSE                   # FALSE = use data/hicp_data.csv (the data used in the report)
                                          # TRUE  = download the latest data from Eurostat
OOS_START      <- as.Date("2015-01-01")   # first out-of-sample forecast origin
H              <- 12                      # forecast horizon: 12 months
P_MAX_YOY      <- 13                      # max AR lags, annual rate (PACF spike at lag 13: base effect)
P_MAX_MM       <- 14                      # max AR lags, monthly change: more than a year of monthly
                                          # lags, so lags 12-13 can capture an evolving seasonal pattern
SEASONAL_BREAK <- as.Date("2011-01-01")   # EU rules on seasonal items (Reg. 330/2009) applied from 2011
ALPHA_BP       <- 0.05                    # Box-Pierce significance level (lecture 9)
dir.create("output", showWarnings = FALSE)
dir.create("data",   showWarnings = FALSE)

# ---------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------
if (DOWNLOAD) {
  library(eurostat)
  raw  <- get_eurostat("prc_hicp_minr", time_format = "date", cache = FALSE)
  hicp <- raw %>% filter(geo == "EA", coicop18 == "TOTAL")   # euro area (changing composition), all items

  data <- inner_join(
    hicp %>% filter(unit == "RCH_A") %>% select(date = TIME_PERIOD, yoy = values),
    hicp %>% filter(unit == "I25")   %>% select(date = TIME_PERIOD, index = values),
    by = "date") %>%
    arrange(date)
  write.csv(data, "data/hicp_data.csv", row.names = FALSE)
} else {
  data <- read.csv("data/hicp_data.csv") %>% mutate(date = as.Date(date))
}

data <- data %>%
  filter(!is.na(index)) %>%
  mutate(mm       = 100 * (log(index) - log(lag(index))),   # monthly change (log, %)
         yoy_calc = 100 * (index / lag(index, 12) - 1),     # annual rate computed from index
         month    = factor(format(date, "%m")),
         period   = ifelse(date >= SEASONAL_BREAK, "2011-2026", "1997-2010"))

# The official annual rate (rounded to 1 decimal) should match the computed one
cat("Max |official - computed| annual rate:",
    round(max(abs(data$yoy - data$yoy_calc), na.rm = TRUE), 2), "pp\n")

data  <- data %>% filter(!is.na(yoy), !is.na(mm))
T_end <- nrow(data)
cat("Sample:", format(min(data$date), "%Y-%m"), "to", format(max(data$date), "%Y-%m"),
    "|", T_end, "obs. | last value:", data$yoy[T_end], "%\n")

# ---------------------------------------------------------------------
# 2. Exploration
# ---------------------------------------------------------------------
# (a) ACF and PACF (lectures 5-9)
png("output/acf_pacf.png", width = 2400, height = 1500, res = 250)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
acf(data$yoy,  lag.max = 36, main = "ACF: annual inflation")
pacf(data$yoy, lag.max = 36, main = "PACF: annual inflation")
acf(data$mm,   lag.max = 36, main = "ACF: monthly change")
pacf(data$mm,  lag.max = 36, main = "PACF: monthly change")
invisible(dev.off())

# (b) Did the seasonal pattern change in 2011? (lecture 11: do not mix regimes)
season_tab <- data %>% group_by(month, period) %>% summarise(mean_mm = mean(mm), .groups = "drop")
g_break <- ggplot(season_tab, aes(month, mean_mm, fill = period)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("grey65", "steelblue4")) +
  labs(title = "Average month-on-month change by calendar month, before and after 2011",
       x = "Month", y = "Percent (log change)", fill = NULL,
       caption = "Source: Eurostat (prc_hicp_minr); author's calculations") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave("output/seasonality_by_period.png", g_break, width = 8, height = 4.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------
# 3. Models
#    Each forecasting function takes the data up to origin T and returns
#    the forecast of the annual inflation rate for T+1, ..., T+H
# ---------------------------------------------------------------------

# AR(p) estimated by OLS: y_t = c + phi_1 y_{t-1} + ... + phi_p y_{t-p} (+ X_t b) + e_t
fit_ar <- function(y, p, X = NULL) {
  n  <- length(y)
  Y  <- embed(y, p + 1)                      # col 1 = y_t, cols 2.. = lags
  df <- data.frame(y = Y[, 1])
  if (p > 0) df <- cbind(df, Y[, -1, drop = FALSE])
  if (!is.null(X)) df <- cbind(df, X[(p + 1):n, , drop = FALSE])
  lm(y ~ ., data = df)
}

# Information criteria as in lecture 10: ln(SSR/T) + penalty. Lower is better.
info_criteria <- function(fit) {
  T  <- nobs(fit)
  k  <- length(coef(fit))
  s2 <- sum(residuals(fit)^2) / T
  c(AIC = log(s2) + 2 * k / T, BIC = log(s2) + k * log(T) / T)
}

# Box-Pierce (lecture 9): Q = T * (r_1^2 + ... + r_R^2) ~ chi2(R), R = floor(sqrt(T)).
# The residuals pass (look like white noise) if Q is below the chi2 critical value.
box_pierce <- function(res) {
  R <- floor(sqrt(length(res)))
  r <- acf(res, lag.max = R, plot = FALSE)$acf[-1]          # [-1] drops lag 0
  c(Q = length(res) * sum(r^2), R = R, crit = qchisq(1 - ALPHA_BP, R))
}

# Choose the lag order p by BIC on a common estimation sample (lectures 10, 11a)
select_p <- function(y, p_max, X = NULL) {
  n <- length(y)
  bics <- sapply(0:p_max, function(p) {
    idx <- (p_max - p + 1):n                 # same effective sample for every p
    info_criteria(fit_ar(y[idx], p, if (is.null(X)) NULL else X[idx, , drop = FALSE]))[["BIC"]]
  })
  (0:p_max)[which.min(bics)]
}

# Chain rule of forecasting for an AR(p) with known future regressors (lecture 10)
forecast_ar <- function(fit, y, p, h, X_future = NULL) {
  b    <- coef(fit)
  path <- y
  out  <- numeric(h)
  for (j in 1:h) {
    x <- b[1]
    if (p > 0) x <- x + sum(b[2:(p + 1)] * rev(tail(path, p)))
    if (!is.null(X_future)) x <- x + sum(b[(p + 2):length(b)] * X_future[j, ])
    out[j] <- x
    path <- c(path, x)
  }
  out
}

# --- Benchmarks and models on the annual rate ---
f_rw     <- function(d, h) rep(tail(d$yoy, 1), h)                 # random walk
f_mean   <- function(d, h) rep(mean(d$yoy), h)                    # sample mean
f_ar_yoy <- function(d, h) {                                      # AR(p) by BIC
  p <- select_p(d$yoy, P_MAX_YOY)
  forecast_ar(fit_ar(d$yoy, p), d$yoy, p, h)
}
f_arma11 <- function(d, h) {                                      # ARMA(1,1)
  fit <- tryCatch(suppressWarnings(arima(d$yoy, order = c(1, 0, 1), method = "ML",
                                         optim.control = list(maxit = 1000))),
                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, h))
  as.numeric(predict(fit, n.ahead = h)$pred)
}

# --- Seasonal AR on the monthly change (lectures 4, 10, 12) ---
# m_t = c + 11 monthly dummies + phi_1 m_{t-1} + ... + phi_p m_{t-p} + e_t
#   brk = FALSE: the monthly dummies are the same for 1997-2026
#   brk = TRUE : the dummies (and the constant) may change from 2011 on (dummies x Post2011)
seasonal_X <- function(dates, brk) {
  m <- factor(format(dates, "%m"), levels = sprintf("%02d", 1:12))
  D <- model.matrix(~ m)[, -1, drop = FALSE]          # 11 dummies (January = base)
  colnames(D) <- month.abb[-1]
  if (!brk) return(D)
  post <- as.numeric(dates >= SEASONAL_BREAK)
  DP   <- D * post
  colnames(DP) <- paste0(month.abb[-1], "_post2011")
  cbind(D, Post2011 = post, DP)
}

# MODEL SELECTION RULE (fixed before looking at which model wins):
#   Candidates: fixed or 2011-break seasonals, p = 0, ..., P_MAX_MM (30 models),
#               all estimated on the same sample so that their BIC is comparable.
#   Step 1 (lecture 9, Box-Jenkins check): keep the models whose residuals pass
#           Box-Pierce. A model that fails leaves predictable information unused.
#   Step 2 (lecture 10): among them, choose the lowest BIC.
#   If no candidate passes (can happen in short early samples), lowest BIC overall.
seasonal_candidates <- function(d) {
  n <- nrow(d)
  do.call(rbind, lapply(c(FALSE, TRUE), function(brk) {
    X <- seasonal_X(d$date, brk)
    do.call(rbind, lapply(0:P_MAX_MM, function(p) {
      idx <- (P_MAX_MM - p + 1):n
      fit <- fit_ar(d$mm[idx], p, X[idx, , drop = FALSE])
      ic  <- info_criteria(fit)
      bp  <- box_pierce(residuals(fit))
      data.frame(seasonals = if (brk) "break 2011" else "fixed", brk = brk, p = p,
                 k = length(coef(fit)), T = nobs(fit),
                 AIC = ic[["AIC"]], BIC = ic[["BIC"]],
                 BP_Q = bp[["Q"]], BP_crit = bp[["crit"]], passes = bp[["Q"]] < bp[["crit"]])
    }))
  }))
}
select_seasonal <- function(d) {
  cand <- seasonal_candidates(d)
  pool <- if (any(cand$passes)) cand[cand$passes, ] else cand
  pool[which.min(pool$BIC), ]
}

fit_seasonal <- function(d, brk, p) {
  X <- seasonal_X(d$date, brk)
  list(brk = brk, p = p, fit = fit_ar(d$mm, p, X),
       labels = c("Constant", if (p > 0) paste0("m(t-", seq_len(p), ")"), colnames(X)))
}

# Forecast the monthly change, then rebuild the annual rate (lecture 12):
#   yoy_{T+h} = 100 * (P_{T+h} / P_{T+h-12} - 1), where P_{T+h-12} is already known for h <= 12
forecast_seasonal <- function(d, s, h) {
  future_dates <- seq(max(d$date), by = "month", length.out = h + 1)[-1]
  mm_future <- forecast_ar(s$fit, d$mm, s$p, h, seasonal_X(future_dates, s$brk))
  logP <- log(tail(d$index, 1)) + cumsum(mm_future) / 100
  P_12 <- d$index[nrow(d) - 12 + (1:h)]               # index level 12 months earlier
  100 * (exp(logP) / P_12 - 1)
}

# The rule is re-applied at every forecast origin, as it would have been in real time
f_seasonal_rule <- function(d, h) {
  sel <- select_seasonal(d)
  out <- forecast_seasonal(d, fit_seasonal(d, sel$brk, sel$p), h)
  attr(out, "choice") <- paste0(sel$seasonals, ", p=", sel$p)
  out
}

# The BIC can only compare models of the same variable (monthly change). The models of
# the annual rate are a different variable, so they are kept as benchmarks and compared
# with the seasonal rule out of sample (part 4).
FINAL_MODEL <- "Seasonal AR (selection rule)"
models <- list("Random walk"           = f_rw,
               "Sample mean"           = f_mean,
               "AR(p) on annual rate"  = f_ar_yoy,
               "ARMA(1,1) annual rate" = f_arma11,
               "Seasonal AR (selection rule)" = f_seasonal_rule)

# ---------------------------------------------------------------------
# 4. Out-of-sample comparison (recursive / expanding window; Handouts 2 and 4)
# ---------------------------------------------------------------------
origins <- which(data$date >= OOS_START & seq_len(T_end) < T_end)
cat("\nOut-of-sample evaluation:", length(origins), "forecast origins. This takes about a minute...\n")

choices <- character(0)
errors <- do.call(rbind, lapply(origins, function(t) {
  d <- data[1:t, ]
  do.call(rbind, lapply(names(models), function(m) {
    f  <- models[[m]](d, H)
    if (m == FINAL_MODEL) choices[as.character(data$date[t])] <<- attr(f, "choice")
    hh <- 1:H
    ok <- t + hh <= T_end
    if (!any(ok)) return(NULL)
    data.frame(model = m, origin = data$date[t], h = hh[ok],
               target = data$date[t + hh[ok]],
               forecast = as.numeric(f[ok]), actual = data$yoy[t + hh[ok]])
  }))
})) %>% mutate(error = actual - forecast)          # e = y - yhat
write.csv(errors, "output/oos_errors.csv", row.names = FALSE)

cat("\nModel chosen by the selection rule at the", length(choices), "origins:\n")
print(table(choices))
write.csv(data.frame(origin = names(choices), choice = choices), "output/oos_rule_choices.csv", row.names = FALSE)

high_inflation <- function(date) date >= as.Date("2021-07-01") & date <= as.Date("2023-12-01")

rmspe <- errors %>%
  group_by(model, h) %>%
  summarise(RMSPE          = sqrt(mean(error^2, na.rm = TRUE)),
            RMSPE_ex_21_23 = sqrt(mean(error[!high_inflation(target)]^2, na.rm = TRUE)),
            .groups = "drop")

rmspe_table <- data.frame(model = names(models))
for (hh in c(1, 3, 6, 12)) {
  rmspe_table[[paste0("h=", hh)]] <- sapply(rmspe_table$model, function(m)
    rmspe$RMSPE[rmspe$model == m & rmspe$h == hh])
}
rmspe_table$`h=1 ex 2021-23` <- sapply(rmspe_table$model, function(m)
  rmspe$RMSPE_ex_21_23[rmspe$model == m & rmspe$h == 1])
rmspe_table <- rmspe_table %>% arrange(`h=1`)
cat("\nOut-of-sample RMSPE (percentage points), sorted by h=1:\n")
print(rmspe_table %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
write.csv(rmspe_table, "output/oos_rmspe.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# 5. Final model: selection table, parameters, forecast and intervals
# ---------------------------------------------------------------------
# (a) Selection table on the full sample
candidates <- seasonal_candidates(data)
cat("\nSelection of the seasonal model (common sample, T =", candidates$T[1],
    "| Box-Pierce R =", floor(sqrt(candidates$T[1])), "| 5% critical value =",
    round(candidates$BP_crit[1], 1), "):\n")
print(candidates %>% select(seasonals, p, k, AIC, BIC, BP_Q, passes) %>%
        mutate(across(c(AIC, BIC), ~ round(.x, 4)), BP_Q = round(BP_Q, 1)), row.names = FALSE)
write.csv(candidates, "output/model_selection.csv", row.names = FALSE)

sel <- select_seasonal(data)
cat("\nSelected model: seasonal AR,", sel$seasonals, "seasonals, p =", sel$p,
    "| BIC =", round(sel$BIC, 4), "| Box-Pierce Q =", round(sel$BP_Q, 1), "\n")
if (sel$p == P_MAX_MM) warning("The selected p equals P_MAX_MM: consider increasing P_MAX_MM.")

# (b) Parameters of the selected model, estimated on the full sample
final_fit <- fit_seasonal(data, sel$brk, sel$p)
param_table <- data.frame(Parameter = final_fit$labels,
                          Estimate  = coef(final_fit$fit),
                          SE        = sqrt(diag(vcov(final_fit$fit))))
rownames(param_table) <- NULL
sigma_hat <- summary(final_fit$fit)$sigma
cat("\nParameters | p =", sel$p, "| sigma =", round(sigma_hat, 3),
    "| R2 =", round(summary(final_fit$fit)$r.squared, 3), "| n =", nobs(final_fit$fit), "\n")
print(param_table %>% mutate(across(Estimate:SE, ~ round(.x, 3))))
write.csv(param_table, "output/parameters_final_model.csv", row.names = FALSE)

# (c) Forecast for the next 12 months
f_final      <- forecast_seasonal(data, final_fit, H)
future_dates <- seq(max(data$date), by = "month", length.out = H + 1)[-1]

# (d) Intervals: forecast +/- z * v_h, z = 0.674 (50%), 1.645 (90%)
#   Main (Handout 4): v_h = standard deviation of the out-of-sample errors at horizon h.
#   Comparison (model based, as in the lecture slides v_h^2 = (1 + phi^2 + ...) sigma^2):
#   the annual-rate error at horizon h is the sum of the monthly errors from T+1 to T+h.
#   With psi_j the AR weights of the monthly model, its variance is
#   sigma^2 * sum_{r=1..h} (psi_0 + ... + psi_{r-1})^2, times (1 + forecast/100) to go
#   from the log index to the annual rate. It ignores parameter uncertainty.
v_oos <- errors %>% filter(model == FINAL_MODEL) %>%
  group_by(h) %>% summarise(v = sd(error, na.rm = TRUE)) %>% pull(v)

phi <- if (sel$p > 0) coef(final_fit$fit)[2:(sel$p + 1)] else numeric(0)
psi <- c(1, numeric(H - 1))
for (j in 2:H) {
  i <- seq_len(min(sel$p, j - 1))
  psi[j] <- sum(phi[i] * psi[j - i])
}
v_model <- sapply(1:H, function(h) sigma_hat * sqrt(sum(cumsum(psi)[1:h]^2))) * (1 + f_final / 100)

z50 <- qnorm(0.75); z90 <- qnorm(0.95)
forecast_table <- data.frame(
  month    = paste(month.abb[as.integer(format(future_dates, "%m"))], format(future_dates, "%Y")),
  h        = 1:H,
  forecast = f_final,
  lo90 = f_final - z90 * v_oos, lo50 = f_final - z50 * v_oos,
  hi50 = f_final + z50 * v_oos, hi90 = f_final + z90 * v_oos,
  v_oos = v_oos, v_model = v_model
)
cat("\nForecast of the annual rate (%), intervals from out-of-sample errors (Handout 4):\n")
print(forecast_table %>% mutate(across(forecast:v_model, ~ round(.x, 2))))
cat("\nOctober 2026 forecast:", round(f_final[1], 1), "% | 90% interval:",
    round(forecast_table$lo90[1], 1), "to", round(forecast_table$hi90[1], 1), "%",
    "\nModel-based 90% interval for comparison:",
    round(f_final[1] - z90 * v_model[1], 1), "to", round(f_final[1] + z90 * v_model[1], 1), "%\n")
write.csv(forecast_table, "output/forecast_table.csv", row.names = FALSE)

# (e) Chart: recent history + forecast with 50% / 90% bands
hist_plot <- data %>% filter(date >= as.Date("2019-01-01"))
fc_plot   <- forecast_table %>% mutate(date = future_dates)
g_fc <- ggplot() +
  geom_ribbon(data = fc_plot, aes(date, ymin = lo90, ymax = hi90), fill = "steelblue", alpha = 0.18) +
  geom_ribbon(data = fc_plot, aes(date, ymin = lo50, ymax = hi50), fill = "steelblue", alpha = 0.35) +
  geom_line(data = hist_plot, aes(date, yoy), color = "grey20", linewidth = 0.6) +
  geom_line(data = rbind(data.frame(date = max(data$date), forecast = tail(data$yoy, 1)),
                         fc_plot[, c("date", "forecast")]),
            aes(date, forecast), color = "steelblue4", linewidth = 0.8) +
  geom_point(data = fc_plot[1, ], aes(date, forecast), color = "firebrick", size = 2.5) +
  geom_hline(yintercept = mean(data$yoy), linetype = "dashed", color = "grey50") +
  labs(title = "Euro area HICP inflation: forecast with 50% and 90% intervals",
       subtitle = paste0("Seasonal AR(", sel$p, "), ", sel$seasonals,
                         " seasonals. Red point: October 2026 (forecast target)"),
       x = NULL, y = "Percent, annual rate",
       caption = "Source: Eurostat (prc_hicp_minr); author's forecasts. Dashed line: sample mean.") +
  theme_minimal(base_size = 12)
ggsave("output/forecast_fan.png", g_fc, width = 8, height = 4.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------
# 6. Diagnostics of the final model (Box-Jenkins step 3, lecture 9)
#    If the model captures the dynamics and the seasonality, its residuals are white
#    noise: no autocorrelation at lag 1, nor at lag 12 (same month a year earlier).
# ---------------------------------------------------------------------
res   <- residuals(final_fit$fit)
r_acf <- acf(res, lag.max = 36, plot = FALSE)$acf[-1]
bp    <- box_pierce(res)
cat("\nFinal model residuals | autocorrelation lag 1:", round(r_acf[1], 3),
    "| lag 12:", round(r_acf[12], 3), "| lag 24:", round(r_acf[24], 3),
    "| 95% band: +/-", round(1.96 / sqrt(length(res)), 3),
    "\n   Box-Pierce Q(", bp[["R"]], ") =", round(bp[["Q"]], 1),
    "| 5% critical value =", round(bp[["crit"]], 1), "\n")

# ---------------------------------------------------------------------
# 7. Robustness checks (reported in the paper; NOT used to choose the model)
#    (a) October 2026 forecast under alternative specifications
#    (b) Selection rule with Box-Pierce degrees of freedom R - p (Box and Pierce, 1970)
#    (c) Empirical coverage of the out-of-sample 50% / 90% intervals
# ---------------------------------------------------------------------
R_bp <- floor(sqrt(candidates$T[1]))
alt  <- data.frame(brk = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE), p = c(6, 12, 13, 14, 6, 13))
alt$oct_2026 <- sapply(seq_len(nrow(alt)), function(i)
  forecast_seasonal(data, fit_seasonal(data, alt$brk[i], alt$p[i]), 1))
alt <- merge(alt, candidates[, c("brk", "p", "seasonals", "BIC", "BP_Q")], by = c("brk", "p"))
alt$passes_dfR  <- alt$BP_Q < qchisq(1 - ALPHA_BP, R_bp)            # course convention
alt$passes_dfRp <- alt$BP_Q < qchisq(1 - ALPHA_BP, R_bp - alt$p)    # R - p degrees of freedom
alt <- alt[order(alt$brk, alt$p), ]
cat("\n(a) October 2026 forecast under alternative specifications:\n")
print(alt %>% select(seasonals, p, BIC, BP_Q, passes_dfR, passes_dfRp, oct_2026) %>%
        mutate(BIC = round(BIC, 3), BP_Q = round(BP_Q, 1), oct_2026 = round(oct_2026, 2)), row.names = FALSE)
write.csv(alt, "output/robustness_specifications.csv", row.names = FALSE)

cand_dfRp <- candidates %>% mutate(passes = BP_Q < qchisq(1 - ALPHA_BP, R_bp - p))
sel_dfRp  <- cand_dfRp %>% filter(passes) %>% slice_min(BIC, n = 1)
cat("\n(b) Rule with R - p degrees of freedom selects:", sel_dfRp$seasonals, "seasonals, p =", sel_dfRp$p,
    "| October 2026 forecast:",
    round(forecast_seasonal(data, fit_seasonal(data, sel_dfRp$brk, sel_dfRp$p), 1), 2), "\n")

coverage <- errors %>% filter(model == FINAL_MODEL) %>%
  group_by(h) %>% mutate(v = sd(error, na.rm = TRUE)) %>% ungroup() %>%
  mutate(in50 = abs(error) <= z50 * v, in90 = abs(error) <= z90 * v)
cov_tab <- rbind(
  coverage %>% filter(h == 1) %>% summarise(sample = "h = 1", in50 = mean(in50), in90 = mean(in90), n = n()),
  coverage %>% summarise(sample = "h = 1 to 12", in50 = mean(in50), in90 = mean(in90), n = n()),
  coverage %>% filter(h == 1, !high_inflation(target)) %>%
    summarise(sample = "h = 1, excl. 2021-23", in50 = mean(in50), in90 = mean(in90), n = n()))
cat("\n(c) Share of out-of-sample errors inside the intervals (nominal 50% and 90%):\n")
print(cov_tab %>% mutate(across(c(in50, in90), ~ round(.x, 2))), row.names = FALSE)
write.csv(cov_tab, "output/interval_coverage.csv", row.names = FALSE)

cat("\nDone. Results saved in the 'output/' folder.\n")
