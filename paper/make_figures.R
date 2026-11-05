# =====================================================================
# Figures of the report that are not produced by forecast.R
# Run from the repository folder, after forecast.R:
#   source("paper/make_figures.R")
# Reads data/ and output/, writes paper/fig/
# =====================================================================

library(dplyr)
library(ggplot2)

dir.create("paper/fig", showWarnings = FALSE)
blue <- "#36648B"; red <- "#B22222"
th <- theme_minimal(base_size = 12) + theme(legend.position = "top")
save_fig <- function(name, g, width, height) {
  ggsave(file.path("paper/fig", name), g, width = width, height = height, dpi = 250, bg = "white")
}

# Figure 1: the series with its sample mean
d <- read.csv("data/hicp_data.csv") %>% mutate(date = as.Date(date)) %>% filter(!is.na(yoy))
m <- mean(d$yoy)
g <- ggplot(d, aes(date, yoy)) +
  geom_hline(yintercept = m, linetype = "dashed", colour = red) +
  geom_line(colour = blue, linewidth = 0.5) +
  annotate("segment", x = as.Date("2015-06-01"), xend = as.Date("2015-06-01"), y = m + 1.6, yend = m + 0.1,
           colour = red, arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("text", x = as.Date("2015-06-01"), y = m + 1.9, label = sprintf("Sample mean = %.2f%%", m),
           colour = red, size = 3) +
  labs(x = NULL, y = "Percent, annual rate") +
  theme_minimal(base_size = 10)
save_fig("fig_series.png", g, 8, 3.8)

# Figure 4: BIC and Box-Pierce statistic of the 30 candidates
sel <- read.csv("output/model_selection.csv") %>%
  mutate(Seasonals = ifelse(brk, "2011 break", "Fixed seasonals"),
         Residuals = ifelse(passes, "passes Box-Pierce", "fails Box-Pierce"))
cols <- c("2011 break" = red, "Fixed seasonals" = blue)
g <- ggplot(sel, aes(p, BIC, colour = Seasonals)) +
  geom_line() + geom_point(aes(shape = Residuals), size = 2.2) +
  scale_shape_manual(values = c("fails Box-Pierce" = 1, "passes Box-Pierce" = 16), name = NULL) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = 0:14) +
  labs(title = "(a) BIC", x = "Number of lags p", y = "BIC (lower is better)") +
  th + theme(legend.box = "vertical")
save_fig("fig_bic.png", g, 5.2, 3.9)

g <- ggplot(sel, aes(p, BP_Q, colour = Seasonals)) +
  geom_hline(yintercept = sel$BP_crit[1], linetype = "dashed") +
  annotate("text", x = 9.5, y = sel$BP_crit[1] * 1.45,
           label = sprintf("5%% critical value = %.1f", sel$BP_crit[1]), size = 3.5) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = 0:14) +
  scale_y_log10(breaks = c(10, 20, 30, 50, 100, 200)) +
  labs(title = sprintf("(b) Box-Pierce Q(%d)", floor(sqrt(sel$T[1]))),
       x = "Number of lags p", y = "Box-Pierce Q (log scale)") +
  th
save_fig("fig_bp.png", g, 5.2, 3.9)

# Figure 5: autocorrelations of the residuals of the selected model
res  <- read.csv("output/residuals_final_model.csv")$residual
band <- 1.96 / sqrt(length(res))
ac   <- data.frame(lag = 1:36, r = acf(res, lag.max = 36, plot = FALSE)$acf[-1])
g <- ggplot(ac, aes(lag, r)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-band, band), linetype = "dashed", colour = "steelblue") +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.8, colour = "grey30") +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  labs(x = "Lag (months)", y = "Residual autocorrelation") +
  theme_minimal(base_size = 10)
save_fig("fig_resid.png", g, 8, 3)

# Figure 6 (left): out-of-sample RMSPE by horizon
e <- read.csv("output/oos_errors.csv")
r <- e %>% filter(model != "Sample mean") %>%
  group_by(model, h) %>% summarise(RMSPE = sqrt(mean(error^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(model = recode(model, "AR(p) on annual rate" = "AR(p), annual rate",
                        "ARMA(1,1) annual rate" = "ARMA(1,1), annual rate"))
g <- ggplot(r, aes(h, RMSPE, colour = model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("Seasonal AR (selection rule)" = red, "AR(p), annual rate" = "darkorange",
                                 "ARMA(1,1), annual rate" = "darkgreen", "Random walk" = "grey30"), name = NULL) +
  scale_x_continuous(breaks = 1:12) +
  labs(x = "Horizon h (months)", y = "RMSPE (percentage points)") +
  th + guides(colour = guide_legend(nrow = 2))
save_fig("fig_rmspe.png", g, 5.2, 3.9)

# Figure 6 (right): model chosen by the selection rule at each forecast origin
ch <- read.csv("output/oos_rule_choices.csv") %>%
  mutate(origin = as.Date(origin),
         Seasonals = ifelse(grepl("break", choice), "2011 break", "Fixed seasonals"),
         p = as.numeric(sub(".*p=", "", choice)))
g <- ggplot(ch, aes(origin, p, colour = Seasonals)) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_y_continuous(breaks = seq(0, 14, 2), limits = c(0, 14)) +
  labs(x = "Forecast origin", y = "Selected p") +
  th
save_fig("fig_rule.png", g, 5.2, 3.9)

cat("Figures saved in paper/fig/\n")
