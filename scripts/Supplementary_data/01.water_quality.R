# ============================================================
# Load libraries
# ============================================================
library(tidyverse)
library(car)
library(emmeans)
library(patchwork)

# Set default ggplot theme
theme_set(theme_bw())

# ============================================================
# Load input data
# ============================================================
dat <- read.csv("data/environmental_data/water_quality.csv", check.names = FALSE, stringsAsFactors = FALSE)

if ("δ13C" %in% names(dat)) dat <- dat %>% rename(d13C = `δ13C`)
if ("δ15N" %in% names(dat)) dat <- dat %>% rename(d15N = `δ15N`)

dat <- dat %>%
  rename(
    Date = Date,
    Water_Depth = Water_Depth,
    site = Area
  ) %>%
  mutate(
    Date = factor(Date, levels = c("Day_0","Day_12","Day_30")),
    day_num = case_when(
      Date == "Day_0"  ~ 0,
      Date == "Day_12" ~ 12,
      Date == "Day_30" ~ 30,
      TRUE ~ NA_real_
    )
  )
  
# ============================================================
# Statistical analysis
# ============================================================
# Calculate time-weighted means using the trapezoidal rule
twa_mean <- function(x, t){
  ok <- !(is.na(x) | is.na(t))
  x <- x[ok]; t <- t[ok]
  if (length(x) <= 1) return(NA_real_)
  o <- order(t); x <- x[o]; t <- t[o]
  if (length(x) == 2) return(mean(x))
  dt <- diff(t)
  auc <- sum( (head(x,-1) + tail(x,-1)) / 2 * dt )
  auc / (max(t) - min(t))
}

vars <- c("DOC", "POC", "HIX", "BIX", "TDN", "TDP", "d13C", "d15N")

env_unit <- dat %>%
  group_by(site, Water_Depth) %>%
  summarise(
    across(all_of(vars),
           ~ twa_mean(.x, day_num), .names = "{.col}"),
    across(all_of(vars),
           ~ sum(!is.na(.x)), .names = "n_{.col}"),
    .groups = "drop"
  )

# Test differences among sites and water depths using linear models
run_anova <- function(response, data = env_unit){
  f <- reformulate(c("site","Water_Depth"), response = response)
  df <- data %>% select(all_of(c("site","Water_Depth", response))) %>% drop_na()
  if (nrow(df) < 5) return(cat(response, ": not enough data\n"))
  fit <- lm(f, data = df)
  cat("\n=== ", response, " ===\n")
  print(Anova(fit, type = 3))
  emm <- emmeans(fit, ~ site)
  print(pairs(emm, adjust = "tukey"))
  invisible(fit)
}

for(v in vars){
  run_anova(v)
}

# ============================================================
# Plot temporal changes in environmental variables
# ============================================================
plot_timeseries <- function(var){
    ylab_var <- switch(var,
    "d13C" = expression(delta^13*C~"["*"\u2030"*"]"),
    "d15N" = expression(delta^15*N~"["*"\u2030"*"]"),
    "TDN"  = expression(TDN~"["*mu*MN*"]"),
    "TDP"  = expression(TDP~"["*mu*MP*"]"),
    "DOC"  = expression(DOC~"["*mu*MC*"]"),
    "POC"  = expression(log[10]~"(POC)"~"["*mu*g~C~L^{-1}*"]"),
    var
  )
    
    p <- ggplot(dat, aes(
        x = day_num,
        y = .data[[var]],
        group = interaction(site, Water_Depth),
        color = Water_Depth,
        shape = as.factor(site)
    )) +
        geom_line(alpha = 0.7, linewidth = 0.8) +
        geom_point(size = 2.5) +
        scale_x_continuous(breaks = c(0, 12, 30)) +
        scale_color_manual(values = c(
            "hyporheic 20cm water" = "#E64B35FF",
            "hyporheic 5cm water"  = "#4DBBD5FF",
            "surface water"        = "#00A087FF"
        )) +
        labs(x = "Day", y = ylab_var, color = "Depth", shape = "Site") +
        theme_bw(base_size = 12) +
        theme(
            text = element_text(family = "Helvetica", size = 13),
            axis.title = element_text(size = 14),
            axis.text  = element_text(size = 12),
            panel.border = element_rect(color = "black"),
            aspect.ratio = 0.5,
            plot.margin = margin(5, 5, 5, 15)
        )
    
    if (var == "POC") {
        p <- p + scale_y_log10(breaks = c(100, 300, 1000, 3000, 10000))
    }
    
    return(p)
}

p1 <- plot_timeseries("d13C")
p2 <- plot_timeseries("d15N")
p3 <- plot_timeseries("HIX")
p4 <- plot_timeseries("BIX")
p5 <- plot_timeseries("POC")
p6 <- plot_timeseries("DOC")
p7 <- plot_timeseries("TDN")
p8 <- plot_timeseries("TDP")

blank <- patchwork::plot_spacer()

supp_fig <- (
    (p1 | p2 | p3) /
        (p4 | p5 | p6) /
        (p7 | p8 | blank)
) +
    plot_layout(
        guides = "collect",
        heights = c(1,1,1)
    ) &
    theme(
        legend.position = "right",
        legend.justification = c("left", "center")
    )

print(supp_fig)