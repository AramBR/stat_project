# ==================================================================
# Bootstrap standard errors and confidence intervals for the Faithful fit
# (resample the data with replacement, refit EM each time, and use the
#  spread of the refitted estimates; see gmm_ci() in R/helpers.R)
# ==================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
source("R/core.R")
source("R/helpers.R")

theme_set(theme_minimal(base_size = 12) + theme(legend.position = "bottom"))

# ------------------------------------------------------------------
# Old Faithful: estimates with bootstrap SEs and 95% CIs
# ------------------------------------------------------------------
set.seed(42)
x <- faithful$waiting
fit_f <- em_gmm(x, init_kmeans(x, 2))
ci_f <- gmm_ci(x, fit_f, B = 500)
kable(
    ci_f,
    digits = 3,
    caption = "Table 9. Bootstrap standard errors and 95% percentile CIs for the Faithful fit (B = 500)."
)

ggplot(ci_f %>% filter(parameter != "pi2"), aes(estimate, parameter)) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
    geom_point(size = 3, colour = "#1b9e77") +
    facet_wrap(~parameter, scales = "free", ncol = 1) +
    labs(
        x = "Estimate with 95% bootstrap CI", y = NULL,
        title = "Faithful: parameter estimates with bootstrap intervals"
    ) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
