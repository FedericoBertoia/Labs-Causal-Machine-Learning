#############################################
##  Instrumental Variable Analysis
##  IGFBP-3 and Dense Breast Tissue (n = 139)
#############################################

# ── Packages ───────────────────────────────────────────────────────────────────
library(foreign)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(SuperLearner)
library(lmtest)
library(sandwich)
library(DiagrammeR)

set.seed(123)


#############################################
##  Data Loading & Preparation
#############################################

data <- read.dta("igf.dta")
N    <- nrow(data)

# ── Log-transform exposure and outcome ────────────────────────────────────────
data <- data %>%
  mutate(
    log_lucent = log(lucent),
    log_igfbp3 = log(igfbp3)
  )

glimpse(data)

# ── Analysis variables ────────────────────────────────────────────────────────
Y <- data$log_lucent    # outcome
A <- data$log_igfbp3   # exposure
R <- data$igfbp3g_num  # instrument (numeric: AA=1, AC=2, CC=3)

# Numeric coding of genotype
data <- data %>%
  mutate(igfbp3g_num = as.numeric(igfbp3g))

R <- data$igfbp3g_num

attach(data)

################################################################################
##  Question 1 — Association of Lucent Area with Genotype
################################################################################

# ── Visualise ────────────────────────────────────────────────────────────────
ggplot(data, aes(x = igfbp3g, y = log_lucent, fill = igfbp3g)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.45) +
  scale_fill_manual(values = c("tomato", "steelblue", "forestgreen")) +
  theme_minimal() +
  labs(x = "Genotype (igfbp3g)", y = "log(lucent area)",
       title = "Log(Lucent Area) by IGFBP-3 Genotype") +
  theme(legend.position = "none")



# ── Linear trend ──────────────────────────────────────────────────────────────
lm_q1 <- lm(log_lucent ~ igfbp3g_num, data = data)
summary(lm_q1)

# ── Factor model (no ordering assumed) ───────────────────────────────────────
lm_q1_cat <- lm(log_lucent ~ factor(igfbp3g), data = data)
summary(lm_q1_cat)

# ── Pearson correlation ───────────────────────────────────────────────────────
cor.test(data$log_lucent, data$igfbp3g_num)

# > Interpretation: The genotype is significantly associated with
# > log(lucent area) (p = 0.018). Each additional C allele corresponds to
# > a decrease of ~0.154 in log-lucent area. Under the IV assumptions,
# > this association propagates entirely through IGFBP-3 levels.


################################################################################
##  Question 2 — Distribution of IGFBP-3 and Relevance of the Instrument
################################################################################

# ── 2.1: Histograms (raw and log-scale) ──────────────────────────────────────
p1 <- ggplot(data, aes(x = igfbp3)) +
  geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7, color = "white") +
  theme_minimal() +
  labs(x = "igfbp3", y = "Frequency", title = "Histogram of igfbp3")

p2 <- ggplot(data, aes(x = log_igfbp3)) +
  geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7, color = "white") +
  theme_minimal() +
  labs(x = "log_igfbp3", y = "Frequency", title = "Histogram of log(igfbp3)")

p1 + p2

# ── 2.2: Distribution by genotype ─────────────────────────────────────────────
data %>%
  group_by(igfbp3g) %>%
  summarise(
    n    = n(),
    mean = round(mean(log_igfbp3), 3),
    sd   = round(sd(log_igfbp3),   3),
    min  = round(min(log_igfbp3),  3),
    max  = round(max(log_igfbp3),  3)
  )

ggplot(data, aes(x = log_igfbp3, fill = igfbp3g, colour = igfbp3g)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 20) +
  scale_fill_manual(values   = c("tomato", "steelblue", "forestgreen")) +
  scale_colour_manual(values = c("tomato", "steelblue", "forestgreen")) +
  theme_minimal() +
  labs(x = "log_igfbp3", y = "Count",
       title = "Distribution of log(IGFBP-3) by Genotype",
       fill = "igfbp3g", colour = "igfbp3g")

# ── 2.3: Testing relevance ────────────────────────────────────────────────────
lm_q2 <- lm(log_igfbp3 ~ igfbp3g_num, data = data)
summary(lm_q2)

lm_q2_cat <- lm(log_igfbp3 ~ factor(igfbp3g), data = data)
anova(lm_q2_cat)

# > Interpretation: The genotype is significantly associated with log(IGFBP-3)
# > (p = 0.009), supporting the relevance assumption. Each additional C allele
# > is associated with a decrease of ~0.063 in log(IGFBP-3). This is the only
# > IV assumption directly testable from data.


################################################################################
##  Question 3 — DAG and IV Assumptions
################################################################################

grViz("
digraph IV_DAG {
  graph [layout = neato, bgcolor = transparent, splines = curved]
  node [shape = plaintext, fontname = 'Helvetica', fontsize = 16]
  R [pos = '0,0!',   label = 'R\n(igfbp3g)']
  A [pos = '3,0!',   label = 'A\n(log igfbp3)']
  Y [pos = '6,0!',   label = 'Y\n(log lucent)']
  U [pos = '4.5,2!', label = 'U']
  L [pos = '3,-2!',  label = 'L\n(age, afbcat,\nbrfcat, parcat, ev_pil)']

  edge [arrowhead = vee, arrowsize = 0.8, penwidth = 1.6, color = '#444444']
  R -> A
  A -> Y
  U -> A
  U -> Y
  L -> A
  L -> Y

  edge [arrowhead = vee, arrowsize = 0.7, penwidth = 1.2,
        style = dashed, color = '#888888']
  R -> Y [style = dashed, constraint = false]
  R -> L [style = dashed]
}
")

# The three IV assumptions:
#
# 1. Relevance (R not-indep A | L):
#    The IGFBP3g genotype must be associated with serum IGFBP-3 after
#    adjusting for L. Verified empirically in Question 2 (p = 0.009).
#
# 2. Exclusion restriction (R indep Y | A, L, U):
#    The genotype should not affect breast tissue density except through
#    IGFBP-3. Potentially violated if the genotype influences parity or age
#    at first birth, which in turn affect breast density.
#
# 3. Randomization (R indep U | L):
#    Mendelian inheritance randomizes allele transmission, so the genotype
#    should be independent of unmeasured confounders after conditioning on L.
#    Not testable from data; relies on biological reasoning.


################################################################################
##  Question 4 — IV Estimate via EIF (No Covariates)
##
##  psi_hat = sum_i (R_i - R_bar)(Y_i - Y_bar) /
##             sum_i (R_i - R_bar)(A_i - A_bar)
################################################################################

# ── Point estimate ────────────────────────────────────────────────────────────
psi <- mean((R - mean(R)) * (Y - mean(Y))) /
       mean((R - mean(R)) * (A - mean(A)))

# ── Standard error via influence function ─────────────────────────────────────
EIC_psi <- (R - mean(R)) * (Y - mean(Y) - psi * (A - mean(A))) /
           mean((R - mean(R)) * (A - mean(A)))
se_psi <- sd(EIC_psi) / sqrt(N)

# ── 95% CI ────────────────────────────────────────────────────────────────────
psi_lb <- psi - 1.96 * se_psi
psi_ub <- psi + 1.96 * se_psi

# ── Results on geometric mean ratio scale (10% increase in A) ─────────────────
results_q4 <- tibble(
  Method   = "EIF (no covariates)",
  Estimate = round(exp(psi    * log(1.1)), 3),
  CI_Lower = round(exp(psi_lb * log(1.1)), 3),
  CI_Upper = round(exp(psi_ub * log(1.1)), 3)
)

knitr::kable(results_q4, digits = 3,
             caption = "Q4 — IV estimate (no covariates): GMR for 10% increase in IGFBP-3")

# > Interpretation: A 10% increase in IGFBP-3 increases the geometric mean of
# > dense breast tissue area by ~26.1% [95% CI: 0.1% to 58.9%]. The CI is wide,
# > reflecting the small sample and modest instrument strength.


################################################################################
##  Question 5 — IV Estimate via EIF (2 Covariates)
##
##  Confounder set: afbcat + parcat (blocks the R -> L pathway)
##  No cross-fitting given n = 139.
################################################################################

L_5 <- data %>% dplyr::select(afbcat, parcat)

SL.library <- c("SL.glm", "SL.ranger")

# ── Nuisance models ───────────────────────────────────────────────────────────
sl_Y <- SuperLearner(Y = Y, X = L_5, SL.library = SL.library)
sl_A <- SuperLearner(Y = A, X = L_5, SL.library = SL.library)
sl_R <- SuperLearner(Y = R, X = L_5, SL.library = SL.library)

# ── Nuisance predictions ──────────────────────────────────────────────────────
m_Y <- predict(sl_Y, newdata = L_5, onlySL = TRUE)$pred
m_A <- predict(sl_A, newdata = L_5, onlySL = TRUE)$pred
m_R <- predict(sl_R, newdata = L_5, onlySL = TRUE)$pred

# ── EIF estimator ─────────────────────────────────────────────────────────────
psi_L <- mean((R - m_R) * (Y - m_Y)) /
         mean((R - m_R) * (A - m_A))

# ── Standard error ────────────────────────────────────────────────────────────
EIC_psi_L <- (R - m_R) * (Y - m_Y - psi_L * (A - m_A)) /
             mean((R - m_R) * (A - m_A))
se_psi_L  <- sd(EIC_psi_L) / sqrt(N)

# ── 95% CI ────────────────────────────────────────────────────────────────────
psi_L_lb <- psi_L - 1.96 * se_psi_L
psi_L_ub <- psi_L + 1.96 * se_psi_L

# ── Results ───────────────────────────────────────────────────────────────────
results_q5 <- tibble(
  Method   = "DML-EIF (afbcat + parcat)",
  Estimate = round(exp(psi_L    * log(1.1)), 3),
  CI_Lower = round(exp(psi_L_lb * log(1.1)), 3),
  CI_Upper = round(exp(psi_L_ub * log(1.1)), 3)
)

knitr::kable(results_q5, digits = 3,
             caption = "Q5 — DML IV estimate with covariates: GMR for 10% increase in IGFBP-3")

# > Interpretation: After conditioning on afbcat and parcat to block the
# > genotype -> confounders pathway, a 10% increase in IGFBP-3 increases
# > the geometric mean of lucent area by ~24.9% [95% CI: 3.4% to 50.9%].
# > The interval no longer includes zero.


################################################################################
##  Question 6 — Assumption-Lean Regression (Full Covariates)
##
##  Model: E[Y | A, L] = beta * A + g(L)  (partially linear)
##  Full confounder set L = (age, afbcat, brfcat, parcat, ev_pil)
##
##  Two approaches:
##   (a) EIF estimator
##   (b) OLS residual-on-residual with sandwich SE
################################################################################

L2 <- data %>% dplyr::select(age, afbcat, brfcat, parcat, ev_pil)

# ── Nuisance models ───────────────────────────────────────────────────────────
sl_Y2 <- SuperLearner(Y = Y, X = L2, SL.library = SL.library)
sl_A2 <- SuperLearner(Y = A, X = L2, SL.library = SL.library)

m_Y2 <- predict(sl_Y2, newdata = L2, onlySL = TRUE)$pred
m_A2 <- predict(sl_A2, newdata = L2, onlySL = TRUE)$pred

# ── (a) EIF estimator ─────────────────────────────────────────────────────────
beta_eif <- mean((A - m_A2) * (Y - m_Y2)) /
            mean((A - m_A2)^2)

EIC_beta <- (A - m_A2) * (Y - m_Y2 - beta_eif * (A - m_A2)) /
            mean((A - m_A2)^2)
se_eif   <- sd(EIC_beta) / sqrt(N)

beta_eif_lb <- beta_eif - 1.96 * se_eif
beta_eif_ub <- beta_eif + 1.96 * se_eif

# ── (b) OLS residual-on-residual with sandwich SE ─────────────────────────────
lm_resid  <- lm(I(Y - m_Y2) ~ -1 + I(A - m_A2))
lm_result <- coeftest(lm_resid, vcov = sandwich)

beta_lm    <- lm_result[1, 1]
se_lm      <- lm_result[1, 2]
beta_lm_lb <- beta_lm - 1.96 * se_lm
beta_lm_ub <- beta_lm + 1.96 * se_lm

# ── Results ───────────────────────────────────────────────────────────────────
results_q6 <- tibble(
  Method   = c("EIF", "LM (sandwich)"),
  Estimate = round(exp(c(beta_eif,    beta_lm)    * log(1.1)), 3),
  CI_Lower = round(exp(c(beta_eif_lb, beta_lm_lb) * log(1.1)), 3),
  CI_Upper = round(exp(c(beta_eif_ub, beta_lm_ub) * log(1.1)), 3)
)

knitr::kable(results_q6, digits = 3,
             caption = "Q6 — Assumption-lean regression: GMR for 10% increase in IGFBP-3")

# > Interpretation: The assumption-lean regression estimates a much smaller
# > association — ~4.5% increase [95% CI: -0.2% to 9.3%], not significant.
# > Compared to the IV estimates (~25%), this suggests negative unmeasured
# > confounding that masks the causal effect in the observational analysis.


#############################################
##  Summary & Forest Plot
#############################################

results_all <- tibble(
  Method   = c(
    "IV — EIF (no Covariates)",
    "IV — EIF with 2 Covariates",
    "IV — EIF with Full Covariates"
  ),
  Estimate = exp(c(psi, psi_L, beta_eif, beta_lm) * log(1.1)),
  Lower    = exp(c(psi_lb, psi_L_lb, beta_eif_lb, beta_lm_lb) * log(1.1)),
  Upper    = exp(c(psi_ub, psi_L_ub, beta_eif_ub, beta_lm_ub) * log(1.1))
)

knitr::kable(
  results_all %>% mutate(across(where(is.numeric), ~round(.x, 3))),
  caption = "Summary — IV and assumption-lean regression results"
)

# ── Forest plot ───────────────────────────────────────────────────────────────
results_all %>%
  mutate(
    Method = factor(Method, levels = rev(Method)),
    Type   = c("IV", "IV", "Regression", "Regression")
  ) %>%
  ggplot(aes(x = Estimate, y = Method, color = Type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1, scale = 100)) +
  scale_color_manual(values = c("IV" = "tomato", "Regression" = "steelblue")) +
  theme_minimal() +
  labs(
    x     = "Geometric mean ratio (10% increase in IGFBP-3)",
    y     = NULL,
    color = NULL,
    title = "IV and Regression Estimates — IGFBP-3 Effect on Breast Density"
  ) +
  theme(legend.position = "bottom")
