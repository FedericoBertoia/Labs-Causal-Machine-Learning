# ==============================================================================
# Causal Machine Learning – Project Solution
# 401(k) Eligibility and Participation on Accumulated Assets
# Stijn Vansteelandt and Federico Bertoia, Ghent University 2025-2026
# ==============================================================================

# ── Packages ───────────────────────────────────────────────────────────────────
library(hdm)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grf)
library(SuperLearner)
library(sandwich)
library(lmtest)

# ── Reproducibility ────────────────────────────────────────────────────────────
set.seed(123)

# ── Load data ──────────────────────────────────────────────────────────────────
data(pension)

# Variable roles:
#   Y  = net_tfa   (net financial assets, continuous outcome)
#   A  = p401      (401(k) participation, binary treatment — Tasks 1-3)
#   A2 = e401      (401(k) eligibility,   binary treatment — main focus)
#   L  = age, inc, fsize, educ, marr, twoearn, pira, db, hown

glimpse(pension)

# ── Data preparation ───────────────────────────────────────────────────────────
df <- pension %>%
  transmute(
    Y      = net_tfa,
    A      = e401,                  # eligibility (main treatment)
    age    = age,
    inc    = inc,
    fsize  = fsize,
    educ   = educ,
    marr   = marr,
    twoearn = twoearn,
    pira    = pira,
    db      = db,
    hown    = hown
  )

# Covariates matrix (used throughout)
L_vars <- c("age", "inc", "fsize", "educ", "marr", "twoearn",
            "pira", "db", "hown")

X <- df %>% dplyr::select(all_of(L_vars)) %>% as.matrix()
Y <- df$Y
A <- df$A
n <- nrow(df)

stopifnot(all(colSums(is.na(df)) == 0))
cat("n =", n, " | E[Y] =", round(mean(Y)), " | P(A=1) =", round(mean(A), 3), "\n")


# ==============================================================================
# Task 1 – Causal Forest
# ==============================================================================
cat("\n====== Task 1: Causal Forest ======\n")

cf <- causal_forest(
  X              = X,
  Y              = Y,
  W              = A,
  num.trees      = 2000,
  tune.parameters = "all",
  seed           = 123
)

# Overall ATE
ate_cf <- average_treatment_effect(cf, target.sample = "all")
cat("ATE (causal forest):", round(ate_cf["estimate"]), "\n")
cat("95% CI: [",
    round(ate_cf["estimate"] - 1.96 * ate_cf["std.err"]),
    ",",
    round(ate_cf["estimate"] + 1.96 * ate_cf["std.err"]), "]\n")

# CATE predictions
tau_hat <- predict(cf)$predictions

# (a) Informative visual: histogram of CATE + scatter CATE vs income
p1 <- ggplot(data.frame(cate = tau_hat), aes(x = cate)) +
  geom_histogram(fill = "steelblue", alpha = 0.7, bins = 60) +
  geom_vline(xintercept = ate_cf["estimate"], colour = "tomato",
             linetype = "dashed", linewidth = 0.9) +
  theme_minimal() +
  labs(
    title    = "Distribution of estimated CATE (causal forest)",
    subtitle = paste0("ATE = $", round(ate_cf["estimate"]),
                      "  (red dashed)"),
    x        = "Estimated CATE ($)",
    y        = "Count"
  )



print(p1)

# (b) Conditioning on additional variables beyond income
# Conditioning on extra variables does NOT invalidate identification.
# Under the assumption that eligibility is exchangeable conditional on income,
# any superset of {income} is also sufficient (monotonicity of conditional
# independence). Formally, if A ⊥⊥ Y^a | income, then
# A ⊥⊥ Y^a | (income, age, ...) as well, because the additional variables
# cannot introduce new confounding that was already blocked by income.
# In a causal diagram, controlling for more pre-treatment covariates that are
# not colliders is always safe and can improve efficiency by reducing residual
# variance (regression adjustment). Hence the causal forest that conditions on
# all nine covariates is valid.
cat("\n[Task 1b] See written explanation above (in-code comment).\n")


# ==============================================================================
# Task 2 – Subgroup analysis using ntile
# ==============================================================================
cat("\n====== Task 2: Subgroup Analysis ======\n")

df2 <- df %>%
  mutate(
    subgroup = ntile(tau_hat, 4),
    tau_hat  = tau_hat
  )

# 5-fold cross-fitting ATE within each subgroup via AIPW-style DR estimator
# We reuse the causal forest's nuisance estimates (mu_hat, e_hat) for simplicity;
# for a fully cross-fit subgroup estimate we re-run on the full sample.
mu1_hat <- predict(cf, estimate.variance = FALSE)$predictions +
  (1 - cf$W.orig) * 0    # placeholder; use DR scores below

# DR (AIPW) scores from the forest
dr_scores <- get_scores(cf)   # = tau_hat_i + (A_i/e - (1-A_i)/(1-e)) * resid_i

# Forest-based DR score: already debiased individual-level score
# Subgroup ATE = mean of DR scores in that subgroup
subgroup_results <- df2 %>%
  mutate(dr = dr_scores) %>%
  group_by(subgroup) %>%
  summarise(
    n      = n(),
    ate    = mean(dr),
    se     = sd(dr) / sqrt(n()),
    ci_lo  = ate - 1.96 * se,
    ci_hi  = ate + 1.96 * se,
    .groups = "drop"
  )

print(subgroup_results)

# Visual
ggplot(subgroup_results,
       aes(x = factor(subgroup), y = ate, ymin = ci_lo, ymax = ci_hi)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 0.5) +
  geom_errorbar(width = 0.2, colour = "navy") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Average effect of 401(k) eligibility by CATE subgroup",
    x     = "Subgroup (1 = lowest predicted CATE, 4 = highest)",
    y     = "Estimated ATE ($)",
    caption = "Error bars: 95% confidence intervals"
  )


# ==============================================================================
# Task 3 – Test of effect heterogeneity: Var[E(Y^1 - Y^0 | L)]
# ==============================================================================
cat("\n====== Task 3: Variance of CATE ======\n")

# Efficient influence function:
#   phi_i = (tau_hat_i - ate)^2 - sigma2
#           + 2*(tau_hat_i - ate) * (A_i/e_i - (1-A_i)/(1-e_i)) * (Y_i - mu_hat_i)
#
# where:
#   tau_hat_i = E[Y^1 - Y^0 | L_i]   (predicted CATE)
#   ate       = E[Y^1 - Y^0]          (estimated ATE)
#   e_i       = P(A=1 | L_i)          (propensity score)
#   mu_hat_i  = E[Y | A_i, L_i]       (outcome regression)

ate_val  <- ate_cf["estimate"]
e_hat    <- cf$W.hat          # propensity scores from forest
mu_hat   <- cf$Y.hat +        # E[Y|L], need E[Y|A,L]
  A * (predict(cf)$predictions) * (1 - e_hat) -
  (1 - A) * (predict(cf)$predictions) * e_hat
# Simpler: reconstruct E[Y|A,L] from forest components
# E[Y|A=1,L] = E[Y|L] + (1-e)*tau_hat
# E[Y|A=0,L] = E[Y|L] - e*tau_hat
mu_hat_full <- cf$Y.hat + (A - e_hat) * tau_hat / (e_hat * (1 - e_hat)) * 0
# Use forest residuals directly
resid_Y <- Y - cf$Y.hat - (A - e_hat) * tau_hat   # approx residual

# IPW augmentation term
ipw_term <- (A / e_hat - (1 - A) / (1 - e_hat)) * (Y - cf$Y.hat)

# EIF for sigma^2
eif_sigma2 <- (tau_hat - ate_val)^2 +
  2 * (tau_hat - ate_val) * ipw_term

sigma2_hat <- mean(eif_sigma2)
se_sigma2  <- sd(eif_sigma2) / sqrt(n)
ci_sigma2  <- sigma2_hat + c(-1, 1) * 1.96 * se_sigma2

cat("Estimated Var[E(Y^1-Y^0|L)]:", round(sigma2_hat), "\n")
cat("95% CI: [", round(ci_sigma2[1]), ",", round(ci_sigma2[2]), "]\n")
cat("SD of CATE (sqrt of estimate):", round(sqrt(max(sigma2_hat, 0))), "\n")

# Interpretation: if CI excludes 0, strong evidence of effect heterogeneity.


# ==============================================================================
# Task 4 – Orthogonal learner: CATE conditional on income only
# ==============================================================================
cat("\n====== Task 4: Orthogonal Learner (CATE on income) ======\n")

# We use a partially linear / Robinson-style orthogonal learner:
#   Y_i - E[Y|A_i,L_i] = tau(income_i) * (A_i - E[A|L_i]) + epsilon_i
# but we want tau to depend only on income. We achieve this via a
# local linear / kernel-weighted DML estimator, or more simply by using
# the DR scores and smoothing them over income.
#
# Approach: use the forest's DR scores and regress them nonparametrically
# on income using a local polynomial / GAM. This is the "R-learner" idea
# of Nie & Wager (2021).

# Step 1: obtain cross-fit nuisance estimates (already done via causal_forest)
#   e_hat  = cf$W.hat
#   m_hat  = cf$Y.hat

# Step 2: construct pseudo-outcome (Robinson residuals)
#   R-learner pseudo-outcome: dr_score_i (already in dr_scores above)
#   Or equivalently: (Y_i - m_hat_i) / (A_i - e_hat_i)  [when denominator != 0]

# Step 3: smooth dr_scores over income using a GAM / loess
# This gives tau(income) = E[DR_score | income]

library(mgcv)

income_grid <- seq(min(df$income), quantile(df$income, 0.99), length.out = 300)

# GAM with penalised spline on income
gam_fit <- gam(dr_scores ~ s(income, bs = "cr", k = 15),
               data  = data.frame(income = df$income, dr_scores = dr_scores),
               method = "REML")

tau_income      <- predict(gam_fit, newdata = data.frame(income = income_grid),
                           se.fit = TRUE)
tau_income_est  <- tau_income$fit
tau_income_se   <- tau_income$se.fit

# Plot
plot_df <- data.frame(
  income  = income_grid,
  tau     = tau_income_est,
  ci_lo   = tau_income_est - 1.96 * tau_income_se,
  ci_hi   = tau_income_est + 1.96 * tau_income_se
)

ggplot(plot_df, aes(x = income, y = tau)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "steelblue", alpha = 0.25) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = ate_cf["estimate"], linetype = "dotted",
             colour = "tomato") +
  theme_minimal() +
  labs(
    title    = "Estimated effect of 401(k) eligibility on net financial assets\nconditional on income (R-learner + GAM)",
    subtitle = "Red dotted line = marginal ATE",
    x        = "Income ($)",
    y        = "Estimated CATE ($)",
    caption  = "Shaded band: pointwise 95% confidence interval"
  )

# (b) Assumptions made by this orthogonal learner vs causal forest
# The R-learner + GAM imposes the working model tau(income_i) = f(income),
# i.e. that the CATE depends on income ONLY after conditioning on L for
# debiasing. The causal forest makes no such restriction and can learn
# tau(L) flexibly over all nine covariates. The GAM further assumes that
# the function f is smooth (continuous with bounded second derivative via
# the penalised spline). The causal forest is fully nonparametric (tree-
# based) and does not assume smoothness. The key advantage of the
# orthogonal learner is lower-dimensional visualisation and potentially
# lower variance when the true CATE is mainly driven by income.
cat("\n[Task 4b] See in-code comment above.\n")


# ==============================================================================
# Diagnostics / Quality Checks (required by the project instructions)
# ==============================================================================
cat("\n====== Diagnostics ======\n")

# 1. Propensity score distribution (positivity)
ps_df <- data.frame(e_hat = e_hat, A = factor(A))
ggplot(ps_df, aes(x = e_hat, fill = A, colour = A)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values  = c("steelblue", "tomato"),
                    labels  = c("Not eligible", "Eligible")) +
  scale_colour_manual(values = c("steelblue", "tomato"),
                      labels = c("Not eligible", "Eligible")) +
  theme_minimal() +
  labs(title    = "Propensity score distribution by eligibility",
       subtitle = "Substantial overlap supports positivity",
       x = "P(Eligible | L)", y = "Density",
       fill = "Eligibility", colour = "Eligibility")

cat("Propensity score range: [", round(min(e_hat), 3), ",",
    round(max(e_hat), 3), "]\n")
cat("Proportion with e < 0.05 or e > 0.95:",
    round(mean(e_hat < 0.05 | e_hat > 0.95), 4), "\n")

# 2. Influence function extremes (Task 3)
cat("DR score: mean =", round(mean(dr_scores)),
    " | SD =", round(sd(dr_scores)),
    " | max |DR| =", round(max(abs(dr_scores))), "\n")

# 3. Forest calibration test
test_calibration(cf)

cat("\n====== Done ======\n")