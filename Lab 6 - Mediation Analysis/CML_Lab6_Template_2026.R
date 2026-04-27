#############################################
###  JOBS II experiment
#############################################

# ── Packages ──────────────────────────────────────────────────────────────────
library(mediation)   # install.packages("mediation")
library(sl3)         # remotes::install_github("tlverse/sl3")
library(medoutcon)   # remotes::install_github("nhejazi/medoutcon")
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DiagrammeR)

# ── Reproducibility ───────────────────────────────────────────────────────────
set.seed(123)

# ── Load data ─────────────────────────────────────────────────────────────────
data(jobs)
jobs <- as.data.frame(jobs)

# ── Initial inspection ────────────────────────────────────────────────────────
glimpse(jobs)


#############################################
###  VISUALIZATION
#############################################

p1 <- jobs %>%
  mutate(treat = factor(treat, labels = c("Control", "Treated"))) %>%
  ggplot(aes(x = treat, fill = treat)) +
  geom_bar(alpha = 0.8, width = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Count", title = "Treatment assignment") +
  theme(legend.position = "none")

p2 <- ggplot(jobs, aes(x = job_seek)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "steelblue", alpha = 0.7) +
  geom_density(color = "navy", linewidth = 0.8) +
  theme_minimal() +
  labs(x = "Job search self-efficacy", y = "Density",
       title = "Mediator: job_seek")

p1 + p2


jobs %>%
  mutate(treat = factor(treat, labels = c("Control", "Treated"))) %>%
  ggplot(aes(x = job_seek, fill = treat, color = treat)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = "Job search self-efficacy", y = "Density",
       fill = "Treatment", color = "Treatment",
       title = "Distribution of job_seek by treatment arm")


p3 <- jobs %>%
  mutate(
    work1 = factor(as.numeric(work1) - 1, labels = c("Unemployed", "Employed")),
    treat = factor(treat, labels = c("Control", "Treated"))
  ) %>%
  count(treat, work1) %>%
  group_by(treat) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = treat, y = prop, fill = work1)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_dodge(width = 0.5), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal() +
  labs(x = NULL, y = "Proportion", fill = "work1",
       title = "Employment status by treatment")

p4 <- jobs %>%
  mutate(treat = factor(treat, labels = c("Control", "Treated"))) %>%
  ggplot(aes(x = depress2, fill = treat, color = treat)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = "Depressive symptoms (depress2)", y = "Density",
       fill = "Treatment", color = "Treatment",
       title = "Depression score by treatment arm")

p3 + p4


#############################################
###  Question 1 — DAG
#############################################

# Draw the DAG manually (check the code in the Solutions)



#############################################
###  Question 2 — NDE and NIE: Effect of Training on Employment
#############################################

# ── Binary Super Learner ──────────────────────────────────────────────────────
# List available binary learners in sl3
sl3_list_learners("binomial")

# Build individual learners for a binary outcome
lrn_glm_bin    <- Lrnr_glm_fast$new(family = binomial())
lrn_gam_bin    <- Lrnr_gam$new(family = binomial())
lrn_earth_bin  <- Lrnr_earth$new(degree = 2)
lrn_ranger_bin <- Lrnr_ranger$new(probability = TRUE)

# Stack them and wrap in a Super Learner
stack_bin <- Stack$new(lrn_glm_bin, lrn_gam_bin, lrn_earth_bin, lrn_ranger_bin)
sl_bin    <- Lrnr_sl$new(learners = stack_bin)

# ── Continuous Super Learner ──────────────────────────────────────────────────
sl3_list_learners("continuous")

lrn_glm_cont    <- Lrnr_glm_fast$new(family = gaussian())
lrn_gam_cont    <- Lrnr_gam$new(family = gaussian())
lrn_earth_cont  <- Lrnr_earth$new(degree = 2)
lrn_ranger_cont <- Lrnr_ranger$new(probability = FALSE)

stack_cont <- Stack$new(lrn_glm_cont, lrn_gam_cont, lrn_earth_cont, lrn_ranger_cont)
sl_cont    <- Lrnr_sl$new(learners = stack_cont)

# ── DAG for Question 2 ────────────────────────────────────────────────────────
# For this question: X = treat, M = job_seek, Z = work1, C = baseline covariates.
# There is no outcome Y with a separate node here — Z plays the role of outcome.
# Draw the relevant DAG.



# ── Define variables ──────────────────────────────────────────────────────────
W <- jobs[, c("econ_hard", "depress1", "sex", "age", "occp",
              "marital", "nonwhite", "educ", "income")]
A <- ...
Y <- ...                # Binary outcome: employed at follow-up (use as.numeric(VAR) -1)
M <- ...                # Mediator: job search self-efficacy

folds <- 5

# ── NDE ───────────────────────────────────────────────────────────────────────
# Estimate the Natural Direct Effect (NDE) of A on Y, blocking the path through M.
# Hint: use medoutcon() with effect = "direct" and Z = NULL (no intermediate confounder here).

nde1 <- medoutcon(
  W = ...,
  A = ...,
  Z = ...,
  M = ...,
  Y = ...,
  g_learners = ...,
  h_learners = ...,
  b_learners = ...,
  effect     = "direct",
  estimator  = "onestep",
  estimator_args = list(cv_folds = folds)
)
nde1

# ── NIE ───────────────────────────────────────────────────────────────────────
# Estimate the Natural Indirect Effect (NIE) of A on Y, operating through M.

nie1 <- medoutcon(
  W = ...,
  A = ...,
  Z = ...,
  M = ...,
  Y = ...,
  g_learners = ...,
  h_learners = ...,
  b_learners = ...,
  effect     = "indirect",
  estimator  = "onestep",
  estimator_args = list(cv_folds = folds)
)
nie1

# Interpret your results.
# Does the job training programme affect employment mainly directly, or mainly
# through improving job search self-efficacy?


#############################################
###  Question 3 — NDE and NIE: Effect of Training on Depression (two mediators)
#############################################

# Now the outcome is depress2, and there are two mediators:
#   M1 = job_seek (continuous)
#   M2 = work1 (binary, employment status)
# Based on the previous results, is there an arrow from M1 to M2?



# ── Define variables ──────────────────────────────────────────────────────────
Y  <- ...
M1 <- ...
M2 <- ...

# ── NDE ───────────────────────────────────────────────────────────────────────
# Estimate the NDE of training on depression, blocking both mediators jointly.
# Hint: pass M = cbind(M1, M2). Which learner family is appropriate for b_learners here?

nde2 <- medoutcon(
  W = W, A = A, Z = NULL,
  M = ...,
  Y = Y,
  g_learners = ...,
  h_learners = ...,
  b_learners = ...,
  effect     = "direct",
  estimator  = "onestep",
  estimator_args = list(cv_folds = folds)
)
nde2

# ── NIE ───────────────────────────────────────────────────────────────────────
# Estimate the NIE of training on depression, flowing through both M1 and M2.

nie2 <- medoutcon(
  W = W, A = A, Z = NULL,
  M = ...,
  Y = Y,
  g_learners = ...,
  h_learners = ...,
  b_learners = ...,
  effect     = "indirect",
  estimator  = "onestep",
  estimator_args = list(cv_folds = folds)
)
nie2

# What does the combined NIE via {M1, M2} represent?
# How does it differ from the NIE you obtained in Question 2?


#############################################
###  Questions 4 & 5 — Path-Specific Effect via Mediator Difference
#############################################

# ── NIE through M1 only (job search self-efficacy) ────────────────────────────
# Estimate the NIE flowing through M1 alone (treating M2 as part of the direct path).

nie_MY <- medoutcon(
  W = W, A = A, Z = NULL,
  M = ...,
  Y = Y,
  g_learners = ...,
  h_learners = ...,
  b_learners = ...,
  effect     = "indirect",
  estimator  = "onestep",
  estimator_args = list(cv_folds = folds)
)
nie_MY

# ── Path-specific NIE through M2 (employment status) ─────────────────────────
N <- length(A)

# Point estimate: difference of NIEs
# Hint: nie_diff = NIE(M1 + M2) - NIE(M1 only)
nie_diff <- ...

# SE via the difference of efficient influence functions (EIFs)
# Hint: the EIF of a difference is the difference of the EIFs
nie_diff_se <- ...

# 95% confidence interval
nie_diff_lb <- round(nie_diff - 1.96 * nie_diff_se, 3)
nie_diff_ub <- round(nie_diff + 1.96 * nie_diff_se, 3)

cat("Path-specific NIE (M2 only):", nie_diff, "\n")
cat("SE:", nie_diff_se, "\n")
cat("95% CI: [", nie_diff_lb, ",", nie_diff_ub, "]\n")

# Is the path-specific NIE via M2 statistically significant?


#############################################
###  Summary of Results
#############################################

results <- tibble(
  Effect = c(
    "NDE (A → work1, blocking job_seek)",
    "NIE (A → work1, via job_seek)",
    "NDE (A → depress2, blocking M1+M2)",
    "NIE (A → depress2, via M1+M2)",
    "NIE (A → depress2, via M1 only)",
    "Path-specific NIE (A → depress2, via M2 only)"
  ),
  Estimate = c(
    nde1$theta,
    nie1$theta,
    nde2$theta,
    nie2$theta,
    nie_MY$theta,
    nie_diff          
  ),
  SE = c(
    sqrt(nde1$var),
    sqrt(nie1$var),
    sqrt(nde2$var),
    sqrt(nie2$var),
    sqrt(nie_MY$var),
    nie_diff_se      
  )
) %>%
  mutate(
    Z       = round(Estimate / SE, 3),
    P_Value = round(2 * pnorm(abs(Z), lower.tail = FALSE), 3),
    Estimate = round(Estimate, 3),
    SE       = round(SE, 3)
  )

print(results)

# ── Forest plot ───────────────────────────────────────────────────────────────
results %>%
  mutate(
    CI_lower = Estimate - 1.96 * SE,
    CI_upper = Estimate + 1.96 * SE,
    Effect   = factor(Effect, levels = rev(Effect))
  ) %>%
  ggplot(aes(x = Estimate, y = Effect)) +
  geom_point(size = 3, color = "tomato") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.25, linewidth = 0.8, color = "tomato") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(x = "Estimate (one-step, Super Learner)",
       y = NULL,
       title = "Natural Direct and Indirect Effects — JOBS II")

# Which effects are statistically significant?
# What is your overall interpretation of the causal mediation story in JOBS II?
# How do the two outcomes (work1 vs. depress2) paint different pictures of the
# programme's mechanisms?