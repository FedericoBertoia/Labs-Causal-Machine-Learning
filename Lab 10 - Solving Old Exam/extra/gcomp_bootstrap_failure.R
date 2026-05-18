# =============================================================================
# G-computation ATE (dose-response slope) estimation:
#   bootstrap failure with adaptive ML nuisances, and DML as the fix.
#
# Directly mirrors Vansteelandt (CML Welcome Lecture) slides 41-44, now
# expressed as a g-computation problem.
#
# Setup (identical DGP to bootstrap_failure.R):
#   E(Y | A, X) = 0   (true psi = 0, true beta = 0)
#   E(A | X)    = 3X  (strong confounding)
#   X ~ N(0,1), A = 3X + N(0,1), Y = N(0,1)
#
# Target estimand: ATE = d/dA E[Y | do(A)] = psi.
# G-computation with a linear outcome model: E[Y|A,X] = psi*A + beta*X.
# For a linear model without interactions, g-computation gives the same
# estimate as the OLS coefficient on A.  The POINT is the post-selection
# mechanism, not the nonlinearity.
#
# Adaptive outcome model procedure (same as slides 41-44):
#   1. Fit full model Y ~ A + X; extract p-value for beta.
#   2. If p < alpha -> g-computation with full model -> ATE = psi_hat.
#   3. Else         -> g-computation with Y ~ A only -> ATE = psi0_hat.
#
# Why bootstrap fails:
#   In the real world beta = 0, so P(select full | truth) = alpha = 5%.
#   Oracle SE reflects mostly psi0_hat (low SE, good use of Var(A)=10).
#   In the bootstrap world beta_hat != 0, so the full model is selected ~16%.
#   Bootstrap distribution reflects more psi_hat (high SE, uses only
#   Var(A|X)=1), inflating bootstrap SE -> OVER-coverage (>> 95%).
#
# Fix: Double Machine Learning (DML) with sample splitting.
#   Split data into halves.  Fit nuisance on fold 1, evaluate on fold 2 (and
#   vice versa).  The cross-fitted residuals are not contaminated by the same
#   selection event, so the bootstrap works.  (Here we use OLS nuisances to
#   isolate the cross-fitting benefit from the selection issue.)
# =============================================================================

set.seed(123)

n     <- 500   # sample size  (same as bootstrap_failure.R)
nsim  <- 200   # simulation replications
B     <- 50    # bootstrap replicates per simulation
alpha <- 0.05  # significance level for the outcome-model selection test
tau   <- 0     # true ATE (dose-response slope)

# ---------------------------------------------------------------------------
# Post-selection g-computation (same mechanism as slides 41-44)
# ---------------------------------------------------------------------------
gcomp_adaptive <- function(Y, A, X, alpha = 0.05) {
  fit_full <- lm(Y ~ A + X)
  pval_X   <- summary(fit_full)$coefficients["X", "Pr(>|t|)"]
  if (pval_X < alpha) {
    # G-comp from full model: (1/n) sum [psi*1 - psi*0] = psi
    coef(fit_full)["A"]
  } else {
    coef(lm(Y ~ A))["A"]
  }
}

# ---------------------------------------------------------------------------
# DML (cross-fitted) g-computation / partial linear regression
#   Split data in two folds.  Residualise Y and A on X using the OTHER fold.
#   ATE = cov(Y_resid, A_resid) / var(A_resid)   (Robinson 1988).
# ---------------------------------------------------------------------------
dml_ate <- function(Y, A, X) {
  n   <- length(Y)
  idx <- sample(n)              # random split
  f1  <- idx[seq_len(n %/% 2)]
  f2  <- idx[(n %/% 2 + 1):n]

  # Fold-1 fits applied to fold-2 data (and vice versa)
  fit_Y1 <- lm(Y[f1] ~ X[f1])
  fit_A1 <- lm(A[f1] ~ X[f1])
  Yr2 <- Y[f2] - predict(fit_Y1, data.frame(X = X[f2]))  # residuals fold 2
  Ar2 <- A[f2] - predict(fit_A1, data.frame(X = X[f2]))

  fit_Y2 <- lm(Y[f2] ~ X[f2])
  fit_A2 <- lm(A[f2] ~ X[f2])
  Yr1 <- Y[f1] - predict(fit_Y2, data.frame(X = X[f1]))
  Ar1 <- A[f1] - predict(fit_A2, data.frame(X = X[f1]))

  Yr_all <- c(Yr1, Yr2)
  Ar_all <- c(Ar1, Ar2)
  coef(lm(Yr_all ~ Ar_all - 1))   # DML estimator (no intercept, residuals)
}

# ---------------------------------------------------------------------------
# Simulation loop
# ---------------------------------------------------------------------------
ate_ada  <- numeric(nsim)
ate_dml  <- numeric(nsim)
se_ada   <- numeric(nsim)
se_dml   <- numeric(nsim)
cover_ada  <- logical(nsim)
cover_dml  <- logical(nsim)
boot_frac  <- numeric(nsim)

cat("Running simulation...\n")
for (s in seq_len(nsim)) {
  X <- rnorm(n)
  A <- 3 * X + rnorm(n)
  Y <- rnorm(n)

  ate_ada[s] <- gcomp_adaptive(Y, A, X, alpha)
  ate_dml[s] <- dml_ate(Y, A, X)

  # Nonparametric bootstrap
  boot_ada <- numeric(B)
  boot_dml <- numeric(B)
  n_full   <- 0L
  for (b in seq_len(B)) {
    idx  <- sample(n, n, replace = TRUE)
    Yb <- Y[idx]; Ab <- A[idx]; Xb <- X[idx]

    # Track full-model selection in bootstrap world
    pval_b <- summary(lm(Yb ~ Ab + Xb))$coefficients["Xb", "Pr(>|t|)"]
    if (pval_b < alpha) n_full <- n_full + 1L

    boot_ada[b] <- gcomp_adaptive(Yb, Ab, Xb, alpha)
    boot_dml[b] <- dml_ate(Yb, Ab, Xb)
  }

  se_ada[s]    <- sd(boot_ada)
  se_dml[s]    <- sd(boot_dml)
  cover_ada[s] <- abs(ate_ada[s] - tau) <= 1.96 * se_ada[s]
  cover_dml[s] <- abs(ate_dml[s] - tau) <= 1.96 * se_dml[s]
  boot_frac[s] <- n_full / B

  if (s %% 50 == 0) cat(sprintf("  sim %d / %d done\n", s, nsim))
}

# ---------------------------------------------------------------------------
# Summary output
# ---------------------------------------------------------------------------
oracle_se_ada <- sd(ate_ada)
oracle_se_dml <- sd(ate_dml)

cat("=================================================================\n")
cat(" G-computation ATE: bootstrap failure (adaptive) vs DML (fix)\n")
cat("=================================================================\n\n")
cat(sprintf(" n = %d,  nsim = %d,  B = %d bootstrap replicates\n\n",
            n, nsim, B))
cat(sprintf(" True ATE: %.3f\n", tau))
cat(sprintf(" DGP: Y = noise, A = 3X + noise  (same as slides 41-44)\n\n"))

cat(" -- Adaptive g-computation (outcome model selected by p-value) --\n")
cat(sprintf("  Mean ATE estimate      :  %.4f\n",  mean(ate_ada)))
cat(sprintf("  Oracle SE              :  %.4f\n",  oracle_se_ada))
cat(sprintf("  Mean bootstrap SE      :  %.4f\n",  mean(se_ada)))
cat(sprintf("  Ratio boot SE / oracle :  %.3f\n",  mean(se_ada) / oracle_se_ada))
cat(sprintf("  Bootstrap coverage     :  %.1f%%  (nominal 95%%)\n",
            100 * mean(cover_ada)))
cat(sprintf("  True P(select full)    :   %.1f%%  (by alpha)\n", 100 * alpha))
cat(sprintf("  Mean P(select full in boot world): %.1f%%\n\n",
            100 * mean(boot_frac)))

cat(" -- DML g-computation (cross-fitted nuisances, no selection) --\n")
cat(sprintf("  Mean ATE estimate      :  %.4f\n",  mean(ate_dml)))
cat(sprintf("  Oracle SE              :  %.4f\n",  oracle_se_dml))
cat(sprintf("  Mean bootstrap SE      :  %.4f\n",  mean(se_dml)))
cat(sprintf("  Ratio boot SE / oracle :  %.3f\n",  mean(se_dml) / oracle_se_dml))
cat(sprintf("  Bootstrap coverage     :  %.1f%%  (nominal 95%%)\n\n",
            100 * mean(cover_dml)))

cat("=================================================================\n")
cat(" Interpretation:\n")
cat("  Adaptive g-computation: the outcome model is selected by a\n")
cat("  p-value test (same post-selection mechanism as slides 41-44).\n")
cat("  In the bootstrap world the fitted beta != 0, so the full model\n")
cat("  is selected ~16%% of the time instead of the nominal 5%%.\n")
cat("  Full model uses Var(A|X)=1 in the denominator; reduced model\n")
cat("  uses Var(A)=10.  Bootstrap over-represents the high-variance\n")
cat("  full model -> bootstrap SE >> oracle SE -> OVER-coverage.\n")
cat("\n")
cat("  DML (cross-fitting): the nuisance is fitted on held-out data,\n")
cat("  decoupling the selection/fitting step from the evaluation step.\n")
cat("  No selection event is shared between bootstrap replicates,\n")
cat("  so bootstrap SE ≈ oracle SE -> ~95%% coverage.\n")
cat("=================================================================\n")

# ---------------------------------------------------------------------------
# Plot: empirical vs bootstrap distribution  (two panels)
# ---------------------------------------------------------------------------
cat("\nGenerating plot...\n")
B_plot <- 800
set.seed(77)
X0 <- rnorm(n); A0 <- 3 * X0 + rnorm(n); Y0 <- rnorm(n)

boot_ada_plot <- numeric(B_plot)
boot_dml_plot <- numeric(B_plot)
for (b in seq_len(B_plot)) {
  idx              <- sample(n, n, replace = TRUE)
  boot_ada_plot[b] <- gcomp_adaptive(Y0[idx], A0[idx], X0[idx], alpha)
  boot_dml_plot[b] <- dml_ate(Y0[idx], A0[idx], X0[idx])
}

de_ada <- density(ate_ada,       bw = "SJ")
db_ada <- density(boot_ada_plot, bw = "SJ")
de_dml <- density(ate_dml,       bw = "SJ")
db_dml <- density(boot_dml_plot, bw = "SJ")

pdf("gcomp_bootstrap_failure.pdf", width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

# Panel 1: adaptive g-comp (bootstrap fails)
xl1 <- range(c(de_ada$x, db_ada$x))
yl1 <- c(0, max(de_ada$y, db_ada$y) * 1.18)
plot(de_ada, col = "black", lwd = 2, xlim = xl1, ylim = yl1,
     main = "Adaptive g-comp (post-selection outcome model)",
     xlab = expression(hat(psi)), ylab = "Density",
     sub = sprintf("Oracle SE = %.4f | Boot SE = %.4f | Coverage = %.1f%%",
                   oracle_se_ada, mean(se_ada), 100 * mean(cover_ada)))
lines(db_ada, col = "red", lwd = 2)
abline(v = tau, lty = 2, col = "grey50")
legend("topright", bty = "n",
       legend = c(sprintf("Empirical  (SD = %.4f)", oracle_se_ada),
                  sprintf("Bootstrap  (SD = %.4f)", sd(boot_ada_plot))),
       col = c("black", "red"), lwd = 2)

# Panel 2: DML (bootstrap works)
xl2 <- range(c(de_dml$x, db_dml$x))
yl2 <- c(0, max(de_dml$y, db_dml$y) * 1.18)
plot(de_dml, col = "black", lwd = 2, xlim = xl2, ylim = yl2,
     main = "DML g-comp (cross-fitted nuisances)",
     xlab = expression(hat(psi)), ylab = "Density",
     sub = sprintf("Oracle SE = %.4f | Boot SE = %.4f | Coverage = %.1f%%",
                   oracle_se_dml, mean(se_dml), 100 * mean(cover_dml)))
lines(db_dml, col = "red", lwd = 2)
abline(v = tau, lty = 2, col = "grey50")
legend("topright", bty = "n",
       legend = c(sprintf("Empirical  (SD = %.4f)", oracle_se_dml),
                  sprintf("Bootstrap  (SD = %.4f)", sd(boot_dml_plot))),
       col = c("black", "red"), lwd = 2)

dev.off()
cat("Plot saved to: gcomp_bootstrap_failure.pdf\n")
