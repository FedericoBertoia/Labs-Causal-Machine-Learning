# =============================================================================
# Post-selection inference: bootstrap failure
# Replicates slides 41-44 from Vansteelandt, CML Welcome Lecture
#
# Setup (slide 42):
#   E(Y | A, X) = 0  (true psi = 0, true beta = 0)
#   E(A | X)    = 3X  (strong confounding)
#
# Procedure:
#   Fit full model E(Y|A,X) = psi*A + beta*X via OLS.
#   If p-value for beta < alpha -> report psi_hat (full model).
#   Else                        -> report psi0_hat (model with A only).
#
# Key finding (slide 43):
#   In the bootstrap world, beta_hat != 0, so the test selects the full model
#   much more often than 5%. Vansteelandt reports ~16.5%.
#   -> Bootstrap distribution is more diffuse than the true sampling distribution.
#   -> Bootstrap SE overestimates oracle SE -> bootstrap CI is miscalibrated.
# =============================================================================

set.seed(123)

n     <- 500    # sample size
nsim  <- 200    # simulation replications (increase for smoother results)
B     <- 50     # bootstrap replicates per simulation (increase for accuracy)
alpha <- 0.05   # significance level for variable selection

# ---------------------------------------------------------------------------
# The post-selection estimator
# ---------------------------------------------------------------------------
procedure <- function(Y, A, X, alpha = 0.05) {
  fit_full <- lm(Y ~ A + X)
  pval_X   <- summary(fit_full)$coefficients["X", "Pr(>|t|)"]
  if (pval_X < alpha) {
    coef(fit_full)["A"]   # X significant: report from full model
  } else {
    coef(lm(Y ~ A))["A"] # X not significant: report from reduced model
  }
}

# ---------------------------------------------------------------------------
# Simulation loop
# ---------------------------------------------------------------------------
psi_emp    <- numeric(nsim)  # post-selection estimate each run
cover_boot <- logical(nsim)  # did 95% bootstrap CI cover the truth (0)?
boot_se    <- numeric(nsim)  # bootstrap SE each run
boot_frac  <- numeric(nsim)  # fraction of boot reps selecting full model

for (s in seq_len(nsim)) {
  
  # --- Generate data under the null (psi = 0, beta = 0) ---
  X <- rnorm(n)
  A <- 3 * X + rnorm(n)  # treatment strongly associated with X (confounding)
  Y <- rnorm(n)           # no treatment effect, no direct X -> Y
  
  psi_s      <- procedure(Y, A, X, alpha)
  psi_emp[s] <- psi_s
  
  # --- Parametric bootstrap ---
  # Y* ~ N(Yhat, sigma^2) where Yhat comes from the FITTED full model.
  # Because beta_hat != 0 in finite samples, the bootstrap world has
  # a non-zero "true" beta -> the test falsely selects X more than alpha%.
  fit_bw <- lm(Y ~ A + X)
  Yhat   <- fitted(fit_bw)
  sg     <- sigma(fit_bw)
  
  # Pre-generate all bootstrap Y matrices at once (faster than looping rnorm)
  Ymat <- matrix(Yhat, n, B) + matrix(rnorm(n * B, 0, sg), n, B)
  
  n_full <- 0L
  pb     <- numeric(B)
  for (b in seq_len(B)) {
    Yb   <- Ymat[, b]
    fb   <- lm(Yb ~ A + X)
    pv   <- summary(fb)$coefficients["X", "Pr(>|t|)"]
    if (pv < alpha) {
      pb[b]  <- coef(fb)["A"]
      n_full <- n_full + 1L
    } else {
      pb[b] <- coef(lm(Yb ~ A))["A"]
    }
  }
  
  boot_se[s]    <- sd(pb)
  boot_frac[s]  <- n_full / B
  cover_boot[s] <- (abs(psi_s) <= 1.96 * boot_se[s])  # normal CI
}

# ---------------------------------------------------------------------------
# Summary output
# ---------------------------------------------------------------------------
oracle_se <- sd(psi_emp)

cat("=================================================================\n")
cat(" Post-selection inference: bootstrap failure\n")
cat("=================================================================\n\n")
cat(sprintf(" n = %d,  nsim = %d,  B = %d bootstrap replicates\n\n", n, nsim, B))
cat(sprintf(" True psi                            :  0.000\n"))
cat(sprintf(" Mean of post-selection estimates    :  %.4f   (should be ~0)\n",
            mean(psi_emp)))
cat(sprintf("\n"))
cat(sprintf(" Oracle SE  (empirical SD)           :  %.4f\n", oracle_se))
cat(sprintf(" Mean bootstrap SE                   :  %.4f\n", mean(boot_se)))
cat(sprintf(" Ratio: mean(boot SE) / oracle SE    :  %.3f\n",
            mean(boot_se) / oracle_se))
cat(sprintf("\n"))
cat(sprintf(" Nominal coverage                    :  95.0%%\n"))
cat(sprintf(" Empirical coverage of bootstrap CI  :  %.1f%%\n\n",
            100 * mean(cover_boot)))
cat(sprintf(" True P(select full | beta=0)        :   %.1f%% (by alpha)\n",
            100 * alpha))
cat(sprintf(" Mean P(select full in boot world)   :  %.1f%%\n",
            100 * mean(boot_frac)))
cat(sprintf("   -> Vansteelandt reports ~16.5%% for similar setup\n\n"))
cat("=================================================================\n")
cat(" Interpretation:\n")
cat("  The bootstrap world is generated from the FITTED model where\n")
cat("  beta_hat != 0.  So the selection test concludes beta != 0 in\n")
cat("  ~15% of boot replicates instead of the nominal 5%.\n")
cat("  This makes the bootstrap distribution much more diffuse\n")
cat("  (heavier tails) than the true sampling distribution:\n")
cat("  bootstrap SE >> oracle SE.\n")
cat("  The 95% bootstrap CI therefore OVER-covers (>> 95%),\n")
cat("  not under-covers -- the bootstrap is MISCALIBRATED regardless.\n")
cat("=================================================================\n")

# ---------------------------------------------------------------------------
# Plot: empirical vs bootstrap distribution  (replicates slide 44)
# ---------------------------------------------------------------------------
# Generate a reference bootstrap distribution with more replicates
B_plot <- 800
X_ref  <- rnorm(n); A_ref <- 3 * X_ref + rnorm(n); Y_ref <- rnorm(n)
fit_ref <- lm(Y_ref ~ A_ref + X_ref)
Ymat_ref <- matrix(fitted(fit_ref), n, B_plot) +
  matrix(rnorm(n * B_plot, 0, sigma(fit_ref)), n, B_plot)

pb_ref <- numeric(B_plot)
for (b in seq_len(B_plot)) {
  Yb  <- Ymat_ref[, b]
  fb  <- lm(Yb ~ A_ref + X_ref)
  pv  <- summary(fb)$coefficients["X_ref", "Pr(>|t|)"]
  pb_ref[b] <- if (pv < alpha) coef(fb)["A_ref"] else coef(lm(Yb ~ A_ref))["A_ref"]
}

de <- density(psi_emp, bw = "SJ")
db <- density(pb_ref,  bw = "SJ")
xl <- range(c(de$x, db$x))
yl <- c(0, max(de$y, db$y) * 1.12)

pdf("post_selection_bootstrap.pdf", width = 7.5, height = 5)
plot(de, col = "black", lwd = 2,
     xlim = xl, ylim = yl,
     main = expression(
       paste("Empirical (black) vs Bootstrap (red) distribution of ", hat(psi))),
     xlab = expression(hat(psi)), ylab = "Density",
     sub  = sprintf(
       "n = %d | Oracle SE = %.4f | Boot SE = %.4f | Boot coverage = %.1f%%",
       n, oracle_se, mean(boot_se), 100 * mean(cover_boot)))
lines(db, col = "red", lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
legend("topright", bty = "n",
       legend = c(
         sprintf("Empirical  (SD = %.4f)", oracle_se),
         sprintf("Bootstrap  (SD = %.4f)", sd(pb_ref))
       ),
       col = c("black", "red"), lwd = 2)
dev.off()

cat("\nPlot saved to: post_selection_bootstrap.pdf\n")