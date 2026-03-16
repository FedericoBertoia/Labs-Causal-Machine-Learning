#' Generate indices for cross validates folds.
#'
#' @param n Number of rows.
#' @param cv_folds Number of folds.
#' @return Vector of indices indicating fold membership.
generate_folds <- function(n, cv_folds) {
  inds <- sample(rep(1:cv_folds, ceiling(n/cv_folds))[1:n])
  return(inds)
}




#' Estimate propensity scores parametrically, via logistic regression.
#'
#' @param A A vector of treatment values, in {0, 1}.
#' @param X A matrix or data frame of covariate values.
#' @return A vector of estimated treatment propensities P(A = 1|X) in [0, 1].
estimate_pi_par <- function(A, X) {
  dat <- as.data.frame(cbind(A, X))
  form_pi <- as.formula(paste("A ~", paste(colnames(X), collapse = " + ")))
  pihat <- glm(form_pi, family = binomial(link = "logit"), dat)$fitted.values
  
  return(pihat)
}

#' Bound estimated propensity scores away from 0 and 1.
#'
#' @param pi Vector of estimated propensity scores in [0, 1].
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @return The bounded scores.
bound_pi <- function(pi, lower=0.025, upper=0.975) {
  pi[pi < lower] <- lower
  pi[pi > upper] <- upper
  return(pi)
}


#' Estimate outcome regression function parametrically, via glm.
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values of length n, in {0, 1}.
#' @param X An n-by-p matrix or data frame of covariate values.
#' @param family  A description of the distribution of Y. Should be either 
#' "binomial" for binary outcomes or "gaussian" for continuous outcomes.
#' @return An n-by-2 data frame of estimated regression values E[Y|A=a, X], for
#' a = 0 (column 1) and a = 1 (column 2).
estimate_mu_par <- function(Y, A, X, family = 'binomial') {
  dat <- as.data.frame(cbind(Y, A, X))
  form_mu <- as.formula(paste("Y ~ A +", paste(colnames(X), collapse = " + ")))
  mumod <- glm(form_mu, family = family, dat)
  dat$A <- 0; mu0 <- predict(mumod, newdata = dat, type = "response")
  dat$A <- 1; mu1 <- predict(mumod, newdata = dat, type = "response")
  mu <- A*mu1 + (1 - A)*mu0
  muhat <- as.data.frame(cbind(mu, mu0, mu1))
  
  return(muhat)
}


#' Estimate propensity scores nonparametrically using SuperLearner.
#'
#' @param A A vector of treatment values, in {0, 1}.
#' @param X A matrix or data frame of covariate values.
#' @param sl_lib A character vector of learner functions created via
#' SuperLearner, as in sl_lib <- SuperLearner::create.Learner(...)$names, where
#' ... are the learner type and learner parameters. These are passed to
#' SuperLearner() to estimate the propensity score function.
#' @param nsplits The number of splits to use for sample splitting, if fold
#' indices aren't provided to the split_inds argument.
#' @param split_inds Fold indices for sample splitting, if nsplits is NULL. One
#' of split_inds or nsplits must be provided. split_inds takes precedence if
#' it's not NULL.
#' @param ... Additional arguments passed to SuperLearner.
#' @return A vector of estimated treatment propensities P(A = 1|X) in [0, 1].
estimate_pi_nonpar <- function(A, X, sl_lib, nsplits=NULL, split_inds=NULL, ...) {
  if (is.null(nsplits) & is.null(split_inds)) {
    stop("One of nsplits or split_inds must be provided.")
  }
  set.seed(123)
  if (is.null(split_inds)) {
    split_inds <- generate_folds(length(A), nsplits)
  }
  else {
    nsplits = length(unique(split_inds))
  }
  X <- as.data.frame(X)
  n <- length(A)
  pihat <- rep(NA, n)
  coefs <- rep(0, length(sl_lib))
  
  for (vfold in 1:nsplits) {
    train <- split_inds != vfold
    test <- split_inds == vfold
    if (nsplits == 1) {
      train <- test
    }
    model <- SuperLearner(A[train], X[train, ], newX = X[test, ],
                          family=binomial(link = "logit"),
                          SL.library = sl_lib, ...)
    pihat[test] <- model$SL.predict
    coefs <- coefs + model$coef
  }
  coefs <- coefs/nsplits
  
  return(list(preds = pihat, coefs = coefs))
}


#' Estimate outcome regression function nonparametrically using SuperLearner.
#'
#' Outcomes are treated as continuous, meaning that the error distribution in
#' SuperLearner() is set to "gaussian."
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param X A matrix or data frame of covariate values.
#' @param sl_lib A character vector of learner functions created via
#' SuperLearner, as in sl_lib <- SuperLearner::create.Learner(...)$names, where
#' ... are the learning type and learner parameters. These are passed to
#' SuperLearner() to estimate the propensity score function.
#' @param family  A description of the distribution of Y. Should be either 
#' "binomial" for binary outcomes or "gaussian" for continuous outcomes.
#' @param nsplits The number of splits to use for sample splitting, if fold
#' indices aren't provided to the split_inds argument.
#' @param split_inds Fold indices for sample splitting, if nsplits is NULL. One
#' of split_inds or nsplits must be provided. split_inds takes precedence if
#' it's not NULL.
#' @param ... Additional arguments passed to SuperLearner.
#' @return An n-by-3 data frame of estimated regression values, where column
#' 1 = E[Y|A, X], column 2 = E[Y|A=0, X], and column 3 = E[Y|A=1, X].
#' a = 0 (column 1) and a = 1 (column 2).
estimate_mu_nonpar <- function(Y, A, X, sl_lib, family = 'binomial', nsplits=NULL, 
                               split_inds=NULL, ...) {
  if (is.null(nsplits) & is.null(split_inds)) {
    stop("One of nsplits or split_inds must be provided.")
  }
  set.seed(123)
  if (is.null(split_inds)) {
    split_inds <- generate_folds(length(A), nsplits)
  }
  else {
    nsplits = length(unique(split_inds))
  }
  XA <- as.data.frame(cbind(X, A))
  X <- as.data.frame(X)
  n <- length(A)
  muhat <- as_tibble(matrix(0, nrow = n, ncol = 3,
                            dimnames = list(c(), c("mu0", "mu1", "mu"))))
  coefs <- rep(0, length(sl_lib))
  
  for (vfold in 1:nsplits) {
    train <- split_inds != vfold
    test <- split_inds == vfold
    if (nsplits == 1) {
      train <- test
    }
    model0 <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ],
                           newX = X[test, ], family = family, SL.library = sl_lib, ...)
    muhat[test, "mu0"] <- model0$SL.predict
    coefs <- coefs + model0$coef
    model1 <- SuperLearner(Y[(A == 1) & train], X[(A == 1) & train, ],
                           newX = X[test, ], family = family, SL.library = sl_lib, ...)
    muhat[test, "mu1"] <- model1$SL.predict
    coefs <- coefs + model1$coef
  }
  coefs <- coefs/(2*nsplits)
  muhat <- muhat %>%
    mutate(mu = A*mu1 + (1 - A)*mu0)
  
  return(list(preds = muhat, coefs = coefs))
}




#' Compute regression, aka g-computation, estimate of the ATE.
#'
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = E[Y|A=0, X] and column 2 = E[Y|A=1, X]. The columns must be named
#' mu0 and mu1.
#' @return A named vector with elements "est": the estimated ATE using the
#' g-computation estimator; and "se", which is always NA and which is a
#' placeholder for the bootstrap-estimated standard error. This element is
#' included to make the output of this function equivalent to the output from
#' other estimator functions.
est_gcomp <- function(muhat) {
  est <- mean(muhat$mu1 - muhat$mu0)
  return("est" = est)
}





#' Compute IPW estimate of ATE.
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param pihat A vector of estimated propensity scores. Estimation can be poor
#' when values are close to 0 or 1, so it may be advisable to bound these in,
#' say, [0.025, 0.975].
#' @return A named vector with elements "est", the estimated ATE, and "se", the
#' standard error, estimated analytically.
est_ipw <- function(Y, A, pihat) {
  weights <- A * (mean(A) / pihat) + (1 - A) * ((1 - mean(A)) / (1 - pihat))
  mod <- lm(Y ~ A, weights = weights)
  est <- coef(mod)[2]
  se <- sqrt(vcovHC(mod, type = "HC")[2,2])
  
  return(c("est" = unname(est), "se" = se))
}






#' Compute tMLE estimate of ATE using tmle::tmle().
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param X n-by-p matrix of covariate data.
#' @param pihat A vector of estimated propensity scores. Estimation can be poor
#' when values are close to 0 or 1, so it may be advisable to bound these in,
#' say, [0.025, 0.975].
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = E[Y|A=0, X] and column 2 = E[Y|A=1, X]. The columns must be named
#' mu0 and mu1.
#' @return A named vector with elements "est", the estimated ATE, and "se", the
#' standard error, estimated analytically.
est_tmle <- function(Y, A, X, muhat, pihat) {
  Q <- cbind(muhat$mu0, muhat$mu1)
  mod <- tmle(Y, A, X, Q = Q, g1W = pihat)
  est <- mod$estimates$ATE$psi
  se <- sqrt(mod$estimates$ATE$var.psi)
  
  return(c("est" = unname(est), "se" = se))
}






#' Compute AIPW estimate of ATE.
#'
#' @param Y A vector of continuous outcomes of length n.
#' @param A A vector of treatment values, in {0, 1}.
#' @param pihat A vector of estimated propensity scores. Estimation can be poor
#' when values are close to 0 or 1, so it may be advisable to bound these in,
#' say, [0.025, 0.975].
#' @param muhat 2-column data frame of estimated regression values, where
#' column 1 = E[Y|A=0, X] and column 2 = E[Y|A=1, X]. The columns must be named
#' mu0 and mu1.
#' @return A named vector with elements "est", the estimated ATE, and "se", the
#' standard error, estimated analytically.
est_aipw <- function(Y, A, pihat, muhat) {
  vec <- (((2*A - 1)*(Y - muhat$mu))/
            ((2*A - 1)*pihat + (1 - A)) + muhat$mu1 - muhat$mu0)
  est <- mean(vec)
  se <- sd(vec)/sqrt(length(Y))
  
  return(c("est" = unname(est), "se" = se))
}










bootstrap_se_gcomp <- function(Y, A, X, family = 'binomial', B = 1000, seed = 123,
                               nonpar = FALSE, sl_lib = NULL, nsplits = 5,
                               stratified = FALSE) {
  set.seed(seed)
  ate_boot <- numeric(B)
  
  if (stratified) {
    strata <- interaction(Y, A)
    stratum_indices <- split(seq_along(Y), strata)
  }
  
  for (b in 1:B) {
    tryCatch({
      if (stratified) {
        boot_idx <- unlist(lapply(stratum_indices, function(idx) 
          sample(idx, length(idx), replace = TRUE)))
      } else {
        boot_idx <- sample(length(Y), replace = TRUE)
      }
      
      Y_boot <- Y[boot_idx]
      A_boot <- A[boot_idx]
      X_boot <- X[boot_idx, , drop = FALSE]
      
      muhat_boot <- if (nonpar) {
        estimate_mu_nonpar(Y_boot, A_boot, X_boot, sl_lib, nsplits, family=family)$preds
      } else {
        estimate_mu_par(Y_boot, A_boot, X_boot, family)
      }
      
      ate_boot[b] <- est_gcomp(muhat_boot)
      
      if (nonpar) {
        cat("Bootstrap iteration", b, "completed.\n")
      }
    })
  }
  
  sd(ate_boot, na.rm = TRUE)
}





















