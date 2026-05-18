######################################
# New Functions for Poisson link
######################################

.SL.require <- function(package, message = paste('loading required package (', package, ') failed', sep = '')) {
  if(!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

.createLibrary <- function(SL.library) {
  if (is.character(SL.library)) { 
    k <- length(SL.library)
    whichScreen <- matrix(1, nrow = 1, ncol = k)
    screenAlgorithm <- "All"
    library <- data.frame(predAlgorithm = SL.library, rowScreen = 1, stringsAsFactors=FALSE)
  } else if (is.list(SL.library)) {
    predNames <- sapply(SL.library, FUN = "[", 1)
    NumberScreen <- (sapply(SL.library, FUN = length) - 1)
    if (sum(NumberScreen == 0) > 0) {
      for(ii in which(NumberScreen == 0)) {
        SL.library[[ii]] <- c(SL.library[[ii]], "All")
        NumberScreen[ii] <- 1
      }
    }
    screenAlgorithmFull <- unlist(lapply(SL.library, FUN="[", -1))
    screenAlgorithm <- unique(screenAlgorithmFull)
    
    library <- data.frame(predAlgorithm = rep(predNames, times=NumberScreen), rowScreen = match(screenAlgorithmFull, screenAlgorithm), stringsAsFactors = FALSE)
  } else {
    stop('format for SL.library is not recognized')
  }
  
  out <- list(library = library, screenAlgorithm = screenAlgorithm)
  return(out)
}

.check.SL.library <- function(library, addPackages = NULL) {
  if("SL.bayesglm" %in% library) .SL.require('arm', message = 'You have selected bayesglm as a library algorithm but either do not have the arm package installed or it can not be loaded')
  if("SL.cforest" %in% library) .SL.require('party', message = 'You have selected cforest as a library algorithm but either do not have the party package installed or it can not be loaded')
  if("SL.DSA" %in% library) .SL.require('DSA', message = 'You have selected DSA as a library algorithm but either do not have the DSA package installed or it can not be loaded')
  if("SL.gam" %in% library) .SL.require('gam', message = 'You have selected gam as a library algorithm but either do not have the gam package installed or it can not be loaded')
  if("SL.gbm" %in% library) .SL.require('gbm', message = 'You have selected gbm as a library algorithm but either do not have the gbm package installed or it can not be loaded')
  if("SL.glmnet" %in% library) .SL.require('glmnet', message = 'You have selected glmnet as a library algorithm but either do not have the glmnet package installed or it can not be loaded')
  if("SL.knn" %in% library) .SL.require('class', message = 'You have selected knn as a library algorithm but either do not have the class package installed or it can not be loaded')
  if("SL.logreg" %in% library) .SL.require('LogicReg', message = 'You have selected logreg as a library algorithm but either do not have the LogicReg package installed or it can not be loaded')
  if("SL.nnet" %in% library) .SL.require('nnet', message = 'You have selected nnet as a library algorithm but either do not have the nnet package installed or it can not be loaded')
  if("SL.polymars" %in% library) .SL.require('polspline', message = 'You have selected polymars or polyclass as a library algorithm but either do not have the polspline package installed or it can not be loaded')
  if("SL.randomForest" %in% library) .SL.require('randomForest', message = 'You have selected randomForest as a library algorithm but either do not have the randomForest package installed or it can not be loaded')
  if("SL.ridge" %in% library) .SL.require('MASS', message = 'You have selected lm.ridge as a library algorithm but either do not have the MASS package installed or it can not be loaded')
  if("SL.spls" %in% library) .SL.require('spls', message = 'You have selected spls as a library algorithm but either do not have the spls package installed or it can not be loaded')
  if("SL.svm" %in% library) .SL.require('e1071', message = 'You have selected svm as a library algorithm but either do not have the e1071 package installed or it can not be loaded')
  if("SL.ipredbagg" %in% library) .SL.require('ipred', message = 'You have selected bagging as a library algorithm but either do not have the ipred package installed or it can not be loaded')
  if("SL.gbm" %in% library) .SL.require('gbm', message = 'You have selected gbm as a library algorithm but either do not have the gbm package installed or it can not be loaded')
  if("SL.mars" %in% library) .SL.require('mda', message = 'You have selected mars as a library algorithm but either do not have the mda package installed or it can not be loaded')
  if("SL.earth" %in% library) .SL.require('earth', message = 'You have selected earth as a library algorithm but either do not have the earth package installed or it can not be loaded')
  if("SL.caret" %in% library) .SL.require('caret', message = 'You have selected caret as a library algorithm but either do not have the caret package installed or it can not be loaded')
  #	#removing this check, need to replace with requireNamespace per writing R extensions 1.1.3.1
  #	if(!is.null(addPackages)) {
  #	  sapply(addPackages, function(x) require(force(x), character.only = TRUE))
  #	}
  invisible(TRUE)
}



######################################
# New Functions for Poisson link
######################################

SL.earth <- function(Y, X, newX,
                     family,
                     obsWeights = NULL,
                     id = NULL,
                     degree = 2,
                     penalty = 3,
                     nk = max(21, 2 * ncol(X) + 1),
                     pmethod = "backward",
                     nfold = 0,
                     ncross = 1,
                     minspan = 0,
                     endspan = 0,
                     ...) {
  
  .SL.require('earth')
  
  # choose appropriate call to earth based on family
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(
      x = X, y = Y,
      degree = degree, nk = nk, penalty = penalty,
      pmethod = pmethod, nfold = nfold, ncross = ncross,
      minspan = minspan, endspan = endspan, ...
    )
  } else if (family$family == "binomial") {
    fit.earth <- earth::earth(
      x = X, y = Y,
      degree = degree, nk = nk, penalty = penalty,
      pmethod = pmethod, nfold = nfold, ncross = ncross,
      minspan = minspan, endspan = endspan,
      glm = list(family = binomial), ...
    )
  } else if (family$family == "poisson") {
    fit.earth <- earth::earth(
      x = X, y = Y,
      degree = degree, nk = nk, penalty = penalty,
      pmethod = pmethod, nfold = nfold, ncross = ncross,
      minspan = minspan, endspan = endspan,
      glm = list(family = poisson), ...
    )
  } else {
    stop("SL.earth currently supports only gaussian, binomial, or poisson families")
  }
  
  pred <- predict(fit.earth, newdata = newX, type = "response")
  
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}



SL.xgboost = function(Y, X, newX, family, obsWeights, id, ntrees = 1000,
                      max_depth = 4, shrinkage = 0.1, minobspernode = 10,
                      params = list(),
                      nthread = 1,
                      verbose = 0,
                      save_period = NULL,
                      ...) {
  .SL.require("xgboost")
  if(packageVersion("xgboost") < "0.6") stop("SL.xgboost requires xgboost version >= 0.6, try help(\'SL.xgboost\') for details")
  
  if(packageVersion("xgboost") > "3.0") {
    if (family$family == "gaussian") {
      model = xgboost::xgboost(x = X, 
                               y = Y,
                               weights = obsWeights,
                               objective = 'reg:squarederror',
                               nrounds = ntrees,
                               max_depth = max_depth, 
                               min_child_weight = minobspernode, 
                               learning_rate = shrinkage,
                               verbosity = verbose, 
                               nthreads = nthread 
                               # params = params  # no longer including params in xgboost, consider xgb.train()
      )
    }
    if (family$family == "poisson") {
      model = xgboost::xgboost(x = X, 
                               y = Y,
                               weights = obsWeights,
                               objective = 'count:poisson',
                               nrounds = ntrees,
                               max_depth = max_depth, 
                               min_child_weight = minobspernode, 
                               learning_rate = shrinkage,
                               verbosity = verbose, 
                               nthreads = nthread 
                               # params = params  # no longer including params in xgboost, consider xgb.train()
      )
    }
    if (family$family == "binomial") {
      model = xgboost::xgboost(x = X,
                               y = as.factor(Y),
                               weights = obsWeights,
                               objective = "binary:logistic", 
                               nrounds = ntrees,
                               max_depth = max_depth, 
                               min_child_weight = minobspernode, 
                               learning_rate = shrinkage,
                               verbosity = verbose, 
                               nthreads = nthread, 
                               # params = params,
                               eval_metric = "logloss")
    }
    pred = predict(model, newdata = newX, type = "response")
    
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
  } else {
    # X needs to be converted to a matrix first, then an xgb.DMatrix.
    if (!is.matrix(X)) {
      X = model.matrix(~ . - 1, X)
    }
    
    # Convert to an xgboost compatible data matrix, using the sample weights.
    xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    
    # TODO: support early stopping, which requires a "watchlist". See ?xgb.train
    
    if (family$family == "gaussian") {
      # reg:linear was deprecated in version 1.1.1.1, changed to reg:squarederror
      if(packageVersion("xgboost") >= "1.1.1.1") {
        objective <- 'reg:squarederror'
      } else {
        objective <- 'reg:linear'
      }
      model = xgboost::xgboost(data = xgmat, objective=objective, nrounds = ntrees,
                               max_depth = max_depth, min_child_weight = minobspernode, eta = shrinkage,
                               verbose = verbose, nthread = nthread, params = params,
                               save_period = save_period)
    }
    if (family$family == "poisson") {
      # reg:linear was deprecated in version 1.1.1.1, changed to reg:squarederror
      if(packageVersion("xgboost") >= "1.1.1.1") {
        objective <- 'count:poisson'
      } else {
        objective <- 'reg:linear'
      }
      model = xgboost::xgboost(data = xgmat, objective=objective, nrounds = ntrees,
                               max_depth = max_depth, min_child_weight = minobspernode, eta = shrinkage,
                               verbose = verbose, nthread = nthread, params = params,
                               save_period = save_period)
    }
    if (family$family == "binomial") {
      model = xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                               max_depth = max_depth, min_child_weight = minobspernode, eta = shrinkage,
                               verbose = verbose, nthread = nthread, params = params,
                               save_period = save_period, eval_metric = "logloss")
    }
    if (family$family == "multinomial") {
      # TODO: test this.
      model = xgboost::xgboost(data = xgmat, objective="multi:softmax", nrounds = ntrees,
                               max_depth = max_depth, min_child_weight = minobspernode, eta = shrinkage,
                               verbose = verbose, num_class = length(unique(Y)), nthread = nthread,
                               params = params,
                               save_period = save_period)
    }
    
    # Newdata needs to be converted to a matrix first, then an xgb.DMatrix.
    if (!is.matrix(newX)) {
      newX = model.matrix(~ . - 1, newX)
    }
    
    pred = predict(model, newdata = newX)
    
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
  }
}


SL.dbarts2 = function(Y, X, newX, family, obsWeights, id,
                      ntree = 200, ndpost = 1000, nskip = 100,
                      verbose = FALSE, ...) {
  
  requireNamespace("dbarts", quietly = FALSE)
  
  # For Poisson: fit on log scale
  Y_fit <- if (family$family == "poisson") log(Y + 0.5) else Y
  
  model <- dbarts::bart(
    x.train   = X,
    y.train   = Y_fit,
    weights   = obsWeights,
    ntree     = ntree,
    ndpost    = ndpost,
    nskip     = nskip,
    verbose   = verbose,
    keeptrees = TRUE
  )
  
  if (family$family == "gaussian") {
    pred <- colMeans(predict(model, newdata = newX))
  } else if (family$family == "binomial") {
    pred <- colMeans(pnorm(predict(model, newdata = newX)))
  } else if (family$family == "poisson") {
    pred <- pmax(colMeans(exp(predict(model, newdata = newX))), 0)
  }
  
  fit <- list(object = model, family = family$family)
  class(fit) <- "SL.dbarts2"
  return(list(pred = pred, fit = fit))
}



SL.ranger <- function(Y, X, newX, family, obsWeights, id,
                      num.trees = 500,
                      mtry = floor(sqrt(ncol(X))),
                      min.node.size = NULL, ...) {
  .SL.require("ranger")
  
  if (family$family == "binomial") {
    fit <- ranger::ranger(
      y            = as.factor(Y),
      x            = X,
      num.trees    = num.trees,
      mtry         = mtry,
      min.node.size = min.node.size,
      case.weights  = obsWeights,
      probability   = TRUE
    )
  } else {
    # gaussian and poisson: standard regression
    fit <- ranger::ranger(
      y             = Y,
      x             = X,
      num.trees     = num.trees,
      mtry          = mtry,
      min.node.size = min.node.size,
      case.weights  = obsWeights
    )
  }
  
  pred <- predict.SL.ranger(
    list(object = fit, verbose = FALSE),
    newdata = newX,
    family  = family
  )
  
  fit <- list(object = fit, verbose = FALSE)
  class(fit) <- "SL.ranger"
  return(list(pred = pred, fit = fit))
}

SL.gam <- function(Y, X, newX, family, obsWeights, id,
                   deg.gam = 2, cts.num = 4, ...) {
  .SL.require("gam")
  
  # identify continuous vs binary predictors
  cts.x <- apply(X, 2, function(x) length(unique(x)) > cts.num)
  
  # build formula: s() for continuous, linear for binary
  gam.model <- as.formula(
    paste("Y ~", paste(
      ifelse(cts.x,
             paste0("gam::s(", colnames(X), ", ", deg.gam, ")"),
             colnames(X)),
      collapse = " + "
    ))
  )
  
  fit.gam <- gam::gam(
    formula  = gam.model,
    data     = data.frame(Y, X),
    family   = family,           # passes gaussian(), binomial(), or poisson() directly
    weights  = obsWeights
  )
  
  pred <- predict(fit.gam, newdata = newX, type = "response")
  
  fit <- list(object = fit.gam)
  class(fit) <- "SL.gam"
  return(list(pred = pred, fit = fit))
}


SL.hal9001 <- function(Y,
                       X,
                       newX,
                       family,
                       obsWeights,
                       id,
                       max_degree = 2,
                       smoothness_orders = 1,
                       num_knots = 5,
                       ...) {
  # create matrix version of X and newX for use with hal9001::fit_hal
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.null(newX) & !is.matrix(newX)) newX <- as.matrix(newX)
  
  # fit hal
  hal_fit <- hal9001::fit_hal(
    Y = Y, X = X, family = family$family, weights = obsWeights, id = id,
    max_degree = max_degree, smoothness_orders = smoothness_orders,
    num_knots = num_knots, ...
  )
  
  # compute predictions based on `newX` or input `X`
  if (!is.null(newX)) {
    pred <- stats::predict(hal_fit, new_data = newX)
  } else {
    pred <- stats::predict(hal_fit, new_data = X)
  }
  
  # build output object
  fit <- list(object = hal_fit)
  class(fit) <- "SL.hal9001"
  out <- list(pred = pred, fit = fit)
  return(out)
}


######################################
# PREDICTION 
######################################

predict.SL.earth <- function(object, newdata, ...) {
  .SL.require('earth')
  pred <- predict(object$object, newdata = newdata, type = "response")
  return(pred)
}


predict.SL.ranger <- function(object, newdata, family,
                              num.threads = 1,
                              verbose = object$verbose,
                              ...) {
  .SL.require("ranger")
  
  # get predictions from ranger
  pred <- predict(object$object, data = newdata,
                  verbose = verbose,
                  num.threads = num.threads)$predictions
  
  # For binomial family, extract P(Y=1|X)
  if (family$family == "binomial") {
    pred <- pred[, "1"]
  }
  
  # For Poisson family, ensure nonnegative predictions
  if (family$family == "poisson") {
    pred <- pmax(pred, 0)
  }
  
  return(pred)
}


predict.SL.xgboost <- function(object, newdata, family, ...) {
  .SL.require("xgboost")
  if(packageVersion("xgboost") < "0.6") stop("SL.xgboost requires xgboost version >= 0.6, try help(\'SL.xgboost\') for details")
  if(packageVersion("xgboost") > "3.0") {
    pred = predict(object$object, newdata = newdata)
    return(pred)
  } else {
    # newdata needs to be converted to a matrix first
    if (!is.matrix(newdata)) {
      newdata = model.matrix(~ . - 1, newdata)
    }
    pred = predict(object$object, newdata = newdata)
    return(pred)
  }
}




predict.SL.dbarts2 <- function(object, newdata, family, ...) {
  model <- object$object
  if (family$family == "gaussian") {
    pred <- colMeans(predict(model, newdata = newdata))
  } else if (family$family == "binomial") {
    pred <- colMeans(pnorm(predict(model, newdata = newdata)))
  } else if (family$family == "poisson") {
    pred <- pmax(colMeans(exp(predict(model, newdata = newdata))), 0)
  }
  return(pred)
}



predict.SL.gam <- function(object, newdata, ...) {
  .SL.require("gam")
  pred <- predict(object$object, newdata = newdata, type = "response")
  return(pred)
}

predict.SL.hal9001 <- function(object, newdata, ...) {
  # coerce newdata to matrix if not already so
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  
  # generate predictions and return
  pred <- stats::predict(object$object, new_data = newdata)
  return(pred)
}
