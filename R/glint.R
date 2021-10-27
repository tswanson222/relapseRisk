### ------------------------------------------------------------------------ ###
### -------------------------- CUSTOM FUNCTIONS ---------------------------- ###
### ------------------------------------------------------------------------ ###
##### glint: reshaping glinternet output, and performs criterion-based selection
glint <- function(x, y = 1, method = 'AIC', m = 'zzzall', type = 'default',
                  gamma = .5, nlam = 50, useSE = TRUE, nfolds = 10,
                  allCoef = FALSE, fit = NULL, lambda = NULL,
                  intOnly = FALSE, lamNum = FALSE, ...){
  if(!missing(x)){if(is(x, 'train')){fit <- x; x <- NULL}}
  suppressMessages(invisible(require(glinternet)))
  if(length(method) > 1){fit <- method; method <- 'aic'}
  suppressMessages(invisible(require(glinternet)))
  suppressMessages(invisible(require(Matrix)))
  dots <- function(k, ff = NULL, safe = TRUE){
    ff <- deparse(substitute(ff))
    ff <- tryCatch({eval(parse(text = ff))}, error = function(e){list()})
    if(length(ff) > 1){if(!isTRUE(attr(ff, 'train'))){return(ff)}}
    k <- paste0(deparse(substitute(k)), collapse = '')
    tryCatch({eval(parse(text = k))}, error = function(e){
      k <- paste0(gsub(', ...)$', '', k), ')')
      if(safe){
        tryCatch({eval(parse(text = k))}, error = function(ee){list()})
      } else {
        eval(parse(text = k))
      }
    })
  }
  if(!is.null(fit)){
    if(is(fit, 'list')){
      if('fitobj' %in% names(fit)){
        fitobj <- fit
        fit <- fit$fitobj
      }
      fobj <- intersect(names(fit), c('fitCV', 'fit'))
      if(length(fobj) == 0){fobj <- 'fit0'}
      stopifnot(fobj %in% names(fit))
      assign(fobj, fit[[fobj]], pos = -1)
    } else if(is(fit, 'train')){
      if(is.null(lambda)){lambda <- fit$bestTune[, 'lambda']}
      if(identical(m, 'zzzall')){m <- fit$bestTune[, 'm']}
      if(ifelse(!missing(x), is.null(x), TRUE)){
        x <- fit$trainingData[, -ncol(fit$trainingData)]
        y <- fit$trainingData[, ncol(fit$trainingData)]
      }
      fit <- fit$finalModel
      attr(fit, 'train') <- TRUE
    } else if('lambdaOpt' %in% names(fit)){
      lambda <- fit$lambdaOpt
      m <- attr(fit, 'm')
      y <- fit$y
      attr(fit, 'train') <- TRUE
    }
    if('nFolds' %in% names(fit) | exists('fitCV', inherits = FALSE)){
      yy <- switch(2 - exists('fitCV', inherits = FALSE), fitCV$fitted, fit$fitted)
      p <- length(switch(2 - exists('fitCV', inherits = FALSE), fitCV$numLevels, fit$numLevels))
      if(missing(x)){x <- matrix(0, nrow = length(yy), ncol = p)}
      if('nFolds' %in% names(fit)){fitCV <- fit}
      family <- fitCV$family
      method <- 'cv'
    } else {
      family <- fit$family
      if(startsWith(tolower(method), 'c')){method <- 'aic'}
      if(missing(x)){x <- matrix(0, nrow = nrow(fit$fitted), ncol = length(fit$numLevels))}
    }
    if(exists('fitobj', inherits = FALSE) & is.null(colnames(x))){
      p <- setdiff(rownames(fitobj$coefs), '(Intercept)')
      p <- p[-grep(':', p)]
      colnames(x) <- p
    }
    if(missing(y)){
      y <- switch(2 - startsWith(family, 'g'), numeric(nrow(x)),
                  sample(0:1, nrow(x), TRUE))
    }
    if('m' %in% names(attributes(fit))){m <- attr(fit, 'm')}
  }
  if(length(y) == 1){
    y0 <- y
    y <- x[, y0]
    if(is.character(y0)){
      x <- as.matrix(x[, setdiff(colnames(x), y0), drop = FALSE])
    } else {
      x <- as.matrix(x[, -y0, drop = FALSE])
    }
  }
  family <- ifelse(dim(table(y)) == 2, 'binomial', 'gaussian')
  if(family == 'binomial' & (is.factor(y) | is.character(y))){
    y <- factor(y)
    levels(y) <- 0:1
    y <- as.numeric(as.character(y))
  }
  method <- match.arg(tolower(method), c('cv', 'aic', 'bic', 'ebic'))
  if(is.null(colnames(x))){colnames(x) <- paste0('X', 1:ncol(x))}
  if(is.data.frame(x)){x <- as.matrix(x)}
  if(identical(type, 'default')){type <- rep(1, ncol(x))}
  if(identical(m, 'zzzall') | identical(m, 0) | identical(m, ncol(x) + 1)){m <- 1:ncol(x)}
  if(is.character(m)){
    stopifnot(all(m %in% colnames(x)))
    m <- which(colnames(x) %in% m)
  }
  stopifnot(length(type) == ncol(x))
  stopifnot(is.numeric(m))
  m <- intersect(1:ncol(x), m)
  if(!identical(lamNum, FALSE)){
    lam <- lambdaGrid(x, y, nlam)
    lambda <- lam[lamNum]
  }
  if(method == 'cv'){
    fitCV <- dots(
      glinternet.cv(x, y, type, interactionCandidates = m,
                    family = family, nFolds = nfolds, lambda = lambda,
                    nLambda = nlam, ...), ff = fitCV
    )
    attr(fitCV, 'm') <- m
    which.lam0 <- which(fitCV$lambda == fitCV$lambdaHat)
    while(is.null(fitCV$glinternetFit$activeSet[[which.lam0]]) & !intOnly){
      which.lam0 <- which.lam0 + 1
    }
    fitCV$lambdaHat <- fitCV$lambda[which.lam0]
    fit0 <- dots(
      glinternet(x, y, type, lambda = fitCV$lambdaHat,
                 interactionCandidates = m,
                 family = family, ...)
    )
    attr(fit0, 'm') <- m
    if(useSE == TRUE){
      lamlist <- fitCV$lambda
      errm <- fitCV$cvErr
      errse <- fitCV$cvErrStd <- fitCV$cvErrStd/sqrt(nfolds)
      o <- which.min(errm)
      lamhat <- lamlist[o]
      oo <- errm <= errm[o] + errse[o]
      fitCV$lambdaHat1Std <- lamlist[oo & lamlist >= lamhat][1]
    }
    which.lam0 <- which(fitCV$lambda == fitCV$lambdaHat)
    while(is.null(fitCV$glinternetFit$activeSet[[which.lam0]]) & !intOnly){
      which.lam0 <- which.lam0 + 1
    }
    fitCV$lambdaHat <- fitCV$lambda[which.lam0]
    which.lam1se <- which(fitCV$lambda == fitCV$lambdaHat1Std)
    while(is.null(fitCV$glinternetFit$activeSet[[which.lam1se]]) & !intOnly){
      which.lam1se <- which.lam1se + 1
    }
    fitCV$lambdaHat1Std <- fitCV$lambda[which.lam1se]
    fit0 <- dots(
      glinternet(x, y, type, lambda = fitCV$lambdaHat,
                 interactionCandidates = m,
                 family = family, ...)
    )
    fit1se <- dots(
      glinternet(x, y, type, lambda = fitCV$lambdaHat1Std,
                 interactionCandidates = m,
                 family = family, ...)
    )
    attr(fit0, 'm') <- m
    attr(fit1se, 'm') <- m
    attributes(fit1se)$useSE <- attributes(fitCV)$useSE <- useSE
    mod0 <- coef(fit0)[[2]]
    mod1se <- coef(fit1se)[[2]]
    mains <- 1:ncol(x)
    ints <- t(combn(mains, 2))
    ints2 <- as.numeric(apply(ints, 1, paste, collapse = ""))
    if(length(mod0$mainEffects$cont) != 0){
      if(any(!mains %in% mod0$mainEffects$cont)){
        mod0miss1 <- mains[!mains %in% mod0$mainEffects$cont]
        mod0coefs1 <- c(mod0$mainEffectsCoef$cont,
                        rep(0, length(mod0miss1)))[order(c(mod0$mainEffects$cont, mod0miss1))]
      } else {
        mod0coefs1 <- mod0$mainEffectsCoef$cont[order(mod0$mainEffects$cont)]
      }
    } else {
      mod0coefs1 <- rep(0, length(mains))
    }
    if(length(mod1se$mainEffects$cont) != 0){
      if(any(!mains %in% mod1se$mainEffects$cont)){
        mod1semiss1 <- mains[!mains %in% mod1se$mainEffects$cont]
        mod1secoefs1 <- c(mod1se$mainEffectsCoef$cont,
                          rep(0, length(mod1semiss1)))[order(c(mod1se$mainEffects$cont, mod1semiss1))]
      } else {
        mod1secoefs1 <- mod1se$mainEffectsCoef$cont[order(mod1se$mainEffects$cont)]
      }
    } else {
      mod1secoefs1 <- rep(0, length(mains))
    }
    if(length(mod0$interactions$contcont) != 0){
      mod0ints1 <- as.numeric(apply(mod0$interactions$contcont, 1, paste, collapse = ""))
      if(nrow(ints) != nrow(mod0$interactions$contcont)){
        mod0coefs2 <- rep(0, nrow(ints))
        mod0coefs2[which(ints2 %in% mod0ints1)] <- mod0$interactionsCoef$contcont
      } else {
        mod0coefs2 <- mod0$interactionsCoef$contcont[match(mod0ints1, ints2)]
      }
    } else {
      mod0coefs2 <- rep(0, nrow(ints))
    }
    if(length(mod1se$interactions$contcont) != 0){
      mod1seints1 <- as.numeric(apply(mod1se$interactions$contcont, 1, paste, collapse = ""))
      if(nrow(ints) != nrow(mod1se$interactions$contcont)){
        mod1secoefs2 <- rep(0, nrow(ints))
        mod1secoefs2[which(ints2 %in% mod1seints1)] <- mod1se$interactionsCoef$contcont
      } else {
        mod1secoefs2 <- mod1se$interactionsCoef$contcont[match(mod1seints1, ints2)]
      }
    } else {
      mod1secoefs2 <- rep(0, nrow(ints))
    }
    mod0coefs1 <- c(fit0$betahat[2][[1]][1], mod0coefs1)
    mod1secoefs1 <- c(fit1se$betahat[2][[1]][1], mod1secoefs1)
    mod0 <- unlist(c(mod0coefs1, mod0coefs2))
    mod1se <- unlist(c(mod1secoefs1, mod1secoefs2))
    coefs <- Matrix(cbind(mod0, mod1se), sparse = TRUE)
    allNames <- c(colnames(x), apply(combn(colnames(x), 2), 2, paste, collapse = ":"))
    rownames(coefs) <- c('(Intercept)', allNames)
    output <- list(mod0 = allNames[mod0 != 0], mod1se = allNames[mod1se != 0], coefs = coefs,
                   fitobj = list(fitCV = fitCV, fit0 = fit0, fit1se = fit1se))
    output$mod0 <- setdiff(output$mod0, '(Intercept)')
    output$mod1se <- setdiff(output$mod1se, '(Intercept)')
  } else {
    criterion <- toupper(method)
    fit <- dots(
      glinternet(x, y, type, interactionCandidates = m, family = family,
                 nLambda = nlam, lambda = lambda, ...), ff = fit
    )
    attr(fit, 'm') <- m
    coefs <- coef(fit)[-1]
    mains <- 1:ncol(x)
    ints <- t(combn(mains, 2))
    ints2 <- as.numeric(apply(ints, 1, paste, collapse = ""))
    vs <- colnames(x)
    vs1 <- vs[m]
    if(length(vs1) > 1){vs1 <- paste0("(", paste(vs1, collapse = " + "), ")")}
    vs2 <- as.formula(paste0("~ . * ", vs1))
    dmnames <- colnames(model.matrix(vs2, data.frame(x)))[-1]
    allCoefs <- lapply(coefs, function(z){
      zmain1 <- z$mainEffects$cont
      zmain2 <- z$mainEffectsCoef$cont
      if(length(zmain1) != 0){
        if(any(!mains %in% zmain1)){
          zmiss1 <- mains[!mains %in% zmain1]
          zcoefs1 <- c(zmain2, rep(0, length(zmiss1)))[order(c(zmain1, zmiss1))]
        } else {
          zcoefs1 <- zmain2[order(zmain1)]
        }
      } else {
        zcoefs1 <- rep(0, length(mains))
      }
      zint1 <- z$interactions$contcont
      zint2 <- z$interactionsCoef$contcont
      if(length(zint1) != 0){
        zints1 <- as.numeric(apply(zint1, 1, paste, collapse = ""))
        if(nrow(ints) != nrow(zint1)){
          zcoefs2 <- rep(0, nrow(ints))
          zcoefs2[which(ints2 %in% zints1)] <- zint2
        } else {
          zcoefs2 <- zint2[match(zints1, ints2)]
        }
      } else {
        zcoefs2 <- rep(0, nrow(ints))
      }
      betas <- unlist(c(zcoefs1, zcoefs2))
      names(betas) <- c(colnames(x), apply(combn(colnames(x), 2), 2, paste, collapse = ":"))
      betas <- betas[which(names(betas) %in% dmnames)]
      return(betas)
    })
    n_neighbors <- sapply(allCoefs, function(z) sum(z != 0))
    if(family == 'gaussian'){
      LL_models <- sapply(1:length(allCoefs), function(z){
        s2 <- sum((y - fit$fitted[, z + 1])^2)/length(y)
        if(identical(s2, 0)){s2 <- 1} # FLAW
        sum(dnorm(y, mean = fit$fitted[, z + 1], sd = sqrt(s2), log = TRUE))
      })
    } else {
      LL_models <- sapply(1:length(allCoefs), function(z){
        prob <- fit$fitted[, z + 1]
        #prob <- exp(prob)/(1 + exp(prob))
        sum(log(prob * y + (1 - prob) * (1 - y)))
      })
    }
    p <- length(mains) + nrow(ints)
    if(all(m == 0)){
      p <- length(mains)
    } else {
      p <- length(c(colnames(x), paste0(colnames(x)[-m], ":", colnames(x)[m])))
    }
    ic_lambda <- -2 * LL_models + n_neighbors * ifelse(
      criterion == "AIC", 2, log(nrow(x))) + ifelse(
        criterion == "EBIC", list(2 * gamma * n_neighbors * log(p)), list(0))[[1]]
    betas <- allCoefs[[which.min(ic_lambda)]]
    lambda_min <- fit$lambda[which.min(ic_lambda) + 1]
    coefs <- Matrix(betas, sparse = TRUE)
    rownames(coefs) <- names(betas)
    fitobj <- list(fit = fit, fit0 = NA, crit = ic_lambda)
    if(length(fit$lambda) > 2){
      fitobj$fit0  <- dots(
        glinternet(x, y, type, interactionCandidates = m,
                   lambda = lambda_min, family = family, ...)
      )
    } else {
      fitobj$fit0 <- fit
    }
    attr(fitobj$fit0, 'm') <- m
    names(fitobj)[3] <- criterion
    output <- list(mod0 = names(betas)[betas != 0], coefs = coefs,
                   fitobj = fitobj, allCoefs = allCoefs)
    output$coefs <- rbind(fitobj$fit0$betahat[2][[1]][1], output$coefs)
    rownames(output$coefs)[1] <- '(Intercept)'
    if(!allCoef){output$allCoefs <- NULL}
    if(length(fitobj[[3]]) == 1){output <- append(output[1:2], fitobj[-1])}
  }
  return(output)
}


##### lambdaGrid: Just returns lambda values from glinternet
lambdaGrid = function(X, Y, nLambda = 50, m = NULL, trim = FALSE, lambda = NULL, numLevels,
                      lambdaMinRatio = 0.01, interactionPairs = NULL, screenLimit = NULL,
                      numToFind = NULL, family = c("gaussian", "binomial"),
                      tol = 1e-5, maxIter = 5000, verbose = FALSE,
                      numCores = 1, lamOnly = TRUE){
  # get call and family
  interactionCandidates <- m
  suppressMessages(invisible(require(glinternet)))
  if(dim(table(Y)) == 2){
    if(is.factor(Y)){
      levels(Y) <- 0:1
      Y <- as.numeric(as.character(Y))
    } else if(is.character(Y)){
      Y <- ifelse(Y == unique(Y)[1], 0, 1)
    }
  }
  family <- ifelse(dim(table(Y)) == 2, 'binomial', 'gaussian')
  thisCall = match.call()
  family = match.arg(family)
  if(missing(numLevels)){numLevels <- rep(1, ncol(X))}
  # make sure inputs are valid
  n = length(Y)
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  stopifnot(n==nrow(X), pCat+pCont==ncol(X), family=="gaussian"||family=="binomial")
  if (family=="binomial" && !all(Y %in% 0:1)) {
    stop("Error:family=binomial but Y not in {0,1}")
  }
  for (i in 1:ncol(X)) {
    if (numLevels[i]>1 && max(X[, i])>=numLevels[i]) {
      stop(sprintf("Column %d of X is categorical, but not coded as {0, 1, ...}. Refer to glinternet help on what the X argument should be.", i))
    }
  }

  contIndices = which(numLevels == 1)
  catIndices = which(numLevels > 1)

  # specific interaction pairs
  if (!is.null(interactionPairs)) {
    # sanity check
    if (!is.matrix(interactionPairs)) {
      stop("interactionPairs must be either NULL or a n x 2 matrix of indices")
    }
    if (!is.null(interactionCandidates)) {
      stop("If interactionPairs is set, interactionCandidates must be NULL.")
    }
    pairs = list(contcont=NULL, catcat=NULL, catcont=NULL)
    for (i in 1:nrow(interactionPairs)) {
      left = interactionPairs[i, 1]
      right = interactionPairs[i, 2]
      if (numLevels[left] == 1 && numLevels[right] == 1) {
        pairs$contcont = c(pairs$contcont, which(contIndices %in% c(left, right)))
      } else if (numLevels[left] == 1) {
        pairs$catcont = c(pairs$catcont, which(catIndices == right), which(contIndices == left))
      } else if (numLevels[right] == 1) {
        pairs$catcont = c(pairs$catcont, which(catIndices == left), which(contIndices == right))
      } else {
        pairs$catcat = c(pairs$catcat, which(catIndices %in% c(left, right)))
      }
    }
    # convert to matrices
    pairs = lapply(pairs, function(x) {
      if (!is.null(x)) {
        return(matrix(x, ncol=2, byrow=TRUE))
      } else {
        return(NULL)
      }
    })
    interactionPairs = pairs
  }

  # separate into categorical and continuous parts
  if (pCont > 0) {
    continuousCandidates = NULL
    Z = as.matrix(apply(as.matrix(X[, contIndices]), 2, glinternet:::standardize))
    if (!is.null(interactionCandidates)) {
      continuousCandidates = which(contIndices %in% interactionCandidates)
    }
  } else {
    Z = NULL
    continuousCandidates = NULL
  }
  if (pCat > 0){
    categoricalCandidates = NULL
    levels = numLevels[catIndices]
    Xcat = as.matrix(X[, catIndices])
    if (!is.null(interactionCandidates)) {
      categoricalCandidates = which(catIndices %in% interactionCandidates)
    }
  } else {
    levels = NULL
    Xcat = NULL
    categoricalCandidates = NULL
  }

  # compute variable norms
  res = Y - mean(Y)
  candidates = glinternet:::get_candidates(Xcat, Z, res, n, pCat, pCont, levels, interactionPairs, categoricalCandidates, continuousCandidates, screenLimit, numCores=numCores)

  # lambda grid if not user provided
  if (is.null(lambda)) {
    nLambda <- nLambda + ifelse(!is.logical(trim), length(trim), ifelse(trim, 2, 0))
    lambda = glinternet:::get_lambda_grid(candidates, nLambda, lambdaMinRatio)
    if(lamOnly){
      if(isTRUE(trim)){trim <- c(1, length(lambda))}
      if(!identical(trim, FALSE)){lambda <- lambda[-trim]}
      return(lambda)
    }
  } else {
    stopifnot(min(lambda) > 0)
    if (any(diff(lambda) > 0)) {
      stop("Error: input lambda sequence is not monotone decreasing.")
    }
    lambdaMax = max(get_lambda_grid(candidates, nLambda, lambdaMinRatio))
    nLambda = length(lambda)
    if (nLambda == 1) {
      lambda = sort(c(lambda, lambdaMax), decreasing=TRUE)
      nLambda = 2
    }
  }

  # initialize storage for results
  fitted = matrix(mean(Y), n, nLambda)
  activeSet = vector("list", nLambda)
  betahat = vector("list", nLambda)
  betahat[[1]] = ifelse(family=="gaussian", mean(Y), -log(1/mean(Y)-1))
  objValue = rep(0, nLambda)
  objValue[1] = ifelse(family=="gaussian", sum(res^2)/(2*n), -mean(Y)*betahat[[1]]+log(1/(1-mean(Y))))

  # ever-active set + sequential strong rules + group lasso
  for (i in 2:nLambda){
    if (verbose) {
      cat("lambda ", i, ": ", lambda[i], "\n")
    }
    activeSet[[i]] = strong_rules(candidates, lambda[i], lambda[i-1])
    betahat[[i]] = initialize_betahat(activeSet[[i]], activeSet[[i-1]], betahat[[i-1]], levels)
    while (TRUE) {
      # group lasso on strong set
      solution = group_lasso(Xcat, Z, Y, activeSet[[i]], betahat[[i]], levels, lambda[i], family, tol, maxIter, verbose)
      activeSet[[i]] = solution$activeSet
      betahat[[i]] = solution$betahat
      res = solution$res
      objValue[i] = solution$objValue
      # check kkt conditions on the rest
      check = check_kkt(Xcat, Z, res, n, pCat, pCont, levels, candidates, activeSet[[i]], lambda[i], numCores)
      candidates$norms = check$norms
      if (check$flag) {
        break
      }
      betahat[[i]] = initialize_betahat(check$activeSet, activeSet[[i]], betahat[[i]], levels)
      activeSet[[i]] = check$activeSet
    }
    # update the candidate set if necessary
    if (!is.null(screenLimit) && (screenLimit<pCat+pCont) && i<nLambda) {
      candidates = get_candidates(Xcat, Z, res, n, pCat, pCont, levels, interactionPairs, categoricalCandidates, continuousCandidates, screenLimit, activeSet[[i]], candidates$norms, numCores)
    }
    # get fitted values
    fitted[, i] = Y - res
    # compute total number of interactions found
    if (!is.null(numToFind)) {
      numFound = sum(sapply(activeSet[[i]][3:5], function(x) ifelse(is.null(x), 0, nrow(x))))
      if (numFound >= numToFind) {
        break
      }
    }
  }

  # rescale betahat
  Z = as.matrix(X[, numLevels==1])
  betahatRescaled = lapply(1:i, function(j) rescale_betahat(activeSet[[j]], betahat[[j]], Xcat, Z, levels, n))

  output = list(call=thisCall, fitted=fitted[, 1:i], lambda=lambda[1:i], objValue=objValue, activeSet=activeSet[1:i], betahat=betahatRescaled[1:i], numLevels=numLevels, family=family)
  class(output) = "glinternet"

  return (output)
}


##### ROCcurve:
ROCcurve <- function(y, X = NULL, model = NULL, plot = FALSE, optPoint = TRUE,
                     grid = FALSE, grid_lty = 3, grid_lwd = 1.5, grid_col = "lightgray",
                     midline = TRUE, midline_lty = 2, midline_lwd = 2, midline_col = "red",
                     pt_pch = 23, pt_border = "black", pt_col = "green", thresh = TRUE,
                     roc_lty = 1, roc_lwd = 2, roc_col = "black", prc = FALSE){
  if(!is.null(model)){if(missing(y)){y <- unname(model$y)}}
  stopifnot(dim(table(y)) == 2)
  if(is.factor(y) | is.character(y)){
    y <- factor(y)
    levels(y) <- 0:1
    y <- as.numeric(as.character(y))
  }
  stopifnot(all(names(table(y)) %in% c("0", "1")))
  if(is.null(model)){
    stopifnot(!is.null(X))
    X <- as.data.frame(X)
    model <- glm(y ~ ., data = X, family = binomial)
    predProbs <- predict(model, type = "response")
  } else if(is(model, 'train')){
    predProbs <- predict(model, X, 'prob')[, 2]
    #lam <- model$finalModel$lambdaOpt
    #predProbs <- predict(model$finalModel, X, lambda = lam, type = 'response')[, 1]
  } else if(is(model, 'glinternet')){
    predProbs <- predict(model, X, type = 'response')[, 2]
  } else if(is(model, 'glm')){
    predProbs <- predict(model, type = "response")
  } else if(is(model, 'numeric')){
    stopifnot(length(model) == length(y))
    predProbs <- model
  }
  p <- unname(sort(predProbs))
  if(thresh | prc){
    t1 <- (c(-Inf, p) + c(p, +Inf))/2
    t2 <- (c(-Inf, p)/2 + c(p, +Inf)/2)
    p <- ifelse(abs(t1) > 1e+100, t2, t1)
    #p <- p[-c(1, length(p))]
  }
  preds <- list()
  tp <- tn <- fp <- fn <- c()
  for(i in 1:length(p)){
    preds[[i]] <- ifelse(predProbs > p[i], 1, 0)
    Y <- cbind(y, preds[[i]])
    tn[i] <- sum(Y[Y[, 1] == 0, 1] == Y[Y[, 1] == 0, 2])
    tp[i] <- sum(Y[Y[, 1] == 1, 1] == Y[Y[, 1] == 1, 2])
    fn[i] <- sum(Y[Y[, 2] == 0, 1] != Y[Y[, 2] == 0, 2])
    fp[i] <- sum(Y[Y[, 2] == 1, 1] != Y[Y[, 2] == 1, 2])
  }
  sens <- tp/sum(y)
  spec <- tn/(length(y) - sum(y))
  npv <- tn/(tn + fn)
  ppv <- tp/(tp + fp)
  #npv <- 1 - (tp/(tp + fp))
  #ppv <- 1 - (tn/(tn + fn))
  #opt <- which.max(sens + spec)
  opt <- ifelse(isTRUE(prc), which.max(ppv + sens), which.max(sens + spec))
  optCut <- p[opt]
  optSens <- sens[opt]
  optSpec <- spec[opt]
  optPPV <- ppv[opt]
  optNPV <- npv[opt]
  if(isTRUE(prc)){
    sx <- c(1, sens)
    sy <- c(0, ppv)
    if(any(is.na(sy))){
      if(sum(is.na(sy)) == 1){
        sy <- na.omit(sy)
      } else {
        sy[which(is.na(sy))] <- sx[which(is.na(sy))]
      }
    }
    if(length(sy) < length(sx)){sy <- c(sy, 1)}
    height <- sy[-1] - sy[-length(sy)]
    width <- (sx[-1] + sx[-length(sx)])/2
  } else {
    sx <- c(0, spec)
    sy <- c(1, sens)
    height <- (sy[-1] + sy[-length(sy)])/2
    width <- sx[-1] - sx[-length(sx)]
  }
  AUC <- sum(height * width)
  if(plot == FALSE){
    return(list(results = data.frame(cutoff = p, sens = sens, spec = spec, ppv = ppv, npv = npv),
                optimal = unlist(list(cutoff = optCut, sensitivity = optSens,
                                      specificity = optSpec, PPV = optPPV,
                                      NPV = optNPV)), AUC = AUC))
  } else {
    xlabel <- ifelse(isTRUE(prc), 'Recall', '1 - Specificity')
    ylabel <- ifelse(isTRUE(prc), 'Precision', 'Sensitivity')
    main <- ifelse(isTRUE(prc), paste0('PR Curve\nAUC = ', round(AUC, 3)),
                   paste0('ROC Curve\nAUC = ', round(AUC, 3)))
    plot(0, 0, type = "n", ylim = c(0, 1), xlim = c(0, 1), axes = FALSE,
         xlab = xlabel, ylab = ylabel, main = main)
    if(grid != FALSE){
      if(grid == TRUE){grid <- 1}
      if(grid == 1){
        grid(NA, 5, lty = grid_lty, lwd = grid_lwd, col = grid_col)
      } else if(grid == 2){
        grid(lty = grid_lty, lwd = grid_lwd, col = grid_col)
      }
    }
    axis(1); axis(2)
    if(isTRUE(prc)){sx <- 1 - sx}
    lines(1 - sx, sy, lty = roc_lty, lwd = roc_lwd, col = roc_col)
    if(midline != FALSE & !isTRUE(prc)){
      if(midline == TRUE){midline <- 1}
      if(midline == "grey" | midline == "gray"){midline_lty <- 1; midline_col = "grey"}
      if(midline == 2){midline_lty <- 1; midline_col = "grey"}
      abline(a = 0, b = 1, lty = midline_lty, lwd = midline_lwd, col = midline_col)
    }
    if(optPoint != FALSE){
      if(optPoint == "black" | optPoint == 2){pt_pch <- 8; pt_col <- "black"}
      if(!pt_pch %in% c(21:25)){pt_border <- pt_col}
      points((1 - sx)[opt + 1], sy[opt + 1], pch = pt_pch, bg = pt_col, col = pt_border)
    }
  }
}


##### coefs
coefs <- function(fit, inds = FALSE, gamma = .5){
  if(!is(fit, 'train')){
    if('coefs' %in% names(fit)){return(fit$coefs)} else {return(list())}
  }
  stopifnot(is(fit, 'train'))
  lam <- fit$bestTune[, 'lambda']
  if('alpha' %in% names(fit$bestTune)){
    coef(fit$finalModel, s = lam)
  } else if(identical(inds, FALSE)){
    glint(fit = fit)$coefs
  } else {
    if(isTRUE(inds) | identical(inds, 'all')){
      inds <- c('AIC', 'BIC', 'EBIC')
    }
    out <- sapply(inds, function(z) glint(fit = fit, method = z, gamma = gamma)[[4]])
    names(out) <- toupper(names(out))
    if('EBIC' %in% names(out)){
      gamma <- gsub('^0[.]', '.', as.character(gamma))
      names(out)[which(names(out) == 'EBIC')] <- paste0('EBIC', gamma)
    }
    return(out)
  }
}


##### simpsim
simpsim <- function(n = 50, p = 3, env = FALSE, binary = TRUE, qt = .2,
                    levs = c('Class1', 'Class2'), b = NULL, intercept = FALSE,
                    mu = 0, sigma = NULL){
  if(isTRUE(mu)){mu <- rnorm(p)}
  if(length(mu) < p){mu <- rep(mu, p)}
  if(is.null(sigma)){sigma <- rcov(p)}
  if(is.null(b)){
    b <- runif(p + as.numeric(intercept), min = -.8, max = .8)
  }
  if(!is.matrix(b)){b <- matrix(b, ncol = 1)}
  dat <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
  y <- switch(2 - isTRUE(intercept), (cbind(1, dat) %*% b)[, 1], (dat %*% b)[, 1])
  if(isTRUE(binary)){
    y <- factor(ifelse(y < quantile(y, probs = qt[1]), levs[1], levs[2]))
  }
  x <- data.frame(dat)
  out <- list(x = x, y = y)
  if(isTRUE(env)){
    list2env(out, envir = .GlobalEnv)
  } else {
    return(out)
  }
}

##### capitalize: capitalize strings
capitalize <- function(x){
  unname(sapply(x, function(z){
    z1 <- substr(z, 1, 1)
    z2 <- substr(z, 2, nchar(z))
    lower <- which(letters == z1)
    if(length(lower) != 0){z1 <- LETTERS[lower]}
    return(paste0(z1, z2))
  }))
}


### ------------------------------------------------------------------------ ###
### --------------------------- CARET WRAPPERS ----------------------------- ###
### ------------------------------------------------------------------------ ###
##### ctrl: easily creates control parameters
ctrl <- function(type = 2, method = 'cv', n = 10, cvreps = 3,
                 sumfun = 'default', res = 'final',
                 preds = 'all', verbose = TRUE, ...){
  k <- n
  suppressMessages(invisible(require(caret)))
  if(is.numeric(type)){
    type <- pmin(2, pmax(1, type))
  } else if(is.character(type)){
    types <- c('regression', 'classification')
    type <- which(types == match.arg(tolower(type), types))
  }
  classProbs <- isTRUE(type == 2)
  args0 <- tryCatch({list(...)}, error = function(e){list()})
  method <- tolower(method)
  if(startsWith(method, 'l')){method <- toupper(method)}
  if(startsWith(method, 'r')){method <- 'repeatedcv'}
  args <- list(method = method, number = k, returnResamp = res, savePredictions = preds)
  if(method == 'repeatedcv'){args$repeats <- cvreps}
  if(type == 1){
    if(identical(sumfun, 'default')){sumfun <- caret::defaultSummary}
    if(is.function(sumfun)){args$summaryFunction <- sumfun}
  } else {
    if(identical(sumfun, 'default')){sumfun <- twoclass2}
    if(is.function(sumfun)){args$summaryFunction <- sumfun}
    args$classProbs <- classProbs
  }
  args <- append(args, args0[setdiff(names(args0), names(args))])
  args$verboseIter <-  isTRUE(verbose)
  do.call(caret::trainControl, args)
}


##### twoclass2
twoclass2 <- function(data, lev = NULL, model = NULL){
  if(length(lev) > 2){
    stop(paste("Your outcome has", length(lev), "levels. The twoClassSummary() function isn't appropriate."))
  }
  caret:::requireNamespaceQuietStop("pROC")
  if(!all(levels(data[, "pred"]) == lev)){
    stop("levels of observed and predicted data do not match")
  }
  suppressMessages(invisible(require('PRROC')))
  ppreds0 <- data[which(data[, "obs"] == lev[2]), lev[2]]
  ppreds1 <- data[which(data[, "obs"] == lev[1]), lev[2]]
  p <- lev[2]; n <- lev[1]
  pred <- data[, "pred"]; obs <- data[, "obs"]
  tp <- sum(pred %in% p & obs %in% p)
  tn <- sum(pred %in% n & obs %in% n)
  fp <- sum(pred %in% p & obs %in% n)
  fn <- sum(pred %in% n & obs %in% p)
  tpr <- tp/(tp + fn)
  tnr <- tn/(tn + fp)
  ppv <- tp/(tp + fp)
  ba <- (tpr + tnr)/2
  f1 <- 2 * ((ppv * tpr)/(ppv + tpr))
  acc <- (tp + tn)/(tp + tn + fp + fn)
  bottom <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
  mcc <- ((tp * tn) - (fp * fn))/sqrt(ifelse(identical(bottom, 0), 1, bottom))
  if(is.na(mcc)){mcc <- 0}
  #kappa <- kappa['kappa']
  caret:::requireNamespaceQuietStop("e1071")
  kappa <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]['kappa']
  rocObject <- try(pROC::roc(data$obs, data[, lev[2]], direction = "<",
                             quiet = TRUE), silent = TRUE)
  rocAUC <- if(inherits(rocObject, "try-error")){NA} else {rocObject$auc}
  prc <- ROCcurve(y = data[, 'obs'], model = data[, lev[2]], prc = TRUE)$AUC
  prAUC <- PRROC::pr.curve(scores.class0 = ppreds0, scores.class1 = ppreds1)[[3]]
  out <- c(rocAUC, prc, caret::sensitivity(data[, "pred"], data[, "obs"], lev[2]),
           caret::specificity(data[, "pred"], data[, "obs"], lev[1]),
           caret::precision(data[, "pred"], data[, "obs"], lev[2]),
           acc, ba, mcc, kappa, prAUC, f1)
  names(out) <- c("ROC", "PRC", "Sens", "Spec", "Prec", "Acc",
                  "BA", "MCC", "Kappoo", "prAUC", "F1")
  out
}


##### trainFit: DMwR package no longer available on CRAN
trainFit <- function(x, y, k = 'default', m = 'zzzall', subsample = 'none',
                     lams = 'default', metric = 'default', pre = NULL, grid = NULL,
                     model = 'glint', seed = NULL, time = TRUE, ...){
  t1 <- Sys.time()
  suppressMessages(invisible(require(caret)))
  args0 <- tryCatch({list(...)}, error = function(e){list()})
  if(any(sapply(m, identical, 'all'))){m <- gsub('all', 'zzzall', m)}
  if(identical(m, 'glmnet') | identical(m, 'glint')){model <- m; m <- 'zzzall'}
  #if(is.character(m)){if(m %in% c('glmnet', 'glint')){model <- m}; m <- 0}
  if(!identical(subsample, 'none') & !'sampling' %in% names(args0)){
    subsample <- match.arg(tolower(subsample), c('down', 'up', 'smote', 'rose'))
    fun <- switch(subsample, down = caret::downSample, up = caret::upSample,
                  #smote = function(x, y){DMwR::SMOTE(Class ~ ., data = data.frame(cbind(x, Class = y)))},
                  rose = function(x, y){ROSE::ROSE(Class ~ ., data = data.frame(cbind(x, Class = y)))$data})
    out <- fun(x = x, y = y)
    y <- out[, 'Class']
    x <- out[, setdiff(colnames(out), 'Class')]
  }
  if(is.null(grid) & identical(model, 'glint')){
    if(identical(lams, 'default')){lams <- 50}
    if(ifelse(length(lams) == 1, identical(round(lams), lams), FALSE)){
      lams <- lambdaGrid(x, y, lams, trim = TRUE)
    }
    grid <- expand.grid(m = m, lambda = lams, stringsAsFactors = FALSE)
  }
  glint_mod <- list(
    label = 'Hierarchical LASSO',
    library = c('glinternet', 'Matrix', 'PRROC'),
    type = c('Regression', 'Classification'),
    parameters = data.frame(parameter = c('lambda', 'm'),
                            class = c('numeric', 'character'),
                            label = c('tuning', 'interaction')),
    grid = function(x, y, len = NULL, search = 'grid'){
      numLev <- ifelse(is.character(y) | is.factor(y), dim(table(y)), NA)
      fam <- ifelse(is.na(numLev), 'gaussian', 'binomial')
      yy <- y
      if(is.factor(yy) & dim(table(yy)) == 2){
        levels(yy) <- 0:1
        yy <- as.numeric(as.character(yy))
      }
      lams <- unique(lambdaGrid(X = x, Y = yy, nLambda = len + 2, family = fam))
      lams <- lams[-c(1, length(lams))]
      lams <- lams[1:min(length(lams), len)]
      expand.grid(m = colnames(x), lambda = lams, stringsAsFactors = FALSE)
    },
    loop = NULL,
    fit = function(x, y, wts, param, lev, last, classProbs, ...){
      mod <- glint(x = x, y = y, m = param$m, lambda = lams)
      mod <- mod$fitobj$fit
      mod$lambdaOpt <- param$lambda[1]
      #mod$lambdaOpt <- mod$lambda[which.min(abs(mod$lambda - param$lambda[1]))]
      if(dim(table(y)) == 2){
        mod$obsLevels <- switch(2 - is.factor(y), levels(y), unique(y))
      }
      mod
    },
    predict = function(modelFit, newdata, submodels = NULL){
      if(!is.matrix(newdata)){newdata <- as.matrix(newdata)}
      obsLevels <- modelFit$obsLevels
      lam <- modelFit$lambdaOpt
      #lam <- modelFit$lambda[which.min(abs(modelFit$lambdaOpt - modelFit$lambda))]
      if(length(obsLevels) < 2){
        predict(modelFit, newdata, lambda = lam)[, 1]
      } else {
        preds <- predict(modelFit, newdata, lambda = lam, type = 'response')[, 1]
        ifelse(preds < .5, obsLevels[1], obsLevels[2])
      }
    },
    prob = function(modelFit, newdata, submodels = NULL){
      if(!is.matrix(newdata)){newdata <- as.matrix(newdata)}
      obsLevels <- switch(2 - ('obsLevels' %in% names(modelFit)), modelFit$obsLevels, NULL)
      #lam <- modelFit$lambdaOpt
      #modelFit <- glint(x = newdata, fit = modelFit)$fit0
      #lam <- modelFit$lambda[which.min(abs(modelFit$lambdaOpt - modelFit$lambda))]
      lam <- modelFit$lambdaOpt
      out <- predict(modelFit, newdata, lambda = lam, type = 'response')
      if(length(obsLevels) == 2){
        out <- out[, 1]
        out <- as.data.frame(cbind(1 - out, out), stringsAsFactors = FALSE)
        colnames(out) <- obsLevels
      } else {
        out <- as.data.frame(out[, , 1, drop = FALSE], stringsAsFactors = FALSE)
        names(out) <- obsLevels
      }
      out
    },
    levels = function(x){if(any(names(x) == 'obsLevels')){x$obsLevels} else {NULL}},
    sort = function(x){x[order(-x$lambda, x$m), ]}
  )
  if(ifelse(is(k, 'list'), all(names(k) %in% formalArgs(ctrl)), FALSE)){
    args0[names(k)] <- k
    k <- 'default'
  }
  if(identical(k, 'default')){
    a1 <- list(type = ifelse(dim(table(y)) == 2, 2, 1))
    a0 <- intersect(names(args0), c(formalArgs(ctrl), formalArgs(caret::trainControl)))
    if(length(a0) > 0){a1 <- append(a1, args0[a0])}
    k <- do.call(ctrl, a1)
  }
  if('method' %in% names(args0)){args0 <- args0[setdiff(names(args0), 'method')]}
  if(identical(metric, 'default')){metric <- ifelse(dim(table(y)) == 2, 'ROC', 'RMSE')}
  model <- switch(2 - identical(model, 'glint'), glint_mod, 'glmnet')
  args <- list(x = x, y = y, method = model, trControl = k, metric = metric)
  if(!is.null(pre)){args$preProcess <- pre}
  if(identical(model, 'glmnet')){
    if(length(lams) == 1){
      nlam <- ifelse(identical(lams, 'default'), 101, lams)
      lams <- glmnet::glmnet(x = as.matrix(x), y = y, alpha = 1,
                             family = ifelse(dim(table(y)) == 2, 'binomial', 'gaussian'),
                             nlambda = nlam + 2)$lambda
      lams <- unique(lams)
      lams <- lams[-c(1, length(lams))]
      lams <- lams[1:min(length(lams), nlam)]
    }
    args$tuneGrid <- expand.grid(alpha = 1, lambda = lams)
  } else {
    args$tuneGrid <- grid
  }
  if(length(args0) > 0){
    args <- append(args, args0[intersect(names(args0), formalArgs(train))])
  }
  if('tuneLength' %in% names(args) & identical(model, 'glmnet')){args$tuneGrid <- NULL}
  output <- suppressWarnings(do.call(caret::train, args))
  if(isTRUE(time)){print(t2 <- Sys.time() - t1)}
  return(output)
}


##### test_roc
test_roc <- function(fit, x = NULL, y = NULL, v = TRUE, ci = TRUE){
  suppressMessages(invisible(library(pROC)))
  if(!isTRUE(ci)){v <- FALSE}
  if(!is.null(x)){if(is.logical(x)){v <- x; x <- NULL}}
  if(!is.null(y)){if(is.logical(y)){v <- y; y <- NULL}}
  if(is(fit, 'list')){
    if(is.null(names(fit))){names(fit) <- paste0('fit', 1:length(fit))}
    args <- as.list(match.call())[-1]
    out <- lapply(fit, function(z){
      args0 <- replace(args, 'fit', list(fit = z))
      do.call(test_roc, args0)
    })
    nn <- names(out)
    out <- setNames(do.call(rbind.data.frame, lapply(out, as.vector)), c('lower', 'ROC', 'upper'))
    rownames(out) <- nn
    return(out)
  }
  if(is.null(x) | is.null(y)){
    xx <- fit$trainingData
    if(is.null(x)){x <- xx[, setdiff(colnames(xx), '.outcome')]}
    if(is.null(y)){y <- xx[, '.outcome']}
  }
  levs <- levels(y)
  roc_obj <- roc(y, predict(fit, x, type = "prob")[, levs[2]],
                 levels = levs)
  if(isTRUE(ci) | isTRUE(v)){roc_obj <- ci(roc_obj)}
  if(isTRUE(v)){roc_obj <- setNames(as.vector(roc_obj), c('lower', 'ROC', 'upper'))}
  roc_obj
}


##### classperf
classperf <- function(pred, obs, metric = 'sens', p = levels(obs)[2]){
  stopifnot(is.factor(pred) & is.factor(obs))
  stopifnot(length(unique(c(levels(obs), levels(pred)))) == 2)
  n <- setdiff(levels(obs), p)
  metric <- match.arg(tolower(metric), c(
    'sensitivity', 'specificity', 'precision', 'ppv', 'kappa', 'f1',
    'fdr', 'mcc', 'accuracy', 'ba'), several.ok = TRUE)
  if('ppv' %in% metric){metric[metric == 'ppv'] <- 'precision'}
  metric <- capitalize(unique(metric))
  if(any(nchar(metric) < 4)){
    metric[nchar(metric) < 4] <- toupper(metric[nchar(metric) < 4])
  }
  caret:::requireNamespaceQuietStop("e1071")
  kappa <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]
  kappa <- kappa['kappa']
  tp <- sum(pred %in% p & obs %in% p)
  tn <- sum(pred %in% n & obs %in% n)
  fp <- sum(pred %in% p & obs %in% n)
  fn <- sum(pred %in% n & obs %in% p)
  tpr <- tp/(tp + fn)
  tnr <- tn/(tn + fp)
  ppv <- tp/(tp + fp)
  fdr <- 1 - ppv
  ba <- (tpr + tnr)/2
  acc <- (tp + tn)/(tp + tn + fp + fn)
  f1 <- 2 * ((ppv * tpr)/(ppv + tpr))
  mcc <- suppressWarnings(((tp * tn) - (fp * fn))/sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
  output <- list(Sensitivity = tpr, Specificity = tnr, Precision = ppv, FDR = fdr,
                 Accuracy = acc, MCC = mcc, BA = ba, Kappa = kappa, F1 = f1)
  return(unlist(output[names(output) %in% metric]))
}


##### sums
sums <- function(res, metric = 'default', means = TRUE){
  if(!is(res, 'resamples') & is(res, 'list')){res <- caret::resamples(res)}
  mets <- res$metrics
  mm <- ifelse(is.character(means), ifelse(
    means %in% mets, means, 'default'), 'default')
  if(any(sapply(c(metric, means), identical, 'all'))){metric <- 'default'; means <- 'all'}
  if(is.logical(metric)){means <- metric; metric <- mm}
  if(identical(metric, 'default')){
    metric <- ifelse('MAE' %in% res$metrics, 'RMSE', 'ROC')
  }
  s <- summary(res, metric = metric)
  if(isTRUE(means)){
    s <- s$statistics[metric]
    if(length(s) == 1){s <- s[[1]]}
    return(s)
  } else if(identical(means, FALSE)){
    if(length(metric) > 1){
      out <- setNames(lapply(metric, function(i){
        z1 <- s$values[, which(gsub('.*.~', '', colnames(s$values)) == i)]
        names(z1) <- gsub('~.*.', '', names(z1))
        z1
      }), metric)
    } else {
      out <- s$values[, which(gsub('.*.~', '', colnames(s$values)) == metric)]
      names(out) <- gsub('~.*.', '', names(out))
    }
    return(out)
  } else {
    return(s$values)
  }
}


##### LL: log-likelihood
LL <- function(fit, y, x = NULL, int = FALSE){
  if(is.factor(y)){
    levels(y) <- 0:1
    y <- as.numeric(as.character(y))
  }
  if(is.null(x)){
    prob <- predict(fit, type = 'response')
  } else {
    prob <- predict(fit, as.matrix(x), type = 'response')
    if(ncol(prob) > 1){prob <- prob[, 2]}
  }
  k <- fit$rank
  if(is.null(k) & is(fit, 'glinternet')){
    k1 <- length(unlist(coef(fit)[[2]][[1]]))
    k2 <- sum(unlist(sapply(coef(fit)[[2]][[3]], nrow)))
    k <- k1 + k2 + as.numeric(int)
  }
  ll <- sum(log(prob * y + (1 - prob) * (1 - y)))
  aic <- -2 * ll + 2 * k
  bic <- -2 * ll + k * log(length(y))
  c(ll = ll, AIC = aic, BIC = bic)
}


##### getMod
getMod <- function(fit, metric = 'PRC', row = NULL, crit = 'AIC'){
  if(!is(fit, 'train') & is(fit, 'list')){
    call <- as.list(match.call())[-1]
    return(lapply(fit, function(z){
      args0 <- replace(call, 'fit', list(fit = z))
      do.call(getMod, args0)
    }))
  }
  stopifnot(is(fit, 'train'))
  parms <- c('lambda', 'm')
  k <- fit$results[, -grep('SD$', colnames(fit$results))]
  k$F1 <- 2 * ((k$Prec * k$Sens)/(k$Prec + k$Sens))
  mets <- setdiff(tolower(colnames(k)), parms)
  metric <- capitalize(match.arg(tolower(metric), mets))
  if(metric %in% c('Roc', 'Prc', 'Ba', 'Mcc')){metric <- toupper(metric)}
  if(metric == 'Prauc'){metric <- 'prAUC'}
  mx <- max(k[, metric], na.rm = TRUE)
  kk <- k[k[, metric] == mx, ]
  if(any(is.na(kk))){
    k0 <- which(apply(kk, 1, function(z) sum(is.na(z))) == ncol(kk))
    if(length(k0) > 0){kk <- kk[-k0, ]}
  }
  if(nrow(kk) > 1 & (is.null(row) | identical(row, FALSE))){
    k2 <- do.call(rbind, lapply(1:nrow(kk), function(i){
      ii <- getMod(fit = fit, metric = metric, row = i)[[1]]
      ints <- ii[grepl(':', ii)]
      mains <- setdiff(ii, ints)
      c(vars = length(mains), ints = length(ints))
    }))
    return(cbind(k2, kk))
  } else if(identical(row, 'auto') & nrow(kk) > 1){
    k1 <- which.min(sapply(1:nrow(kk), function(z) getMod(fit, metric, row = z)[[4]]))
    kk <- kk[k1, , drop = FALSE]
  } else if(!is.null(row) & nrow(kk) > 1){
    kk <- kk[row, , drop = FALSE]
  }
  attr(fit$finalModel, 'm') <- kk[, parms[2]]
  out <- glint(fit, method = crit, lambda = kk[, parms[1]], m = kk[, parms[2]])
  out$parms <- kk[, parms]
  return(out)
}
