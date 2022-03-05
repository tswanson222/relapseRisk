#' Reshaping glinternet output using criterion-based selection
#'
#' Description. May not even be needed as global function.
#'
#' @param x TBD
#' @param y TBD
#' @param method TBD
#' @param m TBD
#' @param type TBD
#' @param gamma TBD
#' @param nlam TBD
#' @param useSE TBD
#' @param nfolds TBD
#' @param allCoef TBD
#' @param fit TBD
#' @param lambda TBD
#' @param intOnly TBD
#' @param lamNum TBD
#' @param ... Additional arguments
#'
#' @return TBD
#' @export
#'
#' @examples
#' 1 + 1
glint <- function(x, y = 1, method = 'AIC', m = 'all', type = 'default',
                  gamma = .5, nlam = 50, useSE = TRUE, nfolds = 10,
                  allCoef = FALSE, fit = NULL, lambda = NULL,
                  intOnly = FALSE, lamNum = FALSE, ...){
  if(!missing(x)){if(is(x, 'train')){fit <- x; x <- NULL}}
  if(length(method) > 1){fit <- method; method <- 'aic'}
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
      if(identical(m, 'all')){m <- fit$bestTune[, 'm']}
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
  if(identical(m, 'all') | identical(m, 0) | identical(m, ncol(x) + 1)){m <- 1:ncol(x)}
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
      glinternet::glinternet.cv(x, y, type, interactionCandidates = m,
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
      glinternet::glinternet(x, y, type, lambda = fitCV$lambdaHat,
                             interactionCandidates = m, family = family, ...)
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
      glinternet::glinternet(x, y, type, lambda = fitCV$lambdaHat,
                             interactionCandidates = m, family = family, ...)
    )
    fit1se <- dots(
      glinternet::glinternet(x, y, type, lambda = fitCV$lambdaHat1Std,
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
    coefs <- Matrix::Matrix(cbind(mod0, mod1se), sparse = TRUE)
    allNames <- c(colnames(x), apply(combn(colnames(x), 2), 2, paste, collapse = ":"))
    rownames(coefs) <- c('(Intercept)', allNames)
    output <- list(mod0 = allNames[mod0 != 0], mod1se = allNames[mod1se != 0], coefs = coefs,
                   fitobj = list(fitCV = fitCV, fit0 = fit0, fit1se = fit1se))
    output$mod0 <- setdiff(output$mod0, '(Intercept)')
    output$mod1se <- setdiff(output$mod1se, '(Intercept)')
  } else {
    criterion <- toupper(method)
    fit <- dots(
      glinternet::glinternet(x, y, type, interactionCandidates = m, family = family,
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
    coefs <- Matrix::Matrix(betas, sparse = TRUE)
    rownames(coefs) <- names(betas)
    fitobj <- list(fit = fit, fit0 = NA, crit = ic_lambda)
    if(length(fit$lambda) > 2){
      fitobj$fit0  <- dots(
        glinternet::glinternet(x, y, type, interactionCandidates = m,
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


##### lambdaGrid: Just returns lambda values from glinternet. Global?
lambdaGrid <- function(X, Y, nLambda = 50, m = NULL, trim = FALSE, lambda = NULL, numLevels,
                       lambdaMinRatio = 0.01, interactionPairs = NULL, screenLimit = NULL,
                       numToFind = NULL, family = c("gaussian", "binomial"),
                       tol = 1e-5, maxIter = 5000, verbose = FALSE,
                       numCores = 1, lamOnly = TRUE){
  # get call and family
  interactionCandidates <- m
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
    Z = as.matrix(apply(as.matrix(X[, contIndices]), 2, standard))
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
    lambda = getLambdaGrid(candidates, nLambda, lambdaMinRatio)
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
    lambdaMax = max(getLambdaGrid(candidates, nLambda, lambdaMinRatio))
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
    activeSet[[i]] = glinternet:::strong_rules(candidates, lambda[i], lambda[i-1])
    betahat[[i]] = glinternet:::initialize_betahat(activeSet[[i]], activeSet[[i-1]], betahat[[i-1]], levels)
    while (TRUE) {
      # group lasso on strong set
      solution = glinternet:::group_lasso(Xcat, Z, Y, activeSet[[i]], betahat[[i]], levels, lambda[i], family, tol, maxIter, verbose)
      activeSet[[i]] = solution$activeSet
      betahat[[i]] = solution$betahat
      res = solution$res
      objValue[i] = solution$objValue
      # check kkt conditions on the rest
      check = glinternet:::check_kkt(Xcat, Z, res, n, pCat, pCont, levels, candidates, activeSet[[i]], lambda[i], numCores)
      candidates$norms = check$norms
      if (check$flag) {
        break
      }
      betahat[[i]] = glinternet:::initialize_betahat(check$activeSet, activeSet[[i]], betahat[[i]], levels)
      activeSet[[i]] = check$activeSet
    }
    # update the candidate set if necessary
    if (!is.null(screenLimit) && (screenLimit<pCat+pCont) && i<nLambda) {
      candidates = glinternet:::get_candidates(Xcat, Z, res, n, pCat, pCont, levels, interactionPairs, categoricalCandidates, continuousCandidates, screenLimit, activeSet[[i]], candidates$norms, numCores)
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
  betahatRescaled = lapply(1:i, function(j) glinternet:::rescale_betahat(activeSet[[j]], betahat[[j]], Xcat, Z, levels, n))

  output = list(call=thisCall, fitted=fitted[, 1:i], lambda=lambda[1:i], objValue=objValue, activeSet=activeSet[1:i], betahat=betahatRescaled[1:i], numLevels=numLevels, family=family)
  class(output) = "glinternet"

  return (output)
}

