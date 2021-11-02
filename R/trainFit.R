##### trainFit: DMwR package no longer available on CRAN
#' Algorithm for iterated resampling and model fitting
#'
#' Description
#'
#' @param x TBD
#' @param y TBD
#' @param k TBD
#' @param m TBD
#' @param subsample TBD
#' @param lams TBD
#' @param metric TBD
#' @param pre TBD
#' @param grid TBD
#' @param model TBD
#' @param seed TBD
#' @param time TBD
#' @param ... Additional arguments
#'
#' @return TBD
#' @export
#'
#' @examples
#' 1 + 1
trainFit <- function(x, y, k = 'default', m = 'zzzall', subsample = 'none',
                     lams = 'default', metric = 'default', pre = NULL, grid = NULL,
                     model = 'glint', seed = NULL, time = TRUE, ...){
  t1 <- Sys.time()
  #suppressMessages(invisible(require(caret)))
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
    args <- append(args, args0[intersect(names(args0), formalArgs(caret::train))])
  }
  if('tuneLength' %in% names(args) & identical(model, 'glmnet')){args$tuneGrid <- NULL}
  output <- suppressWarnings(do.call(caret::train, args))
  if(isTRUE(time)){print(t2 <- Sys.time() - t1)}
  return(output)
}

#' Create control parameters for \code{trainFit}
#'
#' Description
#'
#' @param type TBD
#' @param method TBD
#' @param n TBD
#' @param cvreps TBD
#' @param sumfun TBD
#' @param res TBD
#' @param preds TBD
#' @param verbose TBD
#' @param ... Additional arguments
#'
#' @return TBD
#' @export
#'
#' @examples
#' ctrl()
ctrl <- function(type = 2, method = 'cv', n = 10, cvreps = 3,
                 sumfun = 'default', res = 'final',
                 preds = 'all', verbose = TRUE, ...){
  k <- n
  #suppressMessages(invisible(require(caret)))
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
