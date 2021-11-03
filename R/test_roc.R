#' Test ROC results
#'
#' Description
#'
#' @param fit TBA
#' @param x TBA
#' @param y TBA
#' @param v TBA
#' @param ci TBA
#'
#' @return TBA
#' @export
#'
#' @examples
#' 1 + 1
test_roc <- function(fit, x = NULL, y = NULL, v = TRUE, ci = TRUE){
  #suppressMessages(invisible(library(pROC)))
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
  roc_obj <- pROC::roc(y, predict(fit, x, type = "prob")[, levs[2]], levels = levs)
  if(isTRUE(ci) | isTRUE(v)){roc_obj <- pROC::ci(roc_obj)}
  if(isTRUE(v)){roc_obj <- setNames(as.vector(roc_obj), c('lower', 'ROC', 'upper'))}
  roc_obj
}


#' Summarize results across multiple \code{trainFit} modls
#'
#' @param res List of \code{\link{trainFit}} outputs, or output from
#'   \code{resamples} function from the \code{caret} package.
#' @param metric Character string or vector indicating which performance metrics
#'   to summarize
#' @param means Logical. Determines whether to return means for the performance
#'   metrics or not.
#'
#' @return TBD
#' @export
#'
#' @examples
#' 1 + 1
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
