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
  if(isTRUE(ci) | isTRUE(v)){roc_obj <- ci(roc_obj)}
  if(isTRUE(v)){roc_obj <- setNames(as.vector(roc_obj), c('lower', 'ROC', 'upper'))}
  roc_obj
}
