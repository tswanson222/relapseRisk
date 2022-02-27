#' Log likelihood
#'
#' @param fit TBA
#' @param y tba
#' @param x tba
#' @param int TBA
#'
#' @return TBA
#' @export
#'
#' @examples
#' 1 + 1
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
