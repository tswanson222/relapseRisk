#' Simulate class imbalance data
#'
#' Description
#'
#' @param n TBD
#' @param p TBD
#' @param env TBD
#' @param binary TBD
#' @param qt TBD
#' @param levs TBD
#' @param b TBD
#' @param intercept TBD
#' @param mu TBD
#' @param sigma TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#' x <- simpsim()
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
