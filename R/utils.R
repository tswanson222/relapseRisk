### ------------------------------------------------------------------------ ###
### -------------------------- ADAPTED FUNCTIONS --------------------------- ###
### ------------------------------------------------------------------------ ###
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

##### getLambdaGrid: from glinternet
getLambdaGrid <- function(candidates, nlam, minRatio){
  lmax <- max(sapply(candidates$norms, function(z) ifelse(is.null(z), 0, max(z))))
  lmin <- minRatio * lmax
  f <- seq(0, 1, 1/(nlam - 1))
  lambda <- lmax^(1 - f) * lmin^f
  return(lambda)
}

##### standard: from glinternet
standard <- function(x){
  if(length(unique(x)) == 1){
    out <- rep(0, length(x))
  } else {
    out <- x - mean(x)
    out <- out/rep(sqrt(t(out) %*% out), length(out))
  }
  return(out)
}

##### rcov: simulate a random covariance matrix
rcov <- function(n, dist = "runif", method = 2, args = NULL){
  inner <- paste0(dist, "(", n^2, ifelse(is.null(args), "", paste0(", ", args)), ")")
  inner <- eval(parse(text = inner))
  if(method %in% c(3, 4)){inner <- inner * 2 - 1}
  S <- matrix(inner, n, n)
  if(method %in% c(1, 3)){
    S <- t(S) %*% S
  } else {
    S <- qr.Q(qr(S))
    S <- crossprod(S, S * (n:1))
  }
  return(S)
}
