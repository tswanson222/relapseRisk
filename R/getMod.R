#' Get the final model(s) from \code{trainFit}
#'
#' Description
#'
#' @param fit Output of \code{\link{trainFit}}
#' @param metric Character string for the performance metric (LIST AVAILABLE)
#' @param row Only relevant if there are multiple models that maximize the
#'   performance metric. Allows the user to choose which model (based on the row
#'   number) to return information for.
#' @param crit Criterion for choosing the final model (unsure)
#'
#' @return TBD
#' @export
#'
#' @examples
#' 1 + 1
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

#' Return coefficient values from \code{trainFit} model
#'
#' Description
#'
#' @param fit Output of \code{\link{trainFit}}
#' @param inds TBD
#' @param gamma Hyperparameter gamma value for EBIC.
#'
#' @return Coefficient values associated with final model
#' @export
#'
#' @examples
#' 1 + 1
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
