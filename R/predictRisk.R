#' Predict risk of relapse from PiLR data
#'
#' @param data PiLR survey data file
#' @param cutoff Leave as "default" to use trained cutoff values. Or provide a
#'   single value for a different cutoff
#' @param coefs Leave as "default" to use trained coefficient values. Or provide
#'   a named vector of coefficients.
#' @param pilr Logical. Determines whether the input is from PiLR
#' @param se Logical. Detemines whether to use the standard error or standard
#'   deviation of the cutoff values.
#' @param ... Additional arguments for pilrdata function
#'
#' @return Risk level associated with each timepoint
#' @export
#'
#' @examples
#' 1 + 1
predictRisk <- function(data, cutoff = 'default', coefs = 'default',
                        pilr = FALSE, se = FALSE, ...){
  if(identical(cutoff, 'default')){cutoff <- prc_cutoff}
  if(identical(coefs, 'default')){coefs <- final_coefs}
  predFun <- function(dat, coefs, cutoff, se){
    form <- as.formula(paste0('~', paste0(names(coefs)[-1], collapse = ' + ')))
    dat <- model.matrix(form, dat)
    b <- (dat %*% coefs)[, 1]
    predprob <- exp(b)/(1 + exp(b))
    if(length(cutoff) == 1){
      outcome <- ifelse(predprob < cutoff, 'LOW', 'HIGH')
    } else {
      upper <- paste0('upper_1', ifelse(se, 'SE', 'SD'))
      lower <- paste0('lower_1', ifelse(se, 'SE', 'SD'))
      if(cutoff[lower] < 0){stop('Lower bound of cutoff goes below 0')}
      outcome <- ifelse(predprob > cutoff[upper], 'HIGH',
                        ifelse(predprob < cutoff[lower], 'LOW', 'MODERATE'))
    }
    out <- data.frame(risk = outcome, predprob = predprob)
    return(out)
  }
  if(isTRUE(pilr)){
    data <- lapply(c('epsi', 'idas'), function(z) pilrdata(data, z))
    data <- do.call(rbind, lapply(data, '[[', 'output'))
    epsi <- c('Bulimia', 'Exercise-Focused Behaviors', 'Restrictive Eating')
    idas <- c('Fear', 'Distress', 'Positive Affect')
    stopifnot(all(c(epsi, idas) %in% data$Thetas))
    data <- split(data, data$Time)
    for(i in 1:length(data)){
      data[[i]]$Thetas <- factor(data[[i]]$Thetas, levels = c(epsi, idas),
                                 labels = c(paste0('EP', 1:3), paste0('I', 1:3)))
      temp <- structure(data.frame(t(data.frame(data[[i]]$Estimate))), row.names = '1')
      colnames(temp) <- data[[i]]$Thetas
      data[[i]] <- predFun(dat = temp, coefs = coefs, cutoff = cutoff, se = se)
    }
    out <- do.call(rbind, data)
    out <- data.frame(out, Time = 1:nrow(out))
  } else {
    out <- predFun(dat = data, coefs = coefs, cutoff = cutoff, se = se)
  }
  return(out)
}

