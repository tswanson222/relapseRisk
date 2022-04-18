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
    form <- as.formula(paste('~', paste0(names(coefs)[-1], collapse = ' + ')))
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
  if(pilr | all(c('survey_code', 'event_type', 'question_code') %in% colnames(data))){
    data <- lapply(c('epsi', 'idas'), function(z) pilrdata(data, z, ...))
    nas <- c(data[[1]]$non_response,data[[2]]$non_response)
    if(all(nas)){
      out <- data.frame(risk = 'MISSING', predprob = NA, Time = NA)
      return(out)
    } else if(any(nas)){
      data <- data[[-which(nas)]]$output
    } else {
      data <- do.call(rbind, lapply(data, '[[', 'output'))
    }
    epsi <- c('Bulimia', 'Exercise-Focused Behaviors', 'Restrictive Eating')
    idas <- c('Fear', 'Distress', 'Positive Affect')
    data <- split(data, data$Time)
    for(i in names(data)){
      if(all(c(epsi, idas) %in% data[[i]]$Thetas)){
        data[[i]]$Thetas <- factor(data[[i]]$Thetas, levels = c(epsi, idas),
                                   labels = c(paste0('EP', 1:3), paste0('I', 1:3)))
        temp <- structure(data.frame(t(data.frame(data[[i]]$Estimate))), row.names = '1')
        colnames(temp) <- data[[i]]$Thetas
        data[[i]] <- predFun(dat = temp, coefs = coefs, cutoff = cutoff, se = se)
        data[[i]] <- data.frame(data[[i]], Time = as.numeric(i))
      } else {
        data[[i]] <- data.frame(risk = 'MISSING', predprob = NA, Time = as.numeric(i))
      }
    }
    out <- data.frame(do.call(rbind, data))
  } else {
    out <- predFun(dat = data, coefs = coefs, cutoff = cutoff, se = se)
  }
  return(out)
}

