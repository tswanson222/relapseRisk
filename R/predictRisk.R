#' Predict risk of relapse from PiLR data
#'
#' @param data PiLR survey data file
#' @param ... Additional arguments for pilrdata function
#'
#' @return Risk level associated with each timepoint
#' @export
#'
#' @examples
#' 1 + 1
predictRisk <- function(data, ...){
  data <- lapply(c('epsi', 'idas'), function(z) pilrdata(data, z, ...))
  data <- do.call(rbind, lapply(data, '[[', 'output'))
  epsi <- c('Bulimia', 'Exercise-Focused Behaviors', 'Restrictive Eating')
  idas <- c('Fear', 'Distress', 'Positive Affect')
  stopifnot(all(c(epsi, idas) %in% data$Thetas))
  data <- split(data, data$Time)
  form <- as.formula(paste0('~', paste0(names(trainModel)[-1], collapse = ' + ')))
  for(i in 1:length(data)){
    data[[i]]$Thetas <- factor(data[[i]]$Thetas, levels = c(epsi, idas),
                               labels = c(paste0('EP', 1:3), paste0('I', 1:3)))
    temp <- structure(data.frame(t(data.frame(data[[i]]$Estimate))), row.names = '1')
    colnames(temp) <- data[[i]]$Thetas
    data[[i]] <- model.matrix(form, temp)
    b <- (data[[i]] %*% trainModel)[, 1]
    data[[i]] <- ifelse(exp(b)/(1 + exp(b)) > prc_cutoff, 'High', 'Low')
  }
  out <- do.call(rbind, data)
  out <- data.frame(outcome = out[, 1], Time = 1:nrow(out))
  return(out)
}
