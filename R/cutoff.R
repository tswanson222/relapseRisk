#' Extract cutoff mean and SD across CV folds
#'
#' @param fit output from trainFit
#' @param lambda set to "best" to use optimal lambda
#'
#' @return PRC cutoff mean and SD taken across folds for specific lambda value
#' @export
#'
#' @examples
#' 1 + 1
cutoff <- function(fit, lambda = 'best'){
  if(is(fit, 'train')){
    if(identical(lambda, 'best')){lambda <- fit$bestTune[, 'lambda']}
    out <- unlist(fit$results[which(fit$results$lambda == lambda), c('cutoff', 'cutoffSD')])
    out[3] <- out[2]/sqrt(nrow(fit$resample))
    names(out)[3] <- 'cutoffSE'
  } else {
    out <- c(cutoff = mean(fit, na.rm = TRUE),
             cutoffSD = sd(fit, na.rm = TRUE),
             cutoffSE = sd(fit, na.rm = TRUE)/sqrt(length(fit)))
  }
  out <- c(out, upper_1SD = unname(out[1] + out[2]),
           lower_1SD = unname(out[1] - out[2]),
           upper_1SE = unname(out[1] + out[3]),
           lower_1SE = unname(out[1] - out[3]))
  return(out)
}


#' Resample data to estimate cutoff sampling distribution
#'
#' @param formula formula for selected model
#' @param data training dataset
#' @param niter number of iterations to resample
#' @param seed numeric value of seed
#' @param prc logical, determines whether to find PRC or ROC cutoffs
#' @param method resampling method
#'
#' @return vector of cutoff values
#' @export
#'
#' @examples
#' 1 + 1
sampleCutoffs <- function(formula, data, niter = 1000, seed = NULL, prc = TRUE,
                          method = c('boot', 'up', 'down', 'rose', 'smote')){
  method <- match.arg(method)
  stopifnot('y' %in% colnames(data))
  if(!is.null(seed)){set.seed(seed)}
  sampFun <- switch(method,
                    boot = function(data){
                      inds <- sample(1:nrow(data), replace = TRUE)
                      out <- data[inds, ]
                      return(out)
                    },
                    up = function(data){
                      out <- caret::upSample(x = data[, setdiff(colnames(data), 'y')], y = data$y)
                      colnames(out)[which(colnames(out) == 'Class')] <- 'y'
                      return(out)
                    },
                    down = function(data){
                      out <- caret::downSample(x = data[, setdiff(colnames(data), 'y')], y = data$y)
                      colnames(out)[which(colnames(out) == 'Class')] <- 'y'
                      return(out)
                    },
                    rose = function(data){
                      colnames(data)[which(colnames(data) == 'y')] <- 'Class'
                      out <- ROSE::ROSE(Class ~ ., data)$data
                      colnames(out)[which(colnames(out) == 'Class')] <- 'y'
                      return(out)
                    },
                    smote = function(data){
                      data <- data.frame(data[, setdiff(colnames(data), 'y')],
                                         Class = data$y)
                      out <- themis::smote(data, 'Class')
                      colnames(out)[which(colnames(out) == 'Class')] <- 'y'
                      return(out)
                    })
  out <- sapply(replicate(niter, list(sampFun(data))), function(z){
    mod <- suppressWarnings(glm(formula, z, family = binomial))
    cutoff <- ROCcurve(mod, prc = prc)$optimal['cutoff']
    if(cutoff == -Inf){cutoff <- NA}
    return(cutoff)
  })
  return(out)
}

