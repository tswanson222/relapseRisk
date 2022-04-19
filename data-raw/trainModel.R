### Train model for predicting risk of relapse

# 1) Load internal data (only needed for last step)
devtools::document()

data <- readRDS('Z:/Current Studies/STAR Grant/Data/TrainingData/trainingData_2022-March-9.RDS')
data <- data[,-1]

# 2) Compute the 90% and 95% quantiles
quantile90 <- apply(data,2,FUN=function(x) quantile(x,probs=0.9))
quantile95 <- apply(data,2,FUN=function(x) quantile(x,probs=0.95))

# 3) Compute the 1, 1.5, and 2 z-scores
means <- apply(data,2,mean)
sds <- apply(data,2,sd)

compute_zscore <- function(data,cutpoint=1){
  result <- rep(NA,dim(data)[2])
  for(i in 1:dim(data)[2]){
    result[i] <- cutpoint*sd(data[,i])+mean(data[,i])
  }
  result
}

zscore1 <- compute_zscore(data,1)
zscore1.5 <- compute_zscore(data,1.5)
zscore2 <- compute_zscore(data,2)


# 2) Load training data
data <- readRDS('Z:/Current Studies/STAR Grant/Data/TrainingData/trainingData_2022-March-9.RDS')
x <- data[, -1]
y <- data$y


# 3) Perform variable selection by cross-validating the hierarchical LASSO
fit <- trainFit(x, y, pre = 'center', seed = 1, metric = 'PRC', sampling = 'rose',
                method = 'repeatedcv', cvreps = 10, n = 10)


# 4) Pull final model with selected variables, and fit to training data
fit_coefs <- as.matrix(coefs(fit))
fit_coefs <- fit_coefs[fit_coefs != 0, 1]
form <- as.formula(paste('y ~', paste0(names(fit_coefs)[-1], collapse = ' + ')))
model <- glm(form, data, family = binomial)
final_coefs <- coef(model)


# 5) Refit final model to resamples of the data to estimate cutoff distribution
cut_values <- sampleCutoffs(form, data, niter = 2000, seed = 123, method = 'rose')


# 6) Pull summary statistics of cutoff sampling distribution
prc_cutoff <- cutoff(cut_values)


# 7) Save results as internal data
usethis::use_data(quantile90,quantile95,zscore1,zscore1.5,zscore2,
                  epsi_idas_questions, final_coefs, prc_cutoff,
                  internal = TRUE, overwrite = TRUE)
