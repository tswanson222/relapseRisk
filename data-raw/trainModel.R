### Train model for predicting risk of relapse

# 1) Load internal data (only needed for last step)
devtools::document()


# 2) Load training data
data <- readRDS('~/../Desktop/LINUX_2022/Documents/CARE/CATapp/CATapp/Files/epsidas_FINAL.RDS')
x <- data[, -1]
y <- data$y


# 3) Find optimal hierarchical LASSO model based on area under the PRC
samp_methods <- c('none', 'up', 'down', 'smote', 'rose')
fit1 <- setNames(lapply(samp_methods, function(z){
  relapseRisk::trainFit(x, y, sampling = switch(z, none = NULL, z), pre = 'center',
                        seed = 1, metric = 'PRC', method = 'cv', n = 10)
}), samp_methods)


fit0 <- relapseRisk::trainFit(x, y, sampling = 'rose', pre = 'center', seed = 1,
                              metric = 'PRC', method = 'cv', n = 10)


b <- as.matrix(coefs(fit0))
b <- setdiff(names(b[b[, 1] != 0, ]), '(Intercept)')
form <- paste('Class ~', paste0(b, collapse = ' + '))

set.seed(1)
seeds <- sample(1:1000000, 1000)
fit2 <- sapply(seeds, function(seed){
  set.seed(seed)
  #out <- ROSE::ROSE(y ~ ., data)$data
  out <- themis::smote(data.frame(x, Class = y), 'Class')
  model <- glm(form, data = out, family = binomial)
  out2 <- ROCcurve(model, prc = TRUE)$optimal['cutoff']
  if(any(out2 == -Inf)){out2 <- out2[-which(out2 == -Inf)]}
  return(out2)
})




################################################################################
# 4) Refit model using selected variables with OLS
b <- as.matrix(coefs(fit0))
b <- setdiff(names(b[b[, 1] != 0, ]), '(Intercept)')
form <- paste('y ~', paste0(b, collapse = ' + '))
trainModel <- glm(form, data = data, family = binomial)


# 5) Find probability threshold for classification using the PRC
prc_cutoff <- relapseRisk::ROCcurve(trainModel, prc = TRUE)$optimal['cutoff']


# 6) Pull final model coefficients
trainModel <- coef(trainModel)


# 7) Save results as internal data
usethis::use_data(epsi_idas_questions, trainModel, prc_cutoff, internal = TRUE, overwrite = TRUE)
