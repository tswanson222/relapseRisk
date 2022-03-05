### Train model for predicting risk of relapse

# 1) Load internal data (only needed for last step)
devtools::document()


# 2) Load training data
data <- readRDS('~/../Desktop/LINUX_2022/Documents/CARE/CATapp/CATapp/Files/epsidas_FINAL.RDS')
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
usethis::use_data(epsi_idas_questions, final_coefs, prc_cutoff,
                  internal = TRUE, overwrite = TRUE)
