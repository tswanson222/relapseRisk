## code to prepare `trainModel` dataset goes here

devtools::document()
data <- readRDS('~/../Desktop/LINUX_2022/Documents/CARE/CATapp/CATapp/Files/epsidas_FINAL.RDS')
x <- data[, -1]
y <- data$y

fit <- relapseRisk::trainFit(x, y, subsample = 'rose', pre = 'center', seed = 123,
                             metric = 'PRC', method = 'cv', n = 10)

b <- as.matrix(coefs(fit))
b <- setdiff(names(b[b[, 1] != 0, ]), '(Intercept)')
form <- paste('y ~', paste0(b, collapse = ' + '))
trainModel <- glm(form, data = data, family = binomial)

prc_cutoff <- relapseRisk::ROCcurve(trainModel, prc = TRUE)$optimal['cutoff']

trainModel <- coef(trainModel)

usethis::use_data(epsi_idas_questions, trainModel, prc_cutoff, internal = TRUE, overwrite = TRUE)
