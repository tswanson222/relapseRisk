devtools::document()

# 1) Load DCEB data for high-risk area computation

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

# 4) save as internal data

usethis::use_data(quantile90,quantile95,zscore1,zscore1.5,zscore2,
                  internal = TRUE, overwrite = FALSE)
