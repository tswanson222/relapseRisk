#' Classperf
#'
#' @param pred TBA
#' @param obs TBA
#' @param metric TBA
#' @param p TBA
#'
#' @return TBA
#' @export
#'
#' @examples
#' 1 + 1
classperf <- function(pred, obs, metric = 'sens', p = levels(obs)[2]){
  stopifnot(is.factor(pred) & is.factor(obs))
  stopifnot(length(unique(c(levels(obs), levels(pred)))) == 2)
  n <- setdiff(levels(obs), p)
  metric <- match.arg(tolower(metric), c(
    'sensitivity', 'specificity', 'precision', 'ppv', 'kappa', 'f1',
    'fdr', 'mcc', 'accuracy', 'ba'), several.ok = TRUE)
  if('ppv' %in% metric){metric[metric == 'ppv'] <- 'precision'}
  metric <- capitalize(unique(metric))
  if(any(nchar(metric) < 4)){
    metric[nchar(metric) < 4] <- toupper(metric[nchar(metric) < 4])
  }
  #caret:::requireNamespaceQuietStop("e1071") # NEW
  kappa <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]
  kappa <- kappa['kappa']
  tp <- sum(pred %in% p & obs %in% p)
  tn <- sum(pred %in% n & obs %in% n)
  fp <- sum(pred %in% p & obs %in% n)
  fn <- sum(pred %in% n & obs %in% p)
  tpr <- tp/(tp + fn)
  tnr <- tn/(tn + fp)
  ppv <- tp/(tp + fp)
  fdr <- 1 - ppv
  ba <- (tpr + tnr)/2
  acc <- (tp + tn)/(tp + tn + fp + fn)
  f1 <- 2 * ((ppv * tpr)/(ppv + tpr))
  mcc <- suppressWarnings(((tp * tn) - (fp * fn))/sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
  output <- list(Sensitivity = tpr, Specificity = tnr, Precision = ppv, FDR = fdr,
                 Accuracy = acc, MCC = mcc, BA = ba, Kappa = kappa, F1 = f1)
  return(unlist(output[names(output) %in% metric]))
}
