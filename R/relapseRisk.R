#' relapseRisk: Predict risk of relapse for anorexia nervosa recovery
#'
#' Implements a machine learning algorithm for predicting recovery outcomes
#' among patients diagnosed with anorexia nervosa. Designed for use with data
#' from the STAR app in PiLR for a research project funded by the NIH.
#'
#' @docType package
#' @name relapseRisk
#'
#' @import stats
#' @import ggplot2
#' @importFrom methods formalArgs is
#' @importFrom graphics abline axis lines points
#' @importFrom utils combn read.csv

utils::globalVariables(
  c('question_type', 'session', 'survey_code', 'event_type',
    'Estimate', 'SE', 'Thetas', 'Time')
)
