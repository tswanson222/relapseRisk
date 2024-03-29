% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROCcurve.R
\name{ROCcurve}
\alias{ROCcurve}
\title{Receiver operating characteristic and precision-recall curves}
\usage{
ROCcurve(
  y,
  X = NULL,
  model = NULL,
  plot = FALSE,
  optPoint = TRUE,
  grid = FALSE,
  grid_lty = 3,
  grid_lwd = 1.5,
  grid_col = "lightgray",
  midline = TRUE,
  midline_lty = 2,
  midline_lwd = 2,
  midline_col = "red",
  pt_pch = 23,
  pt_border = "black",
  pt_col = "green",
  thresh = TRUE,
  roc_lty = 1,
  roc_lwd = 2,
  roc_col = "black",
  prc = FALSE,
  plot_na = FALSE,
  cutoff = ">"
)
}
\arguments{
\item{y}{Outcome variable vector}

\item{X}{Matrix of predictors}

\item{model}{Optional. Supply logistic regression model.}

\item{plot}{Logical. Determines whether to plot the ROC or PRC results.}

\item{optPoint}{Logical. Determines whether to plot the optimal threshold
value on the plot.}

\item{grid}{TBD}

\item{grid_lty}{TBD}

\item{grid_lwd}{TBD}

\item{grid_col}{TBD}

\item{midline}{TBD}

\item{midline_lty}{TBD}

\item{midline_lwd}{TBD}

\item{midline_col}{TBD}

\item{pt_pch}{TBD}

\item{pt_border}{TBD}

\item{pt_col}{TBD}

\item{thresh}{TBD}

\item{roc_lty}{TBD}

\item{roc_lwd}{TBD}

\item{roc_col}{TBD}

\item{prc}{Logical. Determines whether to create a PRC (if \code{TRUE}) or an
ROC (if \code{FALSE}).}

\item{plot_na}{Logical}

\item{cutoff}{Character string. Define how to classify cases relative to the
predicted probability cutoff}
}
\value{
TBD
}
\description{
Construct receiver operating characteristic (ROC) and precision-recall curves
(PRC).
}
\examples{
1 + 1
}
