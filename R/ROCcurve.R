#' Receiver operating characteristic and precision-recall curves
#'
#' Construct receiver operating characteristic (ROC) and precision-recall curves
#' (PRC).
#'
#' @param y Outcome variable vector
#' @param X Matrix of predictors
#' @param model Optional. Supply logistic regression model.
#' @param plot Logical. Determines whether to plot the ROC or PRC results.
#' @param optPoint Logical. Determines whether to plot the optimal threshold
#'   value on the plot.
#' @param grid TBD
#' @param grid_lty TBD
#' @param grid_lwd TBD
#' @param grid_col TBD
#' @param midline TBD
#' @param midline_lty TBD
#' @param midline_lwd TBD
#' @param midline_col TBD
#' @param pt_pch TBD
#' @param pt_border TBD
#' @param pt_col TBD
#' @param thresh TBD
#' @param roc_lty TBD
#' @param roc_lwd TBD
#' @param roc_col TBD
#' @param prc Logical. Determines whether to create a PRC (if \code{TRUE}) or an
#'   ROC (if \code{FALSE}).
#' @param plot_na Logical
#'
#' @return TBD
#' @export
#'
#' @examples
#' 1 + 1
ROCcurve <- function(y, X = NULL, model = NULL, plot = FALSE, optPoint = TRUE,
                     grid = FALSE, grid_lty = 3, grid_lwd = 1.5, grid_col = "lightgray",
                     midline = TRUE, midline_lty = 2, midline_lwd = 2, midline_col = "red",
                     pt_pch = 23, pt_border = "black", pt_col = "green", thresh = TRUE,
                     roc_lty = 1, roc_lwd = 2, roc_col = "black", prc = FALSE, plot_na = FALSE){
  stopifnot(all(sapply(list(prc, plot, plot_na), is.logical)))
  if(is(y, 'ROCcurve')){
    sens <- y$results$sens
    spec <- y$results$spec
    ppv <- y$results$ppv
    opt <- ifelse(prc, which.max(ppv + sens), which.max(sens + spec))
  } else {
    if(!missing(y) & is.null(model)){
      if(is(y, 'glm')){
        model <- y
        y <- unname(model$y)
      }
    } else if(!is.null(model) & missing(y)){
      y <- unname(model$y)
    }
    stopifnot(dim(table(y)) == 2)
    if(is.factor(y) | is.character(y)){
      y <- factor(y)
      levels(y) <- 0:1
      y <- as.numeric(as.character(y))
    }
    stopifnot(all(names(table(y)) %in% c("0", "1")))
    if(is.null(model)){
      stopifnot(!is.null(X))
      X <- as.data.frame(X)
      model <- glm(y ~ ., data = X, family = binomial)
      predProbs <- predict(model, type = "response")
    } else if(is(model, 'train')){
      predProbs <- predict(model, X, 'prob')[, 2]
    } else if(is(model, 'glinternet')){
      predProbs <- predict(model, X, type = 'response')[, 2]
    } else if(is(model, 'glm') | is(model, 'glmboost') | is(model, 'gbm')){
      predProbs <- predict(model, type = "response")
    } else if(is(model, 'numeric')){
      stopifnot(length(model) == length(y))
      predProbs <- model
    }
    p <- unname(sort(predProbs))
    if(thresh | prc){
      t1 <- (c(-Inf, p) + c(p, +Inf))/2
      t2 <- (c(-Inf, p)/2 + c(p, +Inf)/2)
      p <- ifelse(abs(t1) > 1e+100, t2, t1)
    }
    preds <- list()
    tp <- tn <- fp <- fn <- c()
    for(i in 1:length(p)){
      preds[[i]] <- ifelse(predProbs > p[i], 1, 0)
      Y <- cbind(y, preds[[i]])
      tn[i] <- sum(Y[Y[, 1] == 0, 1] == Y[Y[, 1] == 0, 2])
      tp[i] <- sum(Y[Y[, 1] == 1, 1] == Y[Y[, 1] == 1, 2])
      fn[i] <- sum(Y[Y[, 2] == 0, 1] != Y[Y[, 2] == 0, 2])
      fp[i] <- sum(Y[Y[, 2] == 1, 1] != Y[Y[, 2] == 1, 2])
    }
    sens <- tp/sum(y)
    spec <- tn/(length(y) - sum(y))
    npv <- tn/(tn + fn)
    ppv <- tp/(tp + fp)
    opt <- ifelse(prc, which.max(ppv + sens), which.max(sens + spec))
    optCut <- p[opt]
    optSens <- sens[opt]
    optSpec <- spec[opt]
    optPPV <- ppv[opt]
    optNPV <- npv[opt]
  }
  if(prc){
    sx <- c(1, sens)
    sy <- c(0, ppv)
    syna <- 0
    if(any(is.na(sy))){
      if(sum(is.na(sy)) == 1){
        sy <- na.omit(sy)
      } else {
        syna <- which(is.na(sy))
        sy[syna] <- sx[syna]
      }
    }
    syextra <- length(sy) < length(sx)
    if(syextra){sy <- c(sy, 1)}
    height <- sy[-1] - sy[-length(sy)]
    width <- (sx[-1] + sx[-length(sx)])/2
  } else {
    sx <- c(0, spec)
    sy <- c(1, sens)
    height <- (sy[-1] + sy[-length(sy)])/2
    width <- sx[-1] - sx[-length(sx)]
  }
  AUC <- sum(height * width)
  if(!plot){
    if(is(y, 'ROCcurve')){
      p <- y$results$cutoff
      npv <- y$results$npv
      optCut <- p[opt]
      optSens <- sens[opt]
      optSpec <- spec[opt]
      optPPV <- ppv[opt]
      optNPV <- npv[opt]
    }
    out <- list(results = data.frame(cutoff = p, sens = sens, spec = spec, ppv = ppv, npv = npv),
                optimal = unlist(list(cutoff = optCut, sensitivity = optSens,
                                      specificity = optSpec, PPV = optPPV,
                                      NPV = optNPV)), AUC = AUC)
    class(out) <- c('ROCcurve', 'list')
    attr(out, 'type') <- ifelse(prc, 'PRC', 'ROC')
    return(out)
  } else {
    xlabel <- ifelse(prc, 'Recall', '1 - Specificity')
    ylabel <- ifelse(prc, 'Precision', 'Sensitivity')
    main <- ifelse(prc, paste0('PR Curve\nAUC = ', round(AUC, 3)),
                   paste0('ROC Curve\nAUC = ', round(AUC, 3)))
    plot(0, 0, type = "n", ylim = c(0, 1), xlim = c(0, 1), axes = FALSE,
         xlab = xlabel, ylab = ylabel, main = main)
    if(grid != FALSE){
      if(grid == TRUE){grid <- 1}
      if(grid == 1){
        grid(NA, 5, lty = grid_lty, lwd = grid_lwd, col = grid_col)
      } else if(grid == 2){
        grid(lty = grid_lty, lwd = grid_lwd, col = grid_col)
      }
    }
    axis(1); axis(2)
    if(!plot_na & prc){
      if(syextra){
        sy <- sy[-length(sy)]
        sx <- sx[-length(sx)]
      }
      if(!identical(syna, 0)){
        sx[syna] <- NA
        sy[syna] <- NA
      }
    }
    if(prc){sx <- 1 - sx}
    lines(1 - sx, sy, lty = roc_lty, lwd = roc_lwd, col = roc_col)
    if(midline != FALSE){
      if(midline == TRUE){midline <- 1}
      if(midline == "grey" | midline == "gray"){midline_lty <- 1; midline_col = "grey"}
      if(midline == 2){midline_lty <- 1; midline_col = "grey"}
      prctrue <- as.numeric(prc)
      abline(a = prctrue, b = 1 - (2 * prctrue), lty = midline_lty,
             lwd = midline_lwd, col = midline_col)
    }
    if(optPoint != FALSE){
      if(optPoint == "black" | optPoint == 2){pt_pch <- 8; pt_col <- "black"}
      if(!pt_pch %in% c(21:25)){pt_border <- pt_col}
      points((1 - sx)[opt + 1], sy[opt + 1], pch = pt_pch, bg = pt_col, col = pt_border)
    }
  }
}
