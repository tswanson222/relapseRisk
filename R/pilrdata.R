#' Pull survey responses and theta values from PiLR data
#'
#' Used for pulling patient data and organizing for therapist output.
#'
#' @param file PiLR dataframe, or character string indicating the path of where
#'   to find the relevant .csv file.
#' @param survey Specify either \code{"epsi"} or \code{"idas"}
#' @param day Indicate a particular day to draw data from, using the format
#'   \code{"YYYY-MM-DD"}
#' @param type Select between \code{"cat", "demographics", "weekly"}. Currently
#'   only works with \code{"cat"}
#' @param process Logical. If \code{TRUE}, re-formats PiLR data and returns
#'   theta values (with standard errors) along with the relevant items and their
#'   associated responses. Setting to \code{FALSE} returns the raw PiLR data for
#'   the CAT surveys.
#' @param questions Determines which set of responses to report. Can be a
#'   numerical value or vector to specify which time points to return items and
#'   responses for. Setting \code{questions = "last"} will only return the set
#'   of responses from the most recent time point. Setting \code{questions =
#'   "all"} will return the set of responses from each time point as a list.
#'
#' @return A list containing theta values along with standard errors, as well as
#'   table of participant items and responses.
#' @export
#'
#' @examples
#' 1 + 1
pilrdata <- function(file = NULL, survey = c('epsi', 'idas'), day = NULL,
                     type = 'cat', process = TRUE, questions = 'last'){
  stopifnot(!is.null(file))
  survey <- match.arg(survey)
  qoptions <- epsi_idas_questions[[which(endsWith(names(epsi_idas_questions), survey))]]
  x <- switch(2 - is.character(file), read.csv(file, stringsAsFactors = FALSE), file)
  type <- match.arg(tolower(type), c('cat', 'demographics', 'weekly'))
  if(type != 'cat'){
    k <- switch(type, demographics = 'cat_demographic_survey',
                weekly = 'weekly_behaviors')
  } else {
    k <- setdiff(unique(x$survey_code), c('cat_demographic_survey', 'weekly_behaviors'))
    if(identical(k, '68158')){ # This is only to make it compatible with "KU CAT"
      survey <- 'epsi'
    } else {
      k <- k[startsWith(tolower(k), survey)]
    }
  }
  x2 <- subset(x, survey_code == k & event_type == 'response' & question_type != 'instruction')
  x2$date <- as.Date(x2[, grep('^metadata.*.timestamp$', colnames(x))])
  if(!is.null(day)){x2 <- subset(x2, date == day)}
  if(isTRUE(process) & type == 'cat'){
    ### Get Thetas
    thetas <- x2[x2$question_code == 'thetas', 'response_values']
    temp <- vector('list', length(thetas))
    for(i in 1:length(temp)){
      temp[[i]] <- trimws(strsplit(thetas[i], ',')[[1]])
      temp[[i]] <- as.numeric(trimws(gsub('[[]|[]]', '', temp[[i]])))
      if(i == length(temp)){
        if(i == 1){
          thetas <- temp[[1]]
        } else {
          thetas <- do.call(rbind, temp)
        }
      }
    }
    ### Get SEs
    ses <- x2[x2$question_code == 'SE_thetas', 'response_values']
    temp <- vector('list', length(ses))
    for(i in 1:length(temp)){
      temp[[i]] <- trimws(strsplit(ses[i], ',')[[1]])
      temp[[i]] <- as.numeric(trimws(gsub('[[]|[]]', '', temp[[i]])))
      if(i == length(temp)){
        if(i == 1){
          ses <- temp[[1]]
        } else {
          ses <- do.call(rbind, temp)
        }
      }
    }
    ### Get responses/items
    sessions <- unique(x2$session)
    values <- items <- vector('list', length(sessions))
    for(i in 1:length(sessions)){
      x0 <- subset(x2, session == sessions[i])
      x0 <- x0[grep('^mc', x0$question_code), ]
      values[[i]] <- as.numeric(x0$response_value)
      items[[i]] <- as.numeric(gsub('mc:', '', x0$question_code))
    }
    responses <- list(values = values, items = items)
    if(survey == 'epsi'){
      resopts <- c('Never', 'Rarely', 'Sometimes', 'Often', 'Very Often')
    } else {
      resopts <- c('Not at all', 'A little bit', 'Moderately', 'Quite a bit', 'Extremely')
    }
    if(identical(questions, 'all')){questions <- 1:length(values)}
    if(identical(questions, 'last')){
      values <- values[[length(values)]]
      for(i in 0:4){values[values == i] <- resopts[i + 1]}
      items <- items[[length(items)]]
      responses <- cbind.data.frame(
        Items = qoptions[items],
        Responses = factor(values, levels = resopts[which(resopts %in% values)])
      )
    } else if(is.numeric(questions)){
      responses <- setNames(vector('list', length(questions)), paste0('time', questions))
      for(i in seq_along(questions)){
        values0 <- values[[questions[i]]]
        for(j in 0:4){values0[values0 == j] <- resopts[j + 1]}
        items0 <- items[[questions[i]]]
        responses[[i]] <- cbind.data.frame(
          Items = qoptions[items0],
          Responses = factor(values0, levels = resopts[which(resopts %in% values0)])
        )
      }
      if(length(responses) == 1){responses <- responses[[1]]}
    }
    if(length(sessions) == 1){
      output <- data.frame(Estimate = thetas, SE = ses,
                           Thetas = factor(paste0('theta', 1:3)),
                           Time = rep(1, 3))
    } else {
      time <- rep(1:length(sessions), 3)
      labs <- rep(paste0('theta', 1:3), each = length(sessions))
      thetas <- c(thetas[, 1], thetas[, 2], thetas[, 3])
      ses <- c(ses[, 1], ses[, 2], ses[, 3])
      output <- data.frame(Estimate = thetas, SE = ses, Thetas = factor(labs), Time = time)
    }
    if(survey == 'epsi'){
      #levels(output$Thetas) <- c('Bulimia', 'Excessive Exercise', 'Restriction')
      levels(output$Thetas) <- c('Bulimia', 'Exercise-Focused Behaviors', 'Restrictive Eating')
    } else {
      levels(output$Thetas) <- c('Fear', 'Distress', 'Positive Affect')
    }
    x2 <- list(output = output, responses = responses)
  }
  x2
}


#' Plot thetas from PiLR app data
#'
#' Description
#'
#' @param data Output from \code{\link{pilrdata}}, or raw PiLR data, or path to
#'   PiLR data .csv file.
#' @param colors Logical. Determines whether to use different colors for each
#'   higher-order dimension.
#' @param zeroLine Logical. Determines whether to draw a dashed line
#' @param pointSize Numerical value for size of points on plot.
#' @param ... Additional arguments for \code{\link{pilrdata}}
#'
#' @return A plot of theta values and their standard errors for each
#'   higher-order dimension.
#' @export
#'
#' @examples
#' 1 + 1
plotThetas <- function(data, colors = FALSE, zeroLine = TRUE, pointSize = 3, ...){
  if(is(data, 'list')){data <- data$output}
  if(ifelse(is(data, 'data.frame'), 'survey_code' %in% colnames(data), is.character(data))){
    data <- pilrdata(data, ...)$output
  }
  thetas <- data
  rects <- 1
  whiskers <- ifelse(max(thetas$Time) == 1, 0, .035)
  background <- data.frame(xlow = c(0, 0), xhigh = rep(max(thetas$Time) + 1, 2),
                           ylow = c(1.5, -2))
  if(FALSE){ # max(thetas$Time) == 1
    if(colors){
      g1 <- ggplot(thetas, aes(x = Thetas, y = Estimate, colour = Thetas)) +
        geom_point(size = pointSize, show.legend = c(colour = FALSE))
    } else {
      g1 <- ggplot(thetas, aes(x = Thetas, y = Estimate)) +
        geom_point(size = pointSize)
    }
    g1 <- g1 + geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = .1, lwd = 1) +
      theme_bw() + scale_x_discrete(labels = c('theta1' = expression(theta[1]), 'theta2' = expression(theta[2]), 'theta3' = expression(theta[3]))) +
      theme(axis.text.x = element_text(size = 18), axis.title.x = element_blank(), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12)) +
      scale_color_discrete(labels = c('theta1' = expression(theta[1]), 'theta2' = expression(theta[2]), 'theta3' = expression(theta[3])))
    if(colors){g1 <- g1 + theme(legend.text = element_text(size = 12))}
    if(zeroLine){g1 <- g1 + geom_hline(yintercept = 0, linetype = 'dashed', alpha = .5)}
  } else {
    if(colors){
      g1 <- ggplot(thetas, aes(x = Time, y = Estimate, colour = Thetas)) +
        geom_point(size = 1.5, show.legend = c(colour = FALSE)) + geom_line(alpha = .4) +
        geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = .035, lwd = .5) + theme_bw() + facet_grid(Thetas ~ .) +
        scale_color_discrete(labels = c('theta1' = expression(theta[1]), 'theta2' = expression(theta[2]), 'theta3' = expression(theta[3])))
    } else {
      g1 <- ggplot(thetas, aes(x = Time, y = Estimate)) + geom_point(size = 1.5) + geom_line(alpha = .4) +
        geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = whiskers, lwd = .5) +
        theme_bw() + facet_grid(Thetas ~ .)
    }
    if(rects > 0){
      g1 <- g1 + geom_rect(mapping = aes(xmin = 1, xmax = max(thetas[, 'Time']), ymin = 1.5, ymax = Inf), alpha = .025, fill = 'red')
      if(rects > 1){
        g1 <- g1 + geom_rect(mapping = aes(xmin = 1, xmax = max(thetas[, 'Time']), ymin = -2, ymax = 1.5), alpha = .025, fill = 'green')
      }
    }
    g1 <- g1 + ylab('Estimate (95% CI)') + scale_x_continuous(breaks = 1:max(thetas$Time), labels = 1:max(thetas$Time))
  }
  return(g1)
}


#' Create table of items and responses from PiLR app
#'
#' Description
#'
#' @param data Output from \code{\link{pilrdata}}, or raw PiLR data, or path to
#'   PiLR data .csv file.
#' @param time Numeric value to indicate which set of responses to use, or set
#'   to \code{"last"} to choose the most recent set of responses.
#' @param ... Additional arguments for \code{\link{pilrdata}}
#'
#' @return A kable table
#' @export
#'
#' @examples
#' 1 + 1
itemTable <- function(data, time = 'last', ...){
  if(ifelse(is(data, 'data.frame'), 'survey_code' %in% colnames(data), is.character(data))){
    data <- pilrdata(data, ...)$responses
  }
  if(ifelse(is(data, 'list'), all(c('output', 'responses') %in% names(data)), FALSE)){
    data <- data$responses
  }
  if(is(data, 'list')){
    if(time == 'last' | !is(time, 'numeric')){time <- length(data)}
    data <- data[[time]]
  }
  resopts_epsi <- c('Never', 'Rarely', 'Sometimes', 'Often', 'Very Often')
  resopts_idas <- c('Not at all', 'A little bit', 'Moderately', 'Quite a bit', 'Extremely')
  resopts <- switch(2 - any(resopts_epsi %in% data$Responses), resopts_epsi, resopts_idas)
  colors <- c('#006D2C', '#41AB5D', '#FDAE61', '#E31A1C', '#800026')
  data <- data[order(data$Responses), ]
  rownames(data) <- 1:nrow(data)
  data <- within(data, {
    Responses <- kableExtra::cell_spec(Responses, color = 'white', bold = TRUE,
                                       background = factor(Responses, resopts, colors))
  })
  out <- kableExtra::kable_styling(kableExtra::kable(data, escape = FALSE),
                                   bootstrap_options = c('striped', 'hover', 'responsive'))
  return(out)
}


#' Create patient report
#'
#' Description
#'
#' @param participant Raw PiLR data or path to .csv file. If path to file is
#'   provided, this must be an absolute path not a relative path.
#' @param id ID number of participant
#' @param output Choose the directory of where to output the report. Defaults to
#'   the current working directory.
#'
#' @return Patient report
#' @export
#'
#' @examples
#' 1 + 1
createReport <- function(participant, id = 0, output = './'){
  date <- Sys.Date()
  input <- system.file('rmd', 'patient_report.Rmd', package = 'relapseRisk')
  colors <- c('green', 'red', 'yellow')
  stoplights <- setNames(sapply(paste0(colors, '_stoplight.png'), function(z){
    system.file('stoplights', z, package = 'relapseRisk')
  }), colors)
  rmarkdown::render(input = input,
                    output_file = paste0('patient', id, '_', date, '.html'),
                    output_dir = output,
                    params = list(participant = participant,
                                  set_title = paste0('Patient ', id),
                                  set_date = date,
                                  stoplights = stoplights))
}


# Predict risk of relapse based on ML model
predictRisk <- function(data, ...){
  data <- lapply(c('epsi', 'idas'), function(z) pilrdata(data, z))
  data <- do.call(rbind, lapply(data, '[[', 'output'))
  epsi <- c('Bulimia', 'Exercise-Focused Behaviors', 'Restrictive Eating')
  idas <- c('Fear', 'Distress', 'Positive Affect')
  stopifnot(all(c(epsi, idas) %in% data$Thetas))
  data <- split(data, data$Time)
  form <- as.formula(paste0('~', paste0(names(trainModel)[-1], collapse = ' + ')))
  for(i in 1:length(data)){
    data[[i]]$Thetas <- factor(data[[i]]$Thetas, levels = c(epsi, idas),
                               labels = c(paste0('EP', 1:3), paste0('I', 1:3)))
    temp <- structure(data.frame(t(data.frame(data[[i]]$Estimate))), row.names = '1')
    colnames(temp) <- data[[i]]$Thetas
    data[[i]] <- model.matrix(form, temp)
    b <- (data[[i]] %*% trainModel)[, 1]
    data[[i]] <- ifelse(exp(b)/(1 + exp(b)) > prc_cutoff, 'Poor', 'Marginal')
  }
  out <- do.call(rbind, data)
  out <- data.frame(outcome = out[, 1], Time = 1:nrow(out))
  return(out)
}

