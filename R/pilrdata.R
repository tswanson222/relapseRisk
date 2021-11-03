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

#x <- catdata('trevor_test2.csv', day = '2020-07-31')




################################################################################

#thetas <- trimws(strsplit(thetas, ',')[[1]])
#thetas <- as.numeric(trimws(gsub('[[]|[]]', '', thetas)))

#ses <- trimws(strsplit(ses, ',')[[1]])
#ses <- as.numeric(trimws(gsub('[[]|[]]', '', ses)))


#values0 <- values[[length(values)]]
#items0 <- items[[length(items)]]

#x <- x[grep('^mc', x$question_code), ]
#values <- as.numeric(x$response_value)
#items <- as.numeric(gsub('mc:', '', x$question_code))
