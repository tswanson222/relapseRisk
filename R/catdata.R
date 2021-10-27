catdata <- function(file = NULL, day = NULL, type = 'cat', process = TRUE,
                    questions = 'last'){
  if(is.null(file)){file <- 'forbushtest1_7_31_20.csv'}
  questions_epsi <- trimws(as.character(readRDS('questions_epsi.RDS')))
  questions_epsi[37] <- gsub('<93>unhealthy<94>', 'unhealthy', questions_epsi[37])
  type <- match.arg(tolower(type), c('cat', 'demographics', 'weekly'))
  x <- read.csv(file, stringsAsFactors = FALSE)
  k <- unique(x$survey_code)
  k <- setdiff(k, c('cat_demographic_survey', 'weekly_behaviors'))
  if(type != 'cat'){
    k <- switch(type, demographics = 'cat_demographic_survey',
                weekly = 'weekly_behaviors')
  }
  x2 <- subset(x, survey_code == k & event_type == 'response' & question_type != 'instruction')
  x2$date <- as.Date(x2$metadata..timestamp)
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
    resopts <- c('Never', 'Rarely', 'Sometimes', 'Often', 'Very Often')
    if(questions == 'last'){
      values <- values[[length(values)]]
      for(i in 0:4){values[values == i] <- resopts[i + 1]}
      items <- items[[length(items)]]
      responses <- cbind.data.frame(Items = questions_epsi[items],
                                    Responses = factor(values, levels = resopts[which(resopts %in% values)]))
    } else if(is.numeric(questions)){
      values <- values[[questions]]
      for(i in 0:4){values[values == i] <- resopts[i + 1]}
      items <- values[[questions]]
      responses <- cbind.data.frame(Items = questions_epsi[items],
                                    Responses = factor(values, levels = resopts[which(resopts %in% values)]))
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
    levels(output$Thetas) <- c('Bulimia', 'Excessive Exercise', 'Restriction')
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
