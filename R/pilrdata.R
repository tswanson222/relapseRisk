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
                     type = 'cat', process = TRUE, questions = 'last',
                     time_start = "2023-01-01"){
  time_now <- lubridate::today()
  Weeks <- floor(difftime(time_now,time_start,units="weeks"))
  if(is.null(file)){
    x2 <- list(non_response =TRUE,
               tech_inter=FALSE)
    return(x2)
  }
  
 else{
  survey <- match.arg(survey)
  qoptions <- epsi_idas_questions[[which(endsWith(names(epsi_idas_questions), survey))]]
  x <- switch(2 - is.character(file), read.csv(file, stringsAsFactors = FALSE), file)
##### extract max time points
  Weeks_extraction <- unique(x$survey_code)[stringr::str_detect(unique(x$survey_code),'EPSI|IDAS')]
  Weeks <- max(as.integer(Weeks),max(as.integer(stringr::str_extract(Weeks_extraction, "\\d+"))))
  
  type <- match.arg(tolower(type), c('cat', 'demographics', 'weekly'))
  if(type == 'cat'){
    if('68158' %in% x$survey_code){
      # This is only to make it compatible with the defunct "KU CAT"
      x$survey_code <- gsub('68158', 'EPSI_survey', x$survey_code)
    }
    k <- setdiff(unique(x$survey_code), c('cat_demographic_survey', 'weekly_behaviors'))
    k <- k[startsWith(tolower(k), survey)]
  } else {
    k <- switch(type, demographics = 'cat_demographic_survey',
                weekly = 'weekly_behaviors')
  }
  k_all <- k
  x2 <- subset(x, survey_code %in% k & event_type == 'response' & question_type != 'instruction')
  k <- x2[x2$question_code == 'thetas', "survey_code"]
  x2 <- subset(x, survey_code %in% k & event_type == 'response' & question_type != 'instruction')
  x2$date <- as.Date(x2[, grep('^metadata.*.timestamp$', colnames(x))])
  if(!is.null(day)){x2 <- subset(x2, date == day)}
  if(isTRUE(process) & type == 'cat'){

    ### Check if any thetas exist
    if(any(x2$question_code == 'thetas')){

      ### Get Thetas
      thetas <- x2[x2$question_code == 'thetas', 'response_values']
      theta_dates <- x2[x2$question_code == 'thetas', 'date']
      theta_surveys <- x2[x2$question_code == 'thetas', 'survey_code']
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
      
      ### Get time points
      survey_weeks <- as.integer(stringr::str_extract(k, "\\d+"))

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
      }
      else if(identical(questions,"current")){
        values <- values[[Weeks]]
        for(i in 0:4){values[values == i] <- resopts[i + 1]}
        items <- items[[length(items)]]
        responses <- cbind.data.frame(
          Items = qoptions[items],
          Responses = factor(values, levels = resopts[which(resopts %in% values)])
        )       
      }
      else if(is.numeric(questions)){
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

        time <- rep(1:Weeks, 3)
        labs <- rep(paste0('theta', 1:3), each = Weeks)
        thetas <- matrix(thetas,ncol=3)
        ses <- matrix(ses,ncol=3)
        thetas_new <- array(NA,dim=c(Weeks,3))
        ses_new <- array(NA,dim=c(Weeks,3))
        theta_dates_new <- as.Date(rep(NA,Weeks))
        
        for(i in 1:length(survey_weeks)){
          thetas_new[survey_weeks[i],] <- thetas[i,]
          ses_new[survey_weeks[i],] <- ses[i,]
          theta_dates_new[survey_weeks[i]] <- as.Date(theta_dates[i])
        }
        
        thetas <- c(thetas_new[, 1], thetas_new[, 2], thetas_new[, 3])
        ses <- c(ses_new[, 1], ses_new[, 2], ses_new[, 3])
        theta_dates <- theta_dates_new
        output <- data.frame(Estimate = thetas, SE = ses, Thetas = factor(labs),
                             Time = time, Date = rep(theta_dates, each = 3))
      
      if(survey == 'epsi'){
        levels(output$Thetas) <- c('Bulimia', 'Exercise-Focused Behaviors', 'Restrictive Eating')
      } else {
        levels(output$Thetas) <- c('Fear', 'Distress', 'Positive Affect')
      }
      tech_inter <- ifelse(max(survey_weeks)<Weeks & 
                             Weeks %in% as.integer(stringr::str_extract(k_all, "\\d+")),TRUE,FALSE)  
      x2 <- list(output = output, responses = responses,
                 non_response =ifelse(max(survey_weeks)==Weeks,FALSE,TRUE),
                 tech_inter=tech_inter)
    }   else {
      x2 <- list(non_response =TRUE,
               tech_inter=FALSE)
    }
  }
  return(x2)
 }
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
  which_survey <- ifelse('Exercise-Focused Behaviors' %in% data$Thetas,1,2)
  if('Exercise-Focused Behaviors' %in% data$Thetas){
    levels(data$Thetas) <- gsub('Exercise-Focused Behaviors',
                                'Exercise-Focused\nBehaviors',
                                levels(data$Thetas))
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
      ##### this is the line for high-rish areas
      g1 <- g1 + geom_rect(mapping = aes(xmin = 0.5, xmax = max(thetas[, 'Time'])+0.5, 
                                         ymin = rep(zscore1.5[(3*which_survey-2):(3*which_survey)],each=max(thetas$Time)), 
                                         ymax = Inf), alpha = 0.2/max(thetas$Time), fill = 'red')
      ##### 
      if(rects > 1){
        g1 <- g1 + geom_rect(mapping = aes(xmin = 0.5, xmax = max(thetas[, 'Time'])+0.5, ymin = -2, ymax = 1.5), alpha = .05, fill = 'green')
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
itemTable <- function(data, ...){
  if(ifelse(is(data, 'data.frame'), 'survey_code' %in% colnames(data), is.character(data))){
    data <- pilrdata(data, ...)$responses
  }
  if(ifelse(is(data, 'list'), all(c('output', 'responses') %in% names(data)), FALSE)){
    data <- data$responses
  }
  if(is(data, 'list')){
    data <- data$responses
  }
  if(!is.null(data)){
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
createReport <- function(participant, therapist=0, epoch=0, id=0, output = './',time_start = "2023-01-01",questions=questions,
                participant_assignment="treat"){
  date <- Sys.Date()
  input <- system.file('rmd', 'patient_report.Rmd', package = 'relapseRisk')
  colors <- c('green', 'red', 'yellow')
  stoplights <- setNames(sapply(paste0(colors, '_stoplight.png'), function(z){
    system.file('stoplights', z, package = 'relapseRisk')
  }), colors)
  rmarkdown::render(input = input,
                    output_file = paste0(therapist, '_', epoch, '_', id, '.html'),
                    output_dir = output,
                    params = list(participant = participant,
                                  set_title = paste0('Patient ', id),
                                  set_date = date,
                                  stoplights = stoplights,
                                  time_start=as.Date(time_start),
                                  epoch=epoch,
                                  questions=questions,
                                  participant_assignment=participant_assignment))
}




questionTable <- function(data, week, questions, participant_assignment){
  if(participant_assignment=="treat"){
    corres <- list()
  for(i in 1:3){
    corres[[i]] <- c(2*i-1,2*i)
  }
  corres[[4]] <- c(7,8,9)
  for(i in 5:6){
    corres[[i]] <- c(2*i,2*i+1)
  }
  corres[[7]] <- c(18,24)
  corres[[8]] <- c(19,20,21)
  corres[[9]] <- c(25,26)
  corres[[10]] <- c(28,29)
  corres[[11]] <- 31
  corres[[12]] <- c(32,33)
  
  sub <- subset(data,substr(question_code,1,2) %in% paste0("m",corres[[week]]))
  
  if(!is.null(sub) & dim(sub)[1]>0){
    tab <- data.frame(array(dim=c(100,1)))
    cols <- rep(NA,100)
    k <- 1
    names(tab) <- "Answers to questions"
    for(i in 1:dim(sub)[1]){
      question_id <- sub$question_code[i]
      subsub <- subset(questions,code==question_id)
      tab[k,1] <- paste0(" Q: ",subsub$text[1],"    ")
      if(dim(subsub)[1]>1){
        tab[(k+1):(k+dim(subsub)[1]-1),1] <- paste0("    ",subsub$text[2:dim(subsub)[1]],"    ")
      }
      cols[k:(k+dim(subsub)[1]-1)] <- "dark blue"
      k <- k+dim(subsub)[1]
      subsub <- subset(sub,question_code==question_id)
      a <- subsub$response[1]
      as <- list()
      s <- 1
      while(nchar(a)>110){
         m <- 95
         while(substr(a,m,m)!=" "){
           m <- m+1
         }
         as[[s]] <- substr(a,1,m)
         a <- substr(a,m+1,nchar(a))
         s <- s+1
      }
      as[[s]] <- a
      for(m in 1:s){
        tab[k,1] <- ifelse(m==1,paste0("A: ",as[[m]]),as[[m]])
        cols[k] <- "black"
        k <- k+1
       }
    }
    tab <- tab[1:(k-1),]
    cols <- matrix(cols[1:(k-1)], ncol=1)
    mytheme <- gridExtra::ttheme_default(base_size=20,padding = grid::unit(c(4, 4), "mm"),
                                         core = list(fg_params = list(hjust=0, x=0.01,
                                                                      fontsize=10,col=cols),
                                                     bg_params = list(fill=rep("grey95",
                                                                               length.out=k-1))),
    colhead = list(fg_params = list(fontsize=10, 
                                    fontface="bold"))
    )
    table_length <- dim(tab)[1]
    L <- floor(table_length/18) + 1
    g1 <- gridExtra::tableGrob(tab, theme = mytheme, rows=NULL)
        g1$widths <- grid::unit(rep(1/ncol(g1), ncol(g1)), "npc")
        out <- grid::grid.draw(g1)
  }
  }
  
  else {
      sub <- subset(data,substr(question_code,1,5)=="diary" & survey_code==paste0("diary_w",week))
      if(!is.null(sub) & dim(sub)[1]>0){
          tab <- data.frame(array(dim=c(100,1)))
          cols <- rep(NA,100)
          k <- 1
          names(tab) <- "Answers to questions"
          for(i in 1:dim(sub)[1]){
               question_id <- sub$question_code[i]
               subsub <- subset(questions,code==question_id)
               tab[k,1] <- paste0(" Q: ",subsub$text[1],"    ")
               if(dim(subsub)[1]>1){
                    tab[(k+1):(k+dim(subsub)[1]-1),1] <- paste0("    ",subsub$text[2:dim(subsub)[1]],"    ")
               }
               cols[k:(k+dim(subsub)[1]-1)] <- "dark blue"
               k <- k+dim(subsub)[1]
               subsub <- subset(sub,question_code==question_id)
               a <- subsub$response[1]
               as <- list()
               s <- 1
               while(nchar(a)>110){
                   m <- 95
                   while(substr(a,m,m)!=" "){
                   m <- m+1
                   }
                   as[[s]] <- substr(a,1,m)
                   a <- substr(a,m+1,nchar(a))
                   s <- s+1
               }
               as[[s]] <- a
               for(m in 1:s){
                   tab[k,1] <- ifelse(m==1,paste0("A: ",as[[m]]),as[[m]])
                   cols[k] <- "black"
                   k <- k+1
               }
          }
          tab <- tab[1:(k-1),]
          cols <- matrix(cols[1:(k-1)], ncol=1)
          mytheme <- gridExtra::ttheme_default(base_size=20,padding = grid::unit(c(4, 4), "mm"),
                                         core = list(fg_params = list(hjust=0, x=0.01,
                                                                      fontsize=10,col=cols),
                                                     bg_params = list(fill=rep("grey95",
                                                                               length.out=k-1))),
                                         colhead = list(fg_params = list(fontsize=10, 
                                                    fontface="bold"))
                      )
          table_length <- dim(tab)[1]
          L <- floor(table_length/18) + 1
          g1 <- gridExtra::tableGrob(tab, theme = mytheme, rows=NULL)
          g1$widths <- grid::unit(rep(1/ncol(g1), ncol(g1)), "npc")
          out <- grid::grid.draw(g1)
      }
   }
}





