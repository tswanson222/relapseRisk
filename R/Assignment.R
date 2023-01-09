# Read from STAR_therapist_startdate.csv file to get the therapist and startdate of each participant
# Create separate files for each individual with the name: [therapist code]_[epoch code]_[participant id].csv
# These csv files are then processed by functions in pilrdata.R to obtain ED factors and generate reports
# Reports are generated under the name of [therapist code]_[epoch code]_[participant id].html
# Convert all reports into the [therapist code]_[epoch code]_[participant id].pdf format to upload
# Assigns new configurations for each participant
# Create Bulk Config Upload file (.csv)
# column headers: participant_id, epoch_code, config_code


# file should be the full data file from PiLR

BESTU_assignment <- function(file = NULL){
  
  data <- read.csv(file)
    
  # reference csv
  refer <- read.csv("STAR_therapist_startdate.csv")
  
  # check if all participants are in the reference csv
  participants <- unique(data$metadata..pt)
  
  # three types: participants in both; participants in file only; participants in log only
  
  p_both <- intersect(participants, refer$ID)
  p_file_only <- participants[!(participants %in% refer$ID)]
  p_ref_only <- refer$ID[!(refer$ID %in% participants)]
  
  # Warning messages
  
  if(length(p_file_only)>0){
    warning("Update the log! Some participants are not logged.")
  }
  
  dir.create(paste0(getwd(),"/",Sys.Date()))
 
  #####################################################
  # create files for p_both
  
  # remove all test subjects
  
  ref <- subset(refer,(ID %in% p_both) & therapist != "TEST")
  p_both <- ref$ID
  
  for(i in 1:length(p_both)){
  }
  
}

