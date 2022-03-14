# Code copied from "full analysis" script from FWB AMD M-N work
# Path to all of my csv's
upperDir <- "C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\FWsurvey\\invert_meas_CJP\\meas_csv\\" 
# vector of all of the paths to csv's within folder
meas <- list.files(upperDir, full.names=T) 
# just the file name of each csv
filename <- list.files(upperDir) 

#### Loop Reading in csv's ####
dataset <- NULL #Empty dataset

for (i in 1:length(meas)){ #length(meas) = number of total files within folder
  # creating temporary dataset for file
  temp_dataset <- read.csv(meas[i],
                           header=T) 
  temp_dataset <- temp_dataset[1:5] # forcing file to only have 4 columns. When I export my measurements into excel, Adobe sometimes adds an additional column
  temp_dataset$file <- filename[i] # adding new colum "file" to temp data set
  dataset <- rbind(dataset,
                   temp_dataset[,1:6]) #row binding 5 columns from temp dataset to full dataset
  rm(temp_dataset) #clearing out temp data for next file
}


# remove "scale" measurements
dataset <- dataset[dataset$Label != "",]#"blank"

# pull out information from "file" column
# pulls out site code from file column
dataset$site <- substr(dataset$file,
                       start = 1,
                       stop = 6)
# pulls out Surber number
dataset$surber <- substr(dataset$file,
                         start = 8,
                         stop = 9)
# pulls out date
dataset$date <- substr(dataset$file,
                       start = 12,
                       stop = 19)
# remove "-" and "_" from site column
dataset$site <- gsub(pattern = "-",
                     replacement = "",
                     x = dataset$site)
dataset$site <- gsub(pattern = "_",
                     replacement = "",
                     x = dataset$site)
# cleaning up names
dataset$Genus <- gsub(
  pattern = "Adult coleoptera",
  replacement = "Adult Coleoptera",
  x = dataset$Genus)
dataset$site <- as.factor(dataset$site)
