# data prep

upperDir <- "data/raw_csv" 
# vector of all of the paths to csv's within folder
meas <- list.files(upperDir, full.names=T) 
# just the file name of each csv
filename <- list.files(upperDir) 



#### Loop Reading in csv's ####
dataset <- NULL #Empty dataset

for (i in 1:length(meas)){ #length(meas) = number of total files within folder
  # adding empty line to csv files
  cat("\n", file = meas[i], append = TRUE)
  # creating temporary dataset for file
  temp_dataset <- read.csv(meas[i],
                           header=T) 
  temp_dataset <- temp_dataset[1:4] # forcing file to only have 4 columns. When I export my measurements into excel, Adobe sometimes adds an additional column
  temp_dataset$file <- filename[i] # adding new colum "file" to temp data set
  dataset <- rbind(dataset,
                   temp_dataset) #row binding from temp dataset to full dataset
  rm(temp_dataset) #clearing out temp data for next file
}

names(dataset)
complete.cases(dataset)

dataset[,1:4]
unique(dataset$Unit)
sort(unique(dataset$Label))
# Cleaning up some names

dataset$Label <- gsub(
  pattern = "Brachycentriidae",
  replacement = "Bracycentridae",
  x = dataset$Label)

# start here
dataset$Label <- gsub(
  pattern = "Bracycentridae",
  replacement = "Brachycentridae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Hrydroptilidae",
  replacement = "Hydroptilidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Pteronidae",
  replacement = "Pteronarcyidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Rhycacophilidae",
  replacement = "Rhyacophilidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Flat_worm",
  replacement = "Planariidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Flatworm",
  replacement = "Planariidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Hydrachnidae",
  replacement = "Acari",
  x = dataset$Label)


sort(unique(dataset$Label))

saveRDS(dataset, "data/ark_csv_stitched.R")
