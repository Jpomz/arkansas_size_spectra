library(tidyverse)

# Download all files from one drive intoa local folder
# 
upperDir <- "cmu/cmu_data" 

# vector of all of the paths to csv's within folder
meas <- list.files(upperDir, full.names=TRUE) 
# just the file name of each csv
filename <- list.files(upperDir) 

# number of files?
length(filename)

#### Loop Reading in csv's ####
dataset <- NULL #Empty dataset


start_time <- Sys.time()

for (i in 1:length(meas)){ 
  #length(meas) = number of total files within folder
  
  # if errors occur when reading in files, you may need to add an empty line to the csv files. 
  # Uncomment the following once 
  
  # adding empty line to csv files
  #cat("\n", file = meas[i], append = TRUE)
  
  
  # creating temporary dataset for file
  temp_dataset <- read.csv(meas[i], header = TRUE) 
  temp_dataset <- temp_dataset[1:4] # forcing file to only have 4 columns. When measurements are exported into excel, Adobe sometimes adds an additional column
  temp_dataset$file <- filename[i] # adding new colum "file" to temp data set
  dataset <- rbind(dataset,
                   temp_dataset) #row binding from temp dataset to full dataset
  rm(temp_dataset) #clearing out temp data for next file
}

# time for code to run?
end_time <- Sys.time()
end_time - start_time
# jpz's work HP= 51 seconds


# check data for errors
dataset
names(dataset)

dataset %>%
  filter(Unit == "pt")

dataset <- dataset %>% 
  separate(file, c("site", "rep", "date", "photo_number"),
           sep = "_", remove = FALSE)

dataset <- dataset %>%
  mutate(
    date = dmy(date),
    year = year(date),
    month = month(date))

distinct(dataset,
         site, 
         date)


sort(unique(dataset$Label))
# lots of typos to fix

dataset %>%
  filter(Label == "",
         Unit == "mm") %>%
  distinct(file)

dataset %>%
  filter(Label == "35.73") %>%
  distinct(file)

### length-weight

# read in csv with length weight equations
lw_coef <- read.csv("data/LW_coeffs.csv")

# which taxa don't have lw coefs?
no_coef <- setdiff(unique(dataset$Label), unique(lw_coef$taxon))

sort(no_coef)

sort(unique(lw_coef$taxon))
