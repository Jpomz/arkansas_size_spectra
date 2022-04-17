# data prep

upperDir <- "data/raw_csv" 
# vector of all of the paths to csv's within folder
meas <- list.files(upperDir, full.names=TRUE) 
# just the file name of each csv
filename <- list.files(upperDir) 



#### Loop Reading in csv's ####

# note, as of 4/17/22, this dataset is getting large
# Takes a while to read it all in
# hopefully don't need to repeat this, or could optimize it if needed later. 

dataset <- NULL #Empty dataset

for (i in 1:length(meas)){ 
  #length(meas) = number of total files within folder
  # adding empty line to csv files
  cat("\n", file = meas[i], append = TRUE)
  # creating temporary dataset for file
  temp_dataset <- read.csv(meas[i],
                           header=T) 
  temp_dataset <- temp_dataset[1:4] # forcing file to only have 4 columns. When measurements are exported into excel, Adobe sometimes adds an additional column
  temp_dataset$file <- filename[i] # adding new colum "file" to temp data set
  dataset <- rbind(dataset,
                   temp_dataset) #row binding from temp dataset to full dataset
  rm(temp_dataset) #clearing out temp data for next file
}

# column names of dataset we just stitched together
names(dataset)

dim(dataset)

# preview the first 5 rows of the data
dataset[1:5,1:4]

# look at all of the Unit values - make sure they're all mm
unique(dataset$Unit)

# look at all of the taxa names in the data
# taxa names are in the "Label" column
sort(unique(dataset$Label))


### Cleaning up some names ###

# there is probably a better way to do this...

# start here
dataset$Label <- gsub(
  pattern = "Brachycentriidae",
  replacement = "Bracycentridae",
  x = dataset$Label)


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
  pattern = "Flat Worm",
  replacement = "Planariidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Flat_Worm",
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

# 2012
dataset$Label <- gsub(
  pattern = "Blepharicid",
  replacement = "Blephariceridae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Ceratropgonidae",
  replacement = "Ceratopogonidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Emphididae",
  replacement = "Empididae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Glossomatidae",
  replacement = "Glossosomatidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Heptagenidae",
  replacement = "Heptageniidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Rhyancophilidae",
  replacement = "Rhyacophilidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Uenoidea",
  replacement = "Uenoidae",
  x = dataset$Label)

dataset$Label <- gsub(
  pattern = "Simulidae",
  replacement = "Simuliidae",
  x = dataset$Label)

# check the updated dataset
# make sure taxa names are correct
sort(unique(dataset$Label))

# check some things out
# No label
dataset[dataset$Label=="",]

# chironomidae adults, how many?
nrow(dataset[dataset$Label=="Chironomidae(adult)",]) + 
  nrow(dataset[dataset$Label=="Chironomidae_(adult)",])
# 24

# horse hair worms, how many?
nrow(dataset[dataset$Label=="Nematomorpha",])


# Elmidae adults, how many?
nrow(dataset[dataset$Label=="Elmidae (adult)",]) + 
  nrow(dataset[dataset$Label=="Elmidae_Adult",])
# 6

# Elmidae adults, how many?
nrow(dataset[dataset$Label=="Perlidae(adult)",])
# 1


# save the data set for further analysis
saveRDS(dataset, "data/ark_csv_stitched.R")
