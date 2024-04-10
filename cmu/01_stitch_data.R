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
# jpz's work HP =  seconds
# hope thinkPad = 1.99 seconds


# check data for errors
dataset
names(dataset)

# measurement units ####
# all units
dataset %>%
  distinct(Unit)
# "pt"? Pixels or something?
dataset %>%
  filter(Unit == "pt")


dataset %>%
  filter(Unit == "in") %>%
  group_by(file) %>%
  count() %>%
  ungroup() %>%
  distinct(n)
  


dataset <- dataset %>% 
  as_tibble() %>%
  separate(file, c("site", "rep", "date", "photo_number"),
           sep = "_", remove = FALSE)

dataset

dataset <- dataset %>%
  mutate(
    date = dmy(date),
    year = year(date),
    month = month(date))

distinct(dataset,
         site, 
         date)

# label typos ####
sort(unique(dataset$Label))
# lots of typos to fix

dataset %>%
  filter(Label == "",
         Unit == "mm") %>%
  distinct(file)

dataset %>%
  filter(Label == "35.73") %>%
  distinct(file)

dataset %>%
  filter(Label == "Hydrophilidae?") %>%
  distinct(file)
# only one, it is a hydroptilidae

dataset %>%
  filter(Label == "Lymnaeidae") %>%
  distinct(file)
# only one, No LW regression, skipping for now (2024-04-09)

# fix Label Typos ####

sort(unique(dataset$Label))

dataset %>%
  mutate(Label = case_when(
    Label == "Akari (aquatic mites)" ~ "Acari",
    Label == "Chironimidae" ~ "Chironomidae",
    Label == "Chloroperilidae" ~"Chloroperlidae",
    Label == "Elmidae (adult)" ~ NA,          
    Label == "Elmidae (larvae)" ~ "Elmidae",
    Label == "Ephemeriliidae" ~ "Ephemerellidae",
    Label == "Hydrophilidae?" ~ "Hydroptilidae",
    Label == "Hydroptilidae (pupae)" ~ NA,
    Label == "Hydropyschidae" ~ "Hydropsychidae",
    Label == "Hydrosychidae" ~ "Hydropsychidae",
    Label == "Nematomorpha " ~ "Nematomorpha", 
    Label == "Nemotomorpha" ~ "Nematomorpha", 
    Label == "Platyhelminthes" ~ "Turbellaria",
    Label == "Platyhelminthes: Turbellaria" ~ "Turbellaria",
    Label == "Pteronarchyidae"~ "Pteronarcyidae",
    Label == "Pteronarcidae"~ "Pteronarcyidae",
    Label == "Pteronarycidae"~ "Pteronarcyidae",
    Label == "Rhyachophilidae" ~ "Rhyacophilidae",
    Label == "Rhyacophildae" ~ "Rhyacophilidae",
    Label == "Simulidae" ~ "Simuliidae",
    Label == "Water Mite " ~ "Acari",
    Label == "water mites" ~ "Acari",
    Label == "Water mites"  ~ "Acari",
    .default = Label
  )) %>%
  distinct(Label) %>%
  arrange(Label) %>%
  pull()

dat_clean <- dataset %>%
  mutate(Label = case_when(
    Label == "Akari (aquatic mites)" ~ "Acari",
    Label == "Chironimidae" ~ "Chironomidae",
    Label == "Chloroperilidae" ~"Chloroperlidae",
    Label == "Elmidae (adult)" ~ NA,          
    Label == "Elmidae (larvae)" ~ "Elmidae",
    Label == "Ephemeriliidae" ~ "Ephemerellidae",
    Label == "Hydrophilidae?" ~ "Hydroptilidae",
    Label == "Hydroptilidae (pupae)" ~ NA,
    Label == "Hydropyschidae" ~ "Hydropsychidae",
    Label == "Hydrosychidae" ~ "Hydropsychidae",
    Label == "Nematomorpha " ~ "Nematomorpha", 
    Label == "Nemotomorpha" ~ "Nematomorpha", 
    Label == "Platyhelminthes" ~ "Turbellaria",
    Label == "Platyhelminthes: Turbellaria" ~ "Turbellaria",
    Label == "Pteronarchyidae"~ "Pteronarcyidae",
    Label == "Pteronarcidae"~ "Pteronarcyidae",
    Label == "Pteronarycidae"~ "Pteronarcyidae",
    Label == "Rhyachophilidae" ~ "Rhyacophilidae",
    Label == "Rhyacophildae" ~ "Rhyacophilidae",
    Label == "Simulidae" ~ "Simuliidae",
    Label == "Water Mite " ~ "Acari",
    Label == "water mites" ~ "Acari",
    Label == "Water mites"  ~ "Acari",
    .default = Label
  )) 

# read in csv with length weight equations
lw_coef <- read.csv("data/LW_coeffs.csv")

# which taxa don't have lw coefs?
no_coef <- setdiff(unique(dat_clean$Label), unique(lw_coef$taxon))

sort(no_coef)

# Using surrogate equations for these taxa
# Morley et al. 2020 Plos One
# Psychodidae --> Empididae
# Uenoidae --> Limnephilidae
# Hydroptilidae --> Rhyacophilidae
dat_clean <- dat_clean %>%
  mutate(Label = case_when(
    Label == "Psychodidae" ~ "Empididae",
    Label == "Uenoidae" ~ "Limnephilidae",
    Label == "Hydroptilidae" ~"Rhyacophilidae",
    .default = Label
  )) 

no_coef2 <- setdiff(unique(dat_clean$Label), unique(lw_coef$taxon))

sort(no_coef2)

# percent of data that does not have lw coeffs?
(nrow(dat_clean[dat_clean$Label %in% no_coef2,])/ nrow(dat_clean))*100
# ~ 9.58% of data does not have lw regression equations

# estimate dry weight (dw) ####

# make string vector of column names we care about
lw_cols <- c("taxon",
             "a",
             "b",
             "L_units",
             "dw_units",
             "formula_type",
             "formula")

# check formulas in LW for taxa in arkansas
equations <- lw_coef[lw_coef$taxon %in% unique(dat_clean$Label),]


# join the ark data with lw coef values
fulldat <- merge(dat_clean, lw_coef[,lw_cols], 
                 by.x = "Label",
                 by.y = "taxon", all.x = FALSE)

# dimensions of the ark and merged data
dim(dataset)
dim(fulldat)

# different lw formulas
unique(fulldat$formula_type)

# estimate dw based on formula type
fulldat <- fulldat %>% 
  mutate(dw = case_when(
    formula_type == 1 ~ a * Value^b,
    formula_type == 2 ~ exp(a + b * log(Value))))


names(fulldat)
# remove any dry weight values < or = to 0
dw <- fulldat %>%
  filter(dw >0) %>%
  select(site, rep, year, month, file, dw)

nrow(fulldat)
nrow(dw)

names(dw)

dw

# filter out data < or = 0.0026 mg
# See SI from Perkins et al. 2018 Ecology Letters
dw <- dw %>%
  filter(dw>0.0026)

min(dw$dw)

# save data with estimated dry weights
saveRDS(dw, "cmu/derived_data/cmu_ark_dw.RDS")
