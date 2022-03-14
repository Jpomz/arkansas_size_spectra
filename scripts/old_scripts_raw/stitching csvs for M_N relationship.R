## Waiuta Mine
##=======================================##
#### reading in .csv measurement files ####
##=======================================##
rm(list=ls()) #start by wiping the slate clean

upperDir <- "C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\Invert CSVs\\Waiuta Test\\"
#Path to all of my csv's
meas <- list.files(upperDir, full.names=T)
#meas2 <- list.files(upperDir, pattern = "*.csv", full.names = T) #code for only reading in files with ".csv" ending
#list of all of the paths to csv's within folder
filename <- list.files(upperDir)
#just the file name of each csv

#### Loop Reading in csv's ####

dataset <- NULL #Empty dataset
for (i in 1:length(meas)){ #length(meas) = number of total files within folder
    temp_dataset <-read.csv(meas[i], header=T) #creating temporary dataset for file
    temp_dataset <- temp_dataset[1:4] #forcing file to only have 4 columns. When I export my measurements into excel, Adobe sometimes adds an additional column
    temp_dataset$file <- filename[i] #adding new colum "file" to temp data set
    dataset<-rbind(dataset, temp_dataset[,1:5]) #row binding 5 columns from temp dataset to full dataset
    rm(temp_dataset) #clearing out temp data for next file
}


head(dataset)
tail(dataset)
dim(dataset)

#alternative code using lapply and dplyr to stitch measurements
#more straight forward, fewer lines, etc
# library(dplyr)
# setwd("C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\Invert CSVs\\Waiuta 
#       Test\\")
# temp <- list.files(pattern = "*.csv")
# files <- lapply(temp, read.csv)
# stitch <- ldply(meas, function (x){rbind(x)})

#### remove "scale" measurements ####
dataset1 <- dataset[dataset$Label != "",]#"blank" labels = scale measurements in Adobe
## Remove "Pupae" and "Unknown" measurements
dataset1 <- dataset1[dataset1$Label != "Chironomidae Pupae",]
dataset1 <- dataset1[dataset1$Label != "Unknown",]

##===================================##
####Removing data from file column ####
##===================================##

#### !!!double check character numbers when switching to other measurement files!!! ####

dataset1$site <- substr(dataset1$file, start = 1, stop = 3)#pulls out site code from file column
dataset1$surber <- substr(dataset1$file, start = 5, stop = 6)#pulls out Surber number
dataset1$date <- substr(dataset1$file, start = 9, stop = 16)#pulls out date

head(dataset1)
dim(dataset1)



dataset1$Label <- gsub(pattern = "Oligachaeta", replacement = "Oligochaeta", x = dataset1$Label)
dataset1$Label <- gsub(pattern = "Oligochaete", replacement = "Oligochaeta", x = dataset1$Label)
dataset1$Label <- gsub(pattern = "Oligachaeta ", replacement = "Oligochaeta", x = dataset1$Label)
dataset1$Label <- gsub(pattern = "Oligachaete", replacement = "Oligochaeta", x = dataset1$Label)
dataset1$Label <- gsub(pattern = "\\bOrthoclad\\b", replacement = "Orthocladinae", x = dataset1$Label)
dataset1$Label <- gsub(pattern = "Elmidae Larvae", replacement = "Elmidae", x = dataset1$Label)

##===============##
#### Write.csv ####
##===============##
### may be alternative way to replace typos, misspellings etc ###
## make csv and then do a find and replace in excel
setwd("C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\")
write.csv(dataset1, file = "Waiuta.Stitched.csv", row.names = F)
list.files() #shows all files
### at this time writes files to my "Documents" folder
### if using this method, name file something meanigful and find it in my "Documents" immediately and move to appropriate place





##=============================================##
#### adding length weight regression formulas ####
##==============================================##
dataset2 <- read.csv("C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\Waiuta.Stitched.csv")#read in stitched data set

range(dataset2$Value)
#removing measurements < 1mm
#dataset2 <- dataset2[dataset2$Value>=1,]
#ran analysis with this, and appears
#to give better results is leave these obs in. 

formula <- read.csv("C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\Files from Helen\\invert_length_weight_JP_modified.csv")#file from Helen containing all variable values
head(formula)

## Merging files ####
#Merging dataset1 with .csv containing LW regression formulas
#merges by matching "Label" columns (taxa names)

#need to add to Helens formula file
#Probably actually want to do this AFTER fixing typos, etc (see gsub() below)

dataset3 <- merge(dataset2, formula, by.x = "Label", by.y = "Label", all.x = T)
#all.x = T includes taxa from my files that do not have corresponding formulas in Helens file
head(dataset3)
tail(dataset3)
dim(dataset3)


#adding dryweight (dw) column
#Rememer, log() is default base e, opposite is exp()
#log10() = log base 10  
dataset3$dw <- exp(dataset3$ln_a+dataset3$b*log(dataset3$Value))
write.csv(dataset3, file = "Waiuta.estDW.csv", row.names = F)

#dataset3 <- read.csv("C:\\Users\\Justin\\Documents\\Data\\Length DW conversion\\Waiuta.estDW.csv")


##================================##
##increasing samples to add fish####
##================================##

dataset3$unique <- 1:dim(dataset3)[1]#adding unique number to each row

dataset3$newsite <-NA #adding "newsite" to combine Wa1+2 and Wa3+4
dataset3[dataset3$site=="Wa1"|dataset3$site=="Wa2","newsite"]<-"Ref"
dataset3[dataset3$site=="Wa3"|dataset3$site=="Wa4","newsite"]<-"Imp"
dataset3[dataset3$site=="Wa5","newsite"]<-"Rec"

ref <- subset(dataset3, newsite=="Ref") #subsetting each site type
imp <- subset(dataset3, newsite=="Imp")
rec <- subset(dataset3, newsite=="Rec")

#ref mean w = 1.25m * 20m = 25m2
#imp mean w = .563m * 20m = 11.26m2
#rec mean w = 1.71m * 20m = 34.2m2
#0.06*6 = 0.36 ##ratio for ref and imp sites (6 surbers)
#0.06*3 = 0.18  ## ratio for rec site (3 surbers)


ref_ratio <- 25/0.36 #total m2 / m2 sampled in cumulative surbers
imp_ratio <- 11.26 / 0.36 
rec_ratio <- 34.2 / 0.18
  
ref_expand <- ref[rep(row.names(ref), ref_ratio),] #multiplying data rows by area ratio
imp_expand <- imp[rep(row.names(imp), imp_ratio),]
rec_expand <- rec[rep(row.names(rec), rec_ratio),]

expand <- rbind(ref_expand, imp_expand, rec_expand) #combining all expanded data sets
#need to add fish/ Koura data to this and then write to .csv 
write.csv(expand, file = "Expanded_Waiuta_estdw.csv", row.names = F) 


dim(expand)
aggregate(dw~unique, expand, FUN=length)


  
##===========================##
#### Random code fragments ####
##===========================##

### code notes from Helens expanding data set
cb_ratio<-21.58676 #this was the ratio of in this case baskets to channel. In your case will be how many surbers fit in your reach
invert_only$ID_orig<-1:dim(invert_only)[1]

invert_expand<-dataset3[rep(row.names(dataset3), cb_ratio),]  #this should replicate up all your data to be at the reach level
head(invert_expand)

#you can sort of check it worked using aggregate
sp <- aggregate(dw~Label+site, invert_expand, FUN=length)


##============##
## scrap code ##
##============##

#### potential if then statement for lw rergression based on base e or 10 ####
#need to test once I have measurements for base 10 taxa (e.g. Oligachaete)
#dataset2$dw <- 
#if dataset2$base == 10{
#  10^(dataset2$ln_a+dataset2$b*log10(dataset2$Value))
#}else{
#  exp(dataset2$ln_a+dataset2$b*log(dataset2$Value)
#}

dataset2[dataset2$base=="10",]
dataset2[dataset2$Label=="Deleatidium",]

#### gsub () Code to rename taxa in label column ####

#### Scrap code ####

###Use below code to replace typos, misidentifications, etc. 

###Useful if you know exactly what the typo is. Or if you want to change all of some taxa into a higher order. e.g. "Deleatidium" --> "Ephemeroptera"

factor(dataset2$Label) #may be useful to show what factors of $Label are present
levels(factor(dataset2$Label))

###need to revise code here to replace whole label. e.g. when trying to replace "Orthoclad" it reaplces both "Orthoclad" and "Orthocladniae" entries bc of shared begining. 
dataset2$Label <- gsub(pattern = "Oligachaeta", replacement = "Oligochaeta", x = dataset2$Label)

dataset2$Label <- gsub(pattern = "Deleatidium", replacement = "Funky butt", x = dataset2$Label) #where label == Deleatidium, gsub replaces with "Funky butt" 
dataset2[dataset2$Label=="Funky butt",] #returns a sunset of the data
dataset2[dataset2$Label=="Deleatidium",] #0 returns

### trying to find all items that begin with "XXX" and replace them with same name
### e.g. Pycnocentroides, Pycnocentrella, Pycnocentria and replace them all with "Pycnocentrella"
### so far it is NOT working. 
dataset2$Label <- gsub(pattern = "Pycno", replacement = "ALOHA", x = dataset2$Label)## only replaces "Pycno" with "ALOHA", not whole word
dataset2[dataset2$Label=="Pycnocentrella",]
dataset2[dataset2$Label=="Pycno",]


