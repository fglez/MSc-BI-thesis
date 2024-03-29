##FIND GOOD BINS

#The following script helps you to do a fast identification 
#of good bins to use as MAGS
#It uses the small new file produced by the script 4mag-tocsv.sh
#It also list these good bins in a list, goodbins.ls

#Required libraries
library(stringr)
library(dplyr)

#Path to your file
samnam <- "PhE"
path <- paste0('/home/mfcg/Descargas/thesis-seq/cult2021/', samnam) #Change this to your path
setwd(path)

#Read file
filnam <- paste0("small", samnam, "check.csv")
smfile <- read.csv(filnam, header = FALSE) #Read the file

#Format the data frame smfile

#Remove words to leave only values
torm <- c("NOTHING", "'Completeness': ", "'Contamination': ") #words to remove must be equal in size to the number of columns
valdf <- t(apply(smfile, 1, str_remove_all, torm)) #Apply by row the str_remove_all function and using t to re-arrange columns and rows

#Name the columns
titulos <- c("bin", "completeness", "contamination") #column names
colnames(valdf) <- titulos

#Make numeric the completeness and contamination columns
valdf <- as.data.frame(valdf)
valdf[ ,2:3] <- apply(valdf[ ,2:3], 2, as.numeric) #Apply by column the as.numeric function
valdf

#Identify the quality bins
comp <- 80 #Mininum completeness, change it as it fits for you
cont <- 10 #Maximum contamination allowed, choose a treshold
good <- valdf %>% 
  filter(completeness > comp) %>% #Leave only the bins above the completeness treshold
  filter(contamination < cont) #Leave only the bins below the contamination treshold
good
goodbins <- good$bin #Retrieve the names of the bins with the desired quality
goodbins
write(goodbins, file = "goodbins.ls") #Make a txt file with the bins' names

