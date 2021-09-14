###TRANSFORM TO PHYLOSEQ AND MERGE###

#The following script shows how to transform to phyloseq objects, 
#the tables produced by the abundance.sh script
#The abundace.sh script can be found here https://github.com/carpentries-incubator/metagenomics/blob/gh-pages/files/abundance.sh

#This script is a modified version from the one found here https://carpentries-incubator.github.io/metagenomics/08-automating_abundance/index.html
#if you wish to learn more, please visit the link above

#This script also shows how to merge the phyloseq object produced by different samples, into one.

#Finally it helps you to prepare your phyloseq object for the following analyses by leaving only bacterial reads, 
#samples with a similar read depth and normalized data

#LIBRARIES
library("readr")
library("phyloseq")
#library("edgeR")


#PATH AND PREFIX TO THE SAMPLE TABLES
sampname <- "PhE"
sampath <- paste0("/home/mfcg/Descargas/thesis-seq/taxonomy/kraken/results", sampname) #write the path to the directory that contains the files produced by the abundance.sh script
setwd(sampath) #set the path as the working directory

#FROM ABUNDANCE.SH TABLES TO PHYLOSEQ OBJECT

psmetagenome <- list() #make an empty list to save the phyloseq object that will be produced by the function ab_to_ps

#The following function automates the steps to transform from the tables done by abundace.sh to a phyloseq object
ab_to_ps <- function(sampname){ #this function uses the sampname object to work, sampname = sampname
  #File names
  rankedwc <- paste0(sampname, ".ranked-wc") #add the name of your ranked-wc file 
  ltablewc <- paste0(sampname, ".lineage_table-wc") #add the name of your lineage_table-wc file 
  
  #Read tables
  OTUS <- read_delim(rankedwc,"\t", escape_double = FALSE, trim_ws = TRUE)
  TAXAS <- read_delim(ltablewc, "\t", escape_double = FALSE,
                      col_types = cols(subspecies = col_character(),
                                       subspecies_2 = col_character()), trim_ws = TRUE)
  #Matrix format
  abundance <- as.matrix(OTUS[ , -1]) #To avoid that the OTU column is taken as a sample the first column is omitted
  lineages <- as.matrix(TAXAS[ , -1]) #To avoid that the OTU column is repeated the first column is omitted
  
  #Retrieving the OTUs' identity as rownames using the the original file as reference
  row.names(abundance) <- OTUS$OTU 
  row.names(lineages) <- TAXAS$OTU
  
  #Phyloseq format
  OTU <- otu_table(abundance, taxa_are_rows = TRUE)
  TAX <- tax_table(lineages)
  psmetagenome <<- phyloseq(OTU, TAX)
  
}

ab_to_ps(sampname = sampname)

#Save phyloseq object to a rds file
filename <- paste0(sampname, ".rds")
write_rds(psmetagenome, file = filename) 

#PHYLOSEQ OBJECT WITH ALL THE SAMPLES
#To make the work easier move all the phyloseq objects that you will merge to the same directory
#The directory must be exclusive for the sample phyloseq objects
psfpath <- "/home/mfcg/Descargas/thesis-seq/taxonomy/phyloseq/indiv" #write the path to the directory that contain the phyloseq objects
setwd(psfpath)
all_samp_ps <- list.files(path = ".", pattern = ".rds") #make a character vector with the names of the phyloseq objects
merged_all_ps <- list()
for(i in 1:length(all_samp_ps)){
  nps <- read_rds(all_samp_ps[i])
  merged_all_ps <<- merge_phyloseq(merged_all_ps, nps)
}

#Save phyloseq object with all the taxa present
write_rds(merged_all_ps, file="allphyla.rds") #The merged phyloseq object can be saved as it is

##PHYLOSEQ OBJECT WITH ONLY BACTERIA TAXA

#Some of these steps were taken from the R script of Diego Garfias which can be found here https://github.com/nselem/mg-covid/blob/master/diego/taxonomy/analysis-covid-240521.R 
bacteria_ps <- subset_taxa(merged_all_ps, superkingdom == "Bacteria")  #Only Bacteria OTUs are left in the merged file

leaveout <- subset_taxa(bacteria_ps, family != "mitocondria" & class != "Chloroplast") #Subset from the bacteria object everything that is not a mitochondrial or chloroplast sequence
bacteria_ps <- prune_taxa(taxa_names(leaveout), bacteria_ps) #Leave only the reads of bacteria without chloroplast and mitochonria ones

write_rds(bacteria_ps, file = "bacthesis.rds") #Save it to the rds file

#To avoid mixing data the rds file with all the samples bacteria taxa was moved to a new directory
abpath <- "/home/mfcg/Descargas/thesis-seq/taxonomy/phyloseq/allbact"
setwd(abpath)
#Read the rds file with the bacteria phyloseq
bacteria_ps <- read_rds("bacthesis.rds") #Read the rds file with the phyloseq object
#Retrieve the tables from the phyloseq object that only contains bacterial OTUs, bacteria_ps
otustab <- otu_table(bacteria_ps) #Table with the OTU abundance matrix
taxatab <- tax_table(bacteria_ps) #Table that contain the taxonomic classification of the OTUs
#Read the metadata table
samptab <- read.csv("metadata-bt.csv", header = TRUE, row.names = 1) #Read a csv file containing a table with the sample data or metadata
#Ensure that the sample names are the same, if they are different follow the following steps
#Change the name of the OTUtable columns, that correspond to the sample names, they are already sorted due to the list.files function
#Sort in the same manner, alphabetically, the SAMtable row or sample names, and use the to replace the names in the OTUtable
#colnames(otustab) <- sort(rownames(samptab), decreasing = FALSE)

#Format the table with the sample data to phyloseq format
samptab <- sample_data(samptab)

#Update the phyloseq object with the metadata
bacteria_ps <- phyloseq(taxatab, otustab, samptab) #Add the metadata from the SAMtable and update the names from the OTUtable in the new bacteria pss
#Save changes
write_rds(bacteria_ps, file = "bacthesis.rds") #This file is ready to be used in the analysis done in phylobject.R
write.csv(otustab, file = "otab-bt.csv") #Save a the OTU table in a CSV
write.csv(taxatab, file = "ttab-bt.csv") #Save a the taxa table in a CSV

##NORMALIZATION
#The normalization will be done with ANCOM-BC

abpath <- "/home/mfcg/Descargas/thesis-seq/taxonomy/phyloseq/allbact" #Set path to your bacteria phyloseq
setwd(abpath) #Set path to your bacteria phyloseq
#Retrieve your bacteria phyloseq files 
bacteria_ps <- read_rds("bacthesis.rds")
samptab <- sample_data(bacteria_ps)
taxatab <-tax_table(bacteria_ps)

#Required libraries
library("ANCOMBC")
library("microbiome")

#Group OTUS to a taxa rank
phyla_ps <- aggregate_taxa(bacteria_ps, "phylum") #Phyla level
genera_ps <- aggregate_taxa(bacteria_ps, "genus") #Genera level

#Using ANCOM BC to normalize
#I used year as this primarily separates the two methods on how the samples were processed method could be used as well

#Normalized OTU data
onormps <- ancombc(phyloseq = bacteria_ps, formula = "Year", 
                   p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000000,
                   group = "Year", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                   max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

nottab <- onormps$feature_table #Save the OTU normalized data as a table
nottab <- otu_table(nottab, taxa_are_rows = TRUE)
normbact <- phyloseq(taxatab, nottab, samptab)
write_rds(normbact, file = "normbact.rds")
