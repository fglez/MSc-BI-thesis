###ANALYSIS AND PLOT OF THE CLUSTERING RESULLT OF BiG-SCAPE

#Set path
path <- "/home/mfcg/Descargas/thesis-seq/bgcs/mtsmr/clstr"
setwd(path)

## VENN DIAGRAM FOR EACH BGC CLASS 
#Requirements
#Libraries
library(nVennR)
#Objects
#Read a list with the sample names, the same used in dvd_clstr_tbl.sh script
smplls <- read.delim(file = "sample.ls", header = FALSE)
#Number of samples
nosam <- 7
#Empty list to save the sst_ll function result
flist <- vector("list", nosam) #Establish its length equal to the sample number
names(flist) <- smplls$V1 #Name each element of the list with the samples in the same order
#Group the sample separated tables in a list of lists
sst_ll <- function(cttbl) {
  for(i in 1:nosam){
    clist <- read.delim(file = cttbl, header = FALSE) #Read the list of the separated tables 
    cuted <- read.delim(file = clist[i,1], header = FALSE) #Read a sample clustering table
    cuted <- as.character(cuted$V2) #Extract the bgcs families present in each sample
    flist[[i]] <- cuted #Save the present families it the corresponding sample element
  }
  return(flist)
}

#nvenn diagrams
#Repeat the plotVenn function on the sets data until obtaining the best result
#Save in a svg file
#NRPS
nrpsll <- sst_ll("NRPScttbl.ls")
nrpsnv<- plotVenn(nrpsll) 
nrpsnv<- plotVenn(nVennObj = nrpsnv)
nrpsnv<- plotVenn(nVennObj = nrpsnv)
showSVG(nVennObj = nrpsnv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "nrpsvenn.svg", 
        fontScale = 2)
#Others
othrll <- sst_ll("OTHERScttbl.ls")
othrnv <- plotVenn(othrll)
othrnv <- plotVenn(nVennObj = othrnv)
showSVG(nVennObj = othrnv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "othrvenn.svg", 
        fontScale = 2)
#RiPPs
rippll <- sst_ll("RIPPcttbl.ls")
rippnv <- plotVenn(rippll)
rippnv <- plotVenn(nVennObj = rippnv)
showSVG(nVennObj = rippnv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "rippvenn.svg", 
        fontScale = 2)
#PKSother
pksoll <- sst_ll("PKSOcttbl.ls")
pksonv <- plotVenn(pksoll)
pksonv <- plotVenn(nVennObj = pksonv)
pksonv <- plotVenn(nVennObj = pksonv)
showSVG(nVennObj = pksonv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "pksovenn.svg", 
        fontScale = 2)
#Terpenes
trpnll <- sst_ll("TRPNcttbl.ls")
trpnnv <- plotVenn(trpnll)
trpnnv <- plotVenn(nVennObj = trpnnv)
trpnnv <- plotVenn(nVennObj = trpnnv)
showSVG(nVennObj = trpnnv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "trpnvenn.svg", 
        fontScale = 2)
#Before running the function in these cases the cttbl was modified
#removing all the empty files and the number of samples was also modified to 2
nosam <- 2
flist <- vector("list", nosam) #Establish its length equal to the sample number
names(flist) <- c("eue", "eub") #Name each element of the list with the samples in the same order
#PKSI
pksill <- sst_ll("PKSIcttbl.ls")
pksinv <- plotVenn(pksill)
showSVG(nVennObj = pksinv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "pksivenn.svg", 
        fontScale = 2)
#PKS-NRPS Hybrid
hbrdll <- sst_ll("HYBRIDcttbl.ls")
hbrdnv <- plotVenn(hbrdll)
showSVG(nVennObj = hbrdnv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "hbrdvenn.svg", 
        fontScale = 2)
#All BGCs classes
allll <- sst_ll("ALLcttbl.ls")
allnv <- plotVenn(allll)
allnv <- plotVenn(nVennObj = allnv)
allnv <- plotVenn(nVennObj = allnv)
showSVG(nVennObj = allnv, 
        opacity = 0.2, 
        borderWidth = 4, 
        outFile =  "alllvenn.svg", 
        fontScale = 2)
#Only samples from the direct samples
diral <- allll[c(5, 6, 7)]
dirnv <- plotVenn(diral)
showSVG(nVennObj = dirnv,
        opacity = 0.2,
        borderWidth = 4,
        outFile = "diralvenn.svg",
        labelRegions = F,
        fontScale = 2)
#Only samples from the co-cultures
culal <- allll[c(1, 2, 3, 4)]
culnv <- plotVenn(culal)
showSVG(nVennObj = culnv,
        opacity = 0.2,
        borderWidth = 4,
        outFile = "culalvenn.svg",
        labelRegions = F,
        fontScale = 2)
#Grouping by microbiome source
metll <- list(direct = c(unlist(diral)), coculture =  c(unlist(culal)))
metnv <- plotVenn(metll)
showSVG(nVennObj = metnv,
        opacity = 0.2,
        borderWidth = 4,
        outFile = "metalvenn.svg",
        labelRegions = F,
        fontScale = 2)

##CONTIGENCY TABLE OF THE BGC FAMILIES PRESENT IN EACH SAMPLE
#Empty data frame to be used in the dd_n_pst function
itbl <- data.frame()
#Column names for the data frame obtained with dd_n_pst function
clmns <- c("Sample", "BGCname","Family")
#Function to add the sample name in the BiG-SCAPE clustering table
#and to paste together each separated sample table
dd_n_pst <- function(cttbl){ #The list of the separated tables obtained with dvd_clstr_tbl.sh
  for(i in 1:nosam){
    clist <- read.delim(file = cttbl, header = FALSE) #Read the list of the separated tables 
    cuted <- read.delim(file = clist[i,1], header = FALSE) #Read a ample clustering table
    nmbr <- dim(cuted)[1] #Number of BGCs foud in each sample, equivalent to row number
    rptnm <- rep(smplls[i,1], nmbr) #Repeat the sample name 
    ddtbl <- cbind(rptnm, cuted) #Paste the sample name column to its clustering table
    colnames(ddtbl) <- clmns #Name the tables
    itbl <- rbind(itbl, ddtbl) #Paste together the sample clustering table
  }
  return(itbl) #Return all the sample clustering tables pasted together
}

#Reformat table of the BGCs present in each samples for a given class 
#with the sample names
#NRPS
nrpsdnp <- dd_n_pst("NRPScttbl.ls")
#Others
othrdnp <- dd_n_pst("OTHERScttbl.ls")
#RiPPs
rippdnp <- dd_n_pst("RIPPcttbl.ls")
#PKS other
pksodnp <- dd_n_pst("PKSOcttbl.ls")
#Terpenes
trpndnp <- dd_n_pst("TRPNcttbl.ls")
#Change the number of samples for these cases before running
nosam <- 2
#PKSI
pksidnp <- dd_n_pst("PKSIcttbl.ls")
#PKS-NRPS Hybrid
hbrddnp <- dd_n_pst("HYBRIDcttbl.ls")

