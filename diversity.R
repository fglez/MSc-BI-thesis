###BACTERIAL DIVERSITY ANALYSIS OF MASTER THESIS
###THE ROLES OF THE MICROBIOME OF THE CYCAD ESPECIALIZED INSECTS

##REQUIRED LIBRARIES
library("readr")
library("phyloseq")
library("ggplot2")
library("microbiome")
library("vegan")
library("dplyr")

##FILES AND OBJECTS

#Open the files with the phyloseq object that contains the normalized bacteria data of all the samples
iepath <- "/home/mfcg/Descargas/thesis-seq/taxonomy/phyloseq/allbact" #Set path to your bacteria phyloseq
setwd(iepath) #Set path to your bacteria phyloseq
normbact <- read_rds("normbact.rds")

#Extract metadata as an object
metadata <- sample_data(normbact)
vrbls <- colnames(metadata)

#Separate the samples by host and source
EuC <- prune_samples(c("EuB", "EuE"), normbact) #Samples of cultures from E. childrenae guts
RhC <- prune_samples(c("PhB", "PhE"), normbact) #Samples of cultutres from Rhopalotria sp. tissues
DEu <- prune_samples(c("S33", "S34", "S35"), normbact) #Samples directly from E. childrenae guts

#Abundance plots
#The following function will help you to plot a relative abundance plot
#to a taxa level you wish to explore, showing only the top 15 taxa in it
relabplot <- function(txrnk,gsmpl){ #Write the taxa rank you wish to explore, txrnk, and the phyloseq object to use, gsmpl
  aggps <- aggregate_taxa(gsmpl, txrnk, verbose = FALSE) #Aggregate normbact by taxa rank
  topnam <- names(sort(taxa_sums(aggps), TRUE)[1:15]) #The name of the ten most abundant phylum in that taxa rank
  topps <-  prune_taxa(topnam, aggps) #Select only the most abundant taxa
  prcntgs <- transform_sample_counts(topps, function(x) x*100 / sum (x)) #Transform your top ten phyloseq to percentages
  perc_plot <- plot_bar(prcntgs, fill = txrnk) + #Plot the information with the percentage data by chosen taxa rank
    geom_bar(aes(), stat = "identity", position = "stack") +
    theme_bw()
  return(perc_plot)}
#Relative abundace plots by taxa rank
EucRab <- relabplot("genus", EuC)
RhcRab <- relabplot("genus", RhC)
DEuRab <- relabplot("genus", DEu)

##ALPHA DIVERSITY
#Alpha diversity plot
alpha_plot <- plot_richness(normbact, 
                            measures = c("Chao1", "Shannon", "Simpson"), #Indexes used 
                            color = "Host", #Color by insect specie and stage
                            shape = "Method") + theme_bw() #Shape showing if the samples was direct, cultured in BG110 with or without cycad extract

#TABLE 9. ALPHA DIVERSITY MEASURES FOR EACH MICROBIOME SAMPLE
riqueza <- alpha_plot$data #Save the richness measures in a data frame
write.csv(riqueza, file = "alphadiv.csv") #Save it as a file

#Separate the alpha diversity measure by their index
chao1 <- riqueza %>% filter(variable == "Chao1")
shannon <- riqueza %>% filter(variable == "Shannon")
simpson <- riqueza %>% filter(variable == "Simpson")

#Kruskal Wallis test

#Test if there are statistical significant differences between the alpha indexes
#when the samples are grouped by a given variable
kt_indata <- function(indata, ipvars) { #Data frame with richness measures, indata, and vector with the variables, ipvars
  klist <- vector("list", length(ipvars)) #list to save the results of the for cycle
  names(klist) <- ipvars #to preserve the identity of each variable tested in the list
  for(i in 1:length(ipvars)){
    groupi <- indata[,ipvars[i]] #index that indicate to what factor of the given variable each sample belongs
    kti <- kruskal.test(indata$value, groupi) #kruskal wallis test by each variable
    klist[[i]] <- as.vector(kti) #saving the test for each variable to the list
  } 
  return(klist)
}

#Kruskal-Wallis for each measure by each variable
chaokw <- kt_indata(chao1, vrbls)
shankw <- kt_indata(shannon, vrbls)
simpkw <- kt_indata(simpson, vrbls)

#BETA DIVERSITY
#Bray-Curtis distances
bradist <- phyloseq::distance(normbact, method = "bray")
#Ordination method Principal Coordinates Analysis
braord <- phyloseq::ordinate(normbact, method = "PCoA", distance = "bray") #Principal Coordinates Analysis as ordination method
#Beta Diversity Plot
braplo <- plot_ordination(normbact, braord,
                          color = "Host", 
                          shape = "Method",
                          title = "Bray-Curtis") + 
  geom_point(size = 4) +
  theme_bw() #Distinguish by color the insect specie and by shape the method
#Function that does a PERMANOVA on distance matrix grouped on different selected variables
pmav_all <- function(distances, samcols){ #distances is distance object created by phyloseq::distance function
  pmavlist <- vector("list", length(samcols)) #empty list to save the PERMANOVA results
  names(pmavlist) <- samcols #sample_data variables
  for(i in 1:length(samcols)){
    rhsi <- c(as.vector(unlist(metadata[,i]))) #create independent variable vector from the columns of the sample_data table
    set.seed(100) #establish a seed so the p-values are the same in each repetition
    pmavi <- adonis(distances ~ rhsi, permutations = 1000)
    pmavlist[[i]] <- pmavi}
  return(pmavlist)
}
#PERMANOVA result for all the independent variables on Bray-Curtis distances matrix
bpmav <- pmav_all(bradist, vrbls)

##GENERA WITH DIFFERENTIAL COMPOSITION

#Phyloseq object at genera level
generaps <- aggregate_taxa(normbact, "genus", verbose=FALSE)

#Using ANCOM BC
#Required libraries
library("ANCOMBC")
library("dplyr")
library("tidyverse")

#The following function sets the parameters used in this work
#It also automates the replacement of the phyloseq object used, taxrank, which is the phyloseq object grouped to a taxa rank of preference
#and of the group variable, ivar, which is could be any variable of interest present in sample data of the phyloseq object 
anco <- function(taxrank, ivar){
  ancombc(phyloseq = taxrank, formula = ivar, 
          p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000000,
          group = ivar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
          max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)
}
#A warning will be displayed for the variables with less than 3 categories, this is not a problem, it only means that the global analysis will not be done, but in those cases it is not needed

#Adjusted log fold change and standard error tables
#These tables will be used to make the plot

#The following function produces a table with the comparisons between one of the categories with the rest of them in the chosen group variable
#Be warned, if your group variable has more than two categories and you want to compare between all of them all the steps are need to be repeated
lf_sd_tabs <- function(ancomres){
  #Part 1 retrieve and adjust the log fold change between the pairwise comparisons
  dfig_pt1 <- data.frame(ancomres$res$beta * ancomres$res$diff_abn, check.names = FALSE) %>%
    rownames_to_column("taxon_id")
  #Part 2 retrieve and adjust the standard error between the pairwise comparisons
  dfig_pt2 <- data.frame(ancomres$res$se * ancomres$res$diff_abn, check.names = FALSE) %>%
    rownames_to_column("taxon_id")
  colnames(dfig_pt2)[-1] <- paste0(colnames(dfig_pt2)[-1], "SD") #Names to distinguish each value
  #Paste together part 1 and 2, filter the comparison you wish to plot
  dfig_c <- dfig_pt1 %>%
    left_join(dfig_pt2, by = "taxon_id") 
  dfig_c$taxon_id <- factor(dfig_c$taxon_id, levels = dfig_c$taxon_id)
  return(dfig_c)
}


#Before doing the log fold change plot the produced table must be filtered and ordered
#If there is not a single taxa with a log fold change different from zero do not plot

#Function to filter and order the log fold change and standard error table, lfsdtab
#to make it work you need a object of vector class with the column names, compar
filtandord <- function(lfsdtab, compar) {
  lfsdtab %>%
    filter(.data[[compar[[2]]]] != 0) %>% #here we keep only the taxa with a log fold change different to zero, for that we use the second element of compar, that refers to coefficient column
    arrange(desc(abs(.data[[compar[[2]]]]))) #arrange the rows in descending order by their absolute value
}

#Log fold change plot
lfc_plot <- function(plotab, compar){
  ggplot(data = head(plotab, n=10L), #Data is adjusted to plot only the top twenty taxa changes
         aes(x = Taxon, y = .data[[compar[[2]]]])) +
    geom_bar(stat = "identity", width = 0.75, position = position_dodge(width = 0.25)) +
    geom_errorbar(aes(ymin = .data[[compar[[2]]]] - SD, 
                      ymax = .data[[compar[[2]]]] + SD), 
                  width = 0.2, position = position_dodge(0.05), color = "black") +
    geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
    labs(x=NULL, y="Log Fold Change") + 
    coord_flip()
}

#ANCOM-BC results

#By genera
#Genera comparisons
igencomp <- c("Taxon", "EumaeusVRhopalotria", "SD")
igencres <- anco(generaps, "Genera")
igentab <- lf_sd_tabs(igencres)
colnames(igentab) <- igencomp
igentab <- filtandord(igentab, igencomp)
igenplot <- lfc_plot(igentab, igencomp)

#By source
sourcomp <- c("Taxon", "CocultureVSDirect", "SD")
sourcres <- anco(generaps, "Source")
sourtab <- lf_sd_tabs(sourcres)
colnames(sourtab) <- sourcomp
sourtab <- filtandord(sourtab, sourcomp)
sourplot <- lfc_plot(sourtab, sourcomp)
