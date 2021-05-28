####INSECT ECOLOGICAL DATA

#Path and read files
#Path
setwd('/home/mfcg/Documentos/Cycad_insects/Ecology-collect-data')
#Table with processes and each individual record
fulldata <- read.csv("All_insects.csv", header = TRUE)
#Table with accumulate count of each insect sample 
counts <- read.csv("eco_insect_count.csv", header = TRUE) 

#Required libraries
library(dplyr)
library(ggplot2)
library(palettesForR)
library(vegan)
library(RColorBrewer)
library(networkD3)

##Working tables with ecological data
#Columns with ecological data
ecocolumns <- c("ind", "individualID", "species", "stage", "plant_part", 
                "location", "sampletype", "year")
#Ecological data of each individual
ecodata <- distinct(fulldata[ ,ecocolumns]) 
#Species found in each location only in the random stratified samplings
ransam <- ecodata %>% 
  filter(sampletype == "random_strat") %>% #Just random stratified
  select(c(location, species)) #Select only the location and species columns
ransam <- as.matrix(table(ransam)) #Make a matrix for diversity analysis

##Insect counts
#Total of collected insects
sum(counts$n)
#Total of collected species
length(unique(counts$TAXON))
#Information about the sampling 
sampling <- list(counts$YEAR, counts$SAMTYPE, counts$LOCATION) 
#Number of collected insects by samplings
tapply(X = counts$n, INDEX = sampling, FUN = sum)
#Insect species by samplings
samplist <- split(x = counts$TAXON, f = sampling, drop = TRUE)
lapply(X = samplist, FUN = unique)

#Estimate diversity indexes
shannon <- vegan::diversity(ransam, index = "shannon")
simpson <- vegan::diversity(ransam, index = "simpson")
spnum <- specnumber(ransam)
pielou <- shannon/log(spnum)
#Data frame with the diversity indexes
vindex <- as.numeric(c(shannon, simpson, spnum, pielou)) #Index values
nindex <- c(rep("H",2), rep("1-S", 2), rep("#species", 2), rep("J", 2)) #Index name
lindex <- rep(c("JC", "MO"), 5) #Location name
dindex_table <- data.frame(Index = nindex,
                           Location = lindex,
                           Value = vindex) #Data frame

#Location and diversity plots
# These are not in the thesis #
#The distribution of insects species and stages by location
ins_loc <- ggplot(counts, aes(fill=STAGE, y=n, x=TAXON)) + #In the X axis is the species
  geom_bar(position = "stack", stat = "identity") + #Stack by stage
  facet_wrap(~LOCATION) + #Divide by location
  theme_bw() + 
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text = element_text(angle= 90,size = 12),
        strip.text.x = element_text(size = 12)) +
  ylab("") +
  xlab("")
#Diversity index comparison between locations
loc_div <- ggplot(dindex_table, aes(y=Value, x=Location, color=Location)) +
  geom_point(size=4) +
  facet_wrap(~Index) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        strip.text.x = element_text(size = 12)) +
  ylab("") +
  xlab("LOCATION")



##Plant part preference plots

##Sankey Network
#Data frame with the links between the insect species, "TAXON", the plant part 
#where they were found, "PLANT_PART", and the number of times, "n"
links <- counts[ , c("TAXON", "PLANT_PART", "n")]
#Data frame for the unique names of "TAXON" and "PLANT_PART" to be used as nodes
nodes <- data.frame(
  name=c(as.character(links$TAXON),
         as.character(links$PLANT_PART)) %>% unique()
)
#Associate the data frame with the links and nodes
links$IDtaxon <- match(links$TAXON, nodes$name)-1
links$IDplpt <- match(links$PLANT_PART, nodes$name)-1
#Vector with the insect development stages to group the links
links$group <- as.factor(c(counts$STAGE))
#Vector that groups the nodes as the insect order or if it is a plant
nodes$group <- as.factor(c("Coleoptera", "Coleoptera", "Lepidoptera", 
                           "Hymenoptera", "Coleoptera", "Coleoptera",
                           "Coleoptera", "Coleoptera", "Coleoptera",
                           "Trichoptera", "Cycad", "Cycad", "Cycad"))
#Associate the link a node groups with a color
colorcito <- 'd3.scaleOrdinal() 
.domain(["adult", "pupa", "larva", "egg", "Lepidoptera", "Coleoptera", "Hymenoptera", "Trichoptera", "Plant"]) 
.range(["#E88D6D", "#F28172", "#DB727D", "#F272BF", "#BD331E", "#123570", "#123570", "#285EBD", "#557012"])' #Colors
#Plot
tvpsn <- sankeyNetwork(Links = links, #Insect-plant part association
                       Nodes = nodes, #Insect taxa and plant part
                       Source = "IDtaxon", #Insect taxa
                       Target = "IDplpt", #Plant part
                       Value = "n", #Each insect individual found
                       NodeID = "name", #Unique values for insect taxa and plant part
                       colourScale = colorcito, #Color association with the groups
                       LinkGroup = "group", #Development stages
                       NodeGroup = "group", #Insect order or cycad 
                       fontSize = 12)
#Legend for the colors and groups
plot(x = c(1:10), y = c(1:10), type = "n") #Blank 10x10 plot
legend(x=3, y=7, legend = c("adult", "pupa", "larva", "egg"), #Developmental stages
       fill = c("#E88D6D", "#F28172", "#DB727D", "#F272BF"),  #Link colors
       title = "Links")
legend(x=7, y=7, legend = c("Lepidoptera", "Coleoptera", "Hymenoptera", "Trichoptera", "Cycad"), #Insect order or cycad
       fill = c("#BD331E", "#123570", "#123570", "#285EBD", "#557012"), #Node colors
       title = "Nodes")

# The following plot was not included in the thesis #
# Replaced by the Sankey plot #

##Network
#but can be visualized with the following code
#Required library
library(igraph)
#Edges table, which insect taxa "TAXON" interacts with each structure "PLANT_PART"
tvp_count <- counts[ , c("TAXON", "PLANT_PART")]
#Attributes of the vertex, "TAXON" and "PLANT_PART"
tvp_atr <- as.data.frame(tibble(
  nodes= c(unique(tvp_count$TAXON), unique(counts$PLANT_PART)), #avoid repetitions
  atr = c(rep("insect", 10), rep("plant", 3))) #to indicate insect or cycad part
)
#The number of individuals are indicated by the edge width 
tvp_ewidth <- c(counts$n)
#The edge color will representate the insect stage
tvp_scol <- gsub("adult", "black", counts$STAGE) 
tvp_scol <- gsub("pupa", "red", tvp_scol)
tvp_scol <- gsub("larva", "pink", tvp_scol)
tvp_scol <- gsub("egg", "orange", tvp_scol)
#igraph graph created from the edge list, tvp_count, and the vertex atributes, tvp_atr
tvp_net <- graph_from_data_frame(d = tvp_count, vertices = tvp_atr, directed = FALSE)
#The vertex color indicates if it corresponds to an insect or a cycad part
scolor <- brewer.pal(4, "Set2") #Color vector
tvp_color <- scolor[as.numeric(as.factor(V(tvp_net)$atr))] #Linking color and factor insect/plant
#Insect and D. edule network
plot(tvp_net,
     vertex.color = tvp_color, #Vertex is insect or plant indicated by color
     vertex.size = 30,
     edge.color = tvp_scol, # individual stage in each interaction
     edge.width = tvp_ewidth, # number of individuals by interaction
     layout = layout.kamada.kawai #network distribution
     )


  
