
setwd("/home/mfcg/Descargas/thesis-seq/bgcs/specific")

library(igraph)

#Terpene network plot
teadma <- read.delim("terp_adj_mat.tsv.csv", header = TRUE, row.names = 1)
teadma <- as.matrix(teadma)
teadma <- 100*(teadma)
gfamte <- graph_from_adjacency_matrix(teadma, 
                                      weighted = TRUE,
                                      mode="undirected", 
                                      diag = FALSE
                                      )
plot.igraph(gfamte,
     vertex.color = c(rep("#B8BCDB", 2), rep("#B8DBA3", 4), "#DBA38C"),
     vertex.frame.color = "white",
     vertex.shape = c(rep(c("square", "rectangle"),2), rep("circle",3)),
     vertex.size = 40,
     vertex.size2 = 35,
     vertex.label = c(rep("Ochrobactrum", 4), rep("Phyllobacterium", 3)),
     vertex.label.color = "black",
     edge.width= 4
     )

#Siderophore network plot
siadma <- read.delim("nrp_adj_mat.tsv.csv", header = TRUE, row.names = 1)
siadma <- as.matrix(siadma)
siadma <- 100*siadma
gfamsi <- graph_from_adjacency_matrix(siadma,
                                      weighted = TRUE,
                                      mode="undirected", 
                                      diag = FALSE
                                      )
plot.igraph(gfamsi,
            vertex.color = c(rep("#B8DBA3", 2), "#DBA38C"),
            vertex.frame.color = "white",
            vertex.shape = c(rep("circle",3)),
            vertex.size = 40,
            vertex.label = c(rep("Phyllobacterium", 3)),
            vertex.label.color = "black",
            edge.width= 4
            )
#Legend
plot(x = 10, y = 10, type = "n")
legend(x = 6, y = 12, legend = c("Rhopalotria adult", "Eumaeus larva", "Eumaeus adult"), fill = c("#B8BCDB","#B8DBA3", "#DBA38C"))
legend(x = 10, y = 12, legend = c("Direct", "Co-culture", "Co-culture with extract"), pch = c(16, 15, 15), pt.cex = 2)
