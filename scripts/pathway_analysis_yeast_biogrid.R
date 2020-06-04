library(readr)
library(dplyr)
library(igraph)
library(bionetdata)
library(ReactomePA)
library(org.Sc.sgd.db)
library(piano)

x <- org.Sc.sgdENTREZID
mapped_genes <- mappedkeys(x)
geneID_to_ORFIF <- as.list(x[mapped_genes])
geneID_to_ORFIF <- as.data.frame(geneID_to_ORFIF)
geneID_to_ORFIF <- as.data.frame(t(geneID_to_ORFIF))
geneID_to_ORFIF$V1 <- as.character(geneID_to_ORFIF$V1)

vec_geneID_to_ORFIF <- row.names(geneID_to_ORFIF)
names(vec_geneID_to_ORFIF) <- geneID_to_ORFIF$V1

pathways <- ReactomePA:::get_Reactome_DATA("yeast")$EXTID2PATHID
pathways_name_mapping <- ReactomePA:::get_Reactome_DATA("yeast")$PATHID2NAME

pathways_df <- list()
i <- 1
for(gene in names(pathways))
{
  sub_path_df <- as.data.frame(matrix(NA,length(pathways[[gene]]),2))
  sub_path_df[,2] <- gene
  sub_path_df[,1] <- unlist(pathways[[gene]])
  pathways_df[[i]] <- sub_path_df
  # j <- 1
  # for(path_name in pathways[[gene]])
  # {
  #   sub_path_df
  #   j <- j + 1
  # }
  i <- i + 1
}

pathways_df <- as.data.frame(do.call(rbind,pathways_df))
for(i in 1:length(pathways_df[,1]))
{
  pathways_df[i,1] <- pathways_name_mapping[pathways_df[i,1]]
  pathways_df[i,2] <- vec_geneID_to_ORFIF[pathways_df[i,2]]
}

pathways_df <- pathways_df[,c(2,1)]
names(pathways_df) <- c("gene","term")

data(Yeast.Biogrid.data)

Yeast.Biogrid.data <- as.data.frame(Yeast.Biogrid.data)

Table1 <- as.data.frame(
  read_delim("Dropbox/conformationomic_yeast_picotti_2020/data/Table1.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

changing_prots <- Table1[Table1$`Qvalue(LiP)` <= 0.05,]
changing_prots <- changing_prots[,c(1,8)]

batches <- changing_prots %>% group_by(Uniprot_ID) %>% summarise_each(funs(min(abs(.))))
batches <- as.data.frame(batches)
batches$Uniprot_ID <- gsub(";.*","",batches$Uniprot_ID) 


mapping <- as.data.frame(
  read_delim("Dropbox/conformationomic_yeast_picotti_2020/supports/uniprot-saccharomyces+cerevisiae-filtered-organism__Saccharomyces+cerevisiae--.tab", 
             "\t", escape_double = FALSE, trim_ws = TRUE))
mapping$`Gene names  (ordered locus )` <- gsub("[;].*","",mapping$`Gene names  (ordered locus )`)
mapping_vec <- mapping$`Gene names  (ordered locus )`
names(mapping_vec) <- mapping$Entry

for(i in 1:length(batches[,1]))
{
  batches[i,1] <- mapping_vec[batches[i,1]]
}
batches <- batches[complete.cases(batches),]

batches <- batches[batches$Uniprot_ID %in% names(Yeast.Biogrid.data),]

igraph_net <- graph_from_adjacency_matrix(as.matrix(Yeast.Biogrid.data))
row.names(batches) <- c(1:length(batches[,1]))


nodes_from_shortest_paths <- list()
k <- 1
for(i in 1:length(batches[,1]))
{
  print(i)
  for(j in i+1:length(batches[,1]))
  {
    if(j != i & j <= length(batches[,1]))
    {
      path_forward <- names(unlist(shortest_paths(igraph_net, from = batches[i,1], to = batches[j,1])[[1]]))
      path_backward <- names(unlist(shortest_paths(igraph_net, from = batches[j,1], to = batches[i,1])[[1]]))
      if(length(path_backward) != 0)
      {
        if(length(path_forward) != 0)
        {
          if(length(path_backward) == length(path_forward))
          {
            nodes_from_shortest_paths[[k]] <- path_backward
            k <- k+1
            
          } else
          {
            nodes_from_shortest_paths[[k]] <- ifelse(length(path_backward) > length(path_forward), path_backward, path_forward)
            k - k+1
          }
        } else
        {
          nodes_from_shortest_paths[[k]] <- path_backward
          k <- k+1
        }
      } else
      {
        if(length(path_forward) != 0)
        {
          nodes_from_shortest_paths[[k]] <- path_forward
          k <- k+1
        }
      }
    }
  }
}

nodes <- unique(unlist(nodes_from_shortest_paths))

biogrid_network <- as.data.frame(as_edgelist(igraph_net))

gsc  <- loadGSC(pathways_df)

ORA_res <- runGSAhyper(genes = nodes, universe = names(Yeast.Biogrid.data), gsc = gsc, adjMethod = "fdr")
ORA_res <- as.data.frame(ORA_res$resTab)
ORA_res$pathway <- row.names(ORA_res)

stres_response_genes <- pathways_df[pathways_df$term == "Cellular responses to stress",1]
stress_response <- biogrid_network[biogrid_network$V1 %in% stres_response_genes & biogrid_network$V2 %in% stres_response_genes,]

cell_cycle_genes <- pathways_df[pathways_df$term == "Cell Cycle",1]
cell_cycle <- biogrid_network[biogrid_network$V1 %in% cell_cycle_genes & biogrid_network$V2 %in% cell_cycle_genes,]

#Network of interaction is so dense...
save.image("~/Dropbox/conformationomic_yeast_picotti_2020/results/shortest_paths/pathway_analysis.RData")
