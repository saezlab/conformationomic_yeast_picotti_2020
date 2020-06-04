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

changing_prots <- Table1[Table1$`Qvalue(LiP)` <= 0.05 & abs(Table1$`Log2FC(LiP_norm)`) > 1,]
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

starts <- c("YIL147C","YER118C")


nodes_from_shortest_paths <- list()
k <- 1
for(start in starts)
{
  print(start)
  for(j in 1:length(batches[,1]))
  {
    path_forward <- names(unlist(shortest_paths(igraph_net, from = start, to = batches[j,1])[[1]]))
    if(length(path_forward) != 0)
    {
      nodes_from_shortest_paths[[k]] <- path_forward
      k <- k+1
    }
  }
}

nodes <- unique(unlist(nodes_from_shortest_paths))

biogrid_network <- as.data.frame(as_edgelist(igraph_net))
biogrid_network$sign <- 1
biogrid_network <- biogrid_network[,c(1,3,2)]
biogrid_network$V1 <- gsub("[-+{},;() ]","______",biogrid_network$V1)
biogrid_network$V2 <- gsub("[-+{},;() ]","______",biogrid_network$V2)
write_tsv(biogrid_network,"~/Dropbox/conformationomic_yeast_picotti_2020/supports/biogrid_yeast_sif.tsv")

carni_meas <- batches
carni_meas <- carni_meas[order(carni_meas$`Pvalue(LiP)`, decreasing = T),]
carni_meas <- as.data.frame(t(carni_meas[,2]))
names(carni_meas) <- batches[,1]
carni_meas <- carni_meas[,names(carni_meas) %in% biogrid_network$V1 | names(carni_meas) %in% biogrid_network$V2,]
names(carni_meas) <-gsub("[-+{},;() ]","______",names(carni_meas))
write_tsv(carni_meas,"~/Dropbox/conformationomic_yeast_picotti_2020/results/carni_meas_ordered_locus_name.tsv")

carni_input <- as.data.frame(matrix(NA,1,2))
names(carni_input) <- c("YIL147C","YER118C")
carni_input[1,] <- c(1,1)
write_tsv(carni_input,"~/Dropbox/conformationomic_yeast_picotti_2020/results/carni_input_ordered_locus_name.tsv")


gsc  <- loadGSC(pathways_df)

full_network_shortest_path <-  biogrid_network[biogrid_network$V1 %in% nodes & biogrid_network$V2 %in% nodes,]

ORA_res <- runGSAhyper(genes = nodes, universe = names(Yeast.Biogrid.data), gsc = gsc, adjMethod = "fdr")
ORA_res <- as.data.frame(ORA_res$resTab)
ORA_res$pathway <- row.names(ORA_res)

write_csv(ORA_res[,c(7,2,3)],"~/Dropbox/conformationomic_yeast_picotti_2020/results/reactome_biogrid_HOG.csv")

ORA_res_top <- ORA_res[order(ORA_res$`Adjusted p-value`, decreasing = F),]
ORA_res_top <- ORA_res_top[1:40,]

top_genes <- batches

gene_mapping <- mapping$`Gene names`
gene_mapping <- gsub(" .*","",gene_mapping)
names(gene_mapping) <- mapping$`Gene names  (ordered locus )`

for(i in 1:length(top_genes[,1]))
{
  if(top_genes[i,1] %in% names(gene_mapping))
  {
    top_genes[i,1] <- gene_mapping[top_genes[i,1]]
  }
}

pathways_df_genename <- pathways_df
for(i in 1:length(pathways_df_genename[,1]))
{
  if(pathways_df_genename[i,1] %in% names(gene_mapping))
  {
    pathways_df_genename[i,1] <- gene_mapping[pathways_df_genename[i,1]]
  }
}

plot_list <- list()
i <- 1
for(pathway in ORA_res_top$pathway)
{
  plot_list[[i]] <- plot_top_genes(pathway,top_genes,pathways_df_genename) #diagonal gene names and pathway name on plot
  i <- i+1
}
names(plot_list) <- ORA_res_top$pathway
names(plot_list) <- gsub("[/ ]","_",names(plot_list))

i <- 1
for(pathway_plot in plot_list)
{
  ggsave(paste0("~/Dropbox/conformationomic_yeast_picotti_2020/visualisation/osmotic_",names(plot_list)[i],".pdf"), plot = pathway_plot, device = "pdf",width = 10, height = 10)
  i <- i+1
}
plot_list[[2]]

full_network_shortest_path_genenames <- full_network_shortest_path

for(i in 1:length(full_network_shortest_path_genenames[,1]))
{
  if(full_network_shortest_path_genenames[i,1] %in% names(gene_mapping))
  {
    full_network_shortest_path_genenames[i,1] <- gene_mapping[full_network_shortest_path_genenames[i,1]]
  }
  if(full_network_shortest_path_genenames[i,3] %in% names(gene_mapping))
  {
    full_network_shortest_path_genenames[i,3] <- gene_mapping[full_network_shortest_path_genenames[i,3]]
  }
}



pathways_df_genename <- pathways_df_genename[complete.cases(pathways_df_genename),]

pathway_name <- "Glycolysis"

# glycolysis_pathway <- full_network_shortest_path_genenames[
#   (full_network_shortest_path_genenames$V1 %in% pathways_df_genename[pathways_df_genename$term == pathway_name,1]) &
#     full_network_shortest_path_genenames$V2 %in% pathways_df_genename[pathways_df_genename$term == pathway_name,1],
# ]

library(visNetwork)
names(top_genes)[1] <- "id"

nodes <- top_genes[top_genes$id %in% pathways_df_genename[pathways_df_genename$term == pathway_name,1],]
names(nodes) <- c("id","size")
nodes$locusname <- nodes$id

gene_mapping_reverse <- mapping$`Gene names  (ordered locus )`
names(gene_mapping_reverse) <- gsub(" .*","",mapping$`Gene names`)

for(i in 1:length(nodes[,1]))
{
  if(nodes[i,"locusname"] %in% names(gene_mapping_reverse))
  {
    nodes[i,"locusname"] <- gene_mapping_reverse[nodes[i,"locusname"]]
  }
}

starts <- c("YIL147C","YER118C")


nodes_from_shortest_paths <- list()
k <- 1
for(start in starts)
{
  print(start)
  for(j in 1:length(nodes[,3]))
  {
    path_forward <- names(unlist(shortest_paths(igraph_net, from = start, to = nodes[j,3])[[1]]))
    if(length(path_forward) != 0)
    {
      nodes_from_shortest_paths[[k]] <- path_forward
      k <- k+1
    }
  }
}
sp_nodes <- unique(unlist(nodes_from_shortest_paths))
for(i in 1:length(sp_nodes))
{
  if(sp_nodes[i] %in% names(gene_mapping))
  {
    sp_nodes[i] <- gene_mapping[sp_nodes[i]]
  }
}

edges <- full_network_shortest_path_genenames[
  full_network_shortest_path_genenames$V1 %in% sp_nodes &
    full_network_shortest_path_genenames$V2 %in% sp_nodes,
]

edges <- edges[,c(1,3,2)]
names(edges) <- c("from","to","action")
edges <- edges[edges$from != edges$to,]
edges$edge_id <- paste(edges$from, edges$to,sep = "_")
edges$edge_id_reverse <- paste(edges$to, edges$from,sep = "_")

for(i in length(edges[,1]):1)
{
  if(edges[i,"edge_id_reverse"] %in% edges$edge_id)
  {
    edges <- edges[-i,]
  }
}

edges <- edges[,c(-4,-5)]
####need to create edge list with these nodes and add them to the node attributes

nodes$pathway <- pathway_name

other_nodes <- unique(c(edges$from, edges$to))
other_nodes <- other_nodes[!other_nodes %in% nodes$id]
other_nodes <- as.data.frame(other_nodes)
names(other_nodes) <- "id"
other_nodes <- merge(other_nodes,top_genes, all.x = T)
other_nodes$pathway <- "other"
other_nodes$locusname <- NA
names(other_nodes)[2] <- "size"
nodes <- as.data.frame(rbind(nodes,other_nodes))

starts_name <- c("SLN1","SHO1")

nodes$pathway <- ifelse(nodes$id %in% starts_name, "starter", nodes$pathway)
nodes$pathway <- ifelse(nodes$id %in% pathways_df_genename[pathways_df_genename$term == pathway_name,"gene"],
                        pathway_name, nodes$pathway)

nodes$color <- "grey"
nodes[nodes$pathway == pathway_name,"color"] <- "lightblue"
nodes[nodes$pathway == "starter","color"] <- "red"

nodes$size <- -log10(nodes$size)*10

nodes$label <- nodes$id

# edges <- glycolysis_pathway[,c(1,3,2)]
# names(edges) <- c("from","to","action")

visNetwork(nodes, edges) 

write_csv(nodes,"~/Dropbox/conformationomic_yeast_picotti_2020/results/shortest_paths/OSMO_glycolysis_att.csv")
write_csv(edges,"~/Dropbox/conformationomic_yeast_picotti_2020/results/shortest_paths/OSMO_glycolysis_sif.csv")


stres_response_genes <- pathways_df[pathways_df$term == "Cellular responses to stress",1]
stress_response <- biogrid_network[biogrid_network$V1 %in% stres_response_genes & biogrid_network$V2 %in% stres_response_genes,]

cell_cycle_genes <- pathways_df[pathways_df$term == "Cell Cycle",1]
cell_cycle <- biogrid_network[biogrid_network$V1 %in% cell_cycle_genes & biogrid_network$V2 %in% cell_cycle_genes,]

#Network of interaction is so dense...
save.image("~/Dropbox/conformationomic_yeast_picotti_2020/results/shortest_paths/pathway_analysis.RData")
