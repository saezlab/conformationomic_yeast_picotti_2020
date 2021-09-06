library(seqinr)
library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viper)

source("~/Dropbox/conformationomic_yeast_picotti_2020/scripts/support_functions.R")

CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

conformationomic <- as.data.frame(
  read_delim("~/Dropbox/conformationomic_yeast_picotti_2020/data/Table1.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

phosphoproteomic <-as.data.frame(read_excel("~/Dropbox/conformationomic_yeast_picotti_2020/data/phospho_osmo.xlsx"))
phosphoproteomic <- phosphoproteomic[phosphoproteomic$Uniprot_ID %in% conformationomic$Uniprot_ID,]
phosphoproteomic$`P-site position` <- as.numeric(phosphoproteomic$`P-site position`)

yeast_prot_seq <- read.fasta("~/Dropbox/conformationomic_yeast_picotti_2020/supports/yeast_full_fasta.fasta")

names(yeast_prot_seq) <- gsub("sp[|]","",names(yeast_prot_seq))
names(yeast_prot_seq) <- gsub("[|].*","",names(yeast_prot_seq))

phosphoproteomic$`P-site position` <- apply(phosphoproteomic, 1, function(x, yeast_prot_seq)
  {
  aa_sequence <-  yeast_prot_seq[[x[1]]]
  psite <- paste0(aa_sequence[[as.numeric(x[2])]],x[2])
  return(psite)
},yeast_prot_seq = yeast_prot_seq)

names(phosphoproteomic)[2] <- "pSite"
phosphoproteomic$pSite <- toupper(gsub(" ","",phosphoproteomic$pSite))


for(i in 1:length(yeast_prot_seq))
{
  yeast_prot_seq[[i]] <- paste(c(unlist(yeast_prot_seq[[i]])), collapse = "")
}

yeast_prot_seq <- as.data.frame(yeast_prot_seq)
yeast_prot_seq <- as.data.frame(t(yeast_prot_seq))
yeast_prot_seq$Uniprot_ID <- row.names(yeast_prot_seq)
names(yeast_prot_seq)[1] <- "full_sequence"
yeast_prot_seq$full_sequence <- as.character(yeast_prot_seq$full_sequence)

conformationomic <- merge(conformationomic,yeast_prot_seq)
conformationomic$meanPos <- apply(conformationomic, 1, function(x) {
  return(mean(str_locate(x[13],tolower(x[5]))))
})

names(phosphoproteomic)[c(3,4)] <- paste0("psite_",names(phosphoproteomic)[c(3,4)])

conformationomic <- merge(conformationomic,phosphoproteomic, by = "Uniprot_ID")
# conformationomic$pSite <- ifelse(gsub("",,conformationomic$pSite)

# conformationomic <- conformationomic[!grepl("[/]",conformationomic$pSite),]

conformationomic$psite_position <- as.numeric(gsub("[A-Z]","",conformationomic$pSite))
conformationomic$delta_psite_lip <- abs(conformationomic$meanPos - conformationomic$psite_position)

plot(density(conformationomic$psite_Log2FC))

conformationomic$psite_UID <- paste(conformationomic$Gene_name, conformationomic$pSite, sep = "_")
conformationomic$lip_psite_pair <- paste(conformationomic$Gene_name, conformationomic$Peptide_sequence, conformationomic$pSite, sep = "_")


conformationomic_top_phosphoFC <- conformationomic
conformationomic_top_phosphoFC <- conformationomic_top_phosphoFC %>% group_by(psite_UID) %>% slice(which.min(delta_psite_lip))
top_phosphoFC_small_delta <- conformationomic_top_phosphoFC[conformationomic_top_phosphoFC$delta_psite_lip < 20,]
top_phosphoFC_big_delta <- conformationomic_top_phosphoFC[conformationomic_top_phosphoFC$delta_psite_lip >= 20,]

mean(top_phosphoFC_small_delta$`Qvalue(P.Abundance)`)
mean(top_phosphoFC_big_delta$`Qvalue(P.Abundance)`)
dim(conformationomic_top_phosphoFC[conformationomic_top_phosphoFC$`Qvalue(P.Abundance)` <= 0.2,1])[1] / dim(conformationomic_top_phosphoFC)[1]

to_plot <- top_phosphoFC_small_delta[,c(21,19,16,9,7)]
to_plot <- to_plot[order(to_plot$`Log2FC(LiP_norm)`),]
to_plot$lip_psite_pair <- factor(to_plot$lip_psite_pair, levels = to_plot$lip_psite_pair)
to_plot <- to_plot[to_plot$`Qvalue(LiP)` <= 0.05,]
to_plot_1 <- to_plot[,c(1,2)]
to_plot_2 <- to_plot[,c(1,3)]
to_plot_3 <- to_plot[,c(1,4)]
to_plot_4 <- to_plot[,c(1,5)]


to_plot_1$variable <- names(to_plot_1)[2]
to_plot_1[,2] <- round(to_plot_1[,2], digits = 3)
to_plot_2$variable <- names(to_plot_2)[2]
to_plot_2[,2] <- round(to_plot_2[,2], digits = 3)
to_plot_3$variable <- names(to_plot_3)[2]
names(to_plot_3)[2] <- "Qvalue_lip"
to_plot_3[,2] <- round(to_plot_3[,2], digits = 3)
to_plot_4$variable <- names(to_plot_4)[2]
names(to_plot_4)[2] <- "log2FC_lip"
to_plot_4[,2] <- round(to_plot_4[,2], digits = 3)

a <- ggplot(to_plot_1, aes(x = variable, y = lip_psite_pair, fill = delta_psite_lip)) + geom_tile() +
  scale_fill_gradient2(low="grey", high="grey", mid = "grey", midpoint = 0) + 
  theme_minimal() + theme(legend.position="none", axis.title.x=element_blank(),) + 
  geom_text(aes(label=delta_psite_lip))

b <- ggplot(to_plot_2, aes(x = variable, y = lip_psite_pair, fill = psite_Log2FC)) + geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint = 0) + 
  theme_minimal() + theme(axis.title.y=element_blank(),
                          axis.title.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(), legend.position="none") + 
  geom_text(aes(label=psite_Log2FC))

c <- ggplot(to_plot_3, aes(x = variable, y = lip_psite_pair, fill = Qvalue_lip)) + geom_tile() +
  scale_fill_gradient2(low="red", high="white", mid = "red", midpoint = 0) +
  theme_minimal() + theme(axis.title.y=element_blank(),
                          axis.title.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(), legend.position="none") + 
  geom_text(aes(label=Qvalue_lip))

d <- ggplot(to_plot_4, aes(x = variable, y = lip_psite_pair, fill = log2FC_lip)) + geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid = "white", midpoint = 0) +
  theme_minimal() + theme(axis.title.y=element_blank(),
                          axis.title.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(), legend.position="none") + 
  geom_text(aes(label=log2FC_lip))

grid.arrange(a,b,c,d, nrow = 1)

###########

phospho_osmo_full <- as.data.frame(
  read_excel("~/Dropbox/conformationomic_yeast_picotti_2020/data/phospho_osmo_full.xlsx"))

new_df <- list()
k <- 1
for(i in 1:length(phospho_osmo_full[,1]))
{
  cur_row <- phospho_osmo_full[i,]
  sites <- unlist(strsplit(as.character(cur_row[4]),"[|]"))
  if(length(sites) > 1)
  {
    for(site in sites)
    {
      new_row <- cur_row
      new_row[4] <- site
      new_df[[k]] <- new_row
      k <- k+1
    }
  } else
  {
    new_df[[k]] <- cur_row
    k <- k+1
  }
}

new_df <- as.data.frame(do.call(rbind,new_df))
new_df <- new_df[grepl("Phospho",new_df$Ptm),]
new_df$Ptm <- gsub("[[]","",new_df$Ptm)
new_df$Ptm <- gsub("[]].*","",new_df$Ptm)
new_df$residue <- apply(new_df,1,function(x)
  {
  return(substr(x[3], x[4], x[4]))
}) 
new_df$relative_psite <- paste(new_df$Uniprot_ID, new_df$residue, sep = "_")
new_df$relative_psite <- paste0(new_df$relative_psite, new_df$Ptm)


new_df$position <- apply(new_df, 1, function(x, yeast_prot_seq)
  {
  start_pos <- str_locate(yeast_prot_seq[x[1],1],tolower(x[3]))[1]
  # print(str_locate(yeast_prot_seq[x[1],1],tolower(x[3]))[1])
  pos <- start_pos + as.numeric(x[4]) - 1
  return(pos)
},yeast_prot_seq = yeast_prot_seq)

new_df$psite <- paste0(new_df$residue, new_df$position)
new_df$psite_ID <- paste(new_df$Uniprot_ID, new_df$psite, sep = "_")

new_df <- new_df %>% group_by(psite_ID) %>% slice(which.max(Log2FC))

phospho_uniprot_to_symbole <- as.data.frame(read_delim("~/Dropbox/conformationomic_yeast_picotti_2020/supports/phospho_uniprot_to_symbole", 
                                                       "\t", escape_double = FALSE, trim_ws = TRUE))

new_df <- merge(new_df, phospho_uniprot_to_symbole, by = "Uniprot_ID")
new_df$pSite <- paste(new_df$gene_symbol, new_df$psite, sep = "_")

phospho_full <-new_df[,c(14,5,6,7)]
phospho_full$z_score <- apply(phospho_full,1,function(x)
  {
  pval <- as.numeric(x[3])
  sign_fc <- sign(as.numeric(x[2]))
  if(sign_fc > 0)
  {
    pval <- 1 - pval
    z_score <- qnorm(pval)
  }
  else
  {
    z_score <- qnorm(pval)
  }
  return(z_score)
})
###########

biogrid_yeast_KSN <- as.data.frame(
  read_csv("~/Dropbox/conformationomic_yeast_picotti_2020/supports/biogrid_yeast_KSN.csv"))

KSN_viper <- unique(biogrid_yeast_KSN[,c(5,6,7)])
KSN_viper <- KSN_viper[KSN_viper$action != "-",]

KSN_viper$action <- ifelse(KSN_viper$action == "kinase",1,-1)

KSN_viper <- df_to_viper_regulon(KSN_viper)

phospho_for_viper <- phospho_full[,5, drop = F]
row.names(phospho_for_viper) <- phospho_full$pSite
# phospho_for_viper <- phospho_for_viper[!grepl("[/]",phospho_for_viper$pSite),]
# row.names(phospho_for_viper) <- paste(phospho_for_viper$`Gene Name`, phospho_for_viper$pSite, sep = "_")

sum(row.names(phospho_for_viper) %in% biogrid_yeast_KSN$psite) # SMALL overlap

# phospho_for_viper <- phospho_for_viper[,c(4),drop = F]

kinact <- as.data.frame(viper(eset = phospho_for_viper, regulon = KSN_viper, nes = T, minsize = 3, eset.filter = F, cores = 3))
kinact <- kinact[abs(kinact$z_score) > 1.7,,drop = F]
##############

kinase_network <- biogrid_yeast_KSN[biogrid_yeast_KSN$kinase %in% row.names(kinact),]
kinase_network <- kinase_network[kinase_network$psite %in% row.names(phospho_for_viper),]
kinase_network <- unique(kinase_network[,c(6,7,5)])

sub_part <- biogrid_yeast_KSN[biogrid_yeast_KSN$kinase %in% row.names(kinact),]
sub_part <- sub_part[sub_part$psite %in% row.names(phospho_for_viper),]
sub_part <- unique(sub_part[,c(5,2)])
sub_part$action <- 1

kinase_network$action <- ifelse(kinase_network$action == "kinase",1,-1)
names(kinase_network) <- c("source","sign","target")

names(sub_part) <- c("source","target","sign")

kinase_network <- as.data.frame(rbind(kinase_network,sub_part))

nodes <- as.data.frame(unique(c(kinase_network$source, kinase_network$target)))
names(nodes) <- "id"
nodes$type <- ifelse(grepl("_",nodes$id),"psite","prot")
nodes$type <- ifelse(nodes$id %in% biogrid_yeast_KSN$kinase, "kinase",nodes$type)
nodes$value <- NA

for(node in nodes$id)
{
  if(node %in% row.names(phospho_for_viper))
  {
    nodes[nodes$id == node,"value"] <- phospho_for_viper[node,]
  }
}
for(node in nodes$id)
{
  if(node %in% row.names(kinact))
  {
    nodes[nodes$id == node,"value"] <- kinact[node,]
  }
}

Table1 <- as.data.frame(
  read_delim("~/Dropbox/conformationomic_yeast_picotti_2020/data/Table1.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

changing_prots <- Table1[Table1$`Qvalue(LiP)` <= 0.1 & abs(Table1$`Log2FC(LiP_norm)`) > 1,]
changing_prots <- changing_prots[,c(2,7)]

batches <- changing_prots %>% group_by(Gene_name) %>% summarise_each(funs(min(abs(.))))
batches <- as.data.frame(batches)
row.names(batches) <- batches$Gene_name


for(node in nodes$id)
{
  if(is.na(nodes[nodes$id == node,3]))
  {
    if(node %in% row.names(batches))
    {
      nodes[nodes$id == node,"value"] <- batches[node,2]
    }
  }
}

nodes$label <- nodes$id
nodes$color <- ifelse(nodes$value > 0, "red","blue")
nodes$shape <- sapply(nodes$type, function(x)
{
  switch(x, "prot" = "circle", "psite" = "star", "kinase" = "diamond")
})

edges <- kinase_network
names(edges) <- c("from","sign","to")
edges$arrows <- "to"


to_filter <- edges
to_filter$source_value <- NA
to_filter$target_value <- NA
for(i in 1:length(to_filter[,1]))
{
  to_filter[i,"source_value"] <- nodes[nodes$id == to_filter[i,1],3]
  to_filter[i,"target_value"] <- nodes[nodes$id == to_filter[i,3],3]
}
to_filter$to_keep <- TRUE
for(i in 1:length(to_filter[,1]))
{
  if(grepl("_",to_filter[i,3]))
  {
    if(sign(to_filter[i,"source_value"]) == sign(to_filter[i,"sign"]) * sign(to_filter[i,"target_value"]))
    {
      to_filter[i,"to_keep"] <- TRUE
    } else
    {
      to_filter[i,"to_keep"] <- FALSE
    }
  }
}

edges <- edges[to_filter$to_keep,]
edges <- edges[!grepl("_",edges$from) | edges$from %in% edges$to,]
edges$color <- ifelse(edges$sign > 0, "red", "blue")
edges$width <- 4

# nodes$value <- abs(nodes$value)

nodes <- nodes[(nodes$id %in% edges$from | nodes$id %in% edges$to),] #need tofilter lonely psites too
nodes <- nodes[complete.cases(nodes),]
edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id,]
edges <- edges[!(grepl("_",edges$to)) | edges$to %in% edges$from,]

nodes <- nodes[nodes$id %in% edges$from | nodes$id %in% edges$to,]

nodes$type <- apply(nodes,1,function(x)
{
  new_type <- switch(x[2],
    "kinase" = ifelse(sum(edges[edges$from == x[1],2]) < 0, "phosphatase", "kinase"),
    "psite" = "psite",
    "prot" = "prot")
  
  return(new_type)
})


##

##

library(visNetwork)


nodes <- nodes[!grepl("_",nodes$id) | abs(nodes$value) > 1.5,]
edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id,]
nodes <- nodes[nodes$id %in% edges$from | nodes$id %in% edges$to,]
visNetwork(nodes = nodes, edges = edges)

nodes$kinase_Act <- ifelse(nodes$type == "kinase" | nodes$type == "phosphatase" & nodes$type, , NA)
nodes$lipFC <- ifelse(nodes$type == "prot", nodes$value, NA)
nodes$phospho_t_val <- ifelse(nodes$type == "psite", nodes$value, NA)

nodes$value_adjusted <- nodes$value
nodes$value_adjusted <- ifelse(nodes$type == "prot", nodes$value + 100, nodes$value_adjusted)
nodes$value_adjusted <- ifelse(nodes$type == "psite", nodes$value - 100, nodes$value_adjusted)



nodes$label <- gsub(".*_","",nodes$label) 
nodes$label_nocap <- tolower(nodes$label)
nodes$label_nocap <- sapply(nodes$label_nocap, CapStr)



write_csv(nodes,"~/Dropbox/conformationomic_yeast_picotti_2020/results/network_kinase_phospho_lip_att.csv")
write_csv(edges,"~/Dropbox/conformationomic_yeast_picotti_2020/results/network_kinase_phospho_lip_sif.csv")
write_csv(nodes[,c(1,7)],"~/Dropbox/conformationomic_yeast_picotti_2020/results/network_kinase_phospho_lip_att_nocaps.csv")

###############

edges_reduced <- edges
nodes_reduced <- nodes

kick_out_kinases <- c("PKP1","OCA1","PPH22","CBK1","TOR1","TPK1")

edges_reduced <- edges_reduced[!(edges_reduced$from %in% kick_out_kinases),]
edges_reduced <- edges_reduced[!grepl("_",edges_reduced$from) | edges_reduced$from %in% edges_reduced$to,]
edges_reduced$color <- ifelse(edges_reduced$sign > 0, "red", "blue")
edges_reduced$width <- 4
nodes_reduced$value <- abs(nodes_reduced$value)

nodes_reduced <- nodes_reduced[(nodes_reduced$id %in% edges_reduced$from | nodes_reduced$id %in% edges_reduced$to),] #need tofilter lonely psites too

nodes_reduced <- nodes_reduced[complete.cases(nodes_reduced),]
edges_reduced <- edges_reduced[edges_reduced$from %in% nodes_reduced$id & edges_reduced$to %in% nodes_reduced$id,]
edges_reduced <- edges_reduced[!(grepl("_",edges_reduced$to)) | edges_reduced$to %in% edges_reduced$from,]

nodes_reduced <- nodes_reduced[nodes_reduced$id %in% edges_reduced$from | nodes_reduced$id %in% edges_reduced$to,]
##

##

library(visNetwork)

visNetwork(nodes = nodes_reduced, edges = edges_reduced)
