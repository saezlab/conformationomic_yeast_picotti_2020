library(seqinr)
library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viper)

source("~/Dropbox/conformationomic_yeast_picotti_2020/scripts/support_functions.R")

conformationomic <- as.data.frame(
  read_delim("Dropbox/conformationomic_yeast_picotti_2020/data/Table1.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

phosphoproteomic <-as.data.frame(read_excel("Dropbox/conformationomic_yeast_picotti_2020/data/Phospho_published.xlsx"))
phosphoproteomic <- phosphoproteomic[phosphoproteomic$`Gene Name` %in% conformationomic$Gene_name,]


yeast_prot_seq <- read.fasta("~/Dropbox/conformationomic_yeast_picotti_2020/supports/yeast_full_fasta.fasta")

names(yeast_prot_seq) <- gsub("sp[|]","",names(yeast_prot_seq))
names(yeast_prot_seq) <- gsub("[|].*","",names(yeast_prot_seq))

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

names(phosphoproteomic)[2] <- "Gene_name"
conformationomic <- merge(conformationomic,phosphoproteomic[,c(2,3,4)], by = "Gene_name")
# conformationomic$pSite <- ifelse(gsub("",,conformationomic$pSite)

conformationomic <- conformationomic[!grepl("[/]",conformationomic$pSite),]

conformationomic$psite_position <- as.numeric(gsub("[A-Z]","",conformationomic$pSite))
conformationomic$delta_psite_lip <- abs(conformationomic$meanPos - conformationomic$psite_position)

plot(density(conformationomic$`(dFC/dT)max`))

conformationomic$psite_UID <- paste(conformationomic$Gene_name, conformationomic$pSite, sep = "_")
conformationomic$lip_psite_pair <- paste(conformationomic$Gene_name, conformationomic$Peptide_sequence, conformationomic$pSite, sep = "_")

conformationomic_top_phosphoFC <- conformationomic[conformationomic$`(dFC/dT)max` > 0.1,]
conformationomic_top_phosphoFC <- conformationomic_top_phosphoFC %>% group_by(psite_UID) %>% slice(which.min(delta_psite_lip))
top_phosphoFC_small_delta <- conformationomic_top_phosphoFC[conformationomic_top_phosphoFC$delta_psite_lip < 20,]
top_phosphoFC_big_delta <- conformationomic_top_phosphoFC[conformationomic_top_phosphoFC$delta_psite_lip >= 20,]

mean(top_phosphoFC_small_delta$`Qvalue(P.Abundance)`)
mean(top_phosphoFC_big_delta$`Qvalue(P.Abundance)`)



to_plot <- top_phosphoFC_small_delta[,c(20,18,16,9,7)]
to_plot$lip_psite_pair <- factor(to_plot$lip_psite_pair, levels = to_plot$lip_psite_pair)
to_plot_1 <- to_plot[,c(1,2)]
to_plot_2 <- to_plot[,c(1,3)]
to_plot_3 <- to_plot[,c(1,4)]
to_plot_4 <- to_plot[,c(1,5)]


to_plot_1$variable <- names(to_plot_1)[2]
to_plot_1[,2] <- round(to_plot_1[,2], digits = 3)
to_plot_2$variable <- names(to_plot_2)[2]
names(to_plot_2)[2] <- "max_delta_fc"
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

b <- ggplot(to_plot_2, aes(x = variable, y = lip_psite_pair, fill = max_delta_fc)) + geom_tile() +
  scale_fill_gradient2(low="white", high="red") + 
  theme_minimal() + theme(axis.title.y=element_blank(),
                          axis.title.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(), legend.position="none") + 
  geom_text(aes(label=max_delta_fc))

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

biogrid_yeast_KSN <- as.data.frame(
  read_csv("Dropbox/conformationomic_yeast_picotti_2020/supports/biogrid_yeast_KSN.csv"))

KSN_viper <- unique(biogrid_yeast_KSN[,c(5,6,7)])
KSN_viper <- KSN_viper[KSN_viper$action != "-",]

KSN_viper$action <- ifelse(KSN_viper$action == "kinase",1,-1)

KSN_viper <- df_to_viper_regulon(KSN_viper)

phospho_for_viper <- as.data.frame(read_excel("Dropbox/conformationomic_yeast_picotti_2020/data/Phospho_published.xlsx"))
phospho_for_viper <- phospho_for_viper[!grepl("[/]",phospho_for_viper$pSite),]
row.names(phospho_for_viper) <- paste(phospho_for_viper$`Gene Name`, phospho_for_viper$pSite, sep = "_")

sum(row.names(phospho_for_viper) %in% biogrid_yeast_KSN$psite) # SMALL overlap

phospho_for_viper <- phospho_for_viper[,c(4),drop = F]

kinact <- viper(eset = phospho_for_viper, regulon = KSN_viper, nes = T, minsize = 3, eset.filter = F, cores = 3)

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
  read_delim("Dropbox/conformationomic_yeast_picotti_2020/data/Table1.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

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

library(visNetwork)

visNetwork(nodes = nodes, edges = edges)

################
library(bio3d)

pdb <- read.pdb("~/Dropbox/conformationomic_yeast_picotti_2020/supports/5nbn.pdb")
pdb$seqres

ca.inds <- atom.select(pdb, "calpha", chain = "D")
ca.inds

##############