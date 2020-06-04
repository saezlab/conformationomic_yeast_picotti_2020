library(seqinr)
library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(gridExtra)


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

### I need to split psites with / to keep HOG
