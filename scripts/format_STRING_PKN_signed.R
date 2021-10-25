library(readr)
library(org.Sc.sgd.db)

X4932_protein_links_detailed_v11_0 <- as.data.frame(
  read_table2("~/Dropbox/conformationomic_yeast_picotti_2020/supports/4932.protein.links.detailed.v11.0.txt")) #need to download from STITCHdb (large file)

X4932_protein_links_detailed_v11_0 <- X4932_protein_links_detailed_v11_0[X4932_protein_links_detailed_v11_0$combined_score >= "70022",]

yeast_PPI_noTM <- X4932_protein_links_detailed_v11_0

yeast_PPI_noTM <- yeast_PPI_noTM[yeast_PPI_noTM$database >= 700,]
# yeast_PPI_noTM <- yeast_PPI_noTM[unlist(apply(yeast_PPI_noTM[,c(3:8)], 1, function(x) {max(x) >= 900})),]
yeast_PPI_noTM$edge_id <- paste(yeast_PPI_noTM$protein1,yeast_PPI_noTM$protein2, sep = "_")


X4932_protein_actions_v11_0 <- as.data.frame(
  read_table2("~/Dropbox/conformationomic_yeast_picotti_2020/supports/4932.protein.actions.v11.0.txt")) #need to download from STITCHdb (large file)

yeast_PKN <- X4932_protein_actions_v11_0[X4932_protein_actions_v11_0$is_directional,]
yeast_PKN$edge_id <- paste(yeast_PKN$item_id_a, yeast_PKN$item_id_b, sep = "_")

yeast_PKN <- yeast_PKN[yeast_PKN$edge_id %in% yeast_PPI_noTM$edge_id,]

yeast_PKN$item_id_a <- gsub(".*[.]","",yeast_PKN$item_id_a)
yeast_PKN$item_id_b <- gsub(".*[.]","",yeast_PKN$item_id_b)

names(yeast_PKN)[c(1,2)] <- c("source","target")
yeast_PKN_signed <- yeast_PKN

yeast_PKN$sign <- 1

yeast_PKN$source <- gsub("[-+{},;() ]","______",yeast_PKN$source)
yeast_PKN$target <- gsub("[-+{},;() ]","______",yeast_PKN$target)
yeast_PKN <- unique(yeast_PKN)

yeast_PKN_signed$sign <- ifelse(yeast_PKN_signed$action == "activation", 1, ifelse(yeast_PKN_signed$action == "inhibition", -1, 0))
yeast_PKN_signed <- yeast_PKN_signed[yeast_PKN_signed$sign != 0,]
yeast_PKN_signed <- yeast_PKN_signed[,c(1,2,8)]
names(yeast_PKN_signed) <- c("source","target","sign")

prots <- unique(c(yeast_PKN_signed$source, yeast_PKN_signed$target))

mapping_ensemblprot_to_symbole <- mapIds(org.Sc.sgd.db, prots, 'GENENAME','ENSEMBLPROT')
yeast_PKN_signed$source <- mapping_ensemblprot_to_symbole[yeast_PKN_signed$source]
yeast_PKN_signed$target <- mapping_ensemblprot_to_symbole[yeast_PKN_signed$target]
yeast_PKN_signed <- unique(yeast_PKN_signed)

save(yeast_PKN_signed, file ="data/STRING_yeast_PKN_700.RData")
