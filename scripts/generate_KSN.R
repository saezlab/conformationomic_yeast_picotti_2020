library(readr)

setwd("~/Dropbox/conformationomic_yeast_picotti_2020/supports/")
# Both dataset downloaded as zip from :
# https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-3.5.186/BIOGRID-PTMS-3.5.186.ptm.zip
BIOGRID_PTM_3_5_186_ptmtab_true_relations <- as.data.frame(
  read_delim("BIOGRID-PTMS-3.5.186.ptm/BIOGRID-PTM-3.5.186.ptmtab.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE))

BIOGRID_PTM_3_5_186_ptmtab_true_relations <- BIOGRID_PTM_3_5_186_ptmtab_true_relations[
BIOGRID_PTM_3_5_186_ptmtab_true_relations$`Has Relationships` == "TRUE",
]


BIOGRID_PTM_RELATIONSHIPS_3_5_186_ptmrel <- as.data.frame(
  read_delim("BIOGRID-PTMS-3.5.186.ptm/BIOGRID-PTM-RELATIONSHIPS-3.5.186.ptmrel.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE))


kinases <- BIOGRID_PTM_RELATIONSHIPS_3_5_186_ptmrel[,c(1,5,7,8,10)]
names(kinases)[c(1,2,3,4)] <- c("psite_ID","kinase","action","type")

psites <- BIOGRID_PTM_3_5_186_ptmtab_true_relations[,c(1,5,11,9)]
psites$psite <- paste(psites$`Official Symbol`,psites$Residue, sep = "_")
psites$psite <- paste0(psites$psite, psites$Position)
names(psites)[1] <- "psite_ID"

psites <- merge(psites,kinases, by = "psite_ID")

write_csv(psites,"biogrid_yeast_KSN.csv")
