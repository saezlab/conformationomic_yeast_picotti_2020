library(readr)
library(dplyr)
library(CARNIVAL)

Table1 <- as.data.frame(
  read_delim("Dropbox/conformationomic_yeast_picotti_2020/data/Table1.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

changing_prots <- Table1[Table1$`Qvalue(LiP)` <= 0.05,]
changing_prots <- changing_prots[,c(1,8)]

batches <- changing_prots %>% group_by(Uniprot_ID) %>% summarise_each(funs(min(abs(.))))
batches <- as.data.frame(batches)
batches$Uniprot_ID <- gsub(";.*","",batches$Uniprot_ID) 

carni_causal_network_yeast <- as.data.frame(
  read_delim("Dropbox/conformationomic_yeast_picotti_2020/supports/carni_causal_network_yeast.tsv", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE))

batches$Uniprot_ID <- gsub("[-+{},;() ]","______",batches$Uniprot_ID)

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

batches <- batches[batches$Uniprot_ID %in% carni_causal_network_yeast$source | batches$Uniprot_ID %in% carni_causal_network_yeast$target,]

batches$zscore <- abs(qnorm(batches$`Pvalue(LiP)`))

carnival_input <- as.data.frame(t(batches[,c(3)]))
names(carnival_input) <- gsub("[-+{},;() ]","______",batches[,1])

carni_causal_network_yeast$source <- gsub("[-+{},;() ]","______",carni_causal_network_yeast$source)
carni_causal_network_yeast$target <- gsub("[-+{},;() ]","______",carni_causal_network_yeast$target)

carni_causal_network_yeast <- unique(carni_causal_network_yeast)
carni_causal_network_yeast <- carni_causal_network_yeast[carni_causal_network_yeast$source != carni_causal_network_yeast$target,]

write_tsv(carnival_input,"~/Dropbox/conformationomic_yeast_picotti_2020/results/changing_prots.tsv")
write_tsv(carni_causal_network_yeast, "~/Dropbox/conformationomic_yeast_picotti_2020/supports/carni_causal_network.tsv")

setwd("~/Dropbox/conformationomic_yeast_picotti_2020/results/carnival/")

CARNIVAL_Result <- runCARNIVAL(CplexPath="~/Documents/cplex",
                               Result_dir=".",
                               CARNIVAL_example=NULL,
                               UP2GS=F,
                               netFile = "../../supports/carni_causal_network.tsv", 
                               measFile = "../changing_prots.tsv", 
                               inverseCR = T,
                               weightFile = NULL, 
                               timelimit = 28800,
                               mipGAP = 0.15
) #11.24
