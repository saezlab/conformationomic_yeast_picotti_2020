library(readr)
library(dplyr)
library(CARNIVAL)

Table1 <- as.data.frame(
  read_delim("Dropbox/conformationomic_yeast_picotti_2020/data/Table1.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ",", grouping_mark = "."), trim_ws = TRUE))

changing_prots <- Table1[Table1$`Qvalue(LiP)` <= 0.05,]
changing_prots <- changing_prots[,c(2,8)]

batches <- changing_prots %>% group_by(Gene_name) %>% summarise_each(funs(min(abs(.))))
batches <- as.data.frame(batches)

##Import and generate a causal network from OMNIPATH

url <- paste0(
  'http://omnipathdb.org/interactions?',
  'fields=sources,references&genesymbols=1'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

omnipath <- download_omnipath()
omnipath <- omnipath[omnipath$is_directed != 0,]
omnipath <- omnipath[,c(3,5,4)]

batches <- batches[batches$Gene_name %in% omnipath$source_genesymbol | batches$Gene_name %in% omnipath$target_genesymbol,]
batches$zscore <- abs(qnorm(batches$`Pvalue(LiP)`))

carnival_input <- as.data.frame(t(batches[,c(3)]))
names(carnival_input) <- gsub("[-+{},;() ]","______",batches[,1])

omnipath$source_genesymbol <- gsub("[-+{},;() ]","______",omnipath$source_genesymbol)
omnipath$target_genesymbol <- gsub("[-+{},;() ]","______",omnipath$target_genesymbol)

omnipath <- omnipath[omnipath$source_genesymbol != omnipath$target_genesymbol,]
omnipath <- unique(omnipath)

write_tsv(carnival_input,"~/Dropbox/conformationomic_yeast_picotti_2020/results/changing_prots.tsv")
write_tsv(omnipath, "~/Dropbox/conformationomic_yeast_picotti_2020/supports/carni_causal_network.tsv")

setwd("~/Dropbox/conformationomic_yeast_picotti_2020/results/carnival/")

CARNIVAL_Result <- runCARNIVAL(CplexPath="~/Documents/cplex",
                               Result_dir=".",
                               CARNIVAL_example=NULL,
                               UP2GS=F,
                               netFile = "../../supports/carni_causal_network.tsv", 
                               measFile = "../changing_prots.tsv", 
                               inverseCR = T,
                               weightFile = NULL, 
                               timelimit = 10000,
                               mipGAP = 0.1
                            
) #11.24
