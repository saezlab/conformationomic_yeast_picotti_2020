library(readr)
library(dplyr)
library(igraph)

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


igraph_net <- graph_from_edgelist(as.matrix(omnipath[,c(1,3)]))
row.names(batches) <- c(1:length(batches[,1]))

batches <- batches[batches$Gene_name %in% omnipath$source_genesymbol | batches$Gene_name %in% omnipath$target_genesymbol,]

nodes_from_shortest_paths <- list()
max_path_length <- 4
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
            if(length(path_backward) <= max_path_length)
            {
              nodes_from_shortest_paths[[k]] <- path_backward
              k <- k+1
              nodes_from_shortest_paths[[k]] <- path_forward
              k <- k+1
            }
          } else
          {
            path_to_add <- ifelse(length(path_backward) > length(path_forward), path_backward, path_forward)
            if(length(path_to_add) <= max_path_length)
            {
              nodes_from_shortest_paths[[k]] <- ifelse(length(path_backward) > length(path_forward), path_backward, path_forward)
              k <- k+1
            }
          }
        } else
        {
          if(length(path_backward) <= max_path_length)
          {
            nodes_from_shortest_paths[[k]] <- path_backward
            k <- k+1
          }
        }
      } else
      {
        if(length(path_forward) != 0)
        {
          if(length(path_forward) <= max_path_length)
          {
            nodes_from_shortest_paths[[k]] <- path_forward
            k <- k+1
          }
        }
      }
    }
  }
}

nodes <- unique(unlist(nodes_from_shortest_paths))

subnet <- omnipath[omnipath$source_genesymbol %in% nodes & omnipath$target_genesymbol %in% nodes,]

#########################
#########################

library(piano)
library(GSEABase)

###SIMPLE MISCALENOUS FUNCTION TO IMPORT GMT FILES INTO AN APPROPRIATE FORMAT FOR PIANO
gmt_to_df <- function(gmtfile, fast = T)
{
  if(fast)
  {
    # genesets = GSEABase::getGmt(con = paste(getwd(),gmtfile, sep = "/"))
    genesets = getGmt(con = gmtfile)
    genesets = unlist(genesets)
    
    gene_to_term =plyr::ldply(genesets,function(geneset){
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
      
    },.progress = plyr::progress_text())
    names(gene_to_term) <- c("gene","term")
    return(gene_to_term[complete.cases(gene_to_term),])
  }
  else
  {
    genesets = getGmt(con = gmtfile)
    genesets = unlist(genesets)
    
    gene_to_term <- data.frame(NA,NA)
    names(gene_to_term) <- c("gene","term")
    for (geneset in genesets)
    {
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
      names(temp3) <- c("gene","term")
      gene_to_term <- rbind(gene_to_term,temp3)
    }
    
    return(gene_to_term[complete.cases(gene_to_term),])
  }
}
### For the hypergeometric test, the success are the nodes of the carnival network output
sucesses <- nodes                   

### Import the pathway collection from GMT file ( downloaded on MsigDB)
pathways <- gmt_to_df("~/Documents/network_tools/data/c2.cp.v7.0.symbols.gmt")
pathways <- pathways[grep("KEGG",pathways$term),]

bg <- unique(c(omnipath$source_genesymbol,omnipath$target_genesymbol))
### Run the hypergeomatric test
kegg_pathways <- runGSAhyper(sucesses, universe = bg, gsc = loadGSC(pathways))
kegg_pathways_df <- as.data.frame(kegg_pathways$resTab)
kegg_pathways_df$pathway <- row.names(kegg_pathways_df)

### Import the pathway collection from GMT file ( downloaded on MsigDB)
GO_BP <- gmt_to_df("~/Dropbox/conformationomic_yeast_picotti_2020/supports/c5.bp.v7.0.symbols.gmt")
GO_BP <- GO_BP[GO_BP$gene %in% bg | GO_BP$gene %in% sucesses,]

### Run the hypergeomatric test
GOBP_pathways <- runGSAhyper(sucesses, universe = bg, gsc = loadGSC(GO_BP))
GOBP_pathways_df <- as.data.frame(GOBP_pathways$resTab)
GOBP_pathways_df$pathway <- row.names(GOBP_pathways_df)

write_csv(kegg_pathways_df[,c(7,2,3)],"~/Dropbox/conformationomic_yeast_picotti_2020/results/KEGG_omnipath.csv")
write_csv(GOBP_pathways_df[,c(7,2,3)],"~/Dropbox/conformationomic_yeast_picotti_2020/results/GOBP_omnipath.csv")
