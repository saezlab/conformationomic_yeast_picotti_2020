library(ggplot2)

plot_top_genes <- function(pathway, top_genes, pathways_df, pval = T)
{
  pathway_genes_conf <- top_genes[top_genes[,1] %in% pathways_df[pathways_df[,2] == pathway,1],]
  
  names(pathway_genes_conf) <- c("ID","stat")
  
  if(pval)
  {
    pathway_genes_conf$stat <- -log10(pathway_genes_conf$stat)
  }
  
  pathway_genes_conf <- pathway_genes_conf[order(pathway_genes_conf$stat),]
  pathway_genes_conf$ID <- factor(pathway_genes_conf$ID , levels = unique(pathway_genes_conf$ID))
  
  gp <- ggplot(pathway_genes_conf, aes(x = ID, y = stat)) + 
    geom_bar(stat= "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45)) +
    ggtitle(pathway)
  
  return(gp)
}

df_to_viper_regulon <- function(df)
{
  names(df) <- c("feature","pathway","sign")
  df <- df[complete.cases(df),]
  
  pathway_regulon <- list(0)
  i <- 1
  for(pathway in unique(df$pathway))
  {
    pathway_feature_list <- list(0)
    features <- df[df$pathway == pathway, 3]
    names(features) <- df[df$pathway == pathway, 1]
    pathway_feature_list[[1]] <- features
    pathway_feature_list[[2]] <- rep(1,length(features))
    names(pathway_feature_list) <- c("tfmode","likelihood")
    
    pathway_regulon[[i]] <- pathway_feature_list
    i <- i+1
  }
  names(pathway_regulon) <- unique(df$pathway)
  return(pathway_regulon)
}
