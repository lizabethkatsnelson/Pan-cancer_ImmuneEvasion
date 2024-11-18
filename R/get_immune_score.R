library(tidyverse)

### function to get cytotoxic immune score from expression table
# df is expression table, column1 = genes, the rest are samples
# genecol = column number where gene names are

get_immune_score <- function(df, genecol="", genes = c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")) {
  
  df_i <- as.data.frame(t(df[df[,genecol] %in% genes, ] %>% remove_rownames %>% column_to_rownames(var=genecol))) # genes of interest
  df_ranks <- as.data.frame(apply(df_i, 2, rank)) # rank each gene per sample
  df_ranks$Ranked_Sum <- rank(rowSums(df_ranks)) # add immune score by taking sum of all ranks, rank the final sum
  
  return(df_ranks)
}