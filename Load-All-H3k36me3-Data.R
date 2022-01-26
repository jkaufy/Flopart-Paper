# Loads in all the H3k36me3 datasets
loadH3K36me3Datasets <- function() {
  library(data.table)
  
  counts <- "~/Flopart-Paper/data/H3K36me3/Counts/"
  labels <- "~/Flopart-Paper/data/H3K36me3/Labels/"
  file_prefixs <- list(count = counts,labels = labels)
  file_locations <- list()
  
  for (file.type in names(file_prefixs)){
    file.prefix <- file_prefixs[[file.type]]
    file_names <- list.files(file.prefix, pattern = "*.rds")
    file_locations[[file.type]] <- paste(file.prefix, 
                                         list.files(file.prefix, pattern = "*.rds"), sep="")
  }
  
  H3K36me3_datasets <- list()
  for (file.type in names(file_locations)){
    H3K36me3_datasets[[file.type]] <- lapply(file_locations[[file.type]], 
                                             function(x){readRDS(file = x)})
  }
  
  return(H3K36me3_datasets)
}