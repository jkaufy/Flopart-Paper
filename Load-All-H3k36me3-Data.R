# Loads in all the H3k36me3 data sets
createH3K36me3DataSetRData <- function() {
  library(data.table)
  
  counts <- "~/R/Flopart-Paper/data/H3K36me3/Counts/"
  labels <- "~/R/Flopart-Paper/data/H3K36me3/Labels/"
  file_prefixs <- list(count = counts,labels = labels)
  file_locations <- list()
  
  for (file.type in names(file_prefixs)){
    file.prefix <- file_prefixs[[file.type]]
    file_names <- paste(file.prefix, 
                        list.files(file.prefix, pattern = "*.csv"), sep = "")
    file_locations[[file.type]] <- file_names
  }
  
  H3K36me3_datasets <- list()
  
  for (file.type in names(file_locations)){
    H3K36me3_datasets[[file.type]] <- lapply(file_locations[[file.type]], 
                                             function(x){read.csv(file = x)})
  }
  saveRDS(H3K36me3_datasets, 
    file = "~/R/Flopart-Paper/data/H3K36me3/H3K36me3-datasets.RData")
}

loadH3K36me3DataSetRData <- function() {
  H3K.loc <- "~/R/Flopart-Paper/data/H3K36me3/H3K36me3-datasets.RData"
  
  if(!file.exists(H3K.loc)){
    createH3K36me3DataSetRData()
  }
    
  H3K36me3_datasets <- readRDS(H3K.loc)
  return(H3K36me3_datasets)
}