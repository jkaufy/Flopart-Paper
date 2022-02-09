# Loads in all the H3K data sets
createH3KData <- function() {
  library(data.table)
  
  counts <- file.path("H3K/Counts")
  labels <- file.path("H3K/Labels")
  file_prefixs <- list(count = counts,labels = labels)
  file_locations <- list()
  
  for (file.type in names(file_prefixs)){
    file.prefix <- file_prefixs[[file.type]]
    file_names <- file.path(file.prefix, 
                            list.files(file.prefix, pattern = "*.csv"))
    file_locations[[file.type]] <- file_names
  }
  
  H3K_data <- list()
  
  for (file.type in names(file_locations)){
    H3K_data[[file.type]] <- lapply(file_locations[[file.type]], 
                                    function(x){read.csv(x)})
  }
  saveRDS(H3K_data, file = "H3K/H3K-Data.RData")
}

loadH3KData <- function() {
  H3K.loc <- "H3K/H3K-Data.RData"
  
  if(!file.exists(H3K.loc)){
    createH3KData()
  }
    
  H3K_data <- readRDS(H3K.loc)
  return(H3K_data)
}
