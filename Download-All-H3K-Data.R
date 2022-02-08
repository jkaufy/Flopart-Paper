library(ggplot2)
library(data.table)
library(dplyr)
library(readr)

datasets <- read.csv("data/H3K-Dataset-List.csv")
data_prefix <- "https://rcdata.nau.edu/genomic-ml/chip-seq-chunk-db/"
count_prefix <- "H3K/Counts"
label_prefix <- "H3K/Labels"
set.seed(1)

for (row in 1:nrow(datasets)) {
  
  count_file <- file.path(count_prefix, datasets$count_destfile[row])
  
  if(!file.exists(count_file)){
    con <- url(paste(data_prefix, datasets$file[row], sep=""))
    load(con)
    names(counts)[names(counts) == 'coverage'] <- 'count'
    dir.create(dirname(count_file), showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(counts, count_file)
    close(con)
  }
  
  label_file <- file.path(label_prefix, datasets$label_destfile[row])
  
  if(!file.exists(label_file)){
    con <- url(paste(data_prefix, datasets$label_file[row], sep=""))
    load(con)
    split.regions <- split(regions, regions$sample.id)
    n.folds <- 2
    
    for (sample.id in names(split.regions)){
      k.fold.vec <- sample(rep(1:n.folds, l=nrow(split.regions[[sample.id]])))
      split.regions[[sample.id]]$random.fold <- k.fold.vec
      split.regions[[sample.id]]$position.fold <- sort(k.fold.vec)
    }
    
    new.regions.dt <- do.call(rbind, split.regions)
    dir.create(dirname(label_file), showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(new.regions.dt, label_file)
    close(con)
  }
}

