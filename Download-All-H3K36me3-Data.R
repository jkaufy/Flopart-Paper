library(ggplot2)
library(data.table)
library(dplyr)
library(readr)

datasets <- read_csv("data/H3K36me3-Dataset-List.csv")
data_prefix <- "https://rcdata.nau.edu/genomic-ml/chip-seq-chunk-db/"
count_prefix <- "~/Flopart-Paper/data/H3K36me3/Counts/"
label_prefix <- "~/Flopart-Paper/data/H3K36me3/Labels/"
set.seed(1)

for (row in 1:nrow(datasets)) {
  
  count_file <- paste(count_prefix, 
        datasets$count_destfile[row], sep="")
  
  if(!file.exists(count_file)){
    con <- url(paste(data_prefix, datasets$file[row], sep=""))
    load(con)
    names(counts)[names(counts) == 'coverage'] <- 'count'
    # write.csv(counts,count_file, row.names = FALSE)
    write.csv(counts, file=count_file)
    close(con)
  }
  
  label_file <- paste(label_prefix, 
        datasets$label_destfile[row], sep="")
  
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
    
    # write.csv(new.regions.dt, label_file, row.names = FALSE)
    write.csv(new.regions.dt, file=label_file)
    close(con)
  }
}


