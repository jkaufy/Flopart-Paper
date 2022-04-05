# need to run Download-All-H3K-Data.R first
# Load in all data sets
library(data.table)
source("Load-All-H3k-Data.R")

H3K_data <- loadH3KData()

# increase above 10 log scale
penalties <- 10^seq(-5, 5, by=0.5)
n.fold <- 2

sets <- list("train","test")
algoritms <- list(Flopart = FLOPART::FLOPART, 
    FPOP =  PeakSegOptimal::PeakSegFPOPchrom)

cache.prefix <- "figure-label-errors-data"
feature.list <- list()

for(dataset in 1:length(H3K_data$count)){
  print(dataset)
  
  one_count <- H3K_data$count[[dataset]]
  one_label <- H3K_data$labels[[dataset]]
  
  sample_split_count <- split(one_count, one_count$sample.id)
  sample_split_label <- split(one_label, one_label$sample.id)
  
  if (!('peaks' %in% one_label$annotation)){
    for (sample.id in names(sample_split_count)){
      sample.err.list <- list()
      
      cache.save <- paste(dataset, "-", sample.id, ".csv", sep = "")
      cache.file <- file.path(cache.prefix, cache.save)
      
      one_sample_count <- sample_split_count[[sample.id]]
      
      bases <- one_sample_count[nrow(one_sample_count), ]$chromEnd - 
        one_sample_count[1,]$chromStart
      
      feature.list[[paste(dataset, sample.id)]] <- data.table(
        dataset,
        sample.id,
        size = bases,
        log.log.data=log(log(bases))
      )
      
      if(!file.exists(cache.file)){
        one_sample_label <- sample_split_label[[sample.id]]
        
        # for each penalty
        for(pen in penalties){
          # for each fold
          for(fold in 1:n.fold){
            segs.list <- list()
            
            setDT(one_sample_label)
            one_sample_label[, set := ifelse(fold==random.fold, "test", "train")]
            
            for (algo in names(algoritms)){
              algo.fun <- algoritms[[algo]]
              
              if(algo == "Flopart"){
                fit <- algo.fun(one_sample_count, 
                                one_sample_label[set == "train"], pen)
                segs <- fit[["segments_dt"]]
              }else{
                fit <- algo.fun(one_sample_count, pen)
                segs <- data.table(fit$segments)
              }
              
              segs.list[[paste(algo)]] <- data.table(
                algo,
                segs)
            }
            
            for (set.i in sets){
              for (model in names(segs.list)){
                pkg.segs <- segs.list[[model]][, .(chromStart, chromEnd, mean, status)]
                pkg.peaks <- pkg.segs[status=="peak"]
                err.dt <- PeakError::PeakErrorChrom(pkg.peaks, 
                                                    one_sample_label[set == set.i])
                
                sample.err.list[[paste(dataset, sample.id, pen, fold, model, set.i)]] <- data.table(
                  dataset,
                  sample.id,
                  pen,
                  fold,
                  model,
                  set.i,
                  possible.fp = sum(err.dt$possible.fp),
                  fp = sum(err.dt$fp),
                  possible.fn = sum(err.dt$possible.tp),
                  fn = sum(err.dt$fn),
                  labels = nrow(one_sample_label[set == set.i])
                )
              }
              
            }
          }
        }
        
        sample.err.dt <- do.call(rbind, sample.err.list)
        dir.create(dirname(cache.file), showWarnings = FALSE, recursive = TRUE)
        data.table::fwrite(sample.err.dt, cache.file)
        gc()
      }
    }
  }
}

feature.dt <- do.call(rbind, feature.list)
data.table::fwrite(feature.dt, "feature-dt.csv")

