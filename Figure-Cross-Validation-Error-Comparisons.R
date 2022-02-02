# Load in all data sets
library(data.table)
source("Load-All-H3k36me3-Data.R")

H3K36me_data <- loadH3K36me3DataSetRData()


# increase above 10 log scale
penalties <- 10^seq(-3, 3, l=20)
penalties <- 1
n.fold <- 2

total_label_error_list <- list()
err_dt_list <- list()

sets <- list("train","test")
algoritms <- list(Flopart = FLOPART::FLOPART, 
    FPOP =  PeakSegOptimal::PeakSegFPOPchrom)

cache.prefix <- "~/R/Flopart-Paper/cache/cross-validation/"

# length(H3K36me_data$count)
for(dataset in 1:1){
  print(dataset)
  
  one_count <- H3K36me_data$count[[dataset]]
  one_label <- H3K36me_data$labels[[dataset]]
  
  sample_split_count <- split(one_count, one_count$sample.id)
  sample_split_label <- split(one_label, one_label$sample.id)
  
  if (!('peaks' %in% one_label$annotation)){
    for (sample.id in names(sample_split_count)){
      one_sample_count <- sample_split_count[[sample.id]]
      one_sample_label <- sample_split_label[[sample.id]]
      
      # for each penalty
      for(pen in penalties){
        # for each fold
        for(fold in 1:n.fold){
          segs.list <- list()
          
          setDT(one_sample_label)
          one_sample_label[, set := ifelse(fold==random.fold, "test", "train")]
          
          cache.save <- paste(dataset, "-", sample.id, "-", pen, "-", fold,".RData", sep = "")
          cache.file <- paste(cache.prefix,cache.save, sep = "")
          
          if(file.exists(cache.file)){
            
            segs.list <- readRDS(cache.file)
            
          }else{
            
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
            saveRDS(segs.list, file = cache.file)
            
          }
          
          for (set in sets){
            for (model in names(segs.list)){
              pkg.segs <- segs.list[[model]][, .(chromStart, chromEnd, mean, status)]
              pkg.peaks <- pkg.segs[status=="peak"]
              err.df <- PeakError::PeakErrorChrom(pkg.peaks, 
                                                  one_sample_label[set == set])
              
              err_dt_list[[paste(dataset, sample.id, pen, fold, model, set)]] <- data.table(
                dataset,
                sample.id,
                pen,
                fold,
                model,
                set,
                err.df)
              
              total_label_error_list[[paste(dataset, sample.id, pen, fold, model, set)]] <- data.table(
                dataset,
                sample.id,
                pen,
                fold,
                model,
                set,
                error = sum(err.df[, 'fp']) + sum(err.df[, 'fn']))
            }
          }
        }
        gc()
      }
    }
  }
}

total.label.error.dt <- do.call(rbind, total_label_error_list)
total.label.error.dt <- total.label.error.dt[order(dataset, sample.id, pen, fold)]

# write.csv(total.label.error.dt, file="total.label.error.dt.csv")


