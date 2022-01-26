# Load in all datasets
source("Load-All-H3k36me3-Data.R")

H3K36me_data <- loadH3K36me3Datasets()

count_data <- H3K36me_data$count
label_data <- H3K36me_data$labels

# lets run with one data set first

penalties <- c(0.001, 0.01, 0.1, 1, 1.5, 2, 2.5, 3, 5, 10)
n.fold <- 2

# for each sample ID
algo = "Flopart"
total_label_error_list <- list()

for( dataset in 5:length(count_data)){
  print(dataset)
  
  one_count <- count_data[[dataset]]
  one_label <- label_data[[dataset]]
  
  sample_split_count <- split(one_count, one_count$sample.id)
  sample_split_label <- split(one_label, one_label$sample.id)
  
  dist
  
  for (sample.id in names(sample_split_count)){
    print(sample.id)
    one_sample_count <- sample_split_count[[sample.id]]
    one_sample_label <- sample_split_label[[sample.id]]
    
    # for each penalty
    for(pen in penalties){
      segs.list <- list()
      # for each fold
      algo = "Flopart"
      for(fold in range(1:n.fold)){
        train_set <- one_sample_label[one_sample_label$random.fold == fold,]
        
        # run flopart in this
        flopart <- FLOPART::FLOPART(one_sample_count, one_sample_label, pen)
        FLOPART.segs <- flopart[["segments_dt"]]
        
        segs.list[[paste(dataset, algo, sample.id, pen, fold)]] <- data.table(
          dataset,
          algo,
          sample.id,
          pen,
          fold,
          flopart[["segments_dt"]])
      }
      
      algo = "FPOP"
      fit <- PeakSegOptimal::PeakSegFPOPchrom(one_sample_count, pen)
      fpop.segs <- data.table(fit$segments)
      
      segs.list[[paste(dataset, algo, sample.id, pen, fold)]] <- data.table(
        dataset,
        algo,
        sample.id,
        pen,
        fold,
        data.table(fit$segments))
      
      for (seg in segs.list){
        pkg.peaks <- seg[status=="peak"]
        err.df <- PeakError::PeakErrorChrom(pkg.peaks, one_sample_label)
        errors <- sum(err.df$fp) + sum(err.df$fn)
        
        algo = seg$algo[1]
        fold = seg$fold[1]
        
        total_label_error_list[[paste(dataset, sample.id, pen, fold, algo)]] <- data.table(
          dataset,
          sample.id,
          pen,
          fold,
          algo,
          errors)
      }
      
    }
  }
}

total.label.error.dt <- do.call(rbind, total_label_error_list)
total.label.error.dt <- total.label.error.dt[order(dataset, sample.id, pen)]

