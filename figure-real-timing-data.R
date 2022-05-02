source("Load-All-H3k-Data.R")
library(microbenchmark)
library(data.table)

H3K_data <- loadH3KData()
pen <- 0.001

timing_list <- list()
# length(H3K_data$count)
for(dataset in 1:length(H3K_data$count)){
  print(dataset)
  
  one_count <- H3K_data$count[[dataset]]
  one_label <- H3K_data$labels[[dataset]]
  
  sample_split_count <- split(one_count, one_count$sample.id)
  sample_split_label <- split(one_label, one_label$sample.id)
  
  if (!('peaks' %in% one_label$annotation)){
    for (sample.id in names(sample_split_count)){
      sample.err.list <- list()
      one_sample_count <- sample_split_count[[sample.id]]
      one_sample_label <- sample_split_label[[sample.id]]
      
      setDT(one_sample_count)
      setDT(one_sample_label)
      
      one_sample_count[, weight.vec := chromEnd-chromStart]
      Lfit <- FLOPART::FLOPART(one_sample_count, one_sample_label, pen)
      
      lopart_labels <- FLOPART::FLOPART_data(one_sample_count, one_sample_label)$label_dt
      
      one_sample_label[, start := lopart_labels$firstRow]
      one_sample_label[, end := ifelse(lopart_labels$lastRow < nrow(one_sample_count),lopart_labels$lastRow , nrow(one_sample_count))]
      one_sample_label[, changes := ifelse(annotation=="noPeaks", 0, 1)]
      
      
      size <- nrow(one_sample_count)
      
      timing <- microbenchmark(
        FLOPART={
          FLOPART::FLOPART(one_sample_count, one_sample_label, pen)
        },
        FPOP={
          PeakSegOptimal::PeakSegFPOPchrom(one_sample_count, pen)
        },
        binsegRcpp={
          binsegRcpp::binseg('poisson', data.vec = one_sample_count$count, 
              weight.vec = one_sample_count$weight.vec, max.segments = nrow(Lfit$segments_dt))
        },
        LOPART={
          LOPART::POISSON_LOPART(one_sample_count$count, one_sample_count$weight.vec, one_sample_label, pen)
        },
        times=3)
      
      timing_list[[paste(dataset, sample.id, size)]] <- data.table(
        dataset,
        sample.id,
        size, 
        timing
      )
      
    }
  }
}

timing.dt <- do.call(rbind, timing_list)
data.table::fwrite(timing.dt, "figure-real-timings-data.csv")
