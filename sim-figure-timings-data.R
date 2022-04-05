library(data.table)
library(microbenchmark)

peakDownCap <- 15
peakDownMin <- 1
peakUpCap <- 50
labelSizeCap <- 6

create_signal <- function (peak_min, peak_cap, label_size_cap){
  signal <- sample(peak_min:peak_cap, label_size_cap)
  return(signal)
}

N <- 500000
between_labels <- 10
label_start <- 0 
which_label <- 0
label_end <- 0
total_signal <- c()
label_vec <- c(0)

while( label_start < N){
  label_size <- sample(3:labelSizeCap, 1)
  label_end <- label_start + (between_labels*label_size)
  label_vec <- c(label_vec, label_end)
  
  if (which_label %% 2 == 0){
    signal <- create_signal(peakDownCap, peakUpCap, label_size)
  }
  else{
    signal <- create_signal(peakDownMin, peakDownCap, label_size)
  }
  
  which_label <- which_label + 1
  label_start <- label_end
  total_signal <- c(total_signal, signal)
}

label.vec.size <- length(label_vec)

all.labels <- data.table(
  chromStart=as.integer(label_vec[1:label.vec.size-1]),
  chromEnd=as.integer(label_vec[2:label.vec.size]),
  annotation = c("peakStart","peakEnd"))

break.vec <- seq(10, label_start, by=10)
all.counts <- data.table(
  chromStart=as.integer(break.vec-between_labels),
  chromEnd=as.integer(break.vec),
  count = as.integer(total_signal))

some.N <- 50000
some.labels <- all.labels[chromEnd <= some.N]
some.counts <- all.counts[chromEnd <= some.N]
some.counts[, weight.vec := chromEnd-chromStart]
size.vec <- unique(ceiling(10^seq(0, log10(nrow(some.labels)-1), l=10)))
size.vec
timing_list <- list()

for(size.i in seq_along(size.vec)){
  size <- size.vec[[size.i]]
  labels <- some.labels[1:size,]
  Lfit <- FLOPART::FLOPART(some.counts, labels, 5)
  timing <- microbenchmark(
    FLOPART={
      FLOPART::FLOPART(some.counts, labels, 5)
    },
    FLOPARTNoLabel={
      FLOPART::FLOPART(some.counts, penalty =  5)
    },
    FPOP={
      PeakSegOptimal::PeakSegFPOPchrom(some.counts, 5)
    },
    binsegRcpp={
      binsegRcpp::binseg('poisson', data.vec = some.counts$count, 
        weight.vec = some.counts$weight.vec, max.segments = nrow(Lfit$segments_dt))
    },
    times=3)
  timing_list[[paste(size)]] <- data.table(size, timing)
}
timing.dt <- do.call(rbind, timing_list)
data.table::fwrite(timing.dt, "figure-timings-data-labels.csv")

## experiment 2: variable number of data, constant number of labels
size.vec <- unique(ceiling(10^seq(0, log10(nrow(all.labels)-1), l=10)))
timing_list <- list()
for(size.i in seq_along(size.vec)){
  size <- size.vec[[size.i]]
  cat(sprintf("%4d / %4d labels=%d\n", size.i, length(size.vec), size))
  labels <- all.labels[1:size]
  n.data <- max(labels$chromEnd)
  some.data <- all.counts[chromEnd <= n.data]
  some.data[, weight.vec := chromEnd-chromStart]
  for(data.over.labels in c(10, 1000)){
    n.labels <- n.data / data.over.labels
    some.labels <- labels[seq(1, .N, l=n.labels)]
    some.labels <- unique(some.labels)
    Lfit <- FLOPART::FLOPART(some.data, some.labels, 5)
    m.args <- list(
      FLOPART=FLOPART::FLOPART(some.data, some.labels, 5),
      FLOPARTNoLabel=FLOPART::FLOPART(some.data, penalty = 5),
      FPOP=PeakSegOptimal::PeakSegFPOPchrom(some.data, 5),
      binsegRcpp=
        binsegRcpp::binseg('poisson', data.vec = some.data$count, 
          weight.vec = some.data$weight.vec)
      ,
      times=3)

    timing <- do.call(microbenchmark, m.args)
    timing_list[[paste(data.over.labels, size)]] <- data.table(
        data.over.labels, size, timing, n.data, n.labels,
        label.density=n.labels/n.data)
    
  }
}
timing.dt <- do.call(rbind, timing_list)
data.table::fwrite(timing.dt, "figure-timings-data.csv")
