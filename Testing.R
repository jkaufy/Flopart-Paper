library(ggplot2)
library(data.table)
library(dplyr)
source("Load-All-H3k36me3-Data.R")


H3K36me_data <- loadH3K36me3DataSetRData()


counts <- H3K36me_data$count[[1]]
regions <- H3K36me_data$labels[[1]]

split.counts <- split(counts, counts$sample.id)
split.regions <- split(regions, regions$sample.id)

test.count <- split.counts[["McGill0001"]]
test.regions <- split.regions[["McGill0001"]]

fold.i <- 1


flopart <- FLOPART::FLOPART(test.count, test.regions, 1400)

FLOPART.segs <- flopart[["segments_dt"]]


seg.dt.list <- list(
  FLOPART=FLOPART.segs)

err.dt.list <- list()
for(model in names(seg.dt.list)){
  pkg.segs <- seg.dt.list[[model]][, .(chromStart, chromEnd, mean, status)]
  pkg.peaks <- pkg.segs[status=="peak"]
  err.df <- PeakError::PeakErrorChrom(pkg.peaks, test.regions)
  err.dt.list[[model]] <- data.table(model, err.df)
}
err.dt <- do.call(rbind, err.dt.list)

FLOPART.peaks <- FLOPART.segs[status == "peak"]
model.color <- "blue"
peak.y <- -2

ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")

gg <- ggplot()+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    fill=annotation,
    ymin=-Inf, ymax=Inf),
    data=err.df,
    color="grey",
    alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=test.count)+
  geom_step(aes(
    chromStart, mean),
    color=model.color,
    data=FLOPART.segs)+
  geom_point(aes(
    chromEnd, peak.y),
    color=model.color,
    shape=21,
    data=FLOPART.peaks)
show(gg)


# now to do test splits


n.fold <- 2
sets <- list("train","test")
fold.labels <- list()
fold.plots <- list()
err.dt.list <- list()
for(fold in 1:n.fold){

  setDT(test.regions)
  test.regions[, set := ifelse(random.fold==fold, "train", "test")]
  
  fold.labels[[fold]] <- data.table(fold, test.regions)
  
  flopart <- FLOPART::FLOPART(test.count,test.regions[set == "train"], 1400)
  FLOPART.segs <- flopart[["segments_dt"]]
  seg.dt.list <- list(
    FLOPART=FLOPART.segs)
  
  
  for (set in sets){
    for(model in names(seg.dt.list)){
      pkg.segs <- seg.dt.list[[model]][, .(chromStart, chromEnd, mean, status)]
      pkg.peaks <- pkg.segs[status=="peak"]
      err.df <- PeakError::PeakErrorChrom(pkg.peaks, test.regions[set == set])
      err.dt.list[[paste(set,fold)]] <- data.table(set,fold, err.df)
    }
  }
  
  FLOPART.peaks <- FLOPART.segs[status == "peak"]
  
  gg <- ggplot()+
    geom_rect(aes(
      xmin=chromStart, xmax=chromEnd,
      fill=annotation,
      ymin=-Inf, ymax=Inf),
      data=err.df,
      color="grey",
      alpha=0.5)+
    theme_bw()+
    scale_fill_manual(values=ann.colors)+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=test.count)+
    geom_step(aes(
      chromStart, mean),
      color=model.color,
      data=FLOPART.segs)+
    geom_point(aes(
      chromEnd, peak.y),
      color=model.color,
      shape=21,
      data=FLOPART.peaks)
  fold.plots[[fold]] <- gg
  
}

fold.labels.dt <- do.call(rbind, fold.labels)
err.dt <- do.call(rbind, err.dt.list)

show(fold.plots[1])
show(fold.plots[2])



