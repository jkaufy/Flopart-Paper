library(ggplot2)
library(data.table)
library(dplyr)

library(data.table)
source("Load-All-H3k-Data.R")

maxJumpRule <- function (seg.dt, count.dt){
  seg.dt[, diff := c(diff(mean),NA)]
  seg.dt[, sign := sign(diff)]
  seg.dt[, new.group := c(TRUE,diff(sign) != 0)]
  seg.dt[, group.i := cumsum(new.group)]
  group.dt <- seg.dt[, .SD[which.max(abs(diff)), .(end, diff, sign, mean)], by=group.i]
  
  chromEnd.at.change <- count.dt$chromEnd[group.dt$end]
  
  group.dt[, .(mean = c(mean,1),
    chromStart=c(count.dt$chromStart[1], chromEnd.at.change),
    chromEnd=c(chromEnd.at.change, count.dt[.N, chromEnd]),
    status=rep(if(sign[1]==1)c("background","peak") else c("peak","background"), l=.N+1)
  )]
}

H3K_data <- loadH3KData()

pen <- 50
dataset <- 1
counts <- H3K_data$count[[dataset]]
regions <- H3K_data$labels[[dataset]]

sampleid <- "McGill0001"

split.counts <- split(counts, counts$sample.id)
split.regions <- split(regions, regions$sample.id)

test.count <- split.counts[[sampleid]]
test.regions <- split.regions[[sampleid]]

colnames(test.regions)


flopart <- FLOPART::FLOPART(test.count, test.regions, pen)

setDT(test.count)
test.count[, weight.vec := chromEnd-chromStart]

setDT(test.regions)
test.regions[, changes := ifelse(annotation=="noPeaks", 0, 1)]

lopart_labels <- FLOPART::FLOPART_data(test.count, test.regions)$label_dt

test.regions[, start := lopart_labels$firstRow]
test.regions[, end := lopart_labels$lastRow]



fit <- LOPART::POISSON_LOPART(test.count$count,test.count$weight.vec, 
                      test.regions, pen)

lopart.segs <- maxJumpRule(data.table(fit$segments), test.count)

FLOPART.segs <- flopart[["segments_dt"]]

model.color <- "blue"
peak.y <- -2

ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")

seg.dt.list <- list(
  FLOPART=FLOPART.segs,
  LOPART=lopart.segs)


err.dt.list <- list()
model.segs.list <- list()
for(model in names(seg.dt.list)){
  pkg.segs <- seg.dt.list[[model]][, .(chromStart, chromEnd, mean, status)]
  pkg.peaks <- pkg.segs[status=="peak"]
  err.df <- PeakError::PeakErrorChrom(pkg.peaks, test.regions)
  err.dt.list[[model]] <- data.table(model, err.df)
  model.segs.list[[model]] <- data.table(model, pkg.segs)
}
err.dt <- do.call(rbind, err.dt.list)
model.segs <- do.call(rbind, model.segs.list)

model.peaks <- model.segs[status=="peak"]

fit.segs <- fit$segments
fit.segs[, chromStart := test.count$chromStart[start]]
fit.segs[, chromEnd := test.count$chromEnd[end]]

gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  scale_fill_manual(values=ann.colors)+
  facet_grid(model ~ ., labeller=label_both)+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    ymin=-Inf, ymax=Inf,
    fill=annotation),
    color="grey",
    alpha=0.5,
    data=test.regions)+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    ymin=-Inf, ymax=Inf,
    linetype=status),
    color="black",
    fill=NA,
    size=1,
    data=err.dt)+
  scale_linetype_manual(
    "error type",
    values=c(
      correct=0,
      "false negative"=3,
      "false positive"=1))+
  geom_point(aes(
    chromEnd, count),
    color="grey50",
    data=test.count)+
  geom_segment(aes(
    chromStart+0.5, mean,
    xend=chromEnd+0.5, yend=mean),
    color=model.color,
    data=model.segs)+
  geom_point(aes(
    chromEnd, peak.y),
    data=model.peaks,
    shape=1,
    color=model.color)+
  geom_segment(aes(
    chromStart+0.5, peak.y,
    xend=chromEnd+0.5, yend=peak.y),
    data=model.peaks,
    color=model.color,
    size=2)+
  scale_y_continuous("Count of aligned DNA sequence reads")+
  scale_x_continuous("Position on chromosome (bases)")
gg
out.png <- "Figure-Lopart-Example.png"
png(out.png, 10, 5, units="in", res=200)
print(gg)
dev.off
