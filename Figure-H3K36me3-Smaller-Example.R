library(ggplot2)
library(data.table)
library(dplyr)


load(url("https://rcdata.nau.edu/genomic-ml/chip-seq-chunk-db/H3K36me3_AM_immune/10/counts.RData"))
load(url("https://rcdata.nau.edu/genomic-ml/chip-seq-chunk-db/H3K36me3_AM_immune/10/regions.RData"))

names(counts)[names(counts) == 'coverage'] <- 'count'

split.counts <- split(counts, counts$sample.id)
split.regions <- split(regions, regions$sample.id)

test.count <- split.counts[["McGill0001"]]
test.regions <- split.regions[["McGill0001"]]


flopart <- FLOPART::FLOPART(test.count, test.regions, 1400)

FLOPART.segs <- flopart[["segments_dt"]]
FLOPART.peaks <- FLOPART.segs[status == "peak"]
model.color <- "blue"
peak.y <- -2

ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")

ggplot()+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    fill=annotation,
    ymin=-Inf, ymax=Inf),
    data=test.regions,
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

max.peaks <- 7L
fit <- PeakSegOptimal::PeakSegPDPAchrom(test.count, max.peaks)
seg.dt <- data.table(fit$segments)
seg.dt.list <- list(
  FLOPART=FLOPART.segs)
for(show.peaks in 5:max.peaks){
  seg.dt.list[[paste(show.peaks, "peaks")]] <- seg.dt[peaks==show.peaks]
}
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
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=test.count)+
  geom_step(aes(
    chromStart, mean),
    color=model.color,
    data=model.segs)+
  geom_segment(aes(
    chromStart, mean,
    xend=chromEnd, yend=mean),
    color=model.color,
    data=model.segs)+
  geom_point(aes(
    chromEnd, peak.y),
    data=model.peaks,
    shape=1,
    color=model.color)+
  geom_segment(aes(
    chromStart, peak.y,
    xend=chromEnd, yend=peak.y),
    data=model.peaks,
    color=model.color,
    size=2)+
  scale_y_continuous("Count of aligned DNA sequence reads")+
  scale_x_continuous("Position on chromosome (bases)")
gg
out.png <- "Figure-H3K36me3-Smaller-Example.png"
png(out.png, 10, 5, units="in", res=200)
print(gg)
dev.off
