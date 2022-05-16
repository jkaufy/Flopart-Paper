library(data.table)
library(ggplot2)

timings.dt <- data.table::fread("figure-real-timings-data.csv")
timings.dt <- timings.dt[expr %in% c("GFPOP", "BINSEG", "LOPART", "FLOPART")]
timings.dt[, seconds := time/1e9]
timing.stats <- timings.dt[, .(
  min=min(seconds),
  max=max(seconds),
  median=median(seconds)
), by=.(Algorithm=expr, size=as.integer(exp(round(log(size)))))]
ref.dt <- data.table(seconds=c(1, 60), label=c("1 second", "1 minute"))

algo.colors <- c(
  GFPOP="deepskyblue",
  FLOPART="darkgrey",
  BINSEG = "orange",
  LOPART = "red")


gg <- ggplot()+
  ggtitle(
    "ChIP-seq Genomic Timings")+
  geom_hline(aes(
    yintercept=seconds),
    data=ref.dt,
    color="grey")+
  geom_text(aes(
    8000, seconds, label=label),
    size=3,
    data=ref.dt,
    vjust=-0.5,
    color="grey50")+
  geom_ribbon(aes(
    size, ymin=min, ymax=max, fill=Algorithm),
    data=timing.stats,
    alpha=0.5)+
  scale_fill_manual(values=algo.colors)+
  scale_color_manual(values=algo.colors)+
  geom_line(aes(
    size, median, color=Algorithm),
    size=2,
    data=timing.stats)+
  scale_x_log10(
    "Size of Dataset",
    limits = c(NA, 1e6))+
  scale_y_log10("Computation time (sec.)")+
  theme_bw()

labelled <- directlabels::direct.label(gg, "right.polygons")

pdf("figure-real-timing.pdf", width=9, height=3.4)
print(labelled)
dev.off()
show(labelled)
