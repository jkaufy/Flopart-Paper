library(data.table)
library(ggplot2)
library(tikzDevice)

timings.dt <- data.table::fread("figure-real-timings-data-small-pen.csv")
timings.dt <- timings.dt[expr %in% c("FPOP", "binsegRcpp", "LOPART", "FLOPART")]
timings.dt[, seconds := time/1e9]
timing.stats <- timings.dt[, .(
  min=min(seconds),
  max=max(seconds),
  median=median(seconds)
), by=.(Algorithm=expr, size=as.integer(exp(round(log(size)))))]
ref.dt <- data.table(seconds=c(1, 60), label=c("1 second", "1 minute"))

algo.colors <- c(
  FPOP="deepskyblue",
  FLOPART="blue",
  binsegRcpp = "#ECAE5E",
  LOPART = "red")


gg <- ggplot()+
  ggtitle(
    "H3K Genomic Data Timings")+
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
    "Size of Dataset $N$")+
  scale_y_log10("Computation time (sec.)")+
  theme_bw()
library(directlabels)
direct.label(gg, method = "Algorithm")
show(gg)

dl <- directlabels::direct.label(gg, list(cex=0.7, "last.qp"))
tikz("figure-timings-labels.tex", width=5, height=3)
print(dl)
dev.off()