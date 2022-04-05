library(data.table)
library(ggplot2)
library(tikzDevice)

timings.dt <- data.table::fread("figure-timings-data-labels.csv")
timings.dt[, seconds := time/1e9]
timing.stats <- timings.dt[, .(
  min=min(seconds),
  max=max(seconds),
  median=median(seconds)
), by=.(Algorithm=expr, size)]
ref.dt <- data.table(seconds=1, label="1 second")

algo.colors <- c(
  FPOP="deepskyblue",
  FLOPART="black",
  FLOPARTNoLabel="blue",
  binsegRcpp = "#ECAE5E")

gg <- ggplot()+
  ggtitle(
    "$N=10,000$ data")+
  geom_hline(aes(
    yintercept=seconds),
    data=ref.dt,
    color="grey")+
  geom_text(aes(
    10, seconds, label=label),
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
    "Number of labels $M$",
    limits=c(1, 8000),
    breaks=10^seq(0, 3))+
  scale_y_log10("Computation time (sec.)")+
  theme_bw()
show(gg)

#dl <- directlabels::direct.label(gg, list(cex=0.7, "last.qp"))
#tikz("figure-timings-labels.tex", width=2.5, height=2)
print(dl)
dev.off()

timings.dt <- data.table::fread("figure-timings-data.csv")
timings.dt[, seconds := time/1e9]
timing.stats <- timings.dt[, .(
  min=min(seconds),
  max=max(seconds),
  median=median(seconds)
), by=.(Algorithm=expr, size, data.over.labels, label.density)]
ref.dt <- data.table(
  seconds=1, label="1 second",
  label.density=timing.stats$label.density[1])
gg <- ggplot()+
  ggtitle("FLOPART is $O(N\\log N)$ time
with $O(N)$ positive labels")+
  theme_bw()+
  theme(
    legend.position="none",
    panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ label.density, labeller=label_both)+
  geom_hline(aes(
    yintercept=seconds),
    data=ref.dt[, .(seconds)],
    color="grey")+
  geom_text(aes(
    100, seconds, label=label),
    data=ref.dt,
    size=3,
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
  directlabels::geom_dl(aes(
    size, median, label=Algorithm, color=Algorithm),
    data=timing.stats[size==max(size)],
    method=list(cex=0.75, "last.qp"))+
  geom_text(aes(
    size, median, label=Algorithm, color=Algorithm),
    hjust=1,
    vjust=0,
    size=3,
    data=timing.stats[Algorithm=="OPART"][size==max(size)])+
  scale_x_log10(
    "Number of data $N$",
    limits=10^c(0.5, 7),
    breaks=10^seq(1, 5, by=2))+
  scale_y_log10("Computation time (sec.)")
#tikz("figure-timings.tex", width=3.3, height=2)
print(gg)
dev.off()

