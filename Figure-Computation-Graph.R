library(ggplot2)
library(data.table)
library(dplyr)

count.vec <- c(3, 4, 6, 5, 
               6, 30, 4, 
               6, 8, 5, 9,
               6, 9, 10, 16, 18, 17,
               21, 19, 22, 18,
               12, 10, 7, 5,
               6, 3, 4, 6,
               12, 18, 20, 17, 16)
num.labels <- 34
label.size <- 1
break.vec <- seq(1, num.labels, by=1)

counts <- data.table(
  chromStart=as.integer(break.vec-label.size),
  chromEnd=as.integer(break.vec),
  count = as.integer(count.vec))

labels = data.table(
  chromStart = c(4, 22, 30),
  chromEnd = c(7, 25, 34),
  annotation = c("noPeaks", "peakEnd", "peakStart"))
  
ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")

penalty <- 5

flopart <- FLOPART::FLOPART(counts, labels, penalty)
FLOPART.segs <- flopart[["segments_dt"]]

fit <- PeakSegOptimal::PeakSegFPOPchrom(counts, penalty)
seg.dt <- data.table(fit$segments)

node.dt <- data.table::CJ(state=c(-1,1), data.num=1:num.labels)
node.dt[, next.num := data.num+1]
edge.dt <- node.dt[node.dt, on=.(next.num=data.num), nomatch=0L]
edge.dt[, type := ifelse(i.state==state, "no change", "change")]

edge.dt[, type := ifelse(i.state==state, "no change", "change")]


showPossibleLabelPaths <- function (edge.dt, labels) {
  for (row in 1:nrow(labels)){
    current.label <- labels[row]
    
    if (current.label$annotation == "noPeaks"){
      edge.dt[data.num >= current.label$chromStart-1 & 
                data.num < current.label$chromEnd  &
                state == -1, 
              delete := ifelse(i.state==-1, "keep", "delete")]
      edge.dt[data.num >= current.label$chromStart & 
                data.num < current.label$chromEnd  &
                state == 1, 
              delete := "delete"]
      edge.dt[data.num == current.label$chromStart-1 &
                state == 1, 
              delete := ifelse(i.state==-1, "keep", "delete")]
      edge.dt[data.num == current.label$chromEnd  &
                state == 1, 
              delete := "delete"]
    } else if (current.label$annotation == "peakStart"){
      edge.dt[data.num == current.label$chromStart-1, 
              delete := ifelse(i.state==-1, "keep", "delete")]
      edge.dt[data.num >= current.label$chromStart & 
                data.num < current.label$chromEnd  &
                state == 1, 
              delete := ifelse(i.state==-1, "delete", "keep")]
      edge.dt[data.num == current.label$chromStart
              & state == 1, 
              delete := "delete"]
      edge.dt[data.num == current.label$chromEnd 
              & state == -1, 
              delete := "delete"]
      edge.dt[data.num == current.label$chromEnd-1 
              & state == -1, 
              delete := ifelse(i.state==-1, "delete", "keep")]
    } else {
      edge.dt[data.num == current.label$chromStart-1,
              delete := ifelse(i.state==-1, "delete", "keep")]
      edge.dt[data.num >= current.label$chromStart & 
                data.num < current.label$chromEnd  &
                state == -1, 
              delete := ifelse(i.state==1, "delete", "keep")]
      edge.dt[data.num == current.label$chromEnd 
              & state == 1, 
              delete := "delete"]
      edge.dt[data.num == current.label$chromEnd-1 
              & state == 1, 
              delete := ifelse(i.state==1, "delete", "keep")]
      edge.dt[data.num == current.label$chromStart & state == -1,
              delete := "delete"]
    }
    
  }
  
  pathing <- edge.dt$delete
  pathing[is.na(pathing)] <- "keep"
  edge.dt$delete <- pathing
  
  possible.edges <- edge.dt[delete == "keep"]
  
  return (possible.edges)
}


showOptimalPath <- function(possible.edges, segs.dt) {
  for (row in 1:nrow(segs.dt)){
    current.seg <- segs.dt[row]
    
    if (current.seg$status == "background"){
      in.state <- -1
    }else {
      in.state <- 1
    }
    
    possible.edges[data.num >= current.seg$chromStart & 
                     data.num < current.seg$chromEnd &
                     state == in.state, 
                   optimal := ifelse(i.state==in.state, "Optimal", "Possible")]
    if (current.seg$chromStart != 0)
    {
      possible.edges[data.num == current.seg$chromStart-1 &
                       state != in.state,
                     optimal := ifelse(i.state==in.state, "Optimal", "Possible")] 
    }
    
  }
  
  pathing <- possible.edges$optimal
  pathing[is.na(pathing)] <- "Possible"
  possible.edges$optimal <- pathing
  
  return(possible.edges)
}

addLast <- function(seg.dt){
  
  last <- seg.dt[nrow(seg.dt)]
  last[, chromStart := chromEnd+0.5]
  seg.dt <- rbind(seg.dt, last)
  
  return (seg.dt)
}

flopart.edges = showPossibleLabelPaths(edge.dt, labels)
flopart.optimal = showOptimalPath(flopart.edges, FLOPART.segs)
FPOP.optimal = showOptimalPath(edge.dt, seg.dt)

FLOPART.segs <- addLast(FLOPART.segs)
FPOP.segs <- addLast(seg.dt)

seg.dt.list <- list(
  FLOPART=FLOPART.segs,
  FPOP= FPOP.segs)

node.dt.list <- list(
  FLOPART=flopart.optimal,
  FPOP= FPOP.optimal
)

model.segs.list = list()
model.node.list = list()
for(model in names(seg.dt.list)){
  pkg.segs <- seg.dt.list[[model]][, .(chromStart, chromEnd, mean, status)]
  model.segs.list[[model]] <- data.table(model, pkg.segs)
  model.node <- node.dt.list[[model]]
  model.node.list[[model]] <- data.table(model, model.node)
}

model.segs <- do.call(rbind, model.segs.list)
model.nodes <- do.call(rbind, model.node.list)

line.colors <- c(
  Optimal = "blue",
  Possible = "grey50"
)
model.color <- "blue"

ggplot()+
  facet_grid(model ~ ., labeller=label_both)+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    fill=annotation,
    ymin=-Inf, ymax=Inf),
    data=labels,
    color="grey",
    alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values=ann.colors)+
  geom_point(aes(
    chromStart, count),
    color="grey50",
    fill="white",
    data=counts)+
  geom_step(aes(
    chromStart-0.5, mean),
    color=model.color,
    data=model.segs)+
  geom_point(aes(
    data.num, state),
    shape=21,
    size=4,
    data = node.dt,
    fill="white",
    color="black")+
  scale_color_manual(values=line.colors)+
  geom_segment(aes(
    data.num, state,
    linetype=type,
    xend=next.num-0.1, yend=i.state+ifelse(i.state!=state, -i.state*0.1, 0), color = optimal),
    arrow=grid::arrow(length=unit(0.1, "in"), type="closed"),
    data=model.nodes)+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=8)) #change legend text font size






