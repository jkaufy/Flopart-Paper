# need to run Download-All-H3K-Data.R first
# Load in all data sets
library(data.table)

maxJumpRule <- function (seg.dt, count.dt){
  seg.dt[, diff := c(diff(mean),NA)]
  seg.dt[, sign := sign(diff)]
  seg.dt[, new.group := c(TRUE,diff(sign) != 0)]
  seg.dt[, group.i := cumsum(new.group)]
  group.dt <- seg.dt[, .SD[which.max(abs(diff)), .(end, diff, sign, mean)], by=group.i]
  
  chromEnd.at.change <- count.dt$chromEnd[group.dt$end]
  
  if (nrow(group.dt)>0){
    group.dt[, .(mean = c(mean,1),
       chromStart=c(count.dt$chromStart[1], chromEnd.at.change),
       chromEnd=c(chromEnd.at.change, count.dt[.N, chromEnd]),
       status=rep(if(sign[1]==1)c("background","peak") else c("peak","background"), l=.N+1)
    )]
  } else
  {
    group.dt <- data.table( mean = numeric(), chromStart=numeric(), chromEnd=numeric(), status=character())
  }
}

source("Load-All-H3k-Data.R")

H3K_data <- loadH3KData()

# increase above 10 log scale
penalties <- 10^seq(-5, 6, by=0.5)
grid.dt <- data.table(penalty=penalties)
grid.dt[, penalty0 := penalty]
n.fold <- 2

sets <- list("train","test")
model <- "BINSEG"

cache.prefix <- "Binseg-figure-label-errors-data"

for(dataset in 1:length(H3K_data$count)){
  print(dataset)
  
  one_count <- H3K_data$count[[dataset]]
  one_label <- H3K_data$labels[[dataset]]
  
  sample_split_count <- split(one_count, one_count$sample.id)
  sample_split_label <- split(one_label, one_label$sample.id)
  
  if (!('peaks' %in% one_label$annotation)){
    for (sample.id in names(sample_split_count)){
      sample.err.list <- list()
      
      cache.save <- paste(dataset, "-", sample.id, ".csv", sep = "")
      cache.file <- file.path(cache.prefix, cache.save)

      one_sample_count <- sample_split_count[[sample.id]]
      
      if(!file.exists(cache.file)){
        print(sample.id)
        
        one_sample_label <- sample_split_label[[sample.id]]
        
        setDT(one_sample_count)
        one_sample_count[, weight.vec := chromEnd-chromStart]
        
        fit <- binsegRcpp::binseg('poisson', one_sample_count$count, 
                             weight.vec = one_sample_count$weight.vec)
        
        
        modelSelection.dt <- data.table(penaltyLearning::modelSelection(
          fit$splits, "loss", "segments"))[, .(min.lambda, max.lambda, segments)]
        
        setkey(grid.dt, penalty, penalty0)
        setkey(modelSelection.dt, min.lambda, max.lambda)
        pen.dt <- foverlaps(grid.dt, modelSelection.dt)
        
        for (pens in 1:nrow(pen.dt)){
          pen <- pen.dt$penalty[pens]
          max.peaks <- pen.dt$segments[pens]
          
          toMaxJump <- coef(fit, max.peaks)
          
          seg.dt <- maxJumpRule(toMaxJump, one_sample_count)
          
          pkg.segs <- seg.dt[, .(chromStart, chromEnd, mean, status)]
          pkg.peaks <- pkg.segs[status=="peak"]

          for (fold in 1:n.fold){
            
            setDT(one_sample_label)
            one_sample_label[, set := ifelse(fold==random.fold, "test", "train")]
            
            for (set.i in sets) {
              
              err.dt <- PeakError::PeakErrorChrom(pkg.peaks, 
                                                  one_sample_label[set == set.i])
              
              sample.err.list[[paste(dataset, sample.id, pen, fold, model, set.i)]] <- data.table(
                dataset,
                sample.id,
                pen,
                fold,
                model,
                set.i,
                possible.fp = sum(err.dt$possible.fp),
                fp = sum(err.dt$fp),
                possible.fn = sum(err.dt$possible.tp),
                fn = sum(err.dt$fn),
                labels = nrow(one_sample_label[set == set.i])
              )
              
            }
          }
          
        }
        
        sample.err.dt <- do.call(rbind, sample.err.list)
        dir.create(dirname(cache.file), showWarnings = FALSE, recursive = TRUE)
        data.table::fwrite(sample.err.dt, cache.file)
        
      }
    }
  }
}
  

