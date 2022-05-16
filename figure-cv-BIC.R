library(data.table)
library(ggplot2)

err.dt <- data.table(csv=Sys.glob("figure-label-errors-data*/*.csv"))[, {data.table::fread(csv)}, by=csv]
binseg_err.dt <- data.table(csv=Sys.glob("Binseg-figure-label-errors-data*/*.csv"))[, {data.table::fread(csv)}, by=csv]
lopart_err.dt <- data.table(csv=Sys.glob("LOPART-figure-label-errors-data*/*.csv"))[, {data.table::fread(csv)}, by=csv]

err.dt <- rbind(err.dt, binseg_err.dt, lopart_err.dt)
err.dt <- err.dt[order(dataset, sample.id, fold, pen)]
err.dt[, errors := fp + fn]
err.dt[, sequenceID := paste(dataset,"-", sample.id, sep = "")]

err.dt[, log.penalty := log(pen)]
err.dt[, min.log.lambda := log.penalty + c(
  -Inf, -diff(log.penalty)/2
), by=.(
  model, fold, set.i, sequenceID)]

err.dt[, max.log.lambda := c(
  min.log.lambda[-1], Inf
), by=.(
  model, fold, set.i, sequenceID)]

err.test <- err.dt[set.i=="test"]

feature.dt <- data.table::fread(
  "feature-dt.csv"
)

feature.dt[, sequenceID := paste(dataset,"-", sample.id, sep = "")]
feature.mat <- feature.dt[, matrix(
  log.log.data,
  ncol=1,
  dimnames=list(sequenceID=sequenceID, feature="log.log.data"))]

err.train <- err.dt[set.i=="train" & model %in% c("GFPOP", "BINSEG", "LOPART")]

pred.dt <- err.train[, {
  best.penalty <- .SD[, .(
    train.errors=sum(errors)
  ), by=pen][which.min(train.errors)]
  target.dt <- penaltyLearning::targetIntervals(.SD, "sequenceID")
  target.mat <- target.dt[
    rownames(feature.mat),
    cbind(min.log.lambda, max.log.lambda),
    on="sequenceID"]
  keep <- -Inf < target.mat[, 1] | target.mat[,2] < Inf
  fit <- penaltyLearning::IntervalRegressionUnregularized(
    feature.mat[keep, , drop=FALSE], target.mat[keep, ])
  rbind(
    data.table(
      sequenceID=rownames(feature.mat),
      Penalty="constant",
      Parameters=1,
      pred.log.lambda=log(best.penalty$pen)),
    data.table(
      sequenceID=rownames(feature.mat),
      Penalty="BIC",
      Parameters=0,
      pred.log.lambda=as.numeric(feature.mat)),
    data.table(
      sequenceID=rownames(feature.mat),
      Penalty="linear",
      Parameters=2,
      pred.log.lambda=as.numeric(fit$predict(feature.mat))))}, by=.(model, fold)]

auc.dt <- err.test[, {
  select.dt <- data.table(
    fold,
    model = "GFPOP")
  pred.fold <- pred.dt[select.dt, on=names(select.dt)]
  model.dt <- .SD[order(sequenceID, min.log.lambda)]
  pred.fold[, {
    roc.list <- penaltyLearning::ROChange(
      model.dt,
      .SD[select.dt, on = "fold"],
      problem.vars="sequenceID")
    with(roc.list, data.table(
      roc=list(roc), auc,
      thresholds[threshold=="predicted"]))
  }, by=.(Penalty, Parameters)]
}, by=.(fold, model)]


roc.dt <- auc.dt[, data.table(
  roc[[1]]
), by=.(fold, model, Penalty, Parameters)]
possible.dt <- unique(auc.dt[, .(
  fold, possible.fp, possible.fn)])
pred.point.dt <- 
    auc.dt[, .(
      FPR, TPR, fp, tp, auc, labels,
      model, fold, Penalty, Parameters
    )][possible.dt, on=.(fold)]
  
pred.point.dt[, label := paste0(model, ifelse(is.na(auc), "", sprintf(" AUC = %.3f", auc)))]
pred.point.dt[, label.x := rep(1, .N)]
pred.point.dt[, label.y := rep(c(0.8,0.75,0.7, 0.65), each = 3, 2)]

algo.colors <- c(
  GFPOP = "blue",
  FLOPART = "grey50",
  BINSEG = "#ECAE5E",
  LOPART = "red")

gg <- ggplot()+
  theme_bw()+
  scale_color_manual(values=algo.colors)+
  scale_size_manual(values=c(
    BINSEG=1.25,
    GFPOP=1.25,
    LOPART=1.25,
    Flopart=1.5))+
  geom_text((aes(
    x = label.x, y = label.y,
    color=model, hjust=1,
    label= label)),
    size = 3,
    data=pred.point.dt)+
  geom_path(aes(
    FPR, TPR,
    color=model,
    size=model,
    group=paste(model, fold)),
    alpha=0.5,
    data=roc.dt)+
  geom_point(aes(
    FPR, TPR,
    color=model),
    size=3,
    shape=21,
    fill="white",
    data=pred.point.dt)+
  theme(
    panel.spacing=grid::unit(0, "lines"),
    legend.position="none"
  )+
  facet_grid(fold ~ Penalty + Parameters, labeller=label_both)+
  coord_equal()+
  scale_x_continuous(
    "False Positive Rate (test set labels)",
    breaks=c(0, 0.5, 1),
    limits=c(0, 1.2),
    labels=c("0", "0.5", "1"))+
  scale_y_continuous(
    "True Positive Rate (test set labels)",
    labels=c("0.5", "1"),
    breaks=c(0.5, 1)) +
  coord_cartesian(ylim=c(0.5,1)) 

expansion <- 2
pdf("figure-cv-BIC-roc.pdf", width=4.5*expansion, height=2*expansion)
print(gg)
dev.off()

auc.wide <- dcast(
  auc.dt,
  fold + Parameters + Penalty ~ model,
  value.var = "auc")
auc.wide[, diff := GFPOP-Flopart]
auc.wide

pred.point.dt[, fn := possible.fn-tp ]
pred.point.dt[, errors := fn + fp ]
pred.point.dt[, percent.error := 100*errors/labels]
pred.point.dt[, Penalty.Params := paste0(Penalty, ".", Parameters)]
pred.point.dt[, percent.accuracy := 100-percent.error]
pred.point.vars <- melt(
  pred.point.dt,
  measure.vars=c("percent.accuracy", "auc"))
pred.point.diff <- dcast(
  pred.point.vars,
  fold + Penalty.Params + variable ~ model,
  value.var="value")
pred.point.compare <- melt(
  pred.point.diff,
  measure.vars=c("GFPOP"),
  variable.name="baseline")
pred.point.compare[, improvement := Flopart - value]
pred.point.compare[, .(
  min=min(improvement),
  max=max(improvement)
), by=.(baseline, variable)]
pred.point.diff[, GFPOP.diff := Flopart - GFPOP ]
pred.point.diff[order(variable, Penalty.Params, fold), .(
  variable, fold, Penalty.Params, GFPOP.diff)]

pred.point.diff


gg.vars <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    value, Penalty.Params, color=model),
    data=pred.point.vars)+
  facet_grid(fold ~ variable, labeller=label_both, scales="free")+
  scale_color_manual(values=algo.colors)
gg.vars.wide <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    value, Penalty.Params, color=model),
    data=pred.point.vars)+
  facet_grid(. ~ variable + fold, labeller=label_both, scales="free")+
  scale_color_manual(values=algo.colors)
pred.point.wide <- dcast(
  pred.point.dt,
  fold + Penalty.Params ~ model,
  value.var="percent.error")
pred.point.tall <- melt(
  pred.point.wide,
  measure.vars=c("GFPOP", "BINSEG", "LOPART"),
  variable.name="competitor",
  value.name="percent.error")
pred.point.tall

to.plot <- pred.point.dt[Penalty != "BIC"]
to.plot[, fold.id := as.character(fold)]

fold.colors <- c("1" = "blue","2" = "red")

gg <- ggplot()+
  theme_bw()+
  theme(
    legend.position="bottom",
    panel.spacing=grid::unit(0, "lines"))+
  geom_point(aes(
    percent.accuracy, model, color = fold.id),
    shape=18,
    size=4,
    data=to.plot)+
  facet_grid(. ~ Penalty.Params, labeller=label_both)+
  scale_color_manual(
    "Fold",
    values=fold.colors,
    breaks=names(fold.colors)
  )+
  scale_x_continuous(
    "Test accuracy (percent correctly predicted labels)",
    ##limits=c(15, 85),
    breaks=seq(20, 80, by=10)) +
  ylab("Algorithms")


pdf("figure-cv-BIC.pdf", width=6, height=1.8)
print(gg)
dev.off()



