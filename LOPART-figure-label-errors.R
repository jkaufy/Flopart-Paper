library(data.table)
library(ggplot2)

err.dt <- data.table(csv=Sys.glob("figure-label-errors-data*/*.csv"))[, {data.table::fread(csv)}, by=csv]
err.dt <- err.dt[model == "Flopart"]

lopart.err.dt <- data.table(csv=Sys.glob("LOPART-figure-label-errors-data*/*.csv"))[, {data.table::fread(csv)}, by=csv]

err.dt <- rbind(err.dt, lopart.err.dt)

err.dt <- err.dt[order(dataset, sample.id, fold, pen)]
err.dt[, errors := fp + fn]

algo.colors <- c(
  LOPART = "red",
  FLOPART = "black"
)

scale.for <- function(algo){
  scale_fill_gradient(
    "log10(seqs)",
    low="white",
    high=algo.colors[[algo]])
}  

total.dt <- dcast(
  data.table(data.frame(err.dt)),
  dataset + sample.id + fold + model + pen ~ set.i,
  value.var="errors")

total.dt[, train.test := test+train]

total.min <- total.dt[, .SD[
  which.min(train.test)], by=.(dataset, sample.id, fold, model)]

total.min.wide <- dcast(
  total.min,
  dataset + sample.id + fold ~ model,
  value.var=c("test", "train", "train.test"))

total.min.wide[, diff := train.test_LOPART-train.test_Flopart]

total.min.wide[, train.test.diff := train.test_LOPART - train.test_Flopart]


mytab <- function(dt, col.name){
  errors <- dt[, .(
    count=.N,
    percent=100*.N/nrow(dt)
  ), by=col.name]
  is.zero <- errors[[col.name]] == 0
  is.pos <- errors[[col.name]] > 0
  is.neg <- errors[[col.name]] < 0
  pos <- errors[is.pos]
  neg <- errors[is.neg]
  sum.wide <- data.table(
    sum.count=sum(errors$count),
    zero.count=errors$count[is.zero],
    pos.count=sum(pos[["count"]]),
    pos.min=min(pos[[col.name]]),
    pos.max=max(pos[[col.name]]),
    neg.count=sum(neg[["count"]]),
    neg.min=min(neg[[col.name]]),
    neg.max=max(neg[[col.name]]))
  sum.tall <- melt(sum.wide, measure.vars=names(sum.wide))
  sum.tall[grepl("count", variable), percent := 100*value/nrow(dt) ]
  list(
    errors=errors,
    summary=sum.tall)
}

mytab(total.min.wide, "train_LOPART")

total.min.wide[, test.diff_LOPART := test_LOPART-test_Flopart]

mytab(total.min.wide, "test.diff_LOPART")

testing <- total.min[train > 0]

train.test.counts <- total.min.wide[, .(
  splits=.N
), by=.(train_LOPART, test.diff_LOPART)]


gg <- ggplot()+
  ggtitle("Best case comparison
with LOPART")+
  geom_hline(yintercept=0, color="grey")+
  geom_vline(xintercept=0, color="grey")+
  geom_tile(aes(
    train_LOPART, test.diff_LOPART, fill=log10(splits)),
    alpha=0.8,
    data=train.test.counts)+
  geom_text(aes(
    train_LOPART, test.diff_LOPART, label=splits),
    data=train.test.counts)+
  scale.for("LOPART")+
  coord_equal()+
  theme_bw()+
  scale_x_continuous(
    "LOPART train label errors
(FLOPART is always=0)")+
  scale_y_continuous(
    "Test label error difference
(LOPART-FLOPART)") + theme(text = element_text(size = 14)) 

pdf("LOPART-figure-label-errors.pdf", width=4, height=4)
print(gg)
dev.off()
show(gg)


# choosing pen which minimizes train first, then train second
total.dt <- total.dt[order(dataset, sample.id, model, fold, train,test)]
rank <- total.dt[, rank := order(train,test), by=.(dataset, sample.id, model, fold)]

total.min <- rank[, .SD[
  which.min(rank)], by=.(dataset, sample.id, fold, model)]

testing <- total.min[train > 0]

total.min.wide <- dcast(
  total.min,
  dataset + sample.id + fold ~ model,
  value.var=c("test", "train", "train.test"))

total.min.wide[, diff := train.test_LOPART-train.test_Flopart]

total.min.wide[, train.test.diff := train.test_LOPART - train.test_Flopart]

mytab(total.min.wide, "train_LOPART")

total.min.wide[, test.diff_LOPART := test_LOPART-test_Flopart]

mytab(total.min.wide, "test.diff_LOPART")

train.test.counts <- total.min.wide[, .(
  splits=.N
), by=.(train_LOPART, test.diff_LOPART)]


gg <- ggplot()+
  ggtitle("Best case comparison
with LOPART")+
  geom_hline(yintercept=0, color="grey")+
  geom_vline(xintercept=0, color="grey")+
  geom_tile(aes(
    train_LOPART, test.diff_LOPART, fill=log10(splits)),
    alpha=0.8,
    data=train.test.counts)+
  geom_text(aes(
    train_LOPART, test.diff_LOPART, label=splits),
    data=train.test.counts)+
  scale.for("LOPART")+
  coord_equal()+
  theme_bw()+
  scale_x_continuous(
    "LOPART train label errors
(FLOPART is always=0)")+
  scale_y_continuous(
    "Test label error difference
(LOPART-FLOPART)")


pdf("LOPART-figure-label-errors-min-train.pdf", width=8, height=8)
print(gg)
dev.off()
show(gg)





