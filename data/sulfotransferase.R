library(aphylo)
library(data.table)

model_aphylo <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)

source("sifter/read_sifter.R")

tree <- read_nhx("sifter/sulfotransferase/reconciled-pf00685.nhx")
ann  <- read_pli("sifter/sulfotransferase/proteinfamily_pf00685.pli")
ann <- ann[moc == c("IDA")]
dupl <- imputate_duplications(tree$tree)

ann[, value:=1L]
ann[, c("moc", "number") := NULL]
anno_i <- dcast(unique(ann), name ~ go, value.var = "value")

# Reading the aphylo tree
ans <- aphylo_from_data_frame(
  tree        = tree$tree,
  annotations = as.data.frame(anno_i),
  types       = data.frame(
    tree$tree$node.label,
    !dupl[-c(1:ape::Ntip(tree$tree))]
  )
)

# Predictions -----------------------------------------------------------------

ids <- which(rowSums(ans$tip.annotation) != Nann(ans) * 9)
pred <- predict(model_aphylo, newdata = ans, ids = ids, loo = TRUE)

pred <- cbind(
  predicted = as.vector(pred[ids,]),
  labels    = as.vector(ans$tip.annotation[ids,])
)

# This is the weird bit
pred[pred==9] <- 0

# Computing AUCs
pred_auc <- auc(pred = pred[,"predicted"], labels = pred[,"labels"])
pred_auc
# Number of observations     : 12656
# Area Under The Curve (AUC) : 0.70
# Rates can be accessed via the $ operator.
plot(pred_auc)