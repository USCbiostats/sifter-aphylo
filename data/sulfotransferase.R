library(aphylo)
library(data.table)

model_aphylo <- readRDS("data-raw/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)

# source("sifter/read_sifter.R")

tree <- read.tree("data-raw/sulfotransferase/reconciled-pf00685.nhx")
ann  <- data.table(read_pli("data-raw/sulfotransferase/proteinfamily_pf00685.pli"))

# Engelhardt (2011):
# We use the annotations with the evidence codes IDA, IMP, and TAS
ann <- ann[moc %in% c("IDA", "IMP", "TAS")]

# They also excluded GO:0000166
ann <- ann[go != "0000166"]

# They excluded the proteins that only had sulfotransferase annotations:
# GO:0008146
ann[, n := .N, by = name]
ann <- ann[n > 1 | (n == 1 & go != "0008146")]


tree$node.label <- sprintf("node%03i", seq_len(ape::Nnode(tree)))

dupl <- imputate_duplications(
  tree,
  gsub(".+[_]|[/].+", "", tree$tip.label)
)

ann[, value:=1L]
ann[, c("moc", "number") := NULL]
anno_i <- dcast(unique(ann), name ~ go, value.var = "value")

# Reading the aphylo tree
ans <- aphylo_from_data_frame(
  tree        = tree,
  annotations = as.data.frame(anno_i),
  types       = data.frame(
    tree$node.label,
    !dupl[-c(1:ape::Ntip(tree))]
  )
)

# Predictions -----------------------------------------------------------------

ids <- which(rowSums(ans$tip.annotation) != Nann(ans) * 9)
pred <- predict(model_aphylo, newdata = ans, ids = ids, loo = TRUE)

pred_table <- cbind(
  predicted = as.vector(pred[ids,]),
  labels    = as.vector(ans$tip.annotation[ids,])
)

# This is the weird bit
pred_table[pred_table==9] <- 0

# Computing AUCs
pred_auc <- auc(pred = pred_table[,"predicted"], labels = pred_table[,"labels"])
pred_auc

# Number of observations     : 1960
# Area Under The Curve (AUC) : 0.75
# Rates can be accessed via the $ operator
plot(pred_auc)

# Computing accuracy a la SIFTER -----------------------------------------------

# Directly from the paper:
# To compute accuracy for a given data set, we counted the number
# of proteins for which the top-ranked prediction exactly matched
# the experimental annotation, and divided by the total number of
# proteins. In the case of multiple experimental annotations or top-
#   ranked predictions, we counted the protein as having an accurate
# prediction when the intersection of the two sets was not empty. 

acc_at_least_one <- lapply(ids, function(i) {
  top_pred <- which(abs(pred[i,] - max(pred[i,])) < 1e-5)
  true_ann <- which(ans$tip.annotation[i,] == 1)
  data.frame(
    pred = paste(top_pred, collapse=","),
    true = paste(true_ann, collapse=","),
    accurate = ifelse(length(intersect(top_pred, true_ann)), 1, 0)
  )
})
acc_at_least_one <- do.call(rbind, acc_at_least_one)
rownames(acc_at_least_one) <- rownames(pred[ids,])
table(acc_at_least_one$accurate)
#  0  1 
# 10 23

acc_at_least_one
