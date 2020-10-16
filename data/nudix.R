library(aphylo)
library(data.table)

model_aphylo <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)

source("sifter/read_sifter.R")

# read.tree("sifter/nudix/PF00293.tree")

tree <- read_nhx("sifter/PF00293.tree")

# Finding the corresponding annotations
ann  <- fread("sifter/nudix/nudix-ann.tsv")

# Are all the leafs in ann included in tree?
(not_in <- which(!ann$protein %in% gsub("[/].+", "", tree$tree$tip.label)))
ann[not_in] # Protein DIPP_ASFB7 with go term 5000000 not present in nudix PFAM 20

ann <- ann[-not_in,]

go <- strsplit(ann$goterms, split=",")
ann <- ann[rep(1:.N, sapply(go, length))]
ann$terms <- unlist(go)
ann[, goterms := NULL]
length(unique(ann$terms)) # I see 68 unique terms, not 66 like mentioned in the paper.

> sort(table(ann$terms))
# 0004170 0004382 0008727 0016887 0019177 0043135 0043262 0047693 0400000 1000000 1020000 
#       1       1       1       1       1       1       1       1       1       1       1 
# 1030000 1500000 1900000 2400000 2500000 3100000 3600000 4100000 4500000 4600000 4800000 
#       1       1       1       1       1       1       1       1       1       1       1 
# 5500000 5600000 6100000 6400000 6800000 7000000 7900000 8400000 8500000 9000000 9400000 
#       1       1       1       1       1       1       1       1       1       1       1 
# 9900000 0047710 1300000 1400000 2000000 2600000 3200000 4000000 4400000 5400000 6000000 
#       1       2       2       2       2       2       2       2       2       2       2 
# 6900000 8200000 8300000 9300000 0047884 3000000 9500000 9700000 3300000 4300000 7300000 
#       2       2       2       2       3       3       3       3       4       4       4 
# 0008796 2300000 3900000 0047840 3800000 5000000 0004081 1600000 0600000 0008828 5900000 
#       5       5       6       7       7       7       8       8       9      10      10 
# 6500000 0047631 
#      11      26 

table(table(ann$terms) > 1)
# Half of the terms have only one annotation. So there's no much we can do
# for cross validation.
# FALSE  TRUE 
# 34    34

dupl <- imputate_duplications(
  tree$tree,
  gsub(".+[_]|[/].+", "", tree$tree$tip.label)
  )

ann[, value:=1L]
anno_i <- dcast(unique(ann), protein ~ terms, value.var = "value")

tree$tree$tip.label <- gsub("[/].+", "", tree$tree$tip.label)

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

pred_table <- cbind(
  predicted = as.vector(pred[ids,]),
  labels    = as.vector(ans$tip.annotation[ids,])
)

# This is the weird bit
pred_table[pred_table==9] <- 0

# Computing AUCs
pred_auc <- auc(pred = pred_table[,"predicted"], labels = pred_table[,"labels"])
pred_auc
# Number of observations     : 12656
# Area Under The Curve (AUC) : 0.70
# Rates can be accessed via the $ operator.
plot(pred_auc)

# Computing accuracy a la SIFTER -----------------------------------------------

threshold <- .5

# Accuracy
tpr <- sapply(ids, function(i) {
  obs <- which(ans$tip.annotation[i,] == 1)
  pre <- which(pred[i,] > threshold)
  length(intersect(obs, pre) ) > 0
}) 

