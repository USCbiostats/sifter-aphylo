library(aphylo)
library(data.table)

model_aphylo <- readRDS("data-raw/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)

tree <- read.tree("data-raw/nudix/PF00293.tree")
tree$node.label <- sprintf("node%03i", seq_len(ape::Nnode(tree)))

# There are duplicates of Q9RXI4_DEIRA and Q9RYE5_DEIRA. Not sure which one is
# SIFTER using. Different regions in the seq.

# Finding the corresponding annotations
ann  <- fread("data-raw/nudix/nudix-ann.tsv")

# Are all the leafs in ann included in tree?
(not_in <- which(!ann$protein %in% gsub("[/].+", "", tree$tip.label)))
# Protein DIPP_ASFB7 with go term 5000000 not present in nudix PFAM 20
# https://www.uniprot.org/uniprot/P32092
ann[not_in]

ann <- ann[-not_in,]

go <- strsplit(ann$goterms, split=",")
ann <- ann[rep(1:.N, sapply(go, length))]
ann$terms <- unlist(go)
ann[, goterms := NULL]
length(unique(ann$terms)) # I see 68 unique terms, not 66 like mentioned in the paper.

sort(table(ann$terms))
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
  tree,
  gsub(".+[_]|[/].+", "", tree$tip.label)
  )

ann[, value:=1L]
anno_i <- dcast(unique(ann), protein ~ terms, value.var = "value")

tree$tip.label <- gsub("[/].+", "", tree$tip.label)

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

ids  <- which(rowSums(ans$tip.annotation) != Nann(ans) * 9)
system.time(
  pred <- predict(model_aphylo, newdata = ans, ids = ids, loo = TRUE)
)
# user  system elapsed 
# 55.408   0.233  55.797 

#  user  system elapsed 
# 2.567   0.028   2.597 

pred_table <- cbind(
  predicted = as.vector(pred[ids,]),
  labels    = as.vector(ans$tip.annotation[ids,])
)

# This is the weird bit
pred_table[pred_table==9] <- 0

# Computing AUCs
pred_auc <- auc(pred = pred_table[,"predicted"], labels = pred_table[,"labels"])
pred_auc
#' Number of observations     : 6664
#' Area Under The Curve (AUC) : 0.57
#' Rates can be accessed via the $ operator.
plot(pred_auc)

# Measuring correct (TPR) at specificity (TNR) ~ 99%
(pred_auc$tpr)[which.min(abs(pred_auc$tnr - .99))] 
# [1] 0.1578947

ans_alt <- ans
ans_alt$tip.annotation[ids,][ans_alt$tip.annotation[ids,] == 9] <- 0
pred_alt <- predict(model_aphylo, newdata = ans_alt, ids = ids, loo = TRUE)

pred_table_alt <- cbind(
  predicted = as.vector(pred_alt[ids,]),
  labels    = as.vector(ans_alt$tip.annotation[ids,])
)

pred_auc_alt <- auc(pred = pred_table_alt[,"predicted"], labels = pred_table_alt[, "labels"])

# Computing accuracy a la SIFTER -----------------------------------------------

# Directly from the paper:
# To compute accuracy for a given data set, we counted the number
# of proteins for which the top-ranked prediction exactly matched
# the experimental annotation, and divided by the total number of
# proteins. In the case of multiple experimental annotations or top-
#   ranked predictions, we counted the protein as having an accurate
# prediction when the intersection of the two sets was not empty. 

acc_at_least_one <- accuracy_sifter(
  pred = pred[ids, ],
  lab  = ans$tip.annotation[ids,],
  tol  = 1e-15,
  highlight = "\\cellcolor{blue!25}\\textbf{%s}"
  )

# Q9RXI4_DEIRA and Q9RYE5_DEIRA are repeated
id <- which(rownames(pred)[ids] == "Q9RXI4_DEIRA")[1]
id <- c(which(rownames(pred)[ids] == "Q9RYE5_DEIRA")[1], id)
acc_at_least_one <- acc_at_least_one[-id,]

table(acc_at_least_one$Accuracy)
#  0  1 
# 63 33 

acc_at_least_one_alt <- accuracy_sifter(
  pred = pred_alt[ids, ],
  lab  = ans_alt$tip.annotation[ids,],
  tol  = 1e-15,
  highlight = "\\cellcolor{blue!25}\\textbf{%s}"
)

acc_at_least_one_alt <- acc_at_least_one_alt[-id,]

# Saving data
saveRDS(
  list(
    aucs       = pred_auc,
    aucs_alt   = pred_auc_alt,
    pred       = pred,
    pred_table = pred_table,
    accuracy   = acc_at_least_one,
    accuracy_alt = acc_at_least_one_alt
  ),
  file = "data/nudix.rds"
)


# Writing table ----------------------------------------------------------------
acc_at_least_one$Gene <- gsub("[_]", "\\\\_", acc_at_least_one$Gene)
acc_at_least_one$Predicted <- gsub(",", ", ", acc_at_least_one$Predicted)
acc_at_least_one$Observed <- gsub(",", ", ", acc_at_least_one$Observed)

fn <- "data/nudix.tex"
if (file.exists(fn))
  file.remove(fn)

cat("\\begin{table}[!htbp]
    \\begin{tabular}{lm{.2\\linewidth}<\\raggedleft m{.2\\linewidth}<\\raggedleft}\\toprule ",
    file = fn)

write.table(
  acc_at_least_one[,-4], sep = "&", eol = "\\\\\n",
  quote = FALSE, row.names = FALSE, file = fn, append = TRUE
  )
cat(
  "\\bottomrule\n\\end{tabular}",
  sprintf(
    "\\caption{\\label{tab:nudix}List of predicted vs experimental annotations. GO terms that are
    consistent with true annotations. In this case, this family included a total
    of %i possible functions. aphylo correctly predicted %i of %i proteins.}",
    ncol(pred), sum(acc_at_least_one$Accuracy), nrow(acc_at_least_one)
    ),
  "\\end{table}",
  file   = fn,
  append = TRUE,
  sep    = "\n"
  )
