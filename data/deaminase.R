library(aphylo)
library(data.table)

model_aphylo <- readRDS("data-raw/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)

tree <- read.tree("data-raw/deaminase/reconciled-pf00962-paup.nex")
ann  <- data.table(read_pli("data-raw/deaminase/proteinfamily_pf00962n.pli"))

# Engelhardt (2011):
# We use the annotations with the evidence codes IDA, IMP, and TAS
ann <- ann[moc %in% c("IDA", "IMP", "TAS")]

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

# Measuring correct (TPR) at specificity (TNR) ~ 99%
(pred_auc$tpr)[which.min(abs(pred_auc$tnr - .99))] # 0.3684211

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

acc_at_least_one_alt <- accuracy_sifter(
  pred = pred_alt[ids, ],
  lab  = ans_alt$tip.annotation[ids,],
  tol  = 1e-15,
  highlight = "\\cellcolor{blue!25}\\textbf{%s}"
)

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
  file = "data/deaminase.rds"
)


# Writing table ----------------------------------------------------------------
acc_at_least_one$Gene <- gsub("[_]", "\\\\_", acc_at_least_one$Gene)
acc_at_least_one$Predicted <- gsub(",", ", ", acc_at_least_one$Predicted)
acc_at_least_one$Observed <- gsub(",", ", ", acc_at_least_one$Observed)

fn <- "data/deaminase.tex"
if (file.exists(fn))
  file.remove(fn)

cat("\\begin{table}[!htbp]
    \\begin{tabular}{lm{.2\\linewidth}<\\raggedleft m{.2\\linewidth}<\\raggedleft}\\toprule\n",
    file = fn)

write.table(
  acc_at_least_one[,-4], sep = "&", eol = "\\\\\n",
  quote = FALSE, row.names = FALSE, file = fn, append = TRUE
)
cat(
  "\\bottomrule\n\\end{tabular}",
  sprintf(
    "\\caption{\\label{tab:deaminase}List of predicted vs experimental annotations. GO terms that are
    consistent with true annotations. In this case, this family included a total
    of %i possible functions. aphylo correctly predicted %i of %i proteins.}",
    ncol(pred), sum(acc_at_least_one$Accuracy), nrow(acc_at_least_one)
  ),
  "\\end{table}",
  file   = fn,
  append = TRUE,
  sep    = "\n"
)
