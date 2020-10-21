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

# Directly from the paper:
# To compute accuracy for a given data set, we counted the number
# of proteins for which the top-ranked prediction exactly matched
# the experimental annotation, and divided by the total number of
# proteins. In the case of multiple experimental annotations or top-
#   ranked predictions, we counted the protein as having an accurate
# prediction when the intersection of the two sets was not empty. 

acc_at_least_one <- lapply(ids, function(i) {
  top_pred <- which(abs(pred[i,] - max(pred[i,])) < 1e-10)
  true_ann <- which(ans$tip.annotation[i,] == 1)
  data.frame(
    pred = paste(top_pred, collapse=","),
    true = paste(true_ann, collapse=","),
    accurate = ifelse(length(intersect(top_pred, true_ann)), 1, 0)
  )
})

acc_at_least_one <- do.call(rbind, acc_at_least_one)

# Q9RXI4_DEIRA is repeated
id <- which(rownames(pred)[ids] == "Q9RXI4_DEIRA")[1]
id <- c(which(rownames(pred)[ids] == "Q9RYE5_DEIRA")[1], id)
acc_at_least_one <- acc_at_least_one[-id,]

# rownames(acc_at_least_one) <- rownames(pred)[ids][-id]

table(acc_at_least_one$accurate)
#  0  1 
# 63 33 

# acc_at_least_one
#               pred                true accurate
# 16              53                  45        0
# 31              50                  45        0
# 32              45               45,50        1
# 41              53            45,50,56        0
# 42              45            25,45,54        1
# 43           25,54                  45        0
# 59              45                  53        0
# 60              53                  45        0
# 115             50                   1        0
# 122             53                   1        0
# 157             45       6,14,17,24,53        0
# 182             33                   6        0
# 193              6               33,37        0
# 223             45               11,49        0
# 224              1               11,49        0
# 227             49               11,49        1
# 231              1                  11        0
# 322             37                   1        0
# 331             37            36,37,44        1
# 350             37       1,33,36,37,59        1
# 386              5                   1        0
# 419             36                1,36        1
# 423             36             1,36,37        1
# 434              1                  36        0
# 495             37               36,37        1
# 565             37          6,29,30,41        0
# 621  3,10,52,61,65                  63        0
# 668            7,9        5,6,14,17,53        0
# 670          47,62                  66        0
# 671             66            47,62,66        1
# 685              5                  66        0
# 717             24                  24        1
# 745              5                  24        0
# 921             42            17,24,53        0
# 971             17                  42        0
# 987             66             2,14,17        0
# 991            7,9                   5        0
# 1000             5      1,7,9,36,37,59        0
# 1114            21               11,26        0
# 1118            53                  11        0
# 1134            26                  11        0
# 1162            53                  11        0
# 1192            26            11,21,26        1
# 1270            53               34,40        0
# 1281            11                  34        0
# 1290            11                  40        0
# 1304            40                  40        1
# 1321            12                5,53        0
# 1421            23                  40        0
# 1433            11                  27        0
# 1439            40               23,46        0
# 1473            11          6,14,17,53        0
# 1505            14                  11        0
# 1530            49       6,12,17,24,53        0
# 1556            11                  34        0
# 1557            34            11,27,34        1
# 1558            11                  34        0
# 1559            34            11,27,34        1
# 1608            34            11,13,15        0
# 1652            68      15,49,55,58,67        0
# 1722            68                  49        0
# 1761            49               49,68        1
# 1790            49            49,55,67        1
# 1800            67                  49        0
# 1841            49                6,53        0
# 1866            24          6,14,17,53        0
# 1881            24                  24        1
# 1883            24                  24        1
# 1886            24         24,30,35,43        1
# 1934            53               14,17        0
# 1940             6         19,32,39,51        0
# 1948   19,32,39,51          6,14,17,53        0
# 2087            63 3,10,41,46,52,61,65        0
# 2108             8                6,53        0
# 2114             6                   8        0
# 2129             5                   5        1
# 2139             5                   5        1
# 2334             8                  57        0
# 2371      16,60,64      22,28,57,60,64        1
# 2443         60,64                  57        0
# 2498      28,60,64         16,57,60,64        1
# 2525             8                  22        0
# 2785            11                  11        1
# 2798            11                  11        1
# 2878            11                  11        1
# 2913            11                  11        1
# 2917            11                  11        1
# 2919            11         11,15,27,49        1
# 2992            11               21,27        0
# 3029            13                  11        0
# 3058            27               11,27        1
# 3086            11            11,31,38        1
# 3088            11               11,31        1
# 3105            11                  11        1
# 3147            13            11,31,38        0
# 3168            38   11,13,20,48,49,67        0
# 3267            45                  11        0
# 3300            13                4,18        0

acc_all <- sapply(ids, function(i) {
  top_pred <- which(pred[i,] == max(pred[i,]))
  true_ann <- which(ans$tip.annotation[i,] == 1)
  if (length(top_pred) != length(true_ann))
    return(0)
  ifelse(all(top_pred %in% true_ann), 1, 0)
})

mean(acc_all)


