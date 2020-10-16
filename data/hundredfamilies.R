library(ape)
library(aphylo)
library(stringi)
library(data.table)

# Functions to read sifter trees
source("sifter/read_sifter.R")

# Parameters
n      <- 100
ncores <- 2

# Reading trees
fn_nhx <- list.files("sifter/hundred/", pattern = "nhx", full.names = TRUE)
fn_pli <- gsub("[.]nhx$", ".pli", fn_nhx)
tn       <- gsub("(.+/)([a-zA-Z0-9]+)\\.nhx", "\\2", fn_nhx, perl = TRUE)

# Getting the sequence by increasing size
read_seq <- order(file.size(fn_nhx), decreasing = FALSE)
trees2read <- read_seq[1:n]

# Step 1: read the trees, including annotations --------------------------------
trees <- parallel::mclapply(fn_nhx[trees2read], read_nhx, mc.cores = ncores)
message("Reading trees, done.")

anno  <- parallel::mclapply(fn_pli[trees2read], function(n) {
  ans <- read_pli(n, dropNAs = FALSE)
  # ans <- ans[!is.na(go)]
  }, mc.cores = ncores)

message("Reading annotations, done.")

rows_per_tree <- sapply(anno, nrow)
anno  <- cbind(
  rbindlist(anno),
  tree_name = rep(tn[trees2read], rows_per_tree)
)

# Step 2: Impute the events ----------------------------------------------------
events <- parallel::mclapply(trees, function(i) {
  imputate_duplications(i$tree)
}, mc.cores = ncores)

message("Imputing events, done.")

# Step 3: Build aphylo objects -------------------------------------------------

aphylo_trees <- parallel::mclapply(seq_len(n), function(i) {
  
  # Preparing the annotations
  anno_i <- anno[tree_name == tn[trees2read][i]][
    , c("tree_name", "moc", "number") := NULL]
  
  # Removing NAs
  anno_i <- anno_i[!is.na(go)]
  
  anno_i[, value:=1L]
  anno_i <- dcast(unique(anno_i), name ~ go, value.var = "value")
  
  # Reading the aphylo tree
  ans <- aphylo_from_data_frame(
    tree        = trees[[i]]$tree,
    annotations = as.data.frame(anno_i),
    types       = data.frame(
      trees[[i]]$tree$node.label,
      !events[[i]][-c(1:ape::Ntip(trees[[i]]$tree))]
    )
  )
  
  message("Tree ", trees2read[i], " (", i, ") done.")
  
  ans
  
}, mc.cores = ncores)

names(aphylo_trees) <- tn[trees2read][1:n]

fwrite(anno, file = "sifter/hundredfamilies.tar.gz", compress = "gzip")
saveRDS(aphylo_trees, "sifter/hundredfamilies.rds", compress = TRUE)

message("All done.")

