library(aphylo)
library(data.table)
library(pfamscanr)

# Loading the trees for the experiment
candidate_trees <- readRDS("data-raw/candidate_trees.rds")

# Extracting database with annotations per tree
dat <- lapply(candidate_trees, function(tree) {
  idx <- which(tree$tip.annotation != 9)
  data.frame(
    go        = colnames(tree$tip.annotation),
    UniProtKB = gsub("UniProtKB=", "", tree$tree$tip.label[idx]),
    qualifier = tree$tip.annotation[idx]
  )
})

dat <- Map(
  function(d, n) cbind(d, tree = n),
  d = dat, n = names(candidate_trees)
  )

dat <- rbindlist(dat)

# ans <- pfamscanr::get_info_uniprot("P00789")

pfam_ids <- unique(dat[,as.character(UniProtKB)])
pfam_ids <- get_info_uniprot(pfam_ids)

# The following IDS were not found in PFAM
# names(pfam_ids)[c(259, 642)]
# [1] "F1QJU1"     "A0A0R4ITP9"

save.image("data/aphylo_families_entries.rda")

# Setting up the data to be used with SIFTER -----------------------------------

ans <- unique(dat[, .(tree, UniProtKB)])[, .(n=.N), by = tree]
ans[order(n),] # Selecting one tree with small number of proteins annotated

# Extracting the genes (with full names) associated the tree
# PTHR10082
idx <- dat[, which(tree == "PTHR24416")]

sample_trees <- pfam_ids # [idx]

cat(names(sample_trees), sep = "\n")

# Creating the fasta sequences for querying pfamscan ---------------------------
counter <- 0L
fasta <- lapply(sample_trees, function(i) {
  
  counter <<- counter + 1L
  
  d <- tryCatch(strsplit(
    i$pfam$entry$sequence[[1]],
    split = "(?<=.{50})",
    perl = TRUE)[[1]],
    error = function(e) e
  )
  
  if (inherits(d, "error")) {
    message("The protein ", counter, " has no sequence.")
    return(NA_character_)
  }
  
  paste0(
    sprintf(">%s\n", attr(i$pfam$entry, "accession")),
    paste(d, collapse="\n"),
    "*")
})

# Removing NAs
fasta <- fasta[!is.na(fasta)]

pfamscan_results <- pfamscan(
  fasta_str     = fasta,
  email         = "vegayon@usc.edu",
  maxchecktime  = 180
)

# Checking which ones were not present
which(!names(fasta) %in% pfamscan_results$seq$name)
# 22 out of 1337 in total
# [1]   53  566  605  606  608  924  925  926  927  928 1066 1067 1070 1071 1072 1073 1074 1075 1076 1077 1078 1079

saveRDS(pfamscan_results, "data/aphylo_families_pfamscan.rds")

# ------------------------------------------------------------------------------

# We will create files in aphylo_families, and the files will be
# indexed by the panther families.
families <- dat[, .N, by = tree][order(N),] # PTHR10438

for (f in as.character(families$tree)) {
  proteins <- dat[tree == f, as.character(UniProtKB)]
  
  # Filtering the results
  proteins <- pfamscan_results[pfamscan_results$seq$name %in% proteins,]
  as_tblout(
    proteins,
    output = sprintf("data/aphylo_families/pfam_%s.txt", f)
    )
}
