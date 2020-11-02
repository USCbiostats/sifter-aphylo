library(aphylo)
library(data.table)
library(pfamscanr)

# Loading the trees for the experiment
candidate_trees <- readRDS("data-raw/candidate_trees_3.rds")

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

save.image("data/aphylo_families_entries_3.rda")

# Setting up the data to be used with SIFTER -----------------------------------

ans <- unique(dat[, .(tree, UniProtKB)])[, .(n=.N), by = tree]
ans[order(n),] # Selecting one tree with small number of proteins annotated

# Creating the fasta sequences for querying pfamscan ---------------------------
counter <- 0L
fasta <- lapply(pfam_ids, function(i) {
  
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

# Gathering taxas
pfam_taxa <- lapply(pfam_ids, function(i) {
  a <- list(
    tax_id = attr(i$pfam$entry$taxonomy, "tax_id"),
    id     = attr(i$pfam$entry, "id")
  )
  if (!length(a$tax_id))
    a$tax_id <- NA_character_
  if (!length(a$id))
    a$id <- NA_character_
  a
})

pfam_taxa <- data.table(
  uniprot = names(pfam_taxa),
  taxa    = sapply(pfam_taxa, "[[", "tax_id"),
  id      = sapply(pfam_taxa, "[[", "id")
)
pfam_taxa <- pfam_taxa[, .(taxa = unique(taxa)), by = .(uniprot, id)]
pfam_taxa[, n := sapply(taxa, length)]
pfam_taxa[n == 0, taxa := NA_character_]
pfam_taxa[, n := NULL]

# When unlisting, this should work OK
stopifnot(length(unlist(pfam_taxa$taxa)) == nrow(pfam_taxa))
pfam_taxa[, taxa := unlist(taxa)]

pfam_taxa <- merge(
  x = dat,
  y = pfam_taxa,
  by.x = "UniProtKB", by.y = "uniprot",
  all.x = TRUE, all.y = FALSE
)

# Making sure everything merged properly
stopifnot(nrow(pfam_taxa) == nrow(dat))
dat <- copy(pfam_taxa)

pfam_taxa <- pfam_taxa[,.(tree, taxa)]
pfam_taxa <- unique(pfam_taxa)
pfam_taxa[is.na(taxa), ]
#         tree taxa
# 1: PTHR11256 <NA>
# 2: PTHR10125 <NA>

pfam_taxa <- pfam_taxa[!is.na(taxa)] # 497 queries

# We will create files in aphylo_families, and the files will be
# indexed by the panther families.
families <- split(pfam_taxa, 1:nrow(pfam_taxa))

for (f in families) {
  proteins <- dat[tree == f$tree & taxa == f$taxa, as.character(UniProtKB)]
  ids      <- dat[tree == f$tree & taxa == f$taxa, as.character(id)]
  
  # Writing the list of HMM scans
  tmp <- pfamscan_results[pfamscan_results$seq$name %in% proteins,]
  if (nrow(tmp) == 0) {
    message(
      sprintf(
        "The family %s (taxid %s) has no rows in pfamscan.",
        f$tree, f$taxa
    ))
    next
  }
  
  # Writing the sequences
  fasta_f <- fasta[names(fasta) %in% proteins]
  fasta_f <- paste(fasta_f, collapse = "\n")
  cat(fasta_f, file = sprintf("data/aphylo_families/%s_taxa-%s.fasta", f$tree, f$taxa))
  

  as_tblout(
    tmp,
    output = sprintf("data/aphylo_families/%s_taxa-%s.txt", f$tree, f$taxa)
    )
  
  # Writing the protein list
  cat(
    ids, sep = "\n",
    file = sprintf("data/aphylo_families/%s_taxa-%s_proteins.txt", f$tree, f$taxa)
    )
}

# The family PTHR12663 (taxid 3702) has no rows in pfamscan.
# The family PTHR11132 (taxid 559292) has no rows in pfamscan.
# The family PTHR12663 (taxid 559292) has no rows in pfamscan.
# The family PTHR11132 (taxid 237561) has no rows in pfamscan.
# The family PTHR23235 (taxid 7955) has no rows in pfamscan.
# The family PTHR10730 (taxid 10090) has no rows in pfamscan.
# The family PTHR11132 (taxid 284812) has no rows in pfamscan.

# Building a super family ------------------------------------------------------
superfam <- dat[!is.na(id)]
proteins <- superfam[, unique(UniProtKB)]
ids      <- superfam[, unique(id)]

# Writing the list of HMM scans
tmp <- pfamscan_results[pfamscan_results$seq$name %in% proteins,]

# Writing the sequences
fasta_f <- fasta[names(fasta) %in% proteins]
fasta_f <- paste(fasta_f, collapse = "\n")
cat(fasta_f, file = "data/aphylo_families/all_families.fasta")

as_tblout(
  tmp,
  output = "data/aphylo_families/all_families.txt"
)

# Writing the protein list
cat(
  ids, sep = "\n",
  file = "data/aphylo_families/all_families_proteins.txt"
)


