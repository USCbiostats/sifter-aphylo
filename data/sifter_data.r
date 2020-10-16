library(aphylo)
library(data.table)
library(pfamscanr)

# Loading the trees for the experiment
candidate_trees <- readRDS("data/candidate_trees.rds")

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

# Obtaining the Pfam name ------------------------------------------------------
# https://www.uniprot.org/help/api_retrieve_entries
get_info_uniprot <- function(uniprotkb, skip_failed = TRUE, db = "pfam") {
  
  # Can be applied to multiple data
  if (length(uniprotkb) > 1L) {
    
    message("Starting multiple queries...")
    
    ans <- vector("list", length(uniprotkb))
    for (i in seq_along(uniprotkb)) {
      
      # Making the query
      ans[[i]] <- get_info_uniprot(
        uniprotkb = uniprotkb[i], skip_failed = skip_failed, db = db
        )
      
      if (!(i %% 10))
        message(sprintf("% 4i/% 4i complete...", i, length(uniprotkb)))
    }
    message("done.")
    names(ans) <- uniprotkb
    return(ans)
  }
  
  # Selecting the source
  query <- if (db == "pfam")
    "http://pfam.xfam.org/protein/%s?output=xml"
  else if (db == "uniprot")
    "https://www.uniprot.org/uniprot/%s.txt"
  
  res <- httr::GET(
    sprintf(query, uniprotkb)
    )
  
  # Checking status
  if (httr::status_code(res) != 200) {
    if (skip_failed) {
      warning("The query failed with status: ", httr::status_code(res))
      return(NULL)
    } else {
      stop("The query failed with status: ", httr::status_code(res))
    }
    
  }
  
  # Returning raw content
  ans <- xml2::as_list(httr::content(res))
  
  if (db == "pfam") {
    return(ans)
  }
  
  return(ans)
  
}

pfam_ids <- unique(dat[,as.character(UniProtKB)])
pfam_ids <- get_info_uniprot(pfam_ids)

# The following IDS were not found in PFAM
# names(pfam_ids)[c(259, 642)]
# [1] "F1QJU1"     "A0A0R4ITP9"

save.image("data/sifter_data.rda")

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

saveRDS(pfamscan_results, "data/sifter_data_pfamscan.rds")
