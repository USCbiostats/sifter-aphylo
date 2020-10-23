library(aphylo)
library(data.table)
library(pfamscanr)
library(rphyloxml)

# Loading the trees for the experiment
candidate_trees <- readRDS("data-raw/candidate_trees.rds")
load("data/aphylo_families_entries.rda")

# What is the distribution of functions per tree
functions_per_tree <- unique(dat[, .(go,tree)])[, .(n=.N), by = tree] 
functions_per_tree[, table(n)]

trees_with_more_than_2 <- functions_per_tree[n >= 2,]

# Building the pli files for sifter --------------------------------------------

# Filtering the annotations:
# - Have at least two functions per tree, and
# - Are all positive
pli_data <- dat[tree %in% trees_with_more_than_2$tree & qualifier == 1]

# Getting the accession and the pfam family

tmp <- lapply(seq_along(pfam_ids), function(i) {
  data.table(
    name = names(pfam_ids)[i],
    fami = attr(pfam_ids[[i]]$pfam$entry$matches[[1]], "accession")
  )
})

tmp <- rbindlist(tmp, fill = TRUE)
tmp <- tmp[!is.na(fami)]

pli_data <- merge(
  x = pli_data,
  y = tmp, by.x = "UniProtKB", by.y = "name",
  all.x = TRUE, all.y = FALSE
)
pli_data <- pli_data[!is.na(fami)]

# For the family name, I will only keep the Panther name so that way
# I make sure that 

# Families
pfam_families <- readLines("data/families_data/annotations_families.txt")
pfam_families <- pfam_families[!grepl("^#", pfam_families)]
pli_data <- pli_data[fami %in% pfam_families]

# # Copying the families' annotations from SIFTER's data
# all_pfam <- list.files(
#   "SIFTER-master/large_scale_v1.0/data/families_data/annotations/",
#   full.names = TRUE,
#   pattern    = "*.pli.gz"
#   )
# 
# # Selecting files to copy
# all_pfam <- data.table(path = all_pfam)
# all_pfam[, file := gsub(".+//(?=PF[0-9]+\\.)", "", path, perl = TRUE)]
# all_pfam[, file := gsub("\\..+", "", file)]
# all_pfam <- all_pfam[file %in% pfam_families]
# 
# file.copy(
#   all_pfam$path, to = "data/families_data/annotations/",
#   overwrite = TRUE,
#   copy.date = TRUE
#   )
# 
# PF00001 <- xml2::read_xml("data/families_data/annotations/PF00001.pli")

for (f in unique(pli_data$fami)) {
  
  # Creating the entire set
  write_pli(
    family_id      = f,
    protein_name   = as.character(pli_data[fami == f, UniProtKB]),
    protein_number = as.character(pli_data[fami == f, UniProtKB]),
    go_number      = gsub("GO:","",as.character(pli_data[fami == f, go])),
    moc            = "EXP",
    file = sprintf("data/families_data/annotations/%s.pli", f)
  )
  
  message("Family ", f, " done.")
  
}

