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

# Getting the accession and the pfam family

tmp <- lapply(seq_along(pfam_ids), function(i) {
  # print(i)
  data.table(
    name   = attr(pfam_ids[[i]]$pfam$entry, "id"),
    number = names(pfam_ids)[i],
    fami   = unlist(sapply(pfam_ids[[i]]$pfam$entry$matches, attr, which = "accession"))
  )
})

tmp <- rbindlist(tmp, fill = TRUE)
tmp <- tmp[!is.na(fami)]

pli_data <- merge(
  x = dat,
  y = tmp, by.x = "UniProtKB", by.y = "number",
  all.x = TRUE, all.y = FALSE, allow.cartesian = TRUE
)

# Only those that matched a family
pli_data <- pli_data[!is.na(fami)]
saveRDS(pli_data, "data/families_data/annotations.rds")

# Filtering the annotations:
# - Have at least two functions per tree (in PFAM, we do that later), and
# - Are all positive
pli_data <- pli_data[qualifier == 1]

unique(pli_data[,.(go,fami)])[, .(n=.N), by = fami][n>1]

# For the family name, I will only keep the Panther name so that way
# I make sure that 

# Families
# pfam_families <- readLines("data/families_data/annotations_families.txt")
# pfam_families <- pfam_families[!grepl("^#", pfam_families)]
# pli_data <- pli_data[fami %in% pfam_families]

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

families_worth <- unique(pli_data[,.(go,fami)])[, .(n=.N), by = fami][n>1, fami]
# 60 families

for (f in families_worth) {
  
  # Creating the entire set
  write_pli(
    family_id      = f,
    protein_name   = as.character(pli_data[fami == f, name]),
    protein_number = as.character(pli_data[fami == f, UniProtKB]),
    go_number      = gsub("GO:","",as.character(pli_data[fami == f, go])),
    moc            = "EXP",
    file = sprintf("data/families_data/annotations/%s.pli", f)
  )
  
  message("Family ", f, " done.")
  
}

