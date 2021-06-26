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
  all.x = TRUE,
  all.y = FALSE,
  # Extending, allowing proteins to show up in multiple families
  allow.cartesian = TRUE
)

# Only those that matched a family
pli_data <- pli_data[!is.na(fami)]
saveRDS(pli_data, "predictions/annotations.rds")

pli_data[qualifier == 1][,.(n = length(unique(go))),by=fami][n > 1]

# Filtering the annotations:
# - Have at least two functions per tree (in PFAM, we do that later), and
# - Are all positive
pli_data <- pli_data[qualifier == 1]

unique(pli_data[,.(go,fami)])[, .(n=.N), by = fami][n>1]

families_worth <- unique(pli_data[,.(go,fami)])[, .(n=.N), by = fami][n>1, fami]
# 60 families

# How many genes
unique(pli_data[fami %in% families_worth, .(UniProtKB)]) # 472

for (f in families_worth) {
  
  # Creating the entire set
  write_pli(
    family_id      = f,
    protein_name   = as.character(pli_data[fami == f, name]),
    protein_number = as.character(pli_data[fami == f, UniProtKB]),
    go_number      = gsub("GO:","",as.character(pli_data[fami == f, go])),
    moc            = "EXP",
    file = sprintf("predictions/annotations/%s.pli", f)
  )
  
  message("Family ", f, " done.")
  
}

