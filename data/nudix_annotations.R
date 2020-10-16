library(data.table)

# Reading annotations
goa <- fread("sifter/gene_association.goa_uniprot.80", header = FALSE)

# Setting colnames
cnames <- readLines("sifter/gaf_1.0_colnames.txt")
colnames(goa) <- cnames

# Keeping only relevant
codes <- c("IDA", "IMP", "TAS")
goa <- goa[Evidence_Code %in% codes]

# ID Synonim
goa[, synm := gsub(".+[|]([a-zA-Z0-9_]+)$", "\\1", DB_Object_Synonym, perl = TRUE)]

