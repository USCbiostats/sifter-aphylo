# Listing missing families
x <- readLines("data/families_data/predictions.Rout")
x <- x[grepl("simpleError.+No such file", x)]
x <- unique(gsub(".+_trees/|_reconciled\\.xml\\.gz>$", "",x))

# The following families were not found in SIFTER's PFAM
cat(x, sep = "\n")
# PF16178
# PF16486
# PF16487
# PF16488
# PF17205
# PF18372

# Downloading directly from pfam
sapply(x, function(f) {
  
})