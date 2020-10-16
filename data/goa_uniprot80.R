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

fwrite(goa, "sifter/goa_uniprot80_experimental.gz", compress = "gzip")

# Reading hundred fams
fams <- fread("sifter/hundredfamilies_annotations.csv")

dat <- merge(x = fams, y = goa, by.y = "synm", by.x = "number",
  all.x = TRUE, all.y = FALSE)

saveRDS(dat, "sifter/goa_uniprot80.rds", compress = TRUE)

goterms <- strsplit(fams[go != ""]$go, split = ",")
goterms <- gsub("\\s+|\\[|\\]", "", unlist(goterms))
length(goterms)
nrow(fams[go != ""])

table(goterms)

dat_wide <- dat[GO_ID %in% paste0("GO:", goterms)]
dat_wide <- dcast(dat_wide, number + famid ~ Qualifier + GO_ID)

library(aphylo)

set.seed(3)
tree0 <- raphylo(100, psi = c(0,0), mu_d=c(.95, .5), mu_s = c(.05, .05), Pi = .01)
# plot(tree0)

tree <- rdrop_annotations(tree0, .9, prob.drop.0 = .9)
plot(tree)

pred <- predict_pre_order(tree, psi = c(0,0), mu_d=c(.95, .5), mu_s = c(.05, .05), Pi = .01, eta = c(1, 1))

pred <- pred[1:Ntip(tree0),,drop=FALSE]

auc(pred, tree0$tip.annotation)
tree0$tip.annotation <- cbind(tree0$tip.annotation, tree$tip.annotation, pred = pred)
plot(tree0)
