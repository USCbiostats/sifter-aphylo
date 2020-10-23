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

sifter_java <- function(
  fam,
  fam_data = "SIFTER-master/large_scale_v1.0/data/families_data",
  fn_tree  = sprintf("%s/reconciled_trees/%s_reconciled.xml", fam_data, fam),
  fn_ann   = sprintf("%s/annotations/%s.pli", fam_data, fam),
  fn_sifter = "/home/vegayon@PREVMED.USC.EDU/sifter-aphylo/SIFTER-master/core/sifter2.1.1.jar",
  fn_godb   = "/home/vegayon@PREVMED.USC.EDU/sifter-aphylo/SIFTER-master/large_scale_v1.0/data/goterms.sqlite"
) {
  
  tmp_d <- tempdir(check = TRUE)
  # on.exit(unlink(tmp_d, recursive = TRUE, force = TRUE))
  
  fn_res <- file.path(tmp_d, "results.txt")
  dir.create(file.path(tmp_d, "results"))
  
  # Moving the files -------------------------------------------------------------
  if (!file.exists(fn_tree)) {
    R.utils::gunzip(
      gsub("\\.xml", ".xml.gz", fn_tree),
      destname = sprintf("%s/%s.xml", tmp_d, fam),
      remove = FALSE, overwrite = TRUE
    )
  } else {
    file.copy(fn_tree, sprintf("%s/%s.xml", tmp_d, fam))
  }
  fn_tree <- sprintf("%s/%s.xml", tmp_d, fam)
  
  if (!file.exists(fn_ann)) {
    R.utils::gunzip(
      gsub("\\.pli", ".pli.gz", fn_ann),
      destname = sprintf("%s/%s.pli", tmp_d, fam),
      remove = FALSE, overwrite = TRUE
    )
  } else {
    file.copy(fn_ann, sprintf("%s/%s.pli", tmp_d, fam))
  }
  fn_ann <- sprintf("%s/%s.pli", tmp_d, fam)
  
  P <- readLines(fn_ann)
  P <- P[grepl("<GONum", P)]
  P <- gsub("</?[a-zA-Z]+>|\\s+|\\[|\\]","",P)
  P <- unique(unlist(strsplit(P, ",")))
  tmp_fun <- as.integer(P)
  P <- length(tmp_fun)
  
  # Creating PxP diagonal and 1 off ---------------------------
  mat <- matrix("1.0", nrow = P, P, dimnames = list(tmp_fun, tmp_fun))
  diag(mat) <- "0.5"
  colnames(mat)[1] <- paste0("\t",colnames(mat)[1])
  
  tmp_infer_fx <- file.path(tmp_d, "results", paste0("infer-", fam, ".fx"))
  cat("# scale parameter: 20.0\n", file = tmp_infer_fx)
  write.table(
    mat, sep="\t", quote=FALSE, col.names = TRUE, file = tmp_infer_fx,
    append = TRUE)
  
  # Creatint the paramters files
  fn_scale <- file.path(tmp_d, "results", sprintf("scale-%s.fx", fam))
  cat("species\t0.03\nduplication\t0.05\n", file = fn_scale)
  fn_alpha <- file.path(tmp_d, "results", sprintf("alpha-%s.fx", fam))
  cat(rep("1.0", P), sep = "\n", file = fn_alpha)
  
  # Trying out sifter for a single family
  cmds <- c(
    "java", sprintf(
      '-jar -Xmx4g %s -v
    --with-exp --with-ida --with-ipi --with-imp --with-igi --with-iep --with-tas --with-nas
    --protein %s
    --reconciled %s
    --ontology %s
    --output %s
    --familyfile %s
    --scale %s
    --alpha %s
    --truncation 1 %s', # --xvalidation --folds 0 
      # SIFTER Path
      fn_sifter,
      # Protein
      fn_ann,
      # Tree
      fn_tree,
      # Ontology
      fn_godb,
      # Results
      fn_res,
      # Family file
      tmp_infer_fx,
      # scale
      fn_scale,
      fn_alpha,
      fam
    ))
  cmds[2] <- gsub("\n", " ", cmds[2])
  system2(cmds[1], cmds[2])
  
}

ans <- sifter_java(
  fam = "PF12474",
  fam_data = "data/families_data/",
  fn_tree = "SIFTER-master/large_scale_v1.0/data/families_data/reconciled_trees/PF12474_reconciled.xml",
  fn_ann = "data/families_data/annotations/PF12474.pli",
)
