library(ape)
library(xml2)
library(data.table)

# Functions to read sifter trees
source("sifter/read_sifter.R")

# Reading in the annotations
fn <- list.files("sifter/hundred", full.names = TRUE, pattern = "pli$")
dat <- parallel::mclapply(fn, read_pli, mc.cores = 4L)

# Adding the corresponding family
nams <- gsub(pattern = ".+[/]|[.]pli$", replacement = "", x = fn)
dat <- cbind(
  do.call(rbind, dat),
  famid = rep(nams, sapply(dat, nrow))
)

# Saving the data
fwrite(dat, "sifter/hundredfamilies_annotations.csv")
