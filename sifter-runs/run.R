# In order for this script to work, the file sifter_db_pass.txt must exists
# at the top level of this project.

#' This is a wrapper that runs SIFTER
#' @param fam A family in the form "PTHR31682_taxa-3702"
#' @param dbpass Password for the database
#' @param path Path to where all the families (fasta, pfamscan, protein list)
#' @param sifter Path to SIFTER's scripts.
#' @return A list with timings and other relevant data.
run_sifter_on_aphylo <- function(
  fam,
  dbpass,
  path   = "data/aphylo_families/",
  sifter = "SIFTER-master/large_scale_v1.0/scripts/"
  ) {
  
  message(paste(rep("-", getOption("width", 80)), collapse = ""))
  message("Starting to work on Family+Tax ", fam)
  
  # Step 1a: Prepare the data
  tmp_dir <- tempdir(check = TRUE)
  
  # Step 1b: Copy the data
  fpaths <- list.files(path = path, pattern = paste(fam, "*"), full.names = TRUE)
  file.copy(from = fpaths, to = tmp_dir)
  
  # Step 2: Check whether sifter has the data ----------------------------------
  message("sifter_find_families on ", fam, "...", appendLF = FALSE)
  time_step2 <- proc.time()

  cmd <- sprintf(
    "%s/sifter_find_families.py --dbpass=%s --ip %s/%s_proteins.txt %s/family_list.txt",
    sifter, dbpass, tmp_dir, fam, tmp_dir
    )
  
  out_step2 <- tryCatch(
    system2("python", cmd, stdout = TRUE, stderr = TRUE, wait = TRUE), 
    error = function(e) e
    )
  
  if (inherits(out_step2, "error")) {
    print(out_step2)
    return(NULL)
  } else {
    message(
      "Done! Families:",
      paste(
        readLines(sprintf("%s/family_list.txt", tmp_dir)),
        collapse=", "
        )
      )
  }
  
  # Step 3: --------------------------------------------------------------------
  message("sifter_prepare ... ", appendLF = FALSE)
  time_step3 <- proc.time()
  cmd <- sprintf(
    "%s/sifter_prepare.py --dbpass=%s --ip %s/%s_proteins.txt %s/../data/families_data/ %s",
    sifter, dbpass, tmp_dir, fam, sifter, tmp_dir
  )
  
  out_step3 <- tryCatch(
    system2("python", cmd, stdout = TRUE, stderr = TRUE, wait = TRUE), 
    error = function(e) e
  )
  
  if (inherits(out_step3, "error")) {
    print(out_step3)
    return(NULL)
  } else {
    message("Done!")
  }
  
  # Step 4 ---------------------------------------------------------------------
  message("sifter_run_cv ... ", appendLF = FALSE)
  time_step4 <- proc.time()
  cmd <- sprintf("%s/sifter_run_cv.py %s %s/results", sifter, tmp_dir, tmp_dir)
  out_step4 <- tryCatch(
    system2("python", cmd, stdout = TRUE, stderr = TRUE, wait = TRUE), 
    error = function(e) e
  )
  
  if (inherits(out_step4, "error")) {
    print(out_step4)
    return(NULL)
  } else {
    message("Done!") 
  }
  
  # Step 5 ---------------------------------------------------------------------
  message("sifter_run_cv ... ", appendLF = FALSE)
  time_step5 <- proc.time()
  cmd <- sprintf(
    "%s/sifter_extract.py --dbpass=%s --ip %s/%s_proteins.txt %s/results/ %s/preds.txt",
    sifter, dbpass, tmp_dir, fam, tmp_dir, tmp_dir
  )
  out_step5 <- tryCatch(
    system2("python", cmd, stdout = TRUE, stderr = TRUE, wait = TRUE), 
    error = function(e) e
  )
  
  if (inherits(out_step5, "error")) {
    print(out_step5)
    return(NULL)
  } else {
    message("Done!")
  }
  
  predictions <- readLines(sprintf("%s/preds.txt", tmp_dir))
  
  # CLosing up
  unlink(tmp_dir, TRUE, TRUE)
  
  # Returning ------------------------------------------------------------------
  ans <- list(
    timings = cbind(time_step2, time_step3, time_step4, time_step5, proc.time()),
    output_cv = out_step4,
    predictions = predictions
  )
  
  ans
  
  
}

SIFTER_PASS <- readLines("sifter_db_pass.txt")

# Listing all the families
families <- list.files(path = "data/aphylo_families/", pattern = "*.fasta")
families <- unique(gsub("[.].+", "", families))

for (f in families) {
  
  ans <- run_sifter_on_aphylo("PTHR31682_taxa-3702", SIFTER_PASS)
  saveRDS(ans, file = sprintf("sifter-runs/%s.rds", f))
  
}

  