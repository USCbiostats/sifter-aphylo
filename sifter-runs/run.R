# In order for this script to work, the file sifter_db_pass.txt must exists
# at the top level of this project.

family_map <- NULL

#' This is a wrapper that runs SIFTER
#' @param fam A family in the form "PTHR31682_taxa-3702"
#' @param dbpass Password for the database
#' @param path Path to where all the families (fasta, pfamscan, protein list)
#' @param sifter Path to SIFTER's scripts.
#' @param timeout 
#' @return A list with timings and other relevant data.
run_sifter_on_aphylo <- function(
  fam,
  dbpass,
  path    = "data/aphylo_families/",
  sifter  = "SIFTER-master/large_scale_v1.0/scripts/",
  fam_data = "SIFTER-master/large_scale_v1.0/data/families_data/",
  timeout = 1
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
    famlist <- readLines(sprintf("%s/family_list.txt", tmp_dir))
    message("Done! Families:", paste(famlist, collapse=", "))
    
    family_map <<- c(family_map, structure(paste(famlist, collapse=", "), names = fam))
  }
  
  # Step 3: --------------------------------------------------------------------
  message("sifter_prepare ... ", appendLF = FALSE)
  time_step3 <- proc.time()
  cmd <- sprintf(
    "%s/sifter_prepare.py --dbpass=%s --ip %s/%s_proteins.txt %s %s",
    sifter, dbpass, tmp_dir, fam, fam_data, tmp_dir
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
  
  # Checking the timing
  which_time <- which(grepl("^Total estimated time", out_step3))
  if (!length(which_time)) {
    message("Couldn't measure estimated time. Trying with the next.")
    cat(out_step3, sep = "\n")
    return(NULL)
  }
  
  # checking time
  time_estimate <- gsub(".+your query is\\s+|\\s+\\(.+", "", out_step3[which_time[1]])
  time_unit <- gsub(".+\\s+", "", time_estimate)
  time_estimate <- as.numeric(gsub("\\s+.+", "", time_estimate))
  if (time_unit == "hours") {
    time_estimate <- time_estimate*60
  } else if (time_unit == "minutes") {
    time_estimate <- time_estimate
  } else if (time_unit == "days")
    time_estimate <- time_estimate*60*24
  
  if (time_estimate >= timeout) {
    message("The estimate was ", time_estimate, "minutes which is above the time out.")
    message(out_step3[which_time[1]], sep = "\n")
    return(NULL)
  } else
    message(out_step3[which_time[1]])
  
  
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
# debug(run_sifter_on_aphylo)
SIFTER_PASS <- readLines("sifter_db_pass.txt")

# Listing all the families
families <- list.files(path = "data/aphylo_families/", pattern = "*.fasta")
families <- setdiff(unique(gsub("[.].+", "", families)), "all_families")

# One example The file PTHR10024_taxa-10090 is mapped to PF00168

PF00168 <- xml2::read_xml("data/families_data/annotations/PF00168.pli")
ids_out <- length(xml2::xml_find_all(PF00168, "//GONumber"))
set.seed(12312)
ids_out <- sample.int(ids_out, ids_out - 10)
xml2::xml_remove(xml2::xml_find_all(PF00168, "//GONumber")[ids_out])
xml2::write_xml(PF00168, file = "data/families_data/annotations/PF00168.pli")
R.utils::gzip("data/families_data/annotations/PF00168.pli", overwrite = TRUE)

ans <- run_sifter_on_aphylo("PTHR10024_taxa-10090", SIFTER_PASS, fam_data = "data/families_data")

failed <- NULL
for (f in families[1:10]) {
  
  ans <- run_sifter_on_aphylo(f, SIFTER_PASS, fam_data = "data/families_data")
  
  if (!length(ans)) {
    failed <- c(failed, f)
    next
  }
  
  saveRDS(ans, file = sprintf("sifter-runs/%s.rds", f))
  
}

# saveRDS(family_map, "sifter-runs/run.rds")
