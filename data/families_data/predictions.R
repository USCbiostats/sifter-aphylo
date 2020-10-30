#' Runs SIFTER
sifter_java <- function(
  fam,
  fam_data  = "SIFTER-master/large_scale_v1.0/data/families_data",
  fn_tree   = sprintf("%s/reconciled_trees/%s_reconciled.xml", fam_data, fam),
  fn_ann    = sprintf("%s/annotations/%s.pli", fam_data, fam),
  fn_sifter = "/home/vegayon@PREVMED.USC.EDU/sifter-aphylo/SIFTER-master/core/sifter2.1.1.jar",
  fn_godb   = "/home/vegayon@PREVMED.USC.EDU/sifter-aphylo/SIFTER-master/large_scale_v1.0/data/goterms.sqlite",
  sifter_args = "--truncation 1 --xvalidation --folds 0",
  out_dir   = "",
  keepfiles = TRUE,
  limit_out = 20,
  timeout   = 60*60*2
) {
  
  message(paste(rep("-", getOption("width", 80)), collapse = ""))
  message("Starting to work on Family ", fam, appendLF = FALSE)
  
  if (out_dir == "")
    tmp_d <- tempfile(paste0("sifter-", fam))
  else
    tmp_d <- out_dir
  
  dir.create(tmp_d, recursive = TRUE)
  
  message(" in directory ", tmp_d)
  
  if (!keepfiles)
    on.exit(unlink(tmp_d, recursive = TRUE, force = TRUE))
  
  fn_res <- file.path(tmp_d, "results.txt")

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
  
  nfuns    <- readLines(fn_ann)
  nfuns    <- nfuns[grepl("<GONum", nfuns)]
  nfuns    <- gsub("</?[a-zA-Z]+>|\\s+|\\[|\\]", "", nfuns)
  nfuns    <- unique(unlist(strsplit(nfuns, ",")))
  funnames <- as.integer(nfuns)
  nfuns    <- length(nfuns)
  
  # Creating PxP diagonal and 1 off ---------------------------
  write_aux_files <- function(
    nfuns,
    funnames,
    fn_familyfile,
    fn_scale,
    fn_alpha
  ) {
    
    mat <- matrix("1.0", nrow = nfuns, nfuns, dimnames = list(funnames, funnames))
    diag(mat) <- "0.5"
    colnames(mat)[1] <- paste0("\t",colnames(mat)[1])
    
    
    cat("# scale parameter: 20.0\n", file = fn_familyfile)
    write.table(
      mat, sep="\t", quote=FALSE, col.names = TRUE, file = fn_familyfile,
      append = TRUE)
    
    # Creatint the paramters files
    cat("species\t0.03\nduplication\t0.05\n", file = fn_scale)
    cat(rep("1.0", nfuns), sep = "\n", file = fn_alpha)
    
  }
  
  # Building filepaths
  fn_familyfile <- file.path(tmp_d, paste0("infer-", fam, ".fx"))
  fn_scale      <- file.path(tmp_d, sprintf("scale-%s.fx", fam))
  fn_alpha      <- file.path(tmp_d, sprintf("alpha-%s.fx", fam))
  
  write_aux_files(
    nfuns = nfuns, funnames = funnames, fn_familyfile = fn_familyfile,
    fn_scale = fn_scale, fn_alpha = fn_alpha
  )
  
  # Trying out sifter for a single family
  cmd_sprintf <- '-jar -Xmx4g %s
    --with-exp --with-ida --with-ipi --with-imp --with-igi --with-iep --with-tas --with-nas
    --protein %s
    --reconciled %s
    --ontology %s
    --output %s
    --familyfile %s
    --scale %s
    --alpha %s
    %s %s'
  
  cmds <- c(
    "java", sprintf(
      cmd_sprintf, # --xvalidation --folds 0 
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
      fn_familyfile,
      # scale
      fn_scale,
      fn_alpha,
      sifter_args,
      fam
    ))
  
  cmds[2] <- gsub("\n", " ", cmds[2])
  
  while (TRUE) {
    
    # Making the call
    time0 <- proc.time()    
    out <- suppressWarnings(system2(
      cmds[1], cmds[2], stdout = TRUE, stderr = TRUE,
      wait = TRUE, timeout  = timeout
    ))
    
    # Wrapping answer
    answer <- list(
      call   = cmds,
      out    = out,
      status = attr(out, "status"),
      time0  = time0,
      time1  = proc.time(),
      tmp_d  = tmp_d
    )
    
    # If status is NULL, then it could be that it finished!
    if (is.null(answer$status)) {
      
      if (file.exists(fn_res)) {
        message("Predictions run suscessfully! Retrieving the results")
        answer$predictions <- data.table::fread(fn_res, fill = TRUE)
        answer$status      <- 0L
        return(answer)
      } else {
        message("It seems that the algorithm finished, but no result file is available.")
        answer$status <- 1L
      }
      
    }
    
    # Trying to find error
    if (!answer$status %in% c(0, 124)) {
      
      message("An error ocurred, we will try to rerun setting a new config, here is the output:")
      cat(tail(answer$out, limit_out), sep = "\n")
      
      msg_line <- which(grepl("^Transition matrix functions do not match", answer$out))
      if (!length(msg_line)) {
        message("We failed to identify the error, returning.")
        return(answer)
      }
      
      # Looking for the Pruned GO DAT list
      funnames <- answer$out[which(grepl("^Pruned GO DAG", answer$out))]
      
      funnames <- strsplit(gsub(".+\\[|\\].*$", "", funnames), ",")[[1L]]
      funnames <- unlist(funnames)
      nfuns    <- length(funnames)
      
      write_aux_files(
        nfuns = nfuns, funnames = funnames, fn_familyfile = fn_familyfile,
        fn_scale = fn_scale, fn_alpha = fn_alpha
      )
      
      next
      
    } else if (answer$status == 124) {
      message("We ran out of time, here is the output:")
      cat(tail(answer$out, limit_out), sep = "\n")
    }
    
    break
  }
  
  return(answer)
  
}

# Running sifter ---------------------------------------------------------------

families <- list.files("data/families_data/annotations/", pattern = "*pli")
families <- gsub("\\..+", "", families)

# Cross validation predictions
parallel::mclapply(families, function(f) {
  
  fn_out <- sprintf("data/families_data/prediction_bis/%s-xval.rds", f)
  if (file.exists(fn_out)) {
    message("Family ", f, " already processed.")
    return(NULL)
  }
  
  tmp_ans <- sifter_java(
    fam         = f,
    fam_data    = "data/families_data/",
    fn_tree     = sprintf("SIFTER-master/large_scale_v1.0/data/families_data/reconciled_trees/%s_reconciled.xml", f),
    fn_ann      = sprintf("data/families_data/annotations/%s.pli", f),
    sifter_args = "--truncation 3 --xvalidation --folds 0",
    out_dir     = sprintf("data/families_data/prediction_bis/%s-xval", f) 
  )
  
  # Saving
  if (tmp_ans$status == 0)
    saveRDS(tmp_ans, fn_out)
  
  return(tmp_ans$out)
  
}, mc.cores = 10L)

