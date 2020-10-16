#' Read New Hampshire eXtended format for trees
#' @param fn Full path to the tree file
#' @return A list with the following elements:
#' - tree An object of class `ape`
#' - edge Edge annotations (length and other annotations)
#' - nhx A list of annotations NHX
read_nhx <- function(fn) {
  text <- paste(readLines(fn), collapse = "")
  
  # This pattern catches most of the information
  pattern <- "([(,)])([a-zA-Z0-9_]*)([:][0-9]+[.]?[0-9]*)(\\[&&NHX[a-zA-Z0-9_:=]+\\])?"
  
  # Capturing the patterns and splitting the data
  x <- gregexpr(text, pattern = pattern, perl = TRUE)
  x <- regmatches(text, x)
  
  # Creating a matrix with the data
  x <- regmatches(x[[1]], regexec(x[[1]], pattern = pattern, perl = TRUE))
  x <- do.call(rbind, x)
  
  # Do all have ids?
  noid <- which(x[,3] == "")
  if (length(noid)) {
    x[noid,3] <- sprintf("unnamed%04i", 1:length(noid))
  }
  
  for (i in noid) {
    text <- sub(
      pattern     = x[i,1],
      replacement = paste0(x[i,2:4], collapse = ""),
      x           = text,
      fixed       = TRUE
    )
  }
  
  # Is there any root?
  text <- sub(
    pattern = "[)][;]$", replacement = ")root;", x = text, perl = TRUE
  )
  
  dat <- x[,-c(1L, 2L)]
  
  # Capturing NHX fields
  nhx <- strsplit(dat[,3], split = "[:]|[=]")
  nhx <- tryCatch(lapply(nhx, function(n) {
    if (length(n) > 0) {
      n <- gsub(pattern = "(^\\[|\\]$)", replacement = "", x = n)
      n <- matrix(n[-1], ncol = 2, byrow = TRUE)
      structure(.Data = n[,2], names = n[,1])
    } else
      n
  }), error = function(e) e)
  
  if (inherits(nhx, "error")) 
    stop(
      "There was a problem when processing the &&NHS blocks.",
      " Possibly, not all the attributes have the right tag. Here is the error",
      ":\n", paste0(nhx, collapse=""), call. = FALSE
    )
  
  list(
    tree = ape::read.tree(text = text),
    edge = dat[,-3L],
    nhx  = nhx
  )
  
}

#' Read PLI files from sifter
#' @param fn Full path to the file
#' @return A data table object including the following columns:
#' - name: Used to match UniProtKB data and GOA,
#' - number,
#' - go: A list of the GO annotations
#' - moc: Evidence code
#' - fam: Name of the family
read_pli <- function(fn, dropNAs = TRUE) {
  
  if (!("data.table" %in% loadedNamespaces()))
    stop("This function needs data.table to be loaded.")
  
  ans <- xml2::as_list(xml2::read_html(fn))
  ans <- ans$html$body$family
  ans <- ans[which(names(ans) == "protein")]
  res <- lapply(ans, function(b) {
    
    # Extracting name
    name   <- b$proteinname
    number <- b$proteinnumber
    
    # And annotations
    go     <- unlist(b$gonumber)
    moc    <- unlist(b$moc)
    
    # We are counting at least one annotation
    nann <- max(1, length(go))
    
    data.table::data.table(
      name   = rep(unname(name), nann),
      number = rep(unname(number), nann),
      go     = if (is.null(go)) NA else unname(as.vector(go)),
      moc    = if (is.null(moc)) NA else unname(as.vector(moc))
    )
    
  })

  res <- data.table::rbindlist(res, fill = FALSE)
  res[, c("name", "number") := list(unlist(name), unlist(number))]
  res[, go := strsplit(gsub("\\[|\\]|\\s", "", go), split = ",")]
  res[, moc := strsplit(gsub("\\[|\\]|\\s", "", moc), split = ",")]
  
  if (dropNAs)
    res <- res[!is.na(go)]
  
  go  <- res[["go"]]
  moc <- res[["moc"]]
  
  nrep <- sapply(go, length)
  nrep[nrep == 0] <- 1L
  
  res[, c("go", "moc") := NULL]
  res <- res[rep(1:.N, nrep)]
  
  cbind(
    res,
    data.table(
      go = unlist(go, recursive = TRUE),
      moc = unlist(moc, recursive = TRUE)
      )
    )
  
  
  
}

get_species <- function(i, offspring, species, ans = NULL) {
  
  # Leaf node?
  if (length(offspring[[i]]) == 0) {
    
    ans <- c(ans, species[i])
    
  } else {
    
    for (j in offspring[[i]])
      ans <- get_species(j, offspring, species, ans)
    
  }
  
  return(ans)
  
}

imput_species <- function(tree, i = 1L, env = NULL) {
  
  if (is.null(env)) {
    if (!inherits(tree, "phylo"))
      stop("-tree- must be a phylo class object.", call. = FALSE)
    
    env <- list()
    env$species     <- gsub(pattern = ".+[_]", replacement = "", tree$tip.label)
    env$offspring   <- aphylo::list_offspring(tree)
    env$checked     <- rep(FALSE, ape::Nnode(tree, internal.only = FALSE))
    env$off_species <- vector("list", length(env$checked))
    env$ntip        <- ape::Ntip(tree)
    env$ntot        <- length(env$checked)
    env$par         <- rbind(tree$edge, c(NA, env$ntip + 1))
    env$par         <- env$par[order(env$par[,2]), 1]
    
    env <- list2env(env)
  }
  
  if (env$checked[i]) {
    return(env)
  }
  env$checked[i] <- TRUE
  
  # Going down
  if (!length(env$offspring[[i]])) { # Is a leaf
    env$off_species[[i]] <- env$species[i]
  } else { # Has offspring
    sapply(env$offspring[[i]], imput_species, tree = tree, env = env)
    # Collecting species
    env$off_species[[i]] <- unique(unlist(env$off_species[env$offspring[[i]]]))
    
  }
  
  # Root node
  if (i == (env$ntip + 1))
    return(as.list(env))
  
  if (is.na(env$par[i])) {
    
    warning("This is fishy")
    
  }
  
  imput_species(tree, env$par[i], env)
  
}

imputate_duplications2 <- function(tree, species = NULL) {
  
  # Checking the input
  if (!inherits(tree, "phylo"))
    stop("This -tree- should be of class 'phylo'", call. = FALSE)
  
  # Getting the species
  if (!length(tree$tip.label))
    stop("This tree has no labels on its tips.", call. = FALSE)
  
  if (!length(species))
    species <- gsub(pattern = ".+[_]", replacement = "", tree$tip.label)
  
  # Retrieving the species
  DAT <- imput_species(tree)
  
  # Vector of answers
  dpl   <- rep(FALSE, length(DAT$par))
  
  for (p in (DAT$ntip + 1):DAT$ntot) {
    
    # Getting the species in p
    species_p <- DAT$off_species[DAT$offspring[[p]]]
    
    # Comparing lists
    for (i in combn(1:length(species_p), 2, simplify = FALSE)) {
      if (any(species_p[[i[1]]] %in% species_p[[i[2]]])) {
        dpl[p] <- TRUE
        break
      }
    }
    
  }
  
  return(dpl)
}

imputate_duplications <- function(tree, species = NULL) {
  
  # Checking the input
  if (!inherits(tree, "phylo"))
    stop("This -tree- should be of class 'phylo'", call. = FALSE)
  
  # Getting the species
  if (!length(tree$tip.label))
    stop("This tree has no labels on its tips.", call. = FALSE)
  
  if (!length(species))
    species <- gsub(pattern = ".+[_]", replacement = "", tree$tip.label)
  
  # Listing offspring
  off <- aphylo::list_offspring(tree)
  size_ <- ape::Nnode(tree, internal.only = FALSE)
  dpl <- logical(size_)
  for (p in (ape::Ntip(tree) + 1):size_) {
    
    # Collecting the species
    species_p <- vector("list", length(off[[p]]))
    
    for (i in seq_along(off[[p]])) {
      species_p[[i]] <- get_species(
        i         = off[[p]][i],
        offspring = aphylo::list_offspring(tree),
        species   = species
      )
    }
    
    # Comparing lists
    for (i in combn(1:length(species_p), 2, simplify = FALSE)) {
      if (any(species_p[[i[1]]] %in% species_p[[i[2]]])) {
        dpl[p] <- TRUE   
        break
      }
      
    }
  }
  
  return(dpl)
  
  
}


if (FALSE) {

  ans <- imput_species(tree2)
  imputate_duplications2(tree2)
  
  microbenchmark::microbenchmark(
    dpl1 <- imputate_duplications(tree2),
    dpl2 <- imputate_duplications2(tree2), times =10, unit = "relative"
  )
  
  # How long to process a large tree
  tree <- read_nhx("sifter/hundred/PF00098.nhx")
  ans <- imputate_duplications2(tree$tree)
}