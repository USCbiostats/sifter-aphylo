library(data.table)
library(aphylo)

analysis <- c(
  trunc1 = "data/families_data/predictions/",
  trunc3 = "data/families_data/predictions_bis/"
  )

analysis_data <- structure(vector(mode = "list", 2), names = names(analysis))
for (a in seq_along(analysis)) {

  # Reading predictions in
  fn_predictions <- list.files(analysis[a], pattern = "*val\\.rds", full.names = TRUE)
  out            <- lapply(fn_predictions, readRDS)
  names(out)     <- gsub(".+/?PF([0-9]+)(-xval)?\\.rds$", "PF\\1", fn_predictions)
  
  # How much time? 110
  message(names(analysis)[a], " took ", sum(sapply(out, function(o) {
    (o$time1 - o$time0)[3]
  }))/60, " minutes.")
  
  # How many were missed?
  all_fams <- gsub("\\..+", "", list.files("data/families_data/annotations/"))
  all_fams[!all_fams %in% names(out)]
  
  # Looking at the pruned terms --------------------------------------------------
  go_terms_info <- fread("data-raw/go_terms_info.csv")
  sifter_log    <- readLines("data/families_data/predictions.Rout")
  sifter_log    <- sifter_log[grepl("^(Transition file|Pruned GO DAG)", sifter_log)]
  sifter_log    <- matrix(gsub(".+\\[|\\]", "", sifter_log), ncol = 2, byrow = TRUE)
  
  # Listing functions that we didn't find
  sifter_log    <- apply(sifter_log, 1, function(a) {
    ans <- strsplit(a, ",")
    setdiff(ans[[1]], ans[[2]])
  })
  
  # Preparing for checking
  sifter_log <- unique(unlist(sifter_log))
  sifter_log <- sprintf("GO:%07i", as.integer(sifter_log))
  
  go_terms_info[id %in% sifter_log, .(id, name, aspect)]
  
  # Was there any cross validation? ----------------------------------------------
  xval <- lapply(out, function(i) i$out[grepl("^x-val", i$out)])
  xval <- lapply(xval, function(i) {
    i <- gsub("x-val\\s*\\(([[:alnum:]_/-]+)\\)\\s*", "\\1 ", i)
    do.call(rbind, strsplit(i, "\\s+"))
  })
  
  valid_xval <- NULL
  for (i in seq_along(xval)) {
    if (ncol(xval[[i]]) == ncol(out[[i]]$predictions))
      colnames(xval[[i]]) <- colnames(out[[i]]$predictions)
    else {
      # Number of columns for PF00754 doesn't match the read.
      message("Number of columns for ", names(out)[i], " doesn't match the read.")
      next
    }
    
    valid_xval <- c(valid_xval, i)
  }
  
  xval <- xval[valid_xval]
  xval <- lapply(names(xval), function(f) {
    data.table(
      family  = f,
      elapsed = (out[[f]]$time1 - out[[f]]$time0)[3],
      xval[[f]]
      )
  })
  xval <- rbindlist(xval, fill = TRUE)
  xval$`#Names` <- toupper(xval$`#Names`)
  
  # Reshaping
  xval <- melt(
    xval,
    id            = 1:3,
    value.name    = "score",
    variable.name = "term"
  )
  xval <- xval[!is.na(score)]
  xval[, name:=gsub("/.+", "", `#Names`)]
  xval[, score := as.double(score)]
  
  xval[, nfuns_sifter     := length(unique(term)), by = family]
  xval[, nproteins_sifter := length(unique(name)), by = family]
  
  # What is under their metric? --------------------------------------------------
  accuracy_sifter <- gsub(
    ".+([0-9]+) out of ([0-9]+).+", "\\1;\\2",
    sapply(out, function(o) o$out[grepl("^Cross-valid", o$out)])
    )
  accuracy_sifter <- strsplit(accuracy_sifter, ";")
  accuracy_sifter <- do.call(rbind, accuracy_sifter)
  accuracy_sifter <- apply(accuracy_sifter, 2, as.integer)
  accuracy_sifter <- data.table(
    family  = names(out),
    correct = accuracy_sifter[,1],
    total   = accuracy_sifter[,2]
  )
  
  # Overall accuracy
  colSums(as.matrix(accuracy_sifter[,-1]))
  # correct   total ~ 0.43 %
  #     185     426 
  
  # Listing genes that it was not able to find -----------------------------------
  
  # Relevant cases
  missing_proteins <- unique(xval$family)
  
  missing_proteins <- unlist(lapply(out[missing_proteins], "[[", "out"))
  missing_proteins <- missing_proteins[grepl("^In hasNode ", missing_proteins)]
  missing_proteins <- unique(gsub(".+\\sID[:]\\s*", "", missing_proteins))
  missing_proteins <- setdiff(
    missing_proteins, 
    xval[,tolower(unique(name))]
  )
  
  cbind(missing_proteins)
  # missing_proteins  
  # [1,] "gper1_danre"     
  # [2,] "gper1_mouse"     
  # [3,] "lirb3_mouse"     
  # [4,] "hp302_arath"     
  # [5,] "svh2_caeel"      
  # [6,] "ddrb_caeel"      
  # [7,] "max2_caeel"      
  # [8,] "a0a0b4kgs4_drome"
  # [9,] "pk1_caeel"       
  # [10,] "utr7_arath"      
  # [11,] "urgt4_arath"     
  # [12,] "urgt1_arath"     
  # [13,] "gons1_arath"     
  # [14,] "g2ox4_arath"     
  # [15,] "g2ox6_arath"     
  # [16,] "a0a0b4jcy1_drome"
  # [17,] "t2r40_chick"
  
  # Storing everything ---------------------------------------------------------
  xval[, folds := names(analysis)[a]]
  analysis_data[[a]]$xval             <- copy(xval)
  analysis_data[[a]]$missing_proteins <- missing_proteins
  analysis_data[[a]]$accuracy         <- copy(accuracy_sifter)

}

# Checking which families
families_analyzed <- lapply(analysis_data, "[[", "xval")
families_analyzed <- lapply(families_analyzed, "[[", "family")

families_analyzed <- lapply(families_analyzed, unique)
families_analyzed <- lapply(families_analyzed, sort)

# Not in two but in one
with(families_analyzed, trunc1[which(!trunc1 %in% trunc3)])
with(families_analyzed, trunc3[which(!trunc3 %in% trunc1)])

# Retrieving the training data -------------------------------------------------

annotations <- readRDS("data/families_data/annotations.rds")
setnames(annotations, "qualifier", "value")
annotations <- unique(annotations[, fami := NULL]) # Don't need this

# Merging the data and computing AUCs! -----------------------------------------

# Since in some cases we have more than tree per protein (multiple hits in pfam)
# we will keep the average
for (a in names(analysis_data)) {
  
  analysis_data[[a]]$xval <- analysis_data[[a]]$xval[nfuns_sifter > 1 & nproteins_sifter > 1]
  
  analysis_data[[a]]$xval <- merge(
    analysis_data[[a]]$xval, 
    annotations[, .(term = go, name, value)],
    by    = c("term", "name"),
    all.x = TRUE,
    all.y = FALSE
  )
  analysis_data[[a]]$xval <- analysis_data[[a]]$xval[!is.na(value)]
  
  analysis_data[[a]]$xval[, smallest := abs(score - value), by = .(term, name)]
  
  # Avaring out scores
  analysis_data[[a]]$xval[, score := mean(score, na.rm = TRUE), by = .(term, name)]

  analysis_data[[a]]$xval[, smallest := which.min(smallest), by = .(term, name)]
  analysis_data[[a]]$xval[, within_id := 1:.N, by = .(term, name)]
  
  analysis_data[[a]]$xval <- analysis_data[[a]]$xval[smallest == within_id]
  analysis_data[[a]]$xval[, c("smallest", "within_id", "value") := NULL]
 
  # Renaming to keep difference between truncation levels 
  setnames(
    analysis_data[[a]]$xval,
    c("nfuns_sifter", "nproteins_sifter", "score"),
    c(
      paste0("nfuns_", a),
      paste0("nproteins_", a),
      paste0("score_", a)
      )
  )
  
  analysis_data[[a]]$xval[, c("elapsed", "folds") := NULL]
}

# Removing just one to keep family and name after merging
analysis_data[[1]]$xval[, c("family", "#Names") := NULL]

annotations_predictions <- merge(
  x = annotations, 
  y = analysis_data$trunc1$xval[nfuns_trunc1 > 1 & nproteins_trunc1 > 1],
  by.x = c("name", "go"),
  by.y = c("name", "term"),
  all.x = FALSE, all.y = FALSE
)

annotations_predictions <- merge(
  x = annotations_predictions, 
  y = analysis_data$trunc3$xval[nfuns_trunc3 > 1 & nproteins_trunc3 > 1],
  by.x = c("name", "go"),
  by.y = c("name", "term"),
  all.x = FALSE, all.y = FALSE
)

# Versus aphylo
aphylo_predictions <- data.table::fread("data-raw/predictions.csv.gz")
aphylo_predictions[, name2 := gsub(".+=", "", name)]
setnames(aphylo_predictions, "score", "score_aphylo")

# Counting how many annotated
aphylo_predictions[, nproteins_aphylo := length(unique(name)), by = .(term, panther)]

annotations_predictions <- merge(
  x = annotations_predictions,
  y = aphylo_predictions,
  by.x = c("UniProtKB", "go"),
  by.y = c("name2", "term"),
  all.x = TRUE, all.y = FALSE
)

set.seed(1231)

graphics.off()
pdf("data/families_data/predictions_assessment_difference_in_input.pdf", width = 6, height = 5)
op <- par(mai = par("mai") * c(1,1,1,1))
with(annotations_predictions[order(nproteins_aphylo - nproteins_trunc1)], {
  
  barplot(
    nproteins_aphylo - nproteins_trunc1,
    ylab   = "aphylo - SIFTER",
    xlab   = paste("All", nrow(annotations_predictions), " annotations sorted by the difference"),
    # main   = "Difference in the number of input proteins\nused for prediction",
    border = "darkgray",
    col    = "darkgray"
    )
  
  test <- t.test(nproteins_aphylo, nproteins_trunc1, paired = TRUE, var.equal = FALSE)
  
  legend(
    "bottomright",
    title = "Paired t-test",
    legend = c(
      "\tH0: #aphylo - #SIFTER = 0",
      sprintf("\tt-statistic: %.2f", test$statistic),
      sprintf("\tp-value: %.4f", test$p.value)
    ),
    bty = "n"
  )
  
  })
par(op)
dev.off()

annotations_predictions[, name.x := NULL]
annotations_predictions[, name.y := NULL]

setnames(annotations_predictions, c("tree", "family", "#Names"), c("panther", "pfam", "name"))

data.table::fwrite(
  annotations_predictions,
  "data/families_data/predictions_assessment.csv"
  )

# Accuracy ---------------------------------------------------------------------
auc_aphylo <- with(annotations_predictions, auc(score_aphylo, value, nc = 100000))
auc_trunc1 <- with(annotations_predictions, auc(score_trunc1, value, nc = 100000))
auc_trunc3 <- with(annotations_predictions, auc(score_trunc3, value, nc = 100000))

cols <- c("black", "steelblue", "tomato")
graphics.off()
pdf("data/families_data/predictions_assessment.pdf", width = 7, height = 6)
plot(auc_aphylo, lty = 1, lwd = 2, col = cols[1])
with(auc_trunc1, lines(x = fpr, y = tpr, lty = 2, lwd = 2, col = cols[2]))
with(auc_trunc3, lines(x = fpr, y = tpr, lty = 4, lwd = 2, col = cols[3]))
legend(
  "bottomright",
  lty = c(1:2,4),
  legend = sprintf(
    "%s AUC: %.2f; MAE: %.2f",
    c("aphylo", "SIFTER (T=1)", "SIFTER (T=3)"),
    c(auc_aphylo$auc, auc_trunc1$auc, auc_trunc3$auc),
    with(annotations_predictions, c(
      mean(abs(value - score_aphylo)),
      mean(abs(value - score_trunc1)),
      mean(abs(value - score_trunc3))
      )
      )
    ),
  col = cols,
  bty = "n",
  lwd = 2,
  title = "Algorithm"
  )
dev.off()


# Why SIFTER trunc level 3 does better in MAE but not in AUC? ------------------
graphics.off()
pdf("data/families_data/predictions_assessment_sifter_t1_vs_t3.pdf", width = 5, height = 5)
set.seed(123123);annotations_predictions[order(-value, score_trunc3),][, {
  
  xcoords <- jitter(score_trunc1, amount = 0.1)
  ycoords <- jitter(score_trunc3, amount = 0.1)
  
  plot(
    x = xcoords,
    y = ycoords,
    col = adjustcolor(c("tomato", "steelblue"), alpha = .7)[value+1],
    # nbin = 30,
    cex  = 1.5,
    pch  = 17,
    xlab = "Prediction score SIFTER (T = 1)",
    ylab = "Prediction score SIFTER (T = 3)",
  )

  points(
    x = xcoords[value == 0 & score_trunc3 > .2 & score_trunc1 < .5],
    y = ycoords[value == 0 & score_trunc3 > .2 & score_trunc1 < .5],
    cex = 2,
    pch = 2,
    lwd = 2
    )
  
  }]

legend(
  "bottom",
  pch = c(17, 17),
  col = c("tomato", "steelblue"),
  cex = c(1.25, 1.25),
  legend = c("Absent", "Present"),
  bty = "n"
  )
dev.off()
# 6/10 not annotations that were annotated with not now have a non-zero annotation
# 