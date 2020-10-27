library(data.table)
library(aphylo)

# Reading predictions in
fn_predictions <- list.files("data/families_data/predictions", pattern = "*val*", full.names = TRUE)
out            <- lapply(fn_predictions, readRDS)
names(out)     <- gsub(".+predictions/PF([0-9]+)-xval\\.rds$", "PF\\1", fn_predictions)

fn_pred_noxval <- list.files("data/families_data/predictions", pattern = "*.rds", full.names = TRUE)
fn_pred_noxval <- fn_pred_noxval[grepl("PF[0-9]+\\.rds", fn_pred_noxval)]
out_no_xval    <- lapply(fn_pred_noxval, readRDS)
names(out_no_xval) <- gsub(".+predictions/PF([0-9]+)\\.rds$", "PF\\1",fn_pred_noxval)

# Are the cross validated too different? ---------------------------------------
common <- intersect(names(out), names(out_no_xval))

# Did the analysis ran with enough information? --------------------------------
proteins <- lapply(names(out), function(i) {
  read_pli(sprintf("data/families_data/annotations/%s.pli", i))
})

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
  data.table(family = f, xval[[f]])
})
xval <- rbindlist(xval, fill = TRUE)
xval$`#Names` <- toupper(xval$`#Names`)

# Reshaping
xval <- melt(
  xval,
  id            = 1:2,
  value.name    = "score",
  variable.name = "term"
)
xval <- xval[!is.na(score)]
xval[, name:=gsub("/.+", "", `#Names`)]
xval[, score := as.double(score)]

xval[, nfuns := length(unique(term)), by = family]
xval[, nproteins := length(unique(name)), by=family]


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

# Retrieving the predictions ---------------------------------------------------
predictions <- lapply(out, "[[", "predictions")
predictions <- lapply(names(predictions), function(f) {
  data.table(family = f, predictions[[f]])
})

names(predictions) <- names(out)

# PF00754 has two columns, but the second column has only missing values
single_column <- which(
  sapply(predictions, function(i) any(is.na(i))) |
    sapply(predictions, function(p) {
      sum(grepl("^GO[:]",colnames(p)))
    }) == 1
  )
predictions <- predictions[-single_column]

predictions <- rbindlist(predictions, fill = TRUE)

# Reshaping
predictions <- melt(
  predictions,
  id            = 1:2,
  value.name    = "score",
  variable.name = "term"
  )
predictions <- predictions[!is.na(score)]
predictions[, name:=gsub("_.+", "", `#Names`)]

# How many per family
predictions[, nfuns := length(unique(term)), by = family]
predictions[, nproteins := length(unique(name)), by=family]

# Listing genes that it was not able to find -----------------------------------
missing_proteins <- unlist(lapply(out, "[[", "out"))
missing_proteins <- missing_proteins[grepl("^In hasNode ", missing_proteins)]
missing_proteins <- unique(gsub(".+\\sID[:]", "", missing_proteins))

cbind(missing_proteins)
#    missing_proteins   
# [1,] " gper1_danre"     
# [2,] " gper1_mouse"     
# [3,] " fgfr1_human"     
# [4,] " vgfr2_human"     
# [5,] " vgfr4_danre"     
# [6,] " lirb3_mouse"     
# [7,] " p3c2g_mouse"     
# [8,] " hp302_arath"     
# [9,] " ddrb_caeel"      
# [10,] " max2_caeel"      
# [11,] " a0a0b4kgs4_drome"
# [12,] " pk1_caeel"       
# [13,] " utr7_arath"      
# [14,] " urgt4_arath"     
# [15,] " urgt1_arath"     
# [16,] " gons1_arath"     
# [17,] " g2ox4_arath"     
# [18,] " g2ox6_arath"     
# [19,] " a0a0b4jcy1_drome"
# [20,] " t2r40_chick"     
# [21,] " tgfr2_human"     
# [22,] " acvr1_mouse"     
# [23,] " acvl1_human"     
# [24,] " bmr1a_mouse"     
# [25,] " tgfr2_mouse"     
# [26,] " svh2_caeel"      
# [27,] " acvr1_human"     
# [28,] " tgfr2_chick"     
# [29,] " piwl1_mouse"     
# [30,] " piwl2_mouse" 

# Retrieving the training data -------------------------------------------------

annotations <- readRDS("data/families_data/annotations.rds")
setnames(annotations, "qualifier", "value")
annotations <- unique(annotations[, fami := NULL]) # Don't need this

# Merging the data and computing AUCs! -----------------------------------------

# Since in some cases we have more than tree per protein (multiple hits in pfam)
# we will keep the best
xval <- xval[nfuns > 1 & nproteins > 1]
xval <- merge(
  xval, 
  annotations[, .(term = go, name, value)],
  by    = c("term", "name"),
  all.x = TRUE,
  all.y = FALSE
)
xval <- xval[!is.na(value)]

xval[, smallest := abs(score - value), by = .(term, name)]
xval[, smallest := which.min(smallest), by = .(term, name)]
xval[, within_id := 1:.N, by = .(term, name)]

xval <- xval[smallest == within_id]
xval[, c("smallest", "within_id", "value") := NULL]

annotations_predictions <- merge(
  x = annotations, 
  y = xval[nfuns > 1 & nproteins > 1],
  by.x = c("name", "go"),
  by.y = c("name", "term"),
)
annotations_predictions

# MAE
with(annotations_predictions, mean((value - score)))

# Versus aphylo
aphylo_predictions <- data.table::fread("data-raw/predictions.csv.gz")
aphylo_predictions[, name2 := gsub(".+=", "", name)]
setnames(aphylo_predictions, "score", "score_aphylo")

annotations_predictions <- merge(
  x = annotations_predictions,
  y = aphylo_predictions,
  by.x = c("UniProtKB", "go"),
  by.y = c("name2", "term"),
  all.x = TRUE, all.y = FALSE
)

with(annotations_predictions, mean((value - score_aphylo)))

set.seed(1231)
with(
  annotations_predictions,
  plot(
    x    = jitter(score, amount = .05),
    y    = jitter(score_aphylo, amount = .05),
    xlim = c(0, 1.1),
    ylim = c(0, 1.1),
    xlab = "Prediction by SIFTER",
    ylab = "Prediction by aphylo",
    main = "Prediction Scores\non all positive annotations",
    sub  = sprintf(
      "This includes %i annotations on %i proteins",
      nrow(annotations_predictions),
      length(unique(annotations_predictions$UniProtKB))
      ),
    col  = c("tomato", "steelblue")[value + 1]
    )
  )

legend(
  "bottomleft",
  legend = sprintf(
    "%s MAE: %.2f; AUC: %.2f",
    c("aphylo", "SIFTER"), 
    with(annotations_predictions, c(mean(abs(value - score_aphylo)), mean(abs(value - score)))),
    with(annotations_predictions, c(auc(score_aphylo,value)$auc, auc(score, value)$auc))
    ),
  bty = "n"
  )


auc_aphylo <- with(annotations_predictions, auc(score_aphylo, value, nc = 100000))
auc_sifter <- with(annotations_predictions, auc(score, value, nc = 100000))

cols <- c("steelblue", "tomato")
graphics.off()
pdf("data/families_data/predictions_assessment.pdf", width = 7, height = 6)
plot(auc_aphylo, lty = 1, lwd = 1.75, col = cols[1])
with(auc_sifter, lines(x = fpr, y = tpr, lty = 2, lwd = 1.75, col = cols[2]))
legend(
  "bottomright",
  lty = 1:2,
  legend = sprintf(
    "%s AUC: %.2f; MAE: %.2f",
    c("aphylo", "SIFTER"),
    c(auc_aphylo$auc, auc_sifter$auc),
    with(annotations_predictions, c(mean(abs(value - score_aphylo)), mean(abs(value - score))))
    ),
  col = cols,
  bty = "n",
  lwd = 1.75,
  title = "Algorithm"
  )
dev.off()
