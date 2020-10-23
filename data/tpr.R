library(aphylo)

# Loading the data
deaminase        <- readRDS("data/deaminase.rds")
sulfotransferase <- readRDS("data/sulfotransferase.rds")
nudix            <- readRDS("data/nudix.rds")

cols <- c(deaminase = "black", sulfotransferase = "steelblue", nudix = "tomato")

# ROC curves -------------------------------------------------------------------
plot(deaminase$aucs, lwd = 1.5, col = cols[1])
with(sulfotransferase$aucs, lines(x = fpr, y = tpr, lty = 2, lwd = 1.5, col = cols[2]))
with(nudix$aucs, lines(x = fpr, y = tpr, lty = 3, lwd = 1.5, col = cols[3]))

legend(
  "bottomright", title = "Family",
  legend = names(cols),
  col    = cols,
  lty    = 1:3,
  lwd    = 1.5,
  bt     = "n"
  )


# True positive rates as a function of the cuttoff -----------------------------
with(deaminase$aucs, plot(x = cutoffs, y = tpr, lwd = 1.5, col = cols[1], 
    main = "True Positive Rate", ylab = "TPR", xlab = "Cutoff threshold",
    type = "l"
  )
)
with(sulfotransferase$aucs, lines(x = cutoffs, y = tpr, lty = 2, lwd = 1.5, col = cols[2]))
with(nudix$aucs, lines(x = cutoffs, y = tpr, lty = 3, lwd = 1.5, col = cols[3]))

legend(
  "bottomleft", title = "Family",
  legend = names(cols),
  col    = cols,
  lty    = 1:3,
  lwd    = 1.5,
  bt     = "n"
)
