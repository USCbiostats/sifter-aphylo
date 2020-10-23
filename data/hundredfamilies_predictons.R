library(aphylo)

# Loading parameter estimates from the previous model
model_aphylo <- readRDS("data-raw/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)
# plot(model_aphylo$hist)

# Loading the data
hundredfamilies <- readRDS("data/hundredfamilies.rds")

ans <- vector("list", length(hundredfamilies))
acc_at_least_one <- ans
for (f in 1:length(ans)) {
  
  # Identifying which are not NAs
  ids <- unname(which(
    rowSums(hundredfamilies[[f]]$tip.annotation) !=
      Nann(hundredfamilies[[f]]) * 9))
  
  # Making the predictions
  tmp <- predict(model_aphylo, newdata = hundredfamilies[[f]], loo = TRUE, ids = ids)
  
  # Computing accuracy
  acc_at_least_one[[f]] <- lapply(ids, function(i) {
    top_pred <- which(abs(tmp[i,] - max(tmp[i,])) < 1e-10)
    true_ann <- which(hundredfamilies[[f]]$tip.annotation[i,] == 1)
    data.frame(
      pred = paste(top_pred, collapse=","),
      true = paste(true_ann, collapse=","),
      accurate = ifelse(length(intersect(top_pred, true_ann)), 1, 0)
    )
  })
  
  # Retrieving the data
  ans[[f]] <- cbind(
    label = as.vector(hundredfamilies[[f]]$tip.annotation[ids,]),
    pred  = as.vector(tmp[ids,])
  )
  message("Tree N ", f, " done.")
}

acc_at_least_one <- lapply(acc_at_least_one, do.call, what = "rbind")

acc <- do.call(rbind, acc_at_least_one)
table(acc$accurate)

# In the SIFTER 2011 paper, they assuming missing annotations as NOTs
# that is how, at least I believe, they calculated the True Negatives,
# which they do feature in the paper.
ans <- do.call(rbind, ans)
ans[ans[,1] == 9,1] <- 0L

plot(print(aphylo::auc(pred = ans[,2], labels=ans[,1])))

# library(AUC)
# 
# accuracy1 <- accuracy(predictions = ans[,2], labels = factor(ans[,1]))
# auc(accuracy1)
