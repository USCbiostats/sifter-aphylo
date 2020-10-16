library(aphylo)

# Loading parameter estimates from the previous model
model_aphylo <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)
plot(model_aphylo$hist)

# Loading the data
hundredfamilies <- readRDS("sifter/hundredfamilies.rds")

ans <- vector("list", length(hundredfamilies))
for (i in 1:length(ans)) {
  
  # Identifying which are not NAs
  ids <- unname(which(
    rowSums(hundredfamilies[[i]]$tip.annotation) !=
      Nann(hundredfamilies[[i]]) * 9))
  
  # Making the predictions
  tmp <- predict(model_aphylo, newdata = hundredfamilies[[i]], loo = TRUE, ids = ids)[ids,,drop=FALSE]
  
  # Retrieving the data
  ans[[i]] <- cbind(
    label = as.vector(hundredfamilies[[i]]$tip.annotation[ids,]),
    pred  = as.vector(tmp)
  )
  message("Tree N ", i, " done.")
}

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
