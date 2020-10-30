# This script generates the raw predictions using leave one out cross validation
library(aphylo)
library(data.table)

estimates <- readRDS("data-raw/mcmc_partially_annotated_no_prior.rds")

# Applying burnin
estimates <- window(
  estimates,
  floor((end(estimates$hist) - start(estimates$hist))/2)
)

# Ids (positions) of the genes that are annotated
ids <- lapply(estimates$dat, function(i) {
  which(i$tip.annotation != 9)
})

# ------------------------------------------------------------------------------
# Proteins included in SIFTER comparison

proteins_in_sifter <- c("A2AHL1", "A3KPS8", "D3XD69", "D4A915", "F1QHM9", "F1R1P2", 
                        "F4JJL3", "G3V6K6", "O00445", "O08835", "O15197", "O22182", "O22183", 
                        "O35681", "O48652", "O59921", "O60568", "O60656", "O61460", "O64692", 
                        "O75310", "O75795", "P00533", "P00914", "P04626", "P04629", "P05066", 
                        "P05708", "P06133", "P06213", "P06494", "P08069", "P08430", "P08541", 
                        "P08542", "P08922", "P09619", "P09759", "P09875", "P11362", "P13368", 
                        "P15127", "P16662", "P17709", "P17710", "P17712", "P18475", "P19488", 
                        "P21521", "P21707", "P22309", "P22607", "P24062", "P29101", "P29323", 
                        "P30530", "P31677", "P31688", "P35503", "P35557", "P35739", "P35968", 
                        "P36511", "P38426", "P38427", "P40387", "P40748", "P40749", "P46096", 
                        "P46097", "P47861", "P50521", "P52792", "P54753", "P54759", "P54760", 
                        "P54762", "P54763", "P54855", "P78875", "P86044", "P93834", "P97523", 
                        "P9WN11", "Q00342", "Q00764", "Q02763", "Q02858", "Q05030", "Q06418", 
                        "Q09756", "Q14AT5", "Q15303", "Q18779", "Q32M45", "Q42525", "Q4KML2", 
                        "Q4WA18", "Q4WCS2", "Q4WHC3", "Q5BGE3", "Q5T7P8", "Q5W676", "Q5XXA6", 
                        "Q61527", "Q62746", "Q62838", "Q63116", "Q64435", "Q64550", "Q64633", 
                        "Q6NUS8", "Q6P9J9", "Q6R6I6", "Q6UWM9", "Q8AXB3", "Q8BHY3", "Q8C5H1", 
                        "Q8CFW1", "Q8LEA2", "Q8LQ68", "Q8VZE9", "Q9BY64", "Q9C521", "Q9FE68", 
                        "Q9FZG4", "Q9HAW7", "Q9HBX9", "Q9LR44", "Q9LSY6", "Q9NQ90", "Q9R0E1", 
                        "Q9R0N4", "Q9R0N5", "Q9R0N7", "Q9R0N8", "Q9SXA6", "Q9SYM4", "Q9U6P7", 
                        "Q9UM73", "Q9V6K3", "Q9VBP0", "Q9VBW3", "Q9VEG4", "Q9XFR9", "Q9Y119", 
                        "Q9ZVY5")

ids_in_sifter <- lapply(estimates$dat, function(e) {
  which(gsub(".+=", "", e$tree$tip.label) %in% proteins_in_sifter)
})
ids_in_sifter <- which(sapply(ids_in_sifter, length) > 0)

microbenchmark::microbenchmark(
  loo_cv = predict(estimates, which.tree = ids_in_sifter, loo = TRUE, ids = ids[ids_in_sifter]),
  times  = 10
)
# Unit: milliseconds
#      min       lq     mean   median      uq      max neval
# 662.9263 675.7675 688.9179 679.9534 701.786 727.2798    10

microbenchmark::microbenchmark(
  predict_all = (predict_all <- predict(estimates, which.tree = ids_in_sifter, loo = FALSE)),
  times = 10
)
# Unit: milliseconds
#        expr      min      lq     mean   median       uq      max neval
# predict_all 147.3905 157.915 175.3287 172.9496 196.1166 209.8919    10

# How many predictions
predict_all <- do.call(rbind, predict_all)
table(grepl("^UniP", rownames(predict_all))) # 18,027
