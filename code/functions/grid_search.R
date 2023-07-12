grid_search <- function(data, response, family){
  trainBRT(data = data,
           preds = vars,
           resp = response,
           family = family,
           treeComplexity = c(1, 5, 10),
           learningRate = c(0.01, 0.05, 0.1),
           bagFraction = c(0.25, 0.5, 0.75),
           minTrees = 1000, # recommended minimum by Elith
           maxTrees = 2500,
           cores = 6, # increase speed
           out = c('model', 'tuning')) # should return model and table with hyperparameters
}
