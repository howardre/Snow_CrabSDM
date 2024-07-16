grid_search <- function(data, response, family){
  trainBRT(data = data,
           preds = c(1:9),
           resp = response,
           family = family,
           treeComplexity = c(1, 5, 10),
           learningRate = c(0.01, 0.05, 0.1),
           bagFraction = c(0.25, 0.5, 0.75),
           step.size = 25,
           minTrees = 1000, # recommended minimum by Elith
           maxTrees = 2500,
           tryBy = c('learningRate', 'treeComplexity', 'bagFraction'),
           cores = 6, # increase speed,
           out = c('tuning', 'model')) # should return model and table with hyperparameters
}
