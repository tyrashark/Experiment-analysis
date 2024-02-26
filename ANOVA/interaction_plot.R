interaction_plot <- function(aov, predictor.scope=1:2, predictors=NULL, drop.blockcom=T,...){
  if(is.null(predictors)){
    groups = list()
    models = list()
    models[[1]] = aov$model[, 1]
    for(i in 1:length(aov$xlevels)){
      groups[[i]] = length(aov$xlevels[[i]])
      models[[i+1]] = aov$model[, (i+1)]
    }
    x_names = names(aov$xlevels)
    x_names = gsub(pattern = "as.factor\\(", replacement = "", x_names)
    x_names = gsub(pattern = "\\)", replacement = "", x_names)
  }
  
  # Generate all combinations of predictors with size 2
  combinations = combn(predictor.scope, 2)
  if(drop.blockcom) combinations = combinations[,1:(length(aov$xlevels)-1)]
  
  seq_const_start = cumsum(c(2, unlist(groups)[-length(groups)]-1))
  seq_const_end = cumsum(unlist(groups)-1)+1
  
  # reparametrization of the population mean
  weight = rep(1/unlist(groups), unlist(groups)-1)
  mu = aov$coefficients[1] + crossprod(aov$coefficients[-1], weight)[1]
  
  for(j in 1:ncol(combinations)){
    combi = combinations[,j]
    group_cols = RColorBrewer::brewer.pal(groups[[combi[2]]], name = "Dark2")
    interaction.plot(models[[(combi[1]+1)]], models[[(combi[2]+1)]], models[[1]], type = "p", 
                     xlab = x_names[combi][1], ylab = names(aov$model)[1],
                     trace.label = x_names[combi][2],
                     col = group_cols, legend = T, pch = 1:groups[[combi[2]]])
    
    ## Predict with reparametrization
    coeff_end = seq_const_start[combi]
    coeff_start = seq_const_end[combi]
    
    repara1 = c(0, aov$coefficients[coeff_start[1]:coeff_end[1]])
    repara1_1 = -sum(repara1)/length(repara1)
    repara1 = repara1 + repara1_1 
    
    repara2 = c(0, aov$coefficients[coeff_start[2]:coeff_end[2]])
    repara2_1 = -sum(repara2)/length(repara2)
    repara2 = repara2 + repara2_1 
    
    for(i in 2:groups[[combi[1]]]){
      y0predict = mu + repara1[i - 1] + repara2
      y1predict = mu + repara1[i] + repara2
      segments(x0 = rep(i - 1, groups[[combi[2]]]), 
               x1 = rep(i, groups[[combi[2]]]),
               y0 = y0predict,
               y1 = y1predict,
               col = group_cols, lty = 1:groups[[combi[2]]])
    }
  }
}
