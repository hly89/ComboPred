# leave-one-out for predicting single compound/combination response
# rm.idx: the indexes for being left out 
LooPred <- function(training.data, response, rm.idx) {
  loo.ctrl <- matrix(NA, nrow = 1, ncol = 3)
  colnames(loo.ctrl) <- c("id", "compound", "DSS.pred")
  for(j in rm.idx) {
    # remove from response data
    tmp.resp <- response[-j, ]
    # remove from training data
    tmp.training <- training.data[-j, ]
    tmp.model <- randomForest(tmp.training, tmp.resp$dss)
    tmp <- matrix(NA, nrow = 1, ncol = 3)
    colnames(tmp) <- c("id", "compound", "DSS.pred")
    tmp[1, 1] <- response$id[j]
    tmp[1, 2] <- response$compound[j]
    set.seed(1)
    tmp[1,3] <- predict(tmp.model, training.data[j,])
    loo.ctrl <- rbind(loo.ctrl, tmp)
  }
  
  loo.ctrl <- loo.ctrl[-1, ]
  loo.ctrl <- data.frame(loo.ctrl, stringsAsFactors = F)
  loo.ctrl$id <- as.numeric(loo.ctrl$id)
  loo.ctrl$DSS.pred <- as.numeric(loo.ctrl$DSS.pred)
  return(loo.ctrl)
}