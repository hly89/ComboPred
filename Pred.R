Pred <- function(dtm, gex, mut, models, model.iteration,
                 response, id, dss.pred){
  comp.list <- response$compound[which(response$id==id)]
  pred.res <- .combospred(models, dtm, gex, mut, comp.list, 
                         id, model.iteration)
  pred.res$mean.pred <- rowMeans(pred.res[, 4:(3+model.iteration)])
  dss.pred <- dss.pred[which(dss.pred$id ==id),]
  pred.res$dss.pred1 <- dss.pred$DSS.pred[match(pred.res$compound1, dss.pred$compound)]
    
  pred.res$dss.pred2 <- dss.pred$DSS.pred[match(pred.res$compound2, dss.pred$compound)]
  pred.res$highest.dss <- ifelse(pred.res$dss.pred1>pred.res$dss.pred2,
                                 pred.res$dss.pred1, pred.res$dss.pred2)
  pred.res$synergy.pred <- pred.res$mean.pred - pred.res$highest.dss
  res.final <- cbind(pred.res$id, pred.res$compound1, pred.res$compound2,
                     pred.res$mean.pred, pred.res$dss.pred1, pred.res$dss.pred2,
                     pred.res$highest.dss, pred.res$synergy.pred)
  colnames(res.final) <- c("ID", "Compound1", "Compound2", "Combo Response Prediction",
                           "Compound1 Response Prediction", "Compound2 Response Prediction",
                           "Highest Single Response", "Synergy Prediction")
  write.csv(res.final, file = paste(id, "predictions.csv", sep = "."),row.names = F)
}

.combospred <- function(models, dtm, gex, mut, comp.list, id, iter){
  num.compound <- length(comp.list)
  res <- matrix(NA, nrow = 1, ncol = (3 + iter))
  iter.names <- paste("iter", c(1:iter), sep = "")
  colnames(res) <- c("id", "compound1", "compound2", iter.names)
  for(i in 1:(num.compound - 1)) {
    for(j in (i+1):num.compound) {
      # assemble gene expression, mutation, and drug-target profiles together
      # dtm features
      comp1 <- comp.list[i]
      comp2 <- comp.list[j]
      #rname <- paste(response$compound, response$id, sep = ".")
      dtm.features <- dtm[match(comp1, rownames(dtm)), ] + 
        dtm[match(comp2, rownames(dtm)), ]
      #rownames(dtm.features) <- rname
      colnames(dtm.features) <- colnames(dtm)
      dtm.features[which(dtm.features>1, arr.ind = T)] <- 1
      
      # gene expression features
      gene.features <- gex[match(id, rownames(gex)),]
      
      # mutation features
      mut.features <- mut[match(id, rownames(mut)),]
      
      training.data <- cbind(dtm.features, gene.features, 
                             mut.features)
      tmp <- matrix(NA, nrow = 1, ncol = (3 + iter))
      colnames(tmp) <- c("id", "compound1", "compound2", iter.names)
      for(t in 1:iter){
        tmp[1, 3+t] <- predict(models[[t]], training.data)
      }
      tmp[1, 1] <- id
      tmp[1, 2] <- comp1
      tmp[1, 3] <- comp2
      res <- rbind(res, tmp)
    }
  }
  res <- res[-1, ]
  res <- data.frame(res, stringsAsFactors = F)
  res$id <- as.numeric(id)
  for(t in 1:iter) {
    res[, 3+t] <- as.numeric(res[, 3+t])
  }
  return(res)
}