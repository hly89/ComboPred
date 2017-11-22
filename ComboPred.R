# gex: gene expression data, each column is for one sample over all the genes.
# mut: mutation data, same sample order as the response data
# dtm: compounds are in the rows while targets are in the columns
# model.iteration: how many iterations for training the model.
# response: sample IDs, sample type (patient/control),compound names, DSS.
# patient.index: specify which patient. 
# control.index: specify which control.
ComboPred <- function(gex, mut, dtm,
                      model.iteration = 6, response, 
                      patient.index = 1, 
                      control.index = 1) {
  
  # patient ID(s)
  p.id <- unique(response$id[which(response$type=="patient")])
  # num.patient: number of patient samples.
  num.patient <- length(p.id)
  if(length(patient.index)>num.patient) {
    stop("Please check patient.index, the maxium number of patient indexes 
         are the number of patients")
  }
  if(anyNA(match(patient.index, 1:length(num.patient)))) {
    stop("Not all the indexes in patient.index exist!")
  }
  # control ID(s)
  ctrl.id <- unique(response$id[which(response$type=="control")])
  # num.ctrl: number of control samples.
  num.ctrl <- length(ctrl.id)
  
  if(length(control.index)>num.ctrl) {
    stop("Please check patient.index, the maxium number of patient indexes 
         are the number of patients")
  }
  if(anyNA(match(control.index, 1:length(num.ctrl)))) {
    stop("Not all the indexes in patient.index exist!")
  }
  
  # num.compound: number of compounds.
  num.compound <- nrow(dtm)
  
  
  # assemble gene expression, mutation, and drug-target profiles together
  # dtm features
  rname <- paste(response$compound, response$id, sep = ".")
  #dtm.features <- matrix(NA, nrow = nrow(response), ncol = ncol(dtm))
  dtm.features <- dtm[match(response$compound, rownames(dtm)), ]
  rownames(dtm.features) <- rname
  colnames(dtm.features) <- colnames(dtm)
  
  # gene expression features
  gene.features <- gex[match(response$id, rownames(gex)),]
  
  # mutation features
  mut.features <- mut[match(response$id, rownames(mut)),]
  
  training.data <- cbind(dtm.features, gene.features, 
                         mut.features)
  # training model
  models <- list()
  for(i in 1:model.iteration) {
    set.seed(i)
    models[[i]] <- randomForest(training.data, response$dss)
  }
  
  # LOO for single compound DSS prediction
  # LOO for control(s)
  comp.idx <- c()
  for(i in ctrl.id){
    # get compounds for this control
    comp.idx <- c(comp.idx, which(response$id==i))
  }
  loo.ctrl <- LooPred(training.data, response, comp.idx)
  # LOO for patient(s)
  comp.idx <- c()
  for(i in p.id){
    # get compounds for this control
    comp.idx <- c(comp.idx, which(response$id==i))
  }
  loo.p <- LooPred(training.data, response, comp.idx)

  
  # combo responses prediction for all possible combos
  # for patients
  
  for (i in p.id) {
    Pred(dtm, gex, mut, models, model.iteration,
                     response, i, loo.p)
    
  }
  # for controls
  for(i in ctrl.id) {
    Pred(dtm, gex, mut, models, model.iteration,
         response, i, loo.ctrl)
  }
}





