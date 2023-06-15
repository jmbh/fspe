# jonashaslbeck@gmail.com; January 19, 2021

fspe <- function(data,
                 maxK,
                 nfold = 10,
                 rep = 1,
                 method = "PE",
                 rotate = "oblimin",
                 pbar = TRUE,
                 ...) {


  # ----- Get basic info ------
  p <- ncol(data)
  n <- nrow(data)


  # ----- Input checks ------

  if(!inherits(data, c("matrix", "data.frame"))) stop("The data has to be provided as a matrix or data.frame.")


  # ----- Check identifiability ------
  # Make sure all considered factor models are in fact identified

  maxK_cell <- NA
  for(q in 1:maxK) {
    k <- q
    df <- p*(p-1)/2 - (k*p - k*(k-1)/2)
    if(df>0) maxK_cell <- q else break
  }

  if(maxK > maxK_cell) warning(paste0("Some of the candidate factor models are not identified with ", p," Variables. \nThe maximum number of factors considered was set to ", maxK_cell,"."))

  maxK_use <- maxK_cell


  # ----- Scale data ------
  data <- apply(data, 2, scale)


  # ------ Set up progress bar ------

  if(pbar==TRUE) if(rep>1) {
    pb <- txtProgressBar(min = 0, max=rep, initial=0, char="-", style = 3)
  } else {
    pb <- txtProgressBar(min = 0, max=maxK_use, initial=0, char="-", style = 3)
  }



  # ----- Method 1: Compute OoS PE + k-fold cross validation ------

  if(method == "PE") {

    a_foldPE <- array(NA, dim=c(maxK_use, nfold, p, rep))

    for(r in 1:rep){

      # Make folds
      times_n <- ceiling(n/nfold)
      id <- rep(1:nfold, times=times_n, each=1)[1:n]
      id_random <- sample(id, size=n, replace=FALSE) # shuffle

      for(k in 1:maxK_use) {
        for(f in 1:nfold) {

          # CV test/train split
          cv_data_train <- data[id_random!=f, ]
          cv_data_test <- data[id_random==f, ]

          # fit factor model
          fit <- fa(cv_data_train, nfactors = k, rotate = rotate)

          # get model implied partial correlations
          mi_cor <- impCov(fit)
          mi_pcor <- cor2pcor(mi_cor)

          # compute prediction error for each variable
          for(i in 1:p) {
            i_pred <- as.matrix(cv_data_test[, -i]) %*% matrix(mi_pcor[i, -i], ncol=1)
            a_foldPE[k, f, i, r] <- mean((cv_data_test[,i]-i_pred)^2)
          } # end for i

          if(pbar) if(rep==1) if(pbar==TRUE) setTxtProgressBar(pb, k)
        } # end for: k

      } # end for: fold

      if(pbar) if(rep>1) if(pbar==TRUE) setTxtProgressBar(pb, r)

    } # end for : rep

    # ----- Obtain Aggregate Prediction ------

    # Aggregation: aggregate within factor, then select
    m_foldPE_2 <- apply(a_foldPE, 1, mean)
    vote2 <- which.min(m_foldPE_2)

    nfactor <- vote2

  } # end if: method == PE


  # ----- Method 2: Compute OoS PE on covariance matrix + k-fold cross validation ------

  if(method == "Cov") {

    a_foldPE <- array(NA, dim=c(maxK_use, nfold, rep))

    for(r in 1:rep){

      # Make folds
      times_n <- ceiling(n/nfold)
      id <- rep(1:nfold, times=times_n, each=1)[1:n]
      id_random <- sample(id, size=n, replace=FALSE) # shuffle

      for(k in 1:maxK_use) {
        for(f in 1:nfold) {

          # cv test/train split
          cv_data_train <- data[id!=f, ]
          cv_data_test <- data[id==f, ]

          # fit factor model
          fit <- fa(cv_data_train, nfactors = k, rotate = rotate, ...)

          # get model-implied cov
          mi_cor <- impCov(fit)
          oos_cor <- cor(cv_data_test)
          a_foldPE[k, f, r] <- mean((oos_cor - mi_cor)^2)

          if(pbar) if(rep==1) if(pbar==TRUE) setTxtProgressBar(pb, k)

        } # end for: k

      } # end for: fold

      if(pbar) if(rep>1) if(pbar==TRUE) setTxtProgressBar(pb, r)

    } # end for : rep

    # Aggregate
    v_foldMSE <- apply(a_foldPE, 1, mean)
    nfactor <- which.min(v_foldMSE)


  } # end if: method == Cov


  # ----- Prepare output ------

  outlist <- list("nfactor" = nfactor,
                  "PEs" = apply(a_foldPE, 1, mean),
                  "PE_array" = a_foldPE)

  return(outlist)

} # eoF


