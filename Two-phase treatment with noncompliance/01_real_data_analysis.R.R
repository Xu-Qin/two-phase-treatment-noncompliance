# Packages
library(mice) # for missing imputation
library(coxed)  # for BCa confidence interval
library(doSNOW) # for parallel computing
library(ggplot2) # for making graphs
library(tictoc) # for recording computational time

# Data Analysis Function
Analyze = function(data.ori, X = c("gender", "race", "gkfreelunch", "age_kentry", "gkrepeat", "gkspeced"), Y = "g1readmath", Z = "Z", D = "D", V = "gkreadmath", nmimp = 10, nboots = 500, ncores = 10, seed = 2025){
  #[input] data.ori: (data.frame) the original data set with missing values
  #[param] nmimp: (integer) number of missing imputations
  #[param] nboots: (integer) number of bootstrap samples
  #[param] ncores: (integer) number of cores used for parallel computing
  #[var] X,Y,Z,D,V: (character) name of the model variables in the data set
  
  # Multiple imputation function for missing data
  Mimputation = function(data.mis, nMimp){
    # Create missing indicators for dichotomous X
    DXmissing = c('gkfreelunch', 'gkrepeat', 'gkspeced')
    for(name in DXmissing){
      data.mis[, paste0(name,'.mis')] = ifelse(is.na(data.mis[, name]), 1, 0)
      data.mis[is.na(data.mis[, name]), name] = 1
    }
    
    # Treat categorical variables as factors
    data.mis$gksurban = as.factor(data.mis$gksurban)
    
    # Set up predictor matrix
    predMatrix = make.predictorMatrix(data = data.mis)
    predMatrix[, c('stdntid', 'gkschid')] = 0
    
    # Impute missing V and Y along with math scores
    imp = mice(data.mis, m = nMimp, predictorMatrix = predMatrix, printFlag = FALSE, seed = 2025)
    
    # Extract the imputed data
    impdatalist = complete(imp, action = 'all')
    return(impdatalist)
  }
  
  # MS2T-IV Estimation Function
  Estimate.MS2T = function(data, ZJk, improper){
    #[input] data: (data.frame) data set for ATE estimation
    #[param] ZJk: (logical) whether to estimate Zitzmann's jackknife variance
    #[param] improper: (logical) whether to estimate the improper variance
    
    # Calculate school-mean-centered covariates
    for(name in X){
      Xkbar = tapply(data[ , name], list(school = data$gkschid), mean)
      data[ , name] = data[ , name] - Xkbar[as.character(data$gkschid)]
    }
    
    # Stage 1
    # Use school-by-school OLS analysis
    K = length(unique(data$gkschid))  # number of schools
    alpha1k <- beta0k <- beta1k <- theta1k <- rep(0, K)
    for(s in 1:K){
      sch = unique(data$gkschid)[s]
      datas = data[which(data$gkschid==sch), ]
      Vmodel = lm(as.formula(paste(V, "~", Z, "*", "(", paste(X, collapse = "+"), ")")), data = datas)
      Dmodel = lm(as.formula(paste(D, "~", Z, "*", "(", paste(X, collapse = "+"), ")")), data = datas)
      Ymodel = lm(as.formula(paste(Y, "~", Z, "*", "(", paste(X, collapse = "+"), ")")), data = datas)
      alpha1k[s] = coef(Vmodel)[2]
      beta0k[s] = coef(Dmodel)[1]
      beta1k[s] = coef(Dmodel)[2]
      theta1k[s] = coef(Ymodel)[2]
    }
    
    # Stage 2
    beta2k = beta0k + beta1k
    model = lm(theta1k ~ beta1k + beta2k + alpha1k)
    ATE = sum(coef(model)*c(1,1,1,mean(alpha1k)))
    
    if(ZJk == T){
      # Zitzmann's jackknife variance
      nout = 4  # creating 19 subgroups for 76 schools
      ATE.ZJk = rep(NA, 19)
      for(j in 1:19){
        kout = (j-1)*nout + (1:nout)
        fit.sub = lm(theta1k[-kout] ~ beta1k[-kout] + beta2k[-kout] + alpha1k[-kout])
        ATE.ZJk[j] = sum(coef(fit.sub)*c(1,1,1,mean(alpha1k[-kout])))
      }
      Var.ZJk = (18/19) * sum((ATE.ZJk - mean(ATE.ZJk))^2)
    }else{
      Var.ZJk = NA
    }
    
    if(improper == T){
      # Improper variance
      Var.improper = as.vector(t(c(1,1,1,mean(alpha1k)))%*%vcov(model)%*%c(1,1,1,mean(alpha1k)))
    }else{
      Var.improper = NA
    }
    
    return(list(ATE = ATE, Var.ZJk = Var.ZJk, Var.improper = Var.improper))
  }
  
  # Safer sampling function when the original sample could have only one element
  resample = function(x, replace){x[sample.int(length(x), replace = replace)]}
  
  # Bootstrap Function
  bootstrap = function(data.ori, nboots){
    #[input] data.ori: (data.frame) the original data set
    #[param] nboots: (numeric) number of bootstrap samples
    
    bsamples = list()
    for(b in 1:nboots){
      sampleid = c()
      samplesch = c()
      # Bootstrap schools
      schsample = sample(unique(data.ori$gkschid), replace = T)
      for(s in 1:length(schsample)){
        # Bootstrap students within each treatment group in each school
        samp = c(resample(which(data.ori$gkschid==schsample[s] & data.ori$Z==0), replace = T), 
                 resample(which(data.ori$gkschid==schsample[s] & data.ori$Z==1), replace = T))
        sampleid = c(sampleid, samp)
        samplesch = c(samplesch, rep(s, length(samp)))
      }
      bsamples[[b]] = list(id = sampleid, sch = samplesch)
    }
    return(bsamples)
  }
  
  set.seed(seed)
  
  # Missing imputation
  tic('Missing imputation')
  
  impdatalist = Mimputation(data.ori, nmimp)
  
  toc()
  
  # Estimation and inference
  tic('Estimation, Zitzmann Jackknife CI and Improper CI')
  
  ATE.imp = rep(NA, nmimp)
  ZJkVar.imp = rep(NA, nmimp)
  ImproperVar.imp = rep(NA, nmimp)
  for(m in 1:nmimp){
    result = Estimate.MS2T(impdatalist[[m]], ZJk = T, improper = T)
    ATE.imp[m] = result$ATE
    ZJkVar.imp[m] = result$Var.ZJk
    ImproperVar.imp[m] = result$Var.improper
  }
  ATE = mean(ATE.imp)
  CI.ZJk = c(ATE - qnorm(0.975)*sqrt(mean(ZJkVar.imp)),
             ATE + qnorm(0.975)*sqrt(mean(ZJkVar.imp)))  # 95% CI
  CI.improper = c(ATE - qnorm(0.975)*sqrt(mean(ImproperVar.imp)),
             ATE + qnorm(0.975)*sqrt(mean(ImproperVar.imp)))
  
  toc()

  tic('Bootstrap Inference')

  bsamples = bootstrap(data.ori, nboots)

  # Prepare for parallel computing
  cl = makeSOCKcluster(ncores)
  registerDoSNOW(cl)
  # Create a progress bar
  pb = txtProgressBar(max = nboots, style = 3)
  progress = function(n){setTxtProgressBar(pb, n)}
  opts = list(progress = progress)

  ATE.boots = foreach(b = 1:nboots, .combine = c, .options.snow = opts) %dopar% {
    ATEb.imp = rep(NA, nmimp)
    for(m in 1:nmimp){
      data.imp = impdatalist[[m]]  # imputed data file
      datab = data.imp[bsamples[[b]]$id, ]  # bootstrap sample b
      datab$gkschid = bsamples[[b]]$sch  # re-assign school id
      ATEb.imp[m] = Estimate.MS2T(datab, ZJk = F, improper = F)$ATE
    }
    mean(ATEb.imp)  # average over 10 missing imputations
  }

  close(pb)
  stopCluster(cl)

  CI = quantile(ATE.boots, probs = c(0.025, 0.975))
  CI.BCa = bca(ATE.boots)

  toc()

  results = c(ATE, CI.ZJk, CI, CI.BCa, CI.improper)
  names(results) = c('ATE', 'CI.ZJk.lower', 'CI.ZJk.upper', 'CI.boots.lower', 'CI.boots.upper', 'CI.BCa.lower', 'CI.BCa.upper', 'CI.Imp.lower', 'CI.Imp.upper')
  
  return(results)
}

# Load the data with missing values
data.ori = read.csv('application_data_2phase.csv')

# Analyze the data
Analyze(data.ori)
