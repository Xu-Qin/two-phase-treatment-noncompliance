# The parameter names are the same as in the paper if not otherwise specified.

# Packages
library(MASS)  # for multivariate Normal distribution
library(coxed)  # for BCa confidence interval
library(doSNOW) # for parallel computing
library(tictoc) # for recording computational time

#### Data Generation ####
# Data Generation Function
gendata = function(K, nk){
  #[param] K: (numeric) number of schools
  #[param] nk: (numeric) number of students per school
  
  # Sample size
  total = K * nk
  
  # School and student ID
  k = rep(1:K, each = nk)  # school ID
  i = rep(1:nk, K)  # (within school) student ID
  
  # Kindergarten treatment assignment Z
  PZk = runif(K, 0.25, 0.35)
  # Z = rbinom(total, size = 1, prob = PZk[k])
  Tsize = round(nk * PZk)
  Z = rep(0, total)
  for (s in 1:K){
    Zk = rep(0, nk)
    Zk[sample(1:nk, Tsize[s])] = 1
    Z[which(k==s)] = Zk
  }
  
  # Unobserved pre-treatment confounder U
  PUk = runif(K, 0.25, 0.45)  # school means
  PU = runif(total, PUk[k] - 0.02, PUk[k] + 0.02)  # individual probabilities
  U = rbinom(total, size = 1, prob = PU)
  Ukbar = tapply(U, list(k), mean)[k]  # observed school averages
  Uc = U - Ukbar  # school-mean-centered U
  
  # Observed pre-treatment confounder X
  PXk = runif(K, 0.3, 0.5)  # school means
  PX = runif(total, PXk[k] - 0.02, PXk[k] + 0.02)  # individual probabilities
  X = rbinom(total, size = 1, prob = PX)
  Xkbar = tapply(X, list(k), mean)[k]  # observed school averages
  Xc = X - Xkbar  # school-mean-centered X
  
  # Post-treatment confounder V (continuous)
  t0k = rnorm(K, 0, 8)  # school-random-effects (z = 0)
  t1k = rnorm(K, 0, 6)  # school-random-effects (1 - 0)
  V0 = 35 + t0k[k] + 10*Xc + 20*Uc
  V1 = 40 + t0k[k] + t1k[k] + 10*Xc + 20*Uc
  V = Z*V1 + (1-Z)*V0
  
  # First-year treatment assignment D (binary)
  s0k = rnorm(K, 0, 1)  # school-random-effects (z = 0)
  s1k = rnorm(K, 0, 1)  # school-random-effects (z = 1)
  u0 = rlogis(total)
  u1 = rlogis(total)
  D0 = 1 * (u0 <= -Xc - Uc - 0.1*V0 + s0k[k])
  D1 = 1 * (u1 <= Xc + Uc + 0.05*V1 + s1k[k])
  D = Z*D1 + (1-Z)*D0
  
  # Potential outcome Y (continuous)
  g0k = rnorm(K, 0, 3)  # school-random-effects (intercept)
  gzk = rnorm(K, 0, 2)  # school-random-effects (z)
  gdk = rnorm(K, 0, 2)  # school-random-effects (d)
  gzdk = rnorm(K, 0, 1)  # school-random-effects (interaction)
  
  Y00 = 80 + g0k[k] + 40*Uc + 20*Xc + 0.2*V0
  Y01 = 95 + g0k[k] + gdk[k] + 40*Uc + 20*Xc + 0.2*V0
  Y10 = 90 + g0k[k] + gzk[k] + 40*Uc + 20*Xc + 0.2*V1
  Y11 = 100 + g0k[k] + gzk[k] + gdk[k] + gzdk[k] + 40*Uc + 20*Xc + 0.2*V1
  
  Y = Z*D*Y11 + (1-Z)*D*Y01 + Z*(1-D)*Y10 + (1-Z)*(1-D)*Y00 + rnorm(total, 0, 6)
  
  # Data output
  data = data.frame(k = k, i = i, U = U, Uc = Uc, X = X, Xc = Xc, Z = Z, V0 = V0, V1 = V1, V = V, D0 = D0, D1 = D1, D = D, Y00 = Y00, Y01 = Y01, Y10 = Y10, Y11 = Y11, Y = Y)
  
  return(data)
}

#### Analytic Methods for Estimation ####
# MS2T-IV Estimation Function
Estimate.MS2T = function(data, adjust.X){
  #[input] data: (data.frame) data set returned by 'gendata'
  #[param] adjust.X: (logical) whether to adjust for X and its interaction with Z in Stage-1 analysis
  
  # Stage 1
  # Use school-by-school OLS analysis
  K = max(data$k)  # number of schools
  alpha1k <- beta0k <- beta1k <- theta1k <- rep(0, K)
  for(s in 1:K){
    datas = data[which(data$k==s), ]
    if(adjust.X == T){
      Vmodel = lm(V ~ Z*Xc, data = datas)
      Dmodel = lm(D ~ Z*Xc, data = datas)
      Ymodel = lm(Y ~ Z*Xc, data = datas)
    }else{
      Vmodel = lm(V ~ Z, data = datas)
      Dmodel = lm(D ~ Z, data = datas)
      Ymodel = lm(Y ~ Z, data = datas)
    }
    alpha1k[s] = coef(Vmodel)[2]
    beta0k[s] = coef(Dmodel)[1]
    beta1k[s] = coef(Dmodel)[2]
    theta1k[s] = coef(Ymodel)[2]
  }
  
  # Stage 2
  beta2k = beta0k + beta1k
  model = lm(theta1k ~ beta1k + beta2k + alpha1k)
  ATE = sum(coef(model)*c(1,1,1,mean(alpha1k)))
  
  # Return Stage-2 parameter estimates and ATE
  return(c(coef(model), mean(alpha1k), ATE))
}

# Naive Estimation Function
Estimate.Naive = function(data, adjust.X){
  #[input] data: (data.frame) data set returned by 'gendata'
  #[param] adjust.X: (logical) whether to adjust for X
  
  # School-by-school naive analysis
  K = max(data$k)
  ATEk = rep(NA, K)
  for(s in 1:K){
    datas = data[which(data$k==s), ]
    if(adjust.X == T){
      index0 = which(datas$X == 0)
      index1 = which(datas$X == 1)
      # Subgroup X = 0
      if(length(index0) > 0){
        datas0 = datas[index0, ]
        group.mean0 = tapply(datas0$Y, list(datas0$Z, datas0$D), mean)
        if(nrow(group.mean0) == 2 & ncol(group.mean0) == 2){
          ATEk0 = group.mean0[2,2] - group.mean0[1,1]
        }else{
          ATEk0 = NA
        }
      }else{
        ATEk0 = NA
      }
      # Subgroup X = 1
      if(length(index1) > 0){
        datas1 = datas[index1, ]
        group.mean1 = tapply(datas1$Y, list(datas1$Z, datas1$D), mean)
        if(nrow(group.mean1) == 2 & ncol(group.mean1) == 2){
          ATEk1 = group.mean1[2,2] - group.mean1[1,1]
        }else{
          ATEk1 = NA
        }
      }else{
        ATEk1 = NA
      }
      # Expectation over X
      if(is.na(ATEk0) | is.na(ATEk1)){
        ATEk[s] = NA
      }else{
        ATEk[s] = (ATEk0*length(index0) + ATEk1*length(index1)) / (length(index0) + length(index1))
      }
    }else{
      group.mean = tapply(datas$Y, list(datas$Z, datas$D), mean)
      if(nrow(group.mean) == 2 & ncol(group.mean) == 2){
        ATEk[s] = group.mean[2,2] - group.mean[1,1]
      }else{
        ATEk[s] = NA
      }
    }
  }
  ATE = mean(ATEk, na.rm = T)
  
  return(ATE)
}

# IPTW Estimation Function
Estimate.IPTW = function(data, adjust.X){
  #[input] data: (data.frame) data set returned by 'gendata'
  #[param] adjust.X: (logical) whether to adjust for X when computing the weights
  
  # School-by-school analysis
  K = max(data$k)
  ATEk = rep(NA, K)
  for(s in 1:K){
    datas = data[which(data$k==s), ]
    datas$weight = NA
    datas$PD.ZXVk = NA
    datas$PD.Zk = NA
    # calculate conditional probabilities
    if(adjust.X == T){
      PD.fit = glm(D ~ Z + X + V, data = datas, family = binomial())
    }else{
      PD.fit = glm(D ~ Z + V, data = datas, family = binomial())
    }
    datas$PD.ZXVk = as.numeric(predict(PD.fit, type = "response"))
    for(z in c(0,1)){
      index0.Zk = which(datas$D==0 & datas$Z==z)
      index1.Zk = which(datas$D==1 & datas$Z==z)
      if(length(index0.Zk) > 0){datas$PD.Zk[index0.Zk] = length(index0.Zk) / (length(index0.Zk)+length(index1.Zk))}
      if(length(index1.Zk) > 0){datas$PD.Zk[index1.Zk] = length(index1.Zk) / (length(index0.Zk)+length(index1.Zk))}
    }
    # Assign weight Pr(D|Z)/Pr(D|Z,X,V) to each individual
    datas$weight = datas$PD.Zk / datas$PD.ZXVk
    ATEk[s] = sum(coef(lm(Y ~ Z*D, weights = weight, data = datas))[-1])
  }
  ATE = mean(ATEk, na.rm = T)
  
  return(ATE)
}

#### Simulations for Comparing Estimation Performance ####
# Computation Time: 5 mins for n = 500, K = 100, max(nk) = 1000
ncores = 10  # number of cores used for parallel computing
n = 500  # number of simulations
K = 100  # number of schools
nkset = c(30, 100, 1000, 5000)  # a set of numbers of students per school
adjust.X = T

results.MS2T <- results.Naive <- results.IPTW <- matrix(nrow = 6*length(nkset), ncol = 3)  # store the results in these matrices
rownames(results.MS2T) <- rownames(results.Naive) <- rownames(results.IPTW) <- rep(c('gamma1', 'gamma2', 'gamma3', 'thetaV', 'alpha1', 'ATE'), length(nkset))
colnames(results.MS2T) <- colnames(results.Naive) <- colnames(results.IPTW) <- c('bias', 'var', 'MSE')

for(j in 1:length(nkset)){
  nk = nkset[j]
  
  # Run MS2T-IV
  tic(paste0('Estimation completed. MS2T-IV #students = ', nkset[j]))
  # Prepare for parallel computing
  cl = makeSOCKcluster(ncores)
  registerDoSNOW(cl)
  # Create a progress bar
  pb = txtProgressBar(max = n, style = 3)
  progress = function(n){setTxtProgressBar(pb, n)}
  opts = list(progress = progress)
  
  est.MS2T = foreach(i = 1:n, .combine = rbind, .packages = c('MASS'), .options.snow = opts) %dopar% {
    set.seed(i)
    data = gendata(K, nk)
    Estimate.MS2T(data, adjust.X = adjust.X)
  }
  
  close(pb)
  stopCluster(cl)
  toc()
  
  # Run Naive
  tic(paste0('Estimation completed. Naive #students = ', nkset[j]))
  # Prepare for parallel computing
  cl = makeSOCKcluster(ncores)
  registerDoSNOW(cl)
  # Create a progress bar
  pb = txtProgressBar(max = n, style = 3)
  progress = function(n){setTxtProgressBar(pb, n)}
  opts = list(progress = progress)
  
  est.Naive = foreach(i = 1:n, .combine = c, .packages = c('MASS'), .options.snow = opts) %dopar% {
    set.seed(i)
    data = gendata(K, nk)
    Estimate.Naive(data, adjust.X = adjust.X)
  }
  
  close(pb)
  stopCluster(cl)
  toc()
  
  # Run IPTW
  tic(paste0('Estimation completed. IPTW #students = ', nkset[j]))
  # Prepare for parallel computing
  cl = makeSOCKcluster(ncores)
  registerDoSNOW(cl)
  # Create a progress bar
  pb = txtProgressBar(max = n, style = 3)
  progress = function(n){setTxtProgressBar(pb, n)}
  opts = list(progress = progress)
  
  est.IPTW = foreach(i = 1:n, .combine = c, .packages = c('MASS'), .options.snow = opts) %dopar% {
    set.seed(i)
    data = gendata(K, nk)
    Estimate.IPTW(data, adjust.X = adjust.X)
  }
  
  close(pb)
  stopCluster(cl)
  toc()
  
  ATE.MS2T = colMeans(est.MS2T)
  ATE.Naive = mean(est.Naive)
  ATE.IPTW = mean(est.IPTW)
  bias.MS2T = ATE.MS2T - c(10, 15, -5, 0.2, 5, 21)
  bias.Naive = ATE.Naive - 21
  bias.IPTW = ATE.IPTW - 21
  var.MS2T = apply(est.MS2T, 2, var)
  var.Naive = var(est.Naive)
  var.IPTW = var(est.IPTW)
  MSE.MS2T = bias.MS2T^2 + var.MS2T
  MSE.Naive = bias.Naive^2 + var.Naive
  MSE.IPTW = bias.IPTW^2 + var.IPTW

  results.MS2T[(1:6)+6*(j-1),1] = bias.MS2T
  results.MS2T[(1:6)+6*(j-1),2] = var.MS2T
  results.MS2T[(1:6)+6*(j-1),3] = MSE.MS2T
  results.Naive[6*j,1] = bias.Naive
  results.Naive[6*j,2] = var.Naive
  results.Naive[6*j,3] = MSE.Naive
  results.IPTW[6*j,1] = bias.IPTW
  results.IPTW[6*j,2] = var.IPTW
  results.IPTW[6*j,3] = MSE.IPTW
}

est.results = cbind(rep(K, nrow(results.MS2T)), rep(nkset, each = 6), results.MS2T, results.IPTW, results.Naive)
colnames(est.results) = c('K', 'nk', 'MS2T bias', 'MS2T var', 'MS2T MSE', 'IPTW bias', 'IPTW var', 'IPTW MSE', 'Naive bias', 'Naive var', 'Naive MSE')





#### Simulations for Comparing Inference Performance ####
# MS2T-IV Estimation Function for Bootstrap
Estimate.MS2T.boots = function(data, adjust.X = T){
  #[input] data: (data.frame) data set returned by 'gendata'
  #[param] adjust.X: (logical) whether to adjust for X and its interaction with Z in stage-1 analysis
  
  # Re-calculate school-mean-centered X
  data$Xkbar = tapply(data$X, list(data$k), mean)[data$k]  # school averages
  data$Xc = data$X - data$Xkbar  # school-mean-centered X
  
  # Stage 1
  # Use school-by-school OLS analysis
  K = max(data$k)
  alpha1k <- beta0k <- beta1k <- theta1k <- rep(0, K)
  for(s in 1:K){
    datas = data[which(data$k==s), ]
    if(adjust.X == T){
      Vmodel = lm(V ~ Z*Xc, data = datas)
      Dmodel = lm(D ~ Z*Xc, data = datas)
      Ymodel = lm(Y ~ Z*Xc, data = datas)
    }else{
      Vmodel = lm(V ~ Z, data = datas)
      Dmodel = lm(D ~ Z, data = datas)
      Ymodel = lm(Y ~ Z, data = datas)
    }
    alpha1k[s] = coef(Vmodel)[2]
    beta0k[s] = coef(Dmodel)[1]
    beta1k[s] = coef(Dmodel)[2]
    theta1k[s] = coef(Ymodel)[2]
  }
  
  # Stage 2
  beta2k = beta0k + beta1k
  model = lm(theta1k ~ beta1k + beta2k + alpha1k)
  ATE = sum(coef(model)*c(1,1,1,mean(alpha1k)))
  
  return(ATE)
}

# MS2T-IV Estimation Function for Improper (naive) CI and Zitzmann's Jackknife
Estimate.MS2T.impp = function(data, adjust.X = T){
  #[input] data: (data.frame) data set returned by 'gendata'
  #[param] adjust.X: (logical) whether to adjust for X and its interaction with Z in stage-1 analysis
  
  # Stage 1
  # Use school-by-school OLS analysis
  K = max(data$k)
  alpha1k <- beta0k <- beta1k <- theta1k <- rep(0, K)
  for(s in 1:K){
    datas = data[which(data$k==s), ]
    if(adjust.X == T){
      Vmodel = lm(V ~ Z*Xc, data = datas)
      Dmodel = lm(D ~ Z*Xc, data = datas)
      Ymodel = lm(Y ~ Z*Xc, data = datas)
    }else{
      Vmodel = lm(V ~ Z, data = datas)
      Dmodel = lm(D ~ Z, data = datas)
      Ymodel = lm(Y ~ Z, data = datas)
    }
    alpha1k[s] = coef(Vmodel)[2]
    beta0k[s] = coef(Dmodel)[1]
    beta1k[s] = coef(Dmodel)[2]
    theta1k[s] = coef(Ymodel)[2]
  }
  
  # Stage 2
  beta2k = beta0k + beta1k
  
  # Improper stage-2 variance
  model = lm(theta1k ~ beta1k + beta2k + alpha1k)
  ATE = sum(coef(model)*c(1,1,1,mean(alpha1k)))
  Var.Naive = t(c(1,1,1,mean(alpha1k)))%*%vcov(model)%*%c(1,1,1,mean(alpha1k))
  
  # Zitzmann Jackknife variance
  nout = K / 25  # the number of subgroups is 25
  ATE.ZJk = rep(NA, 25)
  for(j in 1:25){
    kout = (j-1)*nout + (1:nout)  # schools left out
    fit.sub = lm(theta1k[-kout] ~ beta1k[-kout] + beta2k[-kout] + alpha1k[-kout])
    ATE.ZJk[j] = sum(coef(fit.sub)*c(1,1,1,mean(alpha1k[-kout])))
  }
  Var.ZJk = (24/25) * sum((ATE.ZJk - mean(ATE.ZJk))^2)
  
  return(c(ATE, as.vector(Var.Naive), Var.ZJk))
}

# Bootstrap Function
Bootstrap = function(data, nboots){
  #[input] data: (data.frame) data set returned by 'gendata'
  #[param] nboots: (numeric) number of bootstrap samples
  
  bsamples = list()
  for(b in 1:nboots){
    sampleid = c()
    samplesch = c()
    # Bootstrap schools
    schsample = sample(1:max(data$k), replace = T)
    for(s in 1:length(schsample)){
      # Bootstrap students within each treatment group in each school
      samp = c(sample(which(data$k==schsample[s] & data$Z==0), replace = T), 
               sample(which(data$k==schsample[s] & data$Z==1), replace = T))
      sampleid = c(sampleid, samp)
      samplesch = c(samplesch, rep(s, length(samp)))
    }
    bsamples[[b]] = list(id = sampleid, sch = samplesch)
  }
  return(bsamples)
}

# Coverage rate Function
Cover = function(bATE, tATE = 21, level = 0.95){
  #[input] bATE: (numeric) bootstrapping ATE estimates
  #[param] tATE: (numeric) true ATE
  #[param] level: (numeric) confidence level (between 0 and 1)
  
  covered = (quantile(bATE, probs = (1-level)/2) <= tATE) & (quantile(bATE, probs = 1-(1-level)/2) >= tATE)
  return(covered)
}

# Computation Time: 15.6 hrs for n = 500, K = 100, max(nk) = 1000
ncores = 10  # number of cores used for parallel computing
n = 500  # number of simulations
nboots = 500  # number of bootstrap samples
K = 76  # number of schools
nkset = c(60)  # a set of numbers of students per school

inf.results = matrix(nrow = length(nkset), ncol = 8)  # store the results
rownames(inf.results) = nkset
colnames(inf.results) = c('avgATE', 'N-V', 'ZJk-V', 'N-CR',
                          'ZJk-CR', 'B-V', 'B-CR', 'BCa-CR')

for(j in 1:length(nkset)){
  nk = nkset[j]
  
  # Run MS2T-IV estimation and inference
  tic(paste0('Inference completed. #schools = ', K, '; #students = ', nkset[j]))
  # Prepare for parallel computing
  cl = makeSOCKcluster(ncores)
  registerDoSNOW(cl)
  # Create a progress bar
  pb = txtProgressBar(max = n, style = 3)
  progress = function(n){setTxtProgressBar(pb, n)}
  opts = list(progress = progress)
  
  results.inference = foreach(i = 1:n, .combine = rbind,
                              .packages = c('MASS', 'coxed'),
                              .options.snow = opts) %dopar% {
    # Data generation
    set.seed(i)
    data = gendata(K, nk)
    
    # Bootstrap
    bsamples = Bootstrap(data, nboots)
    bATE = rep(NA, nboots)
    for(b in 1:nboots){
      datab = data[bsamples[[b]]$id, ]  # bootstrap sample b
      datab$k = bsamples[[b]]$sch  # re-assign school id
      bATE[b] = Estimate.MS2T.boots(datab)
    }
    bvar = var(bATE)
    cover.boots = Cover(bATE)
    bBCaCI = bca(bATE)  # BCa confidence interval
    cover.BCa = (bBCaCI[1]<=21) & (bBCaCI[2]>=21)
    
    # Improper (naive) CI and Zitzmann's jackknife
    est = Estimate.MS2T.impp(data)
    cover.Naive = (est[1]-qnorm(0.975)*sqrt(est[2]) <= 21) & (est[1]+qnorm(0.975)*sqrt(est[2]) >= 21)
    cover.ZJk = (est[1]-qnorm(0.975)*sqrt(est[3]) <= 21) & (est[1]+qnorm(0.975)*sqrt(est[3]) >= 21)
    
    c(est, cover.Naive, cover.ZJk, bvar, cover.boots, cover.BCa)
  }
  
  close(pb)
  stopCluster(cl)
  toc()
  
  inf.results[j, ] = colMeans(results.inference)  # average over n simulations
}