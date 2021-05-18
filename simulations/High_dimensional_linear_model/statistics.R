# columnwise standardization
colStand   = function(m) {
  m_std = m
  len = sqrt(colSums(m^2))
  m_std[, which(len!=0)] = t(t(m[, which(len!=0)])/len[len!=0])
  m_std
}

# proportion of the correctly identified varying coefficients.
identifyVC = function(m, set_VC) {
  sum(apply(m[, set_VC], 2, "var") != 0)/length(set_VC)
}

# proportion of the correctly identified constant coefficients.
identifyCC = function(m, set_CC) {
  sum((apply(m[, set_CC], 2, "var") == 0) & 
        (m[1, set_CC] != 0))/length(set_CC)
}

# proportion of the correctly identified zero coefficients.
identifyZC = function(m, set_ZC) {
  sum((apply(m[, set_ZC], 2, "var") == 0) & 
        (m[1, set_ZC] == 0))/length(set_ZC)
}

# classification error between zero and nonzero groups.
classifyER = function(m, set_NZ) {
  NZ = NZ_hat = numeric(ncol(m))
  NZ[set_NZ] = 1
  
  set_NZ_hat = which(colSums(abs(m)) !=0)
  NZ_hat[set_NZ_hat] = 1
  
  mean(abs(NZ - NZ_hat))
}

# Matthews Correlation Coefficient
MCC_vc = function(m, Tset, Fset) {
  beta_hat = (apply(m, 2, "var") != 0)
  TP = sum(beta_hat[Tset])
  FN = length(Tset) - TP
  FP = sum(beta_hat[Fset])
  TN = length(Fset) - FP

  MCC = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  if (MCC != 0) { MCC = (TP*TN - FP*FN)/sqrt(MCC) }
  
  MCC
}

MCC_nz = function(m, Tset, Fset) {
  beta_hat = vector("numeric",ncol(m))
  for (g in 1:ncol(m))  
    beta_hat[g] = 1 - all(m[,g] == 0) 
  TP = sum(beta_hat[Tset])
  FN = length(Tset) - TP
  FP = sum(beta_hat[Fset])
  TN = length(Fset) - FP

  MCC = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  if (MCC != 0) { MCC = (TP*TN - FP*FN)/sqrt(MCC) }
  
  MCC
}

identifyall = function(m, p) {
  var = 0
  const = 0
  zero = 0
  for (i in 1:p) {
    if (var(m[,i]) != 0) {
      var = var + 1
    } else {
      if (all(m[,i] == 0)) {
        zero = zero+1
      } else {
        const = const + 1
      }
    }
  }
  c(var, const, zero)
}

# pick which
thetapick = function(m, group, opt, p) {
  thetapick <- rep(0, 2*p)
  for (i in 1:(2*p)) {
    idx = (group==i)
    thetapick[i] = 1-all(m[idx, opt] == 0)
  }
  thetapick
}

assig = function(n_args) {
  # example: n_args = c(2, 3, 4)
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
  t(expand.grid(cargs))
}

fits = function(case, vc, cc, n,  p,  rho, error, tran) {
  set_VC <- 1:vc
  set_CC <- (vc+1):(vc+cc)
  set_ZC <- (vc+cc+1):p
  set_NZ <- 1:(vc+cc)
  set_C  <- (vc+1):p
  B      <- length(case)
  res0 = readRDS(case[1])
  lfile  <- length(res0) 
  res0 = res0[[1]]$fit
  nmodel = length(res0)
  #numberbic = length(res0[[4]]$optimal)
  
  group <- sort(c(2*(1:p)-1, rep(2*(1:p), each = 6)))
  IMSE = VC = CC = ZC = ER = Mvc = Mnz= matrix(0, B*lfile, nmodel)
  idtfall = matrix(0, B*lfile, nmodel*3)
  AUC = vector("list", nmodel)
  for (j in 1:nmodel) AUC[[j]] = matrix(0, B*lfile, 1)
  AUC_testing = vector("list", nmodel)
  for (j in 1:nmodel) AUC_testing[[j]] = matrix(0, B*lfile, 3)
  
  for (j in 0:(B-1)) {
    res0 = readRDS(case[j+1])
    for (l in 1:lfile) {
      res = res0[[l]]
      set.seed(res$seed)
      res$b.std = colStand( generator(n,  p,  rho, error, tran)$beta )
      Beta.std = vector("list", nmodel)
      for (i in 1:nmodel) {
        Beta.std[[i]] = colStand(res$fit[[i]]$Beta)
        #peason and kendall tau correlation coefficient
        AUC[[i]][j*lfile+l,] = res$fit[[i]]$AUC
        AUC_testing[[i]][j*lfile+l,] = res$fit[[i]]$testing_AUC
        # integrated MSE.
        IMSE[j*lfile+l, i] <- sum(colMeans((Beta.std[[i]] - res$b.std)^2))
        # Varying coefficients' correctly identification rate.
        VC[j*lfile+l,   i] <- identifyVC(Beta.std[[i]], set_VC) 
        # Constant coefficients' correctly identification rate.
        CC[j*lfile+l, i]   <- identifyCC(Beta.std[[i]], set_CC)
        # Zero coefficients' correclty identification rate
        ZC[j*lfile+l, i]   <- identifyZC(Beta.std[[i]], set_ZC) 
        # classification error between zero and nonzero groups.
        ER[j*lfile+l, i]   <- classifyER(Beta.std[[i]], set_NZ)
        Mvc[j*lfile+l, i] <- MCC_vc(Beta.std[[i]], set_VC, set_C)
        Mnz[j*lfile+l, i] <- MCC_nz(Beta.std[[i]], set_NZ, set_ZC)
        
        idtfall[j*lfile+l,(3*i-2):(3*i)] <- identifyall(Beta.std[[i]], p)
      }
    }
  }
  
  AUC = cbind(AUC[[1]], AUC[[2]], AUC[[3]], AUC[[4]])
  AUC_testing = cbind(AUC_testing[[1]], AUC_testing[[2]], AUC_testing[[3]], AUC_testing[[4]])
  
  result = matrix(0, 26, nmodel)
  result[1,]  = round(colMeans(VC), 4)
  result[2,]  = round(apply(VC, 2, sd), 4)
  result[3,]  = round(colMeans(CC), 4)
  result[4,]  = round(apply(CC, 2, sd), 4)
  result[5,]  = round(colMeans(ZC), 4)
  result[6,]  = round(apply(ZC, 2, sd), 4)
  result[7,]  = round(colMeans(ER), 4)
  result[8,]  = round(apply(ER, 2, sd), 4)
  result[9,]  = round(colMeans(IMSE), 4)
  result[10,] = round(apply(IMSE, 2, sd), 4)
  result[11,] = round(colMeans(Mvc), 4)
  result[12,] = round(apply(Mvc, 2, sd), 4)
  result[13,] = round(colMeans(Mnz), 4)
  result[14,] = round(apply(Mnz, 2, sd), 4)
  result[15,] = round(colMeans(AUC[,c(1,4,7,10)]), 4)
  result[16,] = round(apply(AUC[,c(1,4,7,10)], 2, sd), 4)
  result[17,] = round(colMeans(AUC[,c(2,5,8,11)]), 4)
  result[18,] = round(apply(AUC[,c(2,5,8,11)], 2, sd), 4)
  result[19,] = round(colMeans(AUC[,c(3,6,9,12)]), 4)
  result[20,] = round(apply(AUC[,c(3,6,9,12)], 2, sd), 4)
  result[21,] = round(colMeans(AUC_testing[,c(1,4,7,10)]), 4)
  result[22,] = round(apply(AUC_testing[,c(1,4,7,10)], 2, sd), 4)
  result[23,] = round(colMeans(AUC_testing[,c(2,5,8,11)]), 4)
  result[24,] = round(apply(AUC_testing[,c(2,5,8,11)], 2, sd), 4)
  result[25,] = round(colMeans(AUC_testing[,c(3,6,9,12)]), 4)
  result[26,] = round(apply(AUC_testing[,c(3,6,9,12)], 2, sd), 4)
  
  result2 = colSums(apply(Mvc, 2, is.na))
  result3 = colSums(apply(Mnz, 2, is.na))
  result4 = colMeans(idtfall)
  
  list(result, result2, result3, result4)
}


generator = function(n, p, rho, error, tran, censor.rate = NULL) {
  sigma = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y)))
  
  d=5
  x = rmvnorm(n, mean = rep(0,p), sigma = sigma)
  u = runif(n)
  b = cbind(5*sin(2*pi*u), 5*cos(2*pi*u), 5*exp(2*u-1)-5, 5, -5, matrix(0, nrow=n, ncol=p-d))
  
  # linear predictor
  eta = numeric(n)
  for (i in 1:n) eta[i] = x[i, 1:d] %*% b[i, 1:d]
  
  # error   ==================================================
  if(error == "norm") e = rnorm(n, 0, 1)
  if(error =="contaminate") e = c(0.7*rnorm(n)+0.3*rcauchy(n))
  if(substr(error, 1, 1) == "t") e = rt(n, as.numeric(substr(error,2,2)))
  
  # transformation ===========================================
  if(tran == "lm") g = function(x) x
  
  y0 = g(eta + e)
  
  # censor  ==================================================
  y = y0
  status = rep(1, n)
  list(x=x, u=u, y=y, beta=b, status=status)                        
}















