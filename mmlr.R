# Multivariate  Multiple-linear regression 
f_mmlr = function(Y, Xk, mytype = 'lsr', stdpred = F, tregG = 0, ridgedelta = 0){
  # Number of predictors and sample responses
  if(is.null(dim(Xk))){
    p = 1
    n = length(Xk)
  }
  else{
    #There are multiple predictors for multi-variate linear regression
    p = dim(Xk)[2]
    n = dim(Xk)[1]
  }
  
  v1s = rep(1,n)
  X = cbind(v1s, Xk)
  
  if(is.null(dim(Y))){
    d = 1
  }
  else{
    d = dim(Y)[2]
  }
  
  # For Variable standardization
  if (stdpred == T){
    for(k in 1:p){
      X[,k+1] = (X[,k+1] -mean(X[,k+1]))/sd(X[,k+1])
    }
  }
  
  ###
  
  myvif <- c()
  if (p == 1){
    myvif <- append(myvif, 1)
  }
  else{
    for(k in 1:p){
      Xp <- X[,-c(k+1)]
      Yp <- X[, k+1]
      
      Sp <- svd(t(Xp)%*%Xp)
      Up <- Sp$u
      Dp <- diag(Sp$d)
      Vp <- Sp$v
      
      betahatp <- Vp%*%solve(Dp)%*%t(Up)%*%t(Xp)%*%Yp
      yhatp <- Xp%*%betahatp
      
      SSEp <- sum((Yp - yhatp)^2)
      SSTp <- var(Yp)*(n-1)
      r2p <- 1- SSEp/SSTp
      
      myvif <- append(myvif, 1/(1-r2p))        
    }
    
    
  }
  
  #Multicoliniarity indicies
  mci <- sqrt(myvif)
  
  # Singular Value Decomposition
  #Least Squares
  if (mytype == 'lsr'){
    S = svd(t(X)%*%X)
  }
  #tikhonov 
  if (mytype == 'treg'){
    S = svd(t(X)%*%X + t(tregG)%*%tregG)
  }
  
  # Ridge
  if (mytype == 'ridge'){
    S = svd(t(X)%*%X +ridgedelta^2*diag(p+1))
  }
  
  #Least Absolute Differences
  U = S$u
  D = diag(S$d)
  V = S$v
  
  # (UDT)
  Xtxinv = V%*%solve(D)%*%t(U)
  betahat = Xtxinv%*%t(X)%*%Y
  H = X%*%Xtxinv%*%t(X)
  lev = diag(H)
  
  sqdetx = sqrt(det(t(X)%*%X)) 
  
  # Condition number
  sqkappa = sqrt(max(S$d)/min(S$d)) # Large is bad
  
  Yhat = X%*%betahat
  resid = Y-Yhat
  
  
  # Primary Metrics
  if (d == 1) {
    SSE = sum((resid^2))
    SST = var(Y)*(n-1)
    SSM = SST -SSE
  }
  if (d > 1) {
    mat1 = matrix(1, n, n)
    SST = t(Y)%*%(diag(n) - (1/n)*mat1)%*%Y
    SSE = t(Y)%*%(diag(n) - H)%*%Y
    SSM = t(Y)%*%(H - (1/n)*mat1)%*%Y
  }
  
  MST = SST/(d*(n-1))
  MSM = SSM/(d*p)
  MSE = SSE/(d*(n-p-1))
  sig = sqrt(MSE/d)
  
  # Secondary Metrics
  if (d == 1) {
    r2 = 1 - SSE/SST
    R2adj = 1 - MSE/MST 
    SEbhat = sig*sqrt(diag(Xtxinv))
    sresid = resid/(sig*(1-lev))
    bhats = betahat/SEbhat
  }
  if (d > 1) {
    r2 = 1 - sum(diag(SSE))/sum(diag(SST))
    R2adj = 1 - sum(diag(MSE))/sum(diag(MST))
    SEbhat = as.matrix(sqrt(diag(kronecker(d*MSE, solve(t(X)%*%X)))), norw = p+1, ncol = d)
    sresid = NULL
    bhats = NULL
  }
  
  # 3rd(Tertiary) Metrics
  if (d == 1) {
    Fstat = MSM/MSE
    pval = pf(Fstat, p, n-p-1, lower.tail = F)
  }
  if (d > 1) {
    # 5 stat g-1 = p, g = (p + 1), MSE= n-g = n- 9 -1 , MSM = g-1 = p, MST = n-1,
    roy_stat     = max(eigen(solve(MSE) %*% MSM)$values)                      #F dist with df, d, n-g-d-1
    pillai_stat  = sum(diag(solve(MST) %*% MSM))                              #F dist with df, d(g-1), n-g-1
    hl_stat      = sum(diag(solve(MSE) %*% MSM))                              #F dist with df, d(g-1), n-g-d-1
    wilks_stat   = ( (n - p - 1 ) - 1/2 * ( d + 2 ) ) * log(det(MST)/det(MSE) )    #chisq dist with df, d(g-1)
    m_wilks_stat = ( (n - p - d )/ (d*(p)) ) * log(det(MST)/det(MSE) )   #F dist with df, d(g-1), n-g-d-1
    
    #Stat in order Roy, Pillai, Hotelling-Lowry, Wilks, Modified-Wilks
    #Fstat = c(roy_stat, pillai_stat, hl_stat, wilks_stat, m_wilks_stat)
    
    # 5 pval if p = g-1
    pval_roy     = pf(roy_stat,       d,   d*(n-p-1),  lower.tail = FALSE)
    pval_pillai  = pf(pillai_stat,    d*p, d*(n-1),    lower.tail = FALSE)
    pval_hl      = pf(hl_stat,        d*p, d*(n-p-1),  lower.tail = FALSE)
    pval_wilks   = pchisq(wilks_stat, d*p,             lower.tail = FALSE)
    pval_m_wilks = pf(m_wilks_stat,   d*p, d*(n-p-1),  lower.tail = FALSE)
    
    #pval in order Roy, Pillai, Hotelling-Lowry, Wilks, Modified-Wilks
    pval = rbind(c('Roy', 'Pillai', 'HL', 'Wilks', 'Modified-Wilks'),
                 c(pval_roy, pval_pillai, pval_hl, pval_wilks, pval_m_wilks))
  }
  

  results =  list("pred" = Yhat, 
                  "residuals" = resid,
                  "SRES" = sresid,
                  "detx"= sqdetx,
                  "lev"= lev,
                  "Fstat" = Fstat,
                  "pval"= pval,
                  "bhat" = betahat,
                  "sebhat" = SEbhat,
                  "r2" = r2,
                  "r2adj"= R2adj, 
                  "sig" = sig,
                  "mse" = MSE,
                  "mst" = MST,
                  "msm" = MSM,
                  "sse" = SSE,
                  "sst" = SST,
                  "ssm" = SSM,
                  "mci" = mci,
                  "sqrtkappa" = sqkappa#,
                  #"roy_stat" = roy_stat,
                  #"pillai_stat" = pillai_stat,
                  #"hl_stat" = hl_stat,
                  #"wilks_stat" = wilks_stat,
                  #"m_wilks_stat" = m_wilks_stat
  )
  
  
  #Multicolineairryt indices
  # Modified regulrization (lasso, )
  return(results)
}