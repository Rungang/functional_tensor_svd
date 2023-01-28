library(rTensor)
source('FTSVD.R')
library(MASS)
library(fdapace)
library(latex2exp)
library('funData')
library('MFPCA')

######################################################
# Simulation 3: rank-r functional tensor decomposition with imbalanced dimensions.
######################################################

set.seed(1)
Q = 10
n = 30
sigma = 0
tau = 1
lambda = 8
r.set = c(1,2,1,3)
p2.set = c(500, 500, 500, 100)
p1.set = c(20, 20, 100, 100)

iter = 100
error_xi = array(0, dim = c(4, 4, iter))
error_a = array(0, dim = c(4, 4, iter))
error_b = array(0, dim = c(4, 4, iter))
exe_time = array(0, dim = c(4, 4, iter))

for (i in 1:iter){
  print(i)
  for (ind in 1:length(r.set)) {
    r = r.set[ind]
    p1 = p1.set[ind]
    p2 = p2.set[ind]
    X = 0
    A = NULL; B = NULL; XI = NULL
    for (s in 1:r) {
      a = as.matrix(rnorm(p1)); a = a / norm(a, type = 'F') 
      b = as.matrix(rnorm(p2)); b = b / norm(b, type = 'F')
      tn = round(sort(runif(n)),6)
      rxi = xi.generate(tn,Q)
      xin = rxi[[1]]
      xi = rxi[[2]]
      A = cbind(A, a); B = cbind(B, b); XI = cbind(XI, xi)
      X = X + ttl(as.tensor(array(lambda*(r-s+1),dim=c(1,1,1))), c(list(a), list(b), list(xin)), 1:3)
    }
    
    Z = array(0,dim = c(p1,p2,n))
    for (j1 in 1:p1) {
      for (j2 in 1:p2) {
        Z[j1,j2,] = xi.generate(tn,Q)[[1]]
      }
    }
    Y = X + sigma * Z + as.tensor(tau * array(rnorm(p1*p2*n),dim=c(p1,p2,n)))
    
    ptm1 = Sys.time()
    # estimate largest singular value
    a.est = svd(k_unfold(Y,1)@data)$u[,1]
    Y2 = ttm(Y, t(as.matrix(a.est)), 1)
    lambda_lower = svd(k_unfold(Y2,2)@data)$d[1] /sqrt(n)
    ftsvd.est = FTSVD(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=r, alpha = 200*lambda_lower)
    ptm2 = Sys.time()
    CP.est = CP(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=r, L=20)
    ptm3 = Sys.time()
    exe_time[ind, 1, i] = ptm2 - ptm1
    exe_time[ind, 2, i] = ptm3 - ptm2
    
    Lt = rep(list(tn), p1*p2)
    Ly = rep(list(NULL), p1*p2)
    count = 1
    for (j1 in 1:p1) {
      for (j2 in 1:p2) {
        Ly[[count]] = Y@data[j1,j2,]
        count = count + 1
      }
    }
    mu = c(list("t"=tn, "mu"=rep(0,n)))
    
    ptm1 = Sys.time()
    fPCA.uni = FPCA(Ly, Lt, 
                    list(dataType='Dense', maxK = r, kernel='epan', userMu = mu, nRegGrid=101, useBinnedData = 'OFF'))
    ptm2 = Sys.time()
    exe_time[ind, 3, i] = ptm2 - ptm1
    fpca.est = fPCA.uni$phi / norm(fPCA.uni$phi, type = 'F') * sqrt(r)
    
    error_a[ind, 1, i] = dist.measure(ftsvd.est[[2]], A)
    error_b[ind, 1, i] = dist.measure(ftsvd.est[[3]], B)
    error_xi[ind, 1, i] = dist.measure(ftsvd.est[[4]], XI)
    
    error_a[ind, 2, i] = dist.measure(CP.est[[2]], A)
    error_b[ind, 2, i] = dist.measure(CP.est[[3]], B)
    error_xi[ind, 2, i] = dist.measure(CP.est[[4]], XI)
    
    error_xi[ind, 3, i] = dist.measure(fpca.est, XI)
    # error_xi[ind, 4, i] = dist.measure(mfpca.est, XI)
    
    
  }
  
}

mean_a = apply(error_a, c(1,2), mean)
mean_b = apply(error_b, c(1,2), mean)
mean_xi = apply(error_xi, c(1,2), mean)
mean_time = apply(exe_time, c(1,2), mean)

sd_a = apply(error_a, c(1,2), sd)
sd_b = apply(error_b, c(1,2), sd)
sd_xi = apply(error_xi, c(1,2), sd)
sd_time = apply(exe_time, c(1,2), sd)


save(mean_a, mean_b, mean_xi, mean_time,
     sd_a, sd_b, sd_xi, sd_time,
     file = "sim3.RData")