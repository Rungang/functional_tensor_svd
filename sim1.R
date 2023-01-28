library(rTensor)
source('FTSVD.R')
library(MASS)
library(fdapace)
library(latex2exp)
library('funData')
library('MFPCA')

######################################################
# Simulation 1: rank-1 functional tensor decomposition with increasing sampling 
# study the asymptotics and only compare with CP power iteration.
######################################################

set.seed(1)
Q = 10
p.set = c(20,30,50)
n.set = c(10,30,50,100,200)
sigma = 0
tau = 1
lambda = 2
iter = 100
error_a = array(0, dim = c(length(p.set), length(n.set), iter))
error_xi = array(0, dim = c(length(p.set), length(n.set), iter))
error_a_CP = array(0, dim = c(length(p.set), length(n.set), iter))
error_xi_CP = array(0, dim = c(length(p.set), length(n.set), iter))
error_xi_FPCA = array(0, dim = c(length(p.set), length(n.set), iter))
error_a_FPCA = error_xi_FPCA
for (i in 1:iter){
  print(i)
  for (p.ind in 1:length(p.set)) {
    for (n.ind in 1:length(n.set)) {
      p = p.set[p.ind]
      n = n.set[n.ind]
      # generate random loadings
      a = as.matrix(rnorm(p)); a = a / norm(a, type = 'F') 
      b = as.matrix(rnorm(p)); b = b / norm(b, type = 'F')
      tn = round(sort(runif(n)),8)
      rxi = xi.generate(tn,Q)
      xin = rxi[[1]]
      xi = rxi[[2]]
      
      X = ttl(as.tensor(array(lambda,dim=c(1,1,1))), c(list(a), list(b), list(xin)), 1:3)
      Y = X + as.tensor(tau * array(rnorm(p^2*n),dim=c(p,p,n)))
      
      # estimate largest singular value
      a.est = svd(k_unfold(Y,1)@data)$u[,1]
      Y2 = ttm(Y, t(as.matrix(a.est)), 1)
      lambda_lower = svd(k_unfold(Y2,2)@data)$d[1] /sqrt(n)
      #print(lambda_est)
      ftsvd.est = FTSVD(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=1, alpha = 200 * lambda_lower)
      error_a[p.ind, n.ind, i] = sqrt(1 - abs(t(ftsvd.est[[2]])%*%a)^2) 
      error_xi[p.ind, n.ind, i] = sqrt(1 - abs(t(ftsvd.est[[4]])%*%xi)^2)
      
      CP.est = CP(f.tensor = Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=1)
      error_a_CP[p.ind, n.ind, i] = sqrt(1 - abs(t(CP.est[[2]])%*%a)^2) 
      error_xi_CP[p.ind, n.ind, i] = sqrt(1 - abs(t(CP.est[[4]])%*%xi)^2)
      
      Lt = rep(list(tn), p^2)
      Ly = rep(list(NULL), p^2)
      count = 1
      for (j1 in 1:p) {
        for (j2 in 1:p) {
          Ly[[count]] = Y@data[j1,j2,]
          count = count + 1
        }
      }
      mu = c(list("t"=tn, "mu"=rep(0,n)))
      fPCA.uni = FPCA(Ly, Lt,
                      list(dataType='Dense', maxK = 1, kernel='epan', userMu = mu, nRegGrid=101))
      fpca.est = fPCA.uni$phi / norm(fPCA.uni$phi, type = 'F')
      error_xi_FPCA[p.ind, n.ind, i] = sqrt(1 - abs(t(fpca.est)%*%xi)^2)
      
    }
    
  }
  
}

save(error_a, error_xi, error_a_CP, error_xi_CP, error_xi_FPCA, error_a_FPCA, n.set, iter, file = "sim1.RData")
