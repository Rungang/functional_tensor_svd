library(rTensor)
source('FTSVD.R')
library(MASS)
library(fdapace)
library(latex2exp)
library('funData')
library('MFPCA')


## simulation 4: low-dimensional simulation that mimics multivariate FPCA
set.seed(100)
p1 = 100
Q = 10
n = 30
p2.set = c(2, 5, 20, 50)
tau.set = c(0, 0.05, 0.1, 0.2, 0.5, 1)
lambda = 2
iter = 100
r = 1

error_xi = array(NA, dim = c(length(p2.set), length(tau.set), iter))
error_xi_mfpca = error_xi
for (i in 1:iter){
  print(i)
  for (p.ind in 1:length(p2.set)) {
    for (tau.ind in 1:length(tau.set)) {
      p2 = p2.set[p.ind]
      tau = tau.set[tau.ind]
      a = as.matrix(rnorm(p1)); a = a / norm(a, type = 'F') 
      b = as.matrix(rnorm(p2)); b = b / norm(b, type = 'F')
      tn = round(sort(runif(n)),6)
      rxi = xi.generate(tn,Q)
      xin = rxi[[1]]
      xi = rxi[[2]]
      
      X = ttl(as.tensor(array(lambda,dim=c(1,1,1))), c(list(a), list(b), list(xin)), 1:3)
      Z = array(0,dim = c(p1,p2,n))
      for (j1 in 1:p1) {
        for (j2 in 1:p2) {
          Z[j1,j2,] = xi.generate(tn,Q)[[1]]
        }
      }
      Y = X + as.tensor(tau * array(rnorm(p1*p2*n),dim=c(p1,p2,n)))
      # estimate largest singular value
      a.est = svd(k_unfold(Y,1)@data)$u[,1]
      Y2 = ttm(Y, t(as.matrix(a.est)), 1)
      lambda_lower = svd(k_unfold(Y2,2)@data)$d[1] /sqrt(n)
      ftsvd.est = FTSVD(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=1, alpha=200*lambda_lower)
      error_xi[p.ind, tau.ind, i] = sqrt(1 - abs(t(ftsvd.est[[4]])%*%xi)^2)
      
      # mfsvd
      Y_array = Y@data
      mf_list = NULL
      for (j in 1:p2){
        temp = funData(argvals = tn, X = Y_array[,j,])
        mf_list = c(mf_list, list(temp))
      }
      mfd = multiFunData(mf_list)
      uniExpansions <- rep(list(list(type = "uFPCA", npc = 3)),p2)
      
      tryCatch(
        expr = {
          mfpca =  MFPCA(mfd, M = r, uniExpansions = uniExpansions)
          Lt = rep(list(tn), p2 * r)
          Ly = rep(list(NULL), p2 * r)
          count = 1
          for (j2 in 1:p2) {
            for (k in 1:r){
              Ly[[count]] = mfpca$functions[[j2]]@X[k,]
              count = count + 1
            }
          }
          mu = c(list("t"=tn, "mu"=rep(0,n)))
          mfPCA.est = FPCA(Ly, Lt,
                           list(dataType='Dense', maxK = r, kernel='epan', userMu = mu, nRegGrid=101, useBinnedData = 'OFF'))
          mfpca.est = mfPCA.est$phi / norm(mfPCA.est$phi, type = 'F') * sqrt(r)
          error_xi_mfpca[p.ind, tau.ind, i] = sqrt(1 - abs(t(mfpca.est)%*%xi)^2)
        },
        error = function(e){
        }
      )
      
    }
  }
}

mean_xi = apply(error_xi, c(1,2), mean)
mean_xi_mfpca = apply(error_xi_mfpca, c(1,2), mean)

save(error_xi, error_xi_mfpca, file = "sim4.RData")
