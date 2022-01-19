library(rTensor)
source('FTSVD.R')
library(MASS)
library(fdapace)
library(latex2exp)

######################################################
# Simulation 1: rank-1 functional tensor decomposition with increasing sampling 
# study the asymtotics and only compare with CP power iteration.
######################################################

set.seed(100)
Q = 10
p.set = c(20,30,50)
n.set = c(10,30,50,100,200,400,800)
sigma = 1
tau = 1
lambda = 20
iter = 100
error_a = array(0, dim = c(length(p.set), length(n.set), iter))
error_xi = array(0, dim = c(length(p.set), length(n.set), iter))
error_a_CP = array(0, dim = c(length(p.set), length(n.set), iter))
error_xi_CP = array(0, dim = c(length(p.set), length(n.set), iter))


for (i in 1:iter){
  print(i)
  for (p.ind in 1:length(p.set)) {
    for (n.ind in 1:length(n.set)) {
      p = p.set[p.ind]
      n = n.set[n.ind]
      # generate random loadings
      a = as.matrix(rnorm(p)); a = a / norm(a, type = 'F') 
      b = as.matrix(rnorm(p)); b = b / norm(b, type = 'F')
      tn = sort(runif(n))
      rxi = xi.generate(tn,Q)
      xin = rxi[[1]]
      xi = rxi[[2]]
      
      X = ttl(as.tensor(array(lambda,dim=c(1,1,1))), c(list(a), list(b), list(xin)), 1:3)
      Y = X + as.tensor(tau * array(rnorm(p^2*n),dim=c(p,p,n)))
      
      ftsvd.est = FTSVD(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=1, penalty = 1e-5)
      error_a[p.ind, n.ind, i] = sqrt(1 - abs(t(ftsvd.est[[2]])%*%a)^2) 
      error_xi[p.ind, n.ind, i] = sqrt(1 - abs(t(ftsvd.est[[4]])%*%xi)^2)
      
      CP.est = CP(f.tensor = Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=1)
      error_a_CP[p.ind, n.ind, i] = sqrt(1 - abs(t(CP.est[[2]])%*%a)^2) 
      error_xi_CP[p.ind, n.ind, i] = sqrt(1 - abs(t(CP.est[[4]])%*%xi)^2)
      
    }
    
  }
  
}

save(error_a, error_xi, error_a_CP, error_xi_CP, file = "sim1.RData")

error_a.mean = apply(error_a, c(1,2), mean)
error_a_CP.mean = apply(error_a_CP, c(1,2), mean)
error_xi.mean = apply(error_xi, c(1,2), mean)
error_xi_CP.mean = apply(error_xi_CP, c(1,2), mean)

error_a.sd = apply(error_a, c(1,2), sd)
error_a_CP.sd = apply(error_a_CP, c(1,2), sd)
error_xi.sd = apply(error_xi, c(1,2), sd)
error_xi_CP.sd = apply(error_xi_CP, c(1,2), sd)

p.set = c(" p = 20", " p = 30", " p = 50")
data.plot = data.frame(error = c(as.vector(error_xi), as.vector(error_xi_CP), 
                                 as.vector(error_a), as.vector(error_a_CP)),
                       mode = rep(c("Functional","Tabular"), each = length(p.set)*length(n.set)*2*iter),
                       method = rep(rep(c("FTSVD","CP"),each = length(p.set)*length(n.set)*iter),2),
                       p = rep(as.factor(p.set),length(n.set)*4*iter),
                       n = as.factor(rep(rep(n.set, each = length(p.set)), 4*iter)))
data.plot$mode = as.factor(data.plot$mode)
data.plot$method = as.factor(data.plot$method)
library(ggplot2)

p1 = ggplot(data=data.plot, aes(x=n, y=error, fill = method)) + geom_boxplot() +
  facet_wrap(~mode:p) + theme(text=element_text(size=16)) + xlab("Grid density") + ylab("Estimation error")

p1


######################################################
# Simulation 2: rank-1 functional tensor decomposition with increasing perturbation 
# study the asymtotics and compare with CP power iteration and FDA.
######################################################

set.seed(100)
Q = 10
p.set = c(20,50)
n = 20
sigma.set = seq(0.2,1,0.2)
tau = 1
lambda = 20
iter = 100
error_xi = array(0, dim = c(length(p.set), length(sigma.set), iter))
error_xi_CP = array(0, dim = c(length(p.set), length(sigma.set), iter))
error_xi_FPCA = array(0, dim = c(length(p.set), length(sigma.set), iter))


for (i in 1:iter){
  print(i)
  for (p.ind in 1:length(p.set)) {
    for (sigma.ind in 1:length(sigma.set)) {
      p = p.set[p.ind]
      sigma = 1/sigma.set[sigma.ind]
      # generate random loadings
      a = as.matrix(rnorm(p)); a = a / norm(a, type = 'F') 
      b = as.matrix(rnorm(p)); b = b / norm(b, type = 'F')
      tn = round(sort(runif(n)),6)
      rxi = xi.generate(tn,Q)
      xin = rxi[[1]]
      xi = rxi[[2]]
      
      X = ttl(as.tensor(array(lambda,dim=c(1,1,1))), c(list(a), list(b), list(xin)), 1:3)
      Z = array(0,dim = c(p,p,n))
      for (j1 in 1:p) {
        for (j2 in 1:p) {
          Z[j1,j2,] = xi.generate(tn,Q)[[1]]
        }
      }
      Y = X + sigma * Z #+ as.tensor(tau * array(rnorm(p^2*n),dim=c(p,p,n)))
      
      ftsvd.est = FTSVD(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=1, penalty = 1e-5)
      error_xi[p.ind, sigma.ind, i] = sqrt(1 - abs(t(ftsvd.est[[4]])%*%xi)^2)
      
      CP.est = CP(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=1)
      error_xi_CP[p.ind, sigma.ind, i] = sqrt(1 - abs(t(CP.est[[4]])%*%xi)^2)
      
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
      error_xi_FPCA[p.ind, sigma.ind, i] = sqrt(1 - abs(t(fpca.est)%*%xi)^2)
      
    }
    
  }
  
}

save(error_a, error_xi, error_a_CP, error_xi_CP, file = "sim2.RData")


p.set = c(" p = 20", " p = 50")
data.plot = data.frame(error = c(as.vector(error_xi), as.vector(error_xi_CP), as.vector(error_xi_FPCA)),
                       method = as.factor(rep(rep(c("FTSVD","CP","FPCA"),each = length(p.set)*length(sigma.set)*iter))),
                       p = rep(as.factor(p.set),length(sigma.set)*3*iter),
                       n = as.factor(rep(rep(sigma.set, each = length(p.set)), 3*iter)))
data.plot$method = factor(data.plot$method, levels=levels(data.plot$method)[c(1,3,2)])
library(ggplot2)

p2 = ggplot(data=data.plot, aes(x=n, y=error, fill = method)) + geom_boxplot() +
  facet_wrap(~p) + theme(text=element_text(size=16)) + xlab(TeX("Inverse of $\\sigma$")) + ylab("Estimation error")

p2



######################################################
# Simulation 3: rank-r functional tensor decomposition with increasing L.
# study several simple settings.
######################################################

set.seed(1)
Q = 10
n = 30
sigma = 1
tau = 1
lambda = 80
r.set = c(1,2,1,3)
p2.set = c(500, 500, 500, 100)
p1.set = c(20, 20, 100, 100)

iter = 100
error_xi = array(0, dim = c(4, 3, iter))
error_a = array(0, dim = c(4, 3, iter))
error_b = array(0, dim = c(4, 3, iter))
exe_time = array(0, dim = c(4, 3, iter))

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
    ftsvd.est = FTSVD(Y, f.grid = c(list(NULL),list(NULL),list(tn)), rank=r, penalty = 1e-5, L=20, initialization = "Spectral")
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

