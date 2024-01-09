
cjs.mcmc <- function(Y, n.mcmc){
  
  #
  # setup variables
  #
  
  n = dim(Y)[1]
  T = dim(Y)[2]
  
  p.save = rep(0, n.mcmc)
  phi.save = rep(0, n.mcmc)
  
  #
  # priors and starting values
  #
  
  alpha.p = 1
  beta.p = 1
  
  alpha.phi = 1
  beta.phi = 1
  
  p = 0.8
  phi = 0.65
  
  Z = matrix(0, n, T)
  Z[,1] = 1
  Z[Y==1] = 1
  
  for(i in 1:n){
    
    if(sum(Y[i,]) > 0){
      
      idx.tmp = max((1:T)[Y[i,]==1])
      Z[i, 1:idx.tmp] = 1
      
    }
    
  }
  
  #
  # MCMC loop
  #
  
  for(k in 1:n.mcmc){
    
    #
    # sample z
    #
    
    for(i in 1:n){
      
      for(t in 2:(T-1)){
        
        if(sum(Y[i, t:T])==0){
          
          if(Z[i, t-1]==1 & Z[i, t+1]==0){
            
            psi.tmp = phi * (1 - p)
            psi.tmp = psi.tmp / (psi.tmp + 1)
            Z[i, t] = rbinom(1, 1, psi.tmp)
            
          }
          
          if(Z[i, t-1]==0){
            Z[i, t] = 0
          }
          
          if(Z[i, t+1]==1){
            Z[i, t] = 1
          }
        
        }
        
      }
      
      if(Z[i, T-1]==1 & Y[i, T]==0){
        
        psi.tmp = phi * (1 - p)
        psi.tmp = psi.tmp / (psi.tmp + 1 - phi)
        Z[i, T] = rbinom(1, 1, psi.tmp)
      }
      
      if(Z[i, T-1]==0 & Y[i, T]==0){
        Z[i, T] = 0
      }
      if(Y[i, T]==1){
        Z[i, T] = 1
      }
    
    }
    
    #
    # sample p
    #
    
    p = rbeta(1, sum(Y[,-1][Z[,-1]==1]) + alpha.p, sum(1-Y[,-1][Z[,-1]==1]) + beta.p)
    
    #
    # sample phi
    #
    
    phi = rbeta(1, sum(Z[,-1][Z[,-T]==1]) + alpha.phi, sum(1-Z[,-1][Z[,-T]==1]) + beta.phi)
    
    #
    # save samples
    #
    
    p.save[k] = p
    phi.save[k] = phi
  
  }
  
  #
  # write output
  #
  
  list(p.save=p.save, phi.save=phi.save, n.mcmc=n.mcmc)
  
}



#
# fit model to data
#

n.mcmc = 50000
set.seed(1)

mcmc.out = cjs.mcmc(Y = CH, n.mcmc = n.mcmc)

layout(matrix(1:2, 1, 2))
hist(mcmc.out$p.save, breaks = 50, col = 8, xlim = c(0,1), xlab = "p")
hist(mcmc.out$phi.save, breaks = 50, col = 8, xlim = c(0,1), xlab = bquote(phi))


apply(cbind(mcmc.out$p.save, mcmc.out$phi.save)[-(1:2000), ], 2, mean)
apply(cbind(mcmc.out$p.save, mcmc.out$phi.save)[-(1:2000), ], 2, quantile, c(0.025, 0.975))


