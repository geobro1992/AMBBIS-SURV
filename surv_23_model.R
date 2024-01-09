sink("surv_23_RE.txt")
cat("
    model {
    
    #--------------------------------------
    # Parameters:
    # phi: survival probability
    # p: capture probability
    #--------------------------------------
    
    sigma.phi~dunif(0,10)
    tau.phi<-pow(sigma.phi,-2)
    sigma2.phi<-pow(sigma.phi,2)  #optional if we want to see variance too
    
    sigma.p~dunif(0,10)
    tau.p<-pow(sigma.p,-2)
    sigma2.p<-pow(sigma.p,2)  #optional if we want to see variance too
    
    # survival priors
    mean.phi ~ dbeta(1,1)   # mean survival with a beta prior
    mu.phi <- log(mean.phi / (1-mean.phi))
    
    mean.p ~ dbeta(1,1)   # mean detection with a beta prior
    mu.p <- log(mean.p / (1-mean.p))
    
    # survival
    for (i in 1:M){
    for(t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu.phi + b1.mu[t] * zuse[i,t] + epsilon.phi[t]  # size dependent survival
    logit(p[i,t]) <- mu.p + epsilon.p[t]  # detection with random effect of year
    }
    }
    
    for (t in 1:n.occasions-1){ 
    epsilon.phi[t] ~ dnorm(0,tau.phi)
    epsilon.p[t] ~ dnorm(0,tau.p)
    b1.mu[t] ~ dnorm(0,1)I(0,)  # coefficient for linear size effect on survival (phi)
    }
    
    #loop through M individuals and T years
    for (i in 1:M){
    for(t in 1:n.occasions){
    zuse[i,t] <- (x[i,t] - mean.svl)/sd.svl
    }
    }
    
    
    #loop through M individuals
    for (i in 1:M){
    #state transitions and likelihood for the first primary period
    z[i,f[i]] <- 1
    
    #state transition and likelihood for primary periods 2:T
    for(t in (f[i]+1):n.occasions){ 
    
    #State process: draw z(t) conditional on  z(t-1)
    mu2[i,t] <- z[i, t-1] * phi[i, t-1]
    z[i,t] ~ dbern(mu2[i,t]) 
    
    #Observation model 
    obs.mu[i,t] <- z[i,t]*p[i,t-1]
    y[i,t] ~ dbern(obs.mu[i,t])
    
    } #t
    } #i
    
    } #end model
    ",fill = TRUE)
sink()