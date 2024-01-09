sink("model_JS.txt")
cat("
    model {
    
    #--------------------------------------
    # Parameters:
    # phi: survival probability of uninfected
    # gamma: removal entry probability of infected
    # p: capture probability of uninfected
    #--------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive
    # 3 dead
    # 4 skipping
    # Observations (O):
    # 1 seen
    # 2 not seen
    # 
    #--------------------------------------
    
    # Priors and constraints
    
    
    beta ~ dnorm(0, 0.001)I(-1,1)
    sigma ~ dunif(0,0.1)
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)
    p ~ dunif(0.5, 1)     # Prior for mean capture

    for (t in 1:n.occasions-1) {
    mean.phi[t] ~ dunif(0, 1)   # Prior for mean survival
    mu[t] <- log(mean.phi[t] / (1-mean.phi[t]))
    }

 
#----------------------------------------    
    # Dirichlet prior for entry probabilities

    for (t in 1:n.occasions) {
    gam[t] ~ dgamma(1,1)
    b[t] <- gam[t] / sum(gam[1:n.occasions])
    }
    
    #-----------------------------------------------
    # Convert entry probs to conditional entry probs
    nu[1] <- b[1]
    
    for (t in 2:n.occasions){
    nu[t] <- b[t] / (1 - sum(b[1:(t-1)]))
    } #t
    

for (i in 1:M){  
    # Define probabilities of state S(t+1) given S(t)
  for (t in 1:(n.occasions-1)){
    logit(phi[i,t]) <- mu[t] + beta*zuse[i,t] + epsilon[i]  
  }
}
    
  for (i in 1:M){ 
    epsilon[i] ~ dnorm(0,tau)
  }

#loop through M individuals and T years
for (i in 1:M){
  for(t in 1:n.occasions){
    zuse[i,t] <- (x[i,t] - mean.svl)/sd.svl
  }
}

#time-variant parameters 
for(t_ in 1:T){ #loop through primary periods
psi[t_]~dunif(0,1) #inclusion probability
} #t_

#first state transition (pure nuisance; strictly from outside-pop to part of marked-pop)
trmat0[1] <- (1-psi[1]) #remains not-yet-in-pop
trmat0[2] <- 0
trmat0[3] <- psi[1]*(1-lambda) #inclusion into pop, goes outside study are
trmat0[4] <- psi[1]*lambda #inclusion into pop, goes inside study


# State-transition and observation matrices 	
  for (i in 1:M){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    ps[1,i,t,1] <- 1 - nu[t]
    ps[1,i,t,2] <- nu[t]
    ps[1,i,t,3] <- 0
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- phi[i,t]
    ps[2,i,t,3] <- 1 - phi[i,t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1

    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- 0
    po[1,i,t,2] <- 1
    po[2,i,t,1] <- p
    po[2,i,t,2] <- 1 - p
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 1
  } #t
} #i
    # Likelihood 
  for (i in 1:M){
  # Define latent state at first occasion
  z[i, 1] <- 1   # Make sure that all M individuals are in state 1 at t = 1
    
  for (t in 2:n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, ])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1, ])
    y.new[i, t] ~ dcat(po[z[i, t], i, t-1, ])
  } #t
} #i
    # Calculate derived population parameters
    for (i in 1:M){
    for (t in 2:n.occasions){
    al[i,t-1] <- equals(z[i,t], 2)
    } #t
    for (t in 1:(n.occasions-1)){
    d[i,t] <- equals(z[i,t] - al[i,t], 0)
    } #t   
    alive[i] <- sum(al[i,])
    } #i
    for (t in 1:(n.occasions-1)){
    N[t] <- sum(al[,t])        # Actual population size
    B[t] <- sum(d[,t])         # Number of entries
    } #t
    for (i in 1:M){
    w[i] <- 1-equals(alive[i],0)
    } #i
    Nsuper <- sum(w[])            # Superpopulation size
   
    }
    ",fill = TRUE)
sink()


sink("model_JS1.txt")
cat("
    model {
    
    #--------------------------------------
    # Parameters:
    # phi: survival probability of uninfected
    # gamma: removal entry probability of infected
    # p: capture probability of uninfected
    #--------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive
    # 3 dead
    # Observations (O):
    # 1 seen
    # 2 not seen
    # 
    #--------------------------------------
    
    # Priors and constraints
    mean.phi ~ dunif(0.3, 1)   # Prior for mean survival
    mu <- log(mean.phi / (1-mean.phi))
    beta ~ dnorm(0, 0.001)I(-1,1)
    sigma ~ dunif(0,0.1)
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)
    p ~ dunif(0, 1)     # Prior for mean capture

   #----------------------------------------    
    # Dirichlet prior for entry probabilities
    
    beta[1] ~ dgamma(1, 1)
    beta[2] ~ dgamma(1, 1)
    beta[3] ~ dgamma(.01, 1)
    beta[4] ~ dgamma(.01, 1)
    beta[5] ~ dgamma(1, 1)
    beta[6] ~ dgamma(.01, 1)
    beta[7] ~ dgamma(.01, 1)
    beta[8] ~ dgamma(1, 1)
    
    for (t in 1:n.occasions) {
    b[t] <- beta[t] / sum(beta[1:n.occasions])
    }
    
    #-----------------------------------------------
    # Convert entry probs to conditional entry probs
    nu[1] <- b[1]
    
    for (t in 2:n.occasions){
    nu[t] <- b[t] / (1 - sum(b[1:(t-1)]))
    } #t
    

for (i in 1:M){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[i,t] + epsilon[i]  
    }
  }

  for (i in 1:M){ 
    epsilon[i] ~ dnorm(0,tau)
    }

# State-transition and observation matrices 	
  for (i in 1:M){  
# Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
    ps[1,i,t,1] <- 1 - nu[t]
    ps[1,i,t,2] <- nu[t]
    ps[1,i,t,3] <- 0
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- phi[i,t]
    ps[2,i,t,3] <- 1 - phi[i,t]
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- 0
    ps[3,i,t,3] <- 1
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- 0
    po[1,i,t,2] <- 1
    po[2,i,t,1] <- p
    po[2,i,t,2] <- 1 - p
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 1
    } #t
    } #i
    # Likelihood 
    for (i in 1:M){
    # Define latent state at first occasion
    z[i, 1] ~ dbern(nu[1])   # Make sure that all M individuals are in state 1 at t = 1
    
    for (t in 2:n.occasions){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, ])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1, ])
    y.new[i, t] ~ dcat(po[z[i, t], i, t-1, ])
    } #t
    } #i
    # Calculate derived population parameters
    for (t in 1:(n.occasions-1)){
    qnu[t] <- 1-nu[t]
    }
    cprob[1] <- nu[1]
    for (t in 2:(n.occasions-1)){
    cprob[t] <- nu[t] * prod(qnu[1:(t-1)])
    } #t
    psi <- sum(cprob[])            # Inclusion probability
    for (t in 1:(n.occasions-1)){
    b[t] <- cprob[t] / psi      # Entry probability
    } #t
    for (i in 1:M){
    for (t in 2:n.occasions){
    al[i,t-1] <- equals(z[i,t], 2)
    } #t
    for (t in 1:(n.occasions-1)){
    d[i,t] <- equals(z[i,t] - al[i,t], 0)
    } #t   
    alive[i] <- sum(al[i,])
    } #i
    for (t in 1:(n.occasions-1)){
    N[t] <- sum(al[,t])        # Actual population size
    B[t] <- sum(d[,t])         # Number of entries
    } #t
    for (i in 1:M){
    w[i] <- 1-equals(alive[i],0)
    } #i
    Nsuper <- sum(w[])            # Superpopulation size
    
    }
    ",fill = TRUE)
sink()


sink("model_JS2.txt")
cat("
    model {
    
    #--------------------------------------
    # Parameters:
    # phi: survival probability of uninfected
    # gamma: removal entry probability of infected
    # p: capture probability of uninfected
    #--------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive
    # 3 dead
    # Observations (O):
    # 1 seen
    # 2 not seen
    # 
    #--------------------------------------
    
#hyperpriors
   for(t_ in 1:(T-1)){
    phi.mu[t_] ~ dunif(0,1) #mean survival with a Uniform prior
    }
    g1.mu ~ dnorm(0,tau.g1) #mean gamma1, temp migration out | out
    logit(gamma1) <- g1.mu #prob migrate out|out (probability)
    g2.mu ~ dnorm(0,tau.g2) #mean gamma2,  temp migration out | in
    logit(gamma2) <- g2.mu #prob migrate out|in (probability)
    lambda <- (1-gamma1)/(gamma2-gamma1+1) #long-term probability inside study area
    p ~ dunif(0,1) detection probability
    b1.mu ~ dnorm(0,1)  # coefficient for linear size effect on survival (phi)
    
    sigma.g1 ~ dunif(0,0.1)
    tau.g1 <- pow(sigma.g1, -2)
    sigma2.g1 <- pow(sigma.g1, 2)

    sigma.g2 ~ dunif(0,0.1)
    tau.g2 <- pow(sigma.g2, -2)
    sigma2.g2 <- pow(sigma.g2, 2)


#time-variant parameters 

for(t_ in 1:T){ #loop through primary periods for inclusion probs
    gam[t_] ~ dgamma(1,1)
    b[t_] <- gam[t_] / sum(gam[1:T])
    } #t_

    # Convert entry probs to conditional entry probs
    psi[1] <- b[1]
    
    for (t_ in 2:T){
    psi[t_] <- b[t_] / (1 - sum(b[1:(t_-1)]))
    } #t_
    


#first state transition (pure nuisance; strictly from outside-pop to part of marked-pop)
trmat0[1] <- (1-psi[1]) #remains not-yet-in-pop
trmat0[2] <- 0
trmat0[3] <- psi[1]*(1-lambda) #inclusion into pop, goes outside study are
trmat0[4] <- psi[1]*lambda #inclusion into pop, goes inside study

#state transitions (2:T)
for(t_ in 1:(T-1)){
  
# size dependent phi (linear and quadratic term)
  for (i in 1:M){
    logit(phi[i,t_]) <- phi.mu[t_] + b1.mu * zuse[i,t_]


   
#state transitions 
#trmat: transition matrix for Markovian latent-states
    #1 =not yet in population; 2=dead;3=offsite;4=onsite (only observable state)
    #transition are from the column --> rows
    #trmat[row,column,time] = [state at time=t_; state at time t_-1; time=t_]
    #notice that the primary periods are offset by 1 (because we already dealt with T=1)
    trmat[1,1,t_,i]<- 1-psi[t_+1] #excluded from pop
    trmat[2,1,t_,i] <- 0 #dead
    trmat[3,1,t_,i] <- psi[t_+1]*(1-lambda) #inclusion into pop,outside study
    trmat[4,1,t_,i] <- psi[t_+1]*lambda #inclusion into pop,inside study
    trmat[1,2,t_,i]<- 0
    trmat[2,2,t_,i]<- 1 #stay dead
    trmat[3,2,t_,i]<- 0
    trmat[4,2,t_,i]<- 0
    trmat[1,3,t_,i]<- 0
    trmat[2,3,t_,i]<- 1-phi[i,t_] #dies outside
    trmat[3,3,t_,i]<- gamma1*phi[i,t_] #stays outside | outside
    trmat[4,3,t_,i]<- (1-gamma1)*phi[i,t_] #reenters study area | outside
    trmat[1,4,t_,i]<- 0 #
    trmat[2,4,t_,i]<- 1-phi[i,t_] #dies inside
    trmat[3,4,t_,i]<- gamma2*phi[i,t_] #leaves study area | inside
    trmat[4,4,t_,i]<- (1-gamma2)*phi[i,t_] #stays inside | inside
    }
}    
    #loop through M individuals
    for (i in 1:M){
    
    #state transitions and likelihood for the first primary period
    z[i,1]~ dcat(trmat0[1:4]) #z at time 0

    #Observation error
    obs.mu[i,1] <- equals(z[i,1],4)*p
    y[i,1] ~ dbern(obs.mu[i,1]) #inverse-logit transform
    
    #state transition and likelihood for primary periods 2:T
    for(t_ in 2:T){ 
    
    #State process: draw z(t_) conditional on  z(t_-1)
    z[i,t_] ~ dcat(trmat[1:4, z[i,t_-1] , t_-1, i])
    
    #Observation error y[i,tt_,t_] ~ Bernoulli condition on being inside z=4
    obs.mu[i,t_] <- equals(z[i,t_],4)*p
    y[i,t_] ~ dbern(obs.mu[i,t_]) #inverse-logit transform
    
    } #t_
    } #i
    
  #loop through M individuals and T years
    for (i in 1:M){
    for(t_ in 1:T){
    zuse[i,t_] <- (x[i,t_] - mean.svl)/sd.svl
    }
    }
    
    } #end model
    ",fill = TRUE)
sink()

