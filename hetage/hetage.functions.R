# hetage.functions.R

# Shirley Pledger, Murray Efford, Ernestynne Walsh.
# June 2010.

# -------------------------- LIST OF FUNCTIONS ----------------------------------------

# Data processing:
# ---------------- 
# condense.data
# hetage.process.data

# Model fitting:
# --------------

# hetage.fit.model
# choose.start

# Summaries of results:
# ---------------------

# hetage.summary
# hetage.stopover
# hetage.superpop


# Functions used internally for calculations:
# -------------------------------------------

# Functions using C code:

# mLL  
# psibuild 

# Mathematical functions:

# logit     
# expit     
# plg.fn    
# lgp.fn    
# logitmean 

# Matrix construction function:

# nlc       
# nlch      
# nlphi     
# nlrv      
# nlth      
# nlov      

# Start functions, used internally in model fitting:

# start.phic.pc     
# start.phic.pt     
# start.phic.pa     
# start.phic.pta    
# start.phic.ph     
# start.phic.pth    
# start.phit.pc     
# start.phit.pt     
# start.phit.pa     
# start.phit.pta    
# start.phit.ph     
# start.phit.pth    
# start.phia.pc     
# start.phia.pt     
# start.phia.pa     
# start.phia.pta    
# start.phia.ph     
# start.phia.pth    
# start.phita.pc    
# start.phita.pt    
# start.phita.pa    
# start.phita.pta   
# start.phita.ph    
# start.phita.pth   
# start.phih.pc     
# start.phih.pt     
# start.phih.pa     
# start.phih.pta    
# start.phih.ph     
# start.phih.pth    
# start.phith.pc    
# start.phith.pt    
# start.phith.pa    
# start.phith.pta   
# start.phith.ph    
# start.phith.pth   
# start.phiW.pc     
# start.phiW.pt     
# start.phiW.pa     
# start.phiW.pta    
# start.phiW.ph     
# start.phiW.pth    

# Functions used internally to build up a vector
# into a full parameter list:

# v.to.f.phic.pc     
# v.to.f.phic.pt     
# v.to.f.phic.pa     
# v.to.f.phic.pta    
# v.to.f.phic.ph     
# v.to.f.phic.pth    
# v.to.f.phit.pc     
# v.to.f.phit.pt     
# v.to.f.phit.pa     
# v.to.f.phit.pta    
# v.to.f.phit.ph     
# v.to.f.phit.pth    
# v.to.f.phia.pc     
# v.to.f.phia.pt     
# v.to.f.phia.pa     
# v.to.f.phia.pta    
# v.to.f.phia.ph     
# v.to.f.phia.pth    
# v.to.f.phita.pc    
# v.to.f.phita.pt    
# v.to.f.phita.pa    
# v.to.f.phita.pta   
# v.to.f.phita.ph    
# v.to.f.phita.pth   
# v.to.f.phih.pc     
# v.to.f.phih.pt     
# v.to.f.phih.pa     
# v.to.f.phih.pta    
# v.to.f.phih.ph     
# v.to.f.phih.pth    
# v.to.f.phith.pc    
# v.to.f.phith.pt    
# v.to.f.phith.pa    
# v.to.f.phith.pta   
# v.to.f.phith.ph    
# v.to.f.phith.pth   
# v.to.f.phiW.pc     
# v.to.f.phiW.pt     
# v.to.f.phiW.pa     
# v.to.f.phiW.pta    
# v.to.f.phiW.ph     
# v.to.f.phiW.pth    



# ---------------------------------INITIAL DATA SETUP----------------------------------

# Function to condense a capture matrix into like capture
# histories, and create the corresponding frequency vector.

condense.data <- function(inmat)
   {
   RR <- nrow(inmat)
   KK <- ncol(inmat)
   for (jj in KK:1)
   inmat <- inmat[sort.list(inmat[,jj]),]
   # Create new condensed matrix and associated y vector:
   cond.mat <- matrix(0,RR,KK)
   cond.y   <- rep(0,RR)
   current.ch <- inmat[1,]
   cond.mat[1,] <- inmat[1,]
   current.y  <- 1
   rowno    <- 1
   for (ii in 2:RR)
      {
      if (sum((inmat[ii,] - current.ch)^2) < 0.5)  # Same ch
         {
         current.y <- current.y + 1
         if (ii == RR)
            { 
            cond.y[rowno] <- current.y
            cond.mat[rowno,] <- current.ch
            }
         }   else           # New ch
         {
         cond.y[rowno] <- current.y
         cond.mat[rowno,] <- current.ch
         rowno <- rowno + 1
         current.y <- 1
         current.ch <- inmat[ii,]
         if (ii == RR)
            { 
            cond.y[rowno] <- current.y
            cond.mat[rowno,] <- current.ch
            }
         }
      }
   # Shorten the resultant objects:
   cond.mat <- cond.mat[apply(cond.mat,1,sum) > 0.5,]
   cond.y   <- cond.y[cond.y > 0.5]
   cbind(cond.mat,cond.y)
   }


hetage.process.data <- function(data, data.ud = NULL, data.nd = NULL, 
   data.ub = NULL, data.nb = NULL, time.intervals = NULL, mxage = NULL)
# Calculate various summaries; also f and l vectors for unknown age animals.
# Inputs 
#     data  - dataframe of capture-history matrix (unique rows) + 'freq' column
#             (number of animals with each history) 
#     ub dataframe has a final 'cohort' column
#             when extra information regarding (un)natural births/deaths is 
#             available additional data frames in the above format can be inputted
# Outputs (assignments to global variables)
#     x.mat
#     xud.mat
#     xnd.mat
#     xub.mat
#     xnb.mat
#     y.vect
#     yud.vect
#     ynd.vect
#     yub.vect
#     ynb.vect
#     Stimes
#     Etimes
#     n 
#     nud
#     nnd
#     nub
#     nnb
#     birthtime
#     K 
#     f.vect
#     fud.vect
#     fnd.vect
#     fub.vect
#     fnb.vect
#     l.vect
#     lud.vect
#     lnd.vect
#     lub.vect
#     lnb.vect
#     maxage
#     age1    matrix of age rows=cohorts, cols=sample times
#     age2    age at next sample time
{
  if (!is.null(data$freq))
      {
      x.mat <<- as.matrix(data[,-ncol(data)])
      y.vect <<- data$freq
      }  else
      {
      xy.mat <- condense.data(as.matrix(data))
      x.mat <<- xy.mat[,-ncol(xy.mat)]
      y.vect <<- xy.mat[,ncol(xy.mat)]
      }
   R      <<- nrow (x.mat)         # Number of distinct CHs
   K      <<- ncol (x.mat)        # Number of occasions
   n      <<- sum(y.vect)         # Number of animals
  
  lud.vect <<- NULL
  lnd.vect <<- NULL
  lub.vect <<- NULL
  lnb.vect <<- NULL

##########################
   ## Unnatural death information if available

  
   if(!is.null(data.ud$freq))
   {
   xud.mat <<- as.matrix(data.ud[,-ncol(data.ud)])
   #   xud.mat  <<- xmat(data.ud)     # see miscfn.R for xmat()
   Rud      <- nrow (xud.mat)         # Number of distinct CHs
   if (is.null(data.ud$freq)) yud.vect <<- rep(1, R)
   else                    yud.vect <<- data.ud$freq
   nud      <<- sum(yud.vect)         # Number of animals
  # print(data.ud)
  # print(c("xud",xud.mat,"rud",Rud,"yud",yud.vect,"nud",nud))
     }

   # if no unnatural death information is available
   else
     {
     xud.mat  <<- 0    
     Rud      <- 0     
     yud.vect <<- 0
     nud      <<-0       
     }
    # print(c(xud.mat,yud.vect,(length(yud.vect))))  added to check why yud =1

   ##########################
   ## Natural death information if available

  if(!is.null(data.nd))
     {
     #xnd.mat  <<- xmat(data.nd)     # see miscfn.R for xmat()
     xnd.mat <<- as.matrix(data.nd[,-ncol(data.nd)])
     Rnd      <- nrow (xnd.mat)       # Number of distinct CHs
     if (is.null(data.nd$freq)) ynd.vect <<- rep(1, Rnd)
     else                    ynd.vect <<- data.nd$freq
     nnd      <<- sum(ynd.vect)         # Number of animals with known deaths
     }
   # if no natural death information is available
   else
     {
     xnd.mat  <<- 0    
     Rnd      <- 0       
     ynd.vect <<- 0
     nnd      <<-0       
     }

   ##########################
   ## Unnatural birth information if available
   if(!is.null(data.ub))
     {
    # xub.mat  <<- xmat(data.ub)     # see miscfn.R for xmat()
     xub.mat <<- as.matrix(data.ub[,-ncol(data.ub)])
     Rub      <- nrow (xub.mat)       # Number of distinct CHs
     if (is.null(data.ub$freq)) yub.vect <<- rep(1, Rub)
     else                    yub.vect <<- data.ub$freq
     nub      <<- sum(yub.vect)         # Number of animals
     birthtime    <<- data.ub$cohort    # Age cohorts each animal belongs to
     }
   # if no unnatural birth information is available
   else
     {
     xub.mat  <<- 0    
     Rub      <- 0     
     yub.vect <<- 0
     nub      <<-0
     birthtime <<- 0
     }

   ##########################
   ## Natural birth information if available
   if(!is.null(data.nb))
     {
    # xnb.mat  <<- xmat(data.nb)     # see miscfn.R for xmat()
     xnb.mat <<- as.matrix(data.nb[,-ncol(data.nb)])
     Rnb      <- nrow (xnb.mat)       # Number of distinct CHs
     if (is.null(data.nb$freq)) ynb.vect <<- rep(1, Rnb)
     else                    ynb.vect <<- data.nb$freq
     nnb      <<- sum(ynb.vect)         # Number of animals
     }
   # if no natural birth information is available
   else
     {
     xnb.mat  <<- 0    
     Rnb      <- 0     
     ynb.vect <<- 0
     nnb      <<-0       
     }

   ##########################
   ## first and last

   f.vect <<- rep (0,R)           # First capture in each ch
   l.vect <<- rep (0,R)           # Last capture in each ch

   for (i in 1:R) 
      {
      marked <- 0
      for (j in 1:K) if (x.mat[i,j]==1)
         {
        if (marked==0)
            {
            marked    <- 1
            f.vect[i] <<- j
            }
         l.vect[i] <<- j
         }
      }
   
   ##########################
   ## first and last when deaths are unnatural
   if(!is.null(data.ud))
   {
   fud.vect <<- rep (0,Rud)           # First capture in each ch
   lud.vect <<- rep (0,Rud)           # Last capture in each ch
   #print(c("fud overwritten ", fud.vect, "lud overwritten ", lud.vect))
   for (i in 1:Rud) 
      {
      marked <- 0
      for (j in 1:K) if (xud.mat[i,j]==1)
         {
         if (marked==0)
            {
            marked    <- 1
            fud.vect[i] <<- j
            }
         lud.vect[i] <<- j
         }
      }
   # print(c("f",fud.vect,"l",lud.vect))
    }

   ##########################
   ## first and last when deaths are natural
   if(!is.null(data.nd))
   {
   fnd.vect <<- rep (0,Rnd)           # First capture in each ch
   lnd.vect <<- rep (0,Rnd)           # Last capture in each ch

   for (i in 1:Rnd) 
      {
      marked <- 0
      for (j in 1:K) if (xnd.mat[i,j]==1)
         {
         if (marked==0)
            {
            marked    <- 1
            fnd.vect[i] <<- j
            }
         lnd.vect[i] <<- j
         }
      }
   }

   ##########################
   ## first and last when births are unnatural
   
   if(!is.null(data.ub))
   {
   fub.vect <<- rep (0,Rub)           # First capture in each ch
   lub.vect <<- rep (0,Rub)           # Last capture in each ch

   for (i in 1:Rub) 
      {
      marked <- 0
      for (j in 1:K) if (xub.mat[i,j]==1)
         {
         if (marked==0)
            {
            marked    <- 1
            fub.vect[i] <<- j
            }
         lub.vect[i] <<- j
         }
      }
    }

   ##########################
   ## first and last when births are natural
   
   if(!is.null(data.nb))
   {
   fnb.vect <<- rep (0,Rnb)           # First capture in each ch
   lnb.vect <<- rep (0,Rnb)           # Last capture in each ch

   for (i in 1:Rnb) 
      {
      marked <- 0
      for (j in 1:K) if (xnb.mat[i,j]==1)
         {
         if (marked==0)
            {
            marked    <- 1
            fnb.vect[i] <<- j
            }
         lnb.vect[i] <<- j
         }
      }
    }

   ##########################
   ## maximum age : optionally override data

  maxage <<- max(l.vect-f.vect+1)
  if (!is.null(lud.vect)) maxage <<- max(maxage,lud.vect-fud.vect+1)
  if (!is.null(lnd.vect)) maxage <<- max(maxage,lnd.vect-fnd.vect+1)
  if (!is.null(lub.vect)) maxage <<- max(maxage,lub.vect-fub.vect+1)
  if (!is.null(lnb.vect)) maxage <<- max(maxage,lnb.vect-fnb.vect+1)

  if (!is.null(mxage))
      maxage <<- min(mxage, maxage)

   ##########################
   ## sampling times and ages

   if (is.null(time.intervals)) Stimes <<- (1:K)-0.5   ### check this later
   else Stimes <<- cumsum (c(0.5, time.intervals))

   Etimes <<- c(0,(Stimes[-length(Stimes)]+Stimes[-1])/2)   # used by silerfn.R
   temp   <- t(outer (Stimes,Etimes,'-'))
   age1   <<- temp[-K,-K]  # matrix of age rows=cohorts, cols=sample times
   age2   <<- temp[-K,-1]  # age at next sample time
   age1[lower.tri(age1)] <<- 0    # for tidiness
   age2[lower.tri(age2)] <<- 0    # for tidiness


   ###########################
   ## summary

    first.capt <- tabulate (rep(f.vect, y.vect), nbins = K)
    last.capt <- tabulate (rep(l.vect, y.vect), nbins = K)
    temp <- x.mat[rep( 1:nrow(x.mat), y.vect), ]
    nj <- apply (temp,2,sum)
 
    out <- rbind (nj, first.capt, last.capt)
    out <- data.frame(out)
    names(out) <- paste('sample', 1:K, sep = '')
    row.names(out) <- c(
        'Number caught in sample j',
        'Number first caught in sample j',
        'Number last caught in sample j')

    list (TotalIndividuals=n, by.sample = out)

   ###########################
   ## summary for unnatural deaths
   if(!is.null(data.ud))
   {
    firstud.capt <- tabulate (rep(fud.vect, yud.vect), nbins = K)
    lastud.capt <- tabulate (rep(lud.vect, yud.vect), nbins = K)
    tempud <- xud.mat[rep( 1:nrow(xud.mat), yud.vect), ]
    njud <- apply (tempud,2,sum)
 
    outud <- rbind (njud, firstud.capt, lastud.capt)
    outud <- data.frame(outud)
    names(outud) <- paste('sample', 1:K, sep = '')
    row.names(outud) <- c(
        'Number caught in sample j with unnatural deaths',
        'Number first caught in sample j with unnatural deaths',
        'Number last caught in sample j with unnatural deaths')

    list (TotalIndividuals=nud, by.sample = outud)
    }

   ###########################
   ## summary for natural deaths
   if(!is.null(data.nd))
   {
    firstnd.capt <- tabulate (rep(fnd.vect, ynd.vect), nbins = K)
    lastnd.capt <- tabulate (rep(lnd.vect, ynd.vect), nbins = K)
    tempnd <- xnd.mat[rep( 1:nrow(xnd.mat), ynd.vect), ]
    njnd <- apply (tempnd,2,sum)
 
    outnd <- rbind (njnd, firstnd.capt, lastnd.capt)
    outnd <- data.frame(outnd)
    names(outnd) <- paste('sample', 1:K, sep = '')
    row.names(outnd) <- c(
        'Number caught in sample j with natural deaths',
        'Number first caught in sample j with natural deaths',
        'Number last caught in sample j with natural deaths')

    list (TotalIndividuals=nnd, by.sample = outnd)
   }

   ###########################
   ## summary for unnatural births

   if(!is.null(data.ub))
   {
    firstub.capt <- tabulate (rep(fub.vect, yub.vect), nbins = K)
    lastub.capt <- tabulate (rep(lub.vect, yub.vect), nbins = K)
    tempub <- xub.mat[rep( 1:nrow(xub.mat), yub.vect), ]
    njub <- apply (tempub,2,sum)
 
    outub <- rbind (njub, firstub.capt, lastub.capt)
    outub <- data.frame(outub)
    names(outub) <- paste('sample', 1:K, sep = '')
    row.names(outub) <- c(
        'Number caught in sample j with unnatural births',
        'Number first caught in sample j with unnatural births',
        'Number last caught in sample j with unnatural births')

    list (TotalIndividuals=nub, by.sample = outub, BirthCohortVector = birthtime)
   
   }

   ###########################
   ## summary for natural births

   if(!is.null(data.nb))
   {
    firstnb.capt <- tabulate (rep(fnb.vect, ynb.vect), nbins = K)
    lastnb.capt <- tabulate (rep(lnb.vect, ynb.vect), nbins = K)
    tempnb <- xnb.mat[rep( 1:nrow(xnb.mat), ynb.vect), ]
    njnb <- apply (tempnb,2,sum)
 
    outnb <- rbind (njnb, firstnb.capt, lastnb.capt)
    outnb <- data.frame(outnb)
    names(outnb) <- paste('sample', 1:K, sep = '')
    row.names(outnb) <- c(
        'Number caught in sample j with natural births',
        'Number first caught in sample j with natural births',
        'Number last caught in sample j with natural births')

    list (TotalIndividuals=nnb, by.sample = outnb)
   
    }
}

# ---------------------------------FIT ONE MODEL---------------------------------------


hetage.fit.model <- function (modelname, G = 1, start = NULL, getVC = TRUE, 
  printstart = FALSE, printeach = FALSE, maxiter = 5000)

# assume all data values are global to this function and available to mLL
# inadequate maxiter may produce obscure error message

# 'start' is a full parameter set (list with N, beta, phi, p & pi) 
# used by startfn()

{
  # making G available to other functions
  G <<- G

  # optional self start 

  if (is.null(start)) 
  {
    startv <- c(log(n*1.05), plg.fn(c(1-(K-1)/K,rep(1/K,K-1))), 
                logit(0.8), logit(0.5))
    start  <- v.to.f.phic.pc(startv)
  }  

  if(G>1){
    start$pi <- rep(1/G,G)
    for (g in 2:G){
      start$phi[,,g] <- start$phi[,,1]*(1-0.1*(g-1)) 
      start$p[,,g] <- start$p[,,1]*(1-0.1*(g-1)) 
      }
   }

  # set functions for this model
  startfn <- eval(as.name(paste('start.', modelname, sep='')))
  v.to.f  <- eval(as.name(paste('v.to.f.', modelname, sep='')))

  pars.start <- startfn (start)
  if (printstart) print (pars.start)

# Check it evaluates at starting vector:
# print(mLL(pars.start,expandfn=v.to.f))

  npar       <- length (pars.start)
  fneval     <<- 0
  
  this.fit   <- optim (
                  par       = pars.start,
                  fn        = mLL,
                  lower     = c(log(n),rep(-Inf,npar-1)),
                  method    = "L-BFGS-B",
                  control   = list(maxit=maxiter),
                  hessian   = getVC,
                  expandfn  = v.to.f,                      # passed to mLL
                  printeach = printeach                    # passed to mLL
                  )

  RD <- 2*this.fit$value

  if (this.fit$convergence>50) warning(this.fit$message)

  parameters <- v.to.f(this.fit$par)
  cohorts    <- paste('Cohort',1:K, sep='')
  samples    <- paste('Sample',1:K, sep='')
  betas      <- paste('Beta',1:npar, sep='')
  groups     <- paste('Group',1:G, sep='')

  if(G==1){
  dimnames (parameters$p)   <- list(cohorts, samples)
  dimnames (parameters$phi) <- list(cohorts[-K], samples[-K])
  }

  else{

    dimnames (parameters$p)   <- list(cohorts, samples, groups)

    if(is.null(dim(parameters$phi[3]))){
    dimnames (parameters$phi) <- list(cohorts[-K], samples[-K])
    } 

    else{
    dimnames (parameters$phi) <- list(cohorts[-K], samples[-K], groups[-G])
    }
  }

  # Coefficients for the retention curves
  if (substring(modelname,1,4) %in% c('phiW'))
      names(parameters$coefs) <- switch (substring(modelname,1,4),
          phiW = c('a', 'b'),  # shape and scale pars
          paste('c', 1:length(parameters$coefs), sep='')     ## stopgap
          )

  # Coefficients for any user models that require them 
  if (substring(modelname,1,5) %in% c('user1')) 
      names(parameters$coefs) <- switch (substring(modelname,1,5),
          user1 = c('a.beta','b.beta','a.phi','b.phi'),
          paste('c', 1:length(parameters$coefs), sep='')     ## stopgap
          )
  # More models can be added if they require coefficients

  VC <- NULL
  if (this.fit$convergence==0 & getVC) {
      VC <- try (solve(this.fit$hessian), silent = TRUE)
      if (inherits(VC, "try-error")) {
          warning ('Could not invert Hessian to compute variance-covariance matrix')
          VC <- NA
      }  
      else dimnames(VC) <- list (betas, betas)
  }

  AIC <- RD + 2*npar
  AICc <- NA
  if (2^K-npar-1 > 0) AICc <- AIC + 2*npar*(npar+1)/(2^K-npar-1)
  
  list (
    modelname   = modelname,
    maxLL       = -this.fit$value,     ## return LL, not negative LL
    RD          = RD,
    K           = K,
    maxage      = maxage,
    number.seen = list(n=sum(y.vect), u.deaths=nud, n.deaths=nnd, u.births=nub, n.births=nnb),
    npar        = npar, 
    AIC         = AIC,
    AICc        = AICc,
    convergence = this.fit$convergence,
    fncalls     = this.fit$counts[1],
    pars        = this.fit$par,                # non-redundant vector
    parameters  = parameters,                  # full N, beta, phi, p, coef
    VC          = VC
  )
}


# ---------------------------------- SUMMARISE MODELS ---------------------------------


hetage.summary <- function (modellist, AICorder=F, AICcorder=F)
{
  sm <- function(model) 
    if (is.null(model)) numeric(9) 
    else round(c(model$maxLL, model$RD, model$AIC, model$AICc, model$npar, model$K, 
                 length(model$parameters$pi),
                 model$maxage, model$conv, model$fncalls, exp(model$pars[1])),3)
  out1 <- data.frame(t(sapply (modellist, sm)))
  out1 <- data.frame(out1[,1:3], out1[,3]-min(out1[out1[,3]>0,3]),out1[,4],
                     out1[,4]-min(out1[out1[,4]>0,4]), out1[,5:11]) 
  names(out1) <- c("maxLL","RD","AIC","relAIC","AICc","relAICc","npar","K","G",
                   "maxage","conv", "fncalls", "Nhat")
  out1[out1$maxLL==0,] <- NA
  if (AICorder)  out1 <- out1[order(out1$AIC),]
  if (AICcorder) out1 <- out1[order(out1$AICc),]
  out1
}

#-------------------------------- STOPOVER DURATION -----------------------------------


hetage.stopover <- function (model, St = Stimes)
# Mean and SD of stopover duration for discrete arrival and departure times
# using parameter estimates for 'model' in 'modellist'.
# We assume all animals present at last sample time die (phi(K)=0)
{
    if (!is.null(model)) {
    
        psi      <- psibuild (model$parameters)$psi         # lower tri == 0
        Etimes   <- (c(St[1]-1,St[-length(St)]) + St)/2     # Entry times
        Xtimes   <- (St+ c(St[-1], St[length(St)]+1 ))/2    # Exit times
        duration <- t(outer (Xtimes, Etimes,'-'))           # Or simply (col(psi)-row(psi)+1)
        EX       <- sum (duration * psi)
        EX2      <- sum (duration^2 * psi)
        data.frame (SoDhat=EX, SE.SoDhat = sqrt(EX2 - EX^2))
    }
}

#----------------------------- SUPERPOPULATION SIZE  ----------------------------------


# Find estimated N and its S.E.

hetage.superpop <- function(model)
   {
   if (!is.null(model)) 
      {
      N  <- model$parameters$N
      VN <- invlogVar(model$pars[1],model$VC[1,1])
      data.frame(Nhat=N,seNhat=sqrt(VN))
      }
   }

#--------------------------------------------------------------------------------------

# --------------
# Curves
# --------------

Weibull <- function (x, aa, bb, a=1, coef = NULL)
# x = age
# kap = aa = shape parameter > 0
# gam = bb = scale parameter > 0
{
  if (!is.null(coef)) {
      aa <- coef[1]
      bb <- coef[2]
      if (length(coef)>2) a <- coef[3]
  }
  if (aa <= 0 | bb <= 0) stop ('Invalid parameters for Weibull curve')
  a * exp ( - ((x+1)/bb)^aa + (x/bb)^aa )

}

# e.g. plot (x=1:15, y=Weibull (1:15, aa=0.5, bb=0.2), ylim=c(0,1), xlab='Age', ylab='Phi')


#---------------------- ------------------------------------------------------
# Reference list of models and the preceding models to use for starting values

masterlist <- list (
   phic.pc    = 0,                   #  1
   phic.pt    = 1,                   #  2
   phic.pa    = 1,                   #  3
   phic.pta   = 3,                   #  4
   phic.ph    = 1,                   #  5
   phic.pth   = 2,                   #  6

   phit.pc    = 1,                   #  7
   phit.pt    = c(2,7),              #  8
   phit.pa    = c(3,9),              #  9
   phit.pta   = c(4,8,9),            # 10
   phit.ph    = c(1,9),              # 11
   phit.pth   = c(2,5,9),            # 12

   phia.pc    = 1,                   # 13
   phia.pt    = c(2,13),             # 14
   phia.pa    = c(3,13),             # 15
   phia.pta   = c(4,14,15),          # 16
   phia.ph    = 1,                   # 17
   phia.pth   = c(2,5,13),           # 18

   phita.pc   = c(7,13),             # 19
   phita.pt   = c(8,14,19),          # 20
   phita.pa   = c(9,15,19),          # 21
   phita.pta  = c(10,16,20,21),      # 22
   phita.ph   = c(1,11,17),          # 23
   phita.pth  = c(12,18,20,23),      # 24

   phih.pc   = 1,                    # 25
   phih.pt   = c(2,25),              # 26
   phih.pa   = c(3,25),              # 27
   phih.pta  = c(4,26,27),           # 28
   phih.ph   = c(5,25),              # 29
   phih.pth  = c(6,26,29),           # 30

   phith.pc  = c(1,7,25),            # 31
   phith.pt  = c(2,8,26),            # 32
   phith.pa  = c(3,9,27),            # 33
   phith.pta = c(4,10,28,32,33),     # 34
   phith.ph  = c(5,11,29,31),        # 35
   phith.pth = c(6,12,30,32,35),     # 36

   phiW.pc    = 1,                   # 37
   phiW.pt    = 1,                   # 38
   phiW.pa    = 1,                   # 39
   phiW.pta   = 1,                   # 40
   phiW.ph    = 1,                   # 41
   phiW.pth   = 1)                   # 42




#-------------------  CHOOSE STARTING MODEL---------------------

choose.start <- function (model, modellist)
# Return the full parameters of the best model from among those in vector 'model' 
{
  extract.maxLL <- function(x)
    if (is.null(x)) stop ('Starting model not fitted') 
    else x$maxLL
  refvector <- masterlist[[model]]
  refvector <- refvector[!sapply(modellist[refvector],is.null)]
  if (length(refvector)==0) NULL
  else if (refvector==0) NULL
  else {
      if (length(refvector)==1) modellist[[refvector]]$parameters                 
      else {
        maxLL     <- sapply(modellist[refvector], extract.maxLL)
        bestmodel <- refvector[match(max(maxLL), maxLL)]                 
        modellist[[bestmodel]]$parameters
      }
  }
}


#-------------------  LOGIT FUNCTION  --------------------------

# If p is really close to 0 or 1, set logit to some max or min.

logit <- function(pvec)
   {
   # Replace values if too near 0 or 1:
   pvec <- apply(cbind(pvec,rep(0.00001,length(pvec))),1,max)
   pvec <- apply(cbind(pvec,rep(0.99999,length(pvec))),1,min)
   # Output vector:
   log(pvec/(1-pvec))
   }

#-------------------  EXPIT FUNCTION  --------------------------

# Inverse of logit function, acts on a logit vector.

expit <- function(lvec)
   1/(1+exp(-lvec))

#-------------------  invlogitSE FUNCTION  --------------------------

# Inverse of SE logit function, acts on logit and logit SE vectors.
# Delta method Lebreton et al 1992 p 77

invlogitSE <- function(lvec, selvec)
   expit(lvec) * (1-expit(lvec)) * selvec


#-------------------  invlogitVar FUNCTION  --------------------------

# Inverse of Var logit function, acts on logit and logit Var vectors.
# Delta method Lebreton et al 1992 p 77

invlogitVar <- function(lvec, varlvec)
   expit(lvec)^2 * (1-expit(lvec))^2 * varlvec


#-------------------  invlogSE FUNCTION  --------------------------

# Inverse of SE log function, acts on log and log SE vectors.
# Delta method Lebreton et al 1992 p 77

invlogSE <- function(lvec, selvec)
  sqrt(exp(selvec^2) - 1) * exp(lvec)

#-------------------  invlogVar FUNCTION  --------------------------

# Inverse of Var log function, acts on log and log Var vectors.
# Delta method Lebreton et al 1992 p 77

invlogVar <- function(lvec, varlvec)
  (exp(varlvec) - 1) * exp(lvec)^2


#-------------------  PLG FUNCTION  ---------------------------

# Proportions to logits (one fewer elements)
# Use on beta vector or pi vector (which have constraint sum = 1).

plg.fn <- function(invec)
   {
   # How many elements?
   nel <- length(invec)
   # Case of 2 els:
   if (nel==2) 
   {
    outvec <- invec[1]
    # Return logit
    logit(outvec)
    }
   # Case of more than 2 els:
   else if (nel>2)
      {
      # Cumulative sum:
      cs <- cumsum(invec)
      # Output vector:
      outvec <- invec[1:(nel-1)]/(1-c(0,cs[1:(nel-2)]))
      # Return logit
      logit(outvec)
	}
   }


#-------------------  LGP FUNCTION  ---------------------------

# Logit to proportions (one extra element)
# Use to obtain beta or pi vector from shorter vector.
# Works for vector of length 1 or more.

lgp.fn <- function(invec)
   {
   # How many elements?
   nel <- length(invec)
   # Case of 1 el:
   if (nel==1) p.v <- expit(invec)
   # Case of more than 1 els:
   if (nel>1)
      {
      # Get inverse-logit vector:
      thisvec <- expit(invec)
      # Convert to proportions vector: need cum product of 1-thisvec
      cp <- cumprod(1-thisvec)
      # Return proportions vector:
      p.v <- thisvec*c(1,cp[1:(nel-1)])
      }
   c(p.v,1-sum(p.v))
   }


#-------------------  LOGITMEAN FUNCTION  ---------------------------

logitmean <- function(pvec)
   {
    pvec<-pvec[!is.na(pvec)]
    if ((min(pvec)> 0.00001) & (max(pvec) <0.99999)) 
      {
      mlpvec <- mean(log(pvec/(1-pvec)))
      out <- exp(mlpvec)/(1+exp(mlpvec))
      }  
   else out <- mean(pvec)
   out
   }

# ------------------------- FUNCTIONS TO ASSEMBLE MATRICES -------------------

# Functions used to fill 'non-lower-triangular' cells of a matrix
# By default lower triangle is zero

# nlc  - constant
# nlch - constant used for heterogeneous models
# nlphi- used to construct phi arrays for the heterogeneous models
# nlrv - replicate a vector in all rows
# nlth - construct the arrays for th models
# nlov - as with nlrv, but progressively offset the rows and truncate 'overhang'

#-------------------  Fill nonlower array with constant  ---------------------------

nlc <- function (const, K, missing=0)
{
 temp <- matrix (rep(missing,K*K), nr=K, nc=K)
 temp[upper.tri(temp,diag=T)] <- const
 temp2<-array(missing,c(K,K,G))
 for (gg in 1:G) temp2[,,gg]<-temp
 temp2
}


#-------------------  Fill nonlower array with constant  ---------------------------

# Version for heterogeneous model, vector of G constants.

nlch <- function(vect, K, G, missing=0)
{
  temp <- array(missing,c(K,K,G))
  for (gg in 1:G)
  {
    temp[,,gg] <- expit(vect[gg])
    temp[,,gg][lower.tri(temp[,,gg])] <- 0
  }
  temp
}

#-------------------  Fill phi array  ---------------------------


nlphi <- function (vect, K, G, missing=0)
{
  if (length(vect)!=G) stop ('Invalid length for vector in nlphi')
  phiarray<-array(missing,c(K,K,G))

  for (g in 1:G)
  {
  temp <- matrix (rep(missing,(K)*(K)), nr=K, nc=K)
  temp[upper.tri(temp,diag=T)] <- expit(vect[g])
  phiarray[,,g]<-temp
  }
  phiarray
}


# -------------------  Fill nonlower matrix with repeated vector  ---------------------------

nlrv <- function (vect, K, missing=0)
{
  if (length(vect)!=K) stop ('Invalid length for vector in nlrv')
  temp <- matrix (rep(vect,K), nr=K, nc=K, byrow=T)
  temp[lower.tri(temp)] <- missing
  temp2<-array(missing,c(K,K,G))
  for (g in 1:G) temp2[,,g]<-temp
  temp2
}

# -------------------  Fill nonlower array for th models  ---------------------------

nlth <- function (tauv, etav, K, G, missing=0)
{
   if (length(tauv)!=K | length(etav)!=G) stop ('Invalid length for vector in nlth')
  temp <- matrix (rep(expit(tauv),K), nr=K, nc=K, byrow=T)
  temp[lower.tri(temp)] <- missing
  temp2<-array(missing,c(K,K,G))
  temp2[,,1] <- temp
  for (g in 2:G)
    {
    temp2[,,g] <- matrix (rep(expit(tauv+etav[g]),K), nr=K, nc=K, byrow=T)
    temp2[,,g][lower.tri(temp2[,,g])] <- missing
  }
  temp2
}

#---Fill nonlower matrix with progressively offset repeated vector  ------------------

nlov <- function (vect, K, missing=0)
{
  if (length(vect)!=K) stop ('Invalid length for vector in nlov')
  temp <- matrix (rep(missing,K*K), nr=K, nc=K)
  temp[lower.tri(temp,d=T)] <- vect[unlist(sapply(K:1,seq))]
  temp2<-array(missing,c(K,K,G))
  for (g in 1:G) temp2[,,g]<-t(temp)
  temp2
}

# ----------------------- Minus Log Likelihood function -----------------


mLL <- function (pars, expandfn, printeach=F)
{
  temp <- expandfn(pars)

  # save for debugging if we should fail
  # currentpars <<- pars

  # check on validity of parameter values - probably slow
  ptemp <- unlist(temp)
  if (any(is.na(ptemp)) | any(is.nan(ptemp)) | any(ptemp==Inf) )
    1000    # return a large positive value
  else {

  LL <- .C("hetageLL", 
     as.integer (length(y.vect)),
     as.integer (y.vect),
     as.integer (length(yud.vect)),
     as.integer (yud.vect),
     as.integer (length(ynd.vect)),
     as.integer (ynd.vect),
     as.integer (length(yub.vect)),
     as.integer (yub.vect),
     as.integer (length(ynb.vect)),
     as.integer (ynb.vect),
     as.integer (K),
     as.integer (G),
     as.integer (birthtime),   
     as.integer (x.mat),
     as.integer (xud.mat),
     as.integer (xnd.mat),
     as.integer (xub.mat),
     as.integer (xnb.mat),
     as.double  (temp$N), 
     as.double  (temp$beta),
     as.double  (temp$pi),
     as.double  (temp$phi),
     as.double  (temp$p),
     LL          = double (1),
     resultcode  = integer(1)
    )$LL
  fneval <<- fneval + 1
  if (printeach)   cat(fneval, round(c(LL,pars),4), '\n')
  -LL

  }
}

# --------------- psibuild function ------------------------------


psibuild <- function (parameters)
# return matrix of psi (b,d)
{
  temp <-
  .C("psibd", 
     as.integer (length(y.vect)),
     as.integer (y.vect),
     as.integer (K),
     as.integer (x.mat),
     as.double  (parameters$N), 
     as.double  (parameters$beta),
     as.double  (parameters$phi),
     as.double  (parameters$p),
     psi         = double(K*K),
     pch         = double(length(y.vect)),
     pch0        = double(1),
     result      = double(1)
    )

  list (psi=matrix(temp$psi,nr=K), pch=temp$pch, pch0=temp$pch0)

}

# --------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
#                Start functions
# -------------------------------------------------------------------------------------

# ----------------
# phi (c) series
# ----------------

start.phic.pc <- function (parameters)
{
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,1]),      # [K+1]
    logit(parameters$p[1,1,1])         # [K+2]
    )
}

start.phic.pt <- function (parameters)
{
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,1]),      # [K+1]
    logit(parameters$p[1,2:K,1])       # [(K+2):(2*K)] (note p[1,1] fixed)
  )
}

start.phic.pa <- function (parameters)
{
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,1]),      # [K+1]
    logit(parameters$p[1,1:maxage,1])  # [(K+2):(K+maxage)] 
  )
}

start.phic.pta <- function (parameters)
{
  tau   <- logit(diag(parameters$p[,,1])[-1]) # K-1 taus, tau[1] fixed
  alpha <- rep(0,maxage-1)                         # alpha[1] = 0  
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,1]),      # [K+1]
    tau,                               # [(K+2):(2*K)]  
    alpha                              # [(2*K+1):(2*K+maxage-1)]
  )
}

start.phic.ph <- function (parameters)
{

  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,1]),      # [K+1]
    logit(parameters$p[1,1,]),         # [(K+2):(K+G+1)]
    plg.fn(parameters$pi)             # [(K+G+2):(K+2*G)]
  )
}


start.phic.pth <- function (parameters)
{  
  tau <-logit(diag(parameters$p[-1,-1,1])) # K-1 taus tau[1] fixed
  eta  <- rep(0,G-1)                       # G-1 etas
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,1]),      # [K+1]
    tau,                               # [(K+2):(2*K)]
    eta,                               # [(2*K+1):(2*K+G-1)]
    plg.fn(parameters$pi)             # [(2*K+G):(2*K+2*G-2)]
  )
}

# ----------------
# phi (t) series
# ----------------


start.phit.pc <- function (parameters)
{
  c (
    log(parameters$N),                  # [1]
    plg.fn(parameters$beta),            # [2:K]
    logit(parameters$phi[1,1:(K-1),1]), # [(K+1):(2*K-1)]
    logit(parameters$p[1,1,1])          # [2*K]
  )
}

start.phit.pt <- function (parameters)
{
  c (
    log(parameters$N),                   # [1]
    plg.fn(parameters$beta),             # [2:K]
    logit(parameters$phi[1, 1:(K-1),1]), # [(K+1):(2*K-1)]
    logit(parameters$p[1, 2:(K-1),1])    # [(2*K):(3*K-3)] (note p[1], p[K] fixed)
  )
}

start.phit.pa <- function (parameters)
# The top matrix row also gives age-specific rates...
{
  c (
    log(parameters$N),                   # [1]
    plg.fn(parameters$beta),             # [2:K]
    logit(parameters$phi[1,1:(K-1),1]),  # [(K+1):(2*K-1)]
    logit(parameters$p[1,2:maxage,1])    # [(2*K):(2*K+maxage-2)] (note p[1] fixed)
  )
}

start.phit.pta <- function (parameters)
{
  tau   <- logit(diag(parameters$p[,,1])[-c(1,K)]) # K-2 taus tau[1], tau[k] fixed 
  alpha <- rep(0, maxage-1)                        # M-1 alphas

  c (
    log(parameters$N),                   # [1]
    plg.fn(parameters$beta),             # [2:K]
    logit(parameters$phi[1,1:(K-1),1]),  # [(K+1):(2*K-1)]
    tau,                                 # [(2*K):(3*K-3)]
    alpha                                # [(3*K-2):(3*K+maxage-4)]
  )
}

start.phit.ph <- function (parameters)
{
  c (
    log(parameters$N),                   # [1]
    plg.fn(parameters$beta),             # [2:K]
    logit(parameters$phi[1,1:(K-1),1]),  # [(K+1):(2*K-1)]
    logit(parameters$p[1,1,]),           # [(2*K):(2*K+G-1)]
    plg.fn(parameters$pi)               # [(2*K+G):(2*K+2*G-2)]
  )
}


start.phit.pth <- function (parameters) 
{
  tau <- logit(diag(parameters$p[,,1])[-c(1,K)]) # K-2 taus tau[1], tau[K] fixed
  eta <- rep(0, G-1)                             # G-1 etas
  c (
    log(parameters$N),                  # [1]
    plg.fn(parameters$beta),            # [2:K]
    logit(parameters$phi[1,1:(K-1),1]), # [(K+1):(2*K-1)]
    tau,                                # [(2*K):(3*K-3)]
    eta,                                # [(3*K-2):(3*K+G-4)]
    plg.fn(parameters$pi)              # [(3*K+G-3):(3*K+2*G-5)]
  )
}

# ----------------
# phi (a) series
# ----------------

start.phia.pc <- function (parameters)
{
  c (
    log(parameters$N),                       # [1]
    plg.fn(parameters$beta),                 # [2:K]
    logit(parameters$phi[1,1:(maxage-1),1]), # [(K+1):(K+maxage-1)]
    logit(parameters$p[1,1,1])               # [K+maxage]
  )
}

start.phia.pt <- function (parameters)
{
  c (
    log(parameters$N),                       # [1]
    plg.fn(parameters$beta),                 # [2:K]
    logit(parameters$phi[1,1:(maxage-1),1]), # [(K+1):(K+maxage-1)]
    logit(parameters$p[1,2:K,1])             # [(K+maxage):(2*K+maxage-2)] (note p[1] fixed)
  )
}

start.phia.pa <- function (parameters)
{
  c (
    log(parameters$N),                       # [1]
    plg.fn(parameters$beta),                 # [2:K]
    logit(parameters$phi[1,1:(maxage-1),1]), # [(K+1):(K+maxage-1)]
    logit(parameters$p[1,2:maxage,1])        # [(K+maxage):(K+2*maxage-2)] (note p[1] fixed)
  )
}

# New version with only M-2 p's:

start.phia.pa <- function (parameters)
{
  c (
    log(parameters$N),                       # [1]
    plg.fn(parameters$beta),                 # [2:K]
    logit(parameters$phi[1,1:(maxage-1),1]), # [(K+1):(K+maxage-1)]
    logit(parameters$p[1,2:(maxage-1),1])        # [(K+maxage):(K+2*maxage-3)] (p[1],p[M] fixed)
  )
}


start.phia.pta <- function (parameters)
{
  tau   <- logit(diag(parameters$p[,,1])[2:(K-1)]) # K-2 taus tau[1], tau[K] fixed
  alpha <- rep(0, maxage-1)                        # M-1 alphas

  c (
    log(parameters$N),                       # [1]
    plg.fn(parameters$beta),                 # [2:K]
    logit(parameters$phi[1,1:(maxage-1),1]), # [(K+1):(K+maxage-1)]
    tau,                                     # [(K+maxage):(2*K+maxage-3)]
    alpha                                    # [(2*K+maxage-2):(2*K+2*maxage-4)]
  )
}

start.phia.ph <- function (parameters)
{
  c (
    log(parameters$N),                       # [1]
    plg.fn(parameters$beta),                 # [2:K]
    logit(parameters$phi[1,1:(maxage-1),1]), # [(K+1):(K+maxage-1)]
    logit(parameters$p[1,1,]),               # [(K+maxage):(K+maxage+G-1)]
    plg.fn(parameters$pi)                   # [(K+maxage+G):(K+maxage+2*G-2)]
  )
}

start.phia.pth <- function (parameters)
{
  tau <- logit(diag(parameters$p[,,1])[-c(1,K)]) # K-2 taus tau[1], tau[K]fixed
  eta <- rep(0,G-1)                              # G-1 etas
  c (
    log(parameters$N),                        # [1]
    plg.fn(parameters$beta),                  # [2:K]
    logit(parameters$phi[1,1:(maxage-1),1]),  # [(K+1):(K+maxage-1)]
    tau,                                      # [(K+maxage):(2*K+maxage-3)] (note p[1], p[K] fixed)
    eta,                                      # [(2*K+maxage-2):(2*K+maxage+G-4)]
    plg.fn(parameters$pi)                    # [(2*K+maxage+G-3):(2*K+maxage+2*G-5)]                  
  )
}

# ----------------
# phi (t+a) series
# ----------------

start.phita.pc <- function (parameters)
{
  tau.phi   <- logit(diag(parameters$phi[,,1])[-1])
  alpha.phi <- rep(0, maxage-1)
  c (
    log(parameters$N),                # [1]
    plg.fn(parameters$beta),          # [2:K]
    tau.phi,                          # [(K+1):(2*K-2)]
    alpha.phi,                        # [(2*K-1):(2*K+maxage-3)] 
    logit(parameters$p[1,1,1])        # [2*K+maxage-2]
  )
}

start.phita.pt <- function (parameters)
{
  tau.phi   <- logit(diag(parameters$phi[,,1])[-1])
  alpha.phi <- rep(0, maxage-1)
  c (
    log(parameters$N),                # [1]
    plg.fn(parameters$beta),          # [2:K]
    tau.phi,                          # [(K+1):(2*K-2)]
    alpha.phi,                        # [(2*K-1):(2*K+maxage-3)] 
    logit(parameters$p[1,2:(K-1),1])  # [(2*K+maxage-2):(3*K+maxage-5)]
                                      # (note p[1],p[K] fixed) 
  )
}

start.phita.pa <- function (parameters)
{
  tau.phi   <- logit(diag(parameters$phi[,,1])[-1])
  alpha.phi <- rep(0, maxage-1)
  c (
    log(parameters$N),                # [1]
    plg.fn(parameters$beta),          # [2:K]
    tau.phi,                          # [(K+1):(2*K-2)]
    alpha.phi,                        # [(2*K-1):(2*K+maxage-3)] 
    logit(parameters$p[1,2:(maxage-1),1]) # [(2*K+maxage-2):(2*K+2*maxage-5)]
                                      # (note p[1],p[maxage..K] fixed) 
  )
}

start.phita.pta <- function (parameters)
{
  tau.phi   <- logit(diag(parameters$phi[,,1])[-1])
  alpha.phi <- rep(0, maxage-1)

  tau.p     <- logit(diag(parameters$p[,,1])[2:(maxage-1)])
  #alpha.p   <- rep(0, maxage-3)
  alpha.p <- rep(0, maxage-1)

  c (
    log(parameters$N),                # [1]
    plg.fn(parameters$beta),          # [2:K]
    tau.phi,                          # [(K+1):(2*K-1)]
    alpha.phi,                        # [(2*K):(2*K+maxage-3)] 
    tau.p,                            # [(2*K+maxage-2):(3*K+maxage-1)] 
    alpha.p                           # [(3*K+maxage):(3*K+2*maxage-6)]
  )
}

start.phita.ph <- function (parameters)
{
  tau.phi   <- logit(diag(parameters$phi[,,1])[-1])  # Length K-2
  alpha.phi <- rep(0, maxage-1)                      # Length M-1
  c (
    log(parameters$N),                # [1]
    plg.fn(parameters$beta),          # [2:K]
    tau.phi,                          # [(K+1):(2*K-2)]
    alpha.phi,                        # [(2*K-1):(2*K+maxage-3)] 
    logit(parameters$p[1,1,]),        # [(2*K+maxage-2):(2*K+maxage+G-3)]
    plg.fn(parameters$pi)            # [(2*K+maxage+G-2):(2*K+maxage+2*G-4)]
  )
}

start.phita.pth <- function (parameters)
{
  tau.phi   <- logit(diag(parameters$phi[,,1])[-1])
  alpha.phi <- rep(0, maxage-1)
  tau.p <- logit(diag(parameters$p[,,1])[-c(1,K)])  # K-2 taus, tau[1], tau[K] fixed
  eta <- rep(0,G-1)                                 # G-1 etas
  
  c (
    log(parameters$N),                # [1]
    plg.fn(parameters$beta),          # [2:K]
    tau.phi,                          # [(K+1):(2*K-2)]
    alpha.phi,                        # [(2*K-1):(2*K+maxage-3)] 
    tau.p,                            # [(2*K+maxage-2):(3*K+maxage-5))
                                      # (note p[1],p[K] fixed) 
    eta,                              # [(3*K+maxage-4):(3*K+maxage+G-6)]
    plg.fn(parameters$pi)            # [(3*K+maxage+G-5):(3*K+maxage+2*G-7)]
  )
}

# ----------------
# phi (h) series
# ----------------

start.phih.pc <- function(parameters)
{
  c(  
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,]),       # [K+1:K+G]
    logit(parameters$p[1,1,1]),        # [K+G+1]
    plg.fn(parameters$pi)             # [(K+G+2):(K+2*G)]
   )
}

start.phih.pt <- function(parameters)
{
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,]),       # [K+1:K+G]
    logit(parameters$p[1,2:K,1]),      # [(K+G+1):(2*K+G-1)] (note p[1,1] fixed)
    plg.fn(parameters$pi)             # [(2*K+G):(2*K+2*G-2)]
  )
}

start.phih.pa <- function(parameters)
{
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,]),       # [K+1:K+G]
    logit(parameters$p[1,1:maxage,1]), # [(K+G+1):(K+G+maxage)] 
    plg.fn(parameters$pi)             # [(K+G+maxage+1):(K+2*G+maxage-1)]
  )
}

start.phih.pta <- function(parameters)
{
  tau   <- logit(diag(parameters$p[,,1])[-1])                 # K-1 taus tau[1] fixed
  alpha <- (logit(parameters$p[1,-1,1]) - tau[1])[2:maxage]   # alpha[1] = 0  
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,]),       # [K+1:K+G]
    tau,                               # [(K+G+1):(2*K+G-1)]  
    alpha,                             # [(2*K+G):(2*K+G+maxage-2)]
    plg.fn(parameters$pi)             # [(2*K+G+maxage-1):(2*K+2*G+maxage-3)] 
  )
}

start.phih.ph <- function(parameters)
{
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,]),       # [K+1:K+G]
    logit(parameters$p[1,1,]),         # [(K+G+1):(K+2*G)]
    plg.fn(parameters$pi)             # [(K+2*G+1):(K+3*G-1)]
  )
}

start.phih.pth <- function(parameters)
{
  tau <-logit(diag(parameters$p[,,1])[-1]) # K-1 taus tau[1] fixed
  eta <- rep(0,G-1)                        # G-1 etas
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    logit(parameters$phi[1,1,]),       # [K+1:K+G]
    tau,                               # [(K+G+1):(2*K+G-1)] 
    eta,                               # [(2*K+G):(2*K+2*G-2)]
    plg.fn(parameters$pi)             # [(2*K+2*G-1):(2*K+3*G-3)]
  )
}


# ----------------
# phi (t+h) series
# ----------------

start.phith.pc <- function(parameters)
{
  tau <-logit(diag(parameters$phi[,,1])) # K-1 taus
  eta <- rep(0,G-1)                      # G-1 etas
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    tau,                               # [(K+1):(2*K-1)]
    eta,                               # [(2*K):(2*K+G-2)]
    logit(parameters$p[1,1,1]),        # [2*K+G-1] 
    plg.fn(parameters$pi)             # [(2*K+G):(2*K+2*G-2)]
  )
}

start.phith.pt <- function(parameters)
{
  tau <-logit(diag(parameters$phi[,,1])) # K-1 taus
  eta <- rep(0,G-1)                      # G-1 etas 
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    tau,                               # [(K+1):(2*K-1]
    eta,                               # [(2*K):(2*K+G-2)]
    logit(parameters$p[1,2:(K-1),1]),  # [(2*K+G-1):(3*K+G-4)] (p1,pK fixed)
    plg.fn(parameters$pi)             # [(3*K+G-3):(3*K+2*G-5)]
  )
}

start.phith.pa <- function(parameters)
{
  tau <-logit(diag(parameters$phi[,,1]))     # K-1 taus
  eta <- rep(0,G-1)                          # G-1 etas 
  c (
    log(parameters$N),                       # [1]
    plg.fn(parameters$beta),                 # [2:K]
    tau,                                     # [(K+1):(2*K-1)]
    eta,                                     # [(2*K):(2*K+G-2)]
    logit(parameters$p[1,2:(maxage-1),1]),   # [(2*K+G-1):(2*K+M+G-4)] 
    plg.fn(parameters$pi)                   # [(2K+G+maxage-3)(2*K+2*G+maxage-5)]
  )
}

start.phith.pta <- function(parameters)
{
  tau.phi   <- logit(diag(parameters$phi[,,1]))
  eta       <- rep(0,G-1)

  tau.p     <- logit(diag(parameters$p[,,1])[-c(1,K)])
  alpha     <- rep(0,maxage-2)
 
 c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    tau.phi,                           # [(K+1):(2*K-1)]
    eta,                               # [(2*K):(2*K+G-2)]
    tau.p,                             # [(2*K+G-1):(3*K+G-4)]  
    alpha,                             # [(3*K+G-3):(3*K+maxage+G-5)]
    plg.fn(parameters$pi)             # [(3*K+G+maxage-4):(3*K+2*G+maxage-6)]
  )
}

start.phith.ph <- function(parameters)
{
  tau <-logit(diag(parameters$phi[,,1])) # K-1 taus
  eta <- rep(0,G-1)                      # G-1 etas
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    tau,                               # [(K+1):(2*K-1)]
    eta,                               # [(2*K):(2*K+G-2)]
    logit(parameters$p[1,1,]),         # [(2*K+G-1):(2*K+2*G-2)]
    plg.fn(parameters$pi)             # [(2*K+2*G-1):(2*K+3*G-3)]

  )
}

start.phith.pth <- function(parameters)
{
  tau.phi <- logit(diag(parameters$phi[,,1]))   # Length K-1
  eta.phi <- rep(0,G-1)
  tau.p   <- logit(diag(parameters$p[,,1])[-c(1,K)]) # Length K-2
  eta.p   <- rep(0,G-1)
  c (
    log(parameters$N),                 # [1]
    plg.fn(parameters$beta),           # [2:K]
    tau.phi,                           # [K+1:2*K-1]
    eta.phi,                           # [2*K:2*K+G-2]
    tau.p,                             # [(2*K+G-1):(3*K+G-4)]
    eta.p,                             # [(3*K+G-3):(3*K+2*G-5)]
    plg.fn(parameters$pi)             # [(3*K+2*G-4):(3*K+3*G-6)]
  )
}

# ----------------------
# phi (W) Weibull series
# ----------------------

start.phiW.pc <- function (parameters)
{
  # use phic.pc model 1 results
  if (is.null(parameters$coefs))  weib <- c(1, -1/log(1-parameters$phi[1,1,1]))
  else weib <- parameters$coefs
  c (
    log(parameters$N),           # [1]
    plg.fn(parameters$beta),     # [2:K]
    log(weib),                   # [(K+1):(K+2)]    kap, gam, a
    logit(parameters$p[1,1,1])   # [K+3]
  )
}

start.phiW.pt <- function (parameters)
{
  # use phic.pc model 1 results
  if (is.null(parameters$coefs))  weib <- c(1, -1/log(1-parameters$phi[1,1,1]))
  else weib <- parameters$coefs
  c (
    log(parameters$N),           # [1]
    plg.fn(parameters$beta),     # [2:K]
    log(weib),                   # [(K+1):(K+2)]
    logit(parameters$p[1,2:K,1]) # [(K+3):(2*K+1)] (note p[1,1] fixed)
  )
}

start.phiW.pa <- function (parameters)
{
  # use phic.pc model 1 results
  if (is.null(parameters$coefs))  weib <- c(1, -1/log(1-parameters$phi[1,1,1]))
  else weib <- parameters$coefs
  c (
    log(parameters$N),           # [1]
    plg.fn(parameters$beta),     # [2:K]
    log(weib),                   # [(K+1):(K+2)]
    logit(parameters$p[1,2:maxage,1]) # [(K+3):(K+maxage+1)] (note p[1,1] fixed)
  )
}

start.phiW.pta <- function (parameters)
{
  # use phic.pc model 1 results
  if (is.null(parameters$coefs))  weib <- c(1, -1/log(1-parameters$phi[1,1,1]))
  else weib <- parameters$coefs
  tau   <- logit(diag(parameters$p[,,1])[-1])
  alpha <- (logit(parameters$p[1,-1,1]) - tau[1])[1:(maxage-1)]   # alpha[1] = 0  
  c (
    log(parameters$N),            # [1]
    plg.fn(parameters$beta),      # [2:K]
    log(weib),                    # [(K+1):(K+2)]
    tau,                          # [(K+3):(2*K+1)]  
    alpha                         # [(2*K+2):(2*K+maxage)]
  )
}


start.phiW.ph <- function (parameters)
{
  # use phic.pc model 1 results
  if (is.null(parameters$coefs))  weib <- c(1, -1/log(1-parameters$phi[1,1,1]))
  else weib <- parameters$coefs
  c (
    log(parameters$N),           # [1]
    plg.fn(parameters$beta),     # [2:K]
    log(weib),                   # [(K+1):(K+2)]
    logit(parameters$p[1,1,]),   # [(K+3):(K+G+2)]
    plg.fn(parameters$pi)       # [(K+G+3):(K+2*G+1)]
  )
}


start.phiW.pth <- function (parameters)
{
  # use phic.pc model 1 results
  if (is.null(parameters$coefs))  weib <- c(1, -1/log(1-parameters$phi[1,1,1]))
  else weib <- parameters$coefs
  tau <- logit(diag(parameters$p[,,1])[-c(1,K)])
  eta <- rep(0,G-1)

  c (
    log(parameters$N),           # [1]
    plg.fn(parameters$beta),     # [2:K]
    log(weib),                   # [(K+1):(K+2)]
    tau,                         # [(K+3):(2*K)]
    eta,                         # [(2*K+1):(2*K+G-1)]
    plg.fn(parameters$pi)       # [(2*K+G):(2*K+2*G-2)]
  )
}



# -------------------------------------------------------------------------------------
# v.to.f functions
# -------------------------------------------------------------------------------------

# functions used to fill 'non-lower-triangular' cells of a matrix
# by default lower triangle is zero

# nlc  - constant
# nlrv - replicate a vector in all rows
# nlov - as with nlrv, but progressively offset the rows and truncate 'overhang'

# ----------------
# phi (c) series
# ----------------

v.to.f.phic.pc <- function (invect)
{
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlc(expit(invect[K+1]), K-1 ),
        p    = nlc(expit(invect[K+2]), K),
        pi   = c(1)
  )
}

v.to.f.phic.pt <- function (invect)
{
  p.in <- expit(invect[(K+2):(2*K)])
  p <- nlrv( c(logitmean(p.in), p.in), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlc(expit(invect[K+1]), K-1),
        p    = p,
        pi   = c(1)
  )
}

v.to.f.phic.pa <- function (invect)
{
  p.in <- expit(invect[(K+2):(K+maxage+1)])
  p    <- nlov(c(p.in,rep(p.in[length(p.in)],K-maxage)),K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlc(expit(invect[K+1]), K-1),
        p    = p,
        pi   = c(1)
  )
}

v.to.f.phic.pta <- function (invect)
{
  tau.p   <- invect[(K+2):(2*K)]
  tau.p  <- c(mean(tau.p), tau.p)
  alpha.p <- invect[(2*K+1) : (2*K+maxage-1)]
  alpha.p <- c(0, alpha.p, rep(mean(alpha.p), K-maxage))
  p       <- expit(nlrv(tau.p,K) + nlov(alpha.p,K))

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlc(expit(invect[K+1]), K-1),
        p    = p, 
        pi   = c(1)
  )
}

v.to.f.phic.ph <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlc(expit(invect[K+1]), K-1),
        p    = nlch(invect[(K+2):(K+G+1)], K, G),
        pi   = lgp.fn(invect[(K+G+2):(K+2*G)])
  )
}


v.to.f.phic.pth <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')
 
  tau <- invect[(K+2):(2*K)]
  eta <- c(invect[(2*K+1):(2*K+G-1)])
  tau <- c(mean(tau),tau)
  eta <- c(0,eta)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlc(expit(invect[K+1]), K-1),
        p    = nlth(tau, eta, K ,G),
        pi   = lgp.fn(invect[(2*K+G):(2*K+2*G-2)])
  )
}

# ----------------
# phi (t) series
# ----------------

v.to.f.phit.pc <- function (invect)
{
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlrv(expit(invect[(K+1):(2*K-1)]), K-1),
        p    = nlc(expit(invect[2*K]), K),
        pi   = c(1)
  )
}

v.to.f.phit.pt <- function (invect)
{
  p.in <- expit(invect[(2*K):(3*K-3)])
  p.in <- c(logitmean(p.in), p.in, logitmean(p.in))
  p <- nlrv(p.in, K)
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlrv(expit(invect[(K+1):(2*K-1)]), K-1),
        p    = p,
        pi   = c(1)
  )
}

v.to.f.phit.pa <- function (invect)
{
  p.in <- expit(invect[(2*K):(2*K+maxage-2)])
  p <- nlov( c(logitmean(p.in), p.in, rep(p.in[length(p.in)],K-maxage)), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlrv(expit(invect[(K+1):(2*K-1)]), K-1),
        p    = p,
        pi   = c(1)
  )
}

v.to.f.phit.pta <- function (invect)
{
  
  tau.p   <- invect[(2*K):(3*K-3)]        
  tau.p  <- c(mean(tau.p), tau.p, mean(tau.p))

  alpha.p <- invect[(3*K-2):(3*K+maxage-4)]
  alpha.p <- c( 0, alpha.p, rep(mean(alpha.p), K-maxage))   # pad to end with mean
  
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlrv(expit(invect[(K+1):(2*K-1)]), K-1),
        p    = expit(nlrv(tau.p,K) + nlov(alpha.p,K)),
        pi   = c(1)
  )
}

v.to.f.phit.ph <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlrv(expit(invect[(K+1):(2*K-1)]), K-1),
        p    = nlch(invect[(2*K):(2*K+G-1)],K,G),
        pi   = lgp.fn(invect[(2*K+G):(2*K+2*G-2)])
  )
}

v.to.f.phit.pth <- function (invect)  
{
  if (G==1) stop ('heterogeneity models require more than one group')
  tau <- invect[(2*K):(3*K-3)]
  eta <- invect[(3*K-2):(3*K+G-4)]
  tau <- c(mean(tau), tau, mean(tau))
  eta <- c(0,eta)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlrv(expit(invect[(K+1):(2*K-1)]), K-1),
        p    = nlth(tau, eta, K ,G),
        pi   = lgp.fn(invect[(3*K+G-3):(3*K+2*G-5)])
  )
}


# ----------------
# phi (a) series
# ----------------

v.to.f.phia.pc <- function (invect)
{
  phi.in <- expit(invect[(K+1):(K+maxage-1)])
  phi.in <- c(phi.in, rep(phi.in[length(phi.in)],K-maxage))
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlov(phi.in, K-1),
        p    = nlc(expit(invect[K+maxage]), K),
        pi   = c(1)
  )
}

v.to.f.phia.pt <- function (invect)
{
  phi.in <- expit(invect[(K+1):(K+maxage-1)])
  phi.in <- c(phi.in, rep(phi.in[length(phi.in)],K-maxage))
  p.in   <- expit(invect[(K+maxage):(2*K+maxage-2)])
  p      <- nlrv( c(logitmean(p.in), p.in), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlov(phi.in, K-1),
        p    = p,
        pi   = c(1)
  )
}

v.to.f.phia.pa <- function (invect)
{
  phi.in <- expit(invect[(K+1):(K+maxage-1)])
  phi.in <- c(phi.in, rep(phi.in[length(phi.in)],K-maxage))
  phi    <- nlov(phi.in, K-1) 

  p.in <- expit(invect[(K+maxage):(K+2*maxage-3)])
  p <- nlov( c(logitmean(p.in), p.in, rep(p.in[length(p.in)],K-maxage+1)), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = p,
        pi   = c(1)
  )
}

v.to.f.phia.pta <- function (invect)
{
  phi.in <- expit(invect[(K+1):(K+maxage-1)])
  phi.in <- c(phi.in, rep(phi.in[length(phi.in)], K-maxage))
  phi    <- nlov(phi.in, K-1)

  tau.p   <- invect[(K+maxage):(2*K+maxage-3)] # Length K-2
  tau.p   <- c(mean(tau.p), tau.p, mean(tau.p))
  alpha.p <- invect[(2*K+maxage-2):(2*K+2*maxage-4)]  # Length M-1
  if (K>maxage)
     alpha.p <- c(mean(alpha.p),0, alpha.p, rep(mean(alpha.p), K-maxage-1))
  if (K==maxage)
     alpha.p <- c(mean(alpha.p),0, alpha.p)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = expit(nlrv(tau.p,K) + nlov(alpha.p,K)),
        pi   = c(1)
  )
}

# New version:

v.to.f.phia.pta <- function (invect)
{
  phi.in <- expit(invect[(K+1):(K+maxage-1)])
  phi.in <- c(phi.in, rep(phi.in[length(phi.in)], K-maxage))
  phi    <- nlov(phi.in, K-1)

  tau.p   <- invect[(K+maxage):(2*K+maxage-3)] # Length K-2
  tau.p   <- c(mean(tau.p), tau.p, mean(tau.p))
  alpha.p <- invect[(2*K+maxage-2):(2*K+2*maxage-4)]  # Length M-1
  if (K>maxage)
     alpha.p <- c(0, alpha.p, rep(mean(alpha.p), K-maxage))
  if (K==maxage)
     alpha.p <- c(0, alpha.p)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = expit(nlrv(tau.p,K) + nlov(alpha.p,K)),
        pi   = c(1)
  )
}

v.to.f.phia.ph <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  phi.in <- expit(invect[(K+1):(K+maxage-1)])
  phi.in <- c(phi.in, rep(phi.in[length(phi.in)],K-maxage))

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlov(phi.in, K-1),
        p    = nlch(invect[(K+maxage):(K+maxage+G-1)],K,G),
        pi   = lgp.fn(invect[(K+maxage+G):(K+maxage+2*G-2)])
  )
}
v.to.f.phia.pth <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  phi.in <- expit(invect[(K+1):(K+maxage-1)])
  phi.in <- c(phi.in, rep(phi.in[length(phi.in)],K-maxage))

  tau <- invect[(K+maxage):(2*K+maxage-3)]
  eta <- invect[(2*K+maxage-2):(2*K+maxage+G-4)]
  tau <- c(mean(tau) ,tau ,mean(tau))
  eta <- c(0,eta)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlov(phi.in, K-1),
        p    = nlth(tau, eta, K ,G),
        pi   = lgp.fn(invect[(2*K+maxage+G-3):(2*K+maxage+2*G-5)])
  )
}


# ----------------
# phi (t+a) series
# ----------------

v.to.f.phita.pc <- function (invect)
{
  tau.phi   <- invect[(K+1):(2*K-2)]                  # Length K-2 
  tau.phi   <- c(mean(tau.phi),tau.phi)               # Length K-1
  alpha.phi <- invect[(2*K-1):(2*K+maxage-3)]         # Length M-1
  if (K>maxage)
     alpha.phi <- c(0,alpha.phi, rep(mean(alpha.phi),K-maxage-1))
                                                      # Length K-1
  if (K==maxage)
     alpha.phi <- c(0,alpha.phi)
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = expit(nlrv(tau.phi,K-1) + nlov(alpha.phi,K-1)),
        p    = nlc(expit(invect[2*K+maxage-2]), K),
        pi   = c(1)
  )
}

# New version:

v.to.f.phita.pc <- function (invect)
{
  tau.phi   <- invect[(K+1):(2*K-1)]                  # Length K-1 
  alpha.phi <- invect[(2*K):(2*K+maxage-3)]           # Length M-2
  alpha.phi <- c(0,alpha.phi, rep(mean(alpha.phi),K-maxage))
                                                      # Length K-1
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = expit(nlrv(tau.phi,K-1) + nlov(alpha.phi,K-1)),
        p    = nlc(expit(invect[2*K+maxage-2]), K),
        pi   = c(1)
  )
}

v.to.f.phita.pt <- function (invect)
{
  tau.phi   <- invect[(K+1):(2*K-1)]                  # Length K-1 
  alpha.phi <- invect[(2*K):(2*K+maxage-3)]           # Length M-2
  alpha.phi <- c(0,alpha.phi, rep(mean(alpha.phi),K-maxage))
                                                      # Length K-1
  p.in <- expit(invect[(2*K+maxage-2):(3*K+maxage-5)])
  p <- nlrv( c(logitmean(p.in), p.in, logitmean(p.in)), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = expit(nlrv(tau.phi,K-1) + nlov(alpha.phi,K-1)),
        p    = p,
        pi   = c(1)
     )
}

v.to.f.phita.pa <- function (invect)
{
  tau.phi   <- invect[(K+1):(2*K-1)]                  # Length K-1 
  alpha.phi <- invect[(2*K):(2*K+maxage-3)]           # Length M-2
  alpha.phi <- c(0,alpha.phi, rep(mean(alpha.phi),K-maxage))
                                                      # Length K-1
  p.in <- expit(invect[(2*K+maxage-2):(2*K+2*maxage-5)])
  p <- nlov( c(logitmean(p.in), p.in, rep(p.in[length(p.in)],K-maxage+1)), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = expit(nlrv(tau.phi,K-1) + nlov(alpha.phi,K-1)),
        p    = p,
        pi   = c(1)
  )
}

v.to.f.phita.pta <- function (invect)
{
  tau.phi   <- invect[(K+1):(2*K-1)]                  # Length K-1 
  alpha.phi <- invect[(2*K):(2*K+maxage-3)]           # Length M-2
  alpha.phi <- c(0,alpha.phi, rep(mean(alpha.phi),K-maxage))
                                                      # Length K-1
  tau.p   <- invect[(2*K+maxage-2):(3*K+maxage-5)]
  tau.p  <- c(mean(tau.p), tau.p, (mean(tau.p)))

  alpha.p <- invect[(3*K+maxage-4):(3*K+2*maxage-6)]
  if (K>maxage)
     alpha.p <- c(mean(alpha.p), 0, alpha.p, rep(mean(alpha.p),K-maxage-1))
  if (K==maxage)
     alpha.p <- c(mean(alpha.p), 0, alpha.p)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = expit(nlrv(tau.phi,K-1) + nlov(alpha.phi,K-1)),
        p    = expit(nlrv(tau.p,K) + nlov(alpha.p,K)),
        pi   = c(1)
  )
}

v.to.f.phita.ph <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau.phi   <- invect[(K+1):(2*K-2)]                  # Length K-2 
  tau.phi   <- c(mean(tau.phi),tau.phi)               # Length K-1
  alpha.phi <- invect[(2*K-1):(2*K+maxage-3)]         # Length M-1
  if (K>maxage)
     alpha.phi <- c(0,alpha.phi, rep(mean(alpha.phi),K-maxage-1))
                                                      # Length K-1
  if (K==maxage)
     alpha.phi <- c(0,alpha.phi)
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = expit(nlrv(tau.phi,K-1) + nlov(alpha.phi,K-1)),
        p    = nlch(invect[(2*K+maxage-2):(2*K+maxage+G-3)],K,G),
        pi   = lgp.fn(invect[(2*K+maxage+G-2):(2*K+maxage+2*G-4)])
  )
}

v.to.f.phita.pth <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau.phi   <- invect[(K+1):(2*K-2)]                  # Length K-2 
  tau.phi   <- c(mean(tau.phi),tau.phi)               # Length K-1
  alpha.phi <- invect[(2*K-1):(2*K+maxage-3)]         # Length M-1
  if (K>maxage)
    alpha.phi <- c(0,alpha.phi, rep(mean(alpha.phi),K-maxage-1))
                                                      # Length K-1
  if (K==maxage)
    alpha.phi <- c(0,alpha.phi)
  tau.p <- invect[(2*K+maxage-2):(3*K+maxage-5)]
  eta   <- invect[(3*K+maxage-4):(3*K+maxage+G-6)]
  tau.p <- c(mean(tau.p) ,tau.p ,mean(tau.p))
  eta   <- c(0,eta)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = expit(nlrv(tau.phi,K-1) + nlov(alpha.phi,K-1)),
        p    = nlth(tau.p, eta, K ,G),
        pi   = lgp.fn(invect[(3*K+maxage+G-5):(3*K+maxage+2*G-7)])
  )
}

# ----------------
# phi (h) series
# ----------------

v.to.f.phih.pc <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlphi(invect[(K+1):(K+G)],K-1,G),
        p    = nlch(expit(rep(invect[K+G+1],G)),K,G),
        pi   = lgp.fn(invect[(K+G+2):(K+2*G)])
  )
}

v.to.f.phih.pt <- function (invect)
{
 if (G==1) stop ('heterogeneity models require more than one group')

 p.in <- expit(invect[(K+G+1):(2*K+G-1)])
 p <- nlrv( c(logitmean(p.in), p.in), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlphi(invect[(K+1):(K+G)],K-1,G),
        p    = p,
        pi   = lgp.fn(invect[(2*K+G):(2*K+2*G-2)])
  )
}

v.to.f.phih.pa <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  p.in <- expit(invect[(K+G+1):(K+G+maxage)])
  p <- nlov( c( p.in, rep(p.in[length(p.in)],K-maxage)), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlphi(invect[(K+1):(K+G)],K-1,G),
        p    = p,
        pi   = lgp.fn(invect[(K+maxage+G+1):(K+2*G+maxage-1)])
  )
}

v.to.f.phih.pta <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau.p   <- invect[(K+G+1):(2*K+G-1)] 
  tau.p   <- c(mean(tau.p), tau.p)

  alpha.p <- invect[(2*K+G):(2*K+G+maxage-2)]
  alpha.p <- c(0, alpha.p, rep(mean(alpha.p), K-maxage))
  p       <- expit(nlrv(tau.p,K) + nlov(alpha.p,K))
  
  for (g in 1:G)  p[,,g][lower.tri(p[,,g])] <- 0   # for tidiness

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlphi(invect[(K+1):(K+G)],K-1,G),
        p    = p, 
        pi   = lgp.fn(invect[(2*K+G+maxage-1):(2*K+2*G+maxage-3)])
  )
}

v.to.f.phih.ph <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlphi(invect[(K+1):(K+G)],K-1,G),
        p    = nlch(invect[(K+G+1):(K+2+G)],K,G),
        pi   = lgp.fn(invect[(K+2*G+1):(K+3*G-1)])
  )
}

v.to.f.phih.pth <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau <- invect[(K+G+1):(2*K+G-1)]
  eta <- invect[(2*K+G):(2*K+2*G-2)]
  tau <- c(mean(tau),tau)
  eta <- c(0,eta)
     
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlphi(invect[(K+1):(K+G)],K-1,G),
        p    = nlth(tau, eta, K ,G),
        pi   = lgp.fn(invect[(2*K+2*G-1):(2*K+3*G-3)])
  )
}

# -----------------
# phi (t+h) series
# -----------------

v.to.f.phith.pc <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau <- invect[(K+1):(2*K-1)]
  eta <- invect[(2*K):(2*K+G-2)]
  eta <- c(0,eta)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlth(tau, eta, K-1 ,G),
        p    = nlch(rep(invect[2*K+G-1],G), K,G),
        pi   = lgp.fn(invect[(2*K+G):(2*K+2*G-2)])
  )
}

v.to.f.phith.pt <- function (invect) 
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau <- invect[(K+1):(2*K-1)]
  eta <- invect[(2*K):(2*K+G-2)]
  eta <- c(0,eta)
  p.in <- invect[(2*K+G-1):(3*K+G-4)]
  p.in <- c(mean(p.in),p.in,mean(p.in))
  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlth(tau, eta, K-1 ,G),
        p    = nlrv(expit(p.in), K),
        pi   = lgp.fn(invect[(3*K+G-3):(3*K+2*G-5)])
  )
}

v.to.f.phith.pa <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau <- invect[(K+1):(2*K-1)]
  eta <- invect[(2*K):(2*K+G-2)]
  eta <- c(0,eta)

  p.in <- expit(invect[(2*K+G-1):(2*K+G+maxage-4)])
  p.in <- c(logitmean(p.in),p.in,logitmean(p.in), rep(p.in[length(p.in)],K-maxage))

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlth(tau, eta, K-1 ,G),
        p    = nlov(p.in, K),
        pi   = lgp.fn(invect[(2*K+G+maxage-3):(2*K+2*G+maxage-5)])
  )
}

v.to.f.phith.pta <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau.phi <- invect[(K+1):(2*K-1)]
  eta <- invect[(2*K):(2*K+G-2)]
  eta <- c(0,eta)

  tau.p <- invect[(2*K+G-1):(3*K+G-4)]          # Length K-2
  tau.p <- c(mean(tau.p),tau.p,mean(tau.p))
  alpha <- invect[(3*K+G-3):(3*K+G+maxage-6)]   # Length M-2
  alpha <- c(0, alpha, mean(alpha), rep(mean(alpha),K-maxage))

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlth(tau.phi, eta, K-1 ,G),
        p    = expit(nlrv(tau.p,K) + nlov(alpha,K)),
        pi   = lgp.fn(invect[(3*K+G+maxage-5):(3*K+2*G+maxage-7)])
  )
}


v.to.f.phith.ph <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau <- invect[(K+1):(2*K-1)]
  eta <- invect[(2*K):(2*K+G-2)]
  eta <- c(0,eta)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlth(tau, eta, K-1 ,G),
        p    = nlch(invect[(2*K+G-1):(2*K+2*G-2)],K,G),
        pi   = lgp.fn(invect[(2*K+2*G-1):(2*K+3*G-3)])
  )
}

v.to.f.phith.pth <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')

  tau.phi <- invect[(K+1):(2*K-1)]
  eta.phi <- invect[(2*K):(2*K+G-2)]
  eta.phi <- c(0,eta.phi)

  tau.p <- invect[(2*K+G-1):(3*K+G-4)]
  eta.p <- invect[(3*K+G-3):(3*K+2*G-5)]
  eta.p <- c(0,eta.p)
  tau.p <- c(mean(tau.p),tau.p,mean(tau.p))

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = nlth(tau.phi, eta.phi, K-1 ,G),
        p    = nlth(tau.p, eta.p, K ,G),
        pi   = lgp.fn(invect[(3*K+2*G-4):(3*K+3*G-6)])
  )
}


# ----------------
# phi (W) series
# ----------------

v.to.f.phiW.pc <- function (invect)
{
  kap <- exp(invect[K+1])
  gam <- exp(invect[K+2])
  phi   <- exp( -(age2/gam)^kap + (age1/gam)^kap )
  phi[lower.tri(phi)] <- 0 # for tidiness

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = nlc(expit(invect[K+3]), K),
        pi   = c(1),
        coefs = c(kap, gam)
  )
}

v.to.f.phiW.pt <- function (invect)
{
  kap <- exp(invect[K+1])
  gam <- exp(invect[K+2])
  phi <- exp( -(age2/gam)^kap + (age1/gam)^kap )
  phi[lower.tri(phi)] <- 0 # for tidiness
  p.in <- expit(invect[(K+3):(2*K+1)])
  p <- nlrv( c(logitmean(p.in), p.in), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = p,
        pi   = c(1),
        coefs = c(kap, gam)
  )
}

v.to.f.phiW.pa <- function (invect)
{
  kap <- exp(invect[K+1])
  gam <- exp(invect[K+2])
  phi <- exp( -(age2/gam)^kap + (age1/gam)^kap )
  phi[lower.tri(phi)] <- 0 # for tidiness
  p.in <- expit(invect[(K+3):(K+maxage+1)])
  p    <- nlov( c(logitmean(p.in), p.in, rep(logitmean(p.in),K-maxage)), K)

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = p,
        pi   = c(1),
        coefs = c(kap, gam)
  )
}

v.to.f.phiW.pta <- function (invect)
{
  kap <- exp(invect[K+1])
  gam <- exp(invect[K+2])
  phi <- exp( -(age2/gam)^kap + (age1/gam)^kap )
  phi[lower.tri(phi)] <- 0 # for tidiness
  tau.p   <- invect[(K+3):(2*K+1)] 
  tau.p   <- c(mean(tau.p), tau.p)
  alpha.p <- invect[(2*K+2) : (2*K+maxage)]
  alpha.p <- c(0, alpha.p, rep(mean(alpha.p), K-maxage))

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = expit(nlrv(tau.p,K) + nlov(alpha.p,K)),
        pi   = c(1),
        coefs = c(kap, gam)
  )
}

v.to.f.phiW.ph <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')
  kap  <- exp(invect[K+1])
  gam  <- exp(invect[K+2])
  phi.in <- exp( -(age2/gam)^kap + (age1/gam)^kap )
  phi.in[lower.tri(phi.in)] <- 0 # for tidiness
  phi <- array(NA,c(K-1,K-1,G))
  for (gg in 1:G) phi[,,gg] <- phi.in
  p.in <- expit(invect[(K+3):(K+G+2)])

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = nlch(p.in,K,G),
        pi   = lgp.fn(invect[(K+G+3):(K+2*G+1)]) ,
        coefs = c(kap, gam)
  )
}

v.to.f.phiW.pth <- function (invect)
{
  if (G==1) stop ('heterogeneity models require more than one group')
  kap  <- exp(invect[K+1])
  gam  <- exp(invect[K+2])
  phi.in <- exp( -(age2/gam)^kap + (age1/gam)^kap )
  phi.in[lower.tri(phi.in)] <- 0 # for tidiness
  phi <- array(NA,c(K-1,K-1,G))
  for (gg in 1:G) phi[,,gg] <- phi.in
  tau <- invect[(K+3):(2*K)]
  eta <- invect[(2*K+1):(2*K+G-1)]
  eta <- c(0,eta)
  tau <- c(mean(tau),tau,mean(tau))
#print(c("K",K,"tau ",length(tau),"eta",length(eta),"g",g))

  list (N    = exp(invect[1]),
        beta = lgp.fn(invect[2:K]),
        phi  = phi,
        p    = nlth(tau, eta, K ,G),
        pi   = lgp.fn(invect[(2*K+G):(2*K+2*G-2)]) ,
        coefs = c(kap, gam)
  )
}
