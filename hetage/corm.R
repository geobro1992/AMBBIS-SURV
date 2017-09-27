# corm.R

# Commands for analysis of cormorant data, 317 birds, 8 months.

# Load the functions and the data set:

source("hetage.functions.R")
dyn.load("./hetageLL.dll")

corm.df <- read.csv("corm_1994_BPall.csv")

# Check out the appearance of the data set:

str(corm.df)
dim(corm.df)          # 317 birds by 9 sampling occasions
apply(corm.df,2,sum)  # Any leading or trailing zeros?
                      # If so, remove those columns.

# There is no column labelled "freq", so it is assumed
# that all frequencies are 1. In this case, data will be 
# condensed if there are repeated rows.

hetage.process.data(corm.df)

# View the working data:

x.mat         # Only 50 distinct capture histories
y.vect        # Frequencies of the capture histories
sum(y.vect)   # Should be 317.

# Name the models to be fitted:

model.names <- c("phic.pc","phic.pt","phic.pa","phic.pta",
                 "phit.pc","phit.pt","phit.pa","phit.pta",
                 "phia.pc","phia.pt","phia.pa","phia.pta",
                 "phita.pc","phita.pt","phita.pa",
                 "phiW.pc","phiW.pt","phiW.pa","phiW.pta") 
 
# Set up the list:

no.models   <- length(model.names)
model.fits  <- vector("list",no.models)
names(model.fits) <- model.names

# Fit the models:

for (modno in 1:no.models)
   {
   model.fits[[modno]] <- hetage.fit.model(model.names[modno])
   print(paste("Model",model.names[modno],"fitted"))
   }

# Trouble with a model fit?
# Error message including  
#   " L-BFGS-B needs finite values of 'fn' "?
# Default starts may be poor.
# Try starting the optimisation from the output parameters
# of a successful model of a similar type. For example:

model.fits[[17]] <- hetage.fit.model(model.names[17],
   start=model.fits[[16]]$parameters)


# Detailed output information may be found using

hetage.summary(model.fits)             # for a detailed summary table
hetage.summary(model.fits,AICorder=T)  # summary table sorted by AIC
hetage.summary(model.fits,AICcorder=T) # summary table sorted by AICc

phita.pc.out <- model.fits[[13]]       # Best model (AICc)

hetage.stopover(phita.pt.out)  # for mean and sd of stopover duration

hetage.superpop(phita.pt.out)  # for mean and sd of superpopulation 

# Save the results:

corm.out <- model.fits


# We may plot some results:

# Plot 1:

plot(1:K,phita.pc.out$parameters$beta,xlab="Sample",ylab="Proportion of N",
     main="Proportions entering just before each sample")

# Plot 2:

this.mat <- phita.pt.out$parameters$phi[,,1]
plot(1:maxage,this.mat[1,],xlab="Sample",
     ylab="Retention probability",pch="1",ylim=c(0,1),
     main="Probability of staying")
lines(1:maxage,this.mat[1,])
for (j in 2:maxage)
   points(j:maxage,this.mat[j,-(1:(j-1))],pch=as.character(j))
for (j in 2:(maxage-1))
   lines(j:maxage,this.mat[j,-(1:(j-1))])
text(3,0,"Numbers show cohort")


# Practical Session 4.
# --------------------

# Last plot and some calculations for a Weibull model.

# Plot 3:

# Even though model 19 was not chosen, we may look at the results:

phiW.pta.out <- model.fits[[19]]
model.pars  <- phiW.pta.out$parameters$coefs
plot (x=1:(K-1), y=Weibull (1:(K-1), aa=model.pars[1], bb=model.pars[2]),
      ylim=c(0,1), xlab='Duration of stay (months)', ylab='Survival, phi',
      main="Weibull Survival Curve")
lines(1:(K-1), Weibull (1:(K-1),aa=model.pars[1], bb=model.pars[2]))

# Find mean and s.d. for Weibull distribution of survival times.

names(model.pars) <- NULL
Weibull.mean <- model.pars[2]*gamma(1+1/model.pars[1])
Weibull.var  <- model.pars[2]^2*gamma(1+2/model.pars[1]) - Weibull.mean^2
c(mean=Weibull.mean,s.d.=sqrt(Weibull.var))


