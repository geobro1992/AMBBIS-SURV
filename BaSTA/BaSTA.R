# PACKAGES:
library(msm)
library(RColorBrewer)
library(snow)
library(BaSTA)

setwd("BaSTA")

# USER INPUT:
# Specify model to be used:
# (GO = Gompertz, GM = Gompertz-Makeham, SI = Siler)
model = "SI"

# Years with different recapture prob. (pi):
# (Specify the start of the intervals)
diffrec = c(1995, 1999, 2002)

# Specify if covariates will be used:
# (TRUE or FALSE)
Covars = TRUE

# Model variables to be adjusted
# DEFINE PRIORS, STARTING PARAMETERS AND JUMP SDs:
# Survival parameters:
thp = c(-5,0.1,-1,0.001,0.005)

# Jump sd's for survival parameters:
thj = c(0.005, 0.005, 0.02, 0.0075, 0.001)

# Starting values for survival parameters:
thg = c(-1, 0.001, 0, -1, 0.001)


# DATA PREP.:
# Import data:
bd = as.matrix(read.csv("bd.csv", sep = ","),
                            header = TRUE, row.names = NULL)

dd = read.csv("raw.csv")                
dd[,4] = as.Date(dd[,4], format = "%m/%d/%Y")
dd = na.omit(dd[, 3:4])
dd = dd[which(dd[,1] != ""),]
Y = CensusToCaptHist(ID = dd[, 1], d = dd[, 2])


n = nrow(Y)
if(Covars){
Z = read.csv("Z.csv", sep=",", header=TRUE)
} else {
    Z = matrix(1, n, 1)
}
nz = ncol(Z)

# Extract times of birth and death:
bi = bd[,1]
di = bd[,2]

# Define study duration:
Ti = 1995
Tf = 2002
st = Ti:Tf
nt = length(Ti:Tf)
Dx = (st[2]-st[1])
Tm = matrix(st, n, nt, byrow=TRUE)

# Calculate first and last time observed:
ytemp = t(t(Y) * st)
li = c(apply(ytemp,1,max))
ytemp[ytemp==0] = 10000
fi = c(apply(ytemp,1,min))
fi[fi==10000] = 0
rm("ytemp")

# Calculate number of times detected:
oi = Y %*% rep(1, nt)

# MODELS
nth = 5         # number of parameters in mortality function
modm = matrix(1, 3, nth)    # matrix to indicate whether the parameters are used
modm[1,1:3] = 0
modm[2,1:2] = 0
dimnames(modm) = list(c("GO", "GM", "SI"),
                       c("alpha1", "beta1","c","alpha2","beta2"))
pname = paste(rep(colnames(modm),each=nz),
                "[",rep(colnames(Z), nth),"]", sep="")
idm = which(rownames(modm)==model)    # choose the mortality function
idth = which(modm[rep(idm, nz),]==1)  # which parameters will be used, given mortality function

# FUNCTIONS:
# Survival and mortality:

# Hazard function
m.g = function(x,th) exp(th[,1]-th[,2]*x) * modm[idm,1]*modm[idm,2] +
      th[,3] + exp(th[,4] + th[,5]*x)
# Survivorship function
S.g = function(x, th){
      Sg = exp((exp(th[,1])/th[,2] *
          (exp(-th[,2]*x)-1)) * modm[idm,1]*modm[idm,2] +
           x * (-th[,3]) +
           exp(th[,4])/th[,5] * (1-exp(th[,5]*x)))
  return(Sg)
}
# PDF of ages at death
f.g = function(x,th) m.g(x,th) * S.g(x,th)
S.x = function(th) S.g(xv, matrix(th,1,nth))
m.x = function(th) m.g(xv, matrix(th,1,nth))

# Lower bounds for parameter c:
c.low = function(th){
if(idm==1) cl = 0
if(idm==2) if(th[5]>0) cl = -exp(th[4]) else if(th[5]<0) cl = 0
if(idm==3){
    x.minf = (th[1]+log(th[2]) - th[4]-log(th[5]))/(th[2] + th[5])
    cl = -exp(th[1]-th[2]*(x.minf)) - exp(th[4]+th[5]*(x.minf))
  }
  return(cl)
}


th.low = matrix(-Inf, nrow(modm),ncol(modm), dimnames=dimnames(modm))
th.low["SI",c("beta1","beta2")] = 0
low = matrix(th.low[idm,], nz, nth, byrow=TRUE)
dimnames(low) = list(colnames(Z), colnames(modm))

# Observation matrices:
ObsMatFun = function(f, l){
    Fm = Tm - f; Fm[Fm>=0] = 1; Fm[Fm<0] = 0
    Lm = Tm - l; Lm[Lm<=0] = -1; Lm[Lm>0] = 0
    return(Fm * (-Lm))
}



inputMat <- as.data.frame(cbind(1:15, bi, di, Y, Z))

# analysis
install.packages("BaSTA")
library(BaSTA)

# check data
newData <- DataCheck(inputMat, studyStart = 1995,
                      studyEnd = 2002, autofix = rep(1, 7),
                      silent = FALSE)


library(snow)
out <- basta(object = inputMat, studyStart = 1995, studyEnd = 2002, nsim = 4, parallel = T)
summary(out, digits = 3)
plot(out)
plot(out, plot.trace = FALSE, fancy = T)


#------------------
# With Real Data!!!
#------------------
# USER INPUT:
# Specify model to be used:
# (GO = Gompertz, GM = Gompertz-Makeham, SI = Siler)
model = "SI"

# Specify if covariates will be used:
# (TRUE or FALSE)
Covars = TRUE

setwd("C:/Users/boa10gb/Documents/R/BaSTA")

# Model variables to be adjusted
# DEFINE PRIORS, STARTING PARAMETERS AND JUMP SDs:
# Survival parameters:
thp = c(-5,0.1,-1,0.001,0.005)

# Jump sd's for survival parameters:
thj = c(0.005, 0.005, 0.02, 0.0075, 0.001)

# Starting values for survival parameters:
thg = c(-1, 0.001, 0, -1, 0.001)

# PACKAGES:
library(msm)
library(RColorBrewer)
library(snow)

# DATA PREP.:
# Import data:
df = read.csv("ByCatch.csv")
df$Date = as.Date(df$Date, format = "%m/%d/%Y")
# get rid of incomplete records
df = na.omit(df)


# create year column
tick = as.Date(c("2010-06-15",
                 "2011-06-15",
                 "2012-06-15", 
                 "2013-06-15", 
                 "2014-06-15", 
                 "2015-06-15",
                 "2016-06-15"))

year = vector()
for(i in 1:(length(tick)-1)){
  x = length(subset(all.data$Date, all.data$Date > tick[i] & all.data$Date < tick[i+1]))
  year = append(year, rep(i, x)) 
}


bd = vector()
for(i in 1:length(dd)){
  if(dd$Age[i] = Y){
    bd[i] = dd$Year - 1
  }
}

dd = rep(0, length(dd))


Y = as.matrix(read.csv("Y.csv", sep=","),
              header=F, row.names=NULL)
colnames(Y) = 1995:2002


n = nrow(Y)
if(Covars){
  Z = read.csv("Z.csv", sep=",", header=TRUE)
} else {
  Z = matrix(1, n, 1)
}
nz = ncol(Z)
