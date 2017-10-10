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
diffrec = c(2014)

# Specify if covariates will be used:
# (TRUE or FALSE)
Covars = FALSE

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
    di = matrix(0, n, 2)
}
nz = ncol(Z)

# Extract times of birth and death:
dd = read.csv("raw.csv")
dd[,4] = as.Date(dd[,4], format = "%m/%d/%Y")
dd = na.omit(dd[, c(3:4, 19)])
dd = dd[which(dd[,1] != ""),]
dd[,1] = as.numeric(dd[,1])
bd = dd[which(dd[,3] == "J" | dd[,3] == "Y"),1:2]
bd = cbind(bd[,1], as.numeric(format(bd[,2], "%Y")) - 1)

temp = vector()
count = 1

for (i in 2:(length(bd[, 1]))) {
  if (bd[i, 1] == bd[i - 1, 1]) {
    temp[count] = i
    count = count + 1  
    }
  else{
    
  }
}

bi = bd[-temp,]
colnames(bi) = c("ID", "bi")
colnames(dd) = c("ID", "Date", "Age")

dd = merge(dd, bi, by = "ID", all = T)[,c(1,4)]
bi = unique(dd)
bi[is.na(bi)] <- 0

# Define study duration:
Ti = 2010
Tf = 2017
st = Ti:Tf
nt = length(Ti:Tf)
Dx = (st[2]-st[1])
Tm = matrix(st, n, nt, byrow=TRUE)

# Calculate first and last time observed:
ytemp = t(t(Y[,2:9]) * st)
li = c(apply(ytemp,1,max))
ytemp[ytemp==0] = 10000
fi = c(apply(ytemp,1,min))
fi[fi==10000] = 0
rm("ytemp")

# Calculate number of times detected:
oi = as.matrix(Y[,2:9]) %*% rep(1, nt)

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



inputMat <- as.data.frame(cbind(ID = 1:754, Birth = bi[,2], Death = di[,1], Y[,-1],  Z))

# analysis

# check data
newData <- DataCheck(inputMat, studyStart = 2010,
                      studyEnd = 2017, 
                      silent = FALSE)

inputMat[325,2] = inputMat[325,2]-1 
inputMat[325,6] = 0 

ni = 60000
nt = 100
nc = 4
nb = 10000

out <- basta(object = inputMat[,-12], studyStart = 2010, studyEnd = 2017,
             model = "GO", shape = "simple",
             nsim = nc, niter = ni, burnin = nb, thinning = nt, ncpus = 4, parallel = T)



summary(out, digits = 3)
plot(out)
plot(out, plot.trace = FALSE, fancy = T)





#--------------------
# mortality functions
#--------------------
# Exponential
par(mfrow = c(3,4))

a = 1:10
b = seq(0.001, 1, length.out = 10)    # 0.45 < b < 0.65 reasonable

for (i in 1:length(b)) {
  plot(a, rep(b[i], length(a)))
  plot(a, exp(-rep(b[i], length(a)) * a))
}

#--------------------
# Gompertz
a = 1:10
b0 = seq(-2, 2, length.out = 10)      # -1.1 reasonable
b1 = seq(0.001, 1, length.out = 10)   # < 0.1 reasonable 

par(mfrow = c(3,4))

for (i in 1:length(b0)) {
  plot(a, exp(b0[i] + b1[i] * a))
  plot(a, exp((exp(b0[i])/b1[i]) * (1 - exp(b1[i] * a))))
}

#-------------------
# Weibull
a = 1:10
b0 = seq(1, 10, length.out = 10)     # 1 < b < 4 reasonable
b1 = seq(0.001, 1, length.out = 10)  # 0.2 reasonable

par(mfrow = c(3,4))

for (i in 1:length(b0)) {
  plot(a, ((b0[i] * b1[i]) * (b1[i] * a) ^ (b0[i] - 1)))
  plot(a, exp(-(b1[i] * a) ^ b0[i]))
}

#------------------
# Logistic
a = 1:10
b0 = seq(0.1, 3, length.out = 10)
b1 = seq(0.1, 3, length.out = 10)
b2 = seq(0.1, 5, length.out = 10)

par(mfrow = c(3,4))

for (i in 1:length(b0)) {
  plot(a, (exp(b0[i] + (b1[i] * a)) / ((1 + (b2[i] * (exp(b0[i]) / b1[i])) * (exp(b1[i] * a) - 1)))))
  plot(a, ((1 + b2[i] * exp(b0[i]) / b1[i]) * (exp(b1[i] * a) - 1)) ^ -(1 / b2[i]))
}








