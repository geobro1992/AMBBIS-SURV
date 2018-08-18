# PACKAGES:
library(msm)
library(RColorBrewer)
library(snow)
library(BaSTA)

setwd("BaSTA")

# USER INPUT:
# Specify model to be used:
# (GO = Gompertz, GM = Gompertz-Makeham, SI = Siler)
model = "EX"

# Years with different recapture prob. (pi):
# (Specify the start of the intervals)
diffrec = c(2010:2017)

# Specify if covariates will be used:
# (TRUE or FALSE)
Covars = F

# Model variables to be adjusted
# DEFINE PRIORS, STARTING PARAMETERS AND JUMP SDs:
# Survival parameters:
thp = c(-5,0.1,-1,0.001,0.005)

# Starting values for survival parameters:
thg = c(-1, 0.001, 0, -1, 0.001)

# Jump sd's for survival parameters:
thj = c(0.005, 0.005, 0.02, 0.0075, 0.001)




# DATA PREP.:
# Import data:
r.data = read.csv("raw.csv")                
r.data[, 4] = as.Date(r.data[, 4], format = "%m/%d/%Y")

# Assign codes to indv without Master IDs
r.data = na.omit(r.data[, c(3:4, 12, 19)])
r.data[, 1] = as.character(r.data[, 1])

count = 2000

for (i in 1:length(r.data[,1])) {
  if (r.data[i, 1] == "") {
    r.data[i, 1] = count
    count = count + 1
  }
}

colnames(r.data) = c("ID", "Date", "Pond","Age")
r.data[,1] = as.factor(r.data[,1])

# Known birth times
kb = r.data[which(r.data[, 4] == "J" | r.data[, 4] == "Y" | r.data[, 4] == "M"), 1:2]
ys = as.numeric(format(kb[, 2], "%Y")) - 1
kb = cbind(kb, ys)

# remove duplicates birth years
temp = vector()
count = 1

for (i in 2:(length(kb[, 1]))) {
  if (kb[i, 1] == kb[i - 1, 1]) {
    temp[count] = i
    count = count + 1  
  }
  else{
  }
}

kb = kb[-temp, ]
colnames(kb) = c("ID", "Date","bi")

# merge data with birth dates
df = merge(r.data, kb, by = "ID", all = T)[, c(1, 2, 3, 6)]

# birth times for all individuals
bi = df[!duplicated(df$ID), ]
bi[is.na(bi)] <- 0

# pond id for all individuals
#pid = unique(df[c(1,3)])
#pid = pid[-which(duplicated(pid[,1]) == T),]

# cap history
Y = CensusToCaptHist(ID = df[,1], d = df[,2])
n = nrow(Y)

#Y = Y[-385,] # remove pond 53 individual to include group effect
#di = di[-385,]
#Z = Z[-385,]
#i = bi[-385,] # remove pond 53 individual to include group effect



# If no covariates create Z column of 1s
#Z = matrix(1, n, 1)

# If unknown deaths create di column on 0s
di = matrix(0, n, 2)


# get pond ids for group effect
pid = df[!duplicated(df$ID), ]
n_occur <- data.frame(table(pid[,1]))
n_occur[n_occur$Freq > 1,] # find duplicates
pid[, 3] = as.factor(as.numeric(pid[, 3]))
pid = data.frame(pid[,3])
colnames(pid) = "Pond"
Z = MakeCovMat(~ Pond, data = pid)


inputMat <- as.data.frame(cbind(ID = 1:1513, Birth = bi[,4], Death = di[,1], Y[,-1],  Z = Z[,2:4]))
inputMat = inputMat[-which(inputMat$Z.Pond4 > 0), ]
# Define study duration:
Ti = 2010
Tf = 2017
st = Ti:Tf
nt = length(Ti:Tf)
Dx = (st[2] - st[1])
Tm = matrix(st, n, nt, byrow = TRUE)

# Calculate first and last time observed:
ytemp = t(t(Y[, 2:9]) * st)
li = c(apply(ytemp, 1, max))
ytemp[ytemp == 0] = 10000
fi = c(apply(ytemp, 1, min))
fi[fi == 10000] = 0
rm("ytemp")

# Calculate number of times detected:
oi = as.matrix(Y[, 2:9]) %*% rep(1, nt)

# MODELS
nth = 5         # number of parameters in mortality function
modm = matrix(1, 3, nth)    # matrix to indicate whether the parameters are used
modm[1, 1:3] = 0
modm[2, 1:2] = 0
dimnames(modm) = list(c("GO", "GM", "SI"),
                       c("alpha1", "beta1", "c", "alpha2", "beta2"))
pname = paste(rep(colnames(modm), each = nz),
                "[", rep(colnames(Z), nth), "]", sep = "")
idm = which(rownames(modm) == model)    # choose the mortality function
idth = which(modm[rep(idm, nz), ] == 1)  # which parameters will be used, given mortality function

# FUNCTIONS:
# Survival and mortality:

# Hazard function
m.g = function(x,th) exp(th[, 1] - th[, 2] * x) * modm[idm, 1] * modm[idm, 2] +
      th[, 3] + exp(th[, 4] + th[, 5] * x)
# Survivorship function
S.g = function(x, th){
      Sg = exp((exp(th[, 1]) / th[, 2] *
          (exp(-th[, 2] * x) - 1)) * modm[idm, 1] * modm[idm, 2] +
           x * (-th[, 3]) +
           exp(th[, 4]) / th[, 5] * (1 - exp(th[, 5] * x)))
  return(Sg)
}
# PDF of ages at death
f.g = function(x, th) m.g(x, th) * S.g(x, th)
S.x = function(th) S.g(xv, matrix(th, 1, nth))
m.x = function(th) m.g(xv, matrix(th, 1, nth))

# Lower bounds for parameter c:
c.low = function(th){
if (idm == 1) cl = 0
if (idm == 2) if (th[5] > 0) cl = -exp(th[4]) else if (th[5] < 0) cl = 0
if (idm == 3) {
    x.minf = (th[1] + log(th[2]) - th[4] - log(th[5])) / (th[2] + th[5])
    cl = -exp(th[1] - th[2] * (x.minf)) - exp(th[4] + th[5] * (x.minf))
  }
  return(cl)
}


th.low = matrix(-Inf, nrow(modm), ncol(modm), dimnames = dimnames(modm))
th.low["SI", c("beta1", "beta2")] = 0
low = matrix(th.low[idm, ], nz, nth, byrow = TRUE)
dimnames(low) = list(colnames(Z), colnames(modm))

# Observation matrices:
ObsMatFun = function(f, l) {
    Fm = Tm - f; Fm[Fm >= 0] = 1; Fm[Fm < 0] = 0
    Lm = Tm - l; Lm[Lm <= 0] = -1; Lm[Lm > 0] = 0
    return(Fm * (-Lm))
}










# analysis

# check data
newData <- DataCheck(inputMat, studyStart = 2010,
                      studyEnd = 2017, 
                      silent = FALSE)

# correcting errors with ID 325
inputMat[1043, 2] = inputMat[1043, 2] - 1 
inputMat[1043, 6] = 0 


ni = 100000
nt = 90
nc = 4
nb = 10000


inputMat = as.data.frame(inputMat)
multiout <- multibasta(object = inputMat, studyStart = 2010, studyEnd = 2017, recaptTrans = c(2010:2017),
                       models = c("EX", "GO", "WE", "LO"), shapes = c("simple", "Makeham", "bathtub"),
                       lifeTable = T, minAge = 0.5,
                       nsim = nc, niter = ni, burnin = nb, thinning = nt, 
                       ncpus = 4, parallel = T)

save(multiout, file = "BaSTA_all.RData")

summary(multiout, digits = 3)

plot(multiout$runs$Lo.Ma, fancy = T)



inputMat$Birth = 0
out <- basta(object = inputMat, studyStart = 2010, studyEnd = 2017, covarsStruct = "fused",
             model = "LO", shape = "bathtub", lifeTable = T, minAge = 0.5,
             nsim = nc, niter = ni, burnin = nb, thinning = nt, ncpus = 4, parallel = T)



summary(out, digits = 3)
plot(out)
plot(out, plot.trace = FALSE, fancy = T)










inputMat = inputMat[,1:11]
out <- basta(object = inputMat, studyStart = 2010, studyEnd = 2017, diffrec = diffrec,
             model = "LO", shape = "bathtub", lifeTable = T, minAge = 0.5,
             nsim = nc, niter = ni, burnin = nb, thinning = nt, ncpus = 4, parallel = T)

save(out, file = "BaSTA_top.RData")

load("BaSTA_top.RData")
summary(out)
plot(out, fancy = T)
#--------------------
# mortality functions
#--------------------
# Exponential
par(mfrow = c(3, 4))

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

par(mfrow = c(3, 4))

for (i in 1:length(b0)) {
  plot(a, exp(b0[i] + b1[i] * a))
  plot(a, exp((exp(b0[i]) / b1[i]) * (1 - exp(b1[i] * a))))
}

#-------------------
# Weibull
a = 1:10
b0 = seq(1, 10, length.out = 10)     # 1 < b < 4 reasonable
b1 = seq(0.001, 1, length.out = 10)  # 0.2 reasonable

par(mfrow = c(3, 4))

for (i in 1:length(b0)) {
  plot(a, ((b0[i] * b1[i]) * (b1[i] * a) ^ (b0[i] - 1)))
  plot(a, exp(-(b1[i] * a) ^ b0[i]))
}

#------------------
# Logistic
a = 1:10
b0 = seq(-10, -.1, length.out = 10)
b1 = seq(0.1, 5, length.out = 10)
b2 = seq(0.1, 5, length.out = 10)

par(mfrow = c(3, 4))

for (i in 1:length(b0)) {
  plot(a, (exp(b0[i] + (b1[i] * a)) / ((1 + (b2[i] * (exp(b0[i]) / b1[i])) * (exp(b1[i] * a) - 1)))))
  plot(a, ((1 + b2[i] * exp(b0[i]) / b1[i]) * (exp(b1[i] * a) - 1)) ^ -(1 / b2[i]))
}








