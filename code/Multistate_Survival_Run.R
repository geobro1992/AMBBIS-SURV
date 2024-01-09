# Analysis of the JS model as a multistate model
rm(list = ls())
library(R2WinBUGS)
library(ggplot2)
library(plyr)
library(reshape)
library(RODBC)
library(tidyverse)
library(lubridate)
library(forcats)

#------------------------------------------
# Connect to Minnow server
channel <- odbcConnectAccess2007(access.file = "//minnow.cc.vt.edu/cnre2/Eglin/Projects/FlatwoodsSalamander/DriftFence/DriftFenceDatabase/Driftfence Database_GB_29Nov17_BACK_END.accdb")
# Read in capture table
df <- sqlFetch(channel,"Flatwoods Capture") %>%
  tbl_df() %>%
  filter(Age == "M") %>%
  filter(`REC/CAP` == "CAP") %>%
  mutate(season = year(Date - (24*60*60*180))) %>%
  glimpse()
#df = read.csv("FW_main.csv", sep = ",")
df = df[, c("MasterID", "Age", "Date", "Sex", "Weight", "TL", "SVL_1stMeasurement")]
df = as.data.frame(df)
dl = split(df, df$MasterID, drop = T)
names(dl) = NULL


# Individual Capture Histories and data augmentation
y = df %>%
  select(MasterID, Date, SVL_1stMeasurement) %>%
  mutate(SVL_1stMeasurement = 1) %>%
  group_by(MasterID) %>%
  mutate(num.capture = 1:length(MasterID)) %>%
  mutate(Year = year(Date)) %>%
  ungroup() %>%
  group_by(MasterID, Year) %>%
  summarise(mean.svl=mean(SVL_1stMeasurement, na.rm=TRUE)) %>%
  ungroup() %>%
  spread(Year, mean.svl) %>%
  select(-MasterID)

y = as.matrix(y)
sings = which(rowSums(y, na.rm = T) == 1)
#y = y[-sings,]
nz = 200
y.ms <- rbind(y, matrix(0, ncol = dim(y)[2], nrow = nz))

# Recode CH matrix: a 0 is not allowed in WinBUGS!
y.ms[is.na(y.ms)] <- 0
y.ms[y.ms==0] <- 2


# continuous covariate matrix
h = df %>%
  select(MasterID, Date, SVL_1stMeasurement) %>%
  group_by(MasterID) %>%
  mutate(num.capture = 1:length(MasterID)) %>%
  mutate(Year = year(Date)) %>%
  ungroup() %>%
  group_by(MasterID, Year) %>%
  summarise(mean.svl=mean(SVL_1stMeasurement, na.rm=TRUE)) %>%
  ungroup() %>%
  spread(Year, mean.svl) %>%
  select(-MasterID)

h = as.matrix(h)
#h = h[-sings,]
h.ms <- rbind(h, matrix(0, ncol = dim(h)[2], nrow = nz))

# Recode CH matrix: a 0 is not allowed in WinBUGS!
h.ms[is.na(h.ms)] <- 0

for (i in 1:dim(h.ms)[1]) {
  if (h.ms[i,1] == 0) {
    h.ms[i,1] <-  runif(1,35,45)
  } else{}
}


for (i in 1:dim(h.ms)[1]) {
  Linf = rnorm(1,59,5)
  for (j in 2:10) {
    if (h.ms[i,j] == 0) {
      h.ms[i,j] <-  h.ms[i,j-1] + (Linf - h.ms[i,j-1])*(1-exp(-0.91))
    } else{}
  }  
}

# Bundle data
bugs.data <- list(y = y.ms, x = log(h.ms),
                  n.occasions = dim(y.ms)[2], 
                  M = dim(y.ms)[1])
n.occasions = 10
# Initial values
#--------------------- Initial values
inits <- function(){list(mean.phi = runif(1, 0.5, 1), 
                         p = runif(1, 0.7, 1), 
                         beta = runif(1, -0.5, .5),
                         sigma = runif(1, 0, .1),
                         z = cbind(rep(1, times = nrow(y.ms)), y.ms[,-1]))}
# Parameters monitored
params <- c("mean.phi",
            "mu",
            "p",
            "beta",
            "sigma2",
            "nu",
            "b", "Nsuper")

# MCMC settings
ni <- 120000
nb <- 20000
nt <- 100
nc <- 1

# Call JAGS from R

js.ms <- bugs(data = bugs.data, inits = bug.inits, parameters.to.save = variable.names, 
              model.file = "model_JS2.txt", n.chains = nchains, debug=TRUE,
              n.thin = nt, n.iter = ni, n.burnin = nb, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

save(js.ms, file = "SURV_OUT_2021.RData")

dat = as.data.frame(js.ms$sims.array)
dat2 = melt(dat[,c(34,37,40,43,46,49,52,55,58)])

###################
# Posteriors
####################

# save posteriors for IPM projections
beta.post = js.ms$sims.list$beta
save(beta.post, file = "SURV_beta.RData")
phi.post = js.ms$sims.list$mean.phi
save(phi.post, file = "SURV_phi.RData")
mu.post = js.ms$sims.list$mu
save(mu.post, file = "SURV_mu.RData")


odbcCloseAll()


