################################################
# CJS model for flatwoods salamanders
# with a continuous predictor and random effects

# required libraries
library(tidyverse)
library(dplyr)
library(lubridate)
library(R2WinBUGS)

##############
# read in data
dat = read.csv("fw_dat.csv")

dat = dat %>%
  mutate(season = year(as.Date(dat$Date, format = "%m/%d/%Y") - (185))) %>%
  mutate(mon = month(as.Date(dat$Date, format = "%m/%d/%Y"))) %>%
  mutate(Date = as.Date(dat$Date, format = "%m/%d/%Y"))

# continuous covariate matrix
h = dat %>%
  select(MasterID, season, SVL_1stMeasurement) %>%
  group_by(MasterID) %>%
  mutate(num.capture = 1:length(MasterID)) %>%
  ungroup() %>%
  group_by(MasterID, season) %>%
  summarise(mean.svl=mean(SVL_1stMeasurement, na.rm=TRUE)) %>%
  ungroup() %>%
  spread(season, mean.svl) %>%
  select(-MasterID, -`2009`)

h.ms = h[rowSums(h, na.rm = T)>0,]

y = h.ms
# Recode CH matrix
y[!is.na(y)] <- 1
y[is.na(y)] <- 0

# Recode size matrix
h.ms[is.na(h.ms)] <- 0

#function to get occasion of marking for each animal
get.first<-function(x) min(which(x!=0))

f<-apply(h.ms,1,get.first)

# impute body sizes on occasions when animal was not observed
for (i in 1:dim(h.ms)[1]) {
  if (h.ms[i,1] == 0) {
    h.ms[i,1] <-  runif(n = 1,min = 30,max = as.numeric(h.ms[i,f[i]]) )
  } else{}
}

for (i in 1:dim(h.ms)[1]) {
  Linf = rnorm(1,70,5)
  Linf = ifelse(Linf < max(h.ms[i,], na.rm = T), max(h.ms[i,], na.rm = T), Linf)
  k = runif(1, 0.5, 1) 
  for (j in 2:11) {
    if (h.ms[i,j] == 0) {
      h.ms[i,j] <-  h.ms[i,j-1] + (Linf - h.ms[i,j-1])*(1-exp(-k))
    } else{}
  }  
}


# remove single captures
sings = which(rowSums(y, na.rm = T) == 1)
y = y[-sings,]
h.ms = h.ms[-sings,]
f<-apply(y,1,get.first)

#add to data known states (where we know z=1)
known.states.cjs<-function(ch){
  state<-ch
  for (i in 1:dim(ch)[1]){
    n1<-min(which(ch[i,]==1))
    n2<-max(which(ch[i,]==1))
    state[i,n1:n2]<-1
    state[i,n1]<-NA
  }
  state[state==0]<-NA
  return(state)
}

## initialize z states have to be NA where observed
cjs.init.z<-function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2<-max(which(ch[i,]==1))
    ch[i,f[i]:n2]<-NA
  }
  for (i in 1:dim(ch)[1])
  { ch[i,1:f[i]]<-NA
  }
  return(ch)
}


######################
# mixed model in BUGS

# Bundle data
bugs.data <- list(y = as.matrix(y), x = as.matrix(h.ms),
                  z = as.matrix(known.states.cjs(y)), f = f,
                  n.occasions = dim(y)[2], 
                  M = dim(y)[1],
                  mean.svl = mean(as.matrix(h.ms)), sd.svl = sd(as.matrix(h.ms)))


# MCMC settings
ni <- 10000
nb <- 1000
nt <- 9
nc <- 1

inits <- function(){list(mean.phi = runif(1, 0, 1), sigma.phi = runif(1, 0, 1),  
                         mean.p = runif(1, 0, 1), sigma.p = runif(1, 0, 1),
                         b1.mu = runif(10, 0.1, .5),
                         z=as.matrix(cjs.init.z(y,f)))}

# Parameters monitored
params <- c("mean.phi",
            "mu.phi",
            "sigma.phi",
            "epsilon.phi",
            "mu.p",
            "sigma.p",
            "epsilon.p",
            "b1.mu")

# Call JAGS from R

js <- bugs(data = bugs.data, inits = inits, parameters.to.save = params, 
              model.file = "surv_23_RE.txt", n.chains = nc, debug=TRUE,
              n.thin = nt, n.iter = ni, n.burnin = nb, bugs.directory = "C:/Program Files (x86)/WinBUGS14", working.directory = getwd())

###################
# MODEL OUTPUT
###################
beta.post = js$sims.list$b1.mu
phi.post = js$sims.list$mean.phi
mu.post = js$sims.list$mu.phi
sigma.post = js$sims.list$sigma.phi
epsilon.post = js$sims.list$epsilon.phi
mu.p.post = js$sims.list$mu.p
epsilon.p = js$sims.list$epsilon.p


#################################
# size-dependent survival by year

sizes = (30:80-bugs.data$mean.svl)/ bugs.data$sd.svl 

x = mean(mu.post) + sizes*mean(beta.post[,1]) + mean(epsilon.post[,1])

xsim = matrix(nrow = 1000, ncol = length(sizes))

for(j in 1:1000){
  xsim[j,] = mu.post[j] + sizes*beta.post[j,1]  + epsilon.post[j,1]
}

x.u = apply(xsim , 2 , quantile , probs = c(0.9) , na.rm = TRUE)
x.l = apply(xsim , 2 , quantile , probs = c(0.1) , na.rm = TRUE)

y = exp(x)/(1+exp(x))
y.u = exp(x.u)/(1+exp(x.u))
y.l = exp(x.l)/(1+exp(x.l))

pred.dat = data.frame(svl = 30:80, y, y.u, y.l)

# loop through subsequent years
for(i in 2:10){
  
  x = mean(mu.post) + sizes*mean(beta.post[,i]) + mean(epsilon.post[,i])
  
  xsim = matrix(nrow = 1000, ncol = length(sizes))
  
  for(j in 1:1000){
    xsim[j,] = mu.post[j] + sizes*beta.post[j,i] + epsilon.post[j,i]
  }
  
  x.u = apply(xsim , 2 , quantile , probs = c(0.9) , na.rm = TRUE)
  x.l = apply(xsim , 2 , quantile , probs = c(0.1) , na.rm = TRUE)
  
  y = exp(x)/(1+exp(x))
  y.u = exp(x.u)/(1+exp(x.u))
  y.l = exp(x.l)/(1+exp(x.l))
  
  pred.dat = rbind(pred.dat, cbind(svl = 30:80, y, y.u, y.l))
  
}

pred.dat$season = factor(rep(1:10, each = 51))
levels(pred.dat$season) = c("2011", "2012", "2013", "2014", "2015",
                            "2016", "2017", "2018", "2019", "2020")

p1 = ggplot(filter(pred.dat, season != "2011", season != "2020"), aes(x = svl, y = y)) +
  geom_line() + ylim(0,1) +
  geom_ribbon(aes(ymin = y.l, ymax = y.u), alpha = 0.3) +
  theme_Publication()+
  facet_wrap(~season, nrow = 2) +
  xlab("Body Size (mm)") + ylab("Survival Probability")


################################
# survival averaged across years

sizes = (30:80-bugs.data$mean.svl)/ bugs.data$sd.svl 

x = mean(mu.post) + sizes*mean(beta.post[,2:9])

xsim = matrix(nrow = 1000, ncol = length(sizes))

for(j in 1:1000){
  xsim[j,] = mu.post[j] + sizes*mean(beta.post[,2:9])
}

x.u = apply(xsim , 2 , quantile , probs = c(0.9) , na.rm = TRUE)
x.l = apply(xsim , 2 , quantile , probs = c(0.1) , na.rm = TRUE)

y = exp(x)/(1+exp(x))
y.u = exp(x.u)/(1+exp(x.u))
y.l = exp(x.l)/(1+exp(x.l))

pred.dat = data.frame(svl = 30:80, y, y.u, y.l)


p1 = ggplot(pred.dat, aes(x = svl, y = y)) +
  geom_line() + ylim(0,1) +
  geom_ribbon(aes(ymin = y.l, ymax = y.u), alpha = 0.3) +
  theme_Publication() +
  xlab("Body Size (mm)") + ylab("Survival Probability")


ggsave(filename = "surv_re_average_fig.pdf", plot = p1, device = cairo_pdf, width = 6, height = 4)


################
# posterior plots
library(ggdist)

##########
# beta
p.post = gather(as.data.frame(js.re$sims.list$b1.mu))
p.post$key = factor(p.post$key, levels = unique(p.post$key))

levels(p.post$key) = c("2011", "2012", "2013", "2014", "2015",
                       "2016", "2017", "2018", "2019", "2020")

p1 = ggplot(p.post, aes(x = key, y = value)) +
  stat_gradientinterval() +
  theme_Publication() +
  xlab("Year") + ylab("Size:Survival Slope") + ylim(0,2)

ggsave(filename = "surv_beta_posteriors.pdf", plot = p1, device = cairo_pdf, width = 8, height = 6)

##########
# phi

p.post = gather(as.data.frame(js.re$sims.list$epsilon.phi))
p.post$key = factor(p.post$key, levels = unique(p.post$key))
p.post$mu = rep(mu.post,10)

p.post$true.phi = exp(p.post$mu + p.post$value)/(1+exp(p.post$mu + p.post$value))

levels(p.post$key) = c("2011", "2012", "2013", "2014", "2015",
                       "2016", "2017", "2018", "2019", "2020")

p2 = ggplot(p.post, aes(x = key, y = true.phi)) +
  stat_eye(width = 2) +
  theme_Publication() +
  xlab("Year") + ylab("Survival Probability")

ggsave(filename = "surv_phi_posteriors.pdf", plot = p2, device = cairo_pdf, width = 8, height = 6)


##########
# p

p.post = gather(as.data.frame(js.re$sims.list$epsilon.p))
p.post$key = factor(p.post$key, levels = unique(p.post$key))
p.post$mu = rep(mu.p.post,10)

x = mean(mu.post) + sizes*mean(beta.post[,2:9])

p.post$true.p = exp(p.post$mu + p.post$value)/(1+exp(p.post$mu + p.post$value))

levels(p.post$key) = c("2011", "2012", "2013", "2014", "2015",
                       "2016", "2017", "2018", "2019", "2020")

p2 = ggplot(p.post, aes(x = key, y = true.p)) +
  stat_gradientinterval() +
  theme_Publication() +
  xlab("Year") + ylab("Detection Probability")

ggsave(filename = "surv_p_posteriors.pdf", plot = p2, device = cairo_pdf, width = 8, height = 6)



############################################
# correlations with environmental conditions
############################################

tmax = read.csv("summer_temp.csv")
tmax = filter(tmax, Group.1 > "2011", Group.1 < "2020")
names(tmax) = c("ID", "Year", "tmax")

tmin = read.csv("winter_temp.csv")
tmin = filter(tmin, Group.1 > "2011", Group.1 < "2020")
names(tmin) = c("ID", "Year", "tmin")

ppt = read.csv("summer_ppt.csv")
ppt = filter(ppt, Group.1 > "2011", Group.1 < "2020")
names(ppt) = c("ID", "Year", "ppt")

wppt = read.csv("winter_ppt.csv")
wppt = filter(wppt, Group.1 > "2011", Group.1 < "2020")
names(wppt) = c("ID", "Year", "wppt")


# mean survival each year
x = mean(mu.post) + colMeans(epsilon.post[,2:9])
mean.phi.y = exp(x)/(1+exp(x))

# mean survival:size slope
mean.beta.y = colMeans(beta.post[,2:9])

# mean detection each year
x.p = mean(mu.p.post) + colMeans(epsilon.p[,2:9])
mean.p.y = exp(x.p)/(1+exp(x.p))

# combine all environmental data and vital rate estimates
corr.dat = data.frame(tmax, tmin = tmin$tmin, ppt = ppt$ppt, wppt = wppt$wppt, phi = mean.phi.y, beta = mean.beta.y, p = mean.p.y)

###########
# correlation coefficients

cor.test(corr.dat$tmax, corr.dat$phi, method = "kendall", alternative = "less")      # 0.00, p = 0.55
cor.test(corr.dat$tmax, corr.dat$p, method = "kendall", alternative = "greater")     # 0.36, p = 0.14

cor.test(corr.dat$ppt, corr.dat$phi, method = "kendall", alternative = "greater")    # -0.07, p = 0.64
cor.test(corr.dat$ppt, corr.dat$p, method = "kendall", alternative = "less")         # -0.57, p = 0.03

cor.test(corr.dat$tmin, corr.dat$phi, method = "kendall", alternative = "greater")   # -0.14, p = 0.73
cor.test(corr.dat$tmin, corr.dat$p, method = "kendall", alternative = "less")        # -0.36, p = 0.14

cor.test(corr.dat$wppt, corr.dat$phi, method = "kendall", alternative = "greater")   # -0.14, p = 0.72
cor.test(corr.dat$wppt, corr.dat$p, method = "kendall", alternative = "greater")     # -0.07, p = 0.64


############################################################
# plot only significant correlation (detection ~ summer ppt)
library(plotrix) # package plotrix is needed for function "ablineclip""

# if the following line and the line containing "dev.off()" are executed, the plot will be saved as a png file in the current working directory
pdf("surv_corr_plot.pdf", width = 6, height = 6)

op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(corr.dat$ppt, corr.dat$p, col = "black", pch = 21, bg = "grey", cex = 2,
     ylim = c(0, 1), ylab = "", xlab = "", axes = FALSE)
axis(1)
axis(2) 
reg1 <- lm(p ~ ppt, data = corr.dat)
ablineclip(reg1, lwd = 2, x1 = 320, x2 = 900) 
par(las = 0)
mtext("Total Summer Precipitation (mm)", side = 1, line = 2.5, cex = 1.5)
mtext("Detection Probability", side = 2, line = 3.7, cex = 1.5)
text(800, .9, expression(paste(tau, " = -0.57"), sep = ""), cex = 1.5)

dev.off()
