#########################
# Survival Kernel for IPM
#########################
rm(list = ls())
library(R2WinBUGS)
library(ggplot2)
library(plyr)
library(reshape)
library(RODBC)
library(tidyverse)
library(lubridate)
library(forcats)
library(ggalt)
library(lme4)
library(gam)
library(voxel)
library(mgcv)
library(gridExtra)
library(ggthemes)

# Connect to Minnow server
channel <- odbcConnectAccess2007("//minnow.cc.vt.edu/cnre2/Eglin/Projects/FlatwoodsSalamander/DriftFence/DriftFenceDatabase/Driftfence Database_GB_29Nov17_BACK_END.accdb")

# Read in capture table
df <- sqlFetch(channel,"Flatwoods Capture") %>%
  tbl_df() %>%
  filter(MasterID != 0) %>%
  mutate(Year = year(Date)) %>%
  mutate(Season = year(Date - 15768000)) %>%
  #  filter(Season > 2009) %>%
  mutate(Month = month(Date)) %>%
  arrange(Date) %>%
 # filter(Age != "M") %>%
  select("MasterID", "Date", "SVL_1stMeasurement", "Season", "Year", "Month") %>%
  rename("z" = "SVL_1stMeasurement", "ID" = "MasterID") %>%
  glimpse()

df = as.data.frame(df)
df = df[order(df[,1]),]
df = df[with(df, order(Year)),]

#df <- df %>%
#  mutate(T2 = Month/12) %>%
#  select("ID", "z", "Season", "T2")

#df$T2 = round(df$T2, 0)
#df$Season = as.factor(df$Season)
#df$T2 = as.factor(df$T2)
#df$period <- with(df, interaction(Season, T2))

#zuse = df %>%
#  select("ID", "period", "z")


#library(reshape2)
#h = dcast(zuse, ID ~ period, value.var = "z", fun.aggregate = mean)


# code whether last seen leaving or arriving (using month column)
CO = vector()
for(i in 1:length(df[,1])) {
  if(df[i,6] < 7 || df[i,6] == 12) {
    CO[i] = 2    
  } else CO[i] = 1
}

df = cbind(df, CO)

# row for every year, season, individual combination
df = df %>% 
  complete(ID, nesting(Season, CO)) %>%
  select("CO", "Season", "ID", "z") 
df = cbind(df, CH = rep(NA, length(df$z)))

df = as.data.frame(df)
df$ID = as.numeric(as.factor(df$ID))


dl = split(df, df$ID, drop = T)
names(dl) = NULL
n.ind = length(dl)
n.caps = sapply(dl, NROW)
table(n.caps)
# capture histories; 1 for everyhting before last capture, last capture is 0
for(i in 1:n.ind){
  if(n.caps[i] == 21){
    dl[[i]]$CH[which(dl[[i]]$z != "NA")] = 0
  } else {
    dl[[i]]$CH[max(which(dl[[i]]$z != "NA"))] = 0
    dl[[i]]$CH[min(which(dl[[i]]$z != "NA")):(max(which(dl[[i]]$z != "NA"))-1)] = 1
  }
}

surv.df = do.call("rbind", dl)
surv.df = surv.df[complete.cases(surv.df),]

surv.df$CO = as.factor(surv.df$CO)
levels(surv.df$CO) = c("Breeding", "Non-breeding")

surv.df$Season = as.factor(surv.df$Season)
#####################
# logistic regression

# BREEDING SEASON
surv.dfB = surv.df[which(surv.df$CO == "l.s.arriving"),]
modSB = glm(CH ~ z, data = surv.dfB, family = binomial)
summary(modSB)

# SUMMER SEASON
surv.dfD = surv.df[which(surv.df$CO == "l.s.leaving"),]
modSD = glm(CH ~ z, data = surv.dfD, family = binomial)
summary(modSD)

# ANNUAL MODEL
modSA = glm(CH ~ z + CO, data = surv.df, family = binomial)
summary(modSA)
mod2 = glm(CH ~ z, data = surv.df, family = binomial)
summary(mod2)
mod3 = glm(CH ~ CO, data = surv.df, family = binomial)
summary(mod3)
mod4 = glm(CH ~ 1, data = surv.df, family = binomial)
summary(mod4)
mod5 = glm(CH ~ z + CO + Season, data = surv.df, family = binomial)
summary(mod5)
mod6 = glm(CH ~ z + Season, data = surv.df, family = binomial)
summary(mod6)
mod7 = glm(CH ~ Season, data = surv.df, family = binomial)
summary(mod7)
mod8 = glm(CH ~ CO + Season, data = surv.df, family = binomial)
summary(mod8)


aics = AIC(modSA, mod2, mod3, mod4, mod5, mod6, mod7, mod8)$AIC
  
daics = aics - aics[5]


round(exp(-0.5*daics)/sum(exp(-0.5*daics)),3) 


#######################################
# plot predicted survival probabilities

ggplot(surv.df, aes(z, CH, color = CO)) +
  facet_wrap(~Season) +
  stat_smooth(method="glm", formula=y~x,
              alpha=0.2, size=2) 
xz = rep(30:80, 9)
ys = rep(levels(surv.df$Season), length.out = 459)
season1 = rep("Breeding", length(xz)) 
season2 = rep("Non-breeding", length(xz))
preds <- predict(mod5, list(z = xz, CO = season1, Season = ys),type="response", se.fit = T)
preds2 <- predict(mod5, list(z = xz, CO = season2, Season = ys),type="response", se.fit = T)
lower <- preds$fit - (1.96*preds$se.fit) # lower bounds
upper <- preds$fit + (1.96*preds$se.fit) # upper bounds
lower2 <- preds2$fit - (1.96*preds2$se.fit) # lower bounds
upper2 <- preds2$fit + (1.96*preds2$se.fit) # upper bounds

df = as.data.frame(cbind(preds$fit, ys, xz, season1))
df$V1 = as.numeric(as.character(df$V1))
df$xz = as.numeric(as.character(df$xz))

df2 = as.data.frame(cbind(preds2$fit, ys, xz, season2))
df2$V1 = as.numeric(as.character(df2$V1))
df2$xz = as.numeric(as.character(df2$xz))

names(df2) = names(df)


df = rbind(df, df2)

tiff(filename = "fig1.tiff", res = 600, width = 8, height = 6, units = "in")
ggplot(df, aes(xz, V1, color = season1)) +
  stat_smooth(method="glm", formula=y~x,
              alpha=0.2, size=2) + facet_wrap(~ys)+
  xlab("Body Size") + ylab("Survival Probability") + labs(color = "Season")+
 theme_Publication()
dev.off()


y1 = aggregate(preds$fit, list(xz), mean)$x
y2 = aggregate(preds2$fit, list(xz), mean)$x

se1 = aggregate(preds$se.fit, list(xz), mean)$x
se2 = aggregate(preds2$se.fit, list(xz), mean)$x

preds <- predict(mod5, list(z = xz, CO = season1, Season = ys),interval = "prediction", type="response", se.fit = T)


par(cex.main = 1.5, mar = c(5, 8, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(xz, rep(-10, length(xz)), type = "p", ylab = "", xlab = " ", cex = 1.5, 
     ylim = c(0, 0.9), xlim = c(40, 80), lwd = 2, pch = 5, axes = F, main = " ")
axis(1, cex.axis=2, at = c(30, 40, 50, 60, 70, 80), labels = c("30", "40", "50", "60", "70", "80"))
mtext("Size / mm", side = 1, line = 3, cex = 2, font = 2)
axis(2,  cex.axis=2)
par(las = 0)
mtext(expression(paste("Survival ", phi)), side = 2, line = 4, cex = 2, font = 2)

lines(30:80, y1, cex = 1.5, lwd = 2)
lines(30:80, y2, cex = 1.5, lwd = 2, pch = "19")

lines(30:80, y1 + se1, lty="dotted", lwd = 2)
lines(30:80, y1 - se1, lty="dotted", lwd = 2)


lines(30:80, y2 + se2, lty="dotted", lwd = 2)



text(60, 0.2, "Non-breeding", cex = 2, font = 1, adj = 0)
text(40, 0.8, "Breeding", cex = 2, font = 1, adj = 0)


###########################
tiff(filename = "fig2.tiff", res = 600, width = 8, height = 6, units = "in")
ggplot(df, aes(xz, V1, color = season1)) +
  stat_smooth(method="glm", formula=y~x,
              alpha=0.2, level = .95)+
  xlab("Body Size") + ylab("Survival Probability") + labs(color = "Season")+
  theme_Publication()
dev.off()


##################################
# gam for nonlinear effect of size
modS2 = mgcv::gam(CH ~ s(z), data = surv.df, family = binomial)
summary(modS2)

fit <- predict(modS2, se = TRUE, type = "response")$fit
se <- predict(modS2, se = TRUE, type = "response")$se.fit
lcl <- fit - 1.96 * se
ucl <- fit + 1.96 * se

pred = as.data.frame(cbind(fit, surv.df$z, se, lcl, ucl))

ggplot(pred, aes(V2, fit)) +
  geom_line(lwd = 1.5, linetype = 2) +
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha = 0.3)

#########################
# define survival kernels
s_z <- function(z, m.par) {
  
  linear.p <- sample(m.par[["mu"]], 1) + sample(m.par[["beta"]], 1) * z
  # linear predictor
  p <- 1/(1+exp(-linear.p))
  # inv-logistic trans
  return(p)
}


#####################

modS = glm(CH ~ z, data = surv.df, family = binomial)
summary(modS)
plot(modS)

# parameter estimates
mtype = "" # choose "B", "D", or "A" based on seasonal vs annual models
m.par = coef(modS)
#m.par <- c(surv = coef(noquote(paste("modS",mtype, sep=""))))
1/(1+exp(.5))
m.par = cbind(beta = seq(from = 0.05 ,to =  0.1, length.out = 1000), mu = seq(from = -.5, to = -.1, length.out = 1000)) 
m.par = c(order(vGBE$summary.list$mu), order(vGBE$summary.list$beta))
m.par = m.par[-c(1:50,951:1000),]
names(m.par) <- c("mu", "beta") 



#####################################
## Combine the survival-growth kernels
P <- function (z1, z, m.par) {
  return(s_z(z, m.par) * G_z1z(z1, z, m.par))
}





df <- sqlFetch(channel,"Flatwoods Capture") %>%
  tbl_df() %>%
  mutate(Year = year(Date)) %>% 
  mutate(Month = month(Date)) %>% 
  filter(Year > 2017 & Month > 6) 

table(df$`Direction of Travel`)  






















theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}






xz = 40:80
season1 = rep("Breeding", length(xz)) 
season2 = rep("Non-breeding", length(xz))

ys = rep(levels(surv.df$Season)[8], length.out = length(xz))
preds <- predict(mod5, list(z = xz, CO = season1, Season = ys),type="response", se.fit = T)
preds2 <- predict(mod5, list(z = xz, CO = season2, Season = ys),type="response", se.fit = T)
lower <- preds$fit - (1.96*preds$se.fit) # lower bounds
lower2 <- preds2$fit - (1.96*preds2$se.fit) # lower bounds


ys = rep(levels(surv.df$Season)[3], length.out = length(xz))
preds <- predict(mod5, list(z = xz, CO = season1, Season = ys),type="response", se.fit = T)
preds2 <- predict(mod5, list(z = xz, CO = season2, Season = ys),type="response", se.fit = T)
upper <- preds$fit + (1.96*preds$se.fit) # lower bounds
upper2 <- preds2$fit + (1.96*preds2$se.fit) # lower bounds

d1 = data.frame(x = xz, fit = (lower+upper)/2, upr = upper, lwr = lower, season = rep(1, length(xz)))
d2 = data.frame(x = xz, fit = (lower2+upper2)/2, upr = upper2, lwr = lower2, season = rep(2, length(xz)))

db = rbind(d1,d2)
db$season = as.factor(db$season)


tiff(filename = "fig2.tiff", res = 600, width = 8, height = 6, units = "in")
ggplot(db, aes(x, fit, color = season)) +
  geom_line(lwd = 2) +
  geom_ribbon(aes(ymin = upr, ymax = lwr, fill = season), alpha = 0.2)+
  xlab("Body Size") + ylab("Survival Probability") +
  theme_Publication()
dev.off()
