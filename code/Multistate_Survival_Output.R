

load("SURV_OUT2.RData")

beta.post = js.ms$sims.list$beta
phi.post = js.ms$sims.list$mean.phi
mu.post = js.ms$sims.list$mu

mcmc = as.data.frame(cbind(beta.post, phi.post))

x = mean(mu.post) + 30:70*mean(beta.post)
x.u = quantile(mu.post, .75) + 30:70*quantile(beta.post, .75)
x.l = quantile(mu.post, .25) + 30:70*quantile(beta.post, .25)

y = exp(x)/(1+exp(x))
y.u = exp(x.u)/(1+exp(x.u))
y.l = exp(x.l)/(1+exp(x.l))

plot(30:70, y, ylim = c(0,1), type = "l", lwd = 2)
lines(30:70, y.u, lty = "dashed")
lines(30:70, y.l, lty = "dashed")

#------------------
pp = js.ms$sims.list[c(1:5)]
dat.melt <- melt(pp)
levels(dat.melt$variable) = c("beta")
ggplot()+
  geom_density(data = dat.melt, aes(value, fill = L1, colour = L1), 
               alpha=I(.1), size =1)  +
  facet_wrap(~ L1, scales = "free", nrow = 3)


pg = js.ms$sims.list[8]
dat.melt <- melt(pg$B)

ggplot()+
  geom_density(data = dat.melt, aes(value, color = "red", fill = "red"), 
               alpha=I(.1), size =1)


library(bayesplot)
library(ggplot2)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(mcmc,
           prob = 0.95)


tiff("Fig1.tiff", width = 12, height = 8, units = 'in', res = 300)
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
h <- hist(beta.post, freq = FALSE, main = "", xlab = "", ylab = " ", xlim = c(-0.20, .20), axes = FALSE, col = "grey", breaks = 50)
axis(1, seq(-.2, .2, by = .1))
mtext("Posterior Estimate of Beta", side = 1, line = 2.5, cex = 1.5, font = 2)
dev.off()

##################
# predicted curves
##################
mu.phi = mean(dat$`1.mu`) + mean(dat$`1.beta`)*seq(3.4, 4.4, length.out = 20)

u.phi = quantile(dat$`1.mu`, .75) + quantile(dat$`1.beta`, .5)*seq(3.4, 4.4, length.out = 20)
l.phi = quantile(dat$`1.mu`, .25) + quantile(dat$`1.beta`, .5)*seq(3.4, 4.4, length.out = 20)
phi.mu = exp(mu.phi)/(1 + exp(mu.phi))
phi.l = exp(l.phi)/(1 + exp(l.phi))
phi.u = exp(u.phi)/(1 + exp(u.phi))
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(exp(seq(3.4, 4.4, length.out = 20)), phi.mu, type = "l", lwd = 2,  ylim = c(0,1), axes = F, xlab = "", ylab = "")
lines(exp(seq(3.4, 4.4, length.out = 20)), phi.u, lty = "dashed")
lines(exp(seq(3.4, 4.4, length.out = 20)), phi.l, lty = "dashed")
axis(1, at = c(30,40,50,60,70,80))
mtext("Size / mm", side = 1, line = 3, cex = 1.5, font = 2)
axis(2)
par(las = 0)
mtext(expression(paste("Survival Probability ", phi)), side = 2, line = 3, cex = 1.5, font = 2)
