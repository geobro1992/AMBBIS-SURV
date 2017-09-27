# newts.R

# Example with a user-defined function.

# Commands for analysis of newts data, 17 newts, 22 samples (2002).

source("hetage.functions.R")
dyn.load("./hetageLL.dll")

newts.df <- read.table("newts2002.txt")

# There is no column labelled "freq", so it is assumed
# that all frequencies are 1.
# Data will be condensed if there are repeated rows.

hetage.process.data(newts.df)

# Now need to write user-defined functions.
# Build a model with beta and phi both on logistic curves over time.
  
# (1) Need a v.to.f function to expand from vector of minimum
# number of independent parameters to a full parameter list.

# Set up beta.l.function to create beta vector from coefficients:

beta.l.fn <- function (coefficients,K)
  {
  ab <- coefficients[1]
  bb <- coefficients[2]
  beta.logit <- ab*(1:K) + bb
  beta.out <- expit(beta.logit)
  beta.out/sum(beta.out)
  }

# Set up phi.l.function to set up a phi array from coefficients:

phi.l.fn <- function (coefficients,K){
  aphi <- coefficients[1]
  bphi <- coefficients[2]
  phi.logit <- aphi*(1:K) + bphi
  nlrv(expit(phi.logit),K)
  }

# Define the v.to.f function:

v.to.f.user1 <- function (invect)
  {
  list (N     = exp(invect[1]),
        beta  = beta.l.fn(invect[2:3],K),
        phi   = phi.l.fn(invect[4:5],K-1),
        p     = nlc(expit(invect[(6)]),K),
        pie   = c(1),
        coefs = (invect[2:5])
       )
  }


# Need a start function to reduce a full set of parameters to a
# vector of independent parameters. May need to invent coefficients.

start.user1 <- function (parameters)
  {
  if (is.null(parameters$coefs)) parameters$coefs<-c(-1,1.5,-0.1,1.8)
  # Override since the coefficients are null
  c(
    log(parameters$N),           # [1]
    parameters$coefs[1:2],       # [2:3]    a.beta, b.beta
    parameters$coefs[3:4],       # [4:5]    a.phi, b.phi
    logit(parameters$p[1,1,1])   # [6]
   )
  }

# Try using phi(t)p(c) to get a good start.

phit.pc.out <- hetage.fit.model("phit.pc")
beta.vect   <- phit.pc.out$parameters$beta
plot(1:22,beta.vect) # Goes through 0.5 at about sample 1.5.
                     # Hence bbeta approx -1.5*abeta
phi.vect    <- phit.pc.out$parameters$phi[1,,1]
plot(1:21,phi.vect)  # Goes through 0.5 at about sample 18.
                     # Hence bphi approx -18*aphi

inpars <- phit.pc.out$parameters
inpars$coefs <- c(-1,1.5,-1,18)

inparvect <- start.user1(inpars)

v.to.f.user1(inparvect)

mLL(inparvect,expandfn=v.to.f.user1)



betal.phil.pc.out <- hetage.fit.model("user1",start=inpars)

betal.phil.pc.out$parameters
