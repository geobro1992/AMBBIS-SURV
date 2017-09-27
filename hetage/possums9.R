# possums9.R

# Commands for analysis of possums9 data, 270 animals, 9 years.

# Load the functions and the data set:

source("hetage.functions.R")
dyn.load("./hetageLL.dll")

possums9.df <- read.table("possums9.txt")

# Check out the appearance of the data set:

str(possums9.df)

# Remove first column, this was just an index:

possums9.df <- possums9.df[,-1]

dim(possums9.df)          # 270 possums by 9 sampling occasions
apply(possums9.df,2,sum)  # Any leading or trailing zeros?
                      # If so, remove those columns.

# There is no column labelled "freq", so it is assumed
# that all frequencies are 1. In this case, data will be 
# condensed if there are repeated rows.

hetage.process.data(possums9.df)

# View the working data:

x.mat         # Only 91 distinct capture histories
y.vect        # Frequencies of the capture histories
sum(y.vect)   # Should be 270.

# EITHER fit 16 models in dividually, assuming two groups for 
# the heterogeneous models:

phic.pc.out <- hetage.fit.model("phic.pc")
phic.pt.out <- hetage.fit.model("phic.pt")
phic.ph.out <- hetage.fit.model("phic.ph",G=2)  
phic.pth.out <- hetage.fit.model("phic.pth",G=2)
phit.pc.out <- hetage.fit.model("phit.pc")
phit.pt.out <- hetage.fit.model("phit.pt")
phit.ph.out <- hetage.fit.model("phit.ph",G=2)
phit.pth.out <- hetage.fit.model("phit.pth",G=2)
phih.pc.out <- hetage.fit.model("phih.pc",G=2)
phih.pt.out <- hetage.fit.model("phih.pt",G=2)
phih.ph.out <- hetage.fit.model("phih.ph",G=2)
phih.pth.out <- hetage.fit.model("phih.pth",G=2)
phith.pc.out <- hetage.fit.model("phith.pc",G=2)
phith.pt.out <- hetage.fit.model("phith.pt",G=2)
phith.ph.out <- hetage.fit.model("phith.ph",G=2)
phith.pth.out <- hetage.fit.model("phith.pth",G=2)


# OR set up a list:

model.names <- c("phic.pc","phic.pt","phic.ph","phic.pth",
                 "phit.pc","phit.pt","phit.ph","phit.pth",
                 "phih.pc","phih.pt","phih.ph","phih.pth",
                 "phith.pc","phith.pt","phith.ph","phith.pth")
no.models   <- length(model.names)
model.fits  <- vector("list",no.models)
names(model.fits) <- model.names
ngroups <- c(1,1,2,2,1,1,rep(2,10))

# Fit the models:

for (modno in 1:no.models)
   {
   model.fits[[modno]] <- hetage.fit.model(model.names[modno],G=ngroups[modno])
   print(paste("Model",model.names[modno],"fitted"))
   }

hetage.summary(model.fits)             # for a detailed summary table
hetage.summary(model.fits,AICorder=T)  # summary table sorted by AIC
hetage.summary(model.fits,AICcorder=T) # summary table sorted by AICc

round(hetage.summary(model.fits)[,4],1)

# Best model found is phit.ph.
# For one model, found a better solution than in table in 
# Pledger et al. 2010.

possums9.out <- model.fits



# Extra models:

model.fits <- possums9.out

phia.pc.out <- hetage.fit.model("phia.pc")
phia.pt.out <- hetage.fit.model("phia.pt")
phia.ph.out <- hetage.fit.model("phia.ph",G=2)
phia.pth.out <- hetage.fit.model("phia.pth",G=2)

model.fits[[17]] <- phia.pc.out
model.fits[[18]] <- phia.pt.out
model.fits[[19]] <- phia.ph.out
model.fits[[20]] <- phia.pth.out
names(model.fits)[17:20] <- 
   c("phia.pc","phia.pt","phia.ph","phia.pth")

hetage.summary(model.fits)

possums9.out <- model.fits
