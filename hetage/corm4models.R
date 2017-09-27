# corm4models.R

# Commands for analysis of cormorant data, 317 birds, 8 months.

# Four models are fitted, the JSSA model and simplifications
# with phi and/or p constant.

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


# Option 1: fit the four models separately:
# -----------------------------------------

phic.pc.out <- hetage.fit.model("phic.pc")
phic.pt.out <- hetage.fit.model("phic.pt")
phit.pc.out <- hetage.fit.model("phit.pc")
phit.pt.out <- hetage.fit.model("phit.pt")

phic.pc.out   # To inspect the output.


# Option 2: set up a list of model fits.
# --------------------------------------

# Name the models:

model.names <- c("phic.pc","phic.pt","phit.pc","phit.pt") 

# Set up the list:
 
no.models   <- length(model.names)
model.fits  <- vector("list",no.models)
names(model.fits) <- model.names

# Fit the models, store in the list:

for (modno in 1:no.models)
   {
   model.fits[[modno]] <- hetage.fit.model(model.names[modno])
   print(paste("Model",model.names[modno],"fitted"))
   }

# Construct a summary table:

summary.table <- matrix(NA,no.models,7)
dimnames(summary.table) <- list(model.names,
   c("MaxLL","RD","npar","AIC","relAIC","AICc","relAICc"))

for (modno in 1:no.models) if (!is.null(model.fits[[modno]]))
   {
   summary.table[modno,1] <- model.fits[[modno]]$maxLL
   summary.table[modno,2] <- model.fits[[modno]]$RD
   summary.table[modno,3] <- model.fits[[modno]]$npar
   summary.table[modno,4] <- model.fits[[modno]]$AIC
   summary.table[modno,6] <- model.fits[[modno]]$AICc
   }
summary.table[,5] <- summary.table[,4] - min(summary.table[,4],na.rm=T)
summary.table[,7] <- summary.table[,6] - min(summary.table[,6],na.rm=T)

print(round(summary.table,2))

# Save the results:

corm4models.out <- model.fits

