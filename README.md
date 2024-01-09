# _Ambystoma bishopi_ Survival Analysis
Supporting data and code for:
Brooks G. C., T. A. Gorman, and C. A. Haas. 2024. Variation in Flatwoods Salamander Survival Is Unrelated to Temperature and Rainfall. Ichthyology and Herpetology.

Using ten years of mark-recapture data from two adjacent wetlands on the Florida Panhandle, we investigated individual and temporal variability in survival rates of Reticulated Flatwoods Salamanders (Ambystoma bishopi). Our objectives were to 1) provide the first estimates of survival for the species, 2) evaluate the relationship between body size and mortality risk, 3) quantify the degree of variability in survival rates across the study period, and 4) discern whether variability in survival or detection correlates with environmental conditions. To address these objectives, we constructed a modified Cormack-Jolly-Seber model that includes body size and year as covariates. Mean annual survival was estimated to be 0.72 and was strongly correlated with body size; survival rates of the smallest individuals in the study were 0.5 and those of the largest individuals were 0.85. Survival also varied considerably across years, but did not correlate with temperature extremes or rainfall. 

## Contents
### Data
**summer_temp.csv** - mean summer temperatures across the study period\
**winter_temp.csv** - mean winter temperatures across the study period\
**summer_ppt.csv** - mean summer rainfall across the study period\
**winter_ppt.csv** - mean winter rainfall across the study period\
**fw_dat.csv** - Flatwoods salamander capture history dataset

### Code
**surv_23_model.R** - JAGS parameterization of the Cormack-Joly-Seber model used to estimate survival\
**surv_23_run.R** - code to format the raw capture data, run the CJS model, and derive relationships between environmntal covariates and vital rates. 
