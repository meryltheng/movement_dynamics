# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Calibration: predict and optimize
# ====================================

' --- Calculate target values (real data) --- '
cats <- read.csv('Data/realdata.csv')
#catsSub$x <- (catsSub$lon - 763421.1)/10 - 100 # convert to xy starting from (0,0) and scale to resolution of simulation landscape
#catsSub$y <- (catsSub$lat - 6029883.9)/10 - 100
load(file='Data/catsLT.RData')

# HR size
cats.sp <- cats[, c('ID','lon','lat')] # don't mess up original df for later analyses
coordinates(cats.sp ) <- ~lon+lat
KUD <- kernelUD(cats.sp)
KUD95 <- kernel.area(KUD, percent = 95, unin = 'm', unout = 'km2')
HRA <- as.vector(apply(KUD95,2, function(x) as.numeric(x))) 
summary(HRA)
skewness(HRA)

# Daily distance
DailyDist <- lapply(cats.lt, function(x) tapply(x$dist, as.Date(x$date), FUN=sum))
DailyDist.median <- unlist(lapply(dd, function(x) summary(x, na.rm = T)['Median'] ))/1000
summary(DailyDist.median)

# Steplength
Steplength <- unlist(lapply(cats.lt, function(x) summary(x[['dist']], na.rm = T)['Median'] ))/1000
summary(Steplength)
Steplength.skew <- unlist(lapply(cats.lt, function(x) skewness(x[['dist']], na.rm = T) ))
summary(Steplength.skew)

# Target statistics:
# Steplength.median = 0.07 
# Steplength.skew = 2.2 
# DailyDist.mean/median ~ 4.8
# HRA = 1.9 (min = 1.3)
# HRAskew = 1.8

' --- Determine initial values with BRT models --- '
#### PREDICT ####
Cands <- c(-3, 5, 0.05, -7, -3, 5e-5, -6, 3, 2) # candidate parameter set
# to help tuning process
predict(brt.Steplength[[1]], newdata = data.frame(ResourceRegenerationRate = Cands[1], # change BRT model for diff statistic
                                           ForagerSpeedSearch = Cands[2],
                                           ForagerSpeedFeeding = Cands[3],
                                           LongDecayRate = Cands[4],
                                           ShortDecayRate = Cands[5],
                                           MemoryValueUninformed = Cands[6],
                                           ScentDecayRate = Cands[7],
                                           ScentResponseSpatialScale = Cands[8],
                                           ScentDepositionRate = Cands[9]))


# OPTIMISE [warning: error-prone, especially if inits are random]
# functions
predMMoutput <- function(input = NULL, models = NULL, trueVals = NULL){
  propDiff <- vector()
  for (i in 1:length(models)){
    predVal <- predict(models[[i]], newdata = data.frame(ResourceRegenerationRate = input[1], 
                                                         ForagerSpeedSearch = input[2],
                                                         ForagerSpeedFeeding = input[3],
                                                         LongDecayRate = input[4],
                                                         ShortDecayRate = input[5],
                                                         MemoryValueUninformed = input[6],
                                                         ScentDecayRate = input[7],
                                                         ScentResponseSpatialScale = input[8],
                                                         ScentDepositionRate = input[9]))
    
    propDiff[i] <- abs((predVal - trueVals[i]) / trueVals[i])
  } # i
  output = sum(propDiff)^2
  output
}

inits <- function(min = c(-3, 1, 0, -7, -4, 0, -6, 0.5, 0), # for random inits, which we did not use
                  max =c(-1, 6, 0.5, -3, -1, 2e-4, -2, 10, 5)){
  output <- apply(cbind(min,max), 1, function(x){runif(1, x[1], x[2])})
  return(output)
}

# run optim (with candidate parameter set)
result <- optim(par = Cands, fn = predMMoutput, models = list(brt.Steplength[[1]], brt.StepSkew[[2]], brt.DailyDist[[3]], brt.HRA[[3]], brt.HRAskew[[1]]), 
              trueVals = c(0.08, 2.4, 4.8, 1.9, 1.8), # target output values
              lower = c(-3, 1, 0, -7, -4, 0, -6, 0.5, 0), # parameter lower bound
              upper = c(-1, 6, 0.5, -3, -1, 2e-4, -2, 10, 5), # parameter upper bound
              method = 'L-BFGS-B' )
