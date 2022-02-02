# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Functions
# ====================================

library(moments)
library(dplyr)
library(sp)
library(adehabitatHR)
library(rgeos)
library(raster)

' --- Generate output for theoretical experiments --- '
# calculate home-range size and mean resource within hr 
GenOutput_Ex1 <- function(tracksfile="Ch1EX1S3(5)", rep_index = 0:100, r = r, D = 5){
  KDE95 <- KDE50 <- MeanRes <- consumption <- vector()
  for (i in rep_index){ 
    df <- read.csv(paste(tracksfile,'/Tracks',i,'.csv', sep=''), header = T)
    df <- df[, c('id','x','y')]
    coordinates(df) <- ~x+y
    KUD <- kernelUD(df)
    # calculate and store home-range estimates
    KUD95 <- kernel.area(KUD, percent = 95, unin = 'm', unout = 'm2')
    KUD95v <- as.vector(apply(KUD95,2, function(x) as.numeric(x)))
    KDE95 <- c(KDE95, KUD95v)
    KUD50 <- kernel.area(KUD, percent = 50, unin = 'm', unout = 'm2')
    KUD50v <- as.vector(apply(KUD50,2, function(x) as.numeric(x)))
    KDE50 <- c(KDE50, KUD50v)
    # extract home-range polygons, calculate MeanResource
    homepolys95 <- getverticeshr(KUD,percent = 95)
    r.vals <- extract(r, homepolys95)
    r.mean <- lapply(r.vals, FUN=mean)
    r.mean<-do.call(rbind,r.mean)
    MeanRes <- c(MeanRes, r.mean)
  }
  output <- data.frame(replicate = rep(rep_index, each = D), HRA95 = KDE95, HRA50 = KDE50, MeanResource = MeanRes, 
                       Density = rep(D, length(KDE95)))
  return(output)
}

HR_poly <- function(df=NA, KUD = 95){
  df <- df[, c('id','x','y')]
  coordinates(df) <- ~x+y
  kud <- kernelUD(df)
  hrpoly <- getverticeshr(kud,percent = KUD)
  return(hrpoly)
}

# create home-range polygons for phase 1 and phase 2 and calculate overlap in area between remaining and removed individuals
GenOutput_Ex2 <- function(tracksfile = 'Expt2', file_index = 1:50, 
                          phase1 = c(5001,20000), phase2 = c(35001,50000), H1 = TRUE){
  require(adehabitatHR)
  require(rgeos)
  P1 <- P2 <- vector() # to store output
  for (i in file_index){
    df <- read.csv(paste(tracksfile,'/Tracks',i,'.csv', sep=''), header = T)
    if (H1 == "TRUE"){  # H1: removals
      survivors <- 
        df %>% 
        group_by(id) %>% 
        filter(row_number() == 50000) %>% 
        pull(id) # who are the survivors?
    } else { # H0: control
      survivors <- sample(unique(df$id), 6, replace = F) # pick random indivs to survive
    }
    
    # arange data into phases and who is removed/not
    before.alive <-
      df %>% 
      filter(id %in% survivors) %>% 
      group_by(id) %>% 
      filter(row_number() %in% phase1[1]:phase1[2])
    
    P1.alive <- HR_poly(before.alive)
    
    before.dead <-
      df %>% 
      filter(!id %in% survivors) %>% 
      group_by(id) %>% 
      filter(row_number() %in% phase1[1]:phase1[2])
    
    P1.dead <- HR_poly(before.dead)
    
    after.alive <-
      df %>% 
      filter(id %in% survivors) %>% 
      group_by(id) %>% 
      filter(row_number() %in% phase2[1]:phase2[2])
    P2.alive <- HR_poly(after.alive)
    
    # calculate overlap between the areas of the removed in P1 (P1.dead) and the remainders in P1 and P2 (alive)
    P1.overlap <- gIntersection(P1.dead, P1.alive, byid=FALSE)
    P1[i+1] <- P1.overlap@polygons[[1]]@area
    
    P2.overlap <- gIntersection(P1.dead, P2.alive, byid=FALSE) 
    P2[i+1] <- P2.overlap@polygons[[1]]@area
  }
  return(list(P1,P2))
}

' --- Calculate individual-level output metrics --- '
# This version with new pairwise HRA overlap metric 
# Warning: still need to enter new metric into Output Matrix
CalcMetrics <- function(tracksfile="SAcats", file_index = 1:1000, format_num = FALSE, n_sims=1, burnin = 0, numsteps=21600, n_foragers=30){
  StudyArea <- Polygon(matrix(c(0, 0, 400, 0, 400,600, 0, 600 , 0, 0), ncol = 2, byrow=T))
  StudyArea <- SpatialPolygons(list(Polygons(list(StudyArea), ID = "a")))
  skips <-  seq(burnin+15, numsteps, by = 15)
  for (k in file_index){
    if (format_num == "TRUE"){
      id <- sprintf("%02d", k)
    } else {
      id <- k
    }
    summary <- read.csv(paste(tracksfile,id,'.csv', sep=''), header = T)
    rep_index <- 1:n_sims 
    for (i in rep_index){ 
      rownum <- i
      skip_to_next <- FALSE
      tryCatch( {
        df <- read.csv(paste(tracksfile, id,'/Tracks',i,'.csv', sep=''), header = T)} , 
        error = function(e) {
          skip_to_next <- TRUE
          print( paste(tracksfile,id,'/Tracks',i, ' error', sep='') ) })
      if(skip_to_next) {next}
      
      df$x<- df$x - 50 # correct for empty buffer
      df$y<- df$y - 50
      
      df %>% 
        group_by(id) %>% 
        filter(row_number() %in% skips) -> df
      df <- df[, c('id','x','y')]
      coordinates(df) <- ~x+y
      KUD <- kernelUD(df)
      # calculate and store home-range estimates
      KUD95 <- kernel.area(KUD, percent = 95, unin = 'm', unout = 'm2')
      KDE95 <- as.vector(apply(KUD95,2, function(x) as.numeric(x))) / 10000
      KUD50 <- kernel.area(KUD, percent = 50, unin = 'm', unout = 'm2')
      KDE50 <- as.vector(apply(KUD50,2, function(x) as.numeric(x))) / 10000
      # HRA overlap
      homerange <- getverticeshr(KUD,percent = 95)
      
      overlap.prop <- list()
      for(poly in 1:n_foragers){
        overlap.id <- gIntersects(homerange[poly,], homerange[-poly,], byid=TRUE)
        overlap.pair <- vector()
        overlap.pair[which(overlap.id == "FALSE")] <- 0
        
        if (any(overlap.id == "TRUE")) {
          for (j in which(overlap.id == "TRUE")){
            overID <- which(homerange@data[["id"]] == rownames(overlap.id)[j])
            overArea <- gIntersection(homerange[poly,], homerange[overID,], byid = FALSE, id = paste(poly))
            mergedArea <- gUnion(homerange[poly,], homerange[overID,], byid = FALSE)
            propArea <- gArea(overArea , byid= FALSE) / gArea(mergedArea , byid= FALSE)
            overlap.pair[j] <- propArea
          } #j 
          overlap.prop[[poly]] <- overlap.pair
        } else{ overlap.prop[[poly]] <- 0} # check if any overlap at all
      } # poly
      
      IndMedian <- lapply(overlap.prop, median)
      IndMax <- lapply(overlap.prop, max)
      OverlapMedian <- unlist(IndMedian)
      OverlapMax <- unlist(IndMax)
      
      # step-length (skewness)
      df$time <- rep(seq(1, length(skips)*15, by=15), n_foragers) 
      class(df$time) <- c("POSIXct", "POSIXt")
      attr(df$time, "tzone") <- ""
      df.lt <- as.ltraj(df@coords, date=df$time, id=df$id, typeII = TRUE, burst = df$id)
      Skewness.step <- unlist(lapply(df.lt, function(x) skewness(x[['dist']], na.rm = T) )) # previous error - divided by 100
      Median.step <- unlist(lapply(df.lt, function(x) summary(x[['dist']], na.rm = T)['Median'] ))/100
      
      # daily distance -- calc stepl for each 15-min interval and sum per day (48)
      fac <- rep(1:n_foragers, each = 48) # Create a "group" index
      dielmvmt <- lapply(df.lt, function(x) tapply(x[['dist']], fac, sum))
      dd.median <- unlist(lapply(dielmvmt, function(x) summary(x, na.rm = T)['Median'] ))/100
      dd.mean <- unlist(lapply(dielmvmt, function(x) mean(x, na.rm = T) ))/100
      dd.sd <- unlist(lapply(dielmvmt, function(x) sd(x, na.rm = T) ))/100
      
      # Replicate benchmarks
      Skewness.HRA <- skewness(KDE95, na.rm = T)
      Median.HRA <- median(KDE95, na.rm = T)
      AllHRA <- gUnionCascaded(homerange, id = NULL)
      AreaCovered <- gIntersection(AllHRA, StudyArea, byid = FALSE, id = 'moo')
      LandCovered <- gArea(AreaCovered, byid = FALSE) / (400*600)
      
      # Independent vars 
      independents <- summary[rownum, c('ResourceRegenerationRate',
                                        'ForagerSpeedSearch',
                                        'ForagerSpeedFeeding' ,
                                        'LongDecayRate',
                                        'ShortDecayRate',
                                        'MemoryValueUninformed',
                                        'ScentDecayRate',
                                        'ScentResponseSpatialScale',
                                        'ScentDepositionRate')]
      output.store <- data.frame(filename = rep(paste(tracksfile,id, sep=''), n_foragers), 
                                 replicate_no = rep(i, n_foragers),
                                 HRA95 = KDE95, HRA50 = KDE50, Skew.HRA = rep(Skewness.HRA, n_foragers), Median.HRA = rep(Median.HRA, n_foragers),
                                 Skew.step = Skewness.step, Median.step = Median.step, 
                                 DailyDist.mean = dd.mean, DailyDist.median = dd.median, DailyDist.sd = dd.sd,
                                 OverlapMedian = OverlapMedian, OverlapMax = OverlapMax, LandCovered = rep(LandCovered, n_foragers), 
                                 ResourceRegenerationRate = rep(independents[,1], n_foragers),
                                 ForagerSpeedSearch = rep(independents[,2], n_foragers),
                                 ForagerSpeedFeeding = rep(independents[,3], n_foragers),
                                 LongDecayRate = rep(independents[,4], n_foragers),
                                 ShortDecayRate = rep(independents[,5], n_foragers),
                                 MemoryValueUninformed = rep(independents[,6], n_foragers),
                                 ScentDecayRate = rep(independents[,7], n_foragers),
                                 ScentResponseSpatialScale = rep(independents[,8], n_foragers),
                                 ScentDepositionRate = rep(independents[,9], n_foragers) )
      if(i == 1 ){
        output.k <- output.store
      } else{
        output.k <- rbind(output.k, output.store) }
      
    } #i replicates
    if(k == file_index[1]){
      output <- output.k
    } else {
      output <- rbind(output, output.k)
    }
    
  } #k
  return(output)
}

# cumulative HR and breakpoint analysis (TAC: total area covered)
TAC <- function(tracksfile="SAcats", n_steps = 21600, interval = 720, iteration = 1){
  df <- read.csv(paste(tracksfile,'/Tracks',iteration,'.csv', sep=''), header = T)
  df %>% 
    group_by(id) %>% 
    filter(row_number() %in% skips) -> df
  output <- matrix(NA, length(skips), length(unique(df$id)))
  i = 1
  subset_int <- seq(interval, n_steps, by=interval)
  for (q in subset_int){ 
    df.trunc <-
      df %>% 
      group_by(id) %>% 
      filter(row_number()<= q)
    df.trunc <- df.trunc[, c('id','x','y')]
    coordinates(df.trunc) <- ~x+y
    MCP <- mcp.area(df.trunc, percent = 95, unin = 'm', unout = 'm2', plotit = F)
    # calculate and store home-range estimates
    MCP95v <- as.vector(apply(MCP[1,],2, function(x) as.numeric(x)))
    output[i,] <- MCP95v
    i = i + 1
  }
  return(output)
}
break_pt_est <- function(Y = realTAC[,1], X = 1:30, init.val = 5){
  require(segmented)
  mylm <- lm(Y ~ X)
  seg.lm <- segmented(mylm, seg.Z = ~X, psi= init.val)
  result <- summary(seg.lm)[["psi"]][2]
  return(result)
}
TACreal <- function(df=df, n_days = 30, int_days = 5, id_col = "ID", x_name = "lon", y_name = "lat"){
  skips <- seq(int_days, n_days, by=int_days)
  output <- matrix(NA, length(skips), length(unique(df$ID)))
  
  i = 1
  for (q in skips){ 
    df.trunc <-
      as.data.frame(df) %>% 
      group_by(ID)%>% 
      filter(time >= min(time) & time <= (min(time) + q*24*60*60) )
    
    df.trunc <- df.trunc[, c( paste(id_col), paste(x_name), paste(y_name))]
    
    coordinates(df.trunc) <- ~ lon + lat

    MCP <- mcp.area(df.trunc, percent = 100, unin = 'm', unout = 'm2', plotit = F)
    MCP95v <- as.vector(apply(MCP[1,],2, function(x) as.numeric(x)))
    
    output[i,] <- MCP95v
    i = i+1
  }
  return(output)
}

' --- Calculate replicate-level summary statistics for Sensitivity Analysis --- '
SummMetrics <- function(data,log.scale = TRUE){
  metrics <- cbind(cbind(aggregate(data[, c("HRA95","HRA50","Skew.HRA","Median.HRA",
                            "Skew.step","Median.step","DailyDist.mean","DailyDist.median","DailyDist.sd",
                            "OverlapMedian","OverlapMax","LandCovered")], list(data$filename), median)),
                   HRA.min = aggregate(data[, "HRA95"], list(data$filename), min),
                   aggregate(data[, c('ResourceRegenerationRate',
                                      'ForagerSpeedSearch',
                                      'ForagerSpeedFeeding',
                                      'LongDecayRate',
                                      'ShortDecayRate',
                                      'MemoryValueUninformed',
                                      'ScentDecayRate',
                                      'ScentResponseSpatialScale',
                                      'ScentDepositionRate')], list(data$filename), mean))
  
  output <- metrics[,-c(14,16)]
  
  if (log.scale==TRUE){
    output$ResourceRegenerationRate <- log10(output$ResourceRegenerationRate)
    output$LongDecayRate <- log10(output$LongDecayRate)
    output$ShortDecayRate <- log10(output$ShortDecayRate)
    output$ScentDecayRate <- log10(output$ScentDecayRate)
  }
  
  return(output)
}

' --- Evaluation metrics --- '
# Function for Root Mean Squared Error
RMSE <- function(error) { sqrt(mean(error^2)) }

# Function for extracting CV stats
CVstats <- function(model) {
  if (is.list(model)){
    cv.mean <- lapply(model, function(x) {x[["cv.statistics"]][["deviance.mean"]]} )
    cv.se <- lapply(model, function(x) {x[["cv.statistics"]][["deviance.se"]]} )
    int.depth <- lapply(model, function(x) {x[["interaction.depth"]]} )
    output <- data.frame(interaction.depth = unlist(int.depth), cv.mean = unlist(cv.mean), cv.se = unlist(cv.se))
    return(output)} else {
      return(paste('error: not a list yo'))
    }
}

' --- Run BRT --- '
# function to iterate BRT over tree complexity 1-10
runBRT <- function(data, predictor.var = 'HRA95'){
  brt <- list()
  x <- which(colnames(data) %in% c('ResourceRegenerationRate',
                                      'ForagerSpeedSearch',
                                      'ForagerSpeedFeeding',
                                      'LongDecayRate',
                                      'ShortDecayRate',
                                      'MemoryValueUninformed',
                                      'ScentDecayRate',
                                      'ScentResponseSpatialScale',
                                      'ScentDepositionRate'))
  y <- which(colnames(data) == predictor.var)
  for (i in 1:10){ # for tree complexity 1-10
    brt[[i]] <- gbm.step(data = data, gbm.x = x, gbm.y = y, n.folds = 5, 
                    tree.complexity = i , family = "gaussian", n.trees = 15)
    if(brt[[i]][["cv.statistics"]][["deviance.mean"]] < 0.01){
      stop
    }
  }
  return(brt)
}

# function to summarise BRT resilts
summariseBRT <- function(model, output.names = c('Steplength', 'StepSkew', "DailyDist", "HRA", "HRASkew")) {
  target <- model[[1]][["var.names"]]
  output <- matrix(NA, length(target), length(model))
  for (i in 1:length(model)){
    x <- model[[i]]
    df <- x[["contributions"]]
    output[,i] <- df[match(target, df$var),][,2]
  }
  rownames(output) <- target
  colnames(output) <- output.names
  return(output)
}


' --- Plotting --- '
# functions to allow plotting of trajecteories with diminishing opacity over time
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}
colorRampPaletteAlpha <- function(colors, n=100, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
  if (interpolate=='linear') {
    l <- approx(a, n=n)
  } else {
    l <- spline(a, n=n)
  }
  l$y[l$y > 255] <- 255 # Clamp if spline is > 255
  cr <- addalpha(cr, l$y/255.0)
  return(cr)
}
# plot trajectories
plotTraject <- function( df = df, colfunc = colorRampPalette(c("white","black")), num_steps = 25000){
  n_foragers = length(unique(df$id))
  k=0
  for(i in unique(df$id)){
    k=k+1
    foragercols.start <- addalpha(colfunc(n_foragers)[k],0)
    foragercols.end <- addalpha(colfunc(n_foragers)[k],1)
    foragecols <- c(foragercols.start, foragercols.end)
    crpa <- colorRampPaletteAlpha(foragecols, num_steps, interpolate = 'linear')
    segments(x0 = head(df$x[df$id==i],-1),
             y0 = head(df$y[df$id==i],-1),
             x1 = tail(df$x[df$id==i],-1),
             y1 = tail(df$y[df$id==i],-1),
             lwd = 1,
             col = crpa)
  }
}

# plot cv deviance
plot_brt_perf <- function(cv, pred.name = 'Home-range size'){
  for (i in 2:10){
    if (cv[i,2] + cv[i,3] < cv[i-1,2] - cv[i-1,3]){
      next
    } else {
      x_intercept <- i - 1
      y_intercept <- cv[i-1,2] + cv[i-1,3]
      break
    }
  }
  p <- ggplot(cv, aes(x = interaction.depth, y = cv.mean)) +
    geom_errorbar(aes(ymin = cv.mean - cv.se, ymax = cv.mean + cv.se), colour = "lightskyblue2", width = .1) +
    geom_line(colour = "darkgoldenrod3") +
    geom_point(size=2, shape=21, colour = "darkgoldenrod3", fill="darkgoldenrod3") +
    geom_vline(xintercept = x_intercept, linetype="dashed", 
               color = "purple4", size = .25) +
    geom_hline(yintercept = y_intercept, linetype="dashed", 
               color = "purple4", size = .25) +
    scale_x_continuous(breaks = 0:9) + labs(title = pred.name) + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}
