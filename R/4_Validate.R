# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Compare real to simulated data
# ====================================

' --- Calculate approximated values (simulated data) --- '

# calculate individual-level stats for final iteration
simMetrics08 <- CalcMetrics(tracksfile="SAvalid/SAoptim08_", file_index = 1:10, format_num = FALSE, n_sims=1, burnin = 0, numsteps=21600, n_foragers=30)

# quick summary statistics calculation across individuals from all replicates
summary(simMetrics08$HRA95) # summarise HR size
skewness(simMetrics08$HRA95) # calculate HR skewness
summary(simMetrics08$DailyDist.median) # summarise daily distance (median)
summary(simMetrics08$Median.step) # summarise step-length (median)

# summary statistics calculation (by replicate)
datValid08 <- SummMetrics(simMetrics08,log.scale = TRUE)
colMeans(datValid08[, c("HRA95","Skew.HRA","Skew.step","Median.step","DailyDist.mean","DailyDist.median","DailyDist.sd",
                   "OverlapMax","LandCovered")]) # calculate mean across all replicates
apply(datValid08[, c("HRA95","Skew.HRA","Skew.step","Median.step","DailyDist.mean","DailyDist.median","DailyDist.sd",
                        "OverlapMax","LandCovered")], 2, sd) # calculate sd across all replicates

' --- plot HR size distribution(Fig. 5a) --- '
# set-up plot first (empty)
simMetrics = simMetrics08; i = unique(simMetrics$filename)[1]; randoCats <- sample(1:30, 11, replace=F) 
hist(simMetrics$HRA95[simMetrics$filename == i][randoCats], #breaks = seq(0,180,5),
     main = '',col = rgb(0,0,0,0), xlim = c(0,17), ylim = c(0,8),
     breaks =seq(0,ceiling(max(simMetrics$HRA95)/2)*2,1), border=rgb(0,0,0,0.0), 
     xlab = expression(paste("Home-range size (km"^2,')', sep='')), ylab = 'Frequency')

# bootstrap over all replicates 10 times (i.e., run the following loop 10x)
for (i in unique(simMetrics$filename)){
  randoCats <- sample(1:30, 11, replace=F) # select 11 random individuals
  hist(simMetrics$HRA95[simMetrics$filename == i][randoCats], #breaks = seq(0,180,5),
       main = '', xlab='',ylab=' ',
       col = rgb(0,0,0,0.05), breaks =seq(0,ceiling(max(simMetrics$HRA95)/2)*2,1), border=rgb(0,0,0,0.0), add=T)
}
# plot real data (note that data has been read in by previous section '3_calibrate.R')
hist(HRA, col=rgb(0,0,0,0.2), breaks =seq(0,ceiling(max(HRA)/2)*2,1), border='red',lwd=1.5, add=T)
box(bty="l")

# do the same for daily distance (for reference, not in main paper)
hist(simMetrics$DailyDist.median[randoCats], #breaks = seq(0,180,5),
     main = '',col = rgb(0,0,0,0), xlim = c(0,10), ylim = c(0,8),
     breaks =seq(0,ceiling(max(simMetrics$DailyDist.median)/2)*2,1), border=rgb(0,0,0,0.0), 
     xlab = "Daily distance (km)", ylab = 'Frequency')
for (i in unique(simMetrics$filename)){
  randoCats <- sample(1:30, 11, replace=F) # select 11 random individuals
  hist(simMetrics$DailyDist.median[simMetrics$filename == i][randoCats], #breaks = seq(0,180,5),
       main = '', xlab='',ylab=' ',
       col = rgb(0,0,0,0.05), breaks =seq(0,ceiling(max(simMetrics$DailyDist.median)/2)*2,1), border=rgb(0,0,0,0.0), add=T)
}
hist(DailyDist.median, col=rgb(0,0,0,0), breaks =seq(0,ceiling(max(DailyDist.median)/2)*2,1), border='red',add=T)
box(bty="l")

' --- plot cumulative HR size (Fig. 5b) --- '
# simulated data
# par(mar=c(5.1, 4.8, 4.1, 2.1))
i <- sample(1:10,1) # pick random iteration to plot
tracksfile = paste('SAvalid/SAoptim08_',i,sep='') # i = 3
simTACo8 <- TAC(tracksfile=tracksfile, n_steps = 21600/15, interval = 720/15, iteration = 1) # MCP95
matplot(simTACo8/10000, type = "l", ylab = expression(paste("Cumulative home-range size (km"^2,')', sep='')), 
        xlab = "Days", main = '', col = rainbow(11), frame=F, lty = 1, ylim = c(0,16), cex.lab=1.3)
box(bty="l")

# real data
cats$time <- as.POSIXct(paste(cats$time), format="%Y-%m-%d %H:%M:%S", tz="Australia/Adelaide")
realTAC <- TACreal(df=cats, n_days = 30, int_days = 1, id_col = "ID", x_name = "lon", y_name = "lat")

end <- # determine start and end to location data 
  as.data.frame(cats) %>% 
  group_by(ID)%>% 
  summarise(max.time = max(time) )

start <-
  as.data.frame(cats) %>% 
  group_by(ID)%>% 
  summarise(min.time = min(time) )

daysalive <- floor(end$max.time - start$min.time)
daysalive[which(daysalive > 30)] <- 30

for (i in 1:11){ # truncate dataset to days tracked/alive
  realTAC[-(1:daysalive[i]),i] <- NA
}

matplot(realTAC, type = "l",  ylab = "", xlab = "Days", main = '', col = rainbow(11), frame=F, lty = 1, ylim = c(0,16))
box(bty="l")

' --- plot individual locations (Fig. 5c) --- '
env <- read.csv('Data/CatLand.csv', header = T)
#colfunc <- colorRampPalette(c("white","black"))
#par(mar=c(1, 1, 1, 1))
image(1:nrow(env), 1:ncol(env), as.matrix(env), col = colfunc(10), axes = F, xlab = "", ylab = "")

# sim
df <- read.csv('SAvalid/SAoptim07/Tracks1.csv', header = T) # same replicate and individuals as cumulative HR plot
df$x <- df$x - 50; df$y <- df$y - 50
skips <-  seq(15, 21600, by = 15)
df %>% 
  group_by(id) %>% 
  filter(row_number() %in% skips) -> df

k = 0
for(i in randoCats){ # points
  k=k+1
  points(x=df$x[df$id==i],
         y=df$y[df$id==i],
         pch=3,
         col = alpha(rainbow(11)[k],0.2))
} 
box(lty=1, col='black')

#randoCats <- sample(unique(df$id), 11, replace=F)
for(i in unique(df$id)){ # points
  k=k+1
  if (i %in% randoCats){
  points(x=df$x[df$id==i],
         y=df$y[df$id==i],
         pch=3,
         col = alpha(rainbow(30)[k],0.2))
  } else {
    next
  }
} 

for(i in unique(df$id)){ # points
  k=k+1
  points(x=df$x[df$id==i],
         y=df$y[df$id==i],
         pch=3,
         col = alpha(rainbow(30)[k],0.2))
} 

# real
image(1:nrow(env), 1:ncol(env), as.matrix(env), col = colfunc(10), axes = F, xlab = "", ylab = "")
k = 0
for(i in unique(catsSub$ID)){ # points
  k=k+1
  points(x=catsSub$x[catsSub$ID==i],
         y=catsSub$y[catsSub$ID==i],
         pch=3,
         col = alpha(rainbow(11)[k],0.2))
} 
box(lty=1, col='black')
