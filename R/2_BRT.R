# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Boosted regression trees (Sensitivity Analysis)
# ====================================

library(dismo)
library(gbm)
source("R/0_functions.R")

' --- Calculate individual-level output metrics --- '
system.time( # takes a while
  Metrics <- CalcMetrics(tracksfile="SAcats", file_index = 1:1000, format_num = FALSE, n_sims=1, burnin = 0, numsteps=21600, n_foragers=30)
)
write.csv(MetricsALL,'Data/SA_indivMetrics.csv')
#Metrics <- read.csv('Data/SA_indivMetrics.csv')

' --- Calculate replicate-level output metrics for Sensitivity Analysis --- '
dat <- SummMetrics(Metrics,log.scale = TRUE)
#dat <- read.csv('Data/SA_repMeans.csv')

' --- Run BRT --- '
brt.Steplength <- runBRT(dat, predictor.var = 'Median.step')
brt.StepSkew <- runBRT(dat, predictor.var = 'Skew.step')
brt.DailyDist <- runBRT(dat, predictor.var = 'DailyDist.mean')
brt.HRA <- runBRT(dat, predictor.var = 'HRA95')
brt.HRAmin <- runBRT(dat, predictor.var = 'HRA.min.x')
brt.HRAskew <- runBRT(dat, predictor.var = 'Skew.HRA')
brt.OverlapMax <- runBRT(dat, predictor.var = 'OverlapMax')
brt.Coverage <- runBRT(dat, predictor.var = 'LandCovered')

#save(brt.Steplength,brt.StepSkew,brt.DailyDist,brt.HRA,brt.HRAskew,brt.OverlapMax,brt.HRAmin,file='R/BRTresults.RData')
load('Data/BRTresults.RData')

' --- Plot CV --- '
library(ggplot2)
library(ggpubr)
#cv.Steplength <- CVstats(model=brt.Steplength) # cv < 0.01, use first model
#c1 <- plot_brt_perf(cv.Steplength, pred.name = 'Step-length')
cv.Stepskew <- CVstats(model=brt.StepSkew)
c2 <- plot_brt_perf(cv.Stepskew, pred.name = 'Step-length skew')
cv.DailyDist <- CVstats(model=brt.DailyDist)
c3 <- plot_brt_perf(cv.DailyDist, pred.name = 'Daily distance')
cv.HRA <- CVstats(model=brt.HRA)
c4 <- plot_brt_perf(cv.HRA, pred.name = 'Home-range size')
cv.HRAskew <- CVstats(model=brt.HRAskew)
c5 <- plot_brt_perf(cv.HRAskew, pred.name = 'Home-range size skew')
#cv.OverlapMax <- CVstats(model=brt.OverlapMax) # cv < 0.01
#c6 <- plot_brt_perf(cv.OverlapMax, pred.name = 'Maximum HR overlap')
#cv.Coverage <- CVstats(model=brt.Coverage) # cv < 0.01
#c7 <- plot_brt_perf(cv.Coverage, pred.name = '% Land covered')
cv.HRAmin <- CVstats(model=brt.HRAmin)
c8 <- plot_brt_perf(cv.HRAskew, pred.name = 'Min. home-range size')
cv.cHR <- CVstats(model=brt.cHR)
c9 <- plot_brt_perf(cv.cHR, pred.name = 'TAC breakpt')

ggarrange(c2, c3, c4, c5,
           ncol = 2, nrow = 2, align = "v")

' --- Plot PDP (Fig. 4 b,c,d & S3-9)--- '
gbm_obj <- brt.StepSkew[[2]]
# simple
# gbm.plot(gbm_obj, rug = T, plot.layout = c(3,3))

# against response and with bootstrapped CIs
source('R/plot_pdp_w_CIs.R')
#brt_boot = gbm.bootstrap.functions(gbm.object = gbm_obj, list.predictors = plot.gbm.4list(gbm_obj),n.divisions =100, n.reps = 500)
#save(brt_boot, file='R/bootstrapped_CIs/CIs_HRAmin.RData')
#load('R/bootstrapped_CIs/CIs_Stepskew.RData')
pred.name = expression(paste('Step-length skew', sep='')) # 'HR size (km'^'2',')'
g = 
  gbm.plot.boot.response_fullrug(gbm.object = gbm_obj, 
                                 plot_rug = F,
                                 booted.preds = brt_boot$function.preds,
                                 use_ggplot2 = TRUE,
                                 plots_on_one = TRUE,
                                 use_log_scale = c(2,3,4,5), # corresponds to plot order
                                 use_scientific_scale = 6,
                                 plots_exclude = NA,
                                 xlab_vec = gbm_obj[["contributions"]][["var"]])

' --- Plot interactions --- '
#### plot interactions (Fig. S7) ----
library(RColorBrewer)
gbm_obj <- brt.HRA[[3]]
gbm_int = gbm.interactions(gbm_obj)
gbm_int[["rank.list"]]

color1 <- colorRampPalette(brewer.pal(9,'YlGnBu'))(16)
p1 = plot.gbm(gbm_obj, i.var = c(7,5), at = seq(0,7.5, length.out = 16),
              col.regions = color1, colorkey=F, 
              xlab="Scent decay rate (10^)", ylab="Short-term memory decay rate (10^)")
p2 = plot.gbm(gbm_obj, i.var = c(7,4), at = seq(0,7.5, length.out = 16),
              col.regions = color1, colorkey=F,
              xlab="Scent decay rate (10^)", ylab="Long-term memory decay rate (10^)")
p3 = plot.gbm(gbm_obj, i.var = c(7,2), at = seq(0,7.5, length.out = 16),
              col.regions = color1, colorkey=F,
              xlab="Scent decay rate (10^)", ylab="Search speed")
p4 = plot.gbm(gbm_obj, i.var = c(9,7), at = seq(0,7.5, length.out = 16),
              col.regions = color1, colorkey=T,
              xlab="Scent deposition rate", ylab="Scent decay rate (10^)")


' --- Plot contributions (Fig. 4a)--- '
dat <- summariseBRT(list(brt.Steplength[[1]], brt.StepSkew[[2]], brt.DailyDist[[3]], brt.HRA[[3]], brt.HRAskew[[1]], brt.OverlapMax[[1]]), 
                    output.names = c('Step-length', 'Step-length skew', "Daily distance", "HR size", "HR size skew",
                                     'HR overlap'))

library(reshape)
dat.melt <- melt(dat)
dat.melt$X1 <- factor(dat.melt$X1, levels = c("ResourceRegenerationRate","ForagerSpeedSearch","ForagerSpeedFeeding",
                                              "LongDecayRate", "ShortDecayRate", "MemoryValueUninformed",
                                              "ScentDecayRate","ScentResponseSpatialScale","ScentDepositionRate"))
dat.melt$X2 <- factor(dat.melt$X2, levels = c('Step-length', 'Step-length skew', "Daily distance", "HR size", "HR size skew",
                                              'HR overlap'))

ggplot(dat.melt, aes(x=X2, y=X1, fill=value)) + 
  geom_tile(stat='identity') + xlab("") + ylab("") +
  geom_text(aes(label = round(value,1))) +
  scale_fill_gradient(low = "grey90", high = "forestgreen", name="Relative \n influence %") + 
  scale_y_discrete(limits = rev(levels(dat.melt$X1)), 
                   labels = c('Deposition rate', 'Response spatial scale', 'Decay rate', 'Memory expectation',
                              'Short decay rate', 'Long decay rate', 'Feed speed', 'Search speed', 'Regeneration rate')) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11)) +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))  +
  #labs(tag = "Scent") +
  theme(legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.tag.position = c(-0.02, 0.40), 
        plot.tag = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11)) 

