# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Theoretical examples
# ====================================
library(reshape)
library(ggpubr)
source('R/0_functions.R')

' --- Experiment 1: scenarios --- '
# convert landscape matrix to raster
envm<-t(as.matrix(env))
r <- raster(envm[nrow(envm):1,], xmn=0, xmx=nrow(envm), ymn=0, ymx=ncol(envm)) # for mean resource calculation

# run analysis (HR size, mean resource within HR)
Expt1_S3_5 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S3_5", rep_index = 0:99, r = r, D = 5)
Expt1_S3_10 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S3_10", rep_index = 0:99, r = r, D = 10)
Expt1_S3_15 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S3_15", rep_index = 0:99, r = r, D = 15)
Expt1_S2_5 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S2_5", rep_index = 0:99, r = r, D = 5)
Expt1_S2_10 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S2_10", rep_index = 0:99, r = r, D = 10)
Expt1_S2_15 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S2_15", rep_index = 0:99, r = r, D = 15)
Expt1_S1_5 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S1_5", rep_index = 0:99, r = r, D = 5)
Expt1_S1_10 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S1_10", rep_index = 0:99, r = r, D = 10)
Expt1_S1_15 <- GenOutput_Ex1(tracksfile="TheoreticalEg/Expt1_S1_15", rep_index = 0:99, r = r, D = 15)

# assemble and save data
Expt1_S1 <- rbind(Expt1_S1_5,Expt1_S1_10,Expt1_S1_15)
Expt1_S2 <- rbind(Expt1_S2_5,Expt1_S2_10,Expt1_S2_15)
Expt1_S3 <- rbind(Expt1_S3_5,Expt1_S3_10,Expt1_S3_15)
Expt1 <- rbind(Expt1_S1, Expt1_S2, Expt1_S3)
#save(Expt1_S1,Expt1_S2,Expt1_S3,Expt1, file = 'TheoreticalEg/Expt1.RData')
#load('TheoreticalEg/Expt1.RData')

# Pearson's R
cor.test(log10(Expt1_S1[,2]), log10(Expt1_S1[,4]), method=c("pearson", "kendall", "spearman"))

# Two-way ANOVA
aov1 <- aov(log10(HRA95) ~ as.factor(Density) + as.factor(Scenario), data = Expt1)
summary(aov1)
plot(aov1) # check residuals
aggregate(log10(HRA95) ~ as.factor(Density) + as.factor(Scenario), # calculate means
          data = Expt1, FUN = mean)

' --- Experiment 2: removal --- '
#load('TheoreticalEg/Expt2.RData')
H0 <- GenOutput_Ex2(tracksfile = 'TheoreticalEg/Expt2_H0', file_index = 0:99, 
                    phase1 = c(5001,20000), phase2 = c(35001,50000), H1 = FALSE)
H1 <- GenOutput_Ex2(tracksfile = 'TheoreticalEg/Expt2_H1', file_index = 0:99, 
                      phase1 = c(5001,20000), phase2 = c(35001,50000), H1 = TRUE)

#H1 <- Map(c, H1_1,list(na.omit(H1_2[[1]]),na.omit(H1_2[[2]])))
#H0 <- Map(c, H0_1,list(na.omit(H0_2[[1]]),na.omit(H0_2[[2]])))

# Paired t-test
t.test(H0[[1]],H0[[2]],paired=TRUE, alternative = 'less')
t.test(H1[[1]],H1[[2]],paired=TRUE, alternative = 'less')

# check residuals
d <- H0[[1]] - H0[[2]]
# d <- H1[[1]] - H1[[2]]
shapiro.test(d)
plot(density(na.omit(d)))
qqnorm(d)
qqline(d, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75))

' --- Fig. 2 --- '
# Fig 2a
Expt1$colour <- paste(Expt1$Scenario,Expt1$Density, sep='')

ggviolin(Expt1, x = "Density", y = 'HRA95', fill = 'colour', facet.by = 'Scenario',
         panel.labs = list(Scenario = c('Territoriality','Memory','Territoriality + Memory')),
         size = 0.4, add = "boxplot", width = 0.5, add.params = list(size=0.5)) + 
  ggtitle('') + xlab('Density') + ylab(expression(paste("HR size (unit"^2,')', sep=''))) + 
  scale_y_log10() +
  scale_fill_manual(values = c("indianred", "indianred4","lightcoral",
                               "skyblue3", "skyblue4","skyblue1",
                               "seagreen3", "seagreen4", "seagreen2")) + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size=12, face = 'bold')) +
  font("xy.text", size = 9)

# Fig 2b
ggplot(Expt1, aes(x = MeanResource, y = HRA95, color = colour)) +
  geom_point(size = 1) +
  ggtitle('') + xlab('Mean resource') + ylab(expression(paste("HR size (unit"^2,')', sep=''))) + 
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(name = "Legend",
                     values = c("15" = "lightcoral", "110" = "indianred","115"  = "indianred4",
                                "25" = "skyblue1", "210" = "skyblue3","215"  = "skyblue4",
                                "35" = "seagreen2", "310" = "seagreen3","315"  = "seagreen4"),
                     breaks = c("15","110","115","25","210","215","35","310","315"),
                     labels = c("T (n = 5)", "T (n = 10)", "T (n = 15)",
                                "M (n = 5)", "M (n = 10)", "M (n = 15)",
                                "T + M (n = 5)", "T + M (n = 10)", "T + M (n = 15)")) + 
  theme_bw() + 
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_text(size=12, face = 'bold'),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  font("xy.text", size = 9)

# Individual plots
par(mar=c(1, 1, 1, 1), mfrow=c(1,3),
    oma = c(4, 4, 2, 2))

colors <- c("lightcoral","indianred", "indianred4")
plot(log10(Expt1_S1[,4]),log10(Expt1_S1[,2]),xlab="",ylab="",
     col = colors, pch=shapes, frame = FALSE) #, xlim = c(-4.5,-3.7), ylim = c(2.3,4.5))
box(bty="l")
#legend(-4.6,3.0, legend='r = -0.44, \n p < 0.001', bty = "n")
mtext("log10(HR size)", side = 2, line = 2.5)

colors <- c("skyblue1","skyblue3", "skyblue4")
plot(log10(Expt1_S2[,4]),log10(Expt1_S2[,2]),xlab="",ylab="",
     col =colors, pch=shapes, frame = FALSE) #, xlim = c(-4.5,-3.7), ylim = c(2.3,4.5))
box(bty="l")
#legend(-4.0,3.70, legend='r = -0.75, \n p < 0.001', bty = "n")
mtext("log10(Mean resource)", side = 1, line = 3)

colors <- c("seagreen2", "seagreen3", "seagreen4")
plot(log10(Expt1_S3[,4]),log10(Expt1_S3[,2]),xlab="",ylab="",
     col = colors, pch=shapes, frame = FALSE) #, xlim = c(-4.5,-3.7), ylim = c(2.3,4.5))
#legend(-4.17,2.85, legend='r = -0.75, \n p < 0.001', bty = "n")
box(bty="l")

# 2d: trajectories (one rep)
par(mar=c(1, 1, 1, 1), mfrow=c(1,3),
    oma = c(4, 4, 2, 2))

i = sample(0:99,1)
df <- read.csv(paste('TheoreticalEg/Expt1_S3_10/Tracks',i,'.csv',sep=''), header = T)
image(1:nrow(env), 1:ncol(env), as.matrix(env), col = colfunc(10), axes = F, xlab = "", ylab = "")
plotTraject(df, colfunc = rainbow, num_steps = 50000)
box(lty=1, col='black')

' --- Fig. 3 --- '
# Fig. 3a
Expt2_melt <- melt(list(H0,H1))
colnames(Expt2_melt) <- c('AreaOverlap','Phase','Treatment')

ggviolin(Expt2_melt, x = "Phase", y = "AreaOverlap", fill = 'grey80', facet.by = 'Treatment',
         panel.labs = list(Treatment = c('Control','Removal')),
         size = 0.4, add = "boxplot", width = 0.5, add.params = list(size=0.5)) + 
  ggtitle('') + xlab('Phase') + ylab(expression(paste("Aggregate HR overlap (unit"^2,')', sep=''))) + #ylim(-0.6,1.2) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  font("xy.text", size = 9)

# Fig. 3b: plot one rep
df <- read.csv('TheoreticalEg/Expt2_H0/Tracks0.csv')

# arange data into phases
survivors <- sample(unique(df$id), 6, replace=F) # H0
survivors <- # for H1
  df %>% 
  group_by(id) %>% 
  filter(row_number() == 50000) %>% 
  pull(id)

before.alive <-
  df %>% 
  filter(id %in% survivors) %>% 
  group_by(id) %>% 
  filter(row_number() %in% 5001:20000)
before.dead <-
  df %>% 
  filter(!id %in% survivors) %>% 
  group_by(id) %>% 
  filter(row_number() %in% 5001:20000)

after.alive <-
  df %>% 
  filter(id %in% survivors) %>% 
  group_by(id) %>% 
  filter(row_number() %in% 35001:50000)
after.dead <- # only for control
  df %>% 
  filter(!id %in% survivors) %>% 
  group_by(id) %>% 
  filter(row_number() %in% 35001:50000)

par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
image(1:nrow(env), 1:ncol(env), as.matrix(env), col = colfunc(10), axes=F, ann = F)
box(lty=1, col='black')
plotTraject(before.alive, colfunc = colorRampPalette(c("turquoise1","blue3")), num_steps = 25000)
plotTraject(before.dead, colfunc = colorRampPalette(c("orange1","red3")), num_steps = 25000)

plotTraject(after.alive, colfunc = colorRampPalette(c("turquoise1","blue3")), num_steps = 25000)
plotTraject(after.dead, colfunc = colorRampPalette(c("orange1","red3")), num_steps = 25000)

' --- Supplemental figures --- '
# 2a: HR size distributions
par(mar=c(1, 1, 1, 1), mfrow=c(1,3),
    oma = c(4, 4, 2, 2))

hist(Expt1_S1_10$HRA95,
     main = '',col = rgb(0,0,0,0), ylim = c(0,6), 
     breaks =seq(0,ceiling(max(Expt1_S1_10$HRA95)/1000)*1000,1000), border=rgb(0,0,0,0.0), 
     xlab = '', ylab = 'Frequency')
box(bty="l")
for (i in 0:99){
  hist(Expt1_S1_10$HRA95[Expt1_S1_10$replicate == i], #breaks = seq(0,180,5),
       main = '', xlab='',ylab=' ',
       col = rgb(0,0,0,0.05), breaks =seq(0,ceiling(max(Expt1_S1_10$HRA95)/1000)*1000,1000), border=rgb(0,0,0,0.0), add=T)
}
mtext("Frequency", side = 2, line = 2.5)

hist(Expt1_S2_10$HRA95,
     main = '',col = rgb(0,0,0,0), ylim = c(0,6), 
     breaks =seq(0,ceiling(max(Expt1_S2_10$HRA95)/1000)*1000,1000), border=rgb(0,0,0,0.0), 
     xlab = '', ylab = '')
box(bty="l")
for (i in 0:99){
  hist(Expt1_S2_10$HRA95[Expt1_S2_10$replicate == i], #breaks = seq(0,180,5),
       main = '', xlab='',ylab=' ',
       col = rgb(0,0,0,0.05), breaks =seq(0,ceiling(max(Expt1_S2_10$HRA95)/1000)*1000,1000), border=rgb(0,0,0,0.0), add=T)
}
mtext(expression(paste("Home-range size (unit"^2,')', sep='')), side = 1, line = 3)

hist(Expt1_S3_10$HRA95,
     main = '',col = rgb(0,0,0,0), ylim = c(0,6),
     breaks =seq(0,ceiling(max(Expt1_S3_10$HRA95)/1000)*1000,1000), border=rgb(0,0,0,0.0), 
     xlab='',ylab='',)

for (i in 0:99){
  hist(Expt1_S3_10$HRA95[Expt1_S3_10$replicate == i], #breaks = seq(0,180,5),
       main = '', xlab='',ylab=' ',
       col = rgb(0,0,0,0.05), breaks =seq(0,ceiling(max(Expt1_S3_10$HRA95)/1000)*1000,1000), border=rgb(0,0,0,0.0), add=T)
}
