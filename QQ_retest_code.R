### Sort data for CalcSNR

# SCRIPT TO SORT DATA FOR META-ANALYSIS
# 1 SORT DATA
# 2 SAVE TO DATAFRAME FOR INPUT INTO CalculateSNRs_20.R script

####### SNR and CNR data for all databases 
rm(list=ls())
# load packages
library(plyr)
library(gridExtra)    # needed for graphs
library(grid)         # needed for graphs 
library(BayesFactor)
library(dplyr)
library(blme)
library(lme4)
library(readxl)       # read/write excel sheets
library(ggplot2)      # plots
library(ggrepel)
library(extrafont)    # change font in plots 
library(doMC) # reuiqred for multicore=TRUE in anovaBF
font_import()
loadfonts(device="win")
fonts

options(max.print=100000)

# set directories
data_dir <- '/Users/scotti/surfdrive/PhD/meta_analysis/Databases/PLoS One submission/POLS ONE submissions/submission_v4/QQ_retest'
save_dir <- '/Users/scotti/Downloads/QQ_retest'

# save things?
save_stuff <-''

# load in data
dat_orig <- read.csv(file.path(data_dir, 'orig/SNRinfo_all_original.csv')) # original measurements
dat_re <- read.csv(file.path(data_dir, 'retest/SNRinfo_QQretest.csv')) # retest measurements

# START SORTING STUFF
dat_orig$Label <- as.factor(dat_orig$Label)                             # make label a factor 
dat_orig <- rename(dat_orig, "Database_sub_ses" = "Database")           # rename column
dat_orig$Database_sub_ses <- as.character(dat_orig$Database_sub_ses)    # make characters

dat_re$Label <- as.factor(dat_re$Label)                             # make label a factor 
dat_re <- rename(dat_re, "Database_sub_ses" = "Database")           # rename column
dat_re$Database_sub_ses <- as.character(dat_re$Database_sub_ses)    # make characters

# add mean of caudate nuclei and standard deviation of this distibution 
mean_caud_orig <- read.csv('/Users/scotti/Downloads/QQ_retest/orig/mean_caud_all.csv')
sd_caud_orig <- read.csv('/Users/scotti/Downloads/QQ_retest/orig/mean_caud_sd.csv')

mean_caud_re <- read.csv('/Users/scotti/Downloads/QQ_retest/retest/mean_caud_all.csv')
sd_caud_re <- read.csv('/Users/scotti/Downloads/QQ_retest/retest/mean_caud_sd.csv')

# Add to dataframe
dat_orig$mean_CAUD_raw <- mean_caud_orig$X0
dat_orig$mean_CAUD_raw_sd <- sd_caud_orig$X0
dat_orig$mean_CAUD_both <- dat_orig$mean_CAUD_raw/dat_orig$mean_CAUD_raw_sd

dat_re$mean_CAUD_raw <- mean_caud_re$X0
dat_re$mean_CAUD_raw_sd <- sd_caud_re$X0
dat_re$mean_CAUD_both <- dat_re$mean_CAUD_raw/dat_re$mean_CAUD_raw_sd

# Define contrasts and name
dat_orig$contrast <- 0           # pre-allocate
dat_orig$name <- 0               # pre-allocate
for (i in 1:nrow(dat_orig)){
  dat_orig$name[i] <- strsplit(dat_orig$Database_sub_ses, '_',fixed=TRUE)[[i]][1]          # Define Databases per row 
  dat_orig$contrast[i] <- strsplit(dat_orig$Database_sub_ses, '_',fixed=TRUE)[[i]][2]      # Define image contrst for each observation
}

dat_re$contrast <- 0           # pre-allocate
dat_re$name <- 0               # pre-allocate
for (i in 1:nrow(dat_re)){
  dat_re$name[i] <- strsplit(dat_re$Database_sub_ses, '_',fixed=TRUE)[[i]][1]          # Define Databases per row 
  dat_re$contrast[i] <- strsplit(dat_re$Database_sub_ses, '_',fixed=TRUE)[[i]][2]      # Define image contrst for each observation
}

dat_orig$contrast[dat_orig$name=='ATAG'&dat_orig$contrast=='T1'] <- 'qT1'        # T1 image = quantitative T1 image
dat_orig$contrast[dat_orig$name=='ATAG'&dat_orig$contrast=='UNI'] <- 'T1'        # UNI images = T1w images

dat_orig$new_CNR_calc <- 0
for (i in 1:nrow(dat_orig)){
  dat_orig$new_CNR_calc[i] <- (dat_orig$CC_meanvox[i] - dat_orig$mean_CAUD_raw[i])/(sqrt(dat_orig$CC_sd[i]^2+dat_orig$mean_CAUD_raw_sd[i]^2))
}

dat_re$new_CNR_calc <- 0
for (i in 1:nrow(dat_re)){
  dat_re$new_CNR_calc[i] <- (dat_re$CC_meanvox[i] - dat_re$mean_CAUD_raw[i])/(sqrt(dat_re$CC_sd[i]^2+dat_re$mean_CAUD_raw_sd[i]^2))
}


#dat$CNR_old <- dat$new_CNR_calc

# remove columns that we won't use so it looks nicer
dat_orig <- dat_orig %>% dplyr::select(Database_sub_ses, Label, CC, CC_sd, CAUDL, CAUDL_sd, CAUDR, CAUDR_sd, 
                             mean_CAUD_raw,mean_CAUD_raw_sd,mean_CAUD_both, contrast, name, new_CNR_calc)

dat_re <- dat_re %>% dplyr::select(Database_sub_ses, Label, CC, CC_sd, CAUDL, CAUDL_sd, CAUDR, CAUDR_sd, 
                                       mean_CAUD_raw,mean_CAUD_raw_sd,mean_CAUD_both, contrast, name, new_CNR_calc)

# Get rid of data we are not comparing
dat_orig <- dat_orig[dat_orig$contrast=='T1' & !dat_orig$name %in% c('250','MAASTRICHT','MASSIVE') & !dat_orig$Label %in% c('HCPT1PRO','AHEADT1YOUSLAB','AHEADT1MIDSLAB','AHEADT1ELDSLAB','ATAGUNIELDSLAB','ATAGUNIMIDSLAB','ATAGUNIYOUSLAB'),]

# Add columns to separate, then concatenate dataframes
dat_orig$sample <- 1
dat_re$sample <- 2

dat <- rbind(dat_orig, dat_re)

# add voxel sizes for normalization
dat$dim[dat$name=='IXI'&dat$contrast=='T1'] <- 0.94*0.94*1.2
dat$dim[dat$name=='ATAG'&dat$contrast=='T1'] <- 0.7*0.7*0.7
dat$dim[dat$name=='KIR21'&dat$contrast=='T1'] <- 1*1*1.2
dat$dim[dat$name=='AGEILITY'&dat$contrast=='T1'] <- 1*1*1
dat$dim[dat$name=='GSP'&dat$contrast=='T1'] <- 1.2*1.2*1.2
dat$dim[dat$name=='MPICBS'&dat$contrast=='T1'] <- 0.5*0.5*0.5
dat$dim[dat$name=='SALD'&dat$contrast=='T1']<-1
dat$dim[dat$name=='MPILMBB'&dat$contrast=='T1']<-1
dat$dim[dat$name=='MSC'&dat$contrast=='T1']<-0.8*0.8*0.8
dat$dim[dat$name=='Forrest'&dat$contrast=='T1']<- 0.67*0.67*0.7
dat$dim[dat$name=='CAMCAN'&dat$contrast=='T1']<-1
dat$dim[dat$name=='NKIRS'&dat$contrast=='T1']<-1
dat$dim[dat$name=='HCPYA'&dat$contrast=='T1']<-0.7*0.7*0.7
dat$dim[dat$name=='AHEAD'&dat$contrast=='T1']<-0.64*0.64*0.7
dat$dim[dat$name=='DLBS'] <- 1
dat$dim[dat$name=='RAIDERS'] <- 0.9375*0.9375*1
dat$dim[dat$name=='PTBP'] <- 1

# Normalise the values by voxel size
dat$normCNR <- dat$new_CNR_calc/dat$dim 
dat$normSNR_CC <- dat$CC/dat$dim
dat$normSNR_CN <- dat$mean_CAUD_both/dat$dim
dat$normSNR_LCN <- dat$CAUDL/dat$dim
dat$normSNR_RCN <- dat$CAUDR/dat$dim

# FINAL changes
dat$contrast <- as.factor(dat$contrast)
dat$name[dat$name=='Forrest'] <- 'FORREST'
dat$sample <- as.factor(dat$sample)
dat$name <- as.factor(dat$name)
dat$normSNR_CC <- as.numeric(dat$normSNR_CC)

# select columns
dat <- dat %>% dplyr::select(Database_sub_ses, Label, normSNR_CC, normSNR_CN, normSNR_LCN, normSNR_RCN, normCNR, contrast, name, dim, CC, sample)

####### DO THE ANALYSIS #######

deets <- aggregate(normSNR_CC ~ name + sample, dat, function(x) cbind(mean(x),sd(x),length(x))) 

##### PLOT new SNR vs old SNR #####

theSNR <- aggregate(normSNR_CC ~ name + sample, dat, mean)
theSNRsd <- aggregate(normSNR_CC ~ name + sample, dat, sd)
theSNRlen <- aggregate(normSNR_CC ~ name + sample, dat, length)
theSNR$sem <- theSNRsd[,3]/sqrt(theSNRlen[,3])
secondSNR <- theSNR[18:34,]
colnames(secondSNR) <- c('name','sample','normSNR_CC2','sem2')
theSNR <- cbind(theSNR[1:17,],secondSNR[,3:4])

if (save_stuff == 1){tiff(file.path(save_dir,"old_vs_new_SNRCC.tiff"), units="in", width=5,height=5,res=300)}

p <- ggplot(theSNR, aes(y=normSNR_CC, x=normSNR_CC2))+
  geom_smooth(method='lm',formula=y~x, se=FALSE,col="black",linetype="dashed", size=0.4) +
  geom_point(size=1.2) +
  geom_errorbar(aes(ymin=normSNR_CC-sem, ymax=normSNR_CC+sem), width=.2,size=0.3) +
  geom_errorbarh(aes(xmin=normSNR_CC2-sem2, xmax=normSNR_CC2+sem2),size=0.3) +
  scale_x_continuous(name="SNR 2", breaks=seq(0,320,50)) +
  scale_y_continuous(name="SNR 1", breaks=seq(0,320,50)) +
  ggtitle("Scan time vs acquired SNR") +

set.seed(1337)
p <- p + geom_label_repel(aes(label = name),
                          #box.padding   = 0.35,
                          point.padding = 0.48,
                          label.padding = 0.1,
                          max.iter =93,
                          size=1.2,
                          fill='white',
                          label.size=0,
                          segment.color = 'grey10',
                          label.r = 0.05,
                          segment.size = 0.25)

print(p)

dev.off()

plot(deets$normSNR_CC[1:17,1],deets$normSNR_CC[18:34,1])
abline(1,1)

retest_model <- lmer(normSNR_CC[,1] ~ (1 | name) + sample, data=deets)

dat$newname <- c(rep(1,15),rep(2,15),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(7,15),rep(8,15),
                 rep(9,5),rep(10,5),rep(11,15),rep(12,15),rep(13,5),rep(14,15),rep(15,5),rep(16,15),
                 rep(17,5),rep(1,15),rep(2,14),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(7,15),rep(8,15),
                 rep(9,5),rep(10,5),rep(11,15),rep(12,15),rep(13,5),rep(14,15),rep(15,5),rep(16,15),
                 rep(17,5))
dat$newname <- as.factor(dat$newname)

##
bfaovCC <- anovaBF(normSNR_CC ~ sample + newname, data=dat, whichRandom="newname",iterations = 1000000, progress=TRUE, whichModels = 'all')
bfaovCC
bfaovCN <- anovaBF(normSNR_CN ~ sample + newname, data=dat, whichRandom="newname",iterations = 1000000, progress=TRUE, whichModels = 'all')
bfaovCN
bfaovCNR <- anovaBF(normCNR ~ sample + newname, data=dat, whichRandom="newname",iterations = 1000000, progress=TRUE, whichModels = 'all')
bfaovCNR

## Check again with linear model
lmaovCC <- lmBF(normSNR_CC ~ sample + newname, data=dat, whichRandom="newname")
lmaovCC

