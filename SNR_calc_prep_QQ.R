### Sort data for CalcSNR

# SCRIPT TO SORT DATA FOR META-ANALYSIS
# 1 SORT DATA
# 2 INPUT INTO CalculateSNRs.R script

####### SNR and CNR data for all databases 
rm(list=ls())
# load packages
library(plyr)
library(gridExtra)    # needed for graphs
library(grid)         # needed for graphs 
library(BayesFactor)
library(dplyr)
library(blme)
library(readxl)       # read/write excel sheets
library(ggplot2)      # plots
library(ggrepel)
library(extrafont)    # change font in plots 
font_import()
loadfonts(device="win")
fonts

options(max.print=100000)

# set directories
data_dir <- '/Users/scotti/surfdrive/PhD/meta_analysis/Databases/R_stuff'
save_dir <- '/Users/scotti/surfdrive/PhD/meta_analysis/Databases/R_stuff'

# load in data
dat <- read.csv(file.path(data_dir, 'SNRinfo_all_20.csv')) # original measurements
#dat <- read.csv(file.path(data_dir, 'SNRinfo_all_equalvol_20.csv')) # measurements adjusted for volume

# START SORTING STUFF
dat$Label <- as.factor(dat$Label)                             # make label a factor 
dat <- rename(dat, "Database_sub_ses" = "Database")           # rename column
dat$Database_sub_ses <- as.character(dat$Database_sub_ses)    # make characters

# add mean of caudate nuclei and standard deviation of this distibution 
mean_caud <- read.csv('/Users/scotti/Downloads/QQ_retest/orig/mean_caud_all.csv')
sd_caud <- read.csv('/Users/scotti/Downloads/QQ_retest/orig/mean_caud_sd.csv')

dat$mean_CAUD_raw <- mean_caud$X0
dat$mean_CAUD_raw_sd <- sd_caud$X0
dat$mean_CAUD_both <- dat$mean_CAUD_raw/dat$mean_CAUD_raw_sd

dat$contrast <- 0           # pre-allocate
dat$name <- 0               # pre-allocate
for (i in 1:nrow(dat)){
  dat$name[i] <- strsplit(dat$Database_sub_ses, '_',fixed=TRUE)[[i]][1]          # Define Databases per row 
  dat$contrast[i] <- strsplit(dat$Database_sub_ses, '_',fixed=TRUE)[[i]][2]      # Define image contrst for each observation
}

#dat$name[dat$name=='STW'] <- 'AHEAD'                              # Keep al names consistent STW = AHEAD
dat$contrast[dat$name=='ATAG'&dat$contrast=='T1'] <- 'qT1'        # T1 image = quantitative T1 image
dat$contrast[dat$name=='ATAG'&dat$contrast=='UNI'] <- 'T1'        # UNI images = T1w images

dat$new_CNR_calc <- 0
for (i in 1:nrow(dat)){
  dat$new_CNR_calc[i] <- (dat$CC_meanvox[i] - dat$mean_CAUD_raw[i])/(sqrt(dat$CC_sd[i]^2+dat$mean_CAUD_raw_sd[i]^2))
}

# remove columns that we won't use so it looks nicer
dat <- dat %>% dplyr::select(Database_sub_ses, Label, CC, CC_sd, CAUDL, CAUDL_sd, CAUDR, CAUDR_sd, 
                      mean_CAUD_raw,mean_CAUD_raw_sd,mean_CAUD_both, contrast, name, new_CNR_calc)

# CC = SNR of CC
# meanCAUD = SNR of CNs
# CNR_old = CNR of image, not yet normalized by voxel dimensions

#####################################
##################################### MANUALLY ADD VOXEL DIMENSIONS

# IXI
dat$dim <- 0
dat$dim[dat$name=='IXI'&dat$contrast=='T1'] <- 0.94*0.94*1.2
dat$dim[dat$name=='IXI'&dat$contrast=='T2'] <- 0.9*0.9*1.2
dat$dim[dat$name=='IXI'&dat$contrast=='PD'] <- 0.9*0.9*1.2 
# ATAG
dat$dim[dat$name=='ATAG'&dat$contrast=='T1'] <- 0.7*0.7*0.7
dat$dim[dat$name=='ATAG'&dat$contrast=='qT1'] <- 0.7*0.7*0.7
dat$dim[dat$name=='ATAG'&dat$contrast=='FLASH'] <- 0.5*0.5*0.5
dat$dim[dat$name=='ATAG'&dat$contrast=='T1'&dat$Label%in%c('ATAGUNIELDSLAB','ATAGUNIMIDSLAB','ATAGUNIYOUSLAB')] <- 0.6*0.6*0.6
dat$dim[dat$name=='ATAG'&dat$contrast=='qT1'&dat$Label%in%c('ATAGT1ELDSLAB','ATAGT1MIDSLAB','ATAGT1YOUSLAB')] <- 0.6*0.6*0.6
# Kirby 21
dat$dim[dat$name=='KIR21'&dat$contrast=='T1'] <- 1*1*1.2
dat$dim[dat$name=='KIR21'&dat$contrast=='FLAIR'] <- 1.1*1.1*1.1
# Age-ility
dat$dim[dat$name=='AGEILITY'&dat$contrast=='T1'] <- 1*1*1
# GSP
dat$dim[dat$name=='GSP'&dat$contrast=='T1'] <- 1.2*1.2*1.2
# MPI-CBS
dat$dim[dat$name=='MPICBS'&dat$contrast=='T1'] <- 0.5*0.5*0.5
dat$dim[dat$name=='MPICBS'&dat$contrast=='qT1']<- 0.5*0.5*0.5
dat$dim[dat$name=='MPICBS'&dat$contrast=='T2']<- 0.5*0.5*0.5
dat$dim[dat$name=='MPICBS'&dat$contrast=='qT2star']<- 0.5*0.5*0.5
# MASSIVE
dat$dim[dat$name=='MASSIVE'&dat$contrast=='FLAIR']<-1
dat$dim[dat$name=='MASSIVE'&dat$contrast=='T1']<-1
dat$dim[dat$name=='MASSIVE'&dat$contrast=='T2']<-1
# SALD
dat$dim[dat$name=='SALD'&dat$contrast=='T1']<-1
# MPI-LMBB
dat$dim[dat$name=='MPILMBB'&dat$contrast=='T1']<-1
dat$dim[dat$name=='MPILMBB'&dat$contrast=='FLAIR']<-0.49*0.49*1
dat$dim[dat$name=='MPILMBB'&dat$contrast=='qT1']<-1
# Midnight scan club
dat$dim[dat$name=='MSC'&dat$contrast=='T1']<-0.8*0.8*0.8
dat$dim[dat$name=='MSC'&dat$contrast=='T2']<-0.8*0.8*0.8
# StudyForrest
dat$dim[dat$name=='Forrest'&dat$contrast=='T1']<- 0.67*0.67*0.7
dat$dim[dat$name=='Forrest'&dat$contrast=='T2']<- 0.67*0.67*0.7
dat$dim[dat$name=='Forrest'&dat$contrast=='SWI']<- 0.67*0.67*0.7
# Cam-Can
dat$dim[dat$name=='CAMCAN'&dat$contrast=='T1']<-1
dat$dim[dat$name=='CAMCAN'&dat$contrast=='T2']<-1
# NKIRS
dat$dim[dat$name=='NKIRS'&dat$contrast=='T1']<-1
# HCP-YA
dat$dim[dat$name=='HCPYA'&dat$contrast=='T1']<-0.7*0.7*0.7
dat$dim[dat$name=='HCPYA'&dat$contrast=='T2']<-0.7*0.7*0.7
# AHEAD
dat$dim[dat$name=='AHEAD'&dat$contrast=='T1']<-0.64*0.64*0.7
dat$dim[dat$name=='AHEAD'&dat$contrast=='QSM']<-0.64*0.64*0.7
dat$dim[dat$name=='AHEAD'&dat$contrast=='PDmap']<-0.64*0.64*0.7
dat$dim[dat$name=='AHEAD'&dat$contrast=='qt2star']<-0.64*0.64*0.7
dat$dim[dat$name=='AHEAD'&dat$contrast=='qT1']<-0.64*0.64*0.7
dat$dim[501:560] <- 0.64*0.64*0.7
dat$dim[dat$name=='AHEAD'&dat$contrast=='T1'&dat$Label%in%c('AHEADT1YOUSLAB','AHEADT1MIDSLAB','AHEADT1ELDSLAB')]<-0.5*0.5*0.5
dat$dim[dat$name=='AHEAD'&dat$contrast=='PDmap'&dat$Label%in%c('AHEADPDYOUSLAB','AHEADPDMIDSLAB','AHEADPDELDSLAB')]<-0.5*0.5*0.5
dat$dim[dat$name=='AHEAD'&dat$contrast=='qt2star'&dat$Label%in%c('AHEADQT2STARYOUSLAB','AHEADQT2STARMIDSLAB','AHEADQT2STARELDSLAB')]<-0.5*0.5*0.5
dat$dim[dat$name=='AHEAD'&dat$contrast=='qT1'&dat$Label%in%c('AHEADQT1YOUSLAB','AHEADQT1MIDSLAB','AHEADQT1ELDSLAB')]<-0.5*0.5*0.5
dat$dim[561:620] <- 0.5*0.5*0.5
# DLBS
dat$dim[dat$name=='DLBS'] <- 1
# MAASTRICHT
dat$dim[dat$name=='MAASTRICHT'&dat$contrast=='T1'] <- 0.7*0.7*0.7 
dat$dim[dat$name=='MAASTRICHT'&dat$contrast=='PD'] <- 0.7*0.7*0.7
dat$dim[dat$name=='MAASTRICHT'&dat$contrast=='t2star'] <- 0.7*0.7*0.7
# 250
dat$dim[dat$name=='250'] <- 0.25*0.25*0.25
# RAIDERS
dat$dim[dat$name=='RAIDERS'] <- 0.9375*0.9375*1
# PTBP
dat$dim[dat$name=='PTBP'] <- 1

#####################################
##################################### MANUALLY ADD SCAN TIME

# IXI
dat$time <- 0
dat$time[dat$name=='IXI'&dat$contrast=='T1'] <- NA
dat$time[dat$name=='IXI'&dat$contrast=='T2'] <- NA
dat$time[dat$name=='IXI'&dat$contrast=='PD'] <- NA
# ATAG
dat$time[dat$name=='ATAG'&dat$contrast=='T1'] <- 657
dat$time[dat$name=='ATAG'&dat$contrast=='qT1'] <- 657
dat$time[dat$name=='ATAG'&dat$contrast=='FLASH'] <- 1038
dat$time[dat$name=='ATAG'&dat$contrast=='T1'&dat$Label%in%c('ATAGUNIELDSLAB','ATAGUNIMIDSLAB','ATAGUNIYOUSLAB')] <- 547
dat$time[dat$name=='ATAG'&dat$contrast=='qT1'&dat$Label%in%c('ATAGT1ELDSLAB','ATAGT1MIDSLAB','ATAGT1YOUSLAB')] <- 547
# Kirby 21
dat$time[dat$name=='KIR21'&dat$contrast=='T1'] <- 356
dat$time[dat$name=='KIR21'&dat$contrast=='FLAIR'] <- 528
# Age-ility
dat$time[dat$name=='AGEILITY'&dat$contrast=='T1'] <- 405
# GSP
dat$time[dat$name=='GSP'&dat$contrast=='T1'] <- 132
# MPI-CBS
dat$time[dat$name=='MPICBS'&dat$contrast=='T1'] <- 1682 # maybe 3364
dat$time[dat$name=='MPICBS'&dat$contrast=='qT1']<- 1682 #maybe 3364
dat$time[dat$name=='MPICBS'&dat$contrast=='T2']<- 1572
dat$time[dat$name=='MPICBS'&dat$contrast=='qT2star']<- 1572
# MPI-CBS
dat$time[dat$name=='MASSIVE'&dat$contrast=='FLAIR']<-225
dat$time[dat$name=='MASSIVE'&dat$contrast=='T1']<-225
dat$time[dat$name=='MASSIVE'&dat$contrast=='T2']<-165
# SALD
dat$time[dat$name=='SALD'&dat$contrast=='T1']<-266
# MPI-LMBB
dat$time[dat$name=='MPILMBB'&dat$contrast=='T1']<-502
dat$time[dat$name=='MPILMBB'&dat$contrast=='FLAIR']<-NA
dat$time[dat$name=='MPILMBB'&dat$contrast=='qT1']<-502
# midnight scan club
dat$time[dat$name=='MSC'&dat$contrast=='T1']<-537.6
dat$time[dat$name=='MSC'&dat$contrast=='T2']<-NA
# StudyForrest
dat$time[dat$name=='Forrest'&dat$contrast=='T1']<- 769
dat$time[dat$name=='Forrest'&dat$contrast=='T2']<- 460
dat$time[dat$name=='Forrest'&dat$contrast=='SWI']<- NA
# Cam-Can
dat$time[dat$name=='CAMCAN'&dat$contrast=='T1']<-252
dat$time[dat$name=='CAMCAN'&dat$contrast=='T2']<-250
# NKIRS
dat$time[dat$name=='NKIRS'&dat$contrast=='T1']<-258
# HCP YA
dat$time[dat$name=='HCPYA'&dat$contrast=='T1']<-460
dat$time[dat$name=='HCPYA'&dat$contrast=='T2']<-504
# AHEAD
dat$time[dat$name=='AHEAD'&dat$contrast=='T1']<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='QSM']<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='PDmap']<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='qt2star']<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='qT1']<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='T1'&dat$Label%in%c('AHEADT1YOUSLAB','AHEADT1MIDSLAB','AHEADT1ELDSLAB')]<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='PDmap'&dat$Label%in%c('AHEADPDYOUSLAB','AHEADPDMIDSLAB','AHEADPDELDSLAB')]<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='qt2star'&dat$Label%in%c('AHEADQT2STARYOUSLAB','AHEADQT2STARMIDSLAB','AHEADQT2STARELDSLAB')]<-1193
dat$time[dat$name=='AHEAD'&dat$contrast=='qT1'&dat$Label%in%c('AHEADQT1YOUSLAB','AHEADQT1MIDSLAB','AHEADQT1ELDSLAB')]<-1193
dat$time[dat$name=='AHEAD'] <- 1193
# DLBS
dat$time[dat$name=='DLBS'] <- 237
# MAASTRICHT
dat$time[dat$name=='MAASTRICHT'&dat$contrast=='T1'] <- 793.6
dat$time[dat$name=='MAASTRICHT'&dat$contrast=='PD'] <- 353.28
dat$time[dat$name=='MAASTRICHT'&dat$contrast=='t2star'] <- 1280 
# 250
dat$time[dat$name=='250'] <- 3180
# RAIDERS
dat$time[dat$name=='RAIDERS'] <- NA
# PTBP
dat$time[dat$name=='PTBP'] <- 488

### sort naming issue
dat <- transform(dat, Label=revalue(Label,c('SLADT1YOU'='SALDT1YOU'))) # cannot just simply rename as it is a factor

############################################################
########## NORMALIZE DATA BY VOXEL DIMENSIONS ##########
############################################################

## Replace original values with normalized values
dat$normCNR <- dat$new_CNR_calc/dat$dim 
dat$normSNR_CC <- dat$CC/dat$dim
dat$normSNR_CN <- dat$mean_CAUD_both/dat$dim
dat$normSNR_LCN <- dat$CAUDL/dat$dim
dat$normSNR_RCN <- dat$CAUDR/dat$dim

#There are now 2 CNR columncs. CNR and CCNR_ol, both coumns are the same 

############################################################
############################################################

# Now, convert back the quantitative images (qT2* and qT1) as they should not be normalized
dat$contrast[dat$contrast=='qt2star'] <- 'qT2star'

dat$normCNR[dat$contrast %in% c('qT2star', 'qT1')] <- dat$normCNR[dat$contrast %in% c('qT2star', 'qT1')]*dat$dim[dat$contrast %in% c('qT2star', 'qT1')]
dat$normSNR_CC[dat$contrast %in% c('qT2star', 'qT1')] <- dat$normSNR_CC[dat$contrast %in% c('qT2star', 'qT1')]*dat$dim[dat$contrast %in% c('qT2star', 'qT1')]
dat$normSNR_CN[dat$contrast %in% c('qT2star', 'qT1')] <- dat$normSNR_CN[dat$contrast %in% c('qT2star', 'qT1')]*dat$dim[dat$contrast %in% c('qT2star', 'qT1')]
dat$normSNR_LCN[dat$contrast %in% c('qT2star', 'qT1')] <- dat$normSNR_LCN[dat$contrast %in% c('qT2star', 'qT1')]*dat$dim[dat$contrast %in% c('qT2star', 'qT1')]
dat$normSNR_RCN[dat$contrast %in% c('qT2star', 'qT1')] <- dat$normSNR_RCN[dat$contrast %in% c('qT2star', 'qT1')]*dat$dim[dat$contrast %in% c('qT2star', 'qT1')]

# FINAL DAT DATAFFRAME, CLEANED
dat$contrast <- as.factor(dat$contrast)
dat$name[dat$name=='Forrest'] <- 'FORREST'

dat <- dat %>% dplyr::select(Database_sub_ses, Label, normSNR_CC, normSNR_CN, normSNR_LCN, normSNR_RCN, normCNR, contrast, name, dim, time, CC)

######################################################################
######################################################################
######################################################################
# NOW CREATE SUMARY TABLE OF DATA, MEANED OVER DATABASES

spat <- dat[1:500,c('Label','normSNR_CC','normSNR_CN','normCNR','name')]                   # take first 500 observations

allinfo <- aggregate(cbind(normSNR_CC, normSNR_CN, normCNR) ~ Label, spat, function(x) cbind(mean(x),sd(x),length(x)))    # mean SNRs and CNRs by label

# create new labels to mean across again 
allinfo$newlab <- c('AGET1','AHEADPDSLAB','AHEADPDWH','AHEADPDSLAB','AHEADPDWH','AHEADPDSLAB','AHEADPDWH','AHEADQSMWH','AHEADQSMWH',
                    'AHEADQSMWH','AHEADQT1SLAB','AHEADQT1WH','AHEADQT1SLAB','AHEADQT1WH','AHEADQT1SLAB','AHEADQT1WH','AHEADQT2STARSLAB',
                    'AHEADQT2STARWH','AHEADQT2STARSLAB','AHEADQT2STARWH','AHEADQT2STARSLAB','AHEADQT2STARWH','AHEADT1SLAB','AHEADT1WH',
                    'AHEADT1SLAB','AHEADT1WH','AHEADT1SLAB','AHEADT1WH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH',
                    'ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH','ATAGFLASH',
                    'ATAGQT1SLAB','ATAGQT1WH','ATAGQT1SLAB','ATAGQT1WH','ATAGQT1SLAB','ATAGQT1WH','ATAGT1SLAB','ATAGT1WH','ATAGT1SLAB','ATAGT1WH',
                    'ATAGT1SLAB','ATAGT1WH','CAMT1','CAMT1','CAMT1','CAMT2','CAMT2','CAMT2','CBSQT1','CBSQT2STAR','CBST1','CBST2STAR','CBST2STAR',
                    'CBST2STAR','CBST2STAR','CBST2STAR','FORRESTSWI','FORRESTT1','FORRESTT2','GSPT1','HCPT1PRO','HCPT1UNPRO','HCPT2PRO','HCPT2UNPRO',
                    'IXIPD','IXIPD','IXIPD','IXIT1','IXIT1','IXIT1','IXIT2','IXIT2','IXIT2','KIRBFLAIR','KIRBT1','LMBBFLAIR','LMBBFLAIR','LMBBQT1',
                    'LMBBQT1','LMBBQT1','LMBBT1','LMBBT1','LMBBT1','MASSFLAIR','MASST1','MASST2','MSCT1','MSCT2','RST1','RST1','RST1','SALDT1',
                    'SALDT1','SALDT1')

#tableinfo <- aggregate(cbind(normSNR[,1], normCNR[,1],SNRsem,CNRsem) ~ newlab, allinfo, mean)
tableinfo <- aggregate(cbind(normSNR_CC[,1],normSNR_CC[,2], normSNR_CN[,1],normSNR_CN[,2], normCNR[,1], normCNR[,2]) ~ newlab, allinfo, mean)
tableinfo$name <- c('AGE',rep('AHEAD',9),rep('ATAG',5),rep('CAM',2),rep('MPICBS',4),rep('FORREST',3),'GSP',rep('HCPYA',4),
                    rep('IXI',3),rep('KIRB',2),rep('MPILMBB',3),rep('MASS',3),rep('MSC',2),'RS','SALD')

# add number of subjects per meaned group
tableinfo$n <- 0
tableinfo$n[tableinfo$name=='AGE'] <- 5
tableinfo$n[tableinfo$name=='AHEAD'] <- 15
tableinfo$n[tableinfo$name=='ATAG'] <- 15
tableinfo$n[tableinfo$name=='CAM'] <- 15
tableinfo$n[tableinfo$name=='MPICBS'] <- 5
tableinfo$n[tableinfo$name=='FORREST'] <- 5
tableinfo$n[tableinfo$name=='GSP'] <- 5
tableinfo$n[tableinfo$name=='HCPYA'] <- 5
tableinfo$n[tableinfo$name=='IXI'] <- 15
tableinfo$n[tableinfo$name=='KIRB'] <- 5
tableinfo$n[tableinfo$name=='MPILMBB'] <- 15
tableinfo$n[tableinfo$name=='MPILMBB'&tableinfo$newlab=='LMBBFLAIR'] <- 10 # no middle aged group, so only 10 subjects instead of 15
tableinfo$n[tableinfo$name=='MASS'] <- 5
tableinfo$n[tableinfo$name=='MSC'] <- 5
tableinfo$n[tableinfo$name=='RS'] <- 15
tableinfo$n[tableinfo$name=='SALD'] <- 15

tableinfo$SNR_CCsem <- tableinfo$V2/sqrt(tableinfo$n)
tableinfo$SNR_CNsem <- tableinfo$V4/sqrt(tableinfo$n)
tableinfo$CNRsem <- tableinfo$V6/sqrt(tableinfo$n)
tableinfo <- tableinfo %>% rename("Label" = "newlab", "normSNR_CC" = "V1", "SNR_CCsd" = "V2", "normSNR_CN" = "V3", "SNR_CNsd" = "V4", "normCNR" = "V5", "CNRsd" = "V6")

###### THEN ADD REMAINING OBSERVATIONS

# add AHEAD data (dat[501:620]) - this data was added after original so created a new frame for ease
toadd <- dat[501:620,]
toadd$normSNR_CC <- toadd$normSNR_CC#/toadd$dim
toadd$normSNR_CN <- toadd$normSNR_CN#/toadd$dim
toadd$normCNR <- toadd$normCNR#/toadd$dim
toaddinfo <- aggregate(cbind(normSNR_CC, normSNR_CN, normCNR) ~ Label, toadd, function(x) cbind(mean(x),sd(x)))
toaddinfo$index <- c(rep('ELDwh',5),rep('MIDwh',5),rep('ELDsb',5),rep('MIDsb',5),rep('YOUsb',5),rep('YOUwh',5))
aveinfo <- aggregate(cbind(normSNR_CC, normSNR_CN, normCNR) ~ index, toaddinfo, mean)
aveinfo$newlab <- rep(c('T2STARsb','T2STARwh'),3)
newinfo <- aggregate(cbind(V1,V2,V3,V4,V5,V6) ~ newlab, aveinfo, mean)
newinfo$name <- 'AHEAD'
newinfo$n <- 15
newinfo$SNR_CCsem <- newinfo$V2/sqrt(newinfo$n)
newinfo$SNR_CNsem <- newinfo$V4/sqrt(newinfo$n)
newinfo$CNRsem <- newinfo$V6/sqrt(newinfo$n)
newinfo <- newinfo %>% rename("Label" = "newlab", "normSNR_CC" = "V1", "SNR_CCsd" = "V2", "normSNR_CN" = "V3", "SNR_CNsd" = "V4", "normCNR" = "V5", "CNRsd" = "V6")

### add second round of new data (RAIDERS, 250, MAASTRICHT, PTBP, DLBS)
toadd2 <- dat[621:665,]
toadd2$normSNR_CC <- toadd2$normSNR_CC#/toadd2$dim
toadd2$normSNR_CN <- toadd2$normSNR_CN#/toadd2$dim
toadd2$normCNR <- toadd2$normCNR#/toadd2$dim
toaddinfo2 <- aggregate(cbind(normSNR_CC, normSNR_CN, normCNR) ~ Label, toadd2, function(x) cbind(mean(x),sd(x)))
toaddinfo2$name <- c("250","DLBS","DLBS","DLBS","MAASTRICHT","MAASTRICHT","MAASTRICHT","PTBP","RAIDER")
toaddinfo2$n <- 5 # even though only 1 subject in 250 database, there are 5 observations
toaddinfo2$SNR_CCsem <- toaddinfo2$normSNR_CC[,2]/sqrt(toaddinfo2$n)
toaddinfo2$SNR_CNsem <- toaddinfo2$normSNR_CN[,2]/sqrt(toaddinfo2$n)
toaddinfo2$CNRsem <- toaddinfo2$normCNR[,2]/sqrt(toaddinfo2$n)


toaddinfo2 <- data.frame("Label" = toaddinfo2$Label, "normSNR_CC" = toaddinfo2$normSNR_CC[,1], "SNR_CCsd" = toaddinfo2$normSNR_CC[,2], "normSNR_CN" = toaddinfo2$normSNR_CN[,1], "SNR_CNsd" = toaddinfo2$normSNR_CN[,2],  "normCNR" = toaddinfo2$normCNR[,1], "CNRsd" = toaddinfo2$normCNR[,2], 
                         "name" = toaddinfo2$name, "n" = toaddinfo2$n, "SNR_CCsem" = toaddinfo2$SNR_CCsem, "SNR_CNsem" = toaddinfo2$SNR_CNsem, "CNRsem" = toaddinfo2$CNRsem) # split columns and change names 


ALLDAT <- rbind(tableinfo, newinfo, toaddinfo2)              # combine all data

DLBSdat <- colMeans(ALLDAT[ALLDAT$name=='DLBS',c(2:7,9:12)])

addDLBS <- data.frame('Label'='DLBST1','normSNR_CC'= DLBSdat[1], 'SNR_CCsd'=DLBSdat[2], 'normSNR_CN'=DLBSdat[3], 'SNR_CNsd'=DLBSdat[4], 'normCNR'=DLBSdat[5],'CNRsd'=DLBSdat[6], 
           'name'='DLBS', 'n'=5,  'SNR_CCsem'=DLBSdat[8],  'SNR_CNsem'=DLBSdat[9],  'CNRsem'=DLBSdat[10])
rownames(addDLBS) <- NULL

ALLDAT<- ALLDAT[c(1:47,51:55),]

ALLDAT <- rbind(ALLDAT,addDLBS)

ALLDAT$SNRplusSEM <- paste(round(ALLDAT$normSNR_CC,1),"±",round(ALLDAT$SNR_CCsem,1), sep=" ")
ALLDAT$SNR_CNplusSEM <- paste(round(ALLDAT$normSNR_CN,1),"±",round(ALLDAT$SNR_CNsem,1), sep=" ")
ALLDAT$CNRplusSEM <- paste(round(abs(ALLDAT$normCNR),1),"±",round(ALLDAT$CNRsem,1), sep=" ")












