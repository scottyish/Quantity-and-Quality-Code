  ####### SNR and CNR data for all databases 
  #rm(list=ls())
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
  library(lmerTest)
  font_import()
  loadfonts(device="win")
  fonts
  
  options(max.print=100000)
  
  # set directories
  data_dir <- '/Users/scotti/surfdrive/PhD/meta_analysis/Databases/R_stuff'

  # save graphs? set value to 0 for no, 1 for yes
  save_stuff <- ''
  #save_stuff <- 1
  
  dat$normCNR <- abs(dat$normCNR)
  
  ####################################################################################
  ################################################ T1 data
  ####################################################################################
  
  
  T1dat <- dat[dat$contrast=='T1',] # dataframe of just T1 images
  T1dat <- T1dat[!grepl("SLAB", T1dat$Label),] # remove slabs for now 
  
  ### remove unprocessed HCP data
  #T1dat <- T1dat[T1dat$Label!='HCPT1UNP',]
  ## remove processed HCP data
  T1dat <- T1dat[T1dat$Label!='HCPT1PRO',]
  T1dat <- T1dat[T1dat$Label!='250proT1',]
  
  ### T1 graphs for SNR and CNR together 
  # CNR info
  deets <- aggregate(normCNR ~ name, T1dat, function(x) cbind(mean(x),sd(x),length(x)))       # get mean
  deets$sem <- deets$normCNR[,2]/sqrt(deets$normCNR[,3])                                      # get SEM
  deets$name <- factor(deets$name, levels=deets$name[order(deets$normCNR[,1])])               # reorder factor levels
  
  deetsCNR <- deets
  
  # SNR info
  deets <- aggregate(normSNR_CC ~ name, T1dat, function(x) cbind(mean(x),sd(x),length(x)))      # get mean
  deets$sem <- deets$normSNR_CC[,2]/sqrt(deets$normSNR_CC[,3])                                     # get SEM

  for (change in c('250','AHEAD','ATAG','MAASTRICHT','MPICBS')){
    deets$name[deets$name==change] <- paste0("'",change)
  }
  
  deets$name <- factor(deets$name, levels=deets$name[order(deets$normSNR_CC[,1])])              # reorder factor levels                                                                          
  
  deetsdouble <- cbind.data.frame(rep(deets$name,2), c(deetsCNR$normCNR[,1],deets$normSNR_CC[,1]), c(deetsCNR$sem,  deets$sem), rep(c('CNR','SNR'),each=20))
  colnames(deetsdouble) <- c('name', 'value','sem','Type')
  deetsdouble$Tesla <- 3
  deetsdouble$Tesla[deetsdouble$name%in%c("'AHEAD","'ATAG","'MAASTRICHT","'MPICBS","'250")] <- 7
  deetsdouble$Tesla <- as.factor(deetsdouble$Tesla)
  
  if (save_stuff == 1){tiff(file.path(save_dir, "T1_SNR_CNR_dotplot.tiff"), units="in", width=5,height=3,res=300)}                        # START image to be saved
  
  p <- ggplot(deetsdouble, aes(x=value,y=name,color=Type,fill=Type)) +
    geom_dotplot(binaxis='y', stackdir="center", dotsize=0.8, binwidth=0.5) +
    geom_errorbarh(aes(xmin=value-sem, xmax=value+sem), height=.3, size=0.3) +
    geom_vline(xintercept=mean(deetsdouble$value[1:15]),colour='grey50',linetype='dashed',size=0.3) + # maybe change black to grey50
    geom_vline(xintercept=mean(deetsdouble$value[16:30]),colour='black',linetype='dashed', size=0.3) + 
    scale_fill_manual(values=c('grey50','black')) +
    scale_color_manual(values=c('grey50','black')) +
    labs(x = "Value", y = "Database") +
    ggtitle("SNR and CNR estimations for T1 weighted images") +
    scale_x_continuous(breaks=seq(0,320,50)) +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7),
          axis.text=element_text(size=4),
          axis.title=element_text(size=8,face="bold"),
          plot.title = element_text(size=10,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"))
  
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ########## Plot sample size of database against T1 SNR-CC #########
  # plot SNR against number of participants
  # label each SNR with the database name 
  
  # add sample sizes to dataframe
  deetsdouble$size <- rep(c(1, # 250
                          131, # Ageility
                          106, #AHEAD
                          53,# ATAG
                          280,#CAMCAN
                          315,#DLBS
                          20,#forrest
                          1570,#GSP
                          1206,#HCPYA
                          619,#IXI
                          21,#kir21
                          5,#maastricht
                          1,#massive
                          28,#mpicbs
                          321,#mpilmbb
                          10,#MSC
                          1000,#NKIRS
                          120,#PTBP
                          11,#RAIDERS
                          494),#SALD
                          2)
  touse <- deetsdouble[deetsdouble$Type=='SNR',]
  
  touse$name <- revalue(touse$name, c("'AHEAD"="AHEAD","'ATAG"="ATAG","'MAASTRICHT"="MAASTRICHT","'MPICBS"="MPICBS","'250"="250"))
  
  # make value and sem log (this is a tad complicated)
  touse$lowersem <- log(touse$value-touse$sem)
  touse$uppersem <- log(touse$value+touse$sem)
  
  touse$value <- log(touse$value)
  
  if (save_stuff == 1){tiff(file.path(save_dir,"samplesize_vs_SNRlog.tiff"), units="in", width=5,height=3,res=300)}
  
  # plot voxel dims against scan time 
  p <- ggplot(touse, aes(y=value, x=size))+
    #geom_smooth(method='lm',formula=y~x, se=FALSE,col="black",linetype='dashed', size=0.4) +
    geom_errorbar(aes(ymin=lowersem, ymax=uppersem),width=.2,size=0.3) +
    geom_point(size=1.2, aes(shape=Tesla)) +
    scale_x_continuous(name="Sample size", breaks=seq(0,1750,250)) +
    scale_y_continuous(name="log(SNR)") + 
    ggtitle("Sample size vs SNR")+
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7), # removed hjust=1
          axis.text=element_text(size=4),
          axis.title=element_text(size=8,face="bold"),
          plot.title = element_text(size=10,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"))+ 
    
    scale_shape_discrete(name="Field\nstrength",
                         breaks=c("3","7"),
                         labels=c("3T","7T")) +
    
    geom_label_repel(aes(label = name),
                     #box.padding   = 0.35, 
                     #min.segment.length = 0.1,
                     force=3,
                     max.iter=200,
                     point.padding = 0.5,
                     label.padding = 0.1,
                     size = 1.2,
                     fill = 'white',
                     label.size=0,
                     segment.color = 'grey10',label.r = 0.05, segment.size = 0.25) 
  
  p
  
  if (save_stuff ==1){dev.off()}   
 
  ######### plot time against voxel dimensions ###########
  # Define fit for legend
  
  timevox <- aggregate(dim ~ time + name, T1dat, mean)
  timevox$Tesla <- 3
  timevox$Tesla[timevox$name%in%c('AHEAD','ATAG','MAASTRICHT','MPICBS','250')] <- 7
  timevox$Tesla <- as.factor(timevox$Tesla)
  
  fit <- lm(dim ~ time, data=timevox)
  sumfit <- summary(fit)
  
  if (save_stuff == 1){tiff(file.path(save_dir,"dim_vs_scantime_new2.tiff"), units="in", width=5.5,height=3,res=300)}
  #if (save_stuff == 1){pdf(file.path(save_dir,"dim_vs_scantime_new.pdf"), width=5.5,height=3)}
  
  # plot voxel dims against scan time 
  p <- ggplot(timevox, aes(y=dim, x=time))+
    geom_smooth(method='lm',formula=y~x, se=FALSE,col="black",linetype='dashed', size=0.4) +
    geom_point(size=1.2,aes(shape=Tesla)) +
    scale_x_continuous(name="Scan time (s)", breaks=seq(0,3250,250)) +
    scale_y_continuous(name="Voxel volume (mm3)", breaks=seq(0,2,0.5)) + 
    ggtitle("Scan time vs acquired voxel dimensions")+
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7), # removed hjust=1
          axis.text=element_text(size=4),
          axis.title=element_text(size=8,face="bold"),
          plot.title = element_text(size=10,face="bold"),
          legend.title=element_text(size=6.5),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm")) +
  
    scale_shape_discrete(name="Field\nstrength",
                         breaks=c("3","7"),
                         labels=c("3T","7T")) +
    
  geom_label_repel(aes(label = name),
                   #box.padding   = 0.35, 
                   point.padding = 0.4,
                   label.padding = 0.1,
                   max.iter=100,
                   size = 1.2,
                   fill = 'white',
                   label.size=0,
                   segment.color = 'grey10',label.r = 0.05, segment.size = 0.25) 
  
  p
  
  # Setup legend for coefficients
  grid.newpage()
  vpa_ <- viewport(width = 1, height = 1)
  print(p, vp = vpa_)

  grid.text(paste0("Adj R-squared = ",round(sumfit$adj.r.squared,2)),x=0.75,y=0.86,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("Intercept = ",round(sumfit$coefficients[1],2)),x=0.75,y=0.82,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("Slope = ",round(sumfit$coefficients[2],5)),x=0.75,y=0.78,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("F statistic = ", round(sumfit$fstatistic[1],2)),x=0.75,y=0.74,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("DF = ", sumfit$df[2]),x=0.75,y=0.70,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("p = ", round(sumfit$coefficients[8],4)),x=0.75,y=0.66,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  
  if (save_stuff ==1){dev.off()}   
  
  ######## plot time against normalized SNR ##########
  # Define dataframe and get SEMs
  timeSNR <- aggregate(cbind(time, normSNR_CC) ~ name, T1dat, mean)
  timeSNRsd <- aggregate(cbind(time, normSNR_CC) ~ name, T1dat, sd)
  timeSNRlen <- aggregate(cbind(time, normSNR_CC) ~ name, T1dat, length)
  SEM <- timeSNRsd[,3]/sqrt(timeSNRlen[,3])
  timeSNR$sem <- SEM
  timeSNR$Tesla <- 3
  timeSNR$Tesla[timeSNR$name%in%c('AHEAD','ATAG','MAASTRICHT','MPICBS','250')] <- 7
  timeSNR$Tesla <- as.factor(timeSNR$Tesla)
  
  # define model fits for legend
  fit <- lm(normSNR_CC ~ time, data=timeSNR)
  sumfit <- summary(fit)
  
  if (save_stuff == 1){tiff(file.path(save_dir,"SNR_vs_scantime_3T_and_7T.tiff"), units="in", width=5.5,height=3,res=300)}
  #if (save_stuff == 1){pdf(file.path(save_dir,"SNR_vs_scantime_new2.pdf"),  width=5.5,height=3)}
  
  set.seed(1337)
  # plot data
  p <- ggplot(timeSNR, aes(y=normSNR_CC, x=time))+
    geom_smooth(method='lm',formula=y~x, se=FALSE,col="black",linetype="dashed", size=0.4) +
    geom_point(size=1.2,aes(shape=Tesla)) +
    geom_errorbar(aes(ymin=normSNR_CC-sem, ymax=normSNR_CC+sem), width=.2,size=0.3, position=position_dodge(.9)) +
    scale_x_continuous(name="Scan time (s)", breaks=seq(0,3250,250)) +
    scale_y_continuous(name="SNR", breaks=seq(0,320,50)) + 
    ggtitle("Scan time vs acquired SNR") +
    theme(text=element_text(family="Times New Roman", size=8),
        axis.text.y = element_text(size=7),
        axis.text.x = element_text(size=7), 
        axis.text=element_text(size=4),
        axis.title=element_text(size=8,face="bold"),
        plot.title = element_text(size=10,face="bold"),
        legend.title=element_text(size=6.5),
        legend.key.size = unit(0.35, "cm"),
        legend.position="right",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0.2, "cm")) +
    
    scale_shape_discrete(name="Field\nstrength",
                         breaks=c("3","7"),
                         labels=c("3T","7T")) 
  
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
  
  # Setup legend for coefficients
  grid.newpage()
  vpa_ <- viewport(width = 1, height = 1)
  print(p, vp = vpa_)
  
  grid.text(paste0("Adj R-squared = ",round(sumfit$adj.r.squared,2)),x=0.1,y=0.86,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("Intercept = ",round(sumfit$coefficients[1],2)),x=0.1,y=0.82,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("Slope = ",round(sumfit$coefficients[2],2)),x=0.1,y=0.78,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("F statistic = ", round(sumfit$fstatistic[1],2)),x=0.1,y=0.74,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("DF = ", sumfit$df[2]),x=0.1,y=0.7,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  grid.text(paste0("p = ", round(sumfit$coefficients[8],7)),x=0.1,y=0.66,just="left", gp=gpar(fontsize=7,fontfamily="Times"))
  
  if (save_stuff ==1){dev.off()}   
  
  ######## plot time against normalized SNR for BOTH 3T and 7T side-by-side ONLY ##########
  # Define dataframe and get SEMs
  timeSNR <- aggregate(cbind(time, normSNR_CC) ~ name, T1dat, mean)
  timeSNRsd <- aggregate(cbind(time, normSNR_CC) ~ name, T1dat, sd)
  timeSNRlen <- aggregate(cbind(time, normSNR_CC) ~ name, T1dat, length)
  SEM <- timeSNRsd[,3]/sqrt(timeSNRlen[,3])
  timeSNR$sem <- SEM
  timeSNR$Tesla <- 3
  timeSNR$Tesla[timeSNR$name%in%c('AHEAD','ATAG','MAASTRICHT','MPICBS','250')] <- 7
  timeSNR$Tesla <- as.factor(timeSNR$Tesla)
  timeSNR3 <- timeSNR[timeSNR$Tesla==3,]
  timeSNR3 <- droplevels(timeSNR3)
  timeSNR7 <- timeSNR[timeSNR$Tesla==7,]
  timeSNR7 <- droplevels(timeSNR7)
  
  # define model fits for legend
  fit3 <- lm(normSNR_CC ~ time, data=timeSNR3)
  sumfit3 <- summary(fit3)
  
  fit7 <- lm(normSNR_CC ~ time, data=timeSNR7)
  sumfit7 <- summary(fit7)
  
  if (save_stuff == 1){tiff(file.path(save_dir,"SNR_vs_scantime_new2.tiff"), units="in", width=5.5,height=3,res=300)}

  set.seed(1337)
  # plot data
  p3 <- ggplot(timeSNR3, aes(y=normSNR_CC, x=time))+
    geom_smooth(method='lm',formula=y~x, se=FALSE,col="black",linetype="dashed", size=0.4) +
    geom_point(size=1.2) +
    geom_errorbar(aes(ymin=normSNR_CC-sem, ymax=normSNR_CC+sem), width=.2,size=0.3, position=position_dodge(.9)) +
    scale_x_continuous(name="Scan time (s)", breaks=seq(0,800,100)) +
    scale_y_continuous(name="SNR", breaks=seq(0,100,20)) + 
    ggtitle("Scan time vs SNR (3T)") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7), 
          axis.text=element_text(size=4),
          axis.title=element_text(size=8,face="bold"),
          plot.title = element_text(size=10,face="bold"),
          legend.title=element_text(size=6.5),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm")) +
  
  set.seed(1337)
  p3 <- p3 + geom_label_repel(aes(label = name),
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

  set.seed(1337)
  # plot data
  p7 <- ggplot(timeSNR7, aes(y=normSNR_CC, x=time))+
    geom_smooth(method='lm',formula=y~x, se=FALSE,col="black",linetype="dashed", size=0.4) +
    geom_point(size=1.2) +
    geom_errorbar(aes(ymin=normSNR_CC-sem, ymax=normSNR_CC+sem), width=.2,size=0.3, position=position_dodge(.9)) +
    scale_x_continuous(name="Scan time (s)", breaks=seq(500,3250,500)) +
    scale_y_continuous(name="SNR", breaks=seq(0,320,50)) + 
    ggtitle("Scan time vs SNR (7T)") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7), 
          axis.text=element_text(size=4),
          axis.title=element_text(size=8,face="bold"),
          plot.title = element_text(size=10,face="bold"),
          legend.title=element_text(size=6.5),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm")) +
  
  set.seed(1337)
  p7 <- p7 + geom_label_repel(aes(label = name),
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
  
  grid.arrange(p3, p7, ncol=2)
  
  grid.text(paste0("Adj R-squared = ",round(sumfit3$adj.r.squared,2)),x=0.1,y=0.86,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("Intercept = ",round(sumfit3$coefficients[1],2)),x=0.1,y=0.83,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("Slope = ",round(sumfit3$coefficients[2],2)),x=0.1,y=0.80,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("F statistic = ", round(sumfit3$fstatistic[1],2)),x=0.1,y=0.77,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("DF = ", sumfit3$df[2]),x=0.1,y=0.74,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("p = ", round(sumfit3$coefficients[8],3)),x=0.1,y=0.71,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  
  grid.text(paste0("Adj R-squared = ",round(sumfit7$adj.r.squared,2)),x=0.6,y=0.86,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("Intercept = ",round(sumfit7$coefficients[1],2)),x=0.6,y=0.83,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("Slope = ",round(sumfit7$coefficients[2],2)),x=0.6,y=0.8,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("F statistic = ", round(sumfit7$fstatistic[1],2)),x=0.6,y=0.77,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("DF = ", sumfit7$df[2]),x=0.6,y=0.74,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  grid.text(paste0("p = ", round(sumfit7$coefficients[8],3)),x=0.6,y=0.71,just="left", gp=gpar(fontsize=4,fontfamily="Times"))
  
  if (save_stuff ==1){dev.off()}   

  ####### multivariate glm
  glmframe <- data.frame('time'=timevox$time,'dim'=timevox$dim,'SNR'=timeSNR$normSNR_CC)
  test.fit <- lm(cbind(dim, SNR) ~ time, data=glmframe)
  summary(test.fit) # bonferroni correctionn = 0.05/2
  
  ###### bar graph SNR data based on mean CN instead of CC
  ######################################################

  deets <- aggregate(normSNR_CN ~ name, T1dat, function(x) cbind(mean(x),sd(x),length(x)))
  deets$sem <- deets$normSNR_CN[,2]/sqrt(deets$normSNR_CN[,3])

  deets$name <- factor(deets$name, levels=deets$name[order(deets$normSNR_CN[,1])])

  p <- ggplot(deets, aes(y=deets$normSNR_CN[,1], x=deets$name)) +
    geom_bar( stat="identity") +
    geom_errorbar(aes(ymin=deets$normSNR_CN[,1]-deets$sem, ymax=deets$normSNR_CN[,1]+deets$sem), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    ggtitle("SNR estimations for mean of caudate nuclei")

  print(p)

  ##### grouped bar plot for CNs
  both <- aggregate(cbind(normSNR_LCN, normSNR_RCN) ~ name, T1dat, function(x) cbind(mean(x),sd(x),length(x)))
  both$semLCN <- both$normSNR_LCN[,2]/sqrt(both$normSNR_LCN[,3])
  both$semRCN <- both$normSNR_RCN[,2]/sqrt(both$normSNR_RCN[,3])
  both <- both[order(both$normSNR_LCN[,1]),]

  for (change in c('250','AHEAD','ATAG','MAASTRICHT','MPICBS')){
    both$name[both$name==change] <- paste0("'",change)
  }

  new <- data.frame('Database' = rep(both$name,2),'CN' = c(both$normSNR_LCN[,1],both$normSNR_RCN[,1]), Location = c(rep('Left',20),rep('Right',20)),
                    'sem' = c(both$semLCN,both$semRCN), meanCN = rep((both$normSNR_LCN[,1]+both$normSNR_RCN[,1])/2,2))

  if (save_stuff == 1){tiff(file.path(save_dir, "compare_SNR_CNs.tiff"), units="in", width=5,height=3,res=300)}

  p <- ggplot(new, aes(fill=Location, y=CN, x=Database)) +
    geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=CN-sem, ymax=CN+sem), width=.2, size=0.3,
                position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    ggtitle("SNR estimations for left and right caudate nucleus") +
    theme(text=element_text(family="Times New Roman"),
          axis.text=element_text(size=7),
          axis.text.y = element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          axis.text.x = element_text(angle=45, hjust=1),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"),
          legend.text = element_text(size=6)) +
          scale_fill_manual(values=c("grey30","grey70"))
  print(p)

  if (save_stuff ==1){dev.off()}
  
  ## check difference between SNR of left and right CN (t tests)
  testL <- aggregate(normSNR_LCN ~ name, T1dat, mean)
  testR <- aggregate(normSNR_RCN ~ name, T1dat, mean)
  t.test(testL$normSNR_LCN,testR$normSNR_RCN,paired=T) # not significantly different 
  
  library(BayesFactor)
  ttestBF(testL$normSNR_LCN,testR$normSNR_RCN,paired=T) # BF = 0.345
  
  ####################################################################################
  ################################################ T2 data
  ####################################################################################
  
  ###### data for T2 
  T2dat <- dat[dat$contrast=='T2',]
  T2dat <- T2dat[!T2dat$name=='MPICBS',]
  
  ### remove unprocessed HCP data
  #T2dat <- T2dat[T2dat$Label!='HCPT2UNP',]
  ## remove processed HCP data
  T2dat <- T2dat[T2dat$Label!='HCPT2PRO',]
  
  ##### plot data for T2 SNR
  deets <- aggregate(normSNR_CC ~ name, T2dat, function(x) cbind(mean(x),sd(x),length(x)))
  deets$sem <- deets$normSNR_CC[,2]/sqrt(deets$normSNR_CC[,3])
  
  deets$name <- factor(deets$name, levels=deets$name[order(deets$normSNR_CC[,1])])
  
  if (save_stuff == 1){tiff(file.path(save_dir, "T2_SNR_comparison.tiff"), units="in", width=5,height=3,res=300)}
  
  p <- ggplot(deets, aes(y=deets$normSNR_CC[,1], x=deets$name)) + 
    geom_bar( stat="identity") +
    geom_errorbar(aes(ymin=deets$normSNR_CC[,1]-deets$sem, ymax=deets$normSNR_CC[,1]+deets$sem), width=.2, size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    ggtitle("SNR estimations for T2 weighted images") +
    theme(text=element_text(family="Times New Roman"),
          axis.text.y = element_text(size=8),
          axis.text=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ##### plot data for T2 CNR
  deetsSNR_CC<-deets
  
  deets <- aggregate(normCNR ~ name, T2dat, function(x) cbind(mean(x),sd(x),length(x)))
  deets$sem <- deets$normCNR[,2]/sqrt(deets$normCNR[,3])
  
  deetsCNR <- deets
  
  deets$name <- factor(deets$name, levels=deets$name[order(deets$normCNR[,1])])
  
  if (save_stuff == 1){tiff(file.path(save_dir, "T2_CNR_comparison.tiff"), units="in", width=5,height=3,res=300)}
  
  p <- ggplot(deets, aes(y=deets$normCNR[,1], x=deets$name)) + 
    geom_bar( stat="identity") +
    geom_errorbar(aes(ymin=deets$normCNR[,1]-deets$sem, ymax=deets$normCNR[,1]+deets$sem), width=.2,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "CNR") +
    ggtitle("CNR estimations for T2 weighted images") +
    theme(text=element_text(family="Times New Roman"),
          axis.text=element_text(size=14),
          axis.text.y = element_text(size=18),
          axis.title=element_text(size=16,face="bold"),
          plot.title = element_text(size=20,face="bold"))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  #### dotplpot for CNR and SNR values together in same plot # use this for paper 
  deetsdouble <- cbind.data.frame(rep(deetsSNR_CC$name,2), c(deetsSNR_CC$normSNR_CC[,1],deetsCNR$normCNR[,1]), c(deetsSNR_CC$sem,  deetsCNR$sem), rep(c('SNR','CNR'),each=6))
  colnames(deetsdouble) <- c('name', 'value','sem','Type')
  
  if (save_stuff == 1){tiff(file.path(save_dir, "T2_CNR_SNR_dotplot.tiff"), units="in", width=5,height=3,res=300)}
  
  p <- ggplot(deetsdouble, aes(x=value,y=name,color=Type,fill=Type)) +
    geom_dotplot(binaxis='y', stackdir="center", dotsize=0.4, binwidth=0.5) +
    geom_errorbarh(aes(xmin=value-sem, xmax=value+sem), height=.3, size=0.3) +
    geom_vline(xintercept=mean(deetsdouble$value[1:6]),colour='black',linetype='dashed', size = 0.3) + 
    geom_vline(xintercept=mean(deetsdouble$value[7:12]),colour='grey50', linetype='dashed', size = 0.3) + 
    scale_fill_manual(values=c('grey50','black')) +
    scale_color_manual(values=c('grey50','black')) +
    labs(x = "Value", y = "Database") +
    ggtitle("SNR and CNR estimations for T2 weighted images") +
    scale_x_continuous(breaks=seq(0,70,10)) +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7),
          axis.text=element_text(size=7),
          axis.title=element_text(size=8,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"))
  
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ####################################################################################
  ################################################ COMPARE AGES - T1 data SNR CC
  ####################################################################################
  
  ageDat <- T1dat[T1dat$name%in%c('AHEAD','ATAG','CAMCAN','DLBS','IXI','MPILMBB','NKIRS','SALD'),]
  
  model.age <- T1dat[!T1dat$name%in%c('MAASTRICHT','RAIDERS','PTBP','MASSIVE', '250'),] # no ages for these databases OR only 1 subject
  model.age$Age <- c(24, 24, 25, 21, 25, 46, 38, 45, 40, 50, 65, 67, 72, 75, 73, # IXI
                     68, 63, 67, 75, 62, 56, 60, 55, 58, 54, 25, 23, 26, 27, 23, #ATAG
                     25, 25, 30, 28, 22, # Kir21
                     27, 19, 18, 21, 17, # Ageility
                     19, 23, 23, 33, 19, # GSP
                     31, 26, 25, 23, 22, # MPI-CBS
                     #25, 25, 25, 25, 25, # MASSIVE
                     24, 25, 23, 26, 25, 44, 45, 45, 47, 49, 72, 74, 73, 80, 77, # SALD  
                     23, 23, 23, 23, 23, 43, 38, 48, 43, 48, 63, 68, 68, 68, 63, # MPILMBB
                     34, 34, 29, 31, 26, # MSC
                     33, 28, 23, 23, 38, # Forrest
                     24, 18, 20, 19, 22, 47, 45, 45, 42, 48, 77, 77, 78, 75, 71, # CAMCAN
                     21, 20, 23, 20, 22, 53, 50, 42, 41, 34, 81, 67, 63, 70, 75, #NKI-RS
                     28, 28, 28, 28, 24, # HCPYA
                     25, 21, 24, 25, 26, 51, 47, 40, 45, 49, 71, 71, 74, 77, 67, # STW - AHEAD
                     24, 24, 22, 23, 20, 42, 51, 36, 37, 41, 75, 86, 81, 74, 79#, # DLBS
                     #7, 13, 15, 15, 12, # PTBP
                     #35, 35, 35, 35, 35 # 250
                     
  )
  
  young <- dplyr::filter(ageDat, grepl("YOU",Database_sub_ses))
  young$Age <- 'Young'
  middle <- dplyr::filter(ageDat, grepl("MID",Database_sub_ses))
  middle$Age <- 'Middle'
  elderly <- dplyr::filter(ageDat, grepl("ELD",Database_sub_ses))
  elderly$Age <- 'Elderly'
  
  ageDat <- rbind(young,middle,elderly)
  ageDat$Age <- as.factor(ageDat$Age)
  
  ageDat$name[ageDat$name=='AHEAD'] <- "'AHEAD"
  ageDat$name[ageDat$name=='ATAG'] <- "'ATAG"
  
  agedeets <-aggregate(normSNR_CC ~ name + Age, ageDat, function(x) cbind(mean(x),sd(x),length(x)))
  agedeets$sem <- agedeets$normSNR_CC[,2]/sqrt(agedeets$normSNR_CC[,3])
  agedeets$Age <- factor(agedeets$Age,levels(agedeets$Age)[c(3,2,1)]) #reorder age gorups to look nicer
  
  if (save_stuff == 1){tiff(file.path(save_dir, "Age_comparison_T1_SNR.tiff"), units="in", width=5,height=3,res=300)}
  
  p <- ggplot(agedeets, aes(fill=Age, y=normSNR_CC[,1], x=name)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=normSNR_CC[,1]-sem, ymax=normSNR_CC[,1]+sem), width=.2,size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    scale_fill_manual(values=c('grey70','grey35','black')) +
    ggtitle("SNR estimations of CC for T1w images across age groups") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text=element_text(size=7),
          axis.text.y=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ######################################################
  ### FINAL MIXED EFEFCTS MODELS
  ### ANOVA MODEL COMPARISON
  
  library(lme4)
  library(dplyr)
  library(blme)
  
  # sort data
  new.mod.age <- model.age %>% dplyr::select(name, normSNR_CC, normSNR_CN, normCNR, Age) # frame with Database, SNR of CC, SNR of CN, 
  new.mod.age$Age <- as.integer(new.mod.age$Age)
  new.mod.age$name <- as.factor(new.mod.age$name)
  new.mod.age$tmp <- new.mod.age$name
  
  #### FREQUENTIST
  ## SNR of CC
  CC.aov <- aov(normSNR_CC ~ Age + Error(name), data=new.mod.age) # aov CC
  CC.model <- lmer(normSNR_CC ~ Age + (1|name), data=new.mod.age) # CC SNR model
  
  CC.model1 <- lmer(normSNR_CC ~ Age + (1|name), data=new.mod.age, REML=F) # full model
  CC.model2 <- lmer(normSNR_CC ~ (1|name), data=new.mod.age, REML=F) # null model
  #CC.model3 <- lmer(normSNR ~ Age, data=new.mod.age, REML=F) # Age model
  CC.model3 <- lmer(normSNR_CC ~ Age*(1|name), data=new.mod.age, REML=F) # interaction model
  #CC.model4 <- lmer(normSNR ~ Age + name )
  
  anova(CC.model1, CC.model2, CC.model3) # significant
  
  ## SNR of CN
  CN.aov <- aov(normSNR_CN ~ Age + Error(name), data=new.mod.age) # aov CN
  CN.model <- lmer(normSNR_CN ~ Age + (1|name), data=new.mod.age) # CC SNR model
  CN.model1 <- lmer(normSNR_CN ~ Age + (1|name), data=new.mod.age, REML=F) # full model
  CN.model2 <- lmer(normSNR_CN ~ (1|name), data=new.mod.age, REML=F) # null model
  
  anova(CN.model1, CN.model2) # signnificant
  
  ## CNR
  CNR.aov <- aov(normCNR ~ Age + Error(name), data=new.mod.age) # aov CNR
  CNR.model <- lmer(normCNR ~ Age + (1|name), data=new.mod.age) # CC SNR model
  CNR.model1 <- lmer(normCNR ~ Age + (1|name), data=new.mod.age, REML=F) # full model
  CNR.model2 <- lmer(normCNR ~ (1|name), data=new.mod.age, REML=F) # null model
  
  anova(CNR.model1, CNR.model2) # not significant
  
  ### BAYESIAN
  ## SNR of CC
  blm1 <- lmBF(normSNR_CC ~ Age*name, data=new.mod.age, whichRandom = "tmp" )
  blm2 <- lmBF(normSNR_CC ~ Age + name, data=new.mod.age, whichRandom= "tmp" )
  blm3 <- lmBF(normSNR_CC ~ Age, data=new.mod.age, whichRandom= "tmp" )
  blm4 <- lmBF(normSNR_CC ~ name, data=new.mod.age, whichRandom= "tmp" )
  
  #blm2/blm1
  #blm3/blm2
  blm2/blm4 # not signnificant 
  
  ## SNR of CN
  blm1 <- lmBF(normSNR_CN ~ Age*name, data=new.mod.age, whichRandom = "tmp" )
  blm2 <- lmBF(normSNR_CN ~ Age + name, data=new.mod.age, whichRandom= "tmp" )
  blm3 <- lmBF(normSNR_CN ~ Age, data=new.mod.age, whichRandom= "tmp" )
  blm4 <- lmBF(normSNR_CN ~ name, data=new.mod.age, whichRandom= "tmp" )
  
  #blm2/blm1
  #blm3/blm2
  blm2/blm4 #significant 
  
  ## CNR 
  blm1 <- lmBF(normCNR ~ Age*name, data=new.mod.age, whichRandom = "tmp" )
  blm2 <- lmBF(normCNR ~ Age + name, data=new.mod.age, whichRandom= "tmp" )
  blm3 <- lmBF(normCNR ~ Age, data=new.mod.age, whichRandom= "tmp" )
  blm4 <- lmBF(normCNR ~ name, data=new.mod.age, whichRandom= "tmp" )
  
  #blm2/blm1
  #blm3/blm2
  blm2/blm4 # not signicant 
  
  ggplot(new.mod.age, aes(x=Age, y=normCNR, color=name)) + 
    geom_point()+
    geom_smooth(method=lm)
  
  ggplot(new.mod.age, aes(x=Age, y=normSNR_CC, color=name)) + 
    geom_point()+
    geom_smooth(method=lm)
  
  ggplot(new.mod.age, aes(x=Age, y=normSNR_CN, color=name)) + 
    geom_point()+
    geom_smooth(method=lm)
  

  ####################################################################################
  ################################################ COMPARE AGES - T1 data CNR
  ####################################################################################
  
  ageDat <- T1dat[T1dat$name%in%c('AHEAD','ATAG','CAMCAN','DLBS','IXI','MPILMBB','NKIRS','SALD'),]
  
  young <- dplyr::filter(ageDat, grepl("YOU",Database_sub_ses))
  young$Age <- 'Young'
  middle <- dplyr::filter(ageDat, grepl("MID",Database_sub_ses))
  middle$Age <- 'Middle'
  elderly <- dplyr::filter(ageDat, grepl("ELD",Database_sub_ses))
  elderly$Age <- 'Elderly'
  
  ageDat <- rbind(young,middle,elderly)
  ageDat$Age <- as.factor(ageDat$Age)
  
  ageDat$name[ageDat$name=='AHEAD'] <- "'AHEAD"
  ageDat$name[ageDat$name=='ATAG'] <- "'ATAG"
  
  agedeets <-aggregate(normCNR ~ name + Age, ageDat, function(x) cbind(mean(x),sd(x),length(x)))
  agedeets$sem <- agedeets$normCNR[,2]/sqrt(agedeets$normCNR[,3])
  agedeets$Age <- factor(agedeets$Age,levels(agedeets$Age)[c(3,2,1)]) #reorder age gorups to look nicer
  
  if (save_stuff == 1){tiff(file.path(save_dir, "age_comparison_T1_CNR.tiff"), units="in", width=5,height=3,res=300)}
  
  
  p <- ggplot(agedeets, aes(fill=Age, y=normCNR[,1], x=name)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=normCNR[,1]-sem, ymax=normCNR[,1]+sem), width=.2, size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "CNR") +
    scale_fill_manual(values=c('grey70','grey35','black')) +
    ggtitle("CNR estimations for T1w images across age groups") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text=element_text(size=7),
          axis.text.y=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm")
          )
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ####################################################################################
  ################################################ COMPARE AGES - T2 data SNR
  ####################################################################################
  
  ageDat <- T2dat[T2dat$name%in%c('CAMCAN','IXI'),]
  
  mod.T2.dat <- T2dat[!T2dat$name%in%c('MASSIVE'),]
  mod.T2.dat$Age <- c(24, 22, 27, 25, 24, 43,48, 43, 42, 37, 70, 69, 73, 75, 70, # IXI
                  #25, 25, 25, 25, 25, # MASSVE
                  34, 34, 29, 31, 26, # MSC
                  33, 28, 23, 23, 38, # Forrest
                  18, 23, 24, 20, 20, 47, 50,43, 47, 45, 63, 77, 85, 79, 81, # CAMCAN
                  28, 28, 28, 28, 24  # HCPYA
                  )
  
  T2.age <- mod.T2.dat %>% dplyr::select(name, normSNR_CC, normSNR_CN, normCNR, Age) # frame with Database, SNR of CC, SNR of CN, raw ages
  
  young <- dplyr::filter(ageDat, grepl("YOU",Database_sub_ses))
  young$Age <- 'Young'
  middle <- dplyr::filter(ageDat, grepl("MID",Database_sub_ses))
  middle$Age <- 'Middle'
  elderly <- dplyr::filter(ageDat, grepl("ELD",Database_sub_ses))
  elderly$Age <- 'Elderly'
  
  ageDat <- rbind(young,middle,elderly)
  ageDat$Age <- as.factor(ageDat$Age)
  
  agedeets <-aggregate(normSNR_CC ~ name + Age, ageDat, function(x) cbind(mean(x),sd(x),length(x)))
  agedeets$sem <- agedeets$normSNR_CC[,2]/sqrt(agedeets$normSNR_CC[,3])
  agedeets$Age <- factor(agedeets$Age,levels(agedeets$Age)[c(3,2,1)]) #reorder age gorups to look nicer
  
  if (save_stuff == 1){tiff(file.path(save_dir, "age_comparison_T2_SNR.tiff"), units="in", width=3,height=4.5,res=300)}
  
  p <- ggplot(agedeets, aes(fill=Age, y=normSNR_CC[,1], x=name)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=normSNR_CC[,1]-sem, ymax=normSNR_CC[,1]+sem), width=.2,size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    scale_fill_manual(values=c('grey70','grey35','black')) +
    ggtitle("SNR estimations of CC for\nT2w images across age groups") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text=element_text(size=7),
          axis.text.y=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"),
          legend.text = element_text(size=6))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ## mixed effects model to compare age groups SNR
  T2.age$name <- as.factor(T2.age$name)
  T2.age$Age <- as.integer(T2.age$Age)
  T2.age$tmp <- T2.age$name
  
  #### FREQUENTIST
  T2.CC.aov <- aov(normSNR_CC ~ Age + Error(name), data=T2.age) # aov CC
  T2.CC.model <- lmer(normSNR_CC ~ Age + (1|name), data=T2.age) # CC SNR model
  T2.CC.model1 <- lmer(normSNR_CC ~ Age + (1|name), data=T2.age, REML=F) # full model
  T2.CC.model2 <- lmer(normSNR_CC ~ (1|name), data=T2.age, REML=F) # null model
  
  anova(T2.CC.model1, T2.CC.model2) # not significant
  
  T2.CN.aov <- aov(normSNR_CN ~ Age + Error(name), data=T2.age) # aov CN
  T2.CN.model <- lmer(normSNR_CN ~ Age + (1|name), data=T2.age) # CC SNR model
  T2.CN.model1 <- lmer(normSNR_CN ~ Age + (1|name), data=T2.age, REML=F) # full model
  T2.CN.model2 <- lmer(normSNR_CN ~ (1|name), data=T2.age, REML=F) # null model
  
  anova(T2.CN.model1, T2.CN.model2) # signnificant
  
  T2.CNR.aov <- aov(normCNR ~ Age + Error(name), data=T2.age) # aov CNR
  T2.CNR.model <- lmer(normCNR ~ Age + (1|name), data=T2.age) # CC SNR model
  T2.CNR.model1 <- lmer(normCNR ~ Age + (1|name), data=T2.age, REML=F) # full model
  T2.CNR.model2 <- lmer(normCNR ~ (1|name), data=T2.age, REML=F) # null model
  
  anova(T2.CNR.model1, T2.CNR.model2) # not significant
  
  ### BAYESIAN
  # SNR of CC
  blm1 <- lmBF(normSNR_CC ~ Age*name, data=T2.age, whichRandom = "tmp" )
  blm2 <- lmBF(normSNR_CC ~ Age + name, data=T2.age, whichRandom= "tmp" )
  blm3 <- lmBF(normSNR_CC ~ Age, data=T2.age, whichRandom= "tmp" )
  blm4 <- lmBF(normSNR_CC ~ name, data=T2.age, whichRandom= "tmp" )
  
  blm2/blm1
  blm3/blm2
  blm2/blm4 # not significant 
  
  # SNR of CN
  blm1 <- lmBF(normSNR_CN ~ Age*name, data=T2.age, whichRandom = "tmp" )
  blm2 <- lmBF(normSNR_CN ~ Age + name, data=T2.age, whichRandom= "tmp" )
  blm3 <- lmBF(normSNR_CN ~ Age, data=T2.age, whichRandom= "tmp" )
  blm4 <- lmBF(normSNR_CN ~ name, data=T2.age, whichRandom= "tmp" )
  
  blm2/blm1
  blm3/blm2
  blm2/blm4 # significant 
  
  T2.age <- T2.age[!T2.age$name %in% c('IXI','CAMCAN'),]
  
  # CNR 
  blm1 <- lmBF(normCNR ~ Age*name, data=T2.age, whichRandom = "tmp" )
  blm2 <- lmBF(normCNR ~ Age + name, data=T2.age, whichRandom= "tmp" )
  blm3 <- lmBF(normCNR ~ Age, data=T2.age, whichRandom= "tmp" )
  blm4 <- lmBF(normCNR ~ name, data=T2.age, whichRandom= "tmp" )
  
  blm2/blm1
  blm3/blm2
  blm2/blm4 # not significant 

  #### visualize relationship of CNR with age
  ggplot(T2.age, aes(x=Age, y=normCNR, color=name)) + 
    geom_point()+
    geom_smooth(method=lm)
  
  ggplot(T2.age, aes(x=Age, y=normSNR_CC, color=name)) + 
    geom_point()+
    geom_smooth(method=lm)
  
  ggplot(T2.age, aes(x=Age, y=normSNR_CN, color=name)) + 
    geom_point()+
    geom_smooth(method=lm)
  
  ####################################################################################
  ################################################ COMPARE AGES - T2 data CNR
  ####################################################################################
  
  ageDat <- T2dat[T2dat$name%in%c('CAMCAN','IXI'),]
  
  young <- dplyr::filter(ageDat, grepl("YOU",Database_sub_ses))
  young$Age <- 'Young'
  middle <- dplyr::filter(ageDat, grepl("MID",Database_sub_ses))
  middle$Age <- 'Middle'
  elderly <- dplyr::filter(ageDat, grepl("ELD",Database_sub_ses))
  elderly$Age <- 'Elderly'
  
  ageDat <- rbind(young,middle,elderly)
  ageDat$Age <- as.factor(ageDat$Age)
  
  agedeets <-aggregate(normCNR ~ name + Age, ageDat, function(x) cbind(mean(x),sd(x),length(x)))
  agedeets$sem <- agedeets$normCNR[,2]/sqrt(agedeets$normCNR[,3])
  agedeets$Age <- factor(agedeets$Age,levels(agedeets$Age)[c(3,2,1)]) #reorder age gorups to look nicer
  
  if (save_stuff == 1){tiff(file.path(save_dir, "age_comparison_T2_CNR.tiff"), units="in", width=3,height=4.5,res=300)}
  
  p <- ggplot(agedeets, aes(fill=Age, y=normCNR[,1], x=name)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=normCNR[,1]-sem, ymax=normCNR[,1]+sem), width=.2,size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "CNR") +
    scale_fill_manual(values=c('grey70','grey35','black')) +
    ggtitle("CNR estimations for T2w\nimages across age groups") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text=element_text(size=7),
          axis.text.y=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ####################################################################################
  ################################################ COMPARE QUANTITATIVE T1 MAPS TO T1 WEIGHTED IMAGES
  ####################################################################################
  
  QT1 <- dat[dat$contrast=='qT1',] # MPICBS, MPILMBB, AHEAD, AHEAD (slab)
  compT1 <- dat[dat$contrast=='T1'&dat$name%in%c('MPICBS','MPILMBB','AHEAD','ATAG'),]
  compT1$contrast <- 'T1w'
  QT1 <- rbind(QT1, compT1)
  QT1$name[QT1$name=='AHEAD'] <- 'AHEADwb'
  QT1$name[QT1$name=='ATAG'] <- 'ATAGwb'
  QT1$name[QT1$dim==0.12500&QT1$name=='AHEADwb'] <- 'AHEADsb'
  QT1$name[QT1$dim==0.21600&QT1$name=='ATAGwb'] <- 'ATAGsb'
  
  QT1$normSNR_CC[QT1$contrast=='qT1'] <- QT1$normSNR_CC[QT1$contrast=='qT1']/QT1$dim[QT1$contrast=='qT1']
  
  qdeets <- aggregate(normSNR_CC ~ contrast + name, QT1, function(x) cbind(mean(x), sd(x), length(x)))
  qdeets$sem <- qdeets$normSNR_CC[,2]/sqrt(qdeets$normSNR_CC[,3])
  
  for (change in c('AHEADsb','AHEADwb','ATAGsb','ATAGwb','MPICBS')){
    qdeets$name[qdeets$name==change] <- paste0("'", change)
  }
  
  qdeets$name <- as.factor(qdeets$name)
  
  qdeets$name <- factor(qdeets$name, levels=c('MPILMBB',"'AHEADwb","'ATAGwb","'ATAGsb","'AHEADsb","'MPICBS"))
  
  #### if you want to extract SNR from CN run the follwing: ###
  QT1$normSNR_CN[QT1$contrast=='qT1'] <- QT1$normSNR_CN[QT1$contrast=='qT1']/QT1$dim[QT1$contrast=='qT1']
  CNdeets <- aggregate(normSNR_CN ~ contrast + name, QT1, function(x) cbind(mean(x), sd(x), length(x)))
  CNdeets$sem <- CNdeets$normSNR_CN[,2]/sqrt(CNdeets$normSNR_CN[,3])
  for (change in c('AHEADsb','AHEADwb','ATAGsb','ATAGwb','MPICBS')){
    CNdeets$name[CNdeets$name==change] <- paste0("'", change)
  }
  CNdeets$name <- as.factor(CNdeets$name)
  CNdeets$name <- factor(CNdeets$name, levels=c('MPILMBB',"'AHEADwb","'ATAGwb","'ATAGsb","'AHEADsb","'MPICBS"))
  
  if (save_stuff == 1){tiff(file.path(save_dir, "qT1_T1w_slab_whole.tiff"), units="in", width=5,height=3,res=300)}
  
  p <- ggplot(qdeets, aes(fill=contrast, y=normSNR_CC[,1], x=name)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=normSNR_CC[,1]-sem, ymax=normSNR_CC[,1]+sem), width=.2,size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    ggtitle("SNR estimations for T1 maps and T1 weighted images") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text=element_text(size=7),
          axis.text.y=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm")) +
    scale_fill_manual(values=c("grey30","grey70")) + 
    guides(fill=guide_legend(title="Contrast"))
  
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  
  ####################################################################################
  ################################################ COMPARE SLABS to WHOLE BRAIN IMAGES
  ####################################################################################
  
  slabby <- dat[dat$name%in%c('AHEAD','ATAG'),]
  
  # if you want to normalize the quantitative SNRs run this
  slabby$normSNR_CC <- slabby$normSNR_CC/slabby$dim
  slabby$normSNR_CN <- slabby$normSNR_CN/slabby$dim
  slabby$normCNR <- slabby$normCNR/slabby$dim
  
  sb <- dplyr::filter(slabby, grepl("SLAB",Label))
  wh <- dplyr::filter(slabby, grepl("WH",Label))
  
  sbinfo <- aggregate(cbind(normSNR_CC, normCNR, normSNR_CN) ~ name + contrast, sb, function(x) cbind(mean(x),sd(x),length(x)))
  whinfo <- aggregate(cbind(normSNR_CC, normCNR, normSNR_CN) ~ name + contrast, wh, function(x) cbind(mean(x),sd(x),length(x)))
  
  sbinfo$SNRsem <- sbinfo$normSNR_CC[,2]/sqrt(sbinfo$normSNR_CC[,3])
  sbinfo$CNRsem <- sbinfo$normCNR[,2]/sqrt(sbinfo$normCNR[,3])
  sbinfo$SNRCNsem <- sbinfo$normSNR_CN[,2]/sqrt(sbinfo$normSNR_CN[,3])
  
  whinfo$SNRsem <- whinfo$normSNR_CC[,2]/sqrt(whinfo$normSNR_CC[,3])
  whinfo$CNRsem <- whinfo$normCNR[,2]/sqrt(whinfo$normCNR[,3])
  whinfo$SNRCNsem <- whinfo$normSNR_CN[,2]/sqrt(whinfo$normSNR_CN[,3])
  
  t.test(sbinfo$normSNR_CC[,1], whinfo$normSNR_CC[c(1,3:7),1], paired=T) # 0.0017
  t.test(sbinfo$normCNR[,1], whinfo$normCNR[c(1,3:7),1], paired=T) # 0.07425
  
  ttestBF(sbinfo$normSNR_CC[,1], whinfo$normSNR_CC[c(1,3:7),1], paired=T) # 25.8
  ttestBF(sbinfo$normCNR[,1], whinfo$normCNR[c(1,3:7),1], paired=T) # 1.57
  
  ###########################################################################################
  ####################################### Age related differences in SNR of grey matter (CN) T1
  ###########################################################################################
  
  ageDat <- T1dat[T1dat$name%in%c('AHEAD','ATAG','CAMCAN','DLBS','IXI','MPILMBB','NKIRS','SALD'),]
  
  young <- dplyr::filter(ageDat, grepl("YOU",Database_sub_ses))
  young$Age <- 'Young'
  middle <- dplyr::filter(ageDat, grepl("MID",Database_sub_ses))
  middle$Age <- 'Middle'
  elderly <- dplyr::filter(ageDat, grepl("ELD",Database_sub_ses))
  elderly$Age <- 'Elderly'
  
  ageDat <- rbind(young,middle,elderly)
  ageDat$Age <- as.factor(ageDat$Age)
  
  ageDat$name[ageDat$name=='AHEAD'] <- "'AHEAD"
  ageDat$name[ageDat$name=='ATAG'] <- "'ATAG"
  
  agedeets <-aggregate(normSNR_CN ~ name + Age, ageDat, function(x) cbind(mean(x),sd(x),length(x)))
  agedeets$sem <- agedeets$normSNR_CN[,2]/sqrt(agedeets$normSNR_CN[,3])
  agedeets$Age <- factor(agedeets$Age,levels(agedeets$Age)[c(3,2,1)]) #reorder age gorups to look nicer
  
  if (save_stuff == 1){tiff(file.path(save_dir, "age_comparison_T1_SNR_CN.tiff"), units="in", width=5,height=3,res=300)}
  
  p <- ggplot(agedeets, aes(fill=Age, y=normSNR_CN[,1], x=name)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=normSNR_CN[,1]-sem, ymax=normSNR_CN[,1]+sem), width=.2,size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    scale_fill_manual(values=c('grey70','grey35','black')) +
    ggtitle("SNR estimations of CN for T1w images across age groups") +
    theme(text=element_text(family="Times New Roman"),
         axis.text=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
         legend.title=element_text(size=6.5, face="bold"),
         legend.key.size = unit(0.35, "cm"),
         legend.position="right",
         legend.margin = margin(0,0,0,0),
         legend.box.spacing = unit(0.2, "cm"),
         legend.text = element_text(size=6.5))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  

  
  ####################################################################################
  ################################################ COMPARE AGES - T2 data SNR for CN (grey matter)
  ####################################################################################
  
  ageDat <- T2dat[T2dat$name%in%c('CAMCAN','IXI'),]
  
  young <- dplyr::filter(ageDat, grepl("YOU",Database_sub_ses))
  young$Age <- 'Young'
  middle <- dplyr::filter(ageDat, grepl("MID",Database_sub_ses))
  middle$Age <- 'Middle'
  elderly <- dplyr::filter(ageDat, grepl("ELD",Database_sub_ses))
  elderly$Age <- 'Elderly'
  
  ageDat <- rbind(young,middle,elderly)
  ageDat$Age <- as.factor(ageDat$Age)
  
  agedeets <-aggregate(normSNR_CN ~ name + Age, ageDat, function(x) cbind(mean(x),sd(x),length(x)))
  agedeets$sem <- agedeets$normSNR_CN[,2]/sqrt(agedeets$normSNR_CN[,3])
  agedeets$Age <- factor(agedeets$Age,levels(agedeets$Age)[c(3,2,1)]) #reorder age gorups to look nicer
  
  if (save_stuff == 1){tiff(file.path(save_dir, "age_comparison_T2_SNR_CN.tiff"), units="in", width=3,height=4.5,res=300)}
  
  p <- ggplot(agedeets, aes(fill=Age, y=normSNR_CN[,1], x=name)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=normSNR_CN[,1]-sem, ymax=normSNR_CN[,1]+sem), width=.2,size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    scale_fill_manual(values=c('grey70','grey35','black')) +
    ggtitle("SNR estimations of CN for\nT2w images across age groups") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text=element_text(size=7),
          axis.text.y=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.title=element_text(size=6.5, face="bold"),
          legend.key.size = unit(0.35, "cm"),
          legend.position="right",
          legend.margin = margin(0,0,0,0),
          legend.box.spacing = unit(0.2, "cm"),
          legend.text = element_text(size=6.5))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  ## mixed effects model to compare age groups SNR
  ageDat$name <- as.factor(ageDat$name)
  ageDat$Age <- as.factor(ageDat$Age)
  library(lme4)
  # compare null model to full model 
  # full model here only includes age 
  # as name(database) is random, input as (1|name)
  # must put REML to false for information criterions
  null.model <- lmer(normSNR_CN ~ name + (1|name), data=ageDat, REML = F)
  full.model <- lmer(normSNR_CN ~ Age + name + (1|name), data=ageDat, REML = F)
  anova(null.model, full.model) # significant 

  age.model <- lmer(normSNR_CN ~ Age + name + (1|name), data=ageDat)
  
  bruh.model <- lmer(normSNR_CN ~ Age + (1|name), data=ageDat, REML = F)
  no.mod <- lmer(normSNR_CN ~ (1|name),data=ageDat, REML=F)
  anova(null.model, full.model, bruh.model, no.mod)
  
  
  #########################################################################
  ######################## Display raw SNRs and normalized SNRs
  #########################################################################
  
  deets <- aggregate(cbind(normSNR_CC, CC) ~ name, T1dat, function(x) cbind(mean(x),sd(x),length(x)))
  deets$NORMsem <- deets$normSNR_CC[,2]/sqrt(deets$normSNR_CC[,3])
  deets$RAWsem <- deets$CC[,2]/sqrt(deets$CC[,3])
  
  for (change in c('250','AHEAD','ATAG','MAASTRICHT','MPICBS')){
    deets$name[deets$name==change] <- paste0("'",change)
  }
  
  deets$Database <- factor(deets$name, levels=deets$name[order(deets$normSNR_CC[,1])])
  
  deetsframe <- data.frame('Database' = rep(deets$Database, 2), 'SNR' = c(deets$normSNR_CC[,1],deets$CC[,1]), 'Type' = rep(c('Normalized','Raw'),each=length(deets$normSNR_CC[,1])),
                           'SEM' = c(deets$NORMsem,deets$RAWsem))
  
  deetsframe$Database <- factor(deetsframe$Database, levels = deetsframe$Database[order(deetsframe$SNR[deetsframe$Type=='Raw'])])
  
  if (save_stuff == 1){tiff(file.path(save_dir, "supplementary_fig_1.tiff"), units="in", width=5,height=3,res=300)}
  
  p <- ggplot(deetsframe, aes(fill=Type, y=SNR, x=Database)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=SNR-SEM, ymax=SNR+SEM), width=.2,size=0.3,
                  position=position_dodge(.9)) +
    labs(x = "Database", y = "SNR") +
    ggtitle("Raw vs Normalized SNR") +
    theme(text=element_text(family="Times New Roman", size=8),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(angle=45, hjust=1, size =7),  
          axis.text=element_text(size=7),
          axis.title=element_text(size=7,face="bold"),
          plot.title = element_text(size=8,face="bold"),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5,face="bold")) +
          scale_fill_manual(values=c("grey30","grey70"))
  print(p)
  
  if (save_stuff ==1){dev.off()}   
  
  
  ##############################################################
  # Compare quantitative T1 parameters across age 
  
  # Databases: Ahead, MPILMBB, ATAG, MPI-CBS
  
  qcomp <- dat[dat$contrast=='qT1',]
  
  qcomp$Age <- c(68, 63, 67, 75, 62, 56, 60, 55, 58, 54, 25, 23, 26, 27, 23, #ATAG WH
                 74, 75, 73, 72, 63, 40, 45, 53, 55, 49, 25, 25, 22, 24, 25, # ATAG SLAB
                 31, 26, 25, 23, 22, # MPI-CBS
                 23, 23, 23, 23, 23, 43, 38, 48, 43, 48, 63, 68, 68, 68, 63, # MPI-LMBB
                 25, 21, 24, 25, 26, 51, 47, 40, 45, 49, 71, 71, 74, 77, 67, # AHEAD WH
                 25, 21, 24, 25, 26, 51, 47, 40, 45, 49, 71, 71, 74, 77, 67) # AHEAD SLAB
  
  qcomp$normSNR <- qcomp$CC/qcomp$dim
  qcomp$normCN <- qcomp$normSNR_CN
  qcomp$normCNR <- qcomp$normCNR
  
  plot(normSNR ~ Age, data=qcomp)
 
  # Set models
  # sort data
  qt1.check <- qcomp %>% dplyr::select(name, normSNR, normCN, normCNR, Age) # frame with Database, SNR of CC, SNR of CN, 
  qt1.check$Age <- as.integer(qt1.check$Age)
  qt1.check$name <- as.factor(qt1.check$name)
  
  #### FREQUENTIST
  CC.aov <- aov(normSNR ~ Age + Error(name), data=qt1.check) # aov CC
  CC.model <- lmer(normSNR ~ Age + (1|name), data=qt1.check) # CC SNR model
  CC.model1 <- lmer(normSNR ~ Age + (1|name), data=qt1.check, REML=F) # full model
  CC.model2 <- lmer(normSNR ~ (1|name), data=qt1.check, REML=F) # null model
  
  anova(CC.model1, CC.model2) # not significant
  
  CN.aov <- aov(normCN ~ Age + Error(name), data=qt1.check) # aov CN
  CN.model <- lmer(normCN ~ Age + (1|name), data=qt1.check) # CC SNR model
  CN.model1 <- lmer(normCN ~ Age + (1|name), data=qt1.check, REML=F) # full model
  CN.model2 <- lmer(normCN ~ (1|name), data=qt1.check, REML=F) # null model 
  
  anova(CN.model1, CN.model2) # signnificant
  
  CNR.aov <- aov(normCNR ~ Age + Error(name), data=qt1.check) # aov CNR
  CNR.model <- lmer(normCNR ~ Age + (1|name), data=qt1.check) # CC SNR model
  CNR.model1 <- lmer(normCNR ~ Age + (1|name), data=qt1.check, REML=F) # full model
  CNR.model2 <- lmer(normCNR ~ (1|name), data=qt1.check, REML=F) # null model
  
  anova(CNR.model1, CNR.model2) # not significant
  
  ### BAYESIAN

  
  
  
  
  
  
  
  
  
  