## Script for a principal response curve without statistics
## Requirements for datafile:
## - populations (type) coded as follows: 0 = T0rb (reference), 1 = T2rb, 2 = T0ina, 3 = T2ina, 4 = T2cur
## Script by Natalie Dallmann, Institute for environmental research, RWTH
## Version 2 (24.02.2021)

## Load required packages
library("vegan")

## set working directory
#setwd(choose.dir())

## load data
data<-read.table("pcr.data+dummy.csv", header = T, sep = ",")
data<-data[2:10] 
data1<-subset(data, ind != "dummy")

## prepare data for prc
treatment<-as.factor(data1$type2)
week<-as.factor(data1$week)
ID<-as.factor(data1$ID)
l.ID<-length(ID)
n.ID<-length(unique(ID))

##Loop to shuffel ID vector: to see if the implied "repeated measurements" bias the results
filename<-paste(Sys.Date(),"_prc_raw.cor.res.txt", sep = "") ##adjust filename
sink(file = filename, append = T)

  ## calculate prc w/statistics
  res.prc<-prc(response = data1[,5:8], treatment, week, correlation = T) ## for correlation matrix: correlation = T
  
  ## Save results
  print(res.prc)
  print(summary(res.prc))

  ## Plot results
  spec<-abs(res.prc$colsum)
  picname<-paste(Sys.Date(),"_prc_raw.cor.png", sep = "") ##adjust filename
  png(picname, width = 1000, height = 1000, res = 100)
  plot(res.prc, scaling = 0, lty = c(1, 2, 1, 1), col = c("royalblue4", "green", "green4", "red"), lwd = 2, legpos = NA)
  legend("topleft", legend = c("RB(ctrl)", "RB(high)", "ina(ctrl)", "ina(high)", "cur(high)"), lty = c(1, 1, 2, 1, 1), col = c("grey","royalblue4", "green", "green4", "red"), lwd = 2)
  dev.off()

sink()
