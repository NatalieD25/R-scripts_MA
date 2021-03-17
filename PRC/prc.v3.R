## Script for a principal response curve
## Requirements for datafile:
## - populations (type) coded as follows: 0 = T0rb, 1 = T2rb, 2 = T0ina, 3 = T2ina, 4 = T2cur
## - same amount of measurements in each week per treatment
## - ID-vector to imply "repeated measurement" for test of significance
## -> use Script "dummydata_generator" on rawdata

## Load required packages
install.packages("vegan")
install.packages("permute")
library("vegan")
library("permute")

## set working directory
#setwd(choose.dir())

## load data
data<-read.table("pcr.data+dummy.csv", header = T, sep = ",")
data<-data[2:10] 

## prepare data for prc
treatment<-as.factor(data$type2)
week<-as.factor(data$week)
ID<-as.factor(data$ID)
l.ID<-length(ID)
n.ID<-length(unique(ID))

##Loop to shuffel ID vector: to see if the implied "repeated measurements" bias the results
u.ID<-unique(ID)
filename<-paste(Sys.Date(),"_prc.cor.res.txt", sep = "") ##adjust filename
sink(file = filename, append = T)
for (i in 1:10){
  
  for (j in unique(week)){
    w.data<-subset(data, week == j)
    wo.data<-subset(data, week != j)
    for (k in 0:4){
      w.t.data<-subset(w.data, type2 == k)
      w.t.ID<-as.vector(w.t.data$ID)
      w.t.ID<-as.data.frame(sample(w.t.ID))
      rownames(w.t.ID)<-rownames(w.t.data)
      w.t.data<-cbind(w.t.ID, w.t.data[,2:9])
      colnames(w.t.data)<-colnames(data)
      o.data<-subset(w.data, type2 != k)
      w.data<-rbind(w.t.data, o.data)
    }
    data1<-rbind(w.data, wo.data)
  }
  
  treatment<-as.factor(data1$type2)
  week<-as.factor(data1$week)
  ID<-as.factor(data1$ID)
  
  ## calculate prc w/statistics
  res.prc<-prc(response = data1[,5:8], treatment, week, correlation = T) ## for correlation matrix: correlation = T
  ctrl <- how(plots = Plots(strata = ID,type = "free"), within = Within(type = "series"), nperm = 99)
  res.anova <- anova(res.prc, permutations=ctrl, first=TRUE)
  
  ## Save results
  print(paste("Permutation", i, sep = " "))
  print(res.prc)
  print(summary(res.prc))
  print(res.anova)
  print(summary(res.anova))
  
  ## Plot results
  spec<-abs(res.prc$colsum)
  picname<-paste(Sys.Date(),"_prc.cor.perm", i, ".png", sep = "") ##adjust filename
  png(picname, width = 1000, height = 1000, res = 100)
  plot(res.prc, scaling = 0, lty = c(1, 2, 1, 1), col = c("royalblue4", "green", "green4", "red"), lwd = 2, legpos = NA)
  legend("topleft", legend = c("RB(ctrl)", "RB(high)", "ina(ctrl)", "ina(high)", "cur(high)"), lty = c(1, 1, 2, 1, 1), col = c("grey","royalblue4", "green", "green4", "red"), lwd = 2)
  dev.off()
}

sink()
