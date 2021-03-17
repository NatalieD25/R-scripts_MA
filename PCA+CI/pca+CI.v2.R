##Script to analyze data via PCA with 95% confidence interval ellipses
##Script by Natalie Dallmann, Institute for environmental research, RWTH
##Version 2

##1. Setup
install.packages("FactoMineR")
install.packages("shape")
install.packages("vegan")
library("FactoMineR")
library("shape")
library("vegan")
#setwd(choose.dir())
data<-read.table("rawdata.csv", sep = ",", header = T)
cur<-as.vector(c("cur_5", "cur_6", "cur_7", "cur_8", "cur_9", "cur_11"))
ina<-as.vector(c("ina_5", "ina_6", "ina_7", "ina_8", "ina_9", "ina_11"))
rb<-as.vector(c("rb_5", "rb_6", "rb_7", "rb_8", "rb_9", "rb_11"))

##2. Create subset (eg. mice with same treatments, one mouse)
data1<-subset(data, ind != 46)

##3. Create environmental variable week+celltype from subset and add to subset
env.data<-matrix(NA,nrow = 1)
for (i in 1:nrow(data1)) {
  wk.type<-paste(data1[i,3], data1[i,1], sep = "_")
  env.data<-rbind(env.data, wk.type)
}
env.data<-as.data.frame(env.data[2:nrow(env.data)])
rownames(env.data)<-rownames(data1)
colnames(env.data)<-paste("wkty")
data2<-cbind(data1[,1:7], env.data)
wk.ty<-as.vector(unique(data2$wkty))

##4.1 PCA based on covariance matrix
## a) Calculate PCA with quali.sup = supplementary qualitative variable (wkty)
res.cov.pca<-PCA(data2[,4:8], scale.unit=F, ncp=12, quali.sup=5, graph = T)

## b) Get species and variable coordinates from PCA
var<-as.data.frame(res.cov.pca$var$coord)
spec<-cbind.data.frame(data2, res.cov.pca$ind$coord)

## Calculate ellipses and centroids for supplementary variable
CI<-coord.ellipse(spec[,8:12], bary = T)
CI.pt<-as.data.frame(CI$res)
env.mn<-matrix(NA, nrow=1, ncol=3)
mn.header<-c("env.var", "ax1", "ax2")
colnames(env.mn)<-mn.header
for (j in 1:length(wk.ty)){
  env.var<-subset(spec[,8:10], wkty == wk.ty[j])
  ax1<-mean(env.var[,2])
  ax2<-mean(env.var[,3])
  mn<-cbind.data.frame(wk.ty[j], ax1, ax2)
  colnames(mn)<-mn.header
  env.mn<-rbind(env.mn, mn)
  }
env.mn<-env.mn[2:nrow(env.mn),]

## c) Save PCA results
filename<-paste(Sys.Date(), ".pca.cov.T2.txt", sep = "") ## Adjust filename if necessary
sink(file = filename)
summary(res.cov.pca)
sink()
sscorename<-paste(Sys.Date(), ".pca.cov.sscores.T2.csv", sep = "") ## Adjust filename if necessary
write.table(spec, file=sscorename, sep = ",")
x<-paste("AX1 (", signif(res.cov.pca$eig[1,2], 4), " %)", sep = "")
y<-paste("AX2 (", signif(res.cov.pca$eig[2,2], 4), " %)", sep = "")

## d) Get mean and CI for each treatment
var.names<-c(rownames(var))
## cur
cur.pt<-matrix(NA, nrow = 1, ncol = 3)
colnames(cur.pt)<-colnames(env.mn)
cur.CI<-matrix(NA, nrow=1, ncol = 3)
colnames(cur.CI)<-colnames(CI.pt)
for (i in 1:length(cur)) {
  if (any(env.mn == cur[i])) {
  cp<-subset(env.mn, env.var == cur[i])
  cur.pt<-rbind.data.frame(cur.pt, cp)
  }
  if (any(CI.pt == cur[i])) {
  ce<-subset(CI.pt, CI.pt[,1] == cur[i])
  cur.CI<-rbind.data.frame(cur.CI, ce)
  }
}
cur.pt<-cur.pt[2:nrow(cur.pt),]
cur.CI<-cur.CI[2:nrow(cur.CI),]
## ina
ina.pt<-matrix(NA, nrow = 1, ncol = 3)
colnames(ina.pt)<-colnames(env.mn)
ina.CI<-matrix(NA, nrow=1, ncol = 3)
colnames(ina.CI)<-colnames(CI.pt)
for (i in 1:length(ina)) {
  if (any(env.mn == ina[i])) {
    ip<-subset(env.mn, env.var == ina[i])
    ina.pt<-rbind.data.frame(ina.pt, ip)
  }
  if (any(CI.pt == ina[i])) {
    ie<-subset(CI.pt, CI.pt[,1] == ina[i])
    ina.CI<-rbind.data.frame(ina.CI, ie)
  }
}
ina.pt<-ina.pt[2:nrow(ina.pt),]
ina.CI<-ina.CI[2:nrow(ina.CI),]
## rb
rb.pt<-matrix(NA, nrow = 1, ncol = 3)
colnames(rb.pt)<-colnames(env.mn)
rb.CI<-matrix(NA, nrow=1, ncol = 3)
colnames(rb.CI)<-colnames(CI.pt)
for (i in 1:length(rb)) {
  if (any(env.mn == rb[i])) {
    rp<-subset(env.mn, env.var == rb[i])
    rb.pt<-rbind.data.frame(rb.pt, rp)
  }
  if (any(CI.pt == rb[i])) {
    re<-subset(CI.pt, CI.pt[,1] == rb[i])
    rb.CI<-rbind.data.frame(rb.CI, re)
  }
}
rb.pt<-rb.pt[2:nrow(rb.pt),]
rb.CI<-rb.CI[2:nrow(rb.CI),]

##e) Plot with variable arrows, centroids for week+celltype and 95%-CI
picname<-paste(Sys.Date(), ".pca.cov.T2.png", sep = "") ## Adjust filename if necessary
png(picname, width = 1000, height = 1000, res = 200)
title<-paste("PCA of T2 mice (covariance matrix)") ## Adjust Plotname
plot(0, xlab = x, ylab = y, xlim = c(min(CI.pt$Dim.1), max(CI.pt$Dim.1)), ylim = c(min(CI.pt$Dim.2), max(CI.pt$Dim.2)), pch = NA, main = title)
points(cur.pt[,2], cur.pt[, 3], pch = 20, col = "red")
text(cur.pt[,2], cur.pt[, 3], labels = cur.pt[,1], col = "red", cex = 0.5, pos = 4, offset = 0.2)
for (i in 1:nrow(cur.pt)) {
  Arrows(cur.pt[i,]$ax1, cur.pt[i,]$ax2, cur.pt[i+1,]$ax1, cur.pt[i+1,]$ax2,  col = "red", code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
}
for (i in 1:length(cur)) {
  if (any(cur.CI == cur[i])) {
    ell<-subset(cur.CI, cur.CI$wkty == cur[i])
    lines(ell$Dim.1, ell$Dim.2, col = "red", lwd = 0.5)
  }
}
points(ina.pt[,2], ina.pt[, 3], pch = 20, col = "green3")
text(ina.pt[,2], ina.pt[, 3], labels = ina.pt[,1], col = "green3", cex = 0.5, pos = 4, offset = 0.2)
for (i in 1:nrow(ina.pt)) {
  Arrows(ina.pt[i,]$ax1, ina.pt[i,]$ax2, ina.pt[i+1,]$ax1, ina.pt[i+1,]$ax2,  col = "green3", code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
}
for (i in 1:length(ina)) {
  if (any(ina.CI == ina[i])) {
    ell<-subset(ina.CI, ina.CI$wkty == ina[i])
    lines(ell$Dim.1, ell$Dim.2, col = "green3", lwd = 0.5)
  }
}
points(rb.pt[,2], rb.pt[, 3], pch = 20, col = "blue")
text(rb.pt[,2], rb.pt[, 3], labels = rb.pt[,1], col = "blue", cex = 0.5, pos = 4, offset = 0.2)
for (i in 1:nrow(rb.pt)) {
  Arrows(rb.pt[i,]$ax1, rb.pt[i,]$ax2, rb.pt[i+1,]$ax1, rb.pt[i+1,]$ax2,  col = "blue", code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
}
for (i in 1:length(rb)) {
  if (any(rb.CI == rb[i])) {
    ell<-subset(rb.CI, rb.CI$wkty == rb[i])
    lines(ell$Dim.1, ell$Dim.2, col = "blue", lwd = 0.5)
  }
}
Arrows(0, 0, var$Dim.1, var$Dim.2, code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
for (i in 1:length(var.names)){
  text(var[i,]$Dim.1, var[i,]$Dim.2, labels = var.names[i], pos = 4, offset = 0.2, cex = 0.7)
}
legend("topleft", c("cur", "ina", "rb"), col = c("red", "green3", "blue"), pch = 20, lwd = 0.5, cex = 0.8)
dev.off()

##4.2 PCA based on correlation matrix
## a) Calculate PCA
res.cor.pca<-PCA(data2[,4:8], scale.unit=T, ncp=4, quali.sup=5, graph = T)

## b) Get species and variable coordinates from PCA
var<-as.data.frame(res.cor.pca$var$coord)
spec<-cbind.data.frame(data2, res.cor.pca$ind$coord)
CI<-coord.ellipse(spec[,8:12], bary = T)
CI.pt<-as.data.frame(CI$res)
env.mn<-matrix(NA, nrow=1, ncol=3)
mn.header<-c("env.var", "ax1", "ax2")
colnames(env.mn)<-mn.header
for (j in 1:length(wk.ty)){
  env.var<-subset(spec[,8:10], wkty == wk.ty[j])
  ax1<-mean(env.var[,2])
  ax2<-mean(env.var[,3])
  mn<-cbind.data.frame(wk.ty[j], ax1, ax2)
  colnames(mn)<-mn.header
  env.mn<-rbind(env.mn, mn)
}
env.mn<-env.mn[2:nrow(env.mn),]

## c) Save results
filename<-paste(Sys.Date(), ".pca.cor.T2.txt", sep = "") ## Adjust filename
sink(filename)
summary(res.cor.pca)
sink()
sscorename<-paste(Sys.Date(), ".pca.cor.sscores.T2.csv", sep = "") ## Adjust filename
write.table(spec, file=sscorename, sep = ",")
x<-paste("AX1 (", signif(res.cor.pca$eig[1,2], 4), " %)", sep = "")
y<-paste("AX2 (", signif(res.cor.pca$eig[2,2], 4), " %)", sep = "")

## d) Get mean and CI for each treatment
var.names<-c(rownames(var))
## cur
cur.pt<-matrix(NA, nrow = 1, ncol = 3)
colnames(cur.pt)<-colnames(env.mn)
cur.CI<-matrix(NA, nrow=1, ncol = 3)
colnames(cur.CI)<-colnames(CI.pt)
for (i in 1:length(cur)) {
  if (any(env.mn == cur[i])) {
    cp<-subset(env.mn, env.var == cur[i])
    cur.pt<-rbind.data.frame(cur.pt, cp)
  }
  if (any(CI.pt == cur[i])) {
    ce<-subset(CI.pt, CI.pt[,1] == cur[i])
    cur.CI<-rbind.data.frame(cur.CI, ce)
  }
}
cur.pt<-cur.pt[2:nrow(cur.pt),]
cur.CI<-cur.CI[2:nrow(cur.CI),]
## ina
ina.pt<-matrix(NA, nrow = 1, ncol = 3)
colnames(ina.pt)<-colnames(env.mn)
ina.CI<-matrix(NA, nrow=1, ncol = 3)
colnames(ina.CI)<-colnames(CI.pt)
for (i in 1:length(ina)) {
  if (any(env.mn == ina[i])) {
    ip<-subset(env.mn, env.var == ina[i])
    ina.pt<-rbind.data.frame(ina.pt, ip)
  }
  if (any(CI.pt == ina[i])) {
    ie<-subset(CI.pt, CI.pt[,1] == ina[i])
    ina.CI<-rbind.data.frame(ina.CI, ie)
  }
}
ina.pt<-ina.pt[2:nrow(ina.pt),]
ina.CI<-ina.CI[2:nrow(ina.CI),]
## rb
rb.pt<-matrix(NA, nrow = 1, ncol = 3)
colnames(rb.pt)<-colnames(env.mn)
rb.CI<-matrix(NA, nrow=1, ncol = 3)
colnames(rb.CI)<-colnames(CI.pt)
for (i in 1:length(rb)) {
  if (any(env.mn == rb[i])) {
    rp<-subset(env.mn, env.var == rb[i])
    rb.pt<-rbind.data.frame(rb.pt, rp)
  }
  if (any(CI.pt == rb[i])) {
    re<-subset(CI.pt, CI.pt[,1] == rb[i])
    rb.CI<-rbind.data.frame(rb.CI, re)
  }
}
rb.pt<-rb.pt[2:nrow(rb.pt),]
rb.CI<-rb.CI[2:nrow(rb.CI),]

## e) Plot
picname<-paste(Sys.Date(), ".pca.cor.T2.png", sep = "") ## Adjust filename
png(picname, width = 1000, height = 1000, res = 200)
plot(0, xlab = x, ylab = y, xlim = c(min(CI.pt$Dim.1), max(CI.pt$Dim.1)), ylim = c(min(CI.pt$Dim.2), max(CI.pt$Dim.2)), pch = NA, main = "PCA of T2 mice (correlation matrix)")
points(cur.pt[,2], cur.pt[, 3], pch = 20, col = "red")
text(cur.pt[,2], cur.pt[, 3], labels = cur.pt[,1], col = "red", cex = 0.5, pos = 4, offset = 0.2)
for (i in 1:nrow(cur.pt)) {
  Arrows(cur.pt[i,]$ax1, cur.pt[i,]$ax2, cur.pt[i+1,]$ax1, cur.pt[i+1,]$ax2,  col = "red", code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
}
for (i in 1:length(cur)) {
  if (any(cur.CI == cur[i])) {
    ell<-subset(cur.CI, cur.CI$wkty == cur[i])
    lines(ell$Dim.1, ell$Dim.2, col = "red", lwd = 0.5)
  }
}
points(ina.pt[,2], ina.pt[, 3], pch = 20, col = "green3")
text(ina.pt[,2], ina.pt[, 3], labels = ina.pt[,1], col = "green3", cex = 0.5, pos = 4, offset = 0.2)
for (i in 1:nrow(ina.pt)) {
  Arrows(ina.pt[i,]$ax1, ina.pt[i,]$ax2, ina.pt[i+1,]$ax1, ina.pt[i+1,]$ax2,  col = "green3", code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
}
for (i in 1:length(ina)) {
  if (any(ina.CI == ina[i])) {
    ell<-subset(ina.CI, ina.CI$wkty == ina[i])
    lines(ell$Dim.1, ell$Dim.2, col = "green3", lwd = 0.5)
  }
}
points(rb.pt[,2], rb.pt[, 3], pch = 20, col = "blue")
text(rb.pt[,2], rb.pt[, 3], labels = rb.pt[,1], col = "blue", cex = 0.5, pos = 4, offset = 0.2)
for (i in 1:nrow(rb.pt)) {
  Arrows(rb.pt[i,]$ax1, rb.pt[i,]$ax2, rb.pt[i+1,]$ax1, rb.pt[i+1,]$ax2,  col = "blue", code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
}
for (i in 1:length(rb)) {
  if (any(rb.CI == rb[i])) {
    ell<-subset(rb.CI, rb.CI$wkty == rb[i])
    lines(ell$Dim.1, ell$Dim.2, col = "blue", lwd = 0.5)
  }
}
Arrows(0, 0, var$Dim.1, var$Dim.2, code = 2, lwd = 1, arr.length = 0.2, arr.type = "triangle", arr.adj = 1)
for (i in 1:length(var.names)){
  text(var[i,]$Dim.1, var[i,]$Dim.2, labels = var.names[i], pos = 4, offset = 0.2, cex = 0.7)
}
legend("topleft", c("cur", "ina", "rb"), col = c("red", "green3", "blue"), pch = 20, lwd = 0.5, cex = 0.8)
dev.off()
