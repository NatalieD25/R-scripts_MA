## umap.subdatasets.v1
## R script to calculate and visualize 2D UMAPs with variable vectors (2.) and to perform a dbscan (3.) and manova (4.)
## on the umap axes 

## 1. Setup
library("uwot")
#setwd(choose.dir())
data <- read.table("Example data.csv", sep = "", header = T)

##Variables
bm<-c("Ki67", "CD69", "CD62L")
treatment<-c(paste(unique(data$treatment)))
mouse<-c(paste(unique(data$mouse)))
organ<-c(paste(unique(data$organ)))
celltype<-c(paste(unique(data$celltype)))

##Color selection of categorical variables
c.col<-c("green", "royalblue1", "red", "darkgreen", "royalblue4")

## 2. UMAP
## Calculation of UMAP with:
## n_components = number of axes
## scale = T -> normalize data by z-transformation
## y = supervision factor

res.umap<-umap(data[,5:7], n_components = 2, scale = T, y = data$celltype)
  
## Save results
res.umap<-as.data.frame(res.umap)
rownames(res.umap)<-rownames(data)
res.umap<-cbind.data.frame(data[,1:8], res.umap)
filename<-paste(Sys.Date(), "_res.2d.umap.csv", sep = "")
write.table(res.umap, file = filename, row.names = F)
  
## Calculate spline regressions for biomarkers
for (i in 5:7) {
  name<-bm[i-4]
  y<-res.umap[,9]
  x<-res.umap[,i]
  s.u1.m1<-smooth.spline(x,y)
  y<-res.umap[,10]
  x<-res.umap[,i]
  s.u2.m1<-smooth.spline(x,y)
  s.reg<-cbind(s.u1.m1$y, s.u2.m1$y)
  filename<-paste(Sys.Date(), ".2d.umap.spl.", name, ".csv", sep = "")
  write.csv(s.reg, file = filename)
}
filename<-paste(Sys.Date(), ".2d.umap.spl.", bm[1], ".csv", sep = "")
Ki67.spl<-read.csv(filename, sep = ",", header = T)
last1<-nrow(Ki67.spl)
filename<-paste(Sys.Date(), ".2d.umap.spl.", bm[2], ".csv", sep = "")
CD69.spl<-read.csv(filename, sep = ",", header = T)
last2<-nrow(CD69.spl)
filename<-paste(Sys.Date(), ".2d.umap.spl.", bm[3], ".csv", sep = "")
CD62L.spl<-read.csv(filename, sep = ",", header = T)
last3<-nrow(CD62L.spl)

## Plot results
picname<-paste(Sys.Date(), "_Example UMAP.png", sep = "")
png(picname)
titel = paste("Example UMAP") ##adjust title  
plot(res.umap[,9:10], pch=20, cex = 0.5, col = c.col[res.umap$celltype], xlab='AX1',ylab='AX2',zlab='AX3', main = titel)
lines(Ki67.spl[,2:3], lwd = 2, col = "gray0")
arrows(x0=Ki67.spl[last1-1,2], y0=Ki67.spl[last1-1,3], x1=Ki67.spl[last1,2], y1=Ki67.spl[last1,3], lwd = 2, col = "gray0")
text(Ki67.spl[last1, 2:3], labels = "Ki67", offset = 5)
lines(CD69.spl[,2:3], lwd = 2, col = "gray40")
arrows(x0=CD69.spl[last2-1,2], y0=CD69.spl[last2-1,3], x1=CD69.spl[last2,2], y1=CD69.spl[last2,3], lwd = 2, col = "gray40")
text(CD69.spl[last2, 2:3], labels = "CD69", offset = 5)
lines(CD62L.spl[,2:3], lwd = 2, col = "gray80")
arrows(x0=CD62L.spl[last3-1,2], y0=CD62L.spl[last3-1,3], x1=CD62L.spl[last3,2], y1=CD62L.spl[last3,3], lwd = 2, col = "gray80")
text(CD62L.spl[last3, 2:3], labels = "CD62L", offset = 5)
legend("topleft", c("Ki67", "CD69", "CD62L", celltype), col = c("gray0", "gray40", "gray80", c.col), lty = c(1, 1, 1, NA, NA, NA, NA, NA), lwd = 3, pch = c(NA, NA, NA, 20, 20, 20, 20, 20))
dev.off()
  
## Plot variable vector with heatmap
for (i in 1:3) {
  name<-bm[i] 
  val <- as.vector(res.umap[,i+4])
  picname<-paste(Sys.Date(), "_Heatmap+var.vec.", name, ".png", sep = "")
  png(picname)
  plot(res.umap[,9:10], pch=20, cex = 0.5, col = rainbow(length(val), start = 0, end = 0.33, rev = T)[rank(val)], main = name, ylab = "U2", xlab = "U1")
  lines(Ki67.spl[,2:3], lwd = 2, col = "gray0")
  arrows(x0=Ki67.spl[last1-1,2], y0=Ki67.spl[last1-1,3], x1=Ki67.spl[last1,2], y1=Ki67.spl[last1,3], lwd = 2, col = "gray0")
  text(Ki67.spl[last1, 2:3], labels = "Ki67", offset = 5)
  lines(CD69.spl[,2:3], lwd = 2, col = "gray0")
  arrows(x0=CD69.spl[last2-1,2], y0=CD69.spl[last2-1,3], x1=CD69.spl[last2,2], y1=CD69.spl[last2,3], lwd = 2, col = "gray0")
  text(CD69.spl[last2, 2:3], labels = "CD69", offset = 5)
  lines(CD62L.spl[,2:3], lwd = 2, col = "gray0")
  arrows(x0=CD62L.spl[last3-1,2], y0=CD62L.spl[last3-1,3], x1=CD62L.spl[last3,2], y1=CD62L.spl[last3,3], lwd = 2, col = "gray0")
  text(CD62L.spl[last3, 2:3], labels = "CD62L", offset = 5)
  dev.off()
}
