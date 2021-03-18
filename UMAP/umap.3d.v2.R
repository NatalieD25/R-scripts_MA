## R script to calculate and visualize 3D UMAPs with variable vectors
## Script by Natalie Dallmann, Institute for environmental research, RWTH
## Version 2 (24.02.2021)

## 1. Setup
install.packages("rgl")
install.packages("uwot")
library("rgl")
library("uwot")
r3dDefaults$windowRect <- c(20,30, 800, 800)
#setwd(choose.dir()) 
data <- read.table("Example data.csv", sep = " ", header = T)

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

res.umap<-umap(data[,5:7], n_components = 3, scale = T)
  
## Save results
res.umap<-as.data.frame(res.umap)
rownames(res.umap)<-rownames(data)
res.umap<-cbind.data.frame(data[,1:8], res.umap)
filename<-paste(Sys.Date(), "_res.3d.umap.csv", sep = "")
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
  y<-res.umap[,11]
  x<-res.umap[,i]
  s.u3.m1<-smooth.spline(x,y)
  s.reg<-cbind(s.u1.m1$y, s.u2.m1$y, s.u3.m1$y)
  filename<-paste(Sys.Date(), ".3d.umap.spl.", name, ".csv", sep = "")
  write.csv(s.reg, file = filename)
}
filename<-paste(Sys.Date(), ".3d.umap.spl.", bm[1], ".csv", sep = "")
Ki67.spl<-read.csv(filename, sep = ",", header = T)
last1<-nrow(Ki67.spl)
filename<-paste(Sys.Date(), ".3d.umap.spl.", bm[2], ".csv", sep = "")
CD69.spl<-read.csv(filename, sep = ",", header = T)
last2<-nrow(CD69.spl)
filename<-paste(Sys.Date(), ".3d.umap.spl.", bm[3], ".csv", sep = "")
CD62L.spl<-read.csv(filename, sep = ",", header = T)
last3<-nrow(CD62L.spl)
  
## Plot results
titel = paste("Example UMAP") ##adjust title
plot3d(res.umap[,9:11], pch=20, cex = 2, col = c.col[as.factor(res.umap$celltype)], axes=FALSE,xlab='AX1',ylab='AX2',zlab='AX3', main = titel)
axes3d(col='black',box=TRUE,tick=FALSE,axes=FALSE)
lines3d(Ki67.spl[,2:4], lwd = 5, col = "gray0")
arrow3d(p0=c(Ki67.spl[last1-1,2:4]), p1=c(Ki67.spl[last1, 2:4]), s=1, col = "gray0", type = "r", theta = pi/12)
text3d(Ki67.spl[last1, 2:4], texts = "Ki67", offset = 1)
lines3d(CD69.spl[,2:4], lwd = 5, col = "gray40")
arrow3d(p0=c(CD69.spl[last2-1,2:4]), p1=c(CD69.spl[last2, 2:4]), s=1, col = "gray40", type = "r", theta = pi/12)
text3d(CD69.spl[last2, 2:4], texts = "CD69", offset = 1)
lines3d(CD62L.spl[,2:4], lwd = 5, col = "gray80")
arrow3d(p0=c(CD62L.spl[last3-1,2:4]), p1=c(CD62L.spl[last3, 2:4]), s=1, col = "gray80", type = "r", theta = pi/12)
text3d(CD62L.spl[last3, 2:4], texts = "CD62L", offset = 1)
legend3d("topleft", c("Ki67", "CD69", "CD62L", celltype), col = c("gray0", "gray40", "gray80", c.col), lty = c(1, 1, 1, NA, NA, NA, NA, NA), lwd = 3, pch = c(NA, NA, NA, 20, 20, 20, 20, 20))

## Save plot
name<-paste(Sys.Date(),"_3d.umap+var.vec.html", sep = "")
htmlwidgets::saveWidget(rglwidget(), name)
