## dummydata.generator.v2.R
## Script to generate equally sized datasets for each timepoint and treatment
## Quality control of new dataset by t-test of old vs. new datasets: there should be no significant difference between the two datasets
## Populations ("type") should be coded as follows ("type2"): 0 = T0rb, 1 = T2rb, 2 = T0ina, 3 = T2ina, 4 = T2cur

##Setup: working directory, datafile
#setwd(choose.dir())
data<-read.table("rawdata.csv", sep = ",", header = T)
header<-colnames(data)
data1<-data

## Create vectors of populations and week
treatments<-as.vector(unique(data$type))
treatments2<-as.vector(unique(data$type2))
weeks<-as.vector(unique(data$week))

## Find largest dataset for treatment
for (h in 1:length(treatments2)) { ## for each celltype
  t<-subset(data, type2 == treatments2[h])
  maxobs=0
  for (i in weeks) { ## for each week
    obs<-sum(t$week==i)
    if (obs>maxobs) maxobs = obs
  }

  ## Create dummy dataset by sampling from a normal distribution with the corresponding mean and st.dev.
  for (k in weeks) { ## for each week
    t.wk<-subset(t, week == k)
    if (maxobs>nrow(t.wk)) {
      n<-maxobs-nrow(t.wk)
      dummy.wk.v<-matrix(NA, nrow = n)
      for (j in 1:4) { ## for each variable (biomarker)
        mean<-mean(t.wk[,j+3])
        std<-sd(t.wk[,j+3])
        dummy<-rnorm(n, mean = mean, sd = std)
        dummy.wk.v<-cbind.data.frame(dummy.wk.v, dummy)
      }
      
      ## put new dataset together
      wk.wk<-c(rep.int(k, n))
      mouse.wk<-as.vector(rep_len("dummy", n))
      pop<-paste(t.wk[1,3])
      ct.wk<-as.vector(rep_len(pop, n))
      ct2.wk<-as.vector(rep_len(treatments2[h], n))
      dummy.wk<-cbind.data.frame(wk.wk, mouse.wk, ct.wk, dummy.wk.v[,2:5], ct2.wk)
      colnames(dummy.wk)<-header
      data1<-rbind.data.frame(data1, dummy.wk)
    }
  }
}

## generate PRC-required IDs for "repeated measurement"
rd.data<-round(data1[,4:7])
data1<-cbind(data1[,1:3], rd.data, data1[8])
wk5<-subset(data1, week == 5)
ID<-as.data.frame(gl(nrow(wk5), 1, length = nrow(data1)))
rownames(ID)<-rownames(data1)
data1<-(data1[order(data1$type),])
data1<-(data1[order(data1$week),])
data1<-cbind.data.frame(ID, data1[,1:8])
header<-c("ID", "week", "ind", "type", "KLRG1", "CD25", "CCR8", "Ki67", "type2")
colnames(data1)<-header
tab.res<-paste(Sys.Date(), "_pcr.data+dummy.csv", sep = "")
write.csv(data1, file = tab.res, row.names=F)

## t-test of original and original+dummy data
t.res<-paste(Sys.Date(), "_res.t.test.raw+dummy.txt")
sink(t.res)
for (i in 1:length(treatments2)) {
  to<-subset(data, type2 == treatments2[i])
  tp<-subset(data1, type2 == treatments2[i])
  for (j in weeks) {
    t1<-subset(to, week == j)
    t2<-subset(tp, week == j)
    ndummies<-nrow(t2)-nrow(t1)
    dnum.nm<-paste("Number of dummies in week ", j, " and treatment ", treatments2[i], ": ", ndummies, sep = "")
    print(dnum.nm)
    if(nrow(t1)>3) {
    for (k in 1:4) {
      res.t.test<-t.test(t1[,k+3], t2[,k+4])
      print(paste(treatments2[i],".wk", j, ".var", k, sep = ""))
      print(paste("p-value = ", signif(res.t.test$p.value, digits = 2), sep = ""))
    }}
  }
}
sink()

