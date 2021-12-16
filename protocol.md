``` r
# knitr::opts_chunk@set(message=FALSE,warning=FALSE,fig.path = "plots/")
rm(list=ls())
setwd("/Users/feiji/projects/Whetstine/STAR_protocol")
options(warn = -1)
# rmarkdown::render('protocol.R', output_format = 'all',clean=TRUE)
bin.bed = read.table("bin.50Kb.bed")
rownames(bin.bed) = paste("bin",1:dim(bin.bed)[1],sep="")
chr.list = unique(as.matrix(bin.bed[,1]))
bin_crd = rowMeans(bin.bed[,2:3]) # bin coordinate
```

# Overview

Load RT count in genomic bins, apply quantile normalization between
samples and LOESS smoothing.

``` r
# Load RT count #
load_CPM <- function(id){
  dat = read.table(sprintf("cnt.%s.50Kb.txt",id),header=T,check.names = F)
  colnames(dat) = gsub(".bam","",colnames(dat))
  colnames(dat) = gsub(".sort.rmdup","",colnames(dat))
  colnames(dat) = gsub("reptiming","",colnames(dat))
  cnt = as.matrix(dat[,-(1:3)]) # First 3 columns are bed coordinates, "Chr Start End", followed by columns for count number in each bin
  lib.size = colSums(cnt)/1e6
  cpm = t(t(cnt)/lib.size)
  return(cpm)
}
cpm.RT = load_CPM("RT")
```

### Quantile normalization

Perform quantile normalization on multiple samples of the same nature
(e.g. all S1 of Repli-seq)

``` r
gnm_qnm <- function(mtx){
  mtx.sort = matrix(NA,dim(mtx)[1],dim(mtx)[2])
  for(i in 1:dim(mtx)[2]){
    mtx.sort[,i] = sort(mtx[,i])
  }
  avg.cnt = rowMeans(mtx.sort)
  mtx.norm = mtx
  for(i in 1:dim(mtx)[2]){
    mtx.norm[,i] = avg.cnt[rank(mtx[,i])]
  }
  return(mtx.norm)
}
cpm = cpm.RT
cpm.norm = matrix(NA,dim(cpm)[1],dim(cpm)[2])
colnames(cpm.norm) = colnames(cpm)
rep = dim(cpm.norm)[2]/4
for(i in 1:4){
  cpm.norm[,(1:rep-1)*4+i] = gnm_qnm(cpm[,(1:rep-1)*4+i]) # same timing
}
```

### LOESS smoothing

Run Loess smoothing of the resulting normalized tracks

``` r
loess_smooth <- function(cpm.qnorm){
  cpm.qnorm.loess = cpm.qnorm #initialization 
  for(i in 1:dim(cpm.qnorm)[2]){
    # message(i)
    for(chr in chr.list){
      idx = which(bin.bed[,1]==chr)
      x = bin_crd[idx]
      y = cpm.qnorm[idx,i]
      lspan=500000/(max(bin_crd[idx])-min(bin_crd[idx])) #500Kb instead of 300Kb
      RPla=loess(y ~ x, span=lspan)
      cpm.qnorm.loess[idx,i] = RPla$fitted
    }
  }
  cpm.qnorm.loess[cpm.qnorm.loess<0]=0
  return(cpm.qnorm.loess)
}
RT.loess = loess_smooth(cpm.norm)
```

Example below shows the RT signal before and after LOESS normalization

``` r
idx = which(bin.bed[,1]=="chr1" & bin_crd>10e6 & bin_crd<15e6)
par(mar=c(4,4,2,2))
par(mgp=c(2,0.5,0)) # x-axis tick closer
plot(bin_crd[idx]/1e6,cpm.norm[idx,1],type="l",col="grey",xlab="",ylab="Normalized RT")
lines(bin_crd[idx]/1e6,RT.loess[idx,1],col="red")
legend(x="topright",legend = c("before LOESS","after LOESS"),lty=1,col=c("grey","red"),cex=0.7,border=NULL)
```

![](protocol_files/figure-gfm/RT.track.example.LOESS-1.png)<!-- -->
