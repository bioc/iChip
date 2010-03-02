
library(iChip)
library(limma)

#oct4 and p53 data are log2 transformed and quantile-normalized intensities

# Analyze the Oct4 data (average resolution is about 280 bps) 
data(oct4)

# calculate the enrichment measurements --- the limma t-statistics
oct4lmt = lmtstat(oct4[,5:6],oct4[,3:4])

# prepare the data used for the Ising model 
oct4Y = cbind(oct4[,1],oct4lmt)

# Apply the second-order Ising model to the ChIP-chip data
oct4res = iChip2(Y=oct4Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=1.0)

# check the enriched regions detected by the Ising model using
# posterior probability (pp) cutoff at 0.9 or FDR cutoff at 0.01 
enrichreg(pos=oct4[,1:2],enrich=oct4lmt,pp=oct4res$pp,ppcut=0.9,isppcut=TRUE,maxgap=500)
enrichreg(pos=oct4[,1:2],enrich=oct4lmt,pp=oct4res$pp,fdrcut=0.01,isppcut=FALSE,maxgap=500)


# Analyze the p53 data (average resolution is about 35 bps)
data(p53)
p53lmt = lmtstat(p53[,9:14],p53[,3:8])
p53Y = cbind(p53[,1],p53lmt)
p53res = iChip2(Y=p53Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=2.5)

enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,ppcut=0.9,isppcut=TRUE,maxgap=500)
enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,fdrcut=0.01,isppcut=FALSE,maxgap=500)

