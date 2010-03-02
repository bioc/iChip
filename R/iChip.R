## The latest R functions for iChip software
# Quincy Mo, moq@mskcc.org
# Department of Epidemiology and Biostatistics,Memorial Sloan-Kettering Cancer center
# New York, NY 10065

# high-order Ising model for ChIP-chip data
## Y[,1] = chromosome,Y[,2]=enrichment measurements
iChip2 = function(Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=2.5,verbose=FALSE){
  if(missing(Y)){
    stop("Argument Y is missing!\n")
  }
  if(!(is.matrix(Y) || is.data.frame(Y))){
    stop("Y must be a matrix or data frame!\n")
  }
  if((burnin < 1) || (sampling < 1)){
    stop("burin or sampling is too small!\n")
  }
  if(winsize < 1){
    stop("winsize must be >= 1\n")
  }
  if(sdcut < 2){
    warning("sdcut may be too small")
  }
  if(beta < 0){
    stop("beta must be greater than 0! \n")
  }
  burn = as.integer(burnin)
  Size = as.integer(sampling)
  rowdim = as.integer(nrow(Y))
  chr = as.integer(as.factor(Y[,1]))
  idata = as.double(Y[,2]) # exact genomic position is not needed
  post = as.double(rep(0,rowdim))
  halfw = as.integer(winsize)
  X = as.integer(rep(0,rowdim))
  sdcut = as.double(sdcut)
  bt = as.double(beta)
  mu0 = as.double(rep(0,(burnin+sampling)))
  mu1 = mu0
  lambda0 = mu0
  lambda1 = mu0
  verb = as.integer(0)
  if(verbose == TRUE){verb = as.integer(1)}
  res = .C("iChip2",burn,Size,rowdim,chr,idata,halfw,sdcut,bt,postX=post,X=X,
    mu0=mu0,lambda0=lambda0,mu1=mu1,lambda1=lambda1,verbose=verb,PACKAGE="iChip")
  list(pp=res$postX,mu0=res$mu0,lambda0=res$lambda0,mu1=res$mu1,lambda1=res$lambda1)
}

## standard one-dimensional Ising model for ChIP-chip data 
iChip1 = function(enrich,burnin=2000,sampling=10000,sdcut=2,beta0=3,minbeta=0,maxbeta=100,normsd=0.1,verbose=FALSE){
  if(missing(enrich)){
    stop("Argument enrich is missing!\n")
  }
  if((burnin < 1) || (sampling < 1)){
    stop("burin or sampling is too small!\n")
  }
  if(beta0 < 2){
    warning("beta0 may be too small. You may set beta0 between 2 and 4.\n")
  }
  
  if(minbeta < 0){
    stop("minbeta must be greater than 0!\n")
  }
  if(minbeta >= maxbeta){
    stop("minbeta must be less than maxbeta! \n")
  }
  if(sdcut < 2){
    warning("sdcut may be too small. \n")
  }
  if(normsd < 0){
    warning("normsd must be greater than 0. \n")
  }
  burn = as.integer(burnin)
  Size = as.integer(sampling)
  rowdim = as.integer(length(enrich))
  idata = as.double(enrich)
  sdcut = as.double(sdcut)
  btStart = as.double(beta0)
  minbt = as.double(minbeta)
  maxbt = as.double(maxbeta)
  ransd = as.double(normsd)
  post = as.double(rep(0,rowdim))
  X = as.integer(rep(0,rowdim))
  mu0 = as.double(rep(0,(burnin+sampling)))
  mu1 = mu0
  lambda = mu0
  pBeta = mu0
  verb = as.integer(0)
  if(verbose == TRUE){verb = as.integer(1)}
  res = .C("iChip1",burn,Size,rowdim,idata,sdcut,btStart,minbt,maxbt,ransd,postX=post,
    X=X,pBeta=pBeta,mu0=mu0,mu1=mu1,lambda=lambda,verb,PACKAGE="iChip")
  list(pp=res$postX,beta=res$pBeta,mu0=res$mu0,mu1=res$mu1,lambda=res$lambda)
}

#A wrapper function that calls the functions in limma packages for calculating the limma t-statistics
lmtstat = function(IP,CON){
#  require(limma)
  lmtstat = NULL
  if(missing(CON)){
    if(ncol(IP) < 2){
      stop("ERROR: IP must has at least two replicates.\n")
    }
    limfit = lmFit(cbind(IP), design=rep(1,ncol(IP)))
    ebtfit = eBayes(limfit)
    lmtstat = ebtfit$t[,1]
  }else{
    n1 = ncol(CON)
    n2 = ncol(IP)
    if(n1 < 2 || n2 < 2){
      stop("ERROR: IP and CON must have at least two replicates.\n ")
    }
    designX = cbind(int=rep(1,n1+n2),a=c(rep(0,n1),rep(1,n2)))
    limfit = lmFit(cbind(CON,IP), design=designX)
    ebtfit = eBayes(limfit)
    lmtstat = ebtfit$t[,2]
  }
  lmtstat
}

#function to find the binding regions based on posterior probability and selected cutoff
#pos[,1] is the chromosome, pos[,2] is the genomic position
#enrich is enrichment measurements
#pp is the posterior probability, cutoff is selection criteria for pp
#anno and pp are the results of the whole genome. 
#probes are merged if the genomic distance between neighboring probes are less than the maxgap
enrichreg = function(pos,enrich,pp,cutoff,method=c("ppcut","fdrcut"),maxgap=500){
  if((nrow(pos) != length(pp)) || (nrow(pos) != length(enrich))){
    stop("nrow(pos), length(enrich) and length(pp) must be equal. \n")
  }
  if(cutoff <= 0 || cutoff >=1){
    stop("Error: pp should be in region (0, 1).\n")
  }
  if(missing(cutoff)){
    stop("Error: cutoff must be specified!")
  }
  if(maxgap < 1){
    stop("maxgap must be greater than 1! \n")
  }
  
  type = match.arg(method)
  ppcut = FALSE
  fdrcut = FALSE
  if(type == "ppcut"){ppcut = TRUE}
  else if(type == "fdrcut"){fdrcut = TRUE}
  
  chrf = as.factor(pos[,1])
  chrn = as.integer(chrf)
  anno = cbind(chr=chrn,pos=pos[,2],rowID=1:nrow(pos))
  id = NULL
  if(ppcut){
    id = (pp > cutoff)
  }else if(fdrcut){
    sig = dirFDR(pp,cutoff)
    id = (sig$id > 0)
  }else{
    stop("Error: arg for method must be one of 'ppcut','fdrcut'.")
  }
  nsig = sum(id)
  if(nsig == 0){
    cat("- No enriched region found - \n")
    br=NA
    return(br)
  }
  x = NULL
  if(sum(id) != 1){
    x = cbind(anno[id,],pp[id])
  }else{
    x = matrix(c(anno[id,],pp[id]),nrow=1)
  }
  xlen = nrow(x) #length
  ############ chrom   start   end   row.s  row.e  max.pp mean.pp
  y = matrix(c(x[1,1], x[1,2], x[1,2],x[1,3],x[1,3]),nrow=1)
  if(xlen == 1){
    br = cbind(y,x[1,2],x[1,4],x[1,4],1)
    colnames(br) = c("chr","gstart","gend","rstart","rend","peakpos","meanpp","maxpp","nprobe")
    br = as.data.frame(br)
    br[1,1] = pos[br[1,4],1] #use the original chromosome label
    return(br)
  }

  Y = matrix(0,nrow=xlen,ncol=5)
  Y[1,] = y[1,]
  res = .C("MergeRegion",as.integer(as.matrix(x[,1:3])),as.integer(xlen),as.integer(3),
    as.integer(5), as.integer(maxgap),Y=as.integer(t(Y)),nregion=as.integer(0),PACKAGE="iChip")

  y = matrix(res$Y,ncol=5)
  y = y[1:res$nregion,]
  
  if(res$nregion==1){
    y = matrix(y,ncol=5)
  }

#  maxID = .C("maxmvID",as.double(enrich),as.integer(nrow(y)),as.integer(y[,4]-1),as.integer(y[,5]-1),
#    maxid=as.integer(rep(0,nrow(y))),PACKAGE="iChip")
#  maxID = maxID$maxid + 1  #note c idex from 0

  maxID = function(y,enrich){
    (y[4]:y[5])[which.max(enrich[y[4]:y[5]])]
  }
  meanFun = function(y,pp){
    mean(pp[y[4]:y[5]])
  }
  maxFun = function(y,pp){
    max(pp[y[4]:y[5]])
  }
  jpp = matrix(NA,nrow=nrow(y),ncol=4)
  jpp[,1] = pos[apply(y,1,maxID,enrich=enrich),2]
  jpp[,2] = round(apply(y,1,meanFun,pp=pp),2)
  jpp[,3] = round(apply(y,1,maxFun,pp=pp),2)
  jpp[,4] =  y[,5]-y[,4]+1 
  br = cbind(y,jpp)
  br = as.data.frame(br)
  br[,1] = pos[br[,4],1] #use the original chromosome label
  colnames(br) = c("chr","gstart","gend","rstart","rend","peakpos","meanpp","maxpp","nprobe")
  return(br)
}

dirFDR <- function(pp, fdrcut){
  beta.g = 1-pp
  k = sort(unique(beta.g))
  res = .C("fdr",as.integer(length(k)),as.double(k),as.integer(length(beta.g)),
    as.double(beta.g),efdr=as.double(rep(0,length(k))),PACKAGE="iChip")
  cutoff = k[sum(res$efdr<=fdrcut)] #find the biggest k
  pcut = 1 - cutoff  ## this pcut should be >= pcut
  id = rep(0,length(pp))
  id[beta.g <= cutoff] = 1
  return(list(id=id,pcut=pcut))
}
