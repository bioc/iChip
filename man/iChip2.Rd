\name{iChip2}
\alias{iChip2}
\title{Bayesian modeling of ChIP-chip data through hidden Ising models}
\description{
  Function iChip2 implements a method of modeling ChIP-chip data through a high-order hidden Ising model. 
}
\usage{
iChip2(Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=2.5,verbose=FALSE)
}
\arguments{
  \item{Y}{A n by 2 matrix or data frame. The first column of Y contains
    the chromosome IDs; the second column of Y contains the probe enrichment
    measurements.  Y must be sorted, firstly by chromosome and then by
    genomic position.  The probe enrichment measurements could be log2 ratios of
    the intensities of IP-enriched and control samples for a single
    replicate, or summary statistics such as t-like statistics or mean
    differences for multiple replicates.  We suggest to use the
    empirical Bayesian t-statistics implemented in the limma package for
    multiple replicates.  Note, binding probes must have
    a larger mean value than non-binding probes.}
  \item{burnin}{The number of MCMC burn-in iterations.}
  \item{sampling}{The number of MCMC sampling iterations.  The
    posterior probability of binding and non-binding state is calculated
    based on the samples generated in the sampling period. }
  \item{winsize}{The parameter to control the order of interactions
    between probes.  For example, winsize = 2, means that probe i
    interacts with probes i-2,i-1,i+1 and i+2. A balance between high
    sensitivity and low FDR could be achieved by setting winsize = 2.}
  \item{sdcut}{A value used to set the initial state for each
    probe. The enrichment measurements of a enriched probe is typically
    several standard deviations higher than the global mean enrichment
    measurements.}
  \item{beta}{The parameter used to control the strength of interaction
    between probes, which must be a positive value.  A larger value of
    beta represents a stronger interaction between probes.  In general,
    high resolution array such as Affymetrix tiling arrays have relatively
    stronger probe interactions than low resolution array such as Agilent
    tiling arrays. For the second order Ising model (winsize = 2), the
    critical value of beta is around 1.0.  For low resolution array data
    (e.g. 280 bp resolution), beta could be set to close to the critical
    value; For high resolution array data (e.g. 35 bp resolution), beta
    could be set to a value between 2 to 4.  In general, choosing a
    large value of beta amounts to using a more stringent criterion for
    detecting enriched regions in ChIP-chip experiments.
  }
  \item{verbose}{A logical variable.  If TRUE, the number of completed MCMC
    iterations is reported.}
}
\seealso{
  \code{\link{iChip1}}, \code{\link{enrichreg}}, \code{\link{lmtstat}}
}
\value{
  A list with the following elements.
  \item{pp}{The posterior probabilities of probes in
    the binding/enriched state.  There is a strong evidence to be a
    binding/enriched probe if the probe has a posterior probability close to 1. }
  \item{mu0}{The posterior samples of the mean measurement of the probes in
    the non-binding/non-enriched state.}
  \item{mu1}{The posterior samples of the mean measurement of the probes in
    the binding/enriched state.}
  \item{lambda0}{The posterior samples of the precision of the
    enrichment measurements of the probes in the non-binding/non-enriched
    state.}
  \item{lambda1}{The posterior samples of the precision of the
    enrichment measurements of the probes in the binding/enriched
    state.}
}

\examples{
library(iChip)
library(limma)

# oct4 and p53 data are log2 transformed and quantile-normalized intensities

# Analyze the Oct4 data (average resolution is about 280 bps)

data(oct4)

### sort oct4 data, first by chromosome then by genomic position
oct4 = oct4[order(oct4[,1],oct4[,2]),]

# calculate the enrichment measurements --- the limma t-statistics

oct4lmt = lmtstat(oct4[,5:6],oct4[,3:4])

# prepare the data used for the Ising model

oct4Y = cbind(oct4[,1],oct4lmt)

# Apply the second-order Ising model to the ChIP-chip data

oct4res=iChip2(Y=oct4Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=1.25)

# check the enriched regions detected by the Ising model using
# posterior probability (pp) cutoff at 0.9 or FDR cutoff at 0.01

enrichreg(pos=oct4[,1:2],enrich=oct4lmt,pp=oct4res$pp,cutoff=0.9,
          method="ppcut",maxgap=500)
enrichreg(pos=oct4[,1:2],enrich=oct4lmt,pp=oct4res$pp,cutoff=0.01,
          method="fdrcut",maxgap=500)


# Analyze the p53 data (average resolution is about 35 bps)
# uncommenting the following code for running

# data(p53)
# must sort the data first
# p53 = p53[order(p53[,1],p53[,2]),]
# p53lmt = lmtstat(p53[,9:14],p53[,3:8])
# p53Y = cbind(p53[,1],p53lmt)
# p53res=iChip2(Y=p53Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=2.5)

# enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,cutoff=0.9,
#           method="ppcut",maxgap=500)
# enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,cutoff=0.01,
#           method="fdrcut",maxgap=500)

}

\author{Qianxing Mo \email{moq@mskcc.org}}

\references{
  Qianxing Mo, Faming Liang. (2010). Bayesian modeling of ChIP-chip data through
  a high-order Ising model. \emph{Biometrics}, 2010 Jan 29 [Epub ahead
  of print]. DOI: 10.1111/j.1541-0420.2009.01379.x
}
\keyword{models}

