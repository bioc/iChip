\name{iChip1}
\alias{iChip1}
\title{Bayesian modeling of ChIP-chip data through hidden Ising models}
\description{
  Function iChip1 implements an algorithm of modeling ChIP-chip data
  through a standard hidden Ising model.
}
\usage{
iChip1(enrich,burnin=2000,sampling=10000,sdcut=2,beta0=3,
       minbeta=0,maxbeta=100,normsd=0.1,verbose=FALSE)
}
\arguments{
  \item{enrich}{A vector containing the probe enrichment
    measurements. The measurements must be sorted, firstly by chromosome and then by
    genomic position.  The measurements could be log2 ratios of
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
  \item{sdcut}{A value used to set the initial state for each
    probe. The enrichment measurements of a enriched probe is typically
    several standard deviations higher than the global mean enrichment
    measurements.}
  \item{beta0}{The initial parameter used to control the strength of interaction
    between probes, which must be a positive value.  A larger value of
    beta represents a stronger interaction between probes.  The value
    for beta0 could not be too small (e.g. < 1.0).  Otherwise, the Ising
    system may not be able to reach a super-paramagnetic state.}
  \item{minbeta}{The minimum value of beta allowed.}
  \item{maxbeta}{The maximum value of beta allowed.}
  \item{normsd}{iChip1 uses a Metropolis random walk proposal for
    sampling from the posterior distributions of the model
    parameters. The proposal distribution is a normal distribution with
    mean 0 and standard deviation specified by normsd.}
  \item{verbose}{A logical variable.  If TRUE, the number of completed MCMC
    iterations is reported.}
}
\seealso{
  \code{\link{iChip2}},\code{\link{enrichreg}}, \code{\link{lmtstat}}
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
  \item{lambda}{The posterior samples of the precision of the
    enrichment measurements of the probes.}
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

# Apply the standard Ising model to the ChIP-chip data

oct4res = iChip1(enrich=oct4lmt,burnin=2000,sampling=10000,sdcut=2,
                 beta0=3,minbeta=0,maxbeta=100,normsd=0.1)

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
# p53res = iChip1(p53lmt,burnin=2000,sampling=10000,sdcut=2,beta0=3,
#                 minbeta=0,maxbeta=100,normsd=0.1)

# enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,cutoff=0.9,
#           method="ppcut",maxgap=500)
# enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,cutoff=0.01,
#           method="fdrcut",maxgap=500)

}

\author{Qianxing Mo \email{moq@mskcc.org}}

\references{
  Qianxing Mo, Faming Liang. (2010). A hidden Ising model for ChIP-chip
  data analysis.  \emph{Bioinformatics}, 2010 Jan 28 [Epub ahead of print]. 
  doi:10.1093/bioinformatics/btq032
}
\keyword{models}

