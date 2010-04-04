\name{enrichreg}
\alias{enrichreg}
\title{Call and merge enriched probes to enriched regions.}
\description{
  A function used to call and merge enriched probes to enriched regions using the posterior
  probability calculated by iChip2 or iChip1 functions at certain
  posterior probability and false discovery rate (FDR) cutoffs. 
}
\usage{
enrichreg(pos,enrich,pp,cutoff,method=c("ppcut","fdrcut"),maxgap=500)
}

\arguments{
  \item{pos}{A n by 2 matrix or data frame. Rows correspond to probes. The first
    column of the matrix contains chromosome IDs; the second column
    contains the genomic positions.}
  \item{enrich}{A vector containing the probe enrichment measurements.}
  \item{pp}{A vector containing the posterior probabilities returned by
    iChip2 or iChip1. }
  \item{cutoff}{The cutoff value (a scalar) used to call enriched
    probes.  If use posterior probability as a criterion (method="ppcut"), a probe
    is said to be enriched if its pp is greater than the cutoff.  If use
    FDR as a criterion (method="fdrcut"), probes are said to be enriched
    if the probe-based FDR is less than the cutoff.  The FDR is
    calculated using a direct posterior probability approach (Newton et al., 2004).}
  \item{method}{'ppcut' or 'fdrcut'.}
  \item{maxgap}{The criterion used to merge enriched probes.  If the
    genomic distance of adjacent probes is less than maxgap,
    the probes will be merged into the same enriched regions. 
  }
}
\seealso{
  \code{\link{iChip2}}, \code{\link{iChip1}}, \code{\link{lmtstat}}
}
\value{
  A data frame with rows corresponding to enriched regions and columns
  corresponding to the following: 
  \item{chr}{Chromosome IDs.  For human genome, 23 and 24 denote X and Y, respectively.}
  \item{gstart}{The start genomic position of the enriched region.}
  \item{gend}{The end genomic position of the enriched region.}
  \item{rstart}{The row number for gstart in the position matrix.}
  \item{rend}{The row number for gend in the position matrix.}
  \item{peakpos}{The peak genomic position of the enriched region where
    the probe has the largest enrichment value.}
  \item{meanpp}{The mean posterior probability of the probes in the
    enriched region.}
  \item{maxpp}{The maximum posterior probability of the probes in the
    enriched region.}
  \item{nprobe}{The number of probes in the enriched regions.  nprobe =
    rend - rstart + 1}
}

\examples{
library(iChip)
library(limma)

#Analyze the p53 data (average resolution is about 35 bps)
#uncommenting the following code for running

#data(p53)
#p53lmt = lmtstat(p53[,9:14],p53[,3:8])
#p53Y = cbind(p53[,1],p53lmt)
#p53res=iChip2(Y=p53Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=2.5)
#enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,cutoff=0.9,
#          method="ppcut",maxgap=500)
#enrichreg(pos=p53[,1:2],enrich=p53lmt,pp=p53res$pp,cutoff=0.01,
#          method="fdrcut",maxgap=500)

}

\author{Qianxing Mo \email{moq@mskcc.org}}

\references{
  Qianxing Mo, Faming Liang.  (2010). Bayesian modeling of ChIP-chip data through
  a high-order Ising model. \emph{Biometrics}, 2010 Jan 29 [Epub ahead
  of print]. DOI: 10.1111/j.1541-0420.2009.01379.x

  Qianxing Mo, Faming Liang. (2010). A hidden Ising model for ChIP-chip
  data analysis.  \emph{Bioinformatics}  26(6), 777-783.
  doi:10.1093/bioinformatics/btq032
  
  Newton, M., Noueiry, A., Sarkar, D., Ahlquist, P. (2004). Detecting
  differential gene expression with a semiparametric hierarchical mixture method.
  \emph{Biostatistics}  5 , 155-176.
}
\keyword{models}

