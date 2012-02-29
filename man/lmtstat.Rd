\name{lmtstat}
\alias{lmtstat}
\title{A wrapper function used to calculated the limma t-statistics}
\description{
  A wrapper function used to calculated the empirical Bayes t-statistics
  (limma t-statistics) using functions in the limma package.
}
\usage{
lmtstat(IP,CON)
}
\arguments{
  \item{IP}{Data matrix for IP-enriched samples, where the rows and
    columns correspond to the probes and sample replicates,
    respectively.  The number of replicates must be greater than one.
    If CON is missing, IP is assumed to be in log-ratio format
    (e.g. log2(IP-enriched/control)).  In this case, paired t-statistics
    are calculated. If CON is NOT missing, IP and CON are assumed to
    be the normalized intensities for the IP-enriched and control
    samples,  respectively.  In this case, two-sample t-statistics are
    calculated.}
  \item{CON}{Data matrix for control samples, where the rows and
    columns correspond to the probes and sample replicates,
    respectively.  The number of replicates must be greater than one.}
}
\seealso{
  \code{\link{enrichreg}}, \code{\link{iChip2}},\code{\link{iChip1}}
}
\value{
  \item{ }{Empirical Bayes t-statistics calculated using functions in the
    limma package.}
}

\examples{

library(limma)

# load the log2 transformed and quantile-normalized Oct4 data 
data(oct4)
oct4[1:3,]

# calculate the enrichment measurements --- two-sample limma t-statistics
oct4lmt1 = lmtstat(oct4[,5:6],oct4[,3:4])

# calculate paired limma t-statistics for the data that are in
# the log-ratio format (e.g., log2(IP-enriched/control))

oct4log2r = oct4[,5:6] - oct4[,3:4]
oct4lmt2 = lmtstat(oct4log2r)

}

\author{Qianxing Mo \email{qmo@bcm.edu}}

\references{
  Smyth, G. K. (2004). Linear models and empirical Bayes methods for
  assessing differential expression in microarray experiments. Statistical
  Applications in Genetics and Molecular Biology 3, No. 1, Article 3. 
}
\keyword{models}

