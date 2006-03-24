\name{cghMCR}
\alias{cghMCR}
\alias{sampleData}
\title{The constructor for the cghMCR class}
\description{
  Instantiates a cghMCR object using the parameters passed.
}
\usage{
cghMCR(segments, margin = 0, gain.threshold = 0.5, loss.threshold = -0.5)
}
\arguments{
  \item{segments}{\code{segments} is a data frame extracted from the
    "output" element of the object returned by segment of the
    package DNAcopy or \code{\link{getSegments}}}
  \item{margin}{\code{margin} is a integer specifying how far apart two
    adjacent segments are allowed to be considered to be the same locus }
  \item{gain.threshold}{\code{gain.threshold} is a positive integer
    specifying the reshold log2 ratio when a segment of chromosome
    can be considered to be a gain}
  \item{loss.threshold}{\code{loss.threshold} is a negative integer
    specifying the reshold log2 ratio when a segment of chromosome can
    be considered to be a loss}
}
\details{
  The function is just a constructor of the \code{\link{cghMCR}} class
  for the instantiation of a cghMCR object
}
\value{
  An object of the cghMCR class
}
\references{References on S4 class}
\author{Jianhua Zhang}
\note{The function is a contribution of The Center for Applied Cancer
  Science of Dana-Farber Cancer Institute}
\seealso{\code{\link{cghMCR}}}
\examples{
  require("cghMCR")
  data("sampleData")
  segments <- getSegments(sampleData)
  cghmcr <- cghMCR(segments, margin = 0, gain.threshold = 0.8,
                 loss.threshold = -0.8)
}
\keyword{methods}
