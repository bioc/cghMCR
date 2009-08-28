\name{cghMCR}
\alias{cghMCR}
\alias{segData}
\title{The constructor for the cghMCR class}
\description{
  Instantiates a cghMCR object using the parameters passed.
}
\usage{
cghMCR(segments, gapAllowed = 500, alteredLow = 0.03, alteredHigh =
0.97, spanLimit = 2e+07, recurrence = 75, thresholdType = c("quantile",
"value")) 
}
\arguments{
  \item{segments}{\code{segments} is a data frame extracted from the
    "output" element of the object returned by segment of the
    package DNAcopy or \code{\link{getSegments}}}
  \item{gapAllowed}{\code{gapAllowed} is an integer specifying low
    threshold of base pair number to separate two adjacent segments,
    belower which the two segments will be joined as an altered span}
  \item{alteredLow}{\code{alteredLow} is a positive number between 0 and
    1 specifying the lower reshold percential value. Only segments with
    values falling below this threshold are considered as altered span}
  \item{alteredHigh}{\code{alteredHigh} is a positive number between 0 and
    1 specifying the upper reshold percential value. Only segments with
    values falling over this threshold are considered as altered span}
  \item{recurrence}{\code{recurrence} is an integer between 1 and 100
    that specifies the rate of occurrence for a gain or loss that are
    observed across sample. Only gains/losses with ocurrence rate grater
    than the threshold values are declared as MCRs}
  \item{spanLimit}{\code{spanLimit} is an integer that defines the
    leangh of altered spans that can be considered as locus. It is not
    of any use at this time}
  \item{thresholdType}{\code{thresholdType} is a character string that
    can be either "quantile" or "value" indicating wether
    \code{alteredLow} or \code{alteredHigh} is quantial or actual value}
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
  data("segData")
  cghmcr <- cghMCR(segData, gapAllowed = 500,
         alteredLow = 0.20, alteredHigh = 0.80, recurrence = 50)
}
\keyword{methods}
