\name{cghMCR-class}
\docType{class}
\alias{cghMCR-class}
\alias{MCR}
\alias{getSegments}
\alias{MCR,cghMCR-method}
\alias{getSegments,marrayRaw-method}
\alias{getSegments,marrayNorm-method}

\title{Class "cghMCR" is a S4 class for the identification of minimum
  common regions of gains or losses across samples}
\description{Objects of this class provides the functionalities to
  detecting chromosome regions that show gains or losses across differnet
  samples}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("cghMCR", ...)}.
  A constructor \code{\link{cghMCR}} may be used to instantiate object
  of this class
}
\section{Slots}{\describe{
  \item{\code{DNASeg}:}{Object of class \code{"data.frame"} containing
    segmentation data derived from segmentation analysis using segment
    or \code{\link{getSegments}}}
  \item{\code{gapAllowed}:}{Object of class \code{gapAllowed} is an
    integer specifying low threshold of base pair number to separate two
    adjacent segments, belower which the two segments will be joined as
    an altered span} 
  \item{\code{alteredLow}:}{Object of class \code{alteredLow} is a
    positive number between 0 and 1 specifying the lower reshold
    percential value. Only segments with values falling below this
    threshold are considered as altered span} 
  \item{\code{alteredHigh}:}{Object of class \code{alteredHigh} is a
    positive number between 0 and 1 specifying the upper reshold
    percential value. Only segments with values falling over this
    threshold are considered as altered span} 
  \item{\code{recurrence}:}{Object of class \code{recurrence} is an
    integer between 1 and 100 that specifies the rate of occurrence for
    a gain or loss that are observed across sample. Only gains or losses
    with ocurrence rate grater than the threshold values are declared as MCRs}
  \item{\code{spanLimit}:}{Object of class \code{spanLimit} is an
    integer that defines the leangh of altered spans that can be
    considered as locus. It is not of any use at this time}
  }
}

\section{Methods}{
  \describe{
    \item{MCR}{\code{signature(object = "cghMCR")}: identifies minimum
      common regions of gains/losses across samples}
	 }
}
\references{}
\author{Jianhua Zhang}
\note{The function is a contribution of The Center for Applied Cancer
  Science of Dana-Farber Cancer Institute}

\seealso{
  \code{\link{cghMCR}}
}
\examples{
  require("CNTools")
  data("sampleData")
  cghmcr <- cghMCR(sampleData[sampleData[, "ID"] \%in\%
         sample(unique(sampleData[, "ID"]), 20), ], gapAllowed = 500,
         alteredLow = 0.20, alteredHigh = 0.80, recurrence = 50)
}
\keyword{classes}
