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
\section{Slots}{
	 \describe{
    \item{\code{DNASeg}:}{Object of class \code{"data.frame"} containing
      segmentation data derived from segmentation analysis using segment
      or \code{\link{getSegments}}}
    \item{\code{margin}:}{Object of class \code{"numeric"} indicating
      how far apart two adjacent segments are allowed to be considered
      to be of the same locus}
    \item{\code{gain.threshold}:}{Object of class \code{"numeric"}
      specificing the value of log2 ratio for chromosome gains}
    \item{\code{loss.threshold}:}{Object of class \code{"numeric"}
      specificing the value of negative log2 ratio for chromosome losses}
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
  require("cghMCR")
  require("marray")
  sampleData <- read.Agilent(list.files(file.path(.path.package("cghMCR"), 
                  "sampledata"), full.name = TRUE, pattern = "sample"), 
                  path = "")
  segments <- getSegments(sampleData)
  cghmcr <- cghMCR(segments, margin = 0, gain.threshold = 0.8,
                 loss.threshold = -0.8)
}
\keyword{classes}
