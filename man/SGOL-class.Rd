\name{SGOL-class}
\docType{class}
\alias{SGOL}
\alias{SGOL-class}
\alias{gol}
\alias{gol,SGOL-method}
\alias{method,SGOL-method}
\alias{plot}
\alias{plot,methods}
\alias{plot,ANY-method}
\alias{plot,SGOL-method}
\alias{plot,SGOL,ANY-method}
\alias{threshold,SGOL-method}

\title{Class "SGOL" reprents segments of gain or loss across samples}
\description{Segments of gains or losses along chromosomes can be
  calculated based on segmentation data derived from the segment
  function of the DNAcopy package}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SGOL", ...)} or
using a constructor SGOL.  
}
\section{Slots}{
	 \describe{
    \item{\code{gol}:}{Object of class \code{"matrix"} holding the gain
      or loss data for chromosomal segments }
    \item{\code{threshold}:}{Object of class \code{"vector"} of length 2
    indicating the lower and upper thresholds below/over which data
    points will be included in the calculation of SGOL score using the
    method defined by method}
    \item{\code{method}:}{Object of class \code{"function"} giving the
      method used to calculate SGOL scores. Common methods include sum,
      median, and mean}
  }
}
\section{Methods}{
  \describe{
    \item{gol}{\code{signature(object = "SGOL")}: extracts SGOL scores}
    \item{method}{\code{signature(object = "SGOL")}: gets the name of
      the function used to calculate the SGOL scores}
    \item{plot}{\code{signature(x = "SGOL", y = "ANY",...)}: plots the data}
    \item{threshold}{\code{signature(object = "SGOL")}: gets the
      thredhold used for the calculation}
	 }
}
\references{The SGOL score is a modified version of the GISTIC score
  published in PNAS 104: 20007-20012}
\author{Jianhua Zhang}
\note{}

\examples{
showClass("SGOL")
require(CNTools)

}
\keyword{classes}
