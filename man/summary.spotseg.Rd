\name{summary.spotseg}
\alias{summary.spotseg}
\title{Microarray Spot Segmentation Summary}
\description{Summary method for the \code{spotseg} function.
Gives the estimates of foreground and background intensity
obtained from microarray spot segmentation via model-based clustering.}
\usage{\method{summary}{spotseg}(object,...)}
\arguments{
\item{object}{An object of class \code{"spotseg"}, which is the
              output of the function \code{spotseg}.}
\item{...}{Unused, but required by generic \code{"summary"} method.}
}
\value{A list with two components, \code{"channel1"} and \code{"channel2"}
       each of which has subcomponents \code{"background"} and
       \code{"foreground"}, each of which in turn has subcomponents
       \code{"mean"} and \code{"median"}, giving the mean and median
        estimates of background and foreground for each channel.
       There will be missing entries (value \code{NA}) 
       whenever no foreground is detected.
}
\references{
Q. Li, C. Fraley, R. Bumgarner, K. Y. Yeung, and A. Raftery\\
Robust model-based segmentation of microarray images,\\ 
Technical Report No.~473, Department of Statistics, University of Washington,
January 2005.
}
\author{}
\seealso{\code{\link{spotseg}}}
\examples{
data(spotSegTest)

# columns of spotSegTest:
#  1 intensities from the Cy3 (green) channel
#  2 intensities from the Cy5 (red) channel

dataTransformation <- function(x) (256*256-1-x)^2*4.71542407E-05 

chan1 <- matrix(dataTransformation(spotSegTest[,1]), 144, 199)
chan2 <- matrix(dataTransformation(spotSegTest[,2]), 144, 199)

hivGrid <- spotgrid( chan1, chan2, rows = 4, cols = 6, show = TRUE)

library(mclust)

hivSeg <- spotseg( chan1, chan2, hivGrid$rowcut, hivGrid$colcut)

hivSummary <- summary(hivSeg)
}
\keyword{methods}
 

