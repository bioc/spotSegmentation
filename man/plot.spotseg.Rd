\name{plot.spotseg}
\alias{plot.spotseg}
\title{Microarray Spot Segmentation Plot}
\description{Plot method for the \code{spotseg} function.
Displays the result 
obtained from microarray spot segmentation via model-based clustering.}
\usage{\method{plot}{spotseg}(x,...)
}
\arguments{
\item{x}{An object of class \code{"spotseg"}, which is the
              output of the function \code{spotseg}.}
\item{...}{Unused but required by generic \code{"plot"} method.
}
}
\value{None, other than the displayed plot.
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

plot(hivSeg)
}
\keyword{methods}

