\name{plotBlockImage}
\alias{plotBlockImage}
\title{Plot Microarray Image Block}
\description{Displays a block of a microarray image.}
\usage{plotBlockImage(z,title,one)}
\arguments{
\item{z}{Intensities of the image pixels, in the form a of a matrix.}
\item{title}{A title for the image plot (optional).}
\item{one}{Sets appropriate graphics parameters for displaying
           individuals spots (default:FALSE).}
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

plotBlockImage(chan1)
plotBlockImage(chan2)
}
\keyword{methods}

