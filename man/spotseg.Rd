\name{spotseg}
\alias{spotseg}
\title{Microarray Spot Segmentation}
\description{Microarray spot segmentation via model-based clustering.}
\usage{spotseg(chan1, chan2, rowcut, colcut, R=NULL, C=NULL, 
                  threshold=100, hc=FALSE, show=FALSE)}
\arguments{
\item{chan1}{matrix of pixel intensities from the first channel.}
\item{chan2}{matrix of pixel intensities from the second channel.}
\item{rowcut}{row delimiters for the spots. Entries are the
              starting row location in the close of each spot, with
              the last entry being one pixel beyond the border of the
              last spot. For example, from the output of \code{spotgrid}.}
\item{colcut}{column delimiters for the spots. Entries are the
              starting column location in the close of each spot, with
              the last entry being one pixel beyond the border of the
              last spot. For example, from the output of \code{spotgrid}.}
\item{R}{rows over which the spots are to be segmented. The default
         is to segment spots in all rows.}
\item{C}{columns over which the spots are to be segmented. The default
         is to segment spots in all columns.}
\item{threshold}{connected components of size smaller than \code{threshold}
         are ignored. Default: \code{threshold=100}.}
\item{hc}{logical variable indicating whether or not EM should be
          initialized by hierarchical clustering or quantiles 
          in model-based clustering.
          The default is to use quantiles \code{hc = FALSE}, which
          is more efficient both in terms of speed and memory usage.}
\item{show}{logical variable indicating whether or not to display
            the segmentation of each individual spot as it is processed.
            The default is not to display the spots \code{show = FALSE}.}
}
\details{There are \code{plot} and \code{summary} methods that can
         be applied to the result.}
\value{An array of the same dimensions as the image in which the 
       pixels are labeled according to their group within the
       spot area: 1=background,2=uncertain,3=sample.}
\references{
Q. Li, C. Fraley, R. Bumgarner, K. Y. Yeung, and A. Raftery\\
Robust model-based segmentation of microarray images,\\ 
Technical Report No.~473, Department of Statistics, University of Washington,
January 2005.
}
\note{The \code{mclust} package is requiredfor clustering.}
\author{}
\seealso{\code{\link{summary.spotseg}}, \code{\link{plot.spotseg}},
         \code{\link{spotgrid}}}
\examples{
data(spotSegTest)

# columns of spotSegTest:
#  1 intensities from the Cy3 (green) channel
#  2 intensities from the Cy5 (red) channel

dataTransformation <- function(x) (256*256-1-x)^2*4.71542407E-05 

chan1 <- matrix(dataTransformation(spotSegTest[,1]), 144, 199)
chan2 <- matrix(dataTransformation(spotSegTest[,2]), 144, 199)

Grid <- spotgrid( chan1, chan2, rows = 4, cols = 6, show = TRUE)

library(mclust)

Seg <- spotseg( chan1, chan2, Grid$rowcut, Grid$colcut)

plot(Seg)

spotSummary <- summary(Seg)

spot11 <- spotseg( chan1, chan2, Grid$rowcut, Grid$colcut, 
                    R = 1, C = 1, show = TRUE)
}
\keyword{manip}
\keyword{robust}
\keyword{cluster}

