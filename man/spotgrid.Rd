\name{spotgrid}
\alias{spotgrid}
\title{Gridding for Blocks of Microarray Spots}
\description{Determines row or column delimiters for spot locations
             from blocks of microarray slide image data.}
\usage{spotgrid(chan1, chan2, rows = NULL, cols = NULL, span = NULL,
                 show = FALSE)}
\arguments{
\item{chan1}{matrix of pixel intensities from the first channel.}
\item{chan2}{matrix of pixel intensities from the second channel.}
\item{rows}{number of spots in a row of the image block.}
\item{cols}{number of spots in a column of the image block.}
\item{span}{Window size for locating peak signals. This can be of length 2,
            in which case the first value is interpreted as a window size for
            the rows and the second as a window size for the columns.
            A default is estimated from the image dimension 
            and number of spots.}
\item{show}{logical variable indicating whether or not to display
            the gridding result.}
}
\details{}
\value{A list with two elements, \code{rowcut} and \code{colcut}
       giving delimiters for the row and/or column gridding
       of the slide. The indexes indicate the start of a segment
        of the grid, except for the last one, which indicates the
        end of the grid.}
\references{
Q. Li, C. Fraley, R. Bumgarner, K. Y. Yeung, and A. Raftery\\
Robust model-based segmentation of microarray images,\\ 
Technical Report No.~473, Department of Statistics, University of Washington,
January 2005.
}
\note{}
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

Grid <- spotgrid( chan1, chan2, rows = 4, cols = 6, show = TRUE)
}

\keyword{manip}
\keyword{robust}
