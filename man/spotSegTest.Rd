\name{spotSegTest}
\docType{data}
\alias{spotSegTest}
\title{Spot Segmentation Test Data}
\description{
The two columns of this data set represent the Cy3 (green) absorption 
intensities for channel 1,
and the Cy5 (red) absorption intensities for channel 2
for part of a dye-swap experiment with replicates.
They measure expression levels of 
cellular RNA transcripts assessed in CD4+ T cell lines at different 
times after infection with HIV-1BRU using DNA microarrays.
}
\usage{data(spotSegTest)}
\format{Each column is a vector of intensities of 24 spots arranged in 
4 rows and 6 columns, encoded for compact (16-bit TIFF) storage.
For processing each column of \code{spotSegTest} 
should first be converted to a 144x199
matrix, then applying the transformation described below.
}
\details{
The intensities can be obtained from this data by first 
subtracting them from 65535 (256*256-1), then squaring, then multiplying by 
a scale factor 4.71542407E-05.
In other words, a number \code{x} in the \code{spotSegTest} data set 
corresponds to intensity  \preformatted{(256*256 - 1 - x)^2*.0000471542407}. \\
} 
\source{
Dr. Angelique van't Wout, Department of Microbiology, 
University of Washington\\
The data is a subset the first block of a 12 block array image
(\samp{001030\_08\_1.GEL})
in the first data set (\samp{2000095918 A})
in the first experiment (\samp{CEM LAI vs HI-LAI 24hr})
of the following data archive:\\
\url{http://expression.microslu.washington.edu/expression/vantwoutjvi2002.html}
}
\references{
{van't Wout} AB, Lehrman GK, Mikheeva SA, O'Keeffe GC, Katze MG, Bumgarner RE, 
Geiss GK, Mullins JI\\
Cellular gene expression upon human immunodeficiency virus type 1 infection of CD4(+)-T-cell lines.\\
\emph{J Virol.} 2003 Jan;77(2):1392-402. 
}
\keyword{datasets}

