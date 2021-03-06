% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histMatchLength.R
\name{histMatchLength}
\alias{histMatchLength}
\title{Histogram of match length distribution of reads.}
\usage{
histMatchLength(aln, log10Transf = 0, titleHist)
}
\arguments{
\item{aln}{A GAlignments object of the BAM mapping file.}

\item{log10Transf}{A boolean. Either 0 (default) or 1 (log10).}

\item{titleHist}{a character. The main title for the histogram.
Default - none.}
}
\value{
A list with 2 elements.
The first element:
a data.frame of the number of counts per match length distribution.
The second element in the list:
a ggplot2 histogram of the match length distribution.
}
\description{
Histogram of match length distribution of reads.
}
\examples{
#starting from a GAlignment object
data(ctrlGAlignments)
aln <- ctrlGAlignments
#no log10 scaling
matchLenDistr <- histMatchLength(aln, 0)
#to plot the histogram
matchLenDistr[[2]]
}

