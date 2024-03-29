% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakassignment.R
\name{check_nearby_peaks}
\alias{check_nearby_peaks}
\title{Close misses checker}
\usage{
check_nearby_peaks(df, rttable = NA, comparisontable = NA)
}
\arguments{
\item{df}{raw peak data, no manipulation}

\item{rttable}{this will be constructed from the raw data if not submitted
(should be equaivalent to calling the function yourself)}

\item{comparisontable}{If you want to call out a specific dataset, otherwise
it will build it from the raw submitted data equivalently to the standards method}
}
\value{
returns the anti-join table, which (mathematically) should be just peaks
that are not labelled in the dataset
}
\description{
Meant to check for peaks that are just outside of the window
}
