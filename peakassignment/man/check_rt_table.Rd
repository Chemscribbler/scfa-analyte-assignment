% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/referencetable.R
\name{check_rt_table}
\alias{check_rt_table}
\title{Error catching method, flags potetial problems with retention time tables}
\usage{
check_rt_table(df)
}
\arguments{
\item{df}{the retention table you are checking}
}
\value{
will return the original table if there are no errors
}
\description{
Three different conditions are checked: That the minimum is always increasing,
the maximumum is always inreasing, and that the maximum of an analyte is less
then the minmum of the following analyte
}
