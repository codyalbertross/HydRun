% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HydRun_functions.R
\name{match.precip.runoff}
\alias{match.precip.runoff}
\title{Match precipitation and runoff events}
\usage{
match.precip.runoff(PrecipEvents, RunoffEvents, n)
}
\arguments{
\item{PrecipEvents}{A list containing precipitation events, where each event is a data frame (format is 2 columns: datetime and precipitation).}

\item{RunoffEvents}{A list containing runoff events, where each event is a data frame (format is 2 columns: datetime and stormflow).}

\item{n}{The search window. The left edge of the window is n hours before the start of the runoff event.}
}
\value{
A list containing matched precipitation and runoff events.
}
\description{
\code{match.precip.runoff} is used to match precipitation-runoff events.
}
\note{
This is an R implementation of a function from the MATLAB toolbox HydRun.
}
\examples{
events <- match.precip.runoff(PrecipEvents, RunoffEvents, n = 96)
}
\references{
Tang, W., & Carey, S. K. (2017). HydRun: a MATLAB toolbox for rainfall runoff analysis. Hydrological Processes, 31(15), 2670-2682.

Fuka, D. R., Walter, M. T., Archibald, J. A., Steenhuis, T. S., Easton, Z. M., Fuka, M. D., & KeepSource, T. R. U. E. (2014). Package EcoHydRology.
}
