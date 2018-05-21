\name{redfitTauest}
\alias{redfitTauest}
\title{Tauest
}
\description{
The \code{\link{redfitTauest}} function is used by \code{\link{bin_cor}}
function to calculate the persistence for irregular climate time 
series under study. \code{\link{redfitTauest}} is included in the 
\emph{redfit} function that come from the \R \emph{dplR} package 
(Bunn et al. 2015). 
}
\usage{
redfitTauest(t, x)
}
\arguments{
  \item{t, x}{
t and x are the times and the variables for an unevenly spaced time series.
}
}
\details{
The \code{\link{redfitTauest}} function estimate the persistence of an 
irregular times series through the method of Mudelsee (2002). 
\code{\link{redfitTauest}} function is used by the \emph{dplR} package 
to estimate the persistence contained in irregular climate time 
series. More information about \code{\link{redfitTauest}} function can 
be found in Bunn et al. (2015) and Mudelsee (2002).
}
\note{
Needs \emph{dplR} to estimate the persistence contained in the irregular 
time series by means of the method of Mudelsee (2002). Please, look at the 
code \emph{tauest_dplR.R} in the directory \R of our BINCOR package. 
}
\references{
Bunn, A., Korpela, M., Biondi, F., Campelo, F., Mérian, P., Qeadan, F.,
Zang, C., Buras, A., Cecile, J., Mudelsee, M., Schulz, M. (2015). 
Dendrochronology Program Library in R. R package version 1.6.3.
URL \url{https://CRAN.R-project.org/package=dplR}. \cr 

Mudelsee, M. (2002). TAUEST: A computer program for estimating persistence 
in unevenly spaced weather/climate time series. Computers & Geosciences 28 
(1), 69--72. \cr
URL \url{http://www.climate-risk-analysis.com/software/}. \cr

Schulz, M., Mudelsee M. (2002). REDFIT: estimating red-noise spectra directly 
from unevenly spaced paleoclimatic time series. Computers & Geosciences 28(3), 
421--426. \cr
URL \url{https://www.marum.de/Michael-Schulz/Michael-Schulz-Software.html}. \cr

Mudelsee, M. (2010). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods. Springer. \cr

Mudelsee, M. (2014). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods, Second Edition. Springer. \cr 
}
\author{
Mikko Korpela .\cr
2013-2015 Aalto University, FINLAND.\cr 
Web: \url{https://github.com/mvkorpel}.\cr
Email: \url{mvkorpel@iki.fi} 
}
\keyword{Tauest}
\keyword{Memory}
\keyword{Persistence}
