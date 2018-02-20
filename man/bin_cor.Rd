\name{bin_cor}
\alias{bin_cor}
\alias{bin_cor}
\title{Binned correlation 
}
\description{
The \code{\link{bin_cor}} function convert an irregular time series to 
a binned one and its parameters are described in the following lines. 
}
\usage{
bin_cor(ts1, ts2, FLAGTAU=3, ofilename)
}
\arguments{
  \item{ts1, ts2}{
ts1 and ts2 are the unevenly spaced time series.
}
  \item{FLAGTAU}{
FLAGTAU defines the method used to estimate the persistence or memory of 
the unevenly spaced time series. Options (by default is 3): \cr
 If 1 then est_tau = tau_x + tau_y     [Eq. 7.44, Mudelsee (2010, 2014)]. \cr
 If 2 then est_tau = max(tau_x, tau_y) [Eq. 7.45, Mudelsee (2010, 2014)]. \cr
 If 3 then est_tau = dist_x_y/ln(a_x_y_est) [Eq. 7.48, Mudelsee (2010, 2014).] \cr
 If 4 then est_tau via Carlstein [Eq. 3.28, Mudelsee (2010, 2014)]. 
}
  \item{ofilename}{
The output filename (ASCII format) containing the binned time series.
}
}
\details{
The \code{\link{bin_cor}} function convert an irregular times series to 
a binned time series and depends on the \R \emph{dplR} package to carry 
out this task. \emph{dplR} (\emph{redfitTauest} function) estimate the 
persistence contained in the irregular climate time series by means of 
the method of Mudelsee (2002). 
}
\note{
Needs \emph{dplR} (\emph{redfitTauest} function) to estimate the 
persistence contained in the irregular time series by means of the 
method of Mudelsee (2002). Please, look at the code \emph{tauest_dplR.R} 
in the directory \R of our BINCOR package. 
}

\value{
A list of 7 elements:
 \item{Binned_time_series}{An object containing the binned time series.} 
 \item{Auto._cor._coef._ts1}{The autocorrelation for the binned time series number 1.}
 \item{Persistence_ts1}{The persistence or memory for the binned time series number 1.}
 \item{Auto._cor._coef._ts2}{The autocorrelation for the binned time series number 2.}
 \item{Persistence_ts2}{The persistence or memory for the binned time series number 2.}
 \item{Mean_persistence}{The new times obtained for the binned time series.} 
 \item{Number_of_bins}{The number of bins.}
}
\references{
Bunn, A., Korpela, M., Biondi, F., Campelo, F., Mérian, P., Qeadan, F.,
Zang, C., Buras, A., Cecile, J., Mudelsee, M., Schulz, M. (2015). 
Dendrochronology Program Library in R. R package version 1.6.3.
URL \url{https://CRAN.R-project.org/package=dplR} \cr 

Mudelsee, M. (2002). TAUEST: A computer program for estimating persistence 
in unevenly spaced weather/climate time series. Computers & Geosciences 28 
(1), 69--72. \cr

Mudelsee, M. (2010). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods. Springer. \cr

Mudelsee, M. (2014). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods, Second Edition. Springer. \cr 

Polanco-Martínez, J.M., Medina-Elizalde, M.A., Sánchez Goñi, M.F., 
M. Mudelsee. (2018). BINCOR: an R package to estimate the correlation 
between two unevenly spaced climate time series. Ms. submitted. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, SPAIN.\cr 
EPOC UMR CNRS 5805 - U. de Bordeaux, Pessac, FRANCE.\cr
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{http://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \url{josue.m.polanco@gmail.com} 
}
\examples{
 #####################################################################
 #::  Figure 1 D (Polanco-Martínez et al. (2018), (mimeo)). 
 #####################################################################
 library("BINCOR") 

 #####################################################################
 #:: Loading the time series under analysis: example 1 (ENSO vs. NHSST) 
 #####################################################################
 data(ENSO) 
 data(NHSST)

 #####################################################################
 # Testing our bin_cor function
 #####################################################################
 bincor.tmp    <- bin_cor(ENSO.dat, NHSST.dat, FLAGTAU=3, "output_ENSO_NHSST.tmp")
 binnedts      <- bincor.tmp$Binned_time_series
}
\keyword{bincor}
\keyword{Binned correlation}
\keyword{bin_cor}
