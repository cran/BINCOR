\name{BINCOR-package}
\alias{BINCOR-package}
\alias{BINCOR}
\docType{package}
\title{Estimate the Correlation Between Two Irregular Time Series
}
\description{
'BINCOR' estimate the correlation between two irregular 
time series that are not necessarily sampled on identical time points. 
This program is also applicable to the situation of two evenly spaced 
time series that are not on the same time grid. 'BINCOR' is based on 
a novel estimation approach proposed by Mudelsee (2010, 2014) to 
estimate the correlation between two climate time series with different 
timescales. The idea is that autocorrelation (AR1 process) allows to 
correlate values obtained on different time points. The outputs (plots) 
can be displayed in the screen or can be saved as PNG, JPG or PDF formats. 
The 'BINCOR' package also provides two examples with real data: instrumental 
(\code{\link{ENSO.dat}} and \code{\link{NHSST.dat}} data sets) and 
paleoclimatic (\code{\link{ID31.dat}} and \code{\link{ID32.dat}} data sets) 
time series to exemplify its use. 
}
\details{
\tabular{ll}{
Package: \tab BINCOR\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2018-05-18\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
BINCOR package contains four functions: the \code{\link{bin_cor}} (the 
main function to build the binned time series), the \code{\link{plot_ts}} 
(to plot and compare the irregular and binned time series, the 
\code{\link{cor_ts}} (to estimate the correlation between the binned 
time series) and the \code{\link{ccf_ts}} (to estimate the 
cross-correlation between the binned time series). 
}
\note{
 Dependencies: \emph{dplR} and \emph{pracma}. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, SPAIN.\cr 
EPOC UMR CNRS 5805 - U. de Bordeaux, Pessac, FRANCE.\cr
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{http://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \url{josue.m.polanco@gmail.com} \cr\cr
Acknowledgement:\cr 
JMPM was funded by a Basque Government post-doctoral fellowship.  
}
\references{
Borchers, H. W. (2015). pracma: Practical Numerical Math Functions. R
package version 1.8.8.
URL \url{https://CRAN.R-project.org/package=pracma} \cr

Bunn, A., Korpela, M., Biondi, F., Campelo, F., Mérian, P., Qeadan, F.,
Zang, C., Buras, A., Cecile, J., Mudelsee, M., Schulz, M. (2015). Den-
drochronology Program Library in R. R package version 1.6.3.
URL \url{https://CRAN.R-project.org/package=dplR} \cr

Mudelsee, M. (2010). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods. Springer. \cr

Mudelsee, M. (2014). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods, Second Edition. Springer. \cr 

Polanco-Martínez, J.M., Medina-Elizalde, M.A., Sánchez Goñi, M.F., 
M. Mudelsee. (2018). BINCOR: an R package to estimate the correlation 
between two unevenly spaced time series. Ms. under review (second round). \cr
}
