\name{BINCOR-package}
\alias{BINCOR-package}
\alias{BINCOR}
\docType{package}
\title{Estimate the Correlation Between Two Irregular Time Series
}
\description{
'BINCOR' estimates the correlation between two irregular 
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
time series to exemplify its use. A description of the method and package 
is provided in \strong{Polanco-Martínez} et al. (2019), <URL: \doi{10.32614/RJ-2019-035}>. 
}
\details{
\tabular{ll}{
Package: \tab BINCOR\cr
Type: \tab Package\cr
Version: \tab 0.2.1\cr
Date: \tab 2026-02-23\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
BINCOR package contains four functions: the \code{\link{bin_cor}} (the 
main function to build the binned time series), the \code{\link{plot_ts}} 
(to plot and compare the irregular and binned time series), the 
\code{\link{cor_ts}} (to estimate the correlation between the binned 
time series) and the \code{\link{ccf_ts}} (to estimate the 
cross-correlation between the binned time series). 
}
\note{
 Dependencies: \emph{dplR} and \emph{pracma}. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
Dept. of Statistics, Fac. of Sciences, Unversity of Salamanca, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com} \cr\cr
Acknowledgement:\cr 
JMPM was funded by a Basque Government post-doctoral fellowship (2016-2018).  
}
\references{
Borchers, H. W. (2015). \emph{pracma: Practical Numerical Math Functions}. 
R package version 1.8.8, The Comprehensive R Archive Network (CRAN),
<URL: \url{https://CRAN.R-project.org/package=pracma}>. \cr

Bunn, A., Korpela, M., Biondi, F., Campelo, F., Mérian, P., Qeadan, F.,
Zang, C., Buras, A., Cecile, J., Mudelsee, M., Schulz, M. (2015). 
\emph{Dendrochronology Program Library in R}, R package version 1.6.3, 
The Comprehensive R Archive Network (CRAN), <URL: \url{https://CRAN.R-project.org/package=dplR}>. \cr

Mudelsee, M. (2010). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods. Springer. \cr

Mudelsee, M. (2014). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods, Second Edition. Springer. \cr 

\strong{Polanco-Martínez, J.M.}, Medina-Elizalde, M.A., Sánchez Goñi, M.F., 
M. Mudelsee. (2019). BINCOR: an R package to estimate the correlation 
between two unevenly spaced time series. \emph{R Journal}, \strong{11(1)}, 
170-184, <URL: \doi{doi:10.32614/RJ-2019-035}>. \cr
}
