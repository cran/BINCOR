\name{cor_ts}
\alias{cor_ts}
\alias{cor_ts}
\title{Bi-variate correlation 
}
\description{
The \code{\link{cor_ts}} function estimates the correlation between the 
binned time series. \code{\link{cor_ts}} estimates three types of correlation 
coefficients: Pearson’s correlation, Spearman’s and Kendall’s rank correlations
by means of the \R native function cor.test (package:stats). The \code{\link{cor_ts}} function has an option to remove the linear trend of the time series under analysis (other pre-processing methods could be used) and its parameters are described in the following lines. 
}
\usage{
cor_ts(bints1, bints2, varnamets1="NULL", varnamets2="NULL", KoCM, rmltrd="N", 
 device="screen", Hfig, Wfig, Hpdf, Wpdf, resfig, ofilename)
}
\arguments{
  \item{bints1, bints2}{
The bints1 and bints2 are the binned time series.
}
  \item{varnamets1, varnamets2}{
varnamets[1][2] are the names of the variables under study.
}
  \item{KoCM}{
KoCM indicates the correlation estimator: pearson for Pearson (the
option by default), spearman for Spearman and kendall for Kendall. 
}
  \item{rmltrd}{
This is the option used to remove the linear trend in the time series 
under study (by default the linear trend is not removed, but it can 
be activated with the option ``Y'' or ``y''). 
}
  \item{device}{
The type of the output device (by default the option is ``screen'', 
and the other options are ``jpg'', ``png'' and ``pdf'') for the 
scatter plot for the binned time series.
}
  \item{Hfig}{
The height for the scatter plot in ``jpg'' or ``png'' format. 
}
  \item{Wfig}{
The width for the scatter plot in ``jpg'' or ``png'' format.  
}
  \item{Hpdf}{
The height for the scatter plot in ``pdf'' format. 
}
  \item{Wpdf}{
The width for the scatter plot in ``pdf'' format. 
}
  \item{resfig}{
\emph{resfig} is the resolution in ``ppi'' (by default \R does not 
record a resolution in the image file, except for BMP) for the scatter 
plot (``jpg'' or ``png'' formats), an adequate value could be 150 ppi. 
}
  \item{ofilename}{
The output filename for the scatter plot of the binned time series. 
}
}
\details{
The \code{\link{cor_ts}} estimate the correlation between two binned 
time series by means of the \R native function \emph{cor.test} (package:stats).
}
\value{
Output: an object of the form \emph{cor.test} containing the correlation 
coefficient and the statistical significance. \cr

Output plot: \emph{screen or 'ofilename + .png, .jpg or .pdf'}. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
Dept. of Statistics, Fac. of Sciences, Unversity of Salamanca, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com} 
}
\references{
Mudelsee, M. (2010). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods. Springer. \cr

Mudelsee, M. (2014). Climate Time Series Analysis: Classical Statistical and
Bootstrap Methods, Second Edition. Springer. \cr 

\strong{Polanco-Martínez, J.M.}, Medina-Elizalde, M.A., Sánchez Goñi, M.F., 
M. Mudelsee. (2019). \emph{BINCOR: an R package to estimate the correlation 
between two unevenly spaced time series}. \emph{R Journal}, \strong{11(1)}, 
170-184, <URL: \doi{doi:10.32614/RJ-2019-035}>. 
}
\examples{
\donttest{
 #####################################################################
 #::  Figure 2 (Polanco-Martínez et al. (2019) (mimeo)). 
 #####################################################################
 library("BINCOR")  
 library("pracma")

 #####################################################################
 #:: Loading the time series under analysis: example 1 (ENSO vs. NHSST) 
 #####################################################################
 data(ENSO) 
 data(NHSST)

 #####################################################################
 # Computing the binned time series though our bin_cor function 
 #####################################################################
 bincor.tmp    <- bin_cor(ENSO.dat, NHSST.dat, FLAGTAU=3, "output_ENSO_NHSST.tmp")
 binnedts      <- bincor.tmp$Binned_time_series

 #####################################################################
 # Testing our cor_ts function: cor_ts.R
 #####################################################################
 # screen (scatterplot) and Pearson  
 cor_ts(binnedts[,1:2], binnedts[,c(1,3)], "ENSO-Nino3", "SST NH Mean", 
  KoCM="pearson", rmltrd="y") 

 # PDF format (scatterplot) and Kendall 
 cor_ts(binnedts[,1:2], binnedts[,c(1,3)], "ENSO-Nino3", "SST NH Mean", 
  KoCM="kendall", rmltrd="y", device="pdf", Hpdf=6, Wpdf=9, resfig=300, 
  ofilename="scatterplot_ENSO_SST") 

 # JPG format (scatterplot) and Spearman  
 cor_ts( binnedts[,1:2], binnedts[,c(1,3)], "ENSO-Nino3", "SST NH Mean", 
  KoCM="spearman", rmltrd="y", device="jpg", Hfig=900, Wfig=1200, 
  resfig=150, ofilename="scatterplot_ENSO_SST") 
}
}
\keyword{cor_ts}
\keyword{Bi-variate correlation}
\keyword{Correlation}
