\name{ccf_ts}
\alias{ccf_ts}
\alias{ccf_ts}
\title{Cross-correlation 
}
\description{
The \code{\link{ccf_ts}} function estimates and plots the cross-correlation 
between the binned time series. \code{\link{ccf_ts}} has an 
option to remove the linear trend of the time series under analysis (other 
pre-processing methods could be used) and contains several parameters that are 
described in the following lines. 
}
\usage{
ccf_ts(bints1, bints2, lagmax=NULL, ylima=-1, ylimb=1, rmltrd="N", RedL=T, 
 device="screen", Hfig, Wfig, Hpdf, Wpdf, resfig, ofilename)
}
\arguments{
  \item{bints1, bints2}{
The bints1 and bints2 are the binned time series.
}
  \item{lagmax}{
This parameter indicates the maximum lag for which the cross-correlation 
is calculated (its value depends on the length of the data set). 
}
  \item{ylima, ylimb}{
This parameters define the extremes of the range in which the CCF will 
be plotted. 
}
  \item{rmltrd}{
This is the option used to remove the linear trend in the time series 
under study (by default the linear trend is not removed, but it can 
be activated with the option ``Y'' or ``y''). 
}
  \item{RedL}{
RedL plots a right red line to highlight the correlation coefficient at 
the lag-0 (the default option is TRUE).
}
  \item{device}{
The type of the output device (by default the option is ``screen'', 
and the other options are ``jpg'', ``png'' and ``pdf'') for the 
scatter plot of the binned time series.
}
  \item{Hfig}{
The height for the CCF plot in ``jpg'' or ``png'' format.  
}
  \item{Wfig}{
The width for the CCF plot in ``jpg'' or ``png''format.  
}
  \item{Hpdf}{
The height for the CCF plot in ``pdf'' format.  
}
  \item{Wpdf}{
The width for the CCF plot in ``pdf'' format.  
}
  \item{resfig}{
\emph{resfig} is the plot resolution in ``ppi'' (by default \R does not 
record a resolution in the image file, except for BMP) for the CCF plot  
(``jpg'' or ``png'' formats), an adequate value could be 150 ppi. 
}
  \item{ofilename}{
The output filename (CCF plot) for the CCF estimated of the binned time series. 
}
}
\details{
The \code{\link{ccf_ts}} estimate the cross-correlation between two binned 
time series by means of the \R native function \emph{ccf} (package:stats).
}
\value{
Output: an object of the form \emph{ccf} containing the correlation 
coefficients for the defined number of lags (lagmax) and the statistical 
significance. 
}
\references{
Polanco-Martínez, J.M., Medina-Elizalde, M.A., Sánchez Goñi, M.F., 
M. Mudelsee. (2018). BINCOR: an R package to estimate the correlation 
between two unevenly spaced time series. Ms. under review (second round). 
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
 #::  Figure 5 (Polanco-Martínez et al. (2018), (mimeo)). 
 #####################################################################
 library("BINCOR") 
 library("pracma")

 #####################################################################
 #:: Loading the time series under analysis: example 2 (pollen ACER)  
 #####################################################################
 data(MD04_2845_siteID31)
 data(MD95_2039_siteID32)

 #####################################################################
 # Computing the binned time series though our bin_cor function
 #####################################################################
 bincor.tmp <- bin_cor(ID31.dat, ID32.dat, FLAGTAU=3, "salida_ACER_ABRUPT.tmp")
 binnedts   <- bincor.tmp$Binned_time_series

 # To avoid NA's values  
 bin_ts1    <- na.omit(bincor.tmp$Binned_time_series[,1:2])
 bin_ts2    <- na.omit(bincor.tmp$Binned_time_series[,c(1,3)]) 

 #####################################################################
 # Testing our ccf_ts function
 #####################################################################
 # Screen 
 ccf_ts(bin_ts1, bin_ts2, RedL=TRUE, rmltrd="y") 

 # PDF format  
 ccf_ts(bin_ts1, bin_ts2, RedL=TRUE, rmltrd="y", device="pdf", Hpdf=6, 
  Wpdf=9, resfig=300, ofilename="ccf_ID31_ID32_res")
 
 # JPG format 
 ccf_ts(bin_ts1, bin_ts2, RedL=TRUE, rmltrd="y", device="jpg", 
  Hfig=900, Wfig=1200, resfig=150, ofilename="ccf_ID31_ID32_res")
}
\keyword{cor_ts}
\keyword{Bi-variate correlation}
\keyword{Correlation}
