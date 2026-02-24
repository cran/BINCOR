\name{plot_ts}
\alias{plot_ts}
\alias{plot_ts}
\title{Plot time series 
}
\description{
The \code{\link{plot_ts}} function plot and compare the irregular 
and the binned time series. \code{\link{plot_ts}} has 
several parameters that are described in the following lines. 
}
\usage{
plot_ts(ts1, ts2, bints1, bints2, varnamets1="", varnamets2="", colts1=1, colts2=1, 
 colbints1=2, colbints2=2, ltyts1=1, ltyts2=1, ltybints1=2, ltybints2=2, 
 device="screen", Hfig, Wfig, Hpdf, Wpdf, resfig, ofilename)
}
\arguments{
  \item{ts1, ts2}{
ts1 and ts2 are the unevenly spaced time series.
}
  \item{bints1, bints2}{
The bints1 and bints2 are the binned time series.
}
  \item{varnamets1, varnamets2}{
varnamets[1][2] are the names of the variables under study.
}
  \item{colts1, colts2}{
colts[1][2] are the colours for the time series (irregular) under 
 study (by default both curves are in black). 
}
  \item{colbints1, colbints2}{
colbints[1][2] are the colours of the binned time series (by default
both curves are in red). 
}
  \item{ltyts1, ltyts2}{
ltyts[1][2] are the type of lines to be plotted for the irregular time series 
(by default is 1, i.e., solid). 1 = solid, 2 = dashed, 3 = dotted, 4 = dot-dashed, 
5 = long-dashed, 6 = double-dashed. 
}
  \item{ltybints1, ltybints2}{
ltybints[1][2] are the type of lines to be plotted for the binned time series
(by default is 2, i.e., dashed).  1 = solid, 2 = dashed, 3 = dotted, 4 = dot-dashed, 
5 = long-dashed, 6 = double-dashed.
}
  \item{device}{
The type of the output device (by default the option is ``screen'', 
 and the other options are ``jpg'', ``png'' and ``pdf'').
}
  \item{Hfig}{
The height for the plot in ``jpg'' or ``png'' format. 

}
  \item{Wfig}{
The width for the plot in ``jpg'' or ``png'' format. 
}
  \item{Hpdf}{
The height for the plot in ``pdf'' format. 
}
  \item{Wpdf}{
The width for the plot in ``pdf'' format. 
}
  \item{resfig}{
\emph{resfig} is the plot resolution in 'ppi' (by default \R does not 
record a resolution in the image file, except for BMP), an adequate value 
could be 150 ppi. 
}
  \item{ofilename}{
The output filename for the plot. 
}
}
\details{
The \code{\link{plot_ts}} function is used to plot the irregular vs. 
the binned time series and this function uses the native R function ``plot''
 (package:graphics).
}
\value{
Output:

Output plot: \emph{screen or 'ofilename + .png, .jpg or .pdf'}. 
}
\author{
Josué M. Polanco-Martínez (a.k.a. jomopo).\cr
Dept. of Statistics, Fac. of Sciences, Unversity of Salamanca, SPAIN.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: \email{josue.m.polanco@gmail.com} \cr
}
\references{
\strong{Polanco-Martínez, J.M.}, Medina-Elizalde, M.A., Sánchez Goñi, M.F., 
M. Mudelsee. (2019). BINCOR: an R package to estimate the correlation 
between two unevenly spaced time series. \emph{R Journal}, \strong{11(1)}, 
170-184, <URL: \doi{doi:10.32614/RJ-2019-035}>.
}
\examples{
\donttest{
 #####################################################################
 #::  Figure 1 (Polanco-Martínez et al. (2019) (mimeo)). 
 #####################################################################
 library("BINCOR") 

 #####################################################################
 #:: Loading the time series under analysis: example 1 (ENSO vs. NHSST) 
 #####################################################################
 data(ENSO) 
 data(NHSST)

 #####################################################################
 # Computing the binned time series though our bin_cor_function.R 
 #####################################################################
 bincor.tmp    <- bin_cor(ENSO.dat, NHSST.dat, FLAGTAU=3, "output_ENSO_NHSST.tmp")
 binnedts      <- bincor.tmp$Binned_time_series

 #####################################################################
 # Testing our plot_ts function
 #####################################################################
 # "Screen" 
 plot_ts(ENSO.dat, NHSST.dat, binnedts[,1:2], binnedts[,c(1,3)], "ENSO-Nino3", 
  "SST NH Mean", colts1=1, colts2=2, colbints1=3, colbints2=4, device="screen")
 
 # PDF format
 plot_ts(ENSO.dat, NHSST.dat, binnedts[,1:2], binnedts[,c(1,3)], "ENSO-Nino3", 
  "SST NH Mean", colts1=1, colts2=2, colbints1=3, colbints2=4, device="pdf", 
  Hpdf=6, Wpdf=9, resfig=300, ofilename="plot_ts_RAW_BIN_enso_sst") 

 # PNG format 
 plot_ts(ENSO.dat, NHSST.dat, binnedts[,1:2], binnedts[,c(1,3)], "ENSO-Nino3", 
  "SST NH Mean", colts1=1, colts2=2, colbints1=3, colbints2=4, device="png", 
  Hfig=900, Wfig=1200, resfig=150, ofilename="plot_ts_RAW_BIN_enso_sst") 

 #####################################################################
 #::  Figure 4 (Polanco-Martínez et al. (2019) (mimeo)). 
 #####################################################################

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
 # Testing our plot_ts function: plot_ts.R
 #####################################################################
 # "Screen" 
 plot_ts(ID31.dat, ID32.dat, bin_ts1, bin_ts2, "MD04-2845 (Temp. forest)", 
  "MD95-2039 (Temp. forest )", colts1=1, colts2=2, colbints1=3, colbints2=4,
  device="screen") 

 # PDF format
 plot_ts(ID31.dat, ID32.dat, bin_ts1, bin_ts2, "MD04-2845 (Temp. forest)", 
  "MD95-2039 (Temp. forest )", colts1=1, colts2=2, colbints1=3, colbints2=4, 
  device="pdf", Hpdf=6, Wpdf=9, resfig=300, ofilename="ts_ACER_ABRUPT") 

 # PNG format 
 plot_ts(ID31.dat, ID32.dat, bin_ts1, bin_ts2, "MD04-2845 (Temp. forest)", 
  "MD95-2039 (Temp. forest )", colts1=1, colts2=2, colbints1=3, colbints2=4, 
  device="png",  Hfig=900, Wfig=1200, resfig=150, ofilename="ts_ACER_ABRUPT") 
}
}
\keyword{plot}
\keyword{plot time series}
\keyword{plot_ts}
