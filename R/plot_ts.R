######################################################################
#:: plot_ts - R package BINCOR                                       #
#:: Programmed by Josué M. Polanco-Martinez a.k.a jomopo             #
#:: Email: josue.m.polanco@gmail.com                                 #
######################################################################
#   Copyright (C) 2017 by Josué M. Polanco-Martínez 	             #
#   This file/code is part of the R package BINCOR 	             #
######################################################################
#								     
#   BINCOR is free software: you can redistribute it and/or modify it
#   it under the terms of the GNU General Public License as published 
#   by the Free Software Foundation, either version 3 of the License, 
#   or (at your option) any later version.
#
#   BINCOR is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with BINCOR If not, see <http://www.gnu.org/licenses/>.
#
#####################################################################

plot_ts <- 
function(ts1, ts2, bints1, bints2, varnamets1="", varnamets2="", 
         colts1=1, colts2=1, colbints1=2, colbints2=2, ltyts1=1, ltyts2=1, 
         ltybints1=2, ltybints2=2, device="screen", Hfig, Wfig, Hpdf, Wpdf, 
         resfig, ofilename) { 

 #:: Checking the input data 
  if( dim(bints1)[1] != dim(bints2)[1] ) 
  stop ("The binned time series under analysis do not have the same 
   number of elements. Thank you for using our BINCOR package.") 

 if( dim(ts1)[2] !=2 | dim(ts2)[2] != 2) 
   stop ("There is a problem with the input data. The input data 
    should be a couple of vectors of dimension N x 2 (rows x columns). 
    Please, use in the R's command line: dim(ts'x')[2] to verify the 
    number of columns. Thank you for using our BINCOR package.") 

 #:: Getting the max./min. values for times 
 range.timesORts1  <- range(ts1[,1])
 range.timesORts2  <- range(ts2[,1])

 #:: Devices options: png, jpg & pdf 
 if (device=="png") {
  fileout <- paste("plot_", ofilename, "%02d.png", sep="")
  png(fileout, height=Hfig, width=Wfig, res=resfig) 
 }

 if (device=="jpeg" || device=="jpg") {
  fileout <- paste("plot_", ofilename, "%02d.jpg", sep="")
  jpeg(fileout, height=Hfig, width=Wfig, res=resfig)
 }

 if (device=="pdf") {
  fileout <- paste("plot_", ofilename, ".pdf", sep="")
  pdf(fileout, height=Hpdf, width=Wpdf)
 }

 par(mar=c(5, 5, 3, 5) + 0.1)
 #if (device=="screen") dev.new()
  ######
  plot(ts1, t="b", col=colts1, xlab="Time", ylab="", lty=ltyts1, 
   main="", xlim=range.timesORts1, las=1, cex=0.75) 
  #par(new=T)
  points(bints1, t="b", col=colbints1, xaxt="n", yaxt="n", 
   xlab="", ylab="", lty=ltybints1, cex=0.75)
  mtext(2, text=varnamets1, line=3.5, cex=0.95)
  legend("topleft", bty="n", legend=c(paste(varnamets1, " (primary).", 
   " N = ", length(ts1[,1]), " elements", sep=""), paste(varnamets1, 
   " (binned).", " N = ", length(bints1[,1]), " elements", sep="")), 
   lty=c(ltyts1, ltybints1), col=c(colts1,colbints1))
  axis(3, ts1[,1], labels=F, col=colts1) 
  axis(3, bints1[,1], labels=F, line=1, col=colbints1) 

 if (device=="screen") dev.new()
 par(mar=c(5, 5, 3, 5) + 0.1)
  ##### 
  plot(ts2, t="b", col=colts2, xlab="Time", ylab="", lty=ltyts2, 
   main="", xlim=range.timesORts2, las=1, cex=0.75) 
  #par(new=T)
  points(bints2, t="b", col=colbints2, xaxt="n", yaxt="n", 
   xlab="", ylab="", lty=ltybints2, cex=0.75)
  mtext(2, text=varnamets2, line=3.5, cex=0.95)
  legend("topleft", bty="n", legend=c(paste(varnamets2, " (primary).",  
   " N = ", length(ts2[,1]), " elements",  sep=""), paste(varnamets2, 
   " (binned).", " N = ", length(bints2[,1]), " elements",  sep="")), 
   lty=c(ltyts2, ltybints2), col=c(colts2,colbints2))
  axis(3, ts2[,1], labels=F, col=colts2) 
  axis(3, bints2[,1], labels=F, line=1, col=colbints2) 

 if (device=="screen") dev.new()
 par(mar=c(5, 5, 3, 5) + 0.1)
  #####  
  plot(ts1, t="b", col=colts1, xlab="Time", ylab="", yaxt="n", 
   lty=ltyts1, main="", xlim=range.timesORts1, cex=0.75) 
  axis(2, pretty(ts1[,2]), lwd=2, las=2, col=colts1)
  mtext(2, text=varnamets1, line=3.5, cex=0.95)
  par(new=T)
  plot(ts2, t="b", col=colts2, xaxt="n", yaxt="n", xlab="", ylab="", 
   lty=ltyts2, main="", xlim=range.timesORts1, cex=0.75)
   #lty=ltyts2, main="", xlim=range.timesORts2, cex=0.75)
  axis(4, pretty(ts2[,2]), lwd=2, las=2, col=colts2)
  mtext(4, text=varnamets2, line=3.5, cex=0.95)
  legend("topleft", bty="n", legend=c(paste(varnamets1, " (primary). ", 
   " N = ", length(ts1[,1]), " elements", sep=""), paste(varnamets2, 
   " (primary). ",  " N = ", length(ts2[,1]), " elements",  sep="")),  
   lty=c(ltyts1, ltyts2), col=c(colts1, colts2))
  axis(3, ts1[,1], labels=F, col=colts1) 
  axis(3, ts2[,1], labels=F, line=1, col=colts2) 

 if (device=="screen") dev.new()
 par(mar=c(5, 5, 3, 5) + 0.1)
  ##### 
  plot(bints1, t="b", col=colbints1, xlab="Time", ylab="", 
   yaxt="n", lty=ltybints2, cex=0.75)
  axis(2, pretty(bints1[,2]), lwd=2, las=2, col=colbints1)
  mtext(2, text=varnamets1, line=3.5, cex=0.95)
  par(new=T) 
  plot(bints2, t="b", col=colbints2, xaxt="n", yaxt="n", 
   xlab="", ylab="", lty=ltybints2, xlim=range(bints1[,1]), cex=0.75)
  axis(4, pretty(bints2[,2]), lwd=2, las=2, col=colbints2)
  mtext(4, text=varnamets2, line=3.5, cex=0.95)
  legend("topleft", bty="n", legend=c(paste(varnamets1, " (binned). ", 
   " N = ", length(bints1[,1]), " elements", sep=""), paste(varnamets2, 
   " (binned). ", " N = ", length(bints2[,1]), " elements", sep="")), 
   lty=c(ltybints1, ltybints2), col=c(colbints1, colbints2))
  axis(3, bints1[,1], labels=F, col=colbints1) 
  axis(3, bints2[,1], labels=F, line=1, col=colbints2) 

 if (device != "screen") 

 dev.off()
 
 return()

} 

 
