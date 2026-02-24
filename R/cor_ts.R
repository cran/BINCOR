######################################################################
#:: cor_ts function - R package BINCOR                               #
#:: Programmed by Josué M. Polanco-Martinez a.k.a jomopo             #
#:: Email: josue.m.polanco@gmail.com                                 #
######################################################################
#   Copyright (C) 2017-2026 by Josué M. Polanco-Martínez 	             #
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

cor_ts <- 
function(bints1, bints2, varnamets1="NULL", varnamets2="NULL", KoCM, rmltrd="N", 
         device="screen", Hfig, Wfig, Hpdf, Wpdf, resfig, ofilename) {

 #:: Checking the input data 
 if( dim(bints1)[1] != dim(bints2)[1] ) 
  stop ("The binned time series under analysis do not have the same 
    number of elements. Thank you for using our BINCOR package.") 

 if( dim(bints1)[2] !=2 | dim(bints2)[2] != 2) 
   stop ("There is a problem with the input data. The input data should 
    be a couple of vectors of dimension N x 2 (rows x columns). Please, 
    use in the R's command line: dim(ts'x')[2] to verify the number of 
    columns. Thank you for using our BINCOR package.") 

 if (rmltrd == "N" || rmltrd == "n") 
  cor.ts  <- cor.test(bints1[,2], bints2[,2], method=KoCM) 
 if (rmltrd == "Y" || rmltrd == "y") 
 #:: The linear trend is removed -the R pack. "pracma" is required! 
  cor.ts  <- cor.test(c(detrend(bints1[,2])), c(detrend(bints2[,2])), method=KoCM) 

 TOPRINT <- paste("The binned ", KoCM, "'s correlation coefficient is ", sep="")

 if (KoCM=="pearson") 
  cat(paste(TOPRINT, round(cor.ts$estimate,4), " [", round(cor.ts$conf.int[1], 4), "; ",  
  round(cor.ts$conf.int[2],4), "]", sep=""), "\n") 
 if (KoCM=="spearman" || KoCM=="kendall")
  cat(paste(TOPRINT, round(cor.ts$estimate,4), "," ," p-value = ", 
  round(cor.ts$p.value, 15), sep=""), "\n")

 if (device=="png") {
  fileout <- file.path(
    tempdir(),
    paste0("scatterplot_", ofilename, ".png")
  )
  png(fileout, height=Hfig, width=Wfig, res=resfig) 
 }

 if (device=="jpeg" || device=="jpg") {
  fileout <- file.path(
    tempdir(),
    paste0("scaterplot_", ofilename, ".jpg")
  )
  jpeg(fileout, height=Hfig, width=Wfig, res=resfig)
 }

 if (device == "pdf") {
  fileout <- file.path(
    tempdir(),
    paste0("scaterplot_", ofilename, ".pdf")
  )
  pdf(fileout, height = Hpdf, width = Wpdf)
 }
  
 plot(bints1[,2], bints2[,2], t="p", xlab=varnamets1, 
  ylab=varnamets2, las=1, pch=16) 

 if (device != "screen")

 dev.off()

 return(cor.ts)

} 
