######################################################################
#:: bin_cor function - R package BINCOR                              #
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

bin_cor <- function(ts1, ts2, FLAGTAU=3, ofilename) { 

 #:: inputs: 
 #:: ts1 and ts2 are the time series under analysis (the first column 
 #:: are the times/ages (in ascending order) & the second column are the 
 #:: elements of the variable under stiudy). "ofilename" is the output's 
 #:: filename, which will contains the binned data.  
 #:: FLAGTAU (the persistence method). 
 #:: Options (by default it is 3): 
 #:: 1 tau_x + tau_y     [Eq. 7.44, Mudelsee (2010, 2014)]
 #:: 2 max(tau_x, tau_y) [Eq. 7.45, Mudelsee (2010, 2014)]
 #:: 3 dist_x_y/ln(a_x_y_est) [Eq. 7.48, Mudelsee (2010, 2014)]
 #:: tau_x, tau_y are the persistence (memory) time for ts1 and ts2, 
 #:: respectively. 
 
 #:: Checking the input data
 if(dim(ts1)[2] !=2 | dim(ts2)[2] != 2) 
  stop ("There is a problem with the dimension in your input data. 
   The input data should be a couple of vectors of dimension N x 2 
   (rows x columns). Thank you for using the BINCOR package. \n") 

 if(length(which(diff(ts1[,1]) <= 0)) | length(which(diff(ts2[,1]) <= 0)))  
  stop ("There are some times/ages that are not strictly monotonic ascending. 
  Please, check your data. Thank you for using the BINCOR package. \n") 
 
 if(all.equal(ts1[,1], ts2[,1]) == "TRUE")
  cat("The time series have the same timescales, it's not necessary to 
   perform the binned procedure, but it will be computed. \n") 

 #:: Getting the names of the variables 
 names.ts1 <- names(ts1) 
 names.ts2 <- names(ts2) 
 if(is.null(colnames(ts1))) colnames(ts1) <- c("Time-ts1", "Variable-ts1") 
 if(is.null(colnames(ts2))) colnames(ts2) <- c("Time-ts2", "Variable-ts2") 

 #:: Getting Nx & Ny (number of elements for each time series)
 Nx     <- length(ts1[,1])
 Ny     <- length(ts2[,1])			

 #:: Computing Eqs. 7.46 (Mudelsee 2010, 2014)   # or smarter in R ;-) 
 dist_X  <- (ts1[Nx,1] - ts1[1,1])/ (Nx - 1)     # mean(diff(ts"x"[,1]))
 dist_Y  <- (ts2[Ny,1] - ts2[1,1])/ (Ny - 1)     
 Tmax_m  <- max(ts1[Nx,1], ts2[Ny,1])            
 Tmin_m  <- min(ts1[1,1],  ts2[1,1])             
 dist_XY <- (Tmax_m - Tmin_m) / (Nx + Ny - 1)

 #:: Getting Tau (persistence time or memory for each t.s.) 
 #:: There are several ways to get 'tau', but we use Mudelsee (2002). 
 #:: This tau estimation can be obtained in R from the REDFIT (Schulz 
 #:: & Mudelsee 2002) function included in the 'dplR' R package 
 #:: (Bunn et al 2015. https://cran.r-project.org/package=dplR). 
 #:: We compute the raw spectrum (dof=2) for each time series (n50=1 & iwin=0) 
 #:: in order to get tau (subroutine redfitTauest), we have modified slightly 
 #:: the version of the redfit subroutine redfitTauest. This piece of 
 #:: code ("tauest_dplR.R") is provided in our BINCOR package. 
 
 tau1 <- unlist(redfitTauest(ts1[,1], ts1[,2])) #If you use tauest_dplR.R
 tau2 <- unlist(redfitTauest(ts2[,1], ts2[,2]))

 #If you use subroutine redfitTauest from 'dplR' R package 
 #redfit.ts1 <- redfit(ts1[,2], ts1[,1], ofac=1, n50=1, iwin=0) 
 #redfit.ts2 <- redfit(ts2[,2], ts2[,1], ofac=1, n50=1, iwin=0)
 #tau1       <- redfit.ts1$tau
 #tau2       <- redfit.ts2$tau

 a_X_est   <- tau1[2]
 a_Y_est   <- tau2[2] 
 tau_X_est <- tau1[3]
 tau_Y_est <- tau2[3]

 #:: FLAGTAU (persistence for both time series) options! 
 if (FLAGTAU == 1) {
  #:: Eq. 7.44 (Mudelsee 2010 & 2014). 
  taub <- tau_X_est + tau_Y_est
  cat("Hi!, option 1: taub <- tau_X_est + tau_Y_est [Eq. 7.44 (Mudelsee 2010 & 2014)] \n") 
 }
  
 if (FLAGTAU == 2) {
  #:: Eq. 7.45 (Mudelsee 2010 & 2014) 
  taub <- max(tau_X_est, tau_Y_est) 
  cat("Hi!, option 2: taub <- max(tau_X_est, tau_Y_est) [Eq. 7.45 (Mudelsee 2010 & 2014)] \n") 
 }

 if (FLAGTAU == 3) {
  #:: Eq. 7.47 (Mudelsee 2010 & 2014)
  a_XY_est <- sqrt(a_X_est*a_Y_est)
  #:: Eq. 7.48 (Mudelsee 2010 & 2014)  
  taub     <- -dist_XY / log( a_XY_est ) 
  cat("Hi!, option 3: taub <- -dist_XY / log(a_XY_est) [Eq. 7.47 & 7.48 (Mudelsee 2010 & 2014)] \n") 
 }

 #:: Inspired from the M. Mudelsee's code "mc-brxy.f90") 
 taub <- min(taub, (Tmax_m - Tmin_m)*0.5)  
 taub <- max(taub, (Tmax_m - Tmin_m)/(Nx - 1)) 
 
 #:: Computing the number of "bins"
 remi <- (Tmax_m - Tmin_m) / taub
 Nb   <- round(remi) 

 cat("Testing the number of bins: taub=", taub," Nb=", Nb,"\n") 
 
 id1        <- rep(9999, Nb)  
 id2        <- rep(9999, Nb)  
 mean.ts1   <- rep(9999, Nb)  
 mean.ts2   <- rep(9999, Nb)  
 tau.t.mean <- rep(9999, Nb)  
 limI       <- Tmin_m

 #:: Here, the binned time series are created! 
 for (N in 1:Nb) { 
  limS  <- Tmin_m + N*taub
  if (N == 1) { 
   id1t <- which(ts1[,1] >= limI & ts1[,1] <= limS) 
   id2t <- which(ts2[,1] >= limI & ts2[,1] <= limS) 
  }
  if (N > 1)  { 
   id1t <- which(ts1[,1] >  limI & ts1[,1] <= limS) 
   id2t <- which(ts2[,1] >  limI & ts2[,1] <= limS) 
  }
  id1[N] <- length(id1t)
  id2[N] <- length(id2t)

  #:: Evaluating IF a "bin" contains BOTH more than zero 
  #:: X (ts1) & Y (ts2) points (pp. 312, Mudelsee 2010)
  if (id1[N] & id2[N] > 0) {
   mean.ts1[N]   <- mean(ts1[id1t,2]) 
   mean.ts2[N]   <- mean(ts2[id2t,2]) 
   tau.t.mean[N] <- mean(c(limI, limS))  #or (limI + limS)/2  
  }
  else { 
  #:: This's a simple way to face this "problem", but you need to remove 
  #:: the NA's to estimate the correlation btw "bin ts1" and "bin ts2". 
   mean.ts1[N]   <- NA 
   mean.ts2[N]   <- NA 
   tau.t.mean[N] <- mean(c(limI, limS))  #(limI + limS)/2
  } 
  limI   <- limS
 } 

 Datin    <- cbind(tau.t.mean, mean.ts1, mean.ts2) 
 # mean.ts1 and mean.ts2 are the binned time series

 #:: Computing some basic statistics 
 id.noNA     <- which(tau.t.mean != "NA") 
 avg.bin     <- round(mean(diff(tau.t.mean[id.noNA])), 2)
 #avg.bin     <- mean(diff(na.omit(tau.t.mean))) 
 VAR.ts1     <- round(cbind(var(ts1[,2]), var(na.omit(mean.ts1))), 2)
 VAR.ts2     <- round(cbind(var(ts2[,2]), var(na.omit(mean.ts2))), 2)
 chg.VARts1  <- round(VAR.ts1[1] - VAR.ts1[2], 2)
 chg.VARts2  <- round(VAR.ts2[1] - VAR.ts2[2], 2)
 per_chg.VARts1 <- round((chg.VARts1 / VAR.ts1[1])*100, 2)
 per_chg.VARts2 <- round((chg.VARts2 / VAR.ts2[1])*100, 2)
 

 outfile <- file.path(tempdir(), ofilename)
 write.table(Datin, file=outfile, col.names=F, row.names=F)
 
 names.ls <- c("Binned_time_series", "Auto._cor._coef._ts1", "Persistence_ts1", 
	       "Auto._cor._coef._ts2", "Persistence_ts2", "bin width", "Number_of_bins", 
    	       "Average spacing", "VAR. ts1", "VAR. bin ts1", "VAR. ts2", "VAR. bin ts2", 
               "VAR. ts1 - VAR bints1", "VAR. ts2 - VAR bints2", "% of VAR. lost ts1", 
	       "% of VAR. lost ts2")

 LIST        <- list(Datin, a_X_est, tau_X_est, a_Y_est, tau_Y_est, taub, Nb,
		    avg.bin, VAR.ts1[1], VAR.ts1[2], VAR.ts2[1], VAR.ts2[2], 
		    chg.VARts1, chg.VARts2, per_chg.VARts1, per_chg.VARts2)
 names(LIST) <- names.ls 
  
 return(LIST) 

}

