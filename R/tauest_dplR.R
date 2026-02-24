#####################################################################
#:: This piece of code comes from the R package dplR (redfit function
#:: based on Schulz & Mudelsee 2002) programmed by  Mikko Korpela. 
#:: Very minor modifications was done by J. M. Polanco-Martínez 
#:: (josue.m.polanco@gmail.com) in order to used the dplR redfit 
#:: function from the R package BINCOR. 
#:: 10/2016, Bordeaux, FR
#####################################################################

#####################################################################
### This part of dplR was (arguably non-trivially) translated and
### adapted from public domain Fortran program REDFIT, version 3.8e
### (Michael Schulz and Manfred Mudelsee). The possibly non-free parts
### of REDFIT derived from Numerical Recipes were not used.
### http://www.geo.uni-bremen.de/geomod/staff/mschulz/
### Author of the dplR version is Mikko Korpela.
###
### Copyright (C) 2013-2015 Aalto University
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
#####################################################################

## dplR: R version based on Mudelsee's code.
## dplR: Introduction copied from REDFIT (some variables removed).
##
## Manfred Mudelsee's code for tau estimation
## ----------------------------------------------------------------------
##  TAUEST: Routine for persistence estimation for unevenly spaced time series
## ----------------------------------------------------------------------
##        Main variables
##
##        t       :       time
##        x       :       time series value
##        np      :       number of points
##       dt       :       average spacing
##    scalt       :       scaling factor (time)
##      rho       :       in the case of equidistance, rho = autocorr. coeff.
##     mult       :       flag (multiple solution)
##     amin       :       estimated value of a = exp(-scalt/tau)

redfitTauest <- function(t, x) {

    np <- length(t)
    ## Correct time direction; assume that ages are input
    ## dplR: Correction of time direction is done by modifying this
    ## function and redfitMinls, not by explicitly reversing (and
    ## multiplying by one)
    ## tscal <- -rev(t)tauest_dplR.R
    ## xscal <- rev(x)
    ## Scaling of x
    xscal <- x / sd(x)
    ## Scaling of t (=> start value of a = 1/e)
    dt    <- (t[np] - t[1]) / (np - 1)
    ## dplR: rhoest() of REDFIT is now an "inline function" of two
    ## lines + comment line:
    ## Autocorrelation coefficient estimation (equidistant data)
    xscalMNP <- xscal[-np]
    rho      <- sum(xscalMNP * xscal[-1]) / sum(xscalMNP * xscalMNP)

    if (rho <= 0) {
        rho <- 0.05
        warning("rho estimation: <= 0") # comentar tauest_dplR.R
    } else if (rho > 1) {
        rho <- 0.95
        warning("rho estimation: > 1") # comentar  
    }
    scalt  <- -log(rho) / dt
    tscal  <- t * scalt
    ## Estimation
    minRes <- redfitMinls(tscal, xscal)
    amin   <- minRes[["amin"]]
    mult   <- minRes[["nmu"]]
    #warnings <- FALSE 
  
  ##################################################################
  #:: Piece of code*** added by Josué M. Polanco-Martínez:  START 
  ##################################################################
  #:: We use a piece of code*** (translated from Fortran to R) from 
  #:: mc-brxy.f90 ((C) M. Mudelsee) in order to avoid some problems 
  #:: in the estimation of the persistence. Please look at the 
  #:: subroutine tauest in mc-brxy.f90. 
  ##################################################################
  rhoavgmax=0.99 
  rhoavgmin=0.01
  zero= 0.0 
  # determines also rho = rhoavg
      if (amin >= rhoavgmax) {
           rhoavg=rhoavgmax
           tau = -dt / log(rhoavgmax)
          }
      else if (amin >= rhoavgmin & amin < rhoavgmax) {
           tau = -1.0 /(scalt*log(amin))
           rhoavg = exp(-dt / tau)
           # Bias correction (unknown mean)
           rhoavg = (rhoavg * (np - 1.0) + 1.0) / (np - 4.0) 
           rhoavg = min(rhoavgmax,rhoavg)
           rhoavg = max(rhoavgmin,rhoavg)
           tau=-dt/log(rhoavg)
          }
      else if (amin < rhoavgmin) {
           rhoavg=rhoavgmin
           tau = zero
	  }

      rhoout=rhoavg

  salida <- c(rho,  rhoout, tau) 
  return(list(salida))
  ##################################################################
  #:: Piece of code added by Josué M. Polanco-Martínez:  END 
  ##################################################################
}

##################################################################
if(0){ # original code by  Mikko Korpela, please note that we are
       # not using. Josué M. Polanco-Martínez 
    warnings <- FALSE
    if (mult) {
        warning("estimation problem: LS function has > 1 minima")
        warnings <- TRUE
    }
    if (amin <= 0) {
        warning("estimation problem: a_min =< 0")
        warnings <- TRUE
    } else if (amin >= 1) {
        warning("estimation problem: a_min >= 1")
        warnings <- TRUE
    }
    if (!warnings) {
        ## determine tau
        tau <- -1 / (scalt * log(amin))
 ## By jomopo, 10/2016  
        cat("tau", tau, "\n") 
 ## By jomopo, 10/2016  
        ## determine rho, corresponding to tau
        exp(-dt / tau)
    } else {
        ## dplR: fail early
        stop("error in tau estimation")
    }
}
##################################################################

## dplR: Minimization of the built-in least-squares function lsfun
redfitMinls <- function(t, x) {
    ## Least-squares function
    lsfun <- function(a, difft, xM1, xMNP) {
        if (a > 0) {
            tmp <- xMNP - xM1 * a^difft
        } else if (a < 0) {
            tmp <- xMNP + xM1 * (-a)^difft
        } else {
            tmp <- xMNP
        }
        sum(tmp * tmp)
    }
    a_ar1 <- exp(-1) # 1 / e
    tol   <- 3e-8    # Brent's search, precision
    tol2  <- 1e-6    # multiple solutions, precision
    difft <- diff(t)
    np <- length(x)
    xM1 <- x[-1]
    xMNP <- x[-np]
    opt1 <- optimize(lsfun, c(-2, 2),     tol = tol, difft = difft,
                     xM1 = xM1, xMNP = xMNP)
    opt2 <- optimize(lsfun, c(a_ar1, 2),  tol = tol, difft = difft,
                     xM1 = xM1, xMNP = xMNP)
    opt3 <- optimize(lsfun, c(-2, a_ar1), tol = tol, difft = difft,
                     xM1 = xM1, xMNP = xMNP)
    a_ar11 <- opt1[["minimum"]]
    a_ar12 <- opt2[["minimum"]]
    a_ar13 <- opt3[["minimum"]]
    dum1 <- opt1[["objective"]]
    dum2 <- opt2[["objective"]]
    dum3 <- opt3[["objective"]]
    list(amin = c(a_ar11, a_ar12, a_ar13)[which.min(c(dum1, dum2, dum3))],
         nmu = ((abs(a_ar12 - a_ar11) > tol2 && abs(a_ar12 - a_ar1) > tol2) ||
                (abs(a_ar13 - a_ar11) > tol2 && abs(a_ar13 - a_ar1) > tol2)))
}
