      Subroutine NCLD(YDAY,ALAT,ALON,SWOBS,ATC,CLD,CUTOFF,RAD)
C*************************************************************************
C     ECOMSED MODEL
C     VERSION 1.3
C     FEBRUARY 2002
C*************************************************************************
*               Copyright (C) 2002, HydroQual, Inc.                      *
*                                                                        *
*  These coded instructions, statements, and computer programs  contain  *
*  unpublished  proprietary  information of HydroQual, Inc., and         *
*  are protected by Federal copyright law.  They  may  not be disclosed  *
*  to  third  parties  or copied or duplicated in any form, in whole or  *
*  in part, without the prior written consent of HydroQual, Inc.         *
*                                                                        *
* Point of Contact: ecomsed-support@hydroquap.com                        *
C*************************************************************************
c
c Subroutine to calculate percent cloud cover based on observed
c shortwave radiation and day of the year.  Converted from the
c original program to a subroutine. -hlj-        June 21, 1993.
c
c ORIGINAL COMMENTS ------------------------------------------------
c n_cldfrac.f		A.Oien/A. Plueddemann		May 19, 1992
c
c This program estimates the daytime cloud fraction based on the
c ratio of observed to theoretical shortwave radiation.
c No estimate can be made at night, and the cloud fraction is set
c to a flag value.  It is assumed that the nighttime values will
c later be interpolated.  The program is similar to n_clearsky
c except that here the timebase and observed shortwave (with
c nighttime bias removed) are read in from the met data file.  
c Output data is handled more nicely, written to user-defined 
c output data file.
c Values for ATC and CUTOFF are determined empirically from
c iterations of n_clearsky and this program.  See n_clearsky and
c n_cloud for further description of ATC and CUTOFF.
c ORIGINAL COMMENTS ------------------------------------------------
c
c Since this subroutine will be run concurrently with the numerical
c circulation model, cloud cover is persisted through the night using
c the last cloud cover value obtained before cutoff was reached,
c rather than interpolated over the night. -hlj-
c
c ORIGINAL COMMENTS ------------------------------------------------
c 6-23-92
c Modified from n_clearsky to use time series from met data 
c file rather than rather than calculating from start time and
c increment.  Renamed TIME(K) array to yday(k) array to avoid 
c confusion with time variable in cloud subroutine. 
c ORIGINAL COMMENTS ------------------------------------------------

c	set lat, lon of observations, compute lon-based time correction
c     ALAT = 59.5935
c     ALON = 20.9642
      COR = ALON / 360.

c	set desired ATC and CUTOFF
c       started with ATC=.83 , CUTOFF=.05 as with n_clearsky,
c	final version uses ATC=0.77, CUTOFF=0.10 
      ATC = 0.60      
      CUTOFF = 0.05  

C   be careful of sign of COR here (alon is positive W for this convention)

      TD = YDAY - (COR)

CMNT  CORRECT LOCAL TIME TO TRUE SOLAR TIME USING
CMNT  FIT TO TIME EQUATION DATA PROVIDED ON
CMNT  SMITHSONIAN MET TABLES, PAGE 495

      TWOPI = 2. * 3.141592
      TEQN = -5.5 - 9. * COS((TD-45.6)*(TWOPI/182.5))
      TEQN = TEQN / (24.*60.)
      TD = TD + TEQN
      Call N_CLOUD(CUTOFF,ALAT,TD,SWOBS,RAD,ATC,CLD)
      Return
      End
