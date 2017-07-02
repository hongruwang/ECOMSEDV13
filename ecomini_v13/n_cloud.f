      Subroutine N_CLOUD(CUTOFF,ALAT,TIME,SWOBS,RAD,ATC,CLD)
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
c	A. Plueddemann 6-12-92
c	Modified from R. Weller's original routine CLOUD c. 6-18-87
c
c	This routine accepts:
c	  ALAT = latitude of  observations, decimal degrees.
c	  TIME = time in decimal yearday (example, day 1.5 is local
c		noon on the 2nd day).
c	  SWOBS = observed shortwave radiation, it is assumed that 
c		SWOBS has been corrected for nighttime bias.
c	  ATC = atmospheric transmission coeff., it is expected that
c		an appropriate ATC for the site has been empirically
c		derived.
c	and returns:
c	  RAD = theoretical estimate of total shortwave radiation
c		reaching the earths surface
c	  CLD = cloud cover fraction (0.0 to 1.0) for use in
c		correcting longwave radiation, nighttime values of
c		CLD are set to a flag value assuming that they will
c		later be interpolated.
c
c	routine to abandon attempts to compute cloud cover just prior
c	to sunset and just after dawn to avoid erroneous estmiates 
c	caused by persistent dawn/dusk differences in theoretical vs.
c	observed short wave.  The value of this cutoff must be "tuned"
c	for particular data sets by trial and error.

c	--------------------------------------------------------------

c	solar constant
      SC = 1395
c	conversion of degrees to radians
      TWOPI = 3.14159 * 2.
      DTOR = TWOPI / 360.

c	compute altitude of the sun (angular elevation above the
c	horizon) from Smithsonian Met. Tables formulation.
c	  SINA = sine of suns altitude
c	  THTA = declination of sun 
c	  HR = hour angle of the sun

      YRPH = (TWOPI/365.) * (357.-TIME)
c       yrph uses 357 instead of 355.90826
      THTA = (-23.45*COS(YRPH)) * DTOR
c       thta uses 23.45 instead of 23.26658
      HR = (TIME+0.5) * TWOPI
      SINA = COS(THTA) * COS(ALAT*DTOR) * COS(HR) + SIN(THTA) * 
     *    SIN(ALAT*DTOR)

c	if low sun angle go to nighttime loop 
      If (SINA.GT.CUTOFF) Then
c       conditional cut off before R calculation (after R in old version)

c	estimate cloud cover using scheme described in Smithsonian
c	Meteorological Tables (R.J. List, 1984).  
c	  R = radius vector of the earth, equal to the distance from
c		the center of the earth to the center of the sun
c	  RSS = direct E solar radiation at the earths surface.
c	  RDIFF = diffuse sky radiation.
c	  RAD = total clear sky radiation at earths surface, can be
c		compared to SWOBS to determine ATC.
c	  CLD = cloud cover fraction (0 to 1) used for correction to
c		longwave radiation (Fung et al, 1984).

        TRAD = TWOPI * ((TIME-186.)/365.)
c       trad uses 186 instead of 186.72252
        R = 1.0002 + 0.01671 * COS(TRAD)
c       r uses 1.0002 instead of 1.0001531
        RSS = (SC/(R**2.)) * SINA * (ATC**(1./SINA))
        RDIFF = .5 * (SC/(R**2.)) * SINA * (0.91-(ATC**(1./SINA)))
        RAD = RSS + RDIFF
        CLD = (1.-(SWOBS/RAD)) * (10./7.)
c
c	dont allow cloud fraction greater than 1 or less than 0
        If (CLD.GT.1.) CLD = 1.
        If (CLD.LE.0.) CLD = 0.
      Else

c	nighttime loop, set theoretical radiation 
c	to zero.
        RAD = 0.
c redundant, but included to emphasize that the CLD passed INTO the 
c subroutine is passed OUT of the subroutine, i.e. cloud cover is
c PERSISTED at night.
        CLD = CLD
      End If

      Return
      End
