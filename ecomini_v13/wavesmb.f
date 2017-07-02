      SUBROUTINE WAVESMB
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
* Point of Contact: ecomsed-support@hydroqual.com                        *
C*************************************************************************
C
C  CALC. BOTTOM SHEAR STRESS DUE TO WIND WAVES AND CURRENTS
C  USES GRANT/MADSEN/GLENN WAVE MODEL
C**************************************************************
C
      INCLUDE 'comdeck'
      SAVE
C
      PI2=2.*PI
      PIHALF=PI/2.
C
C  CALM WIND ( < 1 MPH = 0.45 m/s) ==> PURE CURRENTS
C
C
C  CALC. SIGNIFICANT WAVE HEIGHT AND PERIOD
C  USE SHALLOW WATER SMB THEORY
C
C  HSIG = SIG. WAVE HEIGHT (m)
C  TSIG = SIG. WAVE PERIOD (s)
C
C  WINDSP = WIND SPEED (m/s)
C  WINDIR = WIND DIRECTION (degrees)
C  HMEAN(N) = MEAN DEPTH ALONG FETCH FOR DIRECTION n (m)
C  FETCH(N) = FETCH FOR DIRECTION n (m)
C
      DO 10 J=2,JM-1
        DO 10 I=2,IM-1
          IF (FSM(I,J).GT.0.0) THEN
C
C
           UWIND=(WU(I,J)*COS(ANG(I,J))+WV(I,J)*SIN(ANG(I,J)))*FSM(I,J)
           VWIND=(-WU(I,J)*SIN(ANG(I,J))+WV(I,J)*COS(ANG(I,J)))*FSM(I,J)
           WSPDSQ=UWIND*UWIND+VWIND*VWIND
           WINDSP=SQRT(WSPDSQ)+1.E-10
C
C     WINDIR is w.r.t. zeta1-direction
C
           WINDIR=ATAN2(VWIND,UWIND)*180./PI
           IF(WINDIR .LT. 0.) WINDIR = WINDIR + 360
           IF(WINDIR .GT. 360.) WINDIR = WINDIR - 360.
C
C     Direction with respect to east (x-direction)
C
           WINDIR = WINDIR + ANG(I,J)*180./PI
           WINDIR = AMOD(WINDIR,360.)
C
C     IN stress.f and in SMB Model WAVE DIR IS IN w.r.t North Clock Wise
C
           WINDIR= 90. - WINDIR
           IF(WINDIR.LE.0.) WINDIR = WINDIR + 360.
           NDIR=NINT(WINDIR/10.)
           IF (NDIR.EQ.0) NDIR=36
           WRAT0=GRAV/(WINDSP*WINDSP)
           HRAT0=0.283/WRAT0
           TRAT0=2.4*PI*WINDSP/GRAV
           WRAT1=WRAT0*HMEAN(NDIR,I,J)
           WRAT2=WRAT0*FETCH(NDIR,I,J)
           HA1=TANH(0.53*(WRAT1**0.75))
           TA1=TANH(0.833*(WRAT1**0.375))
           HSIG(I,J)=HRAT0*HA1*TANH(0.0125*(WRAT2**0.42)/HA1)
           TSIG(I,J)=TRAT0*TA1*TANH(0.077*(WRAT2**0.25)/TA1)
           WAVEDIR(I,J)=WINDIR
C
C  LIMIT WAVE HEIGHT TO BREAKING HEIGHT
C
            TBREAK=TIME*24.0
            HBREAK=0.78*DT(I,J)
            HSIG(I,J)=AMIN1(HBREAK,HSIG(I,J))
          ENDIF
 10   CONTINUE
C
C  OUTPUT WAVE PARAMETERS IF WAVEDYN IS SMBMODEL
C
      WRITE(1001)TIME
      WRITE(1001)HSIG
      WRITE(1001)TSIG
      WRITE(1001)WAVEDIR
C
      RETURN
      END
