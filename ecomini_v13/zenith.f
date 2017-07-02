      SUBROUTINE ZENITH(FJD,DLT,SLAG,XLAT,HA,DHR,COSZ,FRAC,XLNG0,
     1  GMT)
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
C    *******************************************************************
C    *                                                                 *
C    *                           Z E N I T H                           *
C    *                                                                 *
C    *******************************************************************
C
      LOGICAL RISE,  SET
      DATA  PID24/.13089969/
C
C    *******************************************************************
C
C     ZENITH COMPUTES EFFECTIVE MEAN COSINE OF ZENITH ANGLE AND DAYLIGHT
C       FRACTION FROM LATITUDE AND PARAMETERS GIVEN BY SUBROUTINE SOLAR
C
C     INPUT ARGUMENTS TO CIRCULAR FUNCTIONS ARE IN RADIANS
C
      PI=ACOS(-1.0)
      TPI=2.0*PI
      CVPR=TPI
C
C  SHIFT TIME BY GMT HOURS TO PUT IT IN GMT
C
      FJDL=FJD-GMT/24.0
C
      IF (FJDL.LT.0.0) FJDL=FJDL+1.0
      GHA=FJDL*TPI+SLAG
      ARG=DHR*PID24
      SINFAC=SIN(ARG)/ARG
      SS=SIN(XLAT)*SIN(DLT)
      CC=COS(XLAT)*COS(DLT)
      IF(HA.GT.0.0) CONS=SS+CC*SIN(HA)/HA
C
C   SET LONGITUDE IN RADIAN
C
      XLNG=(XLNG0+180.0)*PI/180.0
      HLOC=GHA+ARG+XLNG
C     LOCAL HOUR ANGLE SHIFTED BY HALF OF THE AVERAGING PERIOD
      RISE=.FALSE.
      SET=.FALSE.
      HLOC=MOD(HLOC,TPI)
      IF(HLOC.GT.PI) HLOC=HLOC-TPI
      HLPAR=HLOC+ARG
      ARMHL=ARG-HLOC
      IF(HLPAR.GT.HA) SET=.TRUE.
      IF(ARMHL.GT.HA) RISE=.TRUE.
      IF(RISE.AND.SET) GO TO 7
      IF(HLPAR.GT.PI) GO TO 8
      IF(ARMHL.GT.PI) GO TO 9
      IF(SET) GO TO 3
      IF(RISE) GO TO 4
      FRAC=1.0
      COSZ=SS+CC*COS(HLOC)*SINFAC
      GO TO 2
  3   DELSH=0.5*(HA+ARMHL)
      GO TO 5
  4   DELSH=0.5*(HA+HLPAR)
  5   IF(DELSH.LE.0.0) GO TO 6
      FRAC=DELSH/ARG
C	write(*,*)'DELSH=',DELSH
      COSZ=SS+CC*COS(HA-DELSH)*SIN(DELSH)/DELSH
      GO TO 2
  7   IF(HA.LE.0.0) GO TO 6
      FRAC=HA/ARG
      COSZ=CONS
      GO TO 2
  8   DELE=0.5*MAX(HLPAR+HA-TPI,0.0)
      DELW=0.5*MAX(HA+ARMHL,0.0)
      GO TO 10
  9   DELE=0.5*MAX(HA+HLPAR,0.0)
      DELW=0.5*MAX(ARMHL+HA-TPI,0.0)
 10   FRAC=(DELE+DELW)/ARG
      IF(FRAC.EQ.0.0) GO TO 11
      COSZ=SS+CC*(COS(HA-DELE)*SIN(DELE)+
     1               COS(HA-DELW)*SIN(DELW))/(DELE+DELW)
      GO TO 2
  6   FRAC=0.0
 11   COSZ=0.0
  2   CONTINUE
      COSZ = MIN(1.0,COSZ)
      COSZ = MAX(0.0,COSZ)
      FRAC = MIN(1.0,FRAC)
      RETURN
      END
