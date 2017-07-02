      SUBROUTINE DATEHR(IDAY,IMON,IYR,IHR,IMIN,THOURS)
C
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
C     THIS ROUTINES GIVES THE HOUR (THOURS) CORRESPONDING TO DAY, MON,
C     YEAR, HOUR AND MINUTES
C
C  OSSM01.FOR
C
      COMMON /DEVCES/ LPRT,LCON,LTPE,LDSK,LWCON
      DIMENSION MON(12,2)
      DIMENSION LYEAR(4)
      DATA(LYEAR(I),I=1,4)/8784,8760,8760,8760/
      DATA (MON(I,1),I=1,12)/0,31,59,90,120,151,181,212,243,273,
     C304,334/
      DATA (MON(I,2),I=1,12)/0,31,60,91,121,152,182,213,244,274,
     C305,335/
      THOURS=0.0
      IF((IMON.LT.1).OR.(IMON.GT.12)) GOTO 800
      IF((IDAY.LT.1).OR.(IDAY.GT.31)) GOTO 800
c change year from 2 digitis to 4 digits 
      IYRS=IYR-1900
      IF(IYRS.EQ.0) GOTO 12
      DO 10 J=1,IYRS
      IJ=MOD(J,4)
      IF(IJ.EQ.0) IJ=4
      THOURS=THOURS+FLOAT(LYEAR(IJ))
 10   CONTINUE
      K=1
 11   IF(IJ.EQ.4) K=2
      IDAYS=MON(IMON,K)+IDAY-1
      THOURS=THOURS+FLOAT(24*IDAYS)+FLOAT(IHR)+FLOAT(IMIN)/60.
      GOTO 999
 12   IJ=4
      GOTO 11
 800  WRITE(LWCON,101) IDAY,IMON,IYR,IHR,IMIN
 101  FORMAT(10X,'ERROR IN DATE FIELD',I2,'/',I2,'/',I2,1X,I2,':',I2)
 999  RETURN
      END
