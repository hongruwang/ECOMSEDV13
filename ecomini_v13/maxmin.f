      SUBROUTINE MAXMIN (F, FMASK, IX, JY, FMAX , FMIN)
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
* Point of Contact: ecomsed-support@hydroquap.com                        *
C*************************************************************************
      SAVE
      DIMENSION F(IX,JY),FMASK(IX,JY)
C
      BIG=1.E20
C
      FMAX=-BIG
      FMIN=BIG
      DO 220 J=1, JY
      DO 220 I=1, IX
      IF(FMASK(I,J).EQ.0.) GO TO 220
      FMAX=AMAX1(FMAX,F(I,J))
      FMIN=AMIN1(FMIN,F(I,J))
  220 CONTINUE
C
      RETURN
      END
