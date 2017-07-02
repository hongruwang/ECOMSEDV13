      SUBROUTINE PRTXY(FD,FMASK,IM,JM,KB,KT,IARRAY,SCALE,IUNIT,DEV)
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
      SAVE
      CHARACTER*3 DEV
      CHARACTER*7 LEVEL(2)
      DIMENSION FD(IM,JM,KB),IARRAY(IM,JM),FMASK(IM,JM)
      DIMENSION KP(2)
      KP(1)=KT
      KP(2)=KB-1
      LEVEL(1)='SURFACE'
      LEVEL(2)='BOTTOM '
C
C-----------------------------------------------------------------------
C        THIS SUBROUTINE WRITES HORIZONTAL LAYERS OF A 3-D FIELD
C-----------------------------------------------------------------------
C
      DO 10 KM=1,2
      K=KP(KM)
      WRITE(IUNIT,20) LEVEL(KM)
 20   FORMAT(//A7,' LAYER')
      CALL PRINT(FD(1,1,K),FMASK,IM,JM,IARRAY,SCALE,IUNIT,DEV)
 10   CONTINUE
C
      RETURN
      END
