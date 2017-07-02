      SUBROUTINE SLICEXZ(FD,AMASK,IM,JM,KB,JROW,ARRAY,SCALE,IUNIT)
      SAVE
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
      DIMENSION FD(IM,JM,KB),ARRAY(IM,KB)
      DIMENSION AMASK(IM,JM)
C
      WRITE(IUNIT,5) JROW
 5    FORMAT(///' VERTICAL SECTION AT ROW = ',I4///)
C
      DO 10 K=1,KB
      DO 10 I=1,IM
 10   ARRAY(I,K)=FD(I,JROW,K)
C
      CALL DISPLY(ARRAY,AMASK,IM,0,KB,SCALE,IUNIT)
C
      RETURN
      END
