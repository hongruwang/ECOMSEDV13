      SUBROUTINE SLICEYZ(FD,AMASK,IM,JM,KB,IROW,ARRAY,SCALE,IUNIT)
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
      DIMENSION FD(IM,JM,KB),ARRAY(JM,KB)
      DIMENSION AMASK(IM,JM)
C
      WRITE(IUNIT,5) IROW
 5    FORMAT(///' VERTICAL SECTION AT COLUMN = ',I4///)
C
      DO 10 K=1,KB
      DO 10 J=1,JM
 10   ARRAY(J,K)=FD(IROW,J,K)
C
      CALL DISPLY(ARRAY,AMASK,JM,0,KB,SCALE,IUNIT)
C
      RETURN
      END
