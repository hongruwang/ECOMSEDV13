      SUBROUTINE DISPLY(ARRAY,AMASK,IM,JM,KK,SCALE,IUNIT)
      SAVE
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
C-----------------------------------------------------------------------
C      ARRAY=2D ARRAY(IM,JM) TO BE PRINTTED ON PAPER
C      IM   =IM DIMENSION OF ARRAY
C      JM   =JM DIMENSION OF ARRAY, IF ARRAY IS VERTICAL SLAB, JM=0
C      KK   =2ND DIMENSION OF ARRAY IF VERTICAL SLAB, OTHERWISE KK=0
C      SCALE=SCALES ALL VALUES BY CONSTANT AMOUNT
C-----------------------------------------------------------------------
C
      DIMENSION ARRAY(IM,1),NUM(20),PLINE(20)
      DIMENSION AMASK(IM,KK)
      IF(SCALE.NE.0) GO TO 500
      DO 251 IS=1,IM,12
      IE=IS+11
      IF(IE.GT.IM) IE=IM
      IDIF=IE-IS+1
      DO 100 I=1,IDIF
 100  NUM(I)=IS+I-1
      WRITE(IUNIT,9990) (NUM(I),I=1,IDIF)
 9990 FORMAT(12I10,/)
      JMORKM=JM+KK
      DO 252 JORK=1,JMORKM
      IF(JM.NE.0) L=JMORKM+1-JORK
      IF(KK.NE.0) L=JORK
      WRITE(IUNIT,9966) L, (ARRAY(I,L),I=IS,IE)
 252  CONTINUE
      WRITE(IUNIT,9984)
 251  CONTINUE
 9966 FORMAT(1X,I2,12(1PE10.2))
 9984 FORMAT(//)
      RETURN
 500  CONTINUE
C
      IDIM=JM
      IF(JM.EQ.0) IDIM=KK
      CALL MAXMIN(ARRAY,AMASK,IM,IDIM,VMAX,VMIN)
      WRITE(IUNIT,110) SCALE,VMAX,VMIN
  110 FORMAT(//2X,'ACTUAL VALUE = PRINTED VALUE/',1PE9.2,
     .      5X,'MAX VALUE = ',1PE9.2,2X,'MIN VALUE = ',1PE9.2/)
      DO 751 IS=1,IM,20
      IE=IS+19
      IF(IE.GT.IM) IE=IM
      IDIF=IE-IS+1
      DO 600 I=1,IDIF
 600  NUM(I)=IS+I-1
      WRITE(IUNIT,9991) (NUM(I),I=1,IDIF)
 9991 FORMAT(3X,20I6,/)
      JMORKM=JM+KK
      DO 752 JORK=1,JMORKM
      IF(JM.NE.0) L=JMORKM+1-JORK
      IF(KK.NE.0) L=JORK
      DO 753 I=1,IDIF
      PLINE(I)=ARRAY(IS+I-1,L)*SCALE
 753  CONTINUE
      WRITE(IUNIT,9967) L,(PLINE(I),I=1,IDIF)
 752  CONTINUE
      WRITE(IUNIT,9984)
 751  CONTINUE
 9967 FORMAT(1X,I2,1X,20F6.2)
      RETURN
      END
