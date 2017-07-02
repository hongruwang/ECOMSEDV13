      SUBROUTINE VERTVL(DTI2)
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
      INCLUDE 'comdeck'
      SAVE
C
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),ZMFL3D(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C),(ZMFL3D,VHP)
C
C-------- CALCULATE NEW VERTICAL VELOCITY ------------------------------
      DO 100 J=1,JM
      DO 100 I=1,IM
 100  W(I,J,1)=0.0  
C
      DO 101 K=1,KB
      DO 101 J=1,JM
      DO 101 I=1,IM
 101  ZMFL3D(I,J,K)=0.0
C-----------------------------------------------------------------------
C         IMPOSE MASS FLUX BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
      DO 120 N=1,NUMDBC
      ID=IDD(N)
      JD=JDD(N)
      DO 120 K=1,KBM1
C
        ZMFL3D(ID,JD,K)=ZMFL3D(ID,JD,K)+
     +          QDIFF(N)*RAMP*VDDIST(N,K)/100.0/DZ(K)
 120  CONTINUE
C-----------------------------------------------------------------------
      DO 110 K=1,KBM1
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
 110  W(I,J,K+1)=W(I,J,K)
     .    +DZ(K)*((-ZMFL3D(I,J,K)+XMFL3D(I+1,J,K)-XMFL3D(I,J,K)
     .            +YMFL3D(I,J+1,K)-YMFL3D(I,J,K))/ART(I,J)
     .                        +(ETF(I,J)-ETB(I,J))/DTI2 )
C
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  W(I,J,KB)=0.0  
C
      RETURN
      END
