      SUBROUTINE ADVQ(QB,Q,DTI2,QF)
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
      INCLUDE 'comdeck'
      SAVE
C
C-----------------------------------------------------------------------
C     THIS SUBROUTINE INTEGRATES CONSERVATIVE CONSTITUENT EQUATIONS
C            FOR Q2 AND Q2L
C-----------------------------------------------------------------------
C
      DIMENSION QB(IM,JM,KB),Q(IM,JM,KB),QF(IM,JM,KB)
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C)
C
      DO 10 K=1,KB
      DO 10 J=1,JM
      DO 10 I=1,IM
      XFLUX(I,J,K)=0.0
      YFLUX(I,J,K)=0.0
   10 QF(I,J,K)=0.0
C
C-------- HORIZONTAL ADVECTION -----------------------------------------
      DO 20 K=2,KBM1
      DO 20 J=2,JM
      DO 20 I=2,IM
      XFLUX(I,J,K)=0.25*(Q(I-1,J,K)+Q(I,J,K))
     .     *(XMFL3D(I,J,K-1)+XMFL3D(I,J,K))
 20   YFLUX(I,J,K)=0.25*(Q(I,J-1,K)+Q(I,J,K))
     .     *(YMFL3D(I,J,K-1)+YMFL3D(I,J,K))
C
C-------- HORIZONTAL DIFFUSION -----------------------------------------
C-------- ADD DIFFUSIVE FLUXES -----------------------------------------
      DO 30 K=2,KBM1
      DO 30 J=2,JM
      DO 30 I=2,IM
      AAMAX(I,J,K)=.5*(AAM(I,J,K)+AAM(I-1,J,K))
      AAMAY(I,J,K)=.5*(AAM(I,J,K)+AAM(I,J-1,K))
   30 CONTINUE
C
      DO 35 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 36 K=1,KBM1
   36         AAMAX(IC,JC,K)=0.0
            ELSE
              DO 37 K=1,KBM1
   37         AAMAX(ID,JD,K)=0.0
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 38 K=1,KBM1
   38         AAMAY(IC,JC,K)=0.0
            ELSE
              DO 39 K=1,KBM1
   39         AAMAY(ID,JD,K)=0.0
            ENDIF
      ENDIF
   35 CONTINUE
C
      DO 40 K=2,KBM1
      DO 40 J=2,JM
      DO 40 I=2,IM
      XFLUX(I,J,K)=XFLUX(I,J,K)
     .    -AAMAX(I,J,K)*(H(I,J)+H(I-1,J))*0.5*(H2(I,J)+H2(I-1,J))
     .    *(QB(I,J,K)-QB(I-1,J,K))*DUM(I,J)/(H1(I,J)+H1(I-1,J))
      YFLUX(I,J,K)=YFLUX(I,J,K)
     .    -AAMAY(I,J,K)*(H(I,J)+H(I,J-1))*0.5*(H1(I,J)+H1(I,J-1))
     .    *(QB(I,J,K)-QB(I,J-1,K))*DVM(I,J)/(H2(I,J)+H2(I,J-1))
   40 CONTINUE
C
C-------- VERTICAL ADVECTION -------------------------------------------
C-------- ADD FLUX TERMS & THEN STEP FORWARD IN TIME -------------------
      DO 50 K=2,KBM1
      DO 50 J=2,JMM1
      DO 50 I=2,IMM1
      QF(I,J,K)=(W(I,J,K-1)*Q(I,J,K-1)-W(I,J,K+1)*Q(I,J,K+1))
     .                     /(DZ(K)+DZ(K-1))*ART(I,J)
     .                      +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .                      +YFLUX(I,J+1,K)-YFLUX(I,J,K)
   50 QF(I,J,K)=((H(I,J)+ETB(I,J))*ART(I,J)*QB(I,J,K)-DTI2*QF(I,J,K))
     .             /((H(I,J)+ETF(I,J))*ART(I,J))
C
      RETURN
      END
