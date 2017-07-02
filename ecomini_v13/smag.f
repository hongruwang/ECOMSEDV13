      SUBROUTINE SMAG
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
      DIMENSION IVAR(IM,JM)
      DIMENSION DX1(IM,JM),DX2(IM,JM),DX3(IM,JM),DX4(IM,JM),
     .          DY1(IM,JM),DY2(IM,JM),DY3(IM,JM),DY4(IM,JM)
      EQUIVALENCE (DX1,C11),(DX2,C12),(DX3,C21),(DX4,C22)
      EQUIVALENCE (DY1,A(1,1,1)),(DY2,A(1,1,2)),
     .            (DY3,A(1,1,3)),(DY4,A(1,1,4))
C
C----------------------------------------------------------------------|
C         0.01 < horcon < 0.20                                         |
C         HOR. VISC. = CONST.*H1*H2*SQRT((2.*DU/H1)**2+(2.*DV/H2)**2   |
C                  +2.*(DU/H2+DV/H1)**2)                               |
C----------------------------------------------------------------------|
C
      DO 10 J=1,JM
      DO 10 I=1,IM
      DX1(I,J)=0.
      DX2(I,J)=0.
      DX3(I,J)=0.
      DX4(I,J)=0.
      DY1(I,J)=0.
      DY2(I,J)=0.
      DY3(I,J)=0.
 10   DY4(I,J)=0.
C
      DO 95 J=2,JMM1
      DO 95 I=2,IMM1
C 
      IF (FSM(I,J+1). EQ. 1.) THEN
      DY1(I,J)=0.5*(H2(I,J) + H2(I,J+1))
      ELSE
      DY1(I,J)=H2(I,J)
      END IF
C
      IF (FSM(I,J-1). EQ. 1.) THEN
      DY2(I,J)=0.5*(H2(I,J) + H2(I,J-1))
      ELSE
      DY2(I,J)=H2(I,J)
      END IF
C
      IF (FSM(I+1,J). EQ. 1. .AND. FSM(I+1,J+1).EQ. 1. ) THEN
      DY3(I,J)=0.5*(H2(I+1,J) + H2(I+1,J+1))
      ELSE IF (FSM(I+1,J). EQ. 1. .AND. FSM(I+1,J+1).EQ. 0.0) THEN
      DY3(I,J)=H2(I+1,J)
      ELSE IF (FSM(I+1,J). EQ. 0.0 .AND. FSM(I+1,J+1).EQ. 1.) THEN
      DY3(I,J)=H2(I+1,J+1)
      ELSE IF (FSM(I+1,J). EQ. 0.0 .AND. FSM(I+1,J+1).EQ. 0.0) THEN
      DY3(I,J)=H2(I,J)
      END IF
C
      IF (FSM(I+1,J). EQ. 1. .AND. FSM(I+1,J-1).EQ. 1.) THEN
      DY4(I,J)=0.5*(H2(I+1,J) + H2(I+1,J-1))
      ELSE IF (FSM(I+1,J). EQ. 1. .AND. FSM(I+1,J-1).EQ.0.0) THEN
      DY4(I,J)=H2(I+1,J)
      ELSE IF (FSM(I+1,J). EQ. 0.0 .AND. FSM(I+1,J-1).EQ. 1.) THEN
      DY4(I,J)=H2(I+1,J-1)
      ELSE IF (FSM(I+1,J). EQ. 0.0 .AND. FSM(I+1,J-1).EQ.0.0) THEN
      DY4(I,J)=H2(I,J)
      END IF
C
      IF (FSM(I,J+1). EQ. 1. .AND. FSM(I-1,J+1).EQ.1.) THEN
      DX1(I,J)=0.5*(H1(I,J+1) + H1(I-1,J+1))
      ELSE IF (FSM(I,J+1). EQ. 1. .AND. FSM(I-1,J+1).EQ.0.0) THEN
      DX1(I,J)=H1(I,J+1)
      ELSE IF (FSM(I,J+1). EQ. 0. .AND. FSM(I-1,J+1).EQ. 1.) THEN
      DX1(I,J)=H1(I-1,J+1)
      ELSE IF (FSM(I,J+1). EQ. 0. .AND. FSM(I-1,J+1).EQ.0.) THEN
      DX1(I,J)=H1(I,J)
      END IF
C
      IF (FSM(I+1,J+1). EQ. 1. .AND. FSM(I,J+1).EQ. 1.) THEN
      DX2(I,J)=0.5*(H1(I,J+1) + H1(I+1,J+1))
      ELSE IF (FSM(I+1,J+1). EQ. 1. .AND. FSM(I,J+1).EQ.0.0) THEN
      DX2(I,J)=H1(I+1,J+1)
      ELSE IF (FSM(I+1,J+1). EQ. 0. .AND. FSM(I,J+1).EQ. 1.) THEN
      DX2(I,J)=H1(I,J+1)
      ELSE IF (FSM(I+1,J+1). EQ. 0. .AND. FSM(I,J+1).EQ. 0.) THEN
      DX2(I,J)=H1(I,J)
      END IF
C
      IF (FSM(I-1,J). EQ. 1. ) THEN
      DX3(I,J)=0.5*(H1(I-1,J) + H1(I,J))
      ELSE
      DX3(I,J)=H1(I,J)
      END IF
C
      IF (FSM(I+1,J). EQ. 1. ) THEN
      DX4(I,J)=0.5*(H1(I+1,J) + H1(I,J))
      ELSE
      DX4(I,J)=H1(I,J)
      END IF
C
  95  CONTINUE
C
      IF(TOR.NE.'BAROTROPIC') THEN
C
      DO 110 K=1,KBM1
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
      AAM(I,J,K)=HORCON*ART(I,J)*FSM(I,J)
     1 *SQRT(4.0*(((U(I+1,J,K)-U(I,J,K))/H1(I,J))**2+((V(I,J+1,K)-
     2 V(I,J,K))/H2(I,J))**2) + 0.125*(ABS(U(I,J+1,K)-U(I,J,K))/DY1(I,J)
     3 +ABS(U(I,J,K)-U(I,J-1,K))/DY2(I,J) + ABS(U(I+1,J+1,K)-
     4 U(I+1,J,K))/DY3(I,J) + ABS(U(I+1,J,K)-U(I+1,J-1,K))/DY4(I,J)
     5 + ABS(V(I,J+1,K)-V(I-1,J+1,K))/DX1(I,J) + ABS(V(I+1,J+1,K)-V(I,
     6 J+1,K))/DX2(I,J) + ABS(V(I-1,J,K)-V(I,J,K))/DX3(I,J) + ABS(V(I+1, 
     7 J,K)-V(I,J,K))/DX4(I,J))**2)
 110  CONTINUE
C
      DO 120 N=1,NUMQBC
      DO 120 K=1,KBM1
      IC=IQC(N)
      JC=JQC(N)
      AAM(IC,JC,K)=0.0
  120 CONTINUE
C
      DO 121 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 121 K=1,KBM1
      AAM(IE,JE,K)=AAM(IC,JC,K)
  121 CONTINUE
C 
      DO 96 J=1,JM
      DO 96 I=1,IM
 96   AAM2D(I,J)=0.0
C
      DO 100 K=1,KBM1
      DO 100 J=1,JM
      DO 100 I=1,IM
      AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
 100  CONTINUE
C
      ELSE
C
      DO 130 J=2,JMM1
      DO 130 I=2,IMM1
 130  AAM2D(I,J)=HORCON*ART(I,J)*FSM(I,J)
     1 *SQRT(4.0*(((UA(I+1,J)-UA(I,J))/H1(I,J))**2+((VA(I,J+1)-
     2 VA(I,J))/H2(I,J))**2) + 0.125*(ABS(UA(I,J+1)-UA(I,J))/DY1(I,J)
     3 +ABS(UA(I,J)-UA(I,J-1))/DY2(I,J) + ABS(UA(I+1,J+1)-
     4 UA(I+1,J))/DY3(I,J) + ABS(UA(I+1,J)-UA(I+1,J-1))/DY4(I,J)
     5 + ABS(VA(I,J+1)-VA(I-1,J+1))/DX1(I,J) + ABS(VA(I+1,J+1)-VA(I,
     6 J+1))/DX2(I,J) + ABS(VA(I-1,J)-VA(I,J))/DX3(I,J) + ABS(VA(I+1,
     7 J)-VA(I,J))/DX4(I,J))**2)
C
      DO 112 N=1,NUMQBC
      IC=IQC(N)
      JC=JQC(N)
      AAM2D(IC,JC)=0.0
 112  CONTINUE
C
      DO 113 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      AAM2D(IE,JE)=AAM2D(IC,JC)
 113  CONTINUE
C
      ENDIF
C
      DO 11 J=1,JM
      DO 11 I=1,IM
      DX1(I,J)=0.
      DX2(I,J)=0.
      DX3(I,J)=0.
      DX4(I,J)=0.
      DY1(I,J)=0.
      DY2(I,J)=0.
      DY3(I,J)=0.
 11   DY4(I,J)=0.
C
      RETURN
      END
