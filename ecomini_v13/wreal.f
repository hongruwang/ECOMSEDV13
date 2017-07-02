      SUBROUTINE WREAL(DTI2)
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
      DIMENSION DXR(IM,JM),DXL(IM,JM),DYT(IM,JM),DYB(IM,JM)
C
C-------- CALCULATE REAL VERTICAL VELOCITY
      DO 100 K=1,KB
      DO 100 J=1,JM
      DO 100 I=1,IM
 100  WR(I,J,K) = 0.
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
      DXR(I,J)=2./(H1(I+1,J)+H1(I,J))
      DXL(I,J)=2./(H1(I,J)+H1(I-1,J))
      DYT(I,J)=2./(H2(I,J+1)+H2(I,J))
 110  DYB(I,J)=2./(H2(I,J)+H2(I,J-1))
C
      DO 710 K=1,KBM1
C
      DO 711 J=1,JM
      DO 711 I=1,IM
 711  TPS(I,J)=ZZ(K)*DT(I,J) + ET(I,J)
C
      DO 712 J=2,JMM1
      DO 712 I=2,IMM1
      WR(I,J,K)=0.5*(W(I,J,K)+W(I,J,K+1)) + 0.5*
     1         ( U(I+1,J,K)*(TPS(I+1,J)-TPS(I,J))*DXR(I,J)  +
     2           U(I,J,K)*(TPS(I,J)-TPS(I-1,J))*DXL(I,J)    +
     3           V(I,J+1,K)*(TPS(I,J+1)-TPS(I,J))*DYT(I,J)  +
     4           V(I,J,K)*(TPS(I,J)-TPS(I,J-1))*DYB(I,J)      )
     5      +     (1.+ZZ(K))*(ETF(I,J)-ETB(I,J))/DTI2
  712 CONTINUE
C
  710 CONTINUE
C
      DO 200 K=1,KB
      DO 200 I=1,IM       
      WR(I,1,K)=WR(I,2,K)
 200  WR(I,JM,K)=WR(I,JMM1,K)
      DO 210 K=1,KB
      DO 210 J=1,JM
      WR(1,J,K)=WR(2,J,K)
 210  WR(IM,J,K)=WR(IMM1,J,K)
C
      DO 800 K=1,KBM1
      DO 800 J=1,JM
      DO 800 I=1,IM
  800 WR(I,J,K)=FSM(I,J)*WR(I,J,K)
C
      RETURN
      END
