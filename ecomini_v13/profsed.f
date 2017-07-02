      SUBROUTINE PROFSED(F,DT2,WSET,WFBOT)

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
      REAL*8 AD(IM,JM,KB),AL(IM,JM,KB),AU(IM,JM,KB)
      REAL*8 FF(IM,JM,KB),VHF(IM,JM,KB),VHPF(IM,JM,KB)
      DIMENSION WFBOT(IM,JM),WSET(IM,JM,KB)
C
      DIMENSION F(IM,JM,KB),DH(IM,JM)
      EQUIVALENCE (A,AD),(PROD,AL),(VH,VHF),(AD,FF)
      EQUIVALENCE (TPS,DH)

      IF (TOR.EQ.'BAROTROPIC') RETURN
C
      UMOLPR=UMOL
C
C-----------------------------------------------------------------------
C
C        THE FOLLOWING SECTION SOLVES THE EQUATION
C         DT2*(KH*F')'-F=-FB
C
C-----------------------------------------------------------------------
C
      DO 90 J=2,JMM1
      DO 90 I=2,IMM1
      DH(I,J)=H(I,J)+ETF(I,J)
  90  CONTINUE
C
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      IF(FSM(I,J).LE.0.) GO TO 100
      AL(I,J,1)=0.0
      AD(I,J,1)=1.0 + DT2*(KH(I,J,2)+UMOLPR)/(DZZ(1)*DH(I,J)
     .     *DZ(1)*DH(I,J)) + DT2*WSET(I,J,2)/(DZ(1)*DH(I,J))
      AU(I,J,1)= - DT2*(KH(I,J,2)+UMOLPR)/(DZZ(1)*DH(I,J)
     .            *DZ(1)*DH(I,J))

      DO 95 K=2,KBM1-1
      AL(I,J,K)= - DT2*(KH(I,J,K)+UMOLPR)/(DZZ(K-1)*DH(I,J)
     .     *DZ(K)*DH(I,J)) - DT2*WSET(I,J,K)/(DZ(K)*DH(I,J))
      AD(I,J,K)=1.0 + DT2*(KH(I,J,K+1)+UMOLPR)/(DZZ(K)*DH(I,J)
     .     *DZ(K)*DH(I,J)) + DT2*WSET(I,J,K+1)/(DZ(K)*DH(I,J))
     .     + DT2*(KH(I,J,K)+UMOLPR)/(DZZ(K-1)*DH(I,J)*DZ(K)*DH(I,J)) 
      AU(I,J,K)= - DT2*(KH(I,J,K+1)+UMOLPR)/(DZZ(K)*DH(I,J)*
     .             DZ(K)*DH(I,J))
  95  CONTINUE

      AL(I,J,KBM1)= - DT2*(KH(I,J,KBM1)+UMOLPR)/(DZZ(KBM1-1)*DH(I,J)
     .     *DZ(KBM1)*DH(I,J)) - DT2*WSET(I,J,KBM1)/(DZ(K)*DH(I,J))
      AD(I,J,KBM1)=1.0 + DT2*(KH(I,J,KBM1)+UMOLPR)/(DZZ(KBM1-1)*DH(I,J)
     .     *DZ(KBM1-1)*DH(I,J)) 
      AU(I,J,KBM1)=0.0

 100  CONTINUE
C
C  FORWARD SWEEP OF TRIDIAGONAL SCHEME

C-------- SURFACE BOUNDARY CONDITIONS - WFSURF -------------------------
C
C  ASSUME ZERO SURFACE FLUX FOR SEDIMENT
C
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
      IF(FSM(I,J).LE.0.) GO TO 110
      VHF(I,J,1)=AD(I,J,1)
      VHPF(I,J,1)=F(I,J,1)/VHF(I,J,1)
  
      DO 101 K=2,KBM1
      VHF(I,J,K)=AD(I,J,K)-AL(I,J,K)*AU(I,J,K-1)/VHF(I,J,K-1)
      IF(K.LT.KBM1)THEN
        VHPF(I,J,K)=(F(I,J,K)-AL(I,J,K)*VHPF(I,J,K-1))/VHF(I,J,K)
      ELSE
        IF(WFBOT(I,J).LT.0.0)THEN
          WFBOTMAX=100.*F(I,J,KBM1)*DH(I,J)*DZ(KBM1)
          WFBOT(I,J)=AMAX1(WFBOT(I,J),-WFBOTMAX)
        ENDIF
C
        VHPF(I,J,K)=(F(I,J,K)+ (2.0*WFBOT(I,J)/(DZ(K)*DH(I,J)*100.))
     .             -AL(I,J,K)*VHPF(I,J,K-1))/VHF(I,J,K)
        IF (VHPF(I,J,K).LE. 0.0 )THEN
          IF (WFBOT(I,J).LE. 0.0 ) WFBOT(I,J)=0.0
          VHPF(I,J,K)=0.0  !! loss all mass
        ENDIF
      ENDIF
 101  CONTINUE
 110  CONTINUE
C
      DO 120 K=1,KB
      DO 120 J=1,JM
      DO 120 I=1,IM
 120  FF(I,J,K)=F(I,J,K)
C
C  BACKWARD SWEEP FOR SOLUTION

      DO 130 J=1,JM
      DO 130 I=1,IM
        IF(FSM(I,J).LE.0.) GO TO 130
        FF(I,J,KBM1)=VHPF(I,J,KBM1)

        DO 105 K=2,KBM1
          KI=KB-K
          FF(I,J,KI)=VHPF(I,J,KI)-AU(I,J,KI)*FF(I,J,KI+1)/VHF(I,J,KI)
 105  CONTINUE
 130  CONTINUE
C
      DO 140 K=1,KB
      DO 140 J=1,JM
      DO 140 I=1,IM
      IF(FSM(I,J).LE.0.) GO TO 140
      F(I,J,K)=FF(I,J,K)
 140  CONTINUE
C
      RETURN
      END
