      SUBROUTINE PROFT(F,WFSURF,DT2,SW)
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
      REAL*8 AF(IM,JM,KB),CF(IM,JM,KB),VHF(IM,JM,KB),VHPF(IM,JM,KB)
      REAL*8 FF(IM,JM,KB)
C
      DIMENSION F(IM,JM,KB),WFSURF(IM,JM),DH(IM,JM)
      DIMENSION SW(IM,JM)
      EQUIVALENCE (A,AF),(PROD,CF),(VH,VHF),(AF,FF)
      EQUIVALENCE (TPS,DH)
C
      UMOLPR=UMOL
C
C       This Version of PROFT is modified to accomodate Shortwave Radiation
C       Penetrated through the Water Column. The Algorithm uses the Classi-
C       fication of Jerlov (1976), but Approximate his Methodology Such that
C       Short Wave Radiation, absorved in the first few meters, is added
C       to the surface boundary condition, and the reminder is attenuated
C       according to RAD = SW*TRC*exp(EXTC*Z)
C
C       Here: RAD is attenuated shortwave radiation
C             SW is shortwave radiation (SW=0.0 when PROFT called for Salinity)
C             TR fraction absorved in surface layer
C             TRC fraction absorved in water column (1.-TR)
C             EXTC Extinction Coefficient
C       TR and EXTC are read from MET DATA segment of run_data input file
C
C
C       FRACTION OF SHORT WAVE RAD PENETRATED THRU WATER COLUMN "TRC"
        TRC = 1. - TR
C
      UMOLPR = UMOL
C
C Note: while computing Salinity SW is zero
C
        DO 25 J = 1, JM
        DO 25 I = 1, IM
            WFSURF(I,J) = WFSURF(I,J)-SW(I,J)
   25   CONTINUE

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
      DO 100 K=2,KBM1
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      AF(I,J,K-1)=-DT2*(KH(I,J,K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*DH(I,J)
     .     *DH(I,J))
      CF(I,J,K)=-DT2*(KH(I,J,K)+UMOLPR)/(DZ(K)*DZZ(K-1)*DH(I,J)
     .     *DH(I,J))
 100  CONTINUE
C
C-------- SURFACE BOUNDARY CONDITIONS - WFSURF -------------------------
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
      VHF(I,J,1)=AF(I,J,1)/(AF(I,J,1)-1.)
C  SW Penetration **********************
            VHPF(I,J,1) = -DT2 * (WFSURF(I,J)+SW(I,J)*(1.-TRC))
     *      /(-DZ(1)*DH(I,J)) - F(I,J,1)
 110  VHPF(I,J,1)=VHPF(I,J,1)/(AF(I,J,1)-1.)
C
C-------- SURFACE BOUNDARY CONDITIONS - FSURF --------------------------
      DO 101 K=2,KBM2
      DO 101 J=2,JMM1
      DO 101 I=2,IMM1
      VHPF(I,J,K)=1./(AF(I,J,K)+CF(I,J,K)*(1.-VHF(I,J,K-1))-1.)
      VHF(I,J,K)=AF(I,J,K)*VHPF(I,J,K)
      VHPF(I,J,K)=(CF(I,J,K)*VHPF(I,J,K-1)-DBLE(F(I,J,K)))*VHPF(I,J,K)
 101  CONTINUE
C
      DO 130 K=1,KB
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  FF(I,J,K)=F(I,J,K)
C
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
 102  FF(I,J,KBM1)=((CF(I,J,KBM1)*VHPF(I,J,KBM2)-FF(I,J,KBM1))
     .   /(CF(I,J,KBM1)*(1.-VHF(I,J,KBM2))-1.))
C
      DO 105 K=2,KBM1
      KI=KB-K
      DO 105 J=2,JMM1
      DO 105 I=2,IMM1
      FF(I,J,KI)=(VHF(I,J,KI)*FF(I,J,KI+1)+VHPF(I,J,KI))
 105  CONTINUE
C
C       Penetrative Radiation Calculation. Any unattenuated is deposited
C       in the Bottom Layer
C
      DO 320 K=1,KBM1
      DO 320 J=1,JM
      DO 320 I=1,IM
          RADP(I,J,K)=SW(I,J)*TRC*
     *               EXP(EXTC*Z(K)*DH(I,J))
  320     CONTINUE
C
      DO 420 J=1,JM
      DO 420 I=1,IM
  420 RADP(I,J,KB)=0.
C
      Do 220 K = 1, KBM1
        Do 210 J = 1, JM
          Do 200 I = 1, IM
            IF(DH(I,J).GT.0.0)THEN
              F(I,J,K) = FF(I,J,K)-DT2*(RADP(I,J,K)-RADP(I,J,K+1))/
     *                   (DZ(K)*DH(I,J))
            ENDIF
  200     Continue
  210   Continue
  220 Continue
C
C
      RETURN
      END
