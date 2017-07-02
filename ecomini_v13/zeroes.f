      SUBROUTINE ZEROES(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
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
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
C
	alat=-999.9
	alon=-999.9
C
      DO 250 K=1,KB
      DO 250 J=1,JM
      DO 250 I=1,IM
      T(I,J,K)=0.0
      TB(I,J,K)=0.0
      S(I,J,K)=0.0
      SB(I,J,K)=0.0
      RHO(I,J,K)=0.0
      RMEAN(I,J,K)=0.0
      TMEAN(I,J,K)=0.0
      SMEAN(I,J,K)=0.0
      XMFL3D(I,J,K)=0.0
      YMFL3D(I,J,K)=0.0
      A(I,J,K)=0.0
      C(I,J,K)=0.0
      VH(I,J,K)=0.0
      VHP(I,J,K)=0.0
      PROD(I,J,K)=0.0
      DTEF(I,J,K)=0.0
      U(I,J,K)=0.0
      UB(I,J,K)=0.0
      V(I,J,K)=0.0
      VB(I,J,K)=0.0
      UF(I,J,K)=0.0
      VF(I,J,K)=0.0
      W(I,J,K)=0.0
      DRHOX(I,J,K)=0.0
      DRHOY(I,J,K)=0.0
      Q2B(I,J,K)=0.0
      Q2(I,J,K)=0.0
      Q2LB(I,J,K)=0.0
      Q2L(I,J,K)=0.0
      L(I,J,K)=0.0
      KH(I,J,K)=0.0
      KM(I,J,K)=0.0
      KQ(I,J,K)=0.0
C
        CONC1(I,J,K)=0.0
        CONC1B(I,J,K)=0.0
        ARCC1(I,J,K)=0.0
        CMEAN1(I,J,K)=0.0
C
        CSED1(I,J,K)=0.0
        CSED2(I,J,K)=0.0
        CHEM1(I,J,K)=0.0
        CHEM2(I,J,K)=0.0
C
        WCT1BOT(I,J)=0.0
        WCT2BOT(I,J)=0.0
        CHEMBOT1(I,J)=0.0
        CHEMBOT2(I,J)=0.0
C
         ARCSED1(I,J,K)=0.0
         ARCSED2(I,J,K)=0.0
         ARCTAU(I,J,K)=0.0
         TAU(I,J,K)=0.0
  250 CONTINUE
      DO 640 J=1,JM
      DO 640 I=1,IM
      UAF(I,J)=0.0
      UA(I,J)=0.0
      UAB(I,J)=0.0
      VAF(I,J)=0.0
      VA(I,J)=0.0
      VAB(I,J)=0.0
      ELF(I,J)=0.0
      EL(I,J)=0.0
      ELB(I,J)=0.0
      TPS(I,J)=0.0
      ETF(I,J)=0.0
      ET (I,J)=0.0
      ETB(I,J)=0.0
      UTF(I,J)=0.0
      UTB(I,J)=0.0
      VTF(I,J)=0.0
      VTB(I,J)=0.0
      EGF(I,J)=0.0
      EGB(I,J)=0.0
      PATM(I,J)=0.0
      DPATM(I,J,1)=0.0
      DPATM(I,J,2)=0.0
      WTSURF(I,J)=0.0
      WSSURF(I,J)=0.0
      WUSURF(I,J)=0.0
      WVSURF(I,J)=0.0
      WUBOT(I,J)=0.0
      WVBOT(I,J)=0.0
      TRNU(I,J)=0.0
      TRNV(I,J)=0.0
      CURV42D(I,J)=0.0
      ADVUA(I,J)=0.0
      ADVVA(I,J)=0.0
      ADVUU(I,J)=0.0
      ADVVV(I,J)=0.0
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
C
        WCSURF(I,J)=0.0
C
  640 CONTINUE
C
      DO 100 K=1,KB
      DO 100 J=1,JM
      DO 100 I=1,IM
      ARCU (I,J,K)=0.0
      ARCV (I,J,K)=0.0
      ARCUX(I,J,K)=0.0
      ARCVX(I,J,K)=0.0
      ARCS (I,J,K)=0.0
      ARCT (I,J,K)=0.0
      ARCW (I,J,K)=0.0
      ARCKH(I,J,K)=0.0
      ARCHEM1(I,J,K)=0.0
      ARCHEM2(I,J,K)=0.0
      ARCSED1(I,J,K)=0.0
      ARCSED2(I,J,K)=0.0
      ARCTAU(I,J,K)=0.0
         TAU(I,J,K)=0.0
 100  CONTINUE
      DO 110 J=1,JM
      DO 110 I=1,IM
 110  ARCET (I,J)=0.0
C
      DO 200 N=1,EPTS
      ESAVE(N)=0.0
 200  CONTINUE
      DO 210 N=1,VPTS
      DZSAVE(N)=0.0
      DO 210 K=1,KB
      UZSAVE(N,K)=0.0
      VZSAVE(N,K)=0.0
      SZSAVE(N,K)=0.0
      TZSAVE(N,K)=0.0
C
        C1SAVE(N,K)=0.0
        C2SAVE(N,K)=0.0
        THSAVE(N)  =0.0
        TAUSAVE(N,K)=0.0
        P1SAVE(N,K)=0.0
        P2SAVE(N,K)=0.0
C
      C1ZSAVE(N,K)=0.0
C
 210  CONTINUE
      DO 220 K=1,KB
      DO 220 N=1,FPTS
      CCFLUX(N,K)=0.0
 220  CONTINUE
      ESUM  =0.0
      TKE   =0.0
      APE   =0.0
      VOLUME=0.0
      VSTOR =0.0
      EM    =0.0
      APEC  =0.0
      TSUM  =0.0
      SSUM  =0.0
C
C ------- SET CONSPLT, CONSTSR & CONSTRN TO TRUE FOR FIRST RUN THROUGH --------
      CONSPLT=.TRUE.
      CONSTSR=.TRUE.
      CONSTRN=.TRUE.
C
C  FOR PARTICLE TRACKING
C
      IF (PARTICLE.EQ.'INCLUDE') THEN
       DO 1110 LL=1,NSOURCE
        DO 1110 MM=1,NPART
          DO 1110 NN=1,NCONV
            XP(LL,MM,NN)=0.0
            YP(LL,MM,NN)=0.0
            ZP(LL,MM,NN)=0.0
            IP(LL,MM,NN)=0
            JP(LL,MM,NN)=0
            KP(LL,MM,NN)=0
            INOUT(LL,MM,NN)=0
 1110   CONTINUE
C
        DO 1120 K=1,KBM1
          DO 1130 N=1,DBCM
            PZERO1(N,K)=0.0
 1130     CONTINUE
          DO 1140 N=1,EBCM
            PZERO3(N,K)=0.0
 1140     CONTINUE
 1120   CONTINUE
C
        DO 1150 N=1,QBCM
          PZERO2(N)=0.0
 1150   CONTINUE
      ENDIF
C
      RETURN
      END
