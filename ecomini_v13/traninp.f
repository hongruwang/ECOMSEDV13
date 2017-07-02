      SUBROUTINE TRANINP(IFLAG)
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
C  READS HYDRO INPUT FROM hqi_tran FILE AND PREPARES FOR SEDIMENT
C  TRANSPORT CALCUATIONS
C
C  FOR COUPLING WITH POM AND WET GRID OUTPUT ONLY
C
C**********************************************************************
C
      INCLUDE 'comdeck'
C
C     DECLARE ARRAYS FOR WET GRID OUTPUT
C
      CHARACTER*10 RESTAR
      COMMON /COAST1/ICNT,INDX(MAXWET),JNDX(MAXWET),RESTAR
C 
      DIMENSION  WETGU(MAXWET,KBM1),  WETGV(MAXWET,KBM1),
     .           WETGW(MAXWET,KB),    
     .           WETGAAMX(MAXWET,KBM1), WETGAAMY(MAXWET,KBM1),
     .           WETGKM(MAXWET,KB),   WETGKH(MAXWET,KB),
     .           WETGES(MAXWET),WETGED(MAXWET), 
     .           WETGS(MAXWET,KBM1),WETGT(MAXWET,KBM1) 
C
C
C
      DIMENSION  KHLPF(IM,JM,KB),
     .   SLPF(IM,JM,KB),    TLPF(IM,JM,KB),
     .   ULPF(IM,JM,KB),    VLPF(IM,JM,KB),   WLPF(IM,JM,KB)
C
       REAL*8 ES(IM,JM),         ED(IM,JM)
C
      REAL KHLPF
C
C   CURRENTLY SEDIMENT TRANSPORT MODEL ONLY USES wet_grid TRANSPORT
C   INPUT STRUCTURE HAS BEEN MODIFIED BASED ON CONVENTIONAL gcm_tran
C
      READ (IUTRN) TMIDDLE
      READ (IUTRN) ((WETGU(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGV(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGW(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGAAMX(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGAAMY(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGKH(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGKM(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) (WETGES(I),I=1,ICNT)
      READ (IUTRN) (WETGED(I),I=1,ICNT)
      READ (IUTRN) ((WETGS(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGT(I,K),I=1,ICNT),K=1,KBM1)
C
      DO 330 K=1,KB
        DO 330 I=1,ICNT
          WLPF(INDX(I),JNDX(I),K)=WETGW(I,K)
          KH(INDX(I),JNDX(I),K)=WETGKH(I,K)
          KM(INDX(I),JNDX(I),K)=WETGKM(I,K)
 330  CONTINUE
C
      DO 340 K=1,KBM1
        DO 340 I=1,ICNT
          ULPF(INDX(I),JNDX(I),K)=WETGU(I,K)
          VLPF(INDX(I),JNDX(I),K)=WETGV(I,K)
          AAMAX(INDX(I),JNDX(I),K)=WETGAAMX(I,K)
          AAMAY(INDX(I),JNDX(I),K)=WETGAAMY(I,K)
          T(INDX(I),JNDX(I),K)    =WETGT(I,K)
          S(INDX(I),JNDX(I),K)    =WETGS(I,K)
 340   CONTINUE
C
      DO 350 I=1,ICNT
        ES(INDX(I),JNDX(I))=WETGES(I)
        ED(INDX(I),JNDX(I))=WETGED(I)
 350  CONTINUE
C
C  SET ELEVATION TIME RATE OF CHANGE & INITIAL DEPTH
C
      DO 40 I=1,IM
        DO 40 J=1,JM
C
C  START OF CALC. FOR IFLAG=0
C
          IF (IFLAG.EQ.0) THEN
            ETB(I,J)=ES(I,J)
            ET(I,J)=ES(I,J)
            ETF(I,J)=ES(I,J)
            REWIND (IUTRN)
          ENDIF
C
          EL(I,J)=ES(I,J)
          DT(I,J)=H(I,J)+ES(I,J)
          DETA(I,J)=ED(I,J)
 40   CONTINUE
C
C-----------------------------------------------------------------------
C         CALCULATE HORIZONTAL MASS FLUXES, (H2*U*D) AND (H1*V*D)
C-----------------------------------------------------------------------
C
      DO 10 K=1,KBM1
        DO 10 J=2,JMM1
          DO 10 I=2,IM
            XMFL3D(I,J,K)=ULPF(I,J,K)/DZ(K)
C
C  CALC. VELOCITIES
C
            U(I,J,K)=0.0
            DBAR=0.25*(DT(I,J)+DT(I-1,J))*(H2(I,J)+H2(I-1,J))
            IF (DBAR.GT.0.0) U(I,J,K)=ULPF(I,J,K)/(DZ(K)*DBAR)
 10   CONTINUE
C
      DO 20 K=1,KBM1
        DO 20 J=2,JM
          DO 20 I=2,IMM1
            YMFL3D(I,J,K)=VLPF(I,J,K)/DZ(K)
C
C  CALC. VELOCITIES
C
            V(I,J,K)=0.0
            DBAR=0.25*(DT(I,J)+DT(I,J-1))*(H1(I,J)+H1(I,J-1))
            IF (DBAR.GT.0.0) V(I,J,K)=VLPF(I,J,K)/(DZ(K)*DBAR)
 20   CONTINUE
C
C-----------------------------------------------------------------------
C
C  SET VERTICAL VELOCITIES
C
      DO 30 K=1,KB
        DO 30 J=1,JM
          DO 30 I=1,IM
            W(I,J,K)=WLPF(I,J,K)/ART(I,J) 
 30   CONTINUE
C
      RETURN
      END
