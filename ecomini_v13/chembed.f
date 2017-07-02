      SUBROUTINE CHEMBED(CHEMDRAT1,CHEMDRAT2)
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
C  CHEM BED MODEL
C
C  NOTE:  MUST USE OPPOSITE SIGN OF FLUXES FOR DEPOSITION TO BED
C
C************************************************************
C
      INCLUDE 'comdeck'
C
      REAL*8 CHEMM,SEDM,SUMTEM,TLAY1,TOTLAY,DM,FRAC,SEDTEMP
C
      CDEC1=EXP(-CHEMDRAT1*DTI*XNSBED)
      CDEC2=EXP(-CHEMDRAT2*DTI*XNSBED)
C
      DO 10 J=2,JMM1
        DO 10 I=2,IMM1
          IF (FSM(I,J).GT.0.0) THEN
C
            CHEMM=-XNSBED*(CHEMBOT1(I,J)+CHEMBOT2(I,J))
            SEDM=-XNSBED*(WCT1BOT(I,J)+WCT2BOT(I,J))
C
            IF (SEDM.EQ.0.0) GOTO 10
C
C  DEPOSITION 
C
            IF (SEDM.GT.0.0) THEN
C
C  ADD MASS OF CHEM TO LAYER 1 (ug CHEM/cm**2)
C
              CHEMMASS(1,I,J)=CHEMMASS(1,I,J)+CHEMM
C
C  ADD MASS OF SOLIDS TO LAYER 1 (g SOLIDS/cm**2)
C
              SEDMASS(1,I,J)=SEDMASS(1,I,J)+SEDM
C
C  CALC. THICKNESS OF LAYER 1 (cm)
C
              TLAY1=SEDMASS(1,I,J)/CBED(I,J)
C
C  MIX ACTIVE LAYER
C
              IF (NACTLAY.GT.1) THEN
C
C  TOTAL CHEM MASS IN ACTIVE LAYER (ug CHEM/cm**2)
C
                SUMTEM=0.0
                DO 20 N=1,NACTLAY
                  SUMTEM=SUMTEM+CHEMMASS(N,I,J)
 20             CONTINUE
C
C  TOTAL THICKNESS OF ACTIVE LAYER (cm)
C
                TOTLAY=TLAY1+FLOAT(NACTLAY-1)*CHEMTHIK
C
C  THICKNESS-WEIGHTED DISTRIBUTION OF MASS CHEM IN ACTIVE LAYER
C
                CHEMMASS(1,I,J)=TLAY1*SUMTEM/TOTLAY
C
                DO 30 N=2,NACTLAY
                  CHEMMASS(N,I,J)=CHEMTHIK*SUMTEM/TOTLAY
 30             CONTINUE
              ENDIF
C
C  CHECK IF LAYER 1 THICKNESS IS GREATER THAN BED LAYER THICKNESS
C
C  REORDER ALL LAYERS
C
              IF (SEDMASS(1,I,J).GT.(CBED(I,J)*CHEMTHIK)) THEN
C
                IF (NCHEMLAY.GT.2) THEN
                  CHEMMASS(NCHEMLAY,I,J)=CHEMMASS(NCHEMLAY,I,J)+
     +                                   CHEMMASS(NCHEMLAY-1,I,J)
                  SEDMASS(NCHEMLAY,I,J)=SEDMASS(NCHEMLAY,I,J)+
     +                                   SEDMASS(NCHEMLAY-1,I,J)
C
                  IF (NCHEMLAY.GT.3) THEN
                    DO 50 N=NCHEMLAY-1,3,-1
                      CHEMMASS(N,I,J)=CHEMMASS(N-1,I,J)
                      SEDMASS(N,I,J)=SEDMASS(N-1,I,J)
 50                 CONTINUE
                  ENDIF
C
C  CALC. FRACTION OF SEDIMENT LEFT IN LAYER 1
C
                  DM=SEDMASS(1,I,J)-(CBED(I,J)*CHEMTHIK)
                  FRAC=DM/SEDMASS(1,I,J)
C
C  RESET MASSES IN LAYERS 1 & 2
C
                  SEDMASS(1,I,J)=DM
                  SEDMASS(2,I,J)=CBED(I,J)*CHEMTHIK
C
                  CHEMMASS(2,I,J)=(1.-FRAC)*CHEMMASS(1,I,J)
                  CHEMMASS(1,I,J)=FRAC*CHEMMASS(1,I,J) 
                ELSE
C
C  TWO LAYER MODEL
C
                  DM=SEDMASS(1,I,J)-(CBED(I,J)*CHEMTHIK)
                  FRAC=DM/SEDMASS(1,I,J)
C
                  SEDMASS(1,I,J)=DM
                  SEDMASS(2,I,J)=SEDMASS(2,I,J)+CBED(I,J)*CHEMTHIK
C
                  CHEMMASS(2,I,J)=CHEMMASS(2,I,J)+(1.-FRAC)*
     +                   CHEMMASS(1,I,J)
                  CHEMMASS(1,I,J)=FRAC*CHEMMASS(1,I,J) 
                ENDIF
              ENDIF
C
C  EROSION
C
            ELSE
C
C  REMOVE MASS OF SOLIDS FROM LAYER 1 (g SOLIDS/cm**2)
C
              SEDTEMP=SEDMASS(1,I,J)+SEDM
C
C  LAYER 1 ERODED
C 
              IF (SEDTEMP.LE.0.0) THEN
C
C  REORDER LAYERS
C
C  FIND DEEPEST LAYER WITH MASS
C
                DO 60 N=2,NCHEMLAY
                  IF (SEDMASS(N,I,J).LE.0.0) THEN
                    NMAX=N-1
                    GOTO 70
                  ENDIF
 60             CONTINUE
                NMAX=NCHEMLAY
C
 70             SEDMASS(1,I,J)=SEDMASS(2,I,J)+SEDTEMP
                CHEMMASS(1,I,J)=CHEMMASS(2,I,J)+SEDTEMP*CBEDCHEM(2,I,J)
C
                IF (NMAX.GT.2) THEN
                  DO 80 N=2,NMAX-1
                    SEDMASS(N,I,J)=SEDMASS(N+1,I,J)
                    CHEMMASS(N,I,J)=CHEMMASS(N+1,I,J)
 80               CONTINUE
                ENDIF
C
                SEDMASS(NMAX,I,J)=0.0
                CHEMMASS(NMAX,I,J)=0.0
C
C  LAYER 1 NOT ERODED
C
              ELSE
                SEDMASS(1,I,J)=SEDTEMP
                CHEMMASS(1,I,J)=CHEMMASS(1,I,J)+CHEMM
              ENDIF
            ENDIF
C
C  RE-CALC. CHEM CONC. (ug CHEM/g SOLIDS)
C
C  INCLUDE FIRST-ORDER DECAY RATE
C
            IF (IBMSK(I,J).EQ.0) THEN
              CDEC=PSED1(LAYMAX,I,J)*CDEC1+PSED2(LAYMAX,I,J)*CDEC2
            ELSE
              CDEC=FPBED(1,I,J)*CDEC1+FPBED(2,I,J)*CDEC2
            ENDIF
C
            DO 90 N=1,NCHEMLAY
              IF (SEDMASS(N,I,J).GT.0.0) THEN
                IF(CHEMMASS(N,I,J).le.0.0) then
                   CHEMMASS(N,I,J)=0.0
                ENDIF
                CBEDCHEM(N,I,J)=CDEC*CHEMMASS(N,I,J)/SEDMASS(N,I,J)
              ELSE
                CBEDCHEM(N,I,J)=0.0
              ENDIF
 90         CONTINUE
          ENDIF
 10   CONTINUE
C
      RETURN
      END
