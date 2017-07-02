      SUBROUTINE SUSLOD
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
C  CALCULATES SUSPENDED LOAD OF NON-COHESIVE SEDIMENT
C  MKS UNITS USED
C  REF: VAN RIJN, 1984, J. HYDR. ENGR., p. 1613
C
C  USES 2 CLASSES OF SUSPENDED SEDIMENT
C    CLASS 1:  SILT/CLAY (MEDIUM CLASS)  (D < 75 microns)
C    CLASS 2:  FINE SAND (COARSE CLASS)  (75 < D < 450 microns)
C
C*********************************************************
C
      INCLUDE 'comdeck'
C
      REAL*8 DVS(3),TDV,EROTOT,DEPK(2)
C
      PIHALF=PI/2.
C
C  FOR CONSISTENCY WITH BOTTOM FRICTION FACTOR
C
      NIKURH=11.*Z0B
C
C  SUSPENDED LOAD
C
      DO 10 I=2,IM-1
        DO 30 J=2,JM-1
          IF (FSM(I,J).EQ.0.0.OR.IBMSK(I,J).EQ.0) GOTO 30
C
          WSED2=WSET2(I,J,KBM1)
C
C  CLASS 1, 2 FLUX
C
          WCT1BOT(I,J)=0.0
          WCT2BOT(I,J)=0.0
C
C  HARD BOTTOM
C
          IF (IBMSK(I,J).LT.0) GOTO 30
C
C  CALC. TOTAL VELOCITY
C
          UBAR(I,J,KBM1)=ABS(0.5*(U(I,J,KBM1)+U(I+1,J,KBM1)))
          VBAR(I,J,KBM1)=ABS(0.5*(V(I,J,KBM1)+V(I,J+1,KBM1)))
          UTOT=SQRT(UBAR(I,J,KBM1)*UBAR(I,J,KBM1)+
     +              VBAR(I,J,KBM1)*VBAR(I,J,KBM1)) 
C
          IF (UBAR(I,J,KBM1).GT.0.0) THEN
            THETAC=ATAN(VBAR(I,J,KBM1)/UBAR(I,J,KBM1))
          ELSE
            THETAC=PIHALF
          ENDIF
C
C  CALC. BED-SHEAR VELOCITY
C
          IF (WAVEDYN.EQ.'NEGLECT ') THEN
            UBED=SQRT(CBC(I,J))*UTOT
          ELSE
            UBED=USHEAR(I,J)
          ENDIF
C
C  CALC. USTAR (EQ. 17)
C
          USTAR=UBED
C
C  TEST FOR NO BED LOAD IN THIS ELEMENT
C
          IF (USTAR.LE.UCRVAR(I,J)) GOTO 40        
C
 35       CONTINUE
C
C  CALC. TRANSPORT STAGE PARAMETER (EQ. 2)
C
          TPAR=((USTAR*USTAR)/(UCRVAR(I,J)*UCRVAR(I,J)))-1.
C
C  CALC. REFERENCE CONCENTRATION (EQ. 38)
C
          AREF=MAX(0.01*DT(I,J),NIKURH)
          CREF=CSED(3)*(TPAR**1.5)/AREF
C
C  CALC. CRITICAL VELOCITY FOR SUSPENSION  (EQ. 5)
C
          UCRS=WSET2(I,J,KBM1)
C
          IF (UBED.LE.UCRS) GOTO 40
C
C  RATIO OF WSED TO UBED
C
          WURAT=WSED2/UBED
C
C  COMPUTE BETA (EQ. 22)
C
          BSED=1.+2.*WURAT*WURAT
C
C  COMPUTE PHI (EQ. 34)
C
          PHI=2.97*(WURAT**0.8)*(CREF**0.4)
C
C  COMPUTE Z' (EQ. 33)
C
          ZP1=2.5*WURAT/BSED+PHI
C
C  COMPUTE F FACTOR (EQ. 44)
C
          AD=AREF/DT(I,J)
          FSED=((AD**ZP1)-(AD**1.2))/(((1.-AD)**ZP1)*(1.2-ZP1))
C
C  CALC. WIDTH OF ELEMENT PERPENDICULAR TO FLOW
C  CURRENT ANGLE
C
          IF (UBAR(I,J,KBM1).GT.0.0) THEN
            THETAC=ATAN(VBAR(I,J,KBM1)/UBAR(I,J,KBM1))
          ELSE
            THETAC=PIHALF
          ENDIF
C
C  ELEMENT ANGLE
C
          PHI2=ATAN(H2(I,J)/H1(I,J))
C
          IF (THETAC.LE.(PIHALF-PHI2)) THEN
            DELTAS=UTOT*H2(I,J)/(UBAR(I,J,KBM1)+1.E-06)
          ELSE
            DELTAS=UTOT*H1(I,J)/(VBAR(I,J,KBM1)+1.E-06)
          ENDIF
C
C  COMPUTE TOTAL SUSPENDED LOAD  (EQ. 43) (m**2/s)
C
          QS=FSED*UTOT*DZ(KBM1)*DT(I,J)*CREF
C
C  CALC. TOTAL RESUSPENSION PER UNIT AREA (g/cm**2)
C  100 MULTIPLIER FOR UNIT CONSISTENCY
C
          CTOT=CSED1(I,J,KBM1)+CSED2(I,J,KBM1)
          CTOT=MAX(CTOT,0.0)
          EROTOT=((2.65*QS-UTOT*DZ(KBM1)*DT(I,J)*CTOT)
     +           *DTI*DELTAS/(H1(I,J)*H2(I,J)))*100.
          WCT1BOT(I,J)=FALAY(1,I,J)*EROTOT
          WCT2BOT(I,J)=FALAY(2,I,J)*EROTOT
          E(1,I,J)=WCT1BOT(I,J)
          E(2,I,J)=WCT2BOT(I,J)
          GOTO 45
C
C***********************************************************
C
C  DEPOSITION ONLY IF NO EROSION
C  DEPOSITION PER UNIT AREA (g/cm**2)
C  SETTLING SPEED MUST BE IN cm/s 
C  CLASS 2  (NON-COHESIVE)
C 
C  NO DEPOSITION WHEN BED THICKNESS > BATHYMETRIC DEPTH
C
 40       CONTINUE
        IF (SEDTHK(I,J).LT.H(I,J)) THEN
          IF (CSED2(I,J,KBM1).GT.0.0.AND.UBED.LT.WSED2) THEN
            DD(2,I,J)=CSED2(I,J,KBM1)*(100.*WSED2)*DTI
            DD(2,I,J)=DMAX1(DD(2,I,J),0.0D0)
            WCT2BOT(I,J)=-DD(2,I,J)
          ENDIF
C
C  CLASS 1 (COHESIVE)
C
C   PUT IN ACTUAL RELATIONSHIP FOR SETTLING INSTEAD OF LOOKUP TABLES
C
          IF (CSED1(I,J,KBM1).GT.0.0) THEN
            WSET1(I,J,KBM1)=ADEP*
     +        (1000000.0*CSED1(I,J,KBM1)*TAU(I,J,KBM1))
     +        **DEPEXP/10000.0/100.0 !SETTLING VEL IN m/s
              WSET1(I,J,KBM1)=AMAX1(WSET1(I,J,KBM1),WS1MIN) !IN m/s
C
            DD(1,I,J)=WSET1(I,J,KBM1)*100.0*CSED1(I,J,KBM1)*DTI
            DD(1,I,J)=DMAX1(DD(1,I,J),0.0D0)
C
C  NO DEPOSITION OF COMPONENT 1 IF TAU > TCRDEP
C  USE PROBABILITY OF DEPOSITION
C
            ITAU=NINT(100.*TAU(I,J,KB))+1.
            IF(ITAU.GT.20000) ITAU=20000
            DD(1,I,J)=PDEP(ITAU)*DD(1,I,J)
            WCT1BOT(I,J)=-DD(1,I,J)
          ENDIF
        ELSE
            DD(1,I,J)=0.0
            DD(2,I,J)=0.0
            WCT1BOT(I,J)=0.0
            WCT2BOT(I,J)=0.0
        END IF
C
C  CHANGE IN BED THICKNESS (m**3) AND VOL./MASS FRACTIONS
C
C  COHESIVE FRACTION
C
 45       DVS(1)=WCT1BOT(I,J)*ART(I,J)*XNSBED/265.
C
C  NON-COHESIVE FRACTION
C
          DVS(2)=WCT2BOT(I,J)*ART(I,J)*XNSBED/265.
C
C  CALC. TOTAL CHANGE IN  BED VOLUME
C  TDV > 0 :  EROSION
C  TDV < 0 :  DEPOSITION
C
          TDV=DVS(1)+DVS(2)
C
C  CHECK FOR MAX. RESUPSENSION
C
          DVS(3)=(FALAY(1,I,J)+FALAY(2,I,J))*ACTLAY(2,I,J)          
C
          IF (TDV.GT.DVS(3)) THEN
            TDV=DVS(3)
            DVS(1)=FALAY(1,I,J)*DVS(3)
            DVS(2)=FALAY(2,I,J)*DVS(3)
            WCT1BOT(I,J)=265.*DVS(1)/(XNSBED*ART(I,J))
            WCT2BOT(I,J)=265.*DVS(2)/(XNSBED*ART(I,J))
          ENDIF
          EROTOT=0.0
          DEPTOT=0.0
          IF (TDV.GT.0.0) EROTOT=TDV
          IF (TDV.LT.0.0) DEPTOT=-TDV
          DO 46 K=1,2
            IF (DVS(K).LT.0.0) DEPK(K)=-DVS(K)
 46       CONTINUE
C
C  CALCULATE NEW BED THICKNESS (VOLUME)
C
          BEDTH(1,I,J)=BEDTH(2,I,J)-TDV
C
C  CALC. NEW ACTIVE LAYER VOLUME:  METHOD 1 (NIEKERK et al., 1992)
C
          ACTLAY(1,I,J)=CARMOR(I,J)*TAU(I,J,KB)*SUSARM
C
C  CALC. CHANGE IN ACTIVE LAYER FRACTIONS
C
C  DEPOSITION TO PARENT BED
C
          IF (ACTLAY(1,I,J).EQ.0.0.AND.ACTLAY(2,I,J).EQ.0.0) THEN
            DO 50 K=1,2
              FPBED(K,I,J)=(FPBED(K,I,J)*(BEDTH(2,I,J)-EROTOT)+DEPK(K))
     +                                 /BEDTH(1,I,J)
 50         CONTINUE
            GOTO 80
          ENDIF
C
C  EROSION FROM PARENT BED
C
          IF (ACTLAY(2,I,J).EQ.0.0) THEN
            DO 55 K=1,2
              FALAY(K,I,J)=FPBED(K,I,J)
              FPBED(K,I,J)=FPBED(K,I,J)*(BEDTH(2,I,J)-ACTLAY(1,I,J)
     +                      -EROTOT)/(BEDTH(1,I,J)-ACTLAY(1,I,J))
 55         CONTINUE
            GOTO 80
          ENDIF
C
C  INCREASING ACTIVE LAYER THICKNESS
C
          IF (ACTLAY(1,I,J).GT.ACTLAY(2,I,J)) THEN
            DO 60 K=1,2
              FTEMP=0.0
              IF (ACTLAY(1,I,J).GT.0.0) THEN
                FTEMP=FPBED(K,I,J)+
     +               ((FALAY(K,I,J)-FPBED(K,I,J))*ACTLAY(2,I,J)+
     +               DEPK(K)-FALAY(K,I,J)*EROTOT)/ACTLAY(1,I,J)
              ENDIF
              FPBED(K,I,J)=FPBED(K,I,J)*(BEDTH(2,I,J)-ACTLAY(1,I,J))/
     +              (BEDTH(1,I,J)-ACTLAY(1,I,J))
              FALAY(K,I,J)=FTEMP
 60         CONTINUE
            GOTO 80
          ENDIF
C
C  CONSTANT OR DECREASING ACTIVE LAYER THICKNESS
C
          IF (ACTLAY(1,I,J).LE.ACTLAY(2,I,J)) THEN
            DO 70 K=1,2
              FTEMP=0.0
              IF (ACTLAY(1,I,J).GT.0.0) THEN
                FTEMP=FALAY(K,I,J)+(DEPK(K)-FALAY(K,I,J)*EROTOT)/
     +                                               ACTLAY(1,I,J)
              ENDIF
              FPBED(K,I,J)=(FPBED(K,I,J)*(BEDTH(2,I,J)-ACTLAY(2,I,J))+
     +            FALAY(K,I,J)*(ACTLAY(2,I,J)-ACTLAY(1,I,J)))/
     +            (BEDTH(1,I,J)-ACTLAY(1,I,J))
              FALAY(K,I,J)=FTEMP
 70         CONTINUE
          ENDIF
C
C  PREPARE FOR NEXT TIMESTEP
C
 80       ACTLAY(2,I,J)=ACTLAY(1,I,J)
          BEDTH(2,I,J)=BEDTH(1,I,J)
C
 30     CONTINUE
 10   CONTINUE
      RETURN
      END
