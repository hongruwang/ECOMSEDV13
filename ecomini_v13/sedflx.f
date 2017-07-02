      SUBROUTINE SEDFLX
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
C  CALCULATES DEPOSITION, ENTRAINMENT, NET FLUX, TOTAL THICKNESS,
C  LAYER ORDER, AND COMPONENT LAYER FRACTION
C
C**********************************************************************
C
      INCLUDE 'comdeck'
      SAVE
C
      DOUBLE PRECISION DEP,ETOT,ETOTL,EB(14),ZERO
      DOUBLE PRECISION TEMP,T1TEMP(IM,JM),DDMAX  !hli
C
      EPSTAU=0.001
C
      ZERO=0.0
C
      IF (TOR.EQ.'BAROTROPIC') GOTO 85
C
C  CALC. SHEAR STRESS AT INTERFACES (k-1/2) (dynes/cm**2)
C
      DO 50 K=2,KBM1
        DO 50 J=2,JMM1
          DO 50 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
              DZLOC=0.5*(DZ(K)+DZ(K-1))*DT(I,J)
C
              TAU(I,J,K)=10000.*(UMOL+KM(I,J,K))*
     +          SQRT((UBAR(I,J,K-1)-UBAR(I,J,K))**2.
     +           +(VBAR(I,J,K-1)-VBAR(I,J,K))**2.)/DZLOC
            ENDIF
 50   CONTINUE
C
C  CALC. SHEAR STRESS AT k POINTS
C
      DO 60 J=2,JMM1
        DO 60 I=2,IMM1
          IF (DT(I,J).GT.0.0) THEN
            DO 70 K=1,KBM1
              TAU(I,J,K)=0.5*(TAU(I,J,K)+TAU(I,J,K+1)) 
 70         CONTINUE
          ENDIF
 60   CONTINUE
C
C  SETTLING SPEED IN m/s
C
      DO 80 K=1,KBM1
        DO 80 J=2,JMM1
          DO 80 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
C 
C   NO DEPOSITION WHEN BED THICKNESS > BATHYMETRIC DEPTH
C   PUT IN ACTUAL RELATIONSHIP FOR SETTLING INSTEAD OF LOOKUP TABLES
C
              IF (CSED1(I,J,K).GT.0.0.AND. SEDTHK(I,J).LT.H(I,J))THEN
                  WSET1(I,J,K)=ADEP*(1000000.0*CSED1(I,J,K)*TAU(I,J,K))
     +             **DEPEXP/10000.0/100.0 !SETTLING VEL IN m/s
                  WSET1(I,J,K)=AMAX1(WSET1(I,J,K),WS1MIN) !IN m/s
              ELSE
	        WSET1(I,J,K)=0.0
              ENDIF
            ENDIF
 80   CONTINUE
C
C  CALC. DEPOSITION
C
C  COMP. 1 :  FINE, FLOCCULATING
C  COMP. 2 :  COARSE, NON-COHESVIE
C
 85   DO 90 J=2,JMM1
        DO 90 I=2,IMM1
          DD(1,I,J)=0.0
          DD(2,I,J)=0.0
          E(1,I,J)=0.0
          E(2,I,J)=0.0
          WCT1BOT(I,J)=0.0
          WCT2BOT(I,J)=0.0
C
C   NO DEPOSITION WHEN BED THICKNESS > BATHYMETRIC DEPTH
C
          IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0
     +       .AND.SEDTHK(I,J).LT.H(I,J)) THEN
C
C  COMP. 1
C
	    DD(1,I,J)=WSET1(I,J,KBM1)*100.0*CSED1(I,J,KBM1)*DTI*XNSBED       
            DD(1,I,J)=DMAX1(DD(1,I,J),ZERO)
C
C  CALC. MAX. DEPOSITION BASED UPON MASS OF SEDIMENT IN ELEMENT AT kb-1
C
C  UNITS ARE g/cm**2
C
            DDMAX=100.*CSED1(I,J,KBM1)*DT(I,J)*DZ(KBM1)
C
            DD(1,I,J)=DMIN1(DD(1,I,J),DDMAX)
C
C  NO DEPOSITION OF COMPONENT 1 IF TAU > TCRDEP
C  USE PROBABILITY OF DEPOSITION
C
            ITAU=NINT(100.*TAU(I,J,KB))
            ITAU=MAX0(ITAU,1)
            ITAU=MIN0(ITAU,20000)
            DD(1,I,J)=PDEP(ITAU)*DD(1,I,J)
C
C  COMP. 2
C
            IF (KSED.GT.1) THEN
              DD(2,I,J)=WS2*CSED2(I,J,KBM1)
C
C  NO DEPOSITION OF COMPONENT 2 IF TAU > TCRDP2
C
              IF (TAU(I,J,KB).GE.TCRDP2) DD(2,I,J)=0.0
C
              DD(2,I,J)=DMAX1(DD(2,I,J),ZERO)
C
C  CALC. MAX. DEPOSITION BASED UPON MASS OF SEDIMENT IN ELEMENT AT kb-1
C
C  UNITS ARE g/cm**2
C
              DDMAX=100.*CSED2(I,J,KBM1)*DT(I,J)*DZ(KBM1)
C
              DD(2,I,J)=DMIN1(DD(2,I,J),DDMAX)
            ENDIF
          ENDIF
 90   CONTINUE
C
C  DETERMINE INITIAL AMOUNT OF COMP.1 IN LAYER 1
C
      DO 100 J=2,JMM1
        DO 100 I=2,IMM1
          IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
            T1TEMP(I,J)=PSED1(1,I,J)*TSED(1,I,J)
          ENDIF
 100  CONTINUE
C 
C  ADD DEPOSITION TO TOP LAYER
C
      IF (KSED.GT.1) THEN
        DO 110 J=2,JMM1
          DO 110 I=2,IMM1
            IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
	      TSED(1,I,J)=TSED(1,I,J)+(DD(1,I,J)+DD(2,I,J))
	      EBTOT(1,I,J)=EBTOT(1,I,J)-(DD(1,I,J)+DD(2,I,J))
              EBTOT(1,I,J)=MAX(EBTOT(1,I,J),ZERO)
            ENDIF
 110    CONTINUE
C
        DO 120 J=2,JMM1
          DO 120 I=2,IMM1
            IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
              IF ((DD(1,I,J)+DD(2,I,J)).GT.ZERO) THEN
                LAYER(1,I,J)=1
	        PSED1(1,I,J)=(T1TEMP(I,J)+DD(1,I,J))/TSED(1,I,J)
	        PSED2(1,I,J)=1.-PSED1(1,I,J)
              ENDIF
            ENDIF
 120    CONTINUE
      ELSE
        DO 130 J=2,JMM1
          DO 130 I=2,IMM1
            IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
              TSED(1,I,J)=TSED(1,I,J)+DD(1,I,J)
	      EBTOT(1,I,J)=EBTOT(1,I,J)-DD(1,I,J)
              EBTOT(1,I,J)=MAX(EBTOT(1,I,J),ZERO)
            ENDIF
 130    CONTINUE
C
        DO 135 J=2,JMM1
          DO 135 I=2,IMM1
            IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
              IF (DD(1,I,J).GT.ZERO) THEN
                LAYER(1,I,J)=1
                PSED1(1,I,J)=1.
              ENDIF
            ENDIF
 135    CONTINUE
      ENDIF
C
C*******************************************************************
C
C  CALC. RESUSPENSION, ERODE LAYERS
C
      DO 140 J=2,JMM1
        DO 140 I=2,IMM1
          E(1,I,J)=0.0
          E(2,I,J)=0.0
C
          IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
            ETOT=0.0
C
C  ASSUME LOWER TAU CRITICAL IN TOP NEW LAYER (L=1)
C
            IF (TAU(I,J,KB).LE.TCRIT(1,I,J)) GO TO 140
C
C  DETERMINE MAX. BASE AMOUNT ENTRAINABLE AT TAU
C  BASED ON MAX. STRESS IN 24 HOUR PERIOD
C
	    ITAU=NINT(TAU(I,J,KB)*10.)
C
            IF (ITAU.LT.1) ITAU=1
C
            IF (ITAU.GT.2000) THEN
              ITAU=2000
	      WRITE (IUPRT,101)I,J,TIME
 101      FORMAT (/5X,'*** BOTTOM SHEAR > 200 dynes/cm**2 AT I,J,TIME ='
     +        ,2I4,F10.4)
	    ENDIF
C
	    IF (LAYER(1,I,J).EQ.1) THEN
	      DTAU = ABS(TAU(I,J,KB)-TAUCUR(1,I,J))
	      IF (DTAU.GT.EPSTAU) THEN
	        TAUCUR(1,I,J) = TAU(I,J,KB)
C
                IF (IA0.EQ.0) THEN
	          EBCUR(1,I,J) = A0(1,I,J)*TAUD2(1,ITAU)
                ELSE
	          EBCUR(1,I,J) = A0(1,I,J)*(((TAU(I,J,KB)-TAUCR(1))/
     +              TAUCR(1))**REXP(I,J))
                ENDIF
              ENDIF
            ELSE
	      GOTO 150
            ENDIF
	    IF (TAU(I,J,KB).GT.TAUMAX(1,I,J).AND.LAYER(1,I,J).EQ.1) 
     +      THEN
	      TAUMAX(1,I,J)=TAU(I,J,KB)
C
              IF (IA0.EQ.0) THEN
	        EBMAX(1,I,J) = A0(1,I,J)*TAUD2(1,ITAU)
              ELSE
    	        EBMAX(1,I,J) = A0(1,I,J)*(((TAU(I,J,KB)-TAUCR(1))/
     +          TAUCR(1))**REXP(I,J))
              ENDIF
	    ENDIF
C
C  DETERMINE BASE AMOUNT ENTRAINED IN THIS TIME STEP
C
            IF (IA0.EQ.0) THEN
	      EB(1) = ASED(1,I,J)*TAUD2(1,ITAU)
            ELSE
	      EB(1) = ASED(1,I,J)*(((TAU(I,J,KB)-TAUCR(1))/
     +          TAUCR(1))**REXP(I,J))
            ENDIF
C
C  CHECK IF 1 DAY OLD LAYER (TOP) HAS BEEN ERODED
C
	    IF (LAYER(1,I,J).EQ.1) THEN
	      IF (EBTOT(1,I,J).GT.EBCUR(1,I,J)) GOTO 140
	      EBTOT(1,I,J)=EBTOT(1,I,J)+EB(1)
	      IF (EBTOT(1,I,J).GE.EBCUR(1,I,J)) THEN
	        IF (EBTOT(1,I,J).LT.EBMAX(1,I,J)) THEN
	          EB(1)=EB(1)-EBTOT(1,I,J)+EBCUR(1,I,J)
	          EBTOT(1,I,J)=EBCUR(1,I,J)
	          ETOT=EB(1)
                ELSE
	          EB(1)=EB(1)-EBTOT(1,I,J)+EBMAX(1,I,J)
	          EBTOT(1,I,J)=EBMAX(1,I,J)
	          ETOT=EB(1)
                ENDIF
              ELSE
	        ETOT=EB(1)
              ENDIF
	      IF (EBTOT(1,I,J).GE.EBMAX(1,I,J)) THEN
	        EB(1)=EB(1)-EBTOT(1,I,J)+EBMAX(1,I,J)
	        EBTOT(1,I,J)=EBMAX(1,I,J)
	        ETOT=EB(1)
              ELSE
	        IF (EBTOT(1,I,J).LT.EBCUR(1,I,J)) ETOT=EB(1)
              ENDIF
	      TEMP=TSED(1,I,J)-EB(1)
	      IF (TEMP.LE.0.0) THEN
	        ETOT=TSED(1,I,J)
	        TSED(1,I,J)=0.0
	        LAYER(1,I,J)=0
	        EBTOT(1,I,J)=0.0
	        TAUMAX(1,I,J)=0.0
	        TAUCUR(1,I,J)=0.0
	        EBMAX(1,I,J)=0.0
	        EBCUR(1,I,J)=0.0
              ELSE
                TSED(1,I,J)=TEMP
              ENDIF
	      E(1,I,J)=PSED1(1,I,J)*ETOT
	      E(2,I,J)=PSED2(1,I,J)*ETOT
            ENDIF
C
C  DIFFERENT TAU CRITICAL FOR DEEPER LAYERS
C
 150        CONTINUE
C
C  DETERMINE MAX. BASE AMOUNT ENTRAINABLE AT TAU
C  BASED ON MAX. STRESS IN 24 HOUR PERIOD
C
	    DO 160 LL=2,LAYMAX
              EB(LL)=0.0
              IF (TAU(I,J,KB).LE.TCRIT(2,I,J)) GOTO 140
              IF (TAU(I,J,KB).LE.TCRIT(LL,I,J)) GOTO 160
C
	      IF(LAYER(LL,I,J).EQ.1) THEN
	        DTAU = ABS(TAU(I,J,KB)-TAUCUR(LL,I,J))
	        IF (DTAU.GT.EPSTAU) THEN
	          TAUCUR(LL,I,J) = TAU(I,J,KB)
C
                  IF (IA0.EQ.0) THEN
	            EBCUR(LL,I,J) = A0(LL,I,J)*TAUD2(LL,ITAU)
                  ELSE
	            EBCUR(LL,I,J) = A0(LL,I,J)*(((TAU(I,J,KB)-TAUCR(LL))
     +                /TAUCR(LL))**REXP(I,J))
                  ENDIF
                ENDIF
C
	        IF (TAU(I,J,KB).GT.TAUMAX(LL,I,J)) THEN
	          TAUMAX(LL,I,J)=TAU(I,J,KB)
C
                  IF (IA0.EQ.0) THEN
	            EBMAX(LL,I,J) = A0(LL,I,J)*TAUD2(LL,ITAU)
                  ELSE
	            EBMAX(LL,I,J) = A0(LL,I,J)*(((TAU(I,J,KB)-
     +                TAUCR(LL))/TAUCR(LL))**REXP(I,J))
                  ENDIF
	        ENDIF
              ENDIF
C
C  DETERMINE BASE AMOUNT ENTRAINED IN THIS TIME STEP
C
              IF (IA0.EQ.0) THEN
	        EB(LL) = ASED(LL,I,J)*TAUD2(LL,ITAU)
              ELSE
	        EB(LL) = ASED(LL,I,J)*(((TAU(I,J,KB)-TAUCR(LL))/
     +            TAUCR(LL))**REXP(I,J))
              ENDIF
 160        CONTINUE
C
C  CHECK OTHER LAYERS FOR EROSION, ERODE ONLY 1 LAYER
C  DURING THIS TIMESTEP
C
            IF (LAYER(1,I,J).EQ.0) THEN
	      DO 170 LL=2,LAYMAX
	        IF (LAYER(LL,I,J).EQ.1) THEN
	          ETOTL=0.0
	          IF (EBTOT(LL,I,J).GT.EBCUR(LL,I,J)) GOTO 140
	          EBTOT(LL,I,J)=EBTOT(LL,I,J)+EB(LL)
	          IF (EBTOT(LL,I,J).GE.EBCUR(LL,I,J)) THEN
		    IF (EBTOT(LL,I,J).LT.EBMAX(LL,I,J)) THEN
		      EB(LL)=EB(LL)-EBTOT(LL,I,J)+EBCUR(LL,I,J)
		      EBTOT(LL,I,J)=EBCUR(LL,I,J)
		      ETOTL=EB(LL)
                    ELSE
		      EB(LL)=EB(LL)-EBTOT(LL,I,J)+EBMAX(LL,I,J)
		      EBTOT(LL,I,J)=EBMAX(LL,I,J)
		      ETOTL=EB(LL)
                    ENDIF
	          ELSE
		    ETOTL=EB(LL)
                  ENDIF
	          TEMP=TSED(LL,I,J)-EB(LL)
	          IF (TEMP.LE.0.0) THEN
		    ETOTL=TSED(LL,I,J)
		    TSED(LL,I,J)=0.0
		    LAYER(LL,I,J)=0
		    EBTOT(LL,I,J)=0.0
		    TAUMAX(LL,I,J)=0.0
		    TAUCUR(LL,I,J)=0.0
		    EBMAX(LL,I,J)=0.0
		    EBCUR(LL,I,J)=0.0
                  ELSE
                    TSED(LL,I,J)=TEMP
                  ENDIF
                  ETOT=ETOT+ETOTL
                  E(1,I,J)=E(1,I,J)+PSED1(LL,I,J)*ETOTL
                  E(2,I,J)=E(2,I,J)+PSED2(LL,I,J)*ETOTL
	          IF (LAYER(LL,I,J).EQ.0) THEN
		    PSED1(LL,I,J)=0.0
		    PSED2(LL,I,J)=0.0
                  ENDIF
	          GOTO 140
                ENDIF
 170          CONTINUE
	    ENDIF
          ENDIF
 140  CONTINUE
C
C*****************************************************************
C
C  DETERMINE TOTAL SEDIMENT FLUX
C  TOTAL AMOUNT DEP./RES. (gm/cm**2)
C
      DO 180 J=2,JMM1
        DO 180 I=2,IMM1
          IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
            WCT1BOT(I,J)=(E(1,I,J)-DD(1,I,J))/XNSBED
          ENDIF
 180  CONTINUE
C
      IF (KSED.GT.1) THEN
        DO 190 J=2,JMM1
          DO 190 I=2,IMM1
            IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
              WCT2BOT(I,J)=(E(2,I,J)-DD(2,I,J))/XNSBED
            ENDIF
 190    CONTINUE
      ENDIF
C
      IF (KSED.GT.1) THEN
        DO 200 J=2,JMM1
          DO 200 I=2,IMM1
            IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
              IF (LAYER(1,I,J).EQ.0) THEN
	        PSED1(1,I,J)=0.0
	        PSED2(1,I,J)=0.0
              ELSE
	        PSED1(1,I,J)=(T1TEMP(I,J)-(E(1,I,J)-DD(1,I,J)))/
     +                                   TSED(1,I,J)
	        PSED2(1,I,J)=1.-PSED1(1,I,J)
              ENDIF
            ENDIF
 200    CONTINUE
      ELSE
        DO 210 J=2,JMM1
          DO 210 I=2,IMM1
            IF (FSM(I,J).GT.0.0.AND.IBMSK(I,J).EQ.0) THEN
              IF (LAYER(1,I,J).EQ.0) THEN
	        PSED1(1,I,J)=0.
              ELSE
	        PSED1(1,I,J)=1.
              ENDIF
            ENDIF
 210    CONTINUE
      ENDIF
C
C  CONVERT E AND D BACK TO gm/cm**2/dti (FOR WASTOX)
C  MAKES COHESIVE SEGMENTS CONSISTENT WITH NON-COHESIVE SEGMENTS
C
      DO 220 J=1,JM
        DO 220 I=1,IM
          E(1,I,J)=E(1,I,J)/XNSBED
          DD(1,I,J)=DD(1,I,J)/XNSBED
          E(2,I,J)=E(2,I,J)/XNSBED
          DD(2,I,J)=DD(2,I,J)/XNSBED
 220  CONTINUE
      RETURN
      END
