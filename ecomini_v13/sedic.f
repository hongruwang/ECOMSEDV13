      SUBROUTINE SEDIC(RESTAR)
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
C  SETS INITIAL CONDITIONS AND PARAMETERS FOR SEDIMENT TRANSPORT
C
C***********************************************************************
C
      INCLUDE 'comdeck'
C
      CHARACTER*10 RESTAR
      REAL PINIT(IM,JM)
      REAL D50VAR(IM,JM)
C
      WS1MIN=6.6E-06  ! MINIMUM SETTLING VELOCITY FOR COHESIVE SEDIMENT IN m/s
C
      PI=ACOS(-1.)
C
C  VARIABLE INITIALIZATION
C
C  BED MASK
C
      IF (RESTAR.EQ.'COLD START') THEN
        DO 4 J=1,JM
          DO 4 I=1,IM
            IBMSK(I,J)=0
            DO 3 K=1,KSMAX
              FALAY(K,I,J)=0.0
 3          CONTINUE
 4      CONTINUE
      ENDIF
C
C  INPUT BED LOAD MASK
C  MASK FOR FINE SEDIMENT AREAS
C
C  IBED = 0 :  NO NON-COHESIVE SUSPENDED LOAD
C  IBED = 1 :  NON-COHESIVE SUSPENDED LOAD
C
      IBED=0
C
        OPEN (UNIT=85,FILE='bed_mask',FORM='FORMATTED')
        DO 5 I=1,IM
          READ (85,6) (IBMSK(I,J),J=1,JM)
C
C  CHECK FOR NON-COHESIVE SUSPENDED LOAD
C
          DO 95 J=1,JM
            IF (IBMSK(I,J).EQ.1) IBED=1
 95       CONTINUE
C
C  SEARCH FOR INCORRECT BED MASK FOR SAND TRANSPORT ONLY
C
          IF (SEDTYPE.EQ.'SAND') THEN
            DO 96 J=1,JM
              IF (IBMSK(I,J).EQ.0) THEN
                WRITE (IUPRT,97)I,J
 97             FORMAT (//5X,'**** PROGRAM EXECUTION STOPPED ****',/
     +     /5X,'INCORRECT BED MASK FOR SAND TRANSPORT ONLY',
     +     /5X,'BED MASK SET TO 0 AT i,j =',2I5,
     +     /5X,'MUST SPECIFY BED MASK VALUE AS 1 (SUSPENDED TRANSPORT)',
     +     /5X,'OR BED MASK VALUE AS -1 (HARD BOTTOM)',
     +     /5X,'PLEASE CORRECT bed_mask FILE AND RESUBMIT'/)
                STOP
              ENDIF
 96         CONTINUE
          ENDIF
 5      CONTINUE
 6      FORMAT (40I3)
        CLOSE (85)
C
C  NO SUSPENDED LOAD TRANSPORT IN COHESIVE ONLY RUN
C
        IF (SEDTYPE.EQ.'MUD ') IBED=0
C
C  VARIABLE BED BULK DENSITY (g/cm**3)
C
C  FOR COMPUTING BED ELEVATION CHANGES AND FOR WASTOX/WASP USE
C
      IF (VARIBULK.EQ.'INCLUDE') THEN
        OPEN (UNIT=85,FILE='bed_bulkden',FORM='FORMATTED')
        DO 9019 I=1,IM
          READ (85,9018) (CBED(I,J),J=1,JM)
 9019   CONTINUE
 9018   FORMAT (20F6.0)
        CLOSE(85)
      ELSE
        DO 9020 I=1,IM
          DO 9020 J=1,JM
C
C  COHESIVE ELEMENTS
C
            IF (IBMSK(I,J).EQ.0) THEN
              CBED(I,J)=DENCOH
            ELSE
              CBED(I,J)=DENNON
            ENDIF
 9020   CONTINUE
      ENDIF
C
C  INPUT SPATIALLY VARIABLE P0 
C
      IF (IP0.EQ.1) THEN
        OPEN (UNIT=86,FILE='p0_init',FORM='FORMATTED')
        DO 7 I=1,IM
          READ (86,8) (PINIT(I,J),J=1,JM)
 7      CONTINUE
 8      FORMAT (20F6.2)
        CLOSE (86)
      ENDIF
C
C  INPUT SPATIALLY VARIABLE A0 AND REXP
C
C  A0 INPUT FOR mg/cm**2 EROSION POTENTIAL
C
      IF (IA0.EQ.1) THEN
        OPEN (UNIT=53,FILE='a0_init',FORM='FORMATTED')
        DO 9000 I=1,IM
          READ (53,9001) (A0(1,I,J),J=1,JM)
C
C  CONVERT FROM mg/cm**2 TO g/cm**2
C
          DO 9000 J=1,JM
            A0(1,I,J)=A0(1,I,J)/1000.
 9000   CONTINUE
 9001   FORMAT (8E10.3)
        DO 9002 LL=2,LAYMAX
          DO 9002 J=1,JM
            DO 9002 I=1,IM
              A0(LL,I,J)=A0(1,I,J)
 9002   CONTINUE
        CLOSE (53)
C
        OPEN (UNIT=53,FILE='exp_init',FORM='FORMATTED')
        DO 9010 I=1,IM
          READ (53,9001) (REXP(I,J),J=1,JM)
 9010   CONTINUE
        CLOSE (53)
      ENDIF
C
      DO 10 J=1,JM
        DO 10 I=1,IM
          IF (RESTAR.EQ.'COLD START') THEN  
          DO 15 K=1,KB
            TAU(I,J,K)=0.0
 15       CONTINUE
          ENDIF                            
	  DO 20 LL=1,LAYMAX
C
C  SET CRITICAL SHEAR STRESS AS SPATIALLY CONSTANT OVER EACH LAYER
C
            TCRIT(LL,I,J)=TAUCR(LL)
C
            IF (IBMSK(I,J).EQ.0) THEN
	      TSED0(LL,I,J)=TSED0IN(LL)*CBED(I,J)
            ELSE
              TSED0(LL,I,J)=0.0
            ENDIF
C
C  BED LOAD THICKNESS INPUT IN m, CONVERT TO VOLUME (m**3)
C  NOTE: V(SED) = (1. - BPOR) * V(TOTAL)
C        WHERE V(TOTAL) = AREA * BEDTHI
C
           IF (RESTAR.EQ.'COLD START') THEN
            N24CNT=0
            IF (IBMSK(I,J).EQ.0) THEN
	      TSED(LL,I,J)=TSED0(LL,I,J)
              BEDTH(1,I,J)=0.0
              BEDTH(2,I,J)=0.0
            ELSE
              IF (IBMSK(I,J).EQ.1) THEN
                TSED(LL,I,J)=0.0
                BEDTH(1,I,J)=(CBED(I,J)/2.65)*BEDTHI*H1(I,J)*H2(I,J)
                BEDTH(2,I,J)=(CBED(I,J)/2.65)*BEDTHI*H1(I,J)*H2(I,J)
              ELSE
                TSED(LL,I,J)=0.0
                BEDTH(1,I,J)=0.0
                BEDTH(2,I,J)=0.0
              ENDIF
            ENDIF
C
C
	    TAUMAX(LL,I,J)=0.0
	    TAUCUR(LL,I,J)=0.0
	    EBMAX(LL,I,J)=0.0
	    EBCUR(LL,I,J)=0.0
C
	    IF (TSED(LL,I,J).GT.0.0) THEN
	      LAYER(LL,I,J)=1
            ELSE
	      LAYER(LL,I,J)=0
            ENDIF
	    IF (KSED.EQ.1) THEN
	      IF (LAYER(LL,I,J).EQ.1) THEN
	        PSED1(LL,I,J)=1.0
	        PSED2(LL,I,J)=0.0
              ELSE
	        PSED1(LL,I,J)=0.0
	        PSED2(LL,I,J)=0.0
	      ENDIF
            ELSE
	      IF (LAYER(LL,I,J).EQ.1) THEN
                IF (IP0.EQ.0) THEN
	          PSED1(LL,I,J)=P0(1)
	          PSED2(LL,I,J)=1.-P0(1)
                ELSE
                  PSED1(LL,I,J)=PINIT(I,J)
                  PSED2(LL,I,J)=1.-PINIT(I,J)
                ENDIF
              ELSE
	        PSED1(LL,I,J)=0.0
                PSED2(LL,I,J)=0.0
              ENDIF
	    ENDIF
C
            IF (IBMSK(I,J).EQ.1.OR.IBMSK(I,J).LT.0) THEN
              PSED1(LL,I,J)=0.0
              PSED2(LL,I,J)=0.0
            ENDIF
           ENDIF
 20       CONTINUE
C
          WCT1BOT(I,J)=0.0
          WCT2BOT(I,J)=0.0
 10   CONTINUE
C
      DO 26 K=1,KB
        DO 26 J=1,JM
          DO 26 I=1,IM
            WSET1(I,J,K)=0.0
            WSET2(I,J,K)=0.0
 26   CONTINUE
C
C      TCRDP2=1.
       TCRDP2=(WS2/10000.)**2.
C
      WRITE (IUPRT,27)TCRDP2
 27   FORMAT (/5X,'CLASS 2 CRITICAL SHEAR STRESS FOR DEP (dyne/cm**2)',
     +            F10.4)
C
C  WATER-COLUMN SETTLING SPEED FOR COARSE SEDIMENT 
C  AND NON-COHESIVE SUSPENDED SEDIMENTS
C  CONVERT FROM microns/s TO m/s
C
      WMIN=1000000.
      WMAX=-1000000.
C
      DO 30 K=1,KBM1
        DO 30 J=1,JM
          DO 30 I=1,IM
            WSET2(I,J,K)=WS2/1000000.
C
C  SET MAX. SETTLING SPEED FOR CLASS 2 EQUAL TO 0.25*CFL(VERT)
C
            IF (FSM(I,J).GT.0.0.AND.K.GE.2) THEN
              CFLZ1=WSET2(I,J,K)*DTI/(H(I,J)*DZ(K-1))
              CFLZ2=WSET2(I,J,K)*DTI/(H(I,J)*DZ(K))
C
              WST1=WSET2(I,J,K)
              WST2=WSET2(I,J,K)
C
              IF (CFLZ1.GT.0.1) THEN
                WST1=0.25*H(I,J)*DZ(K-1)/DTI
              ENDIF
C
              IF (CFLZ2.GT.0.1) THEN
                WST2=0.25*H(I,J)*DZ(K)/DTI
              ENDIF
C
              WSET2(I,J,K)=AMIN1(WST1,WST2)
C
              IF (WSET2(I,J,K).LT.WMIN) THEN
                IMIN=I
                JMIN=J
                KMIN=K
                WMIN=WSET2(I,J,K)
              ENDIF
C
              IF (WSET2(I,J,K).GT.WMAX) THEN
                IMAX=I
                JMAX=J
                KMAX=K
                WMAX=WSET2(I,J,K)
              ENDIF
            ENDIF
 30   CONTINUE
C
C  SET WSET2 EQUAL TO MINIMUM SETTLING SPEED
C
      DO 301 K=1,KBM1
        DO 301 J=1,JM
          DO 301 I=1,IM
            WSET2(I,J,K)=WMIN
 301  CONTINUE
        WMIN=1000000.*WMIN
        WMAX=1000000.*WMAX
        WRITE (IUPRT,31)IMIN,JMIN,KMIN,WMIN,H(IMIN,JMIN)
 31     FORMAT (/5X,'**** MIN WSET2 & H: ',3I5,F10.1,F10.2)
        WRITE (IUPRT,32)IMAX,JMAX,KMAX,WMAX,H(IMAX,JMAX)
 32     FORMAT (/5X,'**** MAX WSET2 & H: ',3I5,F10.1,F10.2)
C
C  WS2 FOR COMP. 2 (COARSE), INPUT IN microns/s
C  FOR SEDIMENT-WATER INTERFACE FLUX
C
C  FOR SPEED IMPROVEMENT
C  MULT. WS2 BY XNSBED
C
      WS2=XNSBED*WS2*DTI/10000.
C
C  CALC. SETTLING SPEEDS FOR FLOCCULATION
C
C  SEDZL SETTLING SPEEDS:  BASED ON REANALYSIS OF LICK DATA
C 
      DO 40 II=1,6250000
	CT=FLOAT(II)/100.
C
C  SETTLING SPEED IN cm/s
C
	WS1(II)=ADEP*(CT**DEPEXP)/10000.
 40   CONTINUE
C
C  CALC. PROBABILITY OF DEPOSITION
C
C  INCREMENT EVERY 0.01 dyne/cm**2
C 
      DO 45 II=1,20000 
        TAUBOT=FLOAT(II-1)*0.01
C
C  KRONE FORMULATION
C
        IF (PDEPFORM.EQ.'KRONE') THEN
          IF (TAUBOT.LT.TCRDEP) THEN
            PDEP(II)=1.-TAUBOT/TCRDEP
          ELSE
            PDEP(II)=0.0
          ENDIF
        ELSE
C
C  PARTHENIADES FORMULATION
C
C  ASSUME:  SIGY = 0.49
C           (TAUB - 1)50 = 0.84
C
          IF (TAUBOT.LE.TCRDEP) THEN
            PDEP(II)=1.
          ELSE
C
C  CALC. Y PARAMETER
C
          YPARM=2.04*ALOG10(0.25*((TAUBOT/TCRDEP)-1.)*EXP(1.27*TCRDEP))
C
            IF (YPARM.LT.0.0) THEN
              YPARM=-YPARM
              IFLAG=1
            ELSE
              IFLAG=0
            ENDIF
C
C  CALC. Z(Y)
C
            ZY=EXP(-YPARM*YPARM/2.)/SQRT(2.*PI)
C
C  APPROXIMATE EVALUATION OF CEQ INTEGRAL
C
            TY=1./(1.+0.3327*YPARM)
C
            CEQ=1.-ZY*(0.4362*TY-0.1202*(TY**2.)+0.9373*(TY**3.))
C
            IF (IFLAG.EQ.1) CEQ=1.-CEQ
C
C  CALC. PROBABILITY OF DEPOSITION
C
            PDEP(II)=1.-CEQ
          ENDIF
        ENDIF
 45   CONTINUE
C
C  CALC. TAUD2
C
      TAU1=0.0
      DO 60 I=1,2000
	TAU1=TAU1+.1
        DO 70 LL=1,LAYMAX 
	  RAT1=(TAU1-TAUCR(LL))/TAUCR(LL)
	  IF (RAT1.LT.0.0) RAT1=0.0
	  TAUD2(LL,I)=RAT1**RESEXP
 70     CONTINUE
 60   CONTINUE
C
C********************************************************************
C
C  COHESIVE SEDIMENT CONSTANTS
C
C  FOR SPEED IMPROVEMENT
C  MULT. S1 BY XNSBED
C
      S1=XNSBED*DTI/3600.
C
C  CONVERT A0IN FROM mg/cm**2 TO g/cm**2
C
      A0IN=A0IN/1000.
C
      DO 90 LL=1,LAYMAX
	FTIME2(LL)=FTIME(LL)**EXPM
C
        DO 90 J=1,JM
	  DO 90 I=1,IM
            IF (IBMSK(I,J).EQ.0) THEN
              IF (IA0.EQ.0) THEN
	        A0(LL,I,J)=A0IN/FTIME2(LL)
              ELSE
                A0(LL,I,J)=A0(LL,I,J)/FTIME2(LL)
              ENDIF
	      ASED(LL,I,J)=S1*A0(LL,I,J)
            ENDIF
 90   CONTINUE
C
C
C  SUSPENDED LOAD CONSTANTS
C
C  REF:  VAN RIJN, 1984, J. HYDR. ENGR., p. 1431-1456
C        VAN RIJN, 1984, J. HYDR. ENGR., p. 1613-1641
C
C  MKS SYSTEM USED
C  ASSUMED THAT SEDIMENT DENSITY IS 2.65 g/cm**3 AND KINEMATIC
C  VISCOSITY IS 0.011 cm**2/s
C
      IF (IBED.GE.1) THEN
C
C************************************************
C
C  INPUT VARIABLE D50
C
        OPEN (UNIT=85,FILE='bed_d50',FORM='FORMATTED')
        REWIND(85)
        DO 104 I=1,IM
          READ (85,103) (D50VAR(I,J),J=1,JM)
 104    CONTINUE
 103    FORMAT (20F6.0)
        CLOSE (85)
C
C  CONVERT FROM MICRONS TO METERS
C
        DO 102 J=1,JM
          DO 102 I=1,IM
            D50VAR(I,J)=D50VAR(I,J)/1000000.
 102    CONTINUE
C
C  COMPUTE EFFECTIVE D50 OF CLASS 2 FROM SETTLING SPEED (IN cm)
C
        D50C2=0.05
        ACON1=133636.
        ACON2=9.091*WSET2(1,1,1)*100.
C
C  USE NEWTON'S METHOD TO SOLVE FOR D50, 100 ITERATIONS
C
        DO 105 N=1,20 
          FD50=SQRT(1.+ACON1*(D50C2**3.))-1.-ACON2*D50C2
          FPD50=((1.5*ACON1*D50C2*D50C2)/SQRT(1.+ACON1*(D50C2**3.)))
     +             -ACON2
          D50C2=D50C2-FD50/FPD50
 105    CONTINUE
        D50C2=D50C2*10000.
        WRITE (IUPRT,106)D50C2
 106    FORMAT (/5X,'NON-COHESIVE EFFECTIVE D50 (microns) =',F8.0
     +          /5X,'(BASED UPON SETTLING SPEED)'/)
        D50C2=D50C2/1000000.
C
C  VARIABLE D50
C
        DO 117 J=1,JM
          DO 117 I=1,IM
            IF (IBMSK(I,J).EQ.1.AND.FSM(I,J).GT.0.0) THEN
C
C
C  CALC. DSTAR (EQ. 1)
C
              DSTAR=23731.*D50VAR(I,J)
C
C  CALC. THETA CRITICAL (FIG. 1)
C
              IF (DSTAR.LE.4.) THEN
                THCR=0.24/DSTAR
                GOTO 110
              ENDIF
C
              IF (DSTAR.LE.10.) THEN
                THCR=0.14/(DSTAR**0.64)
                GOTO 110
              ENDIF
C
              IF (DSTAR.LE.20) THEN
                THCR=0.04/(DSTAR**0.10)
                GOTO 110
              ENDIF
C
              IF (DSTAR.LE.150) THEN 
                THCR=0.013*(DSTAR**0.29)
              ELSE
                THCR=0.055
              ENDIF
C
C  CALC. USTAR CRITICAL (FIG. 1)
C
 110          UCRVAR(I,J)=SQRT(16.17*D50VAR(I,J)*THCR)
            ENDIF
 117    CONTINUE
C
        CSED(8)=DSTAR
C
C  CLASS 2
C
C
C  CALC. DSTAR (EQ. 1)
C
        DSTAR2=23731.*D50C2
C
C  CALC. THETA CRITICAL (FIG. 1)
C
        IF (DSTAR2.LE.4.) THEN
          THCR=0.24/DSTAR2
          GOTO 111
        ENDIF
C
        IF (DSTAR2.LE.10.) THEN
          THCR=0.14/(DSTAR2**0.64)
          GOTO 111
        ENDIF
C
        IF (DSTAR2.LE.20) THEN
          THCR=0.04/(DSTAR2**0.10)
          GOTO 111
        ENDIF
C
        IF (DSTAR2.LE.150) THEN 
          THCR=0.013*(DSTAR2**0.29)
        ELSE
          THCR=0.055
        ENDIF
C
C  CALC. USTAR CRITICAL (FIG. 1)
C
 111    UCR2=SQRT(16.17*D50C2*THCR)
C
C  SUSPENDED LOAD CONSTANTS
C
C  NOTE:  CSED(1) AND CSED(2) ARE USED IN BEDLOAD CALC.
C
        CSED(1)=SQRT(9.8)
C
        CSED(3)=0.015*D50C2/(DSTAR2**0.3)
C
C  INPUT NON-COHESIVE MASS FRACTIONS
C
       IF (RESTAR.EQ.'COLD START') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (UNIT=85,FILE='bed_frac.mud',FORM='FORMATTED')
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (UNIT=86,FILE='bed_frac.sand',FORM='FORMATTED')
        ENDIF
C
        DO 120 K=1,2
          DO 120 I=1,IM
            DO 125 J=1,JM
              FPBED(K,I,J)=0.0
 125        CONTINUE
C
            IF((SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH').AND.K.EQ.1) THEN
              READ (85,130) (FPBED(K,I,J),J=1,JM)
            ENDIF
C
            IF((SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH').AND.K.EQ.2) THEN
              READ (86,130) (FPBED(K,I,J),J=1,JM)
            ENDIF
C
              DO 140 J=1,JM
                IF (IBMSK(I,J).EQ.1.AND.FSM(I,J).GT.0.0) THEN
                  FALAY(K,I,J)=0.0
C
                  ACTLAY(2,I,J)=0.0
C
                  CARMOR(I,J)=(2.*D50VAR(I,J)*(CBED(I,J)/2.65)*
     +              H1(I,J)*H2(I,J))/(10000.*UCRVAR(I,J)*UCRVAR(I,J))
                ENDIF
 140          CONTINUE
 120    CONTINUE
 130    FORMAT (20F6.3)
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') CLOSE (85)
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') CLOSE (86)
       ENDIF
      ENDIF
      RETURN
      END
