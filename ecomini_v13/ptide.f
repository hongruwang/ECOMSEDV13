C TIDAL PREDICTION PROGRAM   FOR G GODIN   OCEANOGRAPHY   DEC/72
C THE MOST RECENT REVISION WAS DONE IN APRIL 1977 BY M FOREMAN.
C
C INPUT
C =====
C  1) THE STANDARD CONSTITUENT DATA PACKAGE IS READ BY ROUTINE ASTRO.
C     SEE THIS SUBROUTINE FOR A DESCRIPTION OF THE INPUT REQUIRED. IN
C     ORDER TO INSURE A CONSISTENT INTERPRETATION OF THE AMPLITUDES AND
C     PHASES, THE DATA PACKAGE SHOULD BE IDENTICAL TO ONE USED BY THE
C     TIDAL HEIGHTS ANALYSIS PROGRAM.
C
C  2) ONE RECORD WITH THE TIDAL STATION INFORMATION ISTN,(NA(J),J=1,5),
C     ITZONE,LAD,LAM,LOD,LOM IN THE FORMAT
C     (6X,I4,1X,5A4,2X,A4,I2,1X,I2,2X,I3,1X,I2), WHERE
C          ISTN    = STATION NUMBER
C          NA      = STATION NAME
C          ITZONE  = TIME ZONE
C          LAD,LAM = LATITUDE IN DEGREES AND MINUTES OF THE STATION
C          LOD,LOM = LONGITUDE IN DEGREES AND MINUTES OF THE STATION.
C     FOLLOWED BY ONE RECORD FOR EACH CONSTITUENT TO BE INCLUDED IN THE
C     PREDICTION, WITH THE CONSTITUENT NAME(KON), AMPLITUDE(AMP), AND
C     PHASE LAG(G) IN THE FORMAT (5X,A4,13X,F6.1,12X,F5.1). THE PHASE LAG
C     UNITS SHOULD BE DEGREES WHILE THE UNITS OF THE PREDICTED TIDAL
C     HEIGHTS WILL BE THE SAME AS THOSE OF THE INPUT AMPLITUDES.
C     THE PHASE LAG (THE GREENWICH PHASE LAG) IS THE LOCAL PHASE.
C     THE LAST CONSTITUENT RECORD SHOULD BE FOLLOWED BY A BLANK RECORD.
C
C  3) ONE RECORD CONTAINING THE FOLLOWING INFORMATION ON THE PERIOD AND
C     TYPE OF PREDICTION DESIRED, THE FORMAT IS (3I3,1X,3I3,1X,A4,F17.1)
C        IDY0,IMO0,IYR0 = DAY, MONTH, AND YEAR OF THE BEGINNING OF THE
C                         PREDICTION.
C        IDYE,IMOE,IYRE = DAY, MONTH, AND YEAR OF THE END OF THE
C                         PREDICTION.
C        ITYPE          = "EQUI" IF EQUALLY SPACED PREDICTIONS ARE
C                         DESIRED,
C                         "EXTR" IF DAILY HIGH LOWS ARE DESIRED.
C        DT             = TIME SPACING OF THE PREDICTED VALUES IF
C                         ITYPE="EQUI",
C                         TIME STEP INCREMENT USED TO INITIALLY BRACKET
C                         A HIGH OR LOW IF ITYPE="EXTR".
C     RECOMMENDED VALUES FOR DT WHEN ITYPE="EXTR" ARE GIVEN IN THE USER
C     MANUAL. EQUALLY SPACED PREDICTIONS BEGIN AT DT HOURS ON THE FIRST
C     DAY AND END AT 2400 HOURS OF THE LAST DAY.
C
C     TYPE 3 DATA MAY BE REPEATED ANY NUMBER OF TIMES . A BLANK RECORD
C     CAUSES THE PROGRAM TO GO BACK TO READ MORE TYPE 2 DATA . TWO
C     BLANKS TERMINATE THE JOB.
C
C OUTPUT
C ======
C     EQUI-SPACED PREDICTIONS ARE PRINTED 8 VALUES PER LINE . HI-LO
C     PREDICTIONS ARE PRINTED IN THE FORMAT OF THE CANADIAN TIDES AND
C     WATER LEVELS DAILY HI-LO CARD.
C     AS WELL THE PREDICTIONS ARE PUT ONTO LOGICAL UNIT 10 IN SIMILAR
C     CARD IMAGE FORMAT , WITH A BLANK CARD AFTER EACH SET OF PREDS.
C     A LIST OF THE CONSTITUENTS USED IS PRINTED AFTER EACH SET OF
C     PREDICTIONS.
C
C METHOD OF COMPUTING HIGH-LOW PREDICTIONS
C ========================================
C     THE FIRST DERIVATIVE OF THE TIDAL HGT IS CALCULATED AS THE SUM OF
C     A NUMBER OF TRIGONOMETRIC TERMS (THE CONSTITUENTS) . THE DERIV-
C     ATIVE IS FOUND AT SUCCESIVE TIMES IN STEPS OF DT UNTIL A CHANGE IN
C     SIGN IS NOTED. THE LAST TIME INTERVAL IS THEN BISECTED REPEATEDLY
C     UNTIL THE ROOT IS BRACKETED TO WITHIN 6 MINUTES AND THEN LINEAR
C     INTERPOLATION IS USED TO GET THE ROOT .  THE CHEBYSHEV ITERATION
C     FORMULA ... T(N+1) = 2*DCOS(SIGMA*DT)*T(N)-T(N-1) ...
C     IS USED TO UPDATE EACH CONSTITUENT IN BOTH THE STEPPING AND
C     BISECTION PROCESSES .  FINALLY THE TIDAL HGT AT EXTREMUM IS
C     FOUND USING TABLE LOOK-UP FOR THE COSINES NEEDED.
C     FOR EQUI-SPACED PREDICTIONS THE CHEBYSHEV FORMULA IS APPLIED
C     DIRECTLY.  THE NODAL CORRECTIONS AND THE CHEBYSHEV  ITERATION ARE
C     RENEWED AT THE START OF EACH MONTH .
C
C ROUTINES NEEDED
C ===============
C     ASTRO ... GETS FREQUENCIES,ASTRO ARGS AND NODAL CORRECTIONS
C     CDAY ... A DAY CALENDAR PROGRAM
C
C TIME OF EXECUTION ON CDC 6400 COMPUTER
C ======================================
C     68 SEC FOR 1 YR OF HI-LOS FOR QUEBEC (STN 3250) 174 CONS,DT=3 HRS
C     54 SEC FOR 1 YR OF HOURLY HEIGHTS FOR QUEBEC
C     44 SEC FOR 1YR OF HI-LOS FOR VICTORIA (STN 7120) 62 CONS,DT=.5HRS
C
C
C
C DESCRIPTION OF IMPORTANT VARIABLES OR ARRAYS
C ============================================
C MTAB                  IS THE NUMBER OF CONSTITUENTS INCLUDED IN THE
C                       STANDARD CONSTITUENT DATA PACKAGE. AT PRESENT
C                       THIS IS 146.
C M                     IS THE NUMBER OF CONSTITUENTS TO BE INCLUDED IN
C                       THE PREDICTION.
C KONTAB(I),SIGTAB(I)   ARE ARRAYS CONTAINING THE CONSTITUENT NAMES AND
C                       FREQUENCIES AS THEY ARE READ IN AND CALCULATED
C                       FROM THE STANDARD CONSTITUENT DATA PACKAGE
C V(I),U(I),F(I)        ARE ARRAYS CONTAINING THE ASTRONOMICAL ARGUMENT,
C                       SATELLITE(NODAL MODULATION) PHASE CORRECTION AND
C                       SATELLITE AMPLITUDE CORRECTION RESPECTIVELY, FOR
C                       CONSTITUENT KONTAB(I).
C KON(I),SIG(I)         ARE ARRAYS CONTAINING THE NAMES AND FREQUENCIES
C                       OF CONSTITUENTS TO BE INCLUDED IN THE
C                       PREDICTION.
C AMP(I),G(I)           ARE ARRAYS CONTAINING THE AMPLITUDE AND PHASE
C                       LAG RESPECTIVELY FOR CONSTITUENT KON(I).
C INDX(I)               IS AN ARRAY WHOSE VALUE IS THE LOCATION OF
C                       CONSTITUENT KON(I) IN THE LIST KONTAB(K).
C TWOC(I),BTWOC(I,K)    ARE ARRAYS FOR STORING THE VALUES
C                       2*DCOS(TWOPI*SIG(I)*DT) AND
C                       2*DCOS(TWOPI*SIG(I)*DT/2**K) RESPECTIVELY.
C CHP(I),CH(I),CHA(I),  ARE ARRAYS USED FOR STORING TIDAL HEIGHTS IN
C CHB(I),CHM(I)         THE CHEBYSHEV ITERATION OR THE BISECTION
C                       SEARCH METHOD.
C ANGO(I),AMPNC(I)      ARE ARRAYS USED FOR STORING THE NODALLY
C                       CORRECTED CONSTITUENT ARGUMENTS AND AMPLITUDES
C                       AT THE BEGINNING OF A CHEBYSHEV ITERATION.
C COSINE(I)             IS AN ARRAY CONTAINING 2002 COSINE VALUES IN
C                       THE RANGE OF ZERO TO TWOPI.
C KONTAB,SIGTAB,V,U, AND F SHOULD HAVE MINIMUM DIMENSION MTAB, WHILE
C KON,SIG,AMP,G,INDX,TWOC,CH,CHP,CHA,CHB,CHM,ANGO, AND AMPNC SHOULD
C HAVE MINIMUM DIMENSION M.  BTWOC SHOULD HAVE A MINIMUM DIMENSION OF M
C BY THE MAXIMUM NUMBER OF ITERATIONS REQUIRED TO LOCATE A HIGH OR LOW
C (PRESENTLY THIS IS SET TO 10).
C
C
      SUBROUTINE PTIDE(NEBCM,NUMEBC,NUMC,ALAT,DAMP,DG,Z0,
     . IDY0,IMO0,IYR0,IHR,NDAY,TZERO)
      PARAMETER(NO1=45,NO2=146,NO3=170,NO4=251)
      PARAMETER(NREC=10000)
      IMPLICIT REAL*8 (A-H,O-Z)
c     IMPLICIT INTEGER*4 (I-N)
      INTEGER*2 KR,LU,LP
      REAL*8 II(NO1),JJ(NO1),KK(NO1),LL(NO1),MM(NO1),NN(NO1),SEMI(NO1)
      REAL*8 EE(NO3),LDEL(NO3),MDEL(NO3),NDEL(NO3),PH(NO3)
      
      DIMENSION KONTAB(NO2),SIGTAB(NO2),V(NO2),U(NO2),F(NO2),INDX(NO2)
     .      ,KON(NO2),AMP(NO2),G(NO2),SIG(NO2)
     .      ,TWOC(NO2),CHP(NO2),CH(NO2)
     .      ,BTWOC(NO2,10),CHA(NO2),CHB(NO2),CHM(NO2)
     .      ,ANG0(NO2),AMPNC(NO2),COSINE(2002)
      DIMENSION NA(5)
      REAL*4 DAMP(NEBCM,NUMC-1),DG(NEBCM,NUMC-1)
     .       ,ETABC(300,NREC),Z0(NEBCM),ALAT(NEBCM),TZERO
      INTEGER DKON(7)
      CHARACTER ITYPE*4, KEQUI*4
      DATA KBLANK/'    '/,KEQUI/'EQUI'/
      DATA DKON/'S2  ','M2  ','N2  ','K1  ','P1  ','O1  ','Z0  '/
C
C
      IUPRT=10
      WRITE(IUPRT,*)'TZERO',TZERO
      DO N=1,NUMEBC
      DO K=1,NUMC-1
      WRITE(IUPRT,222)K,DKON(K),DAMP(N,K),DG(N,K),Z0(N),ALAT(N)
      END DO
      END DO
 222  FORMAT(I5,A4,3X,4F10.3)
      IF(NUMEBC.GT.300)THEN
        WRITE(IUPRT,*)'NUMBER OF STATIONS ARE LARGER THAN',300,' !!!!'
        WRITE(IUPRT,*)'INCREASE DIMENSION OF ETABC IN ptide.f ' !!!!'
        STOP
      END IF
      ITYPE='EQUI'
      DT=1.0
      PI=3.14159265358979D0
      TWOPI=2.D0*PI
      DO 20 I=1,2002
      ANGLE=TWOPI*(DBLE(FLOAT(I))-1.D0)/2000.D0
20    COSINE(I)=DCOS(ANGLE)
      IF(NDAY*24.GT.NREC)THEN
        WRITE(IUPRT,*)'NUMBER OF DATES EXCEEDING MAX SET IN ptide.f '
        WRITE(IUPRT,*)
     . 'Reduce NSTEP in run_data or increase NREC in ptide.f'
        STOP
      END IF
C
C READ ASTRO ARG , DOODSON AND NODAL CORRECTION DATA .
C ======================================================================
C
      CALL ASTRO(KONTAB,SIGTAB,V,U,F,MTAB,HOUR,XLAT,1)
      IF(MTAB.EQ.0) GO TO 2100
C
C READ STATION NUMBER, NAME, AND LOCATION.
C ======================================================================
C
      DO 1900 IBB=1,NUMEBC
      ISTN=IBB
      
50    FORMAT(6X,I4,1X,5A4,2X,A4,I2,1X,I2,2X,I3,1X,I2)
      IF(ISTN.EQ.0) GO TO 2100
C
C UNLESS SPECIFIED OTHERWISE, A STATION LATITUDE OF 50 DEGREES IS
C ASSUMED FOR THE PURPOSE OF CALCULATING SATELLITE TO MAIN CONSTITUENT
C AMPLITUDE RATIOS.
C=======================================================================
C
      XLAT=ALAT(IBB)
      IF(XLAT.LT.1.D0) XLAT=50.D0
C
C READ THE NAME, AMPLITUDE, AND GREENWICH PHASE LAG OF THE CONSTITUENTS
C TO BE USED IN THE PREDICTION.  THE CORRESPONDING FREQUENCIES ARE
C FOUND BY SCANNING THE LIST SIGTAB.  THE SEARCH METHOD WILL BE MOST
C EFFICIENT IF ON INPUT, THESE CONSTITUENTS ARE LISTED IN ORDER OF
C INCREASING FREQUENCY.
C=======================================================================
C
      KT1=1
      DO 110 K=1,NUMC
        KON(K)=DKON(K)
       IF(K.EQ.NUMC)THEN
        AMP(K)=Z0(IBB)
        G(K)=0.0
       ELSE
        AMP(K)=DAMP(IBB,K)
        G(K)=  DG(IBB,K)
       END IF
c70    FORMAT(5X,A4,13X,F6.1,12X,F5.1)
70    FORMAT(5X,A4,29X,F8.0,F7.0)
571    FORMAT(5X,A4,29X,F10.2,F10.2)
      IF(KON(K).EQ.KBLANK) GO TO 120
      G(K)=G(K)/360.D0
      KONK=KON(K)
      KTE=KT1-1+MTAB
      DO 80 KKT=KT1,KTE
      KT=KKT
      IF(KT.GT.MTAB) KT=KT-MTAB
      IF(KONTAB(KT).EQ.KONK) GO TO 100
80    CONTINUE
      WRITE(IUPRT,90) KONK
90    FORMAT(/' ROUTINE ASTRO CANNOT FIND CONSTITUENT ',A5 /)
      STOP 1
100   INDX(K)=KT
      SIG(K)=SIGTAB(KT)
      KT1=KT+1
110   CONTINUE
120   M=K-1
C
C READ START AND END DATES , TYPE OF PREDICTIONS AND TIME INCREMENT .
C ======================================================================
C
      DO 1800 ICC=1,1
130   FORMAT(3I3,1X,3I3,1X,A4,F17.1)
      IF(IDY0.EQ.0) GO TO 1900
      DT=DMIN1(6.0D0,DT)
C
C SET UP COEFFICIENTS FOR CHEBYSHEV ITERATION AND INITIALIZE .
C ======================================================================
C
135   CONTINUE
      LMAX=IDINT(0.99+DLOG(DT/0.1)*1.442695)   !*   ****************
C     LMAX=IFIX(0.99+DLOG(DT/0.1)*1.442695)   !*   ****************
      DO 150 K=1,M
      DRAD=TWOPI*SIG(K)*DT
      TWOC(K)=2.0D0*DCOS(DRAD)
      IF(ITYPE.EQ.KEQUI)  GO TO 150
      DO 140 L=1,LMAX
      DRAD=DRAD/2.0D0
      TPY=2.D0*DCOS(DRAD)
      BTWOC(K,L)=TPY
      IF(DABS(TPY).GT. 0.01D0) GO TO 140
      DT=DT*.99D0
      GO TO 135
140   CONTINUE
150   CONTINUE
C
C
      IMOP=0
      IMO=IMO0
      IYR=IYR0
      T=0.0D0
c
       CALL CDAY(IDY0,IMO0,IYR0,KD0,1)
c      ELSE
c       CALL CDAY(31,12,1999,KD1,iidint(1.D0))
c       CALL CDAY(IDY0,IMO0,IYR0,KD2,iidint(1.D0))
c       KD0=KD1+KD2
c      END IF
c
c     CALL CDAY(IDYE,IMOE,IYRE,KDE,iidint(1.D0))
      KDE=KD0+NDAY
      HOUR0=DBLE(FLOAT(KD0))*24.0D0
C
C HERE CALCULATIONS AND OUTPUT OF EQUI-SPACED PREDICTIONS ARE
C PERFORMED. THEY BEGIN AT DT HOURS INTO THE FIRST DAY AND ARE
C SPECIFIED EVERY DT HOURS THEREAFTER UP TO AND INCLUDING THE LAST
C DAY.  THE ASTRONOMICAL ARGUMENT V, AND NODAL MODULATION FACTORS,
C F AND U, AND RECALCULATED AT THE 16-TH OF EVERY MONTH ARE ASSUMED
C TO BE CONSTANT THROUGHOUT THE MONTH.
C=======================================================================
C
800   WRITE(IUPRT,805) ISTN
805   FORMAT('1', 8X,'EQUALLY SPACED PREDICTIONS FOR STATION',I5,' , ',5
     1A4)
c     WRITE(IUPRT,171) LAD,LAM,LOD,LOM,ITZONE,DT
      DO 900 IDD=1,NREC
      IF(IMO.EQ.IMOP) GO TO 830
      CALL CDAY(16,IMO,IYR,KDM,1)
      HOURM=DBLE(FLOAT(KDM))*24.0D0
      CALL ASTRO(KONTAB,SIGTAB,V,U,F,MTAB,HOURM,XLAT,2)
      TM=HOURM-HOUR0
      DO 810 K=1,M
      INDK=INDX(K)
      DTPY=F(INDK)*AMP(K)
      ANG=V(INDK)-(TM-T)*SIG(K)+U(INDK)-G(K)
      ANG=ANG-DINT(ANG)
      CHP(K)=DTPY*DCOS((ANG-DT*SIG(K))*TWOPI)
 810  CH(K)=DTPY*DCOS(ANG*TWOPI)
C
C
830   T=T+DT
      HGT=0.0D0
      DO 840 K=1,M
      DTPY=CH(K)
      CH(K)=TWOC(K)*DTPY-CHP(K)
      CHP(K)=DTPY
840   HGT=HGT+CH(K)
      TTOT=HOUR0+T
      KD=IDINT(TTOT/24.)
c     KD=IFIX(TTOT/24.)
      HR=TTOT-DBLE(FLOAT(KD))*24.D0
      IMOP=IMO
      CALL CDAY(IDY,IMO,IYR,KD,2)
      ETABC(IBB,IDD)=HGT
 3300 FORMAT(I5,4I3,5X,F10.5)
      IF(KD.GT.KDE) GO TO 910
900   CONTINUE
910   NCOUNT=IDD-1
C
C WRITE OUT THE LIST OF CONSTITUENTS USED.
C ======================================================================
C
1000  CONTINUE
      DO 1010 K=1,M
1010  G(K)=G(K)*360.D0
      WRITE(IUPRT,1020)
1020  FORMAT(/// ' NAME , AMPLITUDE AND GREENWICH PHASE LAG OF CONSTITUE
     $NTS USED TO GET THE ABOVE PREDICTIONS ' /)
      WRITE(IUPRT,1030) (K,KON(K),AMP(K),G(K),K=1,M)
1030  FORMAT(2(I5,1X,A4,F8.4,F7.2,5X))
      DO 1040 K=1,M
1040  G(K)=G(K)/360.D0
C
C
1800  CONTINUE
1900  CONTINUE
      DO N=IHR+1,NCOUNT
       WRITE(90,77)FLOAT(N)-IHR+TZERO
       WRITE(90,77)(ETABC(I,N),I=1,NUMEBC)
      END DO
2100  CONTINUE
 77   FORMAT(8E14.7)
      RETURN
      END
      SUBROUTINE ASTRO(KONT,FREQ,V,U,F,NTOTAL,HOUR,XLAT,MOD)
C     THIS ROUTINE CALCULATES THE FREQUENCY FREQ IN CPH , THE ASTRO ARG
C     V IN CYCLES , THE NODAL CORRECTION PHASE U AND AMP F FOR THE
C     CONSTITUENTS KON(1 ... NTOTAL).
C     THE INPUT PARAMETER HOUR GIVES THE TIME AT WHICH THE RESULTS ARE
C     TO BE CALCULATED. HOUR HAS ITS ORIGIN AROUND 1900 AND IS OBTAINED
C     USING ROUTINE CDAY.
C
C DATA INPUT
C ==========
C     THIS ROUTINE MUST BE INITIALIZED BY A CALL TO OPNAST WHICH READS
C     THE FOLLOWING FROM LOGICAL UNIT KR ...
C
C  1) TWO RECORDS CONTAINING ASTRONOMICAL ARGUMENTS AND THEIR RATES OF
C     CHANGE, S0,H0,P0,ENP0,PP0,DS,DH,DP,DNP,DPP, IN THE FORMAT
C     (5F13.10).
C      S0  = MEAN LONGITUDE OF THE MOON (CYCLES) AT 000 ET 1/1/1976.
C      H0  = MEAN LONGITUDE OF THE SUN.
C      P0  = MEAN LONGITUDE OF THE LUNAR PERIGEE.
C      ENP0= NEGATIVE OF THE MEAN LONGITUDE OF THE ASCENDING NODE.
C      PP0 = MEAN LONGITUDE OF THE SOLAR PERIGEE (PERIHELION).
C      DS,DH,DP,DNP,DPP ARE THEIR RESPECTIVE RATES OF CHANGE OVER A 365
C      DAY PERIOD AS OF 000 ET 1/1/1976.
C
C  2) ONE OR MORE RECORDS FOR EACH MAIN CONSTITUENT CONTAINING THE
C     FOLLOWING INFORMATION IN FORMAT (6X,A5,1X,6I3,F5.2,I4)
C      KON    = CONSTITUENT NAME
C   II,JJ,KK,LL,MM,NN = THE SIX DOODSON NUMBERS
C      SEMI   = PHASE CORRECTION
C      NJ     = THE NUMBER OF SATELLITES FOR THIS CONSTITUENT.
C   IF NJ%0, INFORMATION ON THE SATELLITE CONSTITUENTS IS READ , THREE
C   SATELLITES PER CARD, IN THE FORMAT (11X,3(3I3,F4.2,F7.4,1X,I1,1X)).
C   FOR EACH SATELLITE THE VALUES READ ARE
C      LDEL,MDEL,NDEL = THE CHANGES IN THE LAST THREE DOODSON NUMBERS
C                       FROM THOSE OF THE MAIN CONSTITUENT.
C      PH     = THE PHASE CORRECTION
C      EE     = THE AMPLITUDE RATIO OF THE SATELLITE TIDAL POTENTIAL TO
C               THAT OF THE MAIN CONSTITUENT.
C      IR     = 1 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
C                 LATITUDE CORRECTION FACTOR FOR DIURNAL CONSTITUENTS
C               2 IF THE AMPLITUDE RATIO HAS TO BE MULTIPLIED BY THE
C                 LATITUDE CORRECTION FACTOR FOR SEMI-DIURNAL CONSTI-
C                 TUENTS.
C               OTHERWISE IF NO CORRECTION IS REQUIRED TO THE AMPLITUDE
C                 RATIO.
C
C  3) ONE BLANK RECORD
C
C  4) ONE RECORD FOR EACH SHALLOW WATER CONSTITUENT SPECIFYING IN FORMAT
C     (11X,3(3I3,F4.2,F7.4,1X,I1,1X)), THE FOLLOWING INFORMATION
C      KON    = NAME  OF THE SHALLOW WATER CONSTITUENT
C      NJ     = NUMBER OF MAIN CONSTITUENTS FROM WHICH IT IS DERIVED.
C      COEF,KONCO = COMBINATION NUMBER AND NAME OF THESE MAIN
C                   CONSTITUENTS.
C
C  5) ONE BLANK RECORD
C
C MINIMAL DIMENSION REQUIREMENTS
C ==============================
C   THE DIMENSION OF KON, V,U, F, AND NJ SHOULD BE AT LEAST EQUAL TO THE
C   TOTAL NUMBER OF POSSIBLE CONSTITUENTS (PRESENTLY 146), THE DIMENSION
C   OF II, JJ, KK, LL, MM, NN AND SEMI SHOULD BE AT LEAST EQUAL TO THE
C   NUMBER OF MAIN CONSTITUENTS (PRESENTLY 45), THE DIMENSION OF EE,
C   LDEL, MDEL, NDEL, IR, AND PH SHOULD BE AT LEAST EQUAL TO THE TOTAL
C   NUMBER OF SATELLITES TO ALL THE MAIN CONSTITUENTS PLUS THE NUMBER
C   OF CONSTITUENTS WITH NO SATELLITES (PRESENTLY 162+8),
C   AND THE DIMENSION OF KONCO, AND COEFF SHOULD BE AT LEAST EQUAL TO
C   THE SUM FOR ALL SHALLOW WATER CONSTITUENTS OF THE NUMBER OF MAIN
C   CONSTITUENTS FROM WHICH EACH IS DERIVED (PRESENTLY 251).
C=======================================================================
C
      PARAMETER(NO1=45,NO2=146,NO3=170,NO4=251)
c     IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE
cjah  INTEGER*2 KR,LU,LP,MOD
      INTEGER*2 KR,LU,LP
      REAL*8 FREQ(NO2),V(NO2),U(NO2),F(NO2)
      REAL*8 II(NO1),JJ(NO1),KK(NO1),LL(NO1),MM(NO1),NN(NO1),SEMI(NO1)
      REAL*8 EE(NO3),LDEL(NO3),MDEL(NO3),NDEL(NO3),PH(NO3)
      COMMON/TIDECONS/
     .   COEF(NO4),II,JJ,KK,LL,MM,NN,
     .   SEMI,EE,LDEL,MDEL,NDEL,PH,
c    .   COEF(NO4),II(NO1),JJ(NO1),KK(NO1),LL(NO1),MM(NO1),NN(NO1),
c    .   SEMI(NO1),EE(NO3),LDEL(NO3),MDEL(NO3),NDEL(NO3),PH(NO3),
     .   KONTAB(NO2),NJ(NO2),KONCO(NO4),IR(NO3)
      DIMENSION KONT(NO2)
      COMMON/CONSNAME/KON(NO2)
      DATA KBLANK/'    '/
      DATA S0,H0,P0,ENP0,PP0,DS,DH,DP,DNP,DPP/
     .00.7428797055,0.7771900329,0.5187051308,0.3631582592,0.7847990160,
     .13.3594019864,0.9993368945,0.1129517942,0.0536893056,0.0000477414/

      IUPRT=10

      GO TO (1000,180),MOD                !*      GO TO 180
C
C CALL TO OPNAST READS REQUIRED DATA
C=======================================================================
C
1000  CONTINUE                 !*      ENTRY OPNAST
      CALL HARMCONS
      PI=3.1415926536D0
      TWOPI=2.D0*PI
C
C HERE INPUT DATA TYPE 1), THE ASTRONOMICAL VARIABLES AND THEIR
C RATES OF CHANGE ARE READ
C=======================================================================
C
 50   FORMAT(5F13.10)
      DTAU=365.D0+DH-DS
C
C HERE INPUT DATA TYPE 2), MAIN AND SATELLITE CONSTITUENT
C INFORMATION IS READ
C=======================================================================
C
      JBASE=0
      DO 90 K=1,NO1
60    FORMAT(6X,A4,2X,6F3.0,F5.2,I4)
      FREQ(K)=(II(K)*DTAU+JJ(K)*DS+KK(K)*DH+LL(K)*DP+MM(K)*DNP+
     1NN(K)*DPP)/(365.D0*24.D0)
      J1=JBASE+1
 75    JL=JBASE+NJ(K)
80    FORMAT((3X,3(3F3.0,F4.2,F7.4,1X,I1,1X)))
90    JBASE=JL
100   NTIDAL=K-1
   
C
C HERE INPUT DATA TYPE 4), SHALLOW WATER CONSTITUENT INFORMATION IS
C READ
C=======================================================================
C
      NJSUM=0
      JBASE=0
      K1=NTIDAL+1
      DO 160 K=K1,NO2
      J1=JBASE+1
      J4=J1+3
130   FORMAT(6X,A4,I2,2X,4(F5.2,A4,6X))
      IF(KON(K).EQ.KBLANK) GO TO 170
      JBASE=JBASE+NJ(K)
      FREQ(K)=0.D0
      DO 160 J=J1,JBASE
      DO 162 L=1,NTIDAL
      IF(KON(L).EQ.KONCO(J)) GO TO 163
162   CONTINUE
      WRITE(IUPRT,241) KONCO(J)
      STOP 2
163   FREQ(K)=FREQ(K)+COEF(J)*FREQ(L)
160   CONTINUE
170   NTOTAL=K-1
      DO N=1,NTIDAL
       KONT(N)=KON(N)
      END DO
      RETURN
C
C GIVEN THE PARAM HOUR , THE TABLES FREQ(K),V(K),U(K),F(K) ,
C K=1,NTOTAL ARE SET UP .
C=======================================================================
C
180   SLAT=DSIN(PI*XLAT/180.D0)
      CALL CDAY(1,1,1976,KD,1)
      YEARS=(HOUR/24.D0-DBLE(FLOAT(KD)))/365.D0
      S=S0+YEARS*DS
      H=H0+YEARS*DH
      P=P0+YEARS*DP
      ENP=ENP0+YEARS*DNP
      PP=PP0+YEARS*DPP
      INTDYS=HOUR/24.
c     INTDYS=jidint(HOUR/24.D0)
      HH=HOUR-DBLE(FLOAT(INTDYS))*24.D0
      TAU=HH/24.D0+H-S
      JBASE=0
      DO 210 K=1,NTIDAL
      VDBL=II(K)*TAU+JJ(K)*S+KK(K)*H+LL(K)*P+MM(K)*ENP+NN(K)*PP+SEMI(K)
      IV=IDINT(VDBL)
c     IV=IFIX(VDBL)
      IV=(IV/2)*2
      V(K)=VDBL-DBLE(FLOAT(IV))
      J1=JBASE+1
      JL=JBASE+NJ(K)
      SUMC=1.D0
      SUMS=0.D0
      DO 200 J=J1,JL
C
C HERE THE SATELLITE AMPLITUDE RATIO ADJUSTMENT FOR LATITUDE IS MADE
C=======================================================================
C
      RR=EE(J)
      L=IR(J)+1
      GO TO (901,902,903),L
902   RR=EE(J)*0.36309D0*(1.D0-5.D0*SLAT*SLAT)/SLAT
      GO TO 901
903   RR=EE(J)*2.59808D0*SLAT
901   CONTINUE
      UUDBL=LDEL(J)*P+MDEL(J)*ENP+NDEL(J)*PP+PH(J)
      IUU=IDINT(UUDBL)
c     IUU=IFIX(UUDBL)
      UU=UUDBL-DBLE(FLOAT(IUU))
      SUMC=SUMC+RR*DCOS(UU*TWOPI)
200   SUMS=SUMS+RR*DSIN(UU*TWOPI)
      F(K)=DSQRT(SUMC*SUMC+SUMS*SUMS)
      U(K)=DATAN2(SUMS,SUMC)/TWOPI
210   JBASE=JL
C
      JBASE=0
      K1=NTIDAL+1
      IF(K1.GT.NTOTAL) RETURN
      DO 270 K=K1,NTOTAL
      F(K)=1.0D0
      V(K)=0.D0
      U(K)=0.D0
      J1=JBASE+1
      JL=JBASE+NJ(K)
      DO 260 J=J1,JL
      DO 240 L=1,NTIDAL
      IF(KON(L)-KONCO(J))240,250,240
240   CONTINUE
      WRITE(IUPRT,241)KONCO(J)
241   FORMAT(/,' ROUTINE ASTRO UNABLE TO FIND CONSTIT ' ,A5 /)
      STOP 3
250   F(K)=F(K)*F(L)**DABS(COEF(J))
 2700 FORMAT('KON(L) and KONCO(J)',3I5,3x,2A4,F5.1)
      V(K)=V(K)+COEF(J)*V(L)
      U(K)=U(K)+COEF(J)*U(L)
260   CONTINUE
270   JBASE=JL
      RETURN
      END
      SUBROUTINE CDAY(IDD,IMM,IYY,KD,KULL)  
C
C******************* A DAY,MONTH AND YEAR(2 DIGITS) ARE GIVEN IN
C******************* ARGUMENTS 1-3.  THE DAY COUNT KD IS RETURNED.
C******************* THE ORIGIN IS SUCH THAT KD=367 ON 1/1/1901.
C******************* THIS ROUTINE SHOULD NOT BE USED OUTSIDE THE
C******************* RANGE 1901 TO 1999.
C
c     IMPLICIT INTEGER*4 (I-N)
      INTEGER*2 KR,LU,LP
      DIMENSION NDM(12),NDP(12)
      DATA NDP/0,31,59,90,120,151,181,212,243,273,304,334/
      DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/
      GO TO (100,200),KULL
100   IQUOT= IYY/4
      IREM = IYY-IQUOT*4
      L=1
      IF(IREM)3,1,3
1     IF(IMM-2)2,2,3
2     L=0
3     KD=IQUOT*1461+L+IREM*365+NDP(IMM)+IDD
      RETURN
C
C
200   CONTINUE                     !*      ENTRY DMY
C
C******************* GIVEN THE DAY COUNT KD , THE DAY,MONTH AND YEAR
C******************* ARE RETURNED .
C
      NDM(2)=28
      IQQ = (KD-1)/1461
      IYY = IQQ*4
      IDD = KD - IQQ*1461 -366
      IF(IDD)20,20,30
20    NDM(2)=29
      IDD =IDD +366
      GO TO 51
30    IYY = IYY +1
      IDD = IDD -365
      IF(IDD)50,50,30
50    IDD =IDD + 365
51    IMM = 0
60    IMM = IMM +1
      IDD = IDD -NDM(IMM)
      IF(IDD)70,70,60
70    IDD = IDD +NDM(IMM)
      RETURN
      END
      SUBROUTINE HARMCONS
      PARAMETER(NO1=45,NO2=146,NO3=170,NO4=251)
c     IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*2 KR,LU,LP,MOD
      REAL*8 II(NO1),JJ(NO1),KK(NO1),LL(NO1),MM(NO1),NN(NO1),SEMI(NO1)
      REAL*8 EE(NO3),LDEL(NO3),MDEL(NO3),NDEL(NO3),PH(NO3)
      COMMON/TIDECONS/
     .   COEF(NO4),II,JJ,KK,LL,MM,NN,
     .   SEMI,EE,LDEL,MDEL,NDEL,PH,
c    .   COEF(NO4),II(NO1),JJ(NO1),KK(NO1),LL(NO1),MM(NO1),NN(NO1),
c    .   SEMI(NO1),EE(NO3),LDEL(NO3),MDEL(NO3),NDEL(NO3),PH(NO3),
     .   KONTAB(NO2),NJ(NO2),KONCO(NO4),IR(NO3)
      COMMON/CONSNAME/KON(NO2)
      DIMENSION KONCO1(160),KONCO2(91)
       DATA II/
     .       0,   0,   0,   0,   0,   0,   0,   1,   1,   1,
     .       1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
     .       1,   1,   1,   1,   1,   1,   1,   2,   2,   2,
     .       2,   2,   2,   2,   2,   2,   2,   2,   2,   2,
     .       2,   2,   2,   2,   3/
       DATA JJ/
     .       0,   0,   0,   1,   1,   2,   2,  -4,  -3,  -3,
     .      -2,  -2,  -1,  -1,   0,   0,   0,   1,   1,   1,
     .       1,   1,   1,   2,   2,   3,   4,  -3,  -3,  -2,
     .      -2,  -1,  -1,   0,   0,   0,   0,   1,   1,   2,
     .       2,   2,   2,   3,   0/
       DATA KK/
     .       0,   1,   2,  -2,   0,  -2,   0,   2,   0,   2,
     .       0,   2,   0,   2,  -2,   0,   2,  -3,  -2,  -1,
     .       0,   1,   2,  -2,   0,   0,   0,   0,   2,   0,
     .       2,   0,   2,  -2,  -1,   0,   1,  -2,   0,  -3,
     .      -2,  -1,   0,   0,   0/
       DATA LL/
     .       0,   0,   0,   1,  -1,   0,   0,   1,   2,   0,
     .       1,  -1,   0,   0,   1,   1,  -1,   0,   0,   0,
     .       0,   0,   0,   1,  -1,   0,  -1,   3,   1,   2,
     .       0,   1,  -1,   2,   0,   0,   0,   1,  -1,   0,
     .       0,   0,   0,  -1,   0/
       DATA MM/
     .       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   0,   0,   0/
       DATA NN/
     .       0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   0,   0,   0,   0,   0,   1,   0,   1,
     .       0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   0,   0,   1,   0,  -1,   0,   0,   1,
     .       0,  -1,   0,   0,   0/
       DATA SEMI/
     .    0.00,0.00,0.00,0.00,0.00,0.00,0.00,-.25,-.25,-.25,
     .    -.25,-.25,-.25,-.75,-.75,-.75,-.75,-.25,-.25,-.75,
     .    -.75,-.75,-.75,-.75,-.75,-.75,-.75,0.00,0.00,0.00,
     .    0.00,0.00,0.00,-.50,-.50,0.00,0.00,-.50,-.50,0.00,
     .    0.00,-.50,0.00,0.00,-.50/
       DATA LDEL/
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-1.0, 0.0,-2.0,
     .    -1.0,-1.0, 0.0, 0.0,-1.0, 0.0, 0.0, 2.0,-2.0,-2.0,
     .    -1.0,-1.0,-1.0, 0.0,-1.0, 0.0, 1.0, 2.0, 0.0, 0.0,
     .     1.0, 2.0, 2.0,-1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0,
     .     2.0,-2.0,-1.0, 0.0, 0.0, 0.0, 0.0,-2.0,-2.0,-2.0,
     .    -1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 1.0, 2.0, 2.0, 0.0, 0.0,-2.0,-1.0,-1.0,
     .    -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,-2.0,-2.0,
     .     0.0, 0.0, 0.0,-2.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0,-2.0,-2.0,-2.0,
     .    -1.0,-1.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0, 1.0, 1.0,
     .    -1.0, 0.0,-1.0,-1.0, 0.0,-2.0,-1.0,-1.0, 0.0,-1.0,
     .    -1.0, 0.0,-2.0,-1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
     .    -2.0,-1.0, 0.0, 0.0, 1.0,-1.0,-1.0, 0.0, 0.0, 1.0,
     .     1.0, 1.0, 2.0, 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0,
     .     2.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0,-1.0,-1.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 0.0/
       DATA MDEL/
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-1.0,-2.0,
     .    -1.0, 0.0,-2.0,-1.0, 0.0,-2.0,-1.0, 0.0,-3.0,-2.0,
     .    -2.0,-1.0, 0.0,-2.0, 0.0,-1.0, 0.0, 0.0,-2.0,-1.0,
     .     0.0, 0.0, 1.0, 0.0,-2.0,-1.0,-1.0, 0.0, 1.0, 0.0,
     .     1.0, 0.0, 0.0,-1.0, 1.0, 2.0,-1.0,-2.0,-1.0, 0.0,
     .    -1.0, 0.0, 1.0,-1.0, 1.0, 2.0,-1.0, 1.0,-1.0,-2.0,
     .    -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,-1.0,-1.0, 0.0,
     .     1.0,-2.0,-1.0, 1.0, 2.0, 0.0, 1.0, 1.0, 0.0, 1.0,
     .     0.0, 1.0, 2.0,-1.0, 0.0,-1.0, 1.0,-1.0, 1.0, 2.0,
     .    -1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0,-1.0, 0.0, 1.0,
     .     0.0, 1.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 0.0, 1.0,
     .     0.0,-1.0,-1.0, 0.0,-1.0,-2.0,-1.0, 0.0,-1.0,-1.0,
     .     0.0,-1.0,-2.0, 0.0,-2.0,-1.0,-1.0, 0.0, 0.0, 1.0,
     .    -2.0, 0.0,-1.0,-1.0, 0.0,-1.0, 0.0,-2.0,-1.0,-1.0,
     .     0.0, 1.0, 0.0, 1.0,-1.0,-1.0,-1.0,-1.0, 0.0, 1.0,
     .     2.0, 0.0,-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,-1.0,
     .     1.0, 2.0,-1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0,-1.0/
       DATA NDEL/
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 2.0, 0.0, 0.0, 0.0,-2.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .    -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0,
     .     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
       DATA PH/
     .    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.75,0.00,0.50,
     .    0.75,0.75,0.50,0.00,0.75,0.50,0.00,0.50,0.50,0.50,
     .    0.75,0.75,0.75,0.50,0.00,0.00,0.75,0.50,0.50,0.00,
     .    0.75,0.50,0.00,0.25,0.50,0.00,0.25,0.75,0.25,0.50,
     .    0.50,0.00,0.25,0.50,0.50,0.50,0.00,0.50,0.00,0.00,
     .    0.75,0.25,0.75,0.50,0.00,0.50,0.50,0.00,0.50,0.00,
     .    0.50,0.50,0.75,0.50,0.50,0.00,0.50,0.00,0.75,0.25,
     .    0.75,0.00,0.50,0.00,0.50,0.25,0.25,0.00,0.00,0.00,
     .    0.00,0.50,0.50,0.00,0.25,0.50,0.00,0.50,0.00,0.50,
     .    0.75,0.25,0.25,0.25,0.50,0.50,0.50,0.50,0.00,0.00,
     .    0.25,0.25,0.00,0.00,0.00,0.00,0.00,0.00,0.25,0.25,
     .    0.25,0.50,0.25,0.25,0.50,0.50,0.25,0.25,0.50,0.25,
     .    0.25,0.50,0.50,0.00,0.00,0.50,0.50,0.75,0.00,0.50,
     .    0.00,0.25,0.50,0.50,0.50,0.75,0.75,0.00,0.50,0.25,
     .    0.75,0.75,0.00,0.00,0.50,0.50,0.50,0.00,0.50,0.50,
     .    0.50,0.00,0.00,0.75,0.00,0.50,0.00,0.75,0.75,0.50,
     .    0.00,0.00,0.50,0.00,0.00,0.75,0.75,0.75,0.50,0.50/
       DATA EE/
     . 0., 0., 0., 0., 0., 0., 0., 0.0360,.1906,.0063,.0241,.0607,.0063,
     ..1885,.0095,.0061,.1884,.0087,.0007,.0039,.0010,.0115,.0292,.0057,
     ..0008,.1884,.0018,.0028,.0058,.1882,.0131,.0576,.0175,.0003,.0058,
     ..1885,.0004,.0029,.0004,.0064,.0010,.0446,.0426,.0284,.2170,.0142,
     ..2266,.0057,.0665,.3596,.0331,.2227,.0290,.0290,.2004,.0054,.0282,
     ..2187,.0078,.0008,.0112,.0004,.0004,.0015,.0003,.3534,.0264,.0002,
     ..0001,.0007,.0001,.0001,.0198,.1356,.0029,.0002,.0001,.0190,.0344,
     ..0106,.0132,.0384, 0.0185, 0.0300, 0.0141, 0.0317, 0.1993, 0.0294,
     . 0.1980, 0.0047, 0.0027, 0.0816, 0.0331, 0.0027, 0.0152, 0.0098,
     . 0.0057, 0.0037, 0.1496, 0.0296, 0.0240, 0.0099, 0.6398, 0.1342,
     . 0.0086, 0.0611, 0.6399, 0.1318, 0.0289, 0.0257, 0.1042, 0.0386,
     . 0.0075, 0.0402, 0.0373, 0.0061, 0.0117, 0.0678, 0.0374, 0.0018,
     . 0.0104, 0.0375, 0.0039, 0.0008, 0.0005, 0.0373, 0.0373, 0.0042,
     . 0.0042, 0.0036, 0.1429, 0.0293, 0.0330, 0.0224, 0.0447, 0.0001,
     . 0.0004, 0.0005, 0.0373, 0.0001, 0.0009, 0.0002, 0.0006, 0.0002,
     . 0.0217, 0.0448, 0.0366, 0.0047, 0.2505, 0.1102, 0.0156, 0.0000,
     . 0.0022, 0.0001, 0.0001, 0.2535, 0.0141, 0.0024, 0.0004, 0.0128,
     . 0.2980, 0.0324, 0.0187, 0.4355, 0.0467, 0.0747, 0.0482, 0.0093,
     . 0.0078, 0.0564/
       DATA IR/
     .       0,   0,   0,   0,   0,   0,   0,   1,   0,   0,
     .       1,   1,   0,   0,   1,   0,   0,   0,   0,   0,
     .       1,   1,   1,   0,   0,   0,   1,   0,   0,   0,
     .       1,   0,   0,   1,   0,   0,   1,   1,   1,   0,
     .       0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
     .       1,   1,   1,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   1,   0,   0,   0,   0,   0,   1,   1,
     .       1,   0,   0,   0,   0,   1,   1,   0,   0,   0,
     .       0,   0,   0,   0,   1,   0,   0,   0,   0,   0,
     .       1,   1,   1,   1,   0,   0,   0,   0,   0,   0,
     .       1,   1,   0,   0,   0,   0,   0,   0,   1,   1,
     .       2,   0,   2,   2,   0,   0,   2,   2,   0,   2,
     .       2,   0,   0,   0,   0,   0,   0,   2,   0,   0,
     .       0,   2,   0,   0,   0,   2,   2,   0,   0,   2,
     .       2,   2,   0,   0,   0,   0,   0,   0,   0,   0,
     .       0,   0,   0,   2,   0,   0,   0,   2,   2,   0,
     .       0,   0,   0,   0,   0,   2,   2,   2,   0,   0/
       DATA KON/
     .    'Z0  ','SA  ','SSA ','MSM ','MM  ','MSF ','MF  ','ALP1',
     .    '2Q1 ','SIG1','Q1  ','RHO1','O1  ','TAU1','BET1','NO1 ',
     .    'CHI1','PI1 ','P1  ','S1  ','K1  ','PSI1','PHI1','THE1',
     .    'J1  ','OO1 ','UPS1','OQ2 ','EPS2','2N2 ','MU2 ','N2  ',
     .    'NU2 ','GAM2','H1  ','M2  ','H2  ','LDA2','L2  ','T2  ',
     .    'S2  ','R2  ','K2  ','ETA2','M3  ','2PO1','SO1 ','ST36',
     .    '2NS2','ST37','ST1 ','ST2 ','ST3 ','O2  ','ST4 ','SNK2',
     .    'OP2 ','MKS2','ST5 ','ST6 ','2SK2','MSN2','ST7 ','2SM2',
     .    'ST38','SKM2','2SN2','NO3 ','MO3 ','NK3 ','SO3 ','MK3 ',
     .    'SP3 ','SK3 ','ST8 ','N4  ','3MS4','ST39','MN4 ','ST40',
     .    'ST9 ','M4  ','ST10','SN4 ','KN4 ','MS4 ','MK4 ','SL4 ',
     .    'S4  ','SK4 ','MNO5','2MO5','3MP5','MNK5','2MP5','2MK5',
     .    'MSK5','3KM5','2SK5','ST11','2NM6','ST12','ST41','2MN6',
     .    'ST13','M6  ','MSN6','MKN6','2MS6','2MK6','NSK6','2SM6',
     .    'MSK6','ST42','S6  ','ST14','ST15','M7  ','ST16','3MK7',
     .    'ST17','ST18','3MN8','ST19','M8  ','ST20','ST21','3MS8',
     .    '3MK8','ST22','ST23','ST24','ST25','ST26','4MK9','ST27',
     .    'ST28','M10 ','ST29','ST30','ST31','ST32','ST33','M12 ',
     .    'ST34','ST35'/
       DATA NJ/
     .       1,   1,   1,   1,   1,   1,   1,   2,   5,   4,
     .      10,   5,   8,   5,   1,   9,   2,   1,   6,   2,
     .      10,   1,   5,   4,  10,   8,   5,   2,   3,   4,
     .       3,   4,   4,   3,   2,   9,   1,   1,   5,   1,
     .       3,   2,   5,   7,   1,   2,   2,   3,   2,   2,
     .       3,   4,   3,   1,   3,   3,   2,   3,   3,   4,
     .       2,   3,   4,   2,   3,   3,   2,   2,   2,   2,
     .       2,   2,   2,   2,   3,   1,   2,   4,   2,   3,
     .       4,   1,   3,   2,   2,   2,   2,   2,   1,   2,
     .       3,   2,   2,   3,   2,   2,   3,   3,   2,   3,
     .       2,   4,   3,   2,   4,   1,   3,   3,   2,   2,
     .       3,   2,   3,   3,   1,   3,   3,   1,   3,   2,
     .       4,   2,   2,   4,   1,   3,   3,   2,   2,   4,
     .       2,   3,   3,   3,   2,   3,   2,   1,   3,   2,
     .       4,   2,   3,   1,   2,   4/
       DATA COEF/
     .     2.,-1., 1.,-1., 2., 1.,-2., 2.,-1., 3.,-2., 2., 1.,-2., 1.,
     .     1., 1.,-2., 2., 1.,-2., 2., 2., 1.,-2., 1., 1.,-1., 1., 1.,
     .     1., 1.,-1., 1., 2.,-2., 2., 1.,-1.,-1., 2.,-1., 1., 1.,-1.,
     .     2., 1.,-1.,-1., 2.,-1., 2., 1.,-2., 1., 1.,-1., 2.,-1., 1.,
     .     1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 2., 1.,
     .    -1., 2., 3.,-1., 1., 1., 1.,-1., 1., 1., 2., 1.,-1., 1., 1.,
     .     1.,-1., 2., 2., 1.,-1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
     .     1., 2., 1., 1., 1., 1., 1., 2., 1., 3.,-1., 1., 1., 1., 2.,
     .     1., 2., 1., 1., 1., 1., 1., 1., 1., 2., 1., 3., 1.,-1., 2.,
     .     1., 2., 1., 1.,-1., 3., 1.,-1., 2., 1., 2., 1., 1.,-1., 3.,
     .     1., 1., 1., 1., 1., 1., 2., 1., 2., 1., 1., 1., 1., 2., 1.,
     .     1., 1., 1., 2., 2.,-1., 3., 2., 1., 1., 2., 1., 1., 3., 2.,
     .     1., 1., 3., 1., 1., 1., 1., 1., 2., 2., 3., 1., 3., 1., 1.,
     .    -1., 4., 2., 1., 1., 2., 1., 1., 3., 1., 3., 1., 1., 1., 1.,
     .     1., 2., 2., 2., 1., 1., 2., 2., 1., 3.,
     .     1., 1., 4., 1., 3., 1., 1., 4., 1., 5.,
     .     3., 1., 1., 4., 1., 2., 1., 1., 1., 3.,
     .     2., 4., 1., 1., 6., 5., 1., 3., 1., 1.,
     .     1./
       DATA KONCO1/'P1  ','O1  ','S2  ','O1  ','M2  ','N2  ','S2  ',
     .'N2  ','S2  ','M2  ','S2  ','N2  ','K2  ','S2  ','M2  ','N2  ',
     .    'K2  ','S2  ','M2  ','S2  ','K2  ','O1  ','K2  ','N2  ',
     .    'S2  ','S2  ','N2  ','K2  ','O1  ','P1  ','M2  ','K2  ',
     .    'S2  ','M2  ','K2  ','S2  ','S2  ','N2  ','M2  ','K2  ',
     .    'S2  ','K2  ','M2  ','S2  ','N2  ','K2  ','M2  ','S2  ',
     .    'N2  ','S2  ','M2  ','M2  ','S2  ','N2  ','S2  ','K2  ',
     .    'M2  ','S2  ','N2  ','N2  ','O1  ','M2  ','O1  ','N2  ',
     .    'K1  ','S2  ','O1  ','M2  ','K1  ','S2  ','P1  ','S2  ',
     .    'K1  ','M2  ','N2  ','S2  ','N2  ','M2  ','S2  ','M2  ',
     .    'S2  ','N2  ','K2  ','M2  ','N2  ','M2  ','S2  ','K2  ',
     .    'M2  ','N2  ','K2  ','S2  ','M2  ','M2  ','K2  ','S2  ',
     .    'S2  ','N2  ','K2  ','N2  ','M2  ','S2  ','M2  ','K2  ',
     .    'S2  ','L2  ','S2  ','S2  ','K2  ','M2  ','N2  ','O1  ',
     .    'M2  ','O1  ','M2  ','P1  ','M2  ','N2  ','K1  ','M2  ',
     .    'P1  ','M2  ','K1  ','M2  ','S2  ','K1  ','K2  ','K1  ',
     .    'M2  ','S2  ','K1  ','N2  ','K2  ','S2  ','N2  ','M2  ',
     .    'N2  ','M2  ','K2  ','S2  ','M2  ','S2  ','K2  ','M2  ',
     .    'N2  ','M2  ','N2  ','K2  ','S2  ','M2  ','M2  ','S2  ',
     .    'N2  ','M2  ','K2  ','N2  ','M2  ','S2  ','M2  ','K2  '/
       DATA KONCO2/
     .    'N2  ','S2  ','K2  ','S2  ','M2  ','M2  ','S2  ','K2  ',
     .    'M2  ','S2  ','K2  ','S2  ','M2  ','N2  ','O1  ','N2  ',
     .    'M2  ','K1  ','M2  ','M2  ','S2  ','O1  ','M2  ','K1  ',
     .    'M2  ','S2  ','K2  ','O1  ','M2  ','N2  ','M2  ','N2  ',
     .    'M2  ','N2  ','K2  ','S2  ','M2  ','M2  ','S2  ','N2  ',
     .    'M2  ','N2  ','K2  ','M2  ','S2  ','M2  ','K2  ','M2  ',
     .    'S2  ','N2  ','K2  ','M2  ','S2  ','M2  ','S2  ','K2  ',
     .    'M2  ','N2  ','K1  ','M2  ','N2  ','K1  ','M2  ','K1  ',
     .    'M2  ','S2  ','K1  ','M2  ','N2  ','M2  ','M2  ','N2  ',
     .    'S2  ','M2  ','S2  ','M2  ','N2  ','S2  ','K2  ','M2  ',
     .    'S2  ','M2  ','S2  ','K1  ','M2  ','M2  ','S2  ','M2  ',
     .    'N2  ','K2  ','S2  '/
       DO I=1,160
        KONCO(I)=KONCO1(I)
       END DO
       DO I=1,91
        KONCO(I+160)=KONCO2(I)
       END DO
      RETURN
      END
