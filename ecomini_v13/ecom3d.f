      PROGRAM ECOM_3D
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
C*************************************************************************
*                                                                          |
*                                                                          |
*  The development of ECOMSED has its origins in the mid 1980's with       |
*  the creation of the Princeton Ocean Model by George Mellor and me       |
*  at Princeton University and its version for shallow water environments  |
*  - rivers, bays, estuaries and the coastal ocean and reservoirs and      |
*  lakes- named ECOM at HydroQual.  In the mid 1990s, concepts for         |
*  cohesive sediment resuspension, settling and consolidation were         |
*  incorporated within the ECOM modeling framework.                        |
*                                                                          |
*  Todas's version of ECOMSED have been made possible by  the dedicated    |
*  efforts of P. L. Shrestha, B. N. Kim, Q. Ahsan and H. Li. They have     |
*  helped conceive, design and implement the model enhancements and        |
*  worked diligently to debug their (and my) changes. Important            |
*  contributions to various aspects of ECOMSED have also been made         |
*  at one time or another by Boris Galperin,  H. James Herring,            |
*  Eugenio Gomez-Reyes, C. Kirk Zeigler, and Richard P. Signell.           |
*  Finally, the seminal contributions of George L. Mellor must be          |
*  acknowledged. It was he who first managed to secure funding which       |
*  made the Princeton Ocean Model a reality.                               |
*                                                                          |
* Point of Contact: ecomsed-support@hydroqual.com                          |
*__________________________________________________________________________|
C
C                                                                      
*__________________________________________________________________________
C                                                                          |
C                   LOOP LIMITS                                            |
C                                                                          |
C              T,S,etc :  J=2,JMM1                                         |
C                         I=2,IMM1                                         |
C                                                                          |
C                 U    :  J=2,JMM1                                         |
C                         I=3,IMM1                                         |
C                                                                          |
C                 V    :  J=3,JMM1                                         |
C                         I=2,IMM1                                         |
C                                                                          |
C__________________________________________________________________________|
C

CW-------------------------------2017-07-02---------------------------------
CW-----------------------------王永桂测试修改-------------------------------



      INCLUDE 'comdeck'
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION COM(80)
      DIMENSION IVAR(IM),PRT(IM,KB)
C
      REAL DUMB1(QBCM),DUMB2(EBCM,KBM1),DUMB3(DBCM)
      DIMENSION DUMP1(IM,JM),DUMP2(IM,JM,KB)
      INTEGER TMP_IYR,TMP_IMO,TMP_IDA,TMP_IHOUR
C
C  LIMIT DEPOSITION TO BATHYMETRIC DEPTH 
C
      REAL COHTHK(IM,JM),NCOHTHK(IM,JM)
C
      EQUIVALENCE (IVAR,TPS),(PRT,A)
C
      CHARACTER*10 RESTAR,OPTION
c   use common block when HOT START and EXTERNAL
      COMMON /COAST1/ICNT,INDX(MAXWET),JNDX(MAXWET),RESTAR
*
       DIMENSION TBDRY2(EBCM,KBM1), SBDRY2(EBCM,KBM1)
*
C 
      IURUN=1                 ! run_data
      IUGRD=3                 ! model_grid (opened & closed in setdom.f)
      IUTAS=5                 ! init_tands (opened & closed in tands.f)
      IUUAV=7                 ! synop_wind (opened & closed in bcdata.f)
      IURRS=9                 ! restart (opened & closed in ecom3d.f)
C
      IUPRT=10                ! gcmprt
      IUPLT=12                ! gcmplt
      IUTSR=14                ! gcmtsr
      IUWRS=16                ! startup (opened & closed in ecom3d.f)
      IUTRN=18                ! gcm_geom (opened & closed in transport.f)
C                             ! gcm_tran (opened & closed in transport.f)
C
      IUT90=90                ! elevation boundary conditions
      IUT91=91                ! river discharges
      IUT92=92                ! diffuser intake/outfall
      IUT93=93                ! meteorological data
      IUT94=94                ! temperature & salinity boundary conditions
      IUT95=95                ! synoptic wind stress comp (opened in bcdata.f)
      IUT96=96                ! diffuser intake/outfall in loop
      IUT191=191              ! synoptic heat flux input file
      IUT192=192              ! synoptic heat flux (temp file)
      IUT193=193              ! synoptic wind velocity components 
c                               (opened in bcdata.f) used by wave model
C
C********************************************************************
C
C  THESE FILES ARE OPENED IN bcdata.f
C
C  DISSOLVED TRACER TRANSPORT
C
      IUT501=501                ! tracer conc. at open b.c.
      IUT601=601                ! tracer conc. at river discharges
      IUT98=98                  ! tracer conc. at diffuser intake/outfall
      IUT99=99                  ! tracer conc. at diffuser intake/outfall loop
      IUT701=701                ! tracer load at point source
C
C  SEDIMENT TRANSPORT
C
      IUT502=502                !  coh sed conc. at open b.c.
      IUT503=503                !  non-coh sed conc. at open b.c.
      IUT602=602                !  coh sed conc. at river discharges
      IUT603=603                !  non-coh sed conc.river discharges
      IUT702=702                !  coh sed conc. at Diffuser discharges
      IUT703=703                !  non-coh sed conc Diffuser discharges
C
C  PARTICLE-BOUND TRACER TRANSPORT
C
      IUT504=504                !  coh sed conc. at open b.c.
      IUT505=505                !  non-coh sed conc. at open b.c.
      IUT604=604                !  coh sed conc. at river discharges
      IUT605=605                !  non-coh sed conc.river discharges
      IUT704=704                !  coh sed conc. at Diffuser discharges
      IUT705=705                !  non-coh sed conc Diffuser discharges
C
C********************************************************************
C
C

      CALL GETARG(1,USERPATH)   !!!传入用户路径



      OPEN (IURUN,FILE=trim(USERPATH)//'\INC\run_data')
      OPEN (IUPRT,FILE=trim(USERPATH)//'\PRT_RES\gcmprt')
      OPEN (IUT90,FILE=trim(USERPATH)//'\TEMP\gcm_temp90')
      OPEN (IUT91,FILE=trim(USERPATH)//'\TEMP\gcm_temp91')
      OPEN (IUT92,FILE=trim(USERPATH)//'\TEMP\gcm_temp92')
      OPEN (IUT93,FILE=trim(USERPATH)//'\TEMP\gcm_temp93')
c      OPEN (IUT193,FILE='gcm_temp193')
      OPEN (IUT94,FILE=trim(USERPATH)//'\TEMP\gcm_temp94')
      OPEN (IUT96,FILE=trim(USERPATH)//'\TEMP\gcm_temp96')
C
      WRITE(IUPRT,7000)
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,11)  (COM(I),I=1,80)
 11   FORMAT(80A1)
 12   FORMAT(/1X,80A1/)
C
C********************************************************************
C
C  HYDTYPE = 'INTERNAL' :  USE ECOM HYDRODYNAMICS
C          = 'EXTERNAL' :  USE POM HYDRODYNAMICS INPUT FROM hqi_tran
C
C  NHYD = NO. OR TIMESTEPS BETWEEN EACH hqi_tran INPUT (ASSUMED CONSTANT)
C
C  WAVEDYN = 'NEGLECT ' ==>  NO EFFECT OF WAVES ON BOTTOM FRICTION
C          = 'SMBMODEL' ==>  WAVES EFFECT BOT FRIC COEFF - SMB THEORY
C          = 'DONMODEL' ==>  WAVES EFFECT BOT FRIC COEFF - Donelan(1977) THEORY
C          = 'EXTERNAL' ==>  WAVES EFFECT BOT FRIC COEFF - INPUT WAM RESULTS
C
C  TRACER = 'NEGLECT' :  NO CONSERVATIVE TRACER CALC.
C         = 'INCLUDE' :  CONSERVATIVE TRACER
C
C  SEDTRAN = 'NEGLECT' :  NO SEDIMENT TRANSPORT
C          = 'INCLUDE' :  SEDIMENT TRANSPORT
C
C  CHEMTRAN = 'NEGLECT' :  NO CHEMICAL TRANSPORT (PARTICLE BOUND)
C           = 'INCLUDE' :  CHEMICAL TRANSPORT (PARTICLE BOUND)
C                           (ONLY USE IF SEDTRAN = 'INCLUDE')
C
C  PARTICLE = 'NEGLECT' :  NO PARTICLE TRACKING
C           = 'INCLUDE' :  PARTICLE TRACKING
C
C  SEDTYPE = 'BOTH'  :  COHESIVE & NON-COHESIVE
C          = 'SAND'  :  NON-COHESIVE ONLY
C          = 'MUD '  :  COHESIVE ONLY
C
C  NOTE:  SEDTYPE APPLIES TO SEDIMENT TRANSPORT AND PARTICLE-BOUND
C         CHEMICAL (TRACER) TRANSPORT
C
      READ(IURUN,1) HYDTYPE,WAVEDYN,TRACER,SEDTRAN,CHEMTRAN,SEDTYPE,
     +                PARTICLE
 1    FORMAT(2X,A8,2X,A8,3X,A7,3X,A7,3X,A7,6X,A4,3X,A7)
C
          IF (WAVEDYN.NE.'SMBMODEL'.AND.WAVEDYN.NE.'DONMODEL'.AND.
     .        WAVEDYN.NE.'EXTERNAL'.AND.WAVEDYN.NE.'NEGLECT ')THEN
              WRITE(IUPRT,6112)WAVEDYN
              !!!&&&CALL SYSTEM ('rm gcm_temp*')
              STOP
          ENDIF
c
C
C  OPEN AUXILIARY INPUT FILES
C
C  FILE                UNIT      DESCRIPTION
C  water_trace.inp     IUT401    DISSOLVED TRACER
C  coh_sed.inp         IUT402    SEDIMENT TRANSPORT, COHESIVE
C  noncoh_sed.inp      IUT403    SEDIMENT TRANSPORT, NON-COHESIVE
C  coh_trace.inp       IUT404    PARTICLE-BOUND TRACER, COHESIVE
C  noncoh_trace.inp    IUT405    PARTICLE-BOUND TRACER, NON-COHESIVE
C  partrack.inp        IUT406    PARTICLE TRACKING
C
      IF (TRACER.EQ.'INCLUDE') THEN 
        IUT401=401
        OPEN (UNIT=IUT401,FILE='water_trace.inp',FORM='FORMATTED')
      ENDIF
C
      IF (SEDTRAN.EQ.'INCLUDE') THEN 
        IF (SEDTYPE.EQ.'MUD ') THEN
          IUT402=402
          OPEN (UNIT=IUT402,FILE='coh_sed.inp',FORM='FORMATTED')
C
          KSED=1
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND') THEN
          IUT403=403
          OPEN (UNIT=IUT403,FILE='noncoh_sed.inp',FORM='FORMATTED')
C
          KSED=1
        ENDIF
C
        IF (SEDTYPE.EQ.'BOTH') THEN
          IUT402=402
          IUT403=403
          OPEN (UNIT=IUT402,FILE='coh_sed.inp',FORM='FORMATTED')
          OPEN (UNIT=IUT403,FILE='noncoh_sed.inp',FORM='FORMATTED')
C
          KSED=2
        ENDIF
      ENDIF
C
      IF (CHEMTRAN.EQ.'INCLUDE') THEN 
        IF (SEDTYPE.EQ.'MUD ') THEN
          IUT404=404
          OPEN (UNIT=IUT404,FILE='coh_trace.inp',FORM='FORMATTED')
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND') THEN
          IUT405=405
          OPEN (UNIT=IUT405,FILE='noncoh_trace.inp',FORM='FORMATTED')
        ENDIF
C
        IF (SEDTYPE.EQ.'BOTH') THEN
          IUT404=404
          IUT405=405
          OPEN (UNIT=IUT404,FILE='coh_trace.inp',FORM='FORMATTED')
          OPEN (UNIT=IUT405,FILE='noncoh_trace.inp',FORM='FORMATTED')
        ENDIF
      ENDIF
C
      IF (PARTICLE.EQ.'INCLUDE') THEN 
        IUT406=406
        OPEN (UNIT=IUT406,FILE='partrack.inp',FORM='FORMATTED')
      ENDIF
C
C********************************************************************
C
C
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,2) DTI,ISPLIT,IRAMP,IYR,IMO,IDA,IHR,NHYD
      IHOUR=IHR
 2    FORMAT(1F10.4,7I5)
C
      CALL CDAY(IDA,IMO,IYR,IJDAY,1)
      SDAY=float(IJDAY)+FLOAT(IHOUR)/24.
      CALL CDAY(1,1,IYR,J1JDAY,1)
      J1YR = IYR
      CALL CDAY(1,1,J1YR+1,J2JDAY,1)   
C
C
      DTE=DTI/FLOAT(ISPLIT)
      DTE2=2.*DTE
      DTI2=2.*DTI
      ISPI=1./FLOAT(ISPLIT)
      ISP2I=.5*ISPI
      DAYI=1./86400.
      GRAV=9.806
C 
      WRITE(IUPRT,21) DTI,DTE,ISPLIT,IRAMP,IYR,IMO,IDA,IHOUR 
 21   FORMAT(
     . ' BAROCLINIC TIME STEP IS                 ',F10.4,' SECONDS',//,
     . ' BAROTROPIC TIME STEP IS                 ',F10.4,' SECONDS',//,
     . ' INTERNAL/EXTERNAL MODE SPLITTING IS     ',I10//,
     . ' NUMBER OF RAMP TIME STEPS               ',I10//,
     & ' THE STARTING TIME IS  ', 4i5//)  
 211   FORMAT(
     & ' HOT START THE STARTING TIME IS  ', 4i5//)
C
      WRITE (IUPRT,5002)TRACER,SEDTRAN,CHEMTRAN,SEDTYPE,PARTICLE,HYDTYPE
     +  ,WAVEDYN
 5002 FORMAT (//5X,'MODEL OPTIONS: TRACER   =',3X,A7,/20X,'SEDTRAN =',
     + 3X,A7,/20X,'CHEMTRAN =',3X,A7/20X,'SEDTYPE =',6X,A4,
     + /20X,'PARTICLE =',3X,A7,/20X,'HYDTYPE =',2X,A8,
     + /20X,'WAVEDYN =',2X,A8//)
C
C
C
C-----------------------------------------------------------------------
C         TYPE OF RUN -
C      BAROTROPIC: 2-D CALCULATION (BOTTOM STRESS CALCULATED IN ADVAVE)
C      PROGNOSTIC: 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
C
C      TEMP_ONLY : 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
C                : ONLY TEMPERATURE IS CALCULATED
C      SALT_ONLY : 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
C                : ONLY SALINITY IS CALCULATED
C      DIAGNOSTIC: 3-D CALCULATION WITH T AND S HELD FIXED
C----------------------------------------------------------------------- 
C         3-D - TYPE OF MOMENTUM ADVECTION AND BOTTOM FRICTION
C      LINEAR    : ALL MOMENTUM ADVECTION NEGLECTED
C      NON-LINEAR: COMPLETE PHYSICS 
C----------------------------------------------------------------------  
C         3-D - TYPE OF ADVECTION INTEGRATION SCHEME
C      CENTRAL   : CENTRAL FINITE DIFFERENCE SCHEME
C      UPWIND    : UPWIND FINITE DIFFERENCE SCHEME
C      SMOLAR_R  : FINITE DIFFERENCE SCHEME DUE TO SMOLARKIEWICZ USING
C                  RECURSIVE FORMULATION FOR THE ANTIDIFFUSIVE VELOCITIES
C                      (MOMENTUM AND TURBULENCE ARE UPWIND)
C                      ( to be implemented, now central)
C      SMOLAR_2  : FINITE DIFFERENCE SCHEME DUE TO SMOLARKIEWICZ USING
C                  TWO PASSES FOR CORRECTIONS OF THE NUMERICAL DIFFUSION
C                      (MOMENTUM AND TURBULENCE ARE UPWIND)
C                      ( to be implemented, now central)
C----------------------------------------------------------------------  
C
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,10) NSTEPS,IPRINT,IPRTSTART,RESTAR,TOR,ADVECT,SCHEME
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,4)  DEV,VSX,JROW,VSY,IROW,
     .               PTU,PTV,PTW,PTAM,PTS,PTT,PRHO,PTQ2,PTL,PTKM,PTKH
      READ(IURUN,11)  (COM(I),I=1,80)

      TNDAYS=FLOAT(NSTEPS)*DAYI*DTI
      EDAY=SDAY+TNDAYS
C
C  WAVEHYD = 'NEGLECT' ==>  NO EFFECT OF WAVES ON BOTTOM FRICTION
C          = 'INLCUDE' ==>  WAVES EFFECT BOTTOM FRICTION COEFF.
C
C  NWAVE = NUMBER OF TIMESTEPS BETWEEN UPDATING NEW BOTTOM FRICTION COEFF.
C

      READ(IURUN,3)  BFRIC,Z0B,NU,THETA,ALPHA,TLAG,NWAVE,BCTYPE
                IF(BCTYPE.EQ.'       ')BCTYPE='CLAMPED'

      IF(BCTYPE.NE.'IRANDB '.AND.BCTYPE.NE.'OCLAMP '.AND.
     1  BCTYPE.NE.'CLAMPED'.AND.BCTYPE.NE.'PCLAMP '.AND.
     2  BCTYPE.NE.'RANDB  '.AND.BCTYPE.NE.'MIXED  ') THEN
                   WRITE(IUPRT,6110) BCTYPE
                   GOTO 9200
                END IF
                IF(BCTYPE.EQ.'PCLAMP '.AND.TLAG.EQ.0.0) THEN
                   WRITE(IUPRT,6111) BCTYPE,TLAG
                   GOTO 9200
                ENDIF
C
      NWAVECNT=0
C
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,23) HORZMIX,HORCON,HPRNU
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,23) VERTMIX,UMOL,VPRNU
 4    FORMAT(2X,A3,2(4X,A1,I5),11(4X,A1))
 3    FORMAT(6E10.3,I10,3X,A7)
 10   FORMAT(3I10,1X,A10,1X,A10,1X,A10,1X,A10)
 23   FORMAT(A10,2E10.3)
C
C*****************************************************************
C
C  CONSERVATIVE TRACER INPUT SECTION
C  CONDRAT = TRACER DECAY RATE (1/day)
C  CONINIT = INITIAL TRACER CONCENTRATION (ASSUMED SPATIALLY CONSTANT)
C
      IF (TRACER.EQ.'INCLUDE') THEN
        READ (IUT401,11) (COM(I),I = 1,80)
        READ (IUT401,5009)CONDRAT,CONINIT
 5009   FORMAT (2F10.0)
C
        WRITE (IUPRT,11) (COM(I),I = 1,80)
        WRITE (IUPRT,5119)CONDRAT,CONINIT
 5119   FORMAT (2F10.2)
C
C  CONVERT DECAY RATE FROM 1/day TO 1/s
C
        CONDRAT=CONDRAT/86400.
      ENDIF
C
C**********************************************************************
C
C  SEDIMENT TRANSPORT INPUT SECTION
C
      IF (SEDTRAN.EQ.'INCLUDE') THEN
C
C  COHESIVE SEDIMENT INPUT
C
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
C
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5010)NSEDBEG,NSBED
 5010     FORMAT (3I10)
C
          XNSBED=FLOAT(NSBED)
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5010)NSEDBEG,NSBED
C
C  WS1 = ADEP * (C*G) ** DEPEXP
C
C  ADEP = CONSTANT, CLASS 1 SETTLING SPEED IN um/s
C
C  PDEPFORM = 'KRONE'
C           = 'PARTH'
C
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5019)ADEP,DEPEXP,TCRDEP,PDEPFORM
 5020     FORMAT (8F10.0)
 5019     FORMAT (3F10.0,5X,A5)
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5019)ADEP,DEPEXP,TCRDEP,PDEPFORM
C
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5018)A0IN,RESEXP,EXPM,VARIA0N
 5018     FORMAT (3F10.0,3X,A7)
C
C
          IA0=0
          IF (VARIA0N.EQ.'INCLUDE') IA0=1
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020)A0IN,RESEXP,EXPM
C
          READ (IUT402,11) (COM(I),I = 1,80)
C
C  ADDED Z0 FOR WAVES
C
          READ (IUT402,5017)DENCOH,VARIBULK,P0(1),VARIP0,BFCOH,Z0BCOH,
     +                     Z0WAVE
 5017     FORMAT (F10.0,3X,A7,F10.0,3X,A7,3F10.0)
C
          IP0=0
          IF (VARIP0.EQ.'INCLUDE') IP0=1
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5017)DENCOH,VARIBULK,P0(1),VARIP0,BFCOH,Z0BCOH,
     +                      Z0WAVE
C
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) (FTIME(LL),LL=1,LAYMAX)
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) (FTIME(LL),LL=1,LAYMAX)
C
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) (TSED0IN(LL),LL=1,LAYMAX)
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) (TSED0IN(LL),LL=1,LAYMAX)
C
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) (TAUCR(LL),LL=1,LAYMAX)
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) (TAUCR(LL),LL=1,LAYMAX)
C
C  INPUT INITIAL SUSPENDED SEDIMENT CONCENTRATIONS
C
C  ASSUMED TO BE SPATIALLY CONSTANT
C
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) CSI(1)
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) CSI(1)
        ENDIF
C
C  NON-COHESIVE SEDIMENT INPUT
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
C
          READ (IUT403,11) (COM(I),I = 1,80)
          READ (IUT403,5010)NSEDBEG,NSBED
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5010)NSEDBEG,NSBED
C
          READ (IUT403,11) (COM(I),I = 1,80)
          READ (IUT403,5021)WS2,DENNON,VARIBULK,SUSARM,BEDTHI
 5021     FORMAT (2F10.0,3X,A7,2F10.0)
C
C  CONVERT INPUT NON-COHESIVE BED THICKNESS FROM cm TO m
C
          BEDTHI=BEDTHI/100.
C
          XNSBED=FLOAT(NSBED)
C
C
C  INPUT INITIAL SUSPENDED SEDIMENT CONCENTRATIONS
C
C  ASSUMED TO BE SPATIALLY CONSTANT
C
          READ (IUT403,11) (COM(I),I = 1,80)
          READ (IUT403,5020) CSI(2)
C
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) CSI(2)
C
          IF (SEDTYPE.EQ.'SAND') CSI(1)=CSI(2)
        ENDIF
C
C
C  CONVERT INITIAL CONCENTRATIONS FROM mg/l TO g/cm**3
C
        DO 5022 K=1,KSED
          CSI(K)=CSI(K)/1000000.
 5022   CONTINUE
C
C  INITIALIZE COUNTERS
C
        NSEDCT=0           !Check for 0 ini. 
        NBLOW=0            !Check for 0 ini.
        N24CNT=0           !Check for 0 ini.
C
        NHR=NINT(3600./DTI)
        N24HR=24*NHR
      ENDIF
C
C**********************************************************************
C
C  CHEM TRANSPORT INPUT
C
      IF (CHEMTRAN.EQ.'INCLUDE') THEN
C
C  CHEMI(1) = INITIAL WATER COL. CHEM CONC., CLASS 1 (COHESIVE) (ug/l)
C  CHEMI(2) = INITIAL WATER COL. CHEM CONC., CLASS 2 (NON-COHESIVE) (ug/l)
C  NCHEMLAY = NO. OF LAYERS IN CHEM BED MODEL (EXCLUDING SURFACE LAYER)
C  CHEMTHIK = THICKNESS OF LAYERS IN CHEM BED MODEL (cm)
C  CHEMACT  = ACTIVE LAYER THICKNESS (cm)
C
C  CHEMDRAT1 =  PARTICLE-BOUND CHEMICAL DECAY RATE (1/day), CLASS 1
C  CHEMDRAT2 =  PARTICLE-BOUND CHEMICAL DECAY RATE (1/day), CLASS 2
C
C
C  COHESIVE PARTICLE-BOUND TRACER
C
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          READ (IUT404,11) (COM(I),I = 1,80)
          READ (IUT404,5030)CHEMI(1),NCHEMLAY,CHEMTHIK,CHEMACT,CHEMDRAT1
 5030     FORMAT (F10.0,I10,4F10.0)
        ENDIF
C
C  NON-COHESIVE PARTICLE-BOUND TRACER
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          READ (IUT405,11) (COM(I),I = 1,80)
          READ (IUT405,5030)CHEMI(2),NCHEMLAY,CHEMTHIK,CHEMACT,CHEMDRAT2
C
          IF (SEDTYPE.EQ.'SAND') CHEMI(1)=CHEMI(2)
        ENDIF
C
C  CONVERT DECAY RATE FROM 1/day TO 1/s
C
        CHEMDRAT1=CHEMDRAT1/86400.
        CHEMDRAT2=CHEMDRAT2/86400.
C
C  CALC. NUMBER OF LAYERS IN ACTIVE LAYER
C
        NACTLAY=NINT(CHEMACT/CHEMTHIK)
C
C  INPUT INITIAL CHEM BED CONCENTRATIONS (ug CHEM/g SOLIDS)
C
        If (RESTAR.EQ.'COLD START') Then
          OPEN (UNIT=85,FILE='bed_chemic',FORM='FORMATTED')
          DO 5029 N=1,NCHEMLAY
            DO 5029 I=1,IM
              READ (85,5028) (CBEDCHEM(N,I,J),J=1,JM)
 5028         FORMAT (10F8.0)
 5029     CONTINUE
          CLOSE (85)
        ENDIF
C
C  CONVERT FROM ug/l TO ug/cm**3
C
        CHEMI(1)=CHEMI(1)/1000.
        CHEMI(2)=CHEMI(2)/1000.
C
      ENDIF
C
C***********************************************************
C
C  FOR PARTICLE TRACKING     
C
C  NFREQ = FREQUENCY OF PARTICLE RELEASE (timesteps)
C  NPART = NUMBER OF PARTICLES PER RELEASE
C  NCONV = TOTAL NUMBER OF RELEASES BEFORE CONVERSION OF PARTICLES
C          TO CONCENTRATION (timesteps)
C
C  NOTE:  TOTAL NO. OF PARTICLES IN SYSTEM = NPART * NCONV
C
C  IRELST = BEGINNING TIMESTEP OF PARTICLE RELEASE
C  NPCLASS = NUMBER OF PARTICLE CLASSES
C
C  NSOURCE = NO. OF SOURCES OF PARTICLES
C
C  isource = i of particle input
C  jsource = j of particle input
C  ksource = k-level of particle input, for testing only
C
      IF (PARTICLE.EQ.'INCLUDE') THEN
C
        Read (IUT406,11) (COM(I),I = 1,80)
        WRITE (IUPRT,11) (COM(I),I = 1,80)
        READ (IUT406,5810)NFREQ,NPART,IRELST,IRELEND,
     +                          NSOURCE
        NCONV=(IRELEND-IRELST)/NFREQ+1 
        WRITE (IUPRT,5810)NFREQ,NPART,NCONV,IRELST,IRELEND,
     +                          NSOURCE
 5810   FORMAT (9I8)
C
        DO 5809 MM=1,NSOURCE
          READ (IUT406,5901)ISOURCE(MM),JSOURCE(MM),KSOURCE(MM)
          WRITE (IUPRT,5901)ISOURCE(MM),JSOURCE(MM),KSOURCE(MM)
 5809   CONTINUE
 5901   FORMAT (9I8)
C
        IF (NPART.GT.NPARTM) THEN
          WRITE (IUPRT,5811)NPART,NPARTM
 5811     FORMAT (/5X,'EXECUTION STOPPED BECAUSE NPART > NPARTM,'
     +     /6X,'NPART =',I6,3X,'NPARTM =',I6,
     +     /6X,'PLEASE SET NPART <= NPARTM AND RESUBMIT')
          STOP
        ENDIF
C
        IF (NCONV.GT.NCONVM) THEN
          WRITE (IUPRT,5812)NCONV,NCONVM
 5812     FORMAT (/5X,'EXECUTION STOPPED BECAUSE NCONV > NCONVM,'
     +     /6X,'NCONV =',I6,3X,'NCONVM =',I6,
     +     /6X,'PLEASE SET NCONV <= NCONVM AND RESUBMIT')
          STOP
        ENDIF
C
C
C  INTCONV = TIMESTEP AT WHICH TO INITIALLY START CONVERTING PARTICLES
C            TO CONSERVATIVE TRACER
C
        INTCONV=NCONV*NFREQ+IRELST
C
C  OPEN OUTPUT FILE FOR PARTICLES (TESTING ONLY)
C
        OPEN (UNIT=32,FILE='part_location',FORM='FORMATTED')
C
      ENDIF
C
C********************************************************************
C
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,7)  JHM,IAVGE
      WRITE(IUPRT,7) JHM,IAVGE
      IF(JHM.EQ.0) GO TO 200
      IF(JHM.GT.IHISTM) THEN
       WRITE(IUPRT,60) JHM,IHISTM
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  60  FORMAT(//' JHM    (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' IHISTM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      DEI=1./FLOAT(IAVGE)
      READ(IURUN,6)  (IHIST(I,2),I=1,JHM)
      WRITE(IUPRT,6) (IHIST(I,2),I=1,JHM)
 7    FORMAT(8I10)
 6    FORMAT(10I8)
C
 200  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,7)  ISKILL
      WRITE(IUPRT,7) ISKILL
      IF(ISKILL.EQ.0) THEN
       SKILLI=1.0
       GO TO 203
      ELSE
       SKILLI=1./FLOAT(ISKILL)
      ENDIF
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  EPTS
      WRITE(IUPRT,5) EPTS
 5    FORMAT(16I5)
      IF(EPTS.EQ.0) GO TO 201
      IF(EPTS.GT.EPTSM) THEN
       WRITE(IUPRT,61) EPTS,EPTSM
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  61  FORMAT(//' EPTS  (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' EPTSM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      READ(IURUN,5)  (INXIE(I),INXJE(I),I=1,EPTS)
      WRITE(IUPRT,5) (INXIE(I),INXJE(I),I=1,EPTS)
201   READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  VPTS
      WRITE(IUPRT,5) VPTS
      IF(VPTS.EQ.0) GO TO 202
      IF(VPTS.GT.VPTSM) THEN
       WRITE(IUPRT,62) VPTS,VPTSM
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  62  FORMAT(//' VPTS  (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' VPTSM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      READ(IURUN,5)  (INXIV(I),INXJV(I),I=1,VPTS)
      WRITE(IUPRT,5) (INXIV(I),INXJV(I),I=1,VPTS)
 202  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  FPTS
      WRITE(IUPRT,5) FPTS
      IF(FPTS.EQ.0) GO TO 203
      IF(FPTS.GT.FPTSM) THEN
       WRITE(IUPRT,63) VPTS,VPTSM
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  63  FORMAT(//' FPTS  (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' FPTSM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      READ(IURUN,49)  (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
 
      WRITE(IUPRT,49) (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
 49   FORMAT(4(2I5,1X,A4,I5))
C
 203  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,7)  JTM,NPLPF,ITRNFORM,IZERO,IWET
      WRITE(IUPRT,7) JTM,NPLPF,ITRNFORM,IZERO,IWET
      IF(IWET.EQ.0)THEN
       WRITE(IUPRT,*)'WQ INFO FOR ENTIRE GRID OPTION SELECTED'
      ELSE
       WRITE(IUPRT,*)'WQ INFO FOR WET-GRID OPTION SELECTED'
      END IF
      IF(ITRNFORM.EQ.0)
     . WRITE(IUPRT,*)'USER SPECIFIED TIME BREAKS FOR DUMP USED'
      IF(ITRNFORM.EQ.1)
     . WRITE(IUPRT,*)'ECOM WILL GENERATE THE TIME BREAKS VIA JTM,NPLPF,
     . AND IWET'
      IF(JTM.EQ.0) GO TO 204
      IF(JTM.GT.ITRACM) THEN
       WRITE(IUPRT,64) JTM,ITRACM
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  64  FORMAT(//' JTM    (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' ITRACM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      FLTWT=1./FLOAT(NPLPF)
      IF(ITRNFORM.EQ.0)THEN
      READ(IURUN,6) (ITRAC(I,2),I=1,JTM)
      ELSE
       DO I=1,JTM
        ITRAC(I,2)=IZERO+I*NPLPF
       END DO
      END IF
      WRITE(IUPRT,6) (ITRAC(I,2),I=1,JTM)
 204  CONTINUE
C
      IF(TOR.NE.'BAROTROPIC' .AND. TOR.NE.'PROGNOSTIC' .AND.
     .   TOR.NE.'TEMP_ONLY '.AND.TOR.NE.'SALT_ONLY '.AND.
     .   TOR.NE.'DIAGNOSTIC') THEN
      WRITE(IUPRT,24) TOR
 24   FORMAT(//' TYPE OF RUN (TOR=',A10,') IS SPECIFIED INCORRECTLY',/
     .         ' PLEASE FIX AND RESUBMIT'//)
      !!!&&&CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
C
      IF (ADVECT.NE.'LINEAR    ' .AND. ADVECT.NE.'NON-LINEAR') THEN
      WRITE(IUPRT,25) ADVECT
  25  FORMAT(//' TYPE OF MOMENTUN ADVECTION (ADVECT=',A10,') IS',
     .'SPECIFIED INCORRECLTY',/' PLEASE FIX AND RESUBMIT'//)
      !!!&&&CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
C
      IF(SCHEME.NE.'CENTRAL   ' .AND. SCHEME.NE.'UPWIND    ' .AND.
     .   SCHEME.NE.'SMOLAR_R  ' .AND. SCHEME.NE.'SMOLAR_2  ') THEN
      WRITE(IUPRT,34) SCHEME
 34   FORMAT(//' T&S ADVECTION SCHEME (SCHEME=',A10,') IS SPECIFIED',
     .   ' INCORRECTLY',/' PLEASE FIX AND RESUBMIT'//)
      !!!&&&CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
C
      IF (HORZMIX.NE.'CLOSURE   ' .AND. HORZMIX.NE.'CONSTANT  ') THEN
      WRITE(IUPRT,26) HORZMIX
  26  FORMAT(//' TYPE OF HORIZONTAL MIXING (HORZMIX=',A10,') IS',
     .         ' SPECIFIED INCORRECTLY',/
     .         ' PLEASE FIX AND RESUBMIT'//)
      !!!&&&CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
C
      IF (VERTMIX.NE.'CLOSURE   ' .AND. VERTMIX.NE.'CONSTANT  ') THEN
      WRITE(IUPRT,27) VERTMIX
  27  FORMAT(//' TYPE OF VERTICAL MIXING (VERTMIX=',A10,') IS',
     .         ' SPECIFIED INCORRECTLY',/
     .         ' PLEASE FIX AND RESUBMIT'//)
      !!!&&&CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
C
C
      IF(TOR.EQ.'BAROTROPIC') THEN
      WRITE(IUPRT,14) TOR
      ELSE
      WRITE(IUPRT,13) TOR
      END IF
 13   FORMAT(/' THIS IS A THREE DIMENSIONAL MODEL RUN',2X,A10/)
 14   FORMAT(/' THIS IS A TWO DIMENSIONAL MODEL RUN',2X,A10/)
C
      WRITE(IUPRT,22) ADVECT
 22   FORMAT(/' THIS SIMULATION HAS ',A10,' MOMENTUN ADVECTION '/)
C
      WRITE(IUPRT,222) SCHEME
 222  FORMAT(/' THIS SIMULATION USES ',A10,' DIFFERENCING FOR T&S'/)
C
      IF(HORZMIX.EQ.'CLOSURE   ') THEN
      WRITE(IUPRT,29) HORZMIX,HORCON,HPRNU
      ELSE
      WRITE(IUPRT,31) HORZMIX,HORCON,HPRNU
      END IF
 29   FORMAT(/' THIS SIMULATION HAS ',A10,' HORIZONTAL MIXING ',
     . ' HORCON = ',1PE10.3,' HPRNU = ',1PE10.3/)
 31   FORMAT(/' THIS SIMULATION HAS ',A10,' HORIZONTAL MIXING ',
     . ' CONSTANT = ',1PE10.3,'m**2/s  HPRNU = ',1PE10.3/)
C
      IF(VERTMIX.EQ.'CLOSURE   ') THEN
      WRITE(IUPRT,32) VERTMIX,UMOL,VPRNU
      ELSE
      WRITE(IUPRT,33) VERTMIX,UMOL,VPRNU
      END IF
 32   FORMAT(/' THIS SIMULATION HAS ',A10,' VERTICAL MIXING ',
     . ' UMOL = ',1PE10.3,' VPRNU = ',1PE10.3/)
 33   FORMAT(/' THIS SIMULATION HAS ',A10,' VERTICAL MIXING ',
     . ' CONSTANT = ',1PE10.3,'m**2/s  VPRNU = ',1PE10.3/)
C
      CALL ZEROES(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
C
      IF(RESTAR.EQ.'COLD START') THEN
       CALL SETDOM(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
c
       ISTART=INT+1
       IEND=INT+NSTEPS
c
       CALL BCDATA
       CALL TANDS
       CALL DENS
      ELSE
       OPEN (IURRS,FORM='unformatted',FILE='restart')
       IF(HYDTYPE.EQ.'EXTERNAL'.AND.IWET.EQ.1)
     . CALL SETDOM(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
C
C
C  NO HOT START CAPABILITY FOR PARTICLE TRACKING 
C
C
       IF (TRACER.EQ.'INCLUDE') THEN
        IF (SEDTRAN.EQ.'NEGLECT') THEN
         READ (IURRS) 
     .      INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +      CONC1,CONC1B
        ELSE
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
            READ (IURRS) 
     .      INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +      CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +        SEDMASS,SEDEP
          ELSE
            READ (IURRS) 
     .      INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +      CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,TAU
          ENDIF
        ENDIF
       ELSE
        IF (SEDTRAN.EQ.'NEGLECT') THEN
         READ (IURRS) 
     .      INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN
        ELSE
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
            READ (IURRS) 
     .      INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +      N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +        SEDMASS,SEDEP
          ELSE
            READ (IURRS) 
     .      INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +      N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,TAU
c          write(*,*)'HOT READ, INT,T,S,U,V=',
c     .T(20,10,1),S(20,10,1),U(20,10,1),V(20,10,1)
          ENDIF
        ENDIF
       ENDIF
C
C********************************************************************
C
C  FOR PARTICLE TRACKING (corner_loc INPUT IN setdom FOR COLD START)
C
      IF (PARTICLE.EQ.'INCLUDE') THEN
        IRELST=IRELST+INT
C
C  OPEN CORNER LOCATIONS FILE
C
C  CORNER LOCATION CONVENTION:  XCOR(I,J) = x(i-1/2, j-1/2) 
C     (LOWER LEFT-HAND CORNER)  YCOR(I,J) = y(i-1/2, j-1/2)
C
      OPEN (UNIT=33,FILE='corner_loc',FORM='FORMATTED',STATUS='OLD')
C
      DO 2300 N=1,1000000
          READ (33,*,END=2301)I,J,XCOR(I,J),YCOR(I,J)
 2300 CONTINUE
 2301 REWIND(33)
C
C  CALC. H1 AND H2 AT ELEMENT INTERFACES
C
 2310   DO 2320 J=1,JM
          DO 2320 I=1,IM
            H1P(I,J)=SQRT((XCOR(I+1,J)-XCOR(I,J))**2.+
     +                        (YCOR(I+1,J)-YCOR(I,J))**2.)
C
            H2P(I,J)=SQRT((XCOR(I,J+1)-XCOR(I,J))**2.+
     +                        (YCOR(I,J+1)-YCOR(I,J))**2.)
 2320   CONTINUE
      ENDIF
C
C*********************************************************************
C
       CLOSE (IURRS)
       ISTART=INT+1
       IEND=INT+NSTEPS
       EDAY=SDAY+TNDAYS+FLOAT(INT)*DAYI*DTI
C RESET IYR,IMO,IDA FOR PTIDE IF HOT START 
       SD_HOTSTART=SDAY+FLOAT(INT)*DAYI*DTI
       HOUR=AMOD(SD_HOTSTART,1.0)*24.0  
       IHOUR1=IFIX(HOUR)
       KD=IFIX(SDAY+FLOAT(INT)*DAYI*DTI)
       CALL CDAY(IDA1,IMO1,IYR1,KD,2)
       TMP_IYR=IYR
       TMP_IMO=IMO
       TMP_IDA=IDA
       TMP_IHOUR=IHOUR
c
       IYR=IYR1
       IMO=IMO1
       IDA=IDA1
       IHOUR=IHOUR1
       WRITE(IUPRT,211) IYR,IMO,IDA,IHOUR
       CALL BCDATA
       IYR=TMP_IYR
       IMO=TMP_IMO
       IDA=TMP_IDA
       IHOUR=TMP_IHOUR
      ENDIF
C
C  INITIALIZE SEDIMENT TRANSPORT VARIABLES 
C
        IF (SEDTRAN.EQ.'INCLUDE') CALL SEDIC(RESTAR)
C
C**********************************************************************
C
      DO 6030 J = 1, JM
        DO 6020 I = 1, IM
          If (H(I,J).GT.0.0) Then
C
C  VARIABLE BOTTOM FRICTION FOR SEDIMENT TRANSPORT
C
            IF (SEDTRAN.EQ.'INCLUDE') THEN
C
C  COHESIVE ELEMENTS
C
              IF (IBMSK(I,J).EQ.0) THEN
                Z0 = Z0BCOH
                CBCMIN = BFCOH
                IF (TOR.EQ.'BAROTROPIC') THEN
                  CBC(I,J)=BFCOH*FSM(I,J)
                ELSE
                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
     +                      H(I,J)/Z0)**2)*FSM(I,J)
                ENDIF
C
C  NON-COHESIVE ELEMENTS
C
              ELSE
                Z0 = Z0B
                CBCMIN = BFRIC
                IF (TOR.EQ.'BAROTROPIC') THEN
                  CBC(I,J)=BFRIC*FSM(I,J)
                ELSE
                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
     +                      H(I,J)/Z0)**2)*FSM(I,J)
                 ENDIF
               ENDIF
            ELSE
              Z0 = Z0B
              CBCMIN = BFRIC
              IF (TOR.EQ.'BAROTROPIC') THEN
                CBC(I,J)=BFRIC*FSM(I,J)
              ELSE
                CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
     +                      H(I,J)/Z0)**2)*FSM(I,J)
              ENDIF
            ENDIF
          END IF
 6020   CONTINUE
 6030 CONTINUE
C
C  SET INITIAL TRACER CONCENTRATIONS
C
      IF (TRACER.EQ.'INCLUDE'.AND.RESTAR.EQ.'COLD START') THEN
        DO 6034 K=1,KBM1	
          DO 6034 J=2,JMM1
            DO 6034 I=2,IMM1
              IF (H(I,J).GT.0.0) THEN
                CONC1(I,J,K)=CONINIT
                CONC1B(I,J,K)=CONINIT
              ENDIF
 6034   CONTINUE
      ENDIF
C
C  SET INITIAL SEDIMENT CONCENTRATIONS
C
      IF (SEDTRAN.EQ.'INCLUDE'.AND.RESTAR.EQ.'COLD START') THEN
        DO 6035 K=1,KBM1	
          DO 6035 J=2,JMM1
            DO 6035 I=2,IMM1
              IF (H(I,J).GT.0.0) THEN
                IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
                  CSED1(I,J,K)=CSI(1)
                  CSED1B(I,J,K)=CSI(1)
                ENDIF
C
                IF (SEDTYPE.EQ.'SAND') THEN
                  CSED2(I,J,K)=CSI(1)
                  CSED2B(I,J,K)=CSI(1)
                ELSE
                  CSED2(I,J,K)=CSI(2)
                  CSED2B(I,J,K)=CSI(2)
                ENDIF
              ENDIF
 6035   CONTINUE
      ENDIF
C
C***********************************************************
C
C  CHEM TRANSPORT 
C
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.RESTAR.EQ.'COLD START') THEN
C
C  SET:  1. INITIAL CHEM WATER COLUMN CONCENTRATIONS (ug/cm**3)
C        2. INITIAL CHEM MASS IN BED LAYERS (ug CHEM/cm**2)
C        3. INITIAL SEDIMENT MASS IN BED LAYERS (g solids/cm**2)
C
        DO 36 K=1,KBM1	
          DO 36 J=2,JMM1
            DO 36 I=2,IMM1
              IF (H(I,J).GT.0.0) THEN
                IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
                  CHEM1(I,J,K)=CHEMI(1)
                  CHEM1B(I,J,K)=CHEMI(1)
                ENDIF
C
                IF (SEDTYPE.EQ.'SAND') THEN
                  CHEM2(I,J,K)=CHEMI(1)
                  CHEM2B(I,J,K)=CHEMI(1)
                ELSE
                  CHEM2(I,J,K)=CHEMI(2)
                  CHEM2B(I,J,K)=CHEMI(2)
                ENDIF
              ENDIF
 36     CONTINUE
C
        DO 37 N=1,NCHEMLAY
          DO 37 J=1,JM
            DO 37 I=1,IM
              IF (FSM(I,J).GT.0.0) THEN
                SEDMASS(N,I,J)=CBED(I,J)*CHEMTHIK
C
C  TOP LAYER IS INITIALLY VERY THIN
C
                IF (N.EQ.1) SEDMASS(N,I,J)=0.0
                CHEMMASS(N,I,J)=SEDMASS(N,I,J)*CBEDCHEM(N,I,J)
              ELSE
                CBEDCHEM(N,I,J)=0.0
                SEDMASS(N,I,J)=0.0
                CHEMMASS(N,I,J)=0.0
              ENDIF
 37     CONTINUE
      ENDIF
C
C***********************************************************
C
C  CALC. FETCH & MEAN DEPTH FOR WIND WAVE MODEL
C
      IF (WAVEDYN.EQ.'SMBMODEL') CALL FHCALC
C
C  INPUT INITIAL FLOW FIELD FROM gcm_tran
C
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
        OPEN (IUTRN,FORM='unformatted',FILE='gcm_tran')
        NHYDCNT=0
        CALL TRANINP(0)
      ENDIF
C
C*******************************************************************
C
C
      IF(RESTAR.EQ.'COLD START') CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
C
      ISTART=INT+1
      IEND=INT+NSTEPS
      TPRT=FLOAT(IPRINT)*DAYI*DTI
      TAVG=FLOAT(IAVGE)*DAYI*DTI
      TSKILL=FLOAT(ISKILL)*DAYI*DTI
      TRACE=FLOAT(NPLPF)*DAYI*DTI
      DO 16 I=1,JHM
      IHIST(I,2)=INT+IHIST(I,2)
 16   IHIST(I,1)=IHIST(I,2)-IAVGE+1
      DO 18 I=1,JTM
      ITRAC(I,2)=INT+ITRAC(I,2)
 18   ITRAC(I,1)=ITRAC(I,2)-NPLPF+1
C 
      WRITE(IUPRT,15) ISTART,IEND 
 15   FORMAT(//30H MODEL STARTING UP...ISTART = ,I6,8H IEND = ,I10/)
      KD=IFIX(EDAY)
      WRITE(IUPRT,*)'SDAY AND EDAY',SDAY,EDAY
      CALL CDAY(IDAE,IMOE,IYRE,KD,2)
      WRITE(IUPRT,210)IYR,IMO,IDA,IYRE,IMOE,IDAE
 210  FORMAT('SIMULATION FROM',I5,2I3,4X,'TO',3X,I5,2I3) 
      WRITE(IUPRT,20) TNDAYS
 20   FORMAT(//32H NUMBER OF DAYS IN SIMULATION = ,F6.2/) 
      WRITE(IUPRT,30) TPRT,IPRINT,IPRTSTART,TAVG,IAVGE,TSKILL,ISKILL
 30   FORMAT(//' TPRT =   ',F10.3,'  IPRINT =  ',I10,
     .        '  IPRTSTART =  ',I10,//
     .         ' TAVG =   ',F10.3,'  IAVGE =   ',I10,//
     .         ' TSKILL = ',F10.3,'  ISKILL =  ',I10,//)
      WRITE(IUPRT,35)
      WRITE(IUPRT,40) (IHIST(I,1),IHIST(I,2),I=1,JHM)
 35   FORMAT(//' HISTORY TAKEN AT TIMESTEPS   START ----- STOP ')
 40   FORMAT(27X,I8,3X,I8)
      WRITE(IUPRT,47)
      WRITE(IUPRT,48) (ITRAC(I,1),ITRAC(I,2),I=1,JTM)
47    FORMAT(//' QUALITY PARAMETERS FOR RCA  INTEGRATED OVER TIMESTEPS
     .      START      STOP')
48    FORMAT(57X,I8,2X,I8)
C
      WRITE(IUPRT,41) BFRIC,Z0B,NU,THETA,ALPHA,TLAG,NWAVE,BCTYPE
 41   FORMAT(//' BFRIC            =   ',F10.4,' nondimensional'/
     .         ' Z0B              =   ',F10.4,' m'/
     .         ' NU               =   ',F10.4,' nondimensional'/
     .         ' THETA            =   ',F10.4,' nondimensional'/
     .         ' ALPHA            =   ',F10.4,' nondimensional'/
     .         ' TLAG             =   ',F10.4,
     .         ' FRICTION TIME SCALE in PCLAMP BC (Hours)'/
     .         ' WAVE MODEL ACTIVATED AFTER  =   ',I10,  ' Time Steps'/
     .         ' BOUNDARY TYPE  =   ',3x,A7/)
C
      TLAG=TLAG*3600+1.0e-10    ! CONVERT TO SECONDS
      TLAG=1./TLAG
      ALPHA=ALPHA*3600          ! CONVERT TO SECONDS
      IF(ALPHA.GT.0)ALPHA=1./ALPHA
c
      WRITE(IUPRT,70)
  70  FORMAT(/1X,' K',6X,'Z',10X,'ZZ',8X,'DZ',/)
      DO 90 K=1,KBM1
      WRITE(IUPRT,80) K,Z(K),ZZ(K),DZ(K)
  90  CONTINUE
      WRITE(IUPRT,80) K,Z(KB)
  80  FORMAT(I3,3F10.3)
C
      THOUR=FLOAT(INT)*DTI/3600.
      CALL FIRST
C 
      AREA=0.0
      EMI=0.0
      APEI=0.0
      DO 280 J=1,JM
      DO 280 I=1,IM
  280 AREA=AREA+FSM(I,J)*ART(I,J)
      DO 285 K=1,KBM1
      DO 285 J=1,JM
      DO 285 I=1,IM
      TRHO=(RMEAN(I,J,K)+1.)*1000.
      VOL=H(I,J)*ART(I,J)*DZ(K)
      EMI=EMI+TRHO*VOL*FSM(I,J)
  285 APEI=APEI+GRAV*TRHO*ZZ(K)*H(I,J)*VOL*FSM(I,J)
      WRITE(IUPRT,85) AREA,EMI,APEI
   85 FORMAT(//
     . ' SURFACE AREA                          ',1PE14.7,' m**2  ',//,
     . ' INITIAL MASS                          ',1PE14.7,' Kg    ',//,
     . ' INITIAL AVAILABLE POTENTIAL ENERGY    ',1PE14.7,' joules',//)
*
*      INITIALIZE TBDRY2 AND SBDRY2  
*******   NEW VARIABLES FOR USE IN antidif.f 
*
         DO N = 1, NUMEBC
           IE = IETA(N)
           JE = JETA(N)
           DO K = 1, KBM1
C             TBDRY2(N,K) = TB(IE,JE,K)
C             SBDRY2(N,K) = SB(IE,JE,K)
             TBDRY2(N,K) = T(IE,JE,K)
             SBDRY2(N,K) = S(IE,JE,K)
           ENDDO
         ENDDO 
*
*

C
C-----------------------------------------------------------------------
C
C                 BEGIN NUMERICAL INTEGRATION
C
C-----------------------------------------------------------------------
C
                   DO 9000 INT=ISTART,IEND
C 
C     RAMP=TANH(FLOAT(INT)/FLOAT(IRAMP+1))
          RAMP = FLOAT(INT)/FLOAT(IRAMP)
          IF(INT.GT.IRAMP)RAMP=1.0

      TIME=FLOAT(INT)*DAYI*DTI
      THOUR=TIME*24.
C
C  WAVE EFFECTS ON BOTTOM FRICTION (HYDRODYNAMICS ONLY)
C
      CALL BCOND(7,DTI2,0)
C
      IF (WAVEDYN.EQ.'SMBMODEL'.OR.
     .WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'EXTERNAL') THEN
        NWAVECNT=NWAVECNT+1
        IF (NWAVECNT.EQ.NWAVE) THEN
          NWAVECNT=0
          CALL BOTTAU   ! compute shear stress using currents+wave
        ENDIF
      ELSE
          CALL BOTTAU   ! compute shear stress using currents only
      ENDIF
C
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
C
        NHYDCNT=NHYDCNT+1
        IF(INT.EQ.1) GO TO 8200 
C  
C  INPUT NEW FLOW FIELD
C
        IF (NHYDCNT.GT.NHYD) THEN
          CALL TRANINP(1)
          NHYDCNT=1
        ENDIF
C
C
C  SKIP HYDRODYNAMICS
C
        GOTO 9900
      ENDIF
C
      CALL BCOND(7,DTI2,0)
C
      IF(TOR.NE.'BAROTROPIC') THEN
C
      CALL BAROPG(DRHOX,DRHOY,TRNU,TRNV)
C
      DO 50 J=1,JM
      DO 50 I=1,IM
      TRNU(I,J)=TRNU(I,J)+ADVUU(I,J)-ADVUA(I,J)
 50   TRNV(I,J)=TRNV(I,J)+ADVVV(I,J)-ADVVA(I,J)
C
      DO 120 J=1,JM
      DO 120 I=1,IM
 120  EGF(I,J)=EL(I,J)*ISPI
      DO 400 J=2,JM
      DO 400 I=2,IM
      UTF(I,J)=UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
      VTF(I,J)=VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I
400   CONTINUE
C
      ENDIF
C
      IF(HORZMIX.EQ.'CLOSURE   ') CALL SMAG
C
C-----------------------------------------------------------------------
C         BEGIN EXTERNAL MODE
C-----------------------------------------------------------------------
C
                   DO 8000 IEXT=1,ISPLIT
C
      CALL EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2)
C
      IF(TOR.EQ.'BAROTROPIC') GO TO 440
      IF(IEXT.LT.(ISPLIT-2)) GO TO 440
      IF(IEXT.EQ.(ISPLIT-2)) THEN
      DO 402 J=1,JM
      DO 402 I=1,IM
 402  ETF(I,J)=.25*NU*ELF(I,J)
      GO TO 440
      ENDIF
      IF(IEXT.EQ.(ISPLIT-1)) THEN
      DO 404 J=1,JM
      DO 404 I=1,IM
 404  ETF(I,J)=ETF(I,J)+.5*(1.-.5*NU)*ELF(I,J)
      GO TO 440
      ENDIF
      IF(IEXT.EQ.(ISPLIT-0)) THEN
      DO 406 J=1,JM
      DO 406 I=1,IM
 406  ETF(I,J)=(ETF(I,J)+.5*ELF(I,J))*FSM(I,J)
      ENDIF
 440  CONTINUE
C
C  TEST FOR MODEL BLOWUP. IF SO, PRINT AND STOP
C
      VAMAX=-1.E10
      DO 442 J=1,JM
      DO 442 I=1,IM
      VMAXTMP=SQRT(UAF(I,J)**2+VAF(I,J)**2)
      IF(VMAXTMP.GE.VAMAX) THEN
        IAMAX=I
        JAMAX=J
        VAMAX=VAF(I,J)
      END IF
  442 CONTINUE
      IF(VAMAX.GT.10.) GO TO 9100
cstart, check elevation minimum
C elevation
      DO 443 J=1,JM
      DO 443 I=1,IM
      elminx=1.0+elf(i,j)*fsm(i,j)/(h(i,j)+1.0e-5)
      IF(elminx.lt.0.1) THEN ! 10 percent !!
        IAMAX=I
        JAMAX=J
        elmin=elf(I,J)
        go to 9100
      END IF
      IF(elminx.lt.0.20) THEN ! 20 percent !! Warning
        IAMAX=I
        JAMAX=J
        elmin=elf(I,J)
        WRITE(IUPRT,7501)IAMAX,JAMAX,elmin,TIME  
      END IF

  443 CONTINUE
C-----------------------------------------------------------------------
C         APPLY FILTER TO REMOVE TIME SPLIT AND RESET TIME SEQUENCE
C-----------------------------------------------------------------------
C
      DO 150 J=1,JM 
      DO 150 I=1,IM 
      UA(I,J)=UA(I,J)+.5*NU*(UAB(I,J)-2.*UA(I,J)+UAF(I,J)) 
      VA(I,J)=VA(I,J)+.5*NU*(VAB(I,J)-2.*VA(I,J)+VAF(I,J)) 
 150  EL(I,J)=EL(I,J)+.5*NU*(ELB(I,J)-2.*EL(I,J)+ELF(I,J)) 
      IF(TOR.EQ.'BAROTROPIC') THEN
      DO 155 JTRAC=1,JTM
      IF(INT.GE.ITRAC(JTRAC,1).AND.INT.LE.ITRAC(JTRAC,2)) CALL TRANSPORT 
 155  CONTINUE
      END IF
      DO 160 J=1,JM 
      DO 160 I=1,IM 
      ELB(I,J)=EL(I,J)
      EL(I,J)=ELF(I,J)
      D(I,J)=H(I,J)+EL(I,J) 
      UAB(I,J)=UA(I,J)
      UA(I,J)=UAF(I,J)
      VAB(I,J)=VA(I,J)
 160  VA(I,J)=VAF(I,J)
C
      IF(TOR.EQ.'BAROTROPIC') GO TO 8000
      IF(IEXT.EQ.ISPLIT) GO TO 8000
      DO 445 J=1,JM
      DO 445 I=1,IM
 445  EGF(I,J)=EGF(I,J)+EL(I,J)*ISPI
      DO 450 J=2,JM
      DO 450 I=2,IM
      UTF(I,J)=UTF(I,J)+UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
 450  VTF(I,J)=VTF(I,J)+VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I
C
 8000              CONTINUE
C
C-----------------------------------------------------------------------
C         END EXTERNAL (2-D) MODE CALCULATION
C      AND CONTINUE WITH INTERNAL (3-D) MODE CALCULATION
C-----------------------------------------------------------------------
C 
      IF(INT.EQ.1) GO TO 8200    
C
      IF (TOR.EQ.'BAROTROPIC') GO TO 9900
C
C-----------------------------------------------------------------------
C         ADJUST U(Z) AND V(Z) SUCH THAT
C      VERTICAL AVERAGE OF (U,V) = (UA,VA)
C-----------------------------------------------------------------------
C
      DO 299 J=1,JM
      DO 299 I=1,IM
 299  TPS(I,J)=0.0
      DO 300 K=1,KBM1
      DO 300 J=2,JM
      DO 300 I=2,IM
 300  TPS(I,J)=TPS(I,J)+U(I,J,K)*DZ(K)
      DO 302 K=1,KBM1
      DO 302 J=2,JM
      DO 302 I=2,IM
 302  U(I,J,K)=(U(I,J,K)-TPS(I,J))+
     .     (UTB(I,J)+UTF(I,J))/(DT(I,J)+DT(I-1,J))
      DO 303 J=1,JM
      DO 303 I=1,IM
 303  TPS(I,J)=0.0
      DO 304 K=1,KBM1
      DO 304 J=2,JM
      DO 304 I=2,IM
 304  TPS(I,J)=TPS(I,J)+V(I,J,K)*DZ(K)
      DO 306 K=1,KBM1
      DO 306 J=2,JM
      DO 306 I=2,IM
 306  V(I,J,K)=(V(I,J,K)-TPS(I,J))+
     .     (VTB(I,J)+VTF(I,J))/(DT(I,J)+DT(I,J-1))
C
C-------- CORRECT FOR RIVER BOUNDARY CONDITIONS ------------------------
      DO 320 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 320 K=1,KBM1
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            IF(VQDIST(N,K).EQ.0.0) U(IC,JC,K)=0.0
         ELSE
            IF(VQDIST(N,K).EQ.0.0) U(ID,JD,K)=0.0
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            IF(VQDIST(N,K).EQ.0.0) V(IC,JC,K)=0.0
         ELSE
            IF(VQDIST(N,K).EQ.0.0) V(ID,JD,K)=0.0
         ENDIF
      ENDIF
  320 CONTINUE
C
C-----------------------------------------------------------------------
C         CALCULATE HORIZONTAL MASS FLUXES, (H2*U*D) AND (H1*V*D)
C-----------------------------------------------------------------------
      DO 310 K=1,KBM1
      DO 311 J=2,JMM1
      DO 311 I=2,IM
      XMFL3D(I,J,K)=0.25*(H2(I-1,J)+H2(I,J))*(DT(I-1,J)+DT(I,J))*
     .   U(I,J,K)
  311 CONTINUE
      DO 312 J=2,JM
      DO 312 I=2,IMM1
      YMFL3D(I,J,K)=0.25*(H1(I,J-1)+H1(I,J))*(DT(I,J-1)+DT(I,J))*
     .   V(I,J,K)
  312 CONTINUE
  310 CONTINUE
C
C-----------------------------------------------------------------------
C         VERTVL INPUT = U,V,DT(=H+ET),ETF,ETB; OUTPUT = W
C----------------------------------------------------------------------- 
      CALL VERTVL(DTI2)
      CALL BCOND(5,DTI2,0)
      DO 307 JTRAC=1,JTM
      IF(INT.GE.ITRAC(JTRAC,1).AND.INT.LE.ITRAC(JTRAC,2)) CALL TRANSPORT 
 307  CONTINUE
C
C-----------------------------------------------------------------------
C         COMPUTE Q2F AND Q2LF USING UF AND VF AS TEMPORARY VARIABLES
C-----------------------------------------------------------------------
C
      IF(VERTMIX.EQ.'CLOSURE   ') THEN
      CALL ADVQ(Q2B,Q2,DTI2,UF)
      CALL ADVQ(Q2LB,Q2L,DTI2,VF)
      CALL PROFQ(DTI2)
      CALL BCOND(6,DTI2,0)
      DO 325 K=1,KB 
      DO 325 J=1,JM 
      DO 325 I=1,IM 
      Q2 (I,J,K)=Q2 (I,J,K)+.5*NU*(UF(I,J,K)+Q2B(I,J,K)-2.*Q2(I,J,K)) 
      Q2B(I,J,K)=Q2(I,J,K)
      Q2(I,J,K)=UF(I,J,K) 
 325  CONTINUE
      DO 335 K=1,KB 
      DO 335 J=1,JM 
      DO 335 I=1,IM 
      Q2L(I,J,K)=Q2L(I,J,K)+.5*NU*(VF(I,J,K)+Q2LB(I,J,K)-
     .  2.*Q2L(I,J,K))
      Q2LB(I,J,K)=Q2L(I,J,K)
      Q2L(I,J,K)=VF(I,J,K)
 335  CONTINUE
      END IF
C
C-----------------------------------------------------------------------
C     COMPUTE UPDATED BOTTOM FRICTION (CBC) EVERY TIMESTEP
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C         COMPUTE TF AND SF USING UF AND VF AS TEMPORARY VARIABLES
C-----------------------------------------------------------------------
C
*
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'TEMP_ONLY ')
     .   CALL ADVT(TB,T,TMEAN,DTI2,UF,TDIF,TDIS,TBDRY2)
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY ')
     .   CALL ADVT(SB,S,SMEAN,DTI2,VF,SDIF,SDIS,SBDRY2) 
*
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'TEMP_ONLY ')
     .   CALL PROFT(UF,WTSURF,DTI2,SWRAD)
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY ')
     .   CALL PROFT(VF,WSSURF,DTI2,ZEROS)
      CALL BCOND(4,DTI2,0)
*
         DO N = 1, NUMEBC
           IE = IETA(N)
           JE = JETA(N)
           DO K = 1, KBM1
             TBDRY2(N,K) = UF(IE,JE,K)
             SBDRY2(N,K) = VF(IE,JE,K)
           ENDDO
         ENDDO 
*
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'TEMP_ONLY ')THEN
      DO 345 K=1,KB 
      DO 345 J=1,JM 
      DO 345 I=1,IM 
      T(I,J,K)=T(I,J,K)+.5*NU*(UF(I,J,K)+TB(I,J,K)-2.*T(I,J,K))
      TB(I,J,K)=T(I,J,K)
      T(I,J,K)=UF(I,J,K)
 345  CONTINUE
      ENDIF
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY ')THEN
      DO 355 K=1,KB 
      DO 355 J=1,JM 
      DO 355 I=1,IM 
      S(I,J,K)=S(I,J,K)+.5*NU*(VF(I,J,K)+SB(I,J,K)-2.*S(I,J,K))
      SB(I,J,K)=S(I,J,K)
      S(I,J,K)=VF(I,J,K)
 355  CONTINUE

      ENDIF
C 
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY '.OR.
     .   TOR.EQ.'TEMP_ONLY ')THEN
         CALL DENS
      END IF
C
C  FOR CONSERVATIVE TRACER
C
 9900 IF (TRACER.EQ.'INCLUDE') THEN
        CALL ADVCON(CONC1B,CONC1,CMEAN1,DTI2,UF,CDIFF1,CDIS1,CBDRY1,
     +                               CPSTR,CONDRAT)
C
        CALL PROFT(UF,WCSURF,DTI2,ZEROS)
        CALL BCOND(8,DTI2,0)
        DO 346 K=1,KB
          DO 346 J=1,JM 
            DO 346 I=1,IM 
              CONC1(I,J,K)=CONC1(I,J,K)+.5*NU*(UF(I,J,K)+CONC1B(I,J,K)-
     +              2.*CONC1(I,J,K))
              CONC1B(I,J,K)=CONC1(I,J,K)
              CONC1(I,J,K)=UF(I,J,K)
 346    CONTINUE
      ENDIF
C
C********************************************************************
C
C  FOR PARTICLE TRACKING 
C
        IF (PARTICLE.EQ.'INCLUDE') THEN
          IF (INT.GE.IRELST) THEN
            IF (INT.EQ.IRELST) IDUM=-1
            CALL PARTRAK(IDUM,IRELEND)
          ENDIF
        ENDIF
C
C********************************************************************
C
C  SEDIMENT TRANSPORT
C
        IF (SEDTRAN.EQ.'INCLUDE'.AND.INT.GT.NSEDBEG) THEN
C
c         WHEN INT=1 HYDRODYNAMIC SKIPS ie THIS SEGMENT SKIPS    
          IF(INT.EQ.2)NSEDCT=1
          IF(INT.EQ.2)N24CNT=1
          NSEDCT=NSEDCT+1
          N24CNT=N24CNT+1

C  COHESIVE CLASS
C
          IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
            DO 6300 N=1,NUMQBCSE
              DUMB1(N)=CDIS(1,N)
 6300       CONTINUE
C
            DO 6310 N=1,NUMEBCSE
              DO 6310 K=1,KBM1
                DUMB2(N,K)=CBDRY(1,N,K)
 6310      CONTINUE
C
C  FOR DIFFUSER INPUTS
C
            DO 6315 N=1,NUMDBCSE
                DUMB3(N)=CSDIFF(1,N)
 6315       CONTINUE
C
            CALL ADVSED(CSED1B,CSED1,DTI2,UF,DUMB1,DUMB2,DUMB3)
C
            CALL PROFSED(UF,DTI2,WSET1,WCT1BOT) 
C
            CALL BCOND(9,DTI2,1)
C
            DO 6320 K=1,KBM1 
              DO 6320 J=1,JM
                DO 6320 I=1,IM
                  IF (UF(I,J,K).LE. 0.0 )THEN
                    WCT1BOT(I,J)=0.0
                    CSED1(I,J,K)=0.0
                  ELSE
                    CSED1(I,J,K)=CSED1(I,J,K)+.5*NU*(UF(I,J,K)+
     +                           CSED1B(I,J,K)-2.*CSED1(I,J,K))
                  ENDIF
                  CSED1B(I,J,K)=CSED1(I,J,K)
                  CSED1(I,J,K)=UF(I,J,K)
 6320       CONTINUE
          ENDIF
C
8900  Format (f8.3,I7,10F10.5)
9007  Format (2X,'        (02,12)          (03,11)
     .          (03,12)          (03,13)')


C
C  NON-COHESIVE CLASS
C
          IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
C
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
C
            DO 6330 N=1,NUMQBCSE
              DUMB1(N)=CDIS(KK,N)
 6330       CONTINUE
C
            DO 6340 N=1,NUMEBCSE
              DO 6340 K=1,KBM1
                DUMB2(N,K)=CBDRY(KK,N,K)
 6340       CONTINUE
C
C  FOR DIFFUSER INPUTS
C
            DO 6345 N=1,NUMDBCSE
                DUMB3(N)=CSDIFF(KK,N)
 6345       CONTINUE
C
            CALL ADVSED(CSED2B,CSED2,DTI2,UF,DUMB1,DUMB2,DUMB3)
C
C
            CALL PROFSED(UF,DTI2,WSET2,WCT2BOT)  
C
            CALL BCOND(9,DTI2,KK)
C
            DO 6350 K=1,KBM1
              DO 6350 J=1,JM
                DO 6350 I=1,IM
                  IF (UF(I,J,K).LE. 0.0 )THEN
                    WCT2BOT(I,J)=0.0
                    CSED2(I,J,K)=0.0
                  ELSE
                    CSED2(I,J,K)=CSED2(I,J,K)+.5*NU*(UF(I,J,K)+
     +                           CSED2B(I,J,K)-2.*CSED2(I,J,K))
                  ENDIF
C
                  CSED2B(I,J,K)=CSED2(I,J,K)
                  CSED2(I,J,K)=UF(I,J,K)

 6350       CONTINUE
          ENDIF
C
C  CHECK IF MODEL IS BLOWING UP
C
C  STOP PROGRAM IF CONC. > 100,000 mg/l (0.10 g/cm**3)
C
          NBLOW=NBLOW+1
          IF (NBLOW.EQ.100) THEN
            NBLOW=0
            CMAX=0.10
            DO 6360 J=3,JM-2
              DO 6360 I=2,IM-1
                IF (H(I,J).GT.0.01) THEN
                  TOTCON=CSED1(I,J,1)+CSED2(I,J,1)
                  IF (ABS(TOTCON).GT.CMAX) THEN
                  WRITE (IUPRT,6370)INT,I,J,CSED1(I,J,1),CSED2(I,J,1)
 6370             FORMAT (/5X,'******* PROGRAM EXECUTION STOPPED AT',
     +   I8,' TIMESTEP',/8X,'CONCENTRATION > 100,000 mg/l AT ',2I3,/8X,
     +   'C1 AND C2 =',2E12.4)
                  WRITE (6,6380)TIME
 6380             FORMAT (/5X,'$$$$ TIME (days) =',F10.4)
C
                  CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
C
                  STOP
                  ENDIF
                ENDIF
 6360       CONTINUE
          ENDIF
C
C***********************************************************
C
C  CHEM TRANSPORT 
C
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
C
C  COHESIVE CLASS
C
            IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
              DO 401 N=1,NUMQBCCH
                DUMB1(N)=PDIS(1,N)
 401          CONTINUE
C
              DO 410 N=1,NUMEBCCH
                DO 410 K=1,KBM1
                  DUMB2(N,K)=PBDRY(1,N,K)
 410          CONTINUE
C
C  FOR DIFFUSER INPUTS 
C
            DO 415 N=1,NUMDBCCH
                DUMB3(N)=PDIFF(1,N)
 415        CONTINUE
C
              CALL ADVCHEM(CHEM1B,CHEM1,DTI2,UF,DUMB1,DUMB2,DUMB3,
     +                     CHEMDRAT1)        
C
              CALL PROFSED(UF,DTI2,WSET1,CHEMBOT1)
C
              CALL BCOND(10,DTI2,1)
C
              DO 420 K=1,KBM1
                DO 420 J=1,JM
                  DO 420 I=1,IM
 
                    IF (UF(I,J,K).LE. 0.0 )THEN
                      CHEMBOT1(I,J)=0.0
                      CHEM1(I,J,K)=0.0
                    ELSE
                      CHEM1(I,J,K)=CHEM1(I,J,K)+.5*NU*(UF(I,J,K)+
     +                  CHEM1B(I,J,K)-2.*CHEM1(I,J,K))
                    ENDIF
                    CHEM1B(I,J,K)=CHEM1(I,J,K)
                    CHEM1(I,J,K)=UF(I,J,K)
 420          CONTINUE
            ENDIF
C
C  NON-COHESIVE CLASS
C
            IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
              IF (SEDTYPE.EQ.'SAND') THEN
                KK=1
              ELSE
                KK=2
              ENDIF
C
              DO 430 N=1,NUMQBCCH
                DUMB1(N)=PDIS(KK,N)
 430          CONTINUE
C
              DO 441 N=1,NUMEBCCH
                DO 441 K=1,KBM1
                  DUMB2(N,K)=PBDRY(KK,N,K)
 441          CONTINUE
C
C  FOR DIFFUSER INPUTS 
C
            DO 439 N=1,NUMDBCCH
                DUMB3(N)=PDIFF(KK,N)
 439        CONTINUE
C
              CALL ADVCHEM(CHEM2B,CHEM2,DTI2,UF,DUMB1,DUMB2,DUMB3,
     +                                      CHEMDRAT2)
C
              CALL PROFSED(UF,DTI2,WSET2,CHEMBOT2)
C
              CALL BCOND(10,DTI2,KK)
C
              DO 451 K=1,KBM1
                DO 451 J=1,JM
                  DO 451 I=1,IM
 
                    IF (UF(I,J,K).LE. 0.0 )THEN
                      CHEMBOT2(I,J)=0.0
                      CHEM2(I,J,K)=0.0
                    ELSE
                      CHEM2(I,J,K)=CHEM2(I,J,K)+.5*NU*(UF(I,J,K)+
     +                  CHEM2B(I,J,K)-2.*CHEM2(I,J,K))
                    ENDIF
                    CHEM2B(I,J,K)=CHEM2(I,J,K)
                    CHEM2(I,J,K)=UF(I,J,K)
 451          CONTINUE
            ENDIF
          ENDIF
C
C  CALCULATE FLUXES AT SEDIMENT-WATER INTERFACE
C
        IF (NSEDCT.EQ.NSBED) THEN
C
C  LIMIT DEPOSITION TO BATHYMETRIC DEPTH 
C
         DO 459 J=1,JM
         DO 459 I=1,IM
          COHTHK(I,J) =0.0
          NCOHTHK(I,J)=0.0
          IF (IBMSK(I,J).EQ.0) THEN
           DO 458 LL=1,LAYMAX
            COHTHK(I,J)=COHTHK(I,J)+TSED(LL,I,J)
 458       CONTINUE
           COHTHK(I,J)=COHTHK(I,J)/CBED(I,J)
          ELSE
           IF (IBMSK(I,J).EQ.1) THEN
            NCOHTHK(I,J)=NCOHTHK(I,J)+BEDTH(1,I,J)
            NCOHTHK(I,J)=100.*(NCOHTHK(I,J)/((CBED(I,J)/2.65)*
     +                   H1(I,J)*H2(I,J)))
           ENDIF
          ENDIF
          SEDTHK(I,J)=(COHTHK(I,J)+NCOHTHK(I,J))/100.
 459     CONTINUE
C
C  COHESIVE BED MODEL
C
          IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') CALL SEDFLX
C
C  NON-COHESIVE BED MODEL
C
          IF (IBED.EQ.1) CALL SUSLOD
C
C  CALCULATE CHEM FLUX AT SEDIMENT-WATER INTERFACE
C
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
            CALL CHEMFLX
C
C  CHEM BED MODEL
C
            CALL CHEMBED(CHEMDRAT1,CHEMDRAT2)
          ENDIF
C
          NSEDCT=0
        ENDIF
C
C***********************************************************
C
C  REORDER COHESIVE SEDIMENT BED LAYERS, IF NECESSARY
C
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (N24CNT.EQ.N24HR) CALL REORDR
        ENDIF
      ENDIF
C
C  MODIFY ETA AND VELOCITIES FOR NEXT TIMESTEP (HYDRO INFO READ FROM hqi_tran)
C
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
        DO 9920 J=2,JMM1
          DO 9920 I=2,IMM1
            ETB(I,J)=ET(I,J)
            ET(I,J)=ETF(I,J)
C
            ETF(I,J)=ETB(I,J)+2.*DTI*DETA(I,J)
C
            DT(I,J)=H(I,J)+ET(I,J)
 9920   CONTINUE
C
C  CALC NEW VELOCITIES
C
        DO 9930 K=1,KBM1
          DO 9930 J=2,JMM1
            DO 9930 I=2,IM
              U(I,J,K)=0.0
              DBAR=0.25*(DT(I,J)+DT(I-1,J))*(H2(I-1,J)+H2(I,J))
            IF (DBAR.GT.0.0) U(I,J,K)=XMFL3D(I,J,K)/DBAR
 9930   CONTINUE
C
        DO 9940 K=1,KBM1
          DO 9940 J=2,JM
            DO 9940 I=2,IMM1
              V(I,J,K)=0.0
              DBAR=0.25*(DT(I,J)+DT(I,J-1))*(H1(I,J-1)+H1(I,J))
              IF (DBAR.GT.0.0) V(I,J,K)=YMFL3D(I,J,K)/DBAR
 9940   CONTINUE
C
C  SKIP HYDRODYNAMICS IF READING FROM gcm_tran
C
        GOTO 9910
      ENDIF
C
      IF (TOR.EQ.'BAROTROPIC') GO TO 9910
C
C********************************************************************
C
C-----------------------------------------------------------------------
C         COMPUTE UF AND VF
C-----------------------------------------------------------------------
C
      CALL ADVU(DRHOX,ADVUU,DTI2)
      CALL ADVV(DRHOY,ADVVV,DTI2)
      CALL PROFU(DTI2)
      CALL PROFV(DTI2)
      CALL BCOND(3,DTI2,0)
C
      IF(BCTYPE.EQ.'MIXED  ') CALL MIXED(2)
C
      DO 369 J=1,JM 
      DO 369 I=1,IM 
 369  TPS(I,J)=0.0
      DO 370 K=1,KBM1
      DO 370 J=1,JM 
      DO 370 I=1,IM 
 370  TPS(I,J)=TPS(I,J)+(UF(I,J,K)+UB(I,J,K)-2.*U(I,J,K))*DZ(K)
      DO 372 K=1,KBM1
      DO 372 J=1,JM 
      DO 372 I=1,IM 
 372  U(I,J,K)=U(I,J,K)+.5*NU*(UF(I,J,K)+UB(I,J,K)-2.*U(I,J,K)
     .        -TPS(I,J))
      DO 373 J=1,JM 
      DO 373 I=1,IM 
 373  TPS(I,J)=0.0
      DO 374 K=1,KBM1
      DO 374 J=1,JM 
      DO 374 I=1,IM 
 374  TPS(I,J)=TPS(I,J)+(VF(I,J,K)+VB(I,J,K)-2.*V(I,J,K))*DZ(K)
      DO 376 K=1,KBM1
      DO 376 J=1,JM 
      DO 376 I=1,IM 
 376  V(I,J,K)=V(I,J,K)+.5*NU*(VF(I,J,K)+VB(I,J,K)-2.*V(I,J,K)
     .        -TPS(I,J))
      DO 377 K=1,KB 
      DO 377 J=1,JM 
      DO 377 I=1,IM 
      UB(I,J,K)=U(I,J,K)
      U(I,J,K)=UF(I,J,K)
      VB(I,J,K)=V(I,J,K)  
 377  V(I,J,K)=VF(I,J,K)
C
      CALL WREAL(DTI2)
C
 8200 CONTINUE
c
      DO 380 J=1,JM 
      DO 380 I=1,IM 
      EGB(I,J)=EGF(I,J)
      ETB(I,J)=ET(I,J)
      ET(I,J)=ETF(I,J)
      DT(I,J)=H(I,J)+ET(I,J)
      UTB(I,J)=UTF(I,J)
 380  VTB(I,J)=VTF(I,J)
C
9910  continue
C
      CALL ARCHIVE        
C
      IF(MOD(INT,IPRINT).EQ.0 .AND. INT.GE.IPRTSTART)  THEN
       CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
      ENDIF
C
 9000              CONTINUE
C
      OPEN (IUWRS,FORM='unformatted',FILE='startup')
C
      IF (TRACER.EQ.'INCLUDE') THEN
       IF (SEDTRAN.EQ.'NEGLECT') THEN
        WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +     CONC1,CONC1B
       ELSE
         IF (CHEMTRAN.EQ.'INCLUDE') THEN
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +     CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +        SEDMASS,SEDEP
         ELSE
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +     CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,TAU
         ENDIF
       ENDIF
      ELSE
       IF (SEDTRAN.EQ.'NEGLECT') THEN
        WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN
       ELSE
         IF (CHEMTRAN.EQ.'INCLUDE') THEN
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +     N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +        SEDMASS,SEDEP
         ELSE
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     +     N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +        TAUMAX,TAUCUR,TSED,
     +        EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +        FR,NCNT,NHRCNT,NDTCNT,TAU
         ENDIF
       ENDIF
      ENDIF
C
C********************************************************************
C
      CLOSE (IUWRS)
C 
 9100 CONTINUE
      CLOSE (IUT90)
      CLOSE (IUT91)
      CLOSE (IUT92)
      CLOSE (IUT93)
      CLOSE (IUT193)
      CLOSE (IUT94)
C
      IF (TRACER.EQ.'INCLUDE') THEN
        CLOSE (IUT96)
        CLOSE (IUT98)
        CLOSE (IUT99)
      ENDIF
C
C
C********************************************************************
C
c     !!!&&&CALL SYSTEM ('rm gcm_temp*')
C
C      IF(VAMAX.LT.10.) THEN
      IF(HYDTYPE.EQ.'EXTERNAL'.OR.(VAMAX.LT.10.and.elminx.gt.0.1)) THEN 
      WRITE(IUPRT,602) TIME
602   FORMAT(/2X,'JOB SUCCESSFULLY COMPLETED; TIME = ',1P1E10.2,' DAYS',
     . //)
C
      ELSE
      CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
      WRITE(IUPRT,'(''********************************************'')')
      WRITE(IUPRT,'(''********** ABNORMAL JOB END ****************'')')
      WRITE(IUPRT,'(''*********** USER TERMINATED ****************'')')
      WRITE(IUPRT,'(''********************************************'')')
      if(vamax.gt.10.0) 
     .WRITE(IUPRT,7500)IAMAX,JAMAX,VAMAX,TIME,FSM(IAMAX,JAMAX)  
 7500 FORMAT(1X,'I=',I3,3X,'J=',I3,3X,'MAX SPEED=',E12.3,'TIME=',
     .F10.4,'FSM= ',f5.0)
      if(elminx.lt.0.1) 
     .WRITE(IUPRT,7501)IAMAX,JAMAX,elmin,TIME,FSM(IAMAX,JAMAX)
 7501 FORMAT(1X,'I=',I3,3X,'J=',I3,3X,'MIN ELEV=',E12.4,'TIME=',
     .F10.4,'FSM= ',f5.0)
      ENDIF
C Copy Right
 7000 format(
     .'**************************************************************',/        
     .'*     ECOMSED MODEL                                          *',/
     .'*     VERSION 1.3                                            *',/
     .'*     FEBRUARY 2002                                          *',/
     .'**************************************************************',/
     .'*               Copyright (C) 2002, HydroQual, Inc.          ',/
     .'*                                                            ',/
     .'*                                                            ',/
     .'* As per the USER AGREEMENT AND DISCLAIMER to which all users',/
     .'* of ECOMSED have agreed,  no recipient of ECOMSED shall:    ',/
     .'*  (i) copyright or patent it in any form;                   ',/
     .'* (ii) redistribute it to others without HydroQual''s prior  ',/
     .'*      written authorization;                                ',/
     .'* (iii) make a monetary charge for it;                       ',/
     .'* (iv) violate or participate in the violation of the        ',/
     .'*      laws or regulations of the United States or other     ',/      
     .'*      governments, foreign or domestic, in regard to its use',/
     .'*      or distribution.                                      ',/
     .'**************************************************************',
     . //)
                                                         
C--------------------------------------------------------------------
 6110 FORMAT(/' BCTYPE = ',A7,', IS INCORRECTLY SPECIFIED.'//
     *    ' PLEASE FIX AND RESUBMIT'//)
 6111 FORMAT(/' BCTYPE = ',A7,' AND TLAG = ',F10.2,
     *    ' ARE INCORRECTLY SPECIFIED.'//
     *    ' PLEASE FIX AND RESUBMIT'//)
 6112 FORMAT(/' WAVEDYN = ',A8,' IS INCORRECTLY SPECIFIED.'//
     *    ' PLEASE FIX AND RESUBMIT'//)

 9200 STOP
      END
