      SUBROUTINE SOLAR(JD,FJD,R,DLT,ALP,SLAG,N,HANG,TAUDA,COSZ,SOLALT,
     1 PHIT)
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
C    *******************************************************************
C    *                                                                 *
C    *                            S O L A R                            *
C    *                                                                 *
C    *******************************************************************
C
C.....SOLAR COMPUTES RADIUS VECTOR, DECLINATION AND RIGHT ASCENSION OF
C.....SUN, EQUATION OF TIME, AND HOUR ANGLE OF SUN AT SUNSET FOR N
C.....EQUALLY SPACED LATITUDES GIVEN JULIAN DAY AND FRACTION.
C
      PARAMETER (JMT=1)
      DIMENSION ALAT(JMT)
C
                              D I M E N S I O N
     1   HANG( JMT ), COSZ( JMT ),  TAUDA( JMT ), SOLALT (JMT)
C
                                   D A T A
     1   TPI/6.2831853/,     HPI/1.5707963/,
     2   RAD/57.29578/,      CYEAR/365.25/,      CCR/1.E-6/
C
C.....TPP = DAYS BETWEEN EPOCH AND PERIHELION PASSAGE OF 1900
C.....SVT6 = DAYS BETWEEN PERIHELION PASSAGE AND MARCH EQUINOX OF 1900
C.....JDOR = JD OF EPOCH WHICH IS JANUARY 0, 1900 AT 12 HOURS UT
C
                                   D A T A
     1   TPP/1.55/,          SVT6/78.035/,       JDOR/2415020/
      DATA ITERS/100/
C
C    *******************************************************************
C
C
      PI=ACOS(-1.0)
      ALAT(1)=PHIT
      DAT=FLOAT(JD-JDOR)-TPP+FJD
C
C    COMPUTES TIME IN JULIAN CENTURIES AFTER EPOCH
C
      T=FLOAT(JD-JDOR)/36525.0
C
C    COMPUTES LENGTH OF ANOMALISTIC AND TROPICAL YEARS (MINUS 365 DAYS)
C
      YEAR=0.25964134+0.304E-5*T
      TYEAR=0.24219879-0.614E-5*T
C
C    COMPUTES ORBIT ECCENTRICITY AND ANGLE OF EARTH'S INCLINATION FROM T
C
      EC=0.01675104-(0.418E-4+0.126E-6*T)*T
      ANGIN=23.452294-(0.0130125+0.164E-5*T)*T
      ADOR=JDOR
      JDOE=ADOR+(SVT6*CYEAR)/(YEAR-TYEAR)
C
C    DELEQN=UPDATED SVT6 FOR CURRENT DATE
C
C
C
      DELEQN=FLOAT(JDOE-JD)*(YEAR-TYEAR)/CYEAR
C
      YEAR=YEAR+365.0
C
      SNI=SIN(ANGIN/RAD)
      TINI=1.0/TAN(ANGIN/RAD)
      ER=SQRT((1.0+EC)/(1.0-EC))
      QQ=DELEQN*TPI/YEAR
C
C    DETERMINE TRUE ANOMALY AT EQUINOX
C
      E=1.0
C
      DO 200 ITER=1,ITERS
      EP=E-(E-EC*SIN(E)-QQ)/(1.0-EC*COS(E))
      CD=ABS(E-EP)
      E=EP
      IF(CD.LE.CCR) GO TO 250
 200  CONTINUE
      PRINT 400
 250  CONTINUE
      HE=0.5*E
      EQ=2.0*ATAN(ER*TAN(HE))
C
C    DATE=DAYS SINCE LAST PERIHELION PASSAGE
C
      DATE=MOD(DAT,YEAR)
C
C    SOLVE ORBIT EQUATIONS BY NEWTON'S METHOD
C
      EM=TPI*DATE/YEAR
      E=1.0
C
      DO 300 ITER=1,ITERS
      EP=E-(E-EC*SIN(E)-EM)/(1.0-EC*COS(E))
      CR=ABS(E-EP)
      E=EP
      IF(CR.LE.CCR) GO TO 350
 300  CONTINUE
      PRINT 400
 350  CONTINUE
C
      R=1.0-EC*COS(E)
      HE=0.5*E
      W=2.0*ATAN(ER*TAN(HE))
      SIND=SNI*SIN(W-EQ)
      DLT=ASIN(SIND)
      ALP=ASIN(TAN(DLT)*TINI)
      TST=COS(W-EQ)
C
      IF(TST.LT.0.0) ALP=PI-ALP
C
      IF(ALP.LT.0.0) ALP=ALP+TPI
C
      SUN=TPI*(DATE-DELEQN)/YEAR
      IF (SUN.LT.0.0)  SUN=SUN+TPI
      SLAG=SUN-ALP-0.03255
C
C    COMPUTE HOUR ANGLE OF SUNSET AT ALL LATITUDES
C
C
      DO 1 I=1,N
      PH=ALAT(I)
      SS=SIN(PH)*SIN(DLT)
      CC=COS(PH)*COS(DLT)
C
      SOLALT(I)= SS + CC
      SOLALT(I) = ASIN(SOLALT(I)) * RAD
C
      IF(DLT.EQ.0.0) GO TO 16
C
      AP=ABS(PH)
      EPS=ABS(AP-HPI)
      IF(EPS.GT.CCR) GO TO 14
C
      H=HPI*ABS(AP/PH+ABS(DLT)/DLT)
      GO TO 5
C
   14 AR=-SS/CC
      AC=ABS(AR)
      IF(AC-1.0+CCR) 4,2,3
C
  2   H=(AC-AR)*HPI
C
      GO TO 5
C
  3   IF(AR.LT.0.0) GO TO 25
C
      H=0.0
      GO TO 5
C
 25   H=PI
      GO TO 5
C
C
 16   H=HPI
      GO TO 5
C
  4   H=ACOS(AR)
  5   HANG(I)=H
      TAUDA(I)=MAX(H/PI,0.0)
C
      IF(H.EQ.0.0) GO TO 100
C
C	write(*,*)'H in Solar=',H
      COSZ(I)=MAX((SS+CC*SIN(H)/H),0.0)
      GO TO 1
C
C
C
  100 COSZ(I)=0.0
    1 CONTINUE
C
      RETURN
 400  FORMAT('0 MAXIMUM ITERATIONS IN SUBROUTINE SOLAR')
      END
