      SUBROUTINE BCDATA
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
      DIMENSION COM(80),YC(EBCM)
      DIMENSION ZM(KBM1),T1(KSL),T2(KBM1),S1(KSL),S2(KBM1)
      DIMENSION IDBC(DBCM),TEMPD(DBCM)
      DIMENSION SHFLX(IM,JM)    ! SYNOPTIC HEAT FLUX (in Watt/M**2)
C
C  FOR CONSERVATIVE TRACER:  OPEN B.C.
C
      REAL C1D(KSL),C2D(KBM1),WTEMP1(IM,JM),WTEMP2(IM,JM),WTEMP3(IM,JM)
      REAL*4 THZERO
C
C-----------------------------------------------------------------------
C-------- READ IN STANDARD LEVELS HERE ---------------------------------
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
 11   FORMAT(80A1)
 12   FORMAT(/1X,80A1/)
C
      READ(IURUN,24)  IKSL
      WRITE(IUPRT,41) IKSL
   41 FORMAT(' KSL = ',I5,/)
      IF(IKSL.NE.KSL) THEN
       WRITE(6,42) IKSL, KSL
       !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
   42 FORMAT(//' NUMBER OF STANDARD LEVELS IN RUN_DATA',I5,' (IKSL)'/
     .         '           DO NOT EQUAL'/
     .         ' NUMBER OF STANDARD LEVELS IN COMDECK ',I5,' (KSL)'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
C
      READ(IURUN,77)  (DPTHSL(K),K=1,KSL)
      WRITE(IUPRT,77) (DPTHSL(K),K=1,KSL)
C
C-----------------------------------------------------------------------
C-------- INITIAL TEMPERATURE AND SALINITY DATA ------------------------
C-------- TEMPERATURE IN C, SALINITY IN PSU ----------------------------
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,13)  OPTTSI
      WRITE(IUPRT,13) OPTTSI
  13  FORMAT(A20)
      IF(OPTTSI(1:5).EQ.'FIXED') THEN
       READ(IURUN,77)  (TSI(K),K=1,KSL)
       WRITE(IUPRT,77) (TSI(K),K=1,KSL)
       READ(IURUN,77)  (SSI(K),K=1,KSL)
       WRITE(IUPRT,77) (SSI(K),K=1,KSL)
      ENDIF
C
C-----------------------------------------------------------------------
C-------- READ IN BOUNDARY CONDITIONS HERE -----------------------------
C-----------------------------------------------------------------------
C-------- ELEVATION BOUNDARY -------------------------------------------
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 128
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,24)  NUMEBC,OPTEBC
      WRITE(IUPRT,24) NUMEBC,OPTEBC
 24   FORMAT(I5,1X,A20)
C
      IF(NUMEBC.EQ.0) GO TO 128
      IF(NUMEBC.GT.EBCM) THEN
       WRITE(IUPRT,60) NUMEBC,EBCM
       !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  60  FORMAT(//' NUMEBC (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' EBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      IF(OPTEBC(1:5).NE.'DATA ') THEN
      DO 110 N=1,NUMEBC
      READ(IURUN,9)  IETA(N),JETA(N),ICON(N),JCON(N),EMEAN(N)
      WRITE(IUPRT,9) IETA(N),JETA(N),ICON(N),JCON(N),EMEAN(N)
      READ(IURUN,7)  (AMP(N,I),I=1,6)
      WRITE(IUPRT,7) (AMP(N,I),I=1,6)
      READ(IURUN,7)  (PHASE(N,I),I=1,6)
      WRITE(IUPRT,7) (PHASE(N,I),I=1,6)
c     DO 105 I=1,6
 110  CONTINUE
c
      NDAY=IFIX((IEND-ISTART)*DAYI*DTI)
      DO N=1,NUMEBC
       YC(N)=YGRID(IETA(N),JETA(N))
      END DO
      THZERO=FLOAT(INT)*DTI/3600.0
      WRITE(IUPRT,*)'THZERO=',THZERO
      CALL PTIDE(EBCM,NUMEBC,7,YC,AMP,PHASE,EMEAN,IDA,IMO,
     .   IYR,IHOUR,NDAY+1,THZERO)
 9    FORMAT(4I5,1F10.5)
 7    FORMAT(6F10.5)
      ELSE
      READ(IURUN,76)  (IETA(N),JETA(N),ICON(N),JCON(N),N=1,NUMEBC)
      WRITE(IUPRT,76) (IETA(N),JETA(N),ICON(N),JCON(N),N=1,NUMEBC)
 76   FORMAT(16I5)
      DO 122 M=1,100000
      READ(IURUN,77,ERR=125) TIME
      WRITE(IUPRT,77) TIME
      OTIME=TIME
      READ(IURUN,77)  (EBDRY(N),N=1,NUMEBC)
      WRITE(IUPRT,77) (EBDRY(N),N=1,NUMEBC)
      WRITE(IUT90,78) TIME
 122  WRITE(IUT90,78) (EBDRY(N),N=1,NUMEBC)
 125  BACKSPACE IURUN
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE ETA BC TO COMPLETE THE SIMULATIONS ' 
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
   
      END IF
C
C-----------------------------------------------------------------------
C-------- TEMPERATURE AND SALINITY BOUNDARY ----------------------------
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      DO 435 M=1,100000
      READ(IURUN,77,ERR=127) TIME
      WRITE(IUPRT,77) TIME
      WRITE(IUT94,78) TIME
      OTIME=TIME
      DO 239 N=1,NUMEBC
        READ(IURUN,79) ITAS(N),JTAS(N),
     .               (TBDRYSL(N,K),K=1,KSL),(SBDRYSL(N,K),K=1,KSL)
C
            IF ((ITAS(N).NE.IETA(N)).OR.(JTAS(N).NE.JETA(N))) THEN
              WRITE (IUPRT,5701)N,ITAS(N),JTAS(N),IETA(N),JETA(N)
 5701         FORMAT (/5X,'******  MODEL EXECUTION STOPPED  *****',
     +         //5X,'ORDER OF TEMPERATURE/SALINITY OPEN B.C.',
     +         /5X,'SPECIFIED INCORRECTLY',
     +         /5X,'B.C. NUMBER ',I4,' HAS ITAS, JTAS =',2I4,
     +         /5X,'SHOULD BE IETA, JETA =',2I4,
     +         //5X,'PLEASE CORRECT AND RESUBMIT')
               STOP
             ENDIF
C
        WRITE(IUPRT,9) ITAS(N),JTAS(N)
        WRITE(IUPRT,80)((TBDRYSL(N,K)),K=1,KSL)
        WRITE(IUPRT,80)((SBDRYSL(N,K)),K=1,KSL)
        II=ITAS(N)
        JJ=JTAS(N)
        DO 243 K=1,KBM1
243       ZM(K)=ZZ(K)*H(II,JJ)
        DO  241 K=1,KSL
          T1(K)=TBDRYSL(N,K)
241       S1(K)=SBDRYSL(N,K)
      CALL SINTER(DPTHSL,T1,ZM,T2,KSL,KBM1)
      CALL SINTER(DPTHSL,S1,ZM,S2,KSL,KBM1)
        DO 246 K=1,KBM1
          TBDRY(N,K)=T2(K)
246       SBDRY(N,K)=S2(K)
        WRITE(IUT94,78) (TBDRY(N,K),K=1,KBM1)
        WRITE(IUT94,78) (SBDRY(N,K),K=1,KBM1)
239   CONTINUE
435   CONTINUE
127   BACKSPACE IURUN
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE T/S BC TO COMPLETE THE SIMULATIONS ' 
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
79    FORMAT(2I5,100F5.0)
80    FORMAT(50F5.0)
C
C  FOR CONSERVATIVE TRACER:  OPEN B.C.
C
 128  IF (TRACER.EQ.'INCLUDE') THEN
        OPEN (IUT501,FILE='gcm_temp501')
C
        READ(IUT401,11)  (COM(I),I=1,80)
        WRITE(IUPRT,12) (COM(I),I=1,80)
        READ (IUT401,24)NUMEBCTR
C
        IF (NUMEBCTR.EQ.0) GOTO 1128
C
        DO 6010 M=1,100000
          READ(IUT401,77,ERR=6020) TIME
          WRITE(IUPRT,77) TIME
          OTIME=TIME
          WRITE (IUT501,78)TIME
          DO 6030 N=1,NUMEBCTR
            READ(IUT401,179) II,JJ,IIC,JJC,(CBDRYSL1(N,K),K=1,KSL)
 179        FORMAT (4I5,100F5.0)
C
            ITRED(N)=II
            JTRED(N)=JJ
            ITREC(N)=IIC
            JTREC(N)=JJC
C
            WRITE(IUPRT,9) II,JJ,IIC,JJC
            WRITE(IUPRT,80)((CBDRYSL1(N,K)),K=1,KSL)
            DO 6050 K=1,KBM1
              ZM(K)=ZZ(K)*H(II,JJ)
 6050       CONTINUE
            DO 6060 K=1,KSL
              C1D(K)=CBDRYSL1(N,K)
 6060       CONTINUE                
            CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
            DO 6070 K=1,KBM1
              CBDRY1(N,K)=C2D(K)
 6070       CONTINUE
C
            WRITE (IUT501,78) (CBDRY1(N,K),K=1,KBM1)
C
 6030     CONTINUE
 6010   CONTINUE
 6020   BACKSPACE IUT401
      ENDIF
C
C
C*******************************************************************
C
C  SEDIMENT TRANSPORT
C
 1128 IF (SEDTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (IUT502,FILE='gcm_temp502')
          READ(IUT402,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT402,24)NUMEBCSE
C
          IF (NUMEBCSE.EQ.0) GOTO 2127
C
          DO 6110 M=1,100000
            READ(IUT402,77,ERR=6120) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT502,78)TIME
            DO 6130 N=1,NUMEBCSE
              READ(IUT402,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ISEED(N)=II
              JSEED(N)=JJ
              ISEEC(N)=IIC
              JSEEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
C
C  CONVERT FROM mg/l TO g/cm**3
C
              DO 6150 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000000.
 6150         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6160 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6160         CONTINUE
C
              WRITE(IUT502,78) (CBDRY(1,N,K),K=1,KBM1)
C
 6130       CONTINUE
 6110     CONTINUE
 6120     BACKSPACE IUT402
        ENDIF
C
 2127   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (IUT503,FILE='gcm_temp503')
C
          READ(IUT403,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT403,24)NUMEBCSE
C
          IF (NUMEBCSE.EQ.0) GOTO 2128
C
          DO 6170 M=1,100000
            READ(IUT403,77,ERR=6180) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT503,78)TIME
C
            DO 6190 N=1,NUMEBCSE
              READ(IUT403,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ISEED(N)=II
              JSEED(N)=JJ
              ISEEC(N)=IIC
              JSEEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
C
C  CONVERT FROM mg/l TO g/cm**3
C
              DO 6200 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000000.
 6200         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6210 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6210         CONTINUE
C
              WRITE(IUT503,78) (CBDRY(1,N,K),K=1,KBM1)
C
 6190       CONTINUE
 6170     CONTINUE
 6180     BACKSPACE IUT403
        ENDIF
      ENDIF
C
C*******************************************************************
C
C  CHEM TRANSPORT
C
 2128 IF (CHEMTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT504=504
          OPEN (IUT504,FILE='gcm_temp504')
C
          READ(IUT404,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT404,24)NUMEBCCH
C
          IF (NUMEBCCH.EQ.0) GOTO 2129
C
          DO 6310 M=1,100000
            READ(IUT404,77,ERR=6320) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT504,78)TIME
C
            DO 6330 N=1,NUMEBCCH
              READ(IUT404,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ICHED(N)=II
              JCHED(N)=JJ
              ICHEC(N)=IIC
              JCHEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
C
C  CONVERT FROM ug/l TO ug/cm**3
C
              DO 6350 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000.
 6350         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6360 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6360         CONTINUE
C
              WRITE(IUT504,78) (CBDRY(1,N,K),K=1,KBM1)
C
 6330       CONTINUE
 6310     CONTINUE
 6320     BACKSPACE IUT404
        ENDIF
C
 2129   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT505=505
          OPEN (IUT505,FILE='gcm_temp505')
C
          READ(IUT405,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT405,24)NUMEBCCH
C
          IF (NUMEBCCH.EQ.0) GOTO 126 
C
          DO 6370 M=1,100000
            READ(IUT405,77,ERR=6380) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT505,78)TIME
C
            DO 6390 N=1,NUMEBCCH
              READ(IUT405,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ICHED(N)=II
              JCHED(N)=JJ
              ICHEC(N)=IIC
              JCHEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
C
C  CONVERT FROM ug/l TO ug/cm**3
C
              DO 6400 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000.
 6400         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6410 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6410         CONTINUE
C
              WRITE(IUT505,78) (CBDRY(1,N,K),K=1,KBM1)
C
 6390       CONTINUE
 6370     CONTINUE
 6380     BACKSPACE IUT405
        ENDIF
      ENDIF
C
C*******************************************************************
C
 126  CONTINUE
C
C-----------------------------------------------------------------------
C-------- RIVER/DAM AND ONSHORE INTAKE/OUTFALL DISCHARGE BOUNDARY ------
C-----------------------------------------------------------------------
C
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
         NUMEBC = 0
         NUMQBC=0
         GOTO 138
      ENDIF
C
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,8)  NUMQBC
      WRITE(IUPRT,8) NUMQBC
C
      IF(NUMQBC.EQ.0) GO TO 138
      IF(NUMQBC.GT.QBCM) THEN
       WRITE(IUPRT,61) NUMQBC,QBCM
       !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  61  FORMAT(//' NUMQBC (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' QBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      DO 321 N=1,NUMQBC
      READ(IURUN,5)  IQD(N),JQD(N),IQC(N),JQC(N),(VQDIST(N,K),K=1,KBM1)
321   WRITE(IUPRT,6) IQD(N),JQD(N),IQC(N),JQC(N),(VQDIST(N,K),K=1,KBM1)
 8    FORMAT(2I5,4F10.5)
 5    FORMAT(4I5,20F5.0)
 6    FORMAT(4I5,/,20F5.0)
      DO 130 M=1,100000
      READ(IURUN,77,ERR=135) TIME
      WRITE(IUPRT,77) TIME
      OTIME=TIME
      READ(IURUN,77)  (QDIS(N),N=1,NUMQBC)
      READ(IURUN,77)  (TDIS(N),N=1,NUMQBC)
      READ(IURUN,77)  (SDIS(N),N=1,NUMQBC)
      WRITE(IUPRT,77) (QDIS(N),N=1,NUMQBC)
      WRITE(IUPRT,77) (TDIS(N),N=1,NUMQBC)
      WRITE(IUPRT,77) (SDIS(N),N=1,NUMQBC)
      WRITE(IUT91,78) TIME
      WRITE(IUT91,78) (QDIS(N),N=1,NUMQBC)
      WRITE(IUT91,78) (TDIS(N),N=1,NUMQBC)
 130  WRITE(IUT91,78) (SDIS(N),N=1,NUMQBC)
 135  BACKSPACE IURUN
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE RIVER BC TO COMPLETE THE SIMULATIONS' 
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
C
C*******************************************************************
C
C  FOR CONSERVATIVE TRACER:  RIVER DISCHARGE
C
 138  IF (TRACER.EQ.'INCLUDE') THEN
        OPEN (IUT601,FILE='gcm_temp601')
C
        READ(IUT401,11)  (COM(I),I=1,80)
        WRITE(IUPRT,12) (COM(I),I=1,80)
        READ (IUT401,24)NUMQBCTR
        WRITE (IUPRT,24)NUMQBCTR
C
        IF (NUMQBCTR.EQ.0) GOTO 1129
C
        DO 6509 N=1,NUMQBCTR
          READ (IUT401,179)ITRQD(N),JTRQD(N),ITRQC(N),JTRQC(N)
          WRITE (IUPRT,179)ITRQD(N),JTRQD(N),ITRQC(N),JTRQC(N)
 6509   CONTINUE
C
        DO 6510 M=1,100000
          READ(IUT401,77,ERR=6520) TIME
          READ(IUT401,77)  (CDIS1(N),N=1,NUMQBCTR)
          WRITE(IUPRT,77) TIME
          WRITE (IUPRT,77)  (CDIS1(N),N=1,NUMQBCTR)
          WRITE (IUT601,78)TIME
          WRITE (IUT601,78)  (CDIS1(N),N=1,NUMQBCTR)
C
 6510   CONTINUE
 6520   BACKSPACE IUT401
      ENDIF
C
C*******************************************************************
C
C  SEDIMENT TRANSPORT
C
 1129 IF (SEDTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT602=602
          OPEN (IUT602,FILE='gcm_temp602')
C
          READ(IUT402,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT402,24)NUMQBCSE
          WRITE (IUPRT,24)NUMQBCSE
C
          IF (NUMQBCSE.EQ.0) GOTO 1229
C
          DO 6529 N=1,NUMQBCSE
            READ (IUT402,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
            WRITE (IUPRT,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
 6529     CONTINUE
C
          DO 6530 M=1,100000
            READ(IUT402,77,ERR=6540) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT402,77)  (CDIS(1,N),N=1,NUMQBCSE)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCSE)
C
C  CONVERT FROM mg/l TO g/cm**3
C
            DO 6550 NN=1,NUMQBCSE
              CDIS(1,NN)=CDIS(1,NN)/1000000.
 6550       CONTINUE
C
            WRITE (IUT602,78)TIME
            WRITE (IUT602,78)  (CDIS(1,N),N=1,NUMQBCSE)
C
 6530     CONTINUE
 6540     BACKSPACE IUT402
        ENDIF
C
 1229   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT603=603
          OPEN (IUT603,FILE='gcm_temp603')
C
          READ(IUT403,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT403,24)NUMQBCSE
          WRITE (IUPRT,24)NUMQBCSE
C
          IF (NUMQBCSE.EQ.0) GOTO 1239
C
          DO 6559 N=1,NUMQBCSE
            READ (IUT403,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
            WRITE (IUPRT,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
 6559     CONTINUE
C
          DO 6560 M=1,100000
            READ(IUT403,77,ERR=6570) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT403,77)  (CDIS(1,N),N=1,NUMQBCSE)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCSE)
C
C  CONVERT FROM mg/l TO g/cm**3
C
            DO 6580 NN=1,NUMQBCSE
              CDIS(1,NN)=CDIS(1,NN)/1000000.
 6580       CONTINUE
C
            WRITE (IUT603,78)TIME
            WRITE (IUT603,78)  (CDIS(1,N),N=1,NUMQBCSE)
C
 6560     CONTINUE
 6570     BACKSPACE IUT403
        ENDIF
      ENDIF
C
C*******************************************************************
C
C  CHEMICAL TRANSPORT
C
 1239 IF (CHEMTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT604=604
          OPEN (IUT604,FILE='gcm_temp604')
C
          READ(IUT404,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT404,24)NUMQBCCH
          WRITE (IUPRT,24)NUMQBCCH
C
          IF (NUMQBCCH.EQ.0) GOTO 1329
C
          DO 6589 N=1,NUMQBCCH
            READ (IUT404,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
            WRITE (IUPRT,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
 6589     CONTINUE
C
          DO 6590 M=1,100000
            READ(IUT404,77,ERR=6600) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT404,77)  (CDIS(1,N),N=1,NUMQBCCH)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCCH)
C
C  CONVERT FROM ug/l TO ug/cm**3
C
            DO 6610 NN=1,NUMQBCCH
              CDIS(1,NN)=CDIS(1,NN)/1000.
 6610       CONTINUE
C
            WRITE (IUT604,78)TIME
            WRITE (IUT604,78)  (CDIS(1,N),N=1,NUMQBCCH)
C
 6590     CONTINUE
 6600     BACKSPACE IUT404
        ENDIF
C
 1329   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT605=605
          OPEN (IUT605,FILE='gcm_temp605')
C
          READ(IUT405,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT405,24)NUMQBCCH
          WRITE (IUPRT,24)NUMQBCCH
C
          IF (NUMQBCCH.EQ.0) GOTO 1330
C
          DO 6599 N=1,NUMQBCCH
            READ (IUT405,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
            WRITE (IUPRT,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
 6599     CONTINUE
C
          DO 6620 M=1,100000
            READ(IUT405,77,ERR=6630) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT405,77)  (CDIS(1,N),N=1,NUMQBCCH)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCCH)
C
C  CONVERT FROM ug/l TO ug/cm**3
C
            DO 6640 NN=1,NUMQBCCH
              CDIS(1,NN)=CDIS(1,NN)/1000.
 6640       CONTINUE
C
            WRITE (IUT605,78)TIME
            WRITE (IUT605,78)  (CDIS(1,N),N=1,NUMQBCCH)
C
 6620     CONTINUE
 6630     BACKSPACE IUT405
        ENDIF
      ENDIF
C
C*******************************************************************
C
 1330 CONTINUE
      DO 140 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            H2(IC,JC)=H2(ID,JD)
            IF(IC.GT.ID) THEN
              DUM(IC,JC)=1.0
            ELSE
              DUM(ID,JD)=1.0
            ENDIF
      ELSE
            H1(IC,JC)=H1(ID,JD)
            IF(JC.GT.JD) THEN
              DVM(IC,JC)=1.0
            ELSE
              DVM(ID,JD)=1.0
            ENDIF
      ENDIF
            H (IC,JC)=H (ID,JD)
            FSM(IC,JC)=FSM(ID,JD)
  140    CONTINUE
C
      IF(HORZMIX.EQ.'CONSTANT  ') THEN
      DO 220 N=1,NUMQBC
      IC=IQC(N)
      JC=JQC(N)
      AAM2D(IC,JC)=0.0
      DO 220 K=1,KBM1
      AAM(IC,JC,K)=0.0
  220 CONTINUE
      ENDIF
C
 136  CONTINUE
C
C-----------------------------------------------------------------------
C-------- OFFSHORE INTAKE/OUTFALL(DIFFUSER) BOUNDARY -------------------
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
         NUMDBC=0
         GOTO 3700
      ENDIF
C
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,8)  NUMDBC1
      WRITE(IUPRT,8) NUMDBC1
      IF(NUMDBC1.EQ.0) GO TO 156
      IF(NUMDBC1.GT.DBCM) THEN
       WRITE(IUPRT,62) NUMDBC1,DBCM
       !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  62  FORMAT(//' NUMDBC1 (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' DBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      DO N=1,NUMDBC1
        READ(IURUN,51)  IDD(N),JDD(N),(VDDIST(N,K),K=1,KBM1)
        WRITE(IUPRT,52) IDD(N),JDD(N),(VDDIST(N,K),K=1,KBM1)
      ENDDO
51    FORMAT(2I5,20F5.0)
52    FORMAT(2I5,/,20F5.0)
      DO M=1,100000
        READ(IURUN,77,ERR=155) TIME
        WRITE(IUPRT,77) TIME
        READ(IURUN,77)  (QDIFF(N),N=1,NUMDBC1)
        READ(IURUN,77)  (TDIFF(N),N=1,NUMDBC1)
        READ(IURUN,77)  (SDIFF(N),N=1,NUMDBC1)
        WRITE(IUPRT,77) (QDIFF(N),N=1,NUMDBC1)
        WRITE(IUPRT,77) (TDIFF(N),N=1,NUMDBC1)
        WRITE(IUPRT,77) (SDIFF(N),N=1,NUMDBC1)
        WRITE(IUT92,78) TIME
        OTIME=TIME
        WRITE(IUT92,78) (QDIFF(N),N=1,NUMDBC1)
        WRITE(IUT92,78) (TDIFF(N),N=1,NUMDBC1)
        WRITE(IUT92,78) (SDIFF(N),N=1,NUMDBC1)
      ENDDO
 155  BACKSPACE IURUN

      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     &  'NEED MORE DIFFUSER BC TO COMPLETE THE SIMULATIONS' 
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
c
 156  CONTINUE
c
C-------- OFFSHORE INTAKE/OUTFALL(DIFFUSER) BOUNDARY IN LOOPS------------------
c
      READ (IURUN,11) (COM(I),I = 1,80)
      WRITE(IUPRT,12) (COM(I),I = 1,80)
      READ (IURUN,8) NUMDBC2
      NUMDBC = NUMDBC1 + NUMDBC2
      WRITE (IUPRT,8) NUMDBC2,NUMDBC
      IF (NUMDBC2.EQ.0) GOTO 3700
      IF ((NUMDBC1+NUMDBC2).GT.DBCM) THEN
        WRITE (IUPRT,63) NUMDBC1+NUMDBC2, DBCM
        !!!&&&CALL SYSTEM('rm gcm_temp9*')
        STOP
      END IF
  63  FORMAT(//' NUMDBC (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' DBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      DO N = NUMDBC1+1, NUMDBC1+NUMDBC2
        READ (IURUN,51) IDD(N), JDD(N), (VDDIST(N,K),K = 1,KBM1)
        WRITE (IUPRT,52) IDD(N), JDD(N), (VDDIST(N,K),K = 1,KBM1)
      ENDDO    
      DO M = 1, 100000
        READ (IURUN,77,ERR=370) TIME
        WRITE (IUPRT,77) TIME
        OTIME=TIME
        READ (IURUN,77) (QDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        READ (IURUN,77) (TDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        READ (IURUN,77) (SDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
C
        WRITE (IUPRT,77) (QDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        WRITE (IUPRT,77) (TDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        WRITE (IUPRT,77) (SDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
C
        WRITE (IUT96,78) TIME
        WRITE (IUT96,78) (QDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        WRITE (IUT96,78) (TDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        WRITE (IUT96,78) (SDIFF(N),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
C
      ENDDO 
  370 BACKSPACE IURUN
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     .  'NEED MORE DIFFUSER IN LOOP BC TO COMPLETE THE SIMULATIONS'
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
c
 3700 CONTINUE
C
      IF (TRACER.NE.'INCLUDE') GO TO 3140  ! GO TO SED/MET INPUT DECK
      READ(IUT401,11) (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IUT401,179) NUMDBCTR1
      WRITE (IUPRT,179)NUMDBCTR1
      IF(NUMDBCTR1.EQ.0) GO TO 3130   ! GO TO TRACER-IN-LOOP DECK
      IF (NUMDBCTR1.GT.DBCM) THEN
        WRITE (IUPRT,64) NUMDBCTR1, DBCM
        !!!&&&CALL SYSTEM('rm gcm_temp9*')
        STOP
      END IF

      DO 3436 N=1,NUMDBCTR1
        READ (IUT401,179)ITRDD,JTRDD
        WRITE (IUPRT,179)ITRDD,JTRDD
C
C   FIND PROPER DIFFUSER FLOW MATCH
C   2 OPTIONS: IF NUMDBCTR1=NUMDBC1, USER WANTS TO USE SAME ORDER OF DIFFUSERS
C               SPECIFIED IN DIFFUSER INFLOW CONDITION  (THIS OPTION MUST BE
C               USED IF THERE ARE MULTIPLE DIFFUSERS ARE ASSIGNED IN A
C               GRID CELL.  (ORIGINAL CODING)
C              IF NUMDBCTR1 .NE. NUMDBC1 THEN, PROGRAM ASSUMES THAT USER IS
C                ASSIGNING DIFFUSER TRACER INPUT SELECTIVELY.
C
        IF(NUMDBCTR1.EQ.NUMDBC1) THEN
          IDBC(N)=N
          IF(ITRDD.NE.IDD(N).OR.JTRDD.NE.JDD(N))THEN
            WRITE(IUPRT,2481) ITRDD,JTRDD
            !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
            STOP
          END IF
        ELSE
          DO NQ=1,NUMDBC
            IF(ITRDD.EQ.IDD(NQ).AND.JTRDD.EQ.JDD(NQ)) THEN
              IDBC(N)=NQ
              GO TO 3436
            END IF
          ENDDO 
       END IF
 3436 CONTINUE
 2481 FORMAT(' ONE OF DIFFUSER TRACER INPUT DOES NOT MATCH WITH
     .DIFFUSER INFLOW LOCATIONS. ITS I,J ARE',2I5)

      DO N=1,NUMDBCTR1
        CDIFF1(N)=0.0
      END DO
      OPEN(IUT98,FILE='gcm_temp98')
C
      DO M=1,100000
        READ(IUT401,77,ERR=3128) TIME
        READ(IUT401,77)  (TEMPD(N),N=1,NUMDBCTR1)
        DO N=1,NUMDBCTR1
          CDIFF1(IDBC(N))=TEMPD(N)
        END DO
C
        WRITE(IUPRT,77) TIME
        OTIME=TIME
        WRITE (IUPRT,77)  (CDIFF1(N),N=1,NUMDBCTR1)
C
        WRITE (IUT98,78) TIME
        WRITE (IUT98,78)  (CDIFF1(N),N=1,NUMDBCTR1)
C
      ENDDO    

 3128 BACKSPACE IUT401
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     .  'NEED MORE DIFFUSER TRACER BC TO COMPLETE THE RUN'
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
 3130 CONTINUE
C
C FOR TRACER LOADINGS IN DIFFUSER-IN-LOOP  
C-------- OFFSHORE INTAKE/OUTFALL(DIFFUSER IN LOOPS) TRACER BOUNDARY------------------
C CAUTION !!! DIFFUSER-IN-LOOP TRACER INPUTS SHOULD MATCH WITH THE
C             DIFFUSER-IN-LOOP
C
      READ (IUT401,11) (COM(I),I = 1,80)
      WRITE(IUPRT,12) (COM(I),I = 1,80)
      READ (IUT401,8) NUMDBCTR2
C
      IF(NUMDBCTR2.NE.NUMDBC2)THEN
        WRITE(IUPRT,*)'THE NUMBER OF DIFFUSER-IN-LOOP TRACER INPUT
     .  IS NOT EQUAL TO NUMDBC2'
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      END IF
C
      NUMDBCTR = NUMDBCTR1 + NUMDBCTR2
      WRITE (IUPRT,8200) NUMDBCTR2,NUMDBCTR
 8200 FORMAT('NUMBER OF DIFFUSER IN LOOP TRACER LOADING IS',
     . I6,/'TOTAL NUMBER OF DIFFUSER TRACER LOADING IS',I6)
      IF (NUMDBCTR2.EQ.0) GOTO 3140
      IF (NUMDBCTR.GT.DBCM) THEN
        WRITE (IUPRT,64) NUMDBCTR, DBCM
        !!!&&&CALL SYSTEM('rm gcm_temp9*')
        STOP
      END IF
  64  FORMAT(//' NUMDBCTR (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' DBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
C
      OPEN(IUT99,FILE='gcm_temp99')
      DO M = 1, 100000
        READ (IUT401,77,ERR=3133) TIME
        WRITE (IUPRT,77) TIME
        OTIME=TIME
        READ (IUT401,77) (CDIFF1(N),N = 1,NUMDBCTR2)
        WRITE(IUPRT,77) (CDIFF1(N),N = 1,NUMDBCTR2)
C
        WRITE (IUT99,78) TIME
        WRITE (IUT99,78) (CDIFF1(N),N = 1,NUMDBCTR2)
C
      ENDDO     
 3133 BACKSPACE IUT401
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     .  'NEED MORE DIFFUSER IN LOOP BC TO COMPLETE THE SIMULATIONS'
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
C
 3140 CONTINUE
C
C
C*******************************************************************
C
C  SEDIMENT TRANSPORT 
C
      IF (SEDTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT702=702
          OPEN (IUT702,FILE='gcm_temp702')
C
          READ(IUT402,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT402,24)NUMDBCSE
          WRITE (IUPRT,24)NUMDBCSE
C
          IF (NUMDBCSE.EQ.0) GOTO 157
C
          IF (NUMDBCSE.GT.DBCM) THEN
           WRITE(IUPRT,5402)NUMDBCSE,DBCM
           !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
C
          DO 7529 N=1,NUMDBCSE
            READ(IUT402,51)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
 7529     CONTINUE
C
          DO 7560 M=1,100000
            READ(IUT402,77,END=7570) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT402,771)  (CDIFF(1,N),N=1,NUMDBCSE)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCSE)
C
C  Diffuser Loadings are in Kg/day. We need to compute sediment concentration
C  in g/cm3. Therefore Loading should be divided by a factor (86.4*1.0e6)
C
            DO 7580 NN=1,NUMDBCSE
              CDIFF(1,NN)=CDIFF(1,NN)/(86.4*1.0e06)
 7580       CONTINUE
C
            WRITE (IUT702,78)TIME
            WRITE (IUT702,78)  (CDIFF(1,N),N=1,NUMDBCSE)
C
 7560     CONTINUE
 7570     CONTINUE
        ENDIF

 157    IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT703=703
          OPEN (IUT703,FILE='gcm_temp703')
C
          READ(IUT403,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT403,24)NUMDBCSE
          WRITE (IUPRT,24)NUMDBCSE
C
          IF (NUMDBCSE.EQ.0) GOTO 158
C
          IF (NUMDBCSE.GT.DBCM) THEN
           WRITE(IUPRT,5402)NUMDBCSE,DBCM
           !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
C
          DO 7559 N=1,NUMDBCSE
            READ(IUT403,51)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
 7559     CONTINUE
C
C
          DO 7561 M=1,100000
            READ(IUT403,77,END=7571) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT403,771)  (CDIFF(1,N),N=1,NUMDBCSE)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCSE)
C
C
C  Diffuser Loadings are in Kg/day. We need to compute sediment concentration
C  in g/cm3. Therefore Loading should be divided by a factor (86.4*1.0e6)
C
            DO 7581 NN=1,NUMDBCSE
              CDIFF(1,NN)=CDIFF(1,NN)/(86.4*1.0e06)
 7581       CONTINUE
C
            WRITE (IUT703,78)TIME
            WRITE (IUT703,78)  (CDIFF(1,N),N=1,NUMDBCSE)
C
 7561     CONTINUE
 7571     CONTINUE
        ENDIF
      ENDIF

C*******************************************************************
C
C  CHEMICAL TRANSPORT
C
  158 IF (CHEMTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT704=704
          OPEN (IUT704,FILE='gcm_temp704')
C
          READ(IUT404,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT404,24)NUMDBCCH
          WRITE (IUPRT,24)NUMDBCCH
C
          IF (NUMDBCCH.EQ.0) GOTO 1588
          IF (NUMDBCCH.GT.DBCM) THEN
           WRITE(IUPRT,7402)NUMDBCCH,DBCM
           !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
C
          DO 7589 N=1,NUMDBCCH
            READ (IUT404,51)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
 7589     CONTINUE
C
C
          DO 7590 M=1,100000
            READ(IUT404,77,END=7600) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT404,771)  (CDIFF(1,N),N=1,NUMDBCCH)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCCH)
C
C  Diffuser Loadings are in Kg/day. We need to compute chemical concentration
C  in ug/cm3. Therefore Loading should be divided by a factor (86.4)
C
            DO 7610 NN=1,NUMDBCCH
              CDIFF(1,NN)=CDIFF(1,NN)/(86.4)
 7610       CONTINUE
C
            WRITE (IUT704,78)TIME
            WRITE (IUT704,78)  (CDIFF(1,N),N=1,NUMDBCCH)
C
 7590     CONTINUE
 7600     CONTINUE
        ENDIF
C
 1588   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT705=705
          OPEN (IUT705,FILE='gcm_temp705')
C
          READ(IUT405,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT405,24)NUMDBCCH
          WRITE (IUPRT,24)NUMDBCCH
C
          IF (NUMDBCCH.EQ.0) GOTO 1589
          IF (NUMDBCCH.GT.DBCM) THEN
           WRITE(IUPRT,7402)NUMDBCCH,DBCM
           !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
C
          DO 7599 N=1,NUMDBCCH
            READ (IUT405,51)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
 7599     CONTINUE
C
C
          DO 7620 M=1,100000
            READ(IUT405,77,END=7630) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT405,771)  (CDIFF(1,N),N=1,NUMDBCCH)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCCH)
C
C  Diffuser Loadings are in Kg/day. We need to compute chemical concentration
C  in ug/cm3. Therefore Loading should be divided by a factor (86.4)
C
            DO 7640 NN=1,NUMDBCCH
              CDIFF(1,NN)=CDIFF(1,NN)/86.4
 7640       CONTINUE
C
            WRITE (IUT705,78)TIME
            WRITE (IUT705,78)  (CDIFF(1,N),N=1,NUMDBCCH)
C
 7620     CONTINUE
 7630     CONTINUE
        ENDIF
      ENDIF

C
C*******************************************************************
C
C  POINT SOURCE LOADS:  CONSERVATIVE TRACER
C
 1589  IF (TRACER.EQ.'INCLUDE') THEN
        OPEN (IUT701,FILE='gcm_temp701')
C
        READ(IUT401,11)  (COM(I),I=1,80)
        WRITE(IUPRT,12) (COM(I),I=1,80)
C
        READ (IUT401,1120)NUMPSTR,OPTPSTR 
        WRITE (IUPRT,1120)NUMPSTR,OPTPSTR 
 1120   FORMAT (I5,6X,A4)
C
        IF (NUMPSTR.EQ.0) GOTO 3129
C
        DO 1121 N=1,NUMPSTR
          READ (IUT401,179)IPSTR(N),JPSTR(N),KPSTR(N)
          WRITE (IUPRT,179)IPSTR(N),JPSTR(N),KPSTR(N)
 1121   CONTINUE
C
        DO 1122 M=1,100000
          READ(IUT401,77,END=1123) TIME
          READ(IUT401,77) (CDIS1(N),N=1,NUMPSTR)
C
          WRITE(IUPRT,77) TIME
          OTIME=TIME
          WRITE (IUPRT,77) (CDIS1(N),N=1,NUMPSTR)
C
          WRITE (IUT701,78)TIME
          WRITE (IUT701,78) (CDIS1(N),N=1,NUMPSTR)
C
 1122   CONTINUE
C
 1123   CONTINUE 
       IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     . 'NEED MORE POINT SOURCE TRACER IN LOOP BC TO COMPLETE THE RUN' 
         !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
         STOP
       ENDIF
C
      ENDIF
C
3129  CONTINUE
C
C-----------------------------------------------------------------------
C-------- METEOROLOGICAL BOUNDARY CONDITIONS ---------------------------
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 1341  
C
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,131)  OPTMBC,ALAT,ALON,TR,WNDSH
      WRITE(IUPRT,131) OPTMBC,ALAT,ALON,TR,WNDSH
 131  FORMAT(A10,4F10.2)
C
      If (OPTMBC(1:8).NE.'AVERAGED'.AND.OPTMBC(1:8).NE.'LANDPFLX'.
     *AND.OPTMBC(1:8).NE.'AANDBFLX'.AND.OPTMBC(1:8).NE.'RANDMFLX'.
     *AND.OPTMBC(1:8).NE.'SYNOPTIC') Then
        Write (IUPRT,5401) OPTMBC
        !!!&&&CALL SYSTEM('rm gcm_temp9*')
        Stop
      End If
C
C-------- AVERAGED METEOROLOGICAL BOUNDARY CONDITIONS ------------------
C-------- PRECIPITATION(m/yr), EVAPORATION(m/yr) & HEAT FLUX(w/m2) -----
C-------- WIND FROM DIRECTION WIND IS BLOWING --------------------------
      IF(OPTMBC(1:8).EQ.'AVERAGED') THEN
       DO 161 M=1,100000
       READ(IURUN,77,END=165) TIME
       WRITE(IUPRT,77) TIME
       OTIME=TIME
       READ(IURUN,77)  WDS,WDD,HFLUX,QPREC,QEVAP
       WRITE(IUPRT,77)  WDS,WDD,HFLUX,QPREC,QEVAP
       QPREC=QPREC/(86400.*365.)
       QEVAP=QEVAP/(86400.*365.)
C-------- WIND STRESS --------------------------------------------------
       WDD=180.+WDD
       WDD=AMOD(WDD,360.)
       VWIND=WDS*COS(6.28319*WDD/360.)
       UWIND=WDS*SIN(6.28319*WDD/360.)
       CD=1.2E-3
       IF(WDS.GE.11.) CD=(0.49+0.065*WDS)*1.E-3
       IF(WDS.GE.25.) CD=(0.49+0.065*25.)*1.E-3
       TX=1.2*CD*UWIND*WDS
       TY=1.2*CD*VWIND*WDS
       WRITE(IUT93,78) TIME
161    WRITE(IUT93,78) HFLUX,TX,TY,UWIND,VWIND,WDS,WDD,QPREC,QEVAP
165    CONTINUE
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE MET DATA TO COMPLETE THE RUN' 
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
C
      ELSE IF(OPTMBC(1:8).EQ.'SYNOPTIC') THEN
C
C-------- SYNOPTIC METEOROLOGICAL BOUNDARY CONDITIONS ------------------
C-------- PRECIPITATION(m/yr), EVAPORATION(m/yr) & HEAT FLUX(w/m2) -----
*
       OPEN(IUT191,FILE='synop_hflx',FORM='unformatted',STATUS='old')
       OPEN(IUT192,FILE='gcm_temp192',FORM='unformatted')
*
       DO 171 M=1,100000
         READ(IUT191,END=175) TIME
         READ(IUT191) ((SHFLX(I,J),I=1,IM),J=1,JM)
         OTIME=TIME
         WRITE(IUT192)TIME
 171     WRITE(IUT192)((SHFLX(I,J),I=1,IM),J=1,JM)
 175   CONTINUE

C-------- CALCULATE WIND STRESS ----------------------------------------
C-------- FROM U AND V WIND COMPONENTS ---------------------------------
       OPEN(IUUAV,FILE='synop_wind',FORM='unformatted')
       OPEN(IUT95,FILE='gcm_temp95',FORM='unformatted')
       IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
         OPEN(IUT193,FILE='gcm_temp193',FORM='unformatted')
       ENDIF
C
       DO 181 M=1,100000
       READ(IUUAV,END=185) TIME
       READ(IUUAV) ((TX2D(I,J),TY2D(I,J),PATM(I,J),I=1,IM),J=1,JM)

C      WIND SPEED for BOTH DIR is USED BY WAVEDYN MODELS
       IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
         WRITE(IUT193) TIME
         WRITE(IUT193) ((TX2D(I,J),TY2D(I,J),I=1,IM),J=1,JM)
       ENDIF
C
       DO 182 J=1,JM
       DO 182 I=1,IM
       WDS=SQRT(TX2D(I,J)**2+TY2D(I,J)**2)
       CD=1.2E-3
       IF(WDS.GE.11.) CD=(0.49+0.065*WDS)*1.E-3
       IF(WDS.GE.25.) CD=(0.49+0.065*25.)*1.E-3
       TX2D(I,J)=1.2*CD*TX2D(I,J)*WDS
       TY2D(I,J)=1.2*CD*TY2D(I,J)*WDS
       PATM(I,J) = PATM(I,J)*1.0E02  ! from mbar to N/m2
182    CONTINUE
C      WIND STRESS for BOTH DIR is USED TO APPLY on WATER SURFACE
       WRITE(IUT95) TIME
       WRITE(IUT95) ((TX2D(I,J),TY2D(I,J),PATM(I,J),I=1,IM),J=1,JM)
C
 181    CONTINUE
C
185    CONTINUE
       CLOSE(IUUAV)
C------ HEATFLUX CALCULATION --------------------------------------
CSTART, MODIFIED FOR HEATFLUX CALCULATION 
C       ELSE IF (OPTMBC(1:8).EQ.'HEATFLUX') THEN
C
       else If (OPTMBC(1:8).EQ.'LANDPFLX'.OR.
     *          OPTMBC(1:8).EQ.'AANDBFLX'.OR.
     *          OPTMBC(1:8).EQ.'RANDMFLX') Then
        IF(ALAT.LE.-999.0) THEN
C
C  For fraction of SWRAD, TR (0.-1.),
C  Absorbed in surface layer of water colum, rest (1-TR) is penetrated
C
c        READ (IURUN,77) ALAT,ALON, TR, WNDSH ! ALON IS NEGATIVE FOR W
         WRITE(IUPRT,7734) ALAT,ALON, TR, WNDSH
7734  FORMAT(/'Latitude= ', f10.4,'     Longitude= ', f10.4/
     +     5x,'Fraction of SWRad Absorbed in Surface = ', F10.3/
     +     5x,'Wind Sheltering Coeffient Used =', F10.3)
         ALON=-ALON  ! TO FOLLOW THE CONVENTION IN NCLD.F
	ENDIF
C
        AIRLOW=-999.9
        AIRHIGH=999.9
        WRITE(IUPRT,7735) AIRLOW,AIRHIGH
7735  FORMAT(' THE LOWEST AIR TEMPERATURE ALLOWED: ',f10.2/,
     &       ' THE HIGHEST AIR TEMPERATURE ALLOWED: ',f10.2)
        DO 170 M = 1, 100000
          READ (IURUN,77,End=180) TIME
          WRITE (IUPRT,77) TIME
          OTIME=TIME
C
          READ (IURUN,7750) WDS,WDD,SWOBS,AIRTMP,RELHUM,BAROP,
     &      CLD,EXTC,QPREC,QEVAP
c---       
	IF(AIRTMP.LT.AIRLOW) AIRTMP=AIRLOW
	IF(AIRTMP.GT.AIRHIGH) AIRTMP=AIRHIGH
          WRITE (IUPRT,7750) AIRTMP, RELHUM, BAROP,
     &      SWOBS, WDS, WDD, CLD, EXTC, QPREC, QEVAP
C CHANGE EVAP AND PRECIP UNITS FROM m/year TO m/s.
          QPREC = QPREC / (86400.*365.)
          QEVAP = QEVAP / (86400.*365.)
C-------- WIND STRESS (large & pond, newtons/m2)
          WDD = 180. + WDD
          WDD = AMOD(WDD,360.)
          VWIND = WDS * COS(6.28319*WDD/360.)
          UWIND = WDS * SIN(6.28319*WDD/360.)
          CD = 1.2E-3
          If (WDS.GE.11.) CD = (0.49+0.065*WDS) * 1.E-3
          If (WDS.GE.25.) CD = (0.49+0.065*25.) * 1.E-3
          TX = 1.2 * CD * UWIND * WDS
          TY = 1.2 * CD * VWIND * WDS
          WRITE (IUT93,78) TIME
          WRITE (IUT93,7850) AIRTMP, RELHUM, BAROP,
     &     SWOBS, TX, TY, UWIND, VWIND, CLD, EXTC, QPREC, QEVAP
C
  170   CONTINUE
  180   CONTINUE
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE MET DATA TO COMPLETE THE RUN' 
        !!!&&&!!!&&&CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
      ENDIF
C
 77    FORMAT(8F10.1)
 771   FORMAT(5E14.7)
 7750  FORMAT(10F10.5)
 78    FORMAT(8E14.7)
 7850  FORMAT(12E14.7)
C
C  FOR INPUTTING WAM MODEL OUTPUT FOR WIND WAVE CALCULATIONS
C
1341      IF (WAVEDYN.EQ.'EXTERNAL') THEN
C
        OPEN (UNIT=110,FILE='wave_input',FORM='UNFORMATTED')
        OPEN (UNIT=111,FILE='gcm_temp111',FORM='UNFORMATTED')
C
        DO 8010 M=1,1000000
          READ (110,END=8100)TIME
          READ (110) WTEMP1   ! WAVE HEIGHT
          READ (110) WTEMP2   ! WAVE PERIOD
          READ (110) WTEMP3   ! WAVE DIRECTION
C
          WRITE (111)TIME
          WRITE (111) ((WTEMP1(I,J),I=1,IM),J=1,JM)
          WRITE (111) ((WTEMP2(I,J),I=1,IM),J=1,JM)
          WRITE (111) ((WTEMP3(I,J),I=1,IM),J=1,JM)
 8010   CONTINUE 
 8100   CLOSE (110)
      ENDIF
C
 5401 Format (//'OPTMBC HAS BEEN SELECTED AS ',A20//'THIS OPTION IS
     * NOT CORRECT, PLEASE FIX AND RESUBMIT'//)
 5402 Format (//'NUMDBCSE IS EQUAL TO ',I20//'THIS VALUE IS
     * NOT LESS OR EQUAL TO DBCM in comdeck',I20,//,
     * 'PLEASE FIX AND RESUBMIT'//)
 7402 Format (//'NUMDBCCH IS GIVEN EQUAL TO ',I20//'THIS VALUE IS
     * NOT LESS OR EQUAL TO DBCM in comdeck',I20,//,
     * 'PLEASE FIX AND RESUBMIT'//)
      RETURN
      END
