      SUBROUTINE FIRST
C
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
* Point of Contact: ecomsed-support@hydroquap.com                        *
C*************************************************************************
      INCLUDE 'comdeck'
      SAVE
C
      IF(NUMEBC.EQ.0) GO TO 140
       REWIND IUT90
       DO 100 M=1,100000
       IFLAG=1
       READ (IUT90,10,ERR=130) T2E
       READ (IUT90,10) (DEBDRY(N,2),N=1,NUMEBC)
       IF(THOUR.LT.T2E) GO TO 110
       T1E=T2E 
       DO 120 N=1,NUMEBC
 120   DEBDRY(N,1)=DEBDRY(N,2) 
 100   CONTINUE
 110   CONTINUE
 140  CONTINUE
C 
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 1135
C
      IF (NUMEBC.EQ.0) GOTO 135
      REWIND IUT94
C
      DO 101 M=1,100000
       IFLAG=2
      READ (IUT94,10,ERR=130) T2TS
      DO 136 N=1,NUMEBC
      READ (IUT94,10) (DTBDRY(N,K,2),K=1,KBM1)
      READ (IUT94,10) (DSBDRY(N,K,2),K=1,KBM1)
136   CONTINUE
      IF(THOUR.LT.T2TS) GOTO 111
      T1TS=T2TS
      DO 121 N=1,NUMEBC
      DO 121 K=1,KBM1
      DTBDRY(N,K,1)=DTBDRY(N,K,2)
      DSBDRY(N,K,1)=DSBDRY(N,K,2)
 121  CONTINUE
 101  CONTINUE
 111  CONTINUE
C
C  DISSOLVED TRACER AT OPEN BOUNDARY
C
 1135 IF (TRACER.EQ.'INCLUDE'.AND.NUMEBCTR.GT.0) THEN
        REWIND (IUT501)
        DO 1010 M=1,100000
          IFLAG=3
          READ (IUT501,10,ERR=130) T2CONOB
C
          DO 1020 N=1,NUMEBCTR
            READ (IUT501,10) (DCBDRY1(N,K,2),K=1,KBM1)
 1020     CONTINUE
C
          IF(THOUR.LT.T2CONOB) GOTO 1030
C
          T1CONOB=T2CONOB
          DO 1040 N=1,NUMEBCTR
            DO 1040 K=1,KBM1
              DCBDRY1(N,K,1)=DCBDRY1(N,K,2)
 1040     CONTINUE
 1010     CONTINUE
 1030   CONTINUE
      ENDIF
C
C  SEDIMENT TRANSPORT AT OPEN BOUNDARY
C
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMEBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT502)
          DO 1110 M=1,100000
            IFLAG=4
            READ (IUT502,10,ERR=130) T2SEDOBM
C
            DO 1120 N=1,NUMEBCSE
              READ (IUT502,10) (DCBDRY(1,N,K,2),K=1,KBM1)
 1120       CONTINUE
C
            IF(THOUR.LT.T2SEDOBM) GOTO 1130
C
            T1SEDOBM=T2SEDOBM
            DO 1140 N=1,NUMEBCSE
              DO 1140 K=1,KBM1
                DCBDRY(1,N,K,1)=DCBDRY(1,N,K,2)
 1140       CONTINUE
 1110       CONTINUE
 1130     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT503)
          DO 1150 M=1,100000
            IFLAG=5
            READ (IUT503,10,ERR=130) T2SEDOBS
C
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
C
            DO 1160 N=1,NUMEBCSE
              READ (IUT503,10) (DCBDRY(KK,N,K,2),K=1,KBM1)
 1160       CONTINUE
C
            IF(THOUR.LT.T2SEDOBS) GOTO 1180
C
            T1SEDOBS=T2SEDOBS
            DO 1170 N=1,NUMEBCSE
              DO 1170 K=1,KBM1
                DCBDRY(KK,N,K,1)=DCBDRY(KK,N,K,2)
 1170       CONTINUE
 1150       CONTINUE
 1180     CONTINUE
        ENDIF
      ENDIF
C
C  PARTICLE-BOUND TRACER AT OPEN BOUNDARY
C
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMEBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT504)
          DO 1210 M=1,100000
            IFLAG=6
            READ (IUT504,10,ERR=130) T2CHMOBM
C
            KK=1
C
            DO 1220 N=1,NUMEBCCH
              READ (IUT504,10) (DPBDRY(KK,N,K,2),K=1,KBM1)
 1220       CONTINUE
C
            IF(THOUR.LT.T2CHMOBM) GOTO 1230
C
            T1CHMOBM=T2CHMOBM
            DO 1240 N=1,NUMEBCCH
              DO 1240 K=1,KBM1
                DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 1240       CONTINUE
 1210     CONTINUE
 1230     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT505)
          DO 1250 M=1,100000
            IFLAG=7
            READ (IUT505,10,ERR=130) T2CHMOBS
C
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
C
            DO 1260 N=1,NUMEBCCH
              READ (IUT505,10) (DPBDRY(KK,N,K,2),K=1,KBM1)
 1260       CONTINUE
C
            IF(THOUR.LT.T2CHMOBS) GOTO 1270
C
            T1CHMOBS=T2CHMOBS
            DO 1280 N=1,NUMEBCCH
              DO 1280 K=1,KBM1
                DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 1280       CONTINUE
 1250     CONTINUE
 1270     CONTINUE
        ENDIF
      ENDIF
C
 135  CONTINUE
C
C  RIVER DISCHARGE BOUNDARY
C
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 210
C
      IF (NUMQBC.EQ.0) GOTO 210
C
      REWIND (IUT91)
C
      DO 200 M=1,100000
            IFLAG=8
      READ (IUT91,10,ERR=130) T2Q
      READ (IUT91,10) (DQDIS(N,2),N=1,NUMQBC)
      READ (IUT91,10) (DTDIS(N,2),N=1,NUMQBC)
      READ (IUT91,10) (DSDIS(N,2),N=1,NUMQBC)
C
      IF(THOUR.LT.T2Q) GO TO 210
      T1Q=T2Q 
      DO 220 N=1,NUMQBC
      DQDIS(N,1)=DQDIS(N,2) 
      DTDIS(N,1)=DTDIS(N,2) 
      DSDIS(N,1)=DSDIS(N,2) 
 220  CONTINUE
 200  CONTINUE
 210  CONTINUE
C
C  DISSOLVED TRACER AT RIVER DISCHARGE BOUNDARY
C
      IF (TRACER.EQ.'INCLUDE'.AND.NUMQBCTR.GT.0) THEN
        REWIND (IUT601)
        DO 1275 M=1,100000
          IFLAG=9
          READ (IUT601,10,ERR=130) T2CON
          READ (IUT601,10) (DCDIS1(N,2),N=1,NUMQBCTR)
C
          IF(THOUR.LT.T2CON) GO TO 1285
          T1CON=T2CON 
C
          DO 1295 N=1,NUMQBCTR
            DCDIS1(N,1)=DCDIS1(N,2) 
 1295     CONTINUE
 1275   CONTINUE
 1285   CONTINUE
      ENDIF
C
C  SEDIMENT TRANSPORT AT RIVER DISCHARGE BOUNDARY
C
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMQBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT602)
          DO 1300 M=1,100000
            IFLAG=10
            READ (IUT602,10,ERR=130) T2SEDM
            KK=1
            READ (IUT602,10) (DCDIS(KK,N,2),N=1,NUMQBCSE)
            IF(THOUR.LT.T2SEDM) GO TO 1310
            T1SEDM=T2SEDM 
C
            DO 1320 N=1,NUMQBCSE
              DCDIS(KK,N,1)=DCDIS(KK,N,2) 
 1320       CONTINUE
 1300     CONTINUE
 1310     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT603)
          DO 1330 M=1,100000
            IFLAG=11
            READ (IUT603,10,ERR=130) T2SEDS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT603,10) (DCDIS(KK,N,2),N=1,NUMQBCSE)
            IF(THOUR.LT.T2SEDS) GO TO 1340
            T1SEDS=T2SEDS 
C
            DO 1350 N=1,NUMQBCSE
              DCDIS(KK,N,1)=DCDIS(KK,N,2) 
 1350       CONTINUE
 1330     CONTINUE
 1340     CONTINUE
        ENDIF
      ENDIF
C
C  PARTICLE-BOUND TRACER AT RIVER DISCHARGE BOUNDARY
C
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMQBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT604)
          DO 1400 M=1,100000
            IFLAG=12
            READ (IUT604,10,ERR=130) T2CHMM
            KK=1
            READ (IUT604,10) (DPDIS(KK,N,2),N=1,NUMQBCCH)
            IF(THOUR.LT.T2CHMM) GO TO 1410
            T1CHMM=T2CHMM 
C
            DO 1420 N=1,NUMQBCCH
              DPDIS(KK,N,1)=DPDIS(KK,N,2) 
 1420       CONTINUE
 1400     CONTINUE
 1410     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT605)
          DO 1430 M=1,100000
            IFLAG=13
            READ (IUT605,10,ERR=130) T2CHMS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT605,10) (DPDIS(KK,N,2),N=1,NUMQBCCH)
            IF(THOUR.LT.T2CHMS) GO TO 1440
            T1CHMS=T2CHMS 
C
            DO 1450 N=1,NUMQBCCH
              DPDIS(KK,N,1)=DPDIS(KK,N,2) 
 1450       CONTINUE
 1430     CONTINUE
 1440     CONTINUE
        ENDIF
      ENDIF
C 
C  ************************************************************
C  DIFFUSER BOUNDARY
c
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 330
      IF(NUMDBC1.EQ.0) GO TO 1510
      REWIND IUT92
C
      DO M=1,100000
        IFLAG=14
        READ (IUT92,10,ERR=130) T2D
        READ (IUT92,10) (DQDIFF(N,2),N=1,NUMDBC1)
        READ (IUT92,10) (DTDIFF(N,2),N=1,NUMDBC1)
        READ (IUT92,10) (DSDIFF(N,2),N=1,NUMDBC1)
        IF(THOUR.LT.T2D) GO TO 330
        T1D=T2D
        DO N=1,NUMDBC1
          DQDIFF(N,1)=DQDIFF(N,2)
          DTDIFF(N,1)=DTDIFF(N,2)
          DSDIFF(N,1)=DSDIFF(N,2)
        ENDDO 
      ENDDO    
C
 330  CONTINUE
C
C  DISSOLVED TRACER AT DIFFUSER
C
      IF (TRACER.NE.'INCLUDE'.OR.NUMDBCTR1.EQ.0) GO TO 1510
      REWIND IUT98
C
      DO M=1,100000
        IFLAG=15
        READ (IUT98,10,ERR=130) T2DCON
        READ (IUT98,10) (DCDIFF1(N,2),N=1,NUMDBCTR1)
C
        IF(THOUR.LT.T2DCON) GO TO 1510
C
        T1DCON=T2DCON
C
        DO N=1,NUMDBCTR1
          DCDIFF1(N,1)=DCDIFF1(N,2)
        ENDDO     
      ENDDO     
 1510 CONTINUE
c
C  DIFFUSER IN LOOP BOUNDARY
c
      IF(NUMDBC2.EQ.0) GO TO 1511
      REWIND IUT96
C
C
      DO M=1,100000
        IFLAG=144
        READ (IUT96,10,ERR=130) T2D2
        READ (IUT96,10) (DQDIFF(N,2),N=NUMDBC1+1,NUMDBC)
        READ (IUT96,10) (DTDIFF(N,2),N=NUMDBC1+1,NUMDBC)
        READ (IUT96,10) (DSDIFF(N,2),N=NUMDBC1+1,NUMDBC)
        IF(THOUR.LT.T2D2) GO TO 311
        T1D2=T2D2
        DO N=NUMDBC1+1,NUMDBC
          DQDIFF(N,1)=DQDIFF(N,2)
          DTDIFF(N,1)=DTDIFF(N,2)
          DSDIFF(N,1)=DSDIFF(N,2)
        ENDDO 
      ENDDO   
 311  CONTINUE
C
C  DISSOLVED TRACER AT DIFFUSER IN LOOPS
      IF (TRACER.NE.'INCLUDE'.OR.NUMDBCTR2.EQ.0) GO TO 1511
      REWIND IUT99
      DO M=1,100000
        IFLAG=150
        READ (IUT99,10,ERR=130) T2D2CON
        READ (IUT99,10) (DCDIFF1(N,2),N=NUMDBCTR1+1,NUMDBCTR)
C
        IF(THOUR.LT.T2D2CON) GO TO 1511
C
        T1D2CON=T2D2CON
C
        DO N=NUMDBCTR1+1,NUMDBCTR
          DCDIFF1(N,1)=DCDIFF1(N,2)
        ENDDO   
      ENDDO   
 1511 CONTINUE
C  POINT SOURCE LOADS:  DISSOLVED TRACER
C
      IF (TRACER.EQ.'INCLUDE'.AND.NUMPSTR.GT.0) THEN
        REWIND (IUT701)
        DO 4275 M=1,100000
          IFLAG=91
          READ (IUT701,10,ERR=130) T2PSTR
          READ (IUT701,10) (DCPSTR(N,2),N=1,NUMPSTR)
C
          IF(THOUR.LT.T2PSTR) GO TO 4285
C
          T1PSTR=T2PSTR 
C
          DO 4295 N=1,NUMPSTR
            DCPSTR(N,1)=DCPSTR(N,2) 
 4295     CONTINUE
 4275   CONTINUE
 4285   CONTINUE
      ENDIF
C
C  SEDIMENT LOADS AT DIFFUSER 
C
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMDBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT702)
          DO 7300 M=1,100000
            IFLAG=702
            READ (IUT702,10,ERR=730) T2SEDDM
            KK=1
            READ (IUT702,10) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
            IF(THOUR.LT.T2SEDDM) GO TO 7310
            T1SEDDM=T2SEDDM
C
            DO 7320 N=1,NUMDBCSE
              DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7320       CONTINUE
 7300     CONTINUE
 7310     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT703)
          DO 7330 M=1,100000
            IFLAG=703
            READ (IUT703,10,ERR=730) T2SEDDS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT703,10) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
            IF(THOUR.LT.T2SEDDS) GO TO 7340
            T1SEDDS=T2SEDDS
C
            DO 7350 N=1,NUMDBCSE
              DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7350       CONTINUE
 7330     CONTINUE
 7340     CONTINUE
        ENDIF
      ENDIF
C
C  PARTICLE-BOUND TRACER AT DIFFUSER BOUNDARY
C
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMDBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT704)
          DO 7400 M=1,100000
            IFLAG=704
            READ (IUT704,10,ERR=130) T2CHMDM
            KK=1
            READ (IUT704,10) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
            IF(THOUR.LT.T2CHMDM) GO TO 7410
            T1CHMDM=T2CHMDM
C
            DO 7420 N=1,NUMDBCCH
              DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7420       CONTINUE
 7400     CONTINUE
 7410     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT705)
          DO 7430 M=1,100000
            IFLAG=705
            READ (IUT705,10,ERR=130) T2CHMDS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT705,10) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
            IF(THOUR.LT.T2CHMDS) GO TO 7440
            T1CHMDS=T2CHMDS
C
            DO 7450 N=1,NUMDBCCH
              DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7450       CONTINUE
 7430     CONTINUE
 7440     CONTINUE
        ENDIF
      ENDIF
C
C 
      IF(OPTMBC(1:8).EQ.'AVERAGED') THEN
       REWIND IUT93
       DO 400 M=1,100000
            IFLAG=93
       READ (IUT93,10,ERR=130) T2M
       READ (IUT93,10) DHFLUX(2),DTX(2),DTY(2),
     +   DWSX(2),DWSY(2),DWDS(2),DWDD(2),DQPREC(2),DQEVAP(2)
       IF(THOUR.LT.T2M) GO TO 410
       T1M=T2M
       DQPREC(1)=DQPREC(2) 
       DQEVAP(1)=DQEVAP(2) 
       DHFLUX(1)=DHFLUX(2)
       DTX   (1)=DTX   (2) 
       DTY   (1)=DTY   (2) 
       DWSX  (1)=DWSX  (2) 
       DWSY  (1)=DWSY  (2) 
       DWDS  (1)=DWDS  (2) 
       DWDD  (1)=DWDD  (2) 
 400   CONTINUE
 410   CONTINUE
      ELSE IF(OPTMBC(1:8).EQ.'SYNOPTIC') THEN  
       REWIND IUT192
       DO 500 M=1,100000
            IFLAG=92
       READ (IUT192,ERR=130) T2M
       READ (IUT192) ((DHFLX2D(I,J,2),I=1,IM),J=1,JM)
       IF(THOUR.LT.T2M) GO TO 510
       T1M=T2M
       DO 511 J=1,JM     
       DO 511 I=1,IM
        DHFLX2D(I,J,1)=DHFLX2D(I,J,2)
 511   CONTINUE
 500   CONTINUE
 510   CONTINUE

C------  WIND STRESS  ----------------------------------------
       REWIND IUT95
       DO 520 M=1,100000
            IFLAG=95
       READ (IUT95,ERR=130) T2W
       READ (IUT95) ((DTX2D(I,J,2),DTY2D(I,J,2),DPATM(I,J,2),
     +                I=1,IM),J=1,JM)
       IF(THOUR.LT.T2W) GO TO 530
       T1W=T2W
       DO 540 J=1,JM
       DO 540 I=1,IM
       DPATM(I,J,1)=DPATM(I,J,2) 
       DTX2D(I,J,1)=DTX2D(I,J,2) 
 540   DTY2D(I,J,1)=DTY2D(I,J,2) 
 520   CONTINUE
 530   CONTINUE

C------  WIND SPEED FOR WAVEDYN MODEL 
       IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
         REWIND IUT193
         DO 521 M=1,100000
            IFLAG=193
         READ (IUT193,ERR=130) T2WV
         READ (IUT193) ((WU2(I,J),WV2(I,J),I=1,IM),J=1,JM)
         IF(THOUR.LT.T2WV) GO TO 531
         T1WV=T2WV
         DO 541 J=1,JM
         DO 541 I=1,IM
         WU1(I,J)=WU2(I,J) 
 541     WV1(I,J)=WV2(I,J) 
 521     CONTINUE
 531     CONTINUE
       ENDIF

C
       else If (OPTMBC(1:8).EQ.'LANDPFLX'.OR.
     *          OPTMBC(1:8).EQ.'RANDMFLX'.OR.
     *          OPTMBC(1:8).EQ.'AANDBFLX') Then
        REWIND IUT93
        DO 150 M = 1, 100000
          READ (IUT93,810,ERR=130) T2M
          READ (IUT93,810) DAIRTM(2), DRELHU(2),
     *        DBAROP(2), DSWOBS(2), DTX(2), DTY(2), DWSX(2), DWSY(2),
     *     CLOUD(2), EXTCOEF(2), DQPREC(2), DQEVAP(2)
810     format(12e14.7)
C
          IF (THOUR.LT.T2M) GO TO 160
          T1M = T2M
          DQPREC(1) = DQPREC(2)
          DQEVAP(1) = DQEVAP(2)
          DAIRTM(1) = DAIRTM(2)
          DRELHU(1) = DRELHU(2)
          DBAROP(1) = DBAROP(2)
          DSWOBS(1) = DSWOBS(2)
          CLOUD(1)  = CLOUD(2)
          EXTCOEF(1)  = EXTCOEF(2)
          DTX(1) = DTX(2)
          DTY(1) = DTY(2)
          DWSX(1) = DWSX(2)
          DWSY(1) = DWSY(2)
          DWDS  (1)=DWDS  (2) 
          DWDD  (1)=DWDD  (2) 
  150   CONTINUE
  160   CONTINUE
      ENDIF
C
C  FOR WIND WAVE INPUT
C
      IF (WAVEDYN.EQ.'EXTERNAL') THEN
        REWIND (111)
        DO 8010 M=1,100000
          READ (111,ERR=8020) T2WAVE
          READ (111) ((DWHT(I,J,2),I=1,IM),J=1,JM)
          READ (111) ((DWPER(I,J,2),I=1,IM),J=1,JM)
          READ (111) ((DWDIR(I,J,2),I=1,IM),J=1,JM)
C
          IF (THOUR.LT.T2WAVE) GO TO 8020
C
          T1WAVE=T2WAVE
C
          DO 8030 J=1,JM
            DO 8030 I=1,IM
              DWHT(I,J,1)=DWHT(I,J,2) 
              DWPER(I,J,1)=DWPER(I,J,2) 
              DWDIR(I,J,1)=DWDIR(I,J,2) 
 8030     CONTINUE
 8010   CONTINUE
 8020   CONTINUE
      ENDIF
C
      RETURN
C
 130  WRITE(6,20) IFLAG
 20   FORMAT(//' THERE IS INSUFFICIENT TEMPORAL DATA FOR THIS RUN'/,
     .         '        REVISE INPUT DECK AND RESUBMIT '
     +     // 'IFLAG =',I4)
C
 730  WRITE(6,70) IFLAG
 70   FORMAT(//' THERE IS INSUFFICIENT TEMPORAL DATA FOR THIS RUN'/,
     .         '        REVISE INPUT DECK AND RESUBMIT '
     +     // 'IFLAG =',I4)
 10   FORMAT(8E14.7)
      STOP
      END
