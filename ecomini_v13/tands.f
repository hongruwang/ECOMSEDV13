      SUBROUTINE TANDS
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
      INCLUDE 'comdeck'
      SAVE
C
      DIMENSION ZM(KBM1),TI(KB),SI(KB)
      DIMENSION TA(KSL),SA(KSL),TS(IM,JM,KSL),SS(IM,JM,KSL)
C
C  CHECK THAT NO DEPTH IS DEEPER THAN BOTTOM STANDARD LEVEL
C
      DO 9000 J=2,JMM1
        DO 9000 I=2,IMM1
          IF (H(I,J).GT.0.0) THEN
            IF (H(I,J).GT.-DPTHSL(KSL)) THEN
              WRITE (IUPRT,9001)I,J,H(I,J),DPTHSL(KSL)
 9001         FORMAT (/5X,'****** PROGRAM EXECUTION STOPPED *****',
     +  /7X,'SPECIFIED DEPTH IN model_grid AT (i,j) =',2I4,
     +  /7X,'IS ',F10.2,3X,'WHICH IS DEEPER THAN THE LOWEST STANDARD',
     +  /7X,'LEVEL, WHICH IS ',F10.2,
     +  /7X,'PLEASE CORRECT AND RESUBMIT')
C
             STOP
           ENDIF
         ENDIF
 9000 CONTINUE
C
C-------- SPREAD STANDARD LEVEL TEMPERATURE AND SALINITY ---------------
      IF(OPTTSI(1:5).EQ.'FIXED') THEN
       DO 120 K=1,KSL
       DO 120 J=1,JM
       DO 120 I=1,IM
       TS(I,J,K)=TSI(K)*FSM(I,J)
       SS(I,J,K)=SSI(K)*FSM(I,J)
  120  CONTINUE
C
      ELSE
C-------- INPUT TEMPERATURE AND SALINITY DATA FROM FILE ----------------
       OPEN (IUTAS,FILE='init_tands',FORM='formatted',status='OLD')
       DO 130 N=1,(IM*JM)
       READ(IUTAS,132,END=1310)
     .            I,J,(TS(I,J,K),K=1,KSL),(SS(I,J,K),K=1,KSL)
  130  CONTINUE
  131  CONTINUE
       CLOSE(IUTAS)
  132  FORMAT (2I5,100F5.0)
      ENDIF
C
      DO 220 K=1,KSL
      TA(K)=0.0
      SA(K)=0.0
C
      COUNT=0.0
      DO 210 J=1,JM
      DO 210 I=1,IM
C
            IF (K.EQ.1 )THEN
             IF(H(I,J).GT.0.0)THEN
              COUNT = COUNT + FSM(I,J)*ART(I,J)
              TA(K) = TA(K) + ART(I,J)*TS(I,J,K) * FSM(I,J)
              SA(K) = SA(K) + ART(I,J)*SS(I,J,K) * FSM(I,J)
             ENDIF
            ELSE
             IF(H(I,J).GT.-DPTHSL(K-1)) THEN 
              COUNT = COUNT + FSM(I,J)*ART(I,J)
              TA(K) = TA(K) + ART(I,J)*TS(I,J,K) * FSM(I,J)
              SA(K) = SA(K) + ART(I,J)*SS(I,J,K) * FSM(I,J)
             ENDIF
            ENDIF
C
  210 CONTINUE
      IF (COUNT.LT.1.0) GO TO 220
      TA(K)=TA(K)/COUNT
      SA(K)=SA(K)/COUNT
C
  220 CONTINUE
      DO 250 J=1,JM
      DO 250 I=1,IM
      IF (FSM(I,J).EQ.0.0) GO TO 250
      DO 230 K=1,KBM1
      ZM(K)=ZZ(K)*H(I,J)
  230 CONTINUE
      CALL SINTER (DPTHSL,TA,ZM,TI,KSL,KBM1)
      CALL SINTER (DPTHSL,SA,ZM,SI,KSL,KBM1)
C
      DO 240 K=1,KBM1
      T(I,J,K)=TI(K)
      S(I,J,K)=SI(K)
C
  240 CONTINUE
  250 CONTINUE
C
      CALL DENS
C
      DO 280 K=1,KBM1
      DO 280 J=1,JM
      DO 280 I=1,IM
      TMEAN(I,J,K)=T(I,J,K)*FSM(I,J)
      SMEAN(I,J,K)=S(I,J,K)*FSM(I,J)
      RMEAN(I,J,K)=RHO(I,J,K)*FSM(I,J)
  280 CONTINUE
C
      DO 350 J=1,JM
      DO 350 I=1,IM
      IF (FSM(I,J).EQ.0.0) GO TO 350
      DO 320 K=1,KSL
      TA(K)=TS(I,J,K)
      SA(K)=SS(I,J,K)
C
  320 CONTINUE
      DO 330 K=1,KBM1
      ZM(K)=ZZ(K)*H(I,J)
  330 CONTINUE
      CALL SINTER (DPTHSL,TA,ZM,TI,KSL,KBM1)
      CALL SINTER (DPTHSL,SA,ZM,SI,KSL,KBM1)
C
      DO 340 K=1,KBM1
      T(I,J,K)=TI(K)*FSM(I,J)
      S(I,J,K)=SI(K)*FSM(I,J)
C
  340 CONTINUE
  350 CONTINUE
      CALL DENS
C
      DO 430 K=1,KBM1
      DO 430 J=1,JM
      DO 430 I=1,IM
      TB(I,J,K)=T(I,J,K)
      SB(I,J,K)=S(I,J,K)
C
  430 CONTINUE
C
      DO 440 K=1,KB
      DO 440 J=1,JM
      DO 440 I=1,IM
      A(I,J,K)=0.0
      C(I,J,K)=0.0
      VH(I,J,K)=0.0
      VHP(I,J,K)=0.0
      PROD(I,J,K)=0.0
      DTEF(I,J,K)=0.0
  440 CONTINUE
C
      DO 450 K=1,KSL
      DO 450 J=1,JM
      DO 450 I=1,IM
      TS(I,J,K)=0.0
      SS(I,J,K)=0.0
C
  450 CONTINUE
C
      RETURN
1310	WRITE (iuprt,*) 'check inti_tands file and resubmit the run'
	stop
      END
