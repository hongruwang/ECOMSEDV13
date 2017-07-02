      SUBROUTINE BCOND(IDX,DTI2,KS)
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
c       UPWIND FLUX AT OPEN SALT/TEMP BOUNDARY
C
      INCLUDE 'comdeck'
      PARAMETER (INTMAX = 100) 
      REAL SSTINT(INTMAX), WT(INTMAX), WQ(INTMAX), RLWC(INTMAX) 
      REAL TH2D(IM,JM) 
C
      SAVE
C
      Data SBC /2.0411E-7/     ! 4*Boltzman Constant 4*5.67E-8 Army
      Data ZERO /0.0/  
c     Following two lines are for ROSATI and Miyakoda (RANDMFLX)
      Data DHR /1.0/           ! Period of Solar radiation Average in Hr
      Data REFL /0.1/          ! REFL - FRACTION SWOBS THAT REFLECTS BACK
C
C
C===============>    S2     M2     N2     K1     P1     O1    
      PI2=6.283185307
      PERIOD(1)=43200.
      PERIOD(2)=44712.
      PERIOD(3)=45570.
      PERIOD(4)=86164.
      PERIOD(5)=86637.
      PERIOD(6)=92950.
C
C-----------------------------------------------------------------------
C        SPECIFICATION OF OPEN BOUNDARY CONDITIONS
C   NOTE THAT AT J=2 U CALCULATION EXCLUDES ADVECTION AND DIFFUSION
C-----------------------------------------------------------------------
C        I-1           I         I+1
C
C       U(JM)=      EL(JM)      U(JM)
C
C                    V(JM)
C
C     *U(JMM1)*   *EL(JMM1)*   *U(JMM1)*
C
C                  *V(JMM1)*
C
C     *U(JMM2)*   *EL(JMM2)*   *U(JMM2)*
C                       "                                         V(IM)
C                       "
C                       "                *U(IMM1)* EL(IMM1) U(M)  EL(IM)
C                       "                              = = = = = = =
C                       "                             V(IMM1)     V(IM)
C                       "
C       *U(3)*       *EL(3)*   *U(3)*                      BC ON EL & T
C
C                    *V(3)*
C                 "
C       *U(2)*    "  *EL(2)*   *U(2)*
C                 "         "         "
C                  =  V(2)  "         "  BC ON V
C                           "         "
C          U(1)       EL(1)=    U(1) =   BC ON T
C-----------------------------------------------------------------------
C         IDX IDENTIFIES WHICH VARIABLES ARE CONSIDERED
C              1=SURFACE ELEVATION  
C              2=EXTERNAL MODE U,V 
C              3=INTERNAL MODE U,V 
C              4=TEMP,SAL
C              5=W VELOCITY
C              6=KM,KH,Q2,Q2L,L 
C              7=SURFACE FORCING AND TEMPORAL CYCLING
C              8= CONS. TRACER
C              9= SEDIMENT TRANSPORT
C             10= CHEMICAL TRANSPORT
C         
C-----------------------------------------------------------------------
C
C
      GOTO (10,20,30,40,50,60,70,80,81,181), IDX
C
C     ******************************************************************
C
 10   CONTINUE
C
C-------- EXTERNAL ELEVATION BOUNDARY CCONDITIONS ----------------------
      DO 100 N=1,NUMEBC  
      II=IETA(N)
      JJ=JETA(N)
      ELF(II,JJ)=EBDRY(N)*RAMP
 100  CONTINUE
C
      DO 120 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      ELF(IC,JC)=ELF(ID,JD)
120   CONTINUE
C
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  ELF(I,J)=ELF(I,J)*FSM(I,J)
      RETURN
C
C-------- EXTERNAL VELOCITY BOUNDARY CONDITIONS-------------------------
  20  CONTINUE
      DO 200 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)  
      FRESH=QDIS(N)/(H(IC,JC)+ELF(IC,JC))*RAMP
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            UAF(IC,JC)=-FRESH/H2(IC,JC)
         ELSE
            UAF(ID,JD)=+FRESH/H2(ID,JD)
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            VAF(IC,JC)=-FRESH/H1(IC,JC)
         ELSE
            VAF(ID,JD)=+FRESH/H1(ID,JD)
         ENDIF
      ENDIF
 200  CONTINUE
C
      DO 210 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        VAF(IE,JE)=VAF(IE-1,JE)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        VAF(IE,JE)=VAF(IE+1,JE)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        UAF(IE,JE)=UAF(IE,JE-1)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        UAF(IE,JE)=UAF(IE,JE+1)
      END IF
 210  CONTINUE
      DO 135 J=1,JM
      DO 135 I=1,IM
      UAF(I,J)=UAF(I,J)*DUM(I,J)
 135  VAF(I,J)=VAF(I,J)*DVM(I,J)
      RETURN
C
C-------- INTERNAL VELOCITY BOUNDARY CONDITIONS ------------------------
  30  CONTINUE
      DO 320 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 320 K=1,KBM1
      FRESH=QDIS(N)/(H(IC,JC)+ELF(IC,JC))*RAMP*VQDIST(N,K)/100.
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            UF(IC,JC,K)=-FRESH/(H2(IC,JC)*DZ(K))
         ELSE
            UF(ID,JD,K)=+FRESH/(H2(ID,JD)*DZ(K))
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            VF(IC,JC,K)=-FRESH/(H1(IC,JC)*DZ(K))
         ELSE
            VF(ID,JD,K)=+FRESH/(H1(ID,JD)*DZ(K))
         ENDIF
      ENDIF
  320 CONTINUE
C
      DO 321 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 321 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        VF(IE,JE,K)=VF(IE-1,JE,K)
        WVBOT(IE,JE)=WVBOT(IE-1,JE)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        VF(IE,JE,K)=VF(IE+1,JE,K)
        WVBOT(IE,JE)=WVBOT(IE+1,JE)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE-1,K)
        WUBOT(IE,JE)=WUBOT(IE,JE-1)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE+1,K)
        WUBOT(IE,JE)=WUBOT(IE,JE+1)
      END IF
  321 CONTINUE
C 
      DO 360 K=1,KBM1
      DO 360 I=1,IM
      DO 360 J=1,JM
      UF(I,J,K)=UF(I,J,K)*DUM(I,J)
      VF(I,J,K)=VF(I,J,K)*DVM(I,J)
 360  CONTINUE
      RETURN
C
C-------- TEMPERATURE & SALINITY BOUNDARY CONDITIONS -------------------
  40  CONTINUE
      DO 235 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 235 K=1,KBM1
      IF (VQDIST(N,K).NE.0.0 .AND. QDIS(N).GT.0.0) THEN
      UF(IC,JC,K)=TDIS(N)
      VF(IC,JC,K)=SDIS(N)
      ELSE
      UF(IC,JC,K)=UF(ID,JD,K)
      VF(IC,JC,K)=VF(ID,JD,K)
      END IF
 235  CONTINUE

C  FOR DIFFUSER 
      Do 115 N = 1, NUMDBC
        ID = IDD(N)
        JD = JDD(N)
        Do 112 K = 1, KBM1
          If (VDDIST(N,K).NE.0.0.AND.QDIFF(N).GT.0.0) Then
            TDIF(N,K) = TDIFF(N)
            SDIF(N,K) = SDIFF(N)
          Else
            TDIF(N,K) = UF(ID,JD,K)
            SDIF(N,K) = VF(ID,JD,K)
          End If
  112   Continue
  115 Continue
C
C     ADVECTIVE B.C. FOR T AND S
C
      IF(ALPHA.GT.0.0) THEN
        CALL ALPHABC(DTI2)
      ELSE
        DO 236 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        DO 170 K=1,KBM1
        UF(IE,JE,K)= TBDRY(N,K)
        VF(IE,JE,K)= SBDRY(N,K)
  170   CONTINUE
  236   CONTINUE
      ENDIF
C
      DO 240 K=1,KBM1
      DO 240 I=1,IM
      DO 240 J=1,JM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)
      VF(I,J,K)=VF(I,J,K)*FSM(I,J)
 240  CONTINUE
      RETURN
C
C  FOR CONSERVATIVE TRACER
C
  80  CONTINUE
      DO 2235 N=1,NUMQBCTR
        ID=ITRQD(N)
        JD=JTRQD(N)
        IC=ITRQC(N)
        JC=JTRQC(N)
        DO 2235 K=1,KBM1
          IF (HYDTYPE.EQ.'INTERNAL') THEN
            IF (VQDIST(N,K).NE.0.0 .AND. QDIS(N).GT.0.0) THEN
              UF(IC,JC,K)=CDIS1(N)
            ELSE
              UF(IC,JC,K)=UF(ID,JD,K)
            END IF
          ELSE
            UF(IC,JC,K)=CDIS1(N)
          ENDIF
 2235 CONTINUE
      DO 2236 N=1,NUMEBCTR
        IE=ITRED(N)
        JE=JTRED(N)
        DO 2170 K=1,KBM1
          UF(IE,JE,K)= CBDRY1(N,K)
 2170  CONTINUE
 2236 CONTINUE
C  FOR DIFFUSER
      Do 2115 N = 1, NUMDBC
        ID = IDD(N)
        JD = JDD(N)
        Do 2112 K = 1, KBM1
          If (VDDIST(N,K).NE.0.0.AND.QDIFF(N).GT.0.0) Then
            CDIF1(N,K) = CDIFF1(N)
          Else
            CDIF1(N,K) = CONC1(ID,JD,K)
          End If
 2112   Continue
 2115 Continue
      DO 2240 K=1,KBM1
        DO 2240 I=1,IM
          DO 2240 J=1,JM
            UF(I,J,K)=UF(I,J,K)*FSM(I,J)
 2240 CONTINUE
      RETURN
C
C******************************************************************
C
C  FOR SEDIMENT TRANSPORT
C
  81  CONTINUE
C
c******** RIVER SEDIMENT LOADINGS *******************
      Do 2111 N = 1, NUMQBCSE
        ID=ISEQD(N)
        JD=JSEQD(N)
        IC=ISEQC(N)
        JC=JSEQC(N)
        Do 2101 K = 1, KBM1
          IF (HYDTYPE.EQ.'INTERNAL') THEN
            If (VQDIST(N,K).NE.0.0.AND.QDIS(N).GT.0.0) Then
              UF(IC,JC,K) = CDIS(KS,N)
            Else
              UF(IC,JC,K) = UF(ID,JD,K)
            End If
          ELSE
            UF(IC,JC,K) = CDIS(KS,N)
c            UF(ID,JD,K) = CDIS(KS,N)
          ENDIF
 2101   Continue
 2111 Continue
C
c   ******** BOUNDARY SEDIMENT SPEC *********************
      DO 5230 N=1,NUMEBCSE
        IE=ISEED(N)
        JE=JSEED(N)
        DO 5230 K=1,KBM1
          UF(IE,JE,K)= CBDRY(KS,N,K)
 5230 CONTINUE
C
      DO 5241 K=1,KBM1
        DO 5241 I=1,IM
          DO 5241 J=1,JM
            UF(I,J,K)=UF(I,J,K)*FSM(I,J)
 5241 CONTINUE
C
      RETURN
C
C******************************************************************
C
C  FOR CHEM TRANSPORT 
C
 181  CONTINUE
C
      Do 3111 N = 1, NUMQBCCH
        ID=ICHQD(N)
        JD=JCHQD(N)
        IC=ICHQC(N)
        JC=JCHQC(N)
        Do 3101 K = 1, KBM1
          IF (HYDTYPE.EQ.'INTERNAL') THEN
            If (VQDIST(N,K).NE.0.0.AND.QDIS(N).GT.0.0) Then
              UF(IC,JC,K) = PDIS(KS,N)
            Else
              UF(IC,JC,K) = UF(ID,JD,K)
            End If
          ELSE
            UF(IC,JC,K) = PDIS(KS,N)
          ENDIF
 3101   Continue
 3111 Continue
C
      DO 6230 N=1,NUMEBCCH
        IE=ICHED(N)
        JE=JCHED(N)
        DO 6230 K=1,KBM1
          UF(IE,JE,K)= PBDRY(KS,N,K)
 6230 CONTINUE
C
C     Diffuser Loading
C
      DO 6240 K=1,KBM1
        DO 6240 I=1,IM
          DO 6240 J=1,JM
            UF(I,J,K)=UF(I,J,K)*FSM(I,J)
 6240 CONTINUE
C
C
      RETURN
C
C******************************************************************
C
C-------- VERTICAL VELOCITY BOUNDARY CONDITIONS ------------------------
C
  50  CONTINUE
      DO 245 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 245 K=1,KBM1
      W(IC,JC,K)=W(ID,JD,K)
  245 CONTINUE
      DO 246 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 246 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        W(IE,JE,K)=W(IE-1,JE,K)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        W(IE,JE,K)=W(IE+1,JE,K)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        W(IE,JE,K)=W(IE,JE-1,K)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        W(IE,JE,K)=W(IE,JE+1,K)
      END IF
  246 CONTINUE
      DO 250 K=1,KBM1
      DO 250 J=1,JM
      DO 250 I=1,IM
      W(I,J,K)=W(I,J,K)*FSM(I,J)
 250  CONTINUE
      RETURN
C
C-------- Q2 & Q2L BOUNDARY CONDITIONS ---------------------------------
  60  CONTINUE
      DO 255 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 255 K=1,KBM1
      UF(IC,JC,K)=UF(ID,JD,K)
      VF(IC,JC,K)=VF(ID,JD,K)
      L(IC,JC,K)=L(ID,JD,K)
      KM(IC,JC,K)=KM(ID,JD,K)
      KH(IC,JC,K)=KH(ID,JD,K)
      KQ(IC,JC,K)=KQ(ID,JD,K)
  255 CONTINUE
      DO 256 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 256 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        UF(IE,JE,K)=UF(IE-1,JE,K)
        VF(IE,JE,K)=VF(IE-1,JE,K)
        L(IE,JE,K)=L(IE-1,JE,K)
        KM(IE,JE,K)=KM(IE-1,JE,K)
        KH(IE,JE,K)=KH(IE-1,JE,K)
        KQ(IE,JE,K)=KQ(IE-1,JE,K)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        UF(IE,JE,K)=UF(IE+1,JE,K)
        VF(IE,JE,K)=VF(IE+1,JE,K)
        L(IE,JE,K)=L(IE+1,JE,K)
        KM(IE,JE,K)=KM(IE+1,JE,K)
        KH(IE,JE,K)=KH(IE+1,JE,K)
        KQ(IE,JE,K)=KQ(IE+1,JE,K)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE-1,K)
        VF(IE,JE,K)=VF(IE,JE-1,K)
        L(IE,JE,K)=L(IE,JE-1,K)
        KM(IE,JE,K)=KM(IE,JE-1,K)
        KH(IE,JE,K)=KH(IE,JE-1,K)
        KQ(IE,JE,K)=KQ(IE,JE-1,K)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE+1,K)
        VF(IE,JE,K)=VF(IE,JE+1,K)
        L(IE,JE,K)=L(IE,JE+1,K)
        KM(IE,JE,K)=KM(IE,JE+1,K)
        KH(IE,JE,K)=KH(IE,JE+1,K)
        KQ(IE,JE,K)=KQ(IE,JE+1,K)
      END IF
  256 CONTINUE
      DO 300 K=1,KB
      DO 300 J=1,JM
      DO 300 I=1,IM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)                                      
      VF(I,J,K)=VF(I,J,K)*FSM(I,J)
      L (I,J,K)=L (I,J,K)*FSM(I,J)
      KM(I,J,K)=KM(I,J,K)*FSM(I,J)
      KH(I,J,K)=KH(I,J,K)*FSM(I,J)
      KQ(I,J,K)=KQ(I,J,K)*FSM(I,J)
 300  CONTINUE
      RETURN
C
C-------- SURFACE FORCING AND TEMPORAL CYCLING -------------------------
  70  CONTINUE
C
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 4044
C 
C-------- ELEVATION BOUNDARY CONDITIONS --------------------------------
      IF(NUMEBC.EQ.0) GOTO 4045
       IF(THOUR.LT.T2E) GOTO 4020
       T1E=T2E 
       DO 4030 N=1,NUMEBC
 4030  DEBDRY(N,1)=DEBDRY(N,2) 
       READ (IUT90,4000,END=4010) T2E
       READ (IUT90,4000) (DEBDRY(N,2),N=1,NUMEBC)
 4020  CONTINUE
       FACT=(THOUR-T1E)/(T2E-T1E)
       DO 4040 N=1,NUMEBC
 4040  EBDRY(N)=DEBDRY(N,1)+FACT*(DEBDRY(N,2)-DEBDRY(N,1)) 
 4045 CONTINUE
C 
C-------- TEMPERATURE AND SALINITY BOUNDARY CONDITIONS -----------------
C 
      IF(NUMEBC.EQ.0) GOTO 4135
      IF(THOUR.LT.T2TS) GOTO 4111
      T1TS=T2TS
      DO 4121 N=1,NUMEBC
      DO 4122 K=1,KBM1
      DTBDRY(N,K,1)=DTBDRY(N,K,2)
      DSBDRY(N,K,1)=DSBDRY(N,K,2)
4122  CONTINUE
4121  CONTINUE
      READ (IUT94,4000,END=4011) T2TS
      DO 4136 N=1,NUMEBC
      READ (IUT94,4000) (DTBDRY(N,K,2),K=1,KBM1)
      READ (IUT94,4000) (DSBDRY(N,K,2),K=1,KBM1)
4136  CONTINUE
4111  CONTINUE
      FACT=(THOUR-T1TS)/(T2TS-T1TS)
      DO 4042 N=1,NUMEBC
      DO 4042 K=1,KBM1
      TBDRY(N,K)=DTBDRY(N,K,1)+FACT*(DTBDRY(N,K,2)-DTBDRY(N,K,1)) 
      SBDRY(N,K)=DSBDRY(N,K,1)+FACT*(DSBDRY(N,K,2)-DSBDRY(N,K,1)) 
 4042  CONTINUE
C 
C  DISSOLVED TRACER AT OPEN BOUNDARY
C 
 4044 IF (TRACER.EQ.'INCLUDE'.AND.NUMEBCTR.GT.0) THEN
        IF(THOUR.LT.T2CONOB) GOTO 6111
C
        T1CONOB=T2CONOB
C
        DO 6121 N=1,NUMEBCTR
          DO 6121 K=1,KBM1
            DCBDRY1(N,K,1)=DCBDRY1(N,K,2)
 6121   CONTINUE
C
        READ (IUT501,4000,END=7011) T2CONOB
        DO 6136 N=1,NUMEBCTR
          READ (IUT501,4000) (DCBDRY1(N,K,2),K=1,KBM1)
 6136   CONTINUE
C
 6111   CONTINUE
C
        FACT=(THOUR-T1CONOB)/(T2CONOB-T1CONOB)
C
        DO 6042 N=1,NUMEBCTR
          DO 6042 K=1,KBM1
            CBDRY1(N,K)=DCBDRY1(N,K,1)+FACT*(DCBDRY1(N,K,2)-
     +                                          DCBDRY1(N,K,1)) 
 6042   CONTINUE
      ENDIF
C 
C  SEDIMENT TRANSPORT AT OPEN BOUNDARY
C 
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMEBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
C
          IF(THOUR.LT.T2SEDOBM) GOTO 6010
C
          T1SEDOBM=T2SEDOBM
C
          DO 6020 N=1,NUMEBCSE
            DO 6020 K=1,KBM1
              DCBDRY(KK,N,K,1)=DCBDRY(KK,N,K,2)
 6020     CONTINUE
C
          READ (IUT502,4000,END=7012) T2SEDOBM
          DO 6030 N=1,NUMEBCSE
            READ (IUT502,4000) (DCBDRY(KK,N,K,2),K=1,KBM1)
 6030     CONTINUE
C
 6010     CONTINUE
C
          FACT=(THOUR-T1SEDOBM)/(T2SEDOBM-T1SEDOBM)
C
          DO 6040 N=1,NUMEBCSE
            DO 6040 K=1,KBM1
              CBDRY(KK,N,K)=DCBDRY(KK,N,K,1)+FACT*(DCBDRY(KK,N,K,2)-
     +                                          DCBDRY(KK,N,K,1)) 
 6040     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
C
          IF(THOUR.LT.T2SEDOBS) GOTO 6050
C
          T1SEDOBS=T2SEDOBS
C
          DO 6060 N=1,NUMEBCSE
            DO 6060 K=1,KBM1
              DCBDRY(KK,N,K,1)=DCBDRY(KK,N,K,2)
 6060     CONTINUE
C
          READ (IUT503,4000,END=7212) T2SEDOBS
          DO 6070 N=1,NUMEBCSE
            READ (IUT503,4000) (DCBDRY(KK,N,K,2),K=1,KBM1)
 6070     CONTINUE
C
 6050     CONTINUE
C
          FACT=(THOUR-T1SEDOBS)/(T2SEDOBS-T1SEDOBS)
C
          DO 6080 N=1,NUMEBCSE
            DO 6080 K=1,KBM1
              CBDRY(KK,N,K)=DCBDRY(KK,N,K,1)+FACT*(DCBDRY(KK,N,K,2)-
     +                                          DCBDRY(KK,N,K,1)) 
 6080     CONTINUE
        ENDIF
      ENDIF
C 
C  PARTICLE-BOUND TRACER AT OPEN BOUNDARY
C 
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMEBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
C
          IF(THOUR.LT.T2CHMOBM) GOTO 6110
C
          T1CHMOBM=T2CHMOBM
C
          DO 6120 N=1,NUMEBCCH
            DO 6120 K=1,KBM1
              DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 6120     CONTINUE
C
          READ (IUT504,4000,END=7013) T2CHMOBM
          DO 6130 N=1,NUMEBCCH
            READ (IUT504,4000) (DPBDRY(KK,N,K,2),K=1,KBM1)
 6130     CONTINUE
C
 6110     CONTINUE
C
          FACT=(THOUR-T1CHMOBM)/(T2CHMOBM-T1CHMOBM)
C
          DO 6140 N=1,NUMEBCCH
            DO 6140 K=1,KBM1
              PBDRY(KK,N,K)=DPBDRY(KK,N,K,1)+FACT*(DPBDRY(KK,N,K,2)-
     +                                          DPBDRY(KK,N,K,1)) 
 6140     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
C
          IF(THOUR.LT.T2CHMOBS) GOTO 6150
C
          T1CHMOBS=T2CHMOBS
C
          DO 6160 N=1,NUMEBCCH
            DO 6160 K=1,KBM1
              DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 6160     CONTINUE
C
          READ (IUT505,4000,END=7213) T2CHMOBS
          DO 6170 N=1,NUMEBCCH
            READ (IUT505,4000) (DPBDRY(KK,N,K,2),K=1,KBM1)
 6170     CONTINUE
C
 6150     CONTINUE
C
          FACT=(THOUR-T1CHMOBS)/(T2CHMOBS-T1CHMOBS)
C
          DO 6180 N=1,NUMEBCCH
            DO 6180 K=1,KBM1
              PBDRY(KK,N,K)=DPBDRY(KK,N,K,1)+FACT*(DPBDRY(KK,N,K,2)-
     +                                          DPBDRY(KK,N,K,1)) 
 6180     CONTINUE
        ENDIF
      ENDIF
C
 4135  CONTINUE
C
C-------- RIVER DISCHARGE BOUNDARY CONDITIONS --------------------------
C
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 4074
C
      IF(NUMQBC.EQ.0) GOTO 4075
C
      IF(THOUR.LT.T2Q) GOTO 4050
      T1Q=T2Q 
      DO 4060 N=1,NUMQBC
      DQDIS(N,1)=DQDIS(N,2) 
      DTDIS(N,1)=DTDIS(N,2) 
      DSDIS(N,1)=DSDIS(N,2) 
 4060 CONTINUE
C
      READ (IUT91,4000,END=4012) T2Q
      READ (IUT91,4000) (DQDIS(N,2),N=1,NUMQBC)
      READ (IUT91,4000) (DTDIS(N,2),N=1,NUMQBC)
      READ (IUT91,4000) (DSDIS(N,2),N=1,NUMQBC)
C
 4050 CONTINUE
      FACT=(THOUR-T1Q)/(T2Q-T1Q)
      DO 4070 N=1,NUMQBC
      QDIS(N)=DQDIS(N,1)+FACT*(DQDIS(N,2)-DQDIS(N,1)) 
      TDIS(N)=DTDIS(N,1)+FACT*(DTDIS(N,2)-DTDIS(N,1)) 
      SDIS(N)=DSDIS(N,1)+FACT*(DSDIS(N,2)-DSDIS(N,1)) 
 4070 CONTINUE
C
C  CONSERVATIVE TRACER
C
 4074 IF (TRACER.EQ.'INCLUDE'.AND.NUMQBCTR.GT.0) THEN
        IF (THOUR.LT.T2CON) GOTO 5050
C
        T1CON=T2CON 
C
        DO 5060 N=1,NUMQBCTR
          DCDIS1(N,1)=DCDIS1(N,2)
 5060   CONTINUE
C
        READ (IUT601,4000,END=5012) T2CON
        READ (IUT601,4000) (DCDIS1(N,2),N=1,NUMQBCTR)
C
 5050   CONTINUE
C
        FACT=(THOUR-T1CON)/(T2CON-T1CON)
C
        DO 5070 N=1,NUMQBCTR
          CDIS1(N)=DCDIS1(N,1)+FACT*(DCDIS1(N,2)-DCDIS1(N,1)) 
 5070   CONTINUE
      ENDIF
C
C  SEDIMENT TRANSPORT   RIVER DISCHARGE
C
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMQBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2SEDM) GOTO 5155
          T1SEDM=T2SEDM 

          DO 5145 N=1,NUMQBCSE
            DCDIS(KK,N,1)=DCDIS(KK,N,2)
 5145     CONTINUE
C
          READ (IUT602,4000,END=5213) T2SEDM
          READ (IUT602,4000) (DCDIS(KK,N,2),N=1,NUMQBCSE)
C
 5155     CONTINUE
C
          FACT=(THOUR-T1SEDM)/(T2SEDM-T1SEDM)
C
          DO 5165 N=1,NUMQBCSE
            CDIS(KK,N)=DCDIS(KK,N,1)+FACT*(DCDIS(KK,N,2)
     +                                  -DCDIS(KK,N,1))
 5165     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
C
          IF (THOUR.LT.T2SEDS) GOTO 5156
          T1SEDS=T2SEDS 
C
          DO 5146 N=1,NUMQBCSE
            DCDIS(KK,N,1)=DCDIS(KK,N,2)
 5146     CONTINUE
C
          READ (IUT603,4000,END=5313) T2SEDS
          READ (IUT603,4000) (DCDIS(KK,N,2),N=1,NUMQBCSE)
C
 5156     CONTINUE
C
          FACT=(THOUR-T1SEDS)/(T2SEDS-T1SEDS)
C
          DO 5166 N=1,NUMQBCSE
            CDIS(KK,N)=DCDIS(KK,N,1)+FACT*(DCDIS(KK,N,2)
     +                                  -DCDIS(KK,N,1))
 5166     CONTINUE
        ENDIF
      ENDIF
C
C  PARTICLE-BOUND TRACER TRANSPORT
C
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMQBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2CHMM) GOTO 5250
          T1CHMM=T2CHMM 
C
          DO 5240 N=1,NUMQBCCH
            DPDIS(KK,N,1)=DPDIS(KK,N,2)
 5240     CONTINUE
C
          READ (IUT604,4000,END=5014) T2CHMM
          READ (IUT604,4000) (DPDIS(KK,N,2),N=1,NUMQBCCH)
C
 5250     CONTINUE
C
          FACT=(THOUR-T1CHMM)/(T2CHMM-T1CHMM)
C
          DO 5260 N=1,NUMQBCCH
            PDIS(KK,N)=DPDIS(KK,N,1)+FACT*(DPDIS(KK,N,2)
     +                                  -DPDIS(KK,N,1))
 5260     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
C
          IF (THOUR.LT.T2CHMS) GOTO 5255
          T1CHMS=T2CHMS 
          DO 5245 N=1,NUMQBCCH
            DPDIS(KK,N,1)=DPDIS(KK,N,2)
 5245     CONTINUE
          READ (IUT605,4000,END=5214) T2CHMS
          READ (IUT605,4000) (DPDIS(KK,N,2),N=1,NUMQBCCH)
 5255     CONTINUE
C
          FACT=(THOUR-T1CHMS)/(T2CHMS-T1CHMS)
C
          DO 5265 N=1,NUMQBCCH
            PDIS(KK,N)=DPDIS(KK,N,1)+FACT*(DPDIS(KK,N,2)
     +                                  -DPDIS(KK,N,1))
 5265     CONTINUE
        ENDIF
      ENDIF
C
 4075 CONTINUE
C 
C-------- DIFFUSER INTAKE/OUTFALL CONDITIONS ---------------------------
C         (REGULAR INTAKE/OUTFALLS)
C 
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 4104
C
      IF(NUMDBC1.EQ.0) GOTO 4105
      IF(THOUR.LT.T2D) GOTO 4080
      T1D=T2D
      DO N=1,NUMDBC1
        DQDIFF(N,1)=DQDIFF(N,2) 
        DTDIFF(N,1)=DTDIFF(N,2) 
        DSDIFF(N,1)=DSDIFF(N,2) 
      ENDDO   
C
      READ (IUT92,4000,END=4014) T2D
      READ (IUT92,4000) (DQDIFF(N,2),N=1,NUMDBC1)
      READ (IUT92,4000) (DTDIFF(N,2),N=1,NUMDBC1)
      READ (IUT92,4000) (DSDIFF(N,2),N=1,NUMDBC1)
C
 4080 CONTINUE
      FACT=(THOUR-T1D)/(T2D-T1D)
      DO N=1,NUMDBC1
        QDIFF(N)=DQDIFF(N,1)+FACT*(DQDIFF(N,2)-DQDIFF(N,1)) 
        TDIFF(N)=DTDIFF(N,1)+FACT*(DTDIFF(N,2)-DTDIFF(N,1)) 
        SDIFF(N)=DSDIFF(N,1)+FACT*(DSDIFF(N,2)-DSDIFF(N,1)) 
      ENDDO     
 4105 CONTINUE
c
C-------- DIFFUSER INTAKE/OUTFALL CONDITIONS IN LOOPS-------------------------
C           (COUPLED)
c
      If (NUMDBC2.EQ.0) GO TO 640
      If (THOUR.GE.T2D2) Then
        T1D2 = T2D2
        Do N = NUMDBC1+1, NUMDBC1+NUMDBC2
          DQDIFF(N,1) = DQDIFF(N,2)
          DTDIFF(N,1) = DTDIFF(N,2)
          DSDIFF(N,1) = DSDIFF(N,2)
        Enddo    
        Read (IUT96,4000,End=4514) T2D2
        Read (IUT96,4000) (DQDIFF(N,2),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        Read (IUT96,4000) (DTDIFF(N,2),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
        Read (IUT96,4000) (DSDIFF(N,2),N = NUMDBC1+1,NUMDBC1+NUMDBC2)
      End If
C
      FACT = (THOUR-T1D2) / (T2D2-T1D2)
      Do N = NUMDBC1+1, NUMDBC1+NUMDBC2,2
        ID= IDD(N)
        JD= JDD(N)
        QDIFF(N)  = DQDIFF(N,1) + FACT * (DQDIFF(N,2)-DQDIFF(N,1))
        QDIFF(N+1)= DQDIFF(N+1,1) +
     .              FACT * (DQDIFF(N+1,2)-DQDIFF(N+1,1))
        AVGTEMP= 0.0
        AVGSAL = 0.0
        VOLUME = 0.0
        DO K=1,KBM1
          WRITE(55,'(20F8.3)')T(ID,JD,K),S(ID,JD,K)
          AVGTEMP=AVGTEMP+T(ID,JD,K)*DZ(K)*VDDIST(N,K)
          AVGSAL =AVGSAL +S(ID,JD,K)*DZ(K)*VDDIST(N,K)
          VOLUME=VOLUME + DZ(K)*VDDIST(N,K)
        END DO
        IF(VOLUME.EQ.0.0)GO TO  4516
        TDIFF(N)=0.0
        SDIFF(N)=0.0
        TDIFF(N+1)=AVGTEMP/VOLUME + DTDIFF(N+1,1)
     .            + FACT * (DTDIFF(N+1,2)-DTDIFF(N+1,1))
        SDIFF(N+1)=AVGSAL/VOLUME  + DSDIFF(N+1,1)
     .            + FACT * (DSDIFF(N+1,2)-DSDIFF(N+1,1))
C
      ENDDO    
C
  640 CONTINUE
C
C  TRACER FOR DIFFUSERS ARE COMBINED FOR REGULAR AND COUPLED DIFFUSER CASES.
C  USER CAN SPECIFY ANY DIFFUSERS FOR TRACER INPUT (EITHER REGULAR OR COUPLED).
C  THE TRACER CONC. SPECIFIED FOR COUPLED DIFFUSERS ARE NOT TRACER INCREMENTS
C  AS USED IN SALINITY AND TEMPERATURE CASES. 
C
4104  IF (TRACER.NE.'INCLUDE') GO TO 650
      IF (NUMDBCTR.EQ.0)GO TO 5102
      IF (NUMDBCTR1.EQ.0)GO TO 5101
      IF (THOUR.LT.T2DCON) GOTO 5080
C
      T1DCON=T2DCON
C
      DO N=1,NUMDBCTR1
        DCDIFF1(N,1)=DCDIFF1(N,2)
      ENDDO    
C
      READ (IUT98,4000) T2DCON
      READ (IUT98,4000) (DCDIFF1(N,2),N=1,NUMDBCTR1)
C
 5080 CONTINUE
C
      FACT=(THOUR-T1DCON)/(T2DCON-T1DCON)
C
      DO N=1,NUMDBCTR1
        CDIFF1(N)=DCDIFF1(N,1)+FACT*(DCDIFF1(N,2)-DCDIFF1(N,1))
      ENDDO    
 5101 CONTINUE
C
C  TRACER FOR DIFFUSER-IN-LOOPS 
      IF(NUMDBCTR2.EQ.0) GO TO 5102
C
      IF (THOUR.LT.T2D2CON) GOTO 580
C
      T1D2CON=T2D2CON
C
      DO N=NUMDBCTR1+1,NUMDBCTR
        DCDIFF1(N,1)=DCDIFF1(N,2)
      ENDDO 
C
      READ (IUT99,4000) T2D2CON
      READ (IUT99,4000) (DCDIFF1(N,2),N=NUMDBCTR1+1,NUMDBCTR)
C
 580  CONTINUE
C
      FACT=(THOUR-T1D2CON)/(T2D2CON-T1D2CON)
      Do N = NUMDBCTR1+1, NUMDBCTR,2
        ID= IDD(N)
        JD= JDD(N)
        AVGTR= 0.0
        VOLUME = 0.0
        DO K=1,KBM1
          AVGTR=AVGTR+CONC1(ID,JD,K)*DZ(K)*VDDIST(N,K)
          VOLUME=VOLUME + DZ(K)*VDDIST(N,K)
        END DO
        IF(VOLUME.EQ.0.0)GO TO  4516
        CDIFF1(N)=0.0    ! INTAKE CONC.
        CDIFF1(N+1)=AVGTR/VOLUME + DCDIFF1(N+1,1)
     .             + FACT * (DCDIFF1(N+1,2)-DCDIFF1(N+1,1))

C
      ENDDO       
      WRITE(55,'(20E12.4)')TIME,CDIFF1(1),CDIFF1(2)
C
 5102 CONTINUE
C
C  ******************************************************************
C  PARTICLE-BOUND TRACER LOADINGS THRU DIFFUSER 
C
 650  IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMDBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2CHMDM) GOTO 7250
          T1CHMDM=T2CHMDM
C
          DO 7240 N=1,NUMDBCCH
            DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7240     CONTINUE
          READ (IUT704,4000,END=7014) T2CHMDM
          READ (IUT704,4000) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
 7250     CONTINUE
C
          FACT=(THOUR-T1CHMDM)/(T2CHMDM-T1CHMDM)
C
          DO 7260 N=1,NUMDBCCH
            PDIFF(KK,N)=DPDIFF(KK,N,1)+FACT*(DPDIFF(KK,N,2)
     +                                  -DPDIFF(KK,N,1))
 7260     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
C
          IF (THOUR.LT.T2CHMDS) GOTO 7255
          T1CHMDS=T2CHMDS
C
          DO 7245 N=1,NUMDBCCH
            DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7245     CONTINUE
          READ (IUT705,4000,END=5214) T2CHMDS
          READ (IUT705,4000) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
 7255     CONTINUE
C
          FACT=(THOUR-T1CHMDS)/(T2CHMDS-T1CHMDS)
          DO 7265 N=1,NUMDBCCH
            PDIFF(KK,N)=DPDIFF(KK,N,1)+FACT*(DPDIFF(KK,N,2)
     +                                  -DPDIFF(KK,N,1))
 7265     CONTINUE
        ENDIF
      ENDIF
      CONTINUE
C
C  POINT SOURCE LOADS:  CONSERVATIVE TRACER
C
      IF (TRACER.EQ.'INCLUDE'.AND.NUMPSTR.GT.0) THEN
        IF (THOUR.LT.T2PSTR) GOTO 5051
C
C  INSTANTANEOUS SPECIFICATION OF PT. SOURCE CONCENTRATION
C
        IF (OPTPSTR.EQ.'CONC') THEN
          DO 5062 N=1,NUMPSTR
            CPSTR(N)=DCPSTR(N,2)
 5062     CONTINUE
        ENDIF
C
        T1PSTR=T2PSTR 
C
        DO 5061 N=1,NUMPSTR
          DCPSTR(N,1)=DCPSTR(N,2)
 5061   CONTINUE
C
        READ (IUT701,4000,END=5012) T2PSTR
        READ (IUT701,4000) (DCPSTR(N,2),N=1,NUMPSTR)
C
 5051   CONTINUE
C
C  INSTANTANEOUS SPECIFICATION OF PT. SOURCE CONCENTRATION
C
        IF (OPTPSTR.EQ.'CONC') GOTO 5081
C
        IF (OPTPSTR.EQ.'MASS') THEN
          FACT=(THOUR-T1PSTR)/(T2PSTR-T1PSTR)
C
          DO 5071 N=1,NUMPSTR
            CPSTR(N)=DCPSTR(N,1)+FACT*(DCPSTR(N,2)-DCPSTR(N,1)) 
 5071     CONTINUE
        ENDIF
      ENDIF
C
C     *********************************************************
C     DIFFUSER SEDIMENT LOADINGS 
C     SEDIMENT TRANSPORT   DIFFUSER DISCHARGE
C
 5081  IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMDBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2SEDDM) GOTO 7155
          T1SEDDM=T2SEDDM

          DO 7145 N=1,NUMDBCSE
            DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7145     CONTINUE
C
          READ (IUT702,4000,END=7313) T2SEDDM
          READ (IUT702,4000) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
C
 7155     CONTINUE
C
          FACT=(THOUR-T1SEDDM)/(T2SEDDM-T1SEDDM)
C
          DO 7165 N=1,NUMDBCSE
            CSDIFF(KK,N)=DCSDIFF(KK,N,1)+FACT*(DCSDIFF(KK,N,2)
     +                                  -DCSDIFF(KK,N,1))
 7165     CONTINUE
        ENDIF
C
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
C
          IF (THOUR.LT.T2SEDDS) GOTO 7156
          T1SEDDS=T2SEDDS
C
          DO 7146 N=1,NUMDBCSE
            DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7146     CONTINUE
C
          READ (IUT703,4000,END=7413) T2SEDDS
          READ (IUT703,4000) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
C
 7156     CONTINUE
C
          FACT=(THOUR-T1SEDDS)/(T2SEDDS-T1SEDDS)
C
          DO 7166 N=1,NUMDBCSE
            CSDIFF(KK,N)=DCSDIFF(KK,N,1)+FACT*(DCSDIFF(KK,N,2)
     +                                  -DCSDIFF(KK,N,1))
 7166     CONTINUE
        ENDIF
       ENDIF
C
C-------- METEOROLOGICAL BOUNDARY CONDITIONS ---------------------------
C-------- AVERAGED -----------------------------------------------------
      IF(OPTMBC(1:8).EQ.'AVERAGED') THEN
       IF(THOUR.LT.T2M) GOTO 4110
       T1M=T2M
       DQPREC(1)=DQPREC(2) 
       DQEVAP(1)=DQEVAP(2) 
       DHFLUX(1)=DHFLUX(2)
       DTX   (1)=DTX   (2) 
       DTY   (1)=DTY   (2) 
       DWSX   (1)=DWSX   (2) 
       DWSY   (1)=DWSY   (2) 
       DWDS   (1)=DWDS   (2) 
       DWDD   (1)=DWDD   (2) 
       READ (IUT93,4000,END=4016) T2M
       READ (IUT93,4000) DHFLUX(2),DTX(2),DTY(2),
     +      DWSX(2),DWSY(2),DWDS(2),DWDD(2),DQPREC(2),DQEVAP(2)
 4110  CONTINUE
       FACT=(THOUR-T1M)/(T2M-T1M)
       QPREC=DQPREC(1)+FACT*(DQPREC(2)-DQPREC(1)) 
       QEVAP=DQEVAP(1)+FACT*(DQEVAP(2)-DQEVAP(1)) 
       HFLUX=DHFLUX(1)+FACT*(DHFLUX(2)-DHFLUX(1))
       TX   =DTX   (1)+FACT*(DTX   (2)-DTX   (1)) 
       TY   =DTY   (1)+FACT*(DTY   (2)-DTY   (1)) 
       UWIND   =DWSX   (1)+FACT*(DWSX   (2)-DWSX   (1)) 
       VWIND   =DWSY   (1)+FACT*(DWSY   (2)-DWSY   (1)) 
        DO 4115 J=1,JM
         DO 4115 I=1,IM
           WU(I,J) = UWIND
           WV(I,J) = VWIND
4115    CONTINUE

C
C  CONVERT WIND DIRECTION BACK TO CORRECT ORIENTATION (180 deg ROTATION)
C
       DO 4120 J=1,JM
         DO 4120 I=1,IM
           WSSURF(I,J)=0.0
           WTSURF(I,J)=RAMP*HFLUX/(-4.186E6)
           WUSURF(I,J)=RAMP*(-1.E-3)*
     +               (+TX*COS(ANG(I,J))+TY*SIN(ANG(I,J)))*FSM(I,J)
           WVSURF(I,J)=RAMP*
     +      (-1.E-3)*(-TX*SIN(ANG(I,J))+TY*COS(ANG(I,J)))*FSM(I,J)
 4120  CONTINUE
C
      ELSE IF(OPTMBC(1:8).EQ.'SYNOPTIC') THEN  
C-------- SYNOPTIC (PRECIPITATION, EVAPORATION & HEAT FLUX) ------------
       IF(THOUR.LT.T2M) GOTO 4210
       T1M=T2M
       DO 4220 J=1,JM
       DO 4220 I=1,IM
 4220  DHFLX2D(I,J,1)=DHFLX2D(I,J,2)
       READ (IUT192,END=4016) T2M
       READ (IUT192) ((DHFLX2D(I,J,2),I=1,IM),J=1,JM)
 4210  CONTINUE
       FACT=(THOUR-T1M)/(T2M-T1M)
       DO 4230 J=1,JM   
       DO 4230 I=1,IM  
 4230  TH2D(I,J)=DHFLX2D(I,J,1)+FACT*(DHFLX2D(I,J,2)-DHFLX2D(I,J,1))
C-------- SYNOPTIC (WIND STRESS) -------------------------------------
C
       IF(THOUR.LT.T2W) GOTO 4310
       T1W=T2W
       DO 4320 J=1,JM
       DO 4320 I=1,IM
       DPATM(I,J,1)=DPATM(I,J,2) 
       DTX2D(I,J,1)=DTX2D(I,J,2) 
 4320  DTY2D(I,J,1)=DTY2D(I,J,2) 
       READ (IUT95,END=4016) T2W
       READ (IUT95) ((DTX2D(I,J,2),DTY2D(I,J,2),DPATM(I,J,2),
     +                I=1,IM),J=1,JM)
 4310  CONTINUE
       FACT=(THOUR-T1W)/(T2W-T1W)
       DO 4330 J=1,JM
       DO 4330 I=1,IM
       PATM(I,J)=DPATM(I,J,1)+FACT*(DPATM(I,J,2)-DPATM(I,J,1))
       TX2D(I,J)=DTX2D(I,J,1)+FACT*(DTX2D(I,J,2)-DTX2D(I,J,1))
 4330  TY2D(I,J)=DTY2D(I,J,1)+FACT*(DTY2D(I,J,2)-DTY2D(I,J,1))

C-------- SYNOPTIC (WIND SPEED FOR WAVEDYN MODEL
       IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
         IF(THOUR.LT.T2WV) GOTO 5310
         T1WV=T2WV
         DO 5320 J=1,JM
         DO 5320 I=1,IM
         WU1(I,J)=WU2(I,J) 
 5320    WV1(I,J)=WV2(I,J) 
         READ (IUT193,END=4016) T2WV
         READ (IUT193) ((WU2(I,J),WV2(I,J),I=1,IM),J=1,JM)
 5310    CONTINUE
         FACT=(THOUR-T1WV)/(T2WV-T1WV)
         DO 5330 J=1,JM
         DO 5330 I=1,IM
         WU(I,J)=WU1(I,J)+FACT*(WU2(I,J)-WU1(I,J))
 5330    WV(I,J)=WV1(I,J)+FACT*(WV2(I,J)-WV1(I,J))
       ENDIF
C
C  ADDED RAMPING TO WIND STRESSES AND HEAT FLUX
C
       DO 4340 J=1,JM
         DO 4340 I=1,IM
           WSSURF(I,J)=0.0
           WTSURF(I,J)=RAMP*TH2D(I,J)/(-4.186E6)
           WUSURF(I,J)=RAMP*(-1.E-3)*(+TX2D(I,J)*COS(ANG(I,J))+
     .                      TY2D(I,J)*SIN(ANG(I,J)))*FSM(I,J)
           WVSURF(I,J)=RAMP*(-1.E-3)*(-TX2D(I,J)*SIN(ANG(I,J))+
     .                      TY2D(I,J)*COS(ANG(I,J)))*FSM(I,J)
           PATM(I,J) = PATM(I,J)*RAMP
 4340  CONTINUE
C
c---------- heatflux calculation ------------------------------------
C
       else If (OPTMBC(1:8).EQ.'LANDPFLX'.OR.
     *          OPTMBC(1:8).EQ.'RANDMFLX'.OR.
     *          OPTMBC(1:8).EQ.'AANDBFLX') Then
        If (THOUR.GE.T2M) Then
c
          T1M = T2M
          DQPREC(1) = DQPREC(2)
          DQEVAP(1) = DQEVAP(2)
          DAIRTM(1) = DAIRTM(2)
          DRELHU(1) = DRELHU(2)
          DBAROP(1) = DBAROP(2)
          DSWOBS(1) = DSWOBS(2)
          DTX(1) = DTX(2)
          DTY(1) = DTY(2)
	  dwsx(1)=dwsx(2)
	  dwsy(1)=dwsy(2)  
          DWDS(1) = DWDS(2)
          DWDD(1) = DWDD(2)
          CLOUD(1)  = CLOUD(2)
          EXTCOEF(1)  = EXTCOEF(2)
C
           Read (IUT93,5000,End=4016) T2M
           Read (IUT93,5000) DAIRTM(2), DRELHU(2),
     *     DBAROP(2), DSWOBS(2), DTX(2), DTY(2), DWSX(2), DWSY(2),
     *     CLOUD(2), EXTCOEF(2),DQPREC(2), DQEVAP(2)
C
        End If
        FACT = (THOUR-T1M) / (T2M-T1M)
        QPREC = DQPREC(1) + FACT * (DQPREC(2)-DQPREC(1))
        QEVAP = DQEVAP(1) + FACT * (DQEVAP(2)-DQEVAP(1))
        AIRTMP = DAIRTM(1) + FACT * (DAIRTM(2)-DAIRTM(1))
        RELHUM = DRELHU(1) + FACT * (DRELHU(2)-DRELHU(1))
        RHBY100 = RELHUM/100.  ! ROSATI and Miyakoda (RANDMFLX)
        BAROP = DBAROP(1) + FACT * (DBAROP(2)-DBAROP(1))
        BAROPNT = BAROP*100.0  !mb to Newton/m2 to use ROSA heatflux
        SWOBS = DSWOBS(1) + FACT * (DSWOBS(2)-DSWOBS(1))
        SWOBS = SWOBS*(1.0-REFL) ! ALLOWING SWOBS TO REFLECT FROM WATER SURFACE (REFL)
        CLDFRC = CLOUD(1) + FACT * (CLOUD(2)-CLOUD(1))
        EXTC   = EXTCOEF(1) + FACT * (EXTCOEF(2)-EXTCOEF(1))
        TX = DTX(1) + FACT * (DTX(2)-DTX(1))
        TY = DTY(1) + FACT * (DTY(2)-DTY(1))
       UWIND   =DWSX   (1)+FACT*(DWSX   (2)-DWSX   (1))
       VWIND   =DWSY   (1)+FACT*(DWSY   (2)-DWSY   (1))
        DO 4116 J=1,JM
         DO 4116 I=1,IM
           WU(I,J) = UWIND
           WV(I,J) = VWIND
4116    CONTINUE

C
C  CONVERT WIND DIRECTION BACK TO CORRECT ORIENTATION (180 deg ROTATION)
C
C
        WSX = DWSX(1) + FACT * (DWSX(2)-DWSX(1))
        WSY = DWSY(1) + FACT * (DWSY(2)-DWSY(1))
c find max and min sea surface temperatures
        SSTMIN = 100.
        SSTMAX = -100.
        Do 335 J = 1, JM
          Do 334 I = 1, IM
            If (FSM(I,J).EQ.1.) Then
              SSTMAX = AMAX1(SSTMAX,T(I,J,1))
              SSTMIN = AMIN1(SSTMIN,T(I,J,1))
            End If
  334     Continue
  335   Continue
c expand the range slightly to avoid divide by zero problems at T=SSTMIN
c and T=SSTMAX
        SSTMIN = SSTMIN - .001
        SSTMAX = SSTMAX + .001
c Set a few constants
        UZ = SQRT(WSX*WSX+WSY*WSY)
        UZ = AMAX1(UZ,0.1)  ! to set minimun value of UZ as 0.1m/sec
        ZU = 10.
        ZT = 10.
        ZQ = 10.
        SPECHT = 1004.  ! Specific Heat (Joules/kg/deg C)
        RLHE = 2.47E3   ! Latent heat of evaporation (Joules/g)
        RHOAIR = (1.25E-03)     ! air density
C       Get Wind Direction for ROSATI and Miyakoda Heatflux
        Call NORTH(WSX,WSY,ANGLE)
C
C  If CLDFRC is negative from Data i.e. missing then we wil calculate CLDFRC
C  from ncld.f, otherwise it will use cloud data ***************************
C
        IF (CLDFRC.LT.0.0) Then
c days since start of run
        TDAY = THOUR / 24.
c Add days since start of run to Julian date of start of run
        CUJDAY = TDAY + SDAY

c subtract out Julian date of Jan 1 from current year after updating
c current year (if necessary)
        If (CUJDAY.LT.FLOAT(J2JDAY)) Then
          YDAY = CUJDAY - FLOAT(J1JDAY)
        Else
          YDAY = CUJDAY - FLOAT(J2JDAY)
          J1JDAY = J2JDAY
          J1YR = J1YR + 1
          J2YR = J1YR + 1
          Call JDAY(1,1,J2YR,J2JDAY)
        End If
c calculate cloud fraction for current day of the year, based on observed
c vs theoretical sw radiation changed CLDFRC TO CLDX IN ARGS 
        Call NCLD(YDAY,ALAT,ALON,SWOBS,ATC,CLDX,CUTOFF,RSS)
        CLDFRC = CLDX
      Endif
C
C  Here We calculate ShortWave Radiation in Advance before go to the loop
C  Under 2 Circumstances 1. If OPTMBC is RANDMFLX or 2. Short Wave Data
C  is not Available i.e. Negative
C
          If (OPTMBC(1:8).EQ.'RANDMFLX'.OR.DSWOBS(1).LT.0.
     +         .OR.DSWOBS(2).LT.0.) Then
            Call DATEHR(IDA,IMO,IYR,IHR,0,SHOUR)
            CHRS = SHOUR + THOUR
            Call GETDATE(IDC,IMC,IYC,IHC,IMINC,CHRS)
            Call HITFLX (ALAT,ALON,IYC,IMC,IDC,IHC,BAROPNT,ZERO,
     +         AIRTMP,UZ,ANGLE,RHBY100,CLDFRC,DHR,SWR,RLW,HNOT,
     +         WENOT)
            SWOBS = SWR
C
C     ALLOWING SHOBS TO REFLECT BACT FROM WATER SURFACE (REFL)
C
            SWOBS = SWOBS*(1.0-REFL)
          Endif

c number of interpolation steps
        NINTPT = 20   ! this value can be changed
        STEP = (SSTMAX-SSTMIN) / FLOAT(NINTPT-1)
        Do 336 I = 1, NINTPT
          SSTINT(I) = SSTMIN + FLOAT(I-1) * STEP
          SSIN = SSTINT(I)
          Call BULK(UZ,ZU,AIRTMP,ZT,RELHUM,ZQ,SSIN,TOL,UW,WTOUT,WQOUT,
     *        CD,BAROP,NERR)

C
C     ********************************************************************
C     NOTE: IN VAPOR ROUTINE  EA=VAPOR PRESSURE AT SSIN (SEA SURF TEMP),
C     ES=SATURATED VAPORT PRESSURE AT SSIN, IF YOU NEED TO GET SATURATED VAPOR
C     PRESSURE OR VAPOR PRESSURE AT AIR TEMPERATURE YOU SHOULD USE AIR TEMP
C     INSTEAD. ALSO NOTE EA = REL HUM*ES
C     For WHOI(LANDPFLX) ROUTINE ALL VAPOR PRESSURE ARE in mbar
C     **********************************************************************
C
          Call VAPOR(EA,ESS,RELHUM,SSIN,BAROP)   !For Sea Surface Temperature
C
          If (OPTMBC(1:8).EQ.'AANDBFLX') Then
C
C******* Atmospheric Longwave radiation ADOPTED FROM ARMY(AANDBFLX) CODE
           T2K = 273.2+AIRTMP
           FAC = 1.0+0.17*CLDFRC**2   !CLDFRC is in 0-1 scale
           RAN = 1000.0/3600.0*9.37E-6*SBC*T2K**6*FAC*(1.0-0.03)
C******* Reflective Longwave radiation
           RB = 1000.0/3600.0*0.97*SBC*(SSIN+273.2)**4
C
           Call VAPOR(EAS,ES,RELHUM,AIRTMP,BAROP)!  For Air Temperature
           ESEADIFF = (ESS-EAS)/1.3333           ! from mb to mmHg
           UZ = UZ* WNDSH                        ! Wind Sheltering Coeff WNDSH
           fw       = 9.2 + 0.46*UZ**2           ! Wind function Edinger et al. 1974
           WT(I)    = -0.47*fw*(SSIN-AIRTMP)     ! sensible heat flux (w/m2)
           WQ(I)    = -fw*ESEADIFF               ! Latent heat flux (w/m2)
           RLWC(I)  = RAN - RB
          Else If (OPTMBC(1:8).EQ.'LANDPFLX') Then
           WT(I)    = -(RHOAIR*1.E3) * SPECHT * WTOUT  !sensible heat flux (w/m2)
           WQ(I)    = -RLHE * WQOUT              !latent heat flux (w/m2)
           Call LONGWAVE(SSIN,AIRTMP,EA,CLDFRC,RL1,RL2,RLW,RLWCO)
           RLWC(I) = -RLWCO
          Else If (OPTMBC(1:8).EQ.'RANDMFLX') Then
            Call HITFLX (ALAT,ALON,IYC,IMC,IDC,IHC,BAROPNT,SSIN,
     +           AIRTMP,UZ,ANGLE,RHBY100,CLDFRC,DHR,SWR,RLW,HNOT,
     +           WENOT)
            WT(I)   = -HNOT     ! Sensible Heat
            WQ(I)   = -WENOT    ! Latent Heat
            RLWC(I) = -RLW      ! Longwave
          Endif

C
C
  336   Continue
c
c       Set one extra value equals to last NINTPT value to make available for
c       Interpolation purposes in next few statements
c
        WT(NINTPT+1) = WT(NINTPT)
        WQ(NINTPT+1) = WQ(NINTPT)
        RLWC(NINTPT+1) = RLWC(NINTPT)
        Do 350 J = 1, JM
          Do 340 I = 1, IM
C**********************************************************************
C   surface heat and salt fluxes
C**********************************************************************
            WSSURF(I,J) = 33. / 1000 * (QPREC-QEVAP) / (ART(I,J)+1.E-16)

c           interpolate heatfluxes from those calculated
c           First, find the interpolation factor

            If (FSM(I,J).EQ.1.) Then
              ILOW = (T(I,J,1)-SSTMIN) / STEP + 1
              IHIGH = ILOW + 1
              FACTWT = (T(I,J,1)-SSTINT(ILOW)) / (SSTINT(IHIGH)-
     *            SSTINT(ILOW))

c Interpolate heat fluxes from those calculated

              WT1 = WT(ILOW) + FACTWT * (WT(IHIGH)-WT(ILOW))
              WT2 = WQ(ILOW) + FACTWT * (WQ(IHIGH)-WQ(ILOW))
              WT3 = RLWC(ILOW) + FACTWT * (RLWC(IHIGH)-RLWC(ILOW))
              HFLUX = WT1 + WT2 + WT3 + SWOBS
              HFLUXmin = Wt(1) + Wq(1) + rlwc(1) + SWOBS
              HFLUXmax = WT(nintpt) + Wq(nintpt) + rlwc(nintpt) + SWOBS

C
C  This is checking purposes writing heatfluxes at first location where skill
C  assessment output is written in gcm_tsr file
      IF(I.EQ.INXIE(1).AND.J.EQ.INXJE(1))THEN
        WT1SAVE = WT1SAVE + WT1*SKILLI
        WT2SAVE = WT2SAVE + WT2*SKILLI
        WT3SAVE = WT3SAVE + WT3*SKILLI
        SWRSAVE = SWRSAVE + SWOBS*SKILLI
        HFSAVE = HFSAVE + HFLUX*SKILLI
        ESSAVE = ESSAVE + ESS*SKILLI
        EASAVE = EASAVE + EAS*SKILLI
        ATSAVE = ATSAVE + AIRTMP*SKILLI
        STSAVE = STSAVE + SSIN*SKILLI
C
       If (ISKILL.NE.0.AND.MOD(INT,ISKILL).EQ.0) Then
        TMIDDLE = TIME - (.5*DTI*DAYI/(SKILLI*2))
        ESSAVE = 0.0
        EASAVE = 0.0
        HFSAVE = 0.0
        WT1SAVE = 0.0
        WT2SAVE = 0.0
        WT3SAVE = 0.0
        SWRSAVE = 0.0
       Endif
      Endif
777     FORMAT(10F12.5)
C
C  ADDED RAMPING TO WIND STRESSES AND HEAT FLUX
C
              WSSURF(I,J) = RAMP* 33. / 1000 * (QPREC-QEVAP) /(ART(I,J)+
     *            1.E-16)
              WTSURF(I,J) = RAMP* HFLUX / (-4.186E6)
            End If
              WUSURF(I,J) = RAMP*(-1.E-3) * 
     +              (+TX*COS(ANG(I,J))+TY*SIN(ANG(I,J)))
              WVSURF(I,J) = RAMP*(-1.E-3) * 
     +              (-TX*SIN(ANG(I,J))+TY*COS(ANG(I,J)))
              SWRAD(I,J)  = RAMP* SWOBS / (-4.186E6)
  340     Continue
  350   Continue
c-----------------------------------------------------------------
      ENDIF
C
C  FOR WIND WAVE INPUT
C
      IF (WAVEDYN.EQ.'EXTERNAL') THEN
        IF (THOUR.LT.T2WAVE) GOTO 8010
C
        T1WAVE=T2WAVE
C
        DO 8020 J=1,JM
          DO 8020 I=1,IM
            DWHT(I,J,1)=DWHT(I,J,2) 
            DWPER(I,J,1)=DWPER(I,J,2) 
            DWDIR(I,J,1)=DWDIR(I,J,2) 
 8020   CONTINUE
C
        READ (111,END=5016) T2WAVE
        READ (111) ((DWHT(I,J,2),I=1,IM),J=1,JM)
        READ (111) ((DWPER(I,J,2),I=1,IM),J=1,JM)
        READ (111) ((DWDIR(I,J,2),I=1,IM),J=1,JM)
C
 8010  CONTINUE
C
        FACT=(THOUR-T1WAVE)/(T2WAVE-T1WAVE)
C
        DO 8030 J=1,JM
          DO 8030 I=1,IM
            HSIG(I,J)=DWHT(I,J,1)+FACT*(DWHT(I,J,2)-DWHT(I,J,1))
            TSIG(I,J)=DWPER(I,J,1)+FACT*(DWPER(I,J,2)-DWPER(I,J,1))
            WAVEDIR(I,J)=DWDIR(I,J,1)+FACT*(DWDIR(I,J,2)-DWDIR(I,J,1))
C
            HSIG(I,J)=RAMP*HSIG(I,J)
 8030   CONTINUE
      ENDIF
C
      RETURN
C
 4000 FORMAT(8E14.7)
 5000 Format (12E14.7)
C-----------------------------------------------------------------------
 4010 WRITE(6,4130) THOUR
 4130 FORMAT(//' THE MODEL HAS RUN OUT OF ELEVATION DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
C
 4011 WRITE(6,4132) THOUR
 4132 FORMAT(//' THE MODEL HAS RUN OUT OF TEMP-SALINITY DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT',//)
      GOTO 4140
C
 4012 WRITE(6,4013) THOUR
 4013 FORMAT(//' THE MODEL HAS RUN OUT OF DISCHARGE DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
C
 4014 WRITE(6,4015) THOUR
 4015 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
C
 4016 WRITE(6,4017) THOUR
 4017 FORMAT(//' THE MODEL HAS RUN OUT OF METEOROLOGICAL DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
c
 4514 WRITE(6,4515) THOUR
 4515 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER DATA IN LOOP AT TIME
     .  ',F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
 4516 WRITE(6,4517) THOUR
 4517 FORMAT(//' THE DIFFUSER IN LOOPS HAS WRONG VDDIST VALUES AT TIME',
     .   F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
C
 5016 WRITE(6,5017) THOUR
 5017 FORMAT(//' THE MODEL HAS RUN OUT OF WIND WAVE DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
C
 5012 WRITE(6,5113) THOUR
 5113 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** DISSOLVED TRACER INPUT ****'//)
      GOTO 4140
C
 5213 WRITE(6,5114) THOUR
 5114 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
C
 5313 WRITE(6,5124) THOUR
 5124 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
C
 5014 WRITE(6,5115) THOUR
 5115 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
C
 5214 WRITE(6,5116) THOUR
 5116 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
C
 7011 WRITE(6,7113) THOUR
 7113 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** DISSOLVED TRACER INPUT ****'//)
      GOTO 4140
C
 7012 WRITE(6,7114) THOUR
 7114 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
C
 7212 WRITE(6,7214) THOUR
 7214 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
C
 7013 WRITE(6,7115) THOUR
 7115 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
C
 7014 WRITE(6,7116) THOUR
 7116 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER TRACER DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
C
 7213 WRITE(6,7117) THOUR
 7117 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
C
 7313 WRITE(6,7124) THOUR
 7124 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER MUD DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
C
 7413 WRITE(6,7424) THOUR
 7424 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER SAND DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
C
 4140 CONTINUE
      CLOSE (IUT90)
      CLOSE (IUT91)
      CLOSE (IUT92)
      CLOSE (IUT93)
      CLOSE (IUT94)
      !!!&&&CALL SYSTEM('rm gcm_temp*')
C
      STOP
      END
