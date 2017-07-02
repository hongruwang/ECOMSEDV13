      SUBROUTINE ADVCON(FB,F,FMEAN,DTI2,FF,FDIF,TSDIS,TSBDRY,PSLOAD, 
     +                                                         DKRATE)
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
C  FOR TRACER, CONSERVATIVE AND NON-CONSERVATIVE
C
      INCLUDE 'comdeck'
      SAVE
C
C-----------------------------------------------------------------------
C     THIS SUBROUTINE INTEGRATES CONSERVATIVE CONSTITUENT EQUATIONS
C-----------------------------------------------------------------------
C
      DIMENSION FB(IM,JM,KB),F(IM,JM,KB),FF(IM,JM,KB)
      DIMENSION FMEAN(IM,JM,KB),FBSTART(IM,JM,KB),ETA(IM,JM)
      DIMENSION XMFLUX(IM,JM,KB),YMFLUX(IM,JM,KB),ZZFLUX(IM,JM,KB)
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),ZFLUX(IM,JM,KB)
      DIMENSION FDIF(DBCM),TSDIS(QBCM),TSBDRY(EBCM,KBM1)
      DIMENSION PSLOAD(PSCM)
      EQUIVALENCE (XFLUX,A),(YFLUX,C),(ZFLUX,PROD)
      EQUIVALENCE (FBSTART,VH),(ETA,VHP)
C
      IF(SCHEME .EQ. 'CENTRAL   ')THEN
      DO  9 J=1,JM
      DO  9 I=1,IM
      F (I,J,KB)=F (I,J,KBM1)
    9 FB(I,J,KB)=FB(I,J,KBM1)
c
C
C-------- DO ADVECTION FLUXES AS CENTER DIFFERENCES --------------------
      DO 22 K=1,KBM1
      DO 22 J=2,JM
      DO 22 I=2,IM
      XFLUX(I,J,K)=0.5*(F(I-1,J,K)+F(I,J,K)) * XMFL3D(I,J,K)
   22 YFLUX(I,J,K)=0.5*(F(I,J-1,K)+F(I,J,K)) * YMFL3D(I,J,K)
C
C
      DO 236 N=1,NUMEBCTR
      IE=ITRED(N)
      JE=JTRED(N)
      IC=ITREC(N)
      JC=JTREC(N)
C
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
C
C-------- EAST SIDE ----------------------------------------------------
      DO 21 K=1,KBM1
      UUP  =0.5*(XMFL3D(IE,JE,K)+ABS(XMFL3D(IE,JE,K)))
      UDOWN=0.5*(XMFL3D(IE,JE,K)-ABS(XMFL3D(IE,JE,K)))
   21 XFLUX(IE,JE,K)=UUP*FB(IC,JC,K)+UDOWN*FB(IE,JE,K)
      GO TO 236
C
      ELSE IF (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
C
C-------- WEST SIDE ----------------------------------------------------
      DO 222 K=1,KBM1
      UUP  =0.5*(XMFL3D(IC,JC,K)+ABS(XMFL3D(IC,JC,K)))
      UDOWN=0.5*(XMFL3D(IC,JC,K)-ABS(XMFL3D(IC,JC,K)))
  222 XFLUX(IC,JC,K)=UUP*FB(IE,JE,K)+UDOWN*FB(IC,JC,K)
      GO TO 236
C
      ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
C
C-------- NORTH SIDE ---------------------------------------------------
      DO 23 K=1,KBM1
      VUP  =0.5*(YMFL3D(IE,JE,K)+ABS(YMFL3D(IE,JE,K)))
      VDOWN=0.5*(YMFL3D(IE,JE,K)-ABS(YMFL3D(IE,JE,K)))
   23 YFLUX(IE,JE,K)=VUP*FB(IC,JC,K)+VDOWN*FB(IE,JE,K)
      GO TO 236
C
      ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
C
C-------- SOUTH SIDE ---------------------------------------------------
      DO 24 K=1,KBM1
      VUP  =0.5*(YMFL3D(IC,JC,K)+ABS(YMFL3D(IC,JC,K)))
      VDOWN=0.5*(YMFL3D(IC,JC,K)-ABS(YMFL3D(IC,JC,K)))
   24 YFLUX(IC,JC,K)=VUP*FB(IE,JE,K)+VDOWN*FB(IC,JC,K)
C
C-------- DONE ---------------------------------------------------------
      END IF
  236 CONTINUE
C
      DO 121 J=2,JMM1
      DO 121 I=2,IMM1
  121 FF(I,J,1)=-.5*DZR(1)*(F(I,J,1)+F(I,J,2))*W(I,J,2)*ART(I,J)
      DO 130 K=2,KBM1
      DO 130 J=2,JMM1
      DO 130 I=2,IMM1
  130 FF(I,J,K)=.5*DZR(K)*((F(I,J,K-1)+F(I,J,K))*W(I,J,K)
     .                    -(F(I,J,K)+F(I,J,K+1))*W(I,J,K+1))*ART(I,J)
C
C************* ADJUST FOR RIVER/WALL FLUXES *********************
      DO 50 N=1,NUMQBCTR
      ID=ITRQD(N)
      JD=JTRQD(N)
      IC=ITRQC(N)
      JC=JTRQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 60 K=1,KBM1
   60         XFLUX(IC,JC,K)=F(IC,JC,K)*XMFL3D(IC,JC,K)
            ELSE
              DO 70 K=1,KBM1
   70         XFLUX(ID,JD,K)=F(IC,JC,K)*XMFL3D(ID,JD,K)
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 80 K=1,KBM1
   80         YFLUX(IC,JC,K)=F(IC,JC,K)*YMFL3D(IC,JC,K)
            ELSE
              DO 90 K=1,KBM1
   90         YFLUX(ID,JD,K)=F(IC,JC,K)*YMFL3D(ID,JD,K)
            ENDIF
      ENDIF
   50 CONTINUE
C
C-------- ADD NET ADVECTION FLUXES & THEN STEP FORWARD IN TIME --------
C
      DO 141 K=1,KBM1
      DO 141 J=2,JMM1
      DO 141 I=2,IMM1
      FF(I,J,K)=FF(I,J,K)
     .              +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .              +YFLUX(I,J+1,K)-YFLUX(I,J,K)
      FF(I,J,K)=(FB(I,J,K)*(H(I,J)+ETB(I,J))*ART(I,J)-
     +             DTI2*FF(I,J,K))/((H(I,J)+ETF(I,J))*ART(I,J))
C
C  NON-CONSERVATIVE TRACER WITH DKRATE > 0 (FIRST-ORDER DECAY)
C
        FF(I,J,K)=FF(I,J,K)-DTI*DKRATE*F(I,J,K)
 141  CONTINUE
C
      ELSE
C
      DO 10 J=1,JM
      DO 10 I=1,IM
      ETA(I,J)=ETB(I,J)
   10 FB(I,J,KB)=FB(I,J,KBM1)
      DO 11 K=1,KB
      DO 11 J=1,JM
      DO 11 I=1,IM
      FBSTART(I,J,K)=FB(I,J,K)
      XMFLUX(I,J,K)=XMFL3D(I,J,K)
      YMFLUX(I,J,K)=YMFL3D(I,J,K)
      ZZFLUX(I,J,K)=W(I,J,K)
   11 CONTINUE
C-----------------------------------------------------------------------
C            SMOLARKIEWICZ'S SCHEME  
C-----------------------------------------------------------------------
      LOOP=0
 5050 LOOP= LOOP+1
C
C-------- DO ADVECTION FLUXES, FIRST AS UPWIND  -------------
      DO 20 K=1,KBM1
      DO 20 J=2,JM
      DO 20 I=2,IM
      XFLUX(I,J,K)=AMAX1(0.,XMFLUX(I,J,K))*FBSTART(I-1,J,K)+
     .             AMIN1(0.,XMFLUX(I,J,K))*FBSTART(I,J,K)
      YFLUX(I,J,K)=AMAX1(0.,YMFLUX(I,J,K))*FBSTART(I,J-1,K)+
     .             AMIN1(0.,YMFLUX(I,J,K))*FBSTART(I,J,K)
   20 CONTINUE
      DO 120 J=2,JMM1
      DO 120 I=2,IMM1
      ZFLUX(I,J,1)=0.0
  120 ZFLUX(I,J,KB)=0.0
      DO 125 K=2,KBM1
      DO 125 J=2,JMM1
      DO 125 I=2,IMM1
      ZFLUX(I,J,K)=AMAX1(0.,ZZFLUX(I,J,K))*FBSTART(I,J,K)+
     .             AMIN1(0.,ZZFLUX(I,J,K))*FBSTART(I,J,K-1)
  125 ZFLUX(I,J,K)=ZFLUX(I,J,K)*ART(I,J)
C
C----- ADD NET ADVECTIVE FLUXES & STEP FORWARD IN TIME 
C
      DO 140 K=1,KBM1
      DO 140 J=2,JMM1
      DO 140 I=2,IMM1
      FF(I,J,K)=     XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .              +YFLUX(I,J+1,K)-YFLUX(I,J,K)
     .      +DZR(K)*(ZFLUX(I,J,K)  -ZFLUX(I,J,K+1))
      FF(I,J,K)=(FBSTART(I,J,K)*(H(I,J)+ETA(I,J))*
     +    ART(I,J)-DTI2*FF(I,J,K))/((H(I,J)+ETF(I,J))*ART(I,J))
C
C  NON-CONSERVATIVE TRACER WITH DKRATE > 0 (FIRST-ORDER DECAY)
C
        FF(I,J,K)=FF(I,J,K)-DTI*DKRATE*F(I,J,K)
 140  CONTINUE
C
C----- ANTIDIFFUSION VELOCITY -----------------------------------------
C
      IF(SCHEME .EQ. 'UPWIND    ' .OR. LOOP .EQ. 2) GO TO 5051
C
      CALL ANTIDIF(XMFLUX,YMFLUX,ZZFLUX,FB,FF,DTI2,TSDIS,TSBDRY)
C
      DO 146 J=1,JM
      DO 146 I=1,IM
      ETA(I,J)=ETF(I,J)
      DO 146 K=1,KB
      FBSTART(I,J,K)=FF(I,J,K)
  146 CONTINUE
      GO TO 5050
C
      END IF
C-----------------------------------------------------------------------
C            SMOLARKIEWICZ'S SCHEME END
C-----------------------------------------------------------------------
 5051 CONTINUE
C
C-------- ADD DIFFUSIVE FLUXES (HORIZONTAL) ----------------------------
C     DO 30 K=1,KB   ! SHALLOW SEAS
C     DO 30 J=1,JM
C     DO 30 I=1,IM
C  30 FB(I,J,K)=FB(I,J,K)-FMEAN(I,J,K)
      DO 40 K=1,KBM1
      DO 40 J=2,JM
      DO 40 I=2,IM
      XMFLUX(I,J,K)=.5*(AAM(I,J,K)+AAM(I-1,J,K))
      YMFLUX(I,J,K)=.5*(AAM(I,J,K)+AAM(I,J-1,K))
   40 CONTINUE
      DO 100 K=1,KBM1
      DO 100 J=2,JM
      DO 100 I=2,IM
      XFLUX(I,J,K)=
     .    -XMFLUX(I,J,K)/HPRNU*(H(I,J)+H(I-1,J))
     .    *(FB(I,J,K)-FB(I-1,J,K))*DUM(I,J)/(H1(I,J)+H1(I-1,J))
     .    *0.5*(H2(I,J)+H2(I-1,J))
      YFLUX(I,J,K)=
     .    -YMFLUX(I,J,K)/HPRNU*(H(I,J)+H(I,J-1))
     .    *(FB(I,J,K)-FB(I,J-1,K))*DVM(I,J)/(H2(I,J)+H2(I,J-1))
     .    *0.5*(H1(I,J)+H1(I,J-1))
  100 CONTINUE
C     DO 110 K=1,KB   ! SHALLOW SEAS
C     DO 110 J=1,JM
C     DO 110 I=1,IM
C 110 FB(I,J,K)=FB(I,J,K)+FMEAN(I,J,K)
C
C------------- ADJUST FOR RIVER/WALL FLUXES ---------------
      DO 51 N=1,NUMQBCTR
      ID=ITRQD(N)
      JD=JTRQD(N)
      IC=ITRQC(N)
      JC=JTRQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 61 K=1,KBM1
   61         XFLUX(IC,JC,K)=0.0
            ELSE
              DO 71 K=1,KBM1
   71         XFLUX(ID,JD,K)=0.0
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 81 K=1,KBM1
   81         YFLUX(IC,JC,K)=0.0
            ELSE
              DO 91 K=1,KBM1
   91         YFLUX(ID,JD,K)=0.0
            ENDIF
      ENDIF
   51 CONTINUE
C
C----- ADD NET HORIZONTAL FLUXES & STEP FORWARD IN TIME (CENTER DIFF.)
      DO 145 K=1,KBM1
      DO 145 J=2,JMM1
      DO 145 I=2,IMM1
  145 FF(I,J,K)=FF(I,J,K)-DTI2*(
     .              +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .              +YFLUX(I,J+1,K)-YFLUX(I,J,K))
     .                    /((H(I,J)+ETF(I,J))*ART(I,J))
C-----------------------------------------------------------------------
C         IMPOSE TRACER FLUX BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
      DO 150 N=1,NUMDBCTR
      ID=IDD(N)
      JD=JDD(N)
      DO 150 K=1,KBM1
      FF(ID,JD,K)=FF(ID,JD,K)+DTI2*FDIF(N)*QDIFF(N)*RAMP*
     .   VDDIST(N,K)/100./DZ(K)/((H(ID,JD)+ETF(ID,JD))*ART(ID,JD))
  150 CONTINUE
C
C  POINT SOURCE LOADS
C
      DO 160 N=1,NUMPSTR
        IPS=IPSTR(N)
        JPS=JPSTR(N)
        KPS=KPSTR(N)
C
C  MASS LOADING (INPUT AS kg/s, CONVERTED HERE FOR CONC IN mg/l) 
C
        IF (OPTPSTR.EQ.'MASS') THEN
          FF(IPS,JPS,KPS)=FF(IPS,JPS,KPS)+DTI2*(1000.*PSLOAD(N))*RAMP 
     .     /DZ(KPS)/((H(IPS,JPS)+ETF(IPS,JPS))*ART(IPS,JPS))
        ELSE
          IF (PSLOAD(N).GT.0.0) THEN
            FF(IPS,JPS,KPS)=PSLOAD(N)
C
C  RESET INPUT CONC. TO ZERO AFTER INSTANTANEOUS SPECIFICATION
C
            PSLOAD(N)=0.0
          ENDIF
        ENDIF
 160  CONTINUE
C
      RETURN
      END
