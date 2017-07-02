      SUBROUTINE ADVU(DRHOX,ADVUA,DTI2)
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
      DIMENSION DRHOX(IM,JM,KB),XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),
     .           ADVUA(IM,JM),CURV4(IM,JM,KB)
      EQUIVALENCE (A,XFLUX),(C,YFLUX),(VH,CURV4)
C
      DO 95 K=1,KB
      DO 95 J=1,JM
      DO 95 I=1,IM
      UF   (I,J,K)=0.0
      CURV4(I,J,K)=0.0
      XFLUX(I,J,K)=0.0
 95   YFLUX(I,J,K)=0.0
      IF (ADVECT.EQ.'NON-LINEAR')  THEN
      DO 60 K=1,KBM1
      DO 60 J=2,JMM1
      DO 60 I=2,IMM1
      CURV4(I,J,K)=
     .  +.125*((V(I,J+1,K)+V(I,J,K))*
     .     ((H2(I+1,J)*FSM(I+1,J)+H2(I,J)*FSM(I,J))/
     .          (FSM(I+1,J)+FSM(I,J)+1.E-30)-
     .      (H2(I,J)*FSM(I,J)+H2(I-1,J)*FSM(I-1,J))/
     .          (FSM(I,J)+FSM(I-1,J)+1.E-30)) 
     .  -(U(I+1,J,K)+U(I,J,K))*
     .     ((H1(I,J+1)*FSM(I,J+1)+H1(I,J)*FSM(I,J))/
     .          (FSM(I,J+1)+FSM(I,J)+1.E-30)-
     .      (H1(I,J)*FSM(I,J)+H1(I,J-1)*FSM(I,J-1))/
     .          (FSM(I,J)+FSM(I,J-1)+1.E-30)) )
     .                /ART(I,J)
     .            +.25*COR(I,J)
 60   CONTINUE
      ELSE
      DO 62 K=1,KBM1
      DO 62 J=2,JMM1
      DO 62 I=2,IMM1
 62   CURV4(I,J,K)=0.25*COR(I,J)
      END IF
      DO 65 J=1,JM
      DO 65 I=1,IM
  65  CURV42D(I,J)=0.0
      DO 70 K=1,KBM1
      DO 70 J=1,JM
      DO 70 I=1,IM
  70  CURV42D(I,J)=CURV42D(I,J)+CURV4(I,J,K)*DZ(K)
C
C-------- HORIZONTAL ADVECTION -----------------------------------------
      IF (ADVECT.EQ.'NON-LINEAR')  THEN
      DO 100 K=1,KBM1
      DO 100 J=2,JM
      DO 100 I=2,IMM1
 100  XFLUX(I,J,K)=0.25*(XMFL3D(I,J,K)+XMFL3D(I+1,J,K))*
     .             (U(I,J,K)+U(I+1,J,K))
      DO 120 K=1,KBM1
      DO 120 J=2,JMM1
      DO 120 I=2,IMM1
 120  YFLUX(I,J,K)=0.25*(YMFL3D(I-1,J,K)+YMFL3D(I,J,K))*
     .             (U(I,J-1,K)+U(I,J,K))
      END IF
C
C-------- HORIZONTAL DIFFUSION -----------------------------------------
      DO 110 K=1,KBM1
      DO 110 J=1,JM
      DO 110 I=1,IM
 110  PROD(I,J,K)=0.0
      DO 500 K=1,KBM1
      DO 500 J=2,JM
      DO 500 I=2,IM
C
C-------- THIS IS USED IN advv.f AS WELL AS HERE -----------------------
      PROD(I,J,K)=.25*(DT(I,J)+DT(I-1,J)+DT(I,J-1)+DT(I-1,J-1))
     .           *(AAM(I,J,K)+AAM(I-1,J,K)+AAM(I,J-1,K)+AAM(I-1,J-1,K))
 500  CONTINUE
C
      DO 700 K=1,KBM1
      DO 700 J=2,JMM1
      DO 700 I=2,IMM1
      XFLUX(I,J,K)=XFLUX(I,J,K)
     .  -DT(I,J)*AAM(I,J,K)*2.*(UB(I+1,J,K)-UB(I,J,K))*H2(I,J)/H1(I,J)
      YFLUX(I,J,K)=YFLUX(I,J,K)
     .  -PROD(I,J,K)*((UB(I,J,K)-UB(I,J-1,K))
     .  /(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .  +(VB(I,J,K)-VB(I-1,J,K))
     .  /(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)))
     .  *0.25*(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1))
     .  *DUM(I,J)*DUM(I,J-1) 
 700  CONTINUE
C
C-------- DO VERTICAL ADVECTION ----------------------------------------
C-------- ADD HORIZONTAL ADVECTION -------------------------------------
      IF (ADVECT.EQ.'NON-LINEAR')  THEN
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  UF(I,J,1)=0.0
      DO 140 K=2,KBM1
      DO 140 J=2,JMM1
      DO 140 I=2,IMM1
 140  UF(I,J,K)=.25*(W(I,J,K)+W(I-1,J,K))*(U(I,J,K)+U(I,J,K-1))
      DO 160 J=1,JM
      DO 160 I=1,IM
 160  UF(I,J,KB)=0.0
      DO 145 K=1,KBM1
      DO 145 J=2,JMM1
      DO 145 I=2,IMM1
 145  UF(I,J,K)=DZR(K)*(UF(I,J,K)-UF(I,J,K+1))*ARU(I,J)
      END IF
C
      DO 146 K=1,KBM1
      DO 146 J=2,JMM1
      DO 146 I=2,IMM1
 146  UF(I,J,K)=UF(I,J,K)
     .            +XFLUX(I,J,K)-XFLUX(I-1,J,K)
     .            +YFLUX(I,J+1,K)-YFLUX(I,J,K)
C
C-------- SAVE HORIZONTAL ADVECTION & DIFFUSION FLUXES -----------------
C-------------------- FOR EXTERNAL MODE --------------------------------
      DO 790 J=1,JM
      DO 790 I=1,IM
 790  ADVUA(I,J)=0.0
      DO 800 K=1,KBM1
      DO 800 J=2,JMM1
      DO 800 I=2,IMM1
 800  ADVUA(I,J)=ADVUA(I,J)+DZ(K)*UF(I,J,K)
C
      DO 310 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        CURV42D(IE,JE)=CURV42D(IE-1,JE)
        ADVUA  (IE,JE)=0.0
        DO 320 K=1,KBM1
        CURV4(IE,JE,K)=CURV4(IE-1,JE,K)
 320    UF   (IE,JE,K)=0.0
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        CURV42D(IE,JE)=CURV42D(IE+1,JE)
        ADVUA  (IE+1,JE)=0.0
        DO 330 K=1,KBM1
        CURV4(IE,JE,K)=CURV4(IE+1,JE,K)
 330    UF     (IE+1,JE,K)=0.0
      ENDIF
 310  CONTINUE
C
C-------- -FVD + GDEG/H1 + BAROCLINIC TERM -----------------------------
      DO 150 K=1,KBM1
      DO 150 J=2,JMM1
      DO 150 I=3,IMM1
 150  UF(I,J,K)=UF(I,J,K)
     .   -ARU(I,J)*(CURV4(I,J,K)*DT(I,J)*(V(I,J+1,K)+V(I,J,K))
     .             +CURV4(I-1,J,K)*DT(I-1,J)*(V(I-1,J+1,K)+V(I-1,J,K)))
     .        +GRAV*.25*(DT(I,J)+DT(I-1,J))
     .        *.5*(EGF(I,J)-EGF(I-1,J)+EGB(I,J)-EGB(I-1,J))
     .        *(H2(I,J)+H2(I-1,J))
     .        +GRAV*.25*(DT(I,J)+DT(I-1,J))
     .        *(DATUM(I,J)-DATUM(I-1,J))*RAMP*(H2(I,J)+H2(I-1,J))
     .        +DRHOX(I,J,K)

c         Add Atmospheric Pressure Terms 

     *            +ARU(I,J)*(DT(I,J)+DT(I-1,J))*(PATM(I,J)-PATM(I-1,J))
     *             /RHO0/(H1(I,J)+H1(I-1,J))
C
C-------- STEP FORWARD IN TIME -----------------------------------------
      DO 190 K=1,KBM1
      DO 190 J=2,JMM1
      DO 190 I=3,IMM1
 190  UF(I,J,K)=
     .      ((H(I,J)+ETB(I,J)+H(I-1,J)+ETB(I-1,J))*ARU(I,J)*UB(I,J,K)
     .         -2.*DTI2*UF(I,J,K))
     .     /((H(I,J)+ETF(I,J)+H(I-1,J)+ETF(I-1,J))*ARU(I,J))
      RETURN
      END
