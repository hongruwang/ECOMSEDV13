      SUBROUTINE ADVV(DRHOY,ADVVA,DTI2)
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
      DIMENSION DRHOY(IM,JM,KB),XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),
     .                ADVVA(IM,JM),CURV4(IM,JM,KB)
      EQUIVALENCE (A,XFLUX),(C,YFLUX),(VH,CURV4)
C
      DO 10 K=1,KB
      DO 10 J=1,JM
      DO 10 I=1,IM
      VF   (I,J,K)=0.0
      XFLUX(I,J,K)=0.0
   10 YFLUX(I,J,K)=0.0
C
C-------- HORIZONTAL ADVECTION -----------------------------------------
      IF (ADVECT.EQ.'NON-LINEAR') THEN
      DO 20 K=1,KBM1
      DO 20 J=2,JMM1
      DO 20 I=2,IM
   20 XFLUX(I,J,K)=0.25*(XMFL3D(I,J-1,K)+XMFL3D(I,J,K))
     .             *(V(I-1,J,K)+V(I,J,K))
      DO 30 K=1,KBM1
      DO 30 J=2,JMM1
      DO 30 I=2,IMM1
   30 YFLUX(I,J,K)=0.25*(YMFL3D(I,J,K)+YMFL3D(I,J+1,K))
     .             *(V(I,J,K)+V(I,J+1,K))
      END IF
C
C-------- HORIZONTAL DIFFUSION -----------------------------------------
C-------- NOTE:  PROD = DT*AAM LEFT OVER FROM ADVU ---------------------
      DO 40 K=1,KBM1
      DO 40 J=2,JMM1
      DO 40 I=2,IMM1
      XFLUX(I,J,K)=XFLUX(I,J,K)
     .  -PROD(I,J,K)*((UB(I,J,K)-UB(I,J-1,K))
     .  /(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .  +(VB(I,J,K)-VB(I-1,J,K))
     .  /(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)))
     .  *0.25*(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .  *DVM(I,J)*DVM(I-1,J)
      YFLUX(I,J,K)=YFLUX(I,J,K)
     .  -DT(I,J)*AAM(I,J,K)*2.0*(VB(I,J+1,K)-VB(I,J,K))*H1(I,J)/H2(I,J)
   40 CONTINUE
C
C-------- DO VERTICAL ADVECTION ----------------------------------------
C-------- ADD HORIZONTAL ADVECTION -------------------------------------
      IF (ADVECT.EQ.'NON-LINEAR') THEN
      DO 50 J=1,JM
      DO 50 I=1,IM
   50 VF(I,J,1)=0.0
      DO 60 K=2,KBM1
      DO 60 J=3,JMM1
      DO 60 I=2,IMM1
   60 VF(I,J,K)=.25*(W(I,J,K)+W(I,J-1,K))*(V(I,J,K)+V(I,J,K-1))
      DO 70 J=1,JM
      DO 70 I=1,IM
   70 VF(I,J,KB)=0.0
      DO 80 K=1,KBM1
      DO 80 J=3,JMM1
      DO 80 I=2,IMM1
   80 VF(I,J,K)=DZR(K)*(VF(I,J,K)-VF(I,J,K+1))*ARV(I,J)
      END IF
C
      DO 90 K=1,KBM1
      DO 90 J=3,JMM1
      DO 90 I=2,IMM1
   90 VF(I,J,K)=VF(I,J,K)
     .            +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .            +YFLUX(I,J,K)-YFLUX(I,J-1,K)
C
C-------- SAVE HORIZONTAL ADVECTION & DIFFUSION FLUXES -----------------
C------------------- FOR EXTERNAL MODE ---------------------------------
      DO 790 J=1,JM
      DO 790 I=1,IM
 790  ADVVA(I,J)=0.0
      DO 800 K=1,KBM1
      DO 800 J=3,JMM1
      DO 800 I=2,IMM1
 800  ADVVA(I,J)=ADVVA(I,J)+DZ(K)*VF(I,J,K)
C
      DO 305 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        CURV42D(IE,JE)=CURV42D(IE,JE-1)
        ADVVA(IE,JE)=0.0
        DO 310 K=1,KBM1
        CURV4(IE,JE,K)=CURV4(IE,JE-1,K)
 310    VF(IE,JE,K)=0.0
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        CURV42D(IE,JE)=CURV42D(IE,JE+1)
        ADVVA(IE,JE+1)=0.0
        DO 315 K=1,KBM1
        CURV4(IE,JE,K)=CURV4(IE,JE+1,K)
 315    VF(IE,JE+1,K)=0.0
      ENDIF
 305  CONTINUE
C
C--------- +FUD + GDEG/H2 + BAROCLINIC TERM ----------------------------
      DO 340 K=1,KBM1
      DO 340 J=3,JMM1
      DO 340 I=2,IMM1 
 340  VF(I,J,K)=VF(I,J,K)
     .   +ARV(I,J)*(CURV4(I,J,K)*DT(I,J)*(U(I+1,J,K)+U(I,J,K))
     .             +CURV4(I,J-1,K)*DT(I,J-1)*(U(I+1,J-1,K)+U(I,J-1,K)))
     .       +GRAV*.25*(DT(I,J)+DT(I,J-1))
     .       *.5*(EGF(I,J)-EGF(I,J-1)+EGB(I,J)-EGB(I,J-1))
     .       *(H1(I,J)+H1(I,J-1))
     .       +GRAV*.25*(DT(I,J)+DT(I,J-1))
     .       *(DATUM(I,J)-DATUM(I,J-1))*RAMP*(H1(I,J)+H1(I,J-1))
     .       +DRHOY(I,J,K)

c         Add Atmospheric Pressure Terms

     *            +ARV(I,J)*(DT(I,J)+DT(I,J-1))*(PATM(I,J)-PATM(I,J-1))
     *             /RHO0/(H2(I,J)+H2(I,J-1))
C
C-------- STEP FORWARD IN TIME -----------------------------------------
      DO 390 K=1,KBM1
      DO 390 J=3,JMM1
      DO 390 I=2,IMM1
 390  VF(I,J,K)=
     .     ((H(I,J)+ETB(I,J)+H(I,J-1)+ETB(I,J-1))*ARV(I,J)*VB(I,J,K)
     .          -2.*DTI2*VF(I,J,K))
     .    /((H(I,J)+ETF(I,J)+H(I,J-1)+ETF(I,J-1))*ARV(I,J))
      RETURN
      END
