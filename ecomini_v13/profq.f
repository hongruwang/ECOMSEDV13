      SUBROUTINE PROFQ(DT2)
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
      REAL KAPPA,KN
      DIMENSION GH(IM,JM,KB),SM(IM,JM,KB),SH(IM,JM,KB)
      DIMENSION DENOM(IM,JM),KN(IM,JM,KB),BOYGR(IM,JM,KB)
      DIMENSION DH(IM,JM)
      EQUIVALENCE (A,SM),(C,SH),(PROD,KN),(TPS,DH)
      EQUIVALENCE (DTEF,BOYGR,GH)
      a1=0.92
      b1=16.6
      a2=0.74
      b2=10.1
      c1=0.08
      e1=1.8
      e2=1.33
      e3=1.0
      kappa=0.4
      sef=1.
      gee=9.806
      small=1.E-12
C
      DO 10 J=2,JMM1
      DO 10 I=2,IMM1
 10   DH(I,J)=H(I,J)+ETF(I,J)
C
      DO 100 K=2,KBM1
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      A(I,J,K)=-DT2*(KQ(I,J,K+1)+KQ(I,J,K)+2.*UMOL)*.5
     .    /(DZZ(K-1)*DZ(K)*DH(I,J)*DH(I,J))
      C(I,J,K)=-DT2*(KQ(I,J,K-1)+KQ(I,J,K)+2.*UMOL)*.5
     .    /(DZZ(K-1)*DZ(K-1)*DH(I,J)*DH(I,J))
 100  CONTINUE
C
C-----------------------------------------------------------------------
C
C        THE FOLLOWING SECTION SOLVES THE EQUATION
C        DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B
C
C-----------------------------------------------------------------------
C
      CONST1=16.6**.6666667*SEF
      DO 90 J=2,JMM1
      DO 90 I=2,IMM1
      VHP(I,J,1)=SQRT(WUSURF(I,J)**2+WVSURF(I,J)**2)*CONST1
      VH(I,J,1)=0.0
      UF(I,J,KB)=SQRT( (.5*(WUBOT(I,J)+WUBOT(I+1,J)))**2
     .                +(.5*(WVBOT(I,J)+WVBOT(I,J+1)))**2 )*CONST1
  90  CONTINUE
C
      DO 105 K=1,KB
      DO 105 J=2,JMM1
      DO 105 I=2,IMM1
      Q2B(I,J,K)=ABS(Q2B(I,J,K))
 105  Q2LB(I,J,K)=ABS(Q2LB(I,J,K))
C
C-------- INPUT INTERNAL WAVE SHEAR ENERGY ACCORDING TO ----------------
C--------------- 200*(B.V.FREQ)**3 -------------------------------------
      DO 115 K=2,KBM1
      DO 115 J=2,JMM1
      DO 115 I=2,IMM1
      BOYGR(I,J,K)=GEE*(RHO(I,J,K-1)-RHO(I,J,K))/(DZZ(K-1)*DH(I,J))
C     PROD(I,J,K)=KM(I,J,K)*
C    .           200.*(.5*(-BOYGR(I,J,K)+ABS(BOYGR(I,J,K))))**1.5
      PROD(I,J,K)=0.0
 115  CONTINUE
C
      DO 120 K=2,KBM1
      DO 120 J=2,JMM1
      DO 120 I=2,IMM1
      PROD(I,J,K)=PROD(I,J,K)+KM(I,J,K)*.25*SEF
     .       *( (U(I,J,K)-U(I,J,K-1)+U(I+1,J,K)-U(I+1,J,K-1))**2
     .         +(V(I,J,K)-V(I,J,K-1)+V(I,J+1,K)-V(I,J+1,K-1))**2 )
     .              /(DZZ(K-1)*DH(I,J))**2
 120  PROD(I,J,K)=PROD(I,J,K)+KH(I,J,K)*BOYGR(I,J,K)
C
      DO 110 K=2,KBM1
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
 110  DTEF(I,J,K)=Q2B(I,J,K)*SQRT(Q2B(I,J,K))/(B1*Q2LB(I,J,K)+SMALL)
C
      DO 140 K=2,KBM1
      DO 140 J=2,JMM1
      DO 140 I=2,IMM1
      VHP(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-VH(I,J,K-1))
     .    -(2.*DT2*DTEF(I,J,K)+1.) )
      VH(I,J,K)=A(I,J,K)*VHP(I,J,K)
      VHP(I,J,K)=(-2.*DT2*PROD(I,J,K)
     .  +C(I,J,K)*VHP(I,J,K-1)-UF(I,J,K))*VHP(I,J,K)
 140  CONTINUE
C
      DO 150 K=1,KBM1
      KI=KB-K
      DO 150 J=2,JMM1
      DO 150 I=2,IMM1
      UF(I,J,KI)=VH(I,J,KI)*UF(I,J,KI+1)+VHP(I,J,KI)
 150  CONTINUE
C
C-----------------------------------------------------------------------
C
C        THE FOLLOWING SECTION SOLVES THE EQUATION
C        DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
C
C-----------------------------------------------------------------------
C
      DO 155 J=2,JMM1
      DO 155 I=2,IMM1
      VH(I,J,1)=0.0
 155  VHP(I,J,1)=0.0
C
      DO 160 K=2,KBM1
      DO 160 J=2,JMM1
      DO 160 I=2,IMM1
      DTEF(I,J,K) =DTEF(I,J,K)*(1.+E2*((1./ABS(Z(K)-Z(1))+
     .    1./ABS(Z(K)-Z(KB))) *L(I,J,K)/(DH(I,J)*KAPPA))**2)
      VHP(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-VH(I,J,K-1))
     .    -(DT2*DTEF(I,J,K)+1.))
      VH(I,J,K)=A(I,J,K)*VHP(I,J,K)
      VHP(I,J,K)=(DT2*(-PROD(I,J,K)
     .   *L(I,J,K)*E1)+C(I,J,K)*VHP(I,J,K-1)-VF(I,J,K))*VHP(I,J,K)
 160  CONTINUE
C
      DO 170 K=1,KBM1
      KI=KB-K
      DO 170 J=2,JMM1
      DO 170 I=2,IMM1
      VF(I,J,KI)=VH(I,J,KI)*VF(I,J,KI+1)+VHP(I,J,KI)
 170  CONTINUE
C
      DO 112 K=2,KBM1
      DO 112 J=2,JMM1
      DO 112 I=2,IMM1
      IF(UF(I,J,K).GT.SMALL.AND.VF(I,J,K).GT.SMALL) GO TO 112
      UF(I,J,K)=SMALL
      VF(I,J,K)=SMALL
 112  CONTINUE
C
C-----------------------------------------------------------------------
C
C               THE FOLLOWING SECTION SOLVES FOR KM AND KH
C
C-----------------------------------------------------------------------
C
      COEF1=A2*(1.-6.*A1/B1)
      COEF2=3.*A2*B2+18.*A1*A2
      COEF3=A1*(1.-3.*C1-6.*A1/B1)
      COEF4=18.*A1*A1+9.*A1*A2
      COEF5=9.*A1*A2
C
      DO 220 K=2,KBM1
      DO 220 J=2,JMM1
      DO 220 I=2,IMM1
      L(I,J,K)=VF(I,J,K)/UF(I,J,K)
 220  GH(I,J,K)=L(I,J,K)**2/UF(I,J,K)*GEE*(RHO(I,J,K)-RHO(I,J,K-1))
     .   /(-DZZ(K-1)*DH(I,J))
C
C-------- NOTE THAT SM,SH LIMITS TO INFINITY WHEN ----------------------
C--------        GH APPROACHES 0.028 -----------------------------------
      DO 230 K=2,KBM1
      DO 230 J=2,JMM1
      DO 230 I=2,IMM1
 230  GH(I,J,K)=AMIN1(GH(I,J,K),.028)
C
      DO 240 K=2,KBM1
      DO 240 J=2,JMM1
      DO 240 I=2,IMM1
      SH(I,J,K)=COEF1/(1.-COEF2*GH(I,J,K))
      SM(I,J,K)=(COEF3+SH(I,J,K)*COEF4*GH(I,J,K))/(1.-COEF5*GH(I,J,K))
 240  CONTINUE
C
      DO 280 K=2,KBM1
      DO 280 J=2,JMM1
      DO 280 I=2,IMM1
      KN(I,J,K)=L(I,J,K)*SQRT(ABS(Q2(I,J,K)))
      KQ(I,J,K)=(KN(I,J,K)*.41*SM(I,J,K)+KQ(I,J,K))*.5
      KM(I,J,K)=(KN(I,J,K)*SM(I,J,K)+KM(I,J,K))*.5
      KH(I,J,K)=(KN(I,J,K)*SH(I,J,K)+KH(I,J,K))*.5
 280  CONTINUE
C
      RETURN
      END

