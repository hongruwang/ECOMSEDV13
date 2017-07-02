      SUBROUTINE STRESS
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
C  CALC. BOTTOM SHEAR STRESS DUE TO WIND WAVES AND CURRENTS
C  USES GRANT/MADSEN/GLENN WAVE MODEL
C  REVISION DATE:  July 2, 1999
c     The routine Calculates effective CD due to current and waves using
c     the lowest layer velocity from ECOM and wave information
c     from WAVESMB or WAVEDON Model (PERIOD = mean wave period of energy 
c     density spectrum, UBOT   = rms value of the maxima of the orbital 
c     velocity near the bottom (m/s))
c     This routine is a simplified form of Grant & Madsen (1979)
c     which assumes colinear waves and currents.  This is described
c     in "Effect of wave-current interaction on steady wind-driven
c     circulation in narrow, shallow embayments", Signell et al,
c     JGR, 1990, v 95, n c6, 9671--9678.
C**************************************************************
c     Common block REFLE needed for Grant & Madsen routines
      Common /REFLE/ ZR, Z0C, USTARC, USTACW, USTARW
c     USTARW = the friction velocity due to waves only
c     USTACW = the friction velocity due to combined wave and currents
      Real OMEGA, KBR
      INCLUDE 'comdeck'
      Data IFIRST /0/
      PI = ACOS(-1.)
      PI2=2.*PI
      PIHALF=PI/2.
C
C  CALC. TOTAL BOTOM SHEAR STRESS
C
      DO 30 J=2,JM-1
        DO 30 I=2,IM-1
          IF (FSM(I,J).GT.0.0) THEN
C
C  CALC. WAVE NUMBER (1/m)
C  ASSUME OMEGAR=OMEGA
C  CALC. SHALLOW WATER WAVE SPEED (m/s)
C
            C0=SQRT(GRAV*DT(I,J))
C
C  NO SHEAR STRESS IF D < 10 cm
C
            IF (DT(I,J).LT.0.10) THEN
              TAU(I,J,KB)=0.0
              GOTO 30
            ENDIF
            IF(TSIG(I,J).LE.0.1.OR.HSIG(I,J).LE.0.01) THEN
              TAU(I,J,KB)=10000.*CBC(I,J)*QBAR(I,J)*QBAR(I,J)
              GOTO 30
            ENDIF
C
C  ESTIMATE WAVE NO.
C
            WAVEK0=PI2/(C0*TSIG(I,J))
C
C  ITERATE FOR WAVE NO.
C
            TRAT=PI2*PI2/(GRAV*TSIG(I,J)*TSIG(I,J))
C
C  USE NEWTON'S METHOD
C
            DO 40 N=1,20
              TANHKH=TANH(WAVEK0*DT(I,J))
              FKH=WAVEK0*TANHKH-TRAT
              FPKH=TANHKH+WAVEK0*DT(I,J)*(1.-TANHKH*TANHKH)
              WAVEK=WAVEK0-FKH/FPKH
C
C  CHECK FOR CONVERGENCE
C
              DELTA=ABS(WAVEK/WAVEK0-1.)
              IF (DELTA.LE.0.0001) GOTO 45
              WAVEK0=WAVEK
 40         CONTINUE
C
C  CALC. BOTTOM ORBITAL AMPLITUDE (m)
C
 45         CONTINUE
            ABM=HSIG(I,J)/(2.*SINH(WAVEK*DT(I,J))) + 0.00001
C
C  CALC. WAVE FREQUENCY (1/s)
C
            OMEGA=PI2/TSIG(I,J)
C
C  CALC. BOTTOM ORBITAL VELOCITY (m/s)
C
            UB2=ABM*OMEGA + 0.00001
C
C  HEIGHT OF VELOCITY ABOVE BOTTOM
C
            ZR=DZ(KBM1)*DT(I,J)/2.
            KBR=30.*Z0B
C
C GET CUREENTS FROM MODEL AT HEIGHT ZR ABOVE BED & CONVERT TO CGS UNITS
C
            URSPD = SQRT(UBAR(I,J,KBM1)*UBAR(I,J,KBM1)+
     +            VBAR(I,J,KBM1)*VBAR(I,J,KBM1)) * 100.
C
C  CALC. SHEAR VELOCITY DUE TO WAVES & CURRENTS
C  NEAR BOTTOM CURRENT VELOCITY
C  SET MINIMUM CURRENT VELOCITY TO 1 mm/s
C
            URSPD=AMAX1(URSPD,0.1)
C
C CALL GRANT & MADSEN ROUTINE (CGS UNITS) 
C TO GET EFFECTIVE DRAG COEFFICIENT DUE TO WAVES AND CURRENTS
C
            CALL DRGCO(UB2*100,KBR*100.,OMEGA,URSPD,0.,CDE)
            CALL CDNEW(UB2*100,KBR*100.,OMEGA,URSPD,CDE)
C
C DO NOT LET DRAG COEFFICIENT BE LESS THAN BFRIC
C
            CBC(I,J) = AMAX1(CDE,BFRIC)
            TAU(I,J,KB)=AMAX1(USTARC,USTACW)
            TAU(I,J,KB)=TAU(I,J,KB)*TAU(I,J,KB)
          ENDIF
 30   CONTINUE
C
C     FORCING TAU VALUE = 0 AT CELLS CONNECTED BY RIVERS (IC,JC)
C
      DO 120 N=1,NUMQBC
        IC=IQC(N)
        JC=JQC(N)
        TAU(IC,JC,KB) = 0.0
120   CONTINUE
      RETURN
      END
C*****************************************************************
      Subroutine DRGCO(UB,KBR,OMEGA,UR,DIA,CD)
        SAVE FW
C     ******************************************
C     SUBROUTINE DRGCO          DRAG COEFFICIENT
C     W. D. GRANT                   W.H.O.I.
C     GENERATED                    FEB 14, 1984
C     LATEST REVISION              APR  1, 1984
C     Determines the drag coefficient for mean
C     flow.
C     SUBROUTINES CALLED:
C        WFRIC
C        MOVKBR
C     ******************************************
c     If you want movable bed effects, specify KBR=100. and call
c     this routine with a realistic sand diameter (DIA).
c     If you DON'T want movable bed effects, call this routine
c     with DIA = 0.
      Real KBRAB, KBR
      Common /REFLE/ ZR,z0c,ustarc,ustacw,ustarw
      If (KBR.EQ.100.0) Then
        Call MOVKBR(DIA,UB,OMEGA,KBRAB)
        KBR = KBRAB * (UB/OMEGA)
      Else
        KBRAB = KBR / (UB/OMEGA)
      End If
      Call WFRIC(KBRAB,FW)
      USTARW = SQRT(0.5*FW) * UB
C                       ------------------------------------------
C                       INITIAL DRAG COEF ESTIMATE
      DELTAW = 2.0 * 0.4 * USTARW / OMEGA
      Z0 = KBR / 30.0
      Z0C = Z0 * ((DELTAW/Z0)**0.8)
      CD = (0.4/ALOG(ZR*100./Z0C)) ** 2.0   !ZR converted to cm
      Return
      Entry CDNEW(UB,KBR,OMEGA,UR,CD)
C                       ------------------------------------------
C                       UPDATE CD
      
C      WRITE(4001,*) 'INSIDE CDNEW---- ZR ',ZR
C      WRITE(4001,*) 'FW ',FW
C      WRITE(4001,*) 'Z0 ',Z0
C       WRITE(4001,*) 'UR ',UR
C       WRITE(4001,*) 'UB ',UB,'  KBR ',KBR,' OMEGA ',OMEGA
C   ----- INSERT SIMPLE VERSION OF DRAG COEFF. SUBROUTINE   -------
      Do 100 I = 1, 10
        CD0 = CD
        USTARC = SQRT(CD) * UR
C        UAUB= (CD/(0.64*FW))*(UR/UB)**2.0
C        ALPHA= (1.0+UAUB)**2.0
C        USTACW= ((0.5*FW*ALPHA)**0.5)*UB
        USTAW2 = 0.5 * FW * UB ** 2
        USTAW = SQRT(USTAW2)
        USTACW = SQRT(USTAW2+USTARC**2)
        DELTAW = 2.0 * 0.4 * USTACW / OMEGA
        B = 1.0 - (USTARC/USTACW)
        If (B.LT.0) Then
          Write (6,*) 'B LESS THAN 0'
          Stop
        End If
        Z0C = Z0 * ((DELTAW/Z0)**B)
        CD = (0.4/ALOG(ZR*100./Z0C)) ** 2.0   !ZR converted to cm
        DIFF = ABS(CD-CD0)
        If (DIFF.LE.1.E-07) Go To 110
  100 Continue
C 110 WRITE(6,*) 'CD ',CD
  110 Continue
      Return
 5000 Format (I2,8(1X,F7.4))
      End
      Subroutine WFRIC(KBRAB,FW)
C     ******************************************
C     SUBROUTINE WFRIC             WAVE FRICTION
C     W. D. GRANT                   W.H.O.I.
C     GENERATED                    FEB 14, 1984
C     LATEST REVISION              AUG 21, 1986
C     Determines wave friction factor, FW, from
C     relative roughness.
C     ******************************************
      Real KBRAB
C     ------------------------------------------
C     REGION 3: KBRAB .GT. 1.0
      If (KBRAB.GT.1.0) Then
        FW = 0.23
C     ------------------------------------------
C     REGION 2: KBRAB .LT. 0.08
      Else If (KBRAB.LT.0.08) Then
        FW = 0.13 * KBRAB ** 0.40
C                       ------------------------------------------
C     REGION 1: 0.08 .LE. KBRAB .LE. 1.0
      Else
        FW = 0.23 * KBRAB ** 0.62
      End If
      Return
      End
      Subroutine MOVKBR(DIA,UB,OMEGA,KBRAB)
C     ******************************************
C     SUBROUTINE MOVKBR    MOVEABLE BED ROUGHNESS
C     W. D. GRANT                   W.H.O.I.
C     GENERATED                    APR  1, 1984
C     LATEST REVISION              APR  1, 1984
C     Calculates the relative roughness over
C     moveable beds approximately after Grant &
C     Madsen 1982.
C     All units are CGS.
C     This subroutine assumes quartz sand.
C     ******************************************
      Real KBRAB, KRIPAB, KBREDAB
C     ------------------------------------------
C     CONSTANTS
      S = 2.65
      G = 980.6
      VISC = 0.011
C     ------------------------------------------
C     COMMON FACTORS
      SUBWT = (S-1.0) * G * DIA
      SSTAR = (DIA/(4.0*VISC)) * SQRT(SUBWT)
      PHICR = 0.054
      AB = UB / OMEGA
      KBRAB = DIA / AB
      RE = UB * AB / VISC
      TRANRO = 0.56 / SQRT(RE)
      If (KBRAB.LT.TRANRO) Then
        FW = 0.09 * (1.0/(RE**0.2))
      Else
        Call WFRIC(KBRAB,FW)
      End If
      TAUWB = 1.0 * 0.5 * FW * (UB**2.0)
      PHIPRI = TAUWB / SUBWT
      RATIO = PHIPRI / PHICR
      BREAK = 1.8 * (SSTAR**0.6)
      If (RATIO.LE.BREAK) Then
        STEEP = 0.16 * (1.0/(RATIO**0.04))
        HTAB = 0.22 * (1.0/(RATIO**0.16))
      Else
        STEEP = 0.28 * ((SSTAR**0.6)*(1.0/RATIO))
        HTAB = 0.48 * ((SSTAR**0.8)*(1.0/(RATIO**1.5)))
      End If
      KRIPAB = 30.0 * HTAB * STEEP
      KBREDAB = 27.2 * KBRAB * ((SQRT(RATIO)-0.7)**2.0)
      KBRAB = KRIPAB + KBREDAB
      Write (6,5000) KRIPAB, KBREDAB, KBRAB
      Return
 5000 Format (/,5X,' KRIPAB =',F10.6,/,5X,' KBREDAB =',F10.6,/,5X,
     *    '   KBRAB =',F10.6)
      End
