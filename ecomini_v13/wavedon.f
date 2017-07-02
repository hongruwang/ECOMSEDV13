      SUBROUTINE WAVEDON
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
C******************************PROGRAM WAVE***************************  
C 
C PURPOSE:     TO PREDICT SURFACE GRAVITY WAVES FROM THE WIND.      
C 
C SUBROUTINES: FUNCTION WINDX - READS METEOROLOGICAL DATA
C
C VARIABLES  :
C     THERE ARE SIX MAJOR ARRAYS USED TO DESCRIBE THE WAVE FIELD.
C     THEY ARE DIMENSIONED (IM,JM). THE SIX ARRAYS ARE:
C
C          XMOM - THE X-COMPONENT OF MOMENTUM
C          YMOM - THE Y-COMPONENT OF MOMENTUM
C          C    - THE MAGNITUDE OF THE MEAN PHASE SPEED
C          S    - THE WAVE HEIGHT STANDARD DEVIATION
C          CS   - THE COSINE OF THE MEAN WAVE DIRECTION
C          SN   - THE SINE OF THE MEAN WAVE DIRECTION
C
C     THE THREE OTHER ARRAYS ARE:
C
C          WU and WV   - THE X AND Y COMPONENTS OF THE WIND
C
C HISTORY :  
C     WRITTEN BY J.R. BENNETT, 1982, GLERL, ANN ARBOR,MI., BY  
C     MODIFYING A PROGRAM WRITTEN BY M.A. DONELAN, NATIONAL WATER  
C     RESEARCH INSTITUTE,BURLINGTON,ONTARIO. 
C 
      INCLUDE 'comdeck'
C
C   PHYSICAL CONSTANTS
C
      SAVE NR
      DATA G,VK/9.806,0.4/
       DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.028,0.0001/
C      DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.028,0.00005/
C     DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.028,0.00/
      DATA NR/0./
C
C   COMPUTE CONSTANTS
C
      FAC2=0.5*GAMMA*RHOAIR/(RHOH2O*G)
      GTPI=G/(2.*ATAN2(0.,-1.))
      PI=4.*ATAN(1.)  
C--------------------------------------------------------------------------
C  SET INITIAL CONDITIONS   THIS SHOULD GO TO first.f
C
      IF(NR.EQ.0)THEN
        DO 10 I=1,IM
        DO 10 J=1,JM
        XMOM(I,J)=0.0
        YMOM(I,J)=0.0
        CO(I,J)=0.0
        CS(I,J)=1.0
        SN(I,J)=0.0
        SIG(I,J)=0.0
        WN(I,J)=1.
   10   CONTINUE
        NR=NR+1
      ENDIF
C
C   THIS PROGRAM USES TWO TIME STEPS. THE LARGE TIME STEP (DTI)
C   IS SPECIFIED BY THE USER AND CONTROLS THE FREQUENCY AT WHICH DATA
C   CAN BE PRINTED OR STORED FOR LATER ANALYSIS BY WOUTP, THE USER SUP-
C   PLIED OUTPUT ROUTINE. THE SMALL TIME STEP (DTT) IS CALCULATED BY
C   THE PROGRAM FROM THE FORMULA FOR NUMERICAL STABILITY OF THE UPWIND
C   FINITE DIFFERENCE SCHEME. THE NUMBER OF SMALL TIME STEPS PER LARGE
C   TIME STEP IS NDTT; IT IS STORED IN COMMON BLOCK /TPARM/.
C   INSIDE EACH SMALL TIME STEP THERE ARE TWO LOOPS OVER THE SPATIAL
C   COORDINATES I AND J. IN THE FIRST LOOP THE WAVE ADVECTION TERMS ARE
C   CALCULATED. IN THE SECOND LOOP THE WIND INPUT OF MOMENTUM IS ADDED
C   AND THE DIAGNOSTIC VARIABLES ARE CALCULATED.
C
C  CALCULATE SMALL TIME STEP BASED ON CONSTANT WIND FROM MIDDLE OF TIME STEP
C
       WMAX=0.
       NMAX=0
       DMAX=99999.
       DO 40 I=2,IMM1
       DO 40 J=2,JMM1
       IF(FSM(I,J).EQ.0.) GO TO 40
       WIN = SQRT(WU(I,J)**2 + WV(I,J)**2)
       DS  = SQRT(H1(I,J)**2 + H2(I,J)**2)
       WMAX=AMAX1(WMAX,WIN,ABS(CO(I,J)))
       DS  = AMIN1(DMAX,DS)
       DMAX = DS
   40  CONTINUE
       NDTT=IFIX((1.414*WMAX*DTI*NWAVE)/DS)+2
       NDTT=MAX0(NMAX,NDTT)
C       NDTT=IFIX((1.414*WMAX*DTI)/DS)+2
       DTT=DTI*FLOAT(NWAVE)/NDTT
       DTTGAM=FAC2*DTT
C
C  BEGIN LOOP OVER SMALL TIME STEPS
C
      DO 60 NDT=1,NDTT
C
C   FIRST LOOP OVER I,J: CALCULATE THE WAVE ADVECTION TERMS
C
      DO 70 I=2,IMM1
      DO 70 J=2,JMM1
      IF(FSM(I,J).EQ.0.) GO TO 70
      DTTDS=DTT*0.25/ART(I,J)
C
C   X-COMPONENT    UPWIND SCHEME
C
      IF(CS(I,J).GE.0.0)
     .  XFLXD=H2(I,J)*((SIG(I,J)*CS(I,J))**2+0.5*SIG(I,J)**2)
     .  -H2(I-1,J)*((SIG(I-1,J)*CS(I-1,J))**2+0.5*SIG(I-1,J)**2)
      IF(CS(I,J).LT.0.0)
     .  XFLXD=H2(I+1,J)*((SIG(I+1,J)*CS(I+1,J))**2+0.5*SIG(I+1,J)**2)
     .  -H2(I,J)*((SIG(I,J)*CS(I,J))**2+0.5*SIG(I,J)**2)
      IF(SN(I,J).LT.0.0)
     .  YFLXD=H1(I,J+1)*(CS(I,J+1)*SN(I,J+1)*SIG(I,J+1)**2)
     .  -H1(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
      IF(SN(I,J).GE.0.0)
     .  YFLXD=H1(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
     .  -H1(I,J-1)*(CS(I,J-1)*SN(I,J-1)*SIG(I,J-1)**2)

      XMOM(I,J)=XMOM(I,J)-DTTDS*(XFLXD+YFLXD)
C
C   Y-COMPONENT    UPWIND SCHEME
C
      IF(CS(I,J).LT.0.0)
     .  XFLXD=H2(I+1,J)*(CS(I+1,J)*SN(I+1,J)*SIG(I+1,J)**2)
     .  -H2(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
      IF(CS(I,J).GE.0.0)
     .  XFLXD=H2(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
     .  -H2(I-1,J)*(CS(I-1,J)*SN(I-1,J)*SIG(I-1,J)**2)
      IF(SN(I,J).GE.0.0)
     .  YFLXD=H1(I,J)*((SN(I,J)*SIG(I,J))**2+0.5*SIG(I,J)**2)
     .  -H1(I,J-1)*((SN(I,J-1)*SIG(I,J-1))**2+0.5*SIG(I,J-1)**2)
      IF(SN(I,J).LT.0.0)
     .  YFLXD=H1(I,J+1)*((SN(I,J+1)*SIG(I,J+1))**2+0.5*SIG(I,J+1)**2)
     .  -H1(I,J)*((SN(I,J)*SIG(I,J))**2+0.5*SIG(I,J)**2)

         YMOM(I,J)=YMOM(I,J)-DTTDS*(XFLXD+YFLXD)
 80   CONTINUE
C
   70 CONTINUE
C
C   SECOND LOOP OVER I,J: CALCULATE THE WIND INPUT OF MOMENTUM
C   AND COMPUTE THE DIAGNOSTIC VARIABLES (CS,SN,S,C)
C
      DO 90 I=2,IMM1
      DO 90 J=2,JMM1
      IF(FSM(I,J).EQ.0.) GO TO 90
C
C   CALCULATE THE VERTICAL MOMENTUM FLUX IN THE 2 DIRECTIONS,
C   THE WIND AND THE WAVE FIELD
C   NOTE WIND IS ROTATED TO E1 AND E2 DIRECTIONS FROM E-W AND N-S
C
      UWIND=(WU(I,J)*COS(ANG(I,J))+WV(I,J)*SIN(ANG(I,J)))*FSM(I,J)
      VWIND=(-WU(I,J)*SIN(ANG(I,J))+WV(I,J)*COS(ANG(I,J)))*FSM(I,J)
      WSPDSQ=UWIND*UWIND+VWIND*VWIND
      WNDSPD=SQRT(WSPDSQ)+1.E-10
      CSWIND=UWIND/WNDSPD
      SNWIND=VWIND/WNDSPD
      CTHE1=CSWIND*CS(I,J)+SNWIND*SN(I,J)
      SG=AMAX1(0.005,SIG(I,J)*ABS(CTHE1))
      DRAG1=(VK/ALOG(50./SG))**2
      SG=AMAX1(0.000005,SIG(I,J))
      DRAG2=(VK/ALOG(50./SG))**2
      B=CO(I,J)*.83/WNDSPD
      AA=1.-B*CTHE1
      WI1=DRAG1*AA*ABS(AA)
      AA=CTHE1-B
      WI2=DRAG2*AA*ABS(AA)
      VFLXA=WI1*CSWIND+WI2*CS(I,J)
      VFLYA=WI1*SNWIND+WI2*SN(I,J)
C
C   COMPUTE BOTTOM FRICTION LOSS FOR SHALLOW WATER USING STANDARD DEVIATION
C    OF WAVEHEIGHT FOR AMPLITUDE AND WAVENUMBER OF PEAK ENERGY FREQUENCY
C
      IF(DT(I,J)*WN(I,J).LT.10.) THEN
       UBOT=SIG(I,J)*WN(I,J)*CO(I,J)/SINH(WN(I,J)*DT(I,J))
      ELSE
       UBOT=0.
      ENDIF
      DTFAC=DTTGAM*WSPDSQ
      XMOM(I,J)=XMOM(I,J)+DTFAC*VFLXA-DTT*CBF*UBOT*UBOT*CS(I,J)
      YMOM(I,J)=YMOM(I,J)+DTFAC*VFLYA-DTT*CBF*UBOT*UBOT*SN(I,J)
C
C   COMPUTE THE NEW COSINES AND SINES; AND COMPUTE THE PHASE
C   SPEED (C) AND WAVE HEIGHT STANDARD DEVIATION (S).
C
      CM=SQRT(XMOM(I,J)**2+YMOM(I,J)**2)+1.E-10
      CS(I,J)=XMOM(I,J)/CM
      SN(I,J)=YMOM(I,J)/CM
      UU=UWIND*CS(I,J)+VWIND*SN(I,J)
      IF (UU.LE.0.0) GO TO 150
      FM=0.01788735*(((UU*UU)/(CM*CM*CM))**0.14285714)
      IF((UU*FM).GT.0.760545) GO TO 160
  150 FM=(CM*14343.09)**(-.3333333333)
  160 CONTINUE
      SIG(I,J)=SQRT(CM * GTPI/FM)
      WN(I,J)=WNUM(FM,DT(I,J))
      IF(DT(I,J)*WN(I,J).LT.3.14) THEN
       CO(I,J)=SQRT(G*TANH(WN(I,J)*DT(I,J))/WN(I,J))
      ELSE
       CO(I,J)=GTPI/FM
      ENDIF
C
   90 CONTINUE
C
C   END LOOP OVER SMALL TIME STEPS
C
   60 CONTINUE
C
C  CALCULATE WAVE PARAMETERS  (Significant Wave Height, Wave Period
C                              and Wave Direction)
C 
      DO 21 I=1,IM  
      DO 21 J=1,JM  
      HSIG(I,J)=-999.0 
      IF(FSM(I,J).EQ.0.) GO TO 21
      HSIG(I,J)=4.*SIG(I,J) 
   21 CONTINUE  
      DO I=1,IM
      DO J=1,JM
       If (FSM(I,J).EQ.1.) then
          HBREAK=0.78*DT(I,J)
          HSIG(I,J)=AMIN1(HSIG(I,J),HBREAK)
       END IF
      END DO
      END DO
C 
C  SAVE WAVE PERIODS  
C 
      DO 22 I=1,IM  
      DO 22 J=1,JM  
      IF(FSM(I,J).EQ.0.) GO TO 22 
      TSIG(I,J)=2.*PI*CO(I,J)/9.806
22    CONTINUE  
C 
C   SAVE WAVE DIRECTIONS  
C 
      DO 23 I=1,IM  
      DO 23 J=1,JM  
      IF(FSM(I,J).EQ.0.) GO TO 23 
      UU=XMOM(I,J) 
      VV=YMOM(I,J) 
      WAVEDIR(I,J) = 0.
      IF(UU .EQ. 0. .AND. VV .EQ. 0.) GOTO 23

C     Donelan and Schwab's convention of wave direction w.r.t. zeta1 axis
c     wave is moving towards the direction not from, counter clock is positive 

      WAVEDIR(I,J)=ATAN2(VV,UU)*180./PI  ! w.r.t. zeta1-direction
      IF(WAVEDIR(I,J) .LT. 0.) WAVEDIR(I,J) = WAVEDIR(I,J) + 360  
      IF(WAVEDIR(I,J) .GT. 360.) WAVEDIR(I,J) = WAVEDIR(I,J) - 360. 

c     Direction with respect to east (x-direction)

      WAVEDIR(I,J) = WAVEDIR(I,J) + ANG(I,J)*180./PI
      WAVEDIR(I,J) = AMOD(WAVEDIR(I,J),360.)

C     IN stress.f and in SMB Model WAVE DIR IS IN w.r.t North Clock Wise
      WAVEDIR(I,J)= 90.-WAVEDIR(I,J)
      IF(WAVEDIR(I,J).LE.0.) WAVEDIR(I,J) = WAVEDIR(I,J) + 360.

23    CONTINUE  
C
      NR=NR+1

C     HARDWARE for Checking Purposes************

      WRITE(1001)TIME
      WRITE(1001)HSIG
      WRITE(1001)TSIG
      WRITE(1001)WAVEDIR
C
  540 FORMAT(1X, 'DEPTH RELATIVE TO MEAN WATER LEVEL')
      RETURN
      END

      FUNCTION WNUM(F,DD)
C
C  APPROXIMATE SOLUTION OF WAVE DISPERSION EQUATION
C
C  F = FREQUENCY (HZ)
C  DD = DEPTH (M)
C  WNUM = WAVENUMBER (1/M)
C
C  REFERENCE: HUNT, J.N. 1979. DIRECT SOLUTION OF WAVE DISPERSION EQUATION.
C   JWPCOD, ASCE, 105(WW4):457-459.
C
      DATA TPI/6.283185308/
      DATA G/9.806/
      Y=DD*(TPI*F)**2/G
      X=Y*(Y+1./(1.+Y*(0.6522+Y*(0.4622+Y*Y*(0.0864+Y*0.0675)))))
      WNUM=SQRT(X)/DD
      RETURN
      END

