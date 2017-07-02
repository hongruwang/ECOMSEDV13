      SUBROUTINE EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2)
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
C-----------------------------------------------------------------------
C     CALCULATE FREE SURFACE HEIGHT
C-----------------------------------------------------------------------
C
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),TRNU(IM,JM),TRNV(IM,JM)
C
      DO 400 J=2,JMM1
      DO 400 I=2,IM
 400  FLUXUA(I,J)=.25*(D(I,J)+D(I-1,J))*(H2(I,J)+H2(I-1,J))*UA(I,J)
      DO 405 J=2,JM
      DO 405 I=2,IMM1
 405  FLUXVA(I,J)=.25*(D(I,J)+D(I,J-1))*(H1(I,J)+H1(I,J-1))*VA(I,J)
C
      DO 410 J=2,JMM1
      DO 410 I=2,IMM1
 410  ELF(I,J)=ELB(I,J)-DTE2*
     . (FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))/ART(I,J)
C-----------------------------------------------------------------------
C         IMPOSE MASS FLUX BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
      DO 406 N=1,NUMDBC
      ID=IDD(N)
      JD=JDD(N)
      ELF(ID,JD)=ELF(ID,JD)+ DTE2*QDIFF(N)*RAMP /ART(ID,JD)
  406 CONTINUE
C
      CALL BCOND(1,DTI2,0)
C
C     Introduced Inverted Reid and Bodine (IRANDB) and Optimized Clamped 
C     and MIXED TYPE BOUNDARY CONDITIONS In Barotropic or external Mode
C

      IF (BCTYPE.EQ.'OCLAMP '.OR.BCTYPE.EQ.'IRANDB '.OR.
     +    BCTYPE.EQ.'MIXED  ') CALL OCLAMP(CLAM)
      IF (BCTYPE.EQ.'PCLAMP ') CALL PCLAMP
C
      CALL ADVAVE(ADVUA,ADVVA)
C
      DO 95 J=1,JM
      DO 95 I=1,IM
 95   CURV42D(I,J)=.25*COR(I,J)
C
C-----------------------------------------------------------------------
C        COMPUTE BOTTOM FRICTION AND CURVATURE TERMS 
C            IF RUN IS BAROTROPIC MODE ONLY
C-----------------------------------------------------------------------
C
      IF(TOR.NE.'BAROTROPIC') GO TO 5000
C
      DO 60 J=2,JMM1
      DO 60 I=2,IMM1
      CURV42D(I,J)= CURV42D(I,J)
     .  +.125*((VA(I,J+1)+VA(I,J))*
     .     ((H2(I+1,J)*FSM(I+1,J)+H2(I,J)*FSM(I,J))/
     .          (FSM(I+1,J)+FSM(I,J)+1.E-30)-
     .      (H2(I,J)*FSM(I,J)+H2(I-1,J)*FSM(I-1,J))/
     .          (FSM(I,J)+FSM(I-1,J)+1.E-30)) 
     .  -(UA(I+1,J)+UA(I,J))*
     .     ((H1(I,J+1)*FSM(I,J+1)+H1(I,J)*FSM(I,J))/
     .          (FSM(I,J+1)+FSM(I,J)+1.E-30)-
     .      (H1(I,J)*FSM(I,J)+H1(I,J-1)*FSM(I,J-1))/
     .          (FSM(I,J)+FSM(I,J-1)+1.E-30)) )
     .                /ART(I,J)
 60   CONTINUE
C
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      WUBOT(I,J)=-BFRIC
     .     *SQRT(UAB(I,J)**2+(.25*(VAB(I,J)
     .     +VAB(I,J+1)+VAB(I-1,J)+VAB(I-1,J+1)))**2)*UAB(I,J)
 100  CONTINUE
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
      WVBOT(I,J)=-BFRIC
     .     *SQRT((.25*(UAB(I,J)+UAB(I+1,J)
     .     +UAB(I,J-1)+UAB(I+1,J-1)))**2+VAB(I,J)**2)*VAB(I,J)
 102  CONTINUE
 5000 CONTINUE
C
C-----------------------------------------------------------------------
C     CALCULATE U COMPONENT OF VELOCITY
C-----------------------------------------------------------------------
C
      DO 420 J=2,JMM1
      DO 420 I=3,IMM1
 420  UAF(I,J)=ADVUA(I,J)
     .    -ARU(I,J)*( CURV42D(I,J)*D(I,J)*(VA(I,J+1)+VA(I,J))
     .              +CURV42D(I-1,J)*D(I-1,J)*(VA(I-1,J+1)+VA(I-1,J)) )
     .         +.25*GRAV*(H2(I,J)+H2(I-1,J))*(D(I,J)+D(I-1,J))
     .             *( (1.-2.*THETA)*(EL(I,J)-EL(I-1,J))
     .            +THETA*(ELB(I,J)-ELB(I-1,J)+ELF(I,J)-ELF(I-1,J)) )
     .         +.25*GRAV*(H2(I,J)+H2(I-1,J))*(D(I,J)+D(I-1,J))
     .             *(DATUM(I,J)-DATUM(I-1,J))*RAMP
     .         +RAMP*TRNU(I,J)
     .      -ARU(I,J)*(-.5*(WUSURF(I,J)+WUSURF(I-1,J))+WUBOT(I,J)  )
     *            +ARU(I,J)*(D(I,J)+D(I-1,J))*(PATM(I,J)-PATM(I-1,J))
     *             /RHO0/(H1(I,J)+H1(I-1,J))

      DO 425 J=2,JMM1
      DO 425 I=3,IMM1
 425  UAF(I,J)=
     .         ((H(I,J)+ELB(I,J)+H(I-1,J)+ELB(I-1,J))*ARU(I,J)*UAB(I,J)
     .                -4.*DTE*UAF(I,J))
     .        /((H(I,J)+ELF(I,J)+H(I-1,J)+ELF(I-1,J))*ARU(I,J))
C
C-----------------------------------------------------------------------
C     CALCULATE V COMPONENT OF VELOCITY
C-----------------------------------------------------------------------
C
      DO 430 J=3,JMM1
      DO 430 I=2,IMM1
 430  VAF(I,J)=ADVVA(I,J)  
     .    +ARV(I,J)*( CURV42D(I,J)*D(I,J)*(UA(I+1,J)+UA(I,J))
     .                +CURV42D(I,J-1)*D(I,J-1)*(UA(I+1,J-1)+UA(I,J-1)) )
     .           +.25*GRAV*(H1(I,J)+H1(I,J-1))*(D(I,J)+D(I,J-1))
     .                *( (1.-2.*THETA)*(EL(I,J)-EL(I,J-1))
     .                +THETA*(ELB(I,J)-ELB(I,J-1)+ELF(I,J)-ELF(I,J-1)) )
     .           +.25*GRAV*(H1(I,J)+H1(I,J-1))*(D(I,J)+D(I,J-1))
     .                *(DATUM(I,J)-DATUM(I,J-1))*RAMP
     .           +RAMP*TRNV(I,J)
     .    + ARV(I,J)*( .5*(WVSURF(I,J)+WVSURF(I,J-1))-WVBOT(I,J)   )
     *            +ARV(I,J)*(D(I,J)+D(I,J-1))*(PATM(I,J)-PATM(I,J-1))
     *             /RHO0/(H2(I,J)+H2(I,J-1))
      DO 435 J=3,JMM1
      DO 435 I=2,IMM1
 435  VAF(I,J)=
     .        ((H(I,J)+ELB(I,J)+H(I,J-1)+ELB(I,J-1))*VAB(I,J)*ARV(I,J)
     .              -4.*DTE*VAF(I,J))
     .       /((H(I,J)+ELF(I,J)+H(I,J-1)+ELF(I,J-1))*ARV(I,J))
C
      IF(BCTYPE.EQ.'RANDB  ') CALL RANDB
      IF(BCTYPE.EQ.'MIXED  ') CALL MIXED(1)

      CALL BCOND(2,DTI2,0)
C
      RETURN
      END
