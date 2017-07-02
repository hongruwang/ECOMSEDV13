      SUBROUTINE BOTTAU
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
C  CALCULATES BOTTOM SHEAR STRESSES FOR USE IN SEDIMENT TRANSPORT
C  (SUBROUTINES sedflx AND suslod)
C
C**********************************************************************
C
      INCLUDE 'comdeck'
      SAVE
C
C  CALC. BOTTOM SHEAR STRESS
C
      IF(TOR.EQ.'BAROTROPIC') THEN
         DO 10 J=2,JMM1
           DO 10 I=2,IMM1
             UBAR(I,J,KBM1)=0.5*(UA(I,J)+UA(I+1,J))
             VBAR(I,J,KBM1)=0.5*(VA(I,J)+VA(I,J+1))
 10      CONTINUE
      ELSE
         DO 11 K=1,KBM1
           DO 11 J=2,JMM1
             DO 11 I=2,IMM1
               UBAR(I,J,K)=0.5*(U(I,J,K)+U(I+1,J,K))
               VBAR(I,J,K)=0.5*(V(I,J,K)+V(I,J+1,K))
 11      CONTINUE
      ENDIF
C
C
      DO 20 J=2,JMM1
        DO 20 I=2,IMM1
           QBAR(I,J)=SQRT(UBAR(I,J,KBM1)*UBAR(I,J,KBM1)+
     +                    VBAR(I,J,KBM1)*VBAR(I,J,KBM1))
 20   CONTINUE
C
C-------- VELOCITY BOUNDARY CONDITION ----------------------------------
C
C  NO WIND WAVES
C
      IF (WAVEDYN.EQ.'NEGLECT ') THEN
        DO 40 J=2,JMM1
          DO 40 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
C
C  BOTTOM STRESS (AT K= KB)
C
              TAU(I,J,KB)=10000.*CBC(I,J)*QBAR(I,J)*QBAR(I,J)
            ENDIF
 40     CONTINUE
C
C     FORCING TAU VALUE = 0 AT CELLS CONNECTED BY RIVERS (IC,JC)
C
        DO 120 N=1,NUMQBC
          IC=IQC(N)
          JC=JQC(N)
          TAU(IC,JC,KB) = 0.0
120     CONTINUE

      ELSE
C
C  INCLUDE WIND WAVE EFFECTS
C
        IF (WAVEDYN.EQ.'SMBMODEL') CALL WAVESMB
        IF (WAVEDYN.EQ.'DONMODEL') CALL WAVEDON
        CALL STRESS
      ENDIF
C
      RETURN
      END
