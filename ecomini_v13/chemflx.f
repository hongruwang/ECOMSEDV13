      SUBROUTINE CHEMFLX
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
C  CALCULATES DEPOSITION AND DEPOSITION FLUXES OF CHEM AT
C  SEDIMENT-WATER INTERFACE
C
C**********************************************************************
C
      INCLUDE 'comdeck'
C
C
C  CLASS 1 & 2 SEDIMENTS:  FLUX IN ug CHEM/cm**2
C
C  CBEDCHEM = BED CONC. IN LAYER 1 (ug CHEM/g SOLIDS)
C
C  COHESIVE ELEMENTS
C
      DO 10 J=2,JMM1
        DO 10 I=2,IMM1
          CHEMBOT1(I,J)=0.0
          CHEMBOT2(I,J)=0.0
C
          IF (FSM(I,J).GT.0.0) THEN
C
C  EROSION FROM LAYER 1
C
            CHEMBOT1(I,J)=CBEDCHEM(1,I,J)*E(1,I,J)
            CHEMBOT2(I,J)=CBEDCHEM(1,I,J)*E(2,I,J)
C
C  DEPOSITION
C
            IF (CSED1(I,J,KBM1).GT.0.0) THEN
              CHEMBOT1(I,J)=CHEMBOT1(I,J)-
     +                   DD(1,I,J)*CHEM1(I,J,KBM1)/CSED1(I,J,KBM1)
            ENDIF
C
            IF (CSED2(I,J,KBM1).GT.0.0) THEN
              CHEMBOT2(I,J)=CHEMBOT2(I,J)-
     +                   DD(2,I,J)*CHEM2(I,J,KBM1)/CSED2(I,J,KBM1)
            ENDIF
          ENDIF
 10   CONTINUE
C
      RETURN
      END
