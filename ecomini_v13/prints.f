      SUBROUTINE PRINTS(DRHOX,DRHOY,TRNU,TRNV)
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
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION IVAR(IM,JM),PRT(IM,KB),BMASK(IM,JM)
      EQUIVALENCE (IVAR,C),(PRT,A)
      kp1=1
      kp2=2
C-------- KP1=SURFACE LAYER, KP2=NEXT LOWER LEVEL ----------------------
C
      IF(INT.EQ.0) THEN
C
C-------- PRINT INITIAL FIELDS -----------------------------------------
       WRITE(IUPRT,900)
       CALL PRINT  (H,FSM,IM,JM,IVAR,10.0,IUPRT,DEV)
C
      DO 2200 I=1,IM
        DO 2200 J=1,JM
          BMASK(I,J)=FLOAT(IBMSK(I,J))
 2200 CONTINUE
C
       WRITE(IUPRT,9011)
       CALL PRINT  (BMASK,FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
C
       WRITE(IUPRT,913)
       CALL PRINT(FSM,FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
       WRITE(IUPRT,914)
       CALL PRINT(DUM,DUM,IM,JM,IVAR,1.0,IUPRT,DEV)
       WRITE(IUPRT,915)
       CALL PRINT(DVM,DVM,IM,JM,IVAR,1.0,IUPRT,DEV)
       WRITE(IUPRT,901)
       CALL PRINT (H1,FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
       WRITE(IUPRT,902)
       CALL PRINT (H2,FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
       WRITE(IUPRT,903)
       CALL PRINT(ANG,FSM,IM,JM,IVAR,1.0E2,IUPRT,DEV)
       WRITE(IUPRT,904)
       CALL PRINT(COR,FSM,IM,JM,IVAR,1.0E8,IUPRT,DEV)
       WRITE(IUPRT,905)
       CALL PRINT(CBC,FSM,IM,JM,IVAR,1.0E4,IUPRT,DEV)
       WRITE(IUPRT,953)
       CALL PRTXY  (S,FSM,IM,JM,KB,KP1,IVAR,1.0E2,IUPRT,DEV)
       WRITE(IUPRT,954)
       CALL PRTXY  (T,FSM,IM,JM,KB,KP1,IVAR,1.0E2,IUPRT,DEV)
C
C  FOR CONSERVATIVE TRACER
C
      IF (TRACER.EQ.'INCLUDE') THEN
        WRITE(IUPRT,1954)
        CALL PRTXY  (CONC1,FSM,IM,JM,KB,KP1,IVAR,1.0E4,IUPRT,DEV)
      ENDIF
C
C-------- COMPUTE TIME STEP RESTRICTIONS -------------------------------
       DO 500 J=1,JM
       DO 500 I=1,IM
 500   A(I,J,1)=.5/SQRT(GRAV*H(I,J)*(1./H1(I,J)**2
     2   +1./H2(I,J)**2))*FSM(I,J)
C
       WRITE(IUPRT,502)
       CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,10.0,IUPRT,DEV)
C
C------- SPLIT INTO TWO DIRECTIONS FOR CHANNEL MODEL -------------------
       DO 510 J=1,JM
       DO 510 I=1,IM
 510   A(I,J,1)=.5/SQRT(GRAV*H(I,J))*H1(I,J)*FSM(I,J)
C
       WRITE(IUPRT,512)
       CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
C
       DO 520 J=1,JM
       DO 520 I=1,IM
 520   A(I,J,1)=.5/SQRT(GRAV*H(I,J))*H2(I,J)*FSM(I,J)
C
       WRITE(IUPRT,522)
       CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
C
C-------- INTERNAL MODE ------------------------------------------------
       DO 530 J=1,JM
       DO 530 I=1,IM
 530   A(I,J,1)=.5/(1.75*SQRT(1./H1(I,J)**2
     2   +1./H2(I,J)**2))*FSM(I,J)
C
       WRITE(IUPRT,532)
       CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,10.0,IUPRT,DEV)
C
       DO 540 J=1,JM
       DO 540 I=1,IM
 540   A(I,J,1)=1./1.75*H1(I,J)*FSM(I,J)
C
       WRITE(IUPRT,542)
       CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
C
       DO 550 J=1,JM
       DO 550 I=1,IM
 550   A(I,J,1)=1./1.75*H2(I,J)*FSM(I,J)
C
       WRITE(IUPRT,552)
       CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,1.0,IUPRT,DEV)
C
C-----------------------------------------------------------------------
      ELSE
C
       WRITE(IUPRT,910) TIME 
       CALL PRINT(EL,FSM,IM,JM,IVAR,1.0E3,IUPRT,DEV)
       WRITE(IUPRT,920) TIME 
       CALL PRINT(UA,DUM,IM,JM,IVAR,1.0E2,IUPRT,DEV) 
       WRITE(IUPRT,930) TIME 
       CALL PRINT(VA,DVM,IM,JM,IVAR,1.0E2,IUPRT,DEV)
C
       IF(TOR.EQ.'BAROTROPIC') RETURN
C
       IF(PTU.EQ.'Y') THEN
        WRITE(IUPRT,935) TIME 
        CALL PRTXY(U,DUM,IM,JM,KB,KP1,IVAR,1.0E2,IUPRT,DEV)
        IF(VSX.EQ.'Y') CALL SLICEXZ(U,DUM,IM,JM,KB,JROW,PRT,1.0E2,IUPRT)
        IF(VSY.EQ.'Y') CALL SLICEYZ(U,DUM,IM,JM,KB,IROW,PRT,1.0E2,IUPRT)
       ENDIF
C
       IF(PTV.EQ.'Y') THEN
        WRITE(IUPRT,940) TIME 
        CALL PRTXY (V,DVM,IM,JM,KB,KP1,IVAR,1.0E2,IUPRT,DEV) 
        IF(VSX.EQ.'Y') CALL SLICEXZ(V,DVM,IM,JM,KB,JROW,PRT,1.0E2,IUPRT)
        IF(VSY.EQ.'Y') CALL SLICEYZ(V,DVM,IM,JM,KB,IROW,PRT,1.0E2,IUPRT)
       ENDIF
C
       IF(PTW.EQ.'Y') THEN
        WRITE(IUPRT,945) TIME
        CALL PRTXY (WR,FSM,IM,JM,KB,KP2,IVAR,1.0E4,IUPRT,DEV) 
        IF(VSX.EQ.'Y')
     .     CALL SLICEXZ(WR,FSM,IM,JM,KB,JROW,PRT,1.0E4,IUPRT)
        IF(VSY.EQ.'Y')
     .     CALL SLICEYZ(WR,FSM,IM,JM,KB,IROW,PRT,1.0E4,IUPRT)
       ENDIF
C
       IF(PTAM.EQ.'Y') THEN
        WRITE(IUPRT,985) TIME 
        CALL PRTXY (AAM,FSM,IM,JM,KB,KP1,IVAR,1.0E2,IUPRT,DEV)
        IF(VSX.EQ.'Y') CALL SLICEXZ(AAM,FSM,IM,JM,KB,JROW,PRT,1.0,IUPRT)
        IF(VSY.EQ.'Y') CALL SLICEYZ(AAM,FSM,IM,JM,KB,IROW,PRT,1.0,IUPRT)
       ENDIF
C
C  FOR CONSERVATIVE TRACER
C
         IF (TRACER.EQ.'INCLUDE') THEN
           WRITE(IUPRT,1951) TIME 
           CALL PRTXY (CONC1,FSM,IM,JM,KB,KP1,IVAR,1.0E4,IUPRT,DEV)
           IF(VSX.EQ.'Y') 
     +           CALL SLICEXZ(CONC1,FSM,IM,JM,KB,JROW,PRT,1.0E4,IUPRT)
           IF(VSY.EQ.'Y') 
     +           CALL SLICEYZ(CONC1,FSM,IM,JM,KB,IROW,PRT,1.0E4,IUPRT)
         ENDIF
C
       IF (SEDTRAN.EQ.'INCLUDE') THEN
         IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
           WRITE (IUPRT,1900)TIME
           CALL PRTXY(CSED1,FSM,IM,JM,KB,KP1,IVAR,1.0E6,IUPRT,DEV)
           IF (VSX.EQ.'Y') THEN
             CALL SLICEXZ(CSED1,FSM,IM,JM,KB,JROW,PRT,1.0E6,IUPRT)
           ENDIF
           IF (VSY.EQ.'Y') THEN 
             CALL SLICEYZ(CSED1,FSM,IM,JM,KB,JROW,PRT,1.0E6,IUPRT)
           ENDIF
         ENDIF
C
         IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
           WRITE (IUPRT,1910)TIME
           CALL PRTXY(CSED2,FSM,IM,JM,KB,KP1,IVAR,1.0E6,IUPRT,DEV)
           IF (VSX.EQ.'Y') THEN
             CALL SLICEXZ(CSED2,FSM,IM,JM,KB,JROW,PRT,1.0E6,IUPRT)
           ENDIF
           IF (VSY.EQ.'Y') THEN
             CALL SLICEYZ(CSED2,FSM,IM,JM,KB,JROW,PRT,1.0E6,IUPRT)
           ENDIF
         ENDIF
C
C  PRINTS OUTPUT FOR SEDIMENT BED THICKNESS
C
C  DETERMINE INITIAL SED. THICKNESS
C 
C
 9500  FORMAT (I3,2X,20F5.2)
 9501  FORMAT (I3,2X,20F6.2)
 9502  FORMAT (I3,2X,20F6.3)
 9510  FORMAT (20F6.3)
 9940  FORMAT (I4,15F5.1)
C
C  CALC. FINAL SEDIMENT THICKNESS, INITIAL SEDIMENT-WATER INTERFACE
C  IS AT ZERO
C
      DO 8120 J=1,JM
       DO 8120 I=1,IM
C
        A(I,J,1)=0.0
C
        IF (D(I,J).LE.0.0) GOTO 8120
C
        TSET0T=0.0
        DO 8110 LL=1,LAYMAX
	  TSET0T=TSET0T+TSED0(LL,I,J)
 8110   CONTINUE
C
	TSEDT(I,J)=0.0
	DO 8130 LL=1,LAYMAX
	  TSEDT(I,J)=TSEDT(I,J)+TSED(LL,I,J)
 8130   CONTINUE
C
        IF (IBMSK(I,J).EQ.0) THEN
	  THICK(I,J)=TSEDT(I,J)-TSET0T
C
C  CONVERT FROM g/cm**2 TO cm
C
          THICK(I,J)=THICK(I,J)/CBED(I,J)
C
          IF (SEDTYPE.EQ.'SAND') THICK(I,J)=0.0
        ELSE
          IF (IBMSK(I,J).EQ.1) THEN
C
	    THICK(I,J)=BEDTH(1,I,J)/((CBED(I,J)/2.65)*H1(I,J)*H2(I,J))
     +                       -BEDTHI
C
C  CONVERT FROM m TO cm
C
            THICK(I,J)=100.*THICK(I,J)
C
            IF (SEDTYPE.EQ.'MUD ') THICK(I,J)=0.0
          ELSE
            THICK(I,J)=0.0
          ENDIF
        ENDIF
C
C  CONVERT THICKNESS FROM cm TO mm
C
         A(I,J,1)=THICK(I,J)*10.
 8120  CONTINUE
C
C  OUTPUT AS mm
C
        WRITE (IUPRT,9140)TIME
        CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,1.0E0,IUPRT,DEV)
C
C  OUTPUT BOTTOM SHEAR STRESS
C
      DO 8140 J=1,JM
       DO 8140 I=1,IM
        A(I,J,1)=0.0
        IF (D(I,J).LE.0.0) GOTO 8140
        A(I,J,1)=TAU(I,J,KB)
 8140 CONTINUE
C
C  OUTPUT AS dyne/cm**2
C
        WRITE (IUPRT,9150)TIME
        CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,1.0E0,IUPRT,DEV)
        ENDIF
C
C******************************************************************
C
       IF (CHEMTRAN.EQ.'INCLUDE') THEN
         IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
            Write (IUPRT,7203) TIME
            Call PRTXY(CHEM1,FSM,IM,JM,KB,KP1,IVAR,1.0E3,IUPRT,DEV)
            If (VSX.EQ.'Y') Call 
     *        SLICEXZ(CHEM1,FSM,IM,JM,KB,JROW,PRT,1.0E3,IUPRT)
            If (VSY.EQ.'Y') Call 
     *        SLICEYZ(CHEM1,FSM,IM,JM,KB,IROW,PRT,1.0E3,IUPRT)
         ENDIF
C
         IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
            Write (IUPRT,7204) TIME
            Call PRTXY(CHEM2,FSM,IM,JM,KB,KP1,IVAR,1.0E3,IUPRT,DEV)
            If (VSX.EQ.'Y') Call 
     *        SLICEXZ(CHEM2,FSM,IM,JM,KB,JROW,PRT,1.0E3,IUPRT)
            If (VSY.EQ.'Y') Call 
     *        SLICEYZ(CHEM2,FSM,IM,JM,KB,IROW,PRT,1.0E3,IUPRT)
         ENDIF
C
C
         DO 8125 J=1,JM
           DO 8125 I=1,IM
             A(I,J,1)=0.0
             IF (D(I,J).LE.0.0) GOTO 8125
             A(I,J,1)=CBEDCHEM(1,I,J)
 8125    CONTINUE
C
C  OUTPUT TOP LAYER CONCENTRATION (UNITS ARE ppm)
C
         WRITE (IUPRT,9141)TIME
         CALL PRINT(A(1,1,1),FSM,IM,JM,IVAR,1.0E0,IUPRT,DEV)
C
 9651   CONTINUE
          ENDIF
C
       IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'TEMP_ONLY '.OR.
     .    TOR.EQ.'SALT_ONLY ') THEN
C
        IF(PTS.EQ.'Y') THEN
         WRITE(IUPRT,951) TIME 
         CALL PRTXY (S,FSM,IM,JM,KB,KP1,IVAR,1.0E2,IUPRT,DEV)
         IF(VSX.EQ.'Y') CALL SLICEXZ(S,FSM,IM,JM,KB,JROW,PRT,1.0,IUPRT)
         IF(VSY.EQ.'Y') CALL SLICEYZ(S,FSM,IM,JM,KB,IROW,PRT,1.0,IUPRT)
        ENDIF
C
        IF(PTT.EQ.'Y') THEN
         WRITE(IUPRT,950) TIME 
         CALL PRTXY (T,FSM,IM,JM,KB,KP1,IVAR,1.0E2,IUPRT,DEV)
         IF(VSX.EQ.'Y') CALL SLICEXZ(T,FSM,IM,JM,KB,JROW,PRT,1.0,IUPRT)
         IF(VSY.EQ.'Y') CALL SLICEYZ(T,FSM,IM,JM,KB,IROW,PRT,1.0,IUPRT)
        ENDIF
C
        IF(PRHO.EQ.'Y') THEN
         WRITE(IUPRT,952) TIME 
         CALL PRTXY (RHO,FSM,IM,JM,KB,KP1,IVAR,1.0E5,IUPRT,DEV)
         IF(VSX.EQ.'Y')
     .      CALL SLICEXZ(RHO,FSM,IM,JM,KB,JROW,PRT,1.0E3,IUPRT)
         IF(VSY.EQ.'Y')
     .      CALL SLICEYZ(RHO,FSM,IM,JM,KB,IROW,PRT,1.0E3,IUPRT)
        ENDIF
       ENDIF
C
       IF(VERTMIX.EQ.'CLOSURE   ') THEN
C
        IF(PTQ2.EQ.'Y') THEN
         WRITE(IUPRT,955) TIME 
         CALL PRTXY (Q2,FSM,IM,JM,KB,KP2,IVAR,1.0E6,IUPRT,DEV)
         IF(VSX.EQ.'Y')
     .      CALL SLICEXZ(Q2,FSM,IM,JM,KB,JROW,PRT,1.0E6,IUPRT)
         IF(VSY.EQ.'Y')
     .      CALL SLICEYZ(Q2,FSM,IM,JM,KB,IROW,PRT,1.0E6,IUPRT)
        ENDIF
C
        IF(PTL.EQ.'Y') THEN
         WRITE(IUPRT,960) TIME 
         CALL PRTXY (L,FSM,IM,JM,KB,KP2,IVAR,1.0E2,IUPRT,DEV)
         IF(VSX.EQ.'Y') CALL SLICEXZ(L,FSM,IM,JM,KB,JROW,PRT,1.0,IUPRT)
         IF(VSY.EQ.'Y') CALL SLICEYZ(L,FSM,IM,JM,KB,IROW,PRT,1.0,IUPRT)
        ENDIF
C
        IF(PTKM.EQ.'Y') THEN
         WRITE(IUPRT,965) TIME 
         CALL PRTXY (KM,FSM,IM,JM,KB,KP2,IVAR,1.0E4,IUPRT,DEV)
         IF(VSX.EQ.'Y')
     .      CALL SLICEXZ(KM,FSM,IM,JM,KB,JROW,PRT,1.0E4,IUPRT)
         IF(VSY.EQ.'Y')
     .      CALL SLICEYZ(KM,FSM,IM,JM,KB,IROW,PRT,1.0E4,IUPRT)
        ENDIF
C
        IF(PTKH.EQ.'Y') THEN
         WRITE(IUPRT,970) TIME 
         CALL PRTXY (KH,FSM,IM,JM,KB,KP2,IVAR,1.0E4,IUPRT,DEV)
         IF(VSX.EQ.'Y')
     .      CALL SLICEXZ(KH,FSM,IM,JM,KB,JROW,PRT,1.0E4,IUPRT)
         IF(VSY.EQ.'Y')
     .      CALL SLICEYZ(KH,FSM,IM,JM,KB,IROW,PRT,1.0E4,IUPRT)
        ENDIF
C
       ENDIF
C 
      ENDIF
C 
 900  FORMAT(///' .......... DEPTH (m) ..........'/)
 9011 FORMAT(///' .......... SEDIMENT BED MASK ..........'/)
 901  FORMAT(///' .......... METRIC H1 (m) ..........'/)
 902  FORMAT(///' .......... METRIC H2 (m) ..........'/)
 903  FORMAT(///' .......... ANGLE (radians) ..........'/)
 904  FORMAT(///' .......... CORIOLIS PARAMETER (1/sec) ..........'/)
 905  FORMAT(///' .......... BOTTOM DRAG COEFFICIENT ..........'/)
 913  FORMAT(///' .......... ELEVATION MASK ..........'/)
 914  FORMAT(///' .......... U  VEL MASK ..........'/)
 915  FORMAT(///' .......... V  VEL MASK ..........'/)
 953  FORMAT(///' ..... SALINITY (ppt) INITIALIZATION .....'/)
 954  FORMAT(///' ..... TEMPERATURE (degC) INITIALIZATION .....'/)
 1954 FORMAT(///' ..... CONSERVATIVE TRACER INITIALIZATION .....'/)
C
 502  FORMAT(///' MAXIMUM TIME STEP(sec) EXTERNAL MODE'/)
 512  FORMAT(///' MAXIMUM TIME STEP(sec) EXTERNAL MODE XI 1 DIRECTION'/)
 522  FORMAT(///' MAXIMUM TIME STEP(sec) EXTERNAL MODE XI 2 DIRECTION'/)
 532  FORMAT(///' MAXIMUM TIME STEP(sec) INTERNAL MODE'/)
 542  FORMAT(///' MAXIMUM TIME STEP(sec) INTERNAL MODE XI 1 DIRECTION'/)
 552  FORMAT(///' MAXIMUM TIME STEP(sec) INTERNAL MODE XI 2 DIRECTION'/)
C
 910  FORMAT(///' SURFACE ELEVATION (m) '    ,F9.4,' DAYS AFTER T=0 '/)
 920  FORMAT(///' AVERAGED U VELOCITY (m/s) ',F9.4,' DAYS AFTER T=0 '/)
 930  FORMAT(///' AVERAGED V VELOCITY (m/s) ',F9.4,' DAYS AFTER T=0 '/)
 933  FORMAT(///' U-INTEGRALS '              ,F9.4,' DAYS AFTER T=0 '/)
 934  FORMAT(///' V-INTEGRALS '              ,F9.4,' DAYS AFTER T=0 '/)
 935  FORMAT(///' U VELOCITY(m/s) '          ,F9.4,' DAYS AFTER T=0 '/) 
 940  FORMAT(///' V VELOCITY(m/s) '          ,F9.4,' DAYS AFTER T=0 '/)
 945  FORMAT(///' W VELOCITY(m/s) '          ,F9.4,' DAYS AFTER T=0 '/)
 950  FORMAT(///' TEMPERATURE (degC) '       ,F9.4,' DAYS AFTER T=0 '/)
 951  FORMAT(///' SALINITY (ppt) '           ,F9.4,' DAYS AFTER T=0 '/)
1951  FORMAT(///' DISSOLVED TRACER '         ,F9.4,' DAYS AFTER T=0 '/)
 952  FORMAT(///' DENSITY - 1.0 (gm/cm**3) ' ,F9.4,' DAYS AFTER T=0 '/)
 955  FORMAT(///' TURBULENT K.E. '           ,F9.4,' DAYS AFTER T=0 '/)
 960  FORMAT(///' MIXING LENGTH (m) '        ,F9.4,' DAYS AFTER T=0 '/)
 965  FORMAT(///' MIXING KM (m**2/s) '       ,F9.4,' DAYS AFTER T=0 '/)
 970  FORMAT(///' MIXING KH (m**2/s) '       ,F9.4,' DAYS AFTER T=0 '/)
 971  FORMAT(///' BAROCLINIC PRESSURE DRHOX ',F9.4,' DAYS AFTER T=0 '/)
 972  FORMAT(///' BAROCLINIC PRESSURE DRHOY ',F9.4,' DAYS AFTER T=0 '/)
 975  FORMAT(///' U BOTTOM STRESS ((m/s)**2)',F9.4,' DAYS AFTER T=0 '/)
 980  FORMAT(///' V BOTTOM STRESS ((m/s)**2)',F9.4,' DAYS AFTER T=0 '/)
 985  FORMAT(///' HORIZONTAL MIXING '        ,F9.4,' DAYS AFTER T=0 '/)
C
 7203 Format (///' COHESIVE PARTICLE-BOUND CONCENTRATION  ',F9.4,
     +             ' DAYS AFTER T=0 '/
     +               10X,'(IN ug/l AS SHOWN)')
 7204 Format (///' NON-COHESIVE PARTICLE-BOUND CONCENTRATION  ',F9.4,
     +               ' DAYS AFTER T=0 '/
     +               10X,'(IN ug/l AS SHOWN)')
 1900 Format (///' COHESIVE SEDIMENT CONCENTRATION  ',F9.4,
     +             ' DAYS AFTER T=0 '/
     +               10X,'(IN mg/l AS SHOWN)')
 1910 Format (///' NON-COHESIVE SEDIMENT CONCENTRATION  ',F9.4,
     +             ' DAYS AFTER T=0 '/
     +               10X,'(IN mg/l AS SHOWN)')
 9140 Format (///' SEDIMENT BED ELEVATION CHANGE  ',F9.4,
     +             ' DAYS AFTER T=0 '/
     +               10X,'(IN mm AS SHOWN)')
 9141 Format (///' PARTICLE-BOUND TRACER, BED CONCENTRATION  ',F9.4,
     +             ' DAYS AFTER T=0 '/
     +               10X,'(IN ppm AS SHOWN, SURFACE LAYER)')
 9150 Format (///' BOTTOM SHEAR STRESS ',F9.4,
     +             ' DAYS AFTER T=0 '/
     +               10X,'(IN dyne/cm**2 AS SHOWN)')
 9503 FORMAT (I3,2X,20F6.1)
C
C******************************************************************
C
      RETURN
      END
