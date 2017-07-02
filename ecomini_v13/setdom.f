      SUBROUTINE SETDOM(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
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
C
      INCLUDE 'comdeck'
C
C  FOR WET GRID INPUT FROM POM
C	
      CHARACTER*10 RESTAR
      COMMON /COAST1/ICNT,INDX(MAXWET),JNDX(MAXWET),RESTAR
C
      SAVE
C
C-----------------------------------------------------------------------
C     INITIAL CONDITION PREPARATION PROGRAM FOR A GENERAL CIRCULATION
C        MODEL IN ORTHOGONAL CURVILINEAR COORDINATES
C-----------------------------------------------------------------------
C
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION IVAR(IM,JM),PRT(IM,KB)
      EQUIVALENCE (IVAR,TPS),(PRT,A)
      DIMENSION COM(80)
C 
C-------- ALL UNITS IN M.K.S. SYSTEM -----------------------------------
      ISTART=1
      IEND=1
      INT=0 
C
C  FOR EXTERNAL HYDRODYNAMICS RUN
C
C  ASSUME POM OUTPUT FOR WET GRID
C
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
        IF(RESTAR.EQ.'HOT START '.AND.IWET.EQ.1) GOTO 153
        OPEN (IUTRN,FILE='gcm_geom',FORM='UNFORMATTED')
          READ (IUTRN)DZ,DZZ
          READ (IUTRN)H,H1,H2,TPS
          READ (IUTRN) ANG,NU
          DO J=1,JM
           DO I=1,IM
             ANG(I,J)=ANG(I,J)*360.0/(2.*3.14159265)
           ENDDO
          ENDDO
        CLOSE(IUTRN)
C
        OPEN (IUTRN,FORM='unformatted',FILE='wet_grid')
          READ (IUTRN)ICNT
        DO I=1,ICNT
          READ(IUTRN)INDX(I),JNDX(I)
        END DO
        CLOSE (IUTRN)
C
        DO 151 K=1,KBM1
          DZR(K)=1./DZ(K)
 151    CONTINUE
C
C  RESET LAND ELEMENTS TO ZERO DEPTH (SET TO 1.000 m IN POM)
C
        DO 152 J=1,JM
          DO 152 I=1,IM
            IF (H(I,J).LT.1.002) H(I,J)=0.0
 152    CONTINUE
        GOTO 999
C
c     Reads wet_grid when "HOT START" and IWET=1 i.e. IWET on
C
 153    CONTINUE
        OPEN (IUTRN,FORM='unformatted',FILE='wet_grid')
          READ (IUTRN)ICNT
        DO I=1,ICNT
          READ(IUTRN)INDX(I),JNDX(I)
        END DO
        CLOSE (IUTRN)
        RETURN
      ENDIF
C
      OPEN(IUGRD,FILE='model_grid')
C
      WRITE(IUPRT,51)            
 51   FORMAT(/'...... MODEL STARTING UP FROM INITAL CONDITIONS .......')
C
C-------- ESTABLISH DEPTH ARRAY ----------------------------------------
      READ(IUGRD,11) (COM(I),I=1,80)
      WRITE(IUPRT,3) (COM(I),I=1,80)
 3    FORMAT(/1X,80A1/)
 11   FORMAT(80A1)
C
C-------- INITIALIZE SIGMA LEVELS --------------------------------------
      READ(IUGRD,11) (COM(I),I=1,80)
      WRITE(IUPRT,3) (COM(I),I=1,80)
C
      READ(IUGRD,4) IKB
    4 FORMAT(I5)
      WRITE(IUPRT,40) IKB
   40 FORMAT(' KB = ',I5,/)
      IF(IKB.GT.21) THEN
       WRITE(6,41) KBM1
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
   41 FORMAT(//' NUMBER OF MODEL LAYERS',I5,' (KBM1)'/
     .         '           GREATER THAN'/
     .         ' MAXIMUM NUMBER OF VERTICAL DISTRIBUTION LAYERS'/
     .         ' SPECIFIED IN run_data DISCHARGE INFORMATION'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
C
      IF(IKB.NE.KB) THEN
       WRITE(6,42) IKB,KB
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
   42 FORMAT(//' NUMBER OF SIGMA LEVELS IN MODEL_GRID',I5,' (IKB)'/
     .         '           NOT EQUAL TO'/
     .         ' NUMBER OF SIGMA LEVELS IN COMDECK   ',I5,' (KB)'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
C
      IF(IKB.LT.3 .AND. TOR.NE.'BAROTROPIC') THEN
       WRITE(6,43) IKB,TOR
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
   43 FORMAT(//' IKB/KB = ',I5,' NOT ALLOWED IN A ',A10,' RUN'/
     .         ' NUMBER OF SIGMA LEVELS MUST BE GREATER THAN 2'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
C
      DO 140 K=1,IKB
      READ(IUGRD,5) Z(K)
  140 CONTINUE
      WRITE(IUPRT,5) (Z(K),K=1,IKB)
    5 FORMAT(8F10.5)
C
      DO 150 K=1,KBM1
      DZ(K)=Z(K)-Z(K+1)
      DZR(K)=1./DZ(K)
      ZZ(K)=.5*(Z(K)+Z(K+1))
  150 CONTINUE
      DO 460 K=1,KBM2
      DZZ(K)=ZZ(K)-ZZ(K+1)
  460 CONTINUE
      DZZ(KBM1)=0.0
      DZ(KB)=0.0
C
C-------- DEFINE THE METRICS OF THE COORDINATE TRANSFORMATION ----------
      READ(IUGRD,11) (COM(I),I=1,80)
      WRITE(IUPRT,3) (COM(I),I=1,80)
C
      READ(IUGRD,7)   IIX, IJY
      WRITE(IUPRT,71) IIX, IJY
   7  FORMAT(2I5,6F10.2)
   91 FORMAT(2I5,4F10.2,2F10.5,f5.1)
c   91 FORMAT(2I4,8f9.0)
   71 FORMAT(' IM = ',I5,/' JM = ',I5)
      IF(IIX.NE.IM) THEN
       WRITE(6,8) IIX,IM
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
    8 FORMAT (//'     MODEL_GRID I-INDEX',I5,' (IIX)',/
     .          '        DOES NOT EQUAL'/
     .          '     COMDECK    I-INDEX',I5,' (IM)'/
     .          ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
      IF(IJY.NE.JM) THEN
       WRITE(6,9) IJY,JM
       !!!&&&CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
    9 FORMAT (//'     MODEL_GRID J-INDEX',I5,' (IJY)',/
     .          '        DOES NOT EQUAL'/
     .          '     COMDECK    J-INDEX',I5,' (JM)'/
     .          ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
C
C..........SET PARAMETERS TO INITIAL VALUES..........
      DO 50 J=1,JM
      DO 50 I=1,IM
      H1(I,J)    = 0.0
      H2(I,J)    = 0.0
      H(I,J)     = 0.0
      ANG(I,J)   = 0.0
      COR(I,J)   = 0.0
      DATUM(I,J) = 0.0
      YGRID(I,J) = 0.0
   50 CONTINUE
C
      NUMB=IIX*IJY
      DO 100 N=1,NUMB
      READ(IUGRD,91,END=101)
     .     I,J,H1(I,J),H2(I,J),H(I,J),ANG(I,J),COR(I,J),XGRID(I,J),
     .     DATUM(I,J)
      YGRID(I,J)=COR(I,J)
100   continue
      
  101 CONTINUE
      CLOSE(IUGRD)
C
C-------- DEFINE MASK FOR FREE SURFACE HEIGHT = FSM --------------------
C-------- DEFINE MASK FOR (U,V) VELOCITY = (DUM,DVM) -------------------
C-------- IF DEPTH <= 0.0, DEPTH AND MASK SET TO 0 ---------------------
 999  DO 400 J=1,JM
      DO 400 I=1,IM
      FSM(I,J)=1.0
      DUM(I,J)=1.0
      DVM(I,J)=1.0
      IF(H(I,J).LE.0.0) THEN
       H(I,J)=0.0+1.0E-10
       FSM(I,J)=0.0
       DUM(I,J)=0.0
       DVM(I,J)=0.0
      ENDIF
 400  CONTINUE 
C
      DO 420 J=1,JMM1
      DO 420 I=1,IM
      IF(FSM(I,J).EQ.0.0.AND.FSM(I,J+1).NE.0.0) DVM(I,J+1)=0.0
  420 CONTINUE
      DO 440 J=1,JM
      DO 440 I=1,IMM1
      IF(FSM(I,J).EQ.0.0.AND.FSM(I+1,J).NE.0.0) DUM(I+1,J)=0.0
  440 CONTINUE
C
      DO 160 J=1,JM
      DO 160 I=1,IM
      H1(I,J)=H1(I,J)+1.E-10
      H2(I,J)=H2(I,J)+1.E-10
      ANG(I,J)=ANG(I,J)*2.*3.14159265/360.
      COR(I,J)=2.*7.292E-5*SIN(COR(I,J)*2.*3.14159/360.)*FSM(I,J)
  160 CONTINUE
C
      DO 260 J=1,JM
      DO 260 I=1,IM
      D(I,J)=H(I,J)
  260 DT(I,J)=H(I,J)
C
      DO 340 J=2,JMM1
      DO 340 I=2,IMM1
      ART(I,J)=H1(I,J)*H2(I,J)
      ARU(I,J)=.25E0*(H1(I,J)+H1(I-1,J))*(H2(I,J)+H2(I-1,J))
      ARV(I,J)=.25E0*(H1(I,J)+H1(I,J-1))*(H2(I,J)+H2(I,J-1))
  340 CONTINUE
      DO 310 J=1,JM
      ARU(IM,J)=ARU(IMM1,J)
  310 ARU(1,J)=ARU(2,J)
      DO 320 I=1,IM
      ARV(I,JM)=ARV(I,JMM1)
  320 ARV(I,1)=ARV(I,2)
C
      DO 410 J=1,JM
      DO 410 I=1,IM
 410  TPS(I,J)=SQRT(H1(I,J)**2+H2(I,J)**2)
C 
      CALL MAXMIN(TPS,FSM,IM,JM,FMAX,FMIN)
C 
      IF(TOR.NE.'BAROTROPIC') THEN
      DO 350 K=1,KBM1
      DO 350 J=1,JM
      DO 350 I=1,IM
      AAM(I,J,K)=HORCON*TPS(I,J)/FMIN*FSM(I,J)
  350 CONTINUE
      DO 375 J=1,JM
      DO 375 I=1,IM
 375  AAM2D(I,J)=0.0
      DO 380 K=1,KBM1
      DO 380 J=1,JM
      DO 380 I=1,IM
      AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
 380  CONTINUE
C
      ELSE
C
      DO 355 J=1,JM
      DO 355 I=1,IM
      AAM2D(I,J)=HORCON*TPS(I,J)/FMIN*FSM(I,J)
  355 CONTINUE
C 
      ENDIF
C
      DO 360 J=1,JM
      DO 360 I=1,IM
      KM(I,J,1)=0.0
      KM(I,J,KB)=0.0
      KH(I,J,1)=0.0
      KH(I,J,KB)=0.0
      KQ(I,J,1)=0.0
      KQ(I,J,KB)=0.0
      L (I,J,1)=0.0
      L (I,J,KB)=0.0
      Q2(I,J,1)=0.0
      Q2(I,J,KB)=0.0
      Q2B(I,J,1)=0.0
      Q2B(I,J,KB)=0.0
      Q2L(I,J,1)=0.0
      Q2L(I,J,KB)=0.0
      Q2LB(I,J,1)=0.0
      Q2LB(I,J,KB)=0.0
  360 CONTINUE
C
      DO 370 K=2,KBM1
      DO 370 J=1,JM
      DO 370 I=1,IM
      L(I,J,K)   =1.*FSM(I,J)
      Q2B(I,J,K) =1.E-12*FSM(I,J)
      Q2(I,J,K)  =Q2B(I,J,K)
      Q2LB(I,J,K)=Q2B(I,J,K)*L(I,J,K)
      Q2L(I,J,K) =Q2LB(I,J,K)
      KM(I,J,K)  =UMOL*FSM(I,J)
      KQ(I,J,K)  =UMOL*FSM(I,J)
      KH(I,J,K)  =UMOL/VPRNU*FSM(I,J)
  370 CONTINUE
C
      IF(VERTMIX.EQ.'CONSTANT  ') UMOL=0.0
C
C*********************************************************************
C
C  FOR PARTICLE TRACKING
C
      IF (PARTICLE.EQ.'INCLUDE') THEN

C
C  OPEN CORNER LOCATIONS FILE
C
C  CORNER LOCATION CONVENTION:  XCOR(I,J) = x(i-1/2, j-1/2) 
C     (LOWER LEFT-HAND CORNER)  YCOR(I,J) = y(i-1/2, j-1/2)
C
        OPEN (UNIT=33,FILE='corner_loc',FORM='FORMATTED',STATUS='OLD')
C
        DO 1300 N=1,1000000
          READ (33,*,END=1301)I,J,XCOR(I,J),YCOR(I,J)
 1300   CONTINUE
 1301   CLOSE(33)
C
C  CALC. H1 AND H2 AT ELEMENT INTERFACES
C
 1310   DO 1320 J=1,JM
          DO 1320 I=1,IM
            H1P(I,J)=SQRT((XCOR(I+1,J)-XCOR(I,J))**2.+
     +                        (YCOR(I+1,J)-YCOR(I,J))**2.)
C
            H2P(I,J)=SQRT((XCOR(I,J+1)-XCOR(I,J))**2.+
     +                        (YCOR(I,J+1)-YCOR(I,J))**2.)
 1320   CONTINUE
      ENDIF
C
C*********************************************************************
C
C
      RETURN
      END 
