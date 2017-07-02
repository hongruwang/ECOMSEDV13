      SUBROUTINE ARCHIVE
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
C
C  FOR PARTICLE TRACKING
C
      REAL XOUTP(NPARTM),YOUTP(NPARTM),ZOUTP(NPARTM)
C
      SAVE
C
C-----------------------------------------------------------------------
C-------- COMPUTATIONAL HISTORY WRITES AND ACCUMULATION ----------------
C-----------------------------------------------------------------------
      DO 8300 JHIST=1,JHM
      IF(INT.LT.IHIST(JHIST,1).OR.INT.GT.IHIST(JHIST,2)) GO TO 8300
      IF(TOR.EQ.'PROGNOSTIC' .OR. TOR.EQ.'DIAGNOSTIC'.OR.
     .   TOR.EQ.'TEMP_ONLY ' .OR. TOR.EQ.'SALT_ONLY ') THEN
       DO 500 K=1,KB
       DO 500 J=2,JMM1
       DO 500 I=2,IMM1
       ARCU (I,J,K)=ARCU (I,J,K)+U (I,J,K)*DEI
       ARCV (I,J,K)=ARCV (I,J,K)+V (I,J,K)*DEI
       ARCUX(I,J,K)=ARCUX(I,J,K)+U (I,J,K)*DEI*
     2    0.5*(H(I,J)+ET(I,J)+H(I-1,J)+ET(I-1,J))*DZ(K)
       ARCVX(I,J,K)=ARCVX(I,J,K)+V (I,J,K)*DEI*
     2    0.5*(H(I,J)+ET(I,J)+H(I,J-1)+ET(I,J-1))*DZ(K)
       ARCS (I,J,K)=ARCS (I,J,K)+S (I,J,K)*DEI
       ARCT (I,J,K)=ARCT (I,J,K)+T (I,J,K)*DEI
       ARCW (I,J,K)=ARCW (I,J,K)+WR(I,J,K)*DEI
       ARCKH(I,J,K)=ARCKH(I,J,K)+KH(I,J,K)*DEI
C
C  CONSERVATIVE TRACER
C
       ARCC1(I,J,K)=ARCC1(I,J,K)+CONC1(I,J,K)*DEI
C
C  FOR SEDIMENT TRANSPORT
C
        ARCSED1(I,J,K)=ARCSED1(I,J,K)+CSED1(I,J,K)*DEI
        ARCSED2(I,J,K)=ARCSED2(I,J,K)+CSED2(I,J,K)*DEI
        ARCTAU(I,J,K)=ARCTAU(I,J,K)+TAU(I,J,K)*DEI
C
C  FOR CHEMICAL TRANSPORT
C
        ARCHEM1(I,J,K)=ARCHEM1(I,J,K)+CHEM1(I,J,K)*DEI
        ARCHEM2(I,J,K)=ARCHEM2(I,J,K)+CHEM2(I,J,K)*DEI
C
 500   CONTINUE
       DO 520 J=1,JM
         DO 520 I=1,IM
           ARCET (I,J)=ARCET (I,J)+ET (I,J)*DEI
 520   CONTINUE
C
        IF (SEDTRAN.EQ.'INCLUDE') THEN
          DO 521 J=1,JM
            DO 521 I=1,IM
              IF (IBMSK(I,J).EQ.0) THEN
                DO 1521 LL=1,LAYMAX
                  ARCTHIK(I,J)=ARCTHIK(I,J)+TSED(LL,I,J)*DEI
 1521           CONTINUE
              ELSE
                IF (IBMSK(I,J).EQ.1) THEN
                  ARCTHIK(I,J)=ARCTHIK(I,J)+BEDTH(1,I,J)*DEI
                ELSE
                  ARCTHIK(I,J)=0.0
                ENDIF
              ENDIF
 521      CONTINUE
C
C
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
            DO 522 J=1,JM
              DO 522 I=1,IM
                DO 522 LL=1,NCHEMLAY
                    ARCPBED(LL,I,J)=ARCPBED(LL,I,J)
     +                                  +CBEDCHEM(LL,I,J)*DEI
 522        CONTINUE
          ENDIF
        ENDIF
C
      ELSE
       DO 530 J=1,JM
         DO 530 I=1,IM
           ARCUX(I,J,1)=ARCUX(I,J,1)+UAF(I,J)*DEI*
     2       0.5*(H(I,J)+EL(I,J)+H(I-1,J)+EL(I-1,J))
           ARCVX(I,J,1)=ARCVX(I,J,1)+VAF(I,J)*DEI*
     2       0.5*(H(I,J)+EL(I,J)+H(I,J-1)+EL(I,J-1))
           ARCU (I,J,1)=ARCU (I,J,1)+UAF(I,J)*DEI 
           ARCV (I,J,1)=ARCV (I,J,1)+VAF(I,J)*DEI 
           ARCET (I,J)=ARCET (I,J)+EL (I,J)*DEI
C
C  FOR SEDIMENT TRANSPORT
C
           ARCSED1(I,J,1)=ARCSED1(I,J,1)+CSED1(I,J,1)*DEI
           ARCSED2(I,J,1)=ARCSED2(I,J,1)+CSED2(I,J,1)*DEI
           ARCTAU(I,J,K)=ARCTAU(I,J,K)+TAU(I,J,K)*DEI
C
           IF (IBMSK(I,J).EQ.0) THEN
             DO 61 LL=1,LAYMAX
               ARCTHIK(I,J)=ARCTHIK(I,J)+TSED(LL,I,J)*DEI
 61          CONTINUE
           ELSE
             IF (IBMSK(I,J).EQ.1) THEN
               ARCTHIK(I,J)=ARCTHIK(I,J)+BEDTH(1,I,J)*DEI
             ELSE
               ARCTHIK(I,J)=0.0
             ENDIF
           ENDIF
C
C  FOR CHEMICAL TRANSPORT
C
           ARCHEM1(I,J,1)=ARCHEM1(I,J,1)+CHEM1(I,J,1)*DEI
           ARCHEM2(I,J,1)=ARCHEM2(I,J,1)+CHEM2(I,J,1)*DEI
           DO 62 LL=1,NCHEMLAY
             ARCPBED(LL,I,J)=ARCPBED(LL,I,J)
     +                           +CBEDCHEM(LL,I,J)*DEI
 62        CONTINUE
C
 530  CONTINUE
      ENDIF
C
      IF(INT.NE.IHIST(JHIST,2)) GO TO 560
C
C-------- WRITE CONSTANTS FIRST TIME THROUGH ---------------------------
      IF (CONSPLT) THEN
       OPEN (IUPLT,FORM='unformatted',FILE='gcmplt')
c       + CONVERT='BIG_ENDIAN')
C
       WRITE(IUPLT) IM,JM,KB
C-------- NOT NECESSARY TO SPECIFY EBCM & QBCM -------------------------
C-------- PER JAH INSTRUCTIONS - LEAVE IN FOR NOW (06-28-91) -----------
       WRITE(IUPLT) EBCM,QBCM,NCHEMLAY
       WRITE(IUPLT) DTI,GRAV,UMOL,TOR,TRACER,SEDTRAN,CHEMTRAN
       WRITE(IUPLT) NUMEBC
       IF(NUMEBC.GT.0)
     . WRITE(IUPLT) (IETA(I),JETA(I),ICON(I),JCON(I),I=1,NUMEBC)
       WRITE(IUPLT) NUMQBC
       IF(NUMQBC.GT.0)
     . WRITE(IUPLT) (IQC(I),JQC(I),I=1,NUMQBC)
       WRITE(IUPLT) H
       WRITE(IUPLT) H1
       WRITE(IUPLT) H2
       WRITE(IUPLT) ANG
       WRITE(IUPLT) DUM
       WRITE(IUPLT) DVM
       WRITE(IUPLT) FSM
       CONSPLT=.FALSE.
      ENDIF
C
C  CONVERT THICKNESS FROM g/cm**2 TO cm
C
            IF (SEDTRAN.EQ.'INCLUDE') THEN
C
              DO 650 J=1,JM
                DO 650 I=1,IM
                  IF (FSM(I,J).GT.0.0) THEN
                    IF (IBMSK(I,J).EQ.0) THEN
                      TSETOT=0.0
                      DO 640 LL=1,LAYMAX
                        TSETOT=TSETOT+TSED0(LL,I,J)
 640                  CONTINUE
C
                      ARCTHIK(I,J)=ARCTHIK(I,J)-TSETOT
C
C  CONVERT FROM g/cm**2 TO cm
C
                      ARCTHIK(I,J)=ARCTHIK(I,J)/CBED(I,J)
C
                      IF (SEDTYPE.EQ.'SAND') ARCTHIK(I,J)=0.0
                    ELSE
                      IF (IBMSK(I,J).EQ.1) THEN
C
C  IN cm
C
                         ARCTHIK(I,J)=100.*(ARCTHIK(I,J)/
     +                          ((CBED(I,J)/2.65)*
     +                          H1(I,J)*H2(I,J))-BEDTHI)
                         IF (SEDTYPE.EQ.'MUD ') ARCTHIK(I,J)=0.0
                      ELSE
                         ARCTHIK(I,J)=0.0
                      ENDIF
                    ENDIF
                  ENDIF
C
C  CONVERT FROM g/cm**3 TO mg/l
C
                  DO 660 K=1,KBM1
                    ARCSED1(I,J,K)=1000000.*ARCSED1(I,J,K)
                    ARCSED2(I,J,K)=1000000.*ARCSED2(I,J,K)
 660              CONTINUE
 650          CONTINUE
            ENDIF
C
            IF (CHEMTRAN.EQ.'INCLUDE') THEN
C
C  CONVERT FROM ug CHEM/cm**3 TO ug CHEM/l
C
              DO 1660 K=1,KBM1
                DO 1660 J=1,JM
                  DO 1660 I=1,IM
                    ARCHEM1(I,J,K)=1000.*ARCHEM1(I,J,K)               
                    ARCHEM2(I,J,K)=1000.*ARCHEM2(I,J,K)               
 1660         CONTINUE
C
           ENDIF
C
      TMIDDLE=TIME-(.5*DTI*DAYI/DEI)
      WRITE(IUPLT) TMIDDLE
      WRITE(IUPLT) ARCET
C
      IF(TOR.EQ.'BAROTROPIC') THEN
       WRITE(IUPLT) ((ARCU  (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT) ((ARCV  (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT) ((ARCUX (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT) ((ARCVX (I,J,1),I=1,IM),J=1,JM)
C
       IF (TRACER.EQ.'INCLUDE') 
     +             WRITE (IUPLT) ((ARCC1(I,J,1),I=1,IM),J=1,JM)
C
       IF (SEDTRAN.EQ.'INCLUDE') THEN
         WRITE (IUPLT) ((ARCSED1(I,J,1),I=1,IM),J=1,JM)
         WRITE (IUPLT) ((ARCSED2(I,J,1),I=1,IM),J=1,JM)
C
C  TRANSFER FROM ARCHTHIK (real*8) TO ARCET (real*4) FOR OUTPUT
C
          DO 1665 I=1,IM
            DO 1665 J=1,JM
              ARCET(I,J)=ARCTHIK(I,J)
 1665     CONTINUE
C
         WRITE (IUPLT) ((ARCET(I,J),I=1,IM),J=1,JM)
       ENDIF
C
       IF (CHEMTRAN.EQ.'INCLUDE') THEN
         WRITE (IUPLT) ((ARCHEM1(I,J,1),I=1,IM),J=1,JM)
         WRITE (IUPLT) ((ARCHEM2(I,J,1),I=1,IM),J=1,JM)
         WRITE (IUPLT) (((ARCPBED(N,I,J),I=1,IM),J=1,JM),
     +                      N=1,NCHEMLAY)
       ENDIF
C
      ELSE
       WRITE(IUPLT) Z
       WRITE(IUPLT) ZZ
       WRITE(IUPLT) DZ
       WRITE(IUPLT) ARCU
       WRITE(IUPLT) ARCV
       WRITE(IUPLT) ARCUX
       WRITE(IUPLT) ARCVX
       WRITE(IUPLT) ARCT
       WRITE(IUPLT) ARCS
       WRITE(IUPLT) ARCW
       WRITE(IUPLT) ARCKH
C
C  FOR CONSERVATIVE TRACER
C
       IF (TRACER.EQ.'INCLUDE') WRITE(IUPLT) ARCC1
C
        IF (SEDTRAN.EQ.'INCLUDE') THEN
          WRITE (IUPLT) ARCSED1
          WRITE (IUPLT) ARCSED2
C
C  TRANSFER FROM ARCHTHIK (real*8) TO ARCET (real*4) FOR OUTPUT
C
          DO 1670 I=1,IM
            DO 1670 J=1,JM
              ARCET(I,J)=ARCTHIK(I,J)
 1670     CONTINUE
C
          WRITE (IUPLT) ARCET
          WRITE (IUPLT) ARCTAU
        ENDIF
C
        IF (CHEMTRAN.EQ.'INCLUDE') THEN
          WRITE (IUPLT) ARCHEM1
          WRITE (IUPLT) ARCHEM2
          WRITE (IUPLT) (((ARCPBED(N,I,J),I=1,IM),J=1,JM),
     +                      N=1,NCHEMLAY)
        ENDIF
C
      ENDIF
C*************************************************************
C
C  PARTICLE TRACKING OUTPUT
C
c---> convert to physical domain
C
        CPI=ACOS(-1.)/180.
C
        do LL=1,nsource
          do NN=1,ngradeloop
C
C  ICHK = 0 :  NO PARTICLES IN THIS CLASS - NO OUTPUT
C  ICHK = 1 :  PARTICLES IN THIS CLASS - OUTPUT
C
            ICHK=0
C
            do MM=1,NPART
C
              XOUTP(MM)=0.0
              YOUTP(MM)=0.0
              ZOUTP(MM)=0.0
C
              iout=inout(LL,MM,NN)
              if (iout.eq.1) goto 4444
C
              ICHK=1
C
              xstart=XP(LL,MM,NN)
              ystart=YP(LL,MM,NN)
              zstart=ZP(LL,MM,NN)
C
              ic=xstart
              jc=ystart
C
              do k=1,kbm1
                if (zstart.le.z(k).and.zstart.ge.z(k+1)) goto 4445
              enddo
              WRITE (*,*) 'trouble !',zstart,int
 4445         kc=k        
C
              xc= ic+0.5
              yc= jc+0.5
              zc= zz(kc)
C
c---> interpolate D
c     find i-range for D
C
              i0h   = ic
              xlocal= xstart-xc
              xdist = abs(xlocal)
              if(xlocal.gt.0.) i1h=i0h+1
              if(xlocal.le.0.) i1h=i0h-1
C
c     use value at one element
C
              if(fsm(i1h,jc).eq.0.)i1h=i0h
C
c     find j-range for D
C
              j0h   = jc
              ylocal= ystart-yc
              ydist = abs(ylocal)
              if(ylocal.gt.0.) j1h=j0h+1
              if(ylocal.le.0.) j1h=j0h-1
C
c     use value at one element
C
              if(fsm(ic,j1h).eq.0.)j1h=j0h
C
              if(fsm(i1h,j1h).eq.0) then
                i1h=i0h
                j1h=j0h
              endif
C
c     interpolate D
C
              Dj0= D(i0h,j0h) + xdist * (D(i1h,j0h)-D(i0h,j0h))
              Dj1= D(i0h,j1h) + xdist * (D(i1h,j1h)-D(i0h,j1h))
              Dp = Dj0        + ydist * (Dj1-Dj0)
C
c     interpolate EL
C
              ELj0= EL(i0h,j0h) + xdist * (EL(i1h,j0h)-EL(i0h,j0h))
              ELj1= EL(i0h,j1h) + xdist * (EL(i1h,j1h)-EL(i0h,j1h))
              ELp = ELj0        + ydist * (ELj1-ELj0)
C      
              XOUTP(MM)=xcor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +xcor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +xcor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +xcor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
C
              YOUTP(MM)=ycor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +ycor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +ycor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +ycor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
C
              ZOUTP(MM)=ELp + zstart * Dp
C
 4444         continue
             enddo
C
C  OUTPUT PARTICLE POSITIONS FOR PARTICLES IN CLASS nn
C
               IF (ICHK.EQ.0) GOTO 900
C
               WRITE (32,930)LL,NN,NPART,TIME   
 930           FORMAT (3I8,F10.4)
C
               WRITE (32,940) (XOUTP(NP),NP=1,NPART)
               WRITE (32,940) (YOUTP(NP),NP=1,NPART)
               WRITE (32,940) (ZOUTP(NP),NP=1,NPART)
 940           FORMAT (20F10.1)
 900           CONTINUE
           enddo
         enddo
C
C*****************************************************************
C
C
      DO 570 J=1,JM
        DO 570 I=1,IM
          ARCET (I,J)=0.0
          ARCTHIK(I,J) = 0.0
          DO 570 LL=1,NCHEMLAY
            ARCPBED(LL,I,J) = 0.0
 570  CONTINUE
C
      DO 580 K=1,KB
      DO 580 J=1,JM
      DO 580 I=1,IM
      ARCU (I,J,K)=0.0
      ARCV (I,J,K)=0.0
      ARCUX(I,J,K)=0.0
      ARCVX(I,J,K)=0.0
      ARCT (I,J,K)=0.0
      ARCS (I,J,K)=0.0
      ARCW (I,J,K)=0.0
      ARCKH(I,J,K)=0.0
C
        ARCC1(I,J,K)=0.0
C
        ARCSED1(I,J,K) = 0.0
        ARCSED2(I,J,K) = 0.0
        ARCTAU(I,J,K)  = 0.0
C
        ARCHEM1(I,J,K) = 0.0
        ARCHEM2(I,J,K) = 0.0
 580  CONTINUE
 560  CONTINUE
C
 8300 CONTINUE
C
C-----------------------------------------------------------------------
C------- TIME SERIES WRITES AND ACCUMULATIONS --------------------------
      IF(TOR.EQ.'BAROTROPIC') THEN
       DO 140 N=1,EPTS
       II=INXIE(N)
       JJ=INXJE(N)
       ESAVE(N)=ESAVE(N)+EL(II,JJ)*SKILLI
140    CONTINUE
       DO 150 N=1,VPTS
         II=INXIV(N)
         JJ=INXJV(N)
         UZSAVE(N,1)=UZSAVE(N,1)+.5*(UA(II,JJ)+UA(II+1,JJ))*SKILLI
         VZSAVE(N,1)=VZSAVE(N,1)+.5*(VA(II,JJ)+VA(II,JJ+1))*SKILLI
C
C  FOR SEDIMENT TRANSPORT
C
         C1SAVE(N,1)=C1SAVE(N,1)+CSED1(II,JJ,1)*SKILLI
         C2SAVE(N,1)=C2SAVE(N,1)+CSED2(II,JJ,1)*SKILLI
         TAUSAVE(N,1)=TAUSAVE(N,1)+TAU(II,JJ,1)*SKILLI
C
         IF (IBMSK(II,JJ).EQ.0) THEN
           DO 151 LL=1,LAYMAX
             THSAVE(N)=THSAVE(N)+TSED(LL,II,JJ)*SKILLI
 151       CONTINUE
         ELSE
           IF (IBMSK(II,JJ).EQ.1) THEN
             THSAVE(N)=THSAVE(N)+BEDTH(1,II,JJ)*SKILLI
           ELSE
             THSAVE(N)=0.0
           ENDIF
         ENDIF
C
C  FOR CHEMICAL TRANSPORT
C
          P1SAVE(N,1)=P1SAVE(N,1)+CHEM1(II,JJ,1)*SKILLI
          P2SAVE(N,1)=P2SAVE(N,1)+CHEM2(II,JJ,1)*SKILLI
          DO 152 LL=1,NCHEMLAY
            PBEDSAVE(N,LL)=PBEDSAVE(N,LL)+CBEDCHEM(LL,I,J)*SKILLI
 152      CONTINUE
C
150    CONTINUE
C
C-------- COMPUTE CROSS SECTIONAL FLUXES -------------------------------
       DO 180 N=1,FPTS
       IS=ISFLX(N)
       JS=JSFLX(N)
       IF(DIRFLX(N).EQ.'IDIR') THEN
c       IE=IS+NFLXE(N)-1
        DO 181 I=IS,IS+NFLXE(N)-1
 181    CCFLUX(N,1)=CCFLUX(N,1)+VAF(I,JS)*SKILLI*0.25*
     2     (H(I,JS)+EL(I,JS)+H(I,JS-1)+EL(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
       ELSE
C        DIRFLX(N).EQ.'JDIR'
c       JE=JS+NFLXE(N)-1
        DO 182 J=JS,JS+NFLXE(N)-1
 182    CCFLUX(N,1)=CCFLUX(N,1)+UAF(IS,J)*SKILLI*0.25*
     2     (H(IS,J)+EL(IS,J)+H(IS-1,J)+EL(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
       ENDIF
 180   CONTINUE
C
C---- TOR = PROGNOSTIC OR DIAGNOSTIC OR TEMP_ONLY OR SALT_ONLY ----
      ELSE
       DO 160 N=1,EPTS
       II=INXIE(N)
       JJ=INXJE(N)
       ESAVE(N)=ESAVE(N)+ET(II,JJ)*SKILLI
160    CONTINUE
       DO 170 N=1,VPTS
       II=INXIV(N)
       JJ=INXJV(N)
       DZSAVE(N)=DZSAVE(N)+DT(II,JJ)*SKILLI
C
C  FOR SEDIMENT TRANSPORT
C
        IF (IBMSK(II,JJ).EQ.0) THEN
          DO 201 LL=1,LAYMAX
            THSAVE(N)=THSAVE(N)+TSED(LL,II,JJ)*SKILLI
 201      CONTINUE
        ELSE 
          IF (IBMSK(II,JJ).EQ.1) THEN
            THSAVE(N)=THSAVE(N)+BEDTH(1,II,JJ)*SKILLI
          ELSE
            THSAVE(N)=0.0
          ENDIF
        ENDIF
C
C  FOR CHEMICAL TRANSPORT
C
          DO 202 LL=1,NCHEMLAY
            PBEDSAVE(N,LL)=PBEDSAVE(N,LL)+CBEDCHEM(LL,II,JJ)*SKILLI
 202      CONTINUE
C
       DO 170 K=1,KB
       UZSAVE(N,K)=UZSAVE(N,K)+.5*(U(II,JJ,K)+U(II+1,JJ,K))*SKILLI
       VZSAVE(N,K)=VZSAVE(N,K)+.5*(V(II,JJ,K)+V(II,JJ+1,K))*SKILLI
       SZSAVE(N,K)=SZSAVE(N,K)+S(II,JJ,K)*SKILLI
       TZSAVE(N,K)=TZSAVE(N,K)+T(II,JJ,K)*SKILLI
C
C  FOR CONSERVATIVE TRACER
C
       C1ZSAVE(N,K)=C1ZSAVE(N,K)+CONC1(II,JJ,K)*SKILLI
C
C  FOR SEDIMENT TRANSPORT
C
         C1SAVE(N,K)=C1SAVE(N,K)+CSED1(II,JJ,K)*SKILLI
         C2SAVE(N,K)=C2SAVE(N,K)+CSED2(II,JJ,K)*SKILLI
         TAUSAVE(N,K)=TAUSAVE(N,K)+TAU(II,JJ,K)*SKILLI
C
C  FOR CHEMICAL TRANSPORT
C
         P1SAVE(N,K)=P1SAVE(N,K)+CHEM1(II,JJ,K)*SKILLI
         P2SAVE(N,K)=P2SAVE(N,K)+CHEM2(II,JJ,K)*SKILLI
C
 170   CONTINUE
C
C-------- COMPUTE CROSS SECTIONAL FLUXES -------------------------------
       DO 250 N=1,FPTS
       IS=ISFLX(N)
       JS=JSFLX(N)
       IF(DIRFLX(N).EQ.'IDIR') THEN
        DO 251 K=1,KBM1
        DO 251 I=IS,IS+NFLXE(N)-1
 251    CCFLUX(N,K)=CCFLUX(N,K)+V(I,JS,K)*DZ(K)*SKILLI*0.25*
     2     (H(I,JS)+ET(I,JS)+H(I,JS-1)+ET(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
       ELSE
C        DIRFLX(N).EQ.'JDIR'
        DO 252 K=1,KBM1
        DO 252 J=JS,JS+NFLXE(N)-1
 252    CCFLUX(N,K)=CCFLUX(N,K)+U(IS,J,K)*DZ(K)*SKILLI*0.25*
     2     (H(IS,J)+ET(IS,J)+H(IS-1,J)+ET(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
       ENDIF
 250   CONTINUE
      ENDIF
C
cjah     Initialize Variables
      VOLUME=0.0
      ESUM=0.0
      APE=0.0
      TKE=0.0
      TSUM=0.0
      SSUM=0.0
      EM=0.0
      IF(TOR.EQ.'BAROTROPIC') THEN
       DO 100 J=1,JM
       DO 100 I=1,IM
       ESUM=ESUM+EL(I,J)*ART(I,J)*FSM(I,J)/AREA*SKILLI
       APE=APE+.5*GRAV*EL(I,J)*EL(I,J)*FSM(I,J)*ART(I,J)/
     2     AREA*SKILLI
       TKE=TKE+0.125*D(I,J)*((UA(I,J)+UA(I+1,J))**2+
     2     (VA(I,J)+VA(I,J+1))**2)*FSM(I,J)*ART(I,J)/
     3     AREA*SKILLI
 100   CONTINUE
      ELSE
       DO 200 J=1,JM
       DO 200 I=1,IM
       VOLUME=VOLUME+DT(I,J)*ART(I,J)*FSM(I,J)*SKILLI
  200  ESUM=ESUM+ET(I,J)*ART(I,J)*FSM(I,J)*SKILLI
       DO 210 K=1,KBM1
       DO 210 J=1,JM
       DO 210 I=1,IM
       TSUM=TSUM+T(I,J,K)*DT(I,J)*DZ(K)*ART(I,J)*FSM(I,J)*SKILLI
  210  SSUM=SSUM+S(I,J,K)*DT(I,J)*DZ(K)*ART(I,J)*FSM(I,J)*SKILLI
       DO 220 K=1,KBM1
       DO 220 J=1,JM
       DO 220 I=1,IM
       TRHO=(RHO(I,J,K)+1.)*1000.
       VOL=DT(I,J)*ART(I,J)*DZ(K)
       EM=EM+TRHO*VOL*FSM(I,J)
  220  APE=APE+GRAV*TRHO*ZZ(K)*DT(I,J)*VOL*FSM(I,J)
      ENDIF
C
      IF(ISKILL.EQ.0.OR.MOD(INT,ISKILL).NE.0) GO TO 8250
C
C-------- WRITE CONSTANTS THE FIRST TIME THROUGH -----------------------
      IF (CONSTSR) THEN
       OPEN (IUTSR,FORM='unformatted',FILE='gcmtsr')
c       + CONVERT='BIG_ENDIAN')
       WRITE(IUTSR) TOR,TRACER,SEDTRAN,CHEMTRAN
       WRITE(IUTSR) KBM1,NCHEMLAY
       WRITE(IUTSR) EPTS
       IF(EPTS.NE.0)
     . WRITE(IUTSR) (INXIE(N),INXJE(N),N=1,EPTS)
       WRITE(IUTSR) VPTS
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (INXIV(N),INXJV(N),N=1,VPTS)
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (ANG(INXIV(N),INXJV(N)),N=1,VPTS)
       WRITE(IUTSR) FPTS
       IF(FPTS.NE.0)
     . WRITE(IUTSR) (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
       CONSTSR=.FALSE.
      ENDIF
C
          IF (SEDTRAN.EQ.'INCLUDE') THEN
C
            DO 680 N = 1, VPTS
              II = INXIV(N)
              JJ = INXJV(N)
C
              IF (IBMSK(II,JJ).EQ.0) THEN
                TSETOT=0.0
                DO 670 LL=1,LAYMAX
                  TSETOT=TSETOT+TSED0(LL,II,JJ)
 670            CONTINUE
C
                THSAVE(N)=THSAVE(N)-TSETOT
C
                IF (SEDTYPE.EQ.'SAND') THSAVE(N)=0.0
C
C  CONVERT FROM g/cm**2 TO cm
C
                THSAVE(N)=THSAVE(N)/CBED(II,JJ)
              ELSE
                THSAVE(N)=THSAVE(N)/((CBED(II,JJ)/2.65)*H1(II,JJ)
     +                                *H2(II,JJ))-BEDTHI
C
C  CONVERT FROM m TO cm
C
                THSAVE(N)=100.*THSAVE(N)
C
                IF (SEDTYPE.EQ.'MUD ') THSAVE(N)=0.0
              ENDIF
C
C  CONVERT FROM g/cm**3 TO mg/l
C
              DO 690 K=1,KBM1
                C1SAVE(N,K)=1000000.*C1SAVE(N,K)
                C2SAVE(N,K)=1000000.*C2SAVE(N,K)
 690          CONTINUE
C
              IF (CHEMTRAN.EQ.'INCLUDE') THEN
C
C  CONVERT FROM ug CHEM/cm**3 TO ug CHEM/l
C
                DO 691 K=1,KBM1
                  P1SAVE(N,K)=1000.*P1SAVE(N,K)
                  P2SAVE(N,K)=1000.*P2SAVE(N,K)
 691            CONTINUE
              ENDIF
C
 680        CONTINUE
          ENDIF
C
      TMIDDLE=TIME-(.5*DTI*DAYI/SKILLI)
      WRITE(IUTSR) TMIDDLE
C
       IF(EPTS.NE.0)
     . WRITE(IUTSR) (ESAVE(N),N=1,EPTS)
      IF(TOR.EQ.'BAROTROPIC') THEN
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (UZSAVE(N,1),VZSAVE(N,1),N=1,VPTS)
C
       IF (TRACER.EQ.'INCLUDE') THEN
         WRITE(IUTSR) (C1ZSAVE(N,1),N=1,VPTS)
       ENDIF
C
       IF (SEDTRAN.EQ.'INCLUDE') THEN
         IF (VPTS.NE.0) THEN
           WRITE (IUTSR) (C1SAVE(N,1),C2SAVE(N,1),N=1,VPTS)
           WRITE (IUTSR) (THSAVE(N),N=1,VPTS)
           WRITE (IUTSR) (TAUSAVE(N,1),N=1,VPTS)
         ENDIF
       ENDIF
C
       IF (CHEMTRAN.EQ.'INCLUDE') THEN
         IF (VPTS.NE.0) THEN
           WRITE (IUTSR) (P1SAVE(N,1),P2SAVE(N,1),N=1,VPTS)
           WRITE (IUTSR) ((PBEDSAVE(N,LL),N=1,VPTS),LL=1,NCHEMLAY)
         ENDIF
       ENDIF
C
       IF(FPTS.NE.0)
     . WRITE(IUTSR) (CCFLUX(N,1),N=1,FPTS)
       WRITE(IUTSR)  ESUM,TKE,APE
      ELSE
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (DZSAVE(N),N=1,VPTS)
       IF(VPTS.NE.0)
     . WRITE(IUTSR) ((UZSAVE(N,K),VZSAVE(N,K),SZSAVE(N,K),
     .                TZSAVE(N,K),N=1,VPTS),K=1,KBM1)
C
       IF (TRACER.EQ.'INCLUDE') THEN
         WRITE(IUTSR) ((C1ZSAVE(N,K),N=1,VPTS),K=1,KBM1)
       ENDIF
C
         IF (SEDTRAN.EQ.'INCLUDE') THEN
           WRITE (IUTSR) ((C1SAVE(N,K),C2SAVE(N,K),N=1,VPTS),
     +                                   K=1,KBM1)
           WRITE (IUTSR) (THSAVE(N),N=1,VPTS)
           WRITE (IUTSR) (TAUSAVE(N,KB),N=1,VPTS)
         ENDIF
C
         IF (CHEMTRAN.EQ.'INCLUDE') THEN
           IF (VPTS.NE.0) THEN
            WRITE (IUTSR) ((P1SAVE(N,K),P2SAVE(N,K),N=1,VPTS),
     +                              K=1,KBM1)
            WRITE (IUTSR) ((PBEDSAVE(N,LL),N=1,VPTS),LL=1,NCHEMLAY)
           ENDIF
         ENDIF
C
C
       IF(FPTS.NE.0)
     . WRITE(IUTSR) ((CCFLUX(N,K),N=1,FPTS),K=1,KBM1)
C
       VSTOR=ESUM
       EM=(EM-EMI)*SKILLI
       APEC=(APE-APEI)*SKILLI
       TSUM=TSUM/VOLUME
       SSUM=SSUM/VOLUME
C
       WRITE(IUTSR) VSTOR,EM,APEC,TSUM,SSUM
      ENDIF
C
      DO 230 N=1,EPTS
 230  ESAVE(N)=0.0
      DO 235 N=1,VPTS
 235  DZSAVE(N)=0.0
      DO 240 K=1,KB
      DO 240 N=1,VPTS
      UZSAVE(N,K)=0.0
      VZSAVE(N,K)=0.0
      SZSAVE(N,K)=0.0
      TZSAVE(N,K)=0.0
C
C  FOR CONSERVATIVE TRACER
C
      C1ZSAVE(N,K)=0.0
C
      C1SAVE(N,K) = 0.0
      C2SAVE(N,K) = 0.0
      THSAVE(N)   = 0.0
      TAUSAVE(N,K)= 0.0
      P1SAVE(N,K) = 0.0
      P2SAVE(N,K) = 0.0
C
 240  CONTINUE
C
      DO 241 LL=1,NCHEMLAY
        DO 241 N=1,VPTS
          PBEDSAVE(N,LL)=0.0
 241  CONTINUE
C
      DO 245 K=1,KB
      DO 245 N=1,FPTS
 245  CCFLUX(N,K)=0.0
      ESUM  =0.0
      TKE   =0.0
      APE   =0.0
      VSTOR =0.0
      VOLUME=0.0
      EM    =0.0
      APEC  =0.0
      TSUM  =0.0
      SSUM  =0.0
 8250 CONTINUE
C
      RETURN
      END

