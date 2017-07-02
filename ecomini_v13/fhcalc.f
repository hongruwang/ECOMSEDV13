      SUBROUTINE FHCALC
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
C
C  CALC. FETCH AND AVERAGE DEPTH FOR WIND WAVE MODEL
C
C
C************************************************************
C
      INCLUDE 'comdeck'
C
      REAL UTEMP(IM,JM),VTEMP(IM,JM),FTEMP(36,IM,JM)
C
C  STORE ORIGINAL VEL. VALUES (IN CASE OF HOT START)
C
      DO 5 J=1,JM
        DO 5 I=1,IM
          UTEMP(I,J)=U(I,J,1)
          VTEMP(I,J)=V(I,J,1)
 5    CONTINUE
C
C  OPEN CORNER LOCATIONS FILE
C
C  CORNER LOCATION CONVENTION:  XWCOR(I,J) = x(i-1/2, j-1/2) 
C     (LOWER LEFT-HAND CORNER)  YWCOR(I,J) = y(i-1/2, j-1/2)
C
        OPEN (UNIT=33,FILE='corner_loc',FORM='FORMATTED',STATUS='OLD')
C
        DO 300 N=1,10000
          READ (33,*,END=301)I,J,XWCOR(I,J),YWCOR(I,J)
 300    CONTINUE
 301    CLOSE(33)
C
C  CALC. FETCH AND AVERAGE DEPTH FOR EACH OF 36 DIRECTIONS (10 deg SEPARATION)
C
 310  DO 10 N=1,36
C
        WANG=2.*PI*FLOAT(N)/36.
C
C  SET VELOCITIES FOR THIS DIRECTION
C
        DO 20 J=2,JMM1
          DO 20 I=2,IMM1
            DANG=2.5*PI-ANG(I,J)
C
            PHIANG=DANG-WANG
C
            U(I,J,1)=COS(PHIANG)
            V(I,J,1)=SIN(PHIANG)
 20     CONTINUE
C
        DO 21 I=2,IM
          U(I,JM,1)=U(I,JMM1,1)
          V(I,JM,1)=V(I,JMM1,1)
 21     CONTINUE
C
        DO 22 J=2,JM
          U(IM,J,1)=U(IMM1,J,1)
          V(IM,J,1)=V(IMM1,J,1)
 22     CONTINUE
 
C
C  CALC. FETCH FOR EACH (i,j) USING THIS WIND DIRECTION
C
        DO 30 J=2,JMM1
          DO 30 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
C
C  INITIAL STARTING POINT
C
              XP1=FLOAT(I)+0.5
              YP1=FLOAT(J)+0.5
C
              IC=XP1
              JC=YP1
C
              NSUM=1
              HSUM=H(IC,JC)
C
C  IOUT = 0 :  INSIDE DOMAIN
C  IOUT = 1 :  OUTSIDE DOMAIN
C
              IOUT=0
C
 40           XSTART=XP1
              YSTART=YP1
C
              CALL SONEPART(XSTART,YSTART,XP1,YP1,IOUT)
              
C
C  DETERMINE (i,j) LOCAL GRID CELL
C
               IC=XP1
               JC=YP1
C
C  SUM DEPTHS ALONG FETCH
C
               IF (IOUT.EQ.0) THEN
                 HSUM=HSUM+H(IC,JC)
                 NSUM=NSUM+1
                 GOTO 40
               ENDIF
C
C  END OF FETCH HAS BEEN REACHED
C
C  DETERMINE LENGTH OF FETCH
C
              IC=XP1    
              JC=YP1  
C
              XC=IC+0.5 
              YC=JC+0.5 
C
              XLOCAL=XP1-XC
              YLOCAL=YP1-YC
C      
C  CALC. LOCATION OF END OF FETCH
C
              XFIN=xwcor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +xwcor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +xwcor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +xwcor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
C
              YFIN=ywcor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +ywcor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +ywcor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +ywcor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
C
C  CALC. LOCATION OF CENTER OF ELEMENT (i,j)
C
              XSTRT=0.25*(XWCOR(I,J)+XWCOR(I+1,J)+XWCOR(I,J+1)+
     +                                          XWCOR(I+1,J+1))
C
              YSTRT=0.25*(YWCOR(I,J)+YWCOR(I+1,J)+YWCOR(I,J+1)+
     +                                          YWCOR(I+1,J+1))
C
C  CALC. LENGTH OF FETCH
C
              FETCH(N,I,J)=SQRT((XFIN-XSTRT)**2.+(YFIN-YSTRT)**2.)
C
C  CALC. MEAN DEPTH ALONG FETCH
C
              HMEAN(N,I,J)=HSUM/FLOAT(NSUM)
            ENDIF
 30     CONTINUE
 10   CONTINUE
C
C  RESET ORIGINAL VEL. VALUES (IN CASE OF HOT START)
C
      DO 105 J=1,JM
        DO 105 I=1,IM
          U(I,J,1)=UTEMP(I,J)
          V(I,J,1)=VTEMP(I,J)
 105  CONTINUE
C
C  MODIFY FETCH LENGTHS TO ACCOUNT FOR FETCH WIDTH EFFECTS
C
C  USE MODIFIED VERSION OF USACOE METHOD:  AVERAGE OVER 20 DEG ARC
C  INSTEAD OF 24 DEG ARC (3 RADIALS +/- 10 DEG INSTEAD OF 9 RADIALS
C  AT +/- 12 DEG)
C
      DO 200 N=1,36
        DO 200 J=2,JMM1
          DO 200 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
              IF (N.EQ.1.OR.N.EQ.36) THEN
                IF (N.EQ.1) THEN
                  FTEMP(1,I,J)=(FETCH(36,I,J)+FETCH(1,I,J)+
     +                              FETCH(2,I,J))/3.
                ELSE
                  FTEMP(36,I,J)=(FETCH(35,I,J)+FETCH(36,I,J)+
     +                              FETCH(1,I,J))/3.
                ENDIF
              ELSE
                FTEMP(N,I,J)=(FETCH(N-1,I,J)+FETCH(N,I,J)+
     +                              FETCH(N+1,I,J))/3.
              ENDIF
            ENDIF
 200  CONTINUE
C
      DO 210 N=1,36
        DO 210 J=2,JMM1
          DO 210 I=2,IMM1
            FETCH(N,I,J)=FTEMP(N,I,J)
 210  CONTINUE
C
      RETURN
      END
C
C**************************************************************************
C
C  PARTICLE TRACKING SUBROUTINE
C
      Subroutine sonepart(xstart,ystart,xend,yend,iout)
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
* Point of Contact: Alan F. Blumberg, HydroQual, Inc.                    *
C*************************************************************************
      include 'comdeck'
C
c     absolute coordinates in tranformed grid
c     for cell i:  x=[i,i+1) left close '[', right open ')' 
c     for cell j:  y=[j,j+1) left close '[', right open ')'
C
c     find cell indieces
C
      ic= xstart
      jc= ystart
C
c     coordinates of cell center
C
      xc= ic+0.5
      yc= jc+0.5
C
      kc=1
C
c---> interpolate u  
c     find i-range for u 
C
      i0u   = ic
      i1u   = ic+1
      xlocal= xstart-xc
      xdist = abs(xlocal)
C
c     find j-range for u
C
      j0u   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1u=j0u+1
      if(ylocal.le.0.) j1u=j0u-1
C
c     free-slip on land boundary
C
      if(fsm(ic,j1u).eq.0.) j1u=j0u
C
c     find k-range for u
C
      k0u   = kc
C
c     interpolate u in k0 plane
C
      ui0k0=u(i0u,j0u,k0u)+ydist*(u(i0u,j1u,k0u)-u(i0u,j0u,k0u))
      ui1k0=u(i1u,j0u,k0u)+ydist*(u(i1u,j1u,k0u)-u(i1u,j0u,k0u))
      UP   =ui0k0*(0.5-xlocal)+ui1k0*(0.5+xlocal)
C
c---> interpolate v
C
c     find i-range for v
      i0v   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1v=i0v+1
      if(xlocal.le.0.) i1v=i0v-1
C
c     free-slip on land
C
      if(fsm(i1v,jc).eq.0.) i1v=i0v
C
c     find j-range for v
C
      j0v   = jc
      j1v   = j0v+1
      ylocal= ystart-yc
      ydist = abs(ylocal)
C
c     find k-range for v
C
      k0v   = kc
C
c     interpolate v in k0 plane
C
      vj0k0=v(i0v,j0v,k0v)+xdist*(v(i1v,j0v,k0v)-v(i0v,j0v,k0v))
      vj1k0=v(i0v,j1v,k0v)+xdist*(v(i1v,j1v,k0v)-v(i0v,j1v,k0v))
      VP   =vj0k0*(0.5-ylocal)+vj1k0*(0.5+ylocal)
C
c---> interpolate h1 and h2
C
c     find i-range for h1 and h2
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
c     find j-range for h1 and h2
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
C
      if(fsm(i1h,j1h).eq.0.)then
      i1h=i0h
      j1h=j0h
      endif
C
c     interpolate h1
C
      h1j0   = h1(i0h,j0h) + xdist * (h1(i1h,j0h)-h1(i0h,j0h))
      h1j1   = h1(i0h,j1h) + xdist * (h1(i1h,j1h)-h1(i0h,j1h))
      h1pxyz = h1j0        + ydist * (h1j1-h1j0)
C
c     interpolate h2
C
      h2j0   = h2(i0h,j0h) + xdist * (h2(i1h,j0h)-h2(i0h,j0h))
      h2j1   = h2(i0h,j1h) + xdist * (h2(i1h,j1h)-h2(i0h,j1h))
      h2pxyz = h2j0        + ydist * (h2j1-h2j0)
C
c---> advection
C
      deltax=DTI * (Up/h1pxyz)
C
      deltay=DTI*Vp/h2pxyz
C
c---> update particle location
C
      xend=xstart + deltax
      yend=ystart + deltay
C
c     indices for new cell
C
      icnew = xend
      jcnew = yend
C
      IF (XLOCAL.GT.0.49.AND.DUM(IC+1,JC).LE.0.0) IOUT=1
      IF (XLOCAL.LT.0.-49.AND.DUM(IC,JC).LE.0.0) IOUT=1
      IF (YLOCAL.GT.0.49.AND.DVM(IC,JC+1).LE.0.0) IOUT=1
      IF (YLOCAL.LT.-0.49.AND.DVM(IC,JC).LE.0.0) IOUT=1
C
      DXNEW=XEND-XC
C
      IF (DXNEW.GT.0.0) THEN
        DXNEW1=0.5-AMOD(DXNEW,0.5)
        DXNEW1=ABS(DXNEW1)
      ELSE
        DXNEW1=ABS(DXNEW)
        DXNEW1=AMOD(DXNEW1,0.5)-0.5
        DXNEW1=ABS(DXNEW1)
      ENDIF
C
      DYNEW=YEND-YC
C
      IF (DYNEW.GT.0.0) THEN
        DYNEW1=0.5-AMOD(DYNEW,0.5)
        DYNEW1=ABS(DYNEW1)
      ELSE
        DYNEW1=ABS(DYNEW)
        DYNEW1=AMOD(DYNEW1,0.5)-0.5
        DYNEW1=ABS(DYNEW1)
      ENDIF
C
      IF (DXNEW1.LE.ABS(DELTAX)) THEN
        IF (DXNEW.GT.0.0) THEN
          ICNEW=XC
          IF (DUM(ICNEW+1,JCNEW).LE.0.0) THEN
            IOUT=1
            XEND=XC+0.5
          ENDIF
        ELSE
          ICNEW=XC
          IF (DUM(ICNEW,JCNEW).LE.0.0) THEN
            IOUT=1
            XEND=XC-0.5
          ENDIF
        ENDIF
      ENDIF
C
      IF (DYNEW1.LE.ABS(DELTAY)) THEN
        IF (DYNEW.GT.0.0) THEN
          JCNEW=YC
          IF (DVM(ICNEW,JCNEW+1).LE.0.0) THEN
            IOUT=1
            YEND=YC+0.5
          ENDIF
        ELSE
          JCNEW=YC
          IF (DVM(ICNEW,JCNEW).LE.0.0) THEN
            IOUT=1
            YEND=YC-0.5
          ENDIF
        ENDIF
      ENDIF
C
      IF (FSM(ICNEW,JCNEW).LE.0.0) IOUT=1
      IF (FSM(ICNEW,JCNEW).LE.0.0) IOUT=1
C
      RETURN
      END
