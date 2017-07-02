      subroutine PARTRAK(idum,IRELEND)
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
      include 'comdeck'
C
      common/uvwpmax/upmax,vpmax,wpmax,deltaxmax,deltaymax,deltazmax
     .              ,udpmax,vdpmax,wdpmax
C
      if (int.lt.IRELST) RETURN
      intip=int-IRELST 
C
c---> calculate how many groups of particles in the system
C
      ngradeloop=intip/NFREQ + 1
      if (ngradeloop.ge.NCONV) ngradeloop=NCONV
C
c---> advance all groups to a higher grade before new group is introduced
C
chli      if (mod(intip,NFREQ).eq.0.and.ngradeloop.gt.1) then
       if (mod(intip,NFREQ).ne.0.or.int.gt.irelend)go to 31
       if(ngradeloop.eq.1) go to 29      ! nkim 061098
        do 10 LL=1,nsource
          do 10 NN=ngradeloop,2,-1
            do 10 MM=1,NPART
              XP(LL,MM,NN)=   XP(LL,MM,NN-1)
              YP(LL,MM,NN)=   YP(LL,MM,NN-1)
              ZP(LL,MM,NN)=   ZP(LL,MM,NN-1)
              inout(LL,MM,NN)=inout(LL,MM,NN-1)
10      CONTINUE
c      endif                       ! nkim 061098
C
  29   continue
c---> introduce a new group every NFREQ timestep
C
chli      if (mod(intip,NFREQ).eq.0) then
C
        do 30 LL=1,nsource
          do 30 MM=1,NPART
            XP(LL,MM,1)= isource(LL)+0.5
            YP(LL,MM,1)= jsource(LL)+0.5
            ZP(LL,MM,1)= ZZ(ksource(LL))
            inout(LL,MM,1)  = 0
 30     CONTINUE
chli      endif
C
c---> advance all particles in all possible groups at the current time
C
 31   continue
      do 40 LL=1,nsource
        do 40 NN=1,ngradeloop
          do 40 MM=1,NPART
C
            if(inout(LL,MM,NN).eq.1) goto 40
C
            xstart=xp(LL,MM,NN)
            ystart=yp(LL,MM,NN)
            zstart=zp(LL,MM,NN)
            iout=inout(LL,MM,NN)
C
            call onepart(xstart,ystart,zstart,idum,dti,xend,yend,
     +                            zend,iout)
C
            xp(LL,MM,NN)=xend
            yp(LL,MM,NN)=yend
            zp(LL,MM,NN)=zend
            inout(LL,MM,NN)=iout
 40   CONTINUE
C
      return
      end
C
C**************************************************************************
C
      function gasdev(idum)
      data iset/0/
      if( iset.eq.0 )then
 1    continue
       v1=2.0*ran2(idum)-1.
       v2=2.0*ran2(idum)-1.
       r=v1*v1+v2*v2
       if(r.ge.1.0)goto 1
       fac=sqrt(-2.0*log(r)/r)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
      else
       gasdev=gset
       iset=0
      endif
C
      return
      end
C
C**************************************************************************
C
      function ran2(idum)
      parameter(m=714025, ia=1366, ic=150889, rm=1.0/m)
      dimension ir(97)
      data iff/0/
      if(idum.lt.0.or.iff.eq.0)then
         iff=1
         idum=mod(ia-idum,m)
         do 11 j=1,97
         idum=mod(ia*idum+ic,m)
         ir(j)=idum
 11      continue
         idum=mod(ia*idum+ic,m)
         iy=idum
      endif
         j=1+(97*iy)/m
         if(j.gt.97.or.j.lt.1)then
         print *, 'stop in ran2 ',j
         stop
         endif
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end
