      SUBROUTINE DENS
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
      INCLUDE 'comdeck'
      SAVE
C
C-----------------------------------------------------------------------
C         THIS FUNCTION COMPUTES DENSITY-1.0
C-----------------------------------------------------------------------
C
      REAL*8 RHOF(IM,JM,KB),TF(IM,JM,KB),SF(IM,JM,KB)
      EQUIVALENCE (A,RHOF),(VH,TF),(UF,SF)
C
      DO 1 K=1,KBM1
      DO 1 J=1,JM
      DO 1 I=1,IM
      TF(I,J,K)=T(I,J,K)
      SF(I,J,K)=S(I,J,K)
c     SF(I,J,K)=20.*fsm(i,j)
      RHOF(I,J,K)=
     . SF(I,J,K)*SF(I,J,K)*SF(I,J,K)*6.76786136E-6-SF(I,J,K)*SF(I,J,K)*
     . 4.8249614E-4+SF(I,J,K)*8.14876577E-1-0.22584586E0
C
      RHOF(I,J,K)=RHOF(I,J,K)*
     . (TF(I,J,K)*TF(I,J,K)*TF(I,J,K)*1.667E-8-TF(I,J,K)*TF(I,J,K)
     . *8.164E-7+TF(I,J,K)*1.803E-5)
C
      RHOF(I,J,K)=RHOF(I,J,K)+1.- 
     . TF(I,J,K)*TF(I,J,K)*TF(I,J,K)*1.0843E-6+TF(I,J,K)*TF(I,J,K)
     . *9.8185E-5-TF(I,J,K)*4.786E-3
C
      RHOF(I,J,K)=RHOF(I,J,K)*
     . (SF(I,J,K)*SF(I,J,K)*SF(I,J,K)*6.76786136E-6-SF(I,J,K)*SF(I,J,K)*
     . 4.8249614E-4+SF(I,J,K)*8.14876577E-1+3.895414E-2)
C
      RHOF(I,J,K)=RHOF(I,J,K)-
     . (TF(I,J,K)-3.98)**2 * (TF(I,J,K)+283.)/(503.57*(TF(I,J,K)+67.26))
    1 CONTINUE
C
      DO 2 K=1,KBM1
      DO 2 J=1,JM
      DO 2 I=1,IM
      RHO(I,J,K)= RHOF(I,J,K) * 1.E-3*FSM(I,J)
    2 CONTINUE
      RETURN
      END
