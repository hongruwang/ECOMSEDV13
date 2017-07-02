      SUBROUTINE OCLAMP(CLAM)
c
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
C     Optimized open boundary conditions for ECOM3D
C     Open boundaries cross the location of the sea surface elevation
C     Version from Shulman et al. (1998)
C     Otimized Boundary Conditions and Data Assimilation with
C     Application to the M2 Tide in the Yellow Sea
C     By Igor Shulman, J.K. Lewis, A.F. Blumberg and N. Kim
C     Journal of Atmospheric and Oceanic Technology
C     *********************************************
C ************************************
C     CLAM is lambda in eqs. (4) -(5)
C     FT* are estimates of energy fluxes through the open boundary (Pt)
C     FOT* are estimates of integrals in numerator of (5)
C     DIFF3* are estimates of integrals in denominator of (5)
C ************************************
      INCLUDE 'comdeck'
      CLAM=1.0
      IF(BCTYPE.EQ.'IRANDB ') GO TO 2000
      CLAM=0.0
      FT=0.0
      DIFF3=0.0
      DO 1000 N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then
c right
          SIGN=1.0
          DEP=(D(IE,JE)+D(IC,JC))/2.
          DEL=(H2(IE,JE)+H2(IC,JC))/2.
          UVM=UA(IC,JC)
C          UVA=UA(IE,JE)+SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=UA(IE,JE)
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then
c left
          SIGN=-1.0
          DEP=(D(IE,JE)+D(IC,JC))/2.
          DEL=(H2(IC,JC)+H2(IE,JE))/2.
          UVM=UA(IC+1,JC)
C          UVA=UA(IC,JC)-SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=UA(IC,JC)
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then
c top
          SIGN=1.0
          DEP=(D(IC,JC)+D(IE,JE))/2.
          DEL=(H1(IE,JE)+H1(IC,JC))/2.
          UVM=VA(IC,JC)
C          UVA=VA(IE,JE)+SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=VA(IE,JE)
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then
c bottom
          SIGN=-1.0
          DEP=(D(IC,JC)+D(IE,JE))/2.
          DEL=(H1(IC,JC)+H1(IE,JE))/2.
          UVM=VA(IC,JC+1)
C          UVA=VA(IC,JC)-SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=VA(IC,JC)
        End If
        FT=FT+SIGN*GRAV*DEP*DEL*ELF(IE,JE)*(UVM-UVA)
        DIFF3=DIFF3+DEP*SQRT(GRAV*DEP)*DEL*UVM*UVM
1000  CONTINUE
      GAMMA=DIFF3
      IF(TIME.GT.0.05) CLAM=-FT/(DIFF3+GAMMA)
2000  CONTINUE
C     UPDATING BOUNDARY ETA FOLLOWING Eq (4)
      DO N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        DEPTH=D(IE,JE)
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
c right
          UVA=UA(IE,JE)*DUM(IE,JE)
          SIGN=1.
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
c left
          UVA=UA(IC,JC)*DUM(IC,JC)
          SIGN=-1.
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then 
c top
          UVA=VA(IE,JE)*DVM(IE,JE)
          SIGN=1.
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then 
c bottom
          UVA=VA(IC,JC)*DVM(IC,JC)
          SIGN=-1.
        End If
        ELF(IE,JE)=ELF(IE,JE)+SIGN*UVA*SQRT(DEPTH/GRAV)*CLAM
      ENDDO
      RETURN
      END
