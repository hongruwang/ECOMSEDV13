      SUBROUTINE PCLAMP
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
      INCLUDE 'comdeck'
C
        DO 100 N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
C
C  SOUTH BOUNDARY (RADIATION & FORCING) --------------------
C
      IF (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then
      CPH=SQRT(9.8*H(IE,JE))*2.*DTE/(H2(IE,JE)+H2(IE,JE+1))
        ELF(IE,JE)=(ELB(IE,JE)+CPH*(2.*EL(IE,JE+1)-ELB(IE,JE))+2.*DTE*
     .        (ELF(IE,JE)-ELB(IE,JE))*TLAG)/(1.+CPH)
C
C  --- EAST BOUNDARY (RADIATION & FORCING) 
C
      ELSE IF (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then
      CPH=SQRT(9.8*H(IE,JE))*2.*DTE/(H1(IE,JE)+H1(IE-1,JE))
         ELF(IE,JE)=(ELB(IE,JE)+CPH*(2.*EL(IE-1,JE)-ELB(IE,JE))+2.*DTE*
     .        (ELF(IE,JE)-ELB(IE,JE))*TLAG)/(1.+CPH)
C
C  --- WEST BOUNDARY (RADIATION & FORCING)
C
      ELSE IF (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then
      CPH=SQRT(9.8*H(IE,JE))*2.*DTE/(H1(IE,JE)+H1(IE+1,JE))
         ELF(IE,JE)=(ELB(IE,JE)+CPH*(2.*EL(IE+1,JE)-ELB(IE,JE))+2.*DTE*
     .        (ELF(IE,JE)-ELB(IE,JE))*TLAG)/(1.+CPH)

C
C  --- NORTH BOUNDARY (RADIATION & FORCING)
C
      ELSE IF (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then
      CPH=SQRT(9.8*H(IE,JE))*2.*DTE/(H2(IE,JE)+H2(IE,JE-1))
        ELF(IE,JE)=(ELB(IE,JE)+CPH*(2.*EL(IE,JE-1)-ELB(IE,JE))+2.*DTE*
     .        (ELF(IE,JE)-ELB(IE,JE))*TLAG)/(1.+CPH)

      ENDIF
 100  CONTINUE
C
      RETURN
      END
