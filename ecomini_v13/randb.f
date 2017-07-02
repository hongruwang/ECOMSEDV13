      SUBROUTINE RANDB
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
      DO 1000 N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then
c right
          DEPTH=(D(IE,JE)+D(IC,JC))/2.
          UAF(IE,JE)=SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then
c left
          DEPTH=(D(IE,JE)+D(IC,JC))/2.
          UAF(IC,JC)=-SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then
c top
c          DEPTH=(D(IC,JC)+D(IE,JE))/2.
          DEPTH=(h(IC,JC)+h(IE,JE))/2.
          VAF(IE,JE)=SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then
c bottom
          DEPTH=(D(IC,JC)+D(IE,JE))/2.
          VAF(IC,JC)=-SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        End If
1000  CONTINUE
      RETURN
      END
