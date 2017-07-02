      SUBROUTINE GETDATE(ID,IM,IY,IH,IMIN,THOURS)
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
C     THIS ROUTINE GIVES THE DAY, MON, YEAR, HOUR AND MINUTES 
C     CORRESPONDING TO THE GIVEN HOURS
C
C  OSSM01.FOR
C
      COMMON /DEVCES/ LPRT,LCON,LTPE,LDSK,LWCON
      DIMENSION MON(12,2)
      DIMENSION LYEAR(4)
      DATA(LYEAR(I),I=1,4)/8784,8760,8760,8760/
      DATA (MON(I,1),I=1,12)/0,31,59,90,120,151,181,212,243,273,
     C304,334/
      DATA (MON(I,2),I=1,12)/0,31,60,91,121,152,182,213,244,274,
     C305,335/
      THOURS=THOURS+.00417
      XHOURS=MOD(THOURS,24.0)
      IH=INT(XHOURS)
      IMIN=INT((XHOURS-FLOAT(IH))*60.)
      IHOURS=INT(THOURS-XHOURS+.001)
      IY=0
      YDSUM=0
 5    CONTINUE
      DO 10 I=1,4
      YDSUM=YDSUM+LYEAR(I)
      IY=IY+1
      IF(YDSUM.GT.IHOURS) GOTO 75
      IF(YDSUM.EQ.IHOURS) GOTO 85
 10   CONTINUE
      GOTO 5
 75   YDSUM=YDSUM-LYEAR(I)
      IY=IY-1
  85   IF(I.EQ.1) GOTO 200
      GOTO 100
 100  CONTINUE
      K=1
 110  IDAYS=(IHOURS-YDSUM)/24+1
      DO 120 JFLAG=1,12
      IF(MON(JFLAG,K).GE.IDAYS) GOTO 130
 120  CONTINUE
      IM=12
      ID=IDAYS-MON(12,K)
      GOTO 999
 130  ID=IDAYS-MON(JFLAG-1,K)
      IM=JFLAG-1
      GOTO 999
 200  K=2
      GOTO 110
 999  CONTINUE
      THOURS=THOURS-.00417
c change year from 2 digitis to 4 digits
      IY=IY+1900
      RETURN
      END
      Subroutine IJDAY(D,M,Y,J)
C     ************************************************************
C**   WHERE...
C        D.... DAY (1-31) COMPONENT OF GREGORIAN DATE
C        M.... MONTH (1-12) COMPONENT
C        Y.... YEAR (E.G., 1979) COMPONENT
C        J.... CORRESPONDING JULIAN DAY NUMBER
C     ************************************************************
      Integer*4 C, J, MO, YR
c     INTEGER*2 D,M,Y
      Integer*4 D, M, Y
C
      YR = Y
      If (M.GT.2) Then
        MO = M - 3
      Else
        MO = M + 9
        YR = YR - 1
      End If
C
      C = YR / 100
      YR = YR - C * 100
      J = (146097*C) / 4 + (1461*YR) / 4 + (153*MO+2) / 5 + D + 1721119
c
c     convert to days since May 23, 1968 
c       (which in True Julian day is 2440000)
      J = J - 2440000
      Return
      End
