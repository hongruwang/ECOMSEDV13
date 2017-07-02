      SUBROUTINE PRINT(VAR,VMASK,IM,JM,IVAR,S,IUNIT,DEV)
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
      SAVE
      CHARACTER*3 DEV
      DIMENSION VAR(IM,JM),IVAR(IM,JM),VMASK(IM,JM)
      DO 100 J=1,JM
      DO 100 I=1,IM
  100 IVAR(I,J)=IFIX(S*VAR(I,J))
      CALL MAXMIN(VAR,VMASK,IM,JM,VMAX,VMIN)
      IF(DEV.EQ.'LPR') THEN
       KN=25
      ELSE
       KN=15
      ENDIF
      IB=1
      IE=KN
      NUM=(IM-1)/KN
      IF(IM.LE.KN) NUM=1
      IF(IM.LE.KN) IE=IM
      DO 200 K=1,NUM
      IF(K.EQ.1) WRITE(IUNIT,110) S,VMAX,VMIN
  110 FORMAT(//2X,'ACTUAL VALUE = PRINTED VALUE/',1PE9.2,
     .      5X,'MAX VALUE = ',1PE9.2,2X,'MIN VALUE = ',1PE9.2)
      IF(DEV.EQ.'LPR') THEN
       WRITE(IUNIT,120) (I,I=IB,IE)
      ELSE
       WRITE(IUNIT,121) (I,I=IB,IE)
      ENDIF
  120 FORMAT(//5X,25I5)
  121 FORMAT(//5X,15I5)
      WRITE(IUNIT,130)
  130 FORMAT(/)
      DO 150 J=1,JM
      JJ=JM+1-J
      IF(DEV.EQ.'LPR') THEN
       WRITE(IUNIT,140) JJ,(IVAR(I,JJ),I=IB,IE),JJ
      ELSE
       WRITE(IUNIT,141) JJ,(IVAR(I,JJ),I=IB,IE)
      ENDIF
  140 FORMAT(' ',I2,1X,25I5,I3)
  141 FORMAT(' ',I3,1X,15I5)
  150 CONTINUE
      IB=IB+KN
      IE=IE+KN
  200 CONTINUE
      IF(IM.LE.KN) GO TO 300
      IE=IM
      IF(DEV.EQ.'LPR') THEN
       WRITE(IUNIT,120) (I,I=IB,IE)
      ELSE
       WRITE(IUNIT,121) (I,I=IB,IE)
      ENDIF
      WRITE(IUNIT,130)
      DO 250 J=1,JM
      JJ=JM+1-J
      IF(DEV.EQ.'LPR') THEN
       WRITE(IUNIT,140) JJ,(IVAR(I,JJ),I=IB,IE),JJ
      ELSE
       WRITE(IUNIT,141) JJ,(IVAR(I,JJ),I=IB,IE)
      ENDIF
  250 CONTINUE
  300 RETURN
      END
