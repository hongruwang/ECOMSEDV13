        SUBROUTINE NORTH(WX,WY,dir)
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
c       This program was updated by Ray Sheppard. After a couple of bugs
c       are removed, incorporated in ECOMSED
c 
C       WX=wind in x-dir, WY=wind in y-dir
C       DIRECTION  GIVES ANGLE WITH X-
C       AXIS (NORTH) POSITIVE- CLOCKWISE, NEGATIVE- ANTI CLOCKWISE
C       this gives the direction of wind blow from (meteorological conv)
C       
        data RAD/57.29577/
        IF (WX .EQ. 0.) THEN
            IF (WY .GE. 0.) THEN
               dir = 180. 
            ELSE
               dir = 0.0
            ENDIF
        ELSE
            GRAD = WY/WX
            dir = ATAN(GRAD) * RAD
            IF (WY .GE. 0.) THEN
                IF (WX .GT. 0.) then
                   dir = 270. - dir
                ELSE
                   dir = abs(dir) + 90.
                ENDIF
            ELSE
                IF (WX .LT. 0.) THEN
                   dir=90-abs(dir)
                ELSE
                   dir=abs(dir)+270
                ENDIF
            ENDIF
        ENDIF
        RETURN
        END
