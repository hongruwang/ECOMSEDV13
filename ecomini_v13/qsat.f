      FUNCTION QSAT(T,P)
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
C
C	T 	temperature in degrees C
C	P 	pressure in millibars
C	ES	saturation vapor pressure in millibars
C	QSAT	saturation humidity in g/kg 
C
      RD = 287.05
      RV = 461.5

      EPS = RD/RV

C*********Saturation vapor pressure**********
      ES = (1.0007+3.46E-6*P) * 6.1121 * EXP(17.502 * T/ (240.97 + T))

C*********Convert to saturation humidity*****
      QSAT = 1000. * EPS * ES / ( P + ES * (EPS - 1.0))

      RETURN
      END          
