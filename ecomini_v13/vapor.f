      SUBROUTINE VAPOR(EA,ES,RH,SST,BP)
C
CMNT  RH IN PERCENT, SST IN DEG C, BP IN MB
CMNT  COMPUTE VAPOR PRESSURE OF WATER FROM RELATIVE
CMNT  HUMIDITY
CMNT
CMNT  R.WELLER, 6/18/87
CMNT
CMNT  COMPUTE SATURATION VAPOR PRESSURE OF PURE WATER
CMNT  USING TABATA FORMULA
C**********************************************************
C	Modified November 6, 1990	M. Park Samelson
C	
C	Previously used Tabata formula (1973).  This
C	version uses Buck formula (1981).  Temperature
C	is in degrees C instead of degrees K.
C
C**********************************************************

C*******Tabata formula**************
C**** TK=273.16
C**** SSTK=SST+TK
C**** A=1000./SSTK
C**** E=8.42926609-(1.82717483*A)-(.071208271*A*A)
C**** ES=10.**E

C********Buck formula***************
      ES = (1.0007+3.46E-6*BP)*6.1121*EXP(17.502 * SST/ (240.97 + SST))

CMNT  CORRECT FOR SALT WATER AT 35 PPT
      EW=.9815*ES
CMNT  COMPUTE SATURATION MIXING RATIO
      B=EW/BP
      RW=(B*.62197)/(1.-B)
CMNT  COMPUTE MIXING RATIO
      R=(RH/100.)*RW
CMNT  COMPUTE VAPOR PRESSURE
      EA=(R/(.62197+R))*BP
      RETURN
      END
