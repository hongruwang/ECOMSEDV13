      Subroutine JDAY(D,M,Y,J)
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
