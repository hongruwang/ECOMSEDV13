      SUBROUTINE ANTIDIF_SED(XMFLUX,YMFLUX,ZZFLUX,FB,FF,DTI2,
     .TSDIS,TSBDRY)
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
* Point of Contact: ecomsed-support@hydroqual.com                        *
C*************************************************************************
      INCLUDE 'comdeck'
      SAVE
C
      DIMENSION FF(IM,JM,KB),FB(IM,JM,KB)
      DIMENSION XMFLUX(IM,JM,KB),YMFLUX(IM,JM,KB),ZZFLUX(IM,JM,KB)
      DIMENSION UP(IM,JM,KB),DOWN(IM,JM,KB)
      DIMENSION FLUXIN(IM,JM,KB),FLUXOUT(IM,JM,KB)
      DIMENSION FMAX(IM,JM,KB),FMIN(IM,JM,KB)
      DIMENSION TSDIS(QBCM),TSBDRY(EBCM,KBM1)
      EQUIVALENCE (UP,A),(DOWN,C),(FMAX,PROD)
      EQUIVALENCE (FLUXIN,VH),(FLUXOUT,VHP)
      BIG=1.E20
C
C-------- TEMPERATURE & SALINITY BOUNDARY CONDITIONS -------------------
      DO 235 N=1,NUMQBCSE
      ID=ISEQD(N)
      JD=JSEQD(N)
      IC=ISEQC(N)
      JC=JSEQC(N)
      DO 235 K=1,KBM1
      FF(IC,JC,K)=TSDIS(N)
 235  CONTINUE
      DO 236 N=1,NUMEBCSE
      IE=ISEED(N)
      JE=JSEED(N)
      DO 170 K=1,KBM1
      FF(IE,JE,K)=TSBDRY(N,K)
  170 CONTINUE
  236 CONTINUE
      DO 240 K=1,KBM1
      DO 240 I=1,IM
      DO 240 J=1,JM
      FF(I,J,K)=FF(I,J,K)*FSM(I,J)
 240  CONTINUE
      DO 101 J=1,JM
      DO 101 I=1,IM
  101 FF(I,J,KB)=FF(I,J,KBM1)
C-----------------------------------------------------------------------
C    RECALCULATE MASS FLUXES USING ANTIDIFFUSION VELOCITY 
C-----------------------------------------------------------------------
      IF(SCHEME .EQ. 'SMOLAR_2  ') THEN
      DO 142 K=1,KBM1
      DO 142 J=2,JMM1
      DO 142 I=2,IM
      IF(FF(I,J,K).LT.1.E-9.OR.FF(I-1,J,K).LT.1.E-9) THEN
         XMFLUX(I,J,K)=0.0 
      ELSE
         AA=(FF(I,J,K)-FF(I-1,J,K))/(FF(I-1,J,K)+FF(I,J,K)+1.0E-15)
         UDELTAX=ABS(XMFLUX(I,J,K))*AA
          U2DELTAT=(    DTI2*XMFLUX(I,J,K)*XMFLUX(I,J,K)/
     #             (ARU(I,J)*(DT(I-1,J)+DT(I,J))) )*AA
          XMFLUX(I,J,K)=UDELTAX - U2DELTAT
          ANTI1=ABS(UDELTAX) 
          ANTI2=ABS(U2DELTAT)
          IF(ANTI1 .LT. ANTI2)XMFLUX(I,J,K)=0.0
       END IF
  142 CONTINUE
      DO 143 K=1,KBM1
      DO 143 J=2,JM
      DO 143 I=2,IMM1
      IF(FF(I,J,K).LT.1.E-9.OR.FF(I,J-1,K).LT.1.E-9) THEN
         YMFLUX(I,J,K)=0.0  
      ELSE
         AA=(FF(I,J,K)-FF(I,J-1,K))/(FF(I,J-1,K)+FF(I,J,K)+1.0E-15)
         VDELTAY=ABS(YMFLUX(I,J,K))*AA
          V2DELTAT=(    DTI2*YMFLUX(I,J,K)*YMFLUX(I,J,K)/
     #             (ARV(I,J)*(DT(I,J-1)+DT(I,J))) )*AA
          YMFLUX(I,J,K)=VDELTAY - V2DELTAT
          ANTI1=ABS(VDELTAY)
          ANTI2=ABS(V2DELTAT)
          IF(ANTI1 .LT. ANTI2)YMFLUX(I,J,K)=0.0
         END IF
  143 CONTINUE
      DO 144 K=2,KBM1
      DO 144 J=2,JMM1
      DO 144 I=2,IMM1
      IF(FF(I,J,K).LT.1.E-9.OR.FF(I,J,K-1).LT.1.E-9) THEN
         ZZFLUX(I,J,K)=0.0 
      ELSE
         AA=(FF(I,J,K-1)-FF(I,J,K))/(FF(I,J,K)+FF(I,J,K-1)+1.0E-15)
         WDELTAZ=ABS(ZZFLUX(I,J,K))*AA
          W2DELTAT=( DTI2*ZZFLUX(I,J,K)*ZZFLUX(I,J,K)/
     #             (DZZ(K-1)*DT(I,J)) )*AA
          ZZFLUX(I,J,K)= WDELTAZ - W2DELTAT
          ANTI1=ABS(WDELTAZ)
          ANTI2=ABS(W2DELTAT)
          IF(ANTI1 .LT. ANTI2)ZZFLUX(I,J,K)=0.0
         END IF
  144 CONTINUE
C
      ELSE
C

      DO 242 K=1,KBM1
      DO 242 J=2,JMM1
      DO 242 I=2,IM
      IF(FF(I,J,K).LT.1.E-9.OR.FF(I-1,J,K).LT.1.E-9) THEN
         XMFLUX(I,J,K)=0.0 
      ELSE
         AA=(FF(I,J,K)-FF(I-1,J,K))/(FF(I-1,J,K)+FF(I,J,K)+1.0E-15)
         UDELTAX=ABS(XMFLUX(I,J,K))*AA
         XMFLUX(I,J,K)=UDELTAX/(1.-ABS(AA)+1.0E-15) 
      END IF
  242 CONTINUE

      DO 243 K=1,KBM1
      DO 243 J=2,JM
      DO 243 I=2,IMM1
      IF(FF(I,J,K).LT.1.E-9.OR.FF(I,J-1,K).LT.1.E-9) THEN
         YMFLUX(I,J,K)=0.0  
      ELSE
         AA=(FF(I,J,K)-FF(I,J-1,K))/(FF(I,J-1,K)+FF(I,J,K)+1.0E-15)
         VDELTAY=ABS(YMFLUX(I,J,K))*AA
         YMFLUX(I,J,K)=VDELTAY/(1.-ABS(AA)+1.0E-15)  
      END IF
  243 CONTINUE
      DO 244 K=2,KBM1
      DO 244 J=2,JMM1
      DO 244 I=2,IMM1
      IF(FF(I,J,K).LT.1.E-9.OR.FF(I,J,K-1).LT.1.E-9) THEN
         ZZFLUX(I,J,K)=0.0 
      ELSE
         AA=(FF(I,J,K-1)-FF(I,J,K))/(FF(I,J,K)+FF(I,J,K-1)+1.0E-15)
         WDELTAZ=ABS(ZZFLUX(I,J,K))*AA
         ZZFLUX(I,J,K)=WDELTAZ/(1.-ABS(AA)+1.0E-15)
      END IF
  244 CONTINUE

      END IF
C
C------------- ADJUST FOR RIVER/WALL FLUXES -------------
      DO 551 N=1,NUMQBCSE
      ID=ISEQD(N)
      JD=JSEQD(N)
      IC=ISEQC(N)
      JC=JSEQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 661 K=1,KBM1
  661         XMFLUX(IC,JC,K)=0.0
            ELSE
              DO 771 K=1,KBM1
  771         XMFLUX(ID,JD,K)=0.0
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 881 K=1,KBM1
  881         YMFLUX(IC,JC,K)=0.0
            ELSE
              DO 991 K=1,KBM1
  991         YMFLUX(ID,JD,K)=0.0
            ENDIF
      ENDIF
  551 CONTINUE
C-----------------------------------------------------------------------
C    LIMITING THE ANTIDIFFUSION VELOCITY 
C-----------------------------------------------------------------------
      IF(SCHEME .NE. 'SMOLAR_R  ') RETURN
      DO 716 K=1,KB
      DO 716 J=1,JM
      DO 716 I=1,IM
        UP(I,J,K)=0.0
        DOWN(I,J,K)=0.0
        FLUXIN(I,J,K)=0.0
        FLUXOUT(I,J,K)=0.0
        FMAX(I,J,K)=0.0
        FMIN(I,J,K)=0.0
  716 CONTINUE

      DO 719 K=1,KBM1
      DO 719 J=2,JMM1
      DO 719 I=2,IM
      DUMMY=2.*XMFLUX(I,J,K)/(ARU(I,J)*(DT(I-1,J)+DT(I,J)))
      UP(I,J,K)=AMAX1(0.,DUMMY)
  719 DOWN(I,J,K)=AMIN1(0.,DUMMY)
      DO 720 K=1,KBM1
      DO 720 J=2,JMM1
      DO 720 I=2,IMM1
      FLUXIN(I,J,K)=UP(I,J,K)*FF(I-1,J,K)-DOWN(I+1,J,K)*FF(I+1,J,K)
  720 FLUXOUT(I,J,K)=UP(I+1,J,K)*FF(I,J,K)-DOWN(I,J,K)*FF(I,J,K)
C
      DO 730 K=1,KBM1 
      DO 730 J=2,JM
      DO 730 I=2,IMM1
      DUMMY=2.*YMFLUX(I,J,K)/(ARV(I,J)*(DT(I,J-1)+DT(I,J)))
      UP(I,J,K)=AMAX1(0.,DUMMY)
  730 DOWN(I,J,K)=AMIN1(0.,DUMMY)
      DO 820 K=1,KBM1
      DO 820 J=2,JMM1
      DO 820 I=2,IMM1
      FLUXIN(I,J,K)=FLUXIN(I,J,K)+
     .              UP(I,J,K)*FF(I,J-1,K)-DOWN(I,J+1,K)*FF(I,J+1,K)
  820 FLUXOUT(I,J,K)=FLUXOUT(I,J,K)+
     .               UP(I,J+1,K)*FF(I,J,K)-DOWN(I,J,K)*FF(I,J,K)
C
      DO 721 K=2,KBM1
      DO 721 J=2,JMM1
      DO 721 I=2,IMM1
      DUMMY=ZZFLUX(I,J,K)/(DZZ(K-1)*DT(I,J))
      UP(I,J,K)=AMAX1(0.,DUMMY)
  721 DOWN(I,J,K)=AMIN1(0.,DUMMY)
      DO 722 J=2,JMM1
      DO 722 I=2,IMM1
      FLUXIN(I,J,1)=FLUXIN(I,J,1)+UP(I,J,2)*FF(I,J,2) 
      FLUXIN(I,J,KBM1)=FLUXIN(I,J,KBM1)-DOWN(I,J,KBM1)*FF(I,J,KBM2)
      FLUXOUT(I,J,1)=FLUXOUT(I,J,1)-DOWN(I,J,2)*FF(I,J,1)
  722 FLUXOUT(I,J,KBM1)=FLUXOUT(I,J,KBM1)+UP(I,J,KBM1)*FF(I,J,KBM1) 
      DO 718 K=2,KBM2
      DO 718 J=2,JMM1
      DO 718 I=2,IMM1
      FLUXIN(I,J,K)=FLUXIN(I,J,K)+
     .              UP(I,J,K+1)*FF(I,J,K+1)-DOWN(I,J,K)*FF(I,J,K-1)
  718 FLUXOUT(I,J,K)=FLUXOUT(I,J,K)+
     .               UP(I,J,K)*FF(I,J,K)-DOWN(I,J,K+1)*FF(I,J,K)
C
C
      DO 734 J=2,JMM1
      DO 734 I=2,IMM1
      FMAX(I,J,1)=AMAX1(FB(I-1,J,1),FB(I,J,1),FB(I+1,J,1),
     .                     FF(I-1,J,1),FF(I,J,1),FF(I+1,J,1),
     .      FB(I,J-1,1),FB(I,J+1,1), FF(I,J-1,1),FF(I,J+1,1),
     .                                FB(I,J,2), FF(I,J,2  ))
  734 CONTINUE
C
      DO 733 K=2,KBM2 
      DO 733 J=2,JMM1
      DO 733 I=2,IMM1
      FMAX(I,J,K)=AMAX1(FB(I-1,J,K),FB(I,J,K),FB(I+1,J,K),
     .                     FF(I-1,J,K),FF(I,J,K),FF(I+1,J,K),
     .      FB(I,J-1,K),FB(I,J+1,K), FF(I,J-1,K),FF(I,J+1,K),
     .      FB(I,J,K-1),FB(I,J,K+1), FF(I,J,K-1),FF(I,J,K+1))
  733 CONTINUE 
C
      DO 743 J=2,JMM1
      DO 743 I=2,IMM1
      FMAX(I,J,KBM1)=AMAX1(FB(I-1,J,KBM1),FB(I,J,KBM1),FB(I+1,J,KBM1),
     .                     FF(I-1,J,KBM1),FF(I,J,KBM1),FF(I+1,J,KBM1),
     .   FB(I,J-1,KBM1),FB(I,J+1,KBM1), FF(I,J-1,KBM1),FF(I,J+1,KBM1),
     .                     FB(I,J,KBM1-1), FF(I,J,KBM1-1))
  743 CONTINUE 
C
C
      DO 220 J=1,JM
      DO 220 I=1,IM
      IF(FSM(I,J).EQ.1.) GO TO 220
      DO 221 K=1,KB 
      FB(I,J,K)=BIG
  221 FF(I,J,K)=BIG
  220 CONTINUE
C
      DO 728 J=2,JMM1
      DO 728 I=2,IMM1
      FMIN(I,J,1)=AMIN1(FB(I-1,J,1),FB(I,J,1),FB(I+1,J,1),
     .                   FF(I-1,J,1),FF(I,J,1),FF(I+1,J,1),
     .    FB(I,J-1,1),FB(I,J+1,1), FF(I,J-1,1),FF(I,J+1,1),
     .                                FB(I,J,2), FF(I,J,2))
  728 CONTINUE  
C
      DO 735 K=2,KBM2
      DO 735 J=2,JMM1
      DO 735 I=2,IMM1
      FMIN(I,J,K)=AMIN1(FB(I-1,J,K),FB(I,J,K),FB(I+1,J,K),
     .                   FF(I-1,J,K),FF(I,J,K),FF(I+1,J,K),
     .    FB(I,J-1,K),FB(I,J+1,K), FF(I,J-1,K),FF(I,J+1,K),
     .    FB(I,J,K-1),FB(I,J,K+1), FF(I,J,K-1),FF(I,J,K+1))
  735 CONTINUE  
C
      DO 745 J=2,JMM1
      DO 745 I=2,IMM1
      FMIN(I,J,KBM1)=AMIN1(FB(I-1,J,KBM1),FB(I,J,KBM1),FB(I+1,J,KBM1),
     .                     FF(I-1,J,KBM1),FF(I,J,KBM1),FF(I+1,J,KBM1),
     .   FB(I,J-1,KBM1),FB(I,J+1,KBM1), FF(I,J-1,KBM1),FF(I,J+1,KBM1),
     .                     FB(I,J,KBM1-1), FF(I,J,KBM1-1))
  745 CONTINUE 
C
      DO 222 K=1,KB 
      DO 222 J=1,JM
      DO 222 I=1,IM
      FB(I,J,K)=FB(I,J,K)*FSM(I,J)
  222 FF(I,J,K)=FF(I,J,K)*FSM(I,J)
C
C
      DO 736 K=1,KBM1 
      DO 736 J=1,JM
      DO 736 I=1,IM
      UP(I,J,K)=(FMAX(I,J,K)-FF(I,J,K))/
     .            (DTI2*FLUXIN(I,J,K)+1.0E-15)
      DOWN(I,J,K)=(FF(I,J,K)-FMIN(I,J,K))/
     .              (DTI2*FLUXOUT(I,J,K)+1.0E-15)
  736 CONTINUE

 
      DO 750 N=1,NUMQBCSE
      IC=ISEQC(N) 
      JC=JSEQC(N)
      DO 750 K=1,KBM1 
      UP(IC,JC,K)=0.0 
      DOWN(IC,JC,K)=0.0
  750 CONTINUE  
C
      DO 737 K=1,KBM1
      DO 737 J=2,JM
      DO 737 I=2,IM
      UPLUS=AMAX1(0.,XMFLUX(I,J,K))
      VPLUS=AMAX1(0.,YMFLUX(I,J,K))
      UMINUS=AMIN1(0.,XMFLUX(I,J,K))
      VMINUS=AMIN1(0.,YMFLUX(I,J,K))
      XMFLUX(I,J,K)=AMIN1(1.,DOWN(I-1,J,K),UP(I,J,K))*UPLUS +
     .              AMIN1(1.,UP(I-1,J,K),DOWN(I,J,K))*UMINUS
      YMFLUX(I,J,K)=AMIN1(1.,DOWN(I,J-1,K),UP(I,J,K))*VPLUS +
     .              AMIN1(1.,UP(I,J-1,K),DOWN(I,J,K))*VMINUS
  737 CONTINUE
C
      DO 739 K=2,KBM1 
      DO 739 J=2,JMM1
      DO 739 I=2,IMM1
      WPLUS=AMAX1(0.,ZZFLUX(I,J,K))
      WMINUS=AMIN1(0.,ZZFLUX(I,J,K))
      ZZFLUX(I,J,K)=AMIN1(1.,DOWN(I,J,K),UP(I,J,K-1))*WPLUS +
     .              AMIN1(1.,UP(I,J,K),DOWN(I,J,K-1))*WMINUS
  739 CONTINUE
 
      RETURN
      END
