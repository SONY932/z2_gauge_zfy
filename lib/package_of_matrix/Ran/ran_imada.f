C======================================================
      SUBROUTINE RAN_IMADA(ANRAN,NGMA1,NGMAX,IQR1,IQR2)
C**** MAKE RANDOM SEQUENCE: LENGTH = NGMAX ******
      IMPLICIT REAL*8(A-G,O-Z)
      REAL*8 ANRAN(NGMA1)
      INTEGER IQR1(250),IQR2(250)

      XNORM = 1.0/2.0D0**31
      DO 10 JSF=0,NGMAX/250-1
      JSR=JSF*250
      CALL RANDN2(IQR1,IQR2)
      DO 20 JK=1,250
      JJSF=JK+JSR
      ANRAN(JJSF)=IQR1(JK)*XNORM ! 0.4656613D-9
   20 CONTINUE
   10 CONTINUE

      RETURN
      END
C======================================================
******* MAKE 250 NEW RANDOM NUMBERS *****
      SUBROUTINE RANDN2(IQR1,IQR2)
      IMPLICIT REAL*8(A-G,O-Z)
      INTEGER IQR1(250),IQR2(250)
*
      DO 10 I=1,147
      IQR2(I)=IEOR(IQR1(I),IQR1(I+103))
 10   CONTINUE
*
      DO 20 I=1,147
      IQR1(I)=IQR2(I)
 20   CONTINUE
*
      DO 30 I=1,103
      I2=I+147
      IQR2(I2)=IEOR(IQR1(I2),IQR1(I))
 30   CONTINUE
*
      DO 40 I=148,250
      IQR1(I)=IQR2(I)
  40  CONTINUE
*
      RETURN
      END
C======================================================
      SUBROUTINE RANDN3(IQ)
*            MAKE NEW RANDOM NUMBER BY THE SEED IQ
      IMPLICIT REAL*8(A-G,O-Z)
      IQ=IQ*48828125
      IF(IQ) 10,20,20
  10  IQ=(IQ+2147483647)+1
  20  RETURN
      END
C======================================================
      real*8 function  ranf(iq)

      implicit none
      integer  iq
      integer  IP,IR
      parameter (IP = 48828125, IR = 2147483647)

      iq=iq* IP
c       print *,'iq = ',iq
      if(iq) 10,20,20
  10  iq=(iq+IR)+1
  20  ranf = dble(iq)/2.0D0**31
      return
      end
C======================================================
      real*8 function  rang(iq)
      
      integer iq
      real*8 pi, ranmod, theta, ranf
      
      PI = 3.1415926536D0
      RANMOD = SQRT(-2.D0 * LOG(RANF(iq)))
      THETA  = 2.D0 * PI * RANF(iq)
      rang = RANMOD * COS(THETA)

      return
      end
C======================================================
