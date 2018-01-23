#include "CPP_OPTIONS.h"

C     !ROUTINE: GCF
C     !INTERFACE:

      SUBROUTINE GCF(gammcf,a,x,gln)

C     !DESCRIPTION:
C     *==========================================================*
C     | GCF
C     | o From numerical recipes book
C     |  units:
C     |  author: amz
C     *==========================================================*

      IMPLICIT NONE

C     == Global variables
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C

      _RL a,gammcf,gln,x

C     !LOCAL VARIABLES:
C     == Local variables ==
      _RL EPS,FPMIN,an,b,c,d,del,h
      INTEGER ITMAX,i
      PARAMETER(ITMAX=100,EPS=3. _d -7,FPMIN=1. _d -30)

C     begin equations
      CALL GAMMLN(a,gln)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      DO i=1,ITMAX
       an=-i*(i-a)
       b=b+2.
       d=an*d+b
       IF(ABS(d).LT.FPMIN) d=FPMIN
       c=b+an/c
       IF(ABS(c).LT.FPMIN) c=FPMIN
       d=1./d
       del=d*c
       h=h*del
       IF(ABS(del-1.).LT.EPS) GOTO 1
      ENDDO
      write(*,*) 'a too large, ITMAX too small in gcf'
    1 gammcf=EXP(-x+a*LOG(x)-gln)*h

      RETURN
      END
