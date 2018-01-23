#include "CPP_OPTIONS.h"

C     !ROUTINE: GSER
C     !INTERFACE:
      
      SUBROUTINE GSER(gamser,a,x,gln)

C     !DESCRIPTION:
C     *==========================================================*
C     | GSER
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
   
      _RL a,gamser,gln,x

C     !LOCAL VARIABLES:
C     == Local variables ==
      _RL EPS,ap,del,sum
      INTEGER ITMAX,n
      PARAMETER(ITMAX=100,EPS=3. _d -7)


C     begin equations
      CALL GAMMLN(a,gln)
      IF(x .LE. 0.) THEN
       IF(x .LT. 0.) write(*,*) 'x < 0 in gser'
       gamser=0.
       RETURN
      ENDIF
      ap=a
      sum=1./a
      del=sum
      DO n=1,ITMAX
       ap=ap+1.
       del=del*x/ap
       sum=sum+del
       IF(ABS(del).LT.ABS(sum)*EPS) GOTO 1
      ENDDO
      write(*,*) 'a too large, ITMAX too small in gser'
    1 gamser=sum*EXP(-x+a*LOG(x)-gln)

      RETURN
      END
