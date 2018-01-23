#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_INTEGRAL
C     !INTERFACE:
      
      SUBROUTINE GET_INTEGRAL(a,b,ss,integral_code)
     
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_INTEGRAL
C     | o calculate the numerical integral in the flux equation.
C     | o QROMB from numerical recipes
C     |  units: 
C     |  author: amz
C     *==========================================================*
      
      IMPLICIT NONE
      
C     == Global variables
#include "EEPARAMS.h"
     
C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     a              :: lower limit on integral
C     b              :: upper limit on integral
C     ss             :: value of integral (output)
C     integral_code  :: 1 for flux integral, 2 for flux integral derivative
       
      INTEGER integral_code
      _RL a,b,ss

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER JMAX,JMAXP,K,KM,j
      _RL EPS
      PARAMETER (EPS=1. _d -2, JMAX=10, JMAXP=JMAX+1, K=2, KM=K-1)
      _RL dss,h(JMAXP),s(JMAXP)

C     begin equations

      h(1)=1.
      DO j=1,JMAX
       CALL TRAPZD(a,b,s(j),j,integral_code)
       IF(j .GE. K) THEN
        CALL POLINT(h(j-KM),s(j-KM),K,0.,ss,dss)
        IF (ABS(dss) .LE. EPS*ABS(ss)) RETURN
       ENDIF
       s(j+1)=s(j)
       h(j+1)=0.25*h(j)
      ENDDO
      write (*,*) 'too many steps in qromb'
      RETURN
      END

