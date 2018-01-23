#include "CPP_OPTIONS.h"
      
C     !ROUTINE: TRAPZD
C     !INTERFACE:

      SUBROUTINE TRAPZD(a,b,s,n,integral_code)
      
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE TRAPZD
C     | o TRAPZD from numerical recipes
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
C     s              :: value of integral (output)
C     integral_code  :: 1 for flux integral, 2 for flux integral derivative

      INTEGER n,integral_code
      _RL a,b,s 

C     !LOCAL VARIABLES:
C     == Local variables ==
      integer it,j
      _RL del,sum,tnm,x,funca,funcb,funcx


C     begin equations

      IF (n .EQ. 1) THEN
        CALL FLUX_ARG(funca,b,a,integral_code)
        CALL FLUX_ARG(funcb,b,b,integral_code)
       s=0.5*(b-a)*(funca+funcb)
      ELSE
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       DO j=1,it
        CALL FLUX_ARG(funcx,b,x,integral_code)
        sum=sum+funcx
        x=x+del
       ENDDO
       s=0.5*(s+(b-a)*sum/tnm)
      ENDIF

      RETURN
      END        
