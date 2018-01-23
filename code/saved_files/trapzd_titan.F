#include "CPP_OPTIONS.h"
      
C     !ROUTINE: TRAPZD
C     !INTERFACE:

      SUBROUTINE TRAPZD_TITAN(a,b,s,n,tauoo,beta,nn,D)
      
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE TRAPZD_TITAN
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

      INTEGER n,nn
      _RL a,b,s,beta,D,tauoo

C     !LOCAL VARIABLES:
C     == Local variables ==
      integer it,j
      _RL del,sum,tnm,x,funca,funcb,funcx


C     begin equations

      IF (n .EQ. 1) THEN
        funca=(a/tauoo)**(4.*beta/n)*exp(-D*(b-a))
        funcb=(b/tauoo)**(4.*beta/n)*exp(-D*(b-b))
       s=0.5*(b-a)*(funca+funcb)
      ELSE
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       DO j=1,it
        funcx=(x/tauoo)**(4.*beta/n)*exp(-D*(b-x))
        sum=sum+funcx
        x=x+del
       ENDDO
       s=0.5*(s+(b-a)*sum/tnm)
      ENDIF

      RETURN
      END        
