#include "CPP_OPTIONS.h"

C     !ROUTINE: FLUX_ARG
C     !INTERFACE:
        
      SUBROUTINE FLUX_ARG(funcx,b,x,integral_code)
      
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE FLUX_ARG
C     | o Argument of terms under the integral in the flux balance
C     | o equation for the RC equilibrium model
C     |  units:
C     |  author: amz
C     *==========================================================*
       
      IMPLICIT NONE
      
C     == Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "AMZVARS.h"
        
C     !INPUT/OUTPUT PARAMETERS: 
C     == Routine arguments ==
C     x              :: value at which to evaluate argument
C     b              :: either taut (gray scheme) or taut/betacab (semi-gray scheme)
C                       This is a convenient way of getting it to this subroutine,
C                       since it so happens that the limit is also used in the argument
C     funcx          :: value of argument at x (output)
C     integral_code  :: 1 for flux integral, 2 for flux integral derivative
      
      INTEGER integral_code
      _RL x,b,funcx
       
C     !LOCAL VARIABLES:
C     == Local variables ==

       
C     begin equations
      IF(integral_code .EQ. 1) THEN
        funcx = EXP(-(x-b)*1.5)*(-(0.5+0.75*b)*(x/b)**(4.*atm_kappa)
     &          +(0.5+0.75*x))
      ELSE
        funcx = EXP(-(x-b)*1.5)*(-4./b**2.*atm_kappa*x*(-0.5 - 0.75*b)*
     &         (x/b)**(4.*atm_kappa-1)-0.75*(x/b)**(4.*atm_kappa))+
     &         1.5*EXP(-(x-b)*1.5)*(0.5+0.75*x+(-0.5 - 0.75*b)
     &         *(x/b)**(4.*atm_kappa))
      ENDIF

      RETURN
      END
