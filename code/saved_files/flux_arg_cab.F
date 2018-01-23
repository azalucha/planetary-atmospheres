#include "CPP_OPTIONS.h"

C     !ROUTINE: FLUX_ARG_CAB
C     !INTERFACE:
        
      SUBROUTINE FLUX_ARG_CAB(funcx,b,x,integral_code)
      
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE FLUX_ARG_CAB
C     | o Argument of terms under the integral in the flux balance
C     |   equation for the RC equilibrium model of Caballero
C     |   et al. 2008
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
       funcx=(1.-x*betacab/tauinfcab)**(4.*atm_kappa)*EXP(x-b)
      ELSE
       funcx=(1.-x*betacab/tauinfcab)*EXP(x-b)*(-1./betacab)
      ENDIF

      RETURN
      END
