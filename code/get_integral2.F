#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_INTEGRAL2
C     !INTERFACE:

      SUBROUTINE GET_INTEGRAL2(tauo,xamz,flux_integral,integral_code)

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_INTEGRAL2
C     | o calculate the numerical integral in the flux equation,
C     |   using analytical expression (calls gamma function)
C     |  units:
C     |  author: amz
C     *==========================================================*

      IMPLICIT NONE

C     == Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     tauo           :: lower limit on integral (optical depth at surface)
C     xamz           :: upper limit on integral (optical depth a top of convective layer)
C     flux_integral  :: value of integral (output)
C     integral_code  :: 1 for flux integral, 2 for flux integral derivative

      INTEGER integral_code
      _RL tauo,xamz,flux_integral

C     !LOCAL VARIABLES:
C     == Local variables ==
      _RL secondargo,firstarg,secondargt,gammao,gammat

C     begin equations

      firstarg=1. + (4.* atm_Rd)/atm_Cp
      secondargo=(3.*tauo)/2.
      secondargt=(3.*xamz)/2.
      CALL gammq(firstarg,secondargo,gammao)
      CALL gammq(firstarg,secondargt,gammat)
      IF(integral_code .EQ. 1) THEN
       flux_integral=(2./3.)**(4.*atm_Rd/atm_Cp)*EXP(3.*xamz/2.)
     &                *xamz**(-4.*atm_Rd/atm_Cp)*
     &                (1./3.+0.5*xamz)*(-gammao+gammat)-
     &                (2./3.+EXP(-3./2.*(tauo - xamz))*
     &                (-2./3. - 0.5*tauo) + 0.5*xamz)
      ELSE
       flux_integral=(2./3.)**(4.*atm_Rd/atm_Cp)*EXP(3.*xamz/2.)*
     &                xamz**(-4.*atm_Rd/atm_Cp)*
     &                (-(2./3.)**(-1.-4.*atm_Rd/atm_Cp)*EXP(-3.*xamz/2.)
     &                *(1./3.+0.5*xamz)*
     &                xamz**(4.*atm_Rd/atm_Cp)-0.5*gammao+0.5*gammat)+
     &                xamz**(-4.*atm_Rd/atm_Cp)
     &                *EXP(3.*xamz/2.)*(1./3.+0.5*xamz)*(-gammao+gammat)
     &                *(3./2.)**(-4.*atm_Rd/atm_Cp)
     &                *(4.*atm_Rd/atm_Cp/xamz+3./2.)-(0.5-3./2.*
     &                EXP(-3./2.*(tauo - xamz))*(2./3.+0.5*tauo))
      ENDIF

      RETURN
      END
