#include "CPP_OPTIONS.h" 

C     !ROUTINE: GET_TAUT
C     !INTERFACE:  

      SUBROUTINE GET_TAUT(taut,i,j,bi,bj)
      
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_TAUT
C     | o Calculate optical depth of the top of the convective 
C     |   layer in radiative-convective equilibrium
C     |  units: none
C     |  author: amz
C     *==========================================================*
      
      IMPLICIT NONE
      
C     == Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "AMZVARS.h"
      
C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     bi,bj     :: Current tile indices
C     i,j       :: loop counters on tile
C     taut      :: Optical depth of the top of the convective layer
      
      INTEGER bi, bj, i, j
      _RL taut
      
C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER niters, maxiters
      _RL guess,tol,xamz,erroramz,currentfuncvalue,currentfuncder
      _RL tauo,flux_integral,flux_integral_der,guess0

C     begin equations

       IF(selectAddFluid.EQ.1) THEN
        tauo=tauoo*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))/atm_po
       ELSE
        tauo=tauoo*Ro_surf(i,j,bi,bj)/atm_po
       ENDIF
        guess0 = tauo/2.   

      tol = 0.00001
      xamz = guess0
      erroramz = 2.*tol
      niters = 1   
      maxiters = 10

      DO WHILE(erroramz .GT. tol .AND. niters .LT. maxiters)
       CALL GET_INTEGRAL(tauo,xamz,flux_integral,1)
C      CALL GET_INTEGRAL2(tauo,xamz,flux_integral,1)
       currentfuncvalue = flux_integral +emissivityOfSs*
     &                     ((0.5 + 0.75*xamz)*(tauo/xamz)
     &                      **(4.*atm_kappa)
     &                     -(1. + 0.75*tauo))*EXP(-(tauo -
     &                       xamz)*1.5)
       CALL GET_INTEGRAL(tauo,xamz,flux_integral_der,2)
C      CALL GET_INTEGRAL2(tauo,xamz,flux_integral_der,2)
       currentfuncder = flux_integral_der+emissivityOfSs*
     &                  1.5*EXP(-(tauo - xamz)*1.5)*
     &                  (-1.-0.75*tauo+(0.5+0.75*xamz)*
     &                  (tauo/xamz)**(4.*atm_kappa))
     &                  +emissivityOfSs*EXP(-(tauo 
     &                  - xamz)*1.5)*
     &                  (0.75*(tauo/xamz)**(4.*atm_kappa)-
     &                   4.*atm_kappa*tauo/
     &                  xamz**2.*(0.5+0.75*xamz)*(tauo/xamz)
     &                  **(4.*atm_kappa-1.))
       xamz = xamz - currentfuncvalue/currentfuncder
       erroramz = ABS(currentfuncvalue/currentfuncder)
       niters = niters + 1
      ENDDO

C     Solution for taut, in radians
      taut = xamz
C     write(*,*) taut

      RETURN
      END
