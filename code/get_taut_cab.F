#include "CPP_OPTIONS.h" 

C     !ROUTINE: GET_TAUT_CAB
C     !INTERFACE:  

      SUBROUTINE GET_TAUT_CAB(taut,i,j,bi,bj)
      
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_TAUT_CAB
C     | o Calculate optical depth of the top of the convective 
C     |   layer in radiative-convective equilibrium using
C     |   Caballero et al. 2008 scheme
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
      _RL flux_integral,flux_integral_der,guess0

C     begin equations

      guess0=tauinfcab/2.
      tol = 0.00001
      xamz = guess0
      erroramz = 2.*tol
      niters = 1   
      maxiters = 10

      DO WHILE(erroramz .GT. tol .AND. niters .LT. maxiters)
       CALL GET_INTEGRAL_CAB(0.,xamz/betacab,flux_integral,1)
       CALL GET_INTEGRAL_CAB(0.,xamz/betacab,flux_integral_der,2) 
       currentfuncvalue=flux_integral+emissivityOfSs*
     &                   EXP(-xamz/betacab)-
     &                   (2.+(tauinfcab-xamz)/betacab)/
     &                   (1.+(tauinfcab-xamz)/betacab)*
     &                   (1.-xamz/tauinfcab)**(4.*atm_kappa)
       currentfuncder=flux_integral_der-emissivityOfSs
     &                   /betacab*EXP(-xamz/betacab)-
     &                 (2.+(tauinfcab-xamz)/betacab)/
     &                 (1.+(tauinfcab-xamz)/betacab)*
     &                 (1.-xamz/tauinfcab)**(4.*atm_kappa-1.)
     &                 *(-1./tauinfcab)*4.*atm_kappa-
     &                 (1.-xamz/tauinfcab)**(4.*atm_kappa)
     &                  /(1.+(tauinfcab-xamz)/betacab)**2.
     &                 *(-1./betacab)*((1.+(tauinfcab-xamz)
     &                 /betacab)-
     &                 (2.+(tauinfcab-xamz)/betacab))
C amz Leibniz rule part
     &                 +(1.-xamz/tauinfcab)/betacab
       xamz = xamz - currentfuncvalue/currentfuncder
       erroramz = ABS(currentfuncvalue/currentfuncder)
       niters = niters + 1
      ENDDO

C     Solution for taut
      taut = xamz
C     write(*,*) taut

      RETURN
      END
