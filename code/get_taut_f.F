#include "CPP_OPTIONS.h" 

C     !ROUTINE: GET_TAUT_F
C     !INTERFACE:  

      SUBROUTINE GET_TAUT_F(taut,i,j,bi,bj)
      
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_TAUT
C     | o Calculate optical depth of the top of the convective 
C     |   layer in radiative-convective equilibrium, using surface
C     |   temperature
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
      _RL tauo,guess0

C     begin equations

      tauo=tauoo*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))/atm_po
      guess0 = tauoo/2.   

      tol = 0.00001
      xamz = guess0
      erroramz = 2.*tol
      niters = 1
      maxiters = 10


      DO WHILE(erroramz .GT. tol .AND. niters .LT. maxiters)
       currentfuncvalue=(emissivity(i,j,bi,bj)*(0.5+
     &                    0.75*xamz)/(1.+0.75*tauo))**0.25
     &                  -(xamz/tauo)**atm_kappa 
       currentfuncder=3./16.*(emissivity(i,j,bi,bj)/
     &                 (0.5+0.75*xamz)**3./(1.+
     &                  0.75*tauo))**0.25
     &                 -atm_kappa/xamz*(xamz/tauo)**atm_kappa
       xamz = xamz - currentfuncvalue/currentfuncder
       erroramz = ABS(currentfuncvalue/currentfuncder)
       niters = niters + 1
      ENDDO

C     Solution for taut, in radians
      taut = xamz
C     write(*,*) taut

      RETURN
      END
