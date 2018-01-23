#include "CPP_OPTIONS.h"

C     !ROUTINE: STEP_SURFACE_T_2
C     !INTERFACE:

      SUBROUTINE STEP_SURFACE_T_2(i,j,bi,bj,irdown,visdown)             

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE STEP_SURFACE_T
C     | o Calculate surface temperature in conjunction with NASA
C     |  Ames radiation scheme (radcode 9)
C     |  author: amz
C     *==========================================================*

      IMPLICIT NONE

C     == Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "AMZVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     bi,bj     :: Current tile indices
C     i,j       :: Indices on tile
C     irdown    :: Downward component of IR radiation at surface
C                  (lowest model level)
C     visdown   :: Downward component of solar radiation at surface
C                  (lowest model level)

      INTEGER bi, bj,i,j
      _RL visdown,irdown,cp,rho,cpOfIce,iceThickness

C     Begin equations



      IF(selectAddFluid .EQ. 0) THEN
        surfaceT(i,j,bi,bj)=surfaceTPrime(i,j,bi,bj)+
     &     radTimeStep/cpOfSoil/soilThickness
     &    /rhoOfSoil*(visdown
     &        +irdown-emissivity(i,j,bi,bj)*
     &     stephanBoltzmannConstant*
     &     surfaceTPrime(i,j,bi,bj)**4.)
      ELSE
       iceThickness=massonground(i,j,bi,bj)/
     &   rA(i,j,bi,bj)*deltaT/rhoOfIce
C no ice on ground (surface layer all soil)
       IF (iceThickness .LE. 0.) THEN
        surfaceT(i,j,bi,bj)=surfaceTPrime(i,j,bi,bj)+
     &     radTimeStep/cpOfSoil/soilThickness
     &    /rhoOfSoil*(visdown
     &        +irdown-emissivity(i,j,bi,bj)*
     &     stephanBoltzmannConstant*
     &     surfaceTPrime(i,j,bi,bj)**4.)
       ELSE
        cpOfIce=349.+ 4.8*surfaceTPrime(i,j,bi,bj)
C lots of ice on ground	(surface layer all ice)
        IF (iceThickness .GE. soilThickness) THEN
         surfaceT(i,j,bi,bj)=surfaceTPrime(i,j,bi,bj)+
     &     radTimeStep/cpOfIce/soilThickness
     &    /rhoOfIce*(visdown
     &        +irdown-emissivity(i,j,bi,bj)*
     &     stephanBoltzmannConstant*
     &     surfaceTPrime(i,j,bi,bj)**4.)
        ELSE
C some ice on top of some soil
         cp=cpOfIce*(iceThickness/soilThickness)+
     &       cpOfSoil*((soilThickness-iceThickness)/soilThickness)
         rho=rhoOfIce*(iceThickness/soilThickness)+
     & 	    rhoOfSoil*((soilThickness-iceThickness)/soilThickness)
         surfaceT(i,j,bi,bj)=surfaceTPrime(i,j,bi,bj)+
     &     radTimeStep/cp/soilThickness
     &     /rho*(visdown
     &        +irdown-emissivity(i,j,bi,bj)*
     &     stephanBoltzmannConstant*
     &     surfaceTPrime(i,j,bi,bj)**4.)
        ENDIF
       ENDIF
      ENDIF 


      RETURN
      END
