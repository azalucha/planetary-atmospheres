#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_THETAEQ
C     !INTERFACE:      

      SUBROUTINE GET_THETAEQ(thetaEq,termP,bi,bj,myTime,kLev)

C     !DESCRIPTION: 
C     *==========================================================*
C     | SUBROUTINE GET_THETAEQ
C     | o Calculate radiative equilibrium potentialtemperature for
C     |  HS94 forcing
C     |  units: K
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
C     termP     :: Current pressure (r units)
C     thetaEq   :: radiative equilibrium temperature
C     myTime    :: Current time in simulation

      INTEGER bi, bj, kLev
      _RL termP(sNx,sNy),thetaEq(sNx,sNy),myTime
        
C     !LOCAL VARIABLES:
C     == Local variables ==

      INTEGER i,j
      _RL lssol
      _RL trad,tau,qo,tauo
      _RL taut,teqsurf,teff

C     Begin equations

      IF (lsIsConstant .EQ. 0) THEN
       CALL GET_LS(myTime,lssol)
      ELSE
       lssol = lsperpetual
      ENDIF

      DO j=1,sNy
       DO i=1,sNx

        tauo=tauoo*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))/atm_po 

C     Calculate optical depth
        tau = tauoo*(termP(i,j)/atm_po)

C For the case of selectAddFluid=1, when there is ice on the ground, useSurfaceTflag
C is set to 1.  In this case, we calculate Teq using the ground as a source of flux.
        IF(useSurfaceTflag(i,j,bi,bj) .EQ. 1.) THEN
         teff=1./(1./referenceTempInFrostEqn - atm_Rd/
     &             atm_L*LOG((Ro_surf(i,j,bi,bj)+
     &             etaN(i,j,bi,bj))/referencePresInFrostEqn))
     &            *emissivity(i,j,bi,bj)
        ELSE
         CALL GET_QO(i,j,bi,bj,lssol,myTime,solarconstant,
     &     albedoany(i,j,bi,bj),qo)      
         teff=(qo/stephanBoltzmannConstant)**0.25 
        ENDIF
C     pure radiative case
        IF(marsPureRad) THEN
         trad=teff*(0.5+0.75*tau)**0.25
C     radiative convective case
        ELSE
C     Thermodynamics loop calls kLev starting at Nr; taut does not depend on the level
         IF(kLev .EQ. Nr) THEN
C     We only want to recalculate at t=0 and every radTimeStep
          IF(tautall(i,j,bi,bj).EQ.-1. .OR.
     &            MOD(myTime,radTimeStep).EQ.0.) THEN
           CALL GET_TAUT_F(taut,i,j,bi,bj)
           tautall(i,j,bi,bj)=taut
          ENDIF
         ENDIF
         IF(tau .LE. tautall(i,j,bi,bj)) THEN
C Above the convective layer
          trad=teff*(0.5+0.75*tau)**0.25
         ELSE
C Below the convective layer
          trad=teff*(0.5 +0.75*tautall(i,j,bi,bj))**0.25
     &            *(tau/tautall(i,j,bi,bj))**atm_kappa
         ENDIF
        ENDIF
        thetaEq(I,J) = trad / ((termP(I,J)/atm_po)**atm_kappa)
        thetaEqSave(i,j,kLev,bi,bj)=thetaEq(I,J)
       ENDDO
      ENDDO

      RETURN
      END
