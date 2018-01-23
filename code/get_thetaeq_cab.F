#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_THETAEQ_CAB
C     !INTERFACE:      

      SUBROUTINE GET_THETAEQ_CAB(thetaEq,termP,bi,bj,myTime,kLev)

C     !DESCRIPTION: 
C     *==========================================================*
C     | SUBROUTINE GET_THETAEQ_CAB
C     | o Calculate radiative equilibrium potential
C     |  temperature for HS94 forcing, using Cabellero et al. 2008
C     |  radiative-convective scheme
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
      _RL trad,qo
      _RL teff,taut
      _RL taucab,tautcab

C     Begin equations

      IF (lsIsConstant .EQ. 0) THEN
       CALL GET_LS(myTime,lssol)
      ELSE
       lssol = lsperpetual
      ENDIF

      DO j=1,sNy
       DO i=1,sNx

C     Calculate the insolation, qo, in W/m^2
        CALL GET_QO(i,j,bi,bj,lssol,myTime,solarconstant,
     &     albedoany(i,j,bi,bj),qo)

C     Calculate optical depth in the column
        IF(selectAddFluid.EQ.0) THEN
         taucab=tauinfcab*(1.-termP(i,j)/Ro_surf(i,j,bi,bj))
        ELSE
         taucab=tauinfcab*(1.-termP(i,j)/(Ro_surf(i,j,bi,bj)
     &                   +etaN(i,j,bi,bj)))
        ENDIF

C For the case of selectAddFluid=1, when there is ice on the ground, useSurfaceTflag
C is set to 1.  In this case, we calculate Teq using the ground as a source of flux.
        IF(useSurfaceTflag(i,j,bi,bj) .EQ. 1.) THEN
         teff=1./(1./referenceTempInFrostEqn - atm_Rd/
     &             atm_L*LOG((Ro_surf(i,j,bi,bj)+
     &             etaN(i,j,bi,bj))/referencePresInFrostEqn))
     &            *emissivity(i,j,bi,bj)
        ELSE
         teff=(qo/stephanBoltzmannConstant)**0.25
        ENDIF
C     Thermodynamics loop calls kLev starting at Nr; taut does not depend on the level
        IF(kLev .EQ. Nr) THEN
C     We only want to recalculate at t=0 and every radTimeStep
         IF(tautall(i,j,bi,bj).EQ.-1. .OR. 
     &            MOD(myTime,radTimeStep).EQ.0.) THEN
          CALL GET_TAUT_F_CAB(taut,i,j,bi,bj)
          tautall(i,j,bi,bj)=taut
         ENDIF
        ENDIF
        IF(taucab .GE. tautall(i,j,bi,bj)) THEN
C Above the convective layer
         trad=teff*((betacab+tauinfcab-taucab)/
     &            (2.*betacab+(1.-betacab)*tauinfcab))**0.25
        ELSE
C Below the convective layer
         trad=teff*((betacab+
     &            tauinfcab-tautall(i,j,bi,bj))/
     &            (2.*betacab+(1.-betacab)*tauinfcab))**0.25
     &          *((1.-taucab/tauinfcab)/(1.-tautall(i,j,bi,bj)
     &               /tauinfcab))**atm_kappa
        ENDIF
        thetaEq(I,J) = trad / ((termP(I,J)/atm_po)**atm_kappa)
        thetaEqSave(i,j,kLev,bi,bj)=thetaEq(I,J)
       ENDDO
      ENDDO

      RETURN
      END
