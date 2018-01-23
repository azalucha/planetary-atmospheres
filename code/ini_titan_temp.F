#include "CPP_OPTIONS.h"

      SUBROUTINE INI_TITAN_TEMP(i,j,bi,bj,titanT,Too)
C     *==========================================================*
C     | S/R INI_TITAN_TEMP
C     | o Calculate radiative equilibrium temperature to initialize
C     |   Titan temperatures following Robinson and Catling (2012)
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

      INTEGER  bi, bj,k,i,j
      _RL titanT(Nr),TitanSurfT,betat,tau
      _RL lssol,tfrostsamz,qofac,D,k1,k2,alpha
      _RL Too,poo,f,fsunch1,fsunch2,taurc,prc
      _RL gamma,n,myTime,inpFlux,albedo

      IF (lsIsConstant .EQ. 0) THEN
       lssol=lsStart
      ELSE
       lssol=lsperpetual
      ENDIF

      myTime=0.E0
      inpFlux=1.E0
      albedo=0.E0

      CALL GET_QO(i,j,bi,bj,lssol,myTime,inpFlux,albedo,qofac)

C we want to apply our own geometry correction, instead of using
C disk averaged solar constant
      fsunch1=1.5*4.*qofac
      fsunch2=1.1*4.*qofac

      D=1.66
      k1=120.
      k2=0.2
      alpha=0.77
      gamma=7./5.
      poo=atm_po
      f=0.15
      n=2.
      tauoo=5.4
      taurc=4.0
      prc=1.3E5
      betat=alpha*(gamma-1.)/gamma

      DO k=1,Nr
       tau=taurc*(f*rC(k)/prc+(1.-f)*(rC(k)/prc)**n)
       IF (tau .LE. taurc) THEN
C radiative regime
        titanT(k)=((fsunch1/2.*(1.+D/k1+(k1/D-D/k1)*EXP(-k1*tau))
     &           +fsunch2/2.*(1.+D/k2+(k2/D-D/k2)*EXP(-k2*tau)))
     &           /stephanBoltzmannConstant)**0.25
       ELSE
C convective regime
        titanT(k)=Too*(tau/tauoo)**(betat/n)
       ENDIF
      ENDDO
       



      RETURN
      END
