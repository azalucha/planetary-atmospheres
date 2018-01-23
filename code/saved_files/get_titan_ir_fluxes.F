#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_TITAN_IR_FLUXES
C     !INTERFACE:

      SUBROUTINE GET_TITAN_IR_FLUXES(i,j,bi,bj,iPlus,iMinus,myTime)

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_TITAN_IR_FLUXES
C     | o Calculate IR fluxes, iPlus and iMinus 
C     |   using Robinson and Catling (2012) analytic RT
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
C     i,j       :: Indices on tile
C     iMinus    :: Downward IR flux
C     iPlus     :: Upward IR flux (atmosphere + surface)


      INTEGER bi, bj, k,i,j,ki,ks
      _RL taucabC(Nr),taucabF(Nr+1),myTime
      _RL termP(Nr),iMinus1,gammq1,gammq2
      _RL iPlus(Nr),ss,iMinustaurc,albedo,inpFlux
      _RL iMinus(Nr),betat,tau(Nr),n,gamma
      _RL lssol,tfrostsamz,qofac,D,k1,k2,alpha
      _RL Too,poo,f,fsunch1,fsunch2,taurc,prc



C     Begin Equations

      ks=kSurfC(i,j,bi,bj)


      IF (lsIsConstant .EQ. 0) THEN
       CALL GET_LS(myTime,lssol)
      ELSE
       lssol = lsperpetual
      ENDIF


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
      Too=surfaceT(i,j,bi,bj)
      poo=Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj)
      f=0.15
      n=2.
      tauoo=5.4*poo/atm_po
      taurc=4.0
      prc=1.3E5
      betat=alpha*(gamma-1.)/gamma

      IF(selectAddFluid.EQ.0) THEN
       DO k=1,Nr
        termP(k)=rC(k)
       ENDDO
      ELSE
       DO k=1,Nr
        termP(k)=rC(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &           /Ro_surf(i,j,bi,bj)
       ENDDO
      ENDIF



      DO k=1,Nr
       tau(k)=taurc*(f*termP(k)/prc+(1-f)*(termP(k)/prc)**n)
      ENDDO

      iMinustaurc=fsunch1/2.*(1.+D/k1-(1.+D/k1)*exp(-k1*taurc))
     &           +fsunch2/2.*(1.+D/k2-(1.+D/k2)*exp(-k2*taurc))

C Calculate Fluxes
      DO k=ks,Nr
C radiative regime
       IF(tau(k) .LE. taurc) THEN
        iPlus(k)=fsunch1/2.*(1.+D/k1+(1.-D/k1)*exp(-k1*tau(k)))
     &          +fsunch2/2.*(1.+D/k2+(1.-D/k2)*exp(-k2*tau(k)))
        iMinus(k)=fsunch1/2.*(1.+D/k1-(1.+D/k1)*exp(-k1*tau(k)))
     &           +fsunch2/2.*(1.+D/k2-(1.+D/k2)*exp(-k2*tau(k)))
       ELSE
        CALL GAMMQ(1.+4.*betat/n,D*tau(k),gammq1)
        CALL GAMMQ(1.+4.*betat/n,D*tauoo,gammq2)
        iPlus(k)=stephanBoltzmannConstant*Too**4.*exp(D*tau(k))
     &   *(exp(-D*tauoo)+1./(D*tauoo)**(4.*betat/n)
     &   *(gammq1-gammq2))
        CALL GET_INTEGRAL_TITAN(taurc,tau(k),tauoo,betat,n,D,ss)
        iMinus=iMinustaurc*exp(-D*(tau(k)-taurc))+
     &          D*stephanBoltzmannConstant*Too**4.*ss
       ENDIF
      ENDDO

      RETURN
      END
