#include "CPP_OPTIONS.h"

      SUBROUTINE TITANNRAD(bi,bj,kLev,myTime,gt_arr)
C     *==========================================================*
C     | S/R TITANNRAD
C     | o Computes radiation for Newtonian relaxation to Robinson and Catling (2012)
C     |   atmosphere
C     *==========================================================*

      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

      INTEGER kLev, bi, bj
      _RL myTime
      INTEGER i, j,k

      _RL thetaLim,kT,ka,ks,sigma_b
      _RL thetaEq,inpFlux,albedo,qofac
      _RL titanT,betat,tau(Nr)
      _RL tfrostamz,lssol,tfrostsamz,Too,poo,f
      _RL fsunch1,fsunch2,taurc,prc,gamma,n
      _RL deltau(Nr),visdown,temp
      _RL alpha,D,k1,k2
      _RL gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)



C--   Forcing term(s)
      DO j=1,sNy
       DO i=1,sNx

        IF(kLev .EQ. Nr) THEN
         CALL DO_SURFACE_CALC(i,j,bi,bj,myTime)
        ENDIF



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
 
 
        DO k=1,Nr
         tau(k)=taurc*(f*rC(k)/prc+(1.-f)*(rC(k)/prc)**n)
         deltau(k)=delR(k)*tauoo/poo
        ENDDO

C  thermodynamics loop counts down from Nr.  Do calculations
C  for entire column at first instance
        IF(kLev .EQ. Nr) THEN
         IF(initrad(i,j,bi,bj) .EQ. 1. .OR.
     &            MOD(myTime,radTimeStep).EQ.0.) THEN
          initrad(i,j,bi,bj)=0.

          DO k=1,Nr
           IF (tau(k) .LE. taurc) THEN
C radiative regime
            titanT=((fsunch1/2.*(1.+D/k1+(k1/D-D/k1)*
     &            EXP(-k1*tau(k)))
     &           +fsunch2/2.*(1.+D/k2+(k2/D-D/k2)*
     &           EXP(-k2*tau(k))))
     &           /stephanBoltzmannConstant)**0.25
           ELSE
C convective regime
            titanT=Too*(rC(k)/atm_po)**(betat/n)
           ENDIF
           thetaEq=titanT*(rC(k)/atm_po)**(-atm_kappa)

           tfrostamz=1./(1./referenceTempInFrostEqn - atm_Rd/atm_L*
     &            LOG(rC(k)/referencePresInFrostEqn))
     &          / ((rC(k)/atm_po)**atm_kappa)
           IF (thetaEq .LT. tfrostamz) THEN
            thetaEq = tfrostamz
           ENDIF
           thetaEqSave(i,j,k,bi,bj)=thetaEq
          ENDDO
         ENDIF
        ENDIF

C--   Forcing term(s)
        ka=1./hstaua(i,j,bi,bj)
        ks=1./hstaus(i,j,bi,bj)
        sigma_b = 0.7 _d 0
        kT=ka+(ks-ka)
     &      *MAX(0. _d 0,
     &       (rC(kLev)*recip_Rcol(i,j,bi,bj)-sigma_b)/
     &                 (1. _d 0-sigma_b) )
     &      *COS((yC(i,j,bi,bj)*deg2rad))**4.
           gT_arr(i,j)=gT_arr(i,j)
     &        - kT*( theta(i,j,kLev,bi,bj)-
     &            thetaEqSave(i,j,kLev,bi,bj) )
     &            *maskC(i,j,kLev,bi,bj)

       ENDDO
      ENDDO


      RETURN
      END
