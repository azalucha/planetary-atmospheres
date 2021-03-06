#include "CPP_OPTIONS.h"

      SUBROUTINE DO_SURFACE_CALC(i,j,bi,bj,myTime)
C     *==========================================================*
C     | S/R DO_SURFACE_CALC
C     | o Surface scheme that includes calculation of surface 
C     |   temperature for all radcodes.  Allows for no volitile or
C     |   volatile cycle, in the case of the latter, calculates
C     |   surface frost mass and other volatile parameters.
C     |  NOTE: radcodes 2,8,11,12 NEED TO BE TESTED WITH THIS SCHEME.
C     |  NOTE: radcode 9 HAS NOT BEEN WRITTEN YET.
C     |  NOTE: TITAN IS ASSUMED TO HAVE FLAT TOPOGRAPHY
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

      INTEGER i,j,k,bi,bj,ks,nn
      _RL lssol,tfrostsurf,adf,termPAtkSurfC
      _RL tauAtkSurfC,termPs,tau,massneeded,myTime,qo
      _RL tfrostsamz,inpflux,albedo,fsunch1,fsunch2
      _RL qofac,k1,k2,poo,f,n,tauootitan,taurc,prc
      _RL tautitan(Nr),deltautitan(Nr),visdown,temp
      _RL termP(Nr),taucabC(Nr),taucabF(Nr+1)
      _RL downwardsolarflux,tauo,delTau(Nr)
      _RL deltaMdot

      nn=0

C For Newtonian relaxation specifications (2, 8, 12), instantaneous
C surface temperature is calculated (i.e. no thermal inertia).  There is
C not a scientific reason why a Newtonian-relaxed atmosphere could not
C have a time-varying or multilayer subsurface; the algorithm has just
C not been developed here.

C Heating rate specified directly (4, 6, 9, 10, 11) has a time-varying
C surface temperature.  Options are subsurface slab or multilayer subsurface.

C Fixed surface temperature option, for Pluto only
      IF(radcode .EQ. 10 .OR. radcode .EQ. 6) THEN
       IF(plutoSurfTFixedToFrost) THEN
        nn=nn+1
        tfrostsamz=1./(1./referenceTempInFrostEqn - atm_Rd/atm_L*
     &            LOG((Ro_Surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &            /referencePresInFrostEqn))
        surfaceT(i,j,bi,bj)=tfrostsamz
       ENDIF
      ENDIF

C Mars (and gray atmosphere)  Newtonian relaxation
      IF(radcode .EQ. 2 .OR. radcode .EQ. 8) THEN
       n=n+1
       IF (lsIsConstant .EQ. 0) THEN
        CALL GET_LS(myTime,lssol)
       ELSE
        lssol = lsperpetual
       ENDIF
       tfrostsurf=1./(1./referenceTempInFrostEqn -
     &                    atm_Rd/atm_L*LOG((Ro_surf(i,j,bi,bj)+
     &           	etaN(i,j,bi,bj))/referencePresInFrostEqn))
       ks=kSurfC(i,j,bi,bj)
       adf=0.
       termPAtkSurfC=rC(ks)*(Ro_surf(i,j,bi,bj)
     &                    +etaN(i,j,bi,bj))/Ro_surf(i,j,bi,bj)
       tauAtkSurfC=termPAtkSurfC/atm_Po*tauoo
       DO k=Nr,ks,-1
        termPs=rC(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &         /Ro_surf(i,j,bi,bj)
        tau=termPs/atm_Po*tauoo
        adf=adf+
     &        stephanBoltzmannConstant*(theta(i,j,k,bi,bj)
     &        *(termPs/atm_po)**atm_kappa)**4.*
     &        EXP((tau-tauAtkSurfC)/(2./3.))
     &        *drF(k)*tauoo/atm_Po*hFacC(i,j,k,bi,bj)
       ENDDO
       CALL GET_QO(i,j,bi,bj,lssol,myTime,solarconstant,
     &        albedoany(i,j,bi,bj),qo)
       surfaceT(i,j,bi,bj)=((adf+qo)
     &        /stephanBoltzmannConstant
     &        /emissivity(i,j,bi,bj))**(1./4.)
#ifndef ALLOW_ADDFLUID
        IF(surfaceT(i,j,bi,bj).LT.tfrostsurf) THEN
         albedoany(i,j,bi,bj)=albedoFrost
         emissivity(i,j,bi,bj)=emissivityOfFrost
         surfaceT(i,j,bi,bj)=tfrostsurf
         useSurfaceTflag(i,j,bi,bj)=1
        ELSE
         albedoany(i,j,bi,bj)=albedoSs
         emissivity(i,j,bi,bj)=emissivityOfSs
         useSurfaceTflag(i,j,bi,bj)=0
        ENDIF
#endif

#ifdef ALLOW_ADDFLUID
        IF(surfaceT(i,j,bi,bj).LT.tfrostsurf) THEN
         addMass(i,j,ks,bi,bj)=
     &            stephanBoltzmannConstant*emissivityOfFrost*
     &            (surfaceT(i,j,bi,bj)**4.-tfrostsurf**4.)*
     &            rA(i,j,bi,bj)/atm_L
         useSurfaceTflag(i,j,bi,bj)=1.
         massonground(i,j,bi,bj)=massonground(i,j,bi,bj)
     &                              -addMass(i,j,ks,bi,bj)
         albedoany(i,j,bi,bj)=albedoFrost
         emissivity(i,j,bi,bj)=emissivityOfFrost
         surfaceT(i,j,bi,bj)=tfrostsurf
        ELSEIF(massonground(i,j,bi,bj) .GT. 0.) THEN
         massneeded=stephanBoltzmannConstant*emissivityOfFrost*
     &         (surfaceT(i,j,bi,bj)**4.-tfrostsurf**4.)*
     &          rA(i,j,bi,bj)/atm_L
         IF(massneeded.LT.massonground(i,j,bi,bj)) THEN
          addMass(i,j,ks,bi,bj)=massneeded
          massonground(i,j,bi,bj)=
     &        massonground(i,j,bi,bj)-massneeded
          surfaceT(i,j,bi,bj)=tfrostsurf
         ELSE
          addMass(i,j,ks,bi,bj)=massonground(i,j,bi,bj)
          surfaceT(i,j,bi,bj)=(1./emissivityOfSs/
     &         stephanBoltzmannConstant*
     &         (emissivityOfFrost*stephanBoltzmannConstant*
     &         surfaceT(i,j,bi,bj)**4.-atm_L*
     &         massonground(i,j,bi,bj)/rA(i,j,bi,bj)))**0.25
          massonground(i,j,bi,bj)=0.
          albedoany(i,j,bi,bj)=albedoSs
          emissivity(i,j,bi,bj)=emissivityOfSs
         ENDIF
        ELSE
         useSurfaceTflag(i,j,bi,bj)=0.
        ENDIF
#endif
      ENDIF

C Titan
      IF(radcode .EQ. 11 .OR. radcode .EQ. 12) THEN
       nn=nn+1

       IF (lsIsConstant .EQ. 0) THEN
        CALL GET_LS(myTime,lssol)
       ELSE
        lssol = lsperpetual
       ENDIF

       inpFlux=1.E0
       albedo=0.E0

       CALL GET_QO(i,j,bi,bj,lssol,myTime,inpFlux,albedo,qofac)
       fsunch1=1.5*4.*qofac
       fsunch2=1.1*4.*qofac

       k1=120.
       k2=0.2
       poo=Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj)
       f=0.15
       n=2.
       tauootitan=5.4*poo/atm_po
       taurc=4.0
       prc=1.3E5

       DO k=1,Nr
        tautitan(k)=taurc*(f*rC(k)/prc+(1.-f)*(rC(k)/prc)**n)
        deltautitan(k)=delR(k)*tauootitan/poo
       ENDDO

       visdown=fsunch1*EXP(-k1*tauootitan)
     &           +fsunch2*EXP(-k2*tauootitan)
       adf=0.
       DO k=1,Nr
        temp=theta(i,j,k,bi,bj)*(tautitan(k)/tauootitan)**atm_kappa
        adf=adf+stephanBoltzmannConstant*
     &          temp**4.*EXP(tautitan(k)-tauootitan)*deltautitan(k)
       ENDDO


       IF(radcode .EQ. 12) THEN
C Newtonian relaxation case
        surfaceT(i,j,bi,bj)=((adf+visdown)
     &        /stephanBoltzmannConstant
     &        /emissivity(i,j,bi,bj))**(1./4.)
       ELSEIF(radcode .EQ. 11) THEN
C slab model
        surfaceT(i,j,bi,bj)=surfaceT(i,j,bi,bj)+
     &     deltaT/cpOfSs/ssThickness
     &    /rhoOfSs*(visdown
     &        +adf-emissivity(i,j,bi,bj)*
     &     stephanBoltzmannConstant*
     &     surfaceT(i,j,bi,bj)**4.)

       ENDIF

       tfrostsamz=1./(1./referenceTempInFrostEqn - atm_Rd/atm_L*
     &            LOG(Ro_Surf(i,j,bi,bj)/referencePresInFrostEqn))

       IF(surfaceT(i,j,bi,bj) .LE. tfrostsamz) THEN
        surfaceT(i,j,bi,bj)=tfrostsamz
       ENDIF
      ENDIF

C Heating rate calculated directly for Mars or Pluto, time-varying
C surface temperature (nonzero thermal inertia)
      IF(radcode .EQ. 4 .OR. radcode .EQ. 6 
     &                     .OR .radcode .EQ. 10) THEN
       nn=nn+1

       IF (lsIsConstant .EQ. 0) THEN
        CALL GET_LS(myTime,lssol)
       ELSE
	lssol = lsperpetual
       ENDIF

       CALL GET_QO(i,j,bi,bj,lssol,myTime,solarconstant,
     &                            albedoany(i,j,bi,bj),qo)
       visdown=qo
C We assume the downward longwave emission for Mars is significant,
C but not for Pluto
       IF(radcode .EQ. 4) THEN
        DO k=1,Nr
         taucabC(k)=tauinfcab*(1.-rC(k)/Ro_surf(i,j,bi,bj))
         taucabF(k)=tauinfcab*(1.-rF(k)/Ro_surf(i,j,bi,bj))
        ENDDO
        taucabF(Nr+1)=tauinfcab*(1.-rF(Nr+1)/Ro_surf(i,j,bi,bj))
        DO k=1,Nr
         deltau(k)=taucabF(k+1)-taucabF(k)
        ENDDO
        DO k=1,Nr
         IF(selectAddFluid.EQ.0) THEN
          termP(k)=rC(k)
         ELSE
	  termP(k)=rC(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &           /Ro_surf(i,j,bi,bj)
         ENDIF
        ENDDO
       ENDIF
C     atmospheric downward IR flux

       ks=kSurfC(i,j,bi,bj)
       adf=0.
       IF(radcode .EQ. 4) THEN
        DO k=ks,Nr
         adf=adf+stephanBoltzmannConstant*
     &            theta(i,j,k,bi,bj)**4.*
     &          ((termP(k)/atm_po)**atm_kappa)**4.*
     &          EXP(-(taucabC(k)-0.)/betacab)
     &           *deltau(k)
        ENDDO
       ENDIF

       IF(selectAddFluid .EQ. 0) THEN
        surfaceT(i,j,bi,bj)=surfaceT(i,j,bi,bj)+
     &       deltaT/cpOfSs/ssThickness
     &       /rhoOfSs*(visdown+adf
     &       -emissivity(i,j,bi,bj)*stephanBoltzmannConstant*
     &       surfaceT(i,j,bi,bj)**4.)
       ENDIF
#ifdef ALLOW_ADDFLUID
       IF(selectAddFluid .EQ. 1) THEN
        tfrostsurf=1./(1./referenceTempInFrostEqn -
     &        atm_Rd/atm_L*LOG((Ro_surf(i,j,bi,bj)+
     &        etaN(i,j,bi,bj))/referencePresInFrostEqn))
C If there is no frost on the ground, step surfaceT like above
        IF(massonground(i,j,bi,bj) .LE. 0) THEN
         surfaceT(i,j,bi,bj)=surfaceT(i,j,bi,bj)+
     &        deltaT/cpOfSs/ssThickness
     &        /rhoOfSs*(visdown+adf
     &        -emissivity(i,j,bi,bj)*stephanBoltzmannConstant*
     &        surfaceT(i,j,bi,bj)**4.)
         albedoany(i,j,bi,bj)=albedoSs
         emissivity(i,j,bi,bj)=emissivityOfSs
C Check to see if the surface temperature is below freezing.
C If yes, take mass out of atmosphere to bring surface temperature
C up to freezing.  Change albedo, emissivity to that of frost.
C Otherwise, nothing else needs to be done
         IF(surfaceT(i,j,bi,bj).LT.tfrostsurf) THEN
          addMass(i,j,ks,bi,bj)=
     &            rA(i,j,bi,bj)/atm_L*drF(ks)/gravity*
     &            hfacC(i,j,ks,bi,bj)*atm_cp/deltaT*
     &            (surfaceT(i,j,bi,bj)-tfrostsurf)
          massonground(i,j,bi,bj)=massonground(i,j,bi,bj)
     &                              -addMass(i,j,ks,bi,bj)
          albedoany(i,j,bi,bj)=albedoFrost
          emissivity(i,j,bi,bj)=emissivityOfFrost
          surfaceT(i,j,bi,bj)=tfrostsurf
         ENDIF
C There is at least some frost on the ground (and therefore the
C surface temperature, albedo, and emissivity are that of frost).
        ELSE
         deltaMdot=(visdown+adf-emissivity(i,j,bi,bj)
     &                 *stephanBoltzmannConstant*tfrostsurf**4.)
     &                 *rA(i,j,bi,bj)/atm_L
C deltaMdot (kg/s) >= 0 : sublimation
C deltaMdot (kg/s) < 0  : deposition
C 
C Deposition
         IF(deltaMdot .LT. 0.) THEN
          addMass(i,j,ks,bi,bj)=deltaMdot
          massonground(i,j,bi,bj)=massonground(i,j,bi,bj)-
     &          deltaMdot
          albedoany(i,j,bi,bj)=albedoFrost
          emissivity(i,j,bi,bj)=emissivityOfFrost
          surfaceT(i,j,bi,bj)=tfrostsurf
         ELSE
C Check if there is the required frost on the ground to sublimate
C If yes put it in the atmosphere and remove from ground
          IF(deltaMdot .LT. massonground(i,j,bi,bj)) THEN
           addMass(i,j,ks,bi,bj)=deltaMdot
           massonground(i,j,bi,bj)=massonground(i,j,bi,bj)-
     &          deltaMdot
           albedoany(i,j,bi,bj)=albedoFrost
           emissivity(i,j,bi,bj)=emissivityOfFrost
           surfaceT(i,j,bi,bj)=tfrostsurf
C If not, put what's left in the atmosphere and warm the surface
C with the remaining energy
          ELSE
           addMass(i,j,ks,bi,bj)=massonground(i,j,bi,bj)
           surfaceT(i,j,bi,bj)=tfrostsurf+
     &       deltaT*atm_L*massonground(i,j,bi,bj)/
     &       rhoOfSs/cpOfSs/rA(i,j,bi,bj)/ssThickness
           massonground(i,j,bi,bj)=0.
           albedoany(i,j,bi,bj)=albedoSs
           emissivity(i,j,bi,bj)=emissivityOfSs            
C            write(*,*) 'amz flag1'
C            write(*,*) 'mog = ',massonground(i,j,bi,bj)
C            write(*,*) 'albedo = ', albedoany(i,j,bi,bj)
C            write(*,*) 'ts = ', surfaceT(i,j,bi,bj)
C            write(*,*) 'tfrost = ',tfrostsurf
C            STOP
          ENDIF
         ENDIF 
        ENDIF
       ENDIF
#endif
      ENDIF


C NASA Ames Mars RT
      IF(radcode .EQ. 9) THEN
       write(*,*) 'oops! radcode 9 not functional yet. This error was
     &       encountered in DO_SURFACE_CALC.F'
       STOP
      ENDIF


C Checks
      IF(nn .GT. 1) THEN
       write(*,*) 'Error: more than one surface 
     &      temperature scheme was called'
       STOP
      ENDIF

      IF(nn .EQ. 0) THEN
       write(*,*) 'Error: no surface calculation was done in 
     &                     DO_SURFACE_CALC.F'

       STOP
      ENDIF

      RETURN
      END
