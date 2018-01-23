#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_ATM_IR_FLUXES
C     !INTERFACE:

      SUBROUTINE GET_ATM_IR_FLUXES(i,j,bi,bj,iPlus,iMinus)

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_ATM_IR_FLUXES
C     | o Calculate IR fluxes, iPlus and iMinus 
C     |   using Caballero et al. 2008
C     |   semigray model for CO2.
C     |  
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
      _RL taucabC(Nr),taucabF(Nr+1)
      _RL termP(Nr),iMinus1,iPlus1,deltau(Nr)
      _RL iPlus(Nr)
      _RL iMinus(Nr)



C     Begin Equations

      ks=kSurfC(i,j,bi,bj)

C     Calculate optical depth, according to the Caballero
C     et al. 2008 convention where the ground is 0 optical depth
C     c.f. amz convention in get_teq 
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



C Calculate Fluxes
      DO k=ks,Nr
       iMinus1=stephanBoltzmannConstant*theta(i,j,k,bi,bj)**4.*
     &            ((termP(k)/atm_po)**atm_kappa)**4.*
     &         EXP(-(taucabF(k+1)-taucabC(k))/betacab)
     &        *(taucabF(k+1)-taucabC(k))
       IF(k .LT. Nr) THEN
        DO ki=k+1,Nr
         iMinus1=iMinus1+stephanBoltzmannConstant*
     &            theta(i,j,ki,bi,bj)**4.*
     &          ((termP(ki)/atm_po)**atm_kappa)**4.*
     &          EXP(-(taucabC(ki)-taucabC(k))/betacab)
     &           *deltau(ki)
        ENDDO
       ENDIF
       iMinus(k)=iMinus1
       iPlus1=stephanBoltzmannConstant*theta(i,j,k,bi,bj)**4.*
     &            ((termP(k)/atm_po)**atm_kappa)**4.*
     &          (taucabC(k)-taucabF(k))
     &        *EXP((taucabF(k)-taucabC(k))/betacab)
     &        +betacab*stephanBoltzmannConstant
     &        *emissivity(i,j,bi,bj)
     &        *surfaceT(i,j,bi,bj)**4.
     &       *EXP(-taucabC(k)/betacab)
       IF(k .GT. ks) THEN
        DO ki=k-1,ks,-1
         iPlus1=iPlus1+stephanBoltzmannConstant*
     &            theta(i,j,ki,bi,bj)**4.*
     &          ((termP(ki)/atm_po)**atm_kappa)**4.*
     &           EXP((taucabC(ki)-taucabC(k))/betacab)
     &           *deltau(ki)
        ENDDO
       ENDIF
       iPlus(k)=iPlus1
      ENDDO


      RETURN
      END
