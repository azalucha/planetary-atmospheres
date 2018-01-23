#include "CPP_OPTIONS.h"

      SUBROUTINE MARSRAD(bi,bj,kLev,myTime,myThid,gt_arr)
C     *==========================================================*
C     | S/R MARSRAD
C     | o Mars radiation scheme with heating rate supplied directly.
C     |   CO2 is semigray scheme from Caballero et al. 2008.
C     |   Optional dust heating (Briegleb 1992) and cooling.
C     |   Optional convective adjustment
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
      INTEGER myThid

      INTEGER i, j,k
      _RL term1,term2
      _RL termP(sNx,sNy)
      _RL iPlus(Nr)
      _RL iMinus(Nr)
      _RL qco2(Nr),tausgridw(Nr+1)
      _RL qe(Nr),tausgrid(Nr)
      _RL opt(Nr+1),optC(Nr)

      _RL tfrostamz,dTdt
      INTEGER ktop
      _RL gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      DO j=1,sNy
       DO i=1,sNx
C  thermodynamics loop counts down from Nr.  Do flux calculations at first
C  instance
        IF(kLev .EQ. Nr) THEN
         CALL DO_SURFACE_CALC(i,j,bi,bj,myTime)
         IF(initrad(i,j,bi,bj) .EQ. 1. .OR.
     &            MOD(myTime,radTimeStep).EQ.0.) THEN
          CALL GET_ATM_IR_FLUXES(i,j,bi,bj,iPlus,iMinus)
          CALL GET_FLUX_DIV(i,j,bi,bj,iPlus,iMinus,qco2)
          initrad(i,j,bi,bj)=0.
          DO k=1,Nr
           qr(i,j,k,bi,bj)=qco2(k)
          ENDDO
         ENDIF
        ENDIF
        dTdt=qr(i,j,kLev,bi,bj)*(rC(klev)/atm_po)**-atm_kappa
        gT_arr(i,j)=gT_arr(i,j)+dTdt*maskC(i,j,kLev,bi,bj)
       ENDDO
      ENDDO

      RETURN
      END
