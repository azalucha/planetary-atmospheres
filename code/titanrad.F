#include "CPP_OPTIONS.h"

      SUBROUTINE TITANRAD(bi,bj,kLev,myTime,gt_arr)
C     *==========================================================*
C     | S/R TITANRAD
C     | o Titan radiation scheme with heating rate supplied directly.
C     |   Taken from Robinson and Catling (2012).
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
      _RL term1,term2
      _RL termP
      _RL iPlus(Nr)
      _RL iMinus(Nr)
      _RL q(Nr),dTdt
      _RL gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      DO j=1,sNy
       DO i=1,sNx
C  thermodynamics loop counts down from Nr.  Do flux calculations at first
C  instance
        IF(kLev .EQ. Nr) THEN
         CALL DO_SURFACE_CALC(i,j,bi,bj,myTime)
         IF(initrad(i,j,bi,bj) .EQ. 1. .OR.
     &            MOD(myTime,radTimeStep).EQ.0.) THEN
          CALL GET_TITAN_IR_FLUXES(i,j,bi,bj,iPlus,iMinus,myTime)
          CALL GET_FLUX_DIV(i,j,bi,bj,iPlus,iMinus,q)
          initrad(i,j,bi,bj)=0.
          DO k=1,Nr
           qr(i,j,k,bi,bj)=q(k)
          ENDDO
         ENDIF
        ENDIF
        IF(selectAddFluid.EQ.0) THEN
         termP=rC(kLev)
        ELSE
         termP=rC(kLev)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &           /Ro_surf(i,j,bi,bj)
        ENDIF
        dTdt=qr(i,j,kLev,bi,bj)*(termP/atm_po)**(-atm_kappa)
        gT_arr(i,j)=gT_arr(i,j)+dTdt*maskC(i,j,kLev,bi,bj)
       ENDDO
      ENDDO

      RETURN
      END
