#include "CPP_OPTIONS.h"

      SUBROUTINE NRAD(bi,bj,kLev,myTime,gt_arr)
C     *==========================================================*
C     | S/R NRAD
C     | o Computes radiation for Newtonian relaxation to a gray
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
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      _RL thetaLim,kT,ka,ks,sigma_b
      _RL thetaEq(sNx,sNy),termP(sNx,sNy)

      _RL tfrostamz


C--   Forcing term(s)
C      ka=1./hstaua
C      ks=1./hstaus
C      sigma_b = 0.7 _d 0

       DO j=1,sNy
        DO i=1,sNx
         IF(kLev .EQ. Nr) THEN
          CALL DO_SURFACE_CALC(i,j,bi,bj,myTime)
         ENDIF
         IF(selectAddFluid.EQ.0) THEN
          termP(i,j)=rC(kLev)
         ELSE
          termP(i,j)=rC(kLev)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &           /Ro_surf(i,j,bi,bj)
         ENDIF
        ENDDO
       ENDDO

C      -------Get equilibrium temperature, teq
       IF(radcode .EQ. 2) THEN
        CALL GET_THETAEQ(thetaEq,termP,bi,bj,myTime,kLev)
       ELSE IF(radcode .EQ. 8) THEN
        CALL GET_THETAEQ_CAB(thetaEq,termP,bi,bj,myTime,kLev)
       ENDIF
C      -----End specification of teq


C      -----Snap Teq to frost temperature (NOTE: TEQ, not T!)
       DO j=1,sNy
        DO i=1,sNx
         tfrostamz=1./(1./referenceTempInFrostEqn - atm_Rd/atm_L*
     &            LOG(termP(i,j)/referencePresInFrostEqn))
     &          / ((termP(i,j)/atm_po)**atm_kappa)
         IF (thetaEq(i,j) .LT. tfrostamz) THEN
          thetaEq(i,j) = tfrostamz
          hstaua(i,j,bi,bj)=hstauacold
          hstaus(i,j,bi,bj)=hstauscold
         ELSE
          hstaua(i,j,bi,bj)=hstauawarm
          hstaus(i,j,bi,bj)=hstauswarm
         ENDIF

C      ----end teq frost condition

         thetaeqsave(i,j,kLev,bi,bj)=thetaEq(i,j)

         IF(kLev .EQ. Nr) THEN
          CALL DO_SURFACE_CALC(i,j,bi,bj,myTime)
         ENDIF

C--   Forcing term(s)
         ka=1./hstaua(i,j,bi,bj)
         ks=1./hstaus(i,j,bi,bj)
         sigma_b = 0.7 _d 0
         kT=ka+(ks-ka)
     &     *MAX(0. _d 0,
     &       (termP(i,j)*recip_Rcol(i,j,bi,bj)-sigma_b)/
     &                 (1. _d 0-sigma_b) )
     &     *COS((yC(i,j,bi,bj)*deg2rad))**4.
          gT_arr(i,j)=gT_arr(i,j)
     &        - kT*( theta(i,j,kLev,bi,bj)-thetaEq(i,j) )
     &            *maskC(i,j,kLev,bi,bj)
        ENDDO
       ENDDO

       RETURN
       END
