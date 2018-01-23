C $Header: /u/gcmpack/MITgcm/verification/hs94.cs-32x32x5/code/apply_forcing.F,v 1.1 2014/08/20 20:23:10 jmc Exp $
C $Name: checkpoint65e $

#include "CPP_OPTIONS.h"

C--  File apply_forcing.F:
C--   Contents
C--   o APPLY_FORCING_U
C--   o APPLY_FORCING_V
C--   o APPLY_FORCING_T
C--   o APPLY_FORCING_S

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: APPLY_FORCING_U
C     !INTERFACE:
      SUBROUTINE APPLY_FORCING_U(
     U                     gU_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R APPLY_FORCING_U
C     | o Contains problem specific forcing for zonal velocity.
C     *==========================================================*
C     | Adds terms to gU for forcing by external sources
C     | e.g. wind stress, bottom friction etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gU_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gU_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
      INTEGER i, j, k1

C     ===amz vars for Rayleigh friction===
      _RL rayleighuamz(Nr)
C     ---a data table of the damping rate, in inverse seconds
C     ---length(rayleighamz)=Nr; rayleighamz(1) <=> rC(1)
C    ====end amz vars
CEOP
      _RL recip_P0g, termP, rFullDepth
      _RL kV, kF, sigma_b

C--   Forcing terms
      DO k1=1,Nr
       rayleighuamz(k1)=0.
      ENDDO
C Old Pluto/Triton (no Newtonian relaxation spinup)
C      IF(radcode .EQ. 6) THEN
C        rayleighuamz(Nr)=0.02083333
C        rayleighuamz(Nr-1)=0.006944444
C        rayleighuamz(Nr-2)=0.00231481
C        rayleighuamz(Nr-3)=0.00071603
C        rayleighuamz(Nr-4)=0.00023868
C Titan
      IF(radcode .EQ. 11 .OR. radcode .EQ. 12) THEN
       rayleighuamz(Nr)=7.71605E-7
       rayleighuamz(Nr-1)=2.5720E-7
       rayleighuamz(Nr-2)=8.5734E-8
C Pluto/Triton
      ELSE IF (radcode .EQ. 10 .OR. radcode .EQ. 13
     & .OR. radcode .EQ. 6) THEN
       rayleighuamz(Nr)=0.001013799
       rayleighuamz(Nr-1)=0.000337933
       rayleighuamz(Nr-2)=0.000112644
C Mars
      ELSEIF(radcode .EQ. 2 .OR. radcode .EQ. 4
     &  .OR. radcode .EQ. 8 .OR. radcode .Eq. 9) THEN
       rayleighuamz(Nr)=0.0001013799
       rayleighuamz(Nr-1)=0.0000337933
       rayleighuamz(Nr-2)=0.0000112644
      ELSE
       write(*,*) 'Error: radcode not supported'
       STOP
      ENDIF
C

C No mass correction here, since it's a sigma coordinate
      kF = 1. _d 0/86400. _d 0
      sigma_b = 0.7 _d 0
      rFullDepth = rF(1)-rF(Nr+1)
c     DO j=1,sNy
C-jmc: Without CD-scheme, this is OK ; but with CD-scheme, needs to cover [0:sNy+1]
      DO j=0,sNy+1
       DO i=1,sNx+1
        IF ( maskW(i,j,k,bi,bj).EQ.oneRS ) THEN
         IF ( selectSigmaCoord.EQ.0 ) THEN
          IF (selectAddFluid .EQ. 1) THEN
           recip_P0g=MAX(oneRL/(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj)),
     &           oneRL/(Ro_surf(i-1,j,bi,bj)+etaN(i-1,j,bi,bj)))
          termP = 0.5 *(MIN(rF(k)*(Ro_surf(i,j,bi,bj)
     &          +etaN(i,j,bi,bj))/
     &        Ro_surf(i,j,bi,bj)*recip_P0g, oneRL)+rF(k+1)
     &            *(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))/
     &               Ro_surf(i,j,bi,bj)
     &            *recip_P0g )
          ELSE
           recip_P0g = 
     &         MAX(recip_Rcol(i,j,bi,bj),recip_Rcol(i-1,j,bi,bj))
           termP = 0.5 _d 0*( MIN( rF(k)*recip_P0g, oneRL )
     &                      +rF(k+1)*recip_P0g )
c          termP = 0.5 _d 0*( rF(k) + rF(k+1) )*recip_P0g
          ENDIF
         ELSE
C-- Pressure at U.point :
c         midP = rLowW(i,j,bi,bj) + aHybSigmC(k)*rFullDepth
c    &         + bHybSigmC(k)
c    &          *(etaHw(i,j,bi,bj)+rSurfW(i,j,bi,bj)-rLowW(i,j,bi,bj))
C-- Sigma at U.point :
c         termP = ( midP - rLowW(i,j,bi,bj))
c    &          /(etaHw(i,j,bi,bj)+rSurfW(i,j,bi,bj)-rLowW(i,j,bi,bj))
C-  which simplifies to:
          termP = aHybSigmC(k)*rFullDepth
#ifdef NONLIN_FRSURF
     &          /(etaHw(i,j,bi,bj)+rSurfW(i,j,bi,bj)-rLowW(i,j,bi,bj))
#else
     &          /(rSurfW(i,j,bi,bj)-rLowW(i,j,bi,bj))
#endif
     &          + bHybSigmC(k)
         ENDIF
         kV = kF*MAX( zeroRL, (termP-sigma_b)/(1. _d 0-sigma_b) )
     &             +rayleighuamz(k)
         gU_arr(i,j) = gU_arr(i,j)
     &               - kV*uVel(i,j,k,bi,bj)
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: APPLY_FORCING_V
C     !INTERFACE:
      SUBROUTINE APPLY_FORCING_V(
     U                     gV_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R APPLY_FORCING_V
C     | o Contains problem specific forcing for merid velocity.
C     *==========================================================*
C     | Adds terms to gV for forcing by external sources
C     | e.g. wind stress, bottom friction etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gV_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gV_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
      INTEGER i, j, k1

C     ===amz vars for Rayleigh friction===
      _RL rayleighvamz(Nr)
C     ---a data table of the damping rate, in inverse seconds
C     ---length(rayleighamz)=Nr; rayleighamz(1) <=> rC(1)
C    ====end amz vars

CEOP
      _RL recip_P0g, termP, rFullDepth
      _RL kV, kF, sigma_b

C--   Forcing terms
      DO k1=1,Nr
       rayleighvamz(k1)=0.
      ENDDO
C Old Pluto/Triton (no Newtonian relaxation spinup)
C      IF(radcode .EQ. 6) THEN
C        rayleighvamz(Nr)=0.02083333
C        rayleighvamz(Nr-1)=0.006944444
C        rayleighvamz(Nr-2)=0.00231481
C        rayleighvamz(Nr-3)=0.00071603
C        rayleighvamz(Nr-4)=0.00023868
C Titan
      IF(radcode .EQ. 11 .OR. radcode .EQ. 12) THEN
       rayleighvamz(Nr)=7.71605E-7
       rayleighvamz(Nr-1)=2.5720E-7
       rayleighvamz(Nr-2)=8.5734E-8
      ELSE IF(radcode .EQ. 10 .OR. radcode .EQ. 13
     &   .OR. radcode .EQ. 6) THEN
       rayleighvamz(Nr)=0.001013799
       rayleighvamz(Nr-1)=0.000337933
       rayleighvamz(Nr-2)=0.000112644
C Mars
      ELSEIF(radcode .EQ. 2 .OR. radcode .EQ. 4
     & 	.OR. radcode .EQ. 8 .OR. radcode .Eq. 9) THEN
       rayleighvamz(Nr)=0.0001013799
       rayleighvamz(Nr-1)=0.0000337933
       rayleighvamz(Nr-2)=0.0000112644
      ELSE
       write(*,*) 'Error: radcode not supported'
       STOP
      ENDIF
C
      kF = 1. _d 0/86400. _d 0
      sigma_b = 0.7 _d 0
      rFullDepth = rF(1)-rF(Nr+1)
      DO j=1,sNy+1
c      DO i=1,sNx
C-jmc: Without CD-scheme, this is OK ; but with CD-scheme, needs to cover [0:sNx+1]
       DO i=0,sNx+1
        IF ( maskS(i,j,k,bi,bj).EQ.oneRS ) THEN
         IF ( selectSigmaCoord.EQ.0 ) THEN
          IF(selectAddFluid .EQ. 1) THEN
           recip_P0g=MAX(oneRL/(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj)),
     &             oneRL/(Ro_surf(i,j-1,bi,bj)+etaN(i,j-1,bi,bj)))
           termP = 0.5 *(MIN(rF(k)*(Ro_surf(i,j,bi,bj)+
     &                   etaN(i,j,bi,bj))/
     &        Ro_surf(i,j,bi,bj)*recip_P0g,oneRL)+rF(k)
     &            *(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))/
     &                      Ro_surf(i,j,bi,bj)
     &            *recip_P0g )
          ELSE
           recip_P0g = 
     &       MAX(recip_Rcol(i,j,bi,bj),recip_Rcol(i,j-1,bi,bj))
           termP = 0.5 _d 0*( MIN( rF(k)*recip_P0g, oneRL )
     &                      +rF(k+1)*recip_P0g )
c          termP = 0.5 _d 0*( rF(k) + rF(k+1) )*recip_P0g
          ENDIF
         ELSE
C-- Pressure at V.point :
c         midP = rLowS(i,j,bi,bj) + aHybSigmC(k)*rFullDepth
c    &         + bHybSigmC(k)
c    &          *(etaHs(i,j,bi,bj)+rSurfS(i,j,bi,bj)-rLowS(i,j,bi,bj))
C-- Sigma at V.point :
c         termP = ( midP - rLowS(i,j,bi,bj))
c    &          /(etaHs(i,j,bi,bj)+rSurfS(i,j,bi,bj)-rLowS(i,j,bi,bj))
C-  which simplifies to:
          termP = aHybSigmC(k)*rFullDepth
#ifdef NONLIN_FRSURF
     &          /(etaHs(i,j,bi,bj)+rSurfS(i,j,bi,bj)-rLowS(i,j,bi,bj))
#else
     &          /(rSurfS(i,j,bi,bj)-rLowS(i,j,bi,bj))
#endif
     &          + bHybSigmC(k)
         ENDIF
         kV = kF*MAX( zeroRL, (termP-sigma_b)/(1. _d 0-sigma_b) )
     &             +rayleighvamz(k)
         gV_arr(i,j) = gV_arr(i,j)
     &               - kV*vVel(i,j,k,bi,bj)
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: APPLY_FORCING_T
C     !INTERFACE:
      SUBROUTINE APPLY_FORCING_T(
     U                     gT_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R APPLY_FORCING_T
C     | o Contains problem specific forcing for temperature.
C     *==========================================================*
C     | Adds terms to gT for forcing by external sources
C     | e.g. heat flux, climatalogical relaxation, etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
      INTEGER i, j
CEOP

C Loop to determine which radiation scheme to call.
C Gray scheme, Netwonian relaxation.  Mars.
      IF (radcode .EQ. 2 .OR. radcode .EQ. 8) THEN
       CALL NRAD(bi,bj,k,myTime,gT_arr)
C Heating rate supplied directly, Mars case: Caballero
C semi-gray scheme for CO2, optional dust heating and cooling,
C optional convective adjustment.
      ELSE IF (radcode .EQ. 4) THEN
       CALL MARSRAD(bi,bj,k,myTime,myThid,gt_arr)
C Heating rate supplied directly, Pluto case: Yelle and
C Lunine (1989) heating at 3.3 um, cooling at 7.6 um.
      ELSE IF (radcode .EQ. 6) THEN
       CALL PLUTORAD(bi,bj,k,myTime,gT_arr)
CC      ELSE IF (radcode .EQ. 9) THEN
C      write(*,*) 'AMZ external_forcing.F begin Mars RT'
C      flush(6)
CC       CALL MARS_RT_DRIVER(bi,bj,k,myTime,gT_arr)
C      write(*,*) 'AMZ external_forcing.F end Mars RT'
C      flush(6)
      ELSEIF(radcode .Eq. 11) THEN
       CALL TITANRAD(bi,bj,k,myTime,gT_arr)
      ELSEIF(radcode .EQ. 12) THEN
       CALL TITANNRAD(bi,bj,k,myTime,gT_arr)
      ELSEIF(radcode .EQ. 13) THEN
       CALL PLUTONRAD(bi,bj,k,gT_arr)
C      ELSEIF(radcode .EQ. 10) THEN
C       CALL STROBELRAD(bi,bj,k,myTime,gT_arr)

      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: APPLY_FORCING_S
C     !INTERFACE:
      SUBROUTINE APPLY_FORCING_S(
     U                     gS_arr,
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R APPLY_FORCING_S
C     | o Contains problem specific forcing for merid velocity.
C     *==========================================================*
C     | Adds terms to gS for forcing by external sources
C     | e.g. fresh-water flux, climatalogical relaxation, etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "SURFACE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
c     INTEGER i, j
CEOP

C--   Forcing term

      RETURN
      END
