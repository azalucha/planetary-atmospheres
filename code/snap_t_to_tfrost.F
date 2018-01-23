#include "CPP_OPTIONS.h"

C     !ROUTINE: SNAP_T_TO_TFROST
C     !INTERFACE:
      SUBROUTINE SNAP_T_TO_TFROST(myThid)     

C     !DESCRIPTION:                 
C     *==========================================================*
C     | S/R SNAP_T_TO_TFROST
C     | o snaps theta to frost temperature if it goes below the
C     |   frost temperature.  For the case of selectAddFluid=0
C     *==========================================================*
C     

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
C     == Routine arguments ==   
C     myThid    :: Thread Id number

      INTEGER myThid                

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j,k,bi,bj      :: Loop counters
C     tfrostamz        :: frost (potential) temperature 

      _RL tfrostamz,tfrostsamz
      INTEGER i, j,k,bi,bj

      DO bj = myByLo(myThid), myByHi(myThid)           
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO J=1,sNy
         DO I=1,sNx
          DO K=1,Nr
           tfrostamz=1./(1./referenceTempInFrostEqn - atm_Rd/atm_L*
     &            LOG(rC(k)/referencePresInFrostEqn))
     &          / ((rC(k)/atm_po)**atm_kappa)
           IF (theta(i,j,k,bi,bj) .LT. tfrostamz) THEN
            dtsnapdt(i,j,k,bi,bj) = (theta(i,j,k,bi,bj)-tfrostamz)
     &                               /deltaT
            theta(i,j,k,bi,bj) = tfrostamz
           ELSE
            dtsnapdt(i,j,k,bi,bj)=0.
           ENDIF
          ENDDO
          tfrostsamz=1./(1./referenceTempInFrostEqn - atm_Rd/atm_L*
     &            LOG(Ro_Surf(i,j,bi,bj)/referencePresInFrostEqn))
          IF(surfaceT(i,j,bi,bj) .LT. tfrostsamz) THEN
           surfaceT(i,j,bi,bj)=tfrostsamz
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO


      RETURN
      END
