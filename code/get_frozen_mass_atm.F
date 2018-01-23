#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_FROZEN_MASS_ATM
C     !INTERFACE:

      SUBROUTINE GET_FROZEN_MASS_ATM(myThid)

C     !DESCRIPTION: 
C     *==========================================================*
C     | SUBROUTINE GET_FROZEN_MASS_ATM
C     | o Calculate the amount of mass/time removed from the atmosphere
C     |  when the temperature drops below the frost point.
C     |  units: kg/s
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
#include "SURFACE.h"
#include "FFIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
      INTEGER myThid

C     !LOCAL VARIABLES:
C     Local variables
C     i,j,k,bi,bj       :: Loop counters
      INTEGER i, j, k, bi, bj
      _RL tfrostamz,termP,tfrost,delp,newp
      _RL etaNsum,totalarea,thetaoffset

C     Begin equations

#ifdef ALLOW_ADDFLUID
      IF(selectAddFluid.NE.0) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO k=Nr,1,-1
C         termP=0.5 _d 0*( rF(k) + rF(k+1) )
          DO j = 1,sNy
           DO i = 1,sNx
            termP=rC(k)*(Ro_surf(i,j,bi,bj)+
     &        etaN(i,j,bi,bj))/Ro_surf(i,j,bi,bj)
             tfrostamz=1./(1./referenceTempInFrostEqn - 
     &        atm_Rd/atm_L*LOG(termP/referencePresInFrostEqn))
     &        / ((termP/atm_po)**atm_kappa)
C If the temperature falls below the frost point, take mass away from the atmosphere,
C put the mass on the ground, and change the albedo to that of frost.  Mass can freeze out from
C anywhere in the atmosphere, and instantaneously precipitates to the ground
            IF(tfrostamz.GT.theta(i,j,k,bi,bj) 
     &           .AND. k .GE. kSurfC(i,j,bi,bj)) THEN                    
             addMass(i,j,k,bi,bj)=drF(k)*rA(i,j,bi,bj)*
     &            hFacC(i,j,k,bi,bj)
     &        /gravity/deltaT*(theta(i,j,k,bi,bj)-tfrostamz)
     &           *((termP/atm_po)**atm_kappa)*atm_Cp/atm_L
C Carefull, addMass is negative if it is being taken out of the atmosphere.  Want to add
C positive mass to ground.
             massonground(i,j,bi,bj)=
     &         massonground(i,j,bi,bj)-addMass(i,j,k,bi,bj)
             albedoany(i,j,bi,bj)=albedoFrost
             emissivity(i,j,bi,bj)=emissivityOfFrost
C T needs to stay constant during condensation. Adjust theta accordingly,
C accounting for the fact that the pressure of the box has changed.
             delp=addMass(i,j,k,bi,bj)*gravity*deltaT/rA(i,j,bi,bj)
             newp=termP+delP
             tfrost=1./(1./referenceTempInFrostEqn -
     &        atm_Rd/atm_L*LOG(newp/referencePresInFrostEqn))
             theta(i,j,k,bi,bj)=tfrost*(newp/atm_po)**(-atm_kappa)
             hstaua(i,j,bi,bj)=hstauacold
             hstaus(i,j,bi,bj)=hstauscold
             useSurfaceTflag(i,j,bi,bj)=1.
            ELSE
C If we're below the surface, or above the frost temperature, can't condense.
             addMass(i,j,k,bi,bj)=0.
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
C Here we calculate the global mean eta
C       etaNsum=0.
C       totalarea=0.
C       etaNmax=0.
C       etaNmin=0.
c       DO bj=myByLo(myThid),myByHi(myThid)
c        DO bi=myBxLo(myThid),myBxHi(myThid)
c         DO j = 1, sNy
C          DO i = 1, sNx
c           etaNsum=etaNsum+etaN(i,j,bi,bj)*rA(i,j,bi,bj)
c           totalarea=totalarea+rA(i,j,bi,bj)
c           IF( etaN(i,j,bi,bj) .GT. etaNmax) THEN
c            etaNmax=etaN(i,j,bi,bj)
c           ENDIF
c           IF( etaN(i,j,bi,bj) .LT. etaNmin) THEN
C            etaNmin=etaN(i,j,bi,bj)
c           ENDIF
C          ENDDO
C         ENDDO
C        ENDDO
C       ENDDO
c       etaNmean=etaNsum/totalarea
c       write(*,*) etaNmean,etaNmax,etaNmin
      ENDIF
#endif


      RETURN
      END

