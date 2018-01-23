#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_FLUX_DIV
C     !INTERFACE:

      SUBROUTINE GET_FLUX_DIV(i,j,bi,bj,fPlus,fMinus,qf)             

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_FLUX_DIV
C     | o flux divergence: g/cp d(I^+-I^-)/dp
C     |  units: K/s
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
C     qr        :: flux divergence
C     fPlus     :: upward flux
C     fMinus    :: downward flux

      INTEGER bi, bj, i,j,k,ks
      _RL fPlus(Nr)
      _RL fMinus(Nr)
      _RL termP(Nr)
      _RL qf(Nr)

C     Begin equations

      DO k=1,Nr
       IF(selectAddFluid.EQ.0) THEN
        termP(k)=rC(k)
       ELSE
        termP(k)=rC(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &     	 /Ro_surf(i,j,bi,bj)
       ENDIF
      ENDDO

      ks=kSurfC(i,j,bi,bj)

      qf(ks)=
     &    ((fPlus(ks+1)-fMinus(ks+1))
     &    -(fPlus(ks)-fMinus(ks)))/
     &    (termP(ks+1)-termP(ks))*gravity/atm_cp

      DO k=ks+1,Nr-1
       qf(k)=
     &    ((fPlus(k+1)-fMinus(k+1))
     & 	  -(fPlus(k-1)-fMinus(k-1)))/
     & 	  (termP(k+1)-termP(k-1))*gravity/atm_cp
      ENDDO

      qf(Nr)=
     &    ((fPlus(Nr)-fMinus(Nr))
     & 	  -(fPlus(Nr-1)-fMinus(Nr-1)))/
     & 	  (termP(Nr)-termP(Nr-1))*gravity/atm_cp

      RETURN
      END
