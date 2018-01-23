#include "CPP_OPTIONS.h"

      SUBROUTINE PLUTONRAD(bi,bj,kLev,gT_arr)

      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"


      INTEGER i,j,kLev,bi,bj
      _RL kt
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)


      kt=1./(30.*86400.)

      DO j=1,sNy
       DO i=1,sNx
        gT_arr(i,j)=gT_arr(i,j)
     &    -kt*(theta(i,j,kLev,bi,bj)-thetaEqSave(i,j,kLev,bi,bj))
C    &     -10.
       ENDDO
      ENDDO

      RETURN
      END

