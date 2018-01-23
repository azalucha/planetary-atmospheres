#include "CPP_OPTIONS.h"
      SUBROUTINE FILLTAUCUM

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C
C  PURPOSE
C      FILLTAUCUM fills the TAUCUM array - the dust column density of each
C      GCM level.  NDUSK(K) is the column density of the sublayer whose
C      bottom is level K, and whose top is at K-1.
C
      implicit none
#include "EEPARAMS.h"
#include "RT_GRID.h"
#include "DEFINES.h"
#include "RADINC.h"
#include "RADCOMMON.h"
#include "STANDARD.h"
#include "COMP3CMN.h"

      external JSRCHGT

      INTEGER Nstar, N

C     implicit none

      integer J, I, K, JSRCHGT
      _RL   PSTAR, PSTAR1
 
C#=====================================================================

C  Dust optical depth at the bottom of each sub-layer.

      DO N=1,3
        TAUCUM(N) = 0.0
        TAUREF(N) = 0.0
      ENDDO

      DO K=4,L_LEVELS
        PSTAR     = PSF*SIGMA(K)+PTROP
        PSTAR1    = MAX(PSTAR,PRDST(1))
        NSTAR     = MIN0(JSRCHGT(NPDST-1,PRDST,1,PSTAR1)-1,NPDST-1)
        TAUCUM(K) = TAUDST(NSTAR)+(PSTAR1-PRDST(NSTAR))*
     &              (TAUDST(NSTAR+1) - TAUDST(NSTAR))/
     &              (PRDST(NSTAR+1)-PRDST(NSTAR))
        TAUREF(K) = TAUCUM(K) - TAUCUM(K-1)
      ENDDO

      TAUREF(L_LEVELS+1) = 0.0

      RETURN
      END
