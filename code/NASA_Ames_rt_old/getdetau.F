#include "CPP_OPTIONS.h"
      SUBROUTINE GETDETAU(DTDEL,TDEL,TAUCUMIN,WDEL,CDEL,UBAR0,DETAU)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C======================================================================!

      implicit none
#include "EEPARAMS.h"
#include "RT_GRID.h"
#include "RADINC.h"

      _RL   W0(L_NLAYRAD), COSBAR(L_NLAYRAD), DTAU(L_NLAYRAD)
      _RL TAU(L_NLEVRAD), WDEL(L_NLAYRAD), CDEL(L_NLAYRAD)
      _RL DTDEL(L_NLAYRAD), TDEL(L_NLEVRAD)
      _RL FACTOR, TAUCUMIN(L_LEVELS), TAUCUM(L_LEVELS)
      _RL detau

      integer  NAYER, L, K
      _RL   ubar0

C======================================================================C
 
      NAYER  = L_NLAYRAD
      
C  Delta-Eddington Scaling

      FACTOR    = 1.0D0 - WDEL(1)*CDEL(1)**2

      TAU(1)    = TDEL(1)*FACTOR
      TAUCUM(1) = 0.0D0
      TAUCUM(2) = TAUCUMIN(2)*FACTOR
      TAUCUM(3) = TAUCUM(2) +(TAUCUMIN(3)-TAUCUMIN(2))*FACTOR

      DO L=1,L_NLAYRAD-1
        FACTOR      = 1.0D0 - WDEL(L)*CDEL(L)**2
        W0(L)       = WDEL(L)*(1.0D0-CDEL(L)**2)/FACTOR
        COSBAR(L)   = CDEL(L)/(1.0D0+CDEL(L))
        DTAU(L)     = DTDEL(L)*FACTOR
        TAU(L+1)    = TAU(L)+DTAU(L)
        K           = 2*(L+1)
        TAUCUM(K)   = TAU(L+1)
        TAUCUM(K+1) = TAUCUM(K) + (TAUCUMIN(K+1)-TAUCUMIN(K))*FACTOR
      ENDDO

C  Bottom layer

      L             = L_NLAYRAD
      FACTOR        = 1.0D0 - WDEL(L)*CDEL(L)**2
      W0(L)         = WDEL(L)*(1.0D0-CDEL(L)**2)/FACTOR
      COSBAR(L)     = CDEL(L)/(1.0D0+CDEL(L))
      DTAU(L)       = DTDEL(L)*FACTOR
      TAU(L+1)      = TAU(L)+DTAU(L)
      TAUCUM(2*L+1) = TAU(L+1)
      detau         = TAUCUM(2*L+1)

      return
C      end subroutine getdetau
      END
