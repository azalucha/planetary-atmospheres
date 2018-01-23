#include "CPP_OPTIONS.h"
      subroutine ini_optdst(QEXTV,QSCATV,GrefV,QEXTI,QSCATI,GrefI,
     &                  Qxv,Qxi,Qsv,Qsi,gv,gi,Qextref)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!  Initialize dust optical constants, both IR and visible

      implicit none
#include "EEPARAMS.h"
#include "RT_GRID.h"
#include "RADINC.h"
#include "CLDCOMMON.h"

!  Arguments
!  ---------

      _RL  QEXTV(L_NSPECTV)
      _RL  QSCATV(L_NSPECTV)
      _RL  GrefV(L_NSPECTV)

      _RL  QEXTI(L_NSPECTI)
      _RL  QSCATI(L_NSPECTI)
      _RL  GrefI(L_NSPECTI)

      _RL  Qxv(L_LEVELS+1,L_NSPECTV)
      _RL  Qsv(L_LEVELS+1,L_NSPECTV)
      _RL  gv(L_LEVELS+1,L_NSPECTV)

      _RL  Qxi(L_LEVELS+1,L_NSPECTI)
      _RL  Qsi(L_LEVELS+1,L_NSPECTI)
      _RL  gi(L_LEVELS+1,L_NSPECTI)

      _RL  Qextref(L_LEVELS+1)

!  Local variables
!  ---------------

      integer i,k

! Initialyze various variables
! ----------------------------

      DO K = 1, L_LEVELS+1
        Qextref(K) = QEXTV(L_NREFV)
      ENDDO

      do i = 1, nlonv
        DO K = 1, L_LEVELS+1
          Qxv(K,i) = QEXTV(i)
          Qsv(K,i) = QSCATV(i)
          gv(K,i)  = GrefV(i)
        ENDDO
      enddo
      do i = 1, nloni
        DO K = 1, L_LEVELS+1
          Qxi(K,i) = QEXTI(i)
          Qsi(K,i) = QSCATI(i)
          gi(K,i)  = GrefI(i)
        ENDDO
      enddo

      return

      end
