#include "CPP_OPTIONS.h"
      subroutine ini_optcld(QEXTVc,QSCATVc,GrefVc,QEXTIc,QSCATIc,GrefIc,
     &                  Qxvc,Qxic,Qsvc,Qsic,gvc,gic,Qextrefc)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C  Initialize cloud optical constants, both IR and visible

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GRID_H.h"
#include "RADINC_H.h"
#include "CLDCOMMON_H.h"

!  Arguments
!  ---------

      real*8  QEXTVc(L_NSPECTV)
      real*8  QSCATVc(L_NSPECTV)
      real*8  GrefVc(L_NSPECTV)

      real*8  QEXTIc(L_NSPECTI)
      real*8  QSCATIc(L_NSPECTI)
      real*8  GrefIc(L_NSPECTI)

      real*8  Qxvc(L_LEVELS+1,L_NSPECTV)
      real*8  Qsvc(L_LEVELS+1,L_NSPECTV)
      real*8  gvc(L_LEVELS+1,L_NSPECTV)

      real*8  Qxic(L_LEVELS+1,L_NSPECTI)
      real*8  Qsic(L_LEVELS+1,L_NSPECTI)
      real*8  gic(L_LEVELS+1,L_NSPECTI)

      real*8  Qextrefc(L_LEVELS+1)

C  Local variables
C  ---------------

      integer i,k

C Initialyze various variables
C ----------------------------

      DO K = 1, L_LEVELS+1
        Qextrefc(K) = QEXTVc(L_NREFV)
      ENDDO

      do i = 1, nlonv
        DO K = 1, L_LEVELS+1
          Qxvc(K,i) = QEXTVc(i)
          Qsvc(K,i) = QSCATVc(i)
          gvc(K,i)  = GrefVc(i)
        ENDDO
      enddo
      do i = 1, nloni
        DO K = 1, L_LEVELS+1
          Qxic(K,i) = QEXTIc(i)
          Qsic(K,i) = QSCATIc(i)
          gic(K,i)  = GrefIc(i)
        ENDDO
      enddo

      return

      end
