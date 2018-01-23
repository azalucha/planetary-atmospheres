#include "CPP_OPTIONS.h"
      subroutine radsetup

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C     PURPOSE:
C        Bundle the new radiation code setup subroutines and call
C     this one subroutine from main, where the three included files
C     are also listed.  Quantities are passed between this driver
C     and the radiation code via modules (eg radcommon_h).
C
C----------------------------------------------------------------------C

      implicit none
#include "EEPARAMS.h"
#include "RT_GRID.h"
#include "RADINC.h"
#include "RADCOMMON.h"

C======================================================================C

      call setspv(WNOV,DWNV,WAVEV,SOLARF,TAURAY)
      call setspi(WNOI,DWNI,WAVEI,PLANCKIR)
      call setrad(TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV,
     &             QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,     
     &             QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV)

      return
C      end subroutine radsetup
      END
