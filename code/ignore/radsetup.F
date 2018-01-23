#include "CPP_OPTIONS.h"
      subroutine radsetup(DWNV,WAVEV,SOLARF,TAURAY,
     & WNOI,DWNI,WAVEI,
     & TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV0,
     & QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,
     & QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV,PLANCKIR)

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

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GRID_H.h"
#include "RADINC_H.h"
#include "RADCOMMON_H.h"

C======================================================================C

      call setspv(WNOV,DWNV,WAVEV,SOLARF,TAURAY)
      call setspi(WNOI,DWNI,WAVEI,PLANCKIR)
      call setrad(TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV0,       
     &            QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,     
     &            QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV)

      return
      end
