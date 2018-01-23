      subroutine radsetup

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!     PURPOSE:
!        Bundle the new radiation code setup subroutines and call
!     this one subroutine from main, where the three included files
!     are also listed.  Quantities are passed between this driver
!     and the radiation code via modules (eg radcommon_h).
!
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h
      use radcommon_h

      implicit none

!======================================================================C

      call setspv(WNOV,DWNV,WAVEV,SOLARF,TAURAY)
      call setspi(WNOI,DWNI,WAVEI)
      call setrad(TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV,       &
                        QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,     &
                        QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV)

      return
      end subroutine radsetup
