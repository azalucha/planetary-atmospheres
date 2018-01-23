#include "CPP_OPTIONS.h"
      SUBROUTINE cldprofile(psf,ptrop,nlev,sigma,pcld,tautotcld,
     &                      taurefcld)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C
C  PURPOSE
C     CLDPROFILE fills TAUREFCLD(K) with the cloud optical depth (at 
C     the reference wavelength) for each layer.
C
C  PSF       - Surface pressure (mbar)
C  PTROP     - Tropopause pressure (mbar)
C  NLEV      - Number of levels
C  SIGMA     - Sigma values at each level
C  PCLD      - Pressure where cloud is located
C  TAUTOTCLD - Total cloud optical depth of the column
C  TAUREFCLD - Cloud optical depth, in each layer, at the reference
C              wavelength.
C  WTC       - Weight of cloud center
C  WTW       - Weight of cloud wing
C
C----------------------------------------------------------------------

      implicit none
#include "EEPARAMS.h"

      integer  nlev
      _RL  psf, ptrop, pcld, sigma(nlev), pl(nlev)
      _RL   TAUREFCLD(nlev+1), tautotcld
      _RL   wtc,wtw
      PARAMETER(wtc=0.5,wtw=0.25)

C     implicit none

      integer K
 
C#=====================================================================

C  Pressures at layer boundaries and midpoints.

      pl(1) = 0.0D0
      pl(2) = ptrop/2.0D0
      pl(3) = ptrop

      do K=3,nlev
        pl(k) = sigma(k)*(psf-ptrop)+ptrop
      end do

C  Zero out all layers

      do K=1,nlev+1
        taurefcld(k) = 0.0D0
      end do

C  Find layer where cloud is located

      if(pcld.le.pl(5)) then
C       put cloud center in the top layer below the tropopause
        taurefcld(4) = (wtc+wtw)*tautotcld/2.0D0
        taurefcld(5) = taurefcld(4)
        taurefcld(6) = wtw*tautotcld/2.0D0
        taurefcld(7) = taurefcld(6)

      elseif(pcld.gt.pl(nlev-2)) then

C      put cloud in the bottom layer

        taurefcld(nlev)   = (wtc+wtw)*tautotcld/2.0D0
        taurefcld(nlev-1) = taurefcld(nlev)
        taurefcld(nlev-2) = wtw*tautotcld/2.0D0
        taurefcld(nlev-3) = taurefcld(nlev-2)

      else
       
        do k=5,nlev-4
          
          if(pcld.gt.pl(k) .and. pcld.le.pl(k+2)) then

C  center layer

            taurefcld(k+1) = wtc*tautotcld/2.0D0
            taurefcld(k+2) = taurefcld(k+1)

C  layers adjacent to the cloud's central layer (wings)

C  layer above
            taurefcld(k)   = wtw*tautotcld/2.0
            taurefcld(k-1) = taurefcld(k)
C  layer below
            taurefcld(k+3) = taurefcld(k)
            taurefcld(k+4) = taurefcld(k)
            exit
          end if
        end do
      end if

      return
      end
