#include "CPP_OPTIONS.h"
      SUBROUTINE SETSPV(WNOV,DWNV,WAVEV,SOLARF,TAURAY)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C     PURPOSE:
C        Set up the spectral intervals in the visible (solar).  Based
C     on Chris McKay's SETSPV code.
C
C     INPUT PARAMETERS
C     L_NSPECTV  - Number of spectral intervals in the visible
C
C     OUTPUT PARAMETERS
C     WNOV       - Array of wavenumbers at the spectral interval
C                  center for the visible.  Array is NSPECTV
C                  elements long.
C     DWNV       - Array of "delta wavenumber", i.e., the width,
C                  in wavenumbers (cm^-1) of each visible spectral
C                  interval.  NSPECTV elements long.
C     WAVEV      - Array (NSPECTV elements long) of the wavelenght
C                  (in microns) at the center of each visible spectral
C                  interval.
C     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
C                  each spectral interval.  Values are for 1 AU, and
C                  are scaled to the Mars distance elsewhere.
C     TAURAY     - Array (NSPECTV elements) of the wavelength independent
C                  part of Rayleigh Scattering.  The pressure dependent 
C                  part is computed elsewhere (OPTCV).
C
C**********************************************************************C

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GRID_H.h"
#include "CONSTANTS_H.h"
#include "RADINC_H.h"


C     BWNV - Bin wavenumber of the edges of the visible spectral bins
C     units are inverse centimeters.  Dimension needs to be changed
C     if the number of visible bins changes.

      REAL*8  WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
      REAL*8  SOLARF(L_NSPECTV), TAURAY(L_NSPECTV)

      REAL*8   SUM, WL
      INTEGER  N, M

C     P0      - Rayleigh scattering reference pressure in pascals.
C     GRAV    - Acceleration due to gravity (g) - MKS

      real*8  :: P0 = 9.423E6

C     Bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
C     spectral bins.  Go from smaller to larger wavenumbers, the same as
C     in the IR.

C     2222.22D0    ->   4.50 microns
C     3087.37D0    ->   3.24 microns
C     4030.63D0    ->   2.48 microns
C     5370.57D0    ->   1.86 microns
C     7651.11D0    ->   1.31 microns
C     12500.00D0   ->   0.80 microns
C     25000.00D0   ->   0.40 microns            
C     41666.67D0   ->   0.24 microns

      REAL*8  BWNV(L_NSPECTV+1) 
      DATA BWNV / 2222.22, 3087.37, 4030.63, 
     &             5370.57, 7651.11, 12500.00, 
     &             25000.00, 41666.67 /

C     Solar flux within each spectral interval, at 1AU (W/M^2)
C     Sum equals 1356 W/m^2 (values from Wehrli, 1985)

      real*8  SOLAR(L_NSPECTV) 
      DATA SOLAR /12.7, 24.2, 54.6, 145.9, 354.9,657.5, 106.3 /

C======================================================================C

C     Set up mean wavenumbers and wavenumber deltas.  Units of 
C     wavenumbers is cm^(-1); units of wavelengths is microns.

      do M=1,L_NSPECTV
        WNOV(M)  = 0.5*(BWNV(M+1)+BWNV(M))
        DWNV(M)  = BWNV(M+1)-BWNV(M)
        WAVEV(M) = 1.0E+4/WNOV(M)
      end do

C     Sum the solar flux, and write out the result.  

      sum = 0.0
      do N=1,L_NSPECTV
        SOLARF(N) = SOLAR(N)
        sum       = sum+SOLARF(N)
      end do
C      write(6,'("Solar flux at 1AU = ",f7.2," W/M^2")') sum

C     Set up the wavelength independent part of Rayleigh Scattering.
C     The pressure dependent part will be computed elsewhere (OPTCV).
C     WAVEV is in microns.  There is no Rayleigh scattering in the IR.

      do N=1,L_NSPECTV
        WL        = WAVEV(N)
        TAURAY(N) = (8.7/grav)*(1.527*(1.0+0.013/wl**2)/wl**4)*
     &                scalep/P0
      end do

      RETURN
      end 
