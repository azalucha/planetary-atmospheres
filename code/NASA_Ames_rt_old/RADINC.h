C NASA Ames RT code
C module radinc_h
C #include RADINC_h

!======================================================================C
!
!     RADINC.H    RADiation INCludes
!     
!     Includes for the radiation code; RADIATION LAYERS, LEVELS,
!     number of spectral intervals. . .
!
!     GCM2.0  Feb 2003
!
!======================================================================C
!     RADIATION parameters

!     In radiation code, layer 1 corresponds to the stratosphere.  Level
!     1 is the top of the stratosphere.  The dummy layer is at the same
!     temperature as the (vertically isothermal) stratosphere, and
!     any time it is explicitly needed, the appropriate quantities will
!     be dealt with (aka "top". . .)

!     L_NLEVRAD corresponds to the surface - i.e., the GCM Level that
!     is at the surface.  PLEV(L_NLEVRAD) = P(J,I)+PTROP,
!     PLEV(2) = PTROP, PLEV(1) = ptrop

!     L_NLAYRAD is the number of radiation code layers
!     L_NLEVRAD is the number of radiation code levels.  Level N is the
!               top of layer N.
!
!     L_NSPECTI is the number of IR spectral intervals
!     L_NSPECTV is the number of visible (or Solar) spectral intervals
!     L_NGAUSS  is the number of Gauss points for K-coefficients
!               GAUSS POINT 9 (aka the last one) is the special case
!     L_NWNGI   is L_NSPECTI*L_NGAUSS;  the total number of "intervals"
!               in the IR
!     L_NWNGV   is L_NSPECTV*L_NGAUSS;  the total number of "intervals"
!               in the VISIBLE
!
!     L_NPREF   is the number of reference pressures that the
!               k-coefficients are calculated on
!     L_PINT    is the number of Lagrange interpolated reference
!               pressures for the CO2 k-coefficients.
!     L_NTREF   is the number of refernce temperatures for the
!               k-coefficients
!     L_TAUMAX  is the largest optical depth - larger ones are set
!               to this value.
!     
!     MAXEXP    The largest value (TAU/MU) used in exp(-TAU/MU), where
!               TAU is the optical depth, and MU is the cosine of the
!               zenith angle.
!
!     L_REFH2O  The number of different water-mixing ratio values for
!               the k-coefficients that are now CO2+H2O.
!     
!     L_NREFI   The spectral interval number of the IR reference
!               wavelength (i.e. the 9 micron band) (8-12 microns)
!     
!     L_NREFV   The spectral interval number of the visible reference
!               wavelength (i.e. the 0.67 micron band)
!     
!----------------------------------------------------------------------C

      integer L_NLAYRAD
      PARAMETER(L_NLAYRAD  = L_LAYERS+1)

      integer L_NLEVRAD 
      PARAMETER(L_NLEVRAD  = L_LAYERS+2)

      integer L_NSPECTI
      PARAMETER(L_NSPECTI =  5)

      integer L_NSPECTV
      PARAMETER(L_NSPECTV =  7)

      integer L_NGAUSS
      PARAMETER(L_NGAUSS  = 17)

      integer L_NPREF   
      PARAMETER(L_NPREF   = 11)

      integer L_NTREF
      PARAMETER(L_NTREF   =  7)

      integer L_TAUMAX
      PARAMETER(L_TAUMAX  = 35)

      _RL  MAXEXP
      PARAMETER( MAXEXP=35.0)

      integer L_PINT,L_REFH2O,L_NREFV,L_NREFI
      PARAMETER(L_PINT    = 51)
      PARAMETER(L_REFH2O  = 10)
      PARAMETER(L_NREFV   = 6)
      PARAMETER(L_NREFI   = 4)
