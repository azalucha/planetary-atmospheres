C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C                            radinc_h
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C======================================================================C
C
C     RADINC.H    RADiation INCludes
C     
C     Includes for the radiation code; RADIATION LAYERS, LEVELS,
C     number of spectral intervals. . .
C     
C     GCM2.0  Feb 2003
C
C======================================================================C

C     RADIATION parameters

C     In radiation code, layer 1 corresponds to the stratosphere.  Level
C     1 is the top of the stratosphere.  The dummy layer is at the same
C     temperature as the (vertically isothermal) stratosphere, and
C     any time it is explicitly needed, the appropriate quantities will
C     be dealt with (aka "top". . .)

C     L_NLEVRAD corresponds to the surface - i.e., the GCM Level that
C     is at the surface.  PLEV(L_NLEVRAD) = P(J,I)+PTROP,
C     PLEV(2) = PTROP, PLEV(1) = ptrop

C     L_NLAYRAD is the number of radiation code layers
C     L_NLEVRAD is the number of radiation code levels.  Level N is the
C               top of layer N.
C     
C     L_NSPECTI is the number of IR spectral intervals
C     L_NSPECTV is the number of visible (or Solar) spectral intervals
C     L_NGAUSS  is the number of Gauss points for K-coefficients
C               GAUSS POINT 9 (aka the last one) is the special case
C     L_NWNGI   is L_NSPECTI*L_NGAUSS;  the total number of "intervals"
C               in the IR
C     L_NWNGV   is L_NSPECTV*L_NGAUSS;  the total number of "intervals"
C               in the VISIBLE
C
C     L_NPREF   is the number of reference pressures that the
C               k-coefficients are calculated on
C     L_PINT    is the number of Lagrange interpolated reference
C               pressures for the CO2 k-coefficients.
C     L_NTREF   is the number of refernce temperatures for the
C               k-coefficients
C     L_TAUMAX  is the largest optical depth - larger ones are set
C               to this value.
C     
C     MAXEXP    The largest value (TAU/MU) used in exp(-TAU/MU), where
C               TAU is the optical depth, and MU is the cosine of the
C               zenith angle.
C               
C     L_REFH2O  The number of different water-mixing ratio values for
C               the k-coefficients that are now CO2+H2O.
C
C     L_NREFI   The spectral interval number of the IR reference
C               wavelength (i.e. the 9 micron band) (8-12 microns)
C
C     L_NREFV   The spectral interval number of the visible reference
C               wavelength (i.e. the 0.67 micron band)
C               
C----------------------------------------------------------------------

      INTEGER L_NLAYRAD
      PARAMETER(L_NLAYRAD  = L_LAYERS+1)

      INTEGER L_NLEVRAD
      PARAMETER( L_NLEVRAD  = L_LAYERS+2)
      INTEGER L_NSPECTI
      PARAMETER(L_NSPECTI=5)

      INTEGER L_NSPECTV
      PARAMETER(L_NSPECTV=7)

      INTEGER L_NGAUSS
      PARAMETER(L_NGAUSS=17)

      INTEGER L_NPREF
      PARAMETER(L_NPREF=11)

      INTEGER L_NTREF
      PARAMETER(L_NTREF=7)

      INTEGER L_TAUMAX
      PARAMETER(L_TAUMAX=35)

      _RL MAXEXP
      PARAMETER(MAXEXP    = 35.0E0)

      INTEGER L_PINT
      PARAMETER(L_PINT=51)

      INTEGER L_REFH2O
      PARAMETER(L_REFH2O=10)

      INTEGER L_NREFV
      PARAMETER(L_NREFV=6)

      INTEGER L_NREFI
      PARAMETER(L_NREFI=4)

