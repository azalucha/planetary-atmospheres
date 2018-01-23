#include "CPP_OPTIONS.h"
      subroutine setrad(TGASREF,PFGASREF,CO2V,CO2I,
     &                QEXTV,QSCATV,WV,GV, 
     &         QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,     
     &            QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C     PURPOSE:
C        Set up values used by the radiation code, such as the CO2 gas
C     absorption coefficients.  True constants are defined, and the 
C     time-independent quantities used by the radiation code are 
C     calculated. 
C
C     TGASREF        - Temperatures of the opacity grid
C     PFGASREF       - Pressures (on fine mesh) of the opacity grid
C     CO2V           - Visible CO2 k-coefficients (CO2 gas opacity)
C     CO2I           - IR CO2 k-coefficients (CO2 gas opacity)
C     FZEROV         - Fraction of zeros in the visible (k-coefficient
C                      off-line calculation)
C     FZEROI         - Fraction of zeros in the IR (k-coefficient
C                      off-line calculation)
C
C     AEROSOL RADIATIVE OPTICAL CONSTANTS
C     Values are at the wavelenght interval center
C
C  Dust
C
C     MIE SCATTERING - Size distribution weighted
C     Qextv    - Extinction efficiency - in the visible.
C     Qscatv   - Scattering efficiency - in the visible.
C     WV       - Single scattering albedo - in the visible.
C     GV       - Asymmetry parameter - in the visible.
C
C     Qexti    - Extinction efficiency - in the infrared.
C     Qscati   - Scattering efficiency - in the infrared.
C     WI       - Single scattering albedo - in the infrared.
C     GI       - Asymmetry parameter - in the infrared.
C     
C  Water ice cloud constants
C
C     MIE SCATTERING - Size distribution weighted
C     Qextvc   - Extinction efficiency - in the visible.
C     Qscatvc  - Scattering efficiency - in the visible.
C     WVc      - Single scattering albedo - in the visible.
C     GVc      - Asymmetry parameter - in the visible.
C
C     Qextic   - Extinction efficiency - in the infrared.
C     Qscatic  - Scattering efficiency - in the infrared.
C     WIc      - Single scattering albedo - in the infrared.
C     GIc      - Asymmetry parameter - in the infrared.
C     
C----------------------------------------------------------------------C

      implicit none
#include "EEPARAMS.h"
#include "RT_GRID.h"
#include "RADINC.h"

      integer  N, NS, ios, nd

      _RL CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      _RL CO2V(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)
      _RL PFGASREF(L_PINT)
      _RL PGASREF(L_NPREF), TGASREF(L_NTREF)

      _RL qextv(L_NSPECTV), qextvc(L_NSPECTV)
      _RL qscatv(L_NSPECTV), qscatvc(L_NSPECTV)
      _RL wv(L_NSPECTV), wvc(L_NSPECTV)
      _RL gv(L_NSPECTV), gvc(L_NSPECTV)

      _RL qexti(L_NSPECTI), qextic(L_NSPECTI)
      _RL qscati(L_NSPECTI), qscatic(L_NSPECTI)
      _RL wi(L_NSPECTI), wic(L_NSPECTI)
      _RL gi(L_NSPECTI), gic(L_NSPECTI)

C     _RL kvis, kir
      integer :: nt, np, nw, ng

      _RL fzeroi(L_NSPECTI)
      _RL fzerov(L_NSPECTV)
 
C=======================================================================

C  Read in the dust and water ice optical parameters:  Qext, Qscat,
C  w0, and g

C  Dust

      open(20,file='QEXT_DUST',status='old',iostat=ios)
      if(ios.ne.0) then
        write(6,'("setrad.f:  Could not open dust Qext file")')
        stop
      end if

C  Visible

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTV 
        read(20,*) nd, Qextv(n), Qscatv(n), wv(n), gv(n)
      end do

C  IR

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTI 
        read(20,*) nd, Qexti(n), Qscati(n), wi(n), gi(n)
      end do

      close(20)

C  Water ice cloud

      open(20,file='QEXT_WATER',status='old',iostat=ios)
      if(ios.ne.0) then
        write(6,'("setrad.f:  Could not open dust Qext file")')
        stop
      end if

C  Visible

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTV 
        read(20,*) nd, Qextvc(n), Qscatvc(n), wvc(n), gvc(n)
      end do

C  IR

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTI 
        read(20,*) nd, Qextic(n), Qscatic(n), wic(n), gic(n)
      end do

      close(20)
      
C     Set the reference pressure and temperature arrays.  These are
C     the pressures and temperatures at which we have k-coefficients.

      pgasref( 1) = 1.0E-6
      pgasref( 2) = 1.0E-5
      pgasref( 3) = 1.0E-4
      pgasref( 4) = 1.0E-3
      pgasref( 5) = 1.0E-2
      pgasref( 6) = 1.0E-1
      pgasref( 7) = 1.0
      pgasref( 8) = 1.0E+1
      pgasref( 9) = 1.0E+2
      pgasref(10) = 1.0E+3
      pgasref(11) = 1.0E+4

      tgasref(1)  =  50.0
      tgasref(2)  = 100.0
      tgasref(3)  = 150.0
      tgasref(4)  = 200.0
      tgasref(5)  = 250.0
      tgasref(6)  = 300.0
      tgasref(7)  = 350.0
 
C     Insure that w0 < 1

      DO N=1,L_NSPECTV
        IF(wv(n).ge.0.9999D0) then
          Qscatv(n) = 0.9999D0*Qextv(n)
          wv(n)     = 0.9999D0
        END IF
      END DO

      DO N=1,L_NSPECTI
        IF(wi(n).ge.0.9999D0) then
          Qscati(n) = 0.9999D0*Qexti(n)
          wi(n)     = 0.9999D0
        END IF
      END DO

C  Fill the water ice cloud variables

C     Insure that w0 < 1

      DO N=1,L_NSPECTV
        IF(wvc(n).ge.0.9999D0) then
          Qscatvc(n) = 0.9999D0*Qextvc(n)
          wvc(n)     = 0.9999D0
        END IF
      END DO

      DO N=1,L_NSPECTI
        IF(wic(n).ge.0.9999D0) then
          Qscatic(n) = 0.9999D0*Qextic(n)
          wic(n)     = 0.9999D0
        END IF
      END DO

C     Interpolate CO2 k coefficients to the finer pressure grid.

      call laginterp(PGASREF,PFGASREF,CO2I,CO2V,FZEROI,FZEROV)

      return
C      end subroutine setrad 
      END
