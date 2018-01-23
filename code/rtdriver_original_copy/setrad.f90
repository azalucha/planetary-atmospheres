      subroutine setrad(TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV, &
                        QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,     &
                        QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!     PURPOSE:
!        Set up values used by the radiation code, such as the CO2 gas
!     absorption coefficients.  True constants are defined, and the 
!     time-independent quantities used by the radiation code are 
!     calculated. 
!
!     TGASREF        - Temperatures of the opacity grid
!     PFGASREF       - Pressures (on fine mesh) of the opacity grid
!     CO2V           - Visible CO2 k-coefficients (CO2 gas opacity)
!     CO2I           - IR CO2 k-coefficients (CO2 gas opacity)
!     FZEROV         - Fraction of zeros in the visible (k-coefficient
!                      off-line calculation)
!     FZEROI         - Fraction of zeros in the IR (k-coefficient
!                      off-line calculation)
!
!     AEROSOL RADIATIVE OPTICAL CONSTANTS
!     Values are at the wavelenght interval center
!
!  Dust
!
!     MIE SCATTERING - Size distribution weighted
!     Qextv    - Extinction efficiency - in the visible.
!     Qscatv   - Scattering efficiency - in the visible.
!     WV       - Single scattering albedo - in the visible.
!     GV       - Asymmetry parameter - in the visible.
!
!     Qexti    - Extinction efficiency - in the infrared.
!     Qscati   - Scattering efficiency - in the infrared.
!     WI       - Single scattering albedo - in the infrared.
!     GI       - Asymmetry parameter - in the infrared.
!     
!  Water ice cloud constants
!
!     MIE SCATTERING - Size distribution weighted
!     Qextvc   - Extinction efficiency - in the visible.
!     Qscatvc  - Scattering efficiency - in the visible.
!     WVc      - Single scattering albedo - in the visible.
!     GVc      - Asymmetry parameter - in the visible.
!
!     Qextic   - Extinction efficiency - in the infrared.
!     Qscatic  - Scattering efficiency - in the infrared.
!     WIc      - Single scattering albedo - in the infrared.
!     GIc      - Asymmetry parameter - in the infrared.
!     
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h

      implicit none

      integer :: N, NS, ios, nd

      real*8  :: CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: CO2V(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8  :: PFGASREF(L_PINT)
      real*8  :: PGASREF(L_NPREF), TGASREF(L_NTREF)

      real*8  :: qextv(L_NSPECTV), qextvc(L_NSPECTV)
      real*8  :: qscatv(L_NSPECTV), qscatvc(L_NSPECTV)
      real*8  :: wv(L_NSPECTV), wvc(L_NSPECTV)
      real*8  :: gv(L_NSPECTV), gvc(L_NSPECTV)

      real*8  :: qexti(L_NSPECTI), qextic(L_NSPECTI)
      real*8  :: qscati(L_NSPECTI), qscatic(L_NSPECTI)
      real*8  :: wi(L_NSPECTI), wic(L_NSPECTI)
      real*8  :: gi(L_NSPECTI), gic(L_NSPECTI)

!     real*8 kvis, kir
      integer :: nt, np, nw, ng

      real*8  :: fzeroi(L_NSPECTI)
      real*8  :: fzerov(L_NSPECTV)
 
!=======================================================================

!  Read in the dust and water ice optical parameters:  Qext, Qscat,
!  w0, and g

!  Dust

      open(20,file='data/QEXT_DUST',status='old',iostat=ios)
      if(ios.ne.0) then
        write(6,'("setrad.f90:  Could not open dust Qext file")')
        stop
      end if

!  Visible

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTV 
        read(20,*) nd, Qextv(n), Qscatv(n), wv(n), gv(n)
      end do

!  IR

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTI 
        read(20,*) nd, Qexti(n), Qscati(n), wi(n), gi(n)
      end do

      close(20)

!  Water ice cloud

      open(20,file='data/QEXT_WATER',status='old',iostat=ios)
      if(ios.ne.0) then
        write(6,'("setrad.f90:  Could not open dust Qext file")')
        stop
      end if

!  Visible

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTV 
        read(20,*) nd, Qextvc(n), Qscatvc(n), wvc(n), gvc(n)
      end do

!  IR

      do n=1,3
        read(20,*) !header
      end do

      do n=1,L_NSPECTI 
        read(20,*) nd, Qextic(n), Qscatic(n), wic(n), gic(n)
      end do

      close(20)
      
!     Set the reference pressure and temperature arrays.  These are
!     the pressures and temperatures at which we have k-coefficients.

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
 
!     Insure that w0 < 1

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

!  Fill the water ice cloud variables

!     Insure that w0 < 1

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

!     Interpolate CO2 k coefficients to the finer pressure grid.

      call laginterp(PGASREF,PFGASREF,CO2I,CO2V,FZEROI,FZEROV)

      return
      end subroutine setrad
