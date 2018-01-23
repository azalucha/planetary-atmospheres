      program driver

!     V23 RT code 2010

      use grid_h
      use defines_h
      use standard_h
      use fccsave_h
      use radinc_h
      use radcommon_h
      use cldcommon_h
      use constants_h

      implicit none

!  PL and TL are the GCM pressures and temperatures at the layer
!  boundaries and midpoints.

      real*8 PL(L_LEVELS), TL(L_LEVELS)

!  PLEV & TLEV are the pressure and temperatures at the GCM levels
!  while PMID & TMID are the pressure and temperatures at the GCM
!  midpoint.

      real*8 PLEV(L_LEVELS), TLEV(L_LEVELS)
      real*8 TMID(L_LEVELS), PMID(L_LEVELS)

!  DIFFVT is the total diffuse visible flux for a given spectral 
!  interval

      real*8 DIFFVT

!  VISUAL

      real*8 DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      real*8 TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      real*8 COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 taugsurf(L_NSPECTV,L_NGAUSS-1)
      real*8 taucump(L_NSPECTV,L_NGAUSS)
      integer ngwv(L_NSPECTV)
      real*8  :: detau(L_NSPECTV,L_NGAUSS)

!  IR

      real*8 DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8 COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 taugsurfi(L_NSPECTI,L_NGAUSS-1)
      integer ngwi(L_NSPECTI)

!  Water mixing

      real*8 QH2O(L_LEVELS)

      real*8 scaleht
      real*8 SOL(L_NSPECTV)

      integer gcmlayers
      integer NLAYRAD, NLEVRAD, NSPECTI, NSPECTV, NPREF, NTREF
      integer K, L, NW, nn, nnn
      real*8  albi, albv
      real*8  acosz
      real*8  ans
      real*8  fluxid(L_NLAYRAD),fluxvd(L_NLAYRAD)
      real*8  heatingv(L_NLAYRAD), heatingir(L_NLAYRAD)
      real*8  total(L_NLAYRAD)
      real*8  fdmax, fdmin
      real*8  firmax, firmin, fvmax, fvmin, df
      real*4  tstrato, psfo, gto, coszo, alspo, tautoto
      real*8  tstratd, gtd

      real*8 FMNETI(L_NLAYRAD), FMNETV(L_NLAYRAD)
      real*8 fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD), NFLUXTOPI
      real*8 fluxupv(L_NLAYRAD), fluxdnv(L_NLAYRAD), NFLUXTOPV
      real*8 fluxdv(L_NLAYRAD), fluxdi(L_NLAYRAD)
      integer :: J, I

!  Clouds

      integer :: nlev
      real*8  :: pcld, tautotcld

!======================================================================C

      NLAYRAD = L_NLAYRAD
      NLEVRAD = L_NLEVRAD
      NSPECTI = L_NSPECTI
      NSPECTV = L_NSPECTV
      NPREF   = L_NPREF
      NTREF   = L_NTREF
      NLEV    = L_LEVELS

!  ALBI is the IR surface albedo, ALBV the surface albedo in the
!  visible.

      ALBI    = 0.00
      ALBV    = 0.24

!     Set up spectral intervals in the Solar (VISUAL) and then IR.
!     Read in the k-coefficients. . .

      call radsetup

      call ini_optdst(QEXTV,QSCATV,GV,QEXTI,QSCATI,GI,
     *                QXVDST,QXIDST,QSVDST,QSIDST,GVDST,GIDST,
     *                QEXTREFDST)

      call ini_optcld(QEXTVc,QSCATVc,GVc,QEXTIc,QSCATIc,GIc,
     *                QXVCLD,QXICLD,QSVCLD,QSICLD,GVCLD,GICLD,
     *                QEXTREFCLD)

!  PTROP is the pressure (mbar) of the troposphere
!  PSF is the surface pressure (mbar)
!  TAUTOT is the dust optical depth at the reference wavelength
!     (0.67 microns) and reference pressure (RPTAU)
!  TAUTOTCLD is the cloud optical depth
!  PCLD is the pressure where the cloud is placed
!  CONRNU is a parameter defining the vertical dust distribution
!  ACOSZ is the cosine of the solar zenith angle

      ptrop     = 2.0e-6
      psf       = 6.1
      TAUTOT    = 0.3
      CONRNU    = 0.03
      acosz     = 0.8
      tautotcld = 0.0
      pcld      = 0.5

      write(6,'("PTROP  = ",1pe10.3)') PTROP
      write(6,'("PSF    = ",f10.3)') psf
      write(6,'("COSZ   = ",f10.3)') ACOSZ
      write(6,'("ALSP   = ",f10.3)') albv
      write(6,'("TAUTOT = ",f10.3)') tautot

!  Make the optical depth reference pressure equal to the surface 
!  pressure.

      rptau = psf

!  QH2O is the mixing ratio of water.  

      do k=1,L_LEVELS
        QH2O(K) = 1.0D-7
      end do

!     rsdist = 2.428        ! Ls =   0   SCOSZ = 557
!     rsdist = 2.745        ! Ls =  90   SCOSZ = 493
!     rsdist = 2.147        ! Ls = 180   SCOSZ = 630
!     rsdist = 1.927        ! Ls = 270   SCOSZ = 702
!     rsdist = 2.255        !SOLAR FLUX AT MARS:   601.330
                            !used for tests with GCM 1-D model

!  RSDIST is the square of the sun-Mars distance, in AU.

      rsdist = 2.255

      gcmlayers = L_LAYERS

!  TSTRATD is the initial temperature; change as needed
     
      tstratd = 200.0
   
      do k=1,L_LEVELS
        tl(K) = tstratd
      end do

!  GTD is the ground temperature
!  psfo, tstrato, gto are output variables

      gtd     = TL(L_LEVELS)
      psfo    = Psf
      tstrato = tstratd
      gto     = gtd

!  Calculate the sigma values

      sigma(3) = 0.0
      do L=1,L_LAYERS
        K = 2*L+3
        sigma(K) = sigma(K-2)+DSIG(L)
      end do

      do K=4,L_LEVELS-1,2
        sigma(K) = 0.5*(SIGMA(K+1)+SIGMA(K-1))
      end do

!     Fill cumulative dust optical depth arrays (cum. dust optical
!     depth from the top of the atmosphere to the bottom of level K).

!  TAUREF is the dust optical depth at the reference wavelength, in
!  each layer.  

      IF(TAUTOT.LE.0.0) THEN
        do K=1,L_LEVELS+1
          tauref(K) = 0.0
        end do
        TAUCUM(L_LEVELS) = 0.0D0
      ELSE
        CALL dustprofile
      END IF

!     Fill local arrays with pressures at layer boundaries and
!     mid-points.

      PL(2) = PTROP/2.0
      DO K=3,L_LEVELS
        PL(K) = SIGMA(K)*(PSF-PTROP)+PTROP
      END DO

!  FILLPT is the interface subroutine that takes the P & T values
!  on the GCM grid, and puts them on the RT vertical grid.  The
!  radiation code uses PMID, TMID, PLEV, and TLEV values.

      call fillpt(pl,psf,ptrop,gtd,tstratd,tl,plev,tlev,pmid,
     *                   tmid)

!  Fill the TAUREF array, the dust column density for each GCM sub-layer

      if(TAUTOT.gt.0.0) then
        call filltaucum
      end if

!  Fill special bottom radiation level to zero.

      TAUREF(L_LEVELS+1) = 0.0

!     And now back to the regular code. . .

!     Calculate solar flux at the current mars distance

      ans = 0.0
      if(acosz.lt.1.0e-4) then
        do NW=1,L_NSPECTV
          SOL(NW) = 0.0
        end do
      else
        do NW=1,L_NSPECTV
          SOL(nw) = SOLARF(NW)/RSDIST
          ans     = ans+sol(NW)
        end do
      end if

      write(6,'("SOLAR FLUX AT MARS:  ",f8.3)') ANS 
 
!     Set up, and solve for, the solar (visual) fluxes, if the sun
!     is up

!  TAUREFCLD(K) is the cloud optical depth at the reference wavelength.
!  The value is referenced at level K, and measured from the next 
!  higher level.  This is not a cumulative value, but just the optical
!  depth in each sub-layer.

!  If the sun is up, calculate the solar fluxes, else set them to zero.

      if(acosz.ge.1.0e-4) then

        call cldprofile(psf,ptrop,nlev,sigma,pcld,tautotcld,taurefcld)

!  Calculate the optical depths in each layer, spectral interval,
!  and Gauss point.

        call OPTCV(DTAUV,TAUV,TAUCUMV,CO2V,PLEV,PFGASREF,
     *             TGASREF,QXVDST,QSVDST,GVDST,WBARV,COSBV,
     *             TAURAY,TAUREF,TMID,PMID,TAUGSURF,QH2O,WREFH2O,
     *             QEXTREFCLD,TAUREFCLD,QXVCLD,QSVCLD,GVCLD)

!  Compute the visible fluxes

        call SFLUXV(DTAUV,TAUV,TAUCUMV,ALBV,WBARV,COSBV,
     *              ACOSZ,SOL,GWEIGHT,NFLUXTOPV,FMNETV,
     *              FLUXUPV,FLUXDNV,DIFFVT,FZEROV,taugsurf,
     *              detau)
     
      else
        NFLUXTOPV = 0.0
        do L=1,L_NLAYRAD
          FMNETV(L)  = 0.0
          FLUXUPV(L) = 0.0
          FLUXDNV(L) = 0.0
        end do
      end if

!     Set up, and solve for, the Infrared fluxes

!  Calculate the optical depths in each layer, spectral interval,
!  and Gauss point.

      call OPTCI(DTAUI,TAUCUMI,CO2I,PLEV,PFGASREF,TGASREF,
     *           QEXTREFDST,QXIDST,QSIDST,GIDST,COSBI,WBARI,TAUREF,
     *           TMID,PMID,TAUGSURFI,QH2O,WREFH2O,
     *           QEXTREFCLD,TAUREFCLD,QXICLD,QSICLD,GICLD)
 
!  Compute the IR fluxes

      call SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,ALBI,DWNI,
     *            COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,
     *            fluxupi,fluxdni,FZEROI,taugsurfI)

!  Fluxes have been computed.  Below is code that outputs the
!  values for graphics/analysis.

!     Output for IDL

      open(60,file='driver_data')

!     Upward and downward flux

      firmax = FLUXUPI(1)
      firmin = FLUXUPI(1)
      fvmax  = FLUXUPV(1)
      fvmin  = FLUXUPV(1)
     
      do L=1,L_NLAYRAD
        if(FLUXUPI(L).gt.firmax) firmax = FLUXUPI(L) 
        if(FLUXUPI(L).lt.firmin) firmin = FLUXUPI(L) 
        if(FLUXUPV(L).gt.fvmax)  fvmax  = FLUXUPV(L) 
        if(FLUXUPV(L).lt.fvmin)  fvmin  = FLUXUPV(L) 
      end do

      df = 0.05*(firmax-firmin)
      firmax = firmax+df
      firmin = firmin-df
      
      df = 0.05*(fvmax-fvmin)
      fvmax = fvmax+df
      fvmin = fvmin-df
      
      write(60,'(i4)') NLAYRAD
      write(60,'(1pe15.5,3(2x,1pe15.5))') firmax, firmin, fvmax, fvmin

      do L=1,L_NLAYRAD
        write(60,'(1pe15.5,2(2x,1pe15.5))') plev(2*L+1), FLUXUPV(L),
     *                                      FLUXUPI(L)
      end do

!     Downward fluxes

      firmax = FLUXDNI(1)
      firmin = FLUXDNI(1)
      fvmax  = FLUXDNV(1)
      fvmin  = FLUXDNV(1)
     
      do L=1,L_NLAYRAD
        if(FLUXDNI(L).gt.firmax) firmax = FLUXDNI(L) 
        if(FLUXDNI(L).lt.firmin) firmin = FLUXDNI(L) 
        if(FLUXDNV(L).gt.fvmax)  fvmax  = FLUXDNV(L) 
        if(FLUXDNV(L).lt.fvmin)  fvmin  = FLUXDNV(L) 
      end do

      df = 0.05*(firmax-firmin)
      firmax = firmax+df
      firmin = firmin-df
      
      df = 0.05*(fvmax-fvmin)
      fvmax = fvmax+df
      fvmin = fvmin-df
      
      write(60,'(1pe15.5,3(2x,1pe15.5))') firmax, firmin, fvmax, fvmin

      do L=1,L_NLAYRAD
        write(60,'(1pe15.8,2(2x,1pe15.8))') plev(2*L+1), FLUXDNV(L),
     *                                      FLUXDNI(L)
      end do

!     Net fluxes (as well as T-profile)

      firmax = FMNETI(1)
      firmin = FMNETI(1)
      fvmax  = FMNETV(1)
      fvmin  = FMNETV(1)
     
      do L=1,L_NLAYRAD
        if(FMNETI(L).gt.firmax) firmax = FMNETI(L) 
        if(FMNETI(L).lt.firmin) firmin = FMNETI(L) 
        if(FMNETV(L).gt.fvmax)  fvmax  = FMNETV(L) 
        if(FMNETV(L).lt.fvmin)  fvmin  = FMNETV(L) 
      end do

      df = 0.05*(firmax-firmin)
      firmax = firmax+df
      firmin = firmin-df
      
      df = 0.05*(fvmax-fvmin)
      fvmax = fvmax+df
      fvmin = fvmin-df
      
      write(60,'(i4)') NLAYRAD
      write(60,'(f10.4,5(2x,f10.4))') acosz, albv, ans, 
     *           taucum(L_LEVELS), conrnu, tautotcld
      write(60,'(1pe15.5,3(2x,1pe15.5))') firmax, firmin, fvmax, fvmin
      do L=1,L_NLAYRAD-1
        write(60,'(1pe15.5,2(2x,1pe17.7))') plev(2*L+1), FMNETV(L),
     *                              FMNETI(L)
      end do

      L = L_NLAYRAD
      scaleht = 0.0
      write(60,'(1pe15.5,2(2x,1pe15.5))') psf, FMNETV(L),
     *                              FMNETI(L)

!  Flux divergence

      fluxdv(1) = FMNETV(1)-NFLUXTOPV
      fluxdi(1) = FMNETI(1)-NFLUXTOPI

      do L=2,L_NLAYRAD
        fluxdv(L) = FMNETV(L)-FMNETV(L-1)
        fluxdi(L) = FMNETI(L)-FMNETI(L-1)
      end do

!  Heating rates
     
      heatingv(1)  = (FMNETV(1)-NFLUXTOPV)*88775.0*grav/
     *                      (cp*scalep*PLEV(3))
      heatingir(1) = (FMNETI(1)-NFLUXTOPI)*88775.0*grav/
     *                      (cp*scalep*PLEV(3))
      total(1) = heatingv(1) + heatingir(1)

      fdmax = -10000.0
      fdmin =  10000.0

      do L=2,L_NLAYRAD
    
        heatingv(L)   = (FMNETV(L)-FMNETV(L-1))*88775.0*grav/(cp*scalep*
     *                  (PLEV(2*L+1)-PLEV(2*L-1)))
        heatingir(L)  = (FMNETI(L)-FMNETI(L-1))*88775.0*grav/(cp*scalep*
     *                (PLEV(2*L+1)-PLEV(2*L-1)))
        total(L) = heatingv(L) + heatingir(L)
      
        if(heatingv(L) .GT. fdmax)  fdmax = heatingv(L)
        if(heatingir(L) .GT. fdmax) fdmax = heatingir(L)
        if(heatingv(L) .LT. fdmin)  fdmin = heatingv(L)
        if(heatingir(L) .LT. fdmin) fdmin = heatingir(L)
      end do

      write(60,'(1pe10.3,3x,1pe10.3)') fdmin, fdmax
      do L=1,L_NLAYRAD
        write(60,'(1pe13.5,3(2x,1pe13.5))') plev(2*L), heatingv(L),
     *                                      heatingir(L), total(L)
      end do

!  Temperature profile

      do L=1,L_NLAYRAD
        write(60,'(1pe15.7,1x,1pe13.5)') plev(2*L), TLEV(2*L)
      end do

      write(60,'(f7.2,f10.3)') GTD, PSF

      end
