#include "CPP_OPTIONS.h"
      subroutine mars_rt_driver(bi,bj,kLev,myTime,myThid,etaN,
     &       theta,dTdt)

C !DESCRIPTION:
C *==========================================================*
C | SUBROUTINE MARS_RT_DRIVER
C | o Driver for NASA AMES Mars RT (V23 RT code 2010), at 
C | both visible and IR wavelengths due to CO2, dust, and H20.
C | Also calculates the transmission of the solar beam
C | at the surface, for use in STEP_SURFACE_T.F
C | author: amz
C *==========================================================*

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "AMZVARS.h"
#include "GRID.h"
#include "RT_GRID.h"
#include "RADINC.h"
#include "DEFINES.h"
#include "STANDARD.h"
#include "FCCSAVE.h"
#include "RADCOMMON.h"
#include "CLDCOMMON.h"
#include "CONSTANTS.h"

C Input variables (MIT GCM)
      INTEGER kLev, bi, bj
      _RL myTime,etaN(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL theta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER myThid


C  PL and TL are the GCM pressures and temperatures at the layer
C  boundaries and midpoints.

      _RL PL(L_LEVELS), TL(L_LEVELS)

C  PLEV & TLEV are the pressure and temperatures at the GCM levels
C  while PMID & TMID are the pressure and temperatures at the GCM
C  midpoint.

      _RL PLEV(L_LEVELS), TLEV(L_LEVELS)
      _RL TMID(L_LEVELS), PMID(L_LEVELS)

C  DIFFVT is the total diffuse visible flux for a given spectral 
C  interval

      _RL DIFFVT

C  VISUAL

      _RL DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      _RL TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      _RL TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      _RL COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      _RL WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      _RL taugsurf(L_NSPECTV,L_NGAUSS-1)
      _RL taucump(L_NSPECTV,L_NGAUSS)
      integer ngwv(L_NSPECTV)
      _RL  :: detau(L_NSPECTV,L_NGAUSS)

C  IR

      _RL DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      _RL COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL taugsurfi(L_NSPECTI,L_NGAUSS-1)
      integer ngwi(L_NSPECTI)

C  Water mixing

      _RL QH2O(L_LEVELS)

      _RL scaleht
      _RL SOL(L_NSPECTV)

      integer gcmlayers
      integer NLAYRAD, NLEVRAD, NSPECTI, NSPECTV, NPREF, NTREF
      integer K, L, NW, nn, nnn
      _RL  albi, albv
      _RL  acosz
      _RL  ans
      _RL  fluxid(L_NLAYRAD),fluxvd(L_NLAYRAD)
      _RL  heatingv(L_NLAYRAD), heatingir(L_NLAYRAD)
      _RL  total(L_NLAYRAD)
      _RL  fdmax, fdmin
      _RL  firmax, firmin, fvmax, fvmin, df
      _RL  psfo, coszo, alspo, tautoto
      _RL  tstratd, gtd

      _RL FMNETI(L_NLAYRAD), FMNETV(L_NLAYRAD)
      _RL fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD), NFLUXTOPI
      _RL fluxupv(L_NLAYRAD), fluxdnv(L_NLAYRAD), NFLUXTOPV
      _RL fluxdv(L_NLAYRAD), fluxdi(L_NLAYRAD)
      integer J, I

C  Clouds

      integer  nlev
      _RL   pcld, tautotcld

C Local variables
      _RL myTheta(Nr),lssol,mu0,myheating(Nr),p_heating(L_NLAYRAD)
      _RL dTdt
      INTEGER ii,jj

C======================================================================C

C      write(*,*) 'AMZ mars_rt_driver end vairable declarations'

      i=ii
      j=jj

C  thermodynamics loop counts down from Nr.  Do calculations
C  for entire column at first instance
      IF(kLev .EQ. Nr) THEN
       IF(initrad(ii,jj,bi,bj) .EQ. 1. .OR.
     &            MOD(myTime,radTimeStep).EQ.0.) THEN

        NLAYRAD = L_NLAYRAD
        NLEVRAD = L_NLEVRAD
        NSPECTI = L_NSPECTI
        NSPECTV = L_NSPECTV
        NPREF   = L_NPREF
        NTREF   = L_NTREF
        NLEV    = L_LEVELS

C  ALBI is the IR surface albedo, ALBV the surface albedo in the
C  visible.

        ALBI    = 0.00
CC      ALBV    = 0.24
        ALBV    = albedoAny(ii,jj,bi,bj)

C     Set up spectral intervals in the Solar (VISUAL) and then IR.
C     Read in the k-coefficients. . .


C       write(*,*) 'AMZ mars_rt_driver call radsetup'
        call radsetup

C       write(*,*) 'AMZ mars_rt_driver call ini_optdst'

        call ini_optdst(QEXTV,QSCATV,GV,QEXTI,QSCATI,GI,
     &                QXVDST,QXIDST,QSVDST,QSIDST,GVDST,GIDST,
     &                QEXTREFDST)


C       write(*,*) 'AMZ mars_rt_driver call ini_optcld'
        call ini_optcld(QEXTVc,QSCATVc,GVc,QEXTIc,QSCATIc,GIc,
     &                QXVCLD,QXICLD,QSVCLD,QSICLD,GVCLD,GICLD,
     &                QEXTREFCLD)

C Model time and solar zenith angle
C       write(*,*) 'AMZ mars_rt_driver call get_ls'
        IF (lsIsConstant .EQ. 0) THEN
         CALL GET_LS(myTime,lssol)
        ELSE
         lssol = lsperpetual
        ENDIF
C       write(*,*) 'AMZ mars_rt_driver call get_mu0'
        CALL GET_MU0(ii,jj,bi,bj,lssol,myTime,mu0)



C  PTROP is the pressure (mbar) of the troposphere
C  PSF is the surface pressure (mbar)
C  TAUTOT is the dust optical depth at the reference wavelength
C     (0.67 microns) and reference pressure (RPTAU)
C  TAUTOTCLD is the cloud optical depth
C  PCLD is the pressure where the cloud is placed
C  CONRNU is a parameter defining the vertical dust distribution
C  ACOSZ is the cosine of the solar zenith angle

        ptrop     = 2.0e-6
        psf       = (Ro_Surf(ii,jj,bi,bj)+etaN(ii,jj,bi,bj))/100.
        TAUTOT    = 0.3
        CONRNU    = 0.03
        acosz     = mu0
        tautotcld = 0.0
        pcld      = 0.5

CC      write(6,'("PTROP  = ",1pe10.3)') PTROP
CC      write(6,'("PSF    = ",f10.3)') psf
CC      write(6,'("COSZ   = ",f10.3)') ACOSZ
CC      write(6,'("ALSP   = ",f10.3)') albv
CC      write(6,'("TAUTOT = ",f10.3)') tautot

C  Make the optical depth reference pressure equal to the surface 
C  pressure.

        rptau = psf

C  QH2O is the mixing ratio of water.  

        do k=1,L_LEVELS
         QH2O(K) = 1.0D-7
        enddo

CC      rsdist = 2.428        ! Ls =   0   SCOSZ = 557
CC      rsdist = 2.745        ! Ls =  90   SCOSZ = 493
Cc      rsdist = 2.147        ! Ls = 180   SCOSZ = 630
Cc      rsdist = 1.927        ! Ls = 270   SCOSZ = 702
Cc      rsdist = 2.255        !SOLAR FLUX AT MARS:   601.330
C                            !used for tests with GCM 1-D model

C  RSDIST is the square of the sun-Mars distance, in AU.
C       rsdist = 2.255
        rsdist=semimajoraxis/1.496E11/
     &        (1. + eccentricity*COS(lssol - lsp))
     &       *(1. - eccentricity**2.)

        gcmlayers = L_LAYERS

C  psfo is an output variable

        psfo    = Psf

C  Calculate the sigma values

        sigma(3) = 0.0
        do L=1,L_LAYERS
         K = 2*L+3
         sigma(K) = sigma(K-2)+DSIG(L)
        enddo

        do K=4,L_LEVELS-1,2
         sigma(K) = 0.5*(SIGMA(K+1)+SIGMA(K-1))
        enddo

C     Fill cumulative dust optical depth arrays (cum. dust optical
C     depth from the top of the atmosphere to the bottom of level K).

C  TAUREF is the dust optical depth at the reference wavelength, in
C  each layer.  

        IF(TAUTOT.LE.0.0) THEN
         do K=1,L_LEVELS+1
          tauref(K) = 0.0
         enddo
         TAUCUM(L_LEVELS) = 0.0D0
        ELSE
C        write(*,*) 'AMZ mars_rt_driver call dustprofile'
         CALL dustprofile
        ENDIF

C     Fill local arrays with pressures at layer boundaries and
C     mid-points.

        PL(2) = PTROP/2.0
        DO K=3,L_LEVELS
         PL(K) = SIGMA(K)*(PSF-PTROP)+PTROP
        ENDDO

C Now interpolate mtgcm grid to tl
C       write(*,*) 'AMZ mars_rt_driver call interp_t_to_tl'
        DO k=1,Nr
         myTheta(k)=theta(ii,jj,k,bi,bj)
        ENDDO
C          write(*,*) 'jj=',jj,', bj=',bj,', ii=',ii,', bi=',bi
        CALL INTERP_T_TO_TL(pl,tl,L_LEVELS,myTheta)
        tl(L_LEVELS)=surfaceTPrime(ii,jj,bi,bj)
        tstratd=surfaceTPrime(ii,jj,bi,bj)
C  GTD is the ground temperature
        gtd     = TL(L_LEVELS)

C  FILLPT is the interface subroutine that takes the P & T values
C  on the GCM grid, and puts them on the RT vertical grid.  The
C  radiation code uses PMID, TMID, PLEV, and TLEV values.

C       write(*,*) 'AMZ mars_rt_driver call fillpt'
        call fillpt(pl,psf,ptrop,gtd,tstratd,tl,plev,tlev,pmid,
     &                   tmid)

C  Fill the TAUREF array, the dust column density for each GCM sub-layer

        if(TAUTOT.gt.0.0) then
C        write(*,*) 'AMZ mars_rt_driver call filltaucum'
         call filltaucum
        endif

C  Fill special bottom radiation level to zero.

        TAUREF(L_LEVELS+1) = 0.0

C     And now back to the regular code. . .

C     Calculate solar flux at the current mars distance

        ans = 0.0
        if(acosz.lt.1.0e-4) then
         do NW=1,L_NSPECTV
          SOL(NW) = 0.0
         enddo
        else
         do NW=1,L_NSPECTV
          SOL(nw) = SOLARF(NW)/RSDIST
          ans     = ans+sol(NW)
         enddo
        endif

C       write(6,'("SOLAR FLUX AT MARS:  ",f8.3)') ANS 
 
C     Set up, and solve for, the solar (visual) fluxes, if the sun
C     is up

C  TAUREFCLD(K) is the cloud optical depth at the reference wavelength.
C  The value is referenced at level K, and measured from the next 
C  higher level.  This is not a cumulative value, but just the optical
C  depth in each sub-layer.

C  If the sun is up, calculate the solar fluxes, else set them to zero.
 
        if(acosz.ge.1.0e-4) then

C        write(*,*) 'AMZ mars_rt_driver call cldprofile'
         call cldprofile(psf,ptrop,nlev,sigma,pcld,
     &                                     tautotcld,taurefcld)

C  Calculate the optical depths in each layer, spectral interval,
C  and Gauss point.



C        write(*,*) 'AMZ mars_rt_driver call optcv'
         call OPTCV(DTAUV,TAUV,TAUCUMV,CO2V,PLEV,PFGASREF,
     &             TGASREF,QXVDST,QSVDST,GVDST,WBARV,COSBV,
     &             TAURAY,TAUREF,TMID,PMID,TAUGSURF,QH2O,WREFH2O,
     &             QEXTREFCLD,TAUREFCLD,QXVCLD,QSVCLD,GVCLD)



C  Compute the visible fluxes
C        write(*,*) 'AMZ mars_rt_driver call sfluxv'
         call SFLUXV(DTAUV,TAUV,TAUCUMV,ALBV,WBARV,COSBV,
     &              ACOSZ,SOL,GWEIGHT,NFLUXTOPV,FMNETV,
     &              FLUXUPV,FLUXDNV,DIFFVT,FZEROV,taugsurf,
     &              detau)
     
        else
         NFLUXTOPV = 0.0
         do L=1,L_NLAYRAD
          FMNETV(L)  = 0.0
          FLUXUPV(L) = 0.0
          FLUXDNV(L) = 0.0
         enddo
        endif

C     Set up, and solve for, the Infrared fluxes

C  Calculate the optical depths in each layer, spectral interval,
C  and Gauss point.

C       write(*,*) 'AMZ mars_rt_driver call optci'
        call OPTCI(DTAUI,TAUCUMI,CO2I,PLEV,PFGASREF,TGASREF,
     &           QEXTREFDST,QXIDST,QSIDST,GIDST,COSBI,WBARI,TAUREF,
     &           TMID,PMID,TAUGSURFI,QH2O,WREFH2O,
     &           QEXTREFCLD,TAUREFCLD,QXICLD,QSICLD,GICLD)
 

C  Compute the IR fluxes
C       write(*,*) 'AMZ mars_rt_driver call sfluxi'
        call SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,ALBI,DWNI,
     &            COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,
     &            fluxupi,fluxdni,FZEROI,taugsurfI,PLANCKIR)


C  Flux divergence

        fluxdv(1) = FMNETV(1)-NFLUXTOPV
        fluxdi(1) = FMNETI(1)-NFLUXTOPI

        do L=2,L_NLAYRAD
         fluxdv(L) = FMNETV(L)-FMNETV(L-1)
         fluxdi(L) = FMNETI(L)-FMNETI(L-1)
        enddo

C  Heating rates
     
        heatingv(1)  = (FMNETV(1)-NFLUXTOPV)*gravity/
     &                     (atm_cp*scalep*PLEV(3))
        heatingir(1) = (FMNETI(1)-NFLUXTOPI)*gravity/
     &                      (atm_cp*scalep*PLEV(3))
        total(1) = heatingv(1) + heatingir(1)

        do L=2,L_NLAYRAD
    
         heatingv(L)   = (FMNETV(L)-FMNETV(L-1))*gravity/
     &                  (atm_cp*scalep*
     &                  (PLEV(2*L+1)-PLEV(2*L-1)))
         heatingir(L)  = (FMNETI(L)-FMNETI(L-1))*gravity
     &                /(atm_cp*scalep*
     &                (PLEV(2*L+1)-PLEV(2*L-1)))
         total(L) = heatingv(L) + heatingir(L)
      
        enddo

C Interpolate total heating rate to MIT GCM rC grid
C       write(*,*) 'AMZ mars_rt_driver call interp_total_to_rC'
        DO L=1,L_NLAYRAD
         p_heating(L)=plev(2*L)
        ENDDO
        call interp_total_to_rC(p_heating,total,
     &                                  L_NLAYRAD,myheating)

        DO k=1,Nr
         qr(ii,jj,k,bi,bj)=myheating(k)
        ENDDO
C Surface heating
C       write(*,*) 'AMZ mars_rt_driver call step_surface_t_2'
        CALL STEP_SURFACE_T_2(ii,jj,bi,bj,FLUXDNI(L_NLAYRAD),
     &                                        FLUXDNV(L_NLAYRAD))
        initrad(ii,jj,bi,bj)=0.
       ENDIF
      ENDIF
      dTdt=qr(ii,jj,kLev,bi,bj)*(rC(klev)/atm_po)**(-atm_kappa)
C         write(*,*) 'k=',klev,', jj=',jj,', bj=',bj,', qr=',dTdt




      return
      end
