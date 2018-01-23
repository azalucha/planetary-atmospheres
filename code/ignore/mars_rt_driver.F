#include "CPP_OPTIONS.h"
      SUBROUTINE MARS_RT_DRIVER(bi,bj,kLev,myTime)

C     *==========================================================*
C     | S/R MARS_RT_DRIVER
C     | o Driver for NASA Ames GCM radiation scheme, V23 RT code 2010
C     |   
C     *==========================================================*

      IMPLICIT NONE

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"
#include "GRID_H.h"
#include "DEFINES_H.h"
#include "STANDARD_H.h"
#include "RADINC_H.h"
#include "FCCSAVE_H.h"
#include "RADCOMMON_H.h"
#include "CLDCOMMON_H.h"
#include "CONSTANTS_H.h"

      INTEGER bi,bj,kLev
      _RL myTime

C  PL and TL are the GCM pressures and temperatures at the layer
C  boundaries and midpoints.

      real*8 PL(L_LEVELS), TL(L_LEVELS)

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
      _RL ngwv(L_NSPECTV)
      _RL detau(L_NSPECTV,L_NGAUSS)

C  IR

      _RL DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      _RL COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL taugsurfi(L_NSPECTI,L_NGAUSS-1)
      INTEGER ngwi(L_NSPECTI)

C  Water mixing

      _RL QH2O(L_LEVELS)

      _RL scaleht
      _RL SOL(L_NSPECTV)

      INTEGER gcmlayers
      INTEGER NLAYRAD, NLEVRAD, NSPECTI, NSPECTV, NPREF, NTREF
      INTEGER K, L, NW, nn, nnn
      _RL  albi, albv
      _RL  acosz
      _RL  ans
      _RL  fluxid(L_NLAYRAD),fluxvd(L_NLAYRAD)
      _RL  heatingv(L_NLAYRAD), heatingir(L_NLAYRAD)
      _RL  total(L_NLAYRAD)
      _RL  fdmax, fdmin
      _RL  firmax, firmin, fvmax, fvmin, df
      _RL   gtd

      _RL FMNETI(L_NLAYRAD), FMNETV(L_NLAYRAD)
      _RL fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD), NFLUXTOPI
      _RL fluxupv(L_NLAYRAD), fluxdnv(L_NLAYRAD), NFLUXTOPV
      _RL fluxdv(L_NLAYRAD), fluxdi(L_NLAYRAD)
      INTEGER  J, I

C  Clouds

      INTEGER  nlev
      _RL   pcld, tautotcld

C Local variables
      _RL lssol,mu0,termP(Nr),temperature(Nr),tout
      _RL termPtmp(Nr),plevstot(L_NLAYRAD),dy,qr1,dTdt
      _RL temperaturetmp(Nr)

C======================================================================C

      DO j=1,sNy
       DO i=1,sNx
      write(*,*) 'i=',i,', j=',j,', bi=',bi,', bj=',bj
C  thermodynamics loop counts down from Nr.  Do calculations
C  for entire column at first instance
        IF(kLev .EQ. Nr) THEN
         IF(initrad(i,j,bi,bj) .EQ. 1. .OR.
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
          ALBV    = albedoany(i,j,bi,bj)

C     Set up spectral intervals in the Solar (VISUAL) and then IR.
C     Read in the k-coefficients. . .

          call radsetup(DWNV,WAVEV,SOLARF,TAURAY,
     & WNOI,DWNI,WAVEI,
     & TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV0,
     & QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,
     & QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV,PLANCKIR)


          call ini_optdst(QEXTV,QSCATV,GV,QEXTI,QSCATI,GI,
     &                QXVDST,QXIDST,QSVDST,QSIDST,GVDST,GIDST,
     &                QEXTREFDST)

          call ini_optcld(QEXTVc,QSCATVc,GVc,QEXTIc,QSCATIc,GIc,
     &                QXVCLD,QXICLD,QSVCLD,QSICLD,GVCLD,GICLD,
     &                QEXTREFCLD)

C  PTROP is the pressure (mbar) of the troposphere
C  PSF is the surface pressure (mbar)
C  TAUTOT is the dust optical depth at the reference wavelength
C     (0.67 microns) and reference pressure (RPTAU)
C  TAUTOTCLD is the cloud optical depth
C  PCLD is the pressure where the cloud is placed
C  CONRNU is a parameter defining the vertical dust distribution
C  ACOSZ is the cosine of the solar zenith angle

          ptrop     = 2.0e-6
          psf       = (Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))/100.
          TAUTOT    = 0.3
          CONRNU    = 0.03

          IF (lsIsConstant .EQ. 0) THEN
           CALL GET_LS(myTime,lssol)
          ELSE
           lssol = lsperpetual
          ENDIF

          CALL GET_MU0(i,j,bi,bj,lssol,myTime,mu0)

          acosz     = mu0
          tautotcld = 0.0
          pcld      = 0.5

C      write(6,'("PTROP  = ",1pe10.3)') PTROP
C      write(6,'("PSF    = ",f10.3)') psf
C      write(6,'("COSZ   = ",f10.3)') ACOSZ
C      write(6,'("ALSP   = ",f10.3)') albv
C      write(6,'("TAUTOT = ",f10.3)') tautot

C  Make the optical depth reference pressure equal to the surface 
C  pressure.

          rptau = psf

C  QH2O is the mixing ratio of water.  

          do k=1,L_LEVELS
            QH2O(K) = 1.0D-7
          end do

C     rsdist = 2.428        ! Ls =   0   SCOSZ = 557
C     rsdist = 2.745        ! Ls =  90   SCOSZ = 493
C     rsdist = 2.147        ! Ls = 180   SCOSZ = 630
C     rsdist = 1.927        ! Ls = 270   SCOSZ = 702
C     rsdist = 2.255        !SOLAR FLUX AT MARS:   601.330
                            !used for tests with GCM 1-D model

C  RSDIST is the square of the sun-Mars distance, in AU.

C         rsdist = 2.255
          rsdist = (semimajoraxis*(1-eccentricity**2.)/
     &              (1.+eccentricity*COS(lssol - lsp)))**2.

          gcmlayers = L_LAYERS

C  GTD is the ground temperature

          gtd= SurfaceTPrime(i,j,bi,bj)

C  Calculate the sigma values

          sigma(3) = 0.0
          do L=1,L_LAYERS
           K = 2*L+3
           sigma(K) = sigma(K-2)+DSIG(L)
          end do

          do K=4,L_LEVELS-1,2
           sigma(K) = 0.5*(SIGMA(K+1)+SIGMA(K-1))
          end do

C     Fill cumulative dust optical depth arrays (cum. dust optical
C     depth from the top of the atmosphere to the bottom of level K).

C  TAUREF is the dust optical depth at the reference wavelength, in
C  each layer.  

          IF(TAUTOT.LE.0.0) THEN
           do K=1,L_LEVELS+1
            tauref(K) = 0.0
           end do
           TAUCUM(L_LEVELS) = 0.0D0
          ELSE
           CALL dustprofile(prdst,taudst,rptau,psf,conrnu,ptrop)
          END IF

C     Fill local arrays with pressures at layer boundaries and
C     mid-points.

          PL(2) = PTROP/2.0
          DO K=3,L_LEVELS
           PL(K) = SIGMA(K)*(PSF-PTROP)+PTROP
          END DO

C amz: must interpolate MIT GCM grid to pl/tl grid
C remember that grids count in opposite directions!
          DO k=1,Nr
           termPtmp(k)=rC(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &           /Ro_surf(i,j,bi,bj)
           temperaturetmp(k)=theta(i,j,k,bi,bj)*
     &                     (termPtmp(k)/atm_po)**atm_kappa
          ENDDO

          DO k=1,Nr
           termP(Nr-k+1)=termPtmp(k)/100.
           temperature(Nr-k+1)=temperaturetmp(k)
          ENDDO

          do k=1,L_LEVELS
           IF(pl(k) .LT. termP(1)) THEN
            tl(k)=temperature(1)
           ELSEIF(pl(k) .GT. termP(Nr)) THEN
            tl(k)=temperature(Nr)
           ELSE
            CALL POLINT(termP,temperature,Nr,pl(k),tout,dy)
            tl(K) = tout
           ENDIF
          end do



C  FILLPT is the interface subroutine that takes the P & T values
C  on the GCM grid, and puts them on the RT vertical grid.  The
C  radiation code uses PMID, TMID, PLEV, and TLEV values.

          call fillpt(pl,psf,ptrop,gtd,tl,plev,tlev,pmid,
     &                   tmid)

C  Fill the TAUREF array, the dust column density for each GCM sub-layer

          if(TAUTOT.gt.0.0) then
           call filltaucum(tauref,taucum,sigma,prdst,taudst,
     &           psf,ptrop)
          end if

C  Fill special bottom radiation level to zero.

          TAUREF(L_LEVELS+1) = 0.0

C     And now back to the regular code. . .

C     Calculate solar flux at the current mars distance

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

 
C     Set up, and solve for, the solar (visual) fluxes, if the sun
C     is up

C  TAUREFCLD(K) is the cloud optical depth at the reference wavelength.
C  The value is referenced at level K, and measured from the next 
C  higher level.  This is not a cumulative value, but just the optical
C  depth in each sub-layer.

C  If the sun is up, calculate the solar fluxes, else set them to zero.

          if(acosz.ge.1.0e-4) then

           call cldprofile(psf,ptrop,nlev,sigma,pcld,
     &                                tautotcld,taurefcld)

C  Calculate the optical depths in each layer, spectral interval,
C  and Gauss point.

           call OPTCV(DTAUV,TAUV,TAUCUMV,CO2V,PLEV,PFGASREF,
     &             TGASREF,QXVDST,QSVDST,GVDST,WBARV,COSBV,
     &             TAURAY,TAUREF,TMID,PMID,TAUGSURF,QH2O,WREFH2O,
     &             QEXTREFCLD,TAUREFCLD,QXVCLD,QSVCLD,GVCLD)

C  Compute the visible fluxes

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
           end do
          end if

C     Set up, and solve for, the Infrared fluxes

C  Calculate the optical depths in each layer, spectral interval,
C  and Gauss point.

          call OPTCI(DTAUI,TAUCUMI,CO2I,PLEV,PFGASREF,TGASREF,
     &           QEXTREFDST,QXIDST,QSIDST,GIDST,COSBI,WBARI,TAUREF,
     &           TMID,PMID,TAUGSURFI,QH2O,WREFH2O,
     &           QEXTREFCLD,TAUREFCLD,QXICLD,QSICLD,GICLD)
 
C  Compute the IR fluxes

          call SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,ALBI,DWNI,
     &            COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,
     &            fluxupi,fluxdni,FZEROI,taugsurfI,PLANCKIR)

C  Fluxes have been computed. 


C     Upward and downward flux

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
     
C     Downward fluxes

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
      
C     Net fluxes (as well as T-profile)

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
      

          L = L_NLAYRAD
          scaleht = 0.0

C  Flux divergence

          fluxdv(1) = FMNETV(1)-NFLUXTOPV
          fluxdi(1) = FMNETI(1)-NFLUXTOPI

          do L=2,L_NLAYRAD
           fluxdv(L) = FMNETV(L)-FMNETV(L-1)
           fluxdi(L) = FMNETI(L)-FMNETI(L-1)
          end do

C  Heating rates (K/s)
C amz: remember to reverse order of pgrid     
          heatingv(1)  = (FMNETV(1)-NFLUXTOPV)*grav/
     &                       (cp*scalep*PLEV(3))
          heatingir(1) = (FMNETI(1)-NFLUXTOPI)*grav/
     &                        (cp*scalep*PLEV(3))
          total(L_NLAYRAD) = heatingv(1) + heatingir(1)
          plevstot(L_NLAYRAD)=PLEV(3)

          fdmax = -10000.0
          fdmin =  10000.0

          do L=2,L_NLAYRAD
    
           heatingv(L)   = (FMNETV(L)-FMNETV(L-1))*grav/
     &                  (cp*scalep*(PLEV(2*L+1)-PLEV(2*L-1)))
           heatingir(L)  = (FMNETI(L)-FMNETI(L-1))*grav/(cp*scalep*
     &                (PLEV(2*L+1)-PLEV(2*L-1)))
           total(L_NLAYRAD-L+1) = heatingv(L) + heatingir(L)
           plevstot(L_NLAYRAD-L+1)=PLEV(2*L)
          enddo
C amz: interpolate heating rate back to MIT GCM grid
          do k=1,Nr
           IF(rC(k) .LE. plevstot(L_NLAYRAD)) THEN
            qr1=total(L_NLAYRAD)
           ELSEIF(rC(k) .GE. plevstot(1)) THEN
            qr1=total(1)
           ELSE
            CALL POLINT(plevstot,total,L_NLAYRAD,
     &            rC(k),qr1,dy)
           ENDIF
           qr(i,j,k,bi,bj)=qr1
          ENDDO           

C amz: step surface temperature
          CALL STEP_SURFACE_T_2(i,j,bi,bj,
     &           FLUXDNI(L_NLAYRAD),FLUXDNV(L_NLAYRAD))
         ENDIF
        ENDIF
        dTdt=qr(i,j,kLev,bi,bj)*(rC(klev)/atm_po)**(-atm_kappa)
        gT(i,j,kLev,bi,bj)=gT(i,j,kLev,bi,bj)+dTdt
     &               *maskC(i,j,kLev,bi,bj)
       ENDDO
      ENDDO

      RETURN
      END
