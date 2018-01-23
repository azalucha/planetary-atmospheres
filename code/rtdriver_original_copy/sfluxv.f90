      SUBROUTINE SFLUXV(DTAUV,TAUV,TAUCUMV,RSFV,WBARV,COSBV,      &
                        UBAR0,SOL,GWEIGHT,NFLUXTOPV,FMNETV,            &
                        FLUXUPV,FLUXDNV,DIFFVT,FZEROV,taugsurf,        &
                        detau)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

      use grid_h
      use radinc_h
      use radcommon_h, only: tlimits

      implicit none

!  DTAUV     - Visible opacity of layer L.
!  TAUV      - Visible opacity from top of the atmosphere to bottom
!              of layer L.
!  TAUCUMV   - Visible opacity from top of atmosphere to the bottom
!              of level K.
!  RSFV      - Surface albedo in the visible.
!  WBARV     - Scattering parameter
!  COSBV     - Scattering asymetry parameter
!  UBARO     - Cosine of the incidence angle
!  SOL       - Solar flux in each spectral band, at the top of the
!              atmosphere
!  GWEIGHT   - Gauss weight
!  NFLUXTOPV - Net flux at the top of the atmosphere 
!  FMNETV    - Net flux at the bottom of layer L
!  FLUXUPV   - Upward flux at layer L
!  FLUXDNV   - Downward flux at layer L
!  DIFFVT    - THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX
!  FZEROV    - Fraction of zeros (gas opacity k-coefficient off-line
!              computations) in each spectral interval.
!  TAUGSURF  - Gas optical depth (top of atmosphere to the surface).
!  DETAU     - Scaled ( delta-Eddington) optical depth at the surface.
!
!----------------------------------------------------------------------!

      real*8  :: FMNETV(L_NLAYRAD)
      real*8  :: TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      real*8  :: TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: FMUPV(L_NLAYRAD), FMDV(L_NLAYRAD)
      real*8  :: COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: SOL(L_NSPECTV)
      real*8  :: FLUXUPV(L_NLAYRAD), FLUXDNV(L_NLAYRAD)
      real*8  :: NFLUXTOPV, FLUXUP, FLUXDN
      real*8  :: GWEIGHT(L_NGAUSS)
 
!  delta-Eddington surface tau

      real*8  :: detau(L_NSPECTV,L_NGAUSS)

      integer L, NG, NW, NG1
      real*8  rsfv, ubar0, f0pi, btop, bsurf, taumax, eterm
      real*8 FZEROV(L_NSPECTV)

      real*8 DIFFV, DIFFVT

      real*8 taugsurf(L_NSPECTV,L_NGAUSS-1), fzero

!======================================================================C

      TAUMAX = L_TAUMAX

!     ZERO THE NET FLUXES

      NFLUXTOPV = 0.0

      DO L=1,L_NLAYRAD
        FMNETV(L)  = 0.0
        FLUXUPV(L) = 0.0
        FLUXDNV(L) = 0.0
      END DO

      DIFFVT = 0.0

!     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE VISIBLE
!     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      DO 500 NW=1,L_NSPECTV
      
        F0PI = SOL(NW)

        FZERO = FZEROV(NW)
        IF(FZERO.ge.0.99) goto 40
        DO NG=1,L_NGAUSS-1

          call GETDETAU(DTAUV(1,NW,NG),TAUV(1,NW,NG),                &
                        TAUCUMV(1,NW,NG),WBARV(1,NW,NG),             &
                        COSBV(1,NW,NG),UBAR0,detau(NW,NG))

          if(TAUGSURF(NW,NG) .lt. TLIMITS) then
            fzero = fzero + (1.0-FZEROV(NW))*GWEIGHT(NG)
            goto 30
          end if

!         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE

          BTOP = 0.0

!         LOOP OVER THE NTERMS BEGINNING HERE

!  detau(NW,NG) is the scaled optical depth at the surface

          ETERM = MIN(detau(NW,NG)/UBAR0,MAXEXP)
          BSURF = RSFV*UBAR0*SOL(NW)*EXP(-ETERM)

!         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
!         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
!         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
! 
!         FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO 
!         RETURN FLUXES FOR A GIVEN NT

          CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),   &
                      WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV,   &
                      BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN,       &
                      detau(nw,ng))
 
!         NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX 

          NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*GWEIGHT(NG)*           &
                                (1.0-FZEROV(NW))
          DO L=1,L_NLAYRAD
            FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*                  &
                                 GWEIGHT(NG)*(1.0-FZEROV(NW))
            FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*GWEIGHT(NG)*            &
                         (1.0-FZEROV(NW))
            FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*GWEIGHT(NG)*             &
                         (1.0-FZEROV(NW))
          END DO

!         THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX

          DIFFVT = DIFFVT + DIFFV*GWEIGHT(NG)*(1.0-FZEROV(NW))

   30     CONTINUE 

        END DO   ! the Gauss loop 

   40   continue 
!       Special 17th Gauss point

        NG = L_NGAUSS

!       SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
 
        BTOP = 0.0

        call GETDETAU(DTAUV(1,NW,NG),TAUV(1,NW,NG),                &
                      TAUCUMV(1,NW,NG),WBARV(1,NW,NG),             &
                      COSBV(1,NW,NG),UBAR0,detau(NW,NG))

!       LOOP OVER THE NTERMS BEGINNING HERE
 
        ETERM = MIN(detau(NW,NG)/UBAR0,MAXEXP)
        BSURF = RSFV*UBAR0*SOL(NW)*EXP(-ETERM)

!       WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
!       CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
!       WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
! 
!       FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO 
!       RETURN FLUXES FOR A GIVEN NT

        CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),     &
                    WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV,     &
                    BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN,         &
                    detau(nw,ng))
 
!       NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX 

        NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*FZERO
        DO L=1,L_NLAYRAD
          FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*FZERO
          FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*FZERO
          FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*FZERO
        END DO

!       THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX

        DIFFVT = DIFFVT + DIFFV*FZERO

  500 CONTINUE

!     *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE VISIBLE*****
 
      RETURN 
      end subroutine sfluxv
