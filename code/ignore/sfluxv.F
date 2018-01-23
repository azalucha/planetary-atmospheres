#include "CPP_OPTIONS.h"
      SUBROUTINE SFLUXV(DTAUV,TAUV,TAUCUMV,RSFV,WBARV,COSBV,      
     &                   UBAR0,SOL,GWEIGHT,NFLUXTOPV,FMNETV,            
     &                   FLUXUPV,FLUXDNV,DIFFVT,FZEROV,taugsurf,        
     &                   detau)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GRID_H.h"
#include "RADINC_H.h"

C  DTAUV     - Visible opacity of layer L.
C  TAUV      - Visible opacity from top of the atmosphere to bottom
C              of layer L.
C  TAUCUMV   - Visible opacity from top of atmosphere to the bottom
C              of level K.
C  RSFV      - Surface albedo in the visible.
C  WBARV     - Scattering parameter
C  COSBV     - Scattering asymetry parameter
C  UBARO     - Cosine of the incidence angle
C  SOL       - Solar flux in each spectral band, at the top of the
C              atmosphere
C  GWEIGHT   - Gauss weight
C  NFLUXTOPV - Net flux at the top of the atmosphere 
C  FMNETV    - Net flux at the bottom of layer L
C  FLUXUPV   - Upward flux at layer L
C  FLUXDNV   - Downward flux at layer L
C  DIFFVT    - THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX
C  FZEROV    - Fraction of zeros (gas opacity k-coefficient off-line
C              computations) in each spectral interval.
C  TAUGSURF  - Gas optical depth (top of atmosphere to the surface).
C  DETAU     - Scaled ( delta-Eddington) optical depth at the surface.
C
C----------------------------------------------------------------------!

      real*8   FMNETV(L_NLAYRAD)
      real*8   TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      real*8   TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      real*8   DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8   FMUPV(L_NLAYRAD), FMDV(L_NLAYRAD)
      real*8   COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8   WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8   SOL(L_NSPECTV)
      real*8   FLUXUPV(L_NLAYRAD), FLUXDNV(L_NLAYRAD)
      real*8   NFLUXTOPV, FLUXUP, FLUXDN
      real*8   GWEIGHT(L_NGAUSS)
 
C  delta-Eddington surface tau

      real*8  :: detau(L_NSPECTV,L_NGAUSS)

      integer L, NG, NW, NG1
      real*8  rsfv, ubar0, f0pi, btop, bsurf, taumax, eterm
      real*8 FZEROV(L_NSPECTV)

      real*8 DIFFV, DIFFVT

      real*8 taugsurf(L_NSPECTV,L_NGAUSS-1), fzero

      real*8 TLIMITS
      PARAMETER(TLIMITS=1.0E-32)

C======================================================================C

      TAUMAX = L_TAUMAX

C     ZERO THE NET FLUXES

      NFLUXTOPV = 0.0

      DO L=1,L_NLAYRAD
        FMNETV(L)  = 0.0
        FLUXUPV(L) = 0.0
        FLUXDNV(L) = 0.0
      END DO

      DIFFVT = 0.0

C     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE VISIBLE
C     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      DO 500 NW=1,L_NSPECTV
      
        F0PI = SOL(NW)

        FZERO = FZEROV(NW)
        IF(FZERO.ge.0.99) goto 40
        DO NG=1,L_NGAUSS-1

          call GETDETAU(DTAUV(1,NW,NG),TAUV(1,NW,NG),                
     &                   TAUCUMV(1,NW,NG),WBARV(1,NW,NG),             
     &                   COSBV(1,NW,NG),UBAR0,detau(NW,NG))

          if(TAUGSURF(NW,NG) .lt. TLIMITS) then
            fzero = fzero + (1.0-FZEROV(NW))*GWEIGHT(NG)
            goto 30
          end if

C         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE

          BTOP = 0.0

C         LOOP OVER THE NTERMS BEGINNING HERE

C  detau(NW,NG) is the scaled optical depth at the surface

          ETERM = MIN(detau(NW,NG)/UBAR0,MAXEXP)
          BSURF = RSFV*UBAR0*SOL(NW)*EXP(-ETERM)

C         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
C 
C         FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO 
C         RETURN FLUXES FOR A GIVEN NT

          CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),
     &                WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV,   
     &                BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN,       
     &                detau(nw,ng))
 
C         NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX 

          NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*GWEIGHT(NG)*           
     &                           (1.0-FZEROV(NW))
          DO L=1,L_NLAYRAD
            FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*                  
     &                            GWEIGHT(NG)*(1.0-FZEROV(NW))
            FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*GWEIGHT(NG)*            
     &                     (1.0-FZEROV(NW))
            FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*GWEIGHT(NG)*             
     &                    (1.0-FZEROV(NW))
          END DO

C         THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX

          DIFFVT = DIFFVT + DIFFV*GWEIGHT(NG)*(1.0-FZEROV(NW))

   30     CONTINUE 

        END DO   ! the Gauss loop 

   40   continue 
C       Special 17th Gauss point

        NG = L_NGAUSS

C       SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE VISIBLE
 
        BTOP = 0.0

        call GETDETAU(DTAUV(1,NW,NG),TAUV(1,NW,NG),                
     &                 TAUCUMV(1,NW,NG),WBARV(1,NW,NG),             
     &                 COSBV(1,NW,NG),UBAR0,detau(NW,NG))

C       LOOP OVER THE NTERMS BEGINNING HERE
 
        ETERM = MIN(detau(NW,NG)/UBAR0,MAXEXP)
        BSURF = RSFV*UBAR0*SOL(NW)*EXP(-ETERM)

C       WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C       CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C       WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER
C 
C       FUW AND FDW ARE WORKING FLUX ARRAYS THAT WILL BE USED TO 
C       RETURN FLUXES FOR A GIVEN NT

        CALL GFLUXV(DTAUV(1,NW,NG),TAUV(1,NW,NG),TAUCUMV(1,NW,NG),     
     &               WBARV(1,NW,NG),COSBV(1,NW,NG),UBAR0,F0PI,RSFV,     
     &               BTOP,BSURF,FMUPV,FMDV,DIFFV,FLUXUP,FLUXDN,         
     &               detau(nw,ng))
 
C       NOW CALCULATE THE CUMULATIVE VISIBLE NET FLUX 

        NFLUXTOPV = NFLUXTOPV+(FLUXUP-FLUXDN)*FZERO
        DO L=1,L_NLAYRAD
          FMNETV(L)=FMNETV(L)+( FMUPV(L)-FMDV(L) )*FZERO
          FLUXUPV(L) = FLUXUPV(L) + FMUPV(L)*FZERO
          FLUXDNV(L) = FLUXDNV(L) + FMDV(L)*FZERO
        END DO

C       THE DIFFUSE COMPONENT OF THE DOWNWARD SOLAR FLUX

        DIFFVT = DIFFVT + DIFFV*FZERO

  500 CONTINUE

C     *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE VISIBLE*****
 
      RETURN 
      end
