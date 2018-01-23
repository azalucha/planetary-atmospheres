#include "CPP_OPTIONS.h"
      SUBROUTINE SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,RSFI,DWNI,  
     &                   COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,          
     &                   fluxupi,fluxdni,FZEROI,TAUGSURF,PLANCKIR)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

      implicit none
#include "EEPARAMS.h"
#include "RT_GRID.h"
#include "RADINC.h"

C  PLEV      - Pressures at levels (i.e. layer boundaries and 
C              mid-point)
C  TLEV      - Temperatures at levels
C  DTAUI     - IR opacity of layer L.
C  TAUCUMI   - IR opacity from top of atmosphere to the bottom
C              of level K.
C  RSFI      - Surface albedo in the IR.
C  DWNI      - Width of each IR spectral interval
C  WBARI     - Scattering parameter
C  COSBI     - Scattering asymetry parameter
C  UBARI     - Cosine of the incidence angle
C  GWEIGHT   - Gauss weight
C  NFLUXTOPI - Net flux at the top of the atmosphere 
C  FMNETI    - Net flux at the bottom of layer L
C  FLUXUPI   - Upward flux at layer L
C  FLUXDNI   - Downward flux at layer L
C  FZEROI    - Fraction of zeros (gas opacity k-coefficient off-line
C              computations) in each spectral interval.
c  TAUGSURF  - Gas optical depth (top of atmosphere to the surface).
C
C----------------------------------------------------------------------!

      _RL planckir(L_NSPECTI,8501)
      _RL TLIMITI
      PARAMETER(TLIMITI = 5.0D-32)


      integer  NLEVRAD, L, NW, NG, NTS, NTT

      _RL TLEV(L_LEVELS), PLEV(L_LEVELS)
      _RL TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      _RL FMNETI(L_NLAYRAD)
      _RL DWNI(L_NSPECTI)
      _RL DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL FMUPI(L_NLEVRAD), FMDI(L_NLEVRAD)
      _RL COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      _RL GWEIGHT(L_NGAUSS), NFLUXTOPI
      _RL FTOPUP

      _RL UBARI, RSFI, TSURF, BSURF, TTOP, BTOP, TAUTOP
      _RL PLANCK, PLTOP
      _RL fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD)
      _RL FZEROI(L_NSPECTI)
      _RL taugsurf(L_NSPECTI,L_NGAUSS-1), fzero

C======================================================================C
 
      NLEVRAD = L_NLEVRAD
 
C     ZERO THE NET FLUXES
    
      NFLUXTOPI = 0.0
 
      DO L=1,L_NLAYRAD
        FMNETI(L)  = 0.0
        FLUXUPI(L) = 0.0
        FLUXDNI(L) = 0.0
      END DO
 
C     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE INFRARED
C     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      TTOP  = TLEV(2)
      TSURF = TLEV(L_LEVELS)

      NTS   = TSURF*10.0D0-499
      NTT   = TTOP *10.0D0-499

      DO 501 NW=1,L_NSPECTI

C       SURFACE EMISSIONS - INDEPENDENT OF GAUSS POINTS

        BSURF = (1.-RSFI)*PLANCKIR(NW,NTS)
        PLTOP = PLANCKIR(NW,NTT)

C  If FZEROI(NW) = 1, then the k-coefficients are zero - skip to the
C  special Gauss point at the end.
 
        FZERO = FZEROI(NW)
        IF(FZERO.ge.0.99) goto 40
 
C       	WRITE(*,*) "FZERO=",FZERO

        DO NG=1,L_NGAUSS-1
         
          if(TAUGSURF(NW,NG).lt. TLIMITI) then
            fzero = fzero + (1.0-FZEROI(NW))*GWEIGHT(NG)
            goto 30
          end if

C         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE IR
C         CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
C         OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY
 
          TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
          BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP
 
C         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 
          
          CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),         
     &                 TAUCUMI(1,NW,NG),                                
     &                 WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,
     &                 RSFI,BTOP,BSURF,FTOPUP,FMUPI,FMDI,PLANCKIR)
 
C         NOW CALCULATE THE CUMULATIVE IR NET FLUX

          NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*GWEIGHT(NG)*           
     &                            (1.0-FZEROI(NW))

          DO L=1,L_NLEVRAD-1

C           CORRECT FOR THE WAVENUMBER INTERVALS

            FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*        
     &                               GWEIGHT(NG)*(1.0-FZEROI(NW))
            FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*   
     &                                 (1.0-FZEROI(NW))
            FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*GWEIGHT(NG)*    
     &                                 (1.0-FZEROI(NW))
          ENDDO

   30     CONTINUE

C End NGAUSS LOOP
       ENDDO 

   40  CONTINUE

C      SPECIAL 17th Gauss point

       NG     = L_NGAUSS

       TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
       BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP

C      WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
C      CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
C      WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 

       CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),            
     &                 TAUCUMI(1,NW,NG),                                
     &                 WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,
     &                 RSFI,BTOP,BSURF,FTOPUP,FMUPI,FMDI,PLANCKIR)
 
C      NOW CALCULATE THE CUMULATIVE IR NET FLUX

       NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*FZERO

       DO L=1,L_NLEVRAD-1

C        CORRECT FOR THE WAVENUMBER INTERVALS

         FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*FZERO
         FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*FZERO
         FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*FZERO
       ENDDO
  
  501 CONTINUE      !End Spectral Interval LOOP

C *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****

      RETURN
C      end subroutine sfluxi
      END
