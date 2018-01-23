#include "CPP_OPTIONS.h"
      SUBROUTINE GFLUXI(NLL,TLEV,NW,DW,DTAU,TAUCUM,W0,COSBAR,    
     &                 RSF,BTOP,BSURF,FTOPUP,FMIDP,FMIDM,PLANCKIR)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
C  FOR THE INFRARED FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
C  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
C  MEASURED FROM THE TOP OF EACH LAYER.  THE TOP OF EACH LAYER HAS  
C  OPTICAL DEPTH ZERO.  IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
C THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY 
C VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
C
C NLL            = NUMBER OF LEVELS (NLAYERS + 1) MUST BE LESS THAT NL (101)
C TLEV(L_LEVELS) = ARRAY OF TEMPERATURES AT GCM LEVELS
C WAVEN          = WAVELENGTH FOR THE COMPUTATION
C DW             = WAVENUMBER INTERVAL
C DTAU(NLAYER)   = ARRAY OPTICAL DEPTH OF THE LAYERS
C W0(NLEVEL)     = SINGLE SCATTERING ALBEDO
C COSBAR(NLEVEL) = ASYMMETRY FACTORS, 0=ISOTROPIC
C UBARI          = AVERAGE ANGLE, MUST BE EQUAL TO 0.5 IN IR
C RSF            = SURFACE REFLECTANCE
C BTOP           = UPPER BOUNDARY CONDITION ON IR INTENSITY (NOT FLUX)
C BSURF          = SURFACE EMISSION = (1-RSFI)*PLANCK, INTENSITY (NOT FLUX)
C FP(NLEVEL)     = UPWARD FLUX AT LEVELS
C FM(NLEVEL)     = DOWNWARD FLUX AT LEVELS
C FMIDP(NLAYER)  = UPWARD FLUX AT LAYER MIDPOINTS
C FMIDM(NLAYER)  = DOWNWARD FLUX AT LAYER MIDPOINTS
C
C----------------------------------------------------------------------C

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GRID_H.h"
#include "RADINC_H.h"

      integer NLP    ! MUST BE LARGER THAN NLEVEL 
      PARAMETER(NLP=101)

      INTEGER  NLL, NLAYER, L, NW, NT, NT2
      REAL*8   TERM, CPMID, CMMID
      REAL*8   PLANCK
      REAL*8   EM,EP
      REAL*8   COSBAR(L_NLAYRAD), W0(L_NLAYRAD), DTAU(L_NLAYRAD)
      REAL*8   TAUCUM(L_LEVELS),DTAUK
      REAL*8   TLEV(L_LEVELS)
      REAL*8   WAVEN, DW, UBARI, RSF
      REAL*8   BTOP, BSURF, FMIDP(L_NLAYRAD), FMIDM(L_NLAYRAD)
      REAL*8   B0(NLP),B1(NLP),ALPHA(NLP),LAMDA(NLP),XK1(NLP),XK2(NLP)
      REAL*8   GAMA(NLP),CP0(NLP),CM(NLP),CPM1(NLP),CMM1(NLP),E1(NLP)
      REAL*8   E2(NLP),E3(NLP),E4(NLP)

      REAL*8   FTOPUP, FLUXUP, FLUXDN
      real*8   TAUMAX
      PARAMETER(TAUMAX = L_TAUMAX) 
      real*8 planckir(L_NSPECTI,8501)

C======================================================================C

      write(*,*) 'amz flag 1'
      flush(6)

C     WE GO WITH THE HEMISPHERIC CONSTANT APPROACH IN THE INFRARED
 
      IF (NLL .GT. NLP) STOP 'PARAMETER NL TOO SMALL IN GLUFV'

      write(*,*) 'amz flag 1a'
      flush(6)


      NLAYER = L_NLAYRAD

      write(*,*) 'amz flag 1b'
      flush(6)


      DO L=1,L_NLAYRAD-1
      write(*,*) 'amz flag 1c, L=',L
      flush(6)
        ALPHA(L) = SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )
      write(*,*) 'amz flag 1d'
      flush(6)
        LAMDA(L) = ALPHA(L)*(1.0-W0(L)*COSBAR(L))/UBARI
      write(*,*) 'amz flag 1e'
      flush(6)

        NT2   = TLEV(2*L+2)*10.0D0-499
      write(*,*) 'amz flag 1f'
      flush(6)

        NT    = TLEV(2*L)*10.0D0-499
      write(*,*) 'amz flag 1g'
      flush(6)
        write(*,*) 'NT2=',NT2,', NT=',NT

        B1(L) = (PLANCKIR(NW,NT2)-PLANCKIR(NW,NT))/DTAU(L)
      write(*,*) 'amz flag 1h'
      flush(6)

        B0(L) = PLANCKIR(NW,NT)
      write(*,*) 'amz flag 1i'
      flush(6)

      END DO

      write(*,*) 'amz flag 2'
      flush(6)


C     Take care of special lower layer

      L        = L_NLAYRAD
      ALPHA(L) = SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )
      LAMDA(L) = ALPHA(L)*(1.0-W0(L)*COSBAR(L))/UBARI

      NT    = TLEV(2*L+1)*10.0D0-499
      NT2   = TLEV(2*L)*10.0D0-499
      B1(L) = (PLANCKIR(NW,NT)-PLANCKIR(NW,NT2))/DTAU(L)
      B0(L) = PLANCKIR(NW,NT2)

      write(*,*) 'amz flag 3'
      flush(6)

      DO L=1,L_NLAYRAD
        GAMA(L) = (1.0-ALPHA(L))/(1.0+ALPHA(L))
        TERM    = UBARI/(1.0-W0(L)*COSBAR(L))

C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C       BOTTOM OF THE LAYER.  THAT IS AT DTAU OPTICAL DEPTH

        CP0(L) = B0(L)+B1(L)*DTAU(L) +B1(L)*TERM
        CM(L) = B0(L)+B1(L)*DTAU(L) -B1(L)*TERM

C       CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C       AT THE TOP OF THE LAYER, THAT IS ZERO OPTICAL DEPTH

        CPM1(L) = B0(L)+B1(L)*TERM
        CMM1(L) = B0(L)-B1(L)*TERM
      END DO

      write(*,*) 'amz flag 4'
      flush(6)

 
C     NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C     FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C     WARNING IF DTAU(J) IS GREATER THAN ABOUT 35 (VAX)
C     WE CLIP IT TO AVOID OVERFLOW.

      DO L=1,L_NLAYRAD

C       CLIP THE EXPONENTIAL HERE.

        EP    = EXP( MIN((LAMDA(L)*DTAU(L)),TAUMAX))
        EM    = 1.0/EP
        E1(L) = EP+GAMA(L)*EM
        E2(L) = EP-GAMA(L)*EM
        E3(L) = GAMA(L)*EP+EM
        E4(L) = GAMA(L)*EP-EM
      END DO
 
C     B81=BTOP  ! RENAME BEFORE CALLING DSOLVER - used to be to set
C     B82=BSURF ! them to real*8 - but now everything is real*8
C     R81=RSF   ! so this may not be necessary

C     Double precision tridiagonal solver

      write(*,*) 'amz flag 5'
      flush(6)


      CALL DSOLVER(NLAYER,GAMA,CP0,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,       
     &              BSURF,RSF,XK1,XK2)
 
      write(*,*) 'amz flag 6'
      flush(6)


C     NOW WE CALCULATE THE FLUXES AT THE MIDPOINTS OF THE LAYERS.
 
      DO L=1,L_NLAYRAD-1
        DTAUK = TAUCUM(2*L+1)-TAUCUM(2*L)
        EP    = EXP(MIN(LAMDA(L)*DTAUK,TAUMAX)) ! CLIPPED EXPONENTIAL 
        EM    = 1.0/EP
        TERM  = UBARI/(1.-W0(L)*COSBAR(L))

C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C       BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

        CPMID    = B0(L)+B1(L)*DTAUK +B1(L)*TERM
        CMMID    = B0(L)+B1(L)*DTAUK -B1(L)*TERM
        FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
        FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID

C       FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

        FMIDP(L) = FMIDP(L)*PI
        FMIDM(L) = FMIDM(L)*PI
      END DO
 
      write(*,*) 'amz flag 7'
      flush(6)


C     And now, for the special bottom layer

      L    = L_NLAYRAD

      EP   = EXP(MIN((LAMDA(L)*DTAU(L)),TAUMAX)) ! CLIPPED EXPONENTIAL 
      EM   = 1.0/EP
      TERM = UBARI/(1.-W0(L)*COSBAR(L))

      write(*,*) 'amz flag 8'
      flush(6)


C     CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C     BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

      CPMID    = B0(L)+B1(L)*DTAU(L) +B1(L)*TERM
      CMMID    = B0(L)+B1(L)*DTAU(L) -B1(L)*TERM
      FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
      FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID
 
      write(*,*) 'amz flag 9'
      flush(6)


C     FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

      FMIDP(L) = FMIDP(L)*PI
      FMIDM(L) = FMIDM(L)*PI

C     FLUX AT THE TOP LEVEL

      EP   = 1.0
      EM   = 1.0
      TERM = UBARI/(1.0-W0(1)*COSBAR(1))

C     CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C     BOTTOM OF THE LAYER.  THAT IS AT DTAU  OPTICAL DEPTH

      CPMID  = B0(1)+B1(1)*TERM
      CMMID  = B0(1)-B1(1)*TERM

      FLUXUP = XK1(1)*EP + GAMA(1)*XK2(1)*EM + CPMID
      FLUXDN = XK1(1)*EP*GAMA(1) + XK2(1)*EM + CMMID
 
      write(*,*) 'amz flag 10'
      flush(6)


C     FOR FLUX WE INTEGRATE OVER THE HEMISPHERE TREATING INTENSITY CONSTANT

      FTOPUP = (FLUXUP-FLUXDN)*PI

      write(*,*) 'amz flag 11'
      flush(6)
 

      RETURN
      end 
