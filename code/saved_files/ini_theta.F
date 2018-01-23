C $Header: /u/gcmpack/MITgcm/verification/hs94.cs-32x32x5/code/ini_theta.F,v 1.8 2009/09/27 23:49:28 jmc Exp $
C $Name:  $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: INI_THETA
C     !INTERFACE:
      SUBROUTINE INI_THETA( myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE INI_THETA
C     | o Set model initial temperature field.
C     *==========================================================*
C     | There are several options for setting the initial
C     | temperature file
C     |  1. Inline code
C     |  2. Vertical profile ( uniform T in X and Y )
C     |  3. Three-dimensional data from a file. For example from
C     |     Levitus or from a checkpoint file from a previous
C     |     integration.
C     | In addition to setting the temperature field we also
C     | set the initial temperature tendency term here.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "AMZVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid -  Number of this instance of INI_THETA
      INTEGER myThid

C     == Functions ==
c     real*8  PORT_RAND
c     real*8  seed

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER bi, bj, myIter
      INTEGER I, J, K, localWarnings
      _RL     term1,term2,thetaLim,thetaEq
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(MAX_LEN_MBUF) suff

C     amz variables
CEOP
      _RL plutoT(Nr)
      _RL initempA(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy) 
      _RL titanT(Nr)


      myIter=0
          WRITE(suff,'(I10.10)') myIter
      

      J = 99+myBxLo(myThid)+nPx*myByLo(myThid)
c     CALL SRAND( J )
c     seed = j

C     =======AMZ radcode check==========
      IF(radcode .EQ. 1 .OR. radcode .EQ. 3
     &   .OR. radcode .EQ. 5. .OR. radcode .EQ. 7) THEN
       write(*,*) 'Error: radcode not supported'
       STOP
      ENDIF
C     ========volatile (frost cycle) checks======
      IF(radcode .EQ. 10 .OR. radcode .EQ. 6) THEN
       IF(selectAddFluid .EQ. 1) THEN 
        IF(plutoSurfTFixedToFrost) THEN
         write(*,*) 'Error: cannot have active volatile cycle 
     &        with Pluto (Triton) surface temperature fixed'
         STOP
        ENDIF
       ENDIF
      ENDIF
      IF(radcode .EQ. 11 .OR. radcode .EQ. 12) THEN
       IF(selectAddFluid .EQ. 1) THEN
        write(*,*) 'Error: volatile cycle not supported for Titan'
        STOP
       ENDIF
      ENDIF
      IF(useSubSfc) THEN
        write(*,*) 'Error: subsurface scheme not developed yet'
        STOP
      ENDIF



C     ========initialize some parameters here=======
C Initialize scalars
      lengthofyear=SQRT(4.*PI**2./(bigGrav*solarMass)
     &             *semimajoraxis**3.)
      lengthOfSidDay=lengthOfDayInSeconds/
     &      (1.+lengthOfDayInSeconds/lengthofyear)
C Initialize surface and subsurface parameters 
C (includes temp, albedo, etc.)
      CALL INI_SFC_PARMS(myThid)
C Initialize atmospheric temperature
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO J=1,sNy
         DO I=1,sNx
          IF(radcode .EQ. 6) THEN
            CALL INI_PLUTO_TEMP(i,j,bi,bj,plutoT)
           DO k=1,Nr
            initempA(I,J,K,bi,bj)=plutoT(K)
           ENDDO
          ELSEIF(radcode .EQ. 10) THEN
           CALL INI_STROBEL_TEMP(i,j,bi,bj,plutoT)                             
           DO k=1,Nr                                   
            initempA(I,J,K,bi,bj)=plutoT(K)            
           ENDDO                                       
          ELSEIF(radcode .EQ. 11 .OR. radcode .EQ. 12) THEN
           CALL INI_TITAN_TEMP(i,j,bi,bj,titanT,surfaceT(i,j,bi,bj))
           DO k=1,Nr
            initempA(I,J,K,bi,bj)=titanT(K)
           ENDDO
          ELSE
           DO K=1,Nr
            initempA(I,J,K,bi,bj)=iniatmtemp
           ENDDO
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF ( hydrogThetaFile .EQ. ' ' ) THEN

       IF (radcode .EQ. 13) THEN
        CALL READ_FLD_XYZ_RL('INITT',' ',thetaEqSave,-1,myThid)
       ENDIF

       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO K=1,Nr
          DO J=1,sNy
           DO I=1,sNx
            theta(I,J,K,bi,bj)=initempA(I,J,K,bi,bj)
     &          /((rC(K)/atm_po)**atm_kappa)
           ENDDO
          ENDDO
         ENDDO
#ifdef ALLOW_ZONAL_FILT
C--   Zonal FFT filter initial conditions
         IF (useZONAL_FILT) THEN
          CALL ZONAL_FILTER(
     U                       theta(1-OLx,1-OLy,1,bi,bj),
     I                       hFacC(1-OLx,1-OLy,1,bi,bj),
     I                       1, sNy, Nr, bi, bj, 1, myThid )
         ENDIF
#endif /* ALLOW_ZONAL_FILT */
        ENDDO
       ENDDO
      ELSE
       CALL READ_FLD_XYZ_RL( hydrogThetaFile, ' ', theta, 0, myThid )
      ENDIF

      IF(radcode .EQ. 6 .OR.
     &   radcode .EQ. 11 .OR. radcode .EQ. 12
     &   .OR. radcode .EQ. 9 .OR. radcode .EQ. 10) THEN
       CALL SNAP_T_TO_TFROST(myThid)
      ENDIF

C--   Apply mask and test consistency
      localWarnings=0
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO K=1,Nr
         DO J=1-Oly,sNy+Oly
          DO I=1-Olx,sNx+Olx
           IF (maskC(I,J,K,bi,bj).EQ.0.) theta(I,J,K,bi,bj) = 0.
          ENDDO
         ENDDO
         IF ( tRef(k).NE.0. ) THEN
          DO J=1,sNy
           DO I=1,sNx
            IF (  maskC(I,J,K,bi,bj).NE.0.
     &      .AND. theta(I,J,K,bi,bj).EQ.0. ) THEN
              localWarnings=localWarnings+1
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      IF (localWarnings.NE.0) THEN
       WRITE(msgBuf,'(A,A)')
     &  'S/R INI_THETA: theta = 0 identically. If this is intentional',
     &  'you will need to edit ini_theta.F to avoid this safety check'
       CALL PRINT_ERROR( msgBuf , myThid)
       STOP 'ABNORMAL END: S/R INI_THETA'
      ENDIF

      _EXCH_XYZ_RL(theta , myThid )

      IF (debugMode) THEN
        CALL PLOT_FIELD_XYZRL( theta, 'Initial Temperature' ,
     &                         Nr, 1, myThid )
      ENDIF



      RETURN
      END

