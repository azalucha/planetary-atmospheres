#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

      SUBROUTINE INI_SFC_PARMS(myThid)
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE INI_SFC_PARMS
C     | o Initialize surface parameters
C     | o Author: Angela Zalucha
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "AMZVARS.h"

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER i,j,bi,bj,k,myThid

      _RL tfrostsamz,ts,qo,lssol

C      _RL hansen96lat(18),hansen96mass(18),hansen96temp(18)
C      _RL h96mass(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL h96temp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      INTEGER myunit

C Hansen and Paige initial conditions. Leave this here for now

C      IF(radcode .EQ. 6) THEN
C	myunit=25
C	open(myunit,file='hansen96.txt',form='formatted')
C	rewind myunit
C
C	DO i=1,18
C	 read(myunit,*) hansen96lat(i),hansen96mass(i),
C     &                 hansen96temp(i)
C	ENDDO
C      ENDIF

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO J=1,sNy
         DO I=1,sNx
          tautall(I,J,bi,bj)=-1.
          albedoany(I,J,bi,bj)=albedoSs
          massonground(I,J,bi,bj)=0.
          useSurfaceTflag(I,J,bi,bj)=0.
          emissivity(I,J,bi,bj)=emissivityOfSs
          hstaua(I,J,bi,bj)=hstauawarm
          hstaus(I,J,bi,bj)=hstauswarm
          initrad(I,J,bi,bj)=1.

C surface
          IF(radcode .EQ. 6 .OR. radcode .EQ. 10) THEN
           tfrostsamz=1./(1./referenceTempInFrostEqn - atm_Rd/atm_L*
     &            LOG(Ro_Surf(i,j,bi,bj)/referencePresInFrostEqn))
           IF(.NOT. plutoSurfTFixedToFrost) THEN
            IF (lsIsConstant .EQ. 0) THEN
             lssol=lsStart
            ELSE
             lssol=lsperpetual
            ENDIF
            CALL GET_QO(i,j,bi,bj,lssol,0.,solarconstant,
     &                  albedoany(i,j,bi,bj),qo)
            ts=(qo/stephanBoltzmannConstant)**0.25
            IF(ts .LT. tfrostsamz) THEN
             surfaceT(i,j,bi,bj)=tfrostsamz
            ELSE
             surfaceT(i,j,bi,bj)=ts
            ENDIF
           ELSE
            surfaceT(i,j,bi,bj)=tfrostsamz
           ENDIF
          ELSEIF(radcode .EQ. 11 .OR. radcode .EQ. 12) THEN
           surfaceT(I,J,bi,bj)=94.
          ELSEIF(radcode .EQ. 13) THEN
           CALL READ_FLD_XY_RS('INITS',' ',surfaceT,-1,myThid)
          ELSE
           DO K=1,Nr
            surfaceT(I,J,bi,bj)=inisfctemp
           ENDDO
          ENDIF


C subsurface (max depth Nr)
C subsurface scheme only tested for Pluto (eventually Mars)
          IF(useSubSfc) THEN
           DO k=1,Nr
            TsubSfc(i,j,k,bi,bj)=inisubsfctemp
           ENDDO
          ENDIF

C Hansen and Paige initial conditions. Leave this here for now
C          IF(yC(i,j,bi,bj) .GE. -90. .AND.
C     &            yC(i,j,bi,bj) .LE. -80.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(1)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(1)
C          ELSEIF(yC(i,j,bi,bj) .GT. -80. .AND.
C     &            yC(i,j,bi,bj) .LE. -70.) THEN
C          h96mass(i,j,bi,bj)= hansen96mass(2)*rA(i,j,bi,bj)
C          h96temp(i,j,bi,bj)= hansen96temp(2)
C          ELSEIF(yC(i,j,bi,bj) .GT. -70. .AND.
C     &            yC(i,j,bi,bj) .LE. -60.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(3)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(3)
C          ELSEIF(yC(i,j,bi,bj) .GT. -60. .AND.
C     &            yC(i,j,bi,bj) .LE. -50.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(4)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(4)
C          ELSEIF(yC(i,j,bi,bj) .GT. -50. .AND.
C     &            yC(i,j,bi,bj) .LE. -40.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(5)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(5)
C          ELSEIF(yC(i,j,bi,bj) .GT. -40. .AND.
C     &            yC(i,j,bi,bj) .LE. -30.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(6)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(6)
C          ELSEIF(yC(i,j,bi,bj) .GT. -30. .AND.
C     &            yC(i,j,bi,bj) .LE. -20.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(7)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(7)
C          ELSEIF(yC(i,j,bi,bj) .GT. -20. .AND.
C     &            yC(i,j,bi,bj) .LE. -10.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(8)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(8)
C          ELSEIF(yC(i,j,bi,bj) .GT. -10. .AND.
C     &            yC(i,j,bi,bj) .LE. 0.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(9)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(9)
C          ELSEIF(yC(i,j,bi,bj) .GT. 0. .AND.
C     &            yC(i,j,bi,bj) .LE. 10.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(10)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(10)
C           ELSEIF(yC(i,j,bi,bj) .GT. 10. .AND.
C     &            yC(i,j,bi,bj) .LE. 20.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(11)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(11)
C          ELSEIF(yC(i,j,bi,bj) .GT. 20. .AND.
C     &            yC(i,j,bi,bj) .LE. 30.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(12)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(12)
C          ELSEIF(yC(i,j,bi,bj) .GT. 30. .AND.
C     &            yC(i,j,bi,bj) .LE. 40.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(13)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(13)
C          ELSEIF(yC(i,j,bi,bj) .GT. 40. .AND.
C     &            yC(i,j,bi,bj) .LE. 50.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(14)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(14)
C          ELSEIF(yC(i,j,bi,bj) .GT. 50. .AND.
C     &            yC(i,j,bi,bj) .LE. 60.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(15)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(15)
C          ELSEIF(yC(i,j,bi,bj) .GT. 60. .AND.
C     &            yC(i,j,bi,bj) .LE. 70.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(16)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(16)
C          ELSEIF(yC(i,j,bi,bj) .GT. 70. .AND.
C     &            yC(i,j,bi,bj) .LE. 80.) THEN
C           h96mass(i,j,bi,bj)= hansen96mass(17)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(17)
C          ELSE
C           h96mass(i,j,bi,bj)= hansen96mass(18)*rA(i,j,bi,bj)
C           h96temp(i,j,bi,bj)= hansen96temp(18)
C          ENDIF



C Hansen and Paige 96
C
C      CALL WRITE_FLD_XY_RL('H96MASS.',suff,h96mass,
C     &          myiter,myThid)
C      CALL WRITE_FLD_XY_RL('H96TEMP.',suff,h96temp,
C     &          myiter,myThid)

         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END
