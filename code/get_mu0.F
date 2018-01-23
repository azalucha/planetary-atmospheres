#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_MU0
C     !INTERFACE:

      SUBROUTINE GET_MU0(i,j,bi,bj,ls,myTime,mu0)
     
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_MU0
C     | o solar zenith angle as a function of latitude and ls (and longitue
C     |   if in diurnal mode)
C     |  units: W/m^2
C     |  author: amz
C     *==========================================================*
      
      IMPLICIT NONE

C     == Global variables  
#include "SIZE.h"
#include "EEPARAMS.h"   
#include "PARAMS.h"
#include "GRID.h"
#include "AMZVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
      INTEGER i,j,bi,bj

      _RL ls,mu0,myTime

C     !LOCAL VARIABLES:
C     == Local variables ==
      _RL solarDec,haossArg,hourAngleOfSunsetSunrise,gst
      _RL cosz,hourangle

C begin equations

       solarDec = ASIN(SIN(axialTilt)*SIN(ls))
C Case where we do not allow a diurnal cycle.  Ls may be constant or varying.
       IF(allowDiurnal .EQ. 0) THEN
        haossArg = -TAN(yC(i,j,bi,bj)*deg2rad)*TAN(solarDec)
        IF(haossArg .LT. -1.) haossArg = -1.
        IF(haossArg .GT. 1.) haossArg = 1.
        hourAngleOfSunsetSunrise = ACOS(haossArg)
        cosz = (SIN(yC(i,j,bi,bj)
     &       *deg2rad)*SIN(solarDec)*
     &       hourAngleOfSunsetSunrise
     &       + COS(yC(i,j,bi,bj)*deg2rad)*COS(solarDec)*
     &       SIN(hourAngleOfSunsetSunrise))/pi
        IF(cosz .LT. 0.) cosz=0.
        mu0=cosz
       ELSE
        IF(lsIsConstant .EQ. 0) THEN
C Case where there is a diurnal cycle and Ls is varying
C gst="Grenwich" sidereal time (gst), i.e. sidereal time at 0 longitude.  Must be in radians.
         gst = MOD(myTime, lengthOfSidDay)/lengthOfSidDay*2.*PI
         hourangle=gst+xC(i,j,bi,bj)*deg2rad-ls+
     &              TAN(axialTilt/2.)**2.*
     &              SIN(2.*ls)-0.5*TAN(axialTilt/2.)**4.
     &              *SIN(4.*ls)
        ELSE
C Case where there is a diurnal cycle and Ls is constant
	 gst=MOD(myTime,lengthOfDayInSeconds)
     &       /lengthOfDayInSeconds*2.*PI
         hourangle=gst+xC(i,j,bi,bj)*deg2rad
        ENDIF
        cosz=SIN(yC(i,j,bi,bj)*deg2rad)*SIN(solarDec)+
     &        COS(yC(i,j,bi,bj)*deg2rad)*
     &        COS(solarDec)*COS(hourangle)
        IF(cosz .LT. 0.) cosz=0.
        mu0=cosz
       ENDIF
      RETURN            
      END
