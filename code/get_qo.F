#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_QO
C     !INTERFACE:

      SUBROUTINE GET_QO(i,j,bi,bj,ls,myTime,
     &                   inpFlux,albedo,qo)
     
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GET_QO
C     | o solar flux as a function of latitude and ls (and longitue
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

      _RL ls,qo,myTime

C     !LOCAL VARIABLES:
C     == Local variables ==
      _RL solarDec,haossArg,hourAngleOfSunsetSunrise,gst
      _RL cosz,hourangle,inpFlux,albedo

C begin equations

       solarDec = ASIN(SIN(axialTilt)*SIN(ls))
C Case where we do not allow a diurnal cycle.  Ls may be constant or varying.
       IF(allowDiurnal .EQ. 0) THEN
        haossArg = -TAN(yC(i,j,bi,bj)*deg2rad)*TAN(solarDec)
        IF(haossArg .LT. -1.) haossArg = -1.
        IF(haossArg .GT. 1.) haossArg = 1.
        hourAngleOfSunsetSunrise = ACOS(haossArg)
        qo = inpFlux/PI*(1. - albedo)*
     &       ((1. + eccentricity*COS(ls - lsp))
     &       /(1. - eccentricity**2.))**2.*(SIN(yC(i,j,bi,bj)
     &       *deg2rad)*SIN(solarDec)*
     &       hourAngleOfSunsetSunrise
     &       + COS(yC(i,j,bi,bj)*deg2rad)*COS(solarDec)*
     &       SIN(hourAngleOfSunsetSunrise))

       ELSE
        IF(lsIsConstant .EQ. 0) THEN
C Case where there is a diurnal cycle and Ls is varying
C gst="Grenwich" sidereal time (gst), i.e. sidereal time at 0 longitude.  Must be in radians.
         gst = MOD(myTime, lengthOfSidDay)/lengthOfSidDay*2.*PI
         hourangle=gst+xC(i,j,bi,bj)*deg2rad-ls+
     &              TAN(axialTilt/2.)**2.*
     &              SIN(2.*ls)-0.5*TAN(axialTilt/2.)**4.
     &              *SIN(4.*ls)
        ELSEIF(lsIsConstant .EQ. 1) THEN
C Case where there is a diurnal cycle and Ls is constant
	 gst=MOD(myTime,lengthOfDayInSeconds)
     &       /lengthOfDayInSeconds*2.*PI
         hourangle=gst+xC(i,j,bi,bj)*deg2rad
        ELSE
C Case where we have synchronus rotation with star
         hourangle=xC(i,j,bi,bj)*deg2rad
        ENDIF
        cosz=SIN(yC(i,j,bi,bj)*deg2rad)*SIN(solarDec)+
     &        COS(yC(i,j,bi,bj)*deg2rad)*
     &        COS(solarDec)*COS(hourangle)
        IF(cosz .LT. 0.) cosz=0.
        qo=inpFlux*(1. - albedo)*
     &        ((1. + eccentricity*COS(ls - lsp))
     &       /(1. - eccentricity**2.))**2.*cosz
       ENDIF
      RETURN            
      END
