#include "CPP_OPTIONS.h"

C     !ROUTINE: GET_LS
C     !INTERFACE:  

      SUBROUTINE GET_LS(myTime,lssol)

C     !DESCRIPTION:  
C     *==========================================================*
C     | SUBROUTINE GET_Ls
C     | o Calculate Ls as a function of time, assuming t=0 is at
C     |   Ls = 0
C     |  units: radians
C     |  author: amz
C     *==========================================================*

      IMPLICIT NONE

C     == Global variables

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "AMZVARS.h"   

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
      _RL myTime,lssol

C     !LOCAL VARIABLES:
C     Local variables 
      _RL ls0,lhs,rhs0,guess0,tol,arg1
      _RL xamz,erroramz,currentfuncvalue,currentfuncder
      _RL myTime1,arg,theCurrentTime,lsStartTime
      INTEGER niters,maxiters

C ***DO NOT CHANGE ls0 *****
       ls0 = 0.0
C     Note: never make ls0-lsp = +/-180 or the equations below break

      lsStartTime=(2.*ATAN((1.-eccentricity)*TAN((lsStart-lsp)/2.)
     &             /SQRT(1.-eccentricity**2.))
     &            /(1.-eccentricity**2.)**1.5 + eccentricity*
     &            SIN(lsStart-lsp)
     &            /(eccentricity**2.-1.)/(1.+eccentricity*
     &            COS(lsStart-lsp))-
     &            2.*ATAN((1.-eccentricity)*TAN((0.-lsp)/2.)
     &            /SQRT(1.-eccentricity**2.))
     &            /(1.-eccentricity**2.)**1.5 - 
     &            eccentricity*SIN(0.-lsp)
     &            /(eccentricity**2.-1.)/(1.+
     &            eccentricity*COS(0.-lsp)))/
     &            SQRT(bigGrav*solarMass/(semimajoraxis*
     &            (1.-eccentricity**2.))**3.)

       IF(lsStartTime .LT. 0.) lsStartTime=lengthofyear+lsStartTime

       theCurrentTime=myTime + lsStartTime

       myTime1=MOD(theCurrentTime,lengthofyear)
       
       rhs0=2.*ATAN((1.-eccentricity)*TAN((ls0-lsp)*0.5)
     &       /SQRT(1.-eccentricity**2.))
     &       /(1.-eccentricity**2.)**(1.5)
     &       +eccentricity*SIN(ls0-lsp)/(eccentricity**2.-1.)
     &       /(1.+eccentricity*COS(ls0-lsp))
       lhs=SQRT(bigGrav*solarMass/(semimajoraxis*
     &           (1.-eccentricity**2.))**3.)*myTime1
       
       guess0 = 0.  
       tol = 0.00001
       
       xamz = guess0
       erroramz = 2.*tol

       niters = 1
       maxiters = 50

       DO WHILE(erroramz .GT. tol .AND. niters .LT. maxiters)
        arg = ATAN((1.-eccentricity)*TAN((xamz-lsp)*0.5)
     &        /SQRT(1.-eccentricity**2.))
        IF(xamz .EQ. MOD(lsp-PI,2.*PI)) THEN
         arg1 = PI/2.
        ELSEIF(xamz .GT. MOD(lsp-PI,2.*PI)) THEN
         arg1 = arg + PI
        ELSE
         arg1 = arg
        ENDIF
     
        currentfuncvalue = 2.*arg1/(1.-eccentricity**2.)**(1.5)
     &       +eccentricity*SIN(xamz-lsp)/(eccentricity**2.-1.)
     &       /(1.+eccentricity*COS(xamz-lsp))-rhs0-lhs
        currentfuncder = 1./(1.+eccentricity*COS(xamz-lsp))
        xamz = xamz - currentfuncvalue/currentfuncder
        erroramz = ABS(currentfuncvalue/currentfuncder)
        niters = niters + 1
       ENDDO

       lssol=xamz

      RETURN
      END
