#include "CPP_OPTIONS.h"

C     !ROUTINE: GAMMQ
C     !INTERFACE:

      SUBROUTINE GAMMQ(a,x,out)

C     !DESCRIPTION:
C     *==========================================================*
C     | GAMMQ
C     | o incomplete gamma function, gamma(a,x).  Note: output differs from 
C     | numerical recipes book in that we have already multiplied Q 
C     | (what the function with this name is supposed to return) by
C     | gamma(a).
C     |  units:
C     |  author: amz
C     *==========================================================*

      IMPLICIT NONE

C     == Global variables
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     

      _RL a,x,out

C     !LOCAL VARIABLES:
C     == Local variables ==

      _RL gammcf,gamser,gln

C     begin equations

      IF(x .LT. 0. .OR. a .LT. 0.) write(*,*) 'bad arguments in gammq'
      IF(x .LT. a+1.) THEN
       CALL GSER(gamser,a,x,gln)
       out=(1.-gamser)*EXP(gln)
      ELSE
       CALL GCF(gammcf,a,x,gln)
       out=gammcf*EXP(gln)
      ENDIF

      RETURN
      END
