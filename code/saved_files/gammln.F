#include "CPP_OPTIONS.h"
 
C     !ROUTINE: GAMMLN
C     !INTERFACE:

      SUBROUTINE GAMMLN(xx,gln)

C     !DESCRIPTION:
C     *==========================================================*
C     | GAMMLN
C     | o From numerical recipes book
C     |  units:
C     |  author: amz
C     *==========================================================*
      
      IMPLICIT NONE
         
C     == Global variables
#include "EEPARAMS.h"
       
C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C
       
      _RL gln,xx
       
C     !LOCAL VARIABLES:
C     == Local variables ==
      _RL ser,stp,tmp,x,y,cof(6)
      INTEGER j
      SAVE cof,stp
  
      DATA cof,stp/76.18009172947146 _d 0,-86.50532032941677 _d 0,
     .  24.01409824083091 _d 0, -1.231739572450155 _d 0,.1208650973866179 _d -2,
     .  -.5395239384953 _d -5, 2.5066282746310005 _d 0/

C     begin equations
      x=xx
      y=x
      tmp=x+5.5 _d 0
      tmp=(x+0.5 _d 0)*LOG(tmp) - tmp
      ser = 1.000000000190015 _d 0
      DO j=1,6
       y=y+1. _d 0
       ser=ser+cof(j)/y
      ENDDO
      gln = tmp + LOG(stp*ser/x)

      RETURN
      END
