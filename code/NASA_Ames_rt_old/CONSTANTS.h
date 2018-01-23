C NASA Ames RT code
C module constants_h
C #include CONSTANTS.h

#include "PARAMS.h"

      _RL SCALEP,Cmk
C     Factor to convert pressures from millibars to Pascals
      PARAMETER(SCALEP = 1.00 _d 2)

C     A radiation code conversion factor.
      PARAMETER(Cmk = 3.51 _d 22)

