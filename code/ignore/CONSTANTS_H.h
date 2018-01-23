C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C                            constants_h
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



C     Acceleration due to gravity (mks)
      _RL GRAV
      PARAMETER(GRAV=3.72)

C     Gas constant for mars.
      _RL RGAS
      PARAMETER(RGAS=189.02)

C     Stefan-Boltzmann constant
      _RL STBO
      PARAMETER(STBO=5.67E-8)

C     Heat capacity (or specific heat) of CO2 gas.
C     ( units of joules per ( kg * degrees kelvin))
      _RL Cp
      PARAMETER(Cp=860.)

C     Factor to convert pressures from millibars to Pascals
      _RL SCALEP
      PARAMETER(SCALEP = 1.00E2)

C     A radiation code conversion factor.
      _RL Cmk
      PARAMETER(Cmk = 3.51E22)

