C NASA Ames RT code
C module standard_h
C #include STANDARD.h


      _RL TAUCUM(L_LEVELS)
      _RL SIGMA(L_LEVELS)
      _RL PRDST(L_NPDST), TAUDST(L_NPDST)

C     Spatially varying dust - TAUSURF(J,I) is the (visible 0.67 micron)
C     dust opacity at the surface for each grid point
C     TAUTOT added 6/28/01  GCM1.7
C     TAUTOT is the dust opacity at the RPTAU reference pressure
C     for each grid point.

      _RL TAUTOT,RPTAU

!     THE C ARRAY IS A SMALLER RECORD-TYPE STRUCTURE USED TO GROUP
!     TOGETHER A NUMBER OF GLOBALLY USED CONSTANTS.

      _RL RSDIST

      _RL DSIG(L_LAYERS)
      data dsig /3.2718D-07, 5.3944D-07, 8.8940D-07,
     .       1.4663D-06, 2.4175D-06, 3.9859D-06,
     .       6.5718D-06, 1.0835D-05, 1.7863D-05,
     .       2.9452D-05, 4.8556D-05, 8.0060D-05,
     .       1.3199D-04, 2.1763D-04, 3.5879D-04,
     .       5.9155D-04, 9.7530D-04, 1.6080D-03,
     .       2.6512D-03, 4.3714D-03, 7.2060D-03,
     .       1.1882D-02, 1.9590D-02, 3.2298D-02,
     .       5.3245D-02, 8.7800D-02, 1.4475D-01,
     .       2.3865D-01, 3.9347017648D-01 /


      _RL MWCO2
      PARAMETER(MWCO2 = 4.41D-2)
      _RL MWH2O
      PARAMETER(MWH2O = 1.80153D-2)
      _RL MWRATIO

      INTEGER NPDST,LAYERS
      PARAMETER(NPDST     = L_NPDST)
      PARAMETER(LAYERS    = L_LAYERS)

      INTEGER NLAY
      PARAMETER(NLAY    = L_LAYERS)
      _RL PSF,PTROP

C  Standard value = 0.03
C     DATA  CONRNU / 0.03  /   ! Standard value  ~25km half-height
C     DATA  CONRNU / 0.003 /   ! ~50 km half-height
C     DATA  CONRNU / 0.5   /   ! ~10 km half-height

      _RL CONRNU

      common /standard/ TAUCUM,SIGMA,PRDST,TAUDST,TAUTOT,RPTAU,
     &   RSDIST,PSF,PTROP,CONRNU

