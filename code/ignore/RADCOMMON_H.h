C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                  
C                            radcommon_h
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C----------------------------------------------------------------------C
C
C                             radcommon.h
C                         FORTRAN PARAMETERS
C                          GCM2.0  Feb 2003
C                  
C----------------------------------------------------------------------C
C  
C     WNOI	 - Array of wavenumbers at the spectral interval
C                  centers for the infrared.  Array is NSPECTI
C                  elements long.
C     DWNI       - Array of "delta wavenumber", i.e., the width,
C                  in wavenumbers (cm^-1) of each IR spectral
C                  interval.  NSPECTI elements long.
C     WAVEI	 - Array (NSPECTI elements long) of the wavelenght
C                  (in microns) at the center of each IR spectral
C                  interval.
C     WNOV	 - Array of wavenumbers at the spectral interval
C                  center for the VISIBLE.  Array is NSPECTV
C                  elements long.
C     DWNV	 - Array of "delta wavenumber", i.e., the width,
C                  in wavenumbers (cm^-1) of each VISIBLE spectral
C                  interval.  NSPECTV elements long.
C     WAVEV      - Array (NSPECTV elements long) of the wavelenght
C                  (in microns) at the center of each VISIBLE spectral
C                  interval.
C     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
C                  each spectral interval.  Values are for 1 AU, and
C                  are scaled to the Mars distance elsewhere.
C     TAURAY     - Array (NSPECTV elements) of the pressure-independent
C                  part of Rayleigh scattering optical depth.
C     PTOP	 - Pressure at the top of the radiation code coordinate;
C                  = smallest k-coefficient pressure (1.0E-6 mbar)
C     FZEROI     - Fraction of zeros in the IR CO2 k-coefficients, for
C         	   each temperature, pressure, and spectral interval
C     FZEROV     - Fraction of zeros in the VISIBLE CO2 k-coefficients, for
C                  each temperature, pressure, and spectral interval
C     
C     AEROSOL RADIATIVE OPTICAL CONSTANTS
C     Values are at the wavelenght interval center
C     
C     MIE SCATTERING - Size distribution weighted
C     Qextv    - Extinction efficiency - in the visible.
C     QextREF  - Reference visible wavelength (.67 micron band)
C     Qscatv   - Scattering efficiency - in the visible.
C     WV       - Single scattering albedo - in the visible.
C     GV0       - Asymmetry parameter - in the visible.
C
C     Qexti    - Extinction efficiency - in the infrared.
C     Qscati   - Scattering efficiency - in the infrared.
C     WI       - Single scattering albedo - in the infrared.
C     GI       - Asymmetry parameter - in the infrared.
C     
C  Water ice clouds
C     
C     MIE SCATTERING - Size distribution weighted
C     Qextvc   - Extinction efficiency - in the visible.
C     QextREFc - Reference visible wavelength (.67 micron band)
C     Qscatvc  - Scattering efficiency - in the visible.
C     WVc      - Single scattering albedo - in the visible.
C     GVc      - Asymmetry parameter - in the visible.
C
C     Qextic   - Extinction efficiency - in the infrared.
C     Qscatic  - Scattering efficiency - in the infrared.
C     WIc      - Single scattering albedo - in the infrared.
C     GIc      - Asymmetry parameter - in the infrared.
C

      _RL WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI)
      _RL WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
      _RL SOLARF(L_NSPECTV), TAURAY(L_NSPECTV)

      _RL CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      _RL CO2V(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)
      _RL FZEROI(L_NSPECTI)
      _RL FZEROV(L_NSPECTV)
      _RL PGASREF(L_NPREF), TGASREF(L_NTREF)

      _RL qextv(L_NSPECTV), qscatv(L_NSPECTV), wv(L_NSPECTV)
      _RL gv0(L_NSPECTV)
      _RL QextREF

      _RL qexti(L_NSPECTI), qscati(L_NSPECTI), wi(L_NSPECTI)
      _RL gi(L_NSPECTI)

      _RL qextvc(L_NSPECTV), qscatvc(L_NSPECTV), wvc(L_NSPECTV)
      _RL gvc(L_NSPECTV)
      _RL QextREFc
      _RL qextic(L_NSPECTI), qscatic(L_NSPECTI), wic(L_NSPECTI)
      _RL gic(L_NSPECTI)

      _RL planckir(L_NSPECTI,8501)

      _RL TAUREF(L_LEVELS+1)
      _RL PFGASREF(L_PINT)

C     Hemispheric mean in the IR
      _RL UBARI
      PARAMETER(UBARI=0.5)

C     These are for the Gauss-split 0.95 case

      _RL GWEIGHT(L_NGAUSS)
      data GWEIGHT /   4.8083554740E-02, 1.0563099137E-01,
     .                 1.4901065679E-01, 1.7227479710E-01,
     .                 1.7227479710E-01, 1.4901065679E-01,
     .                 1.0563099137E-01, 4.8083554740E-02,
     .                 2.5307134073E-03, 5.5595258613E-03,
     .                 7.8426661469E-03, 9.0670945845E-03,
     .                 9.0670945845E-03, 7.8426661469E-03,
     .                 5.5595258613E-03, 2.5307134073E-03,
     .                 0.0E0 /

C     These are for the CO2+H2O k-coefficients

      _RL WREFCO2(L_REFH2O)
      data WREFCO2 / 9.999999E-1, 9.99999E-1, 9.9999E-1, 9.999E-1,
     .            9.99E-1, 9.9E-1, 9.0E-1, 8.0E-1, 7.0E-1, 6.0E-1/

      _RL WREFH2O(L_REFH2O)
      data WREFH2O /1.0E-7, 1.0E-6, 1.0E-5, 1.0E-4, 1.0E-3, 1.0E-2,
     .                1.0E-1, 2.0E-1, 3.0E-1, 4.0E-1/

C     If the CO2 optical depth (top to the surface) is less than
C     this value, we place that Gauss-point into the "zeros"
C     channel.  NRC parameter.

C  TLIMITS - TLIMIT for solar part of the spectrum
C  TLIMITI - TLIMIT for the IR

      _RL TLIMITS
      PARAMETER(TLIMITS = 1.0E-32)

