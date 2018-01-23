C     !ROUTINE: AMZVARS.h
C     !INTERFACE:
C     #include AMZVARS.h

C     !DESCRIPTION:
C     Header file defining amz's variables not present in the
C     full MIT GCM


      common /amzvars/ teqamz,atm_L,lengthOfDayInSeconds,
     &  albedoSs,albedoFrost,massonground,albedoany,
     &  tautall,etaNmean,lsperpetual,
     &  eccentricity,dTdtsave,thetaEqSave,subsfccond,
     &  semimajoraxis,lsp,lengthofyear,axialTilt,
     &  solarconstant,tauoo,radTimeStep,bigGrav,solarMass,
     &  stephanBoltzmannConstant,lsStart,hstaua,hstaus,surfaceT,
     &  Tsubsfc,cond,cond1,q1,l1,
     &  emissivity,emissivityOfSs,emissivityOfFrost,
     &  lengthOfSidDay,plutoMass,condsurf,
     &  betacab,tauinfcab,rhoOfFrost,q2,l2,cond2,
     &  hstauawarm,hstauswarm,useSurfaceTflag,cpOfFrost,
     &  hstauacold,hstauscold,etaNmax,etaNmin,
     &  dumpFreqFastStart,dumpFreqFast,tdownsurf,
     &  referenceTempInFrostEqn,referencePresInFrostEqn,
     &  topoTextFile,iniatmtemp,ssThickness,taus0,
     &  taudustir0,dtsnapdt,QT76SAVE,HT23SAVE,
     &  HT33SAVE,QTCOSAVE,CONDSAVE,QOSAVE,
     &  cpOfSs,rhoOfSs,iMinus0,qr,initrad,hc,
     &  nlosch,gamma0,kb,clight,amu,hplanck,piF0,comr,
     &  thetabefore,thetaafter,gtadvectsave,
     &  supersfccond,inisfctemp,inisubsfctemp,
     &  radcode,j_planet,lsIsConstant,allowDiurnal,
     &  plutoSurfTFixedToFrost,iniPlutoDiskAvg,readtopofromfile,
     &  calcGradWind,initradread,h96,
     &  plutoTrop,marsPureRad,useSubSfc,switchToPlutoVolatile

      _RL teqamz(sNx,sNy,Nr,nSx),solarconstant,tauoo
      _RL atm_L,lengthOfDayInSeconds,bigGrav,solarMass
      _RL massonground(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL albedoany(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL radTimeStep,cpOfFrost
      _RL surfaceT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL albedoSs,albedoFrost,lsp,lengthofyear
      _RL tautall(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etaNmean,lsperpetual,eccentricity,semimajoraxis
      _RL stephanBoltzmannConstant,lsStart,axialTilt
      _RL emissivity(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL emissivityOfSs,emissivityOfFrost,lengthOfSidDay
      _RL hstaua(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hstaus(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL hstauawarm,hstauswarm,hstauacold,hstauscold
      _RL etaNmax,etaNmin,dumpFreqFast,dumpFreqFastStart
      _RL useSurfaceTflag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL betacab,tauinfcab,iniatmtemp,inisfctemp,inisubsfctemp
      _RL referenceTempInFrostEqn,referencePresInFrostEqn
      _RL cpOfSs,rhoOfSs,ssThickness,iMinus0
      _RL qr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL taus0,rhoOfFrost,taudustir0,piF0,hc,comr
      _RL initrad(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tdownsurf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL nlosch,gamma0,kb,clight,amu,hplanck,plutoMass
      _RL dtsnapdt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL cond(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL cond1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL q1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL cond2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL q2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL QOSAVE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL QT76SAVE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HT23SAVE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HT33SAVE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL QTCOSAVE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL CONDSAVE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL dTdtsave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL thetaEqSave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL thetabefore(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL thetaafter(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL gtadvectsave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL condsurf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Tsubsfc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sfccond(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL subsfccond,supersfccond
      INTEGER radcode,lsIsConstant,allowDiurnal
      INTEGER j_planet

      CHARACTER*(MAX_LEN_FNAM) topoTextFile

      LOGICAL readtopofromfile,plutoSurfTFixedToFrost
      LOGICAL iniPlutoDiskAvg,calcGradWind,useSubSfc
      LOGICAL initradread,h96,plutoTrop,marsPureRad
      LOGICAL switchToPlutoVolatile

C teqamz : Radiative equilibrium temperature, if read in from file (K)

C atm_L : latent heat of sublimation (J/kg)

C lengthOfDayInSeconds : length of solar day (s)

C massonground : Mass/time of accumulated frost on the ground (same units as addMass, kg/s)

C albedoany : surface albedo (dimensionless)

C albedoSs : constant surface albedo for frost-free regions (dimensionless).

C albedoFrost : surface albedo for frost that has accumulated on surface
C (and subject to volatile cycle) (dimensionless).

C tautall : optical depth of the top of the convective layer for the
C radiative-convective equilibrium.  Note:
C tautall may be in either AMZ's gray scheme specification of
C tau or Caballero et al. 2008's semi-gray scheme.
C o/w, tau and tauo always refer to AMZ's specification,
C taucab (things with "cab" in them) are semi-gray scheme. (dimensionless)

C etaNmean : global mean eta (eta units)

C radcode :
C    1. No longer used.
C    2. radiative-convective equilibrium, Newtonian relaxation
C       Calculates instantaneous surface temperature.
C	(Mars)
C    3. No longer used.
C    4. Caballero et al. (2008) CO2 semigray scheme with heating
C       rate supplied directly.  Convective adjustment. Time-
C       varying surface temperature. (Mars)
C    5. No longer used.
C    6. Yelle & Lunine (1989).  Heating rate supplied directly.
C       Non-LTE CH4 heating at 3.3 um and cooling at 7.6 um.
C       No 2.3 um. Time-varying surface temperature. (Pluto)
C    7. No longer used.
C    8. Caballero et al. (2008) CO2 semigray scheme, Newtonian
C       relaxation. Calculates instantaneous surface temperature. (Mars)
C    9. NASA Ames radiation scheme. Calculates time-varying surface temperature. (Mars)
C   10. Strobel & Zhu (2017) radiation scheme. Calculates time-varying surface temperature. (Pluto)
C   11. Robinson & Catling (2012) scheme.
C       Calculates time-varying surface temperature. (Titan)
C   12. Robinson & Catling (2012) scheme, with Newtonian relaxation.
C       Calculates instantaneous surface temperature. (Titan)
C   13. Newtonian relaxation to Pluto radiative-conductive equilibrium.  Used as an
C       intermediate step in spinning up GCM with radcodes 6 or 10.  Not
C       meant as a stand alone experiment. (Pluto)

C lsIsConstant : 1 for for perpetual season; 0 for changing season; -1 for synchronus rotation

C eccentricity : orbital eccentricity (dimensionless)

C lsp : Ls of perihelion (radians).  Warning: never set equal to 180.

C lengthofyear : length of year (s)

C axialTilt : axialTilt (radians).  Warning: is used in radiation code only.  Be
C careful in how you set omega and axialTilt in tricky bodies like Pluto that rotate backwards.

C solarconstant : mean solar constant (W/m^2)

C tauoo : Optical depth at reference pressure, atm_po (dimensionless)

C radTimeStep: Radiation time stepping for RC equilibrium (seconds).

C bigGrav: universal gravitational constant (m^3/kg/s^2)

C solarMass: mass of the sun in kg

C stephanBoltzmannConstant : Stephan-Boltzmann constant (W/m^2/K^4)

C lsStart : value of Ls at initial time (radians), if Ls is not constant.
C Otherwise, set lsperpetual.

C hstaua: Held-Suarez 1994 1/ka (seconds), as a function of lat and lon

C hstaus: Held-Suarez 1994 1/ks (seconds), as a function of lat and lon

C surfaceT: surface temperature (not potential tempearture)

C useSurfaceTflag: used when there is frost on the ground (radcode 2 and 8 only).
C If 1, tells radiation code to calculate Teq
C from surface flux; If 0, tells radiation code to calulate Teq from solar flux only.

C emissivityOfSs: emissivity of Subsurface (dimensionless)

C emissivityOfFrost: emissivity of accumulated frost (dimensionless)

C lengthOfSidDay: length of sidereal day (s).

C allowDiurnal: 1 for diurnally varying solar radiation, 0 for daily averaged solar radiation.

C hstauawarm,hstauacold: value for hstaua if the atmosphere is
C warm (i.e. not at frost temperature)
C or cold (i.e. at frost temperature).  Only applies if selectAddFluid=1.  (seconds)

C hstauswarm,hstauscold: same as hstauawarm,hstauacold, but for hstaus.

C etaNmax, etaNmin: maximum and minimum values for etaN.

C dumpFreqFast, dumpFreqFastStart: for when we want to reset
C the dump frequency to higher time
C resolution midway through the run.  dumpFreqFast sets this
C new frequency (in s), dumpFreqFastStart
C says what time this occurs (in s).

C betacab: beta for Caballero et al. 2008 semi-gray radiation scheme

C tauinfcab: tau_infinity for Caballero et al. 2008 semi-gray radiation scheme

C referenceTempInFrostEqn,referencePresInFrostEqn: reference
C temperature and pressure in frost equation (K, Pa).

C topoTextFile: I don't get how to make binary files, so
CI make text files of topography and read them in.  Make sure
C you also read in a dummy binary file topo.cs.bin

C cpOfSs: specific heat of subsurface material (J/K/kg)

C rhoOfSs: density of subsurface matrial (kg/m^3)

C ssThickness: thickness of surface slab (m)

C initrad: flag to call radiation code on first time step, do not edit

C tdownsurf: transmission of solar beam at the surface

C rhoOfFrost: density of surface frost (kg/m^3)

C taudustir0: IR optical depth of dust, given as a fraction of the solar optical depth of dust

C nlosch: Loschmidt number (m^-3)

C gamma0: mixing ratio of CH4

C kb: Boltzmann constant (J/K)

C clight: speed of light (m/s)

C amu: atmomic mass unit (kg)

C hplanck: Planck's constant (J/s)

C hc: troposhphere critical height (plutoTrop=.TRUE. only) (km)

C piF0: solar constant at 3.3 um (radcode 6 only) (W/m^2/s/um)

C plutoMass: Pluto or Triton's mass (kg)

C readtopofromfile: true, read topography from topo_Hin* type file.  False,
C use user-specified analytical expression ini_depths.F.  Note with this
C flag we have overwritten the default MITGCM configuration for reading in
C topography, so that topoFile is ignored.

C plutoSurfTFixedToFrost: true, Pluto's (Triton's) surface temperature is fixed to
C frost temperaure.  False: temperature is alowed to vary (ignores diffusive IR
C flux from atmosphere).  Only called for radcode 6.  Assumes ground is made
C of Mars-like materials, for now.

C comr: CO mixing ratio (Pluto) (dimensionless)

C dtsnapdt: temperature tendency due to snapping temperature to frost temperature (K/s).

C cond: conduction term (Pluto) (J/s/m^3)

C dTdtsave: qr with potential temperature correction applied (K)

C thetaEqSave: equilibrium potential temperature, radcode 2, 8, 12, and 13 only (K)

C cond1: conduction term, in K/s (radcode 6 only)

C q1: 3.3um heating term, in K/s (radcode 6 only)

c l1: 7.6um cooling term, in K/s (radcode 6 only)

C iniatmtemp: for certain radcodes, sets atmospheric temperature to this value 
C  everywhere (K)

C inisfctemp: for certain radcodes/conditions, sets surface temperature to
C  this value everywhere (K)

C inisubsfctemp: if using subsurface scheme, sets entire subsurface to this
C  value everywhere (K)

C plutoTrop: true includes troposphere for pluto schemes (Yelle only, strobel untested)

C useSubSfc: true uses the subsurface scheme (certain radcodes)

C QT76SAVE: CH4 cooling at 7.6 um, in K/s (radcode 10 only)

C HT23SAVE: CH4 heating at 2.3 um, in K/s (radcode 10 only)

C HT33SAVE: CH4 heating at 3.3 um, in K/s (radcode 10 only)

C QTCOSAVE: CO cooling, in K/s (radcode 10 only)

C CONDSAVE: Atmospheric conduction, in K/s (radcode 10 only)

C switchToPlutoVolatile: short circuits read_pickup when going from Pluto spin up stages (no volatile cycle) to main part of run (with volatile cycle)

C j_planet: Used only with radcodes 10 and 13.  Pluto = 1, Triton = 2. 

C iniPlutoDiskAvg: not used, currently set to FALSE.
