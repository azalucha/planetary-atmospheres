#include "CPP_OPTIONS.h"

      SUBROUTINE INI_STROBEL_TEMP(ii,j,bi,bj,plutoT,myThid)
C     *==========================================================*
C     | S/R STROBELRAD
C     | o Calculate radiative equilibrium temperature to initialize
C     |   Pluto temperatures using Strobel and Zhu (2017) and
C     |   Zalucha et al. (2011) troposphere specification
C     *==========================================================*

C Notes from Xin Zhu:
CCC   Radiative equilibrium temperature for Pluto/Triton's atmosphere. In MKS units.
C
C     pluto_triton_t,        qcgs_sum,  qcgs_sum_j,
C     f_cap_j,   dtdz2,      ZTPR,      ZTPRN_w,  
C     zt_w_parker            IMDIFF3,   IMDIFF2,  IMDIFF1,  IMDjFF2,  
C     tridag,    HTCH4tr1,   HTCH4tr2,  htch4nt_r,  
C     HTCH5t,    htch5nt_r,             INPUT4,   INPUT5, 
C     QTCH4Z,    QJCH4X,     GAMM12,    INPUT6,   str_t_int,
C     zkg_tp_int,            QTCOZ,     QTHCNZ,   QTH2OZ,     
c     QJ_line_co,            QJ_line_hcn,         QJ_line_h2o,  
c     alpd_co,   alpd_hcn,   alpd_h2o, 
c     strt7_co,  strt7_hcn,  strt7_h2o,    
c     ENZD2,     BLACM,      POLINT,    INVERT,   
c     gamma_fels,            gamma_zhu, gamma_zhu_g,        voigt_fs85,
c     erfz1,     erfcz2
c 
CCCCC 
c
c Parameters need to be manually adjusted/changed in order to run the code efficiently:
c
c (1) dt=50.0d0*DAYL     !  50.0d0*DAYL    ! integration time interval in seconds
c
c (2) KMX=175,ITIME=98   ! (vertical levels numbers, time steps)
c     KMX needs to be set to have the top boundary a few grids above the exobase where eta=1
c
c (3) xfac22=1.00d0 !=1.0 ~ [0.95,1.05] to eliminate the oscillating exobase in time integration)
c
c (4) jeuvfuv = (0,1,2,3) = (use f107 parameterization, case 1, case 2, case 3)
c        =(1,2*,3)=(Triton_2009_0704, Pluto_2014paper*, Triton_1989_0808) with f107=(70,130,220)
c     f107=130.0d0   ! used when jeuvfuv=0
c
c (5) j_fels=(1,2,3) = select gamma_fels, gamma_zhu_g, or gamma_zhu for accuracy and speeding consideration.
c
c
C  To calculate rediative equilibrium temperature in Pluto/Triton's atmosphere. In MKS units.
c
C end notes from Xin Zhu

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"


      INTEGER  bi,bj,k,i,j,nlevs,pindex,itime,km,myThid
      INTEGER pinterpindexstart,skipinterp,kmx,kmm,ii
      _RL plutoT(Nr)

      PARAMETER(KMX=175,ITIME=3600,km=kmx,kmm=km-1)    ! (vertical levels numbers, time steps)
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      _RL ZM(KM),GRV(KM),TEM(KM),TEMM(KM),TEMX(KM),rcpr,rrcp,xf2r
     &    ,PRE(KM),RHO(KM),DN(KM),BNS(KM),HSCAL(KM),wtmol,rair                  ! Basic state
      _RL RHO1(KM),RHO2(KM),RHO3(KM),rho4(km),rho5(km),cvair,g00
     & ,DN1(KM),BNS1(KM),FXTOP,DELZ0,FAC1,hc_ch4,rh_ch4,rch4st,xf1r
      _RL RM(KM),ZKM(KM),RKM(KM),eta_z(km),xf3r,xfac1r,xfac2r,Rh2o1
      _RL HT1(KM),HALL(KM),QT1(KM),QT2(KM),QT3(KM),QT4(KM),Rh2o2
     & ,QT5(KM),qt7(km),QALL(KM),QNETAMZ(KM),Rh2o3,Rh2o4,Rh2o5,Rh2o6,
     & qnet_r(km),hdiff(km),ht_mag(km),xf4r
      _RL HT2(KM),HT2CGS(KM),QDIFF(KM),QDFCGS(KM),xf5r,xf6r,xfac1
      _RL QS(KM),DS(KM),CS(KM),bs(km),TALP(KM),TALPM(KM),GAMN4(KM)
      _RL dzz1(km),dzz2(km)  ! background and convective eddy diffusion coefficient
     & ,wrk1(km),wrk2(km),wrk3(km),wrk4(km),wrk5(km),wrk6(km),xfac2
     & ,wrk7(km),wrk8(km),wsk1(km),wsk2(km),wsk3(km),xfact,cgamma
      _RL ae1(km),weigt1(km),w_in(km),w_out(km),zx1,tk1,rsunp0
c
      INTEGER nm3,j_fels,jeuvfuv,j1_efuv,jj_hcn_sat,ix1,jlowbd
      parameter(nm3=3)
      _RL ht3(km),ht_ch4(km,nm3),ht_n2(km,nm3),epstcr,r_gravity
      _RL ht4(km),ht_mag_cgs0(km),qt7cgs(km),f107,f1_mag,mu0,dt00
c
      _RL tsecd(itime),tyear(itime),rsunp(itime),rnormf(itime)
      _RL tsur_t(itime),p0_t(itime),tem_t(km,itime),TSUR0,P00
     & ,pre_t(km,itime),rho_t(km,itime),dn_t(km,itime),qo,capm
      _RL tem0(km),rch4x0,RCO,RHCN1,RHCN2,RHCN3,r_c2h2_0,r00
      _RL frho_t(itime),xx4(4),yy4(4),rsunup_fix,surface_albedo
      INTEGER jeans_t(itime),jcrit,IYELLE,MQCO,MQHCN,MQh2o,MQ_c2h2
      _RL zm_ly(kmm),pre_ly(kmm),dpre_ly(kmm),tem_ly(kmm),bk,capg
      _RL lssol,DAYL,factx1,cpair,dt,ZKVT,XA3,XA4,XA2,XP2,XP3,XP4
      _RL frho,TINF1,frho_in0,XKM,RHOK,AKK0,ALPH,ALPH1,ALPHI,xfac22
      _RL fcapj3,hscal3,eta3,xlambda3,capu3,vesc3,capf_rho_g3,frho3
      _RL fcapj4,hscal4,eta4,PREK,xlambda4,capu4,vesc4,capf_rho_g4
      _RL frho4,fcapj2,hscal2,xlambda2,capu2,vesc,capf_rho_g2
      _RL fcapj5,hscal5,xlambda5,capu5,vesc5,capf_rho_g5,DY11,BN(KM)
      _RL pa_hcn_sat,f_super_sat,rho3_comp,rho3_low,ww_high,ww_low
      _RL rho3_high,vesc2,frho_in,FAC11,xfac11,efficy3,TALP1,RADS1
      _RL ae1c,be1,factor_f,term1,term2,term3,term4,term5
      _RL fxkm,bound_a,xlam3lg,dtpt1,summ1,tem_top,xi_top
      _RL fx00,dtdz0,EPST,DIFTM,temp(km),p(km),tempsub(4),psub(4)
      _RL tempout,dy,xfac3
c
      INTEGER mq_adb4,IZHU,M23,M33,MODE,ISTEP,kxx1,iit,iit_out
      INTEGER KINV,jformula,jjregime,jupperbd
      INTEGER jloop,jj1,LOUT,jkmax1,k1
      DATA DAYL/8.64D4/
c
      nlevs = km

      dt=500.0d0*DAYL  !  50.0d0*DAYL    ! integration time interval in seconds
      j_fels=3         !  j_fels=(1,2,3) to call (gamma_fels,gamma_zhu_g,gamma_zhu) for CO cooling
c
      epstcr=1.0d-3*(dt/(50.0d0*DAYL))  ! (or 1.0d-5*...) criterion for converged temperature variation
c
      r_gravity=1.00d0  ! ={0.75,1.30,+0.05} = effect of surface gravity on the planetary escape rates
c
ccccccccccccccccc
c
CAMZ j_planet is now in AMZVARS.h
C      j_planet=1     ! (1,2) = (Pluto, Triton)
      jeuvfuv=2      ! jeuvfuv=(0,1,2*,3)
      f107=130.0d0   ! used when jeuvfuv=0
      f1_mag=0.0d0   ! factor for magnetic heating input
      j1_efuv=1      ! (0,1) = (off,on) of euv-fuv heating: j1_efuv=1 if f1_mag=0; j1_efuv=0 if f1_mag>0
c
      TSUR0=surfaceT(ii,j,bi,bj)      !  fixed surface temperature
      P00=Ro_Surf(ii,j,bi,bj)+etaN(ii,j,bi,bj)          !  fixed surface pressure in Pascal (1 Pa = 10 micro bar) 
      rch4x0=gamma0      ! 22/3670=5.99455d-3   ! CH4 volume mixing ratio near surface      ************
      RCO=comr       ! 1.0d-3     ! CO volume mixing ratio

      dtdz0=-0.1d-3     ! fixed surface temperature gradient in K m^-1, near moisture adabatic for convective adjustment
c  HCN observed by ALMA:
      RHCN1=4.0d-5        ! 4.0d-5     ! HCN volume mixing ratio (2017_0109)
      RHCN2=3.5d-8        ! 3.5d-8     ! HCN volume mixing ratio (2017_0109)
      RHCN3=7.0d-11       ! 7.0d-11    ! HCN volume mixing ratio (2017_0109)
c
      r_c2h2_0=3.3d-3     ! 1.0d-4     ! C2H2 base volume mixing ratio
c
      jj_hcn_sat=0          ! (0,1,2) = (prescribed, saturated, composite) HCN vapor pressure
c

      IF (lsIsConstant .EQ. 0) THEN
       lssol=lsStart
      ELSE
       lssol=lsperpetual
      ENDIF

C AMZ  We allow orbital distance to vary and use flux that varies with lat/lon
C     rsunp_fix=32.9d0      ! a fixed Run-planet distance
      rsunup_fix=semimajoraxis*(1. - eccentricity**2.)/
     &         (1. + eccentricity*COS(lssol - lsp))/1.496e11 
C     surface_albedo=0.8d0  ! surface albedo of icy planets
      surface_albedo=albedoany(ii,j,bi,bj)
c
      CALL GET_QO(ii,j,bi,bj,lssol,0.,solarconstant,
     &                  albedoany(ii,j,bi,bj),qo)

C AMZ we only do this calculation if it is day
      IF (qo .GT. 0.) THEN


ccccccccccccccccc
c

       r00=rSphere     ! Pluto or Triton's radius in m without K_zz; (1152,1120)km = (without,with) K_zz
       capm=plutoMass    ! Pluto or Triton's  mass (M): kg
       if(j_planet.eq.1) then  ! Pluto
c       open(42,file='pluto_ztp_all_out.txt',status='replace')
c       open(44,file='pluto_parameters.txt',status='replace')
       endif
c
       if(j_planet.eq.2) then  ! Triton
c       open(42,file='triton_ztp_all_out.txt',status='replace')
c       open(44,file='triton_parameters.txt',status='replace')
        open(12,file='pluto_triton_dat/mag_heating2a.dat',status='old')
        do k=1,km
         read(12,'(i5,f10.3,f9.2,e14.5)') ix1,zx1,tk1,ht_mag_cgs0(k)
        enddo
        close(unit=12)
       endif
c
C      rsunp0=30.0d0     ! nominal/reference distance to calculate solar flux
       rsunp0=semimajoraxis/1.496E11 ! in general, in AU
c
       do k=1,km
        ht_mag(k)=0.0d0
       enddo
c
       jlowbd=1          ! =(1,2) = fixed (T*, dT/dz) at the lower boundary
cc     kmm=km-1  ! kmm = vertical layers with mid-level altitude and pressure
       jcrit=km
cc
ccccc
c

       CALL GET_MU0(ii,j,bi,bj,lssol,0.,mu0)
       do i=1,itime
        rsunp(i)=rsunup_fix
        rnormf(i)=rsunp(i)/rsunp0  ! normalized distance for calculating the time-dependent solar flux
       enddo
c
       dt00=0.0d0   ! initial temperature perturbation
c
c      
       IYELLE=2  ! IYELLE=1: thin limit; IYELLE=2: opaque effect;
       MQCO=1    ! MQCO=0 => QCO=0;  MQCO=1 => CO cooling included
       MQHCN=1   ! MQHCN=0 => QHCN=0;  MQHCN=1 => HCN cooling included
       MQh2o=1   ! MQh2o=0 => Qh2o=0;  MQh2o=1 => h2o cooling included
       MQ_c2h2=1 ! MQ_c2h2=0 => Q_c2h2=0;  MQ_c2h2=1 => c2h2 cooling included
       IZHU=1    ! IZHU=0 => no V-V;  IZHU=1 => with V-V transitions; for {3.3, 2.3, 1.6} mu bands
       M23=1     ! M23=0 => H(2.3)=0;  M23=1 => H(2.3) .ne. 0
       M33=1     ! M33=0 => H(1.7)=0;  M33=1 => H(1.7) .ne. 0
       ZKVT=1.5D-16      ! collision rate in [m**3 s**(-1)]
       XA3=25.2d0        ! Einstein coefficient in 1/s for 3.3 micron band
       XA4=2.12d0        ! Einstein coefficient in 1/s for 7.6 micron band
       XA2=XA3+XA4       ! Einstein coefficient in 1/s for 2.3 micron band  ! testingxz
       XP2=7.0D-3
       XP3=7.0D-3
       XP4=3.0D-6
       IF(IZHU.EQ.0) THEN
        XA3=4.24D0        ! Roger Yelle's value
        XP3=3.0D-6        ! Roger Yelle's value
        XA4=2.56D0        ! Roger Yelle's value
       ENDIF
c
       bk=kb         ! Boltzmann constant in J K^-1
       capg=bigGrav      ! gravitational constant (G): m^3 s^-2 kg^-1
c
       rcpr=1./atm_kappa
       rrcp=1.0d0/RCPR        ! R/c_p
       wtmol=28.0d0*(1.0d0-rch4x0)+16.0d0*rch4x0  ! mean molecular weight
       rair=atm_Rd          ! gas constant
       cpair=atm_cp              ! specific heat at constant pressure
       cvair=cpair-rair             ! specific heat at constant volume
       g00=gravity      ! Surface gravity: m s**(-2)
c
       FXTOP=1.0D-3*1.0D-3*0.0D0      ! energy flux at the top boundary in [J m^-2 s^-1]
       MODE=2       ! MODE=1: start from a constant T; MODE=2 read in from 'pluto_restart.dat' or 'triton_restart.dat'
       ISTEP=ITIME  ! ITIME/3     ! output after ISTEP
       DELZ0=0.30D3                ! delta z of the lowest layer
       FAC1=1.03D0
       ZM(1)=0.0D3
       ZM(2)=ZM(1)+DELZ0
       DO 108 K=3,KM
        ZM(K)=ZM(K-1)+(ZM(K-1)-ZM(K-2))*FAC1
 108   CONTINUE
c
       DO 109 K=1,KM
        ZKM(K)=ZM(K)/1.0D3
        RKM(K)=R00/1.0D3+ZKM(K)
        RM(K)=RKM(K)*1.0D3
c
c
        hc_ch4=67.0d0
        rh_ch4=1570.0d0   ! 1665.0d0
        rch4st=rch4x0*(1.0d0     ! based on Strobel's 2015.12.07 notes
     &   +dexp((rkm(k)-rh_ch4)/(hc_ch4*(rkm(k)/rh_ch4))))**0.42857143  ! CH4 volume mixing ratio
        rho1(K)=rch4st*(16.0D0/wtmol)         ! base mass mixing ratio of CH4
c
        RHO2(K)=RCO*(28.0D0/wtmol)            ! base mass mixing ratio of CO
c
c  HCN observed by ALMA:
        RHCN1=4.0d-5        ! 4.0d-5     ! HCN volume mixing ratio
        RHCN2=3.5d-8        ! 3.5d-8     ! HCN volume mixing ratio
        RHCN3=7.0d-11       ! 7.0d-11    ! HCN volume mixing ratio
        xf1r=dexp(-dmin1(59.0d0,((rkm(k)-2150.0d0)/490.0d0)**6))
        xf2r=dexp(-dmin1(59.0d0,((rkm(k)-1245.0d0)/25.0d0)**2))
        xf3r=dexp( (rkm(k)-1300.0d0)/(40.0d0*(rkm(k)/1190.0d0)**2) )
        RHO3(K)=(RHCN1*xf1r+RHCN2*xf2r+RHCN3*xf3r)*(27.0d0/wtmol)   ! base mass mixing ratio of HCN
c
c
        xfac1r=dexp( - dmin1(59.0d0, ((rkm(k)-2200.0d0)/525.0d0)**4) )
        xfac2r=dexp( - dmin1(59.0d0, ((rkm(k)-2200.0d0)/655.0d0)**4) )
        rho4(k)=(r_c2h2_0*xfac1r+1.0d-4*xfac2r)*(26.0D0/wtmol)     ! base mass mixing ratio of C2H2
c
c
        Rh2o1=1.30d-5        ! 1.20d-5     ! h2o volume mixing ratio 
        Rh2o2=1.20d-6        ! 1.10d-6     ! h2o volume mixing ratio 
        Rh2o3=4.00d-9        ! 4.00d-9     ! h2o volume mixing ratio 
        Rh2o4=1.50d-8        ! 1.50d-8     ! h2o volume mixing ratio 
        Rh2o5=2.30d-8        ! 2.50d-8     ! h2o volume mixing ratio 
        Rh2o6=1.50d-8        ! 1.50d-8     ! h2o volume mixing ratio 
        xf1r=dexp(-dmin1(59.0d0,((zkm(k)-1060.0d0)/400.0d0)**6))
        xf2r=dexp(-dmin1(59.0d0,((zkm(k)-490.0d0)/120.0d0)**2))
        xf3r=dexp(-dmin1(59.0d0,((zkm(k)-280.0d0)/80.0d0)**2))
        xf4r=dexp(-dmin1(59.0d0,((zkm(k)-210.0d0)/40.0d0)**2))
        xf5r=dexp(-dmin1(59.0d0,((zkm(k)-130.0d0)/30.0d0)**2))
        xf6r=dexp(-dmin1(59.0d0,((zkm(k)-70.0d0)/20.0d0)**2))
        rho5(k)=(Rh2o1*xf1r+Rh2o2*xf2r+Rh2o3*xf3r
     &        +Rh2o4*xf4r+Rh2o5*xf5r+Rh2o6*xf6r)*(18.0d0/wtmol)   ! base mass mixing ratio of h2o
c
c
        IF(plutoTrop) THEN
         weigt1(k)=( 1.0d0-dtanh((zkm(k)-hc)/5.0d0) )/2.0d0
c        xfac1=0.5d0+7.0d0*weigt1(k)
         xfac1=1.5d0+40.0d0*weigt1(k)
         xfac2=1.5d0+30.0d0*dexp(-dmin1(zkm(k)/10.0d0,200.0d0))
        ELSE
         weigt1(k)=0.
         xfac1=0
         xfac2=0
        ENDIF
c       xfact=xfac1*0.3d0+xfac2*0.7d0
        xfact=(xfac1*1.0d0+xfac2*1.0d0)*0.0d0
        dzz1(k)=xfact*0.01d0    ! thermal eddy diffusion coefficient in m^2 s^-1
        dzz2(k)=xfact*0.99d0    ! convective eddy diffusion coefficient in m^2 s^-1
         
 109   CONTINUE
c
       if(j_planet.eq.1) then
        IF(MODE.EQ.1) THEN
         open(25,FILE='pluto_restart.dat',status='replace')
        ENDIF
        IF(MODE.EQ.2) THEN 
         open(25,FILE='pluto_restart.dat',status='old')
        ENDIF
       endif
c
       if(j_planet.eq.2) then
        IF(MODE.EQ.1) THEN
         open(25,FILE='triton_restart.dat',status='replace')
        ENDIF
        IF(MODE.EQ.2) THEN
         open(25,FILE='triton_restart.dat',status='old')
        ENDIF
       endif
c
ccc
c      
       cgamma=3.0d0   ! initial setting of enhancement factor GAMMA
C AMZ no mass flux (also needs to be set to zero again below)
C      frho_in0=12.31d0  ! prescribed escape mass flux
       frho_in0=0.
c
       IF(MODE.EQ.1) THEN            ! start from initial temperature
        TINF1=100.0D0
        do K=1,KM
         XFAC1=1.0D0-(tsur_t(iit)/TINF1)
         TEM(K)=TINF1*(1.0D0-XFAC1*DEXP(-ZM(K)/30.0D3))
         w_in(k)=0.0d0
        enddo
       ENDIF
c
       IF(MODE.EQ.2) THEN            ! start from previous output pluto/triton.restart
        read(25,*) cgamma,frho_in0
        do K=1,KM
         READ(25,*) kxx1,XKM,TEM(K),PREK,RHOK,dn(k),w_in(k)
        enddo
        CLOSE(UNIT=25)
       ENDIF
C AMZ no mass flux, and no vertical velocity
       frho_in0=0.
       DO k=1,KM
        w_in(k)=0.
       ENDDO
cc
       do k=1,km
        tem0(k)=tem(k)
       enddo
c
       do k=2,km
        tem(k)=tem(k)+dt00   ! introduce a perturbation to study the damping rate
       enddo
cc
ccc
       do i=1,itime
        do k=1,km
         tem_t(k,i)=tem0(k)
        enddo
       enddo
c
        write(*,*) 'AMZ flag 01a: FXTOP=',FXTOP
        FLUSH(6)
       DO 700 iit=1,ITIME
c
        write(*,*) 'AMZ flag 01b: FXTOP=',FXTOP
        FLUSH(6)

        iit_out=iit

c
        if(jlowbd.eq.1) tsur_t(iit)=TSUR0  ! *dsqrt(rsunp(1)/rsunp(iit))

        p0_t(iit)=P00
c
c
        if(jlowbd.eq.1) TEM(1)=tsur_t(iit)
        do K=1,KM
         TEMM(K)=TEM(K)
         TEMX(K)=TEM(K)
        enddo
c
        AKK0=5.63D-5     ! thermal conductivity coefficient in J m^-1 s^-1 K^-1
        ALPH=1.12D0      ! thermal conductivity coefficient index (alpha or s)
        ALPH1=ALPH+1.0D0
        ALPHI=1.0D0/ALPH1
c
C        if(iit.eq.1 .and. MODE.EQ.1) then   ! get dn(km) 
C         CALL ZTPRN(ZM,TEM,PRE,RHO,GRV,DN,BNS,HSCAL,KM,p0_t(iit),bk)
C        endif
        write(*,*) 'AMZ flag 02: FXTOP=',FXTOP
        FLUSH(6)
c
C AMZ not doing jeans escape.
C AMZ we call the function that assumes w is zero
        write(*,*) 'p0_t(iit)=',p0_t(iit),', frho_in0=',frho_in0
        do k=1,km
         write(*,*) 'ZM=',ZM(K),', TEM=',TEM(K)
        ENDDO
        CALL ZTPRN(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,p0_t(iit),bk)
        frho=frho_in0
        write(*,*) 'p0_t(iit)=',p0_t(iit),', frho_in0=',frho_in0,
     &                 ', frho=',frho
        write(*,*) 'AMZ flag 03: FXTOP=',FXTOP
        FLUSH(6)
        do k=1,km
         write(*,*) 'ZM=',ZM(K),', TEM=',TEM(K),', PRE=',PRE(K),
     &    ', RHO=',RHO(K),', GRV=',GRV(K),'DN=',DN(K),', BN=',
     &    BN(K),', HSCAL=',HSCAL(K)
        ENDDO

        do k=1,km
         w_out(k)=0. !AMZ
         tem_t(k,iit)=tem(k)
         pre_t(k,iit)=pre(k)
         rho_t(k,iit)=rho(k)
         dn_t(k,iit)=dn(k)
        enddo
c
c  prescribe absorbers profiles, modified on 2015.10.14 based on Strobel's notes
c
        do k=1,km
c
c
         hc_ch4=67.0d0
         rh_ch4=1570.0d0   ! 1665.0d0
         rch4st=rch4x0*(1.0d0     ! based on Strobel's 2015.12.07 notes
     &   +dexp((rkm(k)-rh_ch4)/(hc_ch4*(rkm(k)/rh_ch4))))**0.42857143  ! CH4 volume mixing ratio
c
         dn1(k)=dn(k)*rch4st                 ! CH4 number density
c
         rho2(k)=rco*(28.0d0/wtmol)       ! CO mass mixing ratio
c
         if(jj_hcn_sat.eq.0) then   ! prescried HCN
c
c  HCN observed by ALMA:
          RHCN1=4.0d-5        ! 4.0d-5     ! HCN volume mixing ratio
          RHCN2=3.5d-8        ! 3.5d-8     ! HCN volume mixing ratio
          RHCN3=7.0d-11       ! 7.0d-11    ! HCN volume mixing ratio
          xf1r=dexp(-dmin1(59.0d0,((rkm(k)-2150.0d0)/490.0d0)**6))
          xf2r=dexp(-dmin1(59.0d0,((rkm(k)-1245.0d0)/25.0d0)**2))
          xf3r=dexp( (rkm(k)-1300.0d0)/(40.0d0*(rkm(k)/1190.0d0)**2))
          RHO3(K)=(RHCN1*xf1r+RHCN2*xf2r+RHCN3*xf3r)*(27.0d0/wtmol)   ! base mass mixing ratio of HCN
c
c
         endif
c
         if(jj_hcn_sat.eq.1) then   ! HCN determined by equilibrium saturation vapor pressure
          pa_hcn_sat=1.0d5*dexp(13.93-3.624d3/tem(k)-1.325d5/tem(k)**2
     &     +6.314d6/tem(k)**3-1.128d8/tem(k)**4)                        ! saturation pressure in Pa
          f_super_sat=1.0d0+0.0d0*(1.0d0+dtanh((zkm(k)-300.0d0)/40.0d0))  ! factor for super saturation vapor pressure
          rho3_comp=(pa_hcn_sat/pre(k)/1.03704d0)*f_super_sat
          if(iit.eq.1) then
           rho3(k)=rho3_comp       ! saturated HCN mass mixing ratio
          else
           rho3(k)=(9.0d0*rho3(k)+rho3_comp)/10.0d0 ! relax to vapor pressure equilibrium
          endif
         endif
c
         if(jj_hcn_sat.eq.2) then  ! composite of 1 and 2 for the upper and lower atmospheres

c  HCN observed by ALMA:
          RHCN1=4.0d-5        ! 4.0d-5     ! HCN volume mixing ratio (2017_0109)
          RHCN2=3.5d-8        ! 3.5d-8     ! HCN volume mixing ratio (2017_0109)
          RHCN3=7.0d-11       ! 7.0d-11    ! HCN volume mixing ratio (2017_0109)
          xf1r=dexp(-dmin1(59.0d0,((rkm(k)-2150.0d0)/490.0d0)**6))
          xf2r=dexp(-dmin1(59.0d0,((rkm(k)-1245.0d0)/25.0d0)**2))
          xf3r=dexp( (rkm(k)-1300.0d0)/(40.0d0*(rkm(k)/1190.0d0)**2) )
          RHO3(K)=(RHCN1*xf1r+RHCN2*xf2r+RHCN3*xf3r)*(27.0d0/wtmol)   ! base mass mixing ratio of HCN
c
          pa_hcn_sat=1.0d5*dexp(13.93-3.624d3/tem(k)-1.325d5/tem(k)**2
     &     +6.314d6/tem(k)**3-1.128d8/tem(k)**4)                ! saturation pressure in Pa
          rho3_low=pa_hcn_sat/pre(k)/1.03704d0                    ! HCN mass mixing ratio
          ww_high=(1.0d0+dtanh((zkm(k)-50.0d0)/10.0d0))/2.0d0
          ww_low=1.0d0-ww_high
          rho3_comp=ww_low*rho3_low+ww_high*rho3_high
          if(iit.eq.1) then
           rho3(k)=rho3_comp
          else
           rho3(k)=(9.0d0*rho3(k)+rho3_comp)/10.0d0       ! HCN mass mixing ratio
          endif
         endif
c
         xfac1r=dexp( - dmin1(59.0d0, ((rkm(k)-2200.0d0)/525.0d0)**4) )
         xfac2r=dexp( - dmin1(59.0d0, ((rkm(k)-2200.0d0)/655.0d0)**4) )
         rho4(k)=(r_c2h2_0*xfac1r+1.0d-4*xfac2r)*(26.0D0/wtmol)     ! base mass mixing ratio of C2H2
c
c
         Rh2o1=1.30d-5        ! 1.20d-5     ! h2o volume mixing ratio 
         Rh2o2=1.20d-6        ! 1.10d-6     ! h2o volume mixing ratio 
         Rh2o3=4.00d-9        ! 4.00d-9     ! h2o volume mixing ratio 
         Rh2o4=1.50d-8        ! 1.50d-8     ! h2o volume mixing ratio 
         Rh2o5=2.30d-8        ! 2.50d-8     ! h2o volume mixing ratio 
         Rh2o6=1.50d-8        ! 1.50d-8     ! h2o volume mixing ratio 
         xf1r=dexp(-dmin1(59.0d0,((zkm(k)-1060.0d0)/400.0d0)**6))
         xf2r=dexp(-dmin1(59.0d0,((zkm(k)-490.0d0)/120.0d0)**2))
         xf3r=dexp(-dmin1(59.0d0,((zkm(k)-280.0d0)/80.0d0)**2))
         xf4r=dexp(-dmin1(59.0d0,((zkm(k)-210.0d0)/40.0d0)**2))
         xf5r=dexp(-dmin1(59.0d0,((zkm(k)-130.0d0)/30.0d0)**2))
         xf6r=dexp(-dmin1(59.0d0,((zkm(k)-70.0d0)/20.0d0)**2))
         rho5(k)=(Rh2o1*xf1r+Rh2o2*xf2r+Rh2o3*xf3r
     &        +Rh2o4*xf4r+Rh2o5*xf5r+Rh2o6*xf6r)*(18.0d0/wtmol)   ! base mass mixing ratio of h2o
c
        enddo
c
        BNS1(KM)=DN1(KM)*20.0D3
        DO 31 K=1,KMM
         KINV=KM-K
         FAC11=(ZM(KINV+1)-ZM(KINV))*(DN1(KINV+1)+DN1(KINV))/2.0D0
         BNS1(KINV)=BNS1(KINV+1)+FAC11
  31    CONTINUE
        CALL QTCH4Z(TEM,PRE,RHO1,QT1,KM,iit,GAMN4,IYELLE,grv,
     &          )   ! [nu 4] cooling, 7.6 micron band
        write(*,*) 'AMZ flag 04: FXTOP=',FXTOP
        FLUSH(6)
        CALL HTCH4tr1(HT1,DN1,BNS1,TEM,PRE,KM,iit
     &   ,IYELLE,GAMN4,IZHU,rnormf(iit),surface_albedo,mu0,kb)         ! [nu 3] heating, 3.3 micron band
        write(*,*) 'AMZ flag 05: FXTOP=',FXTOP
        FLUSH(6)
        if(m33.eq.1) THEN
         call HTCH4tr6(HT4,DN1,BNS1,TEM,PRE,KM,iit
     &   ,IYELLE,GAMN4,IZHU,rnormf(iit),surface_albedo,mu0,kb)         !  2 * [nu 3] heating, 1.67 micron band
         write(*,*) 'AMZ flag 06: FXTOP=',FXTOP
         FLUSH(6)
        ELSE
         DO K=1,KM
          HT4(K)=0.
         ENDDO
        ENDIF
        CALL HTCH4tr2(HT3,ht_ch4,ht_n2,DN1,BNS1,TEM,PRE,KM,nm3
     &   ,rnormf(iit),jeuvfuv,f107,dn,bns,j1_efuv)             ! FUV-EUV heating
        write(*,*) 'AMZ flag 07: FXTOP=',FXTOP
        FLUSH(6)
        IF(M23.EQ.1) THEN
         CALL HTCH5t(HT2,DN1,BNS1,TEM,PRE,KM,iit,GAMN4
     &   ,rnormf(iit),surface_albedo)          ! [nu 3] + [nu 4] heating, 2.3 micron band
         write(*,*) 'AMZ flag 08: FXTOP=',FXTOP
         FLUSH(6)
        ELSE
         DO K=1,KM
          HT2(K)=0.0D0         ! setting CH4 2.3 micron heating equal to 0.
         ENDDO
        ENDIF

        IF(MQCO.EQ.1) THEN
         CALL QTCOZ(TEM,PRE,RHO2,QT2,KM,iit,grv,j_fels)    ! CO cooling 
         write(*,*) 'AMZ flag 09: FXTOP=',FXTOP
         FLUSH(6)
        ELSE
         DO K=1,KM
          QT2(K)=0.0D0
         ENDDO
        ENDIF
        IF(MQHCN.EQ.1) THEN
         CALL QTHCNZ(TEM,PRE,RHO3,QT3,KM,iit,grv,j_fels)  ! HCN cooling
         write(*,*) 'AMZ flag 10: FXTOP=',FXTOP
         FLUSH(6)
        ELSE
         DO K=1,KM
          QT3(K)=0.0D0        ! setting the HCN cooling equal to 0.
         ENDDO
        ENDIF
        IF(MQH2O.EQ.1) THEN
         CALL QTH2OZ(TEM,PRE,RHO5,QT5,KM,iit,grv,j_fels)  ! H2O cooling
         write(*,*) 'AMZ flag 11: FXTOP=',FXTOP
         FLUSH(6)
        ELSE
         DO K=1,KM
          QT5(K)=0.0D0        ! setting the H2O cooling equal to 0.
         ENDDO
        ENDIF
        if(MQ_c2h2.eq.1) THEN
         CALL QTC2H2Z(TEM,PRE,RHO4,QT4,KM,iit,IYELLE,grv)   ! C2H2 [nu 5] cooling, 13.7 micron
         write(*,*) 'AMZ flag 12: FXTOP=',FXTOP
         FLUSH(6)
        ELSE
         DO K=1,KM
          QT4(K)=0.0D0        ! setting the C2H2 cooling equal to 0.
         ENDDO
        ENDIF
        DO K=1,KM
         QT7(K)=0.0D0        ! setting CH4 escape induced adiabatic cooling equal to 0.
        ENDDO

c
c
        DO 50 K=1,KM
         HALL(K)=HT1(K)+HT2(K)+ht3(k)+ht4(k)+ht_mag(k)
         QALL(K)=QT1(K)+QT2(K)+QT3(K)+QT4(K)+QT5(K)+qt7(k)
         qnet_r(k)=HALL(K)-QALL(K)          ! Net radiative heating rate in K/sec
         QNETAMZ(K)=qnet_r(k)      ! Net heating rate in K/sec
  50    CONTINUE
  
        DO 120 K=1,KM
         TALP(K)=TEM(K)**ALPH1      ! T^(alpha+1)
         TALPM(K)=TEMM(K)**ALPH1
         TALP1=TALP(K)/TEM(K)       ! T^alpha
         QS(K)=ALPH1*TALP1*QNETAMZ(K)
         DS(K)=AKK0*TALP1/(cpair*RHO(K))    ! K_0*T^alpha/(c_p*rho) ~ K_molecule
         RADS1=R00+ZM(K)
         CS(K)=2.0D0*AKK0*TALP1/(RADS1*cpair*RHO(K))
         bs(k)=0.0d0
 120    CONTINUE
c
c  Add the background thermal eddy diffusion
c
        IF(plutoTrop) THEN
         do k=2,kmm
          xfac1=(zm(k)-zm(k-1))**2
          xfac2=(zm(k+1)-zm(k))**2
          xfac3=(zm(k+1)-zm(k))*(zm(k)-zm(k-1))*(zm(k+1)-zm(k-1))
          wrk1(k)=((dzz1(k+1)-dzz1(k))*xfac1
     &    +(dzz1(k)-dzz1(k-1))*xfac2)/xfac3  ! dK_zz/dz
          wrk2(k)=((tem(k+1)-tem(k))*xfac1+(tem(k)-tem(k-1))*xfac2)
     &               /xfac3  ! dT/dz
         enddo
         wrk1(1)=(dzz1(2)-dzz1(1))/(zm(2)-zm(1))
         wrk1(km)=(dzz1(km)-dzz1(kmm))/(zm(km)-zm(kmm))
         wrk2(1)=(tem(2)-tem(1))/(zm(2)-zm(1))
         wrk2(km)=(tem(km)-tem(kmm))/(zm(km)-zm(kmm))
         do k=1,km
          qs(k)=qs(k)-alph*alph1*dzz1(k)*tem(k)**(alph-1.0d0)*wrk2(k)**2
          ds(k)=ds(k)+dzz1(k)
          cs(k)=cs(k)+wrk1(k)
         enddo
c
c  Add a moisture convective eddy diffusion to smooth the potential temperature
c
c   A=74.17572370,  B=49.67211044  for the factor exp(-A*exp(-B*T)/(p*T))
c
         ae1c=74.17572370d0*1.0d0  ! *0.0d0 => recover/degenerate to the case of dry convetion when A=0
         be1=49.67211044d0
c 
         do k=1,km
          ae1(k)=ae1c*weigt1(k)
         enddo
c 
         do k=1,km
          wrk1(k)=(pre(1)/pre(k))*rrcp  ! PI
          xfac1=tem(k)**alph            ! T^alpha
c
          xfac2=ae1(k)*dexp(-be1/tem(k))/(pre(k)*tem(k))       ! PSI
          wsk1(k)=dexp(-xfac2)                                 ! LAMBDA    
          wsk2(k)=wsk1(k)*(1.0d0+xfac2*(1.0d0-be1/tem(k)))     ! PHI
C
          wrk2(k)=alph1*xfac1/(wrk1(k)*wsk2(k))    ! [(alpha+1)*T^alpha]/(PI*PHI)
          wrk3(k)=dzz2(k)*wsk1(k)/xfac1
          wrk4(k)=dzz2(k)*wrk1(k)*wsk2(k)/(alph1*xfac1)   ! YY = [related to 2nd term before d[zeta]/dz]
         enddo
c 
         do k=2,kmm
          xfac1=(zm(k)-zm(k-1))**2
          xfac2=(zm(k+1)-zm(k))**2
          xfac3=(zm(k+1)-zm(k))*(zm(k)-zm(k-1))*(zm(k+1)-zm(k-1))
          wrk5(k)=((wrk1(k+1)-wrk1(k))*xfac1
     &    +(wrk1(k)-wrk1(k-1))*xfac2)/xfac3  ! d PI/dz
          wrk6(k)=((wrk4(k+1)-wrk4(k))*xfac1
     &    +(wrk4(k)-wrk4(k-1))*xfac2)/xfac3  ! d YY/dz
         enddo
         wrk5(1)=(wrk1(2)-wrk1(1))/(zm(2)-zm(1))
         wrk5(km)=(wrk1(km)-wrk1(kmm))/(zm(km)-zm(kmm))
         wrk6(1)=(wrk4(2)-wrk4(1))/(zm(2)-zm(1))
         wrk6(km)=(wrk4(km)-wrk4(kmm))/(zm(km)-zm(kmm))
         do k=1,km
          wrk7(k)=wrk3(k)*wrk5(k)
          wsk3(k)=wrk2(k)*(wrk7(k)+wrk6(k))  ! coefficient for the second term
         enddo
c
         do k=2,kmm
          xfac1=(zm(k)-zm(k-1))**2
          xfac2=(zm(k+1)-zm(k))**2
          xfac3=(zm(k+1)-zm(k))*(zm(k)-zm(k-1))*(zm(k+1)-zm(k-1))
          wrk8(k)=((wrk7(k+1)-wrk7(k))*xfac1
     &    +(wrk7(k)-wrk7(k-1))*xfac2)/xfac3
         enddo
         wrk8(1)=(wrk7(2)-wrk7(1))/(zm(2)-zm(1))
         wrk8(km)=(wrk7(km)-wrk7(kmm))/(zm(km)-zm(kmm))
c
         do k=1,km
          ds(k)=ds(k)+dzz2(k)
          cs(k)=cs(k)+wsk3(k)
          bs(k)=bs(k)+wrk2(k)*wrk8(k)
         enddo
        ENDIF !AMZ end troposphere section
c

c
C AMZ
        FXKM=FXTOP*ALPH1/AKK0    ! temperature gradient [d(T^(1+alpha))/dz] at top
        write(*,*) 'FXKM=',FXKM
        DO K=1,KM
C         write(*,*) 'QS=',QS(k),', DS=',DS(k),', CS=',
C     &                                   CS(k),', BS=',BS(k)
C         flush(6)
        ENDDO

        if(jlowbd.eq.1) then
         CALL IMDIFF2(TALP,QS,DS,CS,bs,KM,ZM,DT,FXKM,QDIFF)
        endif
c
        if(jlowbd.eq.2) then
         fx00=dtdz0*tem(1)**ALPH*ALPH1
         CALL IMDjFF2(TALP,QS,DS,CS,bs,KM,ZM,DT,fx00,FXKM,QDIFF)
        endif

        DO 180 K=1,KM
         TEM(K)=TALP(K)**ALPHI
         TEMM(K)=TALPM(K)**ALPHI
 180    CONTINUE
c
        DO 200 K=1,KM
         TALP(K)=TEM(K)**ALPH1
         TALP1=TALP(K)/TEM(K)
         QDIFF(K)=QDIFF(K)/(ALPH1*TALP1)
 200    CONTINUE      

        LOUT=0
        EPST=0.0D0
        jkmax1=0
        DO 210 K=1,KM
         DIFTM=(TEM(K)-TEMX(K))
         IF(ABS(DIFTM).GT.ABS(EPST)) then
          EPST=DIFTM
          jkmax1=k
         endif
         TEMX(K)=TEM(K)
 210    CONTINUE
        IF(ABS(EPST).LE.EPSTCR) LOUT=999
        IF(iit.LE.4) LOUT=0

        IF(LOUT.EQ.999) GO TO 701
c
 700   CONTINUE
c
 701   CONTINUE
c
C Interpolate to GCM rC grid

       DO k=1,km
        temp(k)=tem(k)
        p(k)=pre(k)
C        write(*,*) 'k=',k,', temp=',temp(k),', p=',p(k)
       ENDDO


       DO k=1,Nr

        CALL LOCATE(p,nlevs,rC(k),pindex)

        IF(pindex .EQ. 1) THEN
         pinterpindexstart=1
         skipinterp=0
        ELSE IF (pindex .EQ. nlevs-1) THEN
         pinterpindexstart=nlevs-3
         skipinterp=0
        ELSE IF (pindex .GE. nlevs) THEN
C         pinterpindexstart=nlevs-3
C         write(*,*) 'Locate: p out of bounds high'
          plutoT(k)=temp(nlevs)
          skipinterp=1
        ELSE IF (pindex .EQ. 0) THEN
         pinterpindexstart=1
C         write(*,*) 'Locate: p out of bounds low'
         skipinterp=0
        ELSE
         pinterpindexstart=pindex-1
        ENDIF

        IF(skipinterp .EQ. 0) THEN
         DO k1=1,4
          tempsub(k1)=temp(pinterpindexstart+k1-1)
          psub(k1)=p(pinterpindexstart+k1-1)
         ENDDO
         CALL POLINT(psub,tempsub,4,rC(k),tempout,dy)
C         write(*,*) 'tempout=',tempout
         plutoT(k)=tempout
        ENDIF
       ENDDO
C goes with qo .eq. 0. loop
      ELSE
       DO k=1,Nr
        plutoT(k)=1. _d -8
       ENDDO
      ENDIF

C      write(*,*) i,',',j,',',bi,',',bj
C      DO k=1,Nr
C       write(*,*) 'k=',k,'plutoT=',plutoT(k)
C      ENDDO
C      flush(6)

      RETURN
      END
