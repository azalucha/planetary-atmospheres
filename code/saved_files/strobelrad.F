#include "CPP_OPTIONS.h"

      SUBROUTINE STROBELRAD(bi,bj,kLev,myTime,gT_arr)
C     *==========================================================*
C     | S/R STROBELRAD
C     | o Heating rates on Pluto based on Strobel et al. 1996
C     |   (with troposphere specification of Zalucha et al. 2011)
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

      INTEGER  bi,bj,k,i,j,nlevs,pindex,lm,itime,km
      INTEGER pinterpindexstart,skipinterp,kLev
      _RL plutoT(Nr),PlutoSurfT,myTime


      PARAMETER(KM=Nr,LM=5,ITIME=2)        ! KM= Nr 
      COMMON /PARAIR/ R00,WTMOL
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      _RL ZM(KM),GRV(KM),TEM(KM),TEMM(KM),TEMX(KM)
     &    ,PRE(KM),RHO(KM),DN(KM),BNS(KM),HSCAL(KM)                  ! Basic state
      _RL RHO1(KM),RHO2(KM),DN1(KM),BNS1(KM)
      _RL RM(KM),ZKM(KM),RKM(KM),PREBAR(KM)
      _RL HT1(KM),HALL(KM),QT1(KM),QT2(KM),QALL(KM),QNETAMZ(KM)
      _RL HT1CGS(KM),QT1CGS(KM),QT2CGS(KM),QNETCGS(KM)
      _RL HT2(KM),HT2CGS(KM)
      _RL QDIFF(KM),QDFCGS(KM)
      _RL QS(KM),DS(KM),CS(KM),bs(km),TALP(KM),TALPM(KM),GAMN4(KM)
      _RL tobs(km,3),tsig(km,3),t1err(itime),t2err(itime)  ! measured temperature and standard deviation
      _RL dzz1(km),dzz2(km)  ! background and convective eddy diffusion coefficient
     & ,wrk1(km),wrk2(km),wrk3(km),wrk4(km),wrk5(km),wrk6(km)
     & ,wrk7(km),wrk8(km),wrk9(km),wsk1(km),wsk2(km),wsk3(km)
      _RL ae1(km),weigt1(km)
      _RL tsecd(itime),tyear(itime),rsunp(itime),rnormf(itime)
      _RL tsur_t(itime),p0_t(itime),tem_t(km,itime)
     & ,pre_t(km,itime),rho_t(km,itime),dn_t(km,itime)
c
      _RL tem0(km),dtem0(km,itime),tau0(km)
c
      _RL p(KM),psub(4),tempsub(4),dy,tempout,lssol,ts
      _RL tfrostsamz,dt,epstcr,tsur0,p00,rch4x0,ch4ne,rco
      _RL dtdz0,zkvt,xa3,xa4,xa2,xp2,xp3,xp4,bk,r00
      _RL rcpr,rrcp,wtmol,rair,cp,fxtop,delz0,fac1,xfac1,xfac2
      _RL xfact,tinf,xkm,prek,rhok,bnn,akk0,alph,alph1,alphi
      _RL hfac,rch4st,factx1,talp1,rads1,xfac3,ae1c,be1,fxkm,fx00
      _RL epst,dift0,diftm,dt00,epst2,fac11,tinf1,temp(Nr)
      _RL gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      INTEGER kmm,jlowbd,jopen1,jobs01,iyelle,mqco,izhu,m23,dayl
      INTEGER mode,iit_out,kinv,lout,iit,istep,k1,kbg2,ked2,kdif2

      _RL termP(Nr),termPf(Nr+1),rdiff,conduction,dtdp(Nr),dtdp2(Nr)
      _RL dz(Nr),z(Nr),dtedz(Nr),dtdz(Nr),dtedz2(Nr),dtdz2(Nr)
      _RL div,dTdt,scalefac,tfrostamz,iceThickness,cpg,rhog
      _RL kt(Nr),kc(Nr),bigG1,bigG2,alpha,ko,psurf,thetae(Nr)
      _RL dtedt(Nr)

      DATA DAYL/8.64D4/

C diffusion coefficient (J/m/s/K^(alpha+1))
      ko=5.63 _d -5
C exponent on temperature in diffusion term
      alpha=1.12
C big G1 in CC equation (Pa K)
      bigG1=35.0
C big G2 in CC equation (K)
      bigG2=39.2


      DO j=1,sNy
       DO i=1,sNx
C  thermodynamics loop counts down from Nr.  Do calculations
C  for entire column and surface at first instance
        IF(kLev .EQ. Nr) THEN
         CALL DO_SURFACE_CALC(i,j,bi,bj,myTime)
C Do radiative transfer calculation every radiation time step
         IF(initrad(i,j,bi,bj) .EQ. 1. .OR.
     &            MOD(myTime,radTimeStep).EQ.0.) THEN

          kmm=km-1

          jlowbd=1   ! =(1,2) = fixed (T, dT/dz) at the lower boundary
c
c         dt=1.0d6     ! time step in seconds determined consistently in tvsf2.par
          DT=500.0d0*DAYL  ! 500.0D0*DAYL           ! time interval in seconds
          EPSTCR=1.0D-5*(dt/(500.0d0*DAYL))       ! criterion


          IF (lsIsConstant .EQ. 0) THEN
           lssol=lsStart
          ELSE
           lssol=lsperpetual
          ENDIF


          IF(selectAddFluid.EQ.0) THEN
           DO k=1,Nr
            termP(k)=rC(k)
            termPf(k)=rF(k)
           ENDDO
           termPf(Nr+1)=rF(Nr+1)
           psurf=Ro_surf(i,j,bi,bj)
          ELSE
           DO k=1,Nr
            termP(k)=rC(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &      /Ro_surf(i,j,bi,bj)
            termPf(k)=rF(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &      /Ro_surf(i,j,bi,bj)
           ENDDO
           termPf(Nr+1)=rF(Nr+1)*
     &     (Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &     /Ro_surf(i,j,bi,bj)
           psurf=Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj)
          ENDIF
          DO k=1,Nr
           temp(k)=theta(i,j,k,bi,bj)*
     &                   (termP(k)/atm_po)**atm_kappa
          ENDDO
C Conduction terms
          dtdp(1)=(temp(2)-surfaceT(i,j,bi,bj))
     &               /(termP(2)-psurf)
          DO k=2,Nr-1
           dtdp(k)=(temp(k+1)-temp(k-1))
     &                /(termP(k+1)-termP(k-1))
          ENDDO
          dtdp(Nr)=(temp(Nr)-temp(Nr-1))
     &                /(termP(Nr)-termP(Nr-1))
          DO k=2,Nr-1
           dtdp2(k)=(dtdp(k+1)-dtdp(k-1))/(termP(k+1)-termP(k-1))
          ENDDO
          dtdp2(Nr)=(dtdp(Nr)-dtdp(Nr-1))/(termP(Nr)-termP(Nr-1))
          dtdp2(1)=(dtdp(2)-dtdp(1))/(termP(2)-termP(1))
          DO k=1,Nr
           conduction=ko*termP(k)*gravity**2.
     &                      /atm_Rd**2./temp(k)*
     &            (temp(k)**(alpha-1.)*termP(k)*dtdp2(k)+
     &             temp(k)**(alpha-1.)*dtdp(k)+
     &     termP(k)*dtdp(k)**2.*(alpha-1.)*temp(k)**(alpha-2.))
           cond(i,j,k,bi,bj)=conduction
          ENDDO
C end conduction terms
          TSUR0=surfaceT(i,j,bi,bj)
          P00=psurf
          RCH4X0=gamma0
          CH4NE=0.
          RCO=comr

          dtdz0=-0.1d-3     ! fixed surface temperature gradient in K m^-1, near moisture adabatic for convective adjustment


          jopen1=0

          dt00=0.0d0   ! initial temperature perturbation
c
c
c
ccccccccccccccccccccccccccccccccccc
c
          jobs01=0   ! = (0,1) = (off,on) for using the measured temperature
c
          if(jobs01.eq.1) then
           open(27,FILE='test_para.out',status='new')
          endif

          IYELLE=2  ! IYELLE=1: thin limit; IYELLE=2: opaque effect;
          MQCO=1    ! MQCO=0 => QCO=0;  MQCO=1 => CO cooling included
          IZHU=1    ! IZHU=0 => no V-V;  IZHU=1 => with V-V transitions; for 3.3 mu band
          M23=1     ! M23=0 => H(2.3)=0;  M23=1 => H(2.3)�0
CCCCCCCCCCCCCCAMZ THESE LINES ARE NOT NECESSARY HERE BECAUSE THEY ARE INPUT FROM COMMON BLOCK
cc        RCH4X0=0.004D0             ! CH4 volume mixing ratio at surface      ************
cc        CH4NE=0.00d0  ! 0.0D0  ! 1.5d0                ! mass mixing ratio of CH4 falls exponentially
cc        RCO=5.0d-4  ! 4.6d-4  ! 3.0D-4                 ! CO volume mixing ratio
c         P00=0.30D0                 !  fixed surface pressure in Pascal (1 Pa = 10 micro bar)
c         P00=1.8D0                 !  fixed surface pressure in Pascal (1 Pa = 10 micro bar)
c         TSUR0=35.0D0               ! fixed surface temperature
CCCCCCCCCCCCCAMZ END
CCCCCCCCCCCCCAMZ NOT REALLY SURE WHY XUN HAS THESE HERE, THEY'RE SET ABOVE
c         DT=500.0D0*DAYL           ! time interval in seconds
c         ITIME=60                   ! total steps in time integration = 60, 180, 360, 450
c         EPSTCR=1.0D-4              ! criterion for converged temperature variation
c         EPSTCR=1.0D-9	  ! criterion to be set at a much smaller value to allow a full integration
CCCCCCCCCCCCCCCAMZ END
          ZKVT=1.5D-16	! collision rate in [m**3 s**(-1)]
          XA3=25.2D0        ! Einstein coefficient in 1/s for 3.3 micron band
          XA4=2.12D0        ! Einstein coefficient in 1/s for 7.6 micron band
          XA2=XA3+XA4    ! Einstein coefficient in 1/s for 2.3 micron band  ! testingxz
          XP2=7.0D-3
          XP3=7.0D-3
          XP4=3.0D-6
          IF(IZHU.EQ.0) THEN
           XA3=4.24D0        ! Roger Yelle's value
           XP3=3.0D-6        ! Roger Yelle's value
           XA4=2.56D0        ! Roger Yelle's value
          ENDIF

          BK=kb            ! Boltzmann constant in J K^-1
          R00=rSphere           ! Pluto's radius in m
          RCPR=1./atm_kappa
          rrcp=1.0d0/RCPR        ! R/c_p
          WTMOL=28.0D0*(1.0D0-RCH4X0)+16.0D0*RCH4X0  ! mean molecular weight
          RAIR=atm_Rd	! gas constant
          CP=atm_cp     ! specific heat at the constant pressure
          FXTOP=1.0D-3*1.0D-3*0.0D0        ! energy flux at the top boundary in [J m^-2 s^-1]
          ISTEP=ITIME  ! ITIME/3     ! output after ISTEP
          rdiff=termPf(2)-termPf(1)
          ZM(1)=-RAIR*temp(1)/termP(1)/gravity*rdiff
          DO K=2,KM
           rdiff=termPf(k+1)-termPf(k)
           dz(k)=-RAIR*temp(K)/termP(k)/gravity*rdiff
           ZM(k)=ZM(k-1)+dz(k)
          ENDDO
          dz(1)=ZM(2)-ZM(1)
          DO 109 K=1,KM
           RHO1(K)=RCH4X0*(16.0D0/28.0D0)   ! mass mixing ratio of CH4
           RHO2(K)=RCO*(28.0D0/28.0D0)   ! mass mixing ratio of CO
           ZKM(K)=ZM(K)/1.0D3
           RKM(K)=R00/1.0D3+ZKM(K)
           RM(K)=RKM(K)*1.0D3
C Troposphere not functional at this time!
           IF(plutoTrop) THEN
            weigt1(k)=( 1.0d0-dtanh((zkm(k)-hc)/5.0d0) )/2.0d0
            xfac1=1.5d0+40.0d0*weigt1(k)
            xfac2=1.5d0+30.0d0*dexp(-zkm(k)/10.0d0)
            xfact=(xfac1*1.0d0+xfac2*1.0d0)*0.0d0
            dzz1(k)=xfact*0.01d0    ! thermal eddy diffusion coefficient in m^2 s^-1
            dzz2(k)=xfact*0.99d0    ! convective eddy diffusion coefficient in m^2 s^-1
           ENDIF
 109      CONTINUE
c

          IF(initrad(i,j,bi,bj) .EQ. 1) THEN
           iit=1
           initrad(i,j,bi,bj)=0
          ELSE
           iit=2
          ENDIF
c
          iit_out=iit

c
          if(jlowbd.eq.1) tsur_t(iit)=TSUR0  ! *dsqrt(rsunp(1)/rsunp(iit))

          p0_t(iit)=P00

          DO k=1,KM
           TEM(K)=temp(K) 

          ENDDO
ccc
          do K=1,KM
           TEMM(K)=TEM(K)
           TEMX(K)=TEM(K)
          enddo
c

          AKK0=5.63D-5     ! thermal conductivity coefficient in J m^-1 s^-1 K^-1
          ALPH=1.12D0      ! thermal conductivity coefficient index (alpha or s)
          ALPH1=ALPH+1.0D0
          ALPHI=1.0D0/ALPH1

          CALL ZTPRN(ZM,TEM,PRE,RHO,GRV,DN,BNS,HSCAL,KM,p0_t(iit),
     &     kb,gravity,atm_Rd)
c


          do k=1,km
           tem_t(k,iit)=tem(k)
           pre_t(k,iit)=pre(k)
           rho_t(k,iit)=rho(k)
           dn_t(k,iit)=dn(k)
          enddo
c
          DO 30 K=1,KM
           HFAC=(DN(K)/DN(1))**CH4NE   ! mass mixing ratio of CH4 falls exponentially
           RCH4ST=RCH4X0*HFAC               ! CH4 mixing ratio given by D. F. Strobel
           DN1(K)=DN(K)*RCH4ST              ! number density of CH4
           RHO1(K)=RCH4ST*(16.0D0/28.0D0)   ! mass mixing ratio of CH4
  30      CONTINUE
          BNS1(KM)=DN1(KM)*20.0D3
          DO 31 K=1,KMM
           KINV=KM-K
           FAC11=(ZM(KINV+1)-ZM(KINV))*(DN1(KINV+1)+DN1(KINV))/2.0D0
           BNS1(KINV)=BNS1(KINV+1)+FAC11
  31      CONTINUE
          CALL QTCH4Z(TEM,PRE,RHO1,QT1,KM,iit,GAMN4,IYELLE,
     &    gravity,atm_Rd,atm_cp)
          CALL HTCH4tr(HT1,DN1,BNS1,TEM,PRE,KM,iit,IYELLE,GAMN4,IZHU
     &    ,i,j,bi,bj,lssol,albedoany(i,j,bi,bj),2,kb,atm_kappa,atm_Rd)
          IF(M23.EQ.1) CALL HTCH5t(HT2,DN1,BNS1,TEM,PRE,KM,iit,GAMN4
     &    ,i,j,bi,bj,lssol,albedoany(i,j,bi,bj),kb,atm_kappa,atm_Rd)
          IF(MQCO.EQ.1) CALL QTCOZ(TEM,PRE,RHO2,QT2,KM,iit,
     &      gravity,atm_Rd,cp)
          IF(MQCO.EQ.0) THEN
           DO 42 K=1,KM
  42        QT2(K)=0.0D0        ! setting the CO cooling equal to 0.
          ENDIF
          IF(M23.EQ.0) THEN
           DO 44 K=1,KM
  44        HT2(K)=0.0D0         ! setting CH4 2.3 micron heating equal to 0.
          ENDIF
          DO 50 K=1,KM
           HALL(K)=HT1(K)+HT2(K)
           QALL(K)=QT1(K)+QT2(K)
           QNETAMZ(K)=HALL(K)-QALL(K)     ! Net heating rate in K/sec
           qr(i,j,k,bi,bj)=QNETAMZ(K)
     &          +cond(i,j,k,bi,bj)/atm_cp/rho(k)
           HT33SAVE(i,j,k,bi,bj)=HT1(K)
           HT23SAVE(i,j,k,bi,bj)=HT2(K)
           QT76SAVE(i,j,k,bi,bj)=QT1(K)
           QTCOSAVE(i,j,k,bi,bj)=QT2(K)
           CONDSAVE(i,j,k,bi,bj)=cond(i,j,k,bi,bj)/atm_cp/rho(k)
  50      CONTINUE
C Add in the troposphere.  Not functional at this time!
          IF(plutoTrop) THEN
           DO k=1,Nr
            z(k)=zm(k)
           ENDDO
           DO k=1,Nr
            weigt1(k)=( 1.0-tanh((z(k)/1000.-hc)/5.0) )/2.0
            xfac1=40.0*weigt1(k)
            xfac2=1.5+30.0*exp(-z(k)/1000./10.0)
            kt(k)=1.0 _d -9
            kc(k)=xfac1
            thetae(k)=theta(i,j,k,bi,bj)*EXP(
     &      -bigG1*EXP(-bigG2/temp(k))/termP(k)/temp(k))
           ENDDO
           DO k=1,Nr-1
            dtedz(k)=(thetae(k+1)-thetae(k))/dz(k)
            dtdz(k)=(temp(k+1)-temp(k))/dz(k)
           ENDDO
           dtedz(Nr)=dtedz(Nr-1)
           dtdz(Nr)=dtdz(Nr-1)
           DO k=1,Nr-1
            dtedz2(k)=(kc(k+1)*dtedz(k+1)-kc(k)*dtedz(k))/dz(k)
            dtdz2(k)=(kt(k+1)*dtdz(k+1)-kt(k)*dtdz(k))/dz(k)
           ENDDO
           dtedz2(Nr)=dtedz2(Nr-1)
           dtdz2(Nr)=dtdz2(Nr-1)
           DO k=1,Nr
            dtedT(k)=thetae(k)*(1./temp(k)+
     &         (-termP(k)*temp(k)*bigG1*EXP(-bigG2/temp(k))
     &         *(bigG2/temp(k)**2.)
     &         +bigG1*EXP(-bigG2/temp(k))*termP(k))
     &         /(termP(k)*temp(k))**2.)
            IF(k .EQ. Nr-1 .OR. k .EQ. Nr) THEN
             div=0.
            ELSE
             div=dtedz2(k)/dtedT(k)
            ENDIF
            qr(i,j,k,bi,bj)=qr(i,j,k,bi,bj)+dtdz2(k)+div
           ENDDO
          ENDIF
C ----------end troposphere specification

         ENDIF
        ENDIF
        dTdt=qr(i,j,kLev,bi,bj)*(rC(klev)/atm_po)**(-atm_kappa)
        gT_arr(i,j)=gT_arr(i,j)+dTdt
     &         *maskC(i,j,kLev,bi,bj)
       ENDDO
      ENDDO



      RETURN
      END
