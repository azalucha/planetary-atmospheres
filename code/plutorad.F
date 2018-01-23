#include "CPP_OPTIONS.h"

      SUBROUTINE PlUTORAD(bi,bj,kLev,myTime,gt_arr)
C     *==========================================================*
C     | S/R PLUTORAD
C     | o Pluto radiation scheme with heating rate supplied directly.
C     |   Non-LTE CH4 heating at 3.3 um and cooling at 7.6 um.
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

      INTEGER kLev, bi, bj,k,i,j
      _RL myTime

      _RL ko,alpha,wl,wn,bs,a10,mun2,much4,mi,molwgt,p10
      _RL z10,g,wll,wnl,nul,a10l,c10,epsln,qo,eps
      _RL epslnl,epsl,piF,na,n,rho,q,l,temp(Nr)
      _RL termP(Nr),termPf(Nr+1),dtdp(Nr),conduction
      _RL dtdp2(Nr),dTdt,iceThickness,lssol
      _RL psurf,piF1,cp,tfrostsamz
      _RL qs,p10l,c10l,a10s,wls,bss,p10s,c10s,epslns,epss
      _RL piF0s,piF1s,wns,piFs,bigG1,littleG1,dtedz(Nr),dtedz2(Nr)
      _RL weigt1(Nr),kt(Nr),kc(Nr),z(Nr),dz(Nr),dtdz2(Nr)
      _RL dtedT(Nr),thetae(Nr),xfac1,xfac2,bigG2,div,dtdz(Nr)
      _RL scalefac,Tv
      _RL gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      DO j=1,sNy
       DO i=1,sNx
C  thermodynamics loop counts down from Nr.  Do calculations
C  for entire column and surface at first instance
        IF(kLev .EQ. Nr) THEN
         CALL DO_SURFACE_CALC(i,j,bi,bj,myTime)
C  Call on first iteration then every radiation time step
         IF(initrad(i,j,bi,bj) .EQ. 1. .OR.
     &            MOD(myTime,radTimeStep).EQ.0.) THEN
            initrad(i,j,bi,bj)=0

C diffusion coefficient (J/m/s/K^(alpha+1))
          ko=5.63 _d -5
C exponent on temperature in diffusion term
          alpha=1.12
C wavelength for heating (microns)
          wl=3.3
C convert to wavenumber for heating (m^-1)
          wn=1.e6/wl
C band strength (m)
          bs=30.*(wn/wl)/nlosch
C          bs=34.4*(wn/wl)/nlosch
C Einstein coefficient at 3.3 microns (s^-1)
          a10=4.24
C          a10=25.2
C molecular weight of N2
          mun2=28.
C molecular weight of CH4
          much4=16.0426
C Molecular weight
          mi=gamma0*much4+(1.-gamma0)*mun2
C         molwgt=(mun2*much4)/(mun2+much4)
C         mi=molwgt
C probability that a V-T transition occurs (3.3	microns)
          p10=1.d-6
C          p10=7.d-3
C probability that a V-T transition occurs (7.6 microns)
          p10l=p10
C          p10l=3.d-6
C collision rate (m^3/s)
          z10=1.5 _d -16
C ratio of statistical weights
          g=3.
C big G1 in CC equation (Pa K)
          bigG1=35.0
C big G2 in CC equation (K)
          bigG2=39.2
C wavelength for cooling (microns)
          wll=7.8
C convert to wavenumber (m^-1)
          wnl=1. _d 6/wll
C convert to frequency (1/s)
          nul=clight*wnl
C Einstein coefficient at 7.6 microns (s^-1)
          a10l=2.56
C          a10l=2.12
          c10=p10*z10*0.1/kb/100.
          epsln=c10/a10
C epsilon at 3.3 micron
          eps=epsln/(1.+epsln)
          c10l=p10l*z10*0.1/kb/100.
          epslnl=c10l/a10l
C epsilon at 7.6 micron
          epsl=epslnl/(1.+epslnl)
c 2.3 micron stuff
C         a10s=25.8
C         wls=2.3
C         wns=1.d6/wls
C         bss=2.89*(wns/wls)/nlosch
C         p10s=p10
C         p10s=7.d-3
C         c10s=p10s*z10*0.1/kb/100.
C         epslns=c10s/a10s
C         epss=epslns/(1.+epslns)
C solar forcing at 3.3 um
          IF (lsIsConstant .EQ. 0) THEN
           CALL GET_LS(myTime,lssol)
          ELSE
           lssol = lsperpetual
          ENDIF

          CALL GET_QO(i,j,bi,bj,lssol,myTime,piF0,0.,piF1)
          piF=piF1*(wl/wn)
C solar	forcing	at 2.3 um
C         piF0s=21.3643 _d -3
C         CALL GET_QO(i,j,bi,bj,lssol,0.,piF0s,0.,piF1s)
C         piFs=piF1s*(wls/wns)



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
     &           /Ro_surf(i,j,bi,bj)
            termPf(k)=rF(k)*(Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj))
     &           /Ro_surf(i,j,bi,bj)
           ENDDO
           termPf(Nr+1)=rF(Nr+1)*(Ro_surf(i,j,bi,bj)+
     &       etaN(i,j,bi,bj))/Ro_surf(i,j,bi,bj)
           psurf=Ro_surf(i,j,bi,bj)+etaN(i,j,bi,bj)
          ENDIF
          DO k=1,Nr
           temp(k)=theta(i,j,k,bi,bj)*
     &                   (termP(k)/atm_po)**atm_kappa
          ENDDO
          dtdp(1)=(temp(2)-surfaceT(i,j,bi,bj)
     &                 )/(termP(2)-psurf)
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
           na=termP(k)/kb/temp(k)
           n=gamma0*na
           rho=mi*amu*na
C 3.3 micron heating (J/m^3/s)
           q=piF*eps*n*bs
C 2.3 micron heating (J/m^3/s)
C           qs=piFs*epss*n*bss
C cooling rate
	   l=g*epsl*n*hplanck*nul*a10l*
     &            DEXP(-hplanck*nul/kb/temp(k))
           qr(i,j,k,bi,bj)=(q-l+conduction)/atm_cp/rho
           cond1(i,j,k,bi,bj)=conduction/atm_cp/rho
     &             *(rC(k)/atm_po)**(-atm_kappa)
           q1(i,j,k,bi,bj)=q/atm_cp/rho
     & 	       	   *(rC(k)/atm_po)**(-atm_kappa)
           l1(i,j,k,bi,bj)=l/atm_cp/rho
     & 	       	   *(rC(k)/atm_po)**(-atm_kappa)
          ENDDO
          IF(plutoTrop) THEN
           DO k=1,Nr
            dz(k)=-(termPf(k+1)-termPf(k))*atm_Rd*temp(k)/
     &                  gravity/termPf(k)
           ENDDO
           z(1)=0.
           DO k=2,Nr
            z(k)=z(k-1)+dz(k)
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
         ENDIF
        ENDIF
        dTdt=qr(i,j,kLev,bi,bj)*(rC(klev)/atm_po)**(-atm_kappa)
        dTdtsave(i,j,kLEv,bi,bj)=dTdt
        gT_arr(i,j)=gT_arr(i,j)+dTdt*maskC(i,j,kLev,bi,bj)
       ENDDO
      ENDDO


      RETURN
      END

