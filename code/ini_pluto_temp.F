#include "CPP_OPTIONS.h"

      SUBROUTINE INI_PLUTO_TEMP(i,j,bi,bj,plutoT,ts)
C     *==========================================================*
C     | S/R INI_PLUTO_TEMP
C     | o Calculate radiative equilibrium temperature to initialize
C     |   Pluto temperatures
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "AMZVARS.h"

      INTEGER  bi, bj,k,i,j,k1
      _RL plutoT(Nr),ts

      INTEGER nlevs,niters,nmax,flag,pindex
      INTEGER pinterpindexstart,skipinterp

C For Pluto
      PARAMETER(nlevs=140)
C For Triton
C     PARAMETER(nlevs=120)
      _RL ko,alpha,wl,wn,bs,a10,mun2,much4,mi,molwgt,p10
      _RL z10,g,wll,wnl,nul,a10l,c10,epsln,eps,epslnl,epsl
      _RL lssol,piF1,piF,rs,deltaTp,error,tempout,dy
      _RL alpham(nlevs),lambdam(nlevs),omegam(nlevs),qm(nlevs)
      _RL um(nlevs),bigG(nlevs),bigC(nlevs),bigS(nlevs)
      _RL rm(nlevs),drmm(nlevs),drmp(nlevs),temp(nlevs)
      _RL xi(nlevs),tempPrev(nlevs),xiPrev(nlevs),pint(nlevs)
      _RL p(nlevs),n(nlevs),na(nlevs),rho(nlevs),q(nlevs)
      _RL l(nlevs),rnet(nlevs),drmb(nlevs)
      _RL tempsub(4),psub(4),p10l,c10l
      _RL weigt1(nlevs),kt(nlevs),kc(nlevs)
      _RL bigPi(nlevs),bigPsi(nlevs),bigTheta(nlevs)
      _RL dktdz(nlevs),dbigPidz(nlevs),dtempdz(nlevs)
      _RL dinphi2(nlevs),dinphi3(nlevs),bigB(nlevs)
      _RL bigG1,bigG2,littleG,zm(nlevs)
      _RL xfac1,xfac2,dkp1,dkp

C diffusion coefficient (J/m/s/K^(alpha+1))
      ko=5.63 _d -5
C exponent on temperature in diffusion term
      alpha=1.12
C wavelength for heating (microns)
      wl=3.3
C convert to wavenumber for heating (m^-1)
      wn=1000000./wl
C band strength (m)
      bs=30.*(wn/wl)/nlosch
C      bs=34.4*(wn/wl)/nlosch
C Einstein coefficient at 3.3 microns (s^-1)
      a10=4.24
C      a10=25.2
C molecular weight of N2
      mun2=28.
C molecular weight of CH4
      much4=16.0426
C Molecular weight
      mi=gamma0*much4+(1.-gamma0)*mun2 
C     molwgt=(mun2*much4)/(mun2+much4)
C     mi=molwgt
C probability that a V-T transition occurs (3.3	microns)
      p10=1.d-6
C       p10=7.d-3
C probability that a V-T transition occurs (7.6 microns)
      p10l=p10
C       p10l=3.d-6
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
C      a10l=2.12
      c10=p10*z10*0.1/kb/100.
      epsln=c10/a10
C epsilon at 3.3 micron
      eps=epsln/(1.+epsln)
      c10l=p10l*z10*0.1/kb/100.
      epslnl=c10l/a10l
C epsilon at 7.6 micron
      epsl=epslnl/(1.+epslnl)

C surface radius (m)
      rs=rSphere

C solar forcing at 3.3 um
      IF (lsIsConstant .EQ. 0) THEN
       lssol=lsStart
      ELSE
       lssol=lsperpetual
      ENDIF
      CALL GET_QO(i,j,bi,bj,lssol,0.,piF0,0.,piF1)
      piF=piF1*(wl/wn)


C timestep in numerical scheme (s)
C      deltaTp=1.E6
C      deltaTp=1.E8
       deltaTp=5.E7
C tolerance for numerical convergence (K)
      error=0.0001
C initialize the number of iterations completed
      niters=0
C maximum number of iterations before stopping
      nmax=500000


C Set up the vertical grid
C rm is the radius in m
C drmm is delta r-minus from Zalucha et al. 2011
C drmp is delta r-plus from Zalucha et al. 2011
C drmb is delta r-bar from Zalucha et al. 2011
      rm(1)=rs
      rm(2)=rs+300.
      drmm(1)=0.
      drmm(2)=300.
      drmp(1)=300.

      DO k=2,nlevs-1
       drmp(k)=1.03*drmm(k)
       rm(k+1)=rm(k)+drmp(k)
       drmm(k+1)=drmp(k)
       drmb(k)=(rm(k+1)-rm(k-1))/2.
      ENDDO

C altitude
      DO k=1,nlevs
       zm(k)=rm(k)-rs
      ENDDO

C initialize temperatures.  Here level 1 is surface temperature,
C everywhere else is a constant
      temp(1)=ts
      DO k=2,nlevs
       temp(k)=80.
C       temp(k)=40.
      ENDDO

      flag=1

C begin implicit time stepping
C      IF (piF .GT. 0.) THEN
       IF (piF .GT. 1.E-10) THEN
       DO WHILE(flag .GT. 0)
        DO k=1,nlevs
         xi(k)=temp(k)**(alpha+1.)
         tempPrev(k)=temp(k)
         xiPrev(k)=xi(k)
        ENDDO

C get pressure from integrating hydrostatic balance equation
C variable temperature
        pint(1)=0.
        DO k=2,nlevs
C         pint(k)=pint(k-1)-gravity/atm_Rd/temp(k)*drmm(k)
         pint(k)=pint(k-1)-bigGrav*plutoMass/rm(k)**2./atm_Rd/
     &            temp(k)*drmm(k)
        ENDDO
        DO k=1,nlevs
         p(k)=Ro_Surf(i,j,bi,bj)*DEXP(pint(k))
C number density of atmosphere
         na(k)=p(k)/kb/temp(k)
C number density of CH4
         n(k)=gamma0*na(k)
C mass density of atmosphere
         rho(k)=mi*amu*na(k)
C heating rate (J/m^3/s)
         q(k)=piF*eps*n(k)*bs
C cooling rate (J/m^3/s)
         l(k)=g*epsl*n(k)*hplanck*nul*a10l*
     &      DEXP(-hplanck*nul/kb/temp(k))
C net heating rate
         rnet(k)=q(k)-l(k)
        ENDDO

        IF(plutoTrop) THEN
         DO k=1,nlevs
C set up for dxi/dr coefficients
          weigt1(k)=( 1.0-tanh((zm(k)/1000.-hc)/5.) )/2.
          xfac1=40.*weigt1(k)
          xfac2=1.5+30.*exp(-zm(k)/1000./10.)
          kt(k)=1.0 D-9
          kc(k)=xfac1
C       	    kt(k)=1.0 _d -50
C       	    kc(k)=1.0 _d -50
          bigPi(k)=(atm_po/p(k))**atm_kappa
          bigPsi(k)=bigG1*EXP(-bigG2/temp(k))/(p(k)*temp(k))
          bigTheta(k)=EXP(-bigPsi(k))*
     &    (1.+bigPsi(k)*(1.-bigG2/temp(k)))
         ENDDO
C derivatives
	 DO k=1,nlevs-1
          dktdz(k)=(kt(k+1)-kt(k))/drmp(k)
          dbigPidz(k)=(bigPi(k+1)-bigPi(k))/drmp(k)
          dtempdz(k)=(temp(k+1)-temp(k))/drmp(k)
         ENDDO
	 dktdz(nlevs)=dktdz(nlevs-1)
         dbigPidz(nlevs)=dbigPidz(nlevs-1)
         dtempdz(nlevs)=dtempdz(nlevs-1)
C more derivatives
         DO k=1,nlevs-1
          dkp1=kc(k+1)*bigPi(k+1)*bigTheta(k+1)/
     &        ((alpha+1.)*temp(k+1)**alpha)
          dkp=kc(k)*bigPi(k)*bigTheta(k)/
     &        ((alpha+1.)*temp(k)**alpha)
          dinphi2(k)=(dkp1-dkp)/drmp(k)
          dkp1=kc(k+1)*EXP(-bigPsi(k+1))*dbigPidz(k+1)/
     &         temp(k+1)**alpha
          dkp=kc(k)*EXP(-bigPsi(k))*dbigPidz(k)/
     &         temp(k)**alpha
          dinphi3(k)=(dkp1-dkp)/drmp(k)
         ENDDO
         dinphi2(nlevs)=dinphi2(nlevs-1)
         dinphi3(nlevs)=dinphi3(nlevs-1)
        ENDIF

CCCCC
        DO k=1,nlevs
C        xi(k)=temp(k)**(alpha+1.)
C leading coefficients in dxi/dr, see Zalucha et al. 2011a
C note Zalucha et al. 2011b uses different notation
C S <=> R*_net+C*
C B <=> Phi3
C C <=> Phi2
C G <=> Phi1
C
         bigG(k)=ko*temp(k)**alpha/atm_cp/rho(k)
         bigC(k)=temp(k)**alpha*2.*ko/atm_cp/rho(k)/rm(k)
         bigS(k)=temp(k)**alpha/atm_cp/rho(k)*(1.+alpha)*rnet(k)
         bigB(k)=0.
         IF(radcode .EQ. 7) THEN
          IF(k .LT. 90) THEN
          bigG(k)=bigG(k)+kt(k)+kc(k)
          bigC(k)=bigC(k)+
     &     dktdz(k)+(alpha+1.)*temp(k)**alpha/bigPi(k)/bigTheta(k)*
     &     (kc(k)*EXP(-bigPsi(k))/(temp(k)**alpha)*dbigPidz(k)+
     &     dinphi2(k))
          bigB(k)=bigB(k)+
     &      (alpha+1.)*temp(k)**alpha/bigPi(k)/bigTheta(k)*
     &      dinphi3(k)
          bigS(k)=bigS(k)
     &     -alpha*(alpha+1.)*temp(k)**(alpha-1.)*kt(k)
     &     *dtempdz(k)**2.
          ENDIF
         ENDIF
	ENDDO


C coefficients in xi, see Zalucha et al. 2011
        alpham(1)=0.
        lambdam(1)=1.
        omegam(1)=0.
        qm(1)=temp(1)**(alpha+1.)

        alpham(nlevs)=-1.
        lambdam(nlevs)=1.
        omegam(nlevs)=0.
        qm(nlevs)=0.

        DO k=2,nlevs-1
         alpham(k)=deltaTp/drmb(k)*(-bigG(k)/drmm(k)+bigC(k)/2.)
         lambdam(k)=1.-deltaTp*bigB(k)+
     &       2.*deltaTp*bigG(k)/drmp(k)/drmm(k)
         omegam(k)=-deltaTp/drmb(k)*(bigG(k)/drmp(k)+bigC(k)/2.)
         qm(k)=xi(k)+deltaTp*bigS(k)
        ENDDO


C Call numerical recipes in fortran tridiagonal matrix solver
        CALL TRIDAG(alpham,lambdam,omegam,qm,um,nlevs)
 
        flag=0

C Convert back to temperature.  Check for convergence.
        DO k=1,nlevs
         temp(k)=um(k)**(1./(alpha+1.))
         IF(DABS(um(k)-xiPrev(k)) .GT. error) flag=flag+1
        ENDDO

        niters=niters+1
        IF(niters .GT. nmax) THEN
         flag=0
         write(*,*) 'maxiters reached'
         write(*,*) 'bj=',bj,', j=',j
        ENDIF
C end of do while loop
       ENDDO

C Interpolate to GCM rC grid


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
         plutoT(k)=tempout
        ENDIF
       ENDDO        
C goes with piF .eq. 0 loop
      ELSE
       DO k=1,Nr
        plutoT(k)=1. _d -8
       ENDDO
      ENDIF

      RETURN
      END       

