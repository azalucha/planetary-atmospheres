
CCC   Rediative equilibrium temperature for Pluto's atmosphere
C
C     PLUTOT,   dtdz2,    ZTPR,     IMDIFF,   tridag,   
C     HTCH4,    HTCH4n,   INPUT4,   HTCH5,    HTCH5n,   INPUT5, 
C     QTCH4Z,   QJCH4X,   GAMM12,   INPUT6,   ZKINTE,  
C     QTCOZ,    QJCOX,    ALPDCO,   STRT7,    GAMX2
C     ENZD2,    BLACM,    POLINT,   INVERT,   
c
c     HTCH4M,   HTCH5M,   
C
CCC   Higher resolution to accurately accomadate the tropospheric K_zz parameterization

 
C ==================================================================
c
C  To calculate the vertical derivative (dtdz) of temperature (t) for uneven grids (zm)
c
      subroutine dtdz2(t,zm,dtdz,km)
      implicit real*8(a-h,o-z)
      dimension t(km),zm(km),dtdz(km)
      kmm=km-1
      do k=2,kmm
      xfac1=(zm(k)-zm(k-1))**2
      xfac2=(zm(k+1)-zm(k))**2
      xfac3=(zm(k+1)-zm(k))*(zm(k)-zm(k-1))*(zm(k+1)-zm(k-1))
      dtdz(k)=((t(k+1)-t(k))*xfac1+(t(k)-t(k-1))*xfac2)/xfac3
      enddo
        dtdz(1)=(t(2)-t(1))/(zm(2)-zm(1))
        dtdz(km)=(t(km)-t(kmm))/(zm(km)-zm(kmm))
      return
      end

C ==================================================================
c
C  To calculate p and [rho] for given temperature profile. ZM(1)=0=furface
C  MODE=1 saturated surface pressure; MODE=2 fixed surface number density
c
      SUBROUTINE ZTPRN(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00)
      IMPLICIT real*8(A-H,O-Z)
      COMMON /PARAIR/ R00,WTMOL
      DIMENSION ZM(KM),TEM(KM),PRE(KM),RHO(KM)
     &  ,GRV(KM),DN(KM),BN(KM),HSCAL(KM)
      BC=1.3804D-23          ! Boltzmann constant: J K
      GG0=0.637D0             ! Surface gravity: m s**(-2)
      RGAS=8314.0D0/WTMOL            ! gas constant for N2: J K**(-1) kg**(-1)
      DO 10 K=1,KM
      GRV(K)=GG0*(R00/(R00+ZM(K)))**2   ! gravity
  10  CONTINUE
      HSCAL(1)=RGAS*TEM(1)/GG0
      PRE(1)=P00
      DO 20 K=2,KM
      HSCAL(K)=RGAS*(TEM(K)+TEM(K-1))/(GRV(K)+GRV(K-1))
      FAC2=(ZM(K)-ZM(K-1))/HSCAL(K)
      PRE(K)=PRE(K-1)*DEXP(-FAC2)
  20  CONTINUE
      DO 30 K=1,KM
      RHO(K)=PRE(K)/(RGAS*TEM(K))         ! mass density: kg m**(-3)
      DN(K)=PRE(K)/(BC*TEM(K))            ! number density: m**(-3)
  30  CONTINUE
      BN(KM)=DN(KM)*HSCAL(KM)             ! column number density: m**(-2)
      DO 40 K=2,KM
      IK=KM-K+1
      FACC1=(ZM(IK+1)-ZM(IK))*(DN(IK+1)+DN(IK))/2.0D0
      BN(IK)=BN(IK+1)+FACC1
  40  CONTINUE
      RETURN
      END


C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz) + BT.
C  Boundary condition: T(1)=lower boundary, T(KM)-T(KM-1)=DZ*FXTOP.
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDIFF2(T,QS,DS,CS,bs,KM,ZM,DT,FXTOP,RESIDT)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(KMX=140)                !  KMX = KM
      DIMENSION T(KM),QS(KM),DS(KM),CS(KM),bs(km),ZM(KM),RESIDT(KM)
      DIMENSION DZP(KMX),DZM(KMX),DZB(KMX)
     & ,AA(KMX),BB(KMX),CC(KMX),RR(KMX),WK1(KMX)
      KMM=KM-1
      DZP(1)=ZM(2)-ZM(1)
      DZM(1)=DZP(1)
      DZB(1)=DZP(1)
      DZM(KM)=ZM(KM)-ZM(KMM)
      DZP(KM)=DZM(KM)
      DZB(KM)=DZM(KM)
      DO 10 K=2,KMM
      DZP(K)=ZM(K+1)-ZM(K)
      DZM(K)=ZM(K)-ZM(K-1)
  10  DZB(K)=(ZM(K+1)-ZM(K-1))/2.0D0
      AA(1)=0.0D0
      BB(1)=1.0D0
      CC(1)=0.0D0
      RR(1)=T(1)
      AA(KM)=-1.0D0
      BB(KM)=1.0D0
      CC(KM)=0.0D0
      RR(KM)=DZM(KM)*FXTOP
      DO 20 K=1,KM
      RESIDT(K)=T(K)
  20  WK1(K)=DT/DZB(K)
      DO 30 K=2,KMM
      AA(K)=WK1(K)*(-DS(K)/DZM(K)+CS(K)/2.0D0)
c      BB(K)=1.0D0+WK1(K)*DS(K)*(1.0D0/DZP(K)+1.0D0/DZM(K))
      BB(K)=1.0D0-bs(k)*dt+WK1(K)*DS(K)*(1.0D0/DZP(K)+1.0D0/DZM(K))
      CC(K)=-WK1(K)*(DS(K)/DZP(K)+CS(K)/2.0D0)
  30  RR(K)=T(K)+DT*QS(K)
      CALL tridag(AA,BB,CC,RR,T,KM)
      DO 50 K=1,KM
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT-QS(K)
      RETURN
      END


C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz) + BT.
C  Boundary condition: T(2)-T(1)=DZ*FXbot, T(KM)-T(KM-1)=DZ*FXTOP.
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDjFF2(T,QS,DS,CS,bs,KM,ZM,DT,fxbot,FXTOP,RESIDT)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(KMX=140)                !  KMX = KM
      DIMENSION T(KM),QS(KM),DS(KM),CS(KM),bs(km),ZM(KM),RESIDT(KM)
      DIMENSION DZP(KMX),DZM(KMX),DZB(KMX)
     & ,AA(KMX),BB(KMX),CC(KMX),RR(KMX),WK1(KMX)
      KMM=KM-1
      DZP(1)=ZM(2)-ZM(1)
      DZM(1)=DZP(1)
      DZB(1)=DZP(1)
      DZM(KM)=ZM(KM)-ZM(KMM)
      DZP(KM)=DZM(KM)
      DZB(KM)=DZM(KM)
      DO 10 K=2,KMM
      DZP(K)=ZM(K+1)-ZM(K)
      DZM(K)=ZM(K)-ZM(K-1)
  10  DZB(K)=(ZM(K+1)-ZM(K-1))/2.0D0
      AA(1)=0.0D0
      BB(1)=-1.0D0
      CC(1)=1.0D0
      RR(1)=DZM(1)*fxbot
      AA(KM)=-1.0D0
      BB(KM)=1.0D0
      CC(KM)=0.0D0
      RR(KM)=DZM(KM)*FXTOP
      DO 20 K=1,KM
      RESIDT(K)=T(K)
  20  WK1(K)=DT/DZB(K)
      DO 30 K=2,KMM
      AA(K)=WK1(K)*(-DS(K)/DZM(K)+CS(K)/2.0D0)
c      BB(K)=1.0D0+WK1(K)*DS(K)*(1.0D0/DZP(K)+1.0D0/DZM(K))
      BB(K)=1.0D0-bs(k)*dt+WK1(K)*DS(K)*(1.0D0/DZP(K)+1.0D0/DZM(K))
      CC(K)=-WK1(K)*(DS(K)/DZP(K)+CS(K)/2.0D0)
  30  RR(K)=T(K)+DT*QS(K)
      CALL tridag(AA,BB,CC,RR,T,KM)
      DO 50 K=1,KM
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT-QS(K)
      RETURN
      END
     
      SUBROUTINE tridag(a,b,c,r,u,n)
C  To solve tridiagonal line set. REVISED FROM "Numerical Recipes" for REAL*8 
      INTEGER n,NMAX
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=200)
      INTEGER j
      REAL*8 bet,gam(NMAX)
      if(b(1).eq.0.0D0)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.0D0)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END

 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate solar near IR heating rate by CH4 3.3 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
C  To calculate the globally averaged solar uv heating rate in [K s^-1] by 
C  CH4 nu3 band at Pluto.  AN=CH4 number density (m^-3), BN=CH4 column 
C  number density (m^-2); TEM=temperature, PRE=pressure in Pa, 
C  assuming N2 dominant atmosphere, rnormf=normalized sun-planet distance
c

c2009.0318   With the lower CH4 mixing ratios, I think we need to include
c2009.0318   CH4 heating by absorption of solar 80-140 nm, principally Lyman alpha 
c2009.0318   at 121.6 nm, UV radiation.   The CH4 cross section at 121.6 nm is 1.8 
c2009.0318   x 10^-17 cm^2, so tau = 1 occurs below our upper boundary.  
c2009.0318   For solar medium conditions at 30 AU is the solar UV power flux is 
c2009.0318   0.0108 erg cm^-2 s^-1.  For CH4, the heating efficiency is 50%.  For a global 
c2009.0318   averaged heating rate, divide the flux by 2 and use solar zenith angle of 60 deg.  
c2009.0318   Thus q_CH4_UV = 0.0108 x 0.5 / 2 = 0.0027 erg cm^-2 s^-1 with tau evaluated 
c2009.0318   at 60 deg solar zenith angle.  For solar minimum conditions, 
c2009.0318   multiply by 0.6, and for solar maximum conditions multiply by 1.7.
c
c
      SUBROUTINE HTCH4tr(HT,AN,BN,TEM,PRE,KM,MODE,IYELLE,GAMN4,IZHU
     &  ,i,j,ii,jj,lssol,albedo,jf107)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4,KMAX=200,MM=4)                !  KMAX  >or= KM
      COMMON /PARAIR/ R00,WTMOL
      COMMON /CH4V3/ SIGX(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM)
     &   ,SIG1(JM,KMAX),GAMSL(KMAX),GAMN4(KM)
      DIMENSION XX4(MM),WW4(MM),WK1(KMAX),WK2(KMAX)
      IF(MODE.EQ.1) CALL INPUT4    ! read in CH4 spectral data the 1st time calling HTCH4
      PI=3.14159265D0
      XX4(1)=0.109063D0
      XX4(2)=0.518378D0
      XX4(3)=1.052419D0
      XX4(4)=1.461733D0
      WW4(1)=0.273205D0
      WW4(2)=0.512194D0
      WW4(3)=0.512194D0
      WW4(4)=0.273205D0
      DO 20 K=1,KM
      WK1(K)=0.0D0
  20  CONTINUE
      DO 40 I=1,MM
      DO 40 J=1,MM
      FAC1=DCOS(XX4(I))
      FAC2=DCOS(XX4(J))
      XMU1=FAC1*FAC2
      CALL HTCH4nt(HT,AN,BN,TEM,PRE,KM,XMU1,IYELLE,GAMN4,IZHU,
     &         i,j,ii,jj,lssol,albedo)
      DO 30 K=1,KM
      WK1(K)=WK1(K)+HT(K)*FAC2*WW4(I)*WW4(J)
  30  CONTINUE
  40  CONTINUE
      DO 60 K=1,KM
      FACT1=(2.0D0/PI)/2.0D0
      HT(K)=FACT1*WK1(K)     ! globally averaged solar heating rate in [K s^-1]
  60  CONTINUE
C      IF(IYELLE.EQ.1) THEN
C      DO 70 K=1,KM
C  70  HT(K)=HT(K)/2.0D0      ! a factor 2 difference between Zhu and Yelle
C      ENDIF
c
c  to add heating rate by Lyman alpha with only one solar zenith angle
c
      sig11=1.8d-21    !  CH4 cross section at 121.6 nm is 1.8 x 10^-17 m^2
      falp=0.0108d-3   !  For solar medium conditions at 30 AU, the solar UV power flux is 0.0108 erg cm^-2 s^-1
      if(jf107.le.0)  f11=0.0d0
      if(jf107.eq.1)  f11=falp*0.6d0
      if(jf107.eq.2)  f11=falp
      if(jf107.eq.3)  f11=falp*1.7d0
      xmu1=0.5d0
      xfac1=sig11*f11
      do k=1,km
      wk1(k)=an(k)*xfac1*dexp(-dmin1(sig11*bn(k),300.0d0)/xmu1)  ! heating rate: rho*c_p*dT/dt
      rhoairn2=pre(k)/(287.0d0*tem(k))
      wk2(k)=wk1(k)/(1004.0d0*rhoairn2)   ! heating rate of dT/dt in K/s
      enddo
        xfac2=0.5d0*0.5d0  ! (heating efficiency)*(global average), included in falp
        do k=1,km
        ht(k)=ht(k)+wk2(k)*xfac2  ! add the Lyman-alpha heating to the basic CH4 3.3 mu_bar heating
        enddo
      RETURN
      END

c
C  To calculate the solar uv heating rate in [K s^-1] by CH4 nu3 band at Pluto.
C  AN=CH4 number density (m^-3), BN=CH4 column number density (m^-3); 
C  TEM=temperature, PRE=pressure in Pa, assuming N2 dominant atmosphere
C  IYELLE=1: thin limit; IYELLE=2: thick limit;
c
      SUBROUTINE htch4nt(HT,AN,BN,TEM,PRE,KM,XMU,IYELLE,GAMN4,IZHU
     & ,i,j,ii,jj,lssol,albedo)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4,KMAX=200)                !  KMAX  >or= KM
      COMMON /PARAIR/ R00,WTMOL
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON /CH4V3/ SIGX(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMAX)
     &   ,EPS2(KMAX),EPS3(KMAX),GAMSL1(KMAX),GAMSL2(KMAX)
      DIMENSION WK1(LM),WK2(LM),GAMN4(KM)
      BK=1.38D-23      ! Boltzmann constant in J K^-1
      RCPR=7.0D0/2.0D0
      RAIR=8314.0D0/WTMOL      ! gas constant
      CP=RAIR*RCPR     ! specific heat at the constant pressure 
      FF0=9.10545D0    ! solar flux [pi B] at nu3 band in [W m^-2 / m^-1]
      EPS1=2.40232D-8  ! geometric factor for the solar constant in steradian
      SIGDNU=2.93D-24*4.0D4   ! band strength [sigma * delta nu] in [m]
      XKVT3=XP3*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
      DO 30 K=1,KM
      DO 20 J=1,JM
      IF(TEM(K).LE.TREF(1)) THEN
      SIG1(J,K)=SIGX(J,1)
      GO TO 20
      ENDIF
      IF(TEM(K).GE.TREF(LM)) THEN
      SIG1(J,K)=SIGX(J,LM)
      GO TO 20
      ENDIF
      DO 10 L=1,LM
  10  WK1(L)=DLOG(SIGX(J,L))
      CALL POLINT(TREF,WK1,LM,TEM(K),Y2,DY)
      SIG1(J,K)=DEXP(Y2)
  20  CONTINUE
  30  CONTINUE
      DO 40 K=1,KM
      IF(IYELLE.EQ.1) THEN
      GAMSL1(K)=1.0D0
      GAMSL2(K)=1.0D0
      GO TO 40
      ENDIF
      FAC1=0.0D0
      FAC2=0.0D0
      FAC3=0.0D0
      DO 35 J=1,JM
      FACT2=DMIN1(SIG1(J,K)*BN(K),300.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN(K)/XMU,300.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      IF(IYELLE.EQ.2) GAMSL1(K)=FAC2/FAC1  ! one-side heating/cooling 0.5D0*FAC2/FAC1
      IF(IYELLE.EQ.3) GAMSL1(K)=1.0D0      ! Direct flux only
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
      RNU34=0.13245D0
      DO 50 K=1,KM
      AIRN=PRE(K)/(BK*TEM(K))      ! air number density in m**(-3)
c
c  revised on 2009.02.18, beging
c
      PHIK3=2.8d-17*AIRN/XA3
      xfack4=1.9d-21*dexp((tem(k)-240.0d0)/105.0d0)
     &    *(1.0d0+10.0d0*AN(k)/airn)
      PHIK4=xfack4*AIRN/XA4
      PHIKS=1.16d-17*AN(K)/(2.0D0*XA4)
c
c  ending
c
ccc      PHIKS=XKVT3*AIRN/(2.0*XA4)
      EPS2(K)=PHIK3/(GAMSL1(K)+PHIK3)        ! non-LTE efficient factor
      FAC11=PHIKS/(1.0D0+PHIKS)   ! correction suggested by Roger Yelle
      FAC22=PHIK4/(GAMN4(K)+PHIK4)
      EPS3(K)=RNU34+(1.0D0-RNU34)*FAC11*FAC22
  50  CONTINUE
      DO 80 K=1,KM
      HT(K)=FF0*EPS1*EPS2(K)*AN(K)*SIGDNU*GAMSL2(K)   ! solar heating rate in [W m^-3]
      IF(IZHU.EQ.1) HT(K)=HT(K)*EPS3(K)        ! vibration-vibration cascade
      HT(K)=HT(K)/(RCPR*PRE(K)/TEM(K))         ! solar heating rate in [K s^-1]
  80  CONTINUE
c
      do k=1,km
      httemp0=ht(k)
      CALL GET_QO(i,j,ii,jj,lssol,0.,httemp0,
     &                            albedo,httemp)
      ht(k)=httemp
      enddo
c
      RETURN
      END

      SUBROUTINE INPUT4
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4)
      COMMON /CH4V3/ SIGX(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      OPEN(11,FILE='pluto_spect_dat/ggww_in1.dat',status='OLD')
      DO 10 J=1,JM
      READ(11,*) GG(J),WW(J)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_spect_dat/lbl_k2.dat',status='OLD')
      DO 20 J=1,JM
      READ(12,*) JX,GX,(SIGX(J,L),L=1,LM)
  20  CONTINUE
      CLOSE(UNIT=12)
      FACT1=2.6768D-27      ! conversion factor from cm^2/g to m^2 = 1.673e-24*16*1.0e-4
      DO 25 J=1,JM
      DO 25 L=1,LM
  25  SIGX(J,L)=SIGX(J,L)*FACT1            ! convert the cross section m**2
      DO 30 L=1,LM
      TREF(L)=40.0D0+20.0D0*DBLE(FLOAT(L-1))
  30  CONTINUE
      RETURN
      END
            
            
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate solar near IR heating rate by CH4 2.3 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  To calculate the globally averaged solar uv heating rate in [K s^-1] by 
C  CH4 nu3 band at Pluto.  AN=CH4 number density (m^-3), BN=CH4 column 
C  number density (m^-3); TEM=temperature, PRE=pressure in Pa, 
C  assuming N2 dominant atmosphere
c
      SUBROUTINE HTCH5t(HT,AN,BN,TEM,PRE,KM,MODE,GAMN4,
     &  i,j,ii,jj,lssol,albedo)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4,KMAX=200,MM=4)                !  KMAX  >or= KM
      COMMON /PARAIR/ R00,WTMOL
      COMMON /CH4V34/ SIGX5(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM)
     &   ,SIG1(JM,KMAX),EPS2(KMAX),GAMSL(KMAX)
      DIMENSION XX4(MM),WW4(MM),WK1(KMAX),WK2(KMAX)
      IF(MODE.EQ.1) CALL INPUT5    ! read in CH4 spectral data the 1st time calling HTCH5
      PI=3.14159265D0
      XX4(1)=0.109063D0
      XX4(2)=0.518378D0
      XX4(3)=1.052419D0
      XX4(4)=1.461733D0
      WW4(1)=0.273205D0
      WW4(2)=0.512194D0
      WW4(3)=0.512194D0
      WW4(4)=0.273205D0
      DO 20 K=1,KM
      WK1(K)=0.0D0
  20  CONTINUE
      DO 40 I=1,MM
      DO 40 J=1,MM
      FAC1=DCOS(XX4(I))
      FAC2=DCOS(XX4(J))
      XMU1=FAC1*FAC2
      CALL HTCH5nt(HT,AN,BN,TEM,PRE,KM,XMU1,GAMN4,
     &                      i,j,ii,jj,lssol,albedo)
      DO 30 K=1,KM
      WK1(K)=WK1(K)+HT(K)*FAC2*WW4(I)*WW4(J)
  30  CONTINUE
  40  CONTINUE
      DO 60 K=1,KM
      FACT1=(2.0D0/PI)/2.0D0
      HT(K)=FACT1*WK1(K)     ! globally averaged solar heating rate in [K s^-1]
  60  CONTINUE
      RETURN
      END


C  To calculate the solar uv heating rate in [K s^-1] by CH4 nu3 band at Pluto.
C  AN=CH4 number density (m^-3), BN=CH4 column number density (m^-3); 
C  TEM=temperature, PRE=pressure in Pa, assuming N2 dominant atmosphere
c
      SUBROUTINE htch5nt(HT,AN,BN,TEM,PRE,KM,XMU,GAMN4,
     &      i,j,ii,jj,lssol,albedo)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4,KMAX=200)                !  KMAX  >or=  KM
      COMMON /PARAIR/ R00,WTMOL
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON /CH4V34/ SIGX5(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMAX)
     &   ,EPS2(KMAX),EPS3(KMAX),GAMSL1(KMAX),GAMSL2(KMAX)
      DIMENSION WK1(LM),WK2(LM),GAMN4(KM)
      BK=1.38D-23      ! Boltzmann constant in J K^-1
      RCPR=7.0D0/2.0D0
      RAIR=8314.0D0/WTMOL      ! gas constant
      CP=RAIR*RCPR     ! specific heat at the constant pressure 
      FF0=15.5335D0    ! solar flux [pi B] at [nu3 + nu4] band in [W m^-2 / m^-1]
      EPS1=2.40232D-8  ! geometric factor for the solar constant in steradian
      SIGDNU=7.58D5*2.677D-26   ! band strength [sigma * delta nu] in [m] at 40 K
      XKVT3=XP2*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
      DO 30 K=1,KM
      DO 20 J=1,JM
      IF(TEM(K).LE.TREF(1)) THEN
      SIG1(J,K)=SIGX5(J,1)
      GO TO 20
      ENDIF
      IF(TEM(K).GE.TREF(LM)) THEN
      SIG1(J,K)=SIGX5(J,LM)
      GO TO 20
      ENDIF
      DO 10 L=1,LM
  10  WK1(L)=DLOG(SIGX5(J,L))
      CALL POLINT(TREF,WK1,LM,TEM(K),Y2,DY)
      SIG1(J,K)=DEXP(Y2)
  20  CONTINUE
  30  CONTINUE
      DO 40 K=1,KM
      FAC1=0.0D0
      FAC2=0.0D0
      FAC3=0.0D0
      DO 35 J=1,JM
      FACT2=DMIN1(SIG1(J,K)*BN(K),300.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN(K)/XMU,300.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      GAMSL1(K)=FAC2/FAC1    ! one-side heating/cooling 0.5D0*FAC2/FAC1
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
      RNU34=0.0903D0     ! = [4320 - 3*1310]/4320
      DO 50 K=1,KM
      AIRN=PRE(K)/(BK*TEM(K))              ! air number density in m**(-3)
c
c  revised on 2009.02.18, beging
c
      PHIK3=2.8d-17*AIRN/XA3
      xfack4=1.9d-21*dexp((tem(k)-240.0d0)/105.0d0)
     &    *(1.0d0+10.0d0*AN(k)/airn)
      PHIK4=xfack4*AIRN/XA4
      PHIKS2=1.16d-17*AN(K)/(2.0D0*XA4)
      PHIKS3=1.16d-17*AN(K)/(3.0D0*XA4)
c
c  ending
c
ccc      PHIKS=XKVT3*AIRN/(2.0*XA4)
c
      XRNUX=(1.0D0-RNU34)/3.0D0                ! testingxz
      FAC11=(1.0D0/(1.0D0+1.0D0/PHIKS3))*
     &  (1.0D0+2.0D0/(1.0D0+1.0D0/PHIKS2))*XRNUX
      EPS2(K)=PHIK3/(GAMSL1(K)+PHIK3)        ! non-LTE efficient factor
      FAC22=PHIK4/(GAMN4(K)+PHIK4)
CC      EPS3(K)=RNU34+(1.0D0-RNU34)*FAC11*FAC22
      EPS3(K)=RNU34+XRNUX*FAC11*FAC22
  50  CONTINUE
      DO 80 K=1,KM
      HT(K)=FF0*EPS1*EPS2(K)*AN(K)*SIGDNU*GAMSL2(K)   ! solar heating rate in [W m^-3]
      HT(K)=HT(K)*EPS3(K)                             ! vibration-vibration cascade
      HT(K)=HT(K)/(RCPR*PRE(K)/TEM(K))                ! solar heating rate in [K s^-1]
  80  CONTINUE
c

      do k=1,km
      httemp0=ht(k)
      CALL GET_QO(i,j,ii,jj,lssol,0.,httemp0,
     &                            albedo,httemp)
      ht(k)=httemp
      enddo
c
      RETURN
      END
  
  
      SUBROUTINE INPUT5
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4)
      COMMON /CH4V34/ SIGX5(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      OPEN(11,FILE='pluto_spect_dat/ggww_in1.dat',status='old')
      DO 10 J=1,JM
      READ(11,*) GG(J),WW(J)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_spect_dat/lbl_ks2.dat',status='old')
      DO 20 J=1,JM
      READ(12,*) JX,GX,(SIGX5(J,L),L=1,LM)
  20  CONTINUE
      CLOSE(UNIT=12)
      FACT1=2.6768D-27      ! conversion factor from cm^2/g to m^2 = 1.673e-24*16*1.0e-4
      DO 25 J=1,JM
      DO 25 L=1,LM
  25  SIGX5(J,L)=SIGX5(J,L)*FACT1            ! convert the cross section m**2
      DO 30 L=1,LM
      TREF(L)=40.0D0+20.0D0*DBLE(FLOAT(L-1))
  30  CONTINUE
      RETURN
      END
            
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate non-LTE cooling rate (K/day) by CH4 7.6 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
       SUBROUTINE QTCH4Z(TEM,PRE,RHO,QT,KM,MODE,GAMN4,IYELLE)
C  To calculate the cooling rate QJ [K s^-1] for CH4 7.6 micron band (K/s)
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(MCH4=2,KMX=140)                     ! KMX = KM  *********************
      COMMON /FRENG4/ VMCH4(MCH4),STRCH4(MCH4),ENECH4(MCH4)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QT(KM)
      DIMENSION WK1(KMX),WK2(KMX),GAMN4(KM)
      IF(MODE.EQ.1) CALL INPUT6    ! read in CH4 spectral data the 1st time calling QTCH4Z
      VMCH4(1)=1157.5D2            ! band center in m^-1
      VMCH4(2)=1292.5D2            ! band center in m^-1
      STRCH4(1)=6.106D1*10.0D0     ! band strength in m^-1 / (kg m^-2)
      STRCH4(2)=1.958D5*10.0D0     ! band strength in m^-1 / (kg m^-2)
      DO 20 K=1,KM
      QT(K)=0.0D0
      GAMN4(K)=0.0D0
  20  CONTINUE
      STRTOT=STRCH4(1)+STRCH4(2)
      DO 50 M=1,MCH4
      VMJ=VMCH4(M)
      STRJ=STRCH4(M)
      CALL QJCH4X(TEM,PRE,RHO,WK1,KM,VMJ,STRJ,M,WK2,IYELLE)
      DO 40 K=1,KM
      QT(K)=QT(K)-WK1(K)      ! total CH4 cooling rate in [K s^-1]
      GAMN4(K)=GAMN4(K)+STRJ*WK2(K)/STRTOT
  40  CONTINUE
  50  CONTINUE  
      RETURN
      END

      SUBROUTINE QJCH4X(TEM,PRE,RHO,QJ,KM,VMJ,STRJ,ICON
     &    ,GAMN4,IYELLE)
C  To calculate the cooling rate QJ [K s^-1] for CH4 7.6 micron band
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
C  VMJ [m^-1] is the frequency, (STRJ) is the band strength
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(KMX=140,JM=30,LM=4,KMP=KMX+1)          ! KMX = KM  *********************
      COMMON /PARAIR/ R00,WTMOL
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON /CH4V4/ SIG4A(JM,LM),SIG4B(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QJ(KM),GAMN4(KM)
      DIMENSION BLACX(KMX),H00(KMX),XKZ(KMX,JM),GAM1(KMX,KMX)
     & ,CURT(KMX,KMX),EE(KMX,KMX),CURTN(KMX,KMX),CINVT(KMX,KMP)
      PI=3.14159265D0
      GRAV=0.637D0               ! gravitational acceleration in m s^-2
      RAIR=8314.0D0/WTMOL      ! gas constant
      CP=RAIR*7.0D0/2.0D0       ! specific heat at the constant pressure       
      DO 10 K=1,KM
      BLACX(K)=BLACM(VMJ,TEM(K))      ! Planck function in [W/(m*m) * (m^-1)^-1 * Sr^-1]
  10  CONTINUE
      DO 20 J=1,JM
      DO 20 K=1,KM
      XKZ(K,J)=ZKINTE(TEM(K),J,ICON)
  20  CONTINUE
      DO 30 K=1,KM
      DO 30 J=1,KM
      GAM1(J,K)=GAMM12(J,K,PRE,RHO,KM,XKZ,WW,JM)
  30  CONTINUE
      DO 31 K=1,KM
      GAMN4(K)=GAM1(K,KM)   ! escape to space probability used in nu3 heating rate calculation.
  31  CONTINUE
C  To calculate the Curtis matrices +++++++++++++++++++++++++++++++++++++++++
      DO 35 L2=1,KM
      DO 35 L1=1,KM
      CURT(L1,L2)=0.0D0
  35  CONTINUE
      DO 40 L1=1,KM
      H00(L1)=2.0D0*PI*RHO(L1)*STRJ/CP
      CURT(L1,L1)=-2.0*GAM1(L1,L1)
  40  CONTINUE
      KMM=KM-1
      DO 45 L1=1,KM
      CURT(L1,1)=CURT(L1,1)+GAM1(L1,1)
      DO 42 L2=1,KMM
  42  CURT(L1,L2)=CURT(L1,L2)+0.5D0*DABS(GAM1(L1,L2)-GAM1(L1,L2+1))
      DO 44 L2=2,KM
  44  CURT(L1,L2)=CURT(L1,L2)+0.5D0*DABS(GAM1(L1,L2-1)-GAM1(L1,L2))
  45  CONTINUE
      DO 48 L2=1,KM
      DO 48 L1=1,KM
      CURT(L1,L2)=H00(L1)*CURT(L1,L2)
  48  CONTINUE
C  To calculate the non-LTE Curtis matrices ++++++++++++++++++++++++++++++++++++++
      DO 143 L2=1,KM
      DO 143 L1=1,KM
 143  EE(L1,L2)=0.0D0
      BK=1.38D-23                  ! Boltzmann constant in J K^-1
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
      DO 145 K=1,KM
      AIRN=PRE(K)/(BK*TEM(K))      ! air number density in m**(-3)
      PHI21=XKVT4*AIRN/XA4
      EE(K,K)=1.0D0/(2.0D0*H00(K)*PHI21)
 145  CONTINUE
      IF(IYELLE.EQ.1) THEN
      DO 240 L1=1,KM
 240  QJ(L1)=-H00(L1)*BLACX(L1)/(1.0D0+H00(L1)*EE(L1,L1))
      RETURN
      ENDIF
      DO 148 L2=1,KM
      DO 148 L1=1,KM
      CINVT(L1,L2)=0.0D0
      IF(L1.EQ.L2) CINVT(L1,L2)=1.0D0
      DO 147 L3=1,KM
 147  CINVT(L1,L2)=CINVT(L1,L2)-CURT(L1,L3)*EE(L3,L2)
 148  CONTINUE
      CALL INVERT(CINVT,KM,KMP)
      DO 154 L2=1,KM
      DO 154 L1=1,KM
      CURTN(L1,L2)=0.0D0
      DO 153 L3=1,KM
 153  CURTN(L1,L2)=CURTN(L1,L2)+CINVT(L1,L3)*CURT(L3,L2)
 154  CONTINUE
      DO 54 L1=1,KM
      XX2=0.0D0
      DO 53 L2=1,KM
  53  XX2=XX2+CURTN(L1,L2)*BLACX(L2)
  54  QJ(L1)=XX2
      RETURN
      END


      FUNCTION GAMM12(I,J,PUS,RUS,KM,XKZ,WW,JM)
cc To calculate the escape function according to Zhu (1992, Eq. (13)). 
      IMPLICIT real*8(A-H,O-Z)
      DIMENSION PUS(KM),RUS(KM),XKZ(KM,JM),WW(JM)
      GRAV=0.637D0               ! gravitational acceleration in m s^-2
      IF(I.EQ.J) GO TO 100
      KA=MIN0(I,J)
      KB=MAX0(I,J)-1
      SUM1=0.0D0
      SUM2=0.0D0
      DO 50 M=1,JM
      DELU=0.0D0
      DO 10 K=KA,KB
      FAC1=(RUS(K)+RUS(K+1))*(PUS(K+1)-PUS(K))/GRAV
      FAC2=XKZ(K,M)+XKZ(K+1,M)
      DELU=DELU+FAC1*FAC2
  10  CONTINUE
      DELU=DABS(DELU)/4.0D0
      FAC5=WW(M)*XKZ(I,M)
      SUM1=SUM1+FAC5
      SUM2=SUM2+FAC5*ENZD2(DELU)
  50  CONTINUE
      GAMM12=SUM2/SUM1
      RETURN
 100  EPSM=1.5D0/8.0D0        ! smooth factor = thickness of isothermal layer
      IF(I.EQ.1) FAC1=EPSM*(PUS(2)-PUS(1))
      IF(I.EQ.KM) FAC1=EPSM*(PUS(KM)-PUS(KM-1))
      IF(I.NE.1.AND.I.NE.KM) FAC1=EPSM*(PUS(I+1)-PUS(I-1))/2.0D0
      FAC1=FAC1*RUS(I)/GRAV
      SUM1=0.0D0
      SUM2=0.0D0
      DO 150 M=1,JM
      DELU=DABS(FAC1)*XKZ(I,M)
      FAC5=WW(M)*XKZ(I,M)
      SUM1=SUM1+FAC5
      SUM2=SUM2+FAC5*ENZD2(DELU)
 150  CONTINUE
      GAMM12=SUM2/SUM1
      RETURN
      END
      
      SUBROUTINE INPUT6
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4)
      COMMON /CH4V4/ SIG4A(JM,LM),SIG4B(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      OPEN(11,FILE='pluto_spect_dat/ggww_in1.dat',status='old')
      DO 10 J=1,JM
      READ(11,*) GG(J),WW(J)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_spect_dat/lbl_ka2.dat',status='old')
      OPEN(13,FILE='pluto_spect_dat/lbl_kb2.dat',status='old')
      DO 20 J=1,JM
      READ(12,*) JX,GX,(SIG4A(J,L),L=1,LM)
      READ(13,*) JX,GX,(SIG4B(J,L),L=1,LM)
  20  CONTINUE
      CLOSE(UNIT=12)
      CLOSE(UNIT=13)
      FACT1=0.1D0      ! conversion factor from cm^2/g to m^2/kg
      DO 25 L=1,LM
      DO 25 J=1,JM
      SIG4A(J,L)=SIG4A(J,L)*FACT1
  25  SIG4B(J,L)=SIG4B(J,L)*FACT1
      DO 30 L=1,LM
      TREF(L)=40.0D0+20.0D0*DBLE(FLOAT(L-1))
  30  CONTINUE
      RETURN
      END
      
      FUNCTION ZKINTE(T,J,ICON)
C  To intepolate the k-coefficients to T. ICON=1: SIG4A; ICON=2: SIG4B.
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4)
      COMMON /CH4V4/ SIG4A(JM,LM),SIG4B(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION WR1(4),WR2(4)
      GO TO (100,200) ICON
 100  CONTINUE
      IF(T.LE.TREF(1)) THEN
      ZKINTE=SIG4A(J,1)
      RETURN
      ENDIF
      IF(T.GE.TREF(LM)) THEN
      ZKINTE=SIG4A(J,LM)
      RETURN
      ENDIF
      DO 20 L=1,LM
  20  WR2(L)=DLOG(SIG4A(J,L))
      CALL POLINT(TREF,WR2,LM,T,XKLOGG,DY)
      ZKINTE=DEXP(XKLOGG)
      RETURN
 200  CONTINUE
      IF(T.LE.TREF(1)) THEN
      ZKINTE=SIG4B(J,1)
      RETURN
      ENDIF
      IF(T.GE.TREF(LM)) THEN
      ZKINTE=SIG4B(J,LM)
      RETURN
      ENDIF
      DO 40 L=1,LM
  40  WR2(L)=DLOG(SIG4B(J,L))
      CALL POLINT(TREF,WR2,LM,T,XKLOGG,DY)
      ZKINTE=DEXP(XKLOGG)
      RETURN
      END
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate cooling rate (K/day) by CO rotational lines (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE QTCOZ(TEM,PRE,RHO,QT,KM,MODE)
C  To calculate the cooling rate QJ [K s^-1] for CO rotational lines
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(MCO=25,KMX=140)                     ! KMX = KM  *********************
      COMMON /GAM212/ X212(212),GAM212(212)
      COMMON /FRENG/ VMCO(MCO),STRCO(MCO),ENECO(MCO)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QT(KM),WK1(KMX)
      IF(MODE.EQ.1) THEN
      OPEN(11,FILE='pluto_spect_dat/gamma1.dat',status='old')
      DO 10 L=1,212
      READ(11,12) X212(L),GAM212(L)
  10  CONTINUE
  12  FORMAT(F8.2,F12.8)
      CLOSE(UNIT=11)
      OPEN(11,FILE='pluto_spect_dat/sanda_co.dat',status='old')
      DO 15 I=1,MCO
      READ(11,17) IX,IX,FREQ,STR,AX,HX,ENERGY,IX,IX,JX
      VMCO(I)=FREQ*100.0D0
      STRCO(I)=STR
      ENECO(I)=ENERGY
  15  CONTINUE
  17  FORMAT(I5,I2,F11.6,E10.2,E11.3,F6.3,F11.4,2X,2I3,I4)
      CLOSE(UNIT=11)
      ENDIF
      DO 20 K=1,KM
      QT(K)=0.0D0
  20  CONTINUE
      DO 50 M=1,MCO
      VMJ=VMCO(M)
      STRJ=STRCO(M)
      ENGJ=ENECO(M)
      CALL QJCOX(TEM,PRE,RHO,WK1,KM,VMJ,STRJ,ENGJ)
      DO 40 K=1,KM
      QT(K)=QT(K)-WK1(K)      ! total CO cooling rate in [K s^-1]
  40  CONTINUE
  50  CONTINUE  
      RETURN
      END

 
      SUBROUTINE QJCOX(TEM,PRE,RHO,QJ,KM,VMJ,STRJ,ENGJ)
C  To calculate the cooling rate QJ [K s^-1] for CO rotational line J --> J-1
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
C  VMJ [m^-1] is the frequency, (STRJ,ENGJ) is used to derive the line strength
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(KMX=140)                     ! KMX = KM  *********************
      COMMON /PARAIR/ R00,WTMOL
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QJ(KM)
      DIMENSION ALDX(KMX),STX(KMX),BLACX(KMX),H00(KMX)
     &  ,GAM1(KMX,KMX),CURT(KMX,KMX)
      PI=3.14159265D0
      GRAV=0.637D0               ! gravitational acceleration in m s^-2
      RAIR=8314.0D0/WTMOL      ! gas constant
      CP=RAIR*7.0D0/2.0D0       ! specific heat at the constant pressure
      VCMJ=VMJ/100.0D0
      DO 10 K=1,KM
      ALDX(K)=ALPDCO(VMJ,TEM(K))      ! Doppler half-width in m^-1
      STX(K)=STRT7(TEM(K),STRJ,ENGJ,VCMJ)  ! line strength in m^-1 / (kg m^-2)
      BLACX(K)=BLACM(VMJ,TEM(K))      ! Planck function in [W/(m*m) * (m^-1)^-1 * Sr^-1]
  10  CONTINUE
      KMM=KM-1
      DO 30 I=1,KM
      DO 30 J=1,KM
      IF(J.GT.I) GO TO 30        ! assuming isothermal for escape function
      IF(I.EQ.J) THEN
      IF(I.EQ.1) XXP=RHO(1)*0.25D0*
     &                (PRE(1)-PRE(2))*STX(1)/ALDX(1)
      IF(I.EQ.KM) XXP=RHO(KM)*0.25D0*
     &                (PRE(KMM)-PRE(KM))*STX(KM)/ALDX(KM)
      IF(I.NE.1.AND.I.NE.KM) XXP=RHO(I)*0.125D0*
     &                (PRE(I-1)-PRE(I+1))*STX(I)/ALDX(I)
      XXP=XXP/(1.7724539D0*GRAV)
      GO TO 25
      ENDIF
      XXP=0.0D0
      K1=MIN0(I,J)
      K2=MAX0(I,J)-1
      DO 20 K=K1,K2
      XXP=XXP+RHO(K)*(PRE(K)-PRE(K+1))*STX(K)/ALDX(K)
  20  CONTINUE
      XXP=XXP/(1.7724539D0*GRAV)
  25  CONTINUE
      GAM1(I,J)=GAMX2(XXP)
  30  CONTINUE
      DO 32 I=1,KMM
      IP1=I+1
      DO 32 J=IP1,KM
  32  GAM1(I,J)=GAM1(J,I)
C  To calculate the Curtis matrices +++++++++++++++++++++++++++++++++++++++++
      DO 40 L1=1,KM
      H00(L1)=2.0D0*PI*RHO(L1)*STX(L1)/CP
      DO 40 L2=1,KM
      CURT(L1,L2)=0.0D0
      IF(L1.EQ.L2) CURT(L1,L2)=-2.0*GAM1(L1,L2)
  40  CONTINUE
      DO 45 L1=1,KM
      CURT(L1,1)=CURT(L1,1)+GAM1(L1,1)
      DO 42 L2=1,KMM
  42  CURT(L1,L2)=CURT(L1,L2)+0.5D0*DABS(GAM1(L1,L2)-GAM1(L1,L2+1))
      DO 44 L2=2,KM
  44  CURT(L1,L2)=CURT(L1,L2)+0.5D0*DABS(GAM1(L1,L2-1)-GAM1(L1,L2))
  45  CONTINUE
      DO 48 L1=1,KM
      DO 48 L2=1,KM
      CURT(L1,L2)=H00(L1)*CURT(L1,L2)
  48  CONTINUE
      DO 54 L1=1,KM
      XX2=0.0
      DO 53 L2=1,KM
  53  XX2=XX2+CURT(L1,L2)*BLACX(L2)
  54  QJ(L1)=XX2
      RETURN
      END

      FUNCTION ALPDCO(VM,T)
C  To calculate the Doppler half-width in m^-1, VCM [m^-1] and
C  T [K] are the reference line center and temperature, respectively
      IMPLICIT real*8(A-H,O-Z)
      VX=VM/25.0D2
      TX=T/40.0D0
      ALPDCO=1.28016D-3*VX*DSQRT(TX)
      RETURN
      END
      
      FUNCTION STRT7(T,STR,ENG,VCMJ)
C  To calculate CO line line strength in m^-1 / (kg m^-2)
C  STR = line strength at 296 K in  cm^-1 /(molec cm^-2)
C  VCMJ = line center frequency in cm^-1
      IMPLICIT real*8(A-H,O-Z)
      TS=296.0D0
      XM=28.0D0                            ! molecular weight
      QR=T/TS                              ! rotational partition function
      FAC1=1.439D0*ENG*(T-TS)/(T*TS)
      FAC7=(1.0D0-DEXP(-1.439D0*VCMJ/T))   ! induced emission becomes important
     &    /(1.0D0-DEXP(-1.439D0*VCMJ/TS))  !   for rotational lines at lower T
      IF(FAC1.LE.0.0D0) FAC2=-DMIN1(200.0D0,-FAC1)
      IF(FAC1.GE.0.0D0) FAC2=DMIN1(200.0D0,FAC1)
      BFAC=DEXP(FAC2)
      STRT7=(STR/QR)*BFAC*FAC7
      FAC2=6.023D23/XM                     ! conversion factor of molec g^-1
      STRT7=STRT7*FAC2                     ! line strength in cm^-1 / (g cm^-2)
      STRT7=STRT7*10.0D0                   ! line strength in m^-1 / (kg m^-2)
      RETURN
      END

      FUNCTION GAMX2(X)
C  To calculate the escape function
      IMPLICIT real*8(A-H,O-Z)
      COMMON /GAM212/ X212(212),GAM212(212)
      DIMENSION XA4(4),YA4(4)
      IF(X.LE.1.0D-4) THEN
      GAMX2=1.0D0
      RETURN
      ENDIF
      IF(X.LE.1.0D0) THEN
      GAMX2=1.0D0+(DLOG(X)-0.672784335D0)*X/1.4142136D0
        IF(X.GT.0.05D0) THEN
        X2=X*X
        X3=X2*X
	GAMX2=GAMX2-X2/3.4641016D0+X3/24.0D0
	ENDIF
        IF(X.GT.0.5D0) THEN
         X4=X3*X
         X5=X4*X
         X6=X5*X
	GAMX2=GAMX2-X4/160.99689D0+X5/1175.7551D0+X6/9524.7047D0
	ENDIF
      RETURN
      ENDIF
      IF(X.GE.104.0D0) THEN
      FAC1=DLOG(X)
      FAC2=FAC1/X
      GAMX2=0.2820947D0/(X*DSQRT(FAC1))
      GAMX2=GAMX2-0.057D0*FAC2*FAC2
      RETURN
      ENDIF
      J=IDINT(X/0.5D0)
      DO 10 K=1,4
      XA4(K)=X212(J+K-2)
  10  YA4(K)=GAM212(J+K-2)
      CALL POLINT(XA4,YA4,4,X,YYX,DY)
      GAMX2=YYX
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Some common subroutines
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION ENZD2(Z)
CC
C  Calcualte exponential integral from polynomial and rational 
C  approximation E2(z)=[exp(-z)-zE1(z)]=exp(-z)[1-exp(z)zE1(z)]
CC
      IMPLICIT real*8(A-H,O-Z)
      IF(Z.LE.1.0D0) THEN
        IF(Z.LE.1.0D-50) THEN
	ENZD2=1.0D0
	RETURN
	ENDIF
      ENZ=-DLOG(Z)-0.57721566D0+Z*(0.99999193D0
     & -Z*(0.24991055D0-Z*(0.05519968D0
     & -Z*(0.00976004D0-Z*0.00107857D0))))
      ENZD2=DEXP(-Z)-Z*ENZ
      RETURN 
      ENDIF
      IF(Z.GE.120.0D0) THEN
      ENZD2=0.0D0
      RETURN
      ENDIF
      ENZ=(0.2677737343D0+Z*(8.6347608925D0
     &  +Z*(18.0590169730D0+Z*(8.5733287401D0
     &  +Z))))/(3.9584969228D0+Z*(21.0996530827D0
     &  +Z*(25.6329561486D0+Z*(9.5733223454D0+Z))))
      ENZD2=DEXP(-Z)*(1.0D0-ENZ)
      RETURN
      END


      FUNCTION BLACM(VM,T)
C  Planck function [W/(m*m) * (m^-1)^-1 * Sr^-1] at T[K] & wavenumber VM [m^-1]
      IMPLICIT real*8(A-H,O-Z)
      BLACM=1.191D-16*VM**3/(DEXP(1.438786D-2*VM/T)-1.0D0)
      RETURN
      END
      
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C  Polynomial interpolation from "Numerical Recipes"
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.0D0)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

  
      SUBROUTINE INVERT(D,N,M)
CC  
CC  To invert the matrix D,  M=N+1. Mth row can be any values. Double precision
CC
      IMPLICIT  real*8(A-H,O-Z)
      DIMENSION D(N,M)
      KK=0
      JJ=0
      DO 10 K=1,N
      DO 11 J=1,N
  11  D(J,M)=0.0D0
      D(K,M)=1.0D0
      JJ=KK+1
      LL=JJ
      KK=KK+1
  20  IF(DABS(D(JJ,KK)-1.0D-30)) 21,21,22
  21  JJ=JJ+1
      IF(JJ-N) 20,20,99
  99  WRITE(*,98)
  98  FORMAT('ERRORNEOUS INPUT')
      RETURN
  22  IF(LL-JJ) 23,24,23
  23  DO 25 MM=1,M
      DTEMP=D(JJ,MM)
      D(LL,MM)=D(JJ,MM)
  25  D(JJ,MM)=DTEMP
  24  DIV=D(K,K)
      DO 30 LJ=1,M
  30  D(K,LJ)=D(K,LJ)/DIV
      DO 12 I=1,N
      FAC=D(I,K)
      IF(I-K) 15,12,15
  15  DO 31 LJ=1,M
  31  D(I,LJ)=D(I,LJ)-FAC*D(K,LJ)
  12  CONTINUE
      DO 40 J=1,N
  40  D(J,K)=D(J,M)
  10  CONTINUE
      RETURN
      END
      
      
C ****************************  END OF THE PROGRAM  **********************
cc
cc
cc       old/backup subroutines for the paper of Strobel et al. (1996).
c
C  To calculate the solar uv heating rate in [K s^-1] by CH4 nu3 band at Pluto.
C  AN=CH4 number density (m^-3), BN=CH4 column number density (m^-3); 
C  TEM=temperature, PRE=pressure in Pa, assuming N2 dominant atmosphere
C  IYELLE=1: thin limit; IYELLE=2: thick limit;
c
      SUBROUTINE HTCH4M(HT,AN,BN,TEM,PRE,KM,XMU,IYELLE,GAMN4,IZHU)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4,KMAX=200)                !  KMAX  >or= KM
      COMMON /PARAIR/ R00,WTMOL
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON /CH4V3/ SIGX(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMAX)
     &   ,EPS2(KMAX),EPS3(KMAX),GAMSL1(KMAX),GAMSL2(KMAX)
      DIMENSION WK1(LM),WK2(LM),GAMN4(KM)
      BK=1.38D-23      ! Boltzmann constant in J K^-1
      RCPR=7.0D0/2.0D0
      RAIR=8314.0D0/WTMOL      ! gas constant
      CP=RAIR*RCPR     ! specific heat at the constant pressure 
      FF0=9.10545D0    ! solar flux [pi B] at nu3 band in [W m^-2 / m^-1]
      EPS1=2.40232D-8  ! geometric factor for the solar constant in steradian
      SIGDNU=2.93D-24*4.0D4   ! band strength [sigma * delta nu] in [m]
      XKVT3=XP3*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
      DO 30 K=1,KM
      DO 20 J=1,JM
      IF(TEM(K).LE.TREF(1)) THEN
      SIG1(J,K)=SIGX(J,1)
      GO TO 20
      ENDIF
      IF(TEM(K).GE.TREF(LM)) THEN
      SIG1(J,K)=SIGX(J,LM)
      GO TO 20
      ENDIF
      DO 10 L=1,LM
  10  WK1(L)=DLOG(SIGX(J,L))
      CALL POLINT(TREF,WK1,LM,TEM(K),Y2,DY)
      SIG1(J,K)=DEXP(Y2)
  20  CONTINUE
  30  CONTINUE
      DO 40 K=1,KM
      IF(IYELLE.EQ.1) THEN
      GAMSL1(K)=1.0D0
      GAMSL2(K)=1.0D0
      GO TO 40
      ENDIF
      FAC1=0.0D0
      FAC2=0.0D0
      FAC3=0.0D0
      DO 35 J=1,JM
      FACT2=DMIN1(SIG1(J,K)*BN(K),300.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN(K)/XMU,300.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      IF(IYELLE.EQ.2) GAMSL1(K)=FAC2/FAC1  ! one-side heating/cooling 0.5D0*FAC2/FAC1
      IF(IYELLE.EQ.3) GAMSL1(K)=1.0D0      ! Direct flux only
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
      RNU34=0.13245D0
      DO 50 K=1,KM
      AIRN=PRE(K)/(BK*TEM(K))      ! air number density in m**(-3)
cxz      XKVT3=(DMIN1(1.5D-2,AN(K)/AIRN))*ZKVT
      PHIK3=XKVT3*AIRN/XA3
      PHIK4=XKVT4*AIRN/XA4
ccc      PHIKS=XKVT3*AIRN/(2.0*XA4)
      PHIKS=5.0D-2*ZKVT*AN(K)/(2.0D0*XA4)     ! testingxz
      EPS2(K)=PHIK3/(GAMSL1(K)+PHIK3)        ! non-LTE efficient factor
      FAC11=PHIKS/(1.0D0+PHIKS)   ! correction suggested by Roger Yelle
      FAC22=PHIK4/(GAMN4(K)+PHIK4)
      EPS3(K)=RNU34+(1.0D0-RNU34)*FAC11*FAC22
  50  CONTINUE
      DO 80 K=1,KM
      HT(K)=FF0*EPS1*EPS2(K)*AN(K)*SIGDNU*GAMSL2(K)   ! solar heating rate in [W m^-3]
      IF(IZHU.EQ.1) HT(K)=HT(K)*EPS3(K)        ! vibration-vibration cascade
      HT(K)=HT(K)/(RCPR*PRE(K)/TEM(K))         ! solar heating rate in [K s^-1]
  80  CONTINUE
      RETURN
      END

 
      SUBROUTINE HTCH5M(HT,AN,BN,TEM,PRE,KM,XMU,GAMN4)
C  To calculate the solar uv heating rate in [K s^-1] by CH4 nu3 band at Pluto.
C  AN=CH4 number density (m^-3), BN=CH4 column number density (m^-3); 
C  TEM=temperature, PRE=pressure in Pa, assuming N2 dominant atmosphere
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(JM=30,LM=4,KMAX=200)                !  KMAX  >or=  KM
      COMMON /PARAIR/ R00,WTMOL
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON /CH4V34/ SIGX5(JM,LM)
      COMMON /CH4ALL/ GG(JM),WW(JM),TREF(LM)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMAX)
     &   ,EPS2(KMAX),EPS3(KMAX),GAMSL1(KMAX),GAMSL2(KMAX)
      DIMENSION WK1(LM),WK2(LM),GAMN4(KM)
      BK=1.38D-23      ! Boltzmann constant in J K^-1
      RCPR=7.0D0/2.0D0
      RAIR=8314.0D0/WTMOL      ! gas constant
      CP=RAIR*RCPR     ! specific heat at the constant pressure 
      FF0=15.5335D0    ! solar flux [pi B] at [nu3 + nu4] band in [W m^-2 / m^-1]
      EPS1=2.40232D-8  ! geometric factor for the solar constant in steradian
      SIGDNU=7.58D5*2.677D-26   ! band strength [sigma * delta nu] in [m] at 40 K
      XKVT3=XP2*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
      DO 30 K=1,KM
      DO 20 J=1,JM
      IF(TEM(K).LE.TREF(1)) THEN
      SIG1(J,K)=SIGX5(J,1)
      GO TO 20
      ENDIF
      IF(TEM(K).GE.TREF(LM)) THEN
      SIG1(J,K)=SIGX5(J,LM)
      GO TO 20
      ENDIF
      DO 10 L=1,LM
  10  WK1(L)=DLOG(SIGX5(J,L))
      CALL POLINT(TREF,WK1,LM,TEM(K),Y2,DY)
      SIG1(J,K)=DEXP(Y2)
  20  CONTINUE
  30  CONTINUE
      DO 40 K=1,KM
      FAC1=0.0D0
      FAC2=0.0D0
      FAC3=0.0D0
      DO 35 J=1,JM
      FACT2=DMIN1(SIG1(J,K)*BN(K),300.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN(K)/XMU,300.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      GAMSL1(K)=FAC2/FAC1    ! one-side heating/cooling 0.5D0*FAC2/FAC1
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
      RNU34=0.0903D0     ! = [4320 - 3*1310]/4320
      DO 50 K=1,KM
      AIRN=PRE(K)/(BK*TEM(K))              ! air number density in m**(-3)
cxz      XKVT3=(DMIN1(1.5D-2,AN(K)/AIRN))*ZKVT
      PHIK3=XKVT3*AIRN/XA2
      PHIK4=XKVT4*AIRN/XA4
CC      PHIKS=XKVT3*AIRN/(3.0*XA4)
CC      FAC11=PHIKS/(1.0D0+PHIKS)   ! correction suggested by Roger Yelle
      PHIKS2=5.0D-2*ZKVT*AN(K)/(2.0D0*XA4)     ! testingxz
      PHIKS3=5.0D-2*ZKVT*AN(K)/(3.0D0*XA4)     ! testingxz
      XRNUX=(1.0D0-RNU34)/3.0D0                ! testingxz
      FAC11=(1.0D0/(1.0D0+1.0D0/PHIKS3))*
     &  (1.0D0+2.0D0/(1.0D0+1.0D0/PHIKS2))*XRNUX
      EPS2(K)=PHIK3/(GAMSL1(K)+PHIK3)        ! non-LTE efficient factor
      FAC22=PHIK4/(GAMN4(K)+PHIK4)
CC      EPS3(K)=RNU34+(1.0D0-RNU34)*FAC11*FAC22
      EPS3(K)=RNU34+XRNUX*FAC11*FAC22
  50  CONTINUE
      DO 80 K=1,KM
      HT(K)=FF0*EPS1*EPS2(K)*AN(K)*SIGDNU*GAMSL2(K)   ! solar heating rate in [W m^-3]
      HT(K)=HT(K)*EPS3(K)                             ! vibration-vibration cascade
      HT(K)=HT(K)/(RCPR*PRE(K)/TEM(K))                ! solar heating rate in [K s^-1]
  80  CONTINUE
      RETURN
      END
 
C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz).
C  Boundary condition: T(1)=lower boundary, T(KM)-T(KM-1)=DZ*FXTOP.
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDIFF(T,QS,DS,CS,KM,ZM,DT,FXTOP,RESIDT)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(KMX=140)                !  KMX = KM
      DIMENSION T(KM),QS(KM),DS(KM),CS(KM),ZM(KM),RESIDT(KM)
      DIMENSION DZP(KMX),DZM(KMX),DZB(KMX)
     & ,AA(KMX),BB(KMX),CC(KMX),RR(KMX),WK1(KMX)
      KMM=KM-1
      DZP(1)=ZM(2)-ZM(1)
      DZM(1)=DZP(1)
      DZB(1)=DZP(1)
      DZM(KM)=ZM(KM)-ZM(KMM)
      DZP(KM)=DZM(KM)
      DZB(KM)=DZM(KM)
      DO 10 K=2,KMM
      DZP(K)=ZM(K+1)-ZM(K)
      DZM(K)=ZM(K)-ZM(K-1)
  10  DZB(K)=(ZM(K+1)-ZM(K-1))/2.0D0
      AA(1)=0.0D0
      BB(1)=1.0D0
      CC(1)=0.0D0
      RR(1)=T(1)
      AA(KM)=-1.0D0
      BB(KM)=1.0D0
      CC(KM)=0.0D0
      RR(KM)=DZM(KM)*FXTOP
      DO 20 K=1,KM
      RESIDT(K)=T(K)
  20  WK1(K)=DT/DZB(K)
      DO 30 K=2,KMM
      AA(K)=WK1(K)*(-DS(K)/DZM(K)+CS(K)/2.0D0)
      BB(K)=1.0D0+WK1(K)*DS(K)*(1.0D0/DZP(K)+1.0D0/DZM(K))
      CC(K)=-WK1(K)*(DS(K)/DZP(K)+CS(K)/2.0D0)
  30  RR(K)=T(K)+DT*QS(K)
      CALL tridag(AA,BB,CC,RR,T,KM)
      DO 50 K=1,KM
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT-QS(K)
      RETURN
      END


C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz).
C  Boundary condition: T(2)-T(1)=DZ*FXbot, T(KM)-T(KM-1)=DZ*FXTOP.
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDjFF(T,QS,DS,CS,KM,ZM,DT,fxbot,FXTOP,RESIDT)
      IMPLICIT real*8(A-H,O-Z)
      PARAMETER(KMX=140)                !  KMX = KM
      DIMENSION T(KM),QS(KM),DS(KM),CS(KM),ZM(KM),RESIDT(KM)
      DIMENSION DZP(KMX),DZM(KMX),DZB(KMX)
     & ,AA(KMX),BB(KMX),CC(KMX),RR(KMX),WK1(KMX)
      KMM=KM-1
      DZP(1)=ZM(2)-ZM(1)
      DZM(1)=DZP(1)
      DZB(1)=DZP(1)
      DZM(KM)=ZM(KM)-ZM(KMM)
      DZP(KM)=DZM(KM)
      DZB(KM)=DZM(KM)
      DO 10 K=2,KMM
      DZP(K)=ZM(K+1)-ZM(K)
      DZM(K)=ZM(K)-ZM(K-1)
  10  DZB(K)=(ZM(K+1)-ZM(K-1))/2.0D0
      AA(1)=0.0D0
      BB(1)=-1.0D0
      CC(1)=1.0D0
      RR(1)=DZM(1)*fxbot
      AA(KM)=-1.0D0
      BB(KM)=1.0D0
      CC(KM)=0.0D0
      RR(KM)=DZM(KM)*FXTOP
      DO 20 K=1,KM
      RESIDT(K)=T(K)
  20  WK1(K)=DT/DZB(K)
      DO 30 K=2,KMM
      AA(K)=WK1(K)*(-DS(K)/DZM(K)+CS(K)/2.0D0)
      BB(K)=1.0D0+WK1(K)*DS(K)*(1.0D0/DZP(K)+1.0D0/DZM(K))
      CC(K)=-WK1(K)*(DS(K)/DZP(K)+CS(K)/2.0D0)
  30  RR(K)=T(K)+DT*QS(K)
      CALL tridag(AA,BB,CC,RR,T,KM)
      DO 50 K=1,KM
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT-QS(K)
      RETURN
      END
     


 
