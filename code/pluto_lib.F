C ==================================================================
c
c  To integrate the column heating rate
c
c  Input:
c      zm = altitude profile in m
c      q  = heating rate in erg cm**(-3) s**(-1)
c
c  Output: qcgs_sum = column integrated heating rate in erg cm**(-2) s**(-1)
c
      real*8 function qcgs_sum(zm,q,km)
      implicit real*8(a-h,o-z), INTEGER(I-N)
      dimension zm(km),q(km)
      qcgs_sum=0.0d0
      kmm=km-1
        do k=1,kmm
        qcgs_sum=qcgs_sum+(q(k)+q(k+1))*(zm(k+1)-zm(k))*0.5d2
        enddo
      return
      end


C ==================================================================
c
c  To integrate the column heating rate
c
c  Input:
c      zm = altitude profile in m
c      q  = heating rate in erg cm**(-3) s**(-1)
c
c  Output: qcgs_sum_j = column integrated heating rate in erg cm**(-2) s**(-1)
c
      real*8 function qcgs_sum_j(zm,q,km,nm,j)
      implicit real*8(a-h,o-z), INTEGER(I-N)
      dimension zm(km),q(km,nm)
      qcgs_sum_j=0.0d0
      kmm=km-1
        do k=1,kmm
        qcgs_sum_j=qcgs_sum_j+(q(k,j)+q(k+1,j))*(zm(k+1)-zm(k))*0.5d2
        enddo
      return
      end

C ==================================================================
c
C  to calculate the Jeans escape flux
c
c  Input:
c      rm = radius: m
c      tem = temperature: K
c      dn = number density: m^-3
c
c  Output:
c      fcapj = F_j: m^-2 s^-1
c      hscal = scale height: m
c      eta = (scale height)/(mean free path)
c      xlambda = r/H = rm/hscal
c      capu = U: m s^-1
c      vesc = v_esc: m s^-1
c      capf_rho_g = global mass flux = F_rho_G: kg s^-1
c      frho = integration constant = F_rho: kg s^-1
c
c      cgamma = GAMMA = enhancement factor
c
      subroutine f_cap_j(rm,tem,dn,fcapj
     & ,hscal,eta,xlambda,capu,vesc,capf_rho_g,frho,cgamma,xk_bc,capg)
      implicit real*8(a-h,o-z), INTEGER(I-N)
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
C        xk_bc=1.3807d-23          ! Boltzmann constant (k): J K
C        capg=6.6726d-11           ! gravitational constant (G): m^3 s^-2 kg^-1
cc        capm=1.305d22             ! Pluto/Triton's mass (M): kg
        xm=4.6513d-26             ! N2 molecular mass (m): kg
        capq=8.791d-19            ! collisional cross section (Q): m^2
      hscal=xk_bc*tem*rm*rm/(capg*capm*xm)
      xlmf=1.0d0/(1.4142136d0*dn*capq)   ! mean free path (1/(sqrt(2)*n*Q): m
      eta=hscal/(1.4142136d0*xlmf)
      xlambda=rm/hscal
      capu=dsqrt(2.0d0*xk_bc*tem/xm)
      vesc=dsqrt(2.0d0*capg*capm/rm)
      fcapj=(dn*capu/3.5449077d0)*(1.0d0+xlambda)  ! 3.5449077=2*sqrt(pi)
     &       *dexp(-dmin1(xlambda,199.0d0))
      frho=rm*rm*fcapj*xm*cgamma
      capf_rho_g=frho*12.5663706d0                 ! 12.5663706=4*pi
      return
      end

 
C ==================================================================
c
C  To calculate the vertical derivative (dtdz) of temperature (t) for uneven grids (zm)
c
      subroutine dtdz2(t,zm,dtdz,km)
      implicit real*8(a-h,o-z), INTEGER(I-N)
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
C  to integrate static pressure equation to get ([rho], pre).
c
c
c  Input:
c      p00 = surface pressure in Pa (Pascal)
c      zm = altitude profile in m
c      tem = temperature profile in K
c
c  Output:
c      pre = pressure profile in Pa
c      rho = density profile in kg m^-3
c      grv = gravity in m s^-2
c      dn = number density profile in m^-3
c      bn = column number density profile in m^-2
c      hscal = scale height in m
c
      SUBROUTINE ZTPRN(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00,BC)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &              rsunp
      DIMENSION ZM(KM),TEM(KM),PRE(KM),RHO(KM)
     &  ,GRV(KM),DN(KM),BN(KM),HSCAL(KM)
C      BC=1.3807D-23          ! Boltzmann constant: J K
      DO 10 K=1,KM
      GRV(K)=g00*(R00/(R00+ZM(K)))**2   ! gravity
  10  CONTINUE
      HSCAL(1)=rair*TEM(1)/g00
      PRE(1)=P00
      DO 20 K=2,KM
      HSCAL(K)=rair*(TEM(K)+TEM(K-1))/(GRV(K)+GRV(K-1))
      FAC2=(ZM(K)-ZM(K-1))/HSCAL(K)
      PRE(K)=PRE(K-1)*DEXP(-FAC2)
  20  CONTINUE
      DO 30 K=1,KM
      RHO(K)=PRE(K)/(rair*TEM(K))         ! mass density: kg m**(-3)
      DN(K)=PRE(K)/(BC*TEM(K))            ! number density: m**(-3)
  30  CONTINUE
       rho(km)=rho(km-1)*(pre(km)/pre(km-1))  ! added on 2012.08.06 to count for
       pre(km)=rho(km)*(rair*tem(km))         ! the fact: pre->0 as z->infinity. 
       dn(km)=pre(km)/(bc*tem(km))   ! Model sets T~0.1 at km for Parker solution
      BN(KM)=DN(KM)*HSCAL(KM)    ! column number density: m**(-2)
      DO 40 K=2,KM
      IK=KM-K+1
      FACC1=(ZM(IK+1)-ZM(IK))*(DN(IK+1)+DN(IK))/2.0D0
      BN(IK)=BN(IK+1)+FACC1
  40  CONTINUE
      RETURN
      END


C ==================================================================
c
C  To calculate p, [rho] and w for given temperature profile. ZM(1)=0=furface
c
c  Input:
c      p00 = surface pressure in Pa (Pascal)
c      zm = altitude profile in m
c      tem = temperature profile in K
c      frho = input escape parameter in kg s^-1
c      w_in = vertical velocity profile at previous step in m s^-1
c
c  Output:
c      pre = pressure profile in Pa
c      rho = density profile in kg m^-3
c      dn = number density profile in m^-3
c      bn = column number density profile in m^-2
c      hscal = scale height in m
c      w_out = vertical velocity profile in m s^-1
c
      SUBROUTINE ZTPRN_w(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00
     & ,frho,w_in,w_out,BC)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      parameter(KMX=175)
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &              rsunp
      DIMENSION ZM(KM),TEM(KM),PRE(KM),RHO(KM)
     &  ,GRV(KM),DN(KM),BN(KM),HSCAL(KM),w_in(km),w_out(km)
      dimension rm(kmx),rmsq(kmx),wk1(kmx)
C      BC=1.3807d-23          ! Boltzmann constant: J K
      do k=1,km
      rm(k)=r00+zm(k)             ! radiu in m
      rmsq(k)=rm(k)**2
      grv(k)=g00*(r00/rm(k))**2   ! gravity
      enddo
c
      call dtdz2(w_in,zm,wk1,km)  ! dw/dr
c
      hscal(1)=rair*tem(1)/g00
      pre(1)=p00
      rho(1)=pre(1)/(rair*tem(1))
      w_out(1)=frho/(r00*r00*rho(1))
c
      do k=2,km
      hscal(k)=rair*(tem(k)+tem(k-1))/(grv(k)+grv(k-1))
      fac2=(zm(k)-zm(k-1))/hscal(k)
      fac3=frho*(wk1(k)+wk1(k-1))/(rmsq(k)+rmsq(k-1))/pre(k-1)
      pre(k)=pre(k-1)*dexp(-fac2-fac3)
      rho(k)=pre(k)/(rair*tem(k))         ! mass density: kg m**(-3)
      w_out(k)=frho/(rmsq(k)*rho(k))      ! vertical velocity in m/s
      enddo
c
c    iteration one or twice:
c
      do iterat=1,2
c
      call dtdz2(w_out,zm,wk1,km)  ! dw/dr
c
      do k=2,km
      hscal(k)=rair*(tem(k)+tem(k-1))/(grv(k)+grv(k-1))
      fac2=(zm(k)-zm(k-1))/hscal(k)
      fac3=frho*(wk1(k)+wk1(k-1))/(rmsq(k)+rmsq(k-1))/pre(k-1)
      pre(k)=pre(k-1)*dexp(-fac2-fac3)
      rho(k)=pre(k)/(rair*tem(k))         ! mass density: kg m**(-3)
      w_out(k)=frho/(rmsq(k)*rho(k))      ! vertical velocity in m/s
      enddo
       rho(km)=rho(km-1)*(pre(km)/pre(km-1))  ! added on 2012.08.06 to count for
       pre(km)=rho(km)*(rair*tem(km))         ! the fact: pre->0 as z->infinity. 
       w_out(km)=frho/(rmsq(km)*rho(km)) ! Model sets T~0.1 at km for Parker solution
c
      enddo
c
      do k=1,km
      dn(k)=pre(k)/(bc*tem(k))            ! number density: m**(-3)
      wk1(k)=rair*tem(k)
      rmsq(k)=w_out(k)**2
      enddo
c
      bn(km)=dn(km)*hscal(km)             ! column number density: m**(-2)
c
      do k=2,km
      IK=KM-K+1
      FACC1=(ZM(IK+1)-ZM(IK))*(DN(IK+1)+DN(IK))/2.0D0
      BN(IK)=BN(IK+1)+FACC1
      enddo
c
      RETURN
      END


C ==================================================================
c
C  To calculate w, [rho] and p for given temperature profile. ZM(1)=0=furface
c
c  Input:
c      p00 = surface pressure in Pa (Pascal)
c      zm = altitude profile in m
c      tem = temperature profile in K
c      frho_in = input escape parameter in kg s^-1
c      w_in = vertical velocity profile in m s^-1 for k_crit=0
c
c  Output:
c      pre = pressure profile in Pa
c      rho = density profile in kg m^-3
c      dn = number density profile in m^-3
c      bn = column number density profile in m^-2
c      hscal = scale height in m
c      frho_out = output escape parameter in kg s^-1
c      w_out = vertical velocity profile in m s^-1
c      k_crit = (0,kcrit) index of critical altitude (0=no critical velovity)
c
c  When jeans1<k_crit, it is a subsonic escape with (p,rho,w) directly solved 
c     from the momentum and continuity equations
c
c  When k_crit>=jeans1, it is a supersonic escape where w is first solved
c     according to Parker's arguments, and then (p,rho) are solved
c     from the momentum and continuity equations.
c
c
      SUBROUTINE zt_w_parker(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00
     & ,frho_in,frho_out,w_in,w_out,jeans1,k_crit,jjregime,bc)
      implicit real*8(a-h,o-z), INTEGER(I-N)
      parameter(KMX=175)
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &              rsunp
      dimension ZM(KM),TEM(KM),PRE(KM),RHO(KM)
     &  ,GRV(KM),DN(KM),BN(KM),HSCAL(KM),w_in(km),w_out(km)
      dimension rm(kmx),rmsq(kmx),wk1(kmx),wk2(kmx),wk3(kmx)
c
      if(jeans1.lt.km .and. jjregime.lt.3) then   ! standard subsonic Jeans escape
        call ZTPRN_w(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00
     &    ,frho_in,w_in,w_out,BC)
      k_crit=km  ! there is no sonic critical level
      return
      endif
c
C      BC=1.3807d-23          ! Boltzmann constant: J K
      do k=1,km
      rm(k)=r00+zm(k)             ! radiu in m
      rmsq(k)=rm(k)**2
      grv(k)=g00*(r00/rm(k))**2   ! gravity
c
      wk1(k)=grv(k)/(rmsq(k)*rair)   ! G*M/(r^4*R)
      wk3(k)=-tem(k)/rmsq(k)         ! -T/r^2
      enddo
c
      call dtdz2(wk3,zm,wk2,km)  ! d(-T/r^2)/dr
c
      k_crit=km
      do k=1,km
        if(wk1(k) .le. wk2(k)) then     ! lambda <= 1, critical level
        k_crit=k  ! weight at k_crit, 1-wt_crit=weight at k_crit-1:
        wt_crit=(wk1(k-1)-wk2(k-1))/(wk1(k-1)+wk2(k)-wk2(k-1)-wk1(k))
        goto 10
        endif
      enddo
c
  10  continue
c
      if(k_crit.eq.km .and. jjregime.lt.3) then   ! there is not critical level in the domain
        call ZTPRN_w(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00
     &    ,frho_in,w_in,w_out,BC)
        return
      endif
c
      if(jeans1.lt.k_crit .and. jjregime.lt.3) then   ! subsonic Jeans escape
        call ZTPRN_w(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00
     &    ,frho_in,w_in,w_out,BC)
      return
      endif
c
      if(jjregime.eq.2) then
      wmax11=-99.0d0
         do k=1,km
           if(wmax11.le.w_in(k)) then
             kmx11=k
             wmax11=w_in(k)
           endif
         enddo
      endif
c
      call ZTPRN(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00,BC)
c
      if(jjregime.eq.2) then
        frho_out=rmsq(kmx11)*rho(kmx11)*wmax11  ! F_rho(at dw/dr=0) = frho_out
cc or    frho_out=rmsq(km)*rho(km)*w_in(km)  ! F_rho(top) = frho_out  ! testingxz
      endif
c
      if(jjregime.eq.3) then
        w_crit=dsqrt(rair*tem(k_crit))*wt_crit    ! critical radial velocity
     &        +dsqrt(rair*tem(k_crit-1))*(1.0d0-wt_crit)
        frho_out=rmsq(k_crit)*rho(k_crit)*w_crit  ! F_rho at k_crit = frho_out
      endif
c
      do k=1,km
      wk3(k)=frho_out/(rmsq(k)*rho(k))    ! initial guess of radial velocity
      enddo
 
      call ZTPRN_w(ZM,TEM,PRE,RHO,GRV,DN,BN,HSCAL,KM,P00
     & ,frho_out,wk3,w_out,BC)
c
      return
      end


C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz) + BT.
C  Lower boundary condition: T(1)=lower boundary, 
C  Upper boundarr condition: dT/dz + a*T = F = FXTOP ==> C2*T(km-1)+C1*T(km)=C3 
C   with  C1=(1+a*DZ/2), C2=(-1+a*DZ/2), C3=DZ*FXTOP
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDIFF3(T,QS,DS,CS,bs,KM,ZM,DT,FXTOP,a,RESIDT)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175)                !  KMX = KM
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
c
      xff1=a*DZM(KM)/2.0d0
      c1=1.0d0+xff1
      c2=-1.0d0+xff1
      c3=DZM(KM)*FXTOP
c
      AA(1)=0.0D0
      BB(1)=1.0D0
      CC(1)=0.0D0
      RR(1)=T(1)
      AA(KM)=c2      ! -1.0D0
      BB(KM)=c1      !  1.0D0
      CC(KM)=0.0D0
      RR(KM)=c3      ! DZM(KM)*FXTOP
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
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT   ! -QS(K)
      RETURN
      END


C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz) + BT.
C  Boundary condition: T(1)=lower boundary, T(KM)-T(KM-1)=DZ*FXTOP.
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDIFF2(T,QS,DS,CS,bs,KM,ZM,DT,FXTOP,RESIDT)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175)                !  KMX = KM
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
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT   ! -QS(K)
      RETURN
      END


C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz) + BT.
C  Boundary condition: T(1)=lower boundary, T(KM)=tem_top.
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDIFF1(T,QS,DS,CS,bs,KM,ZM,DT,tem_top,RESIDT)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175)                !  KMX = KM
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
      AA(KM)=0.0D0
      BB(KM)=1.0D0
      CC(KM)=0.0D0
      RR(KM)=tem_top
      DO 20 K=1,KM
      RESIDT(K)=T(K)
  20  WK1(K)=DT/DZB(K)
      DO 30 K=2,KMM
      AA(K)=WK1(K)*(-DS(K)/DZM(K)+CS(K)/2.0D0)
      BB(K)=1.0D0-bs(k)*dt+WK1(K)*DS(K)*(1.0D0/DZP(K)+1.0D0/DZM(K))
      CC(K)=-WK1(K)*(DS(K)/DZP(K)+CS(K)/2.0D0)
  30  RR(K)=T(K)+DT*QS(K)
      CALL tridag(AA,BB,CC,RR,T,KM)
      DO 50 K=1,KM
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT   ! -QS(K)
      RETURN
      END


C  To solve the diffusion equation dT/dt = Q + D(d/dz)(dT/dz) + C(dT/dz) + BT.
C  Boundary condition: T(2)-T(1)=DZ*FXbot, T(KM)-T(KM-1)=DZ*FXTOP.
C  T is used for the temperatures at both t and t+dt.
C  FXTOP=temperature gradient at the top boundary, 
C  RESIDT=residual due to the diffusion processes: dT/dt - Q
c
      SUBROUTINE IMDjFF2(T,QS,DS,CS,bs,KM,ZM,DT,fxbot,FXTOP,RESIDT)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175)                !  KMX = KM
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
  50  RESIDT(K)=(T(K)-RESIDT(K))/DT   ! -QS(K)
      RETURN
      END
     
C  To solve tridiagonal line set. REVISED FROM "Numerical Recipes" for REAL*8 
c
C      SUBROUTINE tridag(a,b,c,r,u,n)
C      INTEGER n,kmx
C      REAL*8 a(n),b(n),c(n),r(n),u(n)
C      PARAMETER (KMX=175)
C      INTEGER j
C      REAL*8 bet,gam(kmx)
C      if(b(1).eq.0.0D0)pause 'tridag: rewrite equations'
C      bet=b(1)
C      u(1)=r(1)/bet
C      do 11 j=2,n
C        gam(j)=c(j-1)/bet
C        bet=b(j)-a(j)*gam(j)
C        if(bet.eq.0.0D0)pause 'tridag failed'
C        u(j)=(r(j)-a(j)*u(j-1))/bet
C11    continue
C      do 12 j=n-1,1,-1
C        u(j)=u(j)-gam(j+1)*u(j+1)
C12    continue
C      return
C      END

 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate solar near IR heating rate by CH4 3.3 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
C  To calculate the globally averaged solar uv heating rate in [K s^-1] by 
C  CH4 nu3 band at Pluto/Triton.  AN=CH4 number density (m^-3), BN=CH4 column 
C  number density (m^-2); TEM=temperature, PRE=pressure in Pa, 
C  assuming N2 dominant atmosphere, rnormf=normalized sun-planet distance
c
      SUBROUTINE HTCH4tr1(HT,AN,BN,TEM,PRE,KM,MODE,IYELLE,GAMN4,IZHU
     & ,rnormf,surface_albedo,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,MM=4,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &              rsunp
      COMMON/CH4V3/sigx(jm,lm,nm_p),strx_b1(lm),strx_b2(lm)     ! [nu 3], [2800, 3200] cm^-1, 3.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM)
     &   ,GAMN4(KM)
      DIMENSION XX4(MM),WW4(MM),WK1(KMX)
      dimension ht_a(kmx),ht_b(kmx)
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
      CALL HTCH4nt_r(ht_a,AN,BN,TEM,PRE,KM,XMU1,IYELLE,GAMN4
     &  ,IZHU,rnormf,xmu0,bk)
      CALL HTCH4nt_r_inv(ht_b,AN,BN,TEM,PRE,KM,XMU1,IYELLE,GAMN4
     &  ,IZHU,rnormf,surface_albedo,xmu0,bk)
      DO 30 K=1,KM
      ht(k)=ht_a(k)+ht_b(k)
      WK1(K)=WK1(K)+HT(K)*FAC2*WW4(I)*WW4(J)
  30  CONTINUE
  40  CONTINUE
      DO 60 K=1,KM
cc      FACT1=(2.0D0/PI)/2.0D0
cc      HT(K)=FACT1*WK1(K)     ! globally averaged solar heating rate in [K s^-1]
C AMZ I don't think this factor of pi should be there for individual columns
C      HT(K)=WK1(K)/pi     ! globally averaged solar heating rate in [K s^-1]
      HT(K)=WK1(K)
  60  CONTINUE
      IF(IYELLE.EQ.1) THEN
      DO 70 K=1,KM
  70  HT(K)=HT(K)/2.0D0      ! a factor 2 difference between Zhu and Yelle
      ENDIF
cc
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate solar near IR heating rate by CH4 1.67 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
C  To calculate the globally averaged solar uv heating rate in [K s^-1] by 
C  CH4 2*nu3 band at Pluto/Triton.  AN=CH4 number density (m^-3), BN=CH4 column 
C  number density (m^-2); TEM=temperature, PRE=pressure in Pa, 
C  assuming N2 dominant atmosphere, rnormf=normalized sun-planet distance
c
      SUBROUTINE HTCH4tr6(HT,AN,BN,TEM,PRE,KM,MODE,IYELLE,GAMN4,IZHU
     & ,rnormf,surface_albedo,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,MM=4,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &              rsunp
      COMMON/CH4V33/sigx6(jm,lm,nm_p),strx6_b1(lm),strx6_b2(lm)  ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM)
     &   ,GAMN4(KM)
      DIMENSION XX4(MM),WW4(MM),WK1(KMX)
      dimension ht_a(kmx),ht_b(kmx)
      IF(MODE.EQ.1) CALL INPUT33    ! read in CH4 spectral data the 1st time calling HTCH4
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
      ht(k)=0.0d0
  20  CONTINUE
      DO 40 I=1,MM
      DO 40 J=1,MM
      FAC1=DCOS(XX4(I))
      FAC2=DCOS(XX4(J))
      XMU1=FAC1*FAC2
      CALL HTCH6t(ht_a,AN,BN,TEM,PRE,KM,XMU1,IYELLE,GAMN4
     &  ,IZHU,rnormf,xmu0,bk)
      CALL HTCH6t_inv(ht_b,AN,BN,TEM,PRE,KM,XMU1,IYELLE,GAMN4
     &  ,IZHU,rnormf,surface_albedo,xmu0,bk)
      DO 30 K=1,KM
      ht(k)=ht_a(k)+ht_b(k)
      WK1(K)=WK1(K)+HT(K)*FAC2*WW4(I)*WW4(J)
  30  CONTINUE
  40  CONTINUE
      DO 60 K=1,KM
cc      FACT1=(2.0D0/PI)/2.0D0
cc      HT(K)=FACT1*WK1(K)     ! globally averaged solar heating rate in [K s^-1]
C      HT(K)=WK1(K)/pi     ! globally averaged solar heating rate in [K s^-1]
C AMZ I think in a single column the factor of pi shouldn't be there
      HT(K)=WK1(K)
  60  CONTINUE
      IF(IYELLE.EQ.1) THEN
      DO 70 K=1,KM
  70  HT(K)=HT(K)/2.0D0      ! a factor 2 difference between Zhu and Yelle
      ENDIF
cc
      RETURN
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate solar EUV and FUV heating rate by CH4 and N2 (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
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
c  jeuvfuv = (0,1,2,3) = (use f107 parameterization, case 1, case 2, case 3)
c
      SUBROUTINE HTCH4tr2(ht,ht_ch4,ht_n2,AN,BN,TEM,PRE,KM,nm3
     & ,rnormf,jeuvfuv,f107,cn,dn,j1_efuv,xmu0)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,MM=4,nm=3,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &              rsunp
      COMMON/CH4V3/sigx(jm,lm,nm_p),strx_b1(lm),strx_b2(lm)     ! [nu 3], [2800, 3200] cm^-1, 3.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      dimension ht_ch4(km,nm3),ht_n2(km,nm3)    ! nm3=nm
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),cn(km),dn(km)
      DIMENSION WK1(KMX),WK2(KMX)
      dimension fuv(nm),sig_ch4(nm),sig_n2(nm),sf_ch4(nm),sf_n2(nm)
      PI=3.14159265D0
      do k=1,km
        ht(k)=0.0d0
        do n=1,nm
        ht_ch4(k,n)=0.0d0
        ht_n2(k,n)=0.0d0
        enddo
      enddo
        if(j1_efuv.eq.0) return
c
c  The following default insertion of adding Lyman alpha heating is
c    replaced by a more accurate update of EUV and FUV spectra of
c    "Pluto_Solar_EUV_FUV_fluxes_xsec.xls" that includes absorption by
c    both CH4 and N2. Two new arrays of (cn,dn) are introduced to
c    represent the N2 (or air) number and column densities, respectively.
c    (2012.04.15)
c
c  Modified on 2012.08.22 to also include the Triton cases:
c  jeuvfuv=(1,2,3)=(Triton_2009_0704, Pluto_2014paper, Triton_1989_0808) with f107=(70,130,220)
c
cc
      xfsq=rsunp0**2  ! correction to reference sun-planet distance, rsunp0**2 = 900
c
      if(jeuvfuv.gt.0) then
c
        if(jeuvfuv.eq.1) then  ! Trion solar EUV and FUV heating for Voyager 2: July 8, 2009, f107=70
        fuv(1)=1.65d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (0-40) nm
        fuv(2)=0.74d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (40-100) nm
        fuv(3)=7.46d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (100-140) nm
          sig_ch4(1)=4.05d-22  ! CH4 cross section in m^2, (0-40) nm
          sig_ch4(2)=4.07d-21  ! CH4 cross section in m^2, (40-100) nm
          sig_ch4(3)=1.88d-21  ! CH4 cross section in m^2, (100-140) nm
            sig_n2(1)=8.27d-22  ! N2 cross section in m^2, (0-40) nm
            sig_n2(2)=2.26d-21  ! N2 cross section in m^2, (40-100) nm
            sig_n2(3)=0.00d-22  ! N2 cross section in m^2, (100-140) nm
        endif
c
        if(jeuvfuv.eq.2) then  ! Pluto solar EUV and FUV heating for Zhu et al (2014) paper, f107=130
        fuv(1)=3.5d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (0-40) nm
        fuv(2)=1.1d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (40-100) nm
        fuv(3)=9.9d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (100-140) nm
          sig_ch4(1)=4.20d-22  ! CH4 cross section in m^2, (0-40) nm
          sig_ch4(2)=4.00d-21  ! CH4 cross section in m^2, (40-100) nm
          sig_ch4(3)=1.90d-21  ! CH4 cross section in m^2, (100-140) nm
            sig_n2(1)=8.10d-22  ! N2 cross section in m^2, (0-40) nm
            sig_n2(2)=2.30d-21  ! N2 cross section in m^2, (40-100) nm
            sig_n2(3)=0.00d-22  ! N2 cross section in m^2, (100-140) nm
        endif
c
        if(jeuvfuv.eq.3) then  ! Trion solar EUV and FUV heating for Voyager 2: August 8, 1989, f107=220
        fuv(1)=5.94d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (0-40) nm
        fuv(2)=1.31d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (40-100) nm
        fuv(3)=13.34d-3/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (100-140) nm
          sig_ch4(1)=3.67d-22  ! CH4 cross section in m^2, (0-40) nm
          sig_ch4(2)=4.06d-21  ! CH4 cross section in m^2, (40-100) nm
          sig_ch4(3)=1.90d-21  ! CH4 cross section in m^2, (100-140) nm
            sig_n2(1)=7.40d-22  ! N2 cross section in m^2, (0-40) nm
            sig_n2(2)=2.23d-21  ! N2 cross section in m^2, (40-100) nm
            sig_n2(3)=0.00d-22  ! N2 cross section in m^2, (100-140) nm
        endif
c
      else 
c       
          xf107=(f107-70.0d0)/150.0d0   ! (220-70)=150
        fuv(1)=(1.65d-3+(5.94d-3-1.65d-3)*xf107)/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (0-40) nm
        fuv(2)=(0.74d-3+(1.31d-3-0.74d-3)*xf107)/xfsq   ! solar flux at 30 AU in [J m^-2 s^-1], (40-100) nm
        fuv(3)=(7.46d-3+(13.34d-3-7.46d-3)*xf107)/xfsq  ! solar flux at 30 AU in [J m^-2 s^-1], (100-140) nm
          sig_ch4(1)=(4.05d-22+4.20d-22+3.67d-22)/3.0d0  ! CH4 cross section in m^2, (0-40) nm
          sig_ch4(2)=(4.07d-21+4.00d-21+4.06d-21)/3.0d0  ! CH4 cross section in m^2, (40-100) nm
          sig_ch4(3)=(1.88d-21+1.90d-21+1.90d-21)/3.0d0  ! CH4 cross section in m^2, (100-140) nm
            sig_n2(1)=(8.27d-22+8.10d-22+7.40d-22)/3.0d0  ! N2 cross section in m^2, (0-40) nm
            sig_n2(2)=(2.26d-21+2.30d-21+2.23d-21)/3.0d0  ! N2 cross section in m^2, (40-100) nm
            sig_n2(3)=0.0d-22                             ! N2 cross section in m^2, (100-140) nm
c
      endif
c
      do n=1,nm
      fuv(n)=fuv(n)/(rnormf*rnormf)*xmu0          ! add effect on the solar flux
cc 2016.0407      sf_ch4(n)=sig_ch4(n)*fuv(n)*(0.5d0*0.5d0)    ! also add (heating efficiency)*(global average)
      sf_ch4(n)=sig_ch4(n)*fuv(n)*(0.5d0*0.5d0)    ! also add (heating efficiency)*(global average)
cc 2016.0502      sf_n2(n)=sig_n2(n)*fuv(n)*(0.25d0*0.5d0)     ! also add (heating efficiency)*(global average)
      sf_n2(n)=sig_n2(n)*fuv(n)*(0.25d0*0.5d0)     ! also add (heating efficiency)*(global average)
      enddo
c
c  to add heating rate by Lyman alpha with only one solar zenith angle
c
      xmu1=0.5d0
C      xmu1=xmu0 !AMZ actual zenith angle
      do k=1,km
      wk1(k)=0.0d0
        do n=1,nm     ! heating rate: rho*c_p*dT/dt
        fac_ch4=an(k)*sf_ch4(n)
     &    *dexp(-dmin1(sig_ch4(n)*bn(k)+sig_n2(n)*dn(k),200.0d0)/xmu1)
        fac_n2=cn(k)*sf_n2(n)
     &    *dexp(-dmin1(sig_ch4(n)*bn(k)+sig_n2(n)*dn(k),200.0d0)/xmu1)
        wk1(k)=wk1(k)+fac_ch4+fac_n2
           ht_ch4(k,n)=fac_ch4
           ht_n2(k,n)=fac_n2
        enddo
      rhoairn2=pre(k)/(rair*tem(k))
      xffcc1=cpair*rhoairn2
      wk2(k)=wk1(k)/xffcc1   ! heating rate of dT/dt in K/s
        do n=1,nm
        ht_ch4(k,n)=ht_ch4(k,n)/xffcc1   ! heating rate of dT/dt in K/s
        ht_n2(k,n)=ht_n2(k,n)/xffcc1   ! heating rate of dT/dt in K/s
        enddo
      enddo
c
        do k=1,km
        ht(k)=ht(k)+wk2(k)  ! add EUV_FUV heating by N2 and CH4 to basic heating profile
        enddo
c
cc
      RETURN
      END

c
C  To calculate the solar uv heating rate in [K s^-1] by CH4 nu3 band at Pluto.
C  AN=CH4 number density (m^-3), BN=CH4 column number density (m^-3); 
C  TEM=temperature, PRE=pressure in Pa, assuming N2 dominant atmosphere
C  IYELLE=1: thin limit; IYELLE=2: thick limit;
c
c  revised on 2014.08.21 and is based on "htch4nt" in the old code (Zhu et al. 2014) 
c    with an updated CH4 parameterization (18 August 2014)
c
c  revised on 2015.10.13 based on the research notes by Strobel: 
c    "Updated CH4 parameterization: 5 October 2015" for including 1.6 micron absorption
c
      SUBROUTINE htch4nt_r(HT,AN,BN,TEM,PRE,KM,XMU,IYELLE,GAMN4,IZHU
     & ,rnormf,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/CH4V3/sigx(jm,lm,nm_p),strx_b1(lm),strx_b2(lm)     ! [nu 3], [2800, 3200] cm^-1, 3.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMX)
     &   ,EPS2(KMX),EPS3(KMX),GAMSL1(KMX),GAMSL2(KMX)
      DIMENSION GAMN4(KM)
C      BK=1.3807D-23      ! Boltzmann constant in J K^-1
      FF0=9.1110D0       ! solar flux [pi B] at nu3 band in [W m^-2 / m^-1]
      EPS1=2.40232D-8    ! geometric factor for the solar constant in steradian
c1      XKVT3=XP3*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
c1      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
c
      do k=1,km
      do j=1,jm
      sig1(j,k)=zkg_tp_int(tem(k),pre(k),j,1)
      enddo
      enddo
c
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
      FACT2=DMIN1(SIG1(J,K)*BN(K),200.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN(K)/XMU,200.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      IF(IYELLE.EQ.2) GAMSL1(K)=FAC2/FAC1  ! one-side heating/cooling 0.5D0*FAC2/FAC1
      IF(IYELLE.EQ.3) GAMSL1(K)=1.0D0      ! Direct flux only
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
cc0504:  [nu 1]=2913, [nu 2]=1533, [nu 3]=3018, [nu 4]=1310
      rnu34=0.131875d0   ! ([nu 3] - 2*[nu 4])/[nu 3] = (3018 - 2*1310)/3018
c
cc  k_3 is taken from Boursier et al (2003). k_4 has Siddles et al (1994)
cc     temperature dependence and is normalized to the Boursier et al value
cc
c  Updated CH4 Parameterization: 3 May 2016 by Strobel
c
      do k=1,km
      airn=pre(k)/(bk*tem(k))    ! air number density in m**(-3)
      phik3=5.7d-17*( dexp(-140.0d0/tem(k)) )*airn/xa3    ! [phi]_3, xa3=25.2d0,  [3018 cm^-1]
      xfaa1=1.0d-21*(dexp((tem(k)-240.0d0)/70.0d0))*(airn+10.0d0*an(k))  ! Siddles et al. (1994)
      phik4=xfaa1/xa4                                     ! [phi]_4, xa4=2.12d0,  [1310 cm^-1]
      phi_2nu4=(1.16d-17*an(k)+xfaa1)/xa4                 ! [phi]_2[nu]4,  xa4=2.12d0,  [1306 cm^-1]
        fac11=phi_2nu4/(1.0d0+phi_2nu4)
        fac22=phik4/(gamn4(k)+phik4)
      eps2(k)=phik3/(gamsl1(K)+phik3)
      eps3(k)=rnu34+(1.0d0-rnu34)*fac11*fac22
      enddo
c
cc
      DO 80 K=1,KM
c
      sigdnu=str_t_int(tem(k),1,2)
c
      HT(K)=FF0*EPS1*EPS2(K)*AN(K)*sigdnu*GAMSL2(K)   ! solar heating rate in [W m^-3]
      IF(IZHU.EQ.1) HT(K)=HT(K)*EPS3(K)        ! vibration-vibration cascade
      HT(K)=HT(K)/(rcpr*PRE(K)/TEM(K))         ! solar heating rate in [K s^-1]
  80  CONTINUE
c
      do k=1,km
      ht(k)=ht(k)/(rnormf*rnormf)*xmu0  ! orbital effect on the solar flux
      enddo
c
      RETURN
      END



cc  calculate the enhanced component due to surface albedo
cc
      SUBROUTINE htch4nt_r_inv(HT,AN,BN,TEM,PRE,KM,XMU,IYELLE,GAMN4,IZHU
     & ,rnormf,surface_albedo,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/CH4V3/sigx(jm,lm,nm_p),strx_b1(lm),strx_b2(lm)     ! [nu 3], [2800, 3200] cm^-1, 3.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMX)
     &   ,EPS2(KMX),EPS3(KMX),GAMSL1(KMX),GAMSL2(KMX)
      DIMENSION GAMN4(KM)
      dimension an_inv(kmx),bn_inv(kmx),tem_inv(kmx),pre_inv(kmx)
     & ,gamn4_inv(kmx)
C      BK=1.3807D-23      ! Boltzmann constant in J K^-1
      FF0=9.1110D0       ! solar flux [pi B] at nu3 band in [W m^-2 / m^-1]
      EPS1=2.40232D-8    ! geometric factor for the solar constant in steradian
c1      XKVT3=XP3*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
c1      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
c
      ff0=ff0*surface_albedo
      do k=1,km
      an_inv(k)=an(km-k+1)
      bn_inv(k)=(bn(1)-bn(k))+bn(km)
      tem_inv(k)=tem(km-k+1)
      pre_inv(k)=pre(km-k+1)
      gamn4_inv(k)=gamn4(km-k+1)
      enddo
c
      do k=1,km
      do j=1,jm
      sig1(j,k)=zkg_tp_int(tem_inv(k),pre_inv(k),j,1)
      enddo
      enddo
c
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
      FACT2=DMIN1(SIG1(J,K)*BN_inv(K),200.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN_inv(K)/XMU,200.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      IF(IYELLE.EQ.2) GAMSL1(K)=FAC2/FAC1  ! one-side heating/cooling 0.5D0*FAC2/FAC1
      IF(IYELLE.EQ.3) GAMSL1(K)=1.0D0      ! Direct flux only
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
c
      trans1=0.0d0
      do j=1,jm
        opt_thick=0.0d0
        do k=2,km
        opt_thick=opt_thick+sig1(j,k)*dabs(bn_inv(k)-bn_inv(k-1))
        enddo
        trans1=trans1+ww(j)*dexp(-dmin1(1.67d0*opt_thick,200.0d0))
      enddo
c
cc0504:  [nu 1]=2913, [nu 2]=1533, [nu 3]=3018, [nu 4]=1310
      rnu34=0.131875d0   ! ([nu 3] - 2*[nu 4])/[nu 3] = (3018 - 2*1310)/3018
c
cc  k_3 is taken from Boursier et al (2003). k_4 has Siddles et al (1994)
cc     temperature dependence and is normalized to the Boursier et al value
cc
c  Updated CH4 Parameterization: 3 May 2016 by Strobel
c
      do k=1,km
      airn=pre_inv(k)/(bk*tem_inv(k))    ! air number density in m**(-3)
      phik3=5.7d-17*( dexp(-140.0d0/tem_inv(k)) )*airn/xa3    ! [phi]_3, xa3=25.2d0,  [3018 cm^-1]
      xfaa1=1.0d-21*(dexp((tem_inv(k)-240.0d0)/70.0d0))
     &    *(airn+10.0d0*an_inv(k))  ! Siddles et al. (1994)
      phik4=xfaa1/xa4                                         ! [phi]_4, xa4=2.12d0,  [1310 cm^-1]
      phi_2nu4=(1.16d-17*an_inv(k)+xfaa1)/xa4                 ! [phi]_2[nu]4,  xa4=2.12d0,  [1306 cm^-1]
        fac11=phi_2nu4/(1.0d0+phi_2nu4)
        fac22=phik4/(gamn4_inv(k)+phik4)
      eps2(k)=phik3/(gamsl1(K)+phik3)
      eps3(k)=rnu34+(1.0d0-rnu34)*fac11*fac22
      enddo
c
cc
      DO 80 K=1,KM
c
      sigdnu=str_t_int(tem_inv(k),1,2)
c
      HT(K)=trans1*FF0*EPS1*EPS2(K)*AN_inv(K)*sigdnu*GAMSL2(K)   ! solar heating rate in [W m^-3]
      IF(IZHU.EQ.1) HT(K)=HT(K)*EPS3(K)                   ! vibration-vibration cascade
      HT(K)=HT(K)/(rcpr*PRE_inv(K)/TEM_inv(K))            ! solar heating rate in [K s^-1]
  80  CONTINUE
c
      do k=1,km
      ht(k)=ht(k)/(rnormf*rnormf)*xmu0  ! orbital effect on the solar flux
      enddo
c
      do k=1,km  ! reverse the vertical coordinate for reflected flux
      ht(k)=ht(km-k+1)
      enddo
c
      RETURN
      END


c
C  To calculate the solar uv heating rate in [K s^-1] by CH4 2*nu3 band at Pluto.
C  AN=CH4 number density (m^-3), BN=CH4 column number density (m^-3); 
C  TEM=temperature, PRE=pressure in Pa, assuming N2 dominant atmosphere
C  IYELLE=1: thin limit; IYELLE=2: thick limit;
c
c  added on 2015.09.28 and is based on "htch4nt_r" for nu3 band 
c
      SUBROUTINE HTCH6t(HT,AN,BN,TEM,PRE,KM,XMU,IYELLE,GAMN4,IZHU
     & ,rnormf,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/CH4V33/sigx6(jm,lm,nm_p),strx6_b1(lm),strx6_b2(lm)  ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMX)
     &   ,EPS2(KMX),EPS3(KMX),GAMSL1(KMX),GAMSL2(KMX)
      DIMENSION GAMN4(KM)
C      BK=1.3807D-23      ! Boltzmann constant in J K^-1
      FF0=23.4454D0      ! solar flux [pi B] at 2*nu3 band in [W m^-2 / m^-1]
      EPS1=2.40232D-8    ! geometric factor for the solar constant in steradian
      XKVT3=XP3*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
c
      do k=1,km
      do j=1,jm
      sig1(j,k)=zkg_tp_int(tem(k),pre(k),j,5)
      enddo
      enddo
c
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
      FACT2=DMIN1(SIG1(J,K)*BN(K),200.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN(K)/XMU,200.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      IF(IYELLE.EQ.2) GAMSL1(K)=FAC2/FAC1  ! one-side heating/cooling 0.5D0*FAC2/FAC1
      IF(IYELLE.EQ.3) GAMSL1(K)=1.0D0      ! Direct flux only
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
c
cc0504:  [nu 1]=2913, [nu 2]=1533, [nu 3]=3018, [nu 4]=1310
      rnu16=0.13188d0   ! ([nu 1.6] - 4*[nu 4])/[nu 1.6] = (6036 - 4*1310)/6036
      rnu46=0.21703d0   ! [nu 4]/[nu 1.6] = 1310/6036
c
cc  k_3 is taken from Boursier et al (2003). k_4 has Siddles et al (1994)
cc     temperature dependence and is normalized to the Boursier et al value
c
c  Updated CH4 Parameterization: 3 May 2016 by Strobel
c
      do k=1,km
      airn=pre(k)/(bk*tem(k))    ! air number density in m**(-3)
      phik3=5.7d-17*( dexp(-140.0d0/tem(k)) )*airn/xa3    ! [phi]_2nu3
      xfaa1=1.0d-21*(dexp((tem(k)-240.0d0)/70.0d0))*(airn+10.0d0*an(k))  ! Siddles et al. (1994)
      phik4=xfaa1/xa4                                     ! [phi]_4, xa4=2.12d0,  [1310 cm^-1]
      phi_3nu4=(1.16d-17*an(k)+xfaa1)/xa4                 ! [phi]_3[nu]4
          phi_2nu4=phi_3nu4
          fac11=phi_3nu4/(1.0d0+phi_3nu4)
          fac11_4=fac11
          fac22=phik4/(gamn4(k)+phik4)
          fac33=1.0d0+2.0d0*phi_2nu4/(1.0d0+phi_2nu4)
      eps2(k)=phik3/(gamsl1(K)+phik3)
      eps3(k)=rnu16+rnu46*fac22*(fac11_4+fac33*fac11)
      enddo
c
cc
      DO 80 K=1,KM
      sigdnu=str_t_int(tem(k),5,2)
      HT(K)=FF0*EPS1*EPS2(K)*AN(K)*sigdnu*GAMSL2(K)   ! solar heating rate in [W m^-3]
      HT(K)=HT(K)/(rcpr*PRE(K)/TEM(K))        ! solar heating rate in [K s^-1]
      IF(IZHU.EQ.1) HT(K)=HT(K)*eps3(k)       ! vibration-vibration cascade
  80  CONTINUE
c
      do k=1,km
      ht(k)=ht(k)/(rnormf*rnormf)*xmu0  ! orbital effect on the solar flux
      enddo
c
      RETURN
      END


      SUBROUTINE HTCH6t_inv(HT,AN,BN,TEM,PRE,KM,XMU,IYELLE,GAMN4,IZHU
     & ,rnormf,surface_albedo,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/CH4V33/sigx6(jm,lm,nm_p),strx6_b1(lm),strx6_b2(lm)  ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMX)
     &   ,EPS2(KMX),EPS3(KMX),GAMSL1(KMX),GAMSL2(KMX)
      DIMENSION GAMN4(KM)
      dimension an_inv(kmx),bn_inv(kmx),tem_inv(kmx),pre_inv(kmx)
     & ,gamn4_inv(kmx)
C      BK=1.3807D-23      ! Boltzmann constant in J K^-1
      FF0=23.4454D0      ! solar flux [pi B] at 2*nu3 band in [W m^-2 / m^-1]
      EPS1=2.40232D-8    ! geometric factor for the solar constant in steradian
      XKVT3=XP3*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
c
c      surface_albedo=0.80d0   !  surface albedo of icy planets
      ff0=ff0*surface_albedo
      do k=1,km
      an_inv(k)=an(km-k+1)
      bn_inv(k)=(bn(1)-bn(k))+bn(km)
      tem_inv(k)=tem(km-k+1)
      pre_inv(k)=pre(km-k+1)
      gamn4_inv(k)=gamn4(km-k+1)
      enddo
c
      do k=1,km
      do j=1,jm
      sig1(j,k)=zkg_tp_int(tem_inv(k),pre_inv(k),j,5)
      enddo
      enddo
c
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
      FACT2=DMIN1(SIG1(J,K)*BN_inv(K),200.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN_inv(K)/XMU,200.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      IF(IYELLE.EQ.2) GAMSL1(K)=FAC2/FAC1  ! one-side heating/cooling 0.5D0*FAC2/FAC1
      IF(IYELLE.EQ.3) GAMSL1(K)=1.0D0      ! Direct flux only
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
c
      trans1=0.0d0
      do j=1,jm
        opt_thick=0.0d0
        do k=2,km
        opt_thick=opt_thick+sig1(j,k)*dabs(bn_inv(k)-bn_inv(k-1))
        enddo
        trans1=trans1+ww(j)*dexp(-dmin1(1.67d0*opt_thick,200.0d0))
      enddo
c
c
cc0504:  [nu 1]=2913, [nu 2]=1533, [nu 3]=3018, [nu 4]=1310
      rnu16=0.13188d0   ! ([nu 1.6] - 4*[nu 4])/[nu 1.6] = (6036 - 4*1310)/6036
      rnu46=0.21703d0   ! [nu 4]/[nu 1.6] = 1310/6036
c
cc  k_3 is taken from Boursier et al (2003). k_4 has Siddles et al (1994)
cc     temperature dependence and is normalized to the Boursier et al value
c
c  Updated CH4 Parameterization: 3 May 2016 by Strobel
c
      do k=1,km
      airn=pre_inv(k)/(bk*tem_inv(k))    ! air number density in m**(-3)
      phik3=5.7d-17*( dexp(-140.0d0/tem_inv(k)) )*airn/xa3    ! [phi]_2nu3
      xfaa1=1.0d-21*(dexp((tem_inv(k)-240.0d0)/70.0d0))
     &   *(airn+10.0d0*an_inv(k))  ! Siddles et al. (1994)
      phik4=xfaa1/xa4                                     ! [phi]_4, xa4=2.12d0,  [1310 cm^-1]
      phi_3nu4=(1.16d-17*an_inv(k)+xfaa1)/xa4                 ! [phi]_3[nu]4
          phi_2nu4=phi_3nu4
          fac11=phi_3nu4/(1.0d0+phi_3nu4)
          fac11_4=fac11
          fac22=phik4/(gamn4_inv(k)+phik4)
          fac33=1.0d0+2.0d0*phi_2nu4/(1.0d0+phi_2nu4)
      eps2(k)=phik3/(gamsl1(K)+phik3)
      eps3(k)=rnu16+rnu46*fac22*(fac11_4+fac33*fac11)
      enddo
c
cc
      DO 80 K=1,KM
      sigdnu=str_t_int(tem_inv(k),5,2)
      HT(K)=trans1*FF0*EPS1*EPS2(K)*AN_inv(K)*sigdnu*GAMSL2(K)   ! solar heating rate in [W m^-3]
      HT(K)=HT(K)/(rcpr*PRE_inv(K)/TEM_inv(K))        ! solar heating rate in [K s^-1]
      IF(IZHU.EQ.1) HT(K)=HT(K)*eps3(k)       ! vibration-vibration cascade
  80  CONTINUE
c
      do k=1,km
      ht(k)=ht(k)/(rnormf*rnormf)  ! orbital effect on the solar flux
      enddo
c
      do k=1,km  ! reverse the vertical coordinate for reflected flux
      ht(k)=ht(km-k+1)
      enddo
c
      RETURN
      END

                        
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate solar near IR heating rate by CH4 2.3 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  To calculate the globally averaged solar uv heating rate in [K s^-1] by 
C  CH4 nu3+nu4 band at Pluto/Triton.  AN=CH4 number density (m^-3), BN=CH4 column 
C  number density (m^-3); TEM=temperature, PRE=pressure in Pa, 
C  assuming N2 dominant atmosphere
c
      SUBROUTINE HTCH5t(HT,AN,BN,TEM,PRE,KM,MODE,GAMN4,rnormf,
     &                      surface_albedo,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,MM=4,nm_p=3)                !  KMX  >or= KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON/CH4V34/sigx5(jm,lm,nm_p),strx5_b1(lm),strx5_b2(lm)  ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),GAMN4(km)
      DIMENSION XX4(MM),WW4(MM),WK1(KMX)
      dimension ht_a(kmx),ht_b(kmx)
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
      CALL HTCH5nt_r(ht_a,AN,BN,TEM,PRE,KM,XMU1,GAMN4,rnormf,xmu0,bk)
      CALL HTCH5nt_r_inv(ht_b,AN,BN,TEM,PRE,KM,XMU1,GAMN4
     &  ,rnormf,surface_albedo,xmu0,bk)
      DO 30 K=1,KM
      ht(k)=ht_a(k)+ht_b(k)
      WK1(K)=WK1(K)+HT(K)*FAC2*WW4(I)*WW4(J)
  30  CONTINUE
  40  CONTINUE
      DO 60 K=1,KM
C      FACT1=(2.0D0/PI)/2.0D0 !AMZ don't pi for individual profiles
      FACT1=1.0
      HT(K)=FACT1*WK1(K)     ! globally averaged solar heating rate in [K s^-1]
  60  CONTINUE
      RETURN
      END


C  To calculate the solar uv heating rate in [K s^-1] by CH4 nu3+nu4 band at Pluto/Triton.
C  AN=CH4 number density (m^-3), BN=CH4 column number density (m^-3); 
C  TEM=temperature, PRE=pressure in Pa, assuming N2 dominant atmosphere
c
c  revised on 2014.08.21 and is based on "htch5nt" in the old code (Zhu et al. 2014) 
c    with an updated CH4 parameterization (18 August 2014)
c
      SUBROUTINE htch5nt_r(HT,AN,BN,TEM,PRE,KM,XMU,GAMN4,rnormf,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,nm_p=3)                !  KMX  >or=  KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/CH4V34/sigx5(jm,lm,nm_p),strx5_b1(lm),strx5_b2(lm)  ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMX)
     &   ,EPS2(KMX),EPS3(KMX),GAMSL1(KMX),GAMSL2(KMX)
      DIMENSION GAMN4(KM)
C      BK=1.3807D-23      ! Boltzmann constant in J K^-1
      FF0=15.5444D0    ! solar flux [pi B] at [nu3 + nu4] band in [W m^-2 / m^-1]
      EPS1=2.40232D-8  ! geometric factor for the solar constant in steradian
c1      XKVT3=XP2*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
c1      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
c
      do k=1,km
      do j=1,jm
      sig1(j,k)=zkg_tp_int(tem(k),pre(k),j,2)
      enddo
      enddo
c
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
c
cc0504:  [nu 1]=2913, [nu 2]=1533, [nu 3]=3018, [nu 4]=1310
      rnu34=0.09112d0   ! ([nu 2.3] - 3*[nu 4])/[nu 2.3] = (4324 - 3*1310)/4324
      rnu41=0.30296d0   ! [nu 4]/[nu 2.3] = 1310/4324
c
cc  k_3 is taken from Boursier et al (2003). k_4 has Siddles et al (1994)
cc     temperature dependence and is normalized to the Boursier et al value
cc
c  Updated CH4 Parameterization: 3 May 2016 by Strobel
c
      do k=1,km
      airn=pre(k)/(bk*tem(k))    ! air number density in m**(-3)
      phik3=5.7d-17*( dexp(-140.0d0/tem(k)) )*airn/xa3    ! [phi]_3, xa3=25.2d0,  [3018 cm^-1]
      xfaa1=1.0d-21*(dexp((tem(k)-240.0d0)/70.0d0))*(airn+10.0d0*an(k))  ! Siddles et al. (1994)
      phik4=xfaa1/xa4                                     ! [phi]_4, xa4=2.12d0,  [1310 cm^-1]
      phi_3nu4=(1.16d-17*an(k)+xfaa1)/xa4                 ! [phi]_3[nu]4
        xfaa1=1.0d-21*(dexp((tem(k)-240.d0)/70.d0))*(airn+10.d0*an(k))  ! Siddles et al. (1994)
        phi_2nu4=(1.16d-17*an(k)+xfaa1)/xa4               ! [phi]_2[nu]4
          fac11=phi_3nu4/(1.0d0+phi_3nu4)
          fac22=phik4/(gamn4(k)+phik4)
          fac33=1.0d0+2.0d0*phi_2nu4/(1.0d0+phi_2nu4)
      eps2(k)=phik3/(gamsl1(K)+phik3)
      eps3(k)=rnu34+rnu41*fac11*fac22*fac33
      enddo
c
cc
      DO 80 K=1,KM
c
      sigdnu=str_t_int(tem(k),2,2)
c
      HT(K)=FF0*EPS1*EPS2(K)*AN(K)*sigdnu*GAMSL2(K)   ! solar heating rate in [W m^-3]
      HT(K)=HT(K)*EPS3(K)                             ! vibration-vibration cascade
      HT(K)=HT(K)/(rcpr*PRE(K)/TEM(K))                ! solar heating rate in [K s^-1]
  80  CONTINUE
c
      do k=1,km
      ht(k)=ht(k)/(rnormf*rnormf)*xmu0  ! orbital effect on the solar flux
      enddo
c
      RETURN
      END


      SUBROUTINE htch5nt_r_inv(HT,AN,BN,TEM,PRE,KM,XMU,GAMN4,rnormf,
     &            surface_albedo,xmu0,bk)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,KMX=175,nm_p=3)                !  KMX  >or=  KM
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/CH4V34/sigx5(jm,lm,nm_p),strx5_b1(lm),strx5_b2(lm)  ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION HT(KM),AN(KM),BN(KM),TEM(KM),PRE(KM),SIG1(JM,KMX)
     &   ,EPS2(KMX),EPS3(KMX),GAMSL1(KMX),GAMSL2(KMX)
      DIMENSION GAMN4(KM)
      dimension an_inv(kmx),bn_inv(kmx),tem_inv(kmx),pre_inv(kmx)
     & ,gamn4_inv(kmx)
C      BK=1.3807D-23      ! Boltzmann constant in J K^-1
      FF0=15.5444D0    ! solar flux [pi B] at [nu3 + nu4] band in [W m^-2 / m^-1]
      EPS1=2.40232D-8  ! geometric factor for the solar constant in steradian
c1      XKVT3=XP2*ZKVT     ! collisional rate coefficient in [m**3 s**(-1)]
c1      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
c
c      surface_albedo=0.80d0   !  surface albedo of icy planets
      ff0=ff0*surface_albedo
      do k=1,km
      an_inv(k)=an(km-k+1)
      bn_inv(k)=(bn(1)-bn(k))+bn(km)
      tem_inv(k)=tem(km-k+1)
      pre_inv(k)=pre(km-k+1)
      gamn4_inv(k)=gamn4(km-k+1)
      enddo
c
      do k=1,km
      do j=1,jm
      sig1(j,k)=zkg_tp_int(tem_inv(k),pre_inv(k),j,2)
      enddo
      enddo
c
      DO 40 K=1,KM
      FAC1=0.0D0
      FAC2=0.0D0
      FAC3=0.0D0
      DO 35 J=1,JM
      FACT2=DMIN1(SIG1(J,K)*BN_inv(K),300.0D0)
      FACT3=DMIN1(SIG1(J,K)*BN_inv(K)/XMU,300.0D0)
      FAC1=FAC1+WW(J)*SIG1(J,K)
      FAC2=FAC2+WW(J)*SIG1(J,K)*ENZD2(FACT2)
      FAC3=FAC3+WW(J)*SIG1(J,K)*DEXP(-FACT3)
  35  CONTINUE
      GAMSL1(K)=FAC2/FAC1    ! one-side heating/cooling 0.5D0*FAC2/FAC1
      GAMSL2(K)=FAC3/FAC1
  40  CONTINUE
c
      trans1=0.0d0
      do j=1,jm
        opt_thick=0.0d0
        do k=2,km
        opt_thick=opt_thick+sig1(j,k)*dabs(bn_inv(k)-bn_inv(k-1))
        enddo
        trans1=trans1+ww(j)*dexp(-dmin1(1.67d0*opt_thick,200.0d0))
      enddo
c
c
cc0504:  [nu 1]=2913, [nu 2]=1533, [nu 3]=3018, [nu 4]=1310
      rnu34=0.09112d0   ! ([nu 2.3] - 3*[nu 4])/[nu 2.3] = (4324 - 3*1310)/4324
      rnu41=0.30296d0   ! [nu 4]/[nu 2.3] = 1310/4324
c
cc  k_3 is taken from Boursier et al (2003). k_4 has Siddles et al (1994)
cc     temperature dependence and is normalized to the Boursier et al value
cc
c  Updated CH4 Parameterization: 3 May 2016 by Strobel
c
      do k=1,km
      airn=pre_inv(k)/(bk*tem_inv(k))    ! air number density in m**(-3)
      phik3=5.7d-17*( dexp(-140.0d0/tem_inv(k)) )*airn/xa3    ! [phi]_3, xa3=25.2d0,  [3018 cm^-1]
      xfaa1=1.0d-21*(dexp((tem_inv(k)-240.0d0)/70.0d0))
     &    *(airn+10.0d0*an_inv(k))  ! Siddles et al. (1994)
      phik4=xfaa1/xa4                                     ! [phi]_4, xa4=2.12d0,  [1310 cm^-1]
      phi_3nu4=(1.16d-17*an_inv(k)+xfaa1)/xa4                 ! [phi]_3[nu]4
        xfaa1=1.0d-21*(dexp((tem_inv(k)-240.d0)/70.d0))
     &   *(airn+10.d0*an_inv(k))  ! Siddles et al. (1994)
        phi_2nu4=(1.16d-17*an_inv(k)+xfaa1)/xa4               ! [phi]_2[nu]4
          fac11=phi_3nu4/(1.0d0+phi_3nu4)
          fac22=phik4/(gamn4_inv(k)+phik4)
          fac33=1.0d0+2.0d0*phi_2nu4/(1.0d0+phi_2nu4)
      eps2(k)=phik3/(gamsl1(K)+phik3)
      eps3(k)=rnu34+rnu41*fac11*fac22*fac33
      enddo
c
cc
      DO 80 K=1,KM
c
      sigdnu=str_t_int(tem_inv(k),2,2)
c
      HT(K)=trans1*FF0*EPS1*EPS2(K)*AN_inv(K)*sigdnu*GAMSL2(K)   ! solar heating rate in [W m^-3]
      HT(K)=HT(K)*EPS3(K)                             ! vibration-vibration cascade
      HT(K)=HT(K)/(rcpr*PRE_inv(K)/TEM_inv(K))                ! solar heating rate in [K s^-1]
  80  CONTINUE
c
      do k=1,km
      ht(k)=ht(k)/(rnormf*rnormf)*xmu0  ! orbital effect on the solar flux
      enddo
c
      do k=1,km  ! reverse the vertical coordinate for reflected flux
      ht(k)=ht(km-k+1)
      enddo
c
      RETURN
      END



c  Read in the k-distribution function derived from the Hitran-2012 dataset
c
      SUBROUTINE INPUT4
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,nm_p=3)
      COMMON/CH4V3/sigx(jm,lm,nm_p),strx_b1(lm),strx_b2(lm)     ! [nu 3], [2800, 3200] cm^-1, 3.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      OPEN(11,FILE='pluto_triton_dat/ggww_30.in',status='OLD')
      DO 10 J=1,JM
      READ(11,*) GG(J),WW(J)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_triton_dat/lbl_nu3_kgtp.dat',status='OLD')
        read(12,'(3i6)') jmgx1,lmgx1,nmgx1
        read(12,'(10f9.6)') (ggjgx1,j=1,jm)         ! g-function
        read(12,'(10f8.2)') (tref(l),l=1,lm)        ! temperaure in K
        read(12,'(10e12.5)') (strx_b1(l),l=1,lm)    ! band strength in m^-1 / (kg m^-2) = m kg^-1
        read(12,'(10e12.5)') (strx_b2(l),l=1,lm)    ! band strength in m
        read(12,'(10e11.4)') (pref(n),n=1,nm_p)     ! pressure in Pascal (1 Pa = 10 micro bar)
      do n=1,nm_p
      pref_log(n)=dlog(pref(n))
      read(12,'(1p10e12.5)') ((sigx(j,l,n),j=1,jm),l=1,lm)  ! k-distribution or k-coefficient in cm^2 / g 
      enddo
      CLOSE(UNIT=12)
      fact1=2.6768d-27  ! convertion factor from cm^2/g to m^2 = 1.673e-24*16*1.0e-4
        do n=1,nm_p
        do l=1,lm
        do j=1,jm
          sigx(j,l,n)=sigx(j,l,n)*fact1    ! convert absorption coefficient to cross section in m**2
        enddo
        enddo
        enddo
      return
      end

c  Read in the k-distribution function derived from the Hitran-2012 dataset
c
      SUBROUTINE INPUT5
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,nm_p=3)
      COMMON/CH4V34/sigx5(jm,lm,nm_p),strx5_b1(lm),strx5_b2(lm)  ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      OPEN(11,FILE='pluto_triton_dat/ggww_30.in',status='old')
      DO 10 J=1,JM
      READ(11,*) gg(j),ww(j)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_triton_dat/lbl_nu34_kgtp.dat',status='old')
        read(12,'(3i6)') jmgx1,lmgx1,nmgx1
        read(12,'(10f9.6)') (ggjgx1,j=1,jm)         ! g-function
        read(12,'(10f8.2)') (tref(l),l=1,lm)        ! temperaure in K
        read(12,'(10e12.5)') (strx5_b1(l),l=1,lm)   ! band strength in m^-1 / (kg m^-2) = m kg^-1
        read(12,'(10e12.5)') (strx5_b2(l),l=1,lm)   ! band strength in m
        read(12,'(10e11.4)') (pref(n),n=1,nm_p)     ! pressure in Pascal (1 Pa = 10 micro bar)
      do n=1,nm_p
      pref_log(n)=dlog(pref(n))
      read(12,'(1p10e12.5)') ((sigx5(j,l,n),j=1,jm),l=1,lm)  ! k-distribution or k-coefficient in cm^2 / g 
      enddo
      CLOSE(UNIT=12)
      fact1=2.6768d-27  ! convertion factor from cm^2/g to m^2 = 1.673e-24*16*1.0e-4
        do n=1,nm_p
        do l=1,lm
        do j=1,jm
          sigx5(j,l,n)=sigx5(j,l,n)*fact1    ! convert absorption coefficient to cross section in m**2
        enddo
        enddo
        enddo
      return
      end            

c  Read in the k-distribution function derived from the Hitran-2012 dataset
c
      SUBROUTINE INPUT33
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,nm_p=3)
      COMMON/CH4V33/sigx6(jm,lm,nm_p),strx6_b1(lm),strx6_b2(lm)  ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      OPEN(11,FILE='pluto_triton_dat/ggww_30.in',status='old')
      DO 10 J=1,JM
      READ(11,*) gg(j),ww(j)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_triton_dat/lbl_2nu3_kgtp.dat',status='old')
        read(12,'(3i6)') jmgx1,lmgx1,nmgx1
        read(12,'(10f9.6)') (ggjgx1,j=1,jm)         ! g-function
        read(12,'(10f8.2)') (tref(l),l=1,lm)        ! temperaure in K
        read(12,'(10e12.5)') (strx6_b1(l),l=1,lm)   ! band strength in m^-1 / (kg m^-2) = m kg^-1
        read(12,'(10e12.5)') (strx6_b2(l),l=1,lm)   ! band strength in m
        read(12,'(10e11.4)') (pref(n),n=1,nm_p)     ! pressure in Pascal (1 Pa = 10 micro bar)
      do n=1,nm_p
      pref_log(n)=dlog(pref(n))
      read(12,'(1p10e12.5)') ((sigx6(j,l,n),j=1,jm),l=1,lm)  ! k-distribution or k-coefficient in cm^2 / g 
      enddo
      CLOSE(UNIT=12)
      fact1=2.6768d-27  ! convertion factor from cm^2/g to m^2 = 1.673e-24*16*1.0e-4
        do n=1,nm_p
        do l=1,lm
        do j=1,jm
          sigx6(j,l,n)=sigx6(j,l,n)*fact1    ! convert absorption coefficient to cross section in m**2
        enddo
        enddo
        enddo
      return
      end            
            
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate non-LTE cooling rate by CH4 7.6 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C  To calculate the cooling rate QJ [K s^-1] for CH4 7.6 micron band (K/s)
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
c
      SUBROUTINE QTCH4Z(TEM,PRE,RHO,QT,KM,MODE,GAMN4,IYELLE,grv)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(MCH4=2,KMX=175)                     ! KMX = KM  *********************
      COMMON /FRENG4/ VMCH4(MCH4),STRCH4(MCH4),ENECH4(MCH4)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QT(KM),grv(km)
      DIMENSION WK1(KMX),WK2(KMX),GAMN4(KM)
      dimension strj_m(kmx,mch4),strtot_k(kmx),strj_k(kmx)
      IF(MODE.EQ.1) CALL INPUT6    ! read in CH4 spectral data the 1st time calling QTCH4Z
      VMCH4(1)=1157.5D2            ! band center in m^-1
      VMCH4(2)=1292.5D2            ! band center in m^-1
      DO 20 K=1,KM
      QT(K)=0.0D0
      GAMN4(K)=0.0D0
  20  CONTINUE
        do k=1,km
        strj_m(k,1)=str_t_int(tem(k),3,1)
        strj_m(k,2)=str_t_int(tem(k),4,1)
        strtot_k(k)=strj_m(k,1)+strj_m(k,2)
        enddo
      DO 50 M=1,MCH4
      VMJ=VMCH4(M)
        do k=1,km
        strj_k(k)=strj_m(k,m)
        enddo
        CALL QJCH4X(TEM,PRE,RHO,WK1,KM,VMJ,strj_k,M,WK2,IYELLE,grv)
      DO 40 K=1,KM
      QT(K)=QT(K)-WK1(K)      ! total CH4 cooling rate in [K s^-1]
      GAMN4(K)=GAMN4(K)+strj_k(k)*WK2(K)/strtot_k(k)
  40  CONTINUE
  50  CONTINUE  
      RETURN
      END


C  To calculate the cooling rate QJ [K s^-1] for CH4 7.6 micron band
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
C  VMJ [m^-1] is the frequency, (STRJ) is the band strength
c
      SUBROUTINE QJCH4X(TEM,PRE,RHO,QJ,KM,VMJ,strj_k,ICON
     &    ,GAMN4,IYELLE,grv)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175,JM=30,LM=7,KMP=KMX+1,nm_p=3)          ! KMX = KM  *********************
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/CH4V4/sig4a(jm,lm,nm_p),str4a_b1(lm),str4a_b2(lm)   ! [nu 4a], [1090,1225] cm^-1, 7.6-a micron
     &            ,sig4b(jm,lm,nm_p),str4b_b1(lm),str4b_b2(lm)   ! [nu 4b], [1225,1360] cm^-1, 7.6-b micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QJ(KM),GAMN4(KM),grv(km)
      DIMENSION BLACX(KMX),H00(KMX),XKZ(KMX,JM),GAM1(KMX,KMX)
     & ,CURT(KMX,KMX),EE(KMX,KMX),CURTN(KMX,KMX),CINVT(KMX,KMP)
      dimension strj_k(km)
      PI=3.14159265D0
      DO 10 K=1,KM
      BLACX(K)=BLACM(VMJ,TEM(K))
  10  CONTINUE
      DO 20 J=1,JM
      DO 20 K=1,KM
      XKZ(K,J)=zkg_tp_int(tem(k),pre(k),j,ICON+2)
  20  CONTINUE
      DO 30 K=1,KM
      DO 30 J=1,KM
      GAM1(J,K)=GAMM12(J,K,PRE,RHO,KM,XKZ,WW,JM,grv)
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
      H00(L1)=2.0D0*PI*RHO(L1)*strj_k(l1)/cpair
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
cc0504:  [nu 1]=2913, [nu 2]=1533, [nu 3]=3018, [nu 4]=1310
      DO 143 L2=1,KM
      DO 143 L1=1,KM
 143  EE(L1,L2)=0.0D0
      BK=1.3807D-23                  ! Boltzmann constant in J K^-1
      XKVT4=XP4*ZKVT     ! collisional rate coefficient at T=80 K in [m**3 s**(-1)]
cc
c  Updated CH4 Parameterization: 3 May 2016 by Strobel
c
      do k=1,km
      airn=pre(k)/(bk*tem(k))        ! air number density in m**(-3)
      an_k4=(rho(k)/(16.0d0/wtmol))*airn    ! CH4 number density in m**(-3), 16.0/28.0=0.57143
      xfaa1=1.0d-21*(dexp((tem(k)-240.0d0)/70.0d0))*(airn+10.0d0*an_k4)  ! Siddles et al. (1994)
      phik4=xfaa1/xa4                       ! [phi]_4, xa4=2.12d0,  [1310 cm^-1]
      ee(k,k)=1.0d0/(2.0d0*h00(k)*phik4)
      enddo
c
cc
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


cc To calculate the escape function according to Zhu (1992, Eq. (13)). 

c

      FUNCTION GAMM12(I,J,PUS,RUS,KM,XKZ,WW,JM,grv)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

      DIMENSION PUS(KM),RUS(KM),XKZ(KM,JM),WW(JM),grv(km)

c      GRAV=0.637D0               ! gravitational acceleration in m s^-2

      IF(I.EQ.J) GO TO 100

      KA=MIN0(I,J)

      KB=MAX0(I,J)-1

      SUM1=0.0D0

      SUM2=0.0D0

      DO 50 M=1,JM

      DELU=0.0D0

      DO 10 K=KA,KB

      FAC1=(RUS(K)+RUS(K+1))*(PUS(K+1)-PUS(K))/(0.5d0*(grv(k)+grv(k+1)))

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

      FAC1=FAC1*RUS(I)/grv(i)

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

c  Read in the k-distribution function derived from the Hitran-2012 dataset
c
      SUBROUTINE INPUT6
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,nm_p=3)
      COMMON/CH4V4/sig4a(jm,lm,nm_p),str4a_b1(lm),str4a_b2(lm)   ! [nu 4a], [1090,1225] cm^-1, 7.6-a micron
     &            ,sig4b(jm,lm,nm_p),str4b_b1(lm),str4b_b2(lm)   ! [nu 4b], [1225,1360] cm^-1, 7.6-b micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      OPEN(11,FILE='pluto_triton_dat/ggww_30.in',status='old')
      DO 10 J=1,JM
      READ(11,*) gg(j),ww(j)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_triton_dat/lbl_nu4a_kgtp.dat',status='old')
      OPEN(13,FILE='pluto_triton_dat/lbl_nu4b_kgtp.dat',status='old')
        read(12,'(3i6)') jmgx1,lmgx1,nmgx1
        read(12,'(10f9.6)') (ggjgx1,j=1,jm)         ! g-function
        read(12,'(10f8.2)') (tref(l),l=1,lm)        ! temperaure in K
        read(12,'(10e12.5)') (str4a_b1(l),l=1,lm)   ! band strength in m^-1 / (kg m^-2) = m kg^-1
        read(12,'(10e12.5)') (str4a_b2(l),l=1,lm)   ! band strength in m
        read(12,'(10e11.4)') (pref(n),n=1,nm_p)     ! pressure in Pascal (1 Pa = 10 micro bar)
      do n=1,nm_p
      pref_log(n)=dlog(pref(n))
      read(12,'(1p10e12.5)') ((sig4a(j,l,n),j=1,jm),l=1,lm)  ! k-distribution or k-coefficient in cm^2 / g 
      enddo
        read(13,'(3i6)') jmgx1,lmgx1,nmgx1
        read(13,'(10f9.6)') (ggjgx1,j=1,jm)         ! g-function
        read(13,'(10f8.2)') (tref(l),l=1,lm)        ! temperaure in K
        read(13,'(10e12.5)') (str4b_b1(l),l=1,lm)   ! band strength in m^-1 / (kg m^-2) = m kg^-1
        read(13,'(10e12.5)') (str4b_b2(l),l=1,lm)   ! band strength in m
        read(13,'(10e11.4)') (pref(n),n=1,nm_p)     ! pressure in Pascal (1 Pa = 10 micro bar)
      do n=1,nm_p
      pref_log(n)=dlog(pref(n))
      read(13,'(1p10e12.5)') ((sig4b(j,l,n),j=1,jm),l=1,lm)  ! k-distribution or k-coefficient in cm^2 / g 
      enddo
       CLOSE(UNIT=12)
       CLOSE(UNIT=13)
cc      fact1=2.6768d-27  ! convertion factor from cm^2/g to m^2 = 1.673e-24*16*1.0e-4
      fact1=0.1d0      ! conversion factor from cm^2/g to m^2/kg
        do n=1,nm_p
        do l=1,lm
        do j=1,jm
          sig4a(j,l,n)=sig4a(j,l,n)*fact1    ! convert absorption coefficient to m^2/kg
          sig4b(j,l,n)=sig4b(j,l,n)*fact1    ! convert absorption coefficient to m^2/kg
        enddo
        enddo
        enddo
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  End of calculating non-LTE cooling rate (K/day) by CH4 7.6 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate non-LTE cooling rate by C2H2 13.7 micron band (K/s)
c     adapted from QTCH4Z for cooling rate by CH4 7.6 micron band
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  To calculate the cooling rate QJ [K s^-1] for C2H2 13.7 micron (K/s)
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
c
      SUBROUTINE QTC2H2Z(TEM,PRE,RHO,QT,KM,MODE,IYELLE,grv)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(mc2h2=1,KMX=175)                     ! KMX = KM  *********************
      COMMON /FRENG5/ vmc2h2(mc2h2),STRc2h2(mc2h2),ENEc2h2(mc2h2)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QT(KM),grv(km)
      DIMENSION WK1(KMX)
      dimension strj_m(kmx,mc2h2),strtot_k(kmx),strj_k(kmx)
      IF(MODE.EQ.1) CALL INPUT6c   ! read in c2h2 spectral data the 1st time calling QTC2H2Z
      vmc2h2(1)=730.0D2            ! band center in m^-1
      DO 20 K=1,KM
      QT(K)=0.0D0
  20  CONTINUE
        do k=1,km
        strj_m(k,1)=str_t_int(tem(k),6,1)
        strtot_k(k)=strj_m(k,1)  ! +strj_m(k,2)
        enddo
      DO 50 M=1,mc2h2
      VMJ=vmc2h2(M)
        do k=1,km
        strj_k(k)=strj_m(k,m)
        enddo
        CALL QJc2h2X(TEM,PRE,RHO,WK1,KM,VMJ,strj_k,IYELLE,grv)
      DO 40 K=1,KM
      QT(K)=QT(K)-WK1(K)      ! total c2h2 cooling rate in [K s^-1]
  40  CONTINUE
  50  CONTINUE  
      RETURN
      END




C  To calculate the cooling rate QJ [K s^-1] for c2h2 13.7 micron band
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
C  VMJ [m^-1] is the frequency, (STRJ) is the band strength
c
      SUBROUTINE QJc2h2X(TEM,PRE,RHO,QJ,KM,VMJ,strj_k,IYELLE,grv)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175,JM=30,LM=7,KMP=KMX+1,nm_p=3)          ! KMX = KM  *********************
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      COMMON /PANDA/ ZKVT,XA2,XA3,XA4,XP2,XP3,XP4
      COMMON/c2h2v5/sig5(jm,lm,nm_p),str5_b1(lm),str5_b2(lm)     ! [nu 5],  [640,820] cm^-1, 13.7 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QJ(KM),grv(km)
      DIMENSION BLACX(KMX),H00(KMX),XKZ(KMX,JM),GAM1(KMX,KMX)
     & ,CURT(KMX,KMX),EE(KMX,KMX),CURTN(KMX,KMX),CINVT(KMX,KMP)
     & ,GAMN4(KMx)
      dimension strj_k(km)
      PI=3.14159265D0
      DO 10 K=1,KM
      BLACX(K)=BLACM(VMJ,TEM(K))      ! Planck function in [W/(m*m) * (m^-1)^-1 * Sr^-1]
  10  CONTINUE
      DO 20 J=1,JM
      DO 20 K=1,KM
      XKZ(K,J)=zkg_tp_int(tem(k),pre(k),j,6)
  20  CONTINUE
      DO 30 K=1,KM
      DO 30 J=1,KM
      GAM1(J,K)=GAMM12(J,K,PRE,RHO,KM,XKZ,WW,JM,grv)
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
      H00(L1)=2.0D0*PI*RHO(L1)*strj_k(l1)/cpair
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
      BK=1.3807D-23      ! Boltzmann constant in J K^-1
      xa5_c2h2=5.09d0    ! Einstein A coefficient in s^-1
      DO 145 K=1,KM
      XKVT4=2.43d-20*dexp( (tem(k)-240.0d0)/105.0d0 )  ! collisional deactivation rate coefficient in [m**3 s**(-1)]
      AIRN=PRE(K)/(BK*TEM(K))      ! air number density in m**(-3)
      phi21_c2h2=XKVT4*AIRN/xa5_c2h2
      EE(K,K)=1.0D0/(2.0D0*H00(K)*phi21_c2h2)
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


c  Read in the k-distribution function derived from the Hitran-2012 dataset
c
      SUBROUTINE INPUT6c
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(JM=30,LM=7,nm_p=3)
      COMMON/c2h2v5/sig5(jm,lm,nm_p),str5_b1(lm),str5_b2(lm)     ! [nu 5],  [640,820] cm^-1, 13.7 micron
      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)
      OPEN(11,FILE='pluto_triton_dat/ggww_30.in',status='old')
      DO 10 J=1,JM
      READ(11,*) gg(j),ww(j)
  10  CONTINUE
      CLOSE(UNIT=11)
      OPEN(12,FILE='pluto_triton_dat/c2h2_lbl_nu5_kgtp.dat'
     & ,status='old')
        read(12,'(3i6)') jmgx1,lmgx1,nmgx1
        read(12,'(10f9.6)') (ggjgx1,j=1,jm)        ! g-function
        read(12,'(10f8.2)') (tref(l),l=1,lm)       ! temperaure in K
        read(12,'(10e12.5)') (str5_b1(l),l=1,lm)   ! band strength in m^-1 / (kg m^-2) = m kg^-1
        read(12,'(10e12.5)') (str5_b2(l),l=1,lm)   ! band strength in m
        read(12,'(10e11.4)') (pref(n),n=1,nm_p)    ! pressure in Pascal (1 Pa = 10 micro bar)
      do n=1,nm_p
      pref_log(n)=dlog(pref(n))
      read(12,'(1p10e12.5)') ((sig5(j,l,n),j=1,jm),l=1,lm)  ! k-distribution or k-coefficient in cm^2 / g 
      enddo
       CLOSE(UNIT=12)
cc      fact1=4.3498d-27  ! conversion factor for C2H2 from cm^2/g to m^2 = 1.673e-24*26*1.0e-4
      fact1=0.1d0      ! conversion factor from cm^2/g to m^2/kg
        do n=1,nm_p
        do l=1,lm
        do j=1,jm
          sig5(j,l,n)=sig5(j,l,n)*fact1    ! convert absorption coefficient to m^2/kg
        enddo
        enddo
        enddo
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  End of calculating non-LTE cooling rate (K/day) by C2H2 13.7 micron band (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  To intepolate the band strength to t grid. jcon=(1,2,3,4,5,6) => (strx_b1,strx5_b1,str4a_b1,str4b_b1,strx6_b1,sig5)

c  icon = (1,2) = (m/kg, m)

c

      function str_t_int(t,jcon,icon)

      implicit real*8(a-h,o-z), INTEGER(I-N)

      parameter(JM=30,LM=7,nm_p=3)

      COMMON/CH4V3/sigx(jm,lm,nm_p),strx_b1(lm),strx_b2(lm)      ! [nu 3], [2800, 3200] cm^-1, 3.3 micron

      COMMON/CH4V33/sigx6(jm,lm,nm_p),strx6_b1(lm),strx6_b2(lm)  ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron

      COMMON/CH4V34/sigx5(jm,lm,nm_p),strx5_b1(lm),strx5_b2(lm)  ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron

      COMMON/CH4V4/sig4a(jm,lm,nm_p),str4a_b1(lm),str4a_b2(lm)   ! [nu 4a], [1090,1225] cm^-1, 7.6-a micron

     &            ,sig4b(jm,lm,nm_p),str4b_b1(lm),str4b_b2(lm)   ! [nu 4b], [1225,1360] cm^-1, 7.6-b micron

      COMMON/c2h2v5/sig5(jm,lm,nm_p),str5_b1(lm),str5_b2(lm)     ! [nu 5],  [640,820] cm^-1, 13.7 micron



      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)

      parameter(lmx=7)  ! lmx >= lm

      dimension wr1(lmx)

c

      if(icon.eq.1) then  !! band strength in m^-1 / (kg m^-2) = m kg^-1

      if(jcon.eq.1) then       ! [nu 3], [2800, 3200] cm^-1, 3.3 micron

        do l=1,lm

        wr1(l)=strx_b1(l)

        enddo

      endif

c

      if(jcon.eq.2) then       ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron

        do l=1,lm

        wr1(l)=strx5_b1(l)

        enddo

      endif

c

      if(jcon.eq.3) then       ! [nu 4a], [1090,1225] cm^-1, 7.6-a micron

        do l=1,lm

        wr1(l)=str4a_b1(l)

        enddo

      endif

c

      if(jcon.eq.4) then       ! [nu 4b], [1225,1360] cm^-1, 7.6-b micron

        do l=1,lm

        wr1(l)=str4b_b1(l)

        enddo

      endif

c

      if(jcon.eq.5) then       ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron

        do l=1,lm

        wr1(l)=strx6_b1(l)

        enddo

      endif

c

      if(jcon.eq.6) then       ! [nu 5],  [640,820] cm^-1, 13.7 micron

        do l=1,lm

        wr1(l)=str5_b1(l)

        enddo

      endif

c

      endif

c

      if(icon.eq.2) then  !! band strength in m

      if(jcon.eq.1) then       ! [nu 3], [2800, 3200] cm^-1, 3.3 micron

        do l=1,lm

        wr1(l)=strx_b2(l)

        enddo

      endif

c

      if(jcon.eq.2) then       ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron

        do l=1,lm

        wr1(l)=strx5_b2(l)

        enddo

      endif

c

      if(jcon.eq.3) then       ! [nu 4a], [1090,1225] cm^-1, 7.6-a micron

        do l=1,lm

        wr1(l)=str4a_b2(l)

        enddo

      endif

c

      if(jcon.eq.4) then       ! [nu 4b], [1225,1360] cm^-1, 7.6-b micron

        do l=1,lm

        wr1(l)=str4b_b2(l)

        enddo

      endif

c

      if(jcon.eq.5) then       ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron

        do l=1,lm

        wr1(l)=strx6_b2(l)

        enddo

      endif

c

      if(jcon.eq.6) then       ! [nu 5],  [640,820] cm^-1, 13.7 micron

        do l=1,lm

        wr1(l)=str5_b2(l)

        enddo

      endif

c

      endif

c

      if(t.le.tref(1)) then

      str_t_int=wr1(1)

      return

      endif

c

      if(t.ge.tref(lm)) then

      str_t_int=wr1(lm)

      return

      endif

c

      do l=2,lm

       if(t.ge.tref(l-1) .and. t.lt.tref(l)) then  ! linear interpolation

       str_t_int=wr1(l-1)

     &  +(t-tref(l-1))*(wr1(l)-wr1(l-1))/(tref(l)-tref(l-1))

       return

       endif

      enddo

c

      return

      end




C  To intepolate the k-coefficients to (t,p) grid. jcon=(1,2,3,4,5,6) => (strx_b1,strx5_b1,str4a_b1,str4b_b1,strx6_b1,sig5)

c

      function zkg_tp_int(t,pre,j,jcon)

      implicit real*8(a-h,o-z), INTEGER(I-N)

      parameter(JM=30,LM=7,nm_p=3)

      COMMON/CH4V3/sigx(jm,lm,nm_p),strx_b1(lm),strx_b2(lm)      ! [nu 3], [2800, 3200] cm^-1, 3.3 micron

      COMMON/CH4V33/sigx6(jm,lm,nm_p),strx6_b1(lm),strx6_b2(lm)  ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron

      COMMON/CH4V34/sigx5(jm,lm,nm_p),strx5_b1(lm),strx5_b2(lm)  ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron

      COMMON/CH4V4/sig4a(jm,lm,nm_p),str4a_b1(lm),str4a_b2(lm)   ! [nu 4a], [1090,1225] cm^-1, 7.6-a micron

     &            ,sig4b(jm,lm,nm_p),str4b_b1(lm),str4b_b2(lm)   ! [nu 4b], [1225,1360] cm^-1, 7.6-b micron

      COMMON/c2h2v5/sig5(jm,lm,nm_p),str5_b1(lm),str5_b2(lm)     ! [nu 5],  [640,820] cm^-1, 13.7 micron

      COMMON/CH4ALL/gg(jm),ww(jm),tref(lm),pref(nm_p),pref_log(nm_p)

      parameter(lmx=7,nmx3=3)  ! lmx >= lm, nmx3 >= nm_p

      dimension wk1(lmx,nmx3),wr1(nmx3)

c

      if(jcon.eq.1) then       ! [nu 3], [2800, 3200] cm^-1, 3.3 micron

        do n=1,nm_p

        do l=1,lm

        wk1(l,n)=dlog(sigx(j,l,n))

        enddo

        enddo

      endif

c

      if(jcon.eq.2) then       ! [nu 3 + nu 4], [4000,4600] cm^-1, 2.3 micron

        do n=1,nm_p

        do l=1,lm

        wk1(l,n)=dlog(sigx5(j,l,n))

        enddo

        enddo

      endif

c

      if(jcon.eq.3) then       ! [nu 4a], [1090,1225] cm^-1, 7.6-a micron

        do n=1,nm_p

        do l=1,lm

        wk1(l,n)=dlog(sig4a(j,l,n))

        enddo

        enddo

      endif

c

      if(jcon.eq.4) then       ! [nu 4b], [1225,1360] cm^-1, 7.6-b micron

        do n=1,nm_p

        do l=1,lm

        wk1(l,n)=dlog(sig4b(j,l,n))

        enddo

        enddo

      endif

c

      if(jcon.eq.5) then       ! [2 * nu 3] band from 5800 to 6200 cm^-1, 1.67 micron

        do n=1,nm_p

        do l=1,lm

        wk1(l,n)=dlog(sigx6(j,l,n))

        enddo

        enddo

      endif

c

      if(jcon.eq.6) then       ! [nu 5],  [640,820] cm^-1, 13.7 micron

        do n=1,nm_p

        do l=1,lm

        wk1(l,n)=dlog(sig5(j,l,n))

        enddo

        enddo

      endif

c

      if(t.le.tref(1)) then

        do n=1,nm_p

        wr1(n)=wk1(1,n)

        enddo

          goto 10

      endif

c

      if(t.ge.tref(lm)) then

        do n=1,nm_p

        wr1(n)=wk1(lm,n)

        enddo

          goto 10

      endif

c

      do l=2,lm

        if(t.ge.tref(l-1) .and. t.lt.tref(l)) then  ! linear interpolation

          do n=1,nm_p

          wr1(n)=wk1(l-1,n)

     &     +(t-tref(l-1))*(wk1(l,n)-wk1(l-1,n))/(tref(l)-tref(l-1))

          enddo

        goto 10

        endif

      enddo

c

  10  continue

c

      p_log=dlog(pre)

      call POLINT(pref_log,wr1,nm_p,p_log,xklog1,DY1)

      zkg_tp_int=dexp(xklog1)

c

      return

      end


      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate cooling rate (K/day) by CO rotational lines (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C  To calculate the cooling rate QJ [K s^-1] for CO rotational lines
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
c
      SUBROUTINE QTCOZ(TEM,PRE,RHO,QT,KM,MODE,grv,j_fels)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(MCO=72,KMX=175)                ! KMX = KM  *********************
      COMMON/FRENG/VMCO(MCO),STRCO(MCO),ENECO(MCO)
     &  ,alp_co(mco),xn_co(mco)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QT(KM),WK1(KMX),grv(km)
      IF(MODE.EQ.1) THEN
        OPEN(11,FILE='pluto_triton_dat/sanda_co_2014.dat',status='old')
        DO 15 I=1,MCO
        READ(11,125) IX,IX,FREQ,STR,alp_Lx,ENERGY,xn_x
        VMCO(I)=FREQ*100.0D0       ! line center wavenumber in m^-1
        STRCO(I)=STR
        ENECO(I)=ENERGY
        alp_co(i)=alp_Lx*100.0d0   ! pressure broadening half-width in m^-1
        xn_co(i)=xn_x              ! temperature dependence power index for alp_co
  15    CONTINUE
 125    format(i7,i2,f10.4,e10.3,f7.4,f10.4,f5.2)
        CLOSE(UNIT=11)
      ENDIF
      DO 20 K=1,KM
      QT(K)=0.0D0
  20  CONTINUE
      DO 50 M=1,MCO
      VMJ=VMCO(M)
      STRJ=STRCO(M)
      ENGJ=ENECO(M)
      alp_j=alp_co(m)
      xn_j=xn_co(m)
      CALL QJ_line_co(TEM,PRE,RHO,WK1,KM,VMJ,STRJ,ENGJ,grv
     &  ,alp_j,xn_j,j_fels)
      DO 40 K=1,KM
      QT(K)=QT(K)-WK1(K)      ! total CO cooling rate in [K s^-1]
  40  CONTINUE
  50  CONTINUE  
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate cooling rate (K/day) by HCN rotational lines (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C  To calculate the cooling rate QJ [K s^-1] for HCN rotational lines
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
c
      SUBROUTINE QTHCNZ(TEM,PRE,RHO,QT,KM,MODE,grv,j_fels)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(MHCN=141,KMX=175)                ! KMX = KM  *********************
      COMMON/FRENG2/VMHCN(MHCN),STRHCN(MHCN),ENEHCN(MHCN)
     &  ,alp_hcn(mhcn),xn_hcn(mhcn)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QT(KM),WK1(KMX),grv(km)
 125    format(i7,i2,f10.4,e10.3,f7.4,f10.4,f5.2)
      IF(MODE.EQ.1) THEN
        OPEN(11,FILE='pluto_triton_dat/sanda_hcn_2014.dat',status='old')
        DO 16 I=1,MHCN
        READ(11,125) IX,IX,FREQ,STR,alp_Lx,ENERGY,xn_x
        VMHCN(I)=FREQ*100.0D0       ! line center wavenumber in m^-1
        STRHCN(I)=STR
        ENEHCN(I)=ENERGY
        alp_hcn(i)=alp_Lx*100.0d0   ! pressure broadening half-width in m^-1
        xn_hcn(i)=xn_x              ! temperature dependence power index for alp_hcn
  16    CONTINUE
        CLOSE(UNIT=11)
      ENDIF
      DO 20 K=1,KM
      QT(K)=0.0D0
  20  CONTINUE
      DO 50 M=1,MHCN
      VMJ=VMHCN(M)
      STRJ=STRHCN(M)
      ENGJ=ENEHCN(M)
      alp_j=alp_hcn(m)
      xn_j=xn_hcn(m)
      CALL QJ_line_hcn(TEM,PRE,RHO,WK1,KM,VMJ,STRJ,ENGJ,grv
     &  ,alp_j,xn_j,j_fels)
      DO 40 K=1,KM
      QT(K)=QT(K)-WK1(K)      ! total HCN cooling rate in [K s^-1]
  40  CONTINUE
  50  CONTINUE  
      RETURN
      END

 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  To calculate cooling rate (K/day) by H2O rotational lines (K/s)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C  To calculate the cooling rate QJ [K s^-1] for H2O rotational lines
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
c
      SUBROUTINE QTH2OZ(TEM,PRE,RHO,QT,KM,MODE,grv,j_fels)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(MH2O=197,KMX=175)                ! KMX = KM  *********************
      COMMON/FRENG3/VMH2O(MH2O),STRH2O(MH2O),ENEH2O(MH2O)
     &  ,alp_h2o(mh2o),xn_h2o(mh2o)
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QT(KM),WK1(KMX),grv(km)
 125    format(i7,i2,f10.4,e10.3,f7.4,f10.4,f5.2)
      IF(MODE.EQ.1) THEN
        OPEN(11,FILE='pluto_triton_dat/sanda_h2o_2014.dat',status='old')
        DO 16 I=1,MH2O
        READ(11,125) IX,IX,FREQ,STR,alp_Lx,ENERGY,xn_x
        VMH2O(I)=FREQ*100.0D0       ! line center wavenumber in m^-1
        STRH2O(I)=STR
        ENEH2O(I)=ENERGY
        alp_h2o(i)=alp_Lx*100.0d0   ! pressure broadening half-width in m^-1
        xn_h2o(i)=xn_x              ! temperature dependence power index for alp_h2o
  16    CONTINUE
        CLOSE(UNIT=11)
      ENDIF
      DO 20 K=1,KM
      QT(K)=0.0D0
  20  CONTINUE
      DO 50 M=1,MH2O
      VMJ=VMH2O(M)
      STRJ=STRH2O(M)
      ENGJ=ENEH2O(M)
      alp_j=alp_h2o(m)
      xn_j=xn_h2o(m)
      CALL QJ_line_h2o(TEM,PRE,RHO,WK1,KM,VMJ,STRJ,ENGJ,grv
     &  ,alp_j,xn_j,j_fels)
      DO 40 K=1,KM
      QT(K)=QT(K)-WK1(K)      ! total H2O cooling rate in [K s^-1]
  40  CONTINUE
  50  CONTINUE  
      RETURN
      END


C  To calculate the cooling rate QJ [K s^-1] for a rotational line J --> J-1
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
C  VMJ [m^-1] is the frequency, (STRJ,ENGJ) is used to derive the line strength
c  Revised on 2014.09.19 to include pressure broadening - voigt line profile
c
      SUBROUTINE QJ_line_co(TEM,PRE,RHO,QJ,KM,VMJ,STRJ,ENGJ
     &  ,grv,alp_j,xn_j,j_fels)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175)                     ! KMX = KM  *********************
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QJ(KM),grv(km)
      DIMENSION ALDX(KMX),STX(KMX),BLACX(KMX),H00(KMX)
     &  ,GAM1(KMX,KMX),CURT(KMX,KMX),alpx(kmx)
      PI=3.14159265D0
      VCMJ=VMJ/100.0D0
      DO 10 K=1,KM
      ALDX(K)=alpd_co(VMJ,TEM(K))      ! Doppler half-width in m^-1
      STX(K)=strt7_co(TEM(K),STRJ,ENGJ,VCMJ)  ! line strength in m^-1 / (kg m^-2)
      BLACX(K)=BLACM(VMJ,TEM(K))      ! Planck function in [W/(m*m) * (m^-1)^-1 * Sr^-1]
      alpx(k)=alp_j*(296.0d0/tem(k))**xn_j*(pre(k)/1.01325d5)
  10  CONTINUE
      KMM=KM-1
      do i=1,km
      do j=1,km
        if(i.eq.j) then  ! assuming isothermal for escape function
          wdx=ALDX(i)
          wdx_b=ALDX(i)
          wlx=alpx(i)
          wlx_b=alpx(i)
          if(i.eq.1) xxm=rho(1)*0.25d0*(pre(1)-pre(2))*stx(1)/grv(1)
          if(i.eq.km) 
     &         xxm=rho(km)*0.25d0*(pre(kmm)-pre(km))*stx(km)/grv(km)
          if(i.ne.1.and.i.ne.km) 
     &         xxm=rho(i)*0.125d0*(pre(i-1)-pre(i+1))*stx(i)/grv(i)
        goto 25
        endif
      wdx=ALDX(i)
      wlx=alpx(i)
        wdx_b=(ALDX(i)+ALDX(j))/2.0d0
        wlx_b=(alpx(i)+alpx(j))/2.0d0
        xxm=0.0d0
        k1=min0(i,j)
        k2=max0(i,j)-1
        do k=k1,k2
        xxm=xxm+rho(k)*(pre(k)-pre(k+1))*stx(k)/grv(k)
        enddo
  25  continue
      if(j_fels.eq.1) gam1(i,j)=gamma_fels(xxm,wdx,wlx,wdx_b,wlx_b)    ! less accurate, fatest
      if(j_fels.eq.2) gam1(i,j)=gamma_zhu_g(xxm,wdx,wlx,wdx_b,wlx_b)   ! more accurate, medium speed
      if(j_fels.eq.3) gam1(i,j)=gamma_zhu(xxm,wdx,wlx,wdx_b,wlx_b)     ! most accurate, slowest
      enddo
      enddo
c
C  To calculate the Curtis matrices +++++++++++++++++++++++++++++++++++++++++
c
      DO 40 L1=1,KM
      H00(L1)=2.0D0*PI*RHO(L1)*STX(L1)/cpair
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

C  To calculate the cooling rate QJ [K s^-1] for a rotational line J --> J-1
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
C  VMJ [m^-1] is the frequency, (STRJ,ENGJ) is used to derive the line strength
c  Revised on 2014.09.19 to include pressure broadening - voigt line profile
c
      SUBROUTINE QJ_line_hcn(TEM,PRE,RHO,QJ,KM,VMJ,STRJ,ENGJ
     &  ,grv,alp_j,xn_j,j_fels)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175)                     ! KMX = KM  *********************
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QJ(KM),grv(km)
      DIMENSION ALDX(KMX),STX(KMX),BLACX(KMX),H00(KMX)
     &  ,GAM1(KMX,KMX),CURT(KMX,KMX),alpx(kmx)
      PI=3.14159265D0
      VCMJ=VMJ/100.0D0
      DO 10 K=1,KM
      ALDX(K)=alpd_hcn(VMJ,TEM(K))      ! Doppler half-width in m^-1
      STX(K)=strt7_hcn(TEM(K),STRJ,ENGJ,VCMJ)  ! line strength in m^-1 / (kg m^-2)
      BLACX(K)=BLACM(VMJ,TEM(K))      ! Planck function in [W/(m*m) * (m^-1)^-1 * Sr^-1]
      alpx(k)=alp_j*(296.0d0/tem(k))**xn_j*(pre(k)/1.01325d5)
  10  CONTINUE
      KMM=KM-1
      do i=1,km
      do j=1,km
        if(i.eq.j) then  ! assuming isothermal for escape function
          wdx=ALDX(i)
          wdx_b=ALDX(i)
          wlx=alpx(i)
          wlx_b=alpx(i)
          if(i.eq.1) xxm=rho(1)*0.25d0*(pre(1)-pre(2))*stx(1)/grv(1)
          if(i.eq.km) 
     &         xxm=rho(km)*0.25d0*(pre(kmm)-pre(km))*stx(km)/grv(km)
          if(i.ne.1.and.i.ne.km) 
     &         xxm=rho(i)*0.125d0*(pre(i-1)-pre(i+1))*stx(i)/grv(i)
        goto 25
        endif
      wdx=ALDX(i)
      wlx=alpx(i)
        wdx_b=(ALDX(i)+ALDX(j))/2.0d0
        wlx_b=(alpx(i)+alpx(j))/2.0d0
        xxm=0.0d0
        k1=min0(i,j)
        k2=max0(i,j)-1
        do k=k1,k2
        xxm=xxm+rho(k)*(pre(k)-pre(k+1))*stx(k)/grv(k)
        enddo
  25  continue
      if(j_fels.eq.1) gam1(i,j)=gamma_fels(xxm,wdx,wlx,wdx_b,wlx_b)    ! less accurate, fatest
      if(j_fels.eq.2) gam1(i,j)=gamma_zhu_g(xxm,wdx,wlx,wdx_b,wlx_b)   ! more accurate, medium speed
      if(j_fels.eq.3) gam1(i,j)=gamma_zhu(xxm,wdx,wlx,wdx_b,wlx_b)     ! most accurate, slowest
      enddo
      enddo
c
C  To calculate the Curtis matrices +++++++++++++++++++++++++++++++++++++++++
c
      DO 40 L1=1,KM
      H00(L1)=2.0D0*PI*RHO(L1)*STX(L1)/cpair
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


C  To calculate the cooling rate QJ [K s^-1] for a rotational line J --> J-1
C  TEM [K], PRE [pa], RHO (mass mixing ratio) = input, ZM(KM) > ZM(1)
C  VMJ [m^-1] is the frequency, (STRJ,ENGJ) is used to derive the line strength
c  Revised on 2014.09.19 to include pressure broadening - voigt line profile
c
      SUBROUTINE QJ_line_h2o(TEM,PRE,RHO,QJ,KM,VMJ,STRJ,ENGJ
     &  ,grv,alp_j,xn_j,j_fels)
      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      PARAMETER(KMX=175)                     ! KMX = KM  *********************
      common/parair/r00,capm,wtmol,rcpr,rair,cvair,cpair,g00,rsunp0,
     &             rsunp
      DIMENSION TEM(KM),PRE(KM),RHO(KM),QJ(KM),grv(km)
      DIMENSION ALDX(KMX),STX(KMX),BLACX(KMX),H00(KMX)
     &  ,GAM1(KMX,KMX),CURT(KMX,KMX),alpx(kmx)
      PI=3.14159265D0
      VCMJ=VMJ/100.0D0
      DO 10 K=1,KM
      ALDX(K)=alpd_h2o(VMJ,TEM(K))      ! Doppler half-width in m^-1
      STX(K)=strt7_h2o(TEM(K),STRJ,ENGJ,VCMJ)  ! line strength in m^-1 / (kg m^-2)
      BLACX(K)=BLACM(VMJ,TEM(K))      ! Planck function in [W/(m*m) * (m^-1)^-1 * Sr^-1]
      alpx(k)=alp_j*(296.0d0/tem(k))**xn_j*(pre(k)/1.01325d5)
  10  CONTINUE
      KMM=KM-1
      do i=1,km
      do j=1,km
        if(i.eq.j) then  ! assuming isothermal for escape function
          wdx=ALDX(i)
          wdx_b=ALDX(i)
          wlx=alpx(i)
          wlx_b=alpx(i)
          if(i.eq.1) xxm=rho(1)*0.25d0*(pre(1)-pre(2))*stx(1)/grv(1)
          if(i.eq.km) 
     &         xxm=rho(km)*0.25d0*(pre(kmm)-pre(km))*stx(km)/grv(km)
          if(i.ne.1.and.i.ne.km) 
     &         xxm=rho(i)*0.125d0*(pre(i-1)-pre(i+1))*stx(i)/grv(i)
        goto 25
        endif
      wdx=ALDX(i)
      wlx=alpx(i)
        wdx_b=(ALDX(i)+ALDX(j))/2.0d0
        wlx_b=(alpx(i)+alpx(j))/2.0d0
        xxm=0.0d0
        k1=min0(i,j)
        k2=max0(i,j)-1
        do k=k1,k2
        xxm=xxm+rho(k)*(pre(k)-pre(k+1))*stx(k)/grv(k)
        enddo
  25  continue
      if(j_fels.eq.1) gam1(i,j)=gamma_fels(xxm,wdx,wlx,wdx_b,wlx_b)    ! less accurate, fatest
      if(j_fels.eq.2) gam1(i,j)=gamma_zhu_g(xxm,wdx,wlx,wdx_b,wlx_b)   ! more accurate, medium speed
      if(j_fels.eq.3) gam1(i,j)=gamma_zhu(xxm,wdx,wlx,wdx_b,wlx_b)     ! most accurate, slowest     
      enddo
      enddo
c
C  To calculate the Curtis matrices +++++++++++++++++++++++++++++++++++++++++
c
      DO 40 L1=1,KM
      H00(L1)=2.0D0*PI*RHO(L1)*STX(L1)/cpair
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
//parair
/
C  To calculate the Doppler half-width in m^-1, VCM [m^-1] and

C  T [K] are the reference line center and temperature, respectively

c

      FUNCTION alpd_co(VM,T)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

      VX=VM/25.0D2

      TX=T/40.0D0

      alpd_co=1.28016D-3*VX*DSQRT(TX)

      RETURN

      END



C  To calculate the Doppler half-width in m^-1, VCM [m^-1] and

C  T [K] are the reference line center and temperature, respectively

c

      FUNCTION alpd_hcn(VM,T)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

      VX=VM/25.0D2

      TX=T/40.0D0

      alpd_hcn=1.30365d-3*VX*DSQRT(TX)

      RETURN

      END

      

C  To calculate the Doppler half-width in m^-1, VCM [m^-1] and

C  T [K] are the reference line center and temperature, respectively

c

      FUNCTION alpd_h2o(VM,T)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
      VX=VM/25.0D2

      TX=T/40.0D0

      alpd_h2o=1.59664d-3*VX*DSQRT(TX)

      RETURN

      END

      
C  To calculate CO line line strength in m^-1 / (kg m^-2)

C  STR = line strength at 296 K in  cm^-1 /(molec cm^-2)

C  VCMJ = line center frequency in cm^-1

c

      FUNCTION strt7_co(T,STR,ENG,VCMJ)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

      TS=296.0D0

      XM=28.0D0                            ! molecular weight, isotops effect has been included in str

      QR=T/TS                              ! rotational partition function

      FAC1=1.439D0*ENG*(T-TS)/(T*TS)

      FAC7=(1.0D0-DEXP(-1.439D0*VCMJ/T))   ! induced emission becomes important

     &    /(1.0D0-DEXP(-1.439D0*VCMJ/TS))  !   for rotational lines at lower T

      IF(FAC1.LE.0.0D0) FAC2=-DMIN1(200.0D0,-FAC1)

      IF(FAC1.GE.0.0D0) FAC2=DMIN1(200.0D0,FAC1)

      BFAC=DEXP(FAC2)

      strt7_co=(STR/QR)*BFAC*FAC7

      FAC2=6.023D23/XM                     ! conversion factor of molec g^-1

      strt7_co=strt7_co*FAC2                     ! line strength in cm^-1 / (g cm^-2)

      strt7_co=strt7_co*10.0D0                   ! line strength in m^-1 / (kg m^-2)

      RETURN

      END



C  To calculate HCN line line strength in m^-1 / (kg m^-2)

C  STR = line strength at 296 K in  cm^-1 /(molec cm^-2)

C  VCMJ = line center frequency in cm^-1

c

      FUNCTION strt7_hcn(T,STR,ENG,VCMJ)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

      TS=296.0D0

      XM=27.0D0                            ! molecular weight, isotops effect has been included in str

      QR=T/TS                              ! rotational partition function

      FAC1=1.439D0*ENG*(T-TS)/(T*TS)

      FAC7=(1.0D0-DEXP(-1.439D0*VCMJ/T))   ! induced emission becomes important

     &    /(1.0D0-DEXP(-1.439D0*VCMJ/TS))  !   for rotational lines at lower T

      IF(FAC1.LE.0.0D0) FAC2=-DMIN1(200.0D0,-FAC1)

      IF(FAC1.GE.0.0D0) FAC2=DMIN1(200.0D0,FAC1)

      BFAC=DEXP(FAC2)

      strt7_hcn=(STR/QR)*BFAC*FAC7

      FAC2=6.023D23/XM                     ! conversion factor of molec g^-1

      strt7_hcn=strt7_hcn*FAC2                     ! line strength in cm^-1 / (g cm^-2)

      strt7_hcn=strt7_hcn*10.0D0                   ! line strength in m^-1 / (kg m^-2)

      RETURN

      END



C  To calculate H2O line line strength in m^-1 / (kg m^-2)

C  STR = line strength at 296 K in  cm^-1 /(molec cm^-2)

C  VCMJ = line center frequency in cm^-1

c

      FUNCTION strt7_h2o(T,STR,ENG,VCMJ)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

      TS=296.0D0

      XM=18.0D0                            ! molecular weight, isotops effect has been included in str

      QR=(T/TS)**1.5                       ! rotational partition function

      FAC1=1.439D0*ENG*(T-TS)/(T*TS)

      FAC7=(1.0D0-DEXP(-1.439D0*VCMJ/T))   ! induced emission becomes important

     &    /(1.0D0-DEXP(-1.439D0*VCMJ/TS))  !   for rotational lines at lower T

      IF(FAC1.LE.0.0D0) FAC2=-DMIN1(200.0D0,-FAC1)

      IF(FAC1.GE.0.0D0) FAC2=DMIN1(200.0D0,FAC1)

      BFAC=DEXP(FAC2)

      strt7_h2o=(STR/QR)*BFAC*FAC7

      FAC2=6.023D23/XM                     ! conversion factor of molec g^-1

      strt7_h2o=strt7_h2o*FAC2                     ! line strength in cm^-1 / (g cm^-2)

      strt7_h2o=strt7_h2o*10.0D0                   ! line strength in m^-1 / (kg m^-2)

      RETURN

      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C

C  Some common subroutines

C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



C  Calcualte exponential integral from polynomial and rational 

C  approximation E2(z)=[exp(-z)-zE1(z)]=exp(-z)[1-exp(z)zE1(z)]

c

      FUNCTION ENZD2(Z)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

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

C  Planck function [W/(m*m) * (m^-1)^-1 * Sr^-1] at T[K] & wavenumber VM [m^-1]

c

      FUNCTION BLACM(VM,T)

      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)

      BLACM=1.191D-16*VM**3/(DEXP(dmin1(1.438786D-2*VM/T,2.0d2))-1.0D0)

      RETURN

      END
      
C  Polynomial interpolation from "Numerical Recipes"
c AMZ already have this subroutine
C      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C      IMPLICIT real*8(A-H,O-Z), INTEGER(I-N)
C      PARAMETER (NMAX=10) 
C      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
C      NS=1
C      DIF=DABS(X-XA(1))
C      DO 11 I=1,N 
C        DIFT=DABS(X-XA(I))
C        IF (DIFT.LT.DIF) THEN
C          NS=I
C          DIF=DIFT
C        ENDIF
C        C(I)=YA(I)
C        D(I)=YA(I)
C11    CONTINUE
C      Y=YA(NS)
C      NS=NS-1
C      DO 13 M=1,N-1
C        DO 12 I=1,N-M
C          HO=XA(I)-X
C          HP=XA(I+M)-X
C          W=C(I+1)-D(I)
C          DEN=HO-HP
C          IF(DEN.EQ.0.0D0)PAUSE
C          DEN=W/DEN
C          D(I)=HP*DEN
C          C(I)=HO*DEN
C12      CONTINUE
C        IF (2*NS.LT.N-M)THEN
C          DY=C(NS+1)
C        ELSE
C          DY=D(NS)
C          NS=NS-1
C        ENDIF
C        Y=Y+DY
C13    CONTINUE
c      RETURN
C      END

  
CC  
CC  To invert the matrix D,  M=N+1. Mth row can be any values. Double precision
c
      SUBROUTINE INVERT(D,N,M)
      IMPLICIT  real*8(A-H,O-Z), INTEGER(I-N)
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



cc subroutines used to calculate the escape function of a voigt profile, 

cc   added by Xun Zhu, on 2014.09.19

c

cccccccccccccccccccccccccc gamma_fels ccccccccccccccccccccccccccccc

CC

CC  To calculate the escape function of a Voigt line based Fels approximation (1979) 

cc    -- default choice for speed

cc

      function gamma_fels(xm,wd,wl,wd_b,wl_b)

      implicit real*8(a-h,o-z), INTEGER(I-N)

      pi=3.14159265359d0

      eta=2.10d0  ! ranging from 1.8 to 2.3, default value = 2.1

c

      v0=4.0d0*wl/pi+1.5d0*wd

      capc=(0.5d0/v0)-wl/(pi*v0*v0)

        v0_b=4.0d0*wl_b/pi+1.5d0*wd_b

        capc_b=(0.5d0/v0_b)-wl_b/(pi*v0_b*v0_b)

      v0_til=(v0+v0_b)/2.0d0

c

        capj1=v0_til*capc*enzd2(xm*capc_b)

        tau=eta*xm*wl_b/pi

        x0=dsqrt(tau)/v0_til

        capj2=wl*erfz1(x0)/(2.0d0*dsqrt(pi*tau))

      gamma_fels=2.0d0*(capj1+capj2)

      return

      end



cccccccccccccccccccccccccc gamma_zhu ccccccccccccccccccccccccccccc

CC

CC  To calculate the escape function of a Voigt line, default choice for accuracy

cc    -- numerical exact voigt in both core and wing 

cc    -- very accurate and much better than gamma_fels but somewhat slower

cc   

      function gamma_zhu(xm,wd,wl,wd_b,wl_b)

      implicit real*8(a-h,o-z), INTEGER(I-N)

      pi=3.14159265359d0

      v1=(wl+wd+wl_b+wd_b)*2.0d0

c

      im=8  ! im=8 is the default setting, leading overall small errors

      dv1=v1/dfloat(im)

      capi12=0.0d0

      do i=1,im

      vx=dv1*(dfloat(i)-0.5d0)

          fx1=voigt_fs85(vx,wl,wd)

          fx2=enzd2(xm*voigt_fs85(vx,wl_b,wd_b))

      capi12=capi12+fx1*fx2

      enddo

      capi12=capi12*dv1

c

      jm=12

      u1_til=1.0d0/v1

      du1=u1_til/dfloat(jm)

      capi3=0.0d0

      do j=1,jm

      ux=du1*(dfloat(j)-0.5d0)

      vx=1.0d0/ux

          fx1=voigt_fs85(vx,wl,wd)

          fx2=enzd2(xm*voigt_fs85(vx,wl_b,wd_b))

      capi3=capi3+fx1*fx2/(ux*ux)

      enddo

      capi3=capi3*du1

c

      gamma_zhu=2.0d0*(capi12+capi3)

c

      return

      end





cccccccccccccccccccccccccc gamma_zhu_g ccccccccccccccccccccccccccccc

CC

CC  To calculate the escape function of a Voigt line,

cc    -- numerical exact voigt in both core and wing 

cc    -- very accurate and much better than gamma_fels but somewhat slower

cc    -- improved in speed based on gamma_zhu by adopting Gaussian quadratures

cc    -- accuracy and speed are between gamma_fels and gamma_zhu

cc   

      function gamma_zhu_g(xm,wd,wl,wd_b,wl_b)

      implicit real*8(a-h,o-z), INTEGER(I-N)
      dimension xa(4),wa(4),xb(6),wb(6)

      pi=3.14159265359d0

        xa(1)=(1.0d0-0.8611363116d0)/2.0d0

        xa(2)=(1.0d0-0.3399810436d0)/2.0d0

        xa(3)=(1.0d0+0.3399810436d0)/2.0d0

        xa(4)=(1.0d0+0.8611363116d0)/2.0d0

          wa(1)=0.3478548451d0/2.0d0

          wa(2)=0.6521451549d0/2.0d0

          wa(3)=0.6521451549d0/2.0d0

          wa(4)=0.3478548451d0/2.0d0

        xb(1)=(1.0d0-0.9324695142d0)/2.0d0

        xb(2)=(1.0d0-0.6612093865d0)/2.0d0

        xb(3)=(1.0d0-0.2386191861d0)/2.0d0

        xb(4)=(1.0d0+0.2386191861d0)/2.0d0

        xb(5)=(1.0d0+0.6612093865d0)/2.0d0

        xb(6)=(1.0d0+0.9324695142d0)/2.0d0

          wb(1)=0.1713244924d0/2.0d0

          wb(2)=0.3607615730d0/2.0d0

          wb(3)=0.4679139346d0/2.0d0

          wb(4)=0.4679139346d0/2.0d0

          wb(5)=0.3607615730d0/2.0d0

          wb(6)=0.1713244924d0/2.0d0

c

      v1=(wl+wd+wl_b+wd_b)*1.5d0

c

      im=4  ! im=4 is the default setting, leading overall small errors

      capi12=0.0d0

      do i=1,im

      vx=v1*xa(i)

          fx1=voigt_fs85(vx,wl,wd)

          fx2=enzd2(xm*voigt_fs85(vx,wl_b,wd_b))

      capi12=capi12+(wa(i)*v1)*fx1*fx2

      enddo

c

      jm=6  ! jm=6 is the default setting, leading overall small errors

      u1=1.0d0/v1

      capi3=0.0d0

      do j=1,jm

      ux=u1*xb(j)

      vx=1.0d0/ux

          fx1=voigt_fs85(vx,wl,wd)

          fx2=enzd2(xm*voigt_fs85(vx,wl_b,wd_b))

      capi3=capi3+(wb(j)*u1)*fx1*fx2/(ux*ux)

      enddo

c

      gamma_zhu_g=2.0d0*(capi12+capi3)

c

      return

      end









cccccccccccccccccccccccccc voigt_fs85 ccccccccccccccccccccccccccccc

CC

CC  To calculate the Voigt function with the approximation by 

CC  Fomichev and Shved (1985) [J. Atmos. Terr. Phys., 47, 1049]. 

CC  

      FUNCTION voigt_fs85(V,WL,WD)

      implicit real*8(a-h,o-z), INTEGER(I-N)

      WDD2=2.772588722D0*WD*WD

      XFAC1=WL+DSQRT(WL*WL+WDD2)

      ALPV=0.5D0*XFAC1+0.05D0*WL*(1.0D0-2.0D0*WL/XFAC1)

      ZETA=WL/ALPV

      ETA=V/ALPV

      ETA2=ETA*ETA

      TERM1=(0.4697186393D0/ALPV)*(1.0D0-ZETA)

     &      *DEXP(-DMIN1(200.0D0,0.6931471806D0*ETA2))

      TERM2=ZETA/(3.141592654D0*ALPV*(1.0D0+ETA2))

      TERM3=(ZETA*(1.0D0-ZETA)/(3.141592654D0*ALPV))

     &       *(3.164042561D0+ZETA)

     &       *(0.066D0*DEXP(-DMIN1(200.0D0,0.4D0*ETA2))

     &         -1.0D0/(40.0D0-5.5D0*ETA2+ETA2*ETA2))

      voigt_fs85=TERM1+TERM2-TERM3

      RETURN

      END





c To calculate the error function according to [JQSRT,1990(43), P511]

c

      function erfz1(xx)

      implicit real*8(a-h,o-z), INTEGER(I-N)

      x=dabs(xx)

      xsq=x*x

      fac1=1.0d0/(1.16358d0*x+dsqrt(1.0d0+0.37072d0*xsq))  ! exp(x^2)*erfc(x)

      erfcz2x=fac1*dexp(-dmin1(xsq,200.0d0))

      if(xx.ge.0.0d0) then

      erfz1=1.0d0-erfcz2x      ! erf(z) = 1 - erfc(z)

      else

      erfz1=erfcz2x-1.0d0      ! erf(-z) = - erf(z)

      endif

      return

      end



c To calculate the complementary error function according to [JQSRT,1990(43), P511]

c

      function erfcz2(xx)

      implicit real*8(a-h,o-z), INTEGER(I-N)

      x=dabs(xx)

      xsq=x*x

      fac1=1.0d0/(1.16358d0*x+dsqrt(1.0d0+0.37072d0*xsq))  ! exp(x^2)*erfc(x)

      erfcz2=fac1*dexp(-dmin1(xsq,200.0d0))

      if(xx.ge.0.0d0) return

      erfcz2=2.0d0-erfcz2   ! erfc(-z) = 2 - erfc(z)

      return

      end



c To calculate acosh(x) for an extremely large x, acosh(x)=ln(x+sqrt(x^2-1))

c

      function dacosh2(x,a,j)   ! x = a*a + 1

      implicit real*8(a-h,o-z), INTEGER(I-N)

c

      if(j.eq.1) then    ! x is the input

      dacosh2=dlog(x+dsqrt(x*x-1.0d0))

      return

      endif

c

      if(j.eq.2) then    ! a is the input

      asq=a*a

      x=asq+1.000001d0

      x2=(asq/x)**2+2.0d0*(a/x)**2

      dacosh2=dlog(x) + dlog(1.0d0+dsqrt(x2))

      return

      endif

c

      end 
