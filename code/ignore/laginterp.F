#include "CPP_OPTIONS.h"
      subroutine laginterp(pgref,pint,co2i,co2v,fzeroi,fzerov)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C  Lagrange interpolation (linear in log pressure) of the CO2 
C  k-coefficients in the pressure domain.  Subsequent use of these
C  values will use a simple linear interpolation in pressure.

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GRID_H.h"
#include "RADINC_H.h"

C  PGREF   - Pressure grid on which k-coefficients were computed.
C  PINT    - Interpolated pressure grid (finer grid spacing than 
C            PGREF)
C  CO2I    - CO2 IR k-coefficients on interpolated grid
C  CO2V    - CO2 visible k-coefficients on interpolated grid
C  FZEROI  - Fraction of zeros in each IR spectral interval
C  FZEROV  - Fraction of zeros in each visible spectral interval
C
C----------------------------------------------------------------------!

      integer  n, nt, np, nh, ng, nw, m, i,j,k,l
      real*8   co2i8(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8   co2v8(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8   pgref(L_NPREF)

      real*8   co2i(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8   co2v(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)

      real*8   fzeroi(L_NSPECTI)
      real*8   fzerov(L_NSPECTV)
 
      real*8   x, xi(4), yi(4), ans
      real*8   pint(L_PINT), pref(L_NPREF), p

      _RL   pin(L_PINT) 
      data pin /  -6.0, -5.8, -5.6, -5.4, -5.2,     
     &  -5.0, -4.8, -4.6, -4.4, -4.2,     
     &  -4.0, -3.8, -3.6, -3.4, -3.2,     
     &  -3.0, -2.8, -2.6, -2.4, -2.2,     
     &  -2.0, -1.8, -1.6, -1.4, -1.2,     
     &  -1.0, -0.8, -0.6, -0.4, -0.2,     
     &   0.0,  0.2,  0.4,  0.6,  0.8,     
     &   1.0,  1.2,  1.4,  1.6,  1.8,     
     &   2.0,  2.2,  2.4,  2.6,  2.8,     
     &   3.0,  3.2,  3.4,  3.6,  3.8,     
     &   4.0/                                  

C======================================================================!

C  Fill pint for output from this subroutine

      do n=1,L_PINT
        PINT(n) = PIN(n)
      end do

C  Take log of the reference pressures

      do n=1,L_NPREF
        pref(n) = LOG10(PGREF(n))
      end do


C amz Fortran runtime error: End of file
C amz when using this method.  Read in as text file instead.
C      open(20,file='CO2H2O_V_12_95_INTEL',form='unformatted')
C      read(20) co2v8
C      read(20) fzerov
C      close(20)
C
C      open(20,file='CO2H2O_IR_12_95_INTEL',form='unformatted')
C      read(20) co2i8
C      read(20) fzeroi
C      close(20)

      open(20,file='CO2H2O_V_12_95_INTEL.txt',form='formatted',
     &       status='old')
      DO i=1,L_NTREF
       DO j=1,L_NPREF
        DO k=1,L_REFH2O
         DO l=1,L_NSPECTV
          DO m=1,L_NGAUSS
           read(20,"(E26.22)") co2v8(i,j,k,l,m)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO l=1,L_NSPECTV
       read(20,"(f26.22)") fzerov(l)
      ENDDO
      close(20)

      open(20,file='CO2H2O_IR_12_95_INTEL.txt',form='formatted',
     &       status='old')
      DO i=1,L_NTREF
       DO j=1,L_NPREF
        DO k=1,L_REFH2O
         DO l=1,L_NSPECTI
          DO m=1,L_NGAUSS
           read(20,"(E26.22)") co2i8(i,j,k,l,m)
          ENDDO
         ENDDO
	ENDDO
       ENDDO
      ENDDO

      DO l=1,L_NSPECTI
       read(20,"(f26.22)") fzeroi(l)
      ENDDO
      close(20)


C  Take Log10 of the values - we interpolate the log10 of the values,
C  not the values themselves.   Smallest value is 1.0E-200.

      do nt=1,L_NTREF
        do np=1,L_NPREF
          do nh=1,L_REFH2O
            do ng = 1,L_NGAUSS

              do nw=1,L_NSPECTV
                if(co2v8(nt,np,nh,nw,ng).gt.1.0d-200) then
                  co2v8(nt,np,nh,nw,ng) = 
     &                  log10(co2v8(nt,np,nh,nw,ng))
                else
                  co2v8(nt,np,nh,nw,ng) = -200.0
                end if
              end do

              do nw=1,L_NSPECTI
                if(co2i8(nt,np,nh,nw,ng).gt.1.0d-200) then
                  co2i8(nt,np,nh,nw,ng) = 
     &              log10(co2i8(nt,np,nh,nw,ng))
                else
                  co2i8(nt,np,nh,nw,ng) = -200.0
                end if
              end do

            end do
          end do
        end do
      end do

C  Interpolate the values:  first the IR

      do nt=1,L_NTREF
        do nh=1,L_REFH2O
        do nw=1,L_NSPECTI
          do ng=1,L_NGAUSS

C  First, the initial interval (P=1e-6 to 1e-5)

            n = 1
            do m=1,5
              x     = pint(m)
              xi(1) = pref(n)
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2i8(nt,n,nh,nw,ng)
              yi(2) = co2i8(nt,n+1,nh,nw,ng)
              yi(3) = co2i8(nt,n+2,nh,nw,ng)
              yi(4) = co2i8(nt,n+3,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2i(nt,m,nh,nw,ng) = 10.0**ans
            end do

            do n=2,L_NPREF-2
              do m=1,5
                i     = (n-1)*5+m
                x     = pint(i)
                xi(1) = pref(n-1)
                xi(2) = pref(n)
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2i8(nt,n-1,nh,nw,ng)
                yi(2) = co2i8(nt,n,nh,nw,ng)
                yi(3) = co2i8(nt,n+1,nh,nw,ng)
                yi(4) = co2i8(nt,n+2,nh,nw,ng)
                call lagrange(x,xi,yi,ans)
                co2i(nt,i,nh,nw,ng) = 10.0**ans
              end do
            end do

C  Now, get the last interval (P=1e+3 to 1e+4)

            n = L_NPREF-1

            do m=1,5
              i     = (n-1)*5+m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n)
              xi(4) = pref(n+1)
              yi(1) = co2i8(nt,n-2,nh,nw,ng)
              yi(2) = co2i8(nt,n-1,nh,nw,ng)
              yi(3) = co2i8(nt,n,nh,nw,ng)
              yi(4) = co2i8(nt,n+1,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2i(nt,i,nh,nw,ng) = 10.0**ans
            end do

C  Fill the last pressure point

            co2i(nt,L_PINT,nh,nw,ng) = 
     &              10.0**co2i8(nt,L_NPREF,nh,nw,ng)

          end do
        end do
	end do
      end do

C  Interpolate the values:  now the visible

      do nt=1,L_NTREF
        do nh=1,L_REFH2O
        do nw=1,L_NSPECTV
          do ng=1,L_NGAUSS

C  First, the initial interval (P=1e-6 to 1e-5)

            n = 1
            do m=1,5
              x     = pint(m)
              xi(1) = pref(n)
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2v8(nt,n,nh,nw,ng)
              yi(2) = co2v8(nt,n+1,nh,nw,ng)
              yi(3) = co2v8(nt,n+2,nh,nw,ng)
              yi(4) = co2v8(nt,n+3,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2v(nt,m,nh,nw,ng) = 10.0**ans
            end do

            do n=2,L_NPREF-2
              do m=1,5
                i     = (n-1)*5+m
                x     = pint(i)
                xi(1) = pref(n-1)
                xi(2) = pref(n)
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2v8(nt,n-1,nh,nw,ng)
                yi(2) = co2v8(nt,n,nh,nw,ng)
                yi(3) = co2v8(nt,n+1,nh,nw,ng)
                yi(4) = co2v8(nt,n+2,nh,nw,ng)
                call lagrange(x,xi,yi,ans)
                co2v(nt,i,nh,nw,ng) = 10.0**ans
              end do
            end do

C  Now, get the last interval (P=1e+3 to 1e+4)

            n = L_NPREF-1

            do m=1,5
              i     = (n-1)*5+m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n)
              xi(4) = pref(n+1)
              yi(1) = co2v8(nt,n-2,nh,nw,ng)
              yi(2) = co2v8(nt,n-1,nh,nw,ng)
              yi(3) = co2v8(nt,n,nh,nw,ng)
              yi(4) = co2v8(nt,n+1,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2v(nt,i,nh,nw,ng) = 10.0**ans
            end do

C  Fill the last pressure point

            co2v(nt,L_PINT,nh,nw,ng) = 
     &            10.0**co2v8(nt,L_NPREF,nh,nw,ng)

          end do
        end do
        end do
      end do

      return
      end

