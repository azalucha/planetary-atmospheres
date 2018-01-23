      SUBROUTINE TRIDAG(a,b,c,r,u,n)

      IMPLICIT NONE

C Inputs: a,b,c,r,n
C Output: u


      INTEGER n,NMAX,j
      DOUBLE PRECISION a(n),b(n),c(n),r(n),u(n)
      PARAMETER(NMAX=500)
      DOUBLE PRECISION bet,gam(NMAX)

      IF(b(1) .EQ. 0.) THEN
       write(*,*) 'TRIDIAG: REWRITE EQUATIONS'
       STOP
      ENDIF

      bet=b(1)
      u(1)=r(1)/bet
      DO j=2,n
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j)*gam(j)
       IF(bet .EQ. 0.) THEN
        write(*,*) 'TRIDIAG FAILED'
        STOP
       ENDIF
       u(j)=(r(j)-a(j)*u(j-1))/bet
      ENDDO

      DO j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
      ENDDO

      RETURN
      END
