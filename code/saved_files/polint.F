
      SUBROUTINE POLINT(xa,ya,n,x,y,dy)

      IMPLICIT NONE


C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     xa            :: array of independent variable
C     ya            :: array of dependent variable
C     n             :: length of xa and ya; n-1 is polynomial degree
C     x             :: value you want to interpolate at
C     y             :: interpolated value (output)
C     dy            :: error estimate

      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=30)

      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)


C     begin equations

      ns=1
      dif=DABS(x-xa(1))
      DO i=1,n
       dift=DABS(x-xa(i))
       IF (dift .LT. dif) THEN
        ns=i
        dif=dift
       ENDIF
       c(i)=ya(i)
       d(i)=ya(i)
      ENDDO
      y=ya(ns)
      ns=ns-1
      DO m=1,n-1
       DO i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        IF(den .EQ. 0.) write(*,*) 'failure in polint'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
       ENDDO
       IF (2*ns .LT. n-m) THEN
        dy=c(ns+1)
       ELSE
        dy=d(ns)
        ns=ns-1
       ENDIF
       y=y+dy
      ENDDO

      RETURN
      END
