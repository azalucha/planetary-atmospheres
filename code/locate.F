      SUBROUTINE LOCATE(xx,n,x,j)

      IMPLICIT NONE

      INTEGER j,n
      DOUBLE PRECISION x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
  10  IF(ju-jl .gt. 1) THEN
       jm=(ju+jl)/2
       IF((xx(n) .ge. xx(1)) .eqv. (x .ge. xx(jm))) THEN
        jl=jm
       ELSE
        ju=jm
       ENDIF
      GOTO 10
      ENDIF
      IF(x .eq. xx(1)) THEN
       j=1
      ELSE IF(x .eq. xx(n)) THEN
       j=n-1
      ELSE
       j=jl
      ENDIF

      RETURN
      END
