#include "CPP_OPTIONS.h"
      SUBROUTINE INTERP_TOTAL_TO_RC(p_heating,total,
     &                              L_NLAYRAD,myheating)

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"

      INTEGER L_NLAYRAD,k,pindex,k1,pinterpindexstart
      _RL p_heating(L_NLAYRAD),tempout,total(L_NLAYRAD)
      _RL psub(4),tempsub(4),dy,myheating(Nr),rCmb(Nr)
      LOGICAL dointerp

C      DO k=1,L_NLAYRAD
C      write(*,*) 'k=',k,'p_heating=',p_heating(k),'total=',total(k)
C      ENDDO

      DO k=1,Nr
       rCmb(k)=rC(k)/100.
      ENDDO

      DO k=1,Nr
       dointerp=.TRUE.
       CALL LOCATE(p_heating,L_NLAYRAD,rCmb(k),pindex)


       IF(pindex .EQ. 1) THEN
        pinterpindexstart=1
       ELSE IF (pindex .EQ. L_NLAYRAD-1) THEN
        pinterpindexstart=L_NLAYRAD-3
       ELSE IF (pindex .EQ. L_NLAYRAD) THEN
C        write(*,*) 'Locate: p out of bounds high'
        myheating(k)=total(L_NLAYRAD)
        dointerp=.FALSE.
       ELSE IF (pindex .EQ. 0) THEN
C        write(*,*) 'Locate: p out of bounds low'
        myheating(k)=total(1)
        dointerp=.FALSE.
       ELSE
        pinterpindexstart=pindex-1
       ENDIF

       IF(dointerp) THEN
        DO k1=1,4
         tempsub(k1)=total(pinterpindexstart+k1-1)
         psub(k1)=p_heating(pinterpindexstart+k1-1)
        ENDDO

        CALL POLINT(psub,tempsub,4,rCmb(k),tempout,dy)
        myheating(k)=tempout
       ENDIF
C       write(*,*) 'rCmb(k)=',rCmb(k),'myheating(k)=',myheating(k)
      ENDDO

      RETURN
      END
