#include "CPP_OPTIONS.h"
      SUBROUTINE INTERP_T_TO_TL(pl,tl,L_LEVELS,myTheta)

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"

      INTEGER L_LEVELS,k,pindex,k1,pinterpindexstart
      _RL pl(L_LEVELS), tl(L_LEVELS),rCmb(Nr)
      _RL psub(4),tempsub(4),tempout,dy,myTheta(Nr)
      LOGICAL dointerp

      DO k=1,Nr
       rCmb(k)=rC(k)/100.
C       write(*,*) 'rCmb(k)=',rCmb(k),
C     &             'myT(k)=',myTheta(k)*(rC(k)/atm_po)**atm_kappa
      ENDDO

      DO k=2,L_LEVELS-1
       dointerp=.TRUE.
       CALL LOCATE(rCmb,Nr,pl(k),pindex)


       IF(pindex .EQ. 1) THEN
        pinterpindexstart=1
       ELSE IF (pindex .EQ. Nr-1) THEN
        pinterpindexstart=Nr-3
       ELSE IF (pindex .EQ. Nr) THEN
C        write(*,*) 'Locate: p out of bounds high'
        tl(k)=myTheta(Nr)*(rC(Nr)/atm_po)**atm_kappa        
        dointerp=.FALSE.
       ELSE IF (pindex .EQ. 0) THEN
C        write(*,*) 'Locate: p out of bounds low'
         tl(k)=myTheta(1)*(rC(1)/atm_po)**atm_kappa
         dointerp=.FALSE.
       ELSE
        pinterpindexstart=pindex-1
       ENDIF

       IF(dointerp) THEN
        DO k1=1,4
         tempsub(k1)=myTheta(pinterpindexstart+k1-1)
         psub(k1)=rCmb(pinterpindexstart+k1-1)
        ENDDO

       CALL POLINT(psub,tempsub,4,pl(k),tempout,dy)
       tl(k)=tempout*(pl(k)/atm_po*100.)**atm_kappa
       ENDIF
C       write(*,*) 'pl(k)=',pl(k),'tl(k)=',tl(k)
      ENDDO

      tl(1)=tl(2)

      RETURN
      END
