#include "CPP_OPTIONS.h"
      subroutine dustprofile

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

C Bob's updates 9/17/99 
C Reference the dust optical depth to RPTAU, and modify the way
C the dust mixing ratio is calculated to more accurately reflect the
C pressure-optical depth relationship.
C  Called by EMISS and GMP.

      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "RT_GRID.h"
#include "DEFINES.h"
#include "CONSTANTS.h"
#include "STANDARD.h"

      _RL  QRDST(L_NPDST)

C     implicit none

      integer i, j, n
      _RL  refpr, conr, sum, pave, qrdst0, pdif1, pdif2
      _RL  deltals, conrtmp, tetalat

C======================================================================C

C Calculate the Reference Pressure Grid (prdst)

      refpr    = (5.0*psf/ptrop)**(1.0/(float(npdst)-1.0))
      prdst(1) = ptrop
   
      do n=2,npdst
        prdst(n) = refpr*prdst(n-1)
      end do

C Calculate the Mixing Ratio at the Reference Pressure Level

C  GCM1.7  6/28/01   spatially varying dust

      conr    = conrnu

      sum = 0.
      do n = 2,npdst
        if (prdst(n).lt.rptau) then
          pave = 0.5*(prdst(n)+prdst(n-1))
          sum  = sum + 
     &           exp(conr*(1.-(rptau/pave)))*
     &              (prdst(n)-prdst(n-1))
        end if
        if (prdst(n).ge.rptau) go to 10 
      end do

10    continue

      pave = 0.5*(rptau+prdst(n-1))
      sum  = sum + exp(conr*(1.-(rptau/pave)))*(rptau-prdst(n-1))

      qrdst0 = tautot/sum

C Now calculate the mixing ratio at all other levels

      do n=1,npdst-1

C Region 1: Mixing ratio changes continuously through the layer

        if (rptau.gt.prdst(n+1)) then
          pave     = 0.5*(prdst(n+1)+prdst(n))
          qrdst(n) = qrdst0*exp(conr*(1.0-(rptau/pave)))
        end if

C Region 2: Reference pressure level within this layer. 

        if (rptau.le.prdst(n+1).and.rptau.ge.prdst(n)) then
          pave     = 0.5*(prdst(n)+rptau)
          pdif1    = rptau-prdst(n)
          pdif2    = prdst(n+1)-rptau
          qrdst(n) = qrdst0*(
     &               exp(conr*(1.0-(rptau/pave)))*pdif1+pdif2) / 
     &               (prdst(n+1)-prdst(n))
        end if

C Region 3: Mixing ratio constant

        if (rptau.lt.prdst(n)) then
          qrdst(n) = qrdst0
        end if

      end do

C Now compute the optical depths (taudst).

      taudst(1) = 0.0

      do n=2,npdst
        taudst(n) = taudst(n-1) +
     &                   qrdst(n-1)*(prdst(n)-prdst(n-1))
      end do

      RETURN
      END
