#include "CPP_OPTIONS.h"
       subroutine settozero(nsize,array)

C     GCM v23   2010
C     Ames Mars GCM group
C     Jeffery Hollingsworth, PI
C     NASA Ames Research Center

       implicit none

#include "EEPARAMS.h"

       integer nsize
       _RL array(nsize)

       integer i

       do i=1,nsize
         array(i) = 0.0D0
       enddo

       return

       end
       subroutine settozero4(nsize,array)

       implicit none

#include "EEPARAMS.h"

       integer nsize
       _RS array(nsize)

       integer i

       do i=1,nsize
         array(i) = 0.0
       enddo

       return

       end
