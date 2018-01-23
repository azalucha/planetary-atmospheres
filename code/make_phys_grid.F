       subroutine make_phys_grid(drF,hfacC,im1,im2,jm1,jm2,Nr,
     . Nsx,Nsy,i1,i2,j1,j2,bi,bj,Nrphys,Lbot,dpphys,numlevphys,nlperdyn)
c***********************************************************************
c subroutine make_phys_grid
c 
c Purpose: Define the grid that the will be used to run the high-end
c          atmospheric physics.
c
c Algorithm: FAKE
c
c Input:   drF         - delta r (p*) edge-to-edge
c          hfacC       - fraction of grid box above topography
c          im1, im2    - beginning and ending i - dimensions
c          jm1, jm2    - beginning and ending j - dimensions
c          Nr          - number of levels in dynamics grid
c          Nsx,Nsy     - number of processes in x and y direction
c          i1, i2      - beginning and ending i - index to fill
c          j1, j2      - beginning and ending j - index to fill
c          bi, bj      - x-dir and y-dir index of process
c          Nrphys      - number of levels in physics grid
c
c Output:  dpphys      - delta r (p*) edge-to-edge of physics grid
c          numlevphys  - number of levels used in the physics
c          nlperdyn    - physics level number atop each dynamics layer
c
c***********************************************************************
       implicit none
c
#include "CPP_OPTIONS.h"
c
       integer im1,im2,jm1,jm2,Nr,Nsx,Nsy,Nrphys
       integer i1,i2,j1,j2,bi,bj
       integer numlevphys
       _RL hfacC(im1:im2,jm1:jm2,Nr,Nsx,Nsy)
       _RL dpphys(im1:im2,jm1:jm2,Nrphys,Nsx,Nsy)
       _RL drF(Nr)
       integer Lbot(im1:im2,jm1:jm2,Nsx,Nsy)
       integer nlperdyn(im1:im2,jm1:jm2,Nr,Nsx,Nsy)
c
       integer i, j, L, Lbotij, Lnew, nlevs, nlphys, ndone
c
       do L = 1,Nr
        do j = j1,j2
        do i = i1,i2+1
         nlperdyn(i,j,L,bi,bj) = 0
        enddo
        enddo
       enddo
       numlevphys = Nrphys
c
c      do for each grid point:
       do j = j1,j2
       do i = i1,i2
        nlphys = 0
        Lbotij = Lbot(i,j,bi,bj)
        do L = Lbotij,Nr
         nlphys = nlphys + 1
         dpphys(i,j,nlphys,bi,bj) = hfacC(i,j,L,bi,bj)*drF(L)
         nlperdyn(i,j,L,bi,bj) = nlphys
        enddo
c
c  All done adding layers - if we need more to make numlevphys, put
c     them as thin (1 mb) layers near the top
        if(nlphys.lt.numlevphys)then
         nlevs = numlevphys-nlphys
c first subtract 100 pa * number of thin levels from dp of top level so far
         dpphys(i,j,nlphys,bi,bj) = dpphys(i,j,nlphys,bi,bj)-100.*nlevs
c now make the new thin levels
         do Lnew = nlphys+1,numlevphys
          dpphys(i,j,Lnew,bi,bj) = 100.
         enddo
         nlperdyn(i,j,Nr,bi,bj) = numlevphys
        endif
c END OF LOOP OVER GRID POINTS

       enddo
       enddo

       return
       end
