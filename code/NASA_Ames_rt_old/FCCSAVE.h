C NASA Ames RT code
C module fccsave_h
C #include FCCSAVE.h

C  Cloud Radiative properties

      _RL  TAUREFCLD(L_LEVELS+1)
      _RL  QEXTREFCLD(L_LEVELS+1)

      _RL  QXVCLD(L_LEVELS+1,L_NSPECTV)
      _RL  QSVCLD(L_LEVELS+1,L_NSPECTV)
      _RL  GVCLD(L_LEVELS+1,L_NSPECTV)

      _RL  QXICLD(L_LEVELS+1,L_NSPECTI)
      _RL  QSICLD(L_LEVELS+1,L_NSPECTI)
      _RL  GICLD(L_LEVELS+1,L_NSPECTI)

c  Dust Radiative properties

      _RL  QEXTREFDST(L_LEVELS+1)

      _RL  QXVDST(L_LEVELS+1,L_NSPECTV)
      _RL  QSVDST(L_LEVELS+1,L_NSPECTV)
      _RL  GVDST(L_LEVELS+1,L_NSPECTV)

      _RL  QXIDST(L_LEVELS+1,L_NSPECTI)
      _RL  QSIDST(L_LEVELS+1,L_NSPECTI)
      _RL  GIDST(L_LEVELS+1,L_NSPECTI)

      common /fccsave/ TAUREFCLD,QEXTREFCLD,QXVCLD,QSVCLD,
     &      GVCLD,QXICLD,QSICLD,GICLD,QEXTREFDST,QXVDST,
     &      QSVDST,GVDST,QXIDST,QSIDST,GIDST
