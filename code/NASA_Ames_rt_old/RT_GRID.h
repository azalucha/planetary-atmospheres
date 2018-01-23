C NASA Ames RT code 
C module grid_h -> RT_GRID.h because already MIT GCM
C grid.h
C #include RT_GRID.h

      INTEGER L_LAYERS,L_LEVELS,NDP
C     Number of atmospheric layers
      PARAMETER(L_LAYERS=29)
C     Number of atmospheric levels:   2 * L_LAYERS + 3
      PARAMETER(L_LEVELS=2*L_LAYERS+3)
C     Number of dust particle sizes
      PARAMETER(NDP=2)

