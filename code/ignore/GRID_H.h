C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C                             grid_h
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C     Number of atmospheric layers
      INTEGER L_LAYERS
      PARAMETER(L_LAYERS=29)

C     Number of atmospheric levels:   2 * L_LAYERS + 3
      INTEGER L_LEVELS
      PARAMETER(L_LEVELS=2*L_LAYERS+3)

C     Number of dust particle sizes
      INTEGER NDP
      PARAMETER(NDP=2)

