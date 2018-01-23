C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C                            defines_h
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C     L_LAYERS + 1
      INTEGER L_LAYERP1
      PARAMETER( L_LAYERP1 = L_LAYERS+1)

C     L_LAYERS - 1
      INTEGER  L_LAYERM1
      PARAMETER(L_LAYERM1 = L_LAYERS-1)

C     L_LEVELS - 1
      INTEGER L_LEVELM1
      PARAMETER( L_LEVELM1 = L_LEVELS-1)

C     L_LEVELS - 2
      INTEGER L_LEVELM2
      PARAMETER(L_LEVELM2 = L_LEVELS-2)

C     L_LEVELS - 3
      INTEGER L_LEVELM3
      PARAMETER(L_LEVELM3 = L_LEVELS-3)

C     L_LEVELS - 4
      INTEGER L_LEVELM4
      PARAMETER(L_LEVELM4 = L_LEVELS-4)

C     Number of fine-mesh layers used in calculating dust profiles
      INTEGER L_NPDST
      PARAMETER(L_NPDST   = 100)

