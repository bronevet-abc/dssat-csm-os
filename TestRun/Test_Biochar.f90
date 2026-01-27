      PROGRAM Test_Biochar
      USE ModuleDefs
      USE Biochar_mod
      IMPLICIT NONE

      TYPE(ControlType) :: CONTROL
      TYPE(SoilType)    :: SOILPROP
      INTEGER :: YR, DOY
      INTEGER :: I, L
      
      ! Environmental Variables (Mocks)
      REAL, DIMENSION(NL) :: SW, ST, NH4, NO3

      ! Initialize Control
      CONTROL % DYNAMIC = RUNINIT
      CONTROL % YRDOY = 2002001
      CONTROL % DAS = 0
      
      ! Initialize Soil Props
      SOILPROP % NLAYR = 3
      SOILPROP % DLAYR = 0.0
      SOILPROP % LL    = 0.10
      SOILPROP % DUL   = 0.25
      SOILPROP % SAT   = 0.40
      
      SOILPROP % DLAYR(1) = 15.0
      SOILPROP % DLAYR(2) = 15.0
      SOILPROP % DLAYR(3) = 30.0
      
      ! Initialize Environmental State
      SW = 0.20   ! Moist
      ST = 25.0   ! 25 C
      NH4 = 5.0   ! kg/ha
      NO3 = 15.0  ! kg/ha
      
      PRINT *, "Initializing Biochar Module..."
      CALL Biochar_Init(CONTROL)

      ! Simulation Loop
      CONTROL % DYNAMIC = INTEGR
      DO I = 1, 365
         CONTROL % YRDOY = 2002000 + I
         CALL YR_DOY(CONTROL % YRDOY, YR, DOY)
         CONTROL % DAS = I
         
         ! Daily Call
         CALL Biochar_Daily(CONTROL, SOILPROP, SW, ST, NH4, NO3)
         
         ! Print status every 50 days
         IF (MOD(I, 50) == 0) THEN
            PRINT *, "Day", I, "processed."
         END IF
      END DO

      PRINT *, "Test Complete. Check BIOCHAR.OUT."

      END PROGRAM Test_Biochar
