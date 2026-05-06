PROGRAM Test_Biochar_Mechanistic
  USE Biochar_mod
  USE ModuleDefs
  IMPLICIT NONE

  TYPE (ControlType) :: CONTROL
  TYPE (SoilType) :: SOILPROP
  REAL, DIMENSION(NL) :: NH4, NO3, SPi_AVAIL
  INTEGER :: L

  PRINT *, "Starting Test_Biochar_Mechanistic..."

  ! Setup mock data
  CONTROL % YRDOY = 2006001
  CONTROL % DYNAMIC = INTEGR ! Or whatever is needed
  
  SOILPROP % NLAYR = 1
  SOILPROP % PH(1) = 6.5
  SOILPROP % CEC(1) = 10.0
  
  NH4(1) = 50.0
  NO3(1) = 20.0
  SPi_AVAIL(1) = 30.0

  ! Initialize
  CALL Biochar_Init(CONTROL)

  PRINT *, "Initial NH4:", NH4(1)
  PRINT *, "Initial NO3:", NO3(1)
  PRINT *, "Initial SPi_AVAIL:", SPi_AVAIL(1)

  ! Run Daily
  CALL Biochar_Daily(CONTROL, SOILPROP, NH4, NO3, SPi_AVAIL)

  PRINT *, "Final NH4:", NH4(1)
  PRINT *, "Final NO3:", NO3(1)
  PRINT *, "Final SPi_AVAIL:", SPi_AVAIL(1)

  ! Check if they changed
  IF (NH4(1) /= 50.0 .OR. NO3(1) /= 20.0 .OR. SPi_AVAIL(1) /= 30.0) THEN
     PRINT *, "SUCCESS: Pools changed."
  ELSE
     PRINT *, "FAILURE: Pools did not change."
  ENDIF

END PROGRAM Test_Biochar_Mechanistic
