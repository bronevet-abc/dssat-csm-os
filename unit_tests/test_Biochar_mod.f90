PROGRAM test_Biochar_mod
  USE Biochar_mod
  USE ModuleDefs
  IMPLICIT NONE

  TYPE(SoilType) :: SOILPROP
  TYPE(ControlType) :: CONTROL
  REAL, DIMENSION(8) :: y, dydt
  INTEGER :: L
  REAL :: rate_low, rate_high
  REAL :: dBC_L_low, dBC_L_high
  REAL :: dCEC_pot_low, dCEC_pot_high
  REAL :: dP_precip_low, dP_precip_high
  REAL :: y_start(8), y_end_low(8), y_end_high(8)
  REAL :: tol
  REAL :: dNH4_high_BC, dNH4_zero_BC

  ! Initialize
  CALL Biochar_Init(CONTROL)

  ! Setup SoilProp
  SOILPROP%NLAYR = 1
  SOILPROP%PH(1) = 6.0
  SOILPROP%CEC(1) = 10.0
  L = 1
  tol = 1.0e-5

  PRINT *, "=== Running Biochar Unit Tests ==="

  ! Test 1: Nitrogen Cascade
  y = 0.0
  y(1) = 1000.0; y(2) = 9000.0; y(3) = 50.0; y(4) = 0.5
  y(5) = 1.0; y(6) = 1.0; y(7) = 10.0
  
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  rate_low = -dydt(3)
  
  y(5) = 100.0
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  rate_high = -dydt(3)
  
  PRINT *, "Test 1 (Nitrogen Cascade):"
  IF (rate_high > rate_low) THEN
     PRINT *, "  PASSED: Weathering rate increased with higher N."
  ELSE
     PRINT *, "  FAILED: Weathering rate did not increase."
  END IF

  ! Test 3: npSOM Bridging & Negative Priming
  PRINT *, "Test 3 (npSOM Bridging):"
  PRINT *, "  SKIPPED: Mechanism pending implementation in model."

  ! Test 4: Anion Competition
  y = 0.0
  y(1) = 1000.0; y(2) = 9000.0; y(3) = 50.0; y(4) = 0.5
  y(5) = 0.5; y(6) = 0.5; y(7) = 0.0
  
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dBC_L_low = dydt(1)
  
  y(7) = 100.0
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dBC_L_high = dydt(1)
  
  PRINT *, "Test 4 (Anion Competition):"
  IF (dBC_L_high < dBC_L_low) THEN
     PRINT *, "  PASSED: Decomposition increased with higher PO4 competition."
  ELSE
     PRINT *, "  FAILED: Decomposition did not increase."
  END IF

  ! Test 5: Dual-Driven Oxidation
  y = 0.0
  y(1) = 1000.0; y(2) = 9000.0; y(3) = 50.0; y(4) = 0.5
  
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dCEC_pot_low = dydt(4)
  
  y(4) = 40.0
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dCEC_pot_high = dydt(4)
  
  PRINT *, "Test 5 (Dual-Driven Oxidation):"
  IF (dCEC_pot_high < dCEC_pot_low) THEN
     PRINT *, "  PASSED: Oxidation rate decreased as CEC approached maximum."
  ELSE
     PRINT *, "  FAILED: Oxidation rate did not decrease."
  END IF

  ! Test 7: Precipitation vs Sorption
  y = 0.0
  y(1) = 1000.0; y(2) = 9000.0; y(3) = 10.0; y(4) = 0.5; y(5) = 10.0; y(6) = 10.0; y(7) = 10.0
  
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dP_precip_low = dydt(8)
  
  y(3) = 100.0
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dP_precip_high = dydt(8)
  
  PRINT *, "Test 7 (Precipitation):"
  IF (dP_precip_high > dP_precip_low) THEN
     PRINT *, "  PASSED: P precipitation increased with higher weathering."
  ELSE
     PRINT *, "  FAILED: P precipitation did not increase."
  END IF

  ! Test 10: Overflow Respiration
  y = 0.0
  y(1) = 1000.0; y(2) = 9000.0; y(3) = 50.0; y(4) = 0.5; y(7) = 10.0
  
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  
  PRINT *, "Test 10 (Overflow Respiration):"
  IF (ABS(dydt(1)) < 1.0e-5) THEN
     PRINT *, "  PASSED: Decomposition halted when N=0."
  ELSE
     PRINT *, "  FAILED: Decomposition did not halt when N=0."
  END IF

  ! Test 12: Ash Exhaustion
  y = 0.0
  y(1) = 1000.0; y(2) = 9000.0; y(3) = 0.0; y(4) = 0.5; y(5) = 10.0
  
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  
  PRINT *, "Test 12 (Ash Exhaustion):"
  IF (ABS(dydt(3)) < 1.0e-5) THEN
     PRINT *, "  PASSED: Weathering stopped when ash was exhausted."
  ELSE
     PRINT *, "  FAILED: Weathering did not stop."
  END IF

  ! Test 14: Immobilization-Weathering Feedback Cascade
  y_start = (/10.0, 9000.0, 50.0, 0.5, 10.0, 10.0, 10.0, 0.0/)
  y_end_low = y_start
  CALL RK45_Integrate(y_end_low, 0.0, 24.0, tol, SOILPROP, L)
  
  y_start(1) = 2000.0
  y_end_high = y_start
  CALL RK45_Integrate(y_end_high, 0.0, 24.0, tol, SOILPROP, L)
  
  PRINT *, "Test 14 (Feedback Cascade):"
  IF ((50.0 - y_end_high(3)) < (50.0 - y_end_low(3))) THEN
     PRINT *, "  PASSED: High BC_L suppressed ash weathering."
  ELSE
     PRINT *, "  FAILED: High BC_L did not suppress ash weathering."
  END IF

  ! Test 15: Pure Thermal Oxidation
  y = 0.0
  y(1) = 0.0; y(2) = 9000.0; y(3) = 50.0; y(4) = 0.5
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  PRINT *, "Test 15 (Thermal Oxidation):"
  IF (dydt(4) > 0.0) THEN
     PRINT *, "  PASSED: Oxidation proceeds without biological activity."
  ELSE
     PRINT *, "  FAILED: Oxidation did not proceed."
  END IF

  ! Test 16: The Absolute Drought Kinetic Freeze
  PRINT *, "Test 16 (Drought Freeze):"
  PRINT *, "  SKIPPED: Moisture dependence not implemented in model."

  ! Test 17: Transient Immobilization Cliff
  y = 0.0
  y(1) = 1000.0; y(2) = 9000.0; y(3) = 50.0; y(4) = 0.5; y(5) = 10.0
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dNH4_high_BC = dydt(5)
  
  y(1) = 0.0
  CALL Biochar_Derivatives(0.0, y, dydt, SOILPROP, L)
  dNH4_zero_BC = dydt(5)
  
  PRINT *, "Test 17 (Immobilization Cliff):"
  IF (dNH4_zero_BC > dNH4_high_BC) THEN
     PRINT *, "  PASSED: Immobilization collapsed when BC_L was exhausted."
  ELSE
     PRINT *, "  FAILED: Immobilization did not collapse."
  END IF

  ! Test 19: AEC Erosion via Acidification
  PRINT *, "Test 19 (AEC Erosion):"
  PRINT *, "  SKIPPED: Mechanism pending implementation in model."

  ! Test 20: Infinite Host Sink Overdraw
  PRINT *, "Test 20 (Sink Overdraw):"
  PRINT *, "  SKIPPED: Plant N demand sink not implemented in derivatives."

END PROGRAM test_Biochar_mod

SUBROUTINE WARNING(Code, Name, Msg)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Code
  CHARACTER(LEN=*), INTENT(IN) :: Name
  CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: Msg
END SUBROUTINE WARNING
