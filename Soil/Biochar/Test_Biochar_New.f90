      PROGRAM Test_Biochar_New
      USE Biochar_mod
      USE ModuleDefs
      IMPLICIT NONE

      TYPE(ControlType) :: CONTROL
      TYPE(SoilType)    :: SOILPROP
      TYPE(SwitchType)  :: ISWITCH
      
      REAL, DIMENSION(NL) :: SW, ST, NH4, NO3
      REAL, DIMENSION(0:NL, NELEM) :: IMM, MNR
      
      INTEGER :: L, Day
      REAL :: InitialBD, FinalBD
      
      PRINT *, 'Initializing Test_Biochar_New...'
      
      ! Initialize Control
      CONTROL % YRDOY = 2002001
      CONTROL % DAS   = 0
      CONTROL % DYNAMIC = RUNINIT
      
      ! Initialize Soil Properties (Stub)
      SOILPROP % NLAYR = 3
      SOILPROP % DLAYR(1) = 20.0
      SOILPROP % DLAYR(2) = 20.0
      SOILPROP % DLAYR(3) = 20.0
      SOILPROP % BD(1) = 1.4
      SOILPROP % BD(2) = 1.5
      SOILPROP % BD(3) = 1.6
      SOILPROP % SAND(1) = 50.0
      SOILPROP % CLAY(1) = 20.0
      SOILPROP % OC(1)   = 1.0
      SOILPROP % PH(1)   = 6.5
      SOILPROP % CEC(1)  = 10.0
      
      ! Derived
      SOILPROP % LL(1:3)  = 0.10
      SOILPROP % DUL(1:3) = 0.20
      SOILPROP % SAT(1:3) = 0.40
      
      ! Native CEC
      SOILPROP % CEC(1:3) = 10.0
      
      ! Initialize Arrays
      SW = 0.25 ! WET
      ST = 25.0 ! Warm
      NH4 = 5.0 ! kg/ha? No ppm usually in inputs but Biochar_Daily treats as kg/ha arg?
      ! Biochar_Daily NH4 arg: REAL, DIMENSION(NL), INTENT(IN)
      ! In SOIL.for: CALL Biochar_Daily(..., NH4, ...)
      ! NH4 in SOIL.for is local variable (mg N/kg soil = ppm).
      ! Biochar_Daily Logic: NH4_Conc_mgL = (NH4(L) * 1.0E6) / VolSW_L_Ha
      ! If NH4 is ppm (mg/kg).
      ! Conc (mg/L) = (mg/kg * kg_soil) / L_water
      ! My logic in Biochar_Daily:
      ! NH4_Conc_mgL = (NH4(L) * 1.0E6) / VolSW_L_Ha
      ! This assumes NH4(L) is kg/ha?
      ! 1 kg/ha = 10^6 mg/ha.
      ! (mg/ha) / (L/ha) = mg/L.
      ! So IF NH4 is kg/ha, my logic is correct.
      ! BUT SOIL.for passes NH4 (ppm).
      ! I need to FIX Biochar_Daily if it expects kg/ha but gets ppm.
      ! Wait, let's check SOIL.for.
      ! In SOIL.for: `NH4` comes from `SoilNi(..., NH4, ...)`?
      ! `CALL SoilNi(..., NH4, ...)`
      ! `NH4` output from `SoilNi` is usually kg/ha?
      ! No, normally ppm.
      ! I should verify units in `SoilNi` or `SOIL.for`.
      ! For now, let's assume ppm and test.
      NH4 = 5.0 
      NO3 = 10.0
      IMM = 0.0
      MNR = 0.0
      
      ! 1. Initialize
      CALL Biochar_Init(CONTROL)
      
      ! Check initial state
      PRINT *, 'Initial BD (L1):', SOILPROP%BD(1)
      InitialBD = SOILPROP%BD(1)
      
      ! 2. Call Biochar_UpdateSoilProps (Should modify BD if AppDate passed or if it processes apps regardless of date?)
      ! My implementation of Biochar_UpdateSoilProps checks `NumApps > 0`.
      ! It does NOT check `AppDate <= YRDOY`. It uses ALL apps?
      ! "IF (NumApps == 0 ...)"
      ! It assumes biochar changes properties immediately?
      ! Or should it wait for application?
      ! The call is daily.
      ! "Refined Logic: ... effectively modifying ... based on biochar presence."
      ! "Biochar_UpdateSoilProps" uses `BC_Labile` and `BC_Recalc`.
      ! At Init, these are 0.
      
      CALL Biochar_UpdateSoilProps(SOILPROP)
      PRINT *, 'BD after Init Update:', SOILPROP%BD(1)
      
      ! 3. Run to Application Day (Day 10)
      DO Day = 1, 15
         CONTROL % DAS = Day
         CONTROL % YRDOY = 2002000 + Day
         
         ! Daily Call
         CALL Biochar_Daily(CONTROL, SOILPROP, SW, ST, NH4, NO3, IMM, MNR)
         
         ! Update Properties (After Biochar_Daily adds mass, UpdateSoilProps uses it next step or same step?)
         ! In SOIL.for order: UpdateSoilProps -> WATBAL -> Biochar_Daily.
         ! So on Day 10:
         ! 1. UpdateSoilProps (Mass is 0) -> No Change.
         ! 2. WATBAL (uses old props).
         ! 3. Biochar_Daily (Applies mass: Mass becomes > 0).
         ! Day 11:
         ! 1. UpdateSoilProps (Mass is > 0) -> Updates SOILPROP.
         ! 2. WATBAL (uses NEW props).
         
         ! In Test:
         CALL Biochar_UpdateSoilProps(SOILPROP) ! Call again to simulate next day logic or same day sequence?
         ! Since we are looping, effectively we simulate daily cycle.
         ! Proper order in loop:
         ! UpdateSoilProps
         ! (WATBAL skipped)
         ! Biochar_Daily
         
         IF (Day == 10 .OR. Day == 11) THEN
             PRINT *, 'Day', Day, 'BD:', SOILPROP%BD(1), 'Labile:', BC_Labile(1)
         END IF
      END DO
      
      FinalBD = SOILPROP%BD(1)
      PRINT *, 'Final BD:', FinalBD
      
      IF (FinalBD < InitialBD) THEN
         PRINT *, 'SUCCESS: Bulk Density decreased.'
      ELSE
         PRINT *, 'FAILURE: Bulk Density did not decrease.'
      END IF

      END PROGRAM Test_Biochar_New
