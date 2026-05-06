MODULE Biochar_mod
  USE ModuleDefs
  IMPLICIT NONE
  SAVE

  ! Biochar state variables per layer
  TYPE BiocharStateType
     REAL, DIMENSION(NL) :: BC_L     ! Labile Biochar C (kg/ha)
     REAL, DIMENSION(NL) :: BC_R     ! Recalcitrant Biochar C (kg/ha)
     REAL, DIMENSION(NL) :: M_BC     ! Unweathered Ash Base Cations (kmol_c/ha)
     REAL, DIMENSION(NL) :: CEC_pot  ! Potential CEC (cmol/kg or similar)
     REAL, DIMENSION(NL) :: P_precip ! Stable Ca-P precipitate (kg/ha)
     REAL, DIMENSION(NL) :: NH4_sorb ! Sorbed NH4 (kg/ha)
     REAL, DIMENSION(NL) :: NO3_sorb ! Sorbed NO3 (kg/ha)
     REAL, DIMENSION(NL) :: PO4_sorb ! Sorbed PO4 (kg/ha)
  END TYPE BiocharStateType

  TYPE(BiocharStateType) :: BCState

  LOGICAL :: Initialized = .FALSE.

  ! Parameters (could be read from file)
  REAL, PARAMETER :: k_labile = 0.05
  REAL, PARAMETER :: k_ox_max = 0.02
  REAL, PARAMETER :: k_weather = 0.1
  REAL, PARAMETER :: k_nitrif = 0.08
  REAL, PARAMETER :: k_volat = 0.15
  REAL, PARAMETER :: CUE_max = 0.6
  REAL, PARAMETER :: CUE_min = 0.2
  REAL, PARAMETER :: CN_mic = 8.0
  REAL, PARAMETER :: CEC_max = 50.0
  REAL, PARAMETER :: pKa_NH4 = 9.24
  REAL, PARAMETER :: pKa_carboxyl = 4.5
  REAL, PARAMETER :: f_poly = 0.6
  REAL, PARAMETER :: k_bridge_eff = 0.8
  REAL, PARAMETER :: K_L_NO3 = 0.2
  REAL, PARAMETER :: K_L_PO4 = 2.5

CONTAINS

!=======================================================================
  SUBROUTINE Biochar_Init(CONTROL)
    USE ModuleDefs
    IMPLICIT NONE
    TYPE (ControlType), INTENT(IN) :: CONTROL

    IF (Initialized) RETURN

    ! Initialize state variables (example values)
    BCState % BC_L = 1000.0
    BCState % BC_R = 9000.0
    BCState % M_BC = 50.0
    BCState % CEC_pot = 0.5
    BCState % P_precip = 0.0
    BCState % NH4_sorb = 0.0
    BCState % NO3_sorb = 0.0
    BCState % PO4_sorb = 0.0

    Initialized = .TRUE.
  END SUBROUTINE Biochar_Init

!=======================================================================
  SUBROUTINE Biochar_Daily(CONTROL, SOILPROP, NH4, NO3, SPi_AVAIL)
    USE ModuleDefs
    IMPLICIT NONE

    TYPE (ControlType), INTENT(IN) :: CONTROL
    TYPE (SoilType), INTENT(INOUT) :: SOILPROP
    REAL, DIMENSION(NL), INTENT(INOUT) :: NH4, NO3, SPi_AVAIL

    INTEGER :: L, Step, NumSteps
    REAL :: dt, t
    REAL :: pH_dynamic, CEC_eff, fraction_NH3, frac_sorbed_NH4
    REAL :: NH4_aq, NH3_aq, AEC_bridge, comp_denom
    REAL :: PO4_aq, NO3_aq
    REAL :: dC_pot, N_req_max, N_req_min, N_avail, N_index, CUE_dyn
    REAL :: R_decomp, R_immob_N, eta_N, R_ox, R_nitrif, H_load_nitric
    REAL :: H_load_total, R_weather, R_volat, R_precip_P
    REAL :: immob_NH4, immob_NO3
    REAL :: dBC_L, dM_BC, dCEC_pot, dTot_NH4, dTot_NO3, dTot_PO4, dP_precip

    REAL :: Tot_NH4, Tot_NO3, Tot_PO4

    IF (.NOT. Initialized) CALL Biochar_Init(CONTROL)

    NumSteps = 24
    dt = 1.0 / REAL(NumSteps) ! Daily step divided into sub-steps

    DO L = 1, SOILPROP%NLAYR
       
       ! Local copies of total pools for this layer
       Tot_NH4 = NH4(L)
       Tot_NO3 = NO3(L)
       Tot_PO4 = SPi_AVAIL(L)

       DO Step = 1, NumSteps
          
          ! 1. pH and CEC Partitioning
          pH_dynamic = SOILPROP%PH(L) + (BCState%M_BC(L) * 0.02) 
          
          ! Henderson-Hasselbalch for Effective CEC
          CEC_eff = BCState%CEC_pot(L) * (1.0 / (1.0 + 10.0**(pKa_carboxyl - pH_dynamic)))
          
          ! Henderson-Hasselbalch for Ammonia Volatilization potential
          fraction_NH3 = 1.0 / (1.0 + 10.0**(pKa_NH4 - pH_dynamic))
          
          ! Cation exchange for NH4
          frac_sorbed_NH4 = MIN(0.9, (CEC_eff * 0.1) / (Tot_NH4 + 1.0e-6))
          BCState%NH4_sorb(L) = Tot_NH4 * frac_sorbed_NH4
          NH4_aq   = Tot_NH4 - BCState%NH4_sorb(L)
          NH3_aq   = NH4_aq * fraction_NH3
          
          ! 2. Cation Bridging & Anion Competition
          AEC_bridge = CEC_eff * f_poly * k_bridge_eff
          comp_denom = 1.0 + (K_L_NO3 * Tot_NO3) + (K_L_PO4 * Tot_PO4)
          
          BCState%NO3_sorb(L) = AEC_bridge * ((K_L_NO3 * Tot_NO3) / comp_denom)
          BCState%PO4_sorb(L) = AEC_bridge * ((K_L_PO4 * Tot_PO4) / comp_denom)
          
          BCState%NO3_sorb(L) = MIN(BCState%NO3_sorb(L), Tot_NO3)
          BCState%PO4_sorb(L) = MIN(BCState%PO4_sorb(L), Tot_PO4)
          
          NO3_aq = Tot_NO3 - BCState%NO3_sorb(L)
          PO4_aq = Tot_PO4 - BCState%PO4_sorb(L)
          
          ! 3. Kinetic Rates
          
          ! Carbon Degradation
          dC_pot = BCState%BC_L(L) * k_labile
          N_req_max = dC_pot * (CUE_max / CN_mic)
          N_req_min = dC_pot * (CUE_min / CN_mic)
          
          N_avail = NH4_aq + NO3_aq
          N_index = MIN(1.0, N_avail / (N_req_max + 1.0e-6))
          CUE_dyn = CUE_min + (CUE_max - CUE_min) * N_index
          
          eta_N = MIN(1.0, N_avail / (N_req_min + 1.0e-6))
          R_decomp = dC_pot * eta_N
          R_immob_N = R_decomp * (CUE_dyn / CN_mic)
          
          ! Surface Oxidation
          R_ox = k_ox_max * ((CEC_max * (BCState%BC_R(L)/10000.0)) - BCState%CEC_pot(L))
          
          ! Nitrification & Weathering
          R_nitrif = k_nitrif * NH4_aq
          H_load_nitric = R_nitrif * 2.0
          H_load_total  = H_load_nitric + 0.01 ! background
          
          R_weather = k_weather * BCState%M_BC(L) * (H_load_total**0.8)
          
          ! Volatilization & Precipitation
          R_volat = k_volat * NH3_aq
          R_precip_P = 0.05 * R_weather * PO4_aq
          
          ! 4. Derivatives (Euler Step)
          immob_NH4 = R_immob_N * 0.8
          immob_NO3 = R_immob_N * 0.2
          
          dBC_L = -R_decomp
          dM_BC = -R_weather
          dCEC_pot = R_ox
          
          dTot_NH4 = -R_nitrif - R_volat - immob_NH4
          dTot_NO3 =  R_nitrif - immob_NO3
          dTot_PO4 = -R_precip_P
          dP_precip = R_precip_P
          
          ! Update state variables
          BCState%BC_L(L) = BCState%BC_L(L) + dBC_L * dt
          BCState%M_BC(L) = BCState%M_BC(L) + dM_BC * dt
          BCState%CEC_pot(L) = BCState%CEC_pot(L) + dCEC_pot * dt
          
          Tot_NH4 = Tot_NH4 + dTot_NH4 * dt
          Tot_NO3 = Tot_NO3 + dTot_NO3 * dt
          Tot_PO4 = Tot_PO4 + dTot_PO4 * dt
          BCState%P_precip(L) = BCState%P_precip(L) + dP_precip * dt
          
       END DO
       
       ! Return updated pools to DSSAT
       NH4(L) = Tot_NH4
       NO3(L) = Tot_NO3
       SPi_AVAIL(L) = Tot_PO4
       
       ! Update soil properties (pH and CEC)
       SOILPROP%PH(L) = SOILPROP%PH(L) + (BCState%M_BC(L) * 0.02) ! simplified
       SOILPROP%CEC(L) = SOILPROP%CEC(L) + CEC_eff ! simplified
       
    END DO

  END SUBROUTINE Biochar_Daily

END MODULE Biochar_mod
