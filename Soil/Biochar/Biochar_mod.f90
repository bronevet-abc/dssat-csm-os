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
  !=======================================================================
  SUBROUTINE Biochar_Derivatives(t, y, dydt, SOILPROP, L)
    USE ModuleDefs
    IMPLICIT NONE
    REAL, INTENT(IN) :: t
    REAL, DIMENSION(8), INTENT(IN) :: y
    REAL, DIMENSION(8), INTENT(OUT) :: dydt
    TYPE (SoilType), INTENT(IN) :: SOILPROP
    INTEGER, INTENT(IN) :: L

    REAL :: BC_L, BC_R, M_BC, CEC_pot, Tot_NH4, Tot_NO3, Tot_PO4, P_precip
    REAL :: pH_dynamic, CEC_eff, fraction_NH3, frac_sorbed_NH4
    REAL :: NH4_aq, NH3_aq, AEC_bridge, comp_denom
    REAL :: PO4_aq, NO3_aq
    REAL :: dC_pot, N_req_max, N_req_min, N_avail, N_index, CUE_dyn
    REAL :: R_decomp, R_immob_N, eta_N, R_ox, R_nitrif, H_load_nitric
    REAL :: H_load_total, R_weather, R_volat, R_precip_P
    REAL :: immob_NH4, immob_NO3
    REAL :: dBC_L, dM_BC, dCEC_pot, dTot_NH4, dTot_NO3, dTot_PO4, dP_precip

    ! Unpack state
    BC_L = MAX(0.0, y(1))
    BC_R = MAX(0.0, y(2))
    M_BC = MAX(0.0, y(3))
    CEC_pot = MAX(0.0, y(4))
    Tot_NH4 = MAX(0.0, y(5))
    Tot_NO3 = MAX(0.0, y(6))
    Tot_PO4 = MAX(0.0, y(7))
    P_precip = MAX(0.0, y(8))

    ! 1. pH and CEC Partitioning
    pH_dynamic = SOILPROP%PH(L) + (M_BC * 0.02) 
    
    ! Henderson-Hasselbalch for Effective CEC
    CEC_eff = CEC_pot * (1.0 / (1.0 + 10.0**(pKa_carboxyl - pH_dynamic)))
    
    ! Henderson-Hasselbalch for Ammonia Volatilization potential
    fraction_NH3 = 1.0 / (1.0 + 10.0**(pKa_NH4 - pH_dynamic))
    
    ! Cation exchange for NH4
    frac_sorbed_NH4 = MIN(0.9, (CEC_eff * 0.1) / (Tot_NH4 + 1.0e-6))
    NH4_aq   = Tot_NH4 - (Tot_NH4 * frac_sorbed_NH4)
    NH3_aq   = NH4_aq * fraction_NH3
    
    ! 2. Cation Bridging & Anion Competition
    AEC_bridge = CEC_eff * f_poly * k_bridge_eff
    comp_denom = 1.0 + (K_L_NO3 * Tot_NO3) + (K_L_PO4 * Tot_PO4)
    
    NO3_aq = Tot_NO3 - MIN(AEC_bridge * ((K_L_NO3 * Tot_NO3) / comp_denom), Tot_NO3)
    PO4_aq = Tot_PO4 - MIN(AEC_bridge * ((K_L_PO4 * Tot_PO4) / comp_denom), Tot_PO4)
    
    ! 3. Kinetic Rates
    
    ! Carbon Degradation
    dC_pot = BC_L * k_labile
    N_req_max = dC_pot * (CUE_max / CN_mic)
    N_req_min = dC_pot * (CUE_min / CN_mic)
    
    N_avail = NH4_aq + NO3_aq
    N_index = MIN(1.0, N_avail / (N_req_max + 1.0e-6))
    CUE_dyn = CUE_min + (CUE_max - CUE_min) * N_index
    
    eta_N = MIN(1.0, N_avail / (N_req_min + 1.0e-6))
    R_decomp = dC_pot * eta_N
    R_immob_N = R_decomp * (CUE_dyn / CN_mic)
    
    ! Surface Oxidation
    R_ox = k_ox_max * ((CEC_max * (BC_R/10000.0)) - CEC_pot)
    
    ! Nitrification & Weathering
    R_nitrif = k_nitrif * NH4_aq
    H_load_nitric = R_nitrif * 2.0
    H_load_total  = H_load_nitric + 0.01 ! background
    
    R_weather = k_weather * M_BC * (H_load_total**0.8)
    
    ! Volatilization & Precipitation
    R_volat = k_volat * NH3_aq
    R_precip_P = 0.05 * R_weather * PO4_aq
    
    ! 4. Derivatives
    immob_NH4 = R_immob_N * 0.8
    immob_NO3 = R_immob_N * 0.2
    
    dBC_L = -R_decomp
    dM_BC = -R_weather
    dCEC_pot = R_ox
    
    dTot_NH4 = -R_nitrif - R_volat - immob_NH4
    dTot_NO3 =  R_nitrif - immob_NO3
    dTot_PO4 = -R_precip_P
    dP_precip = R_precip_P

    ! Pack derivatives
    dydt(1) = dBC_L
    dydt(2) = 0.0 ! BC_R does not change
    dydt(3) = dM_BC
    dydt(4) = dCEC_pot
    dydt(5) = dTot_NH4
    dydt(6) = dTot_NO3
    dydt(7) = dTot_PO4
    dydt(8) = dP_precip

  END SUBROUTINE Biochar_Derivatives

!=======================================================================
  SUBROUTINE RK45_Integrate(y, t_start, t_end, tol, SOILPROP, L)
    USE ModuleDefs
    IMPLICIT NONE
    REAL, DIMENSION(8), INTENT(INOUT) :: y
    REAL, INTENT(IN) :: t_start, t_end, tol
    TYPE (SoilType), INTENT(IN) :: SOILPROP
    INTEGER, INTENT(IN) :: L

    REAL :: t, dt, scale
    REAL, DIMENSION(8) :: y_new, y_temp, error
    REAL, DIMENSION(8) :: k1, k2, k3, k4, k5, k6
    REAL :: err_max

    ! Fehlberg coefficients
    REAL, PARAMETER :: b21 = 0.25
    REAL, PARAMETER :: b31 = 3.0/32.0, b32 = 9.0/32.0
    REAL, PARAMETER :: b41 = 1932.0/2197.0, b42 = -7200.0/2197.0, b43 = 7296.0/2197.0
    REAL, PARAMETER :: b51 = 439.0/216.0, b52 = -8.0, b53 = 3680.0/513.0, b54 = -845.0/4104.0
    REAL, PARAMETER :: b61 = -8.0/27.0, b62 = 2.0, b63 = -3544.0/2565.0, b64 = 1859.0/4104.0, b65 = -11.0/40.0

    REAL, PARAMETER :: c1 = 25.0/216.0, c3 = 1408.0/2565.0, c4 = 2197.0/4104.0, c5 = -0.15
    REAL, PARAMETER :: c2 = 0.0, c6 = 0.0 ! for 4th order

    REAL, PARAMETER :: chat1 = 16.0/135.0, chat3 = 6656.0/12825.0, chat4 = 28561.0/56430.0, chat5 = -9.0/50.0, chat6 = 2.0/55.0
    REAL, PARAMETER :: chat2 = 0.0 ! for 5th order

    t = t_start
    dt = 1.0 ! Initial guess: 1 hour

    DO WHILE (t < t_end)
       IF (t + dt > t_end) dt = t_end - t

       CALL Biochar_Derivatives(t, y, k1, SOILPROP, L)
       
       y_temp = y + dt * b21 * k1
       CALL Biochar_Derivatives(t + 0.25*dt, y_temp, k2, SOILPROP, L)
       
       y_temp = y + dt * (b31*k1 + b32*k2)
       CALL Biochar_Derivatives(t + 0.375*dt, y_temp, k3, SOILPROP, L)
       
       y_temp = y + dt * (b41*k1 + b42*k2 + b43*k3)
       CALL Biochar_Derivatives(t + (12.0/13.0)*dt, y_temp, k4, SOILPROP, L)
       
       y_temp = y + dt * (b51*k1 + b52*k2 + b53*k3 + b54*k4)
       CALL Biochar_Derivatives(t + dt, y_temp, k5, SOILPROP, L)
       
       y_temp = y + dt * (b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
       CALL Biochar_Derivatives(t + 0.5*dt, y_temp, k6, SOILPROP, L)
       
       ! 4th order estimate
       y_new = y + dt * (c1*k1 + c3*k3 + c4*k4 + c5*k5)
       
       ! 5th order estimate
       y_temp = y + dt * (chat1*k1 + chat3*k3 + chat4*k4 + chat5*k5 + chat6*k6)
       
       ! Error estimate
       error = abs(y_temp - y_new)
       err_max = maxval(error)
       
       ! Scale for step size adjustment
       scale = 0.84 * (tol / (err_max + 1.0e-10))**0.25
       scale = min(max(scale, 0.1), 4.0) ! Limit change
       
       IF (err_max <= tol) THEN
          y = y_temp ! Accept 5th order result
          t = t + dt
          dt = dt * scale
       ELSE
          dt = dt * scale ! Reject and reduce dt
       ENDIF
       
    END DO

  END SUBROUTINE RK45_Integrate

!=======================================================================
  SUBROUTINE Biochar_Daily(CONTROL, SOILPROP, NH4, NO3, SPi_AVAIL)
    USE ModuleDefs
    IMPLICIT NONE

    TYPE (ControlType), INTENT(IN) :: CONTROL
    TYPE (SoilType), INTENT(INOUT) :: SOILPROP
    REAL, DIMENSION(NL), INTENT(INOUT) :: NH4, NO3, SPi_AVAIL

    INTEGER :: L
    REAL, DIMENSION(8) :: y
    REAL :: tol

    IF (.NOT. Initialized) CALL Biochar_Init(CONTROL)

    tol = 1.0e-5

    DO L = 1, SOILPROP%NLAYR
       
       ! Pack state
       y(1) = BCState%BC_L(L)
       y(2) = BCState%BC_R(L)
       y(3) = BCState%M_BC(L)
       y(4) = BCState%CEC_pot(L)
       y(5) = NH4(L)
       y(6) = NO3(L)
       y(7) = SPi_AVAIL(L)
       y(8) = BCState%P_precip(L)

       ! Integrate from 0 to 24 hours
       CALL RK45_Integrate(y, 0.0, 24.0, tol, SOILPROP, L)
       
       ! Unpack state
       BCState%BC_L(L) = y(1)
       BCState%BC_R(L) = y(2)
       BCState%M_BC(L) = y(3)
       BCState%CEC_pot(L) = y(4)
       NH4(L) = y(5)
       NO3(L) = y(6)
       SPi_AVAIL(L) = y(7)
       BCState%P_precip(L) = y(8)
       
       ! Update soil properties (pH and CEC)
       SOILPROP%PH(L) = SOILPROP%PH(L) + (BCState%M_BC(L) * 0.02) ! simplified
       SOILPROP%CEC(L) = SOILPROP%CEC(L) + BCState%CEC_pot(L) * (1.0 / (1.0 + 10.0**(pKa_carboxyl - SOILPROP%PH(L)))) ! simplified
       
    END DO

  END SUBROUTINE Biochar_Daily

END MODULE Biochar_mod
