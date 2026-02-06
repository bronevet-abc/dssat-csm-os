      MODULE Biochar_mod
!=======================================================================
!  MODULE Biochar_mod
!  Purpose: Simulate Biochar application and decay in soil.
!  Tracks biochar pools and partitions decay into CO2, Biomass, and Humic pools.
!  Calculates N Mineralization/Immobilization fluxes.
!=======================================================================
      USE ModuleDefs
      IMPLICIT NONE
      SAVE

!     Biochar Application Type
      TYPE BiocharAppType
        INTEGER :: AppDate       ! Date of application (YrDoy)
        REAL    :: Amount        ! Amount applied (kg/ha)
        REAL    :: Depth         ! Depth of application (cm)
        REAL    :: FLoss         ! Fraction lost during application (0-1)
        REAL    :: FCarbon       ! Fraction of Carbon (0-1)
        REAL    :: FLabile       ! Fraction of Labile pool (0-1)
        REAL    :: MRT_Labile    ! Mean Residence Time Labile (years)
        REAL    :: MRT_Recalc    ! Mean Residence Time Recalcitrant (years)
        REAL    :: CN_BC         ! C:N Ratio of Biochar
      END TYPE BiocharAppType

!     Max number of applications
      INTEGER, PARAMETER :: MaxApp = 20
      TYPE(BiocharAppType) :: BC_Apps(MaxApp)
      INTEGER :: NumApps = 0
      
!     Parameters
      REAL :: CNRF_BC = 0.693 ! C:N Ratio Factor (Archontoulis 2015)
      REAL :: Opt_bc  = 25.0  ! Optimal C:N 
      REAL :: EF_BC   = 0.4   ! Carbon retention efficiency (0-1)
      REAL :: FR_BCBIOM = 0.05 ! Fraction of retained C going to BIOM (0-1)
      
      ! Assumed C:N ratios for Soil Organic Matter pools (defaults)

      REAL :: CN_BIOM = 8.0
      REAL :: CN_HUM  = 11.0
      
!     Priming Parameters (Archontoulis et al., 2015)
      REAL :: P_FOM = 0.0   ! Positive priming on FOM decomposition rate
      REAL :: P_E   = 0.0   ! Negative priming on Carbon Efficiency
      REAL :: P_F   = 0.0   ! Negative priming on FOM->BIOM transfer

!     State Variables
      REAL, DIMENSION(NL) :: BC_Labile   ! Labile Biochar C (kg/ha)
      REAL, DIMENSION(NL) :: BC_Recalc   ! Recalcitrant Biochar C (kg/ha)
      
!     Flux Variables (Daily)
      REAL :: Daily_CO2_Gross   ! Total daily CO2 emission from Biochar (kg C/ha/d)
      REAL :: Daily_Biom_Gross  ! Total daily flux to Biomass (kg C/ha/d)
      REAL :: Daily_Hum_Gross   ! Total daily flux to Humic (kg C/ha/d)
      REAL :: Daily_N_Net       ! Net N mineralization (+)/immobilization (-) (kg N/ha/d)

!     Output file unit
      INTEGER :: LUN_BC
      LOGICAL :: FirstOutput = .TRUE.

      CONTAINS

!=======================================================================
      SUBROUTINE Biochar_Init(CONTROL)
        TYPE(ControlType), INTENT(IN) :: CONTROL
        INTEGER :: ERRNUM, LUN
        CHARACTER(LEN=120) :: LINE
        LOGICAL :: FEXIST
        
        ! Initialize State
        BC_Labile = 0.0
        BC_Recalc = 0.0
        NumApps = 0
        CNRF_BC = 0.693
        Opt_bc  = 25.0
        EF_BC   = 0.4
        FR_BCBIOM = 0.05
        CN_BIOM = 8.0
        CN_HUM  = 11.0
        P_FOM   = 0.0
        P_E     = 0.0
        P_F     = 0.0
        
        Daily_CO2_Gross = 0.0
        Daily_Biom_Gross = 0.0
        Daily_Hum_Gross = 0.0
        Daily_N_Net = 0.0

        ! Open and Read BIOCHAR.INP if it exists
        INQUIRE(FILE='BIOCHAR.INP', EXIST=FEXIST)
        IF (FEXIST) THEN
          OPEN(NEWUNIT=LUN, FILE='BIOCHAR.INP', STATUS='OLD', &
               ACTION='READ', IOSTAT=ERRNUM)
          IF (ERRNUM == 0) THEN
            DO WHILE (.TRUE.)
              READ(LUN, '(A)', IOSTAT=ERRNUM) LINE
              IF (ERRNUM /= 0) EXIT
              IF (LINE(1:1) == '!' .OR. TRIM(LINE) == '') CYCLE
              
              IF (LINE(1:6) == '@PARAM') THEN
                 READ(LINE(7:), *, IOSTAT=ERRNUM) CNRF_BC, Opt_bc, &
                      P_FOM, P_E, P_F
                 CYCLE
              END IF
              
              ! Expected Format: AppDate Amount Depth FLoss FCarbon FLabile MRT1 MRT2 CN_BC
              NumApps = NumApps + 1
              IF (NumApps > MaxApp) EXIT
              
              READ(LINE, *, IOSTAT=ERRNUM) &
                BC_Apps(NumApps)%AppDate, &
                BC_Apps(NumApps)%Amount, &
                BC_Apps(NumApps)%Depth, &
                BC_Apps(NumApps)%FLoss, &
                BC_Apps(NumApps)%FCarbon, &
                BC_Apps(NumApps)%FLabile, &
                BC_Apps(NumApps)%MRT_Labile, &
                BC_Apps(NumApps)%MRT_Recalc, &
                BC_Apps(NumApps)%CN_BC
            END DO
            CLOSE(LUN)
          END IF
        END IF
        
        ! Initialize Output
        CALL GetLun('BIOCHAR.OUT', LUN_BC)
        OPEN(UNIT=LUN_BC, FILE='BIOCHAR.OUT', STATUS='REPLACE')
        WRITE(LUN_BC, '(A)') '*BIOCHAR SIMULATION OUTPUT'
        WRITE(LUN_BC, '(A)') '@YEAR DOY   DAS   BC_Labile   BC_Recalc' // &
                             '      dlt_CO2     dlt_Biom    dlt_Hum     dlt_N_Net       TF     WF     NF'
        
      END SUBROUTINE Biochar_Init

!=======================================================================
      SUBROUTINE Biochar_Daily(CONTROL, SOILPROP, SW, ST, NH4, NO3, &
                               IMM, MNR)
        TYPE(ControlType), INTENT(IN) :: CONTROL
        TYPE(SoilType),    INTENT(IN) :: SOILPROP
        REAL, DIMENSION(NL), INTENT(IN) :: SW, ST, NH4, NO3
        REAL, DIMENSION(0:NL, NELEM), INTENT(INOUT) :: IMM, MNR
        
        INTEGER :: YRDOY, DAS, YEAR, DOY, L, iApp
        REAL :: DecayRate1, DecayRate2
        REAL :: dltBC1, dltBC2, dltBC_Total
        REAL :: dlt_bc_CO2, dlt_bc_biom, dlt_bc_hum
        REAL :: AppliedLabile, AppliedRecalc
        REAL :: TotalLabile, TotalRecalc
        REAL :: TF, WF, NF, MF ! Environmental factors
        REAL :: Navail, SoilBCL
        REAL :: Ln2
        ! N-Balance variables
        REAL :: dlt_nbc_need, dlt_nbc_released, dlt_nbc
        REAL :: CN_BC_App

        Ln2 = LOG(2.0)
        Daily_CO2_Gross = 0.0
        Daily_Biom_Gross = 0.0
        Daily_Hum_Gross = 0.0
        Daily_N_Net = 0.0

        YRDOY = CONTROL % YRDOY
        DAS   = CONTROL % DAS
        CALL YR_DOY(YRDOY, YEAR, DOY)

        ! 1. Check for Applications
        DO iApp = 1, NumApps
          IF (BC_Apps(iApp)%AppDate == YRDOY) THEN
             AppliedLabile = BC_Apps(iApp)%Amount * (1.0 - BC_Apps(iApp)%FLoss) * &
                             BC_Apps(iApp)%FCarbon * BC_Apps(iApp)%FLabile
             
             AppliedRecalc = BC_Apps(iApp)%Amount * (1.0 - BC_Apps(iApp)%FLoss) * &
                             BC_Apps(iApp)%FCarbon * (1.0 - BC_Apps(iApp)%FLabile)

             CALL DistributeBiochar(AppliedLabile, AppliedRecalc, BC_Apps(iApp)%Depth, SOILPROP)
          END IF
        END DO

        ! 2. Decay
        DO L = 1, SOILPROP%NLAYR
           IF (SIZE(BC_Apps) > 0 .AND. NumApps > 0) THEN
              
              ! --- Environmental Modifiers ---
              ! WF (Water Factor) - Standard DSSAT SWFAC logic
              WF = 0.0
              IF (SW(L) > SOILPROP%LL(L)) THEN
                 WF = (SW(L) - SOILPROP%LL(L)) / (SOILPROP%DUL(L) - SOILPROP%LL(L))
                 WF = MIN(1.0, WF)
              END IF
              
              ! TF (Temperature Factor) - Lloyd & Taylor
              TF = 0.0
              IF (ST(L) > -10.0) THEN
                 TF = EXP(308.56 * (1.0/56.02 - 1.0/(ST(L) + 46.02)))
              END IF
              TF = MAX(0.0, TF)

              ! NF (Nitrogen Factor)
              SoilBCL = BC_Labile(L)
              Navail = NH4(L) + NO3(L)
              
              NF = 1.0
              IF (Navail > 0.001) THEN
                 IF (Opt_bc > 0.0 .AND. SoilBCL > 1.E-6) THEN
                    NF = MIN(1.0, EXP(-CNRF_BC * ((SoilBCL/Navail) - Opt_bc)/Opt_bc))
                 END IF
              ELSE
                 ! If no N available, decomposition inhibited? 
                 ! Or minimal rate? Archontoulis implies limitation.
                 IF (SoilBCL > 1.E-6) NF = 0.0
              END IF
              
              ! Combined Factor
              MF = WF * TF * NF
              
              ! --- Apply Decay ---
              dltBC1 = 0.0
              dltBC2 = 0.0
              CN_BC_App = BC_Apps(1)%CN_BC

              ! Labile Pool Decay
              IF (BC_Apps(1)%MRT_Labile > 0. .AND. BC_Labile(L) > 1.E-6) THEN
                 DecayRate1 = (Ln2 / (BC_Apps(1)%MRT_Labile * 365.0)) * MF
                 dltBC1 = BC_Labile(L) * (1.0 - EXP(-DecayRate1))
                 BC_Labile(L) = BC_Labile(L) - dltBC1
              END IF
              
              ! Recalcitrant Pool Decay
              IF (BC_Apps(1)%MRT_Recalc > 0. .AND. BC_Recalc(L) > 1.E-6) THEN
                 DecayRate2 = (Ln2 / (BC_Apps(1)%MRT_Recalc * 365.0)) * MF
                 dltBC2 = BC_Recalc(L) * (1.0 - EXP(-DecayRate2))
                 BC_Recalc(L) = BC_Recalc(L) - dltBC2
              END IF
              
              dltBC_Total = dltBC1 + dltBC2
              
              ! --- Partitioning ---
              dlt_bc_CO2  = dltBC_Total * (1.0 - EF_BC)
              dlt_bc_biom = dltBC_Total * EF_BC * FR_BCBIOM
              dlt_bc_hum  = dltBC_Total * EF_BC * (1.0 - FR_BCBIOM)
              
              ! --- N Balance (Mineralization / Immobilization) ---
              ! Eq 7: N Need
              dlt_nbc_need = (dlt_bc_biom / CN_BIOM) + (dlt_bc_hum / CN_HUM)
              
              ! Eq 8: N Released
              dlt_nbc_released = 0.0
              IF (CN_BC_App > 0.0) THEN
                dlt_nbc_released = dltBC_Total / CN_BC_App
              END IF
              
              ! Eq 9: Net Flux
              dlt_nbc = dlt_nbc_released - dlt_nbc_need
              
              ! Update IMM/MNR arrays for SoilNi
              IF (dlt_nbc > 0.0) THEN
                 ! Net Mineralization
                 MNR(L, 1) = MNR(L, 1) + dlt_nbc
              ELSE
                 ! Net Immobilization
                 IMM(L, 1) = IMM(L, 1) + ABS(dlt_nbc)
              END IF
              
              ! Accumulate Profile Totals
              Daily_CO2_Gross  = Daily_CO2_Gross + dlt_bc_CO2
              Daily_Biom_Gross = Daily_Biom_Gross + dlt_bc_biom
              Daily_Hum_Gross  = Daily_Hum_Gross + dlt_bc_hum
              Daily_N_Net      = Daily_N_Net + dlt_nbc

           END IF
        END DO

        ! 3. Output
        TotalLabile = SUM(BC_Labile)
        TotalRecalc = SUM(BC_Recalc)
        
        IF (CONTROL%DYNAMIC == INTEGR .OR. CONTROL%DYNAMIC == OUTPUT) THEN
           WRITE(LUN_BC, '(I5, 1X, I3.3, 1X, I5, 2(F12.2), 4(F12.4), 3(1X, F6.3))') &
                 YEAR, DOY, DAS, TotalLabile, TotalRecalc, &
                 Daily_CO2_Gross, Daily_Biom_Gross, Daily_Hum_Gross, &
                 Daily_N_Net, TF, WF, NF
        END IF

      END SUBROUTINE Biochar_Daily

!=======================================================================
      SUBROUTINE DistributeBiochar(Labile, Recalc, Depth, SOILPROP)
        REAL, INTENT(IN) :: Labile, Recalc, Depth
        TYPE(SoilType), INTENT(IN) :: SOILPROP
        INTEGER :: L
        REAL :: Thickness, ProfileDepth, LayerDepth
        REAL :: Fraction
        REAL :: DistDepth

        ProfileDepth = 0.0
        DO L = 1, SOILPROP%NLAYR
           ProfileDepth = ProfileDepth + SOILPROP%DLAYR(L)
        END DO

        LayerDepth = 0.0
        DO L = 1, SOILPROP%NLAYR
           Thickness = SOILPROP%DLAYR(L)
           IF (LayerDepth < Depth) THEN
            DistDepth = MIN(LayerDepth + Thickness, Depth) - LayerDepth
              IF (DistDepth > 0) THEN
                 Fraction = DistDepth / Depth
                 BC_Labile(L) = BC_Labile(L) + Labile * Fraction
                 BC_Recalc(L) = BC_Recalc(L) + Recalc * Fraction
              END IF
           END IF
           LayerDepth = LayerDepth + Thickness
        END DO
      
      END SUBROUTINE DistributeBiochar

      SUBROUTINE GetBiocharPriming(Layer, BC_Tot, P_Rate, P_Eff, P_Biom)
        INTEGER, INTENT(IN) :: Layer
        REAL, INTENT(OUT)   :: BC_Tot   ! Biochar Total (kg C/ha)
        REAL, INTENT(OUT)   :: P_Rate   ! Rate Modifier (>= 1.0)
        REAL, INTENT(OUT)   :: P_Eff    ! Efficiency Modifier (<= 1.0)
        REAL, INTENT(OUT)   :: P_Biom   ! Partitioning Modifier (<= 1.0)
        
        BC_Tot = BC_Labile(Layer) + BC_Recalc(Layer)
        
        ! Eq 4: Xbc = X * (1 + PFOM * SoilBC / 10000)
        P_Rate = 1.0 + P_FOM * (BC_Tot / 10000.0)
        
        ! Eq 5: ef_fombc = ef_fom * (1 + Pe * SoilBC / 10000)
        ! Assuming P_E is positive coefficient for reduction
        P_Eff  = 1.0 + P_E * (BC_Tot / 10000.0)
        P_Eff  = MAX(0.0, P_Eff) 
        
        ! Eq 6: fr_fom_biom
        ! Assuming P_F is positive coefficient for reduction (or change)
        ! Paper implies negative priming on efficiency/transfer
        P_Biom = 1.0 + P_F * (BC_Tot / 10000.0)
        P_Biom = MAX(0.0, P_Biom)

      END SUBROUTINE GetBiocharPriming

      END MODULE Biochar_mod
