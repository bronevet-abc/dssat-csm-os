      MODULE Biochar_mod
!=======================================================================
!  MODULE Biochar_mod
!  Purpose: Simulate Biochar application and decay in soil.
!  Tracks biochar pools and partitions decay into CO2, Biomass, and Humic pools.
!  Calculates N Mineralization/Immobilization fluxes.
!=======================================================================
      USE ModuleDefs
      USE ModuleData
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
        REAL    :: CEC_INIT      ! Initial CEC (cmol/kg biochar)
        REAL    :: BCLV          ! Biochar Liming Value (cmol/kg biochar)
        REAL    :: Kads          ! Langmuir Adsorption Coeff (L/mg)
        REAL    :: Kdes          ! Langmuir Desorption Coeff (L/mg)
        REAL    :: QLL           ! Quality mod for Lower Limit
        REAL    :: KDUL          ! Quality mod for DUL
        REAL    :: KBD           ! Quality mod for Bulk Density
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
      
      ! CEC Parameters
      REAL :: CEC_MAX = 100.0 ! Maximum CEC after aging (cmol/kg)
      REAL :: K_CEC   = 0.001 ! CEC aging rate constant (1/day)
      
      ! pH Parameters
      REAL :: BCLV_Default = 50.0 ! Default Liming Value if not specified
      REAL :: UpH     = 8.3
      REAL :: LpH     = 3.5
      REAL :: P1_pH   = 10.0
      
      ! NH4 Adsorption Parameters
      REAL :: Kads_Default = 0.006
      REAL :: Kdes_Default = 0.006
      
!     Priming Parameters (Archontoulis et al., 2015)
      REAL :: P_FOM = 0.0   ! Positive priming on FOM decomposition rate
      REAL :: P_E   = 0.0   ! Negative priming on Carbon Efficiency
      REAL :: P_F   = 0.0   ! Negative priming on FOM->BIOM transfer

!     State Variables
      REAL, DIMENSION(NL) :: dlt_nbc_released_yesterday
      REAL, DIMENSION(NL) :: BC_Labile   ! Labile Biochar C (kg/ha)
      REAL, DIMENSION(NL) :: BC_Recalc   ! Recalcitrant Biochar C (kg/ha)
      REAL, DIMENSION(NL) :: BC_NH4_Ads  ! Adsorbed NH4 (kg N/ha)
      
!     Flux Variables (Daily)
      REAL :: Daily_CO2_Gross   ! Total daily CO2 emission from Biochar (kg C/ha/d)
      REAL :: Daily_Biom_Gross  ! Total daily flux to Biomass (kg C/ha/d)
      REAL :: Daily_Hum_Gross   ! Total daily flux to Humic (kg C/ha/d)
      REAL :: Daily_N_Net       ! Net N mineralization (+)/immobilization (-) (kg N/ha/d)

!     Output file unit
      INTEGER :: LUN_BC
      LOGICAL :: FirstOutput = .TRUE.
      LOGICAL :: FirstRun = .TRUE.

      CONTAINS

!=======================================================================
      SUBROUTINE Biochar_Init(CONTROL)
        TYPE(ControlType), INTENT(IN) :: CONTROL
        INTEGER :: ERRNUM, LUN
        CHARACTER(LEN=120) :: LINE
        LOGICAL :: FEXIST
        
        ! Initialize State
        dlt_nbc_released_yesterday = 0.0
        BC_Labile = 0.0
        BC_Recalc = 0.0
        BC_NH4_Ads = 0.0
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
              LINE = ADJUSTL(LINE)
              IF (LINE(1:1) == '!' .OR. TRIM(LINE) == '') CYCLE


              
              IF (LINE(1:6) == '@PARAM') THEN
                 READ(LINE(7:), *, IOSTAT=ERRNUM) CNRF_BC, Opt_bc, &
                      P_FOM, P_E, P_F, CEC_MAX, K_CEC
                 CYCLE
              END IF
              
              ! Expected Format: AppDate Amount Depth FLoss FCarbon FLabile MRT1 MRT2 CN_BC CEC_INIT BCLV Kads Kdes QLL KDUL KBD
              NumApps = NumApps + 1
              IF (NumApps > MaxApp) EXIT
              
              READ(LINE, *, IOSTAT=ERRNUM) BC_Apps(NumApps)%AppDate, &
                  BC_Apps(NumApps)%Amount, BC_Apps(NumApps)%Depth, &
                  BC_Apps(NumApps)%FLoss, BC_Apps(NumApps)%FCarbon, &
                  BC_Apps(NumApps)%FLabile, BC_Apps(NumApps)%MRT_Labile,&
                  BC_Apps(NumApps)%MRT_Recalc, BC_Apps(NumApps)%CN_BC, &
                  BC_Apps(NumApps)%CEC_INIT, BC_Apps(NumApps)%BCLV, &
                  BC_Apps(NumApps)%Kads, BC_Apps(NumApps)%Kdes, &
                   BC_Apps(NumApps)%QLL, BC_Apps(NumApps)%KDUL, &
                   BC_Apps(NumApps)%KBD
              
               IF (ERRNUM /= 0) THEN
                   NumApps = NumApps - 1
                   CYCLE
               END IF

             ! Default Kads/Kdes if missing
             IF (BC_Apps(NumApps)%Kads < 1.E-9) BC_Apps(NumApps)%Kads = Kads_Default
             IF (BC_Apps(NumApps)%Kdes < 1.E-9) BC_Apps(NumApps)%Kdes = Kdes_Default

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
        
        FirstRun = .TRUE.
      END SUBROUTINE Biochar_Init

!=======================================================================
      SUBROUTINE Biochar_Daily(CONTROL, SOILPROP, SW, ST, NH4, NO3, &
                               IMM, MNR)
        TYPE(ControlType), INTENT(IN) :: CONTROL
        TYPE(SoilType),    INTENT(INOUT) :: SOILPROP
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
        
        ! CEC Local Vars
        REAL :: SoilMass, CurrentCEC_Soil, TotalCEC_Soil
        REAL :: TotalCEC_BC, TotalMass_BC, TotalMass_BC_InLayer
        REAL :: WeightedCEC_BC, AvgCEC_BC, CurrentMass_BC
        REAL :: AppDepth, LayerTop, LayerBottom, Fraction, DistDepth
        REAL :: MassApplied, MassInLayer, CEC_t
        INTEGER :: TimeSinceApp, TIMDIF, L2
        REAL :: Age
        ! pH Local Vars
        REAL :: SoilpH, SoilCECBC_Val, Term1, Term2, dpH, AppBCLV
        
        ! NH4 Adsorption Local Vars
        REAL :: NH4_Conc_mgL, NH4_Ads_Target_mgL, NH4_Ads_Target_kgHa
        REAL :: VolSW_L_Ha, CEC_Ratio
        REAL :: AppKads, AppKdes, Ads_Diff
        
        REAL, DIMENSION(NL) :: NativeCEC ! To store initial soil CEC

        Ln2 = LOG(2.0)
        Daily_CO2_Gross = 0.0
        Daily_Biom_Gross = 0.0
        Daily_Hum_Gross = 0.0
        Daily_N_Net = 0.0

        YRDOY = CONTROL % YRDOY
        DAS   = CONTROL % DAS
        CALL YR_DOY(YRDOY, YEAR, DOY)
        
        ! Initialize NativeCEC on first run
        IF (FirstRun) THEN
           NativeCEC = SOILPROP%CEC
           FirstRun = .FALSE.
        END IF
        ! Biochar Daily Simulation
      
      ! 1. Check for new applications
      IF (NumApps > 0) THEN
         DO iApp = 1, NumApps
            IF (TIMDIF(BC_Apps(iApp)%AppDate, YRDOY) == 0) THEN

             AppliedLabile = BC_Apps(iApp)%Amount * (1.0 - BC_Apps(iApp)%FLoss) * &
                             BC_Apps(iApp)%FCarbon * BC_Apps(iApp)%FLabile
             
             AppliedRecalc = BC_Apps(iApp)%Amount * (1.0 - BC_Apps(iApp)%FLoss) * &
                             BC_Apps(iApp)%FCarbon * (1.0 - BC_Apps(iApp)%FLabile)

             CALL DistributeBiochar(AppliedLabile, AppliedRecalc, BC_Apps(iApp)%Depth, SOILPROP)
             
             ! Biochar effects on Soil pH (Eq 12) - Apply ONLY on application day
             ! Calculate incremental effect of THIS application
             AppBCLV = BC_Apps(iApp)%BCLV
             IF (AppBCLV < 1.E-6) AppBCLV = BCLV_Default
             
             AppDepth = BC_Apps(iApp)%Depth
             
             DO L = 1, SOILPROP%NLAYR
                ! Calculate Fraction of this App in this Layer
                LayerTop = 0.0
                DO L2 = 1, L-1
                   LayerTop = LayerTop + SOILPROP%DLAYR(L2)
                END DO
                LayerBottom = LayerTop + SOILPROP%DLAYR(L)
                
                IF (LayerTop < AppDepth) THEN
                   DistDepth = MIN(LayerBottom, AppDepth) - LayerTop
                   IF (DistDepth > 0) THEN
                      Fraction = DistDepth / AppDepth
                      
                      ! Mass of THIS application in this layer (kg/ha)
                      ! used for Massfr in Eq 12
                      MassApplied = BC_Apps(iApp)%Amount
                      MassInLayer = MassApplied * Fraction
                      
                      ! Soil Mass (kg/ha)
                      SoilMass = SOILPROP%BD(L) * SOILPROP%DLAYR(L) * 100000.0
                      
                      ! Mass Fraction (g/g) = BC Mass / Soil Mass
                      ! MassInLayer (kg/ha) / SoilMass (kg/ha) -> g/g
                      IF (SoilMass > 0.0) THEN
                          Fraction = MassInLayer / SoilMass
                          
                          SoilpH = SOILPROP%PH(L)
                          SoilCECBC_Val = SOILPROP%CEC(L)
                          
                          IF (SoilCECBC_Val > 1.E-4 .AND. (UpH - LpH) > 1.E-4) THEN
                             Term1 = (Fraction * AppBCLV) / SoilCECBC_Val
                             
                             IF (SoilpH > LpH .AND. SoilpH < UpH) THEN
                                Term2 = ((UpH - SoilpH) * (SoilpH - LpH)) / (UpH - LpH)
                                dpH = P1_pH * Term1 * Term2
                                SOILPROP%PH(L) = SoilpH + dpH
                             END IF
                          END IF
                      END IF
                   END IF
                END IF
             END DO
          END IF
        END DO
      END IF

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
              TF = (MAX(0.0, ST(L)) / 32.0) ** 2
              TF = MIN(1.0, TF)
              TF = MAX(0.0, TF)

              ! NF (Nitrogen Factor)
              SoilBCL = BC_Labile(L)
              Navail = NH4(L) + NO3(L) + dlt_nbc_released_yesterday(L)
              
              NF = 1.0
              IF (Navail > 0.001) THEN
                 IF (Opt_bc > 0.0 .AND. SoilBCL > 1.E-6) THEN
                    NF = MIN(1.0, EXP(-CNRF_BC * ((SoilBCL/Navail) - Opt_bc)/Opt_bc))
                 END IF
              ELSE
                 IF (SoilBCL > 1.E-6) NF = 0.0
              END IF
              
              ! Combine Factor
              MF = WF * TF * NF
              
              ! --- Apply Decay ---
              dltBC1 = 0.0
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
              
              dlt_nbc_released_yesterday(L) = dlt_nbc_released

           END IF
        END DO

        ! 3. Update Soil CEC (Eq 10 & 11)
        IF (NumApps > 0) THEN
           DO L = 1, SOILPROP%NLAYR
              ! Soil Mass (kg/ha)
              SoilMass = SOILPROP%BD(L) * SOILPROP%DLAYR(L) * 100000.0
              
              TotalCEC_BC = 0.0
              TotalMass_BC_InLayer = 0.0
              WeightedCEC_BC = 0.0
              
              IF (BC_Labile(L) + BC_Recalc(L) > 1.E-6) THEN
                 DO iApp = 1, NumApps
                    AppDepth = BC_Apps(iApp)%Depth
                    
                    ! Calculate Fraction of this App in this Layer
                    LayerTop = 0.0
                    DO L2 = 1, L-1
                       LayerTop = LayerTop + SOILPROP%DLAYR(L2)
                    END DO
                    LayerBottom = LayerTop + SOILPROP%DLAYR(L)
                    
                    IF (LayerTop < AppDepth) THEN
                       DistDepth = MIN(LayerBottom, AppDepth) - LayerTop
                       IF (DistDepth > 0) THEN
                          Fraction = DistDepth / AppDepth
                          
                          ! Mass Applied to this layer (Initial)
                          MassApplied = BC_Apps(iApp)%Amount
                          MassInLayer = MassApplied * Fraction
                          
                          ! Time Since Application
                          ! TIMDIF returns Diff in Days. AppDate is YRDOY.
                          TimeSinceApp = TIMDIF(BC_Apps(iApp)%AppDate, YRDOY)
                          
                          IF (TimeSinceApp >= 0) THEN
                              ! Eq 10: Aging
                              ! CEC_t = CEC_min + (CEC_max - CEC_min) * (1 - exp(-k * t))
                              Age = REAL(TimeSinceApp)
                              CEC_t = BC_Apps(iApp)%CEC_INIT + &
                                      (CEC_MAX - BC_Apps(iApp)%CEC_INIT) * &
                                      (1.0 - EXP(-K_CEC * Age))
                              
                              WeightedCEC_BC = WeightedCEC_BC + CEC_t * MassInLayer
                              TotalMass_BC_InLayer = TotalMass_BC_InLayer + MassInLayer
                          END IF
                       END IF
                    END IF
                 END DO
                 
                 IF (TotalMass_BC_InLayer > 0.0) THEN
                    AvgCEC_BC = WeightedCEC_BC / TotalMass_BC_InLayer
                    CurrentMass_BC = (BC_Labile(L) + BC_Recalc(L)) / MAX(0.1, BC_Apps(1)%FCarbon)
                    
                    TotalCEC_BC = AvgCEC_BC * CurrentMass_BC
                    
                    ! Eq 11: Mixing (Mass Weighted)
                    SOILPROP%CEC(L) = (NativeCEC(L) * SoilMass + TotalCEC_BC) / (SoilMass + CurrentMass_BC)
                 END IF
              END IF
           END DO
        END IF

        ! 4. Update NH4 Adsorption (Eq 13)
        DO L = 1, SOILPROP%NLAYR
           IF (BC_Labile(L) + BC_Recalc(L) > 1.E-6) THEN
              VolSW_L_Ha = SW(L) * SOILPROP%DLAYR(L) * 10.0 * 10000.0
              
              IF (VolSW_L_Ha > 1.0) THEN
                 ! NH4 Conc (mg/L)
                 NH4_Conc_mgL = (NH4(L) * 1.0E6) / VolSW_L_Ha
                 
                 AppKads = BC_Apps(1)%Kads
                 AppKdes = BC_Apps(1)%Kdes
                 
                 ! Ratio of CECbc / CECsoil (Using Ratio of Current Soil CEC to Native?)
                 ! Paper implies amplification factor.
                 IF (NativeCEC(L) > 1.E-6) THEN
                    CEC_Ratio = SOILPROP%CEC(L) / NativeCEC(L)
                    
                    ! Adsorption Potential
                    Term1 = AppKads * CEC_Ratio
                    NH4_Ads_Target_mgL = NH4_Conc_mgL * Term1 / (1.0 + Term1)
                    
                    NH4_Ads_Target_kgHa = (NH4_Ads_Target_mgL * VolSW_L_Ha) / 1.0E6
                    
                    Ads_Diff = NH4_Ads_Target_kgHa - BC_NH4_Ads(L)
                    
                    IF (Ads_Diff > 0.0) THEN
                       ! Adsorption
                       Ads_Diff = MIN(Ads_Diff, NH4(L))
                       IMM(L, 1) = IMM(L, 1) + Ads_Diff
                       BC_NH4_Ads(L) = BC_NH4_Ads(L) + Ads_Diff
                    ELSEIF (Ads_Diff < 0.0) THEN
                       ! Desorption (Mineralization)
                        Term1 = AppKdes * CEC_Ratio
                        NH4_Ads_Target_mgL = NH4_Conc_mgL * Term1 / (1.0 + Term1)
                        NH4_Ads_Target_kgHa = (NH4_Ads_Target_mgL * VolSW_L_Ha) / 1.0E6
                        
                        Ads_Diff = NH4_Ads_Target_kgHa - BC_NH4_Ads(L)
                        
                        IF (Ads_Diff < 0.0) THEN
                            Ads_Diff = MAX(Ads_Diff, -BC_NH4_Ads(L))
                            MNR(L, 1) = MNR(L, 1) + ABS(Ads_Diff)
                            BC_NH4_Ads(L) = BC_NH4_Ads(L) + Ads_Diff
                        END IF
                    END IF
                 END IF
              END IF
           END IF
        END DO

        ! 5. Output
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

!=======================================================================
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
        ! Assuming P_F is positive coefficient for reduction
        P_Biom = 1.0 + P_F * (BC_Tot / 10000.0)
        P_Biom = MAX(0.0, P_Biom)

      END SUBROUTINE GetBiocharPriming

!=======================================================================
!  Biochar_UpdateSoilProps
!  Updates Soil Hydraulic Properties (BD, LL, DUL, SAT) based on 
!  Biochar applications using Saxton & Rawls (2006).
!=======================================================================
      SUBROUTINE Biochar_UpdateSoilProps(SOILPROP)
          TYPE(SoilType), INTENT(INOUT) :: SOILPROP
          
          INTEGER :: L
          REAL, DIMENSION(NL) :: BC_Mass_g_g ! Biochar mass fraction (g/g soil)
          REAL :: Sand, Clay, OrgC_Native, OrgC_Eff
          REAL :: BD_Nat, LL_Nat, DUL_Nat, SAT_Nat
          REAL :: BD_BC,  LL_BC,  DUL_BC,  SAT_BC
          REAL :: SoilMass, BC_Total_kgHa, BC_Total_C_kgHa
          REAL :: QLL_Avg, KDUL_Avg, KBD_Avg
          REAL :: Biochar_OM_Eff
          
          IF (NumApps == 0 .OR. SUM(BC_Labile) + SUM(BC_Recalc) < 1.E-6) RETURN
          
          ! Use weighted average Q parameters if multiple apps?
          ! For simplicity, use First Application parameters for now
          QLL_Avg  = BC_Apps(1)%QLL
          KDUL_Avg = BC_Apps(1)%KDUL
          KBD_Avg  = BC_Apps(1)%KBD
          
          DO L = 1, SOILPROP%NLAYR
             ! 1. Calculate Biochar Fraction

             SoilMass = SOILPROP%BD(L) * SOILPROP%DLAYR(L) * 100000.0 ! kg/ha
             BC_Total_C_kgHa = BC_Labile(L) + BC_Recalc(L)
             
             ! Estimate Total Mass of Biochar
             ! Mass = C / FCarbon. Use Avg FCarbon.
             IF (BC_Apps(1)%FCarbon > 0.0) THEN
                 BC_Total_kgHa = BC_Total_C_kgHa / BC_Apps(1)%FCarbon
             ELSE
                 BC_Total_kgHa = BC_Total_C_kgHa / 0.7 ! Default
             END IF
             
             IF (SoilMass > 1.0) THEN
                 BC_Mass_g_g(L) = BC_Total_kgHa / SoilMass
                 
                 ! 2. Get Native Properties
                  Sand = SOILPROP%SAND(L) / 100.0
                  Clay = SOILPROP%CLAY(L) / 100.0
                  OrgC_Native = SOILPROP%OC(L) / 100.0

                 
                 ! Saxton outputs
                 CALL SaxtonRawls(Sand, Clay, OrgC_Native, &
                                  BD_Nat, LL_Nat, DUL_Nat, SAT_Nat)
                 
                 ! 3. Calculate Effective Properties with Biochar
                 ! Biochar OM% added = BC_Mass_g_g * 100.0
                 ! Apply KBD modifier for density
                  ! Biochar OM% added = BC_Mass_g_g * 100.0
                  ! Apply KBD modifier for density
                  ! NOTE: SaxtonRawls expects fractions (0-1), so keep as fraction.
                  Biochar_OM_Eff = BC_Mass_g_g(L) * EXP(KBD_Avg)
 
                 
                 ! Recalculate BD with Eff OM
                 CALL SaxtonRawls(Sand, Clay, OrgC_Native + Biochar_OM_Eff, &
                                  BD_BC, LL_BC, DUL_BC, SAT_BC)
                 
                 ! Update BD
                 ! Delta approach:
                 SOILPROP%BD(L) = SOILPROP%BD(L) + (BD_BC - BD_Nat)
                 SOILPROP%BD(L) = MAX(0.5, MIN(SOILPROP%BD(L), 2.0))
                 
                 ! Update Water Limits
                 ! DUL
                  ! Update Water Limits
                  ! DUL
                  Biochar_OM_Eff = BC_Mass_g_g(L) * EXP(KDUL_Avg)

                 CALL SaxtonRawls(Sand, Clay, OrgC_Native + Biochar_OM_Eff, &
                                  BD_BC, LL_BC, DUL_BC, SAT_BC)
                 SOILPROP%DUL(L) = SOILPROP%DUL(L) + (DUL_BC - DUL_Nat)
                 
                 ! LL
                 Biochar_OM_Eff = BC_Mass_g_g(L) * EXP(QLL_Avg)
                 CALL SaxtonRawls(Sand, Clay, OrgC_Native + Biochar_OM_Eff, &
                                  BD_BC, LL_BC, DUL_BC, SAT_BC)
                 SOILPROP%LL(L) = SOILPROP%LL(L) + (LL_BC - LL_Nat)
                 
                 ! SAT (Driven by BD primarily)
                 SOILPROP%SAT(L) = 1.0 - (SOILPROP%BD(L)/2.65)
                 SOILPROP%SAT(L) = MAX(SOILPROP%DUL(L)+0.01, SOILPROP%SAT(L))
                 
                 ! Update Porosity
                 SOILPROP%POROS(L) = SOILPROP%SAT(L)
                 
                 ! Recalculate KG2PPM
                 SOILPROP%KG2PPM(L) = 10.0 / (SOILPROP%BD(L) * SOILPROP%DLAYR(L)) 
             END IF
          END DO
      END SUBROUTINE Biochar_UpdateSoilProps

!=======================================================================
!     Saxton & Rawls (2006) Pedotransfer Functions
!=======================================================================
      SUBROUTINE SaxtonRawls(Sand, Clay, OM, BD, LL, DUL, SAT)
          REAL, INTENT(IN)  :: Sand, Clay, OM
          REAL, INTENT(OUT) :: BD, LL, DUL, SAT
          
          REAL :: S, C, O
          REAL :: Theta_1500, Theta_1500t
          REAL :: Theta_33, Theta_33t
          REAL :: Theta_S_33, Theta_S_33t
          
          S = Sand

          C = Clay
          O = OM
          
          ! 1. Permanent Wilting Point (LL) -> Theta_1500 (1500 kPa)
          Theta_1500t = -0.024 * S + 0.487 * C + 0.006 * O + &
                        0.005 * (S * O) - 0.013 * (C * O) + &
                        0.068 * (S * C) + 0.031
          Theta_1500 = Theta_1500t + (0.14 * Theta_1500t - 0.02)
          LL = Theta_1500
          
          ! 2. Field Capacity (DUL) -> Theta_33 (33 kPa)
          Theta_33t = -0.251 * S + 0.195 * C + 0.011 * O + &
                      0.006 * (S * O) - 0.027 * (C * O) + &
                      0.452 * (S * C) + 0.299
          Theta_33 = Theta_33t + (1.283 * Theta_33t * Theta_33t - &
                     0.374 * Theta_33t - 0.015)
          DUL = Theta_33
          
          ! 3. Saturation (SAT)
          Theta_S_33t = 0.278 * S + 0.034 * C + 0.022 * O - &
                        0.018 * (S * O) - 0.027 * (C * O) - &
                        0.584 * (S * C) + 0.078
          Theta_S_33 = Theta_S_33t + (0.636 * Theta_S_33t - 0.107)
          
          SAT = DUL + Theta_S_33 - 0.097 * S + 0.043
          
          ! Now BD
          BD = (1.0 - SAT) * 2.65
          
          ! Sanity Checks

          SAT = MIN(0.9, MAX(0.2, SAT))
          BD  = (1.0 - SAT) * 2.65
          BD  = MIN(2.0, MAX(0.8, BD))
          DUL = MIN(SAT-0.01, MAX(0.05, DUL))
          LL  = MIN(DUL-0.01, MAX(0.01, LL))
          
      END SUBROUTINE SaxtonRawls

      END MODULE Biochar_mod
