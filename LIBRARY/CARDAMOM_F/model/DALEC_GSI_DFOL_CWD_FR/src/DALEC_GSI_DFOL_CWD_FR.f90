
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL        &
         ,dble_one,dble_zero  &
         ,vsmall              &
         ,arrhenious          &
         ,acm_gpp             &
         ,acm_albedo_gc       &
         ,daylength_in_hours  &
         ,daylength_hours     &
         ,daylength_seconds   &
         ,freeze              &
         ,co2comp_saturation  &
         ,co2comp_half_sat_conc &
         ,co2_half_saturation &
         ,co2_compensation_point &
         ,kc_saturation       &
         ,kc_half_sat_conc    &
         ,calculate_Rtot      &
         ,calculate_aerodynamic_conductance &
         ,linear_model_gradient &
         ,seconds_per_day     &
         ,seconds_per_hour    &
         ,seconds_per_step    &
         ,root_biomass        &
         ,root_reach          &
         ,minlwp              &
         ,min_root            &
         ,max_depth           &
         ,root_k              &
         ,top_soil_depth      &
         ,soil_depth          &
         ,previous_depth      &
         ,nos_root_layers     &
         ,deltat_1            &
         ,water_flux          &
         ,layer_thickness     &
         ,min_layer        &
         ,nos_soil_layers  &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,meant            &
         ,meant_K          &
         ,stomatal_conductance &
         ,avN              &
         ,iWUE             &
         ,NUE              &
         ,pn_max_temp      &
         ,pn_opt_temp      &
         ,pn_kurtosis      &
         ,e0               &
         ,co2_half_sat     &
         ,co2_comp_point   &
         ,max_lai_lwrad_absorption  &
         ,lai_half_lwrad_absorption &
         ,max_lai_nir_absorption    &
         ,lai_half_nir_absorption   &
         ,max_lai_par_absorption    &
         ,lai_half_par_absorption   &
         ,max_lai_swrad_reflected   &
         ,lai_half_swrad_reflected  &
         ,lai_half_lwrad_to_sky     &
         ,soil_swrad_absorption     &
         ,max_lai_lwrad_release     &
         ,lai_half_lwrad_release    &
         ,soilevap_rad_intercept    &
         ,soilevap_rad_coef         &
         ,mint             &
         ,maxt             &
         ,swrad            &
         ,co2              &
         ,doy              &
         ,rainfall         &
         ,wind_spd         &
         ,vpd_pa           &
         ,lai              &
         ,days_per_step    &
         ,days_per_step_1  &
         ,dayl_seconds     &
         ,dayl_hours       &
         ,canopy_storage   &
         ,intercepted_rainfall          &
         ,disturbance_residue_to_litter &
         ,disturbance_residue_to_cwd    &
         ,disturbance_residue_to_som    &
         ,disturbance_loss_from_litter  &
         ,disturbance_loss_from_cwd     &
         ,disturbance_loss_from_som     &
         ,itemp,ivpd,iphoto&
         ,extracted_C      &
         ,dim_1,dim_2      &
         ,nos_trees        &
         ,nos_inputs       &
         ,leftDaughter     &
         ,rightDaughter    &
         ,nodestatus       &
         ,xbestsplit       &
         ,nodepred         &
         ,bestvar

!!!!!!!!!!
! Random Forest GPP emulator
!!!!!!!!!!

! arrays for the emulator, just so we load them once and that is it cos they be
! massive
integer ::    dim_1, & ! dimension 1 of response surface
              dim_2, & ! dimension 2 of response surface
          nos_trees, & ! number of trees in randomForest
         nos_inputs    ! number of driver inputs
double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                    rightDaughter, & ! right daughter for forets
                                                       nodestatus, & ! nodestatus for forests
                                                       xbestsplit, & ! for forest
                                                         nodepred, & ! prediction value for each tree
                                                          bestvar    ! for randomForests
!!!!!!!!!
! Parameters
!!!!!!!!!

! useful technical parameters
double precision, parameter :: xacc = 1d-4        & ! accuracy parameter for zbrent bisection proceedure ! 0.0001
                              ,dble_zero = 0d0    &
                              ,dble_one = 1d0     &
                              ,vsmall = tiny(0d0)

integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, parameter :: pi = 3.1415927d0,    &
                             pi_1 = pi**(-1d0),     &
                              pi2 = pi**2d0,        &
                           two_pi = pi*2d0,         &
                       deg_to_rad = pi/180d0,       &
              sin_dayl_deg_to_rad = sin( 23.45d0 * deg_to_rad ), & ! repeated function in acm
                          gravity = 9.8067d0,       & ! acceleration due to gravity, ms-1
                            boltz = 5.670400d-8,    & ! Boltzmann constant (W.m-2.K-4)
                       emissivity = 0.96d0,         &
                      emiss_boltz = emissivity * boltz, &
                           freeze = 273.15d0,       &
                       gs_H2O_CO2 = 1.646259d0,     & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
                     gs_H2O_CO2_1 = gs_H2O_CO2 ** (-dble_one), &
                       gb_H2O_CO2 = 1.37d0,         & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
          partial_molar_vol_water = 18.05d-6,       & ! partial molar volume of water, m3 mol-1 at 20C
                       umol_to_gC = 1d-6*12d0,      & ! conversion of umolC -> gC
                 mmol_to_kg_water = 1.8d-5,         & ! milli mole conversion to kg
                   mol_to_g_water = 18d0,           & ! molecular mass of water
                     mol_to_g_co2 = 12d0,           & ! molecular mass of CO2 (g)
                     g_to_mol_co2 = 1d0/12d0,       &
!snowscheme       density_of_water = 998.9d0,         & ! density of !water kg.m-3
                   gas_constant_d = 287.04d0,       & ! gas constant for dry air (J.K-1.mol-1)
                             Rcon = 8.3144d0,       & ! Universal gas constant (J.K-1.mol-1)
                        vonkarman = 0.41d0,         & ! von Karman's constant
                      vonkarman_2 = vonkarman**2d0, & ! von Karman's constant^2
                            cpair = 1004.6d0          ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

! photosynthesis / respiration parameters
double precision, parameter :: &
                    kc_saturation = 310d0,        & ! CO2 half saturation, saturation value
                 kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, half sat
               co2comp_saturation = 36.5d0,       & ! CO2 compensation point, saturation
            co2comp_half_sat_conc = 9.46d0,       & ! CO2 comp point, half sat
                                                    ! Each of these are temperature
                                                    ! sensitivty
              leaf_life_weighting = 1d0/2d0,      & ! inverse of averaging period of lagged effects
                                                    ! probably should be an actual parmeter
                      Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                    ! lost as growth respiration
                                                    ! (i.e. 0.28 .eq. xNPP)
                  one_Rg_fraction = dble_one - Rg_fraction
! hydraulic parameters
double precision, parameter :: &
                       tortuosity = 2.5d0,        & ! tortuosity
                           gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                        max_depth = 2d0,          & ! max root depth (m)
                           root_k = 100d0,        & ! root biomass needed to reach 50% depth (gbiomass/m2)
                      root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                    !                 Williams et al 1996 = 0.0001
                    root_radius_1 = root_radius**(-dble_one), &
              root_cross_sec_area = 3.141593d-08, & ! root cross sectional area (m2)
                                                    ! = pi * root_radius * root_radius
                     root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                    ! 0.5e6 Williams et al 1996
                                                    ! 0.31e6 Bonan et al 2014
          root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-1d0), &
               const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                             head = 0.009807d0,   & ! head of pressure (MPa/m)
                           head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

! structural parameters
double precision, parameter :: &
                    canopy_height = 9d0,          & ! canopy height assumed to be 9 m
                     tower_height = canopy_height + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                         min_wind = 0.1d0,        & ! minimum wind speed at canopy top
                        min_layer = 0.01d0,       & ! minimum thickness of the second rooting layer (m)
                      soil_roughl = 0.05d0,       & ! soil roughness length (m)
                   top_soil_depth = 0.3d0,        & ! depth to which we consider the top soil to extend (m)
                         min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                  min_throughfall = 0.2d0           ! minimum fraction of precipitation which
                                                    ! is through fall

! timing parameters
double precision, parameter :: &
                 seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                  seconds_per_day = 86400d0,        & ! Number of seconds per day
                seconds_per_day_1 = 1d0/seconds_per_day       ! Inverse of seconds per day

!!!!!!!!!
! Module level variables
!!!!!!!!!

! management and gsi related values
integer :: gsi_lag_remembered
! local variables for GSI phenology model
double precision :: Tfac,Photofac,VPDfac & ! oC, seconds, Pa
                   ,Rg_from_labile       &
                   ,delta_gsi,tmp,tmp_lai,gradient         &
                   ,fol_turn_crit,lab_turn_crit    &
                   ,gsi_history(22),just_grown

double precision, allocatable, dimension(:) :: extracted_C,itemp,ivpd,iphoto, &
                                               disturbance_residue_to_litter, &
                                               disturbance_residue_to_som,    &
                                               disturbance_residue_to_cwd,    &
                                               disturbance_loss_from_litter,  &
                                               disturbance_loss_from_cwd,     &
                                               disturbance_loss_from_som,     &
                                               tmp_x, tmp_m

! Phenological choices
! See source below for details of these variables
double precision :: root_cost,root_life,       &
                    leaf_cost,leaf_life

! hydraulic model variables
integer :: water_retention_pass, soil_layer
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand, &
                                                layer_thickness    ! thickness of soil layers (m)
double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                         demand, & ! maximum potential canopy hydraulic demand
                                                     water_flux    ! potential transpiration flux (mmol.m-2.s-1)

double precision :: root_reach, root_biomass,soil_depth, &
                    soilRT, &
  new_depth,previous_depth, & ! depth of bottom of soil profile
               canopy_wind, & ! wind speed (m.s-1) at canopy top
                     ustar, & ! friction velocity (m.s-1)
            air_density_kg, & ! air density kg/m3
                    roughl, & ! roughness length (m)
              displacement, & ! zero plane displacement (m)
                max_supply, & ! maximum water supply (mmolH2O/m2/day)
                     meant, & ! mean air temperature (oC)
                   meant_K, & ! mean air temperature (K)
        canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
          soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
            soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
      stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.day-1)
   aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
          soil_conductance, & ! soil surface conductance (m.s-1)
         convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                    lambda, & ! latent heat of vapourisation (J.kg-1)
                     psych, & ! psychrometric constant (kPa K-1)
                     slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
            canopy_storage, & ! water storage on canopy (kg.m-2)
      intercepted_rainfall    ! intercepted rainfall rate equivalent (kg.m-2.s-1)

! Module level variables for ACM_GPP_ET parameters
double precision :: delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                      avN, & ! average foliar N (gN/m2)
                     iWUE, & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
                      NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                             ! ,unlimited by CO2, light and photoperiod
                             ! (gC/gN/m2leaf/day)
              pn_max_temp, & ! Maximum temperature for photosynthesis (oC)
              pn_opt_temp, & ! Optimum temperature fpr photosynthesis (oC)
              pn_kurtosis, & ! Kurtosis of photosynthesis temperature response
                       e0, & ! Quantum yield gC/MJ/m2/day PAR
             co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
           co2_comp_point, & ! CO2 at which photosynthesis > 0 (ppm)
                   minlwp, & ! min leaf water potential (MPa)
 max_lai_lwrad_absorption, & ! Max fraction of LW from sky absorbed by canopy
lai_half_lwrad_absorption, & ! LAI at which canopy LW absorption = 50 %
   max_lai_nir_absorption, & ! Max fraction of NIR absorbed by canopy
  lai_half_nir_absorption, & ! LAI at which canopy NIR absorption = 50 %
   max_lai_par_absorption, & ! Max fraction of PAR absorbed by canopy
  lai_half_par_absorption, & ! LAI at which canopy PAR absorption = 50 %
  max_lai_swrad_reflected, & ! Max fraction of SW rad reflected back to sky
 lai_half_swrad_reflected, & ! LAI at which SW reflection = 50 %
    lai_half_lwrad_to_sky, & ! LAI at which 50 % LW is reflected back to sky
    soil_swrad_absorption, & ! Fraction of SW rad absorbed by soil
    max_lai_lwrad_release, & ! Max fraction of LW emitted from canopy to be
   lai_half_lwrad_release, & ! LAI at which LW emitted from canopy to be released at 50 %
   soilevap_rad_intercept, & ! Intercept (kgH2O/m2/day) on linear adjustment to soil evaporation
                             ! to account for non-calculation of energy balance
       soilevap_rad_coef     ! Coefficient on linear adjustment to
                             ! soil evaporation to account for non-calculation of energy balance

! Module level variables for step specific met drivers
double precision :: mint, & ! minimum temperature (oC)
                    maxt, & ! maximum temperature (oC)
                   swrad, & ! incoming short wave radiation (MJ/m2/day)
                     co2, & ! CO2 (ppm)
                     doy, & ! Day of year
                rainfall, & ! rainfall (kgH2O/m2/s)
                wind_spd, & ! wind speed (m/s)
                  vpd_pa, & ! Vapour pressure deficit (Pa)
                     lai    ! leaf area index (m2/m2)

! Module level varoables for step specific timing information
double precision :: seconds_per_step, & !
                       days_per_step, & !
                     days_per_step_1, & !
                        dayl_seconds, & ! day length in seconds
                          dayl_hours    ! day length in hours

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days
                                                meant_time, &
                                           daylength_hours, &
                                         daylength_seconds, &
                                             rainfall_time, &
                                       co2_half_saturation, & ! (ppm)
                                    co2_compensation_point    ! (ppm)

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai_out,NEE_out,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP_out)

    ! The Data Assimilation Linked Ecosystem Carbon - Growing Season
    ! Index using DeltaGSI to determine FOLiar phenology- Forest Rotation (DALEC_GSI_DFOL_FR) model.
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP and
    ! partitions between various ecosystem carbon pools. These pools are
    ! subject to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2

    implicit none

    ! declare input variables
    integer, intent(in) :: start &
                          ,finish &
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat(nodays)    & ! time step in decimal days
                         ,pars(nopars)      & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout) :: lai_out & ! leaf area index
                                               ,GPP_out & ! Gross primary productivity
                                               ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    double precision :: Rtot    ! Total hydraulic resistance (MPa.s-1.m-2.mmol-1)
    integer :: p,f,nxp,n,test,m

    ! local fire related variables
    double precision :: burnt_area &
                       ,CFF(7) = dble_zero, CFF_res(4) = dble_zero    & ! combusted and non-combustion fluxes
                       ,NCFF(7) = dble_zero, NCFF_res(4) = dble_zero  & ! with residue and non-residue seperates
                       ,combust_eff(5)                                & ! combustion efficiency
                       ,rfac                                            ! resilience factor

    ! local deforestation related variables
    double precision, dimension(4) :: post_harvest_burn   & ! how much burning to occur after
                                     ,foliage_frac_res    &
                                     ,roots_frac_res      &
                                     ,rootcr_frac_res     &
                                     ,stem_frac_res       &
                                     ,branch_frac_res     &
                                     ,Cbranch_part        &
                                     ,Crootcr_part        &
                                     ,soil_loss_frac

    double precision :: labile_loss,foliar_loss      &
                       ,roots_loss,wood_loss         &
                       ,labile_residue,foliar_residue&
                       ,roots_residue,wood_residue   &
                       ,wood_pellets,C_total         &
                       ,labile_frac_res              &
                       ,Cstem,Cbranch,Crootcr        &
                       ,stem_residue,branch_residue  &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots

    integer :: reforest_day, harvest_management, restocking_lag, gsi_lag

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th precipitation (kgH2O.m-2.s-1)
    ! 8th deforestation fraction
    ! 9th burnt area fraction
    ! 10th 21 day average min temperature (oC)
    ! 11th 21 day average photoperiod (seconds)
    ! 12th 21 day average VPD (kPa)
    ! 13th Forest management practice to accompany any clearing
    ! 14th avg daily temperature (oC)
    ! 15th avg daily wind speed (m.s-1)
    ! 16th vapour pressure deficit (Pa)

    ! POOLS are:
    ! 1 = labile (p18)
    ! 2 = foliar (p19)
    ! 3 = root   (p20)
    ! 4 = wood   (p21)
    ! 5 = litter (p22)
    ! 6 = som    (p23)
    ! 7 = cwd    (p37)
    ! 8 = soil water content (currently assumed to field capacity)

    ! p(30) = labile replanting
    ! p(31) = foliar replanting
    ! p(32) = fine root replanting
    ! p(33) = wood replanting

    ! FLUXES are:
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile production
    ! 9 = leaffall factor
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som
    ! 16 = labrelease factor
    ! 17 = carbon flux due to fire
    ! 18 = growing season index
    ! 19 = CWD turnover to litter
    ! 20 = C extracted as harvest
    ! 21 = labile loss due to disturbance
    ! 22 = foliage loss due to disturbance
    ! 23 = root loss due to disturbance
    ! 24 = wood loss due to disturbance

    ! PARAMETERS
    ! 31 process parameters; 7 C pool initial conditions

    ! p(1) = Litter to SOM conversion rate (fraction)
    ! p(2) = RmGPP fraction
    ! p(3) = GSI sensitivity for leaf growth
    ! p(4) = Max labile turnover to roots (fraction)
    ! p(5) = Max leaf turnover (GSI; fraction)
    ! p(6) = Turnover rate of wood (fraction)
    ! p(7) = Turnover rate of roots (fraction)
    ! p(8) = Litter turnover rate to heterotrophic respiration (fraction)
    ! p(9) = SOM turnover rate to heterotrophic respiration (fraction)
    ! p(10) = Exponential coefficient for temperature response for heterotrophic respiration
    ! p(11) = Average foliar nitrogen content (log10(gN/m2))
    ! p(12) = Max labile turnover to leaves (GSI; fraction)
    ! p(13) = Max labile turnover to wood (fraction)
    ! p(14) = Min temp threshold (GSI; Kelvin)
    ! p(15) = Max temp threshold (GSI; Kelvin)
    ! p(16) = Min photoperiod threshold (GSI; seconds)
    ! p(17) = Leaf Mass per unit Area (gC.m-2)
    ! p(18) = Initial labile pool (gC/m2)
    ! p(19) = Initial foliage pool (gC/m2)
    ! p(20) = Initial root pool (gC/m2)
    ! p(21) = Initial wood pool (gC/m2)
    ! p(22) = Initial litter pool (gC/m2)
    ! p(23) = Initial som pool (gC/m2)
    ! p(24) = Max photoperiod threshold (GSI; seconds)
    ! p(25) = Min vapour pressure deficit threshold (GSI; Pa)
    ! p(26) = Max vapour pressure deficit (GSI; Pa)
    ! p(27) = Fraction return threshold on GPP for LAI growth
    ! p(28) = Fraction of Cwood which is Cbranch
    ! p(29) = Fraction of Cwood which is Ccoarseroot
    ! p(37) = Initial CWD pool (gC/m2)
    ! p(38) = CWD turnover fraction (fraction)

    ! variables related to deforestation
    ! labile_loss = total loss from labile pool from deforestation
    ! foliar_loss = total loss form foliar pool from deforestation
    ! roots_loss = total loss from root pool from deforestation
    ! wood_loss = total loss from wood pool from deforestation
    ! labile_residue = harvested labile remaining in system as residue
    ! foliar_residue = harested foliar remaining in system as residue
    ! roots_residue = harvested roots remaining in system as residue
    ! wood_residue = harvested wood remaining in system as residue
    ! coarse_root_residue = expected coarse woody root left in system as residue

    ! parameters related to deforestation
    ! labile_frac_res = fraction of labile harvest left as residue
    ! foliage_frac_res = fraction of foliage harvest left as residue
    ! roots_frac_res = fraction of roots harvest left as residue
    ! wood_frac_res = fraction of wood harvest left as residue
    ! Crootcr_part = fraction of wood pool expected to be coarse root
    ! Crootcr_frac_res = fraction of coarse root left as residue
    ! soil_loss_frac = fraction determining Csom expected to be physically
    ! removed along with coarse roots

! profiling example
!real :: begin, done,f1=0,f2=0,f3=0,f4=0,f5=0
!real :: Rtot_times=0, aero_time=0 , soilwater_time=0 , acm_et_time = 0
!call cpu_time(begin)
!call cpu_time(finish)

    ! load ACM-GPP-ET parameters
    NUE                       = 1.850535d+01  ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                              ! ,unlimited by CO2, light and
                                              ! photoperiod
                                              ! (gC/gN/m2leaf/day)
    pn_max_temp               = 6.982614d+01  ! Maximum temperature for photosynthesis (oC)
    pn_opt_temp               = 3.798068d+01  ! Optimum temperature for photosynthesis (oC)
    pn_kurtosis               = 1.723531d-01  ! Kurtosis of photosynthesis temperature response
    e0                        = 4.489652d+00  ! Quantum yield gC/MJ/m2/day PAR
    max_lai_lwrad_absorption  = 9.282892d-01  ! Max fraction of LW from sky absorbed by canopy
    lai_half_lwrad_absorption = 5.941333d-01  ! LAI at which canopy LW absorption = 50 %
    max_lai_nir_absorption    = 8.333743d-01  ! Max fraction of NIR absorbed by canopy
    lai_half_nir_absorption   = 2.148633d+00  ! LAI at which canopy NIR absorption = 50 %
    minlwp                    = -1.990154d+00 ! minimum leaf water potential (MPa)
    max_lai_par_absorption    = 8.737539d-01  ! Max fraction of PAR absorbed by canopy
    lai_half_par_absorption   = 1.804925d+00  ! LAI at which canopy PAR absorption = 50 %
    lai_half_lwrad_to_sky     = 2.489314d+00  ! LAI at which 50 % LW is reflected back to sky
    iWUE                      = 1.722579d-02  ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
    soil_swrad_absorption     = 7.375071d-01  ! Fraction of SW rad absorbed by soil
    max_lai_swrad_reflected   = 2.796492d-01  ! Max fraction of SW reflected back to sky
    lai_half_swrad_reflected  = (lai_half_nir_absorption+lai_half_par_absorption) * 0.5d0
    max_lai_lwrad_release     = 2.481599d-01  ! Max fraction of LW emitted from canopy to be released
    lai_half_lwrad_release    = 5.020443d-01  ! LAI at which LW emitted from canopy to be released at 50 %
    soilevap_rad_intercept    = 1.122969d-02  ! Intercept (kgH2O/m2/day) on linear adjustment to soil evaporation
                                              ! to account for non-calculation
                                              ! of energy balance
    soilevap_rad_coef         = 1.748044d+00  ! Coefficient on linear adjustment to
                                              ! soil evaporation to account for
                                              ! non-calculation of energy
                                              ! balance
    ! load some values
    avN = 10d0**pars(11)  ! foliar N gN/m2

    ! initial values for deforestation variables
    labile_loss = dble_zero    ; foliar_loss = dble_zero
    roots_loss = dble_zero     ; wood_loss = dble_zero
    labile_residue = dble_zero ; foliar_residue = dble_zero
    roots_residue = dble_zero  ; wood_residue = dble_zero
    stem_residue = dble_zero   ; branch_residue = dble_zero
    reforest_day = 0
    soil_loss_with_roots = dble_zero
    coarse_root_residue = dble_zero
    post_harvest_burn = dble_zero

    ! now load the hardcoded forest management parameters into their locations

    ! Parameter values for deforestation variables
    ! scenario 1
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = dble_one
    roots_frac_res(1)   = dble_one
    rootcr_frac_res(1) = dble_one
    branch_frac_res(1) = dble_one
    stem_frac_res(1)   = dble_zero !
    ! wood partitioning (fraction)
    Crootcr_part(1) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(1) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(1) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(1) = dble_one

    !## scen 2
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = dble_one
    roots_frac_res(2)   = dble_one
    rootcr_frac_res(2) = dble_one
    branch_frac_res(2) = dble_one
    stem_frac_res(2)   = dble_zero !
    ! wood partitioning (fraction)
    Crootcr_part(2) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(2) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(2) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(2) = dble_zero

    !## scen 3
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5d0
    roots_frac_res(3)   = dble_one
    rootcr_frac_res(3) = dble_one
    branch_frac_res(3) = dble_zero
    stem_frac_res(3)   = dble_zero !
    ! wood partitioning (fraction)
    Crootcr_part(3) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(3) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(3) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(3) = dble_zero

    !## scen 4
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5d0
    roots_frac_res(4)   = dble_one
    rootcr_frac_res(4) = dble_zero
    branch_frac_res(4) = dble_zero
    stem_frac_res(4)   = dble_zero
    ! wood partitioning (fraction)
    Crootcr_part(4) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(4) =  0.20d0 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(4) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(4) = dble_zero

    ! for the moment override all paritioning parameters with those coming from
    ! CARDAMOM
    Cbranch_part = pars(28)
    Crootcr_part = pars(29)

    ! declare fire constants (labile, foliar, roots, wood, litter)
    combust_eff(1) = 0.1d0 ; combust_eff(2) = 0.9d0
    combust_eff(3) = 0.1d0 ; combust_eff(4) = 0.5d0
    combust_eff(5) = 0.3d0 ; rfac = 0.5d0

    ! in preparation for EnKF development DALEC is begin made to run 1 time step
    ! at a time without resetting all variables
    if (start == 1) then

        ! assigning initial conditions
        POOLS(1,1)=pars(18)
        POOLS(1,2)=pars(19)
        POOLS(1,3)=pars(20)
        POOLS(1,4)=pars(21)
        POOLS(1,5)=pars(22)
        POOLS(1,6)=pars(23)
        POOLS(1,7)=pars(37)

        if (.not.allocated(disturbance_residue_to_som)) then
            allocate(disturbance_residue_to_litter(nodays), &
                     disturbance_residue_to_som(nodays),    &
                     disturbance_residue_to_cwd(nodays),    &
                     disturbance_loss_from_litter(nodays),  &
                     disturbance_loss_from_som(nodays),     &
                     disturbance_loss_from_cwd(nodays))
        endif
        disturbance_residue_to_litter = dble_zero ; disturbance_loss_from_litter = dble_zero
        disturbance_residue_to_som = dble_zero ; disturbance_loss_from_som = dble_zero
        disturbance_residue_to_cwd = dble_zero ; disturbance_loss_from_cwd = dble_zero

        ! calculate some values once as these are invarient between DALEC runs
        if (.not.allocated(tmp_x)) then
            ! 21 days is the maximum potential so we will fill the maximum potential
            ! + 1 for safety
            allocate(tmp_x(22),tmp_m(nodays))
            do f = 1, 22
               tmp_x(f) = f
            end do
            do n = 1, nodays
              ! calculate the gradient / trend of GSI
              if (sum(deltat(1:n)) < 21d0) then
                  tmp_m(n) = n-1
              else
                 ! else we will try and work out the gradient to see what is
                 ! happening
                 ! to the system over all. The default assumption will be to
                 ! consider
                 ! the averaging period of GSI model (i.e. 21 days). If this is not
                 ! possible either the time step of the system is used (if step
                 ! greater
                 ! than 21 days) or all available steps (if n < 21).
                 m = 0 ; test = 0
                 do while (test < 21d0)
                    m = m + 1 ; test = nint(sum(deltat((n-m):n)))
                    if (m > (n-1)) test = 21
                 end do
                 tmp_m(n) = m
              endif ! for calculating gradient
            end do ! calc daily values once
            ! allocate GSI history dimension
            gsi_lag_remembered=nint(max(2d0,maxval(tmp_m)))
        end if ! .not.allocated(tmp_x)
        ! assign our starting value
        gsi_history = pars(36)-dble_one
         just_grown = pars(35)

        ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
        if (.not.allocated(deltat_1)) then
           allocate(deltat_1(nodays), &
                    daylength_hours(nodays),daylength_seconds(nodays), &
                    meant_time(nodays),rainfall_time(nodays), &
                    co2_half_saturation(nodays),co2_compensation_point(nodays))
           ! inverse of time step (days-1) to avoid divisions
           deltat_1 = deltat**(-dble_one)
           ! meant time step temperature
           meant_time = (met(2,1:nodays)+met(3,1:nodays)) * 0.5d0
           do n = 1, nodays
              ! calculate day length as invarient between iterations
              daylength_hours(n)=daylength_in_hours(met(6,n)-(deltat(n)*0.5d0),lat)
              ! check positive values only for rainfall input
              rainfall_time(n)=max(dble_zero,met(7,n))
              ! Temperature adjustments for Michaelis-Menten coefficients
              ! for CO2 (kc) and O2 (ko) and CO2 compensation point.
              co2_compensation_point(n) = arrhenious(co2comp_saturation,co2comp_half_sat_conc,met(3,n))
              co2_half_saturation(n) = arrhenious(kc_saturation,kc_half_sat_conc,met(3,n))
           end do
           ! generate daylength per seconds
           daylength_seconds = daylength_hours * seconds_per_hour

        else

        endif ! has SWP already been determined?

    endif ! start == 1

    ! assign climate sensitivities
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory
    fol_turn_crit=pars(34)-1d0
    lab_turn_crit=pars(3)-1d0

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = met(6,n)   ! Day of year
      rainfall = rainfall_time(n) ! rainfall (kgH2O/m2/s)
      meant = meant_time(n)  ! mean air temperature (oC)
      meant_K = meant + freeze
      wind_spd = met(15,n) ! wind speed (m/s)
      vpd_pa = met(16,n)  ! Vapour pressure deficit (Pa)

      ! states needed for module variables
      lai_out(n) = POOLS(n,2)/pars(17)
      lai = lai_out(n) ! leaf area index (m2/m2)

      ! Temperature adjustments for Michaelis-Menten coefficients
      ! for CO2 (kc) and O2 (ko) and CO2 compensation point
      ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
      co2_half_sat   = co2_half_saturation(n)
      co2_comp_point = co2_compensation_point(n)

      ! extract timing related values
      dayl_hours = daylength_hours(n)
      dayl_seconds = daylength_seconds(n)
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      !!!!!!!!!!
      ! calculate soil water potential and total hydraulic resistance
      !!!!!!!!!!

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,POOLS(n,3)*2d0)
      call calculate_Rtot(Rtot)

      !!!!!!!!!!
      ! Calculate surface exchange coefficients
      !!!!!!!!!!

      ! calculate aerodynamic using consistent approach with SPA
      call calculate_aerodynamic_conductance
      ! calculate variables used commonly between ACM_GPP and ACM_ET
      call acm_albedo_gc(abs(minlwp),Rtot)

      !!!!!!!!!!
      ! GPP (gC.m-2.day-1)
      !!!!!!!!!!

      if (stomatal_conductance > dble_zero) then
          FLUXES(n,1) = max(dble_zero,acm_gpp(stomatal_conductance))
      else
          FLUXES(n,1) = dble_zero
      endif

      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(10)*meant)
      ! (maintenance) autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = dble_zero !(FLUXES(n,1)-FLUXES(n,3))*pars(3)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))*pars(13)
      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5))*pars(4)
      ! wood production
      FLUXES(n,7) = FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5)-FLUXES(n,6)

      ! GSI added to fortran version by TLS 24/11/2014
      ! /* 25/09/14 - JFE
      ! Here we calculate the Growing Season Index based on
      ! Jolly et al. A generalized, bioclimatic index to predict foliar
      ! phenology in response to climate Global Change Biology, Volume 11, page 619-632,
      ! 2005 (doi: 10.1111/j.1365-2486.2005.00930.x)
      ! Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S.
      ! Denning (2011), A global reanalysis of vegetation phenology, J. Geophys. Res.,
      ! 116, G03020, doi:10.1029/2010JG001545.

      ! It is the product of 3 limiting factors for temperature, photoperiod and
      ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated
      ! min and max value. Photoperiod, VPD and avgTmin are direct input

      ! temperature limitation, then restrict to 0-1; correction for k-> oC
      Tfac = (met(10,n)-(pars(14)-freeze)) / (pars(15)-pars(14))
      Tfac = min(dble_one,max(dble_zero,Tfac))
      ! photoperiod limitation
      Photofac = (met(11,n)-pars(16)) / (pars(24)-pars(16))
      Photofac = min(dble_one,max(dble_zero,Photofac))
      ! VPD limitation
      VPDfac = 1.0 - ( (met(12,n)-pars(25)) / (pars(26)-pars(25)) )
      VPDfac = min(dble_one,max(dble_zero,VPDfac))

      ! calculate and store the GSI index
      FLUXES(n,18) = Tfac*Photofac*VPDfac

      ! we will load up some needed variables
      m = tmp_m(n)
      ! update gsi_history for the calculation
      if (n == 1) then
          ! in first step only we want to take the initial GSI value only
          gsi_history(gsi_lag) = FLUXES(n,18)
      else
          gsi_history((gsi_lag-m):gsi_lag) = FLUXES((n-m):n,18)
      endif
      ! calculate gradient
      gradient = linear_model_gradient(tmp_x(1:(gsi_lag)),gsi_history(1:gsi_lag),gsi_lag)
      ! adjust gradient to daily rate
      gradient = gradient /  nint((sum(deltat((n-m+1):n))) / (gsi_lag-1))
      gsi_lag_remembered = gsi_lag

      ! first assume that nothing is happening
      FLUXES(n,9) = dble_zero  ! leaf turnover
      FLUXES(n,16) = dble_zero ! leaf growth

      ! now update foliage and labile conditions based on gradient calculations
      if (gradient < fol_turn_crit .or. FLUXES(n,18) == dble_zero) then
         ! we are in a decending condition so foliar turnover
         FLUXES(n,9) = pars(5)*(dble_one-FLUXES(n,18))
         just_grown = 0.5d0
      else if (gradient > lab_turn_crit) then
         ! we are in a assending condition so labile turnover
         FLUXES(n,16) = pars(12)*FLUXES(n,18)
         just_grown = 1.5d0
         ! check carbon return
         tmp = POOLS(n,1)*min(dble_one,dble_one-(dble_one-FLUXES(n,16))**deltat(n))/deltat(n)
         tmp = (POOLS(n,2)+tmp)/pars(17)
         tmp_lai = lai ; lai = tmp
         tmp = max(dble_zero,acm_gpp(stomatal_conductance))
         lai = tmp_lai
         ! determine if increase in LAI leads to an improvement in GPP greater
         ! than
         ! critical value, if not then no labile turnover allowed
         if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(27) ) then
             FLUXES(n,16) = dble_zero
         endif
      else
         ! probably we want nothing to happen, however if we are at the seasonal
         ! maximum we will consider further growth still
         if (just_grown >= dble_one) then
            ! we are between so definitely not losing foliage and we have
            ! previously been growing so maybe we still have a marginal return on
            ! doing so again
            FLUXES(n,16) = pars(12)*FLUXES(n,18)
            ! but possibly gaining some?
            ! determine if this is a good idea based on GPP increment
            tmp = POOLS(n,1)*min(dble_one,dble_one-(dble_one-FLUXES(n,16))**deltat(n))/deltat(n)
            tmp = (POOLS(n,2)+tmp)/pars(17)
            tmp_lai = lai ; lai = tmp
            tmp = max(dble_zero,acm_gpp(stomatal_conductance))
            lai = tmp_lai            ! determine if increase in LAI leads to an improvement in GPP greater
            ! than
            ! critical value, if not then no labile turnover allowed
            if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(27) ) then
                FLUXES(n,16) = dble_zero
            endif
         end if ! Just grown?
      endif ! gradient choice

      ! these allocated if post-processing
      if (allocated(itemp)) then
         itemp(n) = Tfac
         ivpd(n) = VPDfac
         iphoto(n) = Photofac
      endif

      !
      ! those with time dependancies
      !

      ! total labile release
      FLUXES(n,8)  = POOLS(n,1)*min(dble_one,dble_one-(dble_one-FLUXES(n,16))**deltat(n))/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*min(dble_one,dble_one-(dble_one-FLUXES(n,9))**deltat(n))/deltat(n)
      ! total wood litter production
      FLUXES(n,11) = POOLS(n,4)*min(dble_one,dble_one-(dble_one-pars(6))**deltat(n))/deltat(n)
      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*min(dble_one,dble_one-(dble_one-pars(7))**deltat(n))/deltat(n)

      !
      ! those with temperature AND time dependancies
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(8)))**deltat(n))/deltat(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(9)))**deltat(n))/deltat(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(1)))**deltat(n))/deltat(n)
      ! CWD to litter
      FLUXES(n,19) = POOLS(n,7)*min(dble_one,dble_one-(dble_one-(FLUXES(n,2)*pars(38)))**deltat(n))/deltat(n)

      ! calculate growth respiration and adjust allocation to pools assuming
      ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
      ! foliage
      Rg_from_labile = FLUXES(n,8)*0.21875d0 ; FLUXES(n,8) = FLUXES(n,8) * 0.78125d0
      ! now update the Ra flux
      FLUXES(n,3) = FLUXES(n,3) + Rg_from_labile
      ! roots
      FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,6)*0.21875d0) ; FLUXES(n,6) = FLUXES(n,6) * 0.78125d0
      ! wood
      FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,7)*0.21875d0) ; FLUXES(n,7) = FLUXES(n,7) * 0.78125d0

      ! calculate the NEE
      NEE_out(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! load GPP
      GPP_out(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      !

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8)-Rg_from_labile)*deltat(n)
      ! foliar pool
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)-FLUXES(n,10) + FLUXES(n,8))*deltat(n)
      ! wood pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
      ! root pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)+FLUXES(n,19)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14))*deltat(n)
      ! cwd pool
      POOLS(n+1,7) = POOLS(n,7) + (FLUXES(n,11)-FLUXES(n,19))*deltat(n)

      !
      ! deal first with deforestation
      !

      if (n == reforest_day) then
          POOLS(n+1,1) = pars(30)
          POOLS(n+1,2) = pars(31)
          POOLS(n+1,3) = pars(32)
          POOLS(n+1,4) = pars(33)
      endif

      ! reset variables
      FLUXES(n,21:24) = dble_zero ;  harvest_management = 0 ; burnt_area = dble_zero

      if (met(8,n) > dble_zero) then

          ! pass harvest management to local integer
          harvest_management = int(met(13,n))

          ! assume that labile is proportionally distributed through the plant
          ! and therefore so is the residual fraction
          C_total = POOLS(n+1,2) + POOLS(n+1,3) + POOLS(n+1,4)
          ! partition wood into its components
          Cbranch = POOLS(n+1,4)*Cbranch_part(harvest_management)
          Crootcr = POOLS(n+1,4)*Crootcr_part(harvest_management)
          Cstem   = POOLS(n+1,4)-(Cbranch + Crootcr)
          ! now calculate the labile fraction of residue
          labile_frac_res = ( (POOLS(n+1,2)/C_total) * foliage_frac_res(harvest_management) ) &
                          + ( (POOLS(n+1,3)/C_total) * roots_frac_res(harvest_management)   ) &
                          + ( (Cbranch/C_total)      * branch_frac_res(harvest_management)  ) &
                          + ( (Cstem/C_total)        * stem_frac_res(harvest_management)    ) &
                          + ( (Crootcr/C_total)      * rootcr_frac_res(harvest_management)  )

          ! loss of carbon from each pools
          labile_loss = POOLS(n+1,1)*met(8,n)
          foliar_loss = POOLS(n+1,2)*met(8,n)
          roots_loss  = POOLS(n+1,3)*met(8,n)
          wood_loss   = POOLS(n+1,4)*met(8,n)

          ! for output / EDC updates
          if (met(8,n) <= 0.99d0) then
              FLUXES(n,21) = labile_loss / deltat(n)
              FLUXES(n,22) = foliar_loss / deltat(n)
              FLUXES(n,23) = roots_loss / deltat(n)
              FLUXES(n,24) = wood_loss / deltat(n)
          endif

          ! transfer fraction of harvest waste to litter or som pools
          ! easy pools first
          labile_residue = POOLS(n+1,1)*met(8,n)*labile_frac_res
          foliar_residue = POOLS(n+1,2)*met(8,n)*foliage_frac_res(harvest_management)
          roots_residue  = POOLS(n+1,3)*met(8,n)*roots_frac_res(harvest_management)
          ! explicit calculation of the residues from each fraction
          coarse_root_residue  = Crootcr*met(8,n)*rootcr_frac_res(harvest_management)
          branch_residue = Cbranch*met(8,n)*branch_frac_res(harvest_management)
          stem_residue = Cstem*met(8,n)*stem_frac_res(harvest_management)
          ! now finally calculate the final wood residue
          wood_residue = stem_residue + branch_residue + coarse_root_residue
          ! mechanical loss of Csom due to coarse root extraction
          soil_loss_with_roots = Crootcr*met(8,n)*(dble_one-rootcr_frac_res(harvest_management)) &
                              * soil_loss_frac(harvest_management)

          ! update living pools directly
          POOLS(n+1,1) = max(dble_zero,POOLS(n+1,1)-labile_loss)
          POOLS(n+1,2) = max(dble_zero,POOLS(n+1,2)-foliar_loss)
          POOLS(n+1,3) = max(dble_zero,POOLS(n+1,3)-roots_loss)
          POOLS(n+1,4) = max(dble_zero,POOLS(n+1,4)-wood_loss)

          ! ensure fire related values are reset
          FLUXES(n,17) = dble_zero
          CFF = dble_zero ; NCFF = dble_zero
          CFF_res = dble_zero ; NCFF_res = dble_zero

          ! update all pools this time
          POOLS(n+1,1) = max(dble_zero, POOLS(n+1,1) - CFF(1) - NCFF(1) )
          POOLS(n+1,2) = max(dble_zero, POOLS(n+1,2) - CFF(2) - NCFF(2) )
          POOLS(n+1,3) = max(dble_zero, POOLS(n+1,3) - CFF(3) - NCFF(3) )
          POOLS(n+1,4) = max(dble_zero, POOLS(n+1,4) - CFF(4) - NCFF(4) )
          POOLS(n+1,5) = max(dble_zero, POOLS(n+1,5) + (labile_residue+foliar_residue+roots_residue) &
                                              + (NCFF(1)+NCFF(2)+NCFF(3))-CFF(5)-NCFF(5) )
          POOLS(n+1,6) = max(dble_zero, POOLS(n+1,6) - soil_loss_with_roots + (NCFF(4)+NCFF(5)+NCFF(7)))
          POOLS(n+1,7) = max(dble_zero, POOLS(n+1,7) + wood_residue - CFF(7) - NCFF(7) )
          ! some variable needed for the EDCs
          ! reallocation fluxes for the residues
          disturbance_residue_to_litter(n) = (labile_residue+foliar_residue+roots_residue) &
                                           + (NCFF(1)+NCFF(2)+NCFF(3))
          disturbance_loss_from_litter(n)  = CFF(5)+NCFF(5)
          disturbance_residue_to_cwd(n)    = wood_residue
          disturbance_loss_from_cwd(n)     = CFF(7) - NCFF(7)
          disturbance_residue_to_som(n)    = NCFF(4)+NCFF(5)+NCFF(7)
          disturbance_loss_from_som(n)     = soil_loss_with_roots
          ! convert all to rates to be consistent with the FLUXES in EDCs
          disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n) / deltat(n)
          disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) / deltat(n)
          disturbance_residue_to_cwd(n)    = disturbance_residue_to_cwd(n) / deltat(n)
          disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) / deltat(n)
          disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) / deltat(n)
          disturbance_loss_from_som(n)     = disturbance_loss_from_som(n) / deltat(n)
          ! harvested carbon from all pools
          FLUXES(n,20) = (wood_loss-(wood_residue+CFF_res(4)+NCFF_res(4))) &
                       + (labile_loss-(labile_residue+CFF_res(1)+NCFF_res(1))) &
                       + (foliar_loss-(foliar_residue+CFF_res(2)+NCFF_res(2))) &
                       + (roots_loss-(roots_residue+CFF_res(3)+NCFF_res(3)))
          ! convert to daily rate
          FLUXES(n,20) = FLUXES(n,20) / deltat(n)
          ! total carbon loss from the system
          C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)) &
                  - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF))

          ! if total clearance occured then we need to ensure some minimum
          ! values and reforestation is assumed one year forward
          if (met(8,n) > 0.99d0) then
              m = 0 ; test = nint(sum(deltat(n:(n+m))))
              ! FC Forest Statistics 2015 lag between harvest and restocking ~ 2 year
              restocking_lag = 365*2
              do while (test < restocking_lag)
                 m = m+1 ; test = nint(sum(deltat(n:(n+m))))
                 !  get out clause for hitting the end of the simulation
                 if (m+n >= nodays) test = restocking_lag
              enddo
              reforest_day = min((n+m), nodays)
          endif ! if total clearance

      endif ! end deforestation info

      !
      ! then deal with fire
      !

      if (met(9,n) > dble_zero .or.(met(8,n) > dble_zero .and. harvest_management > 0)) then

          burnt_area = met(9,n)
          if (met(8,n) > dble_zero .and. burnt_area > dble_zero) then
              ! pass harvest management to local integer
              burnt_area = min(dble_one,burnt_area + post_harvest_burn(harvest_management))
          else if (met(8,n) > dble_zero .and. burnt_area <= dble_zero) then
              burnt_area = post_harvest_burn(harvest_management)
          endif

         if (burnt_area > dble_zero) then

             !/*first fluxes*/
             !/*LABILE*/
             CFF(1) = POOLS(n+1,1)*met(9,n)*combust_eff(1)
             NCFF(1) = POOLS(n+1,1)*met(9,n)*(dble_one-combust_eff(1))*(1-rfac)
             !/*foliar*/
             CFF(2) = POOLS(n+1,2)*met(9,n)*combust_eff(2)
             NCFF(2) = POOLS(n+1,2)*met(9,n)*(dble_one-combust_eff(2))*(1-rfac)
             !/*root*/
             CFF(3) = dble_zero ! POOLS(n+1,3)*met(9,n)*combust_eff(3)
             NCFF(3) = dble_zero ! POOLS(n+1,3)*met(9,n)*(dble_one-combust_eff(3))*(dble_one-rfac)
             !/*wood*/
             CFF(4) = POOLS(n+1,4)*met(9,n)*combust_eff(4)
             NCFF(4) = POOLS(n+1,4)*met(9,n)*(dble_one-combust_eff(4))*(dble_one-rfac)
             !/*litter*/
             CFF(5) = POOLS(n+1,5)*met(9,n)*combust_eff(5)
             NCFF(5) = POOLS(n+1,5)*met(9,n)*(dble_one-combust_eff(5))*(1-rfac)
             !/*CWD*/ Using Combustion factors for wood
             CFF(7) = POOLS(n+1,7)*met(9,n)*combust_eff(4)
             NCFF(7) = POOLS(n+1,7)*met(9,n)*(dble_one-combust_eff(4))*(1-rfac)

             !/*fires as daily averages to comply with units*/
             FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)+CFF(7)) / deltat(n)
             !/*all fluxes are at a daily timestep*/
             NEE_out(n)=NEE_out(n)+FLUXES(n,17)
             ! determine the as daily rate impact on live tissues for use in EDC and
             ! MTT calculations
             FLUXES(n,21) = FLUXES(n,21) + ((CFF(1) + NCFF(1)) / deltat(n))
             FLUXES(n,22) = FLUXES(n,22) + ((CFF(2) + NCFF(2)) / deltat(n))
             FLUXES(n,23) = FLUXES(n,23) + ((CFF(3) + NCFF(3)) / deltat(n))
             FLUXES(n,24) = FLUXES(n,24) + ((CFF(4) + NCFF(4)) / deltat(n))

             !// update pools
             !/*Adding all fire pool transfers here*/
             POOLS(n+1,1)=POOLS(n+1,1)-CFF(1)-NCFF(1)
             POOLS(n+1,2)=POOLS(n+1,2)-CFF(2)-NCFF(2)
             POOLS(n+1,3)=POOLS(n+1,3)-CFF(3)-NCFF(3)
             POOLS(n+1,4)=POOLS(n+1,4)-CFF(4)-NCFF(4)
             POOLS(n+1,5)=POOLS(n+1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3)
             POOLS(n+1,6)=POOLS(n+1,6)+NCFF(4)+NCFF(5)+NCFF(7)
             POOLS(n+1,7)=POOLS(n+1,7)-CFF(7)-NCFF(7)
             ! some variable needed for the EDCs
             ! reallocation fluxes for the residues
             disturbance_residue_to_litter(n) = (NCFF(1)+NCFF(2)+NCFF(3))
             disturbance_residue_to_som(n)    = (NCFF(4)+NCFF(5)+NCFF(7))
             disturbance_loss_from_litter(n)  = CFF(5)+NCFF(5)
             disturbance_loss_from_cwd(n)     = CFF(7) - NCFF(7)
             ! convert all to rates to be consistent with the FLUXES in EDCs
             disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n) / deltat(n)
             disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) / deltat(n)
             disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) / deltat(n)
             disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) / deltat(n)

         end if ! burnt_area  > 0

      endif ! end burnst area issues

      do nxp = 1, nopools
         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < dble_zero) then
            print*,"step",n,"POOL",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
            print*,"POOLS+1",POOLS(n+1,:)
            stop
         endif
      enddo

    end do ! nodays loop
!call cpu_time(done)
!print*,done-begin
  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp(gs)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pn, pd, pp, qq, ci, mult, pl &
                       ,gc ,gs_mol, gb_mol

    ! Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
!    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,maxt)
!    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,maxt)

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
    pn = lai*avN*NUE*opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,maxt)

    !
    ! Diffusion limited photosynthesis
    !

    ! daily canopy conductance (mmolH2O.m-2.day-1-> molCO2.m-2.day-1)
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    ! i.e. gcH2O*1.646259 = gcCO2
    gs_mol = gs * 1d-3 * seconds_per_day * gs_H2O_CO2
    ! canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2).
    gb_mol = aerodynamic_conductance * seconds_per_day * convert_ms1_mol_1 * gb_H2O_CO2
    ! Combining in series the stomatal and boundary layer conductances
    gc = (gs_mol ** (-dble_one) + gb_mol ** (-dble_one)) ** (-dble_one)

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = (pn/umol_to_gC)/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mult = co2+qq-pp
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))
    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
    pd = (gc * (co2-ci)) * umol_to_gC
    ! scale to day light period as this is then consistent with the light
    ! capture period
    pd = pd * (dayl_hours / 24d0)

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    pl = e0 * canopy_par_MJday

    !
    ! CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis
    acm_gpp = pl*pd/(pl+pd)

    ! don't forget to return
    return

  end function acm_gpp
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE (iWUE).

    ! Arguments
    double precision, intent(in) :: gs_in

    ! Local variables
    double precision :: gs_high, gpp_high, gpp_low

    ! Estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp(gs_in)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! Estimate photosynthesis with incremented gs
    gpp_high = acm_gpp(gs_high)

    ! Determine impact of gs increment on pd and how far we are from iWUE
    find_gs = iWUE - ((gpp_high - gpp_low)/lai)

  end function find_gs
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_aerodynamic_conductance

    !
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are with SPA at the time
    ! of the calibration
    !

    implicit none

    ! Local variables
    double precision :: ustar_Uh

    ! Calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh)
    ! Calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
!    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd * ustar_Uh

    ! Based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! calculate bulk conductance (Jones p68)
    aerodynamic_conductance = (canopy_wind * vonkarman_2) &
                            / (log((canopy_height-displacement)/roughl))**2

  end subroutine calculate_aerodynamic_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine acm_albedo_gc(deltaWP,Rtot)

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI

    implicit none

    ! Arguments
    double precision, intent(in) :: deltaWP, & ! minlwp-wSWP (MPa)
                                       Rtot    ! total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    ! Local variables
    double precision :: s, denom, mult, tmp
    double precision, parameter :: max_gs = 5000d0, & ! mmolH2O.m-2.s-1
                                   min_gs = 5d-12

    !!!!!!!!!!
    ! Determine some multiple use constants
    !!!!!!!!!

    ! Density of air (kg/m3)
    air_density_kg = 353d0/(maxt+freeze)
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / ((maxt+freeze)*Rcon)
    ! Latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    if (maxt < dble_one) then
        lambda = 2.835d6
    else
        lambda = 2501000d0-2364d0*maxt
    endif
    ! Psychrometric constant (kPa K-1)
    psych = (0.0646d0*exp(0.00097d0*maxt))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = maxt+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    s = 2502.935945d0*exp(17.269d0*maxt/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)

    !!!!!!!!!!
    ! Determine net shortwave and isothermal longwave energy balance
    !!!!!!!!!!

    call calculate_shortwave_balance
    call calculate_longwave_isothermal

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (deltaWP > vsmall) then
        ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
        max_supply = (deltaWP/Rtot) * dayl_seconds
    else
        ! set minimum (computer) precision level flow
        max_supply = vsmall
    end if

    ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
    ! maximum possible evaporation for the day.
    ! This will then be reduced based on CO2 limits for diffusion based
    ! photosynthesis
    denom = slope * ((canopy_swrad_MJday * 1d6 * seconds_per_day_1) + canopy_lwrad_Wm2) &
          + (air_density_kg*cpair*vpd_pa*1d-3*aerodynamic_conductance)
    denom = (denom / (lambda * max_supply * mmol_to_kg_water * seconds_per_day_1)) - slope
    denom = denom / psych
    stomatal_conductance = aerodynamic_conductance / denom

    ! Convert m.s-1 to mmolH2O.m-2.s-1
    stomatal_conductance = stomatal_conductance * 1d3 * convert_ms1_mol_1
    if (stomatal_conductance < dble_zero .or. stomatal_conductance > max_gs) stomatal_conductance = max_gs

    ! Solve for photosynthesis limits on gs through iterative solution
    delta_gs = 1d-3*dayl_seconds ! mmolH2O/m2/dayl
    stomatal_conductance = zbrent('acm_albedo_gc:find_gs',find_gs,min_gs,stomatal_conductance,delta_gs)

  end subroutine acm_albedo_gc
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

    ! Local variables
    double precision :: lwrad, & ! Downward long wave radiation from sky (W.m-2)
             longwave_release, & ! Emission of long wave radiation from surfaces per m2
                                 ! assuming isothermal condition (W.m-2)
     canopy_absorbed_fraction, & ! Fraction of incoming longwave absorbed by canopy
        sky_returned_fraction, & ! Fraction of incoming longwave reflected back into sky
      canopy_release_fraction, & ! Fraction of longwave emitted from within the canopy to ultimately be released
   canopy_absorption_from_sky, & ! Canopy absorbed radiation from downward LW (W.m-2)
  canopy_absorption_from_soil, & ! Canopy absorbed radiation from soil surface (W.m-2)
                  canopy_loss, & ! Longwave radiation released from canopy surface (W.m-2).
                                 ! i.e. this value is released from the top and the bottom
     soil_absorption_from_sky, & ! Soil absorbed radiation from sky (W.m-2)
  soil_absorption_from_canopy    ! Soil absorbed radiation emitted from canopy (W.m-2)

    ! Estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (maxt+freeze-20d0) ** 4
    ! Estimate isothermal long wave emission per unit area
    longwave_release = emiss_boltz * (maxt+freeze) ** 4

    !!!!!!!!!!
    ! Determine fraction of longwave absorbed by canopy and returned to the sky
    !!!!!!!!!!

    ! Calculate potential canopy absorbed fraction of longwave radiation
    ! incoming from outside of the canopy
    ! as a function of LAI
    canopy_absorbed_fraction = max_lai_lwrad_absorption*lai/(lai+lai_half_lwrad_absorption)
    ! Calculate the fraction of longwave radiation from sky which is reflected
    ! back into the sky
    sky_returned_fraction = (dble_one-emissivity) - &
                            (((dble_one-emissivity)*lai) / (lai+lai_half_lwrad_to_sky))
    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky
    canopy_release_fraction = max_lai_lwrad_release - &
                              ((max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release))

    !!!!!!!!!!
    ! Calculate fluxes for canopy long wave absorption
    !!!!!!!!!!

    ! Long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad * canopy_absorbed_fraction
    ! Long wave absorbed by the canopy from soil
    canopy_absorption_from_soil = longwave_release * canopy_absorbed_fraction
    ! Calculate two-sided long wave radiation emitted from canopy which is
    ! ultimately lost from to soil or sky (i.e. this value is used twice, once
    ! to soil once to sky)
    canopy_loss = 2d0 * longwave_release * lai * canopy_release_fraction

    !!!!!!!!!!
    ! Calculate fluxes for soil long wave absorption
    !!!!!!!!!!

    ! Long wave absorbed by soil from the sky.
    ! Accounting for fraction absorbed by the canopy and returned to the sky.
    ! We assume that long wave absorption is equivalent the the emissivity of
    ! the surface
    soil_absorption_from_sky = dble_one - (canopy_absorbed_fraction + sky_returned_fraction)
    soil_absorption_from_sky = soil_absorption_from_sky * lwrad * emissivity
    ! Calculate longwave absorbed by soil which si released by the canopy itself
    soil_absorption_from_canopy = canopy_loss * emissivity

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! Determine isothermal net canopy. Note two canopy_loss used to account for
    ! upwards and downwards emissions
    canopy_lwrad_Wm2 = (canopy_absorption_from_sky + canopy_absorption_from_soil) - (canopy_loss + canopy_loss)
    ! Determine isothermal net soil
    soil_lwrad_Wm2 = (soil_absorption_from_sky + soil_absorption_from_canopy) - longwave_release

  end subroutine calculate_longwave_isothermal
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation (MJ/m2/day).
    ! Radiation absorption is paritioned into NIR and PAR for canopy, and NIR +
    ! PAR for soil.

    ! SPA uses a complex multi-layer radiative transfer scheme including
    ! reflectance, transmittance any absorption. However, for a given
    ! canopy vertical profiles, the LAI absorption relationship is readily
    ! predicted via Michaelis-Menten or non-rectangular hyperbola as done here.

    implicit none

    ! Parameters
    double precision, parameter :: sw_diffuse_fraction = 0.5d0

    ! Local variables
    double precision :: absorbed_nir_fraction & !
                       ,absorbed_par_fraction & !
                       ,sky_absorped_fraction   !

    !!!!!!!!!!
    ! Determine canopy absorption / reflectance as function of LAI
    !!!!!!!!!!

    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = max_lai_nir_absorption*lai/(lai+lai_half_nir_absorption)
    absorbed_par_fraction = max_lai_par_absorption*lai/(lai+lai_half_par_absorption)
    ! Canopy reflectance of SW radiation back into the sky
    sky_absorped_fraction = max_lai_swrad_reflected*lai/(lai+lai_half_swrad_reflected)

    !!!!!!!!!!
    ! Determine SW balance
    !!!!!!!!!!

    ! now determine shortwave radiation absorbed by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = (sw_diffuse_fraction * swrad * absorbed_par_fraction)
    canopy_swrad_MJday = canopy_par_MJday + (sw_diffuse_fraction * swrad * absorbed_nir_fraction)
    ! absorption of shortwave radiation by soil (MJ.m-2.day-1)
    soil_swrad_MJday = dble_one - (sw_diffuse_fraction*absorbed_nir_fraction + &
                                   sw_diffuse_fraction*absorbed_par_fraction + sky_absorped_fraction)
    soil_swrad_MJday = swrad * soil_swrad_MJday * soil_swrad_absorption

  end subroutine calculate_shortwave_balance
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

    ! Local parameters
    double precision, parameter :: min_wind = 0.01d0 ! minimum wind speed at canopy top

    ! Log law decay
    canopy_wind = (ustar / vonkarman) * log(displacement / roughl)

    ! Set minimum value for wind speed at canopy top (m.s-1)
    canopy_wind = max(min_wind,canopy_wind)

  end subroutine log_law_decay
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot(Rtot)

    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! Declare inputs
    double precision,intent(inout) :: Rtot ! MPa.s-1.m-2.mmol-1

    ! Local variables
    integer :: i
    double precision :: slpa, cumdepth, prev, curr, sum_water_flux, &
        soilR1,soilR2,transpiration_resistance,root_reach_local
    double precision, dimension(nos_root_layers) :: root_mass    &
                                                   ,soilRT_local &
                                                   ,root_length  &
                                                   ,ratio

    ! reset water flux
    water_flux = dble_zero ; soilRT_local = dble_zero ; soilRT = dble_zero
    ratio = dble_zero ; ratio(1) = dble_one
    ! Calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! explicitly update the soil profile if there has been rooting depth
    ! changes
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = max(min_layer,root_reach-layer_thickness(1))
    layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
    ! Calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
!    transpiration_resistance = (gplant * lai)**(-dble_one)
    transpiration_resistance = canopy_height / (gplant * lai)

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 2 root layer system this can be estimates more simply

    ! Top 25 % of root profile
    slpa = (root_reach * 0.25d0) - layer_thickness(1)
    if (slpa <= dble_zero) then
        ! > 50 % of root is in top layer
        root_mass(1) = root_biomass * 0.5d0
        root_mass(1) = root_mass(1) + ((root_biomass-root_mass(1)) * (abs(slpa)/root_reach))
    else
        ! < 50 % of root is in bottom layer
        root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/(abs(slpa)+layer_thickness(1)))
    endif
    root_mass(2) = max(dble_zero,root_biomass - root_mass(1))
    root_length = root_mass * root_mass_length_coef_1
!    root_length = root_mass / (root_density * root_cross_sec_area)

    ! Soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
    root_reach_local = min(root_reach,layer_thickness(1))
    soilR2=root_resistance(root_mass(1),root_reach_local)
    soilRT_local(1) = soilR2 + transpiration_resistance
    ! Calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    ! also note that the head correction is now added rather than
    ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
    ! soilWP prior to application of minlwp
    demand = abs(minlwp)+head*canopy_height
    water_flux(1) = demand(1)/(transpiration_resistance + soilR2)
    ! Bottom root layer
    if (root_mass(2) > dble_zero ) then
        ! Soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
        soilR2=root_resistance(root_mass(2),layer_thickness(2))
        soilRT_local(2) = soilR2 + transpiration_resistance
        ! Calculate and accumulate steady state water flux in mmol.m-2.s-1
        water_flux(2) = demand(2)/(transpiration_resistance + soilR2)
        ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif ! roots present in second layer?

    ! if freezing then assume soil surface is frozen
    if (meant < dble_one) then
        water_flux(1) = dble_zero ; ratio(1) = dble_zero ; ratio(2) = dble_one
    endif
    ! Calculate sum value
    sum_water_flux = sum(water_flux)

    ! Calculate weighted SWP and uptake fraction
    soilRT = sum(soilRT_local(1:nos_root_layers) * water_flux(1:nos_root_layers))
    uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux
    soilRT = soilRT / sum_water_flux

    ! Sanity check in case of zero flux
    if (sum_water_flux == dble_zero) then
        soilRT = sum(soilRT_local)*0.5d0
        uptake_fraction = dble_zero ; uptake_fraction(1) = dble_one
    endif

    ! Determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = sum(demand) / sum(water_flux)
    ! Finally convert transpiration flux (mmol.m-2.s-1)
    ! into kg.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux * mmol_to_kg_water * seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(out) :: ustar_Uh ! ratio of friction velocity over wind speed at canopy top
    ! local variables
    double precision  sqrt_cd1_lai &
                     ,local_lai &
                     ,phi_h       ! roughness sublayer influence function
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
!                          ustar_Uh_max = 0.3,   & ! Maximum observed ratio of
                                                   ! (friction velocity / canopy top wind speed) (m.s-1)
                          ustar_Uh_max = 1d0, ustar_Uh_min = 0.2d0, &
                               min_lai = 1d0,   & ! Minimum LAI parameter as height does not vary with growth
                                    Cw = 2d0      ! Characterises roughness sublayer depth (m)

    ! assign new value to min_lai to avoid max min calls
    local_lai = max(min_lai,lai)
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (dble_one-((dble_one-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
    ! calculate roughness sublayer influence function;
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    phi_h = 0.19314718056d0
!    phi_h = log(Cw)-dble_one+Cw**(-dble_one) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    roughl = ((dble_one-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

    ! sanity check
!    if (roughl /= roughl) then
!        write(*,*)"TLS:  ERROR roughness length calculations"
!        write(*,*)"Roughness lenght", roughl, "Displacement", displacement
!        write(*,*)"canopy height", canopy_height, "lai", lai
!    endif

  end subroutine z0_displacement
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )

    ! The equation is simply...                        !
    !    a * exp( b * ( t - 25.0 ) / ( t + 273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split & d0 has been used. !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    ! local variables..
    double precision :: answer, denominator, numerator

    numerator   = t - 25d0
    denominator = t + freeze
    answer      = a * exp( b * dble_one * numerator / denominator )
    arrhenious  = answer

  end function arrhenious
  !
  !------------------------------------------------------------------
  !
  double precision function daylength_in_hours(doy,lat)

    ! Function uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate
    ! day length in hours

    implicit none

    ! arguments
    double precision, intent(in) :: doy, lat

    ! local variables
    double precision :: dec, mult, sinld, cosld, aob

!    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
    mult = lat * deg_to_rad
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-dble_one,min(dble_one,sinld / cosld))

    ! define output
    daylength_in_hours = 12d0 * ( dble_one + 2d0 * asin( aob ) * pi_1 )

    ! return to user
    return

  end function daylength_in_hours
  !
  !------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x,y,interval)

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer :: interval
    double precision, dimension(interval) :: x,y

    ! declare local variables
    double precision :: sum_x, sum_y, sumsq_x,sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    sum_product_xy = sum(x*y)
    ! calculate the gradient
    linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val , optimum , kurtosis , current )

    ! estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis

    implicit none

    ! arguments..
    double precision,intent(in) :: max_val, optimum, kurtosis, current

    ! local variables..
    double precision :: dummy

    if ( current >= max_val ) then
       opt_max_scaling = dble_zero
    else
       dummy     = ( max_val - current ) / ( max_val - optimum )
       dummy     = exp( log( dummy ) * kurtosis * ( max_val - optimum ) )
       opt_max_scaling = dummy * exp( kurtosis * ( current - optimum ) )
    end if

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  double precision function root_resistance (root_mass,thickness)

   !
   ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
   !

   implicit none

   ! arguments
   double precision :: root_mass, & ! root biomass in layer (gbiomass)
                       thickness    ! thickness of soil zone roots are in

   ! calculate root hydraulic resistance
   root_resistance = root_resist / (root_mass*thickness)
   ! return
   return

  end function root_resistance
  !
  !------------------------------------------------------------------
  !
  double precision function zbrent( called_from , func , x1 , x2 , tol )

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    !  reference to a particular function and also two values which bound !
    !  the arguments for the function of interest. ZBRENT finds a root of !
    !  the function (i.e. the point where the function equals zero), that !
    !  lies between the two bounds.                                       !
    ! For a full description see Press et al. (1986).                     !

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    double precision,intent(in) :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       double precision function func( xval )
         double precision ,intent(in) :: xval
       end function func
    end interface

    ! local variables..
    integer            :: iter
    integer, parameter :: ITMAX = 30
    double precision   :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    double precision, parameter :: EPS = 3d-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0d0 ) then
      zbrent = a
      return
    elseif ( fb .eq. 0d0 ) then
      zbrent = b
      return
    end if
    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
    if ( sign(1d0,fa) .eq. sign(1d0,fb) ) then
       fa = func( a )
       fb = func( b )
       ! tell me otherwise what is going on
!       print*,"Supplied values must bracket the root of the function.",new_line('x'),  &
!         "     ","You supplied x1:",x1,new_line('x'),                     &
!         "     "," and x2:",x2,new_line('x'),                             &
!         "     "," which give function values of fa :",fa,new_line('x'),  &
!         "     "," and fb:",fb," .",new_line('x'),                        &
!         " zbrent was called by: ",trim(called_from)
!       fa = func( a )
!       fb = func( b )
    end if
    c = b
    fc = fb

    do iter = 1 , ITMAX

       ! If the new value (f(c)) doesn't bracket
       ! the root with f(b) then adjust it..
       if ( sign(1d0,fb) .eq. sign(1d0,fc) ) then
          c  = a
          fc = fa
          d  = b - a
          e  = d
       end if
       if ( abs(fc) .lt. abs(fb) ) then
          a  = b
          b  = c
          c  = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2d0 * EPS * abs(b) + 0.5d0 * tol
       xm   = 0.5d0 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
          zbrent = b
          return
       end if
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
          s = fb / fa
          if ( a .eq. c ) then
             p = 2d0 * xm * s
             q = 1d0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * ( 2d0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1d0 ) )
             q = ( q - 1d0 ) * ( r - 1d0 ) * ( s - 1d0 )
          end if
          if ( p .gt. 0d0 ) q = -q
          p = abs( p )
          if ( (2d0*p) .lt. min( 3d0*xm*q-abs(tol1*q) , abs(e*q) ) ) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a  = b
       fa = fb
       if ( abs(d) .gt. tol1 ) then
          b = b + d
       else
          b = b + sign( tol1 , xm )
       end if
       fb = func(b)
    enddo

!    print*,"zbrent has exceeded maximum iterations",new_line('x'),&
!           "zbrent was called by: ",trim(called_from)

    zbrent = b

  end function zbrent
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEl_MOD
