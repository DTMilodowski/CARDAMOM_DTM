
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL        &
         ,acm                 &
         ,acm_et              &
         ,calculate_update_soil_water &
         ,calculate_Rtot      &
         ,calculate_aerodynamic_conductance &
         ,saxton_parameters   &
         ,initialise_soils    &
         ,seconds_per_day,seconds_per_step &
         ,root_biomass        &
         ,root_reach          &
         ,minlwp              &
         ,min_root            &
         ,max_depth           &
         ,root_k              &
         ,top_soil_depth      &
         ,soil_depth,previous_depth &
         ,nos_root_layers     &
         ,wSWP,SWP,SWP_initial&
         ,deltat_1            &
         ,water_flux          &
         ,layer_thickness     &
         ,waterloss,watergain &
         ,potA,potB           &
         ,cond1,cond2,cond3   &
         ,soil_conductivity   &
         ,soil_waterfrac,soil_waterfrac_initial &
         ,porosity,porosity_initial &
         ,field_capacity,field_capacity_initial &
         ,soilwatermm      &
         ,Rtot_time        &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
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
! ACM related parameters
double precision, parameter :: pi = 3.1415927
double precision, parameter :: pi2 = pi*2.0
double precision, parameter :: deg_to_rad = pi/180.0

! management and gsi related values
integer :: gsi_lag_remembered 
! local variables for GSI phenology model
double precision :: Tfac,Photofac,VPDfac & ! oC, seconds, Pa
                   ,SLA & ! Specific leaf area
                   ,avail_labile,Rg_from_labile    &
                   ,Cwood_labile_release_gradient  &
                   ,Cwood_labile_half_saturation   &
                   ,Croot_labile_release_gradient  &
                   ,Croot_labile_half_saturation   &
                   ,Cwood_hydraulic_gradient       &
                   ,Cwood_hydraulic_half_saturation&
                   ,Cwood_hydraulic_limit          &
                   ,delta_gsi,tmp,gradient         &
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

! Variables related to integration by ordinary differenctial
! equations
integer             :: kmax   !
double precision    :: dxsav  !
integer,parameter :: fname_length = 100   & ! length of filename variables
              ,max_nos_iterations = 2     & ! number of iterations in math-loops
                          ,kmaxx  = 200   & ! descriptions
                          ,maxstp = 10000   !   would
double precision,parameter :: tiny   = 1.e-30  !        nice!

! hydraulic model parameters
integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, parameter :: gravity = 9.8067,       & ! acceleration due to gravity, ms-1
                                minlwp = -2.060814,    & ! min leaf water potential (MPa)
                             vonkarman = 0.41,         & ! von Karman's constant
                           vonkarman_2 = 0.41**2,      & ! von Karman's constant^2
                                 cpair = 1004.6,       & ! Specific heat capacity of air; used in energy balance J.kg-1.K-1
                       seconds_per_day = 86400,        & ! number of seconds per day
                     seconds_per_day_1 = 1d0/86400,    & ! inverse of seconds per day
                                  xacc = 0.0001,       & ! accuracy parameter for zbrent bisection proceedure ! 0.0001
                        mol_to_g_water = 18.0,         & ! molecular mass of water
                      density_of_water = 998.9,        & ! density of water kg.m-3 
                                gplant = 5.0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                                  head = 0.009807,     & ! head of pressure  (MPa/m)
                                head_1 = 101.968,      & ! inverse head of pressure (m/MPa)
                           root_radius = 0.0001,       & ! root radius (m) Bonen et al 2014 = 0.00029
                   root_cross_sec_area = 3.141593e-08, & ! root cross sectional area (m2)
                                                         ! = pi * root_radius * root_radius 
                          root_density = 0.5e6,        & ! root density (g biomass m-3 root) 
                                                         ! 0.5e6 Williams et al 1996                                       
                                                         ! 0.31e6 Bonan et al 2014
                 root_mass_length_coef = root_cross_sec_area * root_density, &
                         canopy_height = 9.0,          & ! canopy height assumed to be 9 m
                        top_soil_depth = 0.3,          & ! depth to which we conider the top soil to extend (m)
                              min_root = 5.0,          & ! minimum root biomass (gBiomass.m-2)
!                             max_depth = 2.0,          & ! maximum possible root depth (m)
!                                root_k = 100,          & ! biomass to reach half max_depth
                           root_resist = 25.0            ! Root resistivity (MPa s g mmol−1 H2O)

! useful technical parameters
double precision, parameter :: dble_zero = 0.0, &
                                dble_one = 1.0

! hydraulic model variables 
integer :: water_retention_pass, soil_layer
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand
double precision, dimension(nos_root_layers) :: uptake_fraction
double precision :: root_reach, root_biomass,soil_depth, &
                    demand, & ! maximum potential canopy hydraulic demand
!              surf_biomass, & ! 
                    soilRT, &
                      wSWP, & ! weighted soil water potential (MPa) used in GSI calculate. 
                              ! Removes / limits the fact that very low root density in young plants
                              ! give values too large for GSI to handle.
                 max_depth, & ! maximum possible root depth (m)
                    root_k, & ! biomass to reach half max_depth
   liquid,drainlayer,unsat, & ! variables used in drainage (m)
                    runoff, & ! runoff (kg.m-2.step)
          seconds_per_step, & !
                        x1, & ! lower boundary condition for zbrent calculation
                        x2, & ! upper boundary condition for zbrent calculation
  new_depth,previous_depth, & ! depth of bottom of soil profile
   aerodynamic_conductance, & ! bulk surface layer conductance
                    roughl, & ! roughness length (m)
              displacement    ! zero plane displacement (m)

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days 
                                 Cwood_labile_release_coef, & ! time series of labile release to wood
                                 Croot_labile_release_coef, & ! time series of labile release to root
                                           layer_thickness, & ! thickness of soil layers
                                         soil_conductivity, & ! soil conductivity
                                                water_flux, & ! potential transpiration flux (mmol.m-2.s-1)
                                                  porosity, & ! soil layer porosity, (fraction)
                                          porosity_initial, &
                                                       SWP, & ! soil water potential (MPa)
                                            field_capacity, & ! soil field capacity (m3.m-3)
                                                 waterloss, & ! water loss from specific soil layers (m)
                                                 watergain, & ! water gained by specfic soil layers (m)
                                            soil_waterfrac, & ! soil water content (m3.m-3)
                                    soil_waterfrac_initial, &
                                               SWP_initial, & 
                                    field_capacity_initial, &
                                               soilwatermm, & ! water content (mm) in rooting zone
                                                 Rtot_time, &
                           cond1, cond2, cond3, potA, potB    ! Saxton equation values


contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai,NEE,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP)

    ! The Data Assimilation Linked Ecosystem Carbon - Growing Season
    ! Index - BUCKET (DALEC_GSI_BUCKET) model. 
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP and 
    ! partitions between various ecosystem carbon pools. These pools are
    ! subject to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2
    ! The ACM_ET simulates the potential evapotranspiration and updates the water balance of a simple 3 layer soil model 
    ! into which roots are distributed.
    ! The purpose of the simple hydraulic (or Bucket) model is to link water deficit in the calculation of GPP.

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   & 
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat(nodays)    & ! time step in decimal days
                         ,pars(nopars)      & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout) :: lai & ! leaf area index
                                               ,GPP & ! Gross primary productivity
                                               ,NEE   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
 
    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes
                                             
    ! declare general local variables
    double precision ::  infi & ! used to calculate infinity for diagnositc
                ,Tfac_range_1 &
            ,Photofac_range_1 &
              ,VPDfac_range_1 &
                   ,root_cost &
                   ,root_life &
                       ,meant & ! mean air temperature (oC)
                     ,deltaWP & ! deltaWP (MPa) minlwp-soilWP
                          ,ET & ! Evapotranspiration (kg.m-2.s-1)
                 ,gpppars(12) & ! ACM inputs (LAI+met)
                 ,evappars(6) & ! ACM_ET parameters
               ,constants(10)   ! parameters for ACM

    integer :: p,f,nxp,n,test,m

    ! local fire related variables
    double precision :: CFF(7) = 0, CFF_res(4) = 0    & ! combusted and non-combustion fluxes
                       ,NCFF(7) = 0, NCFF_res(4) = 0  & ! with residue and non-residue seperates
                       ,combust_eff(5)                & ! combustion efficiency
                       ,rfac                            ! resilience factor

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

    integer :: reforest_day, harvest_management,restocking_lag, gsi_lag
   
    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th precipitation (kg.m-2.s-1)
    ! 8th deforestation fraction
    ! 9th burnt area fraction
    ! 10th 21 day average min temperature
    ! 11th 21 day average photoperiod
    ! 12th 21 day average VPD 
    ! 13th Forest management practice to accompany any clearing
    ! 14th avg daily temperature
    ! 15th avg daily wind speed (ms-1)
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
    ! 19 = Evapotranspiration (kgH2O.m-2.day-1)

    ! PARAMETERS
    ! 17+4(GSI)+2(BUCKET)

    ! p(1) Litter to SOM conversion rate  - m_r
    ! p(2) Fraction of GPP respired - f_a
    ! p(3) Fraction of NPP allocated to foliage - f_f 
    ! p(4) Fraction of NPP allocated to roots - f_r
    ! p(5) max leaf turnover (GSI) ! Leaf lifespan - L_f (CDEA)
    ! p(6) Turnover rate of wood - t_w
    ! p(7) Turnover rate of roots - t_r
    ! p(8) Litter turnover rate - t_l
    ! p(9) SOM turnover rate  - t_S
    ! p(10) Parameter in exponential term of temperature - \theta
    ! p(11) Canopy efficiency parameter - C_eff (part of ACM)
    ! p(12) = max labile turnover(GSI) ! date of Clab release - B_day (CDEA)
    ! p(13) = Fraction allocated to Clab - f_l
    ! p(14) = min temp threshold (GSI) ! lab release duration period - R_l (CDEA)
    ! p(15) = max temp threshold (GSI)! date of leaf fall - F_day
    ! p(16) = min photoperiod threshold (GIS) 
    ! p(17) = LMA
    ! p(24) = max photoperiod threshold (GSI)
    ! p(25) = min VPD threshold (GSI)
    ! p(26) = max VPD threshold (GSI)
    ! p(27) = minimum GPP benefit of increased LAI for labile allocation to be allowed
    ! p(28) = fraction of Cwood which is Cbranch
    ! p(29) = fraction of Cwood which is Ccoarseroot
    ! p(37) = Initial CWD pool
    ! p(38) = CWD turnover fraction
    ! p(39) = Fine root (gbiomass.m-2) needed to reach 50% of max depth
    ! p(40) = Maximum rooting depth (m)

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
!call cpu_time(start)
!call cpu_time(finish)

    ! load some values
    gpppars(4) = 10.0**pars(11) !TLS 1 ! foliar N
    gpppars(7) = lat
    gpppars(9) = abs(minlwp) ! leafWP-soilWP (i.e. -2-0) ! p11 from ACM recal
                             ! NOTE: sign is forced positive for use in varous
                             ! equations
    gpppars(10) = 1.0 ! totaly hydraulic resistance ! p12 from ACM recal (updated)
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)=24.24129      ! pars(11) ! p1  from ACM recal
    constants(2)=7.798524e-03  ! 0.0156935! p2  from ACM recal
    constants(3)=154.7495      ! 4.22273  ! p3  from ACM recal
    constants(4)=465.7482      ! 208.868  ! p4  from ACM recal
    constants(5)=7.817923e-02  ! 0.0453194! p5  from ACM recal
    constants(6)=5.674312e-01  ! 0.37836  ! p6  from ACM recal
    constants(7)=1.076729e+01  ! 7.19298  ! p7  from ACM recal
    constants(8)=5.577107e-03  ! 0.011136 ! p8  from ACM recal
    constants(9)=3.154374e+00  ! 2.1001   ! p9  from ACM recal
    constants(10)=3.959395e-01 ! 0.789798 ! p10 from ACM recal

    ! assign acm_et parameters
    evappars(1) = 0.3251466    ! p13 from ACM recal
    evappars(2) = 9.468052e-04 ! p14 from ACM recal
    evappars(3) = 4.700843e-02 ! p15 from ACM recal
    evappars(4) = 94.83173     ! p16 from ACM recal
    evappars(5) = 2.038884e+00 ! p17 from ACM recal
    evappars(6) = 1.275885e+00 ! p18 from ACM recal
    ! plus ones being calibrated
    root_k = pars(39) ; max_depth = pars(40)

    ! initial values for deforestation variables
    labile_loss = 0.    ; foliar_loss = 0.
    roots_loss = 0.     ; wood_loss = 0.
    labile_residue = 0. ; foliar_residue = 0.
    roots_residue = 0.  ; wood_residue = 0.
    stem_residue = 0.   ; branch_residue = 0.
    reforest_day = 0
    soil_loss_with_roots = 0.
    coarse_root_residue = 0.
    post_harvest_burn = 0.

    ! now load the hardcoded forest management parameters into their locations

    ! Parameter values for deforestation variables
    ! scenario 1
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = 1.0
    roots_frac_res(1)   = 1.0
    rootcr_frac_res(1) = 1.0
    branch_frac_res(1) = 1.0
    stem_frac_res(1)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(1) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(1) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(1) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(1) = 1. 

    !## scen 2
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = 1.0
    roots_frac_res(2)   = 1.0
    rootcr_frac_res(2) = 1.0
    branch_frac_res(2) = 1.0
    stem_frac_res(2)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(2) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(2) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(2) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(2) = 0.

    !## scen 3
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5
    roots_frac_res(3)   = 1.0
    rootcr_frac_res(3) = 1.0
    branch_frac_res(3) = 0.
    stem_frac_res(3)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(3) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(3) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(3) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(3) = 0.

    !## scen 4
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5
    roots_frac_res(4)   = 1.0
    rootcr_frac_res(4) = 0.
    branch_frac_res(4) = 0.
    stem_frac_res(4)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(4) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(4) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(4) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(4) = 0.

    ! for the moment override all paritioning parameters with those coming from
    ! CARDAMOM
    Cbranch_part = pars(28)
    Crootcr_part = pars(29)

    ! declare fire constants (labile, foliar, roots, wood, litter)
    combust_eff(1) = 0.1 ; combust_eff(2) = 0.9
    combust_eff(3) = 0.1 ; combust_eff(4) = 0.5
    combust_eff(5) = 0.3 ; rfac = 0.5

    if (start == 1) then

        ! assigning initial conditions
        POOLS(1,1)=pars(18)
        POOLS(1,2)=pars(19)
        POOLS(1,3)=pars(20)
        POOLS(1,4)=pars(21)
        POOLS(1,5)=pars(22)
        POOLS(1,6)=pars(23)
        POOLS(1,7)=pars(37)
        ! POOL(1,8) assigned later

        if (.not.allocated(disturbance_residue_to_som)) then
            allocate(disturbance_residue_to_litter(nodays), &
                     disturbance_residue_to_cwd(nodays),    &
                     disturbance_residue_to_som(nodays),    &
                     disturbance_loss_from_litter(nodays),  &
                     disturbance_loss_from_cwd(nodays),     &
                     disturbance_loss_from_som(nodays))
        endif
        disturbance_residue_to_litter = 0.0 ; disturbance_loss_from_litter = 0.0
        disturbance_residue_to_som = 0.0 ; disturbance_loss_from_som = 0.0
        disturbance_residue_to_cwd = 0.0 ; disturbance_loss_from_cwd = 0.0

        if (.not.allocated(Cwood_labile_release_coef)) then
            allocate(Cwood_labile_release_coef(nodays),Croot_labile_release_coef(nodays))
            ! Wood/fine root labile turnover parameters
            ! parmeters generated on the assumption of 5 % / 95 % activation at key
            ! temperature values. Roots 1oC/30oC, wood 5oC/30oC.
            Croot_labile_release_gradient = 0.2998069 ; Croot_labile_half_saturation = 15.28207
            Cwood_labile_release_gradient = 0.2995754 !0.2995754 ! 0.25
            Cwood_labile_half_saturation  = 17.49752  !17.49752  ! 20.0
            ! calculate temperature limitation on potential wood/root growth
            Cwood_labile_release_coef = (1.0+exp(-Cwood_labile_release_gradient* &
                                      (((met(3,:)+met(2,:))*0.5)-Cwood_labile_half_saturation)))**(-1.0)
            Croot_labile_release_coef = (1.0+exp(-Croot_labile_release_gradient* &
                                      (((met(3,:)+met(2,:))*0.5)-Croot_labile_half_saturation)))**(-1.0)
        endif
        ! hydraulic limitation parameters for wood cell expansion, i.e. growth
        Cwood_hydraulic_gradient = 10.0 ; Cwood_hydraulic_half_saturation = -1.5

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
              if (sum(deltat(1:n)) < 21) then
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
                 do while (test < 21)
                    m=m+1 ; test = sum(deltat((n-m):n))
                    if (m > (n-1)) test = 21
                 end do
                 tmp_m(n) = m
              endif ! for calculating gradient
            end do ! calc daily values once
            ! allocate GSI history dimension
            gsi_lag_remembered=max(2,maxval(nint(tmp_m)))
        end if ! .not.allocated(tmp_x)
        ! assign our starting value
        gsi_history = pars(36)-1d0
         just_grown = pars(35)

        ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
        if (.not.allocated(SWP)) then
           allocate(deltat_1(nodays),soilwatermm(nodays),Rtot_time(nodays))
           deltat_1 = deltat**(-1.0) 
           allocate(water_flux(nos_root_layers),SWP(nos_soil_layers)                           &
                   ,layer_thickness(nos_soil_layers),soil_conductivity(nos_soil_layers+1)      &
                   ,waterloss(nos_soil_layers+1),watergain(nos_soil_layers+1)                  &    
                   ,field_capacity(nos_soil_layers+1),soil_waterfrac(nos_soil_layers+1)        &
                   ,porosity(nos_soil_layers+1),porosity_initial(nos_soil_layers+1)            &
                   ,SWP_initial(nos_soil_layers+1),field_capacity_initial(nos_soil_layers+1)   &
                   ,soil_waterfrac_initial(nos_soil_layers+1)                                  &
                   ,cond1(nos_soil_layers+1),cond2(nos_soil_layers+1),cond3(nos_soil_layers+1) &
                   ,potA(nos_soil_layers+1),potB(nos_soil_layers+1))
           ! zero variables not done elsewhere
           water_flux = 0.0 
           ! initialise some time invarient parameters
           call saxton_parameters(soil_frac_clay,soil_frac_sand)
           call initialise_soils(soil_frac_clay,soil_frac_sand)
           soil_waterfrac_initial = soil_waterfrac
           SWP_initial = SWP
           field_capacity_initial = field_capacity
           porosity_initial = porosity
        else
           water_flux = 0.0
           soil_waterfrac = soil_waterfrac_initial
           SWP = SWP_initial
           field_capacity = field_capacity_initial
           porosity = porosity_initial
        endif

        ! zero evapotranspiration for beginning
        ET = 0.0 ; seconds_per_step = seconds_per_day * deltat(1)
        ! initialise root reach based on initial conditions
        root_biomass = max(min_root,POOLS(1,3)*2)
        root_reach = max_depth * root_biomass / (root_k + root_biomass)
        ! determine initial soil layer thickness
        layer_thickness(1) = top_soil_depth ; layer_thickness(2)=max(0.1,root_reach-layer_thickness(1))
        layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
        previous_depth = max(top_soil_depth,root_reach)
        soil_depth = 0.0 ; previous_depth = 0.0
        ! needed to initialise soils 
        call calculate_Rtot(gpppars(9),gpppars(10),met(15,1) &
                           ,deltat(1),((met(3,1)+met(2,1))*0.5))
        ! used to initialise soils
        ET = calculate_update_soil_water(ET*deltat(1),dble_zero,((met(3,1)+met(2,1))*0.5))
        ! store soil water content of the rooting zone (mm)
        POOLS(1,8) = 1e3*sum(soil_waterfrac(1:nos_root_layers)*layer_thickness(1:nos_root_layers))
    
    else
        ! load ET from memory
        ET = FLUXES(start-1,19)
    endif !  start == 1

    ! assign climate sensitivities
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory
    fol_turn_crit=pars(34)-dble_one
    lab_turn_crit=pars(3)-dble_one
    Tfac_range_1 = (pars(15)-pars(14))**(-1.0)
    Photofac_range_1 = (pars(24)-pars(16))**(-1.0)
    VPDfac_range_1 = abs(pars(26)-pars(25))**(-1.0)
    SLA = pars(17)**(-1.0)
    root_cost = dble_zero
    root_life = pars(7)**(-1.0)

    infi=0d0

    ! 
    ! Begin looping through each time step
    ! 

    do n = start, finish

      ! timing variable
      seconds_per_step = seconds_per_day * deltat(n)
      meant = (met(3,n)+met(2,n))*0.5
      ! calculate LAI value
      lai(n)=POOLS(n,2)*SLA
      ! load next met / lai values for ACM and acm_et
      gpppars(1)=lai(n)
      gpppars(2)=met(3,n) ! max temp
      gpppars(3)=met(2,n) ! min temp
      gpppars(5)=met(5,n)!+200 ! co2
      gpppars(6)=ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8)=met(4,n) ! radiation
      tmp = acm(gpppars,constants)

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,POOLS(n,3)*2)
      deltaWP = min(dble_zero,minlwp-wSWP)
!      gpppars(9) = abs(deltaWP) ! update deltaWP

      call calculate_Rtot(gpppars(9),gpppars(10),lai(n) &
                         ,deltat(n),meant)

      ! pass Rtot to output variable and update deltaWP between minlwp and
      ! current weighted soil WP 
      Rtot_time(n) = wSWP !gpppars(10) ! soilRT ! gpppars(10)
      ! calculate aerodynamic resistance (1/conductance) using consistent
      ! approach with SPA
      call calculate_aerodynamic_conductance(lai(n),met(15,n))      

      ! Potential latent energy (kg.m-2.day-1)
      if (deltaWP < 0.0) then
          FLUXES(n,19) = acm_et(meant,met(4,n),met(16,n),lai(n) &
                               ,evappars,gpppars(9),gpppars(10))
      else
          FLUXES(n,19) = 0.0
      endif

      ! do mass balance (i.e. is there enough water to support ET)
      FLUXES(n,19) = calculate_update_soil_water(FLUXES(n,19)*deltat(n),met(7,n)*seconds_per_step,meant)
      ! now reverse the time correction (step -> sec)
      FLUXES(n,19) = FLUXES(n,19) * deltat_1(n)
      ! pass to local variable for soil mass balance
      ET = FLUXES(n,19)
      ! store soil water content of the rooting zone (mm)
      POOLS(n,8) = 1e3*sum(soil_waterfrac(1:nos_root_layers)*layer_thickness(1:nos_root_layers))

      ! GPP (gC.m-2.day-1)
      if (lai(n) > 1e-10 .and. deltaWP < 0.0 .and. abs(gpppars(10)) /= abs(log(infi))) then
         FLUXES(n,1) = acm(gpppars,constants)
      else
         FLUXES(n,1) = 0.0
      endif

      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(10)*meant)
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = 0.0 !(FLUXES(n,1)-FLUXES(n,3))*pars(3)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))

      ! GSI added to fortran version by TLS 24/11/2014
      ! /* 25/09/14 - JFE
      ! Here we calculate the Growing Season Index based on 
      ! Jolly et al., Global Change Biology, Volume 11, page 619-632, 2005 (doi: 10.1111/j.1365-2486.2005.00930.x) 
      ! Stoeckli et al., J. Geophys. Res., 116, G03020, doi:10.1029/2010JG001545.
        
      ! It is the product of 3 limiting factors for temperature, photoperiod and
      ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated 
      ! min and max value. Photoperiod, VPD and avgTmin are direct input

      ! temperature limitation, then restrict to 0-1; correction for k-> oC
      Tfac = (met(10,n)-(pars(14)-273.15)) * Tfac_range_1
      Tfac = min(dble_one,max(dble_zero,Tfac))
      ! photoperiod limitation
      Photofac = (met(11,n)-pars(16)) * Photofac_range_1
      Photofac = min(dble_one,max(dble_zero,Photofac))
      ! VPD limitation
!      ! soil hydraulic limitation (but not root) 
!      VPDfac = 1.0 - ( (soilRT-pars(25)) * VPDfac_range_1 )
!      VPDfac = min(dble_one,max(dble_zero,VPDfac))
      ! deltaWP = minlwp-wSWP (MPa)
      VPDfac = 1.0 - ((deltaWP-pars(25)) * VPDfac_range_1)
      VPDfac = min(dble_one,max(dble_zero,VPDfac))

      ! calculate and store the GSI index (DEFAULT)
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
      FLUXES(n,9) = 0.0  ! leaf turnover
      FLUXES(n,16) = 0.0 ! leaf growth

      ! now update foliage and labile conditions based on gradient calculations
      if (gradient < fol_turn_crit .or. FLUXES(n,18) == 0) then
         ! we are in a decending condition so foliar turnover
         FLUXES(n,9) = pars(5)*(1.0-FLUXES(n,18))
         just_grown = 0.5
      else if (gradient > lab_turn_crit .and. deltaWP < 0.0 .and. abs(gpppars(10)) /= abs(log(infi))) then
         ! we are in a assending condition so labile turnover
         FLUXES(n,16) = pars(12)*FLUXES(n,18)
         just_grown = 1.5
         ! check carbon return
         tmp = POOLS(n,1)*(1.0-(1.0-FLUXES(n,16))**deltat(n))*deltat_1(n)
         tmp = (POOLS(n,2)+tmp)*SLA
         gpppars(1)=tmp
         tmp = acm(gpppars,constants)
         ! determine if increase in LAI leads to an improvement in GPP greater
         ! than critical value, if not then no labile turnover allowed
         if ( (tmp - FLUXES(n,1)) < (pars(27)*FLUXES(n,1)) ) then
             FLUXES(n,16) = 0.0
         endif
      else if (deltaWP < 0.0 .and. abs(gpppars(10)) /= abs(log(infi))) then
         ! probaly we want nothing to happen, however if we are at the seasonal
         ! maximum we will consider further growth still
         if (just_grown >= 1.0) then
            ! we are between so definitely not losing foliage and we have
            ! previously been growing so maybe we still have a marginal return on
            ! doing so again
            FLUXES(n,16) = pars(12)*FLUXES(n,18)
            ! but possibly gaining some?
            ! determine if this is a good idea based on GPP increment
            tmp = POOLS(n,1)*(1.0-(1.0-FLUXES(n,16))**deltat(n))*deltat_1(n)
            tmp = (POOLS(n,2)+tmp)*SLA
            gpppars(1)=tmp
            tmp = acm(gpppars,constants)
            ! determine if increase in LAI leads to an improvement in GPP greater
            ! than
            ! critical value, if not then no labile turnover allowed
            if ( (tmp - FLUXES(n,1)) < (pars(27)*FLUXES(n,1)) ) then
                FLUXES(n,16) = 0.0
            endif
         end if ! Just grown?
      endif ! gradient choice

      ! total labile release to foliage
      FLUXES(n,8) = POOLS(n,1)*(1.0-(1.0-FLUXES(n,16))**deltat(n))*deltat_1(n)
      avail_labile = max(dble_zero,POOLS(n,1) - (FLUXES(n,8)*deltat(n)))

      ! these allocated if post-processing
      if (allocated(itemp)) then
         itemp(n) = Tfac
         ivpd(n) = VPDfac
         iphoto(n) = Photofac
      endif

      ! Premise of wood and root phenological controls

      ! Assumption 1:
      ! Based on plant physiology all cell expansion can only occur if there is 
      ! sufficient water pressure available to drive the desired expansion. 
      ! Moreover, as there is ample evidence that shows wood and root growth 
      ! do not follow the same phenology as leaves or GPP availability.
      ! Therefore, their phenological constrols should be separate from both that of 
      ! the GSI model driving canopy phenology or GPP. Wood growth is limited by a 
      ! logisitic temperature response assuming <5 % growth potential at 5oC and
      ! >95 % growth potential at 30 oC. Wood growth is also limited by a
      ! logistic response to water availability. When deltaWP (i.e. minleaf-wSWP) 
      ! is less than -1 MPa wood growth is restricted to <5 % of potential. 

      ! Assumption 2:
      ! As with wood, root growth biologically speaking is independent of
      ! observed foliar phenological dynamics and GPP availabilty and thus has 
      ! a separate phenology model. Similar to wood, a logistic temperature
      ! response is applied such that root growth is <5 % of potential at 0oC and 
      ! >95 % of potential at 30oC. The different temperature minimua between
      ! wood and root growth is due to observed root growth when ever the soil
      ! is not frozen. We also assume that root growth is less sensitive to
      ! available hydraulic pressure, see assumption 4.

      ! Assumption 3:
      ! Actual biological theory suggests that roots support demands for resources
      ! made by the rest of the plant in this current model this is water only.
      ! Therefore there is an implicit assumption that roots should grow so long as 
      ! growth is environmentally possible as growth leads to an improvement in 
      ! C balance over their life time greater than their construction cost. 
      ! NOTE: here we accept any positive impact on instanteous GPP.

      ! Assumption 4:
      ! Determining when root growth should stop is poorly constrained. 
      ! Similar to wood growth, here we assume root expansion is also dependent on water availability,
      ! but is less sensitive than wood. Root growth is assumed to stop when deltaWP approaches 0.
      ! determine root growth based on marginal returns, temperature and hydraulic limits.

      ! Assumption 5: 
      ! The fact that NPP allocation fractions vary both with plant age and
      ! environmental conditions is well established though poorly modelled from
      ! a mechanistic perspective. Here, we argue that if the GSI model
      ! indicates that the canopy is undergoing senescence due to hydraulic
      ! limitations, that the plant will respond by directing an increased
      ! fraction of NPP to roots assuming, assuming both Assumption 3 and 4 are
      ! met. Given that wood growth is also hampered by hydraulic limits, we
      ! assume that NPP normally allocated to wood growth is re-directed towards
      ! root growth in this scenario.

      ! reset allocation to roots and wood
      FLUXES(n,6) = 0.0 ; FLUXES(n,7) = 0.0

      ! Based on simplified assumptions regarding wood and root growth and
      ! allocation. Assumes that both wood and root growth is primarily based on hydraulic
      ! and temperature limitations. Root growth has a second constraint that
      ! growth requires that growth leads to a marginal return on GPP which over
      ! the life of the root will be greater than the rNPP and Rg_root costs. 
      ! Breaking the link to canopy hydraulic limtations from the GSI model.
      if ( deltaWP < dble_zero ) then
          ! Assume potential root growth is dependent on hydraulic and temperature conditions. 
          ! Actual allocation is only allowed if the marginal return on GPP,
          ! averaged across the life span of the root is greater than the rNPP and Rg_root.

          ! Temperature limited turnover rate of labile -> roots
          FLUXES(n,6) = pars(4)*Croot_labile_release_coef(n)
          ! Estimate potential root allocation over time for potential root allocation
          tmp = avail_labile*(1.0-(1.0-FLUXES(n,6))**deltat(n))*deltat_1(n)
          root_cost = tmp*deltat(n)
          ! adjust to extressed biomass (i.e. less Rg_root)
          tmp = (POOLS(n,3)+(root_cost*0.78125))*2
          ! estimate new Rtot and then GPP 
          tmp = calc_pot_root_alloc_Rtot(gpppars(9),meant,lai(n),tmp)
          gpppars(1) = lai(n) ; gpppars(10) = tmp ; tmp = acm(gpppars,constants)
          ! calculate marginal return on new root growth, scaled over life span
          ! of new root.
          tmp = (tmp-FLUXES(n,1))*root_life
          ! if return is negative or canopy under no stress don't grow new roots
!          if (tmp < root_cost) FLUXES(n,6) = dble_zero
          if (tmp < dble_zero) FLUXES(n,6) = dble_zero

          ! calculate hydraulic limits on wood growth. 
          ! NOTE: PARAMETERS NEED TO BE CALIBRATRED
          Cwood_hydraulic_limit = (1.0+exp(Cwood_hydraulic_gradient*(deltaWP-Cwood_hydraulic_half_saturation)))**(-1.0)
          ! determine wood growth based on temperature and hydraulic limits
          FLUXES(n,7) = pars(13)*Cwood_labile_release_coef(n)*Cwood_hydraulic_limit

      endif ! grow root and wood?

      ! track labile reserves to ensure that fractional losses are applied
      ! sequencially in assumed order of importance (leaf->root->wood)

      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = avail_labile*(1.0-(1.0-FLUXES(n,6))**deltat(n))*deltat_1(n)
      avail_labile = max(dble_zero,avail_labile - (FLUXES(n,6)*deltat(n)))
      ! wood production (gC.m-2.day-1)
      FLUXES(n,7) = avail_labile*(1.0-(1.0-FLUXES(n,7))**deltat(n))*deltat_1(n)

      !
      ! litter creation with time dependancies
      !

      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1.0-(1.0-FLUXES(n,9))**deltat(n))*deltat_1(n)
      ! total wood litter production
      FLUXES(n,11) = POOLS(n,4)*(1.0-(1.0-pars(6))**deltat(n))*deltat_1(n)
      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*(1.0-(1.0-pars(7))**deltat(n))*deltat_1(n)

      ! 
      ! those with temperature AND time dependancies
      ! 

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*(1.0-(1.-FLUXES(n,2)*pars(8))**deltat(n))*deltat_1(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*(1.0-(1.-FLUXES(n,2)*pars(9))**deltat(n))*deltat_1(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*(1.0-(1.-FLUXES(n,2)*pars(1))**deltat(n))*deltat_1(n)
      ! CWD to litter
      FLUXES(n,20) = POOLS(n,7)*(1.0-(1.-FLUXES(n,2)*pars(38))**deltat(n))*deltat_1(n)

      ! calculate growth respiration and adjust allocation to pools assuming
      ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)

      ! foliage 
      Rg_from_labile = FLUXES(n,8)*0.21875 ; FLUXES(n,8) = FLUXES(n,8) * 0.78125
      ! roots
      Rg_from_labile = Rg_from_labile + (FLUXES(n,6)*0.21875) ; FLUXES(n,6) = FLUXES(n,6) * 0.78125
      ! wood
      Rg_from_labile = Rg_from_labile + (FLUXES(n,7)*0.21875) ; FLUXES(n,7) = FLUXES(n,7) * 0.78125
      ! now update the Ra flux
      FLUXES(n,3) = FLUXES(n,3) + Rg_from_labile

      ! calculate the NEE 
      NEE(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! load GPP
      GPP(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      ! 

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8)-FLUXES(n,6)-FLUXES(n,7)-Rg_from_labile)*deltat(n)
      ! foliar pool
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)-FLUXES(n,10)+FLUXES(n,8))*deltat(n)
      ! wood pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
      ! root pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)+FLUXES(n,20)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14))*deltat(n)
      ! cwd pool
      POOLS(n+1,7) = POOLS(n,7) + (FLUXES(n,11)-FLUXES(n,20))*deltat(n)

      ! 
      ! deal first with deforestation
      ! 

      if (n == reforest_day) then
          POOLS(n+1,1) = pars(30) 
          POOLS(n+1,2) = pars(31) 
          POOLS(n+1,3) = pars(32) 
          POOLS(n+1,4) = pars(33) 
      endif 

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
          soil_loss_with_roots = Crootcr*met(8,n)*(1.-rootcr_frac_res(harvest_management)) &
                              * soil_loss_frac(harvest_management)

          ! update living pools directly
          POOLS(n+1,1) = max(dble_zero,POOLS(n+1,1)-labile_loss)
          POOLS(n+1,2) = max(dble_zero,POOLS(n+1,2)-foliar_loss)
          POOLS(n+1,3) = max(dble_zero,POOLS(n+1,3)-roots_loss)
          POOLS(n+1,4) = max(dble_zero,POOLS(n+1,4)-wood_loss)
          ! then work out the adjustment due to burning if there is any
          if (post_harvest_burn(harvest_management) > 0.) then
              !/*first fluxes*/
              !/*LABILE*/
              CFF(1) = POOLS(n+1,1)*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF(1) = POOLS(n+1,1)*post_harvest_burn(harvest_management)*(1-combust_eff(1))*(1-rfac)
              CFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*(1-combust_eff(1))*(1-rfac)
              !/*foliar*/
              CFF(2) = POOLS(n+1,2)*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF(2) = POOLS(n+1,2)*post_harvest_burn(harvest_management)*(1-combust_eff(2))*(1-rfac)
              CFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*(1-combust_eff(2))*(1-rfac)
              !/*root*/
              CFF(3) = 0. !POOLS(n+1,3)*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF(3) = 0. !POOLS(n+1,3)*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
              CFF_res(3) = 0. !roots_residue*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF_res(3) = 0. !roots_residue*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
              !/*wood*/
              CFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*(1-combust_eff(4))*(1-rfac)
              CFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*(1-combust_eff(4))*(1-rfac)
              !/*litter*/
              CFF(5) = POOLS(n+1,5)*post_harvest_burn(harvest_management)*combust_eff(5)
              NCFF(5) = POOLS(n+1,5)*post_harvest_burn(harvest_management)*(1-combust_eff(5))*(1-rfac)
              !/*CWD*/ Using Combustion factors for wood
              CFF(7) = POOLS(n+1,7)*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF(7) = POOLS(n+1,7)*post_harvest_burn(harvest_management)*(1-combust_eff(4))*(1-rfac)
              !/*fires as daily averages to comply with units*/
              FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)+CFF(7) & 
                           +CFF_res(1)+CFF_res(2)+CFF_res(3)+CFF_res(4))/deltat(n)
              ! update the residue terms
              labile_residue = labile_residue - CFF_res(1) - NCFF_res(1)
              foliar_residue = foliar_residue - CFF_res(2) - NCFF_res(2)
              roots_residue  = roots_residue  - CFF_res(3) - NCFF_res(3)
              wood_residue   = wood_residue   - CFF_res(4) - NCFF_res(4)
              ! now update NEE
              NEE(n)=NEE(n)+FLUXES(n,17)
          else
              FLUXES(n,17) = 0.
              CFF = 0. ; NCFF = 0.
              CFF_res = 0. ; NCFF_res = 0.
          end if
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
          ! this is intended for use with the R interface for subsequent post
          ! processing
          FLUXES(n,21) =  (wood_loss-(wood_residue+CFF_res(4)+NCFF_res(4))) &
                           + (labile_loss-(labile_residue+CFF_res(1)+NCFF_res(1))) &
                           + (foliar_loss-(foliar_residue+CFF_res(2)+NCFF_res(2))) &
                           + (roots_loss-(roots_residue+CFF_res(3)+NCFF_res(3)))
          ! total carbon loss from the system
          C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)) &
                  - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF))

          ! if total clearance occured then we need to ensure some minimum
          ! values and reforestation is assumed one year forward
          if (met(8,n) > 0.99) then
              m=0 ; test=sum(deltat(n:(n+m)))
              ! FC Forest Statistics 2015 lag between harvest and restocking ~ 2 year
              restocking_lag = 365*2
              do while (test < restocking_lag)
                 m=m+1 ; test = sum(deltat(n:(n+m)))
                 !  get out clause for hitting the end of the simulation
                 if (m+n >= nodays) test = restocking_lag
              enddo
              reforest_day = min((n+m), nodays)
          endif ! if total clearance

      endif ! end deforestation info

      ! 
      ! then deal with fire
      ! 

      if (met(9,n) > 0.) then

         !/*first fluxes*/
         !/*LABILE*/
         CFF(1) = POOLS(n+1,1)*met(9,n)*combust_eff(1)
         NCFF(1) = POOLS(n+1,1)*met(9,n)*(1-combust_eff(1))*(1-rfac)
         !/*foliar*/
         CFF(2) = POOLS(n+1,2)*met(9,n)*combust_eff(2)
         NCFF(2) = POOLS(n+1,2)*met(9,n)*(1-combust_eff(2))*(1-rfac)
         !/*root*/
         CFF(3) = 0. ! POOLS(n+1,3)*met(9,n)*combust_eff(3)
         NCFF(3) = 0. ! POOLS(n+1,3)*met(9,n)*(1-combust_eff(3))*(1-rfac)
         !/*wood*/
         CFF(4) = POOLS(n+1,4)*met(9,n)*combust_eff(4)
         NCFF(4) = POOLS(n+1,4)*met(9,n)*(1-combust_eff(4))*(1-rfac)
         !/*litter*/
         CFF(5) = POOLS(n+1,5)*met(9,n)*combust_eff(5)
         NCFF(5) = POOLS(n+1,5)*met(9,n)*(1-combust_eff(5))*(1-rfac)
         !/*fires as daily averages to comply with units*/
         FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5))/deltat(n)

         !/*all fluxes are at a daily timestep*/
         NEE(n)=NEE(n)+FLUXES(n,17)

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

      endif ! end burnst area issues

      do nxp = 1, nopools
         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp)) then
            print*,"step",n, nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
!            print*,"Rtot",Rtot_time(n)
            print*,"wSWP",wSWP
            print*,"waterfrac",soil_waterfrac
            stop
         endif
      enddo

    end do ! nodays loop

    ! do mass balance
    ET = calculate_update_soil_water(FLUXES(n,19)*deltat(n),met(7,n)*seconds_per_step,meant)

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm(drivers,constants)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: drivers(12) & ! acm input requirements
                         ,constants(10) ! ACM parameters

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, ci, e0, dayl, cps, dec, nit &
             ,sinld, cosld,aob, mult &
             ,mint,maxt,radiation,co2,lai,doy,lat &
             ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
             ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
             ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! load driver values to correct local vars
    lai = drivers(1)
    maxt = drivers(2)
    mint = drivers(3)
    nit = drivers(4)   
    co2 = drivers(5)
    doy = drivers(6)
    radiation = drivers(8)
    lat = drivers(7)

    ! load parameters into correct local vars
    deltaWP = drivers(9)
    Rtot = drivers(10)
    NUE = constants(1)
    dayl_coef = constants(2)
    co2_comp_point = constants(3) 
    co2_half_sat = constants(4)
    dayl_const = constants(5)
    hydraulic_temp_coef = constants(6)
    lai_coef = constants(7)
    temp_exponent = constants(8)
    lai_const = constants(9)
    hydraulic_exponent = constants(10)

    ! daily canopy conductance, of CO2 or H2O? 
    gc=deltaWP**(hydraulic_exponent)/(hydraulic_temp_coef*Rtot)
    ! daily canopy conductance but now consistent with the demand term of acm_et 
!    gc=(deltaWP+head*canopy_height)**(hydraulic_exponent)/(hydraulic_temp_coef*Rtot)
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn=lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp=pn/gc ; qq=co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    mult = co2+qq-pp
    ci=0.5*(mult+sqrt(((mult)*(mult))-4.0*(co2*qq-pp*co2_comp_point)))
    ! limit maximum quantium efficiency by leaf area, hyperbola
    mult = lai*lai
    e0=lai_coef*(mult)/((mult)+lai_const)
    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    mult = lat*deg_to_rad
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    dayl = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )
    ! calculate CO2 limited rate of photosynthesis
    pd=gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps=e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm=cps*(dayl_coef*dayl+dayl_const)

    ! don't forget to return
    return

  end function acm
  !
  !------------------------------------------------------------------
  !
  double precision function acm_et(meant,swrad_MJ,vpd_pa,lai & 
                                  ,pars,deltaWP,Rtot)

    ! simple response function based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kg.m-2.day-1)
    
    implicit none

    ! arguments
    double precision, intent(in) :: meant    & ! daily mean temperature (oC)
                                   ,swrad_MJ & ! daily short wave rad (MJ.m-2.day-1)
                                   ,vpd_pa   & ! daily mean vapoure pressure deficit (Pa)
                                   ,lai      & ! daily leaf area index (m2/m2)
                                   ,deltaWP  & ! soil-minleafwp (MPa)
                                   ,Rtot       ! total soil root hydraulic resistance (MPa)

    double precision, intent(in), dimension(6) :: pars

    ! local variables
    double precision :: mult     & 
                       ,swrad    & ! daily mean radiation 
                       ,fun_rad  &
                       ,fun_Rtot &
                       ,s,slope  &
                       ,psych    &
                       ,lambda   &
                       ,rho      &                     
                       ,gc 
 
    ! first determine mean shortwave radiation as (J.m-2.day-1)
    swrad = swrad_MJ*1e6
    ! calculate canopy conductance of evaporation. Assumes logitistic functions
    ! to radiation and hydraulic resistance
    fun_rad  = (1.0 + exp(-pars(3)*((swrad*seconds_per_day_1)-pars(4))))**(-1.0)
    fun_Rtot = (1.0 + exp(pars(5)*(sqrt(Rtot)-pars(6))))**(-1.0)
    gc = lai * deltaWP * fun_rad * fun_Rtot

    ! calculate coefficient for Penman Montieth
    ! density of air (kg.m-3)
    rho = 353.0/(meant+273.15)
    if (meant < 1.0) then
        lambda = 2.835e6
    else 
       ! latent heat of vapourisation (J.kg-1)
       lambda = 2501000.0-2364.0*meant
    endif
    ! psychrometric constant (kPa K-1)
    psych = (0.0646*exp(0.00097*meant))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult=meant+237.3
    s = 0.61078*17.269*237.3*exp(17.269*meant/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)
    ! calculate numerator of Penman Montheith (kg.m-2.day-1); note vpd Pa->kPa
    acm_et = (pars(1)*(slope*(swrad-pars(2))) + (rho*cpair*aerodynamic_conductance*vpd_pa*1e-3))
    !  then demoninator
    acm_et = acm_et / (lambda*(slope+(psych*(aerodynamic_conductance/gc))))

    return

  end function acm_et
  !
  !------------------------------------------------------------------
  !
  double precision function calc_pot_root_alloc_Rtot(deltaWP,meant,lai,root_biomass)

    ! 
    ! Description
    !

    implicit none

    ! declare arguments
    double precision,intent(in) :: root_biomass & ! potential root biomass (g.m-2)
                                  ,deltaWP      &
                                  ,meant        &
                                  ,lai

    ! declare local variables
    double precision, dimension(nos_root_layers) :: water_flux_local &
                                                   ,root_mass    &
                                                   ,root_length  &
                                                   ,ratio
    double precision, dimension(nos_soil_layers) :: layer_thickness_save(nos_soil_layers)
    double precision, dimension(nos_soil_layers+1) :: soil_waterfrac_save 
    double precision :: slpa,transpiration_resistance,root_reach_local &
                       ,demand,soilR1,soilR2,depth_change

    ! estimate rooting depth with potential root growth
    root_reach_local = max_depth * root_biomass / (root_k + root_biomass)   
   
    ! save soil water information
    soil_waterfrac_save = soil_waterfrac

    ! if roots extent down into the bucket 
    if (root_reach_local > layer_thickness(1)) then
       ! how much has root depth extended since last step?
       depth_change=root_reach_local-previous_depth
       ! if there has been an increase
       if (depth_change > 0.0) then
           ! calculate weighting between current lowest root layer and new soil 
           depth_change = depth_change / (depth_change+layer_thickness(nos_root_layers))
           ! add to bottom root layer
           soil_waterfrac(nos_root_layers) = (soil_waterfrac(nos_root_layers)*(1.0-depth_change)) &
                                           + (soil_waterfrac(nos_soil_layers)*depth_change)
       else
           ! calculate weighting between bottom soil layer and new bit coming
           ! from lowest root 
           depth_change = abs(depth_change) / (abs(depth_change)+layer_thickness(nos_root_layers))
           ! and add back to the bottom soil layer
           soil_waterfrac(nos_soil_layers) = (soil_waterfrac(nos_soil_layers)*(1.0-depth_change)) &
                                           + (soil_waterfrac(nos_root_layers)*depth_change)
       end if ! depth change 
 
    end if ! root reach beyond top layer
    ! determine soil layer thickness
    layer_thickness_save = layer_thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2)=max(0.1,root_reach-layer_thickness(1))
    layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

    ! estimate water flux based on soil and root hydraulic resistances with potential growth. 
    ! See subroutine calculate_Rtot for further details

    ! calculate the plant hydraulic resistance component
    transpiration_resistance = (gplant * lai)**(-1.0)
    ! top 25 % of root profile
    slpa = (root_reach_local * 0.25) - layer_thickness(1)
    if (slpa <= dble_zero) then
        ! > 50 % of root is in top layer
        root_mass(1) = root_biomass * 0.5
        root_mass(1) = root_mass(1) + ((root_biomass-root_mass(1)) * (abs(slpa)/root_reach_local))
        root_mass(2) = max(dble_zero,root_biomass - root_mass(1))
    else
        ! < 50 % of root is in bottom layer
        root_mass(1) = root_biomass * 0.5 * (layer_thickness(1)/(abs(slpa)+layer_thickness(1)))
        root_mass(2) = max(dble_zero,root_biomass - root_mass(1))
    endif
    root_length = root_mass / root_mass_length_coef !(root_density * root_cross_sec_area)
    !! Top root layer.
    ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
    soilR1=soil_resistance(root_length(1),min(root_reach_local,layer_thickness(1)) &
                          ,calculate_soil_conductivity(soil_waterfrac(1))*head_1)
    soilR2=root_resistance(root_mass(1),min(root_reach_local,layer_thickness(1)))
    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    demand = deltaWP+head*canopy_height
!    water_flux_local(1) = max(dble_zero,demand/(transpiration_resistance + soilR1 + soilR2))
    water_flux_local(1) = demand/(transpiration_resistance + soilR1 + soilR2)
    ! Bottom root layer
    if (root_mass(2) > 0.0 ) then
       ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
       soilR1=soil_resistance(root_length(2),layer_thickness(2) &
                             ,calculate_soil_conductivity(soil_waterfrac(2))*head_1)
       soilR2=root_resistance(root_mass(2),layer_thickness(2))
       ! calculate and accumulate steady state water flux in mmol.m-2.s-1
!       water_flux_local(2) = max(dble_zero,demand/(transpiration_resistance + soilR1 + soilR2))
       water_flux_local(2) = demand/(transpiration_resistance + soilR1 + soilR2)
    endif ! roots present in second layer?

    ! if freezing then assume soil surface is frozen
    if (meant < 1.0) water_flux_local(1) = dble_zero

    ! WARNING: should probably have updated the wSWP here as well...do this
    ! later I thinks...

    ! determine effective resistance (MUST THINK ON WHETHER THIS IS RIGHT)
    ratio = layer_thickness/sum(layer_thickness(1:nos_root_layers))
    calc_pot_root_alloc_Rtot = demand / sum(water_flux_local*ratio)

    ! return layer_thickness and soil_waterfrac back to
    ! orginal values
    layer_thickness = layer_thickness_save 
    soil_waterfrac = soil_waterfrac_save

    return

  end function calc_pot_root_alloc_Rtot
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_aerodynamic_conductance(lai,wind_spd)

    ! 
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are with SPA at the time
    ! of the calibration
    ! 

    implicit none

    ! arguments
    double precision, intent(in) :: lai, wind_spd

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(lai)

    ! calculate bulk conductance (Jones p68)
    aerodynamic_conductance = (wind_spd * vonkarman_2) / (log((canopy_height-displacement)/roughl))**2

  end subroutine calculate_aerodynamic_conductance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !

    implicit none

    ! local variables..
    integer        :: i

    x1 = 0.1 ; x2 = 0.7 ! low/high guess
    do i = 1 , nos_soil_layers+1
      water_retention_pass = i
      ! field capacity is water content at which SWP = -10 kPa
      field_capacity( i ) = zbrent( 'water_retention:water_retention_saxton_eqns' , water_retention_saxton_eqns , x1 , x2 , xacc )
    enddo

  end subroutine calculate_field_capacity
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot(deltaWP,Rtot,lai,deltat,meant)
  
    ! purpose of this function is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The minimum is assumed to be the same and the
    ! soil layer with the greated root content. Here we use the same assumption
    ! as used in SPA to calcule the root hydraulic resistance in the top soil
    ! layer only. 

    ! This could be extended later for include a 2 soil layer bucket model

    ! declare inputs
    double precision,intent(in) :: deltat,deltaWP,lai,meant
    double precision,intent(inout) :: Rtot

    ! local variables
    integer :: i
    double precision :: slpa, cumdepth, prev, curr, &
                        soilR1,soilR2,transpiration_resistance
    double precision, dimension(nos_root_layers) :: root_mass    &
                                                   ,soilRT_local &
                                                   ,root_length  &
                                                   ,ratio

    ! reset water flux
    water_flux = 0.0 ; wSWP = 0.0 ; soilRT_local = 0.0 ; soilRT = 0.0
    ratio = 0.0 ; ratio(1) = 1.0
    ! calculate soil depth to which roots reach    
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate scale factor needed to ensure 50 % of root mass in in top 25 %
    ! of rooted layers. See SPA src code for this relationship
!    mult = min(10., max(2.0, 11.*exp(-0.006 * root_biomass)))
!    ! assume surface root biomass density
!    surf_biomass = root_biomass * mult
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    transpiration_resistance = (gplant * lai)**(-1.0)
!    transpiration_resistance = canopy_height / (gplant * lai)

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 2 root layer system this can be estimates more simply

    ! top 25 % of root profile
    slpa = (root_reach * 0.25) - layer_thickness(1) 
    if (slpa <= dble_zero) then
        ! > 50 % of root is in top layer
        root_mass(1) = root_biomass * 0.5
        root_mass(1) = root_mass(1) + ((root_biomass-root_mass(1)) * (abs(slpa)/root_reach))
        root_mass(2) = max(dble_zero,root_biomass - root_mass(1))
    else
        ! < 50 % of root is in bottom layer
        root_mass(1) = root_biomass * 0.5 * (layer_thickness(1)/(abs(slpa)+layer_thickness(1)))
        root_mass(2) = max(dble_zero,root_biomass - root_mass(1))
    endif
    root_length = root_mass / root_mass_length_coef !(root_density * root_cross_sec_area)
    !! Top root layer.
    ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
    soilR1=soil_resistance(root_length(1),min(root_reach,layer_thickness(1)),soil_conductivity(1)*head_1)
    soilR2=root_resistance(root_mass(1),min(root_reach,layer_thickness(1)))
    soilRT_local(1)=soilR1 + soilR2 + transpiration_resistance
    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    ! also note that the head correction is now added rather than
    ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
    ! soilWP prior to application of minlwp
    demand = deltaWP+head*canopy_height
!    water_flux(1) = max(dble_zero,demand/(transpiration_resistance + soilR1 + soilR2))
    water_flux(1) = demand/(transpiration_resistance + soilR1 + soilR2)
    ! Bottom root layer
    if (root_mass(2) > dble_zero ) then
       ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
       soilR1=soil_resistance(root_length(2),layer_thickness(2),soil_conductivity(2)*head_1)
       soilR2=root_resistance(root_mass(2),layer_thickness(2))
       soilRT_local(2)=soilR1 + soilR2 + transpiration_resistance
       ! calculate and accumulate steady state water flux in mmol.m-2.s-1
!       water_flux(2) = max(dble_zero,demand/(transpiration_resistance + soilR1 + soilR2))
       water_flux(2) = demand/(transpiration_resistance + soilR1 + soilR2)
       ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif ! roots present in second layer?
!!----------------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------------
!!-------ORIGINAL SRC BASED ON SPA'S ACTUAL SUBROUTINE------------------------------------
!!----------------------------------------------------------------------------------------
!!       ! use zbrent bisection proceedure to calculate the exponential decay
!!       ! coefficient
!!       x1 = 0.1 ; x2 = 10.0
!!       slpa = zbrent('calculate_Rtot:root_dist', root_dist, x1, x2, xacc)
!!       ! calculate the root mass (gbiomass) in the top 0.1 m of soil (because this is what was
!!       ! done in SPA during training)
!!       prev = 1.0 / slpa
!!       cumdepth = 0.0
!!       do i = 1, nos_root_layers !
!!          cumdepth       = cumdepth + layer_thickness(i)
!!          curr           = 1.0 / slpa * exp(-slpa*cumdepth)
!!          ! don't forget root mass here is biomass. i.e. 2*Croot
!!          root_mass(i)   = (prev - curr) * surf_biomass
!!          ! calculate root length in soil layer (m m-3 soil)
!!          root_length(i) = root_mass(i) / ( root_density * root_cross_sec_area )
!!          prev           = curr
!!          ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
!!          soilC=soil_conductivity(i) / head
!!          soilR1=soil_resistance(root_length(i),layer_thickness(i),soilC)
!!          soilR2=root_resistance(root_mass(i),layer_thickness(i))
!!          ! calculate total resistance to water flow
!!          transpiration_resistance = transpiration_resistance + soilR1 + soilR2 
!!          ! calculate and accumulate steady state water flux in mmol.m-2.s-1
!!          ! NOTE: Depth correction already accounted for in soil resistance
!!          ! calculations and this is the maximum potential rate of transpiration
!!          ! assuming saturated soil and leaves at their minimum water potential.
!!          ! also note that the head correction is now added rather than
!!          ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
!!          ! soilWP prior to application of minlwp
!!          water_flux(i) = max(0.,(deltaWP+head*canopy_height)/(transpiration_resistance))
!!       end do ! both soil layers or just one
!!----------------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------------

    ! if freezing then assume soil surface is frozen
    if (meant < 1.0) then 
        water_flux(1) = 0.0 ; ratio(1) = 0.0
    endif

    ! calculate weighted SWP and uptake fraction
    do soil_layer = 1 , nos_root_layers
       wSWP = wSWP + SWP(soil_layer) * water_flux(soil_layer)
       soilRT = soilRT + soilRT_local(soil_layer) * water_flux(soil_layer)
       ! fraction of total et taken from layer i...
       uptake_fraction(soil_layer) = water_flux(soil_layer) / sum(water_flux)
    enddo
    wSWP = wSWP / sum(water_flux)
    soilRT = soilRT / sum(water_flux)
    if (sum(water_flux) == dble_zero) then
        wSWP = minlwp ; soilRT = sum(soilRT_local)*0.5
        uptake_fraction = 0.0 ; uptake_fraction(1) = 1.0
    endif

    ! determine effective resistance
    Rtot = demand / sum(water_flux*ratio)
    ! finally convert transpiration flux into kg.m-2.step-1 for consistency with
    ! ET in calculate_update_soil_water""
    water_flux = water_flux * 1e-6 * mol_to_g_water * seconds_per_step 

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall)

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    implicit none

    ! arguments 
    double precision, intent(in) :: rainfall ! rainfall (kg.m-2.step-1)

    ! local argumemts
    integer :: i
    double precision    :: add   & ! surface water available for infiltration (m)
                          ,wdiff   ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kg.m-2.step-1 -> Mg.m-2.step-1)
    add = rainfall * 1e-3
    do i = 1 , nos_soil_layers
       ! determine the available pore space in current soil layer
       wdiff=max(dble_zero,(porosity(i)-soil_waterfrac(i))*layer_thickness(i)-watergain(i)+waterloss(i))
       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer
       if (add .gt. wdiff) then
          ! if so fill and subtract from input and move on to the next layer
          watergain(i)=watergain(i)+wdiff
          add=add-wdiff
       else
          ! otherwise infiltate all in the current layer
          watergain(i)=watergain(i)+add
          add = 0.0
       end if
       ! if we have added all available water we are done
       if (add <= 0.0) exit
    end do

    ! if after all of this we have some water left assume it is runoff
    if (add > 0.0) then
       runoff=add * 1e3
    else
       runoff=0.0
    end if

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage(meant)

    ! integrator for soil gravitational drainage !

    implicit none

    ! arguments
    double precision, intent(in) :: meant ! daily mean temperature (oC)

    ! local variables..
    integer,parameter :: nvar = 1
    integer           :: nbad, nok
    double precision  :: change, eps, h1, hmin, newwf, ystart(nvar) &
                        ,iceprop(nos_soil_layers),drainage

    ! --calculations begin below--
!    eps        = 1.0e-4 ! accuracy
!    h1         = 0.001  ! first guess at step size
!    hmin       = 0.0    ! minimum step size
!    kmax       = 100    ! maximum number of iterations (?)
!    x1         = 1.0    ! x1 and x2 define the integrator range for process to
!                        ! occur. i.e. 2-1 = 1, 1 step for integrator
!    x2         = 2.0
!    dxsav      = ( x2 - x1 ) * 0.05

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = 0.0
    ! except the surface layer in the mean daily temperature is < 1oC
    if (meant < 1.0) iceprop(1) = 1.0
    do soil_layer = 1, nos_soil_layers

       ! liquid content of the soil layer, i.e. fraction avaiable for drainage
       liquid     = soil_waterfrac( soil_layer ) * ( 1.0 - iceprop( soil_layer ) )     ! liquid fraction
       ! soil water capacity of the current layer
       drainlayer = field_capacity( soil_layer )
       ! initial conditions; i.e. is there liquid water and more water than
       ! layer can hold
       if ( (liquid > 0.0)  .and. (soil_waterfrac( soil_layer ) > drainlayer) ) then

          ! unsaturated volume of layer below (m3 m-2)..
          unsat      = max( dble_zero , ( porosity( soil_layer+1 ) - soil_waterfrac( soil_layer+1 ) ) &
                                 * layer_thickness( soil_layer+1 ) / layer_thickness( soil_layer ) )

!       ! initial conditions; i.e. is there liquid water and more water than layer
!       ! can hold
!       if ( (liquid .gt. 0.0)  .and. (soil_waterfrac( soil_layer ) .gt. drainlayer) ) then

          ! there is liquid water..
!-------------------------------
!--modified version without ode
!          drainage = calculate_soil_conductivity( soil_waterfrac(soil_layer) ) * seconds_per_step
          drainage = soil_conductivity( soil_layer ) * seconds_per_step
          if ( soil_waterfrac(soil_layer) <= drainlayer ) then  ! gravitational drainage above field_capacity 
             drainage = 0.0
          end if
          if ( drainage >= liquid ) then  ! ice does not drain
             drainage = liquid
          end if
          if ( drainage >= unsat ) then   ! layer below cannot accept more water than unsat
             drainage = unsat
          end if
          change = drainage * layer_thickness(soil_layer)    ! waterloss from this layer
!-------------------------------
!--Original code
!          ystart(1) = soil_waterfrac( soil_layer )    ! total layer   
!          ! call integrator
!          call ode_int( 'gravitational_drainage:soil_water_store' , ystart , nvar , x1 , &
!                         x2 , eps , h1 , hmin , nok , nbad , soil_water_store )
!          newwf  = ystart(1)
!          ! convert from waterfraction to absolute amount (m)
!          change = ( soil_waterfrac( soil_layer ) - newwf ) * layer_thickness( soil_layer )
!-------------------------------
          ! update soil layer below with drained liquid
          watergain( soil_layer + 1 ) = watergain( soil_layer + 1 ) + change
          waterloss( soil_layer     ) = waterloss( soil_layer     ) + change
       end if

    end do ! soil layers

  end subroutine gravitational_drainage
  !
  !-----------------------------------------------------------------
  !
  subroutine soil_porosity(soil_frac_clay,soil_frac_sand)

   ! Porosity is estimated from Saxton equations. !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand
    ! local variables..
    integer :: i
    double precision    :: H, J, K

    ! saxton params relevant to porosity..
    H = 0.332  ;  J = -7.251e-4  ;  K = 0.1276

    ! loop over soil layers..
    do i = 1 , nos_soil_layers
      porosity(i) = H + J * soil_frac_sand(i) + K * log10( soil_frac_clay(i) )
    enddo
    porosity(nos_soil_layers+1) = H + J * soil_frac_sand(nos_soil_layers) + K * log10(soil_frac_clay(nos_soil_layers))

  end subroutine soil_porosity
  !
  !---------------------------------------------------------------------
  !
  subroutine soil_water_store( time_dummy , y , dydt , max_iter)

    ! determines gravitational water drainage !

    implicit none

    ! arguments..
    integer, intent(in) :: max_iter
    double precision,intent(in)  :: y(max_iter)
    double precision,intent(in)  :: time_dummy ! dummy argument, provided for ode_int
    double precision,intent(out) :: dydt(max_iter)

    ! local variables..
    double precision    :: drainage

    drainage = calculate_soil_conductivity( y(1) ) * seconds_per_step
    if ( y(1) .le. drainlayer ) then  ! gravitational drainage above field_capacity 
      drainage = 0.
    end if
    if ( drainage .gt. liquid ) then  ! ice does not drain
      drainage = liquid
    end if
    if ( drainage .gt. unsat ) then   ! layer below cannot accept more water than unsat
      drainage = unsat
    end if
    dydt(1) = -drainage               ! waterloss from this layer

  end subroutine soil_water_store
  !
  !---------------------------------------------------------------------
  !
  subroutine initialise_soils(soil_frac_clay,soil_frac_sand)

    ! 
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none
   
    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay,soil_frac_sand)
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity
    ! calculate initial soil water fraction
    soil_waterfrac = field_capacity
    ! calculate initial soil water potential
    SWP = 0.0
    ! seperately calculate the soil conductivity as this applies to each layer
    do soil_layer = 1, nos_soil_layers
       soil_conductivity(soil_layer) = calculate_soil_conductivity( soil_waterfrac(soil_layer) )
    end do

  end subroutine initialise_soils
  !
  !------------------------------------------------------------------
  !
  subroutine saxton_parameters(soil_frac_clay,soil_frac_sand)

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer :: i
    double precision :: A, B, CC, D, E, F, G, H, J, K, P, Q, R, T, U, V, &
                        mult1, mult2!, mult3

    mult1 = 100.
    mult2 = 2.778e-6
!    mult3 = 1000.0
    A = -4.396    ;  B = -0.0715   ; CC = -4.880e-4 ; D = -4.285e-5
    E = -3.140    ;  F = -2.22e-3  ;  G = -3.484e-5 ; H =  0.332
    J = -7.251e-4 ;  K =  0.1276   ;  P = 12.012    ; Q = -7.551e-2
    R = -3.895    ;  T =  3.671e-2 ;  U = -0.1103   ; V =  8.7546e-4

    do i = 1 , nos_soil_layers
       potA(i)  = exp( A + B * soil_frac_clay(i) + CC * soil_frac_sand(i) * soil_frac_sand(i) + &
                    D * soil_frac_sand(i) * soil_frac_sand(i) * soil_frac_clay(i) ) * mult1
       potB(i)  = E + F * soil_frac_clay(i) * soil_frac_clay(i) + G * soil_frac_sand(i) * soil_frac_sand(i) * soil_frac_clay(i)
       cond1(i) = mult2
       cond2(i) = P + Q * soil_frac_sand(i)
       cond3(i) = R + T * soil_frac_sand(i) + U * soil_frac_clay(i) + V * soil_frac_clay(i) * soil_frac_clay(i)
    enddo
    potA(nos_soil_layers+1)  = exp( A + B * soil_frac_clay(nos_soil_layers) + &
                                   CC * soil_frac_sand(nos_soil_layers) * soil_frac_sand(nos_soil_layers) + &
                                    D * soil_frac_sand(nos_soil_layers) * soil_frac_sand(nos_soil_layers) * &
                                    soil_frac_clay(nos_soil_layers) ) * mult1
    potB(nos_soil_layers+1)  = E + F * soil_frac_clay(nos_soil_layers) * soil_frac_clay(nos_soil_layers) + &
                               G * soil_frac_sand(nos_soil_layers) * soil_frac_sand(nos_soil_layers) * &
                               soil_frac_clay(nos_soil_layers)
    cond1(nos_soil_layers+1) = mult2
    cond2(nos_soil_layers+1) = P + Q * soil_frac_sand(nos_soil_layers)
    cond3(nos_soil_layers+1) = R + T * soil_frac_sand(nos_soil_layers) + U * soil_frac_clay(nos_soil_layers) + & 
                               V * soil_frac_clay(nos_soil_layers) * soil_frac_clay(nos_soil_layers)

  end subroutine saxton_parameters
  !
  !------------------------------------------------------------------
  !
  subroutine soil_water_potential

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    implicit none

    ! local variables..
    integer :: i

    do i = 1 , nos_soil_layers
      if ( soil_waterfrac(i) >= 0.005 ) then
        SWP(i) = -0.001 * potA(i) * soil_waterfrac(i)**potB(i)   !  Soil water potential (MPa)
      else
        SWP(i) = -9999.0
      end if
    enddo

  end subroutine soil_water_potential
  ! 
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(lai)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none
    double precision, intent(in) :: lai
    double precision  min_lai       & ! minimum LAI parameter as height does not vary with growth  
                     ,cd1           & ! canopy drag parameter; fitted to data
                     ,Cr            & ! Roughness element drag coefficient
                     ,Cs            & ! Substrate drag coefficient
                     ,ustar_Uh_max  & ! Maximum observed ratio of (friction velocity / canopy top wind speed) (m.s-1)
                     ,ustar_Uh      & ! ratio of fricition velocity / canopy top wind speed (m.s-1)
!                     ,Cw            & ! characterises roughness sublayer depth (m)
                     ,phi_h           ! roughness sublayer inflence function

    ! set parameters
    cd1=7.5 ; Cs=0.003 ; Cr=0.3 ; ustar_Uh_max=0.3 ; min_lai = 1.0! ; Cw=2.0
    ! assign new value to min_lai to avoid max min calls
    min_lai = max(min_lai,lai)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement=(1.0-((1.0-exp(-sqrt(cd1*min_lai)))/sqrt(cd1*min_lai)))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh=min(sqrt(Cs+Cr*min_lai*0.5),ustar_Uh_max)
    ! calculate roughness sublayer influence function; 
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    phi_h=0.1931472 ! log(Cw)-1.0+Cw**(-1.0) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    ! NOTE that the more empircal 0.13*canopy_height +/-30 % is a fudge but I
    ! can't remember why. Probably due to canopy height not varying with growth
    roughl=((1.0-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

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
  !------------------------------------------------------------------
  ! Functions other than the primary ACM and ACM ET are stored 
  ! below this line.
  !------------------------------------------------------------------
  !
  !------------------------------------------------------------------
  !
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
    linear_model_gradient = ( (interval*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (interval*sumsq_x) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (interval*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !------------------------------------------------------------------
  !
  double precision function calculate_update_soil_water(ET,rainfall,meant)

   !
   ! 1) Limits ET by available water in the soil
   ! 2) Updates soil water balance based on ET, drainage and rainfall
   !

   implicit none

   ! arguments
   double precision, intent(in) :: ET & ! evapotranspiration estimate (kg.m-2.step-1)
                            ,rainfall & ! rainfall (kg.m-2.step-1)
                               ,meant   ! daily mean temperature (oC)

   ! local variables
   double precision ::  depth_change, avail_flux
   double precision, dimension(nos_root_layers) :: evaporation_losses

   ! seperately calculate the soil conductivity as this applies to each layer
   do soil_layer = 1, nos_soil_layers
      soil_conductivity(soil_layer) = calculate_soil_conductivity( soil_waterfrac(soil_layer) )
   end do

   ! for simplicity assume that all evaporation occurs in same distribution as
   ! transpiration. As the surface layer will also have the bulk of the roots
   ! too and steady state flux from surface is linked to water availability too.
   evaporation_losses = ET * uptake_fraction
   ! simplified version of drythick correction on potential evaporation
   ! TLS: removed 25/05/2017 for testing purposes...
!   evaporation_losses(1) = evaporation_losses(1) * min(dble_one,(soil_waterfrac(1) / field_capacity(1)))
   do soil_layer = 1, nos_root_layers
      avail_flux = (soil_waterfrac(soil_layer)*layer_thickness(soil_layer)*1e3)
      if (evaporation_losses(soil_layer) > avail_flux) then 
         ! just to give a buffer against numerical precision error make
         ! extraction slightly less than available
         evaporation_losses(soil_layer) = avail_flux*0.99
      endif
   end do
   where (evaporation_losses < 0) evaporation_losses = 0.0

   ! this will update the ET estimate outside of the function
   ! unit / time correction also occurs outside of this function
   calculate_update_soil_water = sum(evaporation_losses)

   ! pass information to waterloss variable and zero watergain 
   ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
   waterloss = 0.0 ; watergain = 0.0
   waterloss(1:nos_root_layers) = evaporation_losses(1:nos_root_layers)*1e-3 
   ! determine drainage flux between surface -> sub surface and sub surface
   call gravitational_drainage(meant)
   ! determine infiltration from rainfall,
   ! if rainfall is probably liquid / soil surface is probably not frozen
   if (meant >= 1.0 .and. rainfall > 0.0) call infiltrate(rainfall)
   ! update soil profiles. Convert fraction into depth specific values (rather than m3/m3) then update fluxes
   soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness) &
                                        + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                     / layer_thickness(1:nos_soil_layers)

   ! if roots extent down into the bucket 
   if (root_reach > layer_thickness(1)) then
      ! how much has root depth extended since last step?
      depth_change=root_reach-previous_depth
      ! if there has been an increase
      if (depth_change > 0.0) then
          ! calculate weighting between current lowest root layer and new soil 
          depth_change = depth_change / (depth_change+layer_thickness(nos_root_layers))
          ! add to bottom root layer
          soil_waterfrac(nos_root_layers) = (soil_waterfrac(nos_root_layers)*(1.0-depth_change)) & 
                                          + (soil_waterfrac(nos_soil_layers)*depth_change)
      else
          ! calculate weighting between bottom soil layer and new bit coming from lowest root 
          depth_change = abs(depth_change) / (abs(depth_change)+layer_thickness(nos_root_layers))
          ! and add back to the bottom soil layer
          soil_waterfrac(nos_soil_layers) = (soil_waterfrac(nos_soil_layers)*(1.0-depth_change)) &
                                          + (soil_waterfrac(nos_root_layers)*depth_change)
      end if ! depth change 

   end if ! root reach beyond top layer
   ! update new soil states
   previous_depth = root_reach
   ! determine soil layer thickness
   layer_thickness(1) = top_soil_depth ; layer_thickness(2)=max(0.1,root_reach-layer_thickness(1))
   layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
   ! finally update soil water potential
   call soil_water_potential

   ! sanity check for catastrophic failure
!   do soil_layer = 1, nos_soil_layers
!      if (soil_waterfrac(soil_layer) < 0.0 .and. soil_waterfrac(soil_layer) > -0.01) then
!          soil_waterfrac(soil_layer) = 0.0
!      endif
!      if (soil_waterfrac(soil_layer) < 0.0 .or. soil_waterfrac(soil_layer) /= soil_waterfrac(soil_layer)) then
!         print*,'ET',ET,"rainfall",rainfall
!         print*,'evaporation_losses',evaporation_losses
!         print*,"watergain",watergain
!         print*,"waterloss",waterloss
!         print*,'depth_change',depth_change
!         print*,"soil_waterfrac",soil_waterfrac
!         print*,"porosity",porosity
!         print*,"layer_thicknes",layer_thickness
!         print*,"Uptake fraction",uptake_fraction
!         print*,"max_depth",max_depth,"root_k",root_k,"root_reach",root_reach
!         print*,"fail" ; stop
!      endif
!   end do

   ! explicit return needed to ensure that function runs all needed code
   return

  end function calculate_update_soil_water
  !
  !------------------------------------------------------------------
  !
!  double precision function root_dist( xin )    ! parameters
!
!    ! function to calculate the exponential gradient of fine root decay through
!    ! the soil profile. This is aimed at maintaining 50 % of the root mass in
!    ! the top 25 % of the soil profile
!
!    implicit none
!
!    ! arguments..
!    double precision, intent(in) :: xin       ! slope parameter
!
!    ! see pipo chronosequence.xls for this relationship
!    root_dist = ( 1.0 - exp( -xin * root_reach ) ) / xin - root_biomass / surf_biomass
!
!  end function root_dist
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
  !-----------------------------------------------------------------
  !
  double precision function calculate_soil_conductivity( wf )

    ! Used in the soil drainage integrator. !
    ! Returns a single-point value.         !
    ! 'slayer' is a module variable that    !
    !  provides the soil-layer number.      !

    implicit none

    ! arguments..
    double precision, intent(in) :: wf ! fraction of water in soils

    if ( wf .lt. 0.05 ) then    ! Avoid floating-underflow fortran error
      calculate_soil_conductivity = 1e-30
    else
      calculate_soil_conductivity = cond1(soil_layer) * exp( cond2(soil_layer) + cond3(soil_layer) / wf )
      ! Soil conductivity (m s-1 )
    end if

  end function calculate_soil_conductivity
  !
  !-----------------------------------------------------------------
  !
  double precision function soil_resistance(root_length,thickness,soilC)

    !
    ! Calculates the soil hydraulic resistance (MPa m2 s mmol-1) for a given
    ! soil-root zone
    ! 

    implicit none

    ! arguments
    double precision :: root_length, & ! root length in soil layer (m)
                          thickness, & ! thickness of soil layer (m)
                              soilC    ! soil conductivity m2.s-1.MPa-1

    ! local variables
    double precision :: rs, rs2

    ! calculate
    rs  = (root_length*pi)**(-0.5) 
    rs2 = log( rs / root_radius ) / (pi2*root_length*thickness*soilC)
    ! soil water resistance
    soil_resistance = rs2*1e-9*mol_to_g_water

    ! return
    return

  end function soil_resistance
  !
  !------------------------------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    implicit none

    ! arguments..
    double precision, intent(in) :: xin

    ! local variables..
    double precision ::soil_wp

    ! calculate the soil water potential (MPa)..
    ! note that some modifications to scaling values have been made compared to
    ! SPA src to reduce computational cost
!    soil_wp = -0.001 * potA( water_retention_pass ) * xin**potB( water_retention_pass )
!    water_retention_saxton_eqns = -1000.0 * soil_wp + 10.0    ! 10 kPa represents air-entry swp
    soil_wp = potA( water_retention_pass ) * xin**potB( water_retention_pass )
    water_retention_saxton_eqns = -1.0 * soil_wp + 10.0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns
  !
  !------------------------------------------------------------------
  !
  ! 
  !------------------------------------------------------------------
  ! Generic mathematical functions such as bisection and intergrator proceedures
  ! are stored below here
  !------------------------------------------------------------------
  !
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
    double precision,intent(in)             :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       double precision function func( xval )
         double precision ,intent(in) :: xval
       end function func
    end interface

    ! local variables..
    integer            :: iter
    integer,parameter  :: ITMAX = 30
    double precision   :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    double precision,parameter :: EPS = 3e-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0. ) then
      zbrent = a
      return
    elseif ( fb .eq. 0. ) then
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
       tol1 = 2. * EPS * abs(b) + 0.5 * tol
       xm   = .5 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0. ) ) then
          zbrent = b
          return
       end if
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
          s = fb / fa
          if ( a .eq. c ) then
             p = 2. * xm * s
             q = 1. - s
          else
             q = fa / fc
             r = fb / fc
             p = s * ( 2. * xm * q * ( q - r ) - ( b - a ) * ( r - 1. ) )
             q = ( q - 1. ) * ( r - 1. ) * ( s - 1. )
          end if
          if ( p .gt. 0. ) q = -q
          p = abs( p )
          if ( (2.*p) .lt. min( 3.*xm*q-abs(tol1*q) , abs(e*q) ) ) then
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
!  subroutine ode_int( called_from , ystart , nvar , x1 , x2 , eps , h1 , hmin , nok , nbad , derivs )
!
!    ! This is an integrator for ordinary differential equations. We use !
!    ! the RUNGE_KUTTA to track the dynamic behaviour of various model   !
!    ! state variables.. uch a leaf water potential, and water stored on !
!    ! leaf surfaces. The RUNGE_KUTTA finds the time step that ensures   !
!    ! dynamics are smooth and that the relevant feedbacks are properly  !
!    ! incorporated. For a full description see Press et al. (1986).     !
!
!    implicit none
!
!    ! arguments..
!    character(len=*),intent(in)    :: called_from    ! name of procedure calling (used to pass through for errors)
!    integer,intent(in)             :: nvar
!    double precision,intent(in)    :: h1, hmin, x1, x2
!    double precision,intent(inout) :: ystart(nvar), eps
!    integer,intent(out)            :: nbad, nok
!
!    ! Interfaces are the correct way to pass procedures as arguments.
!    ! (note 'derivs' is actually one of
!    ! canopy_water_store,soil_water_store,lwp_diff_eqn)
!!    external :: derivs
!    interface
!       subroutine derivs(time,y,dydt, max_iter)
!         integer, intent(in) :: max_iter
!         double precision,intent(in)    :: time
!         double precision,intent(in)    :: y(max_iter)
!         double precision,intent(out)   :: dydt(max_iter)
!       end subroutine derivs
!    end interface
!
!    ! local variables..
!    integer :: i, kount, nstp
!    double precision :: h, hdid, hnext, x, xsav ,         &
!                        dydx(max_nos_iterations), y(max_nos_iterations), &
!                        yscal(max_nos_iterations), xp(kmaxx), yp(max_nos_iterations,kmaxx)
!
!    ! calculations..
!    x = x1
!    h = sign( h1 , x2-x1 )
!    nok = 0
!    nbad = 0
!    kount = 0
!    do i = 1 , nvar
!       y(i) = ystart(i)
!    enddo
!    if ( kmax .gt. 0 ) xsav = x - 2. * dxsav
!    do nstp = 1 , MAXSTP
!       call derivs( x , y , dydx , max_nos_iterations)
!       do i = 1 , nvar
!          yscal(i) = abs(y(i))+abs(h*dydx(i))+TINY
!       enddo
!       if ( kmax .gt. 0 ) then
!          if ( abs( x - xsav ) .gt. abs( dxsav ) ) then
!             if ( kount .lt. kmax - 1 ) then
!                kount = kount + 1
!                xp(kount) = x
!                do i = 1 , nvar
!                   yp(i,kount) = y(i)
!                enddo
!                xsav = x
!             end if
!          end if
!       end if
!       if ( (x+h-x2) * (x+h-x1) .gt. 0. ) h = x2 - x
!
!       call runge_kutta_q_step( trim(called_from)//":ode_int" , y , dydx , nvar , x , h , eps , yscal , hdid , hnext , derivs )
!       if (hdid.eq.h) then
!          nok = nok+1
!       else
!          nbad = nbad+1
!       end if
!       if ( ( x-x2 ) * ( x2-x1 ) .ge. 0. ) then
!          do i = 1 , nvar
!            ystart(i) = y(i)
!          enddo
!          if ( kmax .ne. 0 ) then
!            kount = kount + 1
!            xp(kount) = x
!            do i = 1 , nvar
!              yp(i,kount) = y(i)
!            enddo
!          end if
!          return
!       end if
!       if ( abs(hnext) .lt. hmin ) then
!         print*,"stepsize smaller than permitted minimum in ode_int",new_line('x'),&
!                          "ode_int was called by: ",trim(called_from)
!       end if
!       h  =  hnext
!    enddo
!
!  end subroutine ode_int
  !
  !------------------------------------------------------------------
  !  
!  subroutine runge_kutta_check( y , dydx , n , x , h , yout , yerr , derivs )
!
!    ! > subroutine summary? < !
!
!    implicit none
!
!    ! arguments..
!    integer,intent(in) :: n
!    double precision,intent(in)    :: h,x,dydx(n),y(n)
!    double precision,intent(out)   :: yerr(n),yout(n)
!
!    ! Interfaces are the correct way to pass procedures as arguments.
!    ! (note 'derivs' is actually one of
!    ! canopy_water_store,soil_water_store,lwp_diff_eqn)
!    external :: derivs
!    interface
!       subroutine derivs( time , y , dydt , max_iter)
!         integer, intent(in) :: max_iter
!         double precision,intent(in)   :: y(max_iter)
!         double precision,intent(in)   :: time
!         double precision,intent(out)  :: dydt(max_iter)
!       end subroutine derivs
!    end interface
!
!    ! local variables..
!    integer        :: i
!    double precision  :: ak2(max_nos_iterations),ak3(max_nos_iterations) &
!                        ,ak4(max_nos_iterations),ak5(max_nos_iterations) &
!                        ,ak6(max_nos_iterations),ytemp(max_nos_iterations)
!    double precision,parameter :: & 
!         A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,B32=9./40., &
!         B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,B53=-70./27.,B54=35./27., &
!         B61=1631./55296.,B62=175./512.,B63=575./13824.,B64=44275./110592., &
!         B65=253./4096.,C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771., &
!         DC1=C1-2825./27648.,DC3=C3-18575./48384.,DC4=C4-13525./55296., &
!         DC5=-277./14336.,DC6=C6-.25
!
!    ! calculations...
!    do i = 1 , n
!       ytemp(i) = y(i) + B21 * h * dydx(i)
!    enddo
!    call derivs( x + A2 * h , ytemp , ak2 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B31 * dydx(i) + B32 * ak2(i) )
!    enddo
!    call derivs( x + A3 * h , ytemp , ak3 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B41 * dydx(i) + B42 * ak2(i) + B43 * ak3(i) )
!    enddo
!    call derivs( x + A4 * h , ytemp , ak4 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B51 * dydx(i) + B52 * ak2(i) + B53 * ak3(i) + B54 * ak4(i) )
!    enddo
!    call derivs( x + A5 * h , ytemp , ak5 , max_nos_iterations )
!    do i = 1 , n
!       ytemp(i) = y(i) + h * ( B61 * dydx(i) + B62 * ak2(i) + B63 * ak3(i) + B64 * ak4(i) + B65 * ak5(i) )
!    enddo
!    call derivs( x + A6 * h , ytemp , ak6 , max_nos_iterations )
!    do i = 1 , n
!       yout(i) = y(i) + h * ( C1 * dydx(i) + C3 * ak3(i) + C4 * ak4(i) + C6 * ak6(i) )
!       yerr(i) = h * ( DC1 * dydx(i) + DC3 * ak3(i) + DC4 * ak4(i) + DC5 * ak5(i) + DC6 * ak6(i) )
!    enddo
!
!  end subroutine runge_kutta_check
  !
  !----------------------------------------------------------------------
  !
!  subroutine runge_kutta_q_step( called_from , y , dydx , n , x , htry , eps , yscal , hdid , hnext , derivs )
!
!    ! > subroutine summary? < !
!
!    implicit none
!
!    ! arguments..
!    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
!    integer,intent(in)          :: n
!    double precision,intent(in)             :: eps, htry, dydx(n), yscal(n)
!    double precision,intent(inout)          :: x
!    double precision,intent(out)            :: hdid, hnext, y(n)
!
!    ! Interfaces are the correct way to pass procedures as arguments.
!    ! canopy_water_store,soil_water_store,lwp_diff_eqn)
!    interface
!       subroutine derivs( time , y , dydt , max_iter)
!         integer, intent(in) :: max_iter
!         double precision,intent(in)  :: y(max_iter)
!         double precision,intent(in)  :: time
!         double precision,intent(out) :: dydt(max_iter)
!       end subroutine derivs
!    end interface
!
!    ! local variables..
!    integer        :: i
!    double precision           :: errmax, h, htemp, xnew, yerr(max_nos_iterations), ytemp(max_nos_iterations)
!    double precision,parameter :: ERRCON = 1.89e-4, & !
!                                  PGROW  = -0.2,    & !
!                                  PSHRNK = -0.25,   & !
!                                  SAFETY = 0.9        !
!
!    ! calculations...
!    h = htry
!1   call runge_kutta_check( y , dydx , n , x , h , ytemp , yerr , derivs )  ! !YES THERE SHOULD BE A ONE AT THE BEGINNING OF THIS LINE
!    errmax = 0.
!    do  i = 1 , n
!       errmax = max( errmax , abs( yerr(i) / yscal(i) ) )
!    enddo
!    errmax = errmax / eps
!    if ( errmax .gt. 1. ) then
!       htemp = SAFETY * h * ( errmax**PSHRNK )
!       h = sign( max( abs(htemp) , 0.1*abs(h) ) , h )
!       xnew = x + h
!       if ( xnew .eq. x ) then
!         print*,"stepsize underflow in runge_kutta_q_step",new_line('x'),&
!                          "runge_kutta_q_step called from: ",trim(called_from)
!       end if
!       goto 1
!    else
!       if ( errmax .gt. ERRCON ) then
!          hnext = SAFETY * h * ( errmax**PGROW )
!       else
!          hnext = 5. * h
!       end if
!       hdid = h
!       x = x + h
!       do i = 1 , n
!          y(i) = ytemp(i)
!       enddo
!
!    end if
!
!  end subroutine runge_kutta_q_step
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEl_MOD
