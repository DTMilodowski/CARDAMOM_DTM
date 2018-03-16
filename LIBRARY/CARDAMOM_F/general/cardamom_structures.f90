module cardamom_structures

implicit none

private

public :: data_type, DATAin, emulator_parameters ,emulator_pars

  !!!!! such as the data type !!!!!
  type DATA_type

      ! drivers
      double precision, allocatable, dimension(:,:) :: MET ! contains our met fields
      double precision ::  meantemp, meanrad ! mean conditions used in some EDCs

      ! OBS: more can obviously be added
      double precision, allocatable, dimension(:) :: GPP     & ! GPP (gC.m-2.day-1)
                                          ,NEE               & ! NEE (gC.m-2.day-1)
                                          ,LAI               & ! LAI (m2/m2)
                                          ,WOO               & ! Wood increment observations (gC.m-2.yr-1)
                                          ,Reco              & ! Ecosystem respiration (gC.m-2.day-1)
                                          ,Cfol_stock        & ! time specific estimate of foliage carbon (gC.m-2)
                                          ,Cwood_stock       & ! time specific estimate of wood carbon (gC.m-2)
                                          ,Croots_stock      & ! time specific estimate of roots carbon (gC.m-2)
                                          ,Csom_stock        & ! time specific estimate if som carbon (gC.m-2)
                                          ,Cagb_stock        & ! time specific agb woody estimate (gC.m-2)
                                          ,Clit_stock        & ! time specific estimate of litter carbon (gC.m-2)
                                          ,Cstem_stock       & ! time specific estimate of stem carbon (gC.m-2)
                                          ,Cbranch_stock     & ! time specific estimate of branch carbon (gC.m-2)
                                          ,Ccoarseroot_stock & ! time specific estimate of coarse root carbon (gC.m-2)
                                          ,Cfolmax_stock     & ! maximum annual foliar stock (gC.m-2)
                                          ,Evap                ! Evapotranspiration (kg.m-2.day-1)

      ! OBS uncertainties: obv these must be pared with OBS above
      double precision, allocatable, dimension(:) :: GPP_unc     & ! (gC.m-2.day-1)
                                          ,NEE_unc               & ! (gC.m-2.day-1)
                                          ,LAI_unc               & ! (log(m2/m2))
                                          ,WOO_unc               & ! (log(gC.m-2.yr-1))
                                          ,Reco_unc              & ! (gC.m-2.day-1)
                                          ,Cfol_stock_unc        & ! (%)
                                          ,Cwood_stock_unc       & ! (%)
                                          ,Croots_stock_unc      & ! (%)
                                          ,Csom_stock_unc        & ! (%)
                                          ,Cagb_stock_unc        & ! (%)
                                          ,Clit_stock_unc        & ! (%)
                                          ,Cstem_stock_unc       & ! (%)
                                          ,Cbranch_stock_unc     & ! (%)
                                          ,Ccoarseroot_stock_unc & ! (%)
                                          ,Cfolmax_stock_unc     & ! (%)
                                          ,Evap_unc                ! (kg.m-2.day-1)

      ! location of observations in the data stream
      integer, allocatable, dimension(:) :: gpppts               & ! gpppts vector used in deriving ngpp
                                           ,neepts               & ! same for nee
                                           ,woopts               & ! same for wood
                                           ,laipts               & ! same for lai
                                           ,recopts              & ! same for ecosystem respiration
                                           ,Cfol_stockpts        & ! same for Cfoliage
                                           ,Cwood_stockpts       & ! smae for Cwood
                                           ,Croots_stockpts      & ! same for Croots
                                           ,Csom_stockpts        & ! same for Csom
                                           ,Cagb_stockpts        & ! same 
                                           ,Clit_stockpts        & ! same for Clitter
                                           ,Cstem_stockpts       & ! same for Csom
                                           ,Cbranch_stockpts     & ! same 
                                           ,Ccoarseroot_stockpts & ! same for Clitter
                                           ,Cfolmax_stockpts     & !
                                           ,Evappts

      ! counters for the number of observations per data stream
      integer :: ngpp               & ! number of GPP observations
                ,nnee               & ! number of NEE observations
                ,nlai               & ! number of LAI observations
                ,nwoo               & ! number of wood increment obervations
                ,nreco              & ! number of Reco observations
                ,nCfol_stock        & ! number of Cfol observations
                ,nCwood_stock       & ! number of Cwood observations
                ,nCroots_stock      & ! number of Croot observations
                ,nCsom_stock        & ! number of Csom obervations
                ,nCagb_stock        & !
                ,nClit_stock        & ! number of Clitter observations
                ,nCstem_stock       & ! 
                ,nCbranch_stock     & !
                ,nCcoarseroot_stock & !
                ,nCfolmax_stock     & !
                ,nEvap

      ! saving computational speed by allocating memory to model output
      double precision, allocatable, dimension(:) :: M_GPP    & ! 
                                          ,M_NEE    & !
                                          ,M_LAI      ! 
      ! timing variable
      double precision, allocatable, dimension(:) :: deltat ! time step (decimal day)

      double precision, allocatable, dimension(:,:) :: M_FLUXES & !
                                                      ,M_POOLS  & !
                                                      ,C_POOLS    !
      ! static data
      integer :: nodays   & ! number of days in simulation
                ,ID       & ! model ID, currently 1=DALEC_CDEA, 2=DALEC_BUCKET
                ,noobs    & ! number of obs fields
                ,nomet    & ! number met drivers
                ,nofluxes & ! number of fluxes
                ,nopools  & ! number of pools
                ,nopars   & ! number of parameters
                ,EDC      & ! Ecological and dynamical contraints on (1) or off (0)
                ,yield    & ! yield class for ecosystem (forest only)
                ,age      & ! time in years since ecosystem established (forest only)
                ,pft        ! plant functional type information used to select appropriate DALEC submodel / ACM

      double precision :: LAT ! site latitude

      ! binary file mcmc options (need to add all options HERE except
      ! inout files)
      integer :: edc_random_search ! 

      ! priors
      double precision, dimension(100) :: parpriors     & ! prior values
                                         ,parpriorunc   & ! prior uncertainties
                                         ,otherpriors   & ! other prior values
                                         ,otherpriorunc   ! other prior uncertainties

  end type ! DATA_type
  type (DATA_type), save :: DATAin

  type emulator_parameters

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

  end type ! emulator parameters
  type (emulator_parameters), save :: emulator_pars

end module cardamom_structures
