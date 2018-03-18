
module model_likelihood_module
  implicit none

  ! make all private
  private

  ! which to make open
  public :: model_likelihood, find_edc_initial_values

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: PASSFAIL(100) ! allow space for 100 possible checks
    integer :: EDC
    integer :: DIAG
    integer :: nedc ! number of edcs being assessed
  end type
  type (EDCDIAGNOSTICS), save :: EDCD

  contains
  !
  !------------------------------------------------------------------
  !
  subroutine find_edc_initial_values (PI)
    use MCMCOPT, only: mcmc_output, parameter_info, mcmc_options, initialise_mcmc_output
    use cardamom_structures, only: DATAin ! will need to change due to circular dependance
    use cardamom_io, only: restart_flag
    use MHMCMC_MODULE, only: MHMCMC

    ! subroutine deals with the determination of initial parameter and initial
    ! conditions which are consistent with EDCs

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI

    ! declare local variables
    type ( mcmc_output ) :: MCOUT
    type ( mcmc_options ) :: MCOPT
    integer :: n, counter_local
    double precision :: PEDC

    ! initialise output for this EDC search
    call initialise_mcmc_output(PI,MCOUT)

    ! set MCMC options needed for EDC run
    MCOPT%APPEND=0
    MCOPT%nADAPT=20
    MCOPT%fADAPT=0.5d0
    MCOPT%nOUT=5000
    MCOPT%nPRINT=0
    MCOPT%nWRITE=0
    ! the next two lines ensure that parameter inputs are either given or
    ! entered as -9999
    MCOPT%randparini=.true.
    MCOPT%returnpars=.true.
    MCOPT%fixedpars=.false.

    do n = 1, PI%npars
       PI%stepsize(n)=0.02d0
       PI%parini(n)=DATAin%parpriors(n)
       ! assume we need to find random parameters
       PI%parfix(n)=0
       ! if the prior is not missing and we have not told the edc to be random
       ! keep the value
       if (PI%parini(n) /= -9999d0 .and. DATAin%edc_random_search < 1) PI%parfix(n)=1
    end do ! parameter loop

    if (.not. restart_flag) then
       ! set up edc log likelihood for MHMCMC initial run
       PEDC = -1
       counter_local = 0
       do while (PEDC < 0)
         write(*,*)"Beginning EDC search attempt"
         ! reset the parameter step size at the beginning of each attempt
         PI%stepsize(1:PI%npars) = 0.02d0 !0.0005
         ! call the MHMCMC directing to the appropriate likelihood
         call MHMCMC(EDC_MODEL_LIKELIHOOD,PI,MCOPT,MCOUT)

         ! store the best parameters from that loop
         PI%parini(1:PI%npars)=MCOUT%best_pars(1:PI%npars)
         ! turn off random selection for initial values
         MCOPT%randparini = .false.

         ! call edc likelihood function to get final edc probability
         call edc_model_likelihood(PI,PI%parini,PEDC)

         ! keep track of attempts
         counter_local=counter_local+1
         ! periodically reset the initial conditions
         if (PEDC < 0 .and. mod(counter_local,3) == 0) then
             PI%parini(1:PI%npars)=DATAin%parpriors(1:PI%npars)
             ! reset to select random starting point
             MCOPT%randparini = .true.
         endif
       end do ! for while condition
    endif ! if for restart

    ! reset
    PI%parfix(1:PI%npars)=0

    ! clean up some memory
    deallocate(MCOUT%best_pars)

  end subroutine find_edc_initial_values
  !
  !------------------------------------------------------------------
  !
  subroutine edc_model_likelihood(PI, PARS, prob_out)
!TLS    use, intrinsic :: ieee_arithmetic
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PARAMETER_INFO
    use CARBON_MODEL_MOD, only: carbon_model
    use CARBON_MODEL_CROP_MOD, only: CARBON_MODEL_CROP

    ! Model likelihood function specifically intended for the determination of
    ! appropriate initial parameter choices, consistent with EDCs for DALEC2 /
    ! DALEC_GSI

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI
    double precision, dimension(PI%npars), intent(inout) :: PARS
    ! output
    double precision, intent(inout) :: prob_out

    ! declare local variables
    integer ::  n
    double precision :: tot_exp, ML, exp_orig, decay_coef, EDC1, EDC2, infini

    ! set initial values
    EDCD%DIAG=1

    if (DATAin%PFT == 1) then
       ! PFT has been provided and is crop! Best try running the crop model
       ! then...

       ! call EDCs which can be evaluated prior to running the model
       call EDC1_CROP(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
       ! next need to run the model itself
       call CARBON_MODEL_CROP(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                             ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE  &
                             ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                             ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                             ,DATAin%nofluxes,DATAin%M_GPP                &
                             ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                             ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                             ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

        ! assess post running EDCs
        call EDC2_CROP(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                      ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                      ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                      ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

    else
       ! aaahhh relax we are using the generic DALEC model this time so things
       ! should move smoothly...
       ! call EDCs which can be evaluated prior to running the model
       call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

       ! next need to run the model itself
       call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat&
                        ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                        ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                        ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                        ,DATAin%M_GPP)

       ! assess post running EDCs
       call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                    ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                    ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                    ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

    end if ! crop or not if

    ! calculate the likelihood
    tot_exp = sum(1d0-EDCD%PASSFAIL(1:EDCD%nedc))
!    tot_exp=0d0
!    do n = 1, EDCD%nedc
!       tot_exp=tot_exp+(1d0-EDCD%PASSFAIL(n))
!       if (EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
!    end do ! checking EDCs
!    ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) then
!        print*,"Found it!" ; stop
!    endif

    ! convert to a probability
    prob_out = -0.5d0*(tot_exp*10d0)*DATAin%EDC

    ! override probability if parameter set gives NaN or near -infinitiy output
    call model_likelihood(PI,PARS,ML)

    infini = 0d0
    if (ML /= ML .or. ML == log(infini) .or. ML == -log(infini) .or. &
        sum(DATAin%M_LAI) /= sum(DATAin%M_LAI) .or. sum(DATAin%M_GPP) /= sum(DATAin%M_GPP)) then
       prob_out=prob_out-0.5d0*10d0
    end if

    ! now add the exponential component
    ! prob_out is the Log-Likelihood
    prob_out=prob_out

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine EDC1_CROP (PARS, npars, meantemp, meanrad, EDC1)
    ! the first of two subroutine to assess current parameters for passing
    ! realism tests for crop
    ! ecosystems

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local variables
    integer :: n, DIAG
    double precision :: torfol,tmp ! yearly leaf loss fraction

    ! set initial value
!    torfol=1./(pars(5)*365.25)
    EDC1 = 1
    DIAG = EDCD%DIAG

    ! set all EDCs to 1 (pass)
    EDCD%nedc = 100
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(10) > pars(9))) then
        EDC1 = 0 ; EDCD%PASSFAIL(1) = 0
    endif

    ! decomposition of litter to SOM greater than SOM to air
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(10) > pars(1))) then
        EDC1 = 0 ; EDCD%PASSFAIL(2) = 0
    endif

    ! turnover of foliage faster than turnover of wood
! TLS: turnover off because foliage and stem turnovers are made same
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(6) > pars(5)) then
!       EDC1=0 ; EDCD%PASSFAIL(3)=0
!    end if

    ! pre_DR should be greater than post_DR
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(4) > pars(3))) then
!        EDC1=0 ; EDCD%PASSFAIL(4)=0
!    endif

    ! for development: Tmin should be < topt and topt should be < tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(26) > pars(28) &
                                     .or. pars(28) > pars(27) &
                                     .or. pars(26) > pars(27))) then
        EDC1=0 ; EDCD%PASSFAIL(5)=0
    endif

    ! for development: the difference between each Tmin,Topt,Tmax > 1.
    if ((EDC1 == 1 .or. DIAG == 1) .and. (abs(pars(26)-pars(28)) < 1d0 &
                                     .or. abs(pars(28)-pars(27)) < 1d0  &
                                     .or. abs(pars(26)-pars(27)) < 1d0)) then
        EDC1=0 ; EDCD%PASSFAIL(6)=0
    endif

   ! for vernalisation: Tmin < Topt < Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(31) &
                                     .or. pars(31) > pars(30) &
                                     .or. pars(29) > pars(30))) then
        EDC1=0 ; EDCD%PASSFAIL(7)=0
    endif

   ! for vernalisation: the difference between each Tmin, Topt, Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. ( abs(pars(29)-pars(31)) < 1d0 &
                                      .or. abs(pars(31)-pars(30)) < 1d0 &
                                      .or. abs(pars(29)-pars(30)) < 1d0 ) ) then
        EDC1=0 ; EDCD%PASSFAIL(8)=0
    endif

   ! development temperature value should be larger corresponding vernalisation
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(26) &
                                     .or. pars(31) > pars(28) &
                                     .or. pars(30) > pars(27))) then
        EDC1=0 ; EDCD%PASSFAIL(9)=0
    endif

!    ! plough must be before sow and after harvest: WINTER ONLY
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(16) > pars(12) ) then
!        EDC1=0 ; EDCD%PASSFAIL(10)=0
!    endif
!
!    ! harvest cannot be more than 345 after harvest
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(15) < pars(12)+345.25) ) then
!        EDC1=0 ; EDCD%PASSFAIL(11)=0
!    endif
!
!    ! plough must be before sow and after harvest: WINTER ONLY
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(16) < pars(15) .or. &
!        pars(16) > pars(12)) ) then
!        EDC1=0 ; EDCD%PASSFAIL(12)=0
!    endif

    ! CN ratio of leaf should also be between 95CI of trait database values
    ! Kattge et al (2011)
    tmp = (pars(17)/(10d0**pars(11)))
    if ((EDC1 == 1 .or. DIAG == 1) .and. (tmp > 43.76895d0 .or. tmp < 10.82105d0)) then
       EDC1 = 0 ; EDCD%PASSFAIL(13) = 0
    endif

    ! could and probably should add some more
  end subroutine EDC1_CROP
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_CROP(npars,nomet,nofluxes,nopools,nodays,deltat &
                      ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                      ,meantemp,EDC2)

    use CARBON_MODEL_CROP_MOD, only: resp_rate_temp_coeff, ts_length 
    use CARBON_MODEL_MOD,      only: seconds_per_hour, seconds_per_day

    ! the second of two subroutines for assessing current parameters for passing
    ! realism tests for crop ecosystems

    implicit none

    ! declare input variables
    integer, intent(in) :: npars    & ! number of model parameters
                          ,nomet    & ! number of met drivers
                          ,nofluxes & ! number of fluxes from model
                          ,nopools  & ! number of pools in model
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_LAI(nodays)               & ! LAI output from current model simulation
                                   ,M_NEE(nodays)               & ! NEE output from current model simulation
                                   ,M_GPP(nodays)               & ! GPP output from current model simulation
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation model
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs

    ! declare local variables
    integer :: n, DIAG, no_years, nn
    double precision :: mean_pools(nopools), G, decay_coef, meangpp, EQF, PEDC, infi
    double precision, dimension(:), allocatable :: mean_annual_pools
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,flit  & !
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    ! set initial value
    fauto=sum(M_FLUXES(:,3))/sum(M_FLUXES(:,1))
    ffol=sum(M_FLUXES(:,4))/(sum(M_FLUXES(:,1))*fauto)
    flab=sum(M_FLUXES(:,5))/(sum(M_FLUXES(:,1))*fauto)
    froot=sum(M_FLUXES(:,6))/(sum(M_FLUXES(:,1))*fauto)
    fwood=sum(M_FLUXES(:,7))/(sum(M_FLUXES(:,1))*fauto)
    fsom=fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(10))
    flit=(froot+flab+ffol)
    ! length of time step in hours..
!    ts_length = ((sum(deltat)/(nodays)) * sec_in_day) / sec_in_hour
    ts_length = 24d0
    ! update initial values
    DIAG=EDCD%DIAG
    ! give EDC2 an initial value
    EDC2 = 1

    ! SOM attractor - must be within a factor of 2 from Csom0
    ! eqiulibrium factor (in comparison with initial conditions)
    EQF=10d0

    ! initialise and then calculate mean gpp values
    meangpp=sum(M_GPP(1:nodays))/dble(nodays)

    ! EDC 11 - SOM steady state within order magnitude of initial conditions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(pars(10)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) > (pars(23)*EQF)) then
       EDC2 = 0 ; EDCD%PASSFAIL(14) = 0
    end if
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(pars(10)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) < (pars(23)/EQF)) then
       EDC2 = 0 ; EDCD%PASSFAIL(15) = 0
    endif

    ! EDC 12 - Litter steady state assumptions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(pars(9)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) > (pars(22)*EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(16) = 0
    endif
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(pars(9)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) < (pars(22)/EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(17) = 0
    endif

    ! EDC 13
    ! assesses the exponential decay/growth of the Csom pool

    !  work out how many completed years there are in the system
    no_years=int(nint(sum(deltat)/365.25d0))

    ! only do this for the Csom pool
    do n = 1, 1 !nopools
       if (EDC2 == 1 .or. DIAG == 1) then
          decay_coef=expdecay2(M_POOLS(1:(nodays+1),6),n,deltat,1,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
          if (abs(-log(2d0)/decay_coef) < (365.25d0*dble(no_years)) .and. decay_coef < 0d0 ) then
             EDC2 = 0 ; EDCD%PASSFAIL(18)=0
          end if ! EDC conditions
       end if ! EDC .or. DIAG condition
    end do ! pools loop

    ! EDC 14
    ! assesses the exponential decay/growth of the Clit pool

    ! only do this for the Clit pool
    do n = 1, 1 !nopools
       if (EDC2 == 1 .or. DIAG == 1) then
          decay_coef=expdecay2(M_POOLS(1:(nodays+1),5),n,deltat,1,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
          if (abs(-log(2d0)/decay_coef) < (365.25d0*dble(no_years)) .and. decay_coef < 0d0 ) then
             EDC2 = 0 ; EDCD%PASSFAIL(19)=0
          end if ! EDC conditions
       end if ! EDC .or. DIAG condition
    end do ! pools loop

    ! we know that the crop areas should produce some yield - therefore reject
    ! parmeter sets which not produce any yield ever!
    if ((EDC2 == 1 .or. DIAG == 1) .and. sum(M_FLUXES(1:nodays,17)) < 1d0) then
        EDC2 = 0 ; EDCD%PASSFAIL(20)=0
    endif

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 50-100 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
    infi = 0d0
       n=1
       do while (n <= nopools .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_POOLS(nn,n) < 0d0 .or. M_POOLS(nn,n) /= M_POOLS(nn,n) .or. &
                 M_POOLS(nn,n) == log(infi) .or. M_POOLS(nn,n) == -log(infi)) then
                 EDC2=0 ; PEDC=0 ; EDCD%PASSFAIL(50+n)=0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
       n = 1
       do while (n <= nofluxes .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_FLUXES(nn,n) < 0d0 .or. M_FLUXES(nn,n) /= M_FLUXES(nn,n) .or. &
                 M_FLUXES(nn,n) == log(infi) .or. M_FLUXES(nn,n) == -log(infi)) then
                 EDC2=0 ; PEDC=0 ; EDCD%PASSFAIL(50+nopools+n)=0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition

    end if ! min pool assessment

  end subroutine EDC2_CROP
  !
  !------------------------------------------------------------------
  !
  subroutine EDC1_GSI(PARS, npars, meantemp, meanrad, EDC1)

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (modified from Bloom et al., 2014).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local variables
    integer :: n, DIAG
    double precision :: tmp

    ! set initial value
    EDC1=1
    DIAG=EDCD%DIAG

    ! set all EDCs to 1 (pass)
    EDCD%nedc=100
    EDCD%PASSFAIL(1:EDCD%nedc)=1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9) > pars(8))) then
        EDC1=0 ; EDCD%PASSFAIL(1)=0
    endif

    ! turnover of cwd should be slower than litter
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(8)+pars(1)) < pars(38) ) then
        EDC1=0 ; EDCD%PASSFAIL(2)=0
    endif

    ! turnover of cwd should be faster than wood
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(6) > (pars(38)*exp(pars(10)*meantemp))) ) then
        EDC1=0 ; EDCD%PASSFAIL(3)=0
    endif

    ! litter2som greater than som to atm rate
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(1) < pars(9))) then
       EDC1=0 ; EDCD%PASSFAIL(4)=0
    endif

    ! root turnover greater than som turnover at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(7) < (pars(9)*exp(pars(10)*meantemp)))) then
       EDC1=0 ; EDCD%PASSFAIL(5)=0
    endif

    ! turnover of roots should be faster than that of wood
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(7) < pars(6)) then
       EDC1=0 ; EDCD%PASSFAIL(6)=0
    endif

    ! replanting 30 = labile ; 31 = foliar ; 32 = roots ; 33 = wood
    ! initial    18 = labile ; 19 = foliar ; 20 = roots ; 21 = wood
    ! initial replanting labile must be consistent with available wood storage
    ! space. Labile storage cannot be greater than 12.5 % of the total ecosystem carbon stock.
    ! Gough et al (2009) Agricultural and Forest Meteorology. Avg 11, 12.5, 3 %
    ! (Max across species for branch, bole and coarse roots). Evidence that
    ! Branches accumulate labile C prior to bud burst from other areas.
    ! Wurth et al (2005) Oecologia, Clab 8 % of living biomass (DM) in tropical
    ! forest Richardson et al (2013), New Phytologist, Clab 2.24 +/- 0.44 % in
    ! temperate (max = 4.2 %)
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(30) > (pars(33)*0.125d0)) then
       EDC1 = 0 ; EDCD%PASSFAIL(7) = 0
    endif
    ! also apply to initial conditions
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(18) > (pars(21)*0.125d0)) then
       EDC1 = 0 ; EDCD%PASSFAIL(8) = 0
    endif

    ! initial replanting foliage and fine roots ratio must be consistent with
    ! ecological ranges. Because this is the initial condition and not the mean
    ! only the upper foliar:fine root bound is applied
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(32)/pars(31) < 0.04d0) ) then
       EDC1 = 0 ; EDCD%PASSFAIL(9) = 0
    endif
    ! also apply to initial conditions
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(20)/pars(19) < 0.04d0) ) then
       EDC1 = 0 ; EDCD%PASSFAIL(10) = 0
    endif

    ! replanting stock of foliage is unlikely to have much lai, thus limit lai
    ! to less than 1.5 m2/m2
    if ((EDC1 == 1 .or. DIAG == 1) .and. ((pars(31)/pars(17)) > 1.5d0)) then
       EDC1 = 0 ; EDCD%PASSFAIL(11) = 0
    endif

    ! initial replanting wood stocks must be sufficient to support intended
    ! foliar stocks. Again as this is the initial values and not the annual mean /
    ! maximum stocks we constrain against the minimum foliar:wood stock only.
    ! We assume that at least as much structural C is needed to hold a given
    ! leaf mass.
    ! NOTE: only half that used from Thomas & Williams to allow for non-forested
    ! systems
    if ((EDC1 == 1 .or. DIAG == 1) .and. ((pars(31) / pars(33)) > 2.0d0) ) then
       EDC1 = 0 ; EDCD%PASSFAIL(12) = 0
    endif
    ! also apply to initial conditions
    if ((EDC1 == 1 .or. DIAG == 1) .and. ((pars(19) / pars(21)) > 2.0d0) ) then
       EDC1 = 0 ; EDCD%PASSFAIL(13) = 0
    endif

    ! --------------------------------------------------------------------
    ! TLS: some added specifically to deal with GSI conditions
    ! Note that the EDC numbers do not run on
    ! avgTmin min threshold should not be larger than max
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(14) > pars(15)) ) then
       EDC1=0 ; EDCD%PASSFAIL(14)=0
    endif

    ! photoperiod, min threshold should not be larger than max
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(16) > pars(24)) ) then
       EDC1=0 ; EDCD%PASSFAIL(15)=0
    endif

    ! VPD min threshold should not be larger than max
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(25) > pars(26)) ) then
       EDC1=0 ; EDCD%PASSFAIL(16)=0
    endif

    ! critical values for labile must be greater than foliage
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(34) > pars(3))) then
       EDC1=0 ; EDCD%PASSFAIL(17)=0
    endif

    ! CN ratio of leaf should also be between 95CI of trait database values
    ! Kattge et al (2011)
    tmp = (pars(17)/(10d0**pars(11)))
    if ((EDC1 == 1 .or. DIAG == 1) .and. (tmp > 43.76895d0 .or. tmp < 10.82105d0)) then
       EDC1 = 0 ; EDCD%PASSFAIL(18) = 0
    endif

    ! --------------------------------------------------------------------
    ! could always add more / remove some

  end subroutine EDC1_GSI
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_GSI(npars,nomet,nofluxes,nopools,nodays,deltat &
                      ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                      ,meantemp,EDC2)

    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: disturbance_residue_to_litter, &
                                disturbance_residue_to_som,    &
                                disturbance_residue_to_cwd,    &
                                disturbance_loss_from_cwd,     &
                                disturbance_loss_from_litter,  &
                                disturbance_loss_from_som

    ! Determines whether the dynamical contraints for the search of the initial
    ! parameters has been successful or whether or not we should abandon the
    ! current set and move on

    implicit none

    ! declare input variables
    integer, intent(in) :: npars    & ! number of model parameters
                          ,nomet    & ! number of met drivers
                          ,nofluxes & ! number of fluxes from model
                          ,nopools  & ! number of pools in model
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_LAI(nodays)               & ! LAI output from current model simulation
                                   ,M_NEE(nodays)               & ! NEE output from current model simulation
                                   ,M_GPP(nodays)               & ! GPP output from current model simulation
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation model
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs

    ! declare local variables
    logical :: found
    integer :: n, DIAG, no_years, y, PEDC, nn, num_EDC, max_location(1) &
              ,i, exp_adjust, no_years_adjust, disturb_year, replant_year &
              ,disturb_begin, disturb_end
    double precision :: mean_pools(nopools), G, decay_coef, meangpp &
                       ,EQF2, EQF5, EQF10, EQF20, sumgpp, sumnpp, model_living_C &
                       ,target_living_C(2),hold, steps_per_year
    double precision, dimension(nodays) :: mean_ratio, resid_fol,resid_lab
    integer, dimension(nodays) :: hak ! variable to determine number of NaN in foliar residence time calculation
    double precision, dimension(:), allocatable :: mean_annual_pools,tmp
    double precision :: max_wood    & !
               ,in_out_root_disturb &
               ,in_out_wood_disturb &
                       ,in_out_root &
                       ,in_out_wood &
                       ,in_out_lit  &
                       ,in_out_cwd  &
                       ,in_out_som  &
                       ,in_out_dead &
                       ,torfol      & ! yearly average turnover
                       ,torlab      & !
                       ,fNPP        & ! fraction of NPP to foliage
                       ,rNPP        & ! fraction of NPP to roots
                       ,wNPP        & ! fraction of NPP to wood
                       ,fauto       & ! fraction of GPP to autotrophic respiration
                       ,ffol        & ! fraction of GPP to foliage
                       ,flab        & ! fraction of GPP to labile pool
                       ,froot       & ! fraction of GPP to root
                       ,fwood       & ! fraction of GPP to wood
                       ,delta_gsi

    ! update initial values
    hak = 0 ; resid_fol = 0d0
    ! calculate mean turnover rate for leaves
    resid_fol(1:nodays) = (M_FLUXES(1:nodays,10)+M_FLUXES(1:nodays,22))/M_POOLS(1:nodays,2)
    ! division by zero results in NaN plus obviously I can't have turned
    ! anything over if there was nothing to start out with...
    where ( M_POOLS(1:nodays,2) == 0d0 )
           hak = 1 ; resid_fol(1:nodays) = 0d0
    end where
    torfol = sum(resid_fol) / dble(nodays-sum(hak))

    hak = 0 ; resid_lab = 0d0
    ! calculate mean turnover rate for labile pool
    resid_lab(1:nodays) = (M_FLUXES(1:nodays,8)+M_FLUXES(1:nodays,21))/M_POOLS(1:nodays,1)
    ! division by zero results in NaN plus obviously I can't have turned
    ! anything over if there was nothing to start out with...
    where ( M_POOLS(1:nodays,1) == 0d0 )
           hak = 1 ; resid_lab(1:nodays) = 0d0
    end where
    torlab = sum(resid_lab) / dble(nodays-sum(hak))

    ! timing information needed
    no_years=int(nint(sum(deltat)/365.25d0))
    steps_per_year = sum(deltat)/dble(no_years)
    no_years_adjust=no_years
    allocate(mean_annual_pools(no_years))

    ! Some EDCs can only be used if the management periods are except in their
    ! analysis timeframe. For example EDC 8 assesses expoential shifts
    found = .false. ; exp_adjust = 1 ; disturb_begin = 1 ; disturb_end = 1
    if (maxval(met(8,:)) > 0.99d0 ) then
       ! so we will find the location of the management
       i = 0
       do while (.not.found)
          i = i + 1
          ! if we find what we are looking for
          if (met(8,i) > 0.99d0 .or. i == nodays) found = .true.
       enddo
       disturb_begin = i-1 ; disturb_end = i + nint((real(nodays)/real(no_years))*2d0)
       ! if the end is more than 1 year away we are good to go. Otherwise bail
       ! on the second half of the EDC by setting disturb_end == nodays
       if ((nodays-disturb_end) < nint((real(nodays)/real(no_years)))) then
          disturb_end = nodays
       endif
       ! check if this is in the first year
       if (sum(deltat(1:i)) < (365.25d0*2d0)) then
          ! if so then we need to calculate the adjustment
          exp_adjust=i
       endif
       ! calculate new number of whole years to assess over
       no_years_adjust=int(nint(sum(deltat(exp_adjust:nodays))/365.25d0))
    endif

    ! initialise and then calculate mean gpp values
    fauto=sum(M_FLUXES(1:nodays,3))/sum(M_FLUXES(1:nodays,1))
!    meangpp=sum(M_GPP(1:nodays))/real(nodays)
    sumgpp=1/sum(M_GPP(1:nodays))
    sumnpp=1/(sum(M_GPP(1:nodays))*(1d0-fauto))

    DIAG=EDCD%DIAG
    EDC2=1

    ! GPP allocation fractions
    ffol = sum(M_FLUXES(1:nodays,4)+M_FLUXES(1:nodays,8)) * sumgpp
    froot = sum(M_FLUXES(1:nodays,6)) * sumgpp
    fwood = sum(M_FLUXES(1:nodays,7)) * sumgpp

    ! NPP allocations; note that because of possible labile accumulation this
    ! might not be equal to 1
    fNPP = sum(M_FLUXES(1:nodays,4)+M_FLUXES(1:nodays,8)) * sumnpp
    wNPP = sum(M_FLUXES(1:nodays,7)) * sumnpp
    rNPP = sum(M_FLUXES(1:nodays,6)) * sumnpp

    ! calculate input and output ratios for all pools
    if (maxval(met(8,:)) > 0.99d0 .and. disturb_end == nodays) then
       ! there has been a replacement level event, but there is less than 2
       ! years before the end so we will assess the beginning of the analysis
       ! only
       in_out_root = sum(M_FLUXES(1:disturb_begin,6)) / sum(M_FLUXES(1:disturb_begin,12)+M_FLUXES(1:disturb_begin,23))
       in_out_root_disturb = 1d0
       in_out_wood = sum(M_FLUXES(1:disturb_begin,7)) / sum(M_FLUXES(1:disturb_begin,11)+M_FLUXES(1:disturb_begin,24))
       in_out_wood_disturb = 1d0
       in_out_lit = sum(M_FLUXES(1:disturb_begin,10)+M_FLUXES(1:disturb_begin,12)+M_FLUXES(1:disturb_begin,19) &
                        +disturbance_residue_to_litter(1:disturb_begin)) &
                  / sum(M_FLUXES(1:disturb_begin,13)+M_FLUXES(1:disturb_begin,15)+disturbance_loss_from_litter(1:disturb_begin))
       in_out_som = sum(M_FLUXES(1:disturb_begin,15)+disturbance_residue_to_som(1:disturb_begin)) &
                  / sum(M_FLUXES(1:disturb_begin,14)+disturbance_loss_from_som(1:disturb_begin))
       in_out_cwd = sum(M_FLUXES(1:disturb_begin,11)+disturbance_residue_to_cwd(1:disturb_begin)) &
                  / sum(M_FLUXES(1:disturb_begin,19)+disturbance_loss_from_cwd(1:disturb_begin))
       in_out_dead = sum(M_FLUXES(1:disturb_begin,10)+M_FLUXES(1:disturb_begin,12) &
                        +M_FLUXES(1:disturb_begin,11)                       &
                        +disturbance_residue_to_som(1:disturb_begin)        &
                        +disturbance_residue_to_cwd(1:disturb_begin)        &
                        +disturbance_residue_to_litter(1:disturb_begin))    &
                   / sum(M_FLUXES(1:disturb_begin,13)+M_FLUXES(1:disturb_begin,14) &
                        +disturbance_loss_from_som(1:disturb_begin)         &
                        +disturbance_loss_from_cwd(1:disturb_begin)         &
                        +disturbance_loss_from_litter(1:disturb_begin))
    else if (maxval(met(8,:)) > 0.99d0 .and. disturb_end /= nodays) then
       ! there has been a replacement level event, we will remove filter out a 2
       ! year period to allow for the most severe non-steady state response
       ! Croot
       in_out_root = sum(M_FLUXES(1:disturb_begin,6)) / sum(M_FLUXES(1:disturb_begin,12)+M_FLUXES(1:disturb_begin,23))
       in_out_root_disturb = sum(M_FLUXES(disturb_end:nodays,6)) / sum(M_FLUXES(disturb_end:nodays,12)+M_FLUXES(disturb_end:nodays,23))
!       in_out_root = in_out_root + (sum(M_FLUXES(disturb_end:nodays,6)) / sum(M_FLUXES(disturb_end:nodays,12)))
!       in_out_root = in_out_root * 0.5
       ! Cwood
       in_out_wood = sum(M_FLUXES(1:disturb_begin,7)) / sum(M_FLUXES(1:disturb_begin,11)+M_FLUXES(1:disturb_begin,24))
       in_out_wood_disturb = sum(M_FLUXES(disturb_end:nodays,7)) / sum(M_FLUXES(disturb_end:nodays,11)+M_FLUXES(disturb_end:nodays,24))
!       in_out_wood = in_out_wood + (sum(M_FLUXES(disturb_end:nodays,7)) / sum(M_FLUXES(disturb_end:nodays,11)))
!       in_out_wood = in_out_wood * 0.5
       ! Clitter
!       in_out_lit = sum(M_FLUXES(1:nodays,10)+M_FLUXES(1:nodays,12)+M_FLUXES(1:nodays,19) &
!                        +disturbance_residue_to_litter(1:nodays)) &
!                  / sum(M_FLUXES(1:nodays,13)+M_FLUXES(1:nodays,15)+disturbance_loss_from_litter(1:nodays))
       in_out_lit = (sum(M_FLUXES(1:disturb_begin,10)+M_FLUXES(1:disturb_begin,12)+ &
                         M_FLUXES(1:disturb_begin,19)+disturbance_residue_to_litter(1:disturb_begin))) &
                  / (sum(M_FLUXES(1:disturb_begin,13)+M_FLUXES(1:disturb_begin,15)+ &
                         disturbance_loss_from_litter(1:disturb_begin)))
       in_out_lit = in_out_lit + &
                   ((sum(M_FLUXES(disturb_end:nodays,10)+M_FLUXES(disturb_end:nodays,12)+ &
                         M_FLUXES(disturb_end:nodays,19)+disturbance_residue_to_litter(disturb_end:nodays))) &
                  / (sum(M_FLUXES(disturb_end:nodays,13)+M_FLUXES(disturb_end:nodays,15)+ &
                         disturbance_loss_from_litter(disturb_end:nodays))))
       in_out_lit = in_out_lit * 0.5d0
       ! Csom
!       in_out_som = sum(M_FLUXES(1:nodays,15)+disturbance_residue_to_som(1:nodays)) &
!                  / sum(M_FLUXES(1:nodays,14)+disturbance_loss_from_som(1:nodays))
       in_out_som = (sum(M_FLUXES(1:disturb_begin,15)+disturbance_residue_to_som(1:disturb_begin))) &
                  / (sum(M_FLUXES(1:disturb_begin,14)+disturbance_loss_from_som(1:disturb_begin)))
       in_out_som = in_out_som + &
                   ((sum(M_FLUXES(disturb_end:nodays,15)+disturbance_residue_to_som(disturb_end:nodays))) &
                  / (sum(M_FLUXES(disturb_end:nodays,14)+disturbance_loss_from_som(disturb_end:nodays))))
       in_out_som = in_out_som * 0.5d0
       ! Ccwd
!       in_out_cwd = sum(M_FLUXES(1:nodays,11)+disturbance_residue_to_cwd(1:nodays)) &
!                  / sum(M_FLUXES(1:nodays,19)+disturbance_loss_from_cwd(1:nodays))
       in_out_cwd = sum(M_FLUXES(1:disturb_begin,11)+disturbance_residue_to_cwd(1:disturb_begin)) &
                   / sum(M_FLUXES(1:disturb_begin,19)+disturbance_loss_from_cwd(1:disturb_begin))
       in_out_cwd = in_out_cwd + &
                    (sum(M_FLUXES(disturb_end:nodays,11)+disturbance_residue_to_cwd(disturb_end:nodays)) &
                  / sum(M_FLUXES(disturb_end:nodays,19)+disturbance_loss_from_cwd(disturb_end:nodays)))
       in_out_cwd = in_out_cwd * 0.5d0
       ! combined dead organic matter pools
!       in_out_dead = sum(M_FLUXES(1:nodays,10)+M_FLUXES(1:nodays,12) &
!                        +M_FLUXES(1:nodays,11)                       &
!                        +disturbance_residue_to_som(1:nodays)        &
!                        +disturbance_residue_to_cwd(1:nodays)        &
!                        +disturbance_residue_to_litter(1:nodays))    &
!                   / sum(M_FLUXES(1:nodays,13)+M_FLUXES(1:nodays,14) &
!                        +disturbance_loss_from_som(1:nodays)         &
!                        +disturbance_loss_from_cwd(1:nodays)         &
!                        +disturbance_loss_from_litter(1:nodays))
       in_out_dead = sum(M_FLUXES(1:disturb_begin,10) &
                        +M_FLUXES(1:disturb_begin,12) &
                        +M_FLUXES(1:disturb_begin,11) &
                        +disturbance_residue_to_som(1:disturb_begin) &
                        +disturbance_residue_to_cwd(1:disturb_begin) &
                        +disturbance_residue_to_litter(1:disturb_begin) ) &
                   / sum(M_FLUXES(1:disturb_begin,14) &
                        +M_FLUXES(1:disturb_begin,13) &
                        +disturbance_loss_from_som(1:disturb_begin) &
                        +disturbance_loss_from_cwd(1:disturb_begin) &
                        +disturbance_loss_from_litter(1:disturb_begin))
       in_out_dead = in_out_dead + &
                    (sum(M_FLUXES(disturb_end:nodays,10)+M_FLUXES(disturb_end:nodays,12) &
                        +M_FLUXES(disturb_end:nodays,11)+disturbance_residue_to_som(disturb_end:nodays) &
                        +disturbance_residue_to_cwd(disturb_end:nodays)        &
                        +disturbance_residue_to_litter(disturb_end:nodays) ) &
                   / sum(M_FLUXES(disturb_end:nodays,14)+M_FLUXES(disturb_end:nodays,13) &
                        +disturbance_loss_from_som(disturb_end:nodays) &
                        +disturbance_loss_from_cwd(disturb_end:nodays)         &
                        +disturbance_loss_from_litter(disturb_end:nodays)))
       in_out_dead = in_out_dead * 0.5d0
    else
       ! no replacement level disturbance so we assume everything must be in
       ! balance
       in_out_root = sum(M_FLUXES(1:nodays,6)) / sum(M_FLUXES(1:nodays,12)+M_FLUXES(1:nodays,23))
       in_out_root_disturb = 1d0
       in_out_wood = sum(M_FLUXES(1:nodays,7)) / sum(M_FLUXES(1:nodays,11)+M_FLUXES(1:nodays,24))
       in_out_wood_disturb = 1d0
       in_out_lit = sum(M_FLUXES(1:nodays,10)+M_FLUXES(1:nodays,12)+M_FLUXES(1:nodays,19) &
                        +disturbance_residue_to_litter(1:nodays)) &
                  / sum(M_FLUXES(1:nodays,13)+M_FLUXES(1:nodays,15)+disturbance_loss_from_litter(1:nodays))
       in_out_som = sum(M_FLUXES(1:nodays,15)+disturbance_residue_to_som(1:nodays)) &
                  / sum(M_FLUXES(1:nodays,14)+disturbance_loss_from_som(1:nodays))
       in_out_cwd = sum(M_FLUXES(1:nodays,11)+disturbance_residue_to_cwd(1:nodays)) &
                  / sum(M_FLUXES(1:nodays,19)+disturbance_loss_from_cwd(1:nodays))
       in_out_dead = sum(M_FLUXES(1:nodays,10)+M_FLUXES(1:nodays,12) &
                        +M_FLUXES(1:nodays,11)                       &
                        +disturbance_residue_to_som(1:nodays)        &
                        +disturbance_residue_to_cwd(1:nodays)        &
                        +disturbance_residue_to_litter(1:nodays))    &
                   / sum(M_FLUXES(1:nodays,13)+M_FLUXES(1:nodays,14) &
                        +disturbance_loss_from_som(1:nodays)         &
                        +disturbance_loss_from_cwd(1:nodays)         &
                        +disturbance_loss_from_litter(1:nodays))
    endif ! calculate in_out ratios, disturbed or not?

    ! derive mean pools
    do n = 1, nopools
       mean_pools(n)=cal_mean_pools(M_POOLS,n,nodays+1,nopools)
    end do

    !
    ! Begin EDCs here
    !

    ! tissue expantion has poorly described physiological limits, however we can
    ! safely constrain foliar(8), root(6) and wood(7) growth the < 20
    ! gCm-2.day-1
    if ((EDC2 == 1 .or. DIAG == 1) .and. (maxval(M_FLUXES(1:nodays,8)) > 20d0) ) then
       EDC2=0 ; EDCD%PASSFAIL(19)=0
    endif

    ! GPP allocation to foliage and labile cannot be 5 orders of magnitude
    ! difference from GPP allocation to roots
    if ((EDC2 == 1 .or. DIAG == 1) .and. (ffol > (5d0*froot) .or. (ffol*5d0) < froot)) then
       EDC2=0 ; EDCD%PASSFAIL(20)=0
    endif

    ! Part of the GSI test, we will assess EDC(3) here
    ! average turnover of foliage should not be less than wood
    ! neither should the average leaf life span be greater than 8 years
    if ((EDC2 == 1 .or. DIAG == 1) .and. (torfol < pars(6) .or. (1d0/(torfol*365.25d0)) > 8d0 .or. (1/(torfol*365.25)) < 0.15d0) ) then
         EDC2=0 ; EDCD%PASSFAIL(21)=0
    endif

    ! in contrast to the leaf longevity labile carbon stocks can be quite long
    ! lived, particularly in forests. (original value 8yr)
    ! Richardson et al (2015) New Phytologist, Clab residence time = 11 +/- 7.4 yrs (95CI = 18yr)
    if ((EDC2 == 1 .or. DIAG == 1) .and. (1d0/(torlab*365.25d0)) > 18d0 ) then
        EDC2=0 ; EDCD%PASSFAIL(22)=0
    endif

    ! finally we would not expect that the mean labile stock is greater than
    ! 12.5 % of the total ecosystem carbon stock (as we need structure to store
    ! labile).
    ! Gough et al (2009) Agricultural and Forest Meteorology. Avg 11, 12.5, 3 %
    ! (Max across species for branch, bole and coarse roots). Evidence that
    ! Branches accumulate labile C prior to bud burst from other areas.
    ! Wurth et al (2005) Oecologia, Clab 8 % of living biomass (DM) in tropical forest
    ! Richardson et al (2013), New Phytologist, Clab 2.24 +/- 0.44 % in temperate (max = 4.2 %)
    mean_ratio=M_POOLS(1:nodays,1)/M_POOLS(1:nodays,4) ; hak = 0
    where ( M_POOLS(1:nodays,4) == 0d0 )
           hak = 1 ; mean_ratio(1:nodays) = 0d0
    end where
    if ((EDC2 == 1 .or. DIAG == 1) .and. (sum(mean_ratio)/dble(nodays-sum(hak))) > 0.08d0 ) then
        EDC2=0 ; EDCD%PASSFAIL(23)=0
    endif

    !---------------------------------------------------
    ! First section will deal with Temperate plantation forest specific issues
    !---------------------------------------------------

    ! EDC 6
    ! ensure fine root : foliage ratio is between 0.1 and 0.45 (Albaugh et al
    ! 2004; Samuelson et al 2004; Vogel et al 2010; Akers et al 2013
    ! Duke ambient plots between 0.1 and 0.55
    ! Black et al 2009 Sitka Spruce chronosquence
    ! Q1 = 0.1278, median = 0.7488, mean = 1.0560 Q3 = 1.242
    ! lower CI = 0.04180938, upper CI = 4.06657167
    if ((EDC2 == 1 .or. DIAG == 1) .and. ( ((mean_pools(3)/mean_pools(2)) < 0.04d0) .or. &
                                           ((mean_pools(3)/mean_pools(2)) > 4.07d0) ) ) then
        EDC2=0 ; EDCD%PASSFAIL(24)=0
    end if

    ! generate mean annual pool values
    mean_annual_pools=0d0

    ! EDC 8
    ! assesses the exponential decay of specific pools

    ! loop vegetation foliar, roots, wood, litter, som and cwd
    if (EDC2 == 1 .or. DIAG == 1) then
        do n = 2, nopools ! previously started at 2 , nopools
           decay_coef=expdecay2(M_POOLS(exp_adjust:(nodays+1),:),n,deltat(exp_adjust:nodays),nopools,(nodays+1-exp_adjust+1))
           ! next assess the decay coefficient for meetings the EDC criterion
           if (abs(-log(2d0)/decay_coef) < (365.25d0*real(no_years_adjust)) .and. decay_coef < 0d0 ) then
              EDC2 = 0 ; EDCD%PASSFAIL(25)=0
           end if ! EDC conditions
        enddo
    endif

    if (EDC2 == 1 .or. DIAG == 1 .and. (maxval(met(8,:)) > 0.99d0 .and. disturb_end < (nodays-steps_per_year-1)) ) then
        do n = 4, 4 !, nopools
           decay_coef=expdecay2(M_POOLS(disturb_end:(nodays+1),:),n,deltat(disturb_end:nodays),nopools,(nodays+1-disturb_end+1))
           ! next assess the decay coefficient for meetings the EDC criterion
           if (abs(-log(2d0)/decay_coef) < (365.25d0*real(no_years_adjust)) .and. decay_coef < 0d0 ) then
              EDC2 = 0 ; EDCD%PASSFAIL(26)=0
           end if ! EDC conditions
        enddo
    endif

!    ! Cfoliar
!    if (EDC2 == 1 .or. DIAG == 1) then
!!        decay_coef=expdecay2(M_POOLS,2,deltat,nopools,nodays+1)
!        decay_coef=expdecay2(M_POOLS(exp_adjust:(nodays+1),:),2,deltat(exp_adjust:nodays),nopools,(nodays+1-exp_adjust+1))
!        ! next assess the decay coefficient for meetings the EDC criterion
!        if (abs(-log(2.0)/decay_coef) < (365.25*real(no_years_adjust)) .and. decay_coef < 0.0 ) then
!           EDC2 = 0 ; EDCD%PASSFAIL(8)=0
!        end if ! EDC conditions
!    end if ! EDC .or. DIAG condition

!    ! Croots
!    if (EDC2 == 1 .or. DIAG == 1) then
!!        decay_coef=expdecay2(M_POOLS,3,deltat,nopools,nodays+1)
!        decay_coef=expdecay2(M_POOLS(exp_adjust:(nodays+1),:),3,deltat(exp_adjust:nodays),nopools,(nodays+1-exp_adjust+1))
!        ! next assess the decay coefficient for meetings the EDC criterion
!        if (abs(-log(2.0)/decay_coef) < (365.25*real(no_years_adjust)) .and. decay_coef < 0.0 ) then
!           EDC2 = 0 ; EDCD%PASSFAIL(8)=0
!        end if ! EDC conditions
!    end if ! EDC .or. DIAG condition
!
!    ! Cwood
!    if (EDC2 == 1 .or. DIAG == 1) then
!!        decay_coef=expdecay2(M_POOLS,4,deltat,nopools,nodays+1)
!        decay_coef=expdecay2(M_POOLS(exp_adjust:(nodays+1),:),4,deltat(exp_adjust:nodays),nopools,(nodays+1-exp_adjust+1))
!        ! next assess the decay coefficient for meeting the EDC criterion
!        if (abs(-log(2.0)/decay_coef) < (365.25*real(no_years_adjust)) .and. decay_coef < 0.0 ) then
!           EDC2 = 0 ; EDCD%PASSFAIL(8)=0
!        endif ! EDC condition
!    endif ! EDC .or. DIAG condition
!
!    ! Clitter
!    if (EDC2 == 1 .or. DIAG == 1) then
!!        decay_coef=expdecay2(M_POOLS,5,deltat,nopools,nodays+1)
!        decay_coef=expdecay2(M_POOLS(exp_adjust:(nodays+1),:),5,deltat(exp_adjust:nodays),nopools,(nodays+1-exp_adjust+1))
!        ! next assess the decay coefficient for meetings the EDC criterion
!        if (abs(-log(2.0)/decay_coef) < (365.25*real(no_years_adjust)) .and. decay_coef < 0.0 ) then
!           EDC2 = 0 ; EDCD%PASSFAIL(8)=0
!        end if ! EDC conditions
!    end if ! EDC .or. DIAG condition
!
!    ! Csom
!    if (EDC2 == 1 .or. DIAG == 1) then
!!        decay_coef=expdecay2(M_POOLS,6,deltat,nopools,nodays+1)
!        decay_coef=expdecay2(M_POOLS(exp_adjust:(nodays+1),:),6,deltat(exp_adjust:nodays),nopools,(nodays+1-exp_adjust+1))
!        ! next assess the decay coefficient for meetings the EDC criterion
!        if (abs(-log(2.0)/decay_coef) < (365.25*real(no_years_adjust)) .and. decay_coef < 0.0 ) then
!           EDC2 = 0 ; EDCD%PASSFAIL(8)=0
!        end if ! EDC conditions
!    end if ! EDC .or. DIAG condition
!
!    ! Ccwd
!    if (EDC2 == 1 .or. DIAG == 1) then
!!        decay_coef=expdecay2(M_POOLS,7,deltat,nopools,nodays+1)
!        decay_coef=expdecay2(M_POOLS(exp_adjust:(nodays+1),:),7,deltat(exp_adjust:nodays),nopools,(nodays+1-exp_adjust+1))
!        ! next assess the decay coefficient for meetings the EDC criterion
!        if (abs(-log(2.0)/decay_coef) < (365.25*real(no_years_adjust)) .and. decay_coef < 0.0 ) then
!           EDC2 = 0 ; EDCD%PASSFAIL(8)=0
!        end if ! EDC conditions
!    end if ! EDC .or. DIAG condition

    ! EDC 9
    ! Mature forest maximum foliar biomass (gC.m-2) can be expected to be
    ! between 430 gC.m-2 and 768 gC.m-2, assume 50 % uncertainty (Loblolly Pine)
    ! Black et al Sitka Spruce estimates (gC.m-2)
    ! Lower CI = 379.2800 median = 477.1640 upper CI = 575.1956
    ! Harwood = 1200 ; Griffin = 960
    do y = 1, no_years
       ! derive mean annual foliar pool
       mean_annual_pools(y)=cal_mean_annual_pools(M_POOLS,y,2,nopools,deltat,nodays+1)
    end do ! year loop
    ! now check
    if ((EDC2 == 1 .or. DIAG == 1) .and. ( (sum(mean_annual_pools)/dble(no_years)) > 1200d0 ) ) then
       EDC2=0 ; EDCD%PASSFAIL(27)=0
    endif

    ! Steady state attractor - log ratio difference between inputs and outputs of
    ! system (maybe this should be pool specific?)
    EQF5=log(5d0) ! 10.0 = order magnitude; 2 = double and half (not accurate any longer)
    EQF2=log(2d0) ; EQF10=log(10d0) ; EQF20=log(20d0)

    if (EDC2 == 1 .or. DIAG == 1) then
       ! roots input / output ratio
       if (abs(log(in_out_root)) > EQF2) then
!   print*,"in_out_root",in_out_root
          EDC2 = 0 ; EDCD%PASSFAIL(28) = 0
       endif
       ! wood input / output ratio
       if (abs(log(in_out_wood)) > EQF5) then
   !print*,"in_out_wood",in_out_wood
          EDC2 = 0 ; EDCD%PASSFAIL(29) = 0
       endif
       ! litter input / output ratio
       if (abs(log(in_out_lit)) > EQF2) then
   !print*,"in_out_lit",in_out_lit
          EDC2 = 0 ; EDCD%PASSFAIL(30) = 0
       endif
       ! som input / output ratio
       if (abs(log(in_out_som)) > EQF2) then
   !print*,"in_out_som",in_out_som
          EDC2 = 0 ; EDCD%PASSFAIL(31) = 0
       endif
       ! cwd input / output ratio
       if (abs(log(in_out_cwd)) > EQF5) then
   !print*,"in_out_cwd",in_out_som
          EDC2 = 0 ; EDCD%PASSFAIL(32) = 0
       endif
       ! combined dead organic matter pools
       if (abs(log(in_out_dead)) > EQF2) then
          EDC2 = 0 ; EDCD%PASSFAIL(33) = 0
       endif
    endif ! EDC2 == 1 .or. DIAG == 1

    ! incase if disturbance
    if (maxval(met(8,:)) > 0.99d0 .and. disturb_end < (nodays-steps_per_year-1)) then
        ! roots input / output ratio
        if ((EDC2 == 1 .or. DIAG == 1) .and. abs(log(in_out_root_disturb)) > EQF5) then
!print*,"in_out_root_disturb",in_out_root_disturb
           EDC2 = 0 ; EDCD%PASSFAIL(34) = 0
        endif
        ! wood input / output ratio
        if ((EDC2 == 1 .or. DIAG == 1) .and. abs(log(in_out_wood_disturb)) > EQF20) then
!print*,"in_out_wood_disturb",in_out_wood_disturb
           EDC2 = 0 ; EDCD%PASSFAIL(35) = 0
        endif
    endif ! been cleared

    !
    ! EDC 14 - Fractional allocation to foliar biomass is well constrained
    ! across dominant ecosystem types (boreal -> temperate evergreen and deciduous),
    ! therefore this information can be used to
    ! contrain the foliar pool further. Through control of the
    ! photosynthetically active compoent of the carbon balance we can enforce
    ! additional contraint on the remainder of the system.
    ! Luyssaert et al (2007)

    if ((EDC2 == 1 .or. DIAG == 1) .and. (fNPP < 0.1d0 .or. fNPP > 0.5d0)) then
        EDC2 = 0 ; EDCD%PASSFAIL(36) = 0
    endif

    if ((EDC2 == 1 .or. DIAG == 1) .and. rNPP < 0.05d0) then
        EDC2 = 0 ; EDCD%PASSFAIL(37) = 0
    endif

    ! EDC 19 - Constrain the initial condition of wood stocks to that consistent
    ! with forestry age~yeild curves. UK forestry commission yield curves for
    ! evergreen species lowest yield and largest yield at year 60 is similar bound to those used above,
    ! so for generality these yield curves will be used here in broadest sense

    ! can only do this is we have age information
    if ((EDC2 == 1 .or. DIAG == 1) .and. DATAin%age > -1d0) then
        ! we will do this for the beginning of the simulation only.
        ! calculate sum pools (gC.m-2)
        model_living_C=M_POOLS(1,4) !M_POOLS(1,2)+M_POOLS(1,3)+M_POOLS(1,4)
        ! find out how many years into the simulation this is
        max_location=1
        ! call for empirical approximation of C accumulation curvies from
        ! forestry commissions
        call UK_forestry_commission_growth_curves(target_living_C,max_location)
        ! yield curve approximations result in unrealistic values early in
        ! the rotation so only assess if these values are sensible.
        ! This assumption means that the lower value may be negative which means its
        ! condition will always be passed but that the upper value must not be
        ! negative and of forest reasonable size
        if (target_living_C(2) > 100d0) then
            if (model_living_C < (target_living_C(1)) .or. model_living_C > (target_living_C(2))) then
                EDC2=0 ; EDCD%PASSFAIL(38)=0
             end if
        end if
    endif ! EDC2 .or. DIAG .and. age

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 50-100 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nopools .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_POOLS(nn,n) < 0d0 .or. M_POOLS(nn,n) /= M_POOLS(nn,n)) then
                 EDC2=0 ; PEDC=0 ; EDCD%PASSFAIL(50+n)=0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
       n = 1
       do while (n <= nofluxes .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_FLUXES(nn,n) < 0d0 .or. M_FLUXES(nn,n) /= M_FLUXES(nn,n)) then
                 EDC2=0 ; PEDC=0 ; EDCD%PASSFAIL(50+nopools+n)=0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

  end subroutine EDC2_GSI
  !
  !------------------------------------------------------------------
  !
  subroutine UK_forestry_commission_growth_curves(target_living_C,max_location)
    use cardamom_structures, only: DATAin

    ! subroutine uses PFT and yield classification to generate an estimate of
    ! expected living C accumulated at a given age. Equation generated Mg.ha-1
    ! we need to correct this to gC.m-2 for the model

    implicit none

    ! declare input / output variables
    double precision, intent(out) :: target_living_C(2) ! (gC.m-2)
    integer, intent(in) :: max_location(1) ! additional years from initial

    ! local variables
    double precision :: adjusted_age, tmp1(2),tmp2(2)
    integer :: i

    ! calculate adjusted age from initial conditions to max point
    adjusted_age=DATAin%age+max_location(1)

    ! set initial value for output
    target_living_C = 0d0

    ! loop through to get the minimum (1) and maximum estimates (2)
    ! which will be passed back to the model

    ! if we have an age (therefore it is a forest but we don't know even
    ! if it is evergreen or deciduos) we will assume the most generous
    ! range of values possible

    ! broadleaf
    tmp1(1) = 2.07956043460835d-05*adjusted_age**3d0 &
            + (-0.0141108480550955d0)*adjusted_age**2d0 &
            + 3.14928740556523d0*adjusted_age
    tmp1(2) = 0.000156065120683174d0*adjusted_age**3d0 &
            + (-0.0629544794948499d0)*adjusted_age**2d0 &
            + 8.30163202577001d0*adjusted_age
    ! evergreen
    tmp2(1) =  8.8519973125961d-06*adjusted_age**3d0 &
            + (-0.00822909089061558d0)*adjusted_age**2d0 &
            + 1.98952585135788d0*adjusted_age
    tmp2(2) = 0.00014916728414466d0*adjusted_age**3d0 &
            + (-0.0662815983372182d0)*adjusted_age**2d0 &
            + 9.55519207729034d0*adjusted_age
    ! work out which to use
    ! use smallest
    if (tmp1(1) < tmp2(1)) then
        target_living_C(1) = tmp1(1)*0.70d0
    else
        target_living_C(1) = tmp2(1)*0.70d0
    endif
    ! use biggest
    if (tmp1(2) > tmp2(2)) then
        target_living_C(2) = tmp1(2)*1.30d0
    else
        target_living_C(2) = tmp2(2)*1.30d0
    endif

    ! correct units from MgC.ha-1 to gC.m-2
    target_living_C=target_living_C*1d2

  end subroutine UK_forestry_commission_growth_curves
  !
  !------------------------------------------------------------------
  !
  double precision function cal_mean_pools(pools,pool_number,averaging_period,nopools)

    ! Function calculate the mean values of model pools / states across the
    ! entire simulation run

    implicit none

    ! declare input variables
    integer, intent(in) :: nopools          & !
                          ,pool_number      & !
                          ,averaging_period   !

    double precision,dimension(averaging_period,nopools), intent (in) :: pools

    ! declare local variables
    integer :: c

    ! initial conditions
    cal_mean_pools=0d0

    ! loop through now
    cal_mean_pools=sum(pools(1:averaging_period,pool_number))/dble(averaging_period)

    ! ensure return command issued
    return

  end function cal_mean_pools
  !
  !------------------------------------------------------------------
  !
  double precision function cal_mean_annual_pools(pools,year,pool_number,nopools,interval,averaging_period)

    ! Function calculates the mean model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: nopools     & ! how many pools in the model
                          ,year        & ! which year are we working on
                          ,averaging_period & ! number of days in analysis period
                          ,pool_number   ! which pool are we currently working on

    double precision, intent(in) :: pools(averaging_period,nopools) & ! input pool state variables
                         ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    integer :: startday, endday, c

    ! initialise the output variable
    cal_mean_annual_pools = 0d0

    ! calculate some constants
    startday=floor(365.25d0*dble(year-1)/(sum(interval)/dble(averaging_period-1)))+1
    endday=floor(365.25d0*dble(year)/(sum(interval)/dble(averaging_period-1)))

    ! pool through and work out the annual mean values
    cal_mean_annual_pools=sum(pools(startday:endday,pool_number))/dble(endday-startday)

    ! ensure function returns
    return

  end function cal_mean_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function cal_max_annual_pools(pools,year,pool_number,nopools,interval,averaging_period)

    ! Function calculates the max model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: nopools     & ! how many pools in the model
                          ,year        & ! which year are we working on
                          ,averaging_period & ! number of days in analysis period
                          ,pool_number   ! which pool are we currently working on

    double precision, intent(in) :: pools(averaging_period,nopools) & ! input pool state variables
                         ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    integer :: startday, endday, c

    ! initialise the output variable
    cal_max_annual_pools=0d0

    ! calculate some constants
    startday=floor(365.25d0*dble(year-1)/(sum(interval)/dble(averaging_period-1)))+1
    endday=floor(365.25d0*dble(year)/(sum(interval)/dble(averaging_period-1)))

    ! pool through and work out the annual max values
    cal_max_annual_pools=maxval(pools(startday:endday,pool_number))

    ! ensure function returns
    return

  end function cal_max_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function expdecay2(pools,pool_number,interval,nopools,averaging_period)

   ! Function to calculate the exponential decay coefficients used several EDCs.
   ! We assumpe the equation Cexp= a + b*exp(c*t)

   implicit none

   ! declare input variables
   integer, intent(in) :: nopools     & ! how many pools in the model
                         ,averaging_period & ! i.e. nodays + 1
                         ,pool_number   ! which pool are we currently working on

   double precision, intent(in) :: pools(averaging_period,nopools) & ! input pool state variables
                                  ,interval((averaging_period-1))      ! model time step in decimal days

   ! declare local variables
   integer :: n
   double precision :: P0    & ! initial pool value
            ,os,aw &
            ,MP0   & ! mean pool (year 1 to year end-2)
            ,MP1   & ! mean pool (year 2 to year end-1)
            ,MP0os & ! mean pool (year 1+os to year end-2+os)
            ,MP1os & ! mean pool (year 2+os to year end-2+os)
            ,dcdt1 & ! gradient of exponential over time in second year
            ,dcdt0   ! gradient of exponential over time in first year

   ! declare initial values / constants
   os = 1d0 ! offset in days
   aw = floor(365.25d0/(sum(interval)/dble(averaging_period-1))) ! averaging window
   MP0 = 0d0 ; MP1 = 0d0 ; MP0os = 0d0 ; MP1os = 0d0

   ! calculate mean pools within defined averaging window
   do n = 1, int(aw)
      MP0=MP0+pools(n,pool_number)
   end do ! for first year
   ! now average
   MP0=MP0/aw

   do n = int(aw)+1, int(aw*2)
      MP1=MP1+pools(n,pool_number)
   end do ! for second year
   ! now average
   MP1=MP1/aw

   do n = (1+int(os)), int(aw+os)
      MP0os=MP0os+pools(n,pool_number)
   end do ! for first year with offset
   ! now average
   MP0os=MP0os/aw

   do n = (int(aw+os)+1), int(aw*2d0+os)
      MP1os=MP1os+pools(n,pool_number)
   end do ! for second year withoffset
   ! now average
   MP1os=MP1os/aw

   ! derive mean gradient ratio (dcdt1/dcdt0)
   ! where dcdt1 is the numeric gradient between n+1 and n+365+1
   ! and dcdt0 os the numeric gradient between n and n+365
   dcdt1 = MP1os-MP0os
   dcdt0 = MP1-MP0

   ! using multiple year mean to determine c
   if ((dcdt1 > 0d0 .and. dcdt0 < 0d0) .or. (dcdt1 < 0d0 .and. dcdt0 > 0d0) &
       .or. dcdt1 == 0d0 .or. dcdt0 == 0d0) then
       ! then return error values
       expdecay2 = 1d0
   else
       expdecay2 = log(dcdt1/dcdt0) / (os*(sum(interval)/dble(averaging_period-1)))
   end if
   ! ensure return
   return

  end function expdecay2
  !
  !------------------------------------------------------------------
  !
  subroutine model_likelihood(PI,PARS,ML_out)
    use MCMCOPT, only:  PARAMETER_INFO
    use CARBON_MODEL_MOD, only: carbon_model
    use CARBON_MODEL_CROP_MOD, only: carbon_model_crop
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI ! parameter information

    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_out ! output variables for log-likelihood

    ! declare local variables
    double precision :: EDC,EDC1,EDC2

    ! initial values
    ML_out=0d0
    EDCD%DIAG=0

    if (DATAin%PFT == 1) then
       ! then we are crops so run these EDCs instead
       ! call EDCs which can be evaluated prior to running the model
       call EDC1_CROP(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
    else
       ! generic DALEC in use therefore relax...a little more
       ! call EDCs which can be evaluated prior to running the model
       call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
    endif ! crop choice

    ! now use the EDCD%EDC flag to determine if effect is kept
    if (DATAin%EDC == 1) then
        EDC = EDC1
    else
        EDC = 1
    end if

    ! update effect to the probabity
    ML_out=ML_out+log(EDC)

    ! if first set of EDCs have been passed
    if (EDC == 1) then
       ! calculate parameter log likelihood (assumed we have estimate of
       ! uncertainty)
       ML_out=ML_out+likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,PARS)

       if (DATAin%PFT == 1) then
          ! then this is a crop run....
           ! run the dalec model
           call CARBON_MODEL_CROP(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat,DATAin%nodays &
                                 ,DATAin%LAT, DATAin%M_LAI, DATAin%M_NEE &
                                 ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft &
                                 ,DATAin%nopars,DATAin%nomet,DATAin%nopools &
                                 ,DATAin%nofluxes,DATAin%M_GPP &
                                 ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                                 ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV &
                                 ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

           ! check edc2
           call EDC2_CROP(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                          ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                          ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                          ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)
       else
          ! yay generic run chill
          ! run the dalec model
          call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat,DATAin%nodays  &
                           ,DATAin%LAT, DATAin%M_LAI, DATAin%M_NEE       &
                           ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                           ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                           ,DATAin%M_GPP)

          ! check edc2
          call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                       ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                       ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                       ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)
       endif ! crop choice

       ! check if EDCs are switched on
       if (DATAin%EDC == 1) then
           EDC = EDC2
       else
           EDC = 1
       end if

       ! extra checks to ensure correct running of the model
       if (sum(DATAin%M_LAI) /= sum(DATAin%M_LAI) .or. sum(DATAin%M_GPP) /= sum(DATAin%M_GPP)) then
           EDC=0
       end if

       ! add EDC2 log-likelihood
       ML_out=ML_out+log(EDC)

       ! calculate final model likelihood when compared to obs
!       if (EDC == 1) ML_out=ML_out+likelihood(PI%npars,PARS)
       ML_out=ML_out+likelihood(PI%npars,PARS)

    end if ! EDC == 1

  end subroutine model_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood_p(npars,parpriors,parpriorunc,pars)
    ! function calculates the parameter based log-likelihood for the current set
    ! of parameters. This assumes that we have any actual priors / prior
    ! uncertainties to be working with. This does include initial states, as we
    ! consider them to be parameters

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars      & ! current parameter vector
                                                     ,parpriors & ! prior values for parameters
                                                     ,parpriorunc ! prior uncertainties

    ! declare local variables
    integer :: n

    ! set initial value
    likelihood_p = 0d0

    ! now loop through defined parameters for their uncertainties
    do n = 1, npars
       ! if there is actually a value
       if (parpriors(n) > -9999d0) then
           if (n == 17 .or. n == 11) then
               likelihood_p=likelihood_p-0.5d0*((pars(n)-parpriors(n))/parpriorunc(n))**2d0
           else if (n == 21) then
               likelihood_p=likelihood_p-0.5d0*((pars(n)-parpriors(n))/(parpriors(n)*parpriorunc(n)))**2d0
           else
               likelihood_p=likelihood_p-0.5d0*(log(pars(n)/parpriors(n))/log(parpriorunc(n)))**2d0
           end if
       end if
    end do

    ! dont for get to return
    return

  end function likelihood_p
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood(npars,pars)
    use cardamom_structures, only: DATAin

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y
    double precision :: tot_exp, pool_dynamics, tmp_var, infini
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    likelihood=0d0 ; infini=0d0

    ! GPP Log-likelihood
    tot_exp = 0d0
    if (DATAin%ngpp > 0) then
       do n = 1, DATAin%ngpp
         dn=DATAin%gpppts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((DATAin%M_GPP(dn)-DATAin%GPP(dn))/DATAin%GPP_unc(dn))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! LAI log-likelihood
    tot_exp = 0d0
    if (DATAin%nlai > 0) then
       ! loop split to allow vectorisation
       do n = 1, DATAin%nlai
         dn=DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (DATAin%M_LAI(dn) >= 0d0) then
             ! note that division is the uncertainty
             tot_exp=tot_exp+(log(max(0.001d0,DATAin%M_LAI(dn))/max(0.001d0,DATAin%LAI(dn)))/log(DATAin%LAI_unc(dn)))**2d0
         endif
       end do
       do n = 1, DATAin%nlai
         dn=DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (DATAin%M_LAI(dn) < 0d0) then
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp=tot_exp+(-log(infini))
         endif
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! NEE likelihood
    tot_exp = 0d0
    if (DATAin%nnee > 0) then
       do n = 1, DATAin%nnee
         dn=DATAin%neepts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((DATAin%M_NEE(dn)-DATAin%NEE(dn))/DATAin%NEE_unc(dn))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Reco likelihood
    tot_exp = 0d0
    if (DATAin%nreco > 0) then
       do n = 1, DATAin%nreco
         dn=DATAin%recopts(n)
         tmp_var=DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp=tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Cwood increment log-likelihood
    tot_exp = 0d0
    if (DATAin%nwoo > 0) then
       do n = 1, DATAin%nwoo
         dn=DATAin%woopts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+(log((DATAin%M_POOLS(dn,4)-DATAin%M_POOLS(dn-365,4)) &
                          / DATAin%WOO(dn))/log(DATAin%WOO_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Cfoliage log-likelihood
    tot_exp = 0d0
    if (DATAin%nCfol_stock > 0) then
       do n = 1, DATAin%nCfol_stock
         dn=DATAin%Cfol_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,2)/DATAin%Cfol_stock(dn))/log(2.))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Cfol_stock(dn)) &
                          / (DATAin%Cfol_stock(dn)*DATAin%Cfol_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Annual foliar maximum
    tot_exp = 0.
    if (DATAin%nCfolmax_stock > 0) then
       no_years=int(nint(sum(DATAin%deltat)/365.25d0))
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(no_years))
       ! determine the annual max for each pool
       do y = 1, no_years
          ! derive mean annual foliar pool
           mean_annual_pools(y)=cal_max_annual_pools(DATAin%M_POOLS,y,2,DATAin%nopools,DATAin%deltat,DATAin%nodays+1)
       end do ! year loop
       ! loop through the observations then
       do n = 1, DATAin%nCfolmax_stock
         ! load the observation position in stream
         dn=DATAin%Cfolmax_stockpts(n)
         ! determine which years this in in for the simulation
         y = ceiling( (dble(dn)*(sum(DATAin%deltat)/(DATAin%nodays))) / 365.25d0 )
         ! load the correct year into the analysis
         tmp_var = mean_annual_pools(y)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) &
                          / (DATAin%Cfolmax_stock(dn)*DATAin%Cfolmax_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    tot_exp = 0d0
    if (DATAin%nCwood_stock > 0) then
       do n = 1, DATAin%nCwood_stock
         dn=DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,4)/DATAin%Cwood_stock(dn))/log(2.))**2.
!         tot_exp=tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/(DATAin%Cwood_stock(dn)*0.20))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/(DATAin%Cwood_stock(dn)*DATAin%Cwood_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Cagb log-likelihood
    tot_exp = 0d0
    if (DATAin%nCagb_stock > 0) then
       do n = 1, DATAin%nCagb_stock
         dn=DATAin%Cagb_stockpts(n)
         ! remove coarse root fraction from wood (pars29)
         tmp_var = DATAin%M_POOLS(dn,4)-(DATAin%M_POOLS(dn,4)*pars(29))
         tot_exp=tot_exp+((tmp_var-DATAin%Cagb_stock(dn))/(DATAin%Cagb_stock(dn)*DATAin%Cagb_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Cstem log-likelihood
    tot_exp = 0d0
    if (DATAin%nCstem_stock > 0) then
       do n = 1, DATAin%nCstem_stock
         dn=DATAin%Cstem_stockpts(n)
         ! remove coarse root and branches from wood (pars29 and pars28)
         tmp_var = DATAin%M_POOLS(dn,4)-( (DATAin%M_POOLS(dn,4)*pars(29))+((DATAin%M_POOLS(dn,4)*pars(28))) )
         tot_exp=tot_exp+((tmp_var-DATAin%Cstem_stock(dn))/(DATAin%Cstem_stock(dn)*DATAin%Cstem_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Cbranch log-likelihood
    tot_exp = 0d0
    if (DATAin%nCbranch_stock > 0) then
       do n = 1, DATAin%nCbranch_stock
         dn=DATAin%Cbranch_stockpts(n)
         ! extract branch component from only
         tmp_var = DATAin%M_POOLS(dn,4)*pars(28)
         tot_exp=tot_exp+((tmp_var-DATAin%Cbranch_stock(dn))/(DATAin%Cbranch_stock(dn)*DATAin%Cbranch_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Ccoarseroot log-likelihood
    tot_exp = 0d0
    if (DATAin%nCcoarseroot_stock > 0) then
       do n = 1, DATAin%nCcoarseroot_stock
         dn=DATAin%Ccoarseroot_stockpts(n)
         ! extract coarse root component from wood only
         tmp_var = DATAin%M_POOLS(dn,4)*pars(29)
         tot_exp=tot_exp+((tmp_var-DATAin%Ccoarseroot_stock(dn))/(DATAin%Ccoarseroot_stock(dn)*DATAin%Ccoarseroot_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Croots log-likelihood
    tot_exp = 0d0
    if (DATAin%nCroots_stock > 0) then
       do n = 1, DATAin%nCroots_stock
         dn=DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,3)/DATAin%Croots_stock(dn))/log(2.))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Croots_stock(dn)) &
                         / (DATAin%Croots_stock(dn)*DATAin%Croots_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    tot_exp = 0d0
    if (DATAin%nClit_stock > 0) then
       do n = 1, DATAin%nClit_stock
         dn=DATAin%Clit_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+((log((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12)))*DATAin%M_POOLS(dn,5))/DATAin%Clit_stock(dn))/log(2.))**2d0
         tot_exp=tot_exp+(((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
                           *(DATAin%M_POOLS(dn,5))-DATAin%Clit_stock(dn))/(DATAin%Clit_stock(dn)*DATAin%Clit_stock_unc(dn)))**2d0
      end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

!    ! Clitter flux log-likelihood
!    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
!    ! originating from surface pools WARNING WARNING WARNING hack to leaf litter
!    ! flux not pool
!    tot_exp = 0.
!    if (DATAin%nClit_stock > 0) then
!       do n = 1, DATAin%nClit_stock
!         dn=DATAin%Clit_stockpts(n)
!         ! note that division is the uncertainty
!         tot_exp=tot_exp+((DATAin%M_FLUXES(dn,10)-DATAin%Clit_stock(dn))/(DATAin%Clit_stock(dn)*DATAin%Clit_stock_unc(dn)))**2
!      end do
!       likelihood=likelihood-0.5*tot_exp
!    endif

    ! Csom log-likelihood
    tot_exp = 0d0
    if (DATAin%nCsom_stock > 0) then
       do n = 1, DATAin%nCsom_stock
         dn=DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,6)/DATAin%Csom_stock(dn))/log(2.))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,6)-DATAin%Csom_stock(dn))/(DATAin%Csom_stock(dn)*DATAin%Csom_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! growth rate contraints ensures that end point pools are within uncertainty
    ! of 2 from the prior value for proportional rate of change
    ! this could / should be altered or extended to allow for time varying, pool
    ! specific values as above
    ! increments as done above for NEE,LAI and GPP
!    if (DATAin%otherpriors(1) > -9999) then
!       do n = 1, DATAin%nopools
!          pool_dynamics=DATAin%M_POOLS(DATAin%nodays+1,n)/DATAin%M_POOLS(1,n)
!          likelihood=likelihood-(0.5*(log(pool_dynamics/DATAin%otherpriors(1))/log(DATAin%otherpriors(1)))**2d0)
!       end do
!    end if

    ! check that log-likelihood is an actual number
    if (likelihood /= likelihood) then
       likelihood=log(infini)
    end if
    ! don't forget to return
    return

  end function likelihood
  !
  !------------------------------------------------------------------
  !
end module model_likelihood_module
