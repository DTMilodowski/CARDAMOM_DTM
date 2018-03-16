
module model_likelihood_module
  implicit none

  ! make all private
  private

  ! which to make open
  public :: model_likelihood, find_edc_initial_values

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: EDC
    integer :: DIAG
    integer :: PASSFAIL(100) ! allow space for 100 possible checks
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
    double precision :: tot_exp, ML, exp_orig, decay_coef, prob_exp, EDC, EDC1, EDC2, infini

    ! set initial values
    EDCD%DIAG=1

    ! call EDCs which can be evaluated prior to running the model
    call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat,DATAin%nodays  &
                   ,DATAin%LAT, DATAin%M_LAI, DATAin%M_NEE       &
                   ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                   ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                   ,DATAin%M_GPP)

    ! assess post running EDCs
    call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                  ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                  ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                  ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

    ! combine results
    if (DATAin%EDC == 1 .and. (EDC1 == 0 .or. EDC2 == 0 .or. &
        sum(DATAin%M_LAI) /= sum(DATAin%M_LAI) .or. sum(DATAin%M_GPP) /= sum(DATAin%M_GPP))) then
        EDC=0
    else
        EDC=1
    end if

    ! calculate the likelihood
    tot_exp=0d0
    do n = 1, EDCD%nedc
       tot_exp=tot_exp+(1d0-EDCD%PASSFAIL(n))
!       if (EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
    end do ! checking EDCs
    ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) stop

    ! convert to a probability
    prob_out=-0.5d0*(tot_exp*10d0)*DATAin%EDC

    ! override probability if parameter set gives NaN or near -infinitiy output
    call model_likelihood(PI,PARS,ML)

    infini=0d0
    if (DATAin%EDC == 0 .and. (ML /= ML .or. ML == log(infini) .or. ML == -log(infini) )) then
       prob_out=prob_out-0.5d0*10d0
    end if

    ! prob_out is the Log-Likelihood
    prob_out=prob_out

  end subroutine edc_model_likelihood
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
    integer :: n, DIAG, lai
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to la`bile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    ! set initial value
    EDC1=1
    DIAG=EDCD%DIAG


    ! set all EDCs to 1 (pass)
    EDCD%nedc=100
    EDCD%PASSFAIL(1:EDCD%nedc)=1

    !
    ! begin checking EDCs
    !

    ! the absorption of NIR should always be less than PAR at all LAI values
    do lai = 1,10
       if ((EDC1 == 1 .or. DIAG == 1) .and. &
           pars(8)*lai/(lai+pars(9)) >= pars(11)*lai/(lai+pars(12)) ) then
           EDC1 = 0 ; EDCD%PASSFAIL(1) = 0
       endif
    enddo


  end subroutine EDC1_GSI
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_GSI(npars,nomet,nofluxes,nopools,nodays,deltat &
                      ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                      ,meantemp,EDC2)

    use cardamom_structures, only: DATAin

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
    integer :: n, DIAG, no_years, y, PEDC, nn, num_EDC, max_location(1),i
    double precision :: mean_pools(nopools), G, decay_coef, meangpp, EQF &
                       ,model_living_C, target_living_C(2),hold
    double precision, dimension(:), allocatable :: mean_annual_pools,tmp
    double precision :: max_wood & !
                       ,torfol   & ! yearly average turnover
                       ,torlab   & !
                       ,fauto    & ! Fractions of GPP to autotrophic respiration
                       ,ffol     & ! Fraction of GPP to foliage
                       ,flab     & ! Fraction of GPP to labile pool
                       ,froot    & ! Fraction of GPP to root
                       ,fwood    & ! Fraction of GPP to wood
                       ,fsom     & ! fraction of GPP som under eqilibrium conditions
                       ,flit     & ! fraction of GPP to litter under equilibrium condition
                       ,delta_gsi

    ! set initial values
    DIAG=EDCD%DIAG
    EDC2=1

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 35 - 40 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nopools .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_POOLS(nn,n) < 0d0 .or. M_POOLS(nn,n) /= M_POOLS(nn,n)) then
                 EDC2=0 ; PEDC=0 ; EDCD%PASSFAIL(35+n)=0
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
  subroutine model_likelihood(PI,PARS,ML_out)
    use MCMCOPT, only:  PARAMETER_INFO
    use CARBON_MODEL_MOD, only: carbon_model
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

    ! call EDCs which can be evaluated prior to running the model
    call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
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
       if (EDC == 1) then
          ML_out=ML_out+likelihood(PI%npars,PARS)
       endif ! EDC still == 1

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
           if (n >= 1 .or. n <= 19) then
               ! Gaussian uncertainties for CO2 compensation and half saturation
               ! points
               likelihood_p=likelihood_p-0.5d0*((pars(n)-parpriors(n))/parpriorunc(n))**2
           else
               likelihood_p=likelihood_p-0.5d0*(log(pars(n)/parpriors(n))/log(parpriorunc(n)))**2
           endif
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
         tot_exp=tot_exp+((DATAin%M_GPP(dn)-DATAin%GPP(dn))/DATAin%GPP_unc(dn))**2
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Evapotranspiration (kg.m-2.day-1) Log-likelihood
    ! in this case transpiration only
    tot_exp = 0d0
    if (DATAin%nEvap > 0) then
       do n = 1, DATAin%nEvap
        dn=DATAin%Evappts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((DATAin%M_FLUXES(dn,2)-DATAin%Evap(dn))/DATAin%Evap_unc(dn))**2
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

    ! Borrowed wood increment to provide soil evaporation for ACM recal  (kg.m-2.day-1) Log-likelihood
    tot_exp = 0d0
    if (DATAin%nwoo > 0) then
       do n = 1, DATAin%nwoo
        dn=DATAin%woopts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((DATAin%M_FLUXES(dn,3)-DATAin%woo(dn))/DATAin%woo_unc(dn))**2
       end do
       likelihood=likelihood-0.5d0*tot_exp
    endif

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
