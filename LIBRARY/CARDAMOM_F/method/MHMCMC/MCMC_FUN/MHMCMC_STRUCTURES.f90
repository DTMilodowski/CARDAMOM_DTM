
module MCMCOPT 

  ! module contains the main derived types required by the MCMC itself. Input met
  ! drivers and observational data are however stored else where in a cardamom
  ! related module
  
  implicit none

  ! assume all private
  private

  ! explicit public statements
  public :: MCMC_OPTIONS, MCMC_OUTPUT ,COUNTERS, PARAMETER_INFO &
           ,PI,MCOUT,MCO, initialise_mcmc_output

  !
  ! module contains the declarable types used in the MCMC
  !

  ! contains MHMCMC options
  type MCMC_OPTIONS

    logical :: returnpars = .true. & ! return best fit parameters or not
              ,randparini = .true. & ! use random initial values parameters
              ,fixedpars  = .true.   ! use fixed initial values where inputs are not = -9999
    character(350) :: outfile   & ! output file name
                     ,stepfile    ! step file name
    double precision :: fADAPT    ! adapt step size for a given fraction of the full run
    integer :: nADAPT     & ! adapt step size after every N iterations
              ,nOUT       & ! number of requested output parameter sets
              ,nPRINT     & ! print info to screen every N solutions (0 to silent)
              ,nWRITE     & ! write to file every N solutions
              ,APPEND       ! append to existing output files (0 = delete existing file or 1 = append)

  end type ! MCMC_OPTIONS 
  ! create options type
  type(MCMC_OPTIONS), save :: MCO

  ! contains output information
  type MCMC_OUTPUT
    integer :: complete ! is MHMCMC completed (1 = yes, 0 = no)
    ! further metrics could be added here
    double precision, allocatable, dimension(:) :: best_pars ! store current best parameter set
  end type ! MCMC_OUTPUT
  ! create output type
  type(MCMC_OUTPUT), save :: MCOUT

  ! information which is needed determine progress
  type COUNTERS
    integer :: ACC    & ! total number of accepted solutions
              ,ACCLOC & ! number of recently accepted solutions
              ,ITER   & ! number of iterations attempted
              ,ACCLOC_ZEROS ! number of consecutive adaption periods to pass
                            ! without any acceptance

    double precision :: ACCRATE   ! local acceptance rate
  end type ! COUNTERS
  ! create counters type

  ! parameter structure defined here
  type PARAMETER_INFO
    double precision, allocatable, dimension(:) :: parmax   & ! maximum parameter values
                                        ,parmin   & ! minimum parameter values
                                        ,parini   & ! initial parameter values
                                        ,parfix   & ! do they need fixing (i.e. randomly generated)
                                        ,stepsize   ! parameter specific stepsize
    integer :: npars ! number of parameters to be solved
    ! crop specific variables
    double precision :: stock_seed_labile
    double precision, allocatable, dimension(:)  ::          DS_shoot, & !
                                                              DS_root, & ! 
                                                             fol_frac, & !
                                                            stem_frac, & !
                                                            root_frac, & !
                                                              DS_LRLV, & !
                                                                 LRLV, & ! 
                                                              DS_LRRT, & !
                                                                 LRRT

  end type ! PARAMETER_INFO
  ! create parameter info type
  type (PARAMETER_INFO), save :: PI

  contains

  !
  !------------------------------------------------------------------
  !
  subroutine initialise_mcmc_output(PI,MCOUT)

    ! subroutine allocated memory to the output arrays

    implicit none
    ! declare input variables
    type ( mcmc_output ), intent(inout) :: MCOUT
    type ( parameter_info ), intent(in) :: PI

    ! define dimensions for output structure
    allocate(MCOUT%best_pars(PI%npars))

    ! set to zero, will become 1 when MCMCM complete
    MCOUT%complete=0

  end subroutine initialise_mcmc_output
  !
  !------------------------------------------------------------------
  !
end module 
