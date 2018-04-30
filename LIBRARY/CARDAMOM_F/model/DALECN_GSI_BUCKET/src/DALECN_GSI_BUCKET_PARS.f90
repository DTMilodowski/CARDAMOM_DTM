module MODEL_PARAMETERS

  implicit none

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains

  !
  !------------------------------------------------------------------
  !
  subroutine pars_info(PI)
    use MCMCOPT, only: parameter_info
    use cardamom_structures, only: DATAin

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or
    ! possibly should go into an alternate file which can be read in.
    ! This may
    ! improve the usability when it comes to reading these information
    ! in for
    ! different PFTs

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI

!    PI%npars=29 ! dont forget to change in cardamom_io.f90

    if (DATAin%PFT == 1) then
       ! crop model will be ran and therefore needs specific parameters to be
       ! called
       call crop_parameters(PI)
       call crop_development_parameters(PI%stock_seed_labile,PI%DS_shoot &
                                       ,PI%DS_root,PI%fol_frac,PI%stem_frac &
                                       ,PI%root_frac,PI%DS_LRLV,PI%LRLV &
                                       ,PI%DS_LRRT,PI%LRRT)

    else

       ! generic model

       !
       ! declare parameters
       !

       ! Decomposition of litter to som (fraction; temperature adjusted)
       PI%parmin(1)=0.0001368925d0 ! 20  years
       PI%parmax(1)=0.005475702d0  ! 0.5 years

       ! CN_Root (gC/gN; Kattge et al., 2011)
       PI%parmin(2)=13.85d0
       PI%parmax(2)=192.31d0

       ! GSI sensitivity for leaf growth
       PI%parmin(3)=1d0
       PI%parmax(3)=1.02d0

       ! Max labile turnover (fraction) to roots
       PI%parmin(4)=0.0001368925d0 ! 20 years
       PI%parmax(4)=0.05d0         ! 20 days

       ! GSI max leaf turnover (fraction)
       PI%parmin(5)=0.0003422313d0 !  8 years
       PI%parmax(5)=0.10d0         ! 10 days

       ! Turnover fraction of wood
       PI%parmin(6)=0.000009d0 ! 300  years
       PI%parmax(6)=0.001d0    ! 2.73 years

       ! Turnover fraction of roots
       PI%parmin(7)=0.0003422313d0 ! 8    years
       PI%parmax(7)=0.01d0         ! 0.27 years

       ! Turnover of litter to Rhet (fraction; temperature adjusted)
       PI%parmin(8)=0.0001368925d0 ! 27 years
       PI%parmax(8)=0.005475702d0  ! 0.50 years

       ! Turnover of som to Rhet (fraction; temperature adjusted)
       PI%parmin(9)=1.368925d-6 ! 2000 years
       PI%parmax(9)=0.0005d0      ! ~ 5 years

       ! Exponential coefficient for Rhet temperature response
       PI%parmin(10)=0.018d0
       PI%parmax(10)=0.06d0

       ! log10 avg foliar N (gN.m-2)
       ! Kattge et al., (2011)
       PI%parmin(11)=-0.2218487d0
       PI%parmax(11)= 0.5563025d0

       ! Max labile turnover fraction to foliage
       PI%parmin(12)=0.0001368925d0 ! 20 years
       PI%parmax(12)=0.05d0         ! 20 days

       ! Max labile turnover fraction to wood
       PI%parmin(13)=0.0001368925d0 ! 20 years
       PI%parmax(13)=0.05d0         ! 20 days

       ! GSI min temperature threshold (oC)
       PI%parmin(14)=225.15d0 !225d0
       PI%parmax(14)=305.15d0 !330d0

       ! GSI max temperature threshold (oC)
       PI%parmin(15)=263.15d0 !225d0
       PI%parmax(15)=330.15d0 !330d0

       ! GSI min photoperiod threshold (sec)
       PI%parmin(16)=3600d0*1d0  ! 1 hours
       PI%parmax(16)=3600d0*12d0 ! 12 hours

       ! Leaf Mass per unit Area (gC/m2)
       ! Kattge et al. 2011,
       PI%parmin(17)=10d0
       PI%parmax(17)=180d0 !200d0

       ! GSI max photoperiod threshold (sec)
       PI%parmin(24)=3600d0*6d0  ! 6 hours
       PI%parmax(24)=3600d0*18d0 ! 18 hours

       ! GSI min wSWP threshold (MPa)
       PI%parmin(25)=-2d0
       PI%parmax(25)=-0.01d0

       ! GSI max wSWP threshold (MPa)
       PI%parmin(26)=-2d0
       PI%parmax(26)=-0.01d0

       ! CN_wood (gC/gN; Kattge et al., 2011)
       PI%parmin(27)=169.5d0
       PI%parmax(27)=909.1d0

!       ! CN_wood baseline value (gC/gN) for logarithmic relation with woody biomass
!       PI%parmin(27)=100
!       PI%parmax(27)=300

       ! fraction of Cwood which is branch
       PI%parmin(28)=0.05d0
       PI%parmax(28)=0.40d0

       ! fraction of Cwood which is coarse root
       PI%parmin(29)=0.15d0
       PI%parmax(29)=0.30d0

       ! GSI senstivity for leaf senescence
       PI%parmin(34)=0.97d0 ! 0.96
       PI%parmax(34)=1.00d0

       ! GSI - have I just left a growing state (>1)
       PI%parmin(35)=0.50d0
       PI%parmax(35)=1.5d0

       ! GSI - initial GSI value
       PI%parmin(36)=1d0
       PI%parmax(36)=2d0

       ! Turnover fraction of CWD to litter (temperature adjusted)
       PI%parmin(38)=0.0001d0 ! 0.00001
       PI%parmax(38)=1d0/365.25d0 !0.005  ! 0.01

       ! BUCKET - root biomass (gbiomass/m2) needed to reach 50 % of max depth (m)
       PI%parmin(39)=50d0
       PI%parmax(39)=500d0

       ! BUCKET - maximum rooting depth (m)
       PI%parmin(40)=0.35d0
       PI%parmax(40)=10d0 !15.0 !20.0

       ! Reich - Leaf N linked respiration exponential coefficient
       PI%parmin(41)=0.935d0 ! 1.639-0.01
       PI%parmax(41)=1.774d0 ! 1.639+0.01

       ! Reich - Leaf N linked respiration intercept
       ! max/min values based on observed ranges from Reich et al (2008)
       ! Figure 1
       PI%parmin(42)=0.01d0 !0.01 !0.645
       PI%parmax(42)=1.65d0 !1.25 !0.911

       ! Reich - root N linked respiration exponential coefficient
       PI%parmin(43)=1.012d0 ! 1.352-0.01
       PI%parmax(43)=1.478d0 ! 1.352+0.01

       ! Reich - root N linked respiration intercept
       ! max/min values based on observed ranges from Reich et al (2008)
       ! Figure 1
       PI%parmin(44)=0.10d0 !0.01 !0.915
       PI%parmax(44)=1.90d0 !1.25 !1.079

       ! Reich - wood N linked respiration exponential coefficient
       PI%parmin(45)=1.170d0 ! 1.344-0.01
       PI%parmax(45)=1.478d0 ! 1.344+0.01

       ! Reich - wood N linked respiration intercept
       ! max/min values based on observed ranges from Reich et al (2008)
       ! Figure 1
       PI%parmin(46)=0.01d0 !0.01 !0.839
       PI%parmax(46)=1.90d0 !1.25 !1.053

       ! Initial leaf life span (days)
       PI%parmin(47)=50d0
       PI%parmax(47)=365.25d0*8d0

       ! baseline NUE (gC/gN/m2/day-1)
       PI%parmin(48)= 3d0 ! 11.0
       PI%parmax(48)=30d0 ! 12.0

!       ! CN_wood coefficient for increase due to C_wood (deltalCN per gC.m-2)
!       ! NOTE: values in log scale
!       PI%parmin(49)=0.01
!       PI%parmax(49)=0.25

       !
       ! INITIAL VALUES DECLARED HERE (gC/m2)
       !

       ! C labile
       PI%parmin(18)=1d0
       PI%parmax(18)=1000d0

       ! C foliar
       PI%parmin(19)=1d0
       PI%parmax(19)=1000d0

       ! C roots
       PI%parmin(20)=1d0
       PI%parmax(20)=1000d0

       ! C_wood
       PI%parmin(21)=1d0
       PI%parmax(21)=50000d0

       ! C litter
       PI%parmin(22)=1d0
       PI%parmax(22)=5000d0

       ! C_som
       PI%parmin(23)=100d0
       PI%parmax(23)=100000d0

       ! C CWD
       PI%parmin(37)=1d0
       PI%parmax(37)=10000d0

       !
       ! Replanting pools (gC/m2) values
       !

       ! C labile
       PI%parmin(30)=1d0
       PI%parmax(30)=500d0

       ! C foliar
       PI%parmin(31)=1d0
       PI%parmax(31)=500d0

       ! C roots
       PI%parmin(32)=1d0
       PI%parmax(32)=500d0

       ! C_wood
       PI%parmin(33)=1d0
       PI%parmax(33)=1000d0

    endif ! crop / default split

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
  subroutine crop_parameters(PI)

    ! Subroutine reads specific parameter ranges for the
    ! generic AT_DALEC model

    use MCMCOPT, only: parameter_info

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI

    !
    ! declare parameters
    !

!    PI%npars=34;

    ! Decomposition rate (frac/day)
    PI%parmin(1)=0.000001d0*24d0 ; PI%parmax(1)=0.0001d0*24d0

    ! Fraction of GPP to autotrophic pool
    PI%parmin(2)=0.2d0 ; PI%parmax(2)=0.7d0

    ! max development rate (day-1) DS (0->1)
    PI%parmin(3)=0.020d0 ; PI%parmax(3)=0.050d0
    ! max development rate (day-1) DS (1->2)
    PI%parmin(4)=0.010d0 ; PI%parmax(4)=0.050d0

    ! turnover rate foliage (frac/day)
    PI%parmin(5)=1.0d-4*24d0 ; PI%parmax(5)=0.02d0*24d0

    ! TOR wood* - 1% loss per year value (day-1)
    PI%parmin(6)=1d-4*24d0 ; PI%parmax(6)=0.01d0*24d0
    ! maximum rate of foliar turnover (hr-1) due to self-shading
    PI%parmin(7)=1d-5*24d0 ; PI%parmax(7)=0.01d0*24d0

    ! effective vernalisation days when plants are 50 % vernalised
    PI%parmin(8)=12d0 ; PI%parmax(8)=32d0

    ! mineralisation rate of litter (hr-1)
    PI%parmin(9)=1d-5*24d0 ; PI%parmax(9)=1d-2*24d0
    ! mineralisation rate of SOM (hr-1)
    PI%parmin(10)=1d-8*24d0 ; PI%parmax(10)=1d-3*24d0

    ! log10 avg foliar N (gN.m-2)
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=-0.50d0 ; PI%parmax(11)=0.7d0

    ! sow day
    PI%parmin(12)=115d0 ; PI%parmax(12)=350d0

    ! respiratory cost of labile transfer (per gC.m-2 labile)
    PI%parmin(13)=0.05d0 ; PI%parmax(13)=0.4d0

    ! phenological heat units required for emergence
    PI%parmin(14)=100d0 ; PI%parmax(14)=150d0

    ! harvest day
    PI%parmin(15)=15d0 ; PI%parmax(15)=350d0
    ! plough day
    PI%parmin(16)=365.25d0 ; PI%parmax(16)=365.25d0*4d0

    ! LMA
    PI%parmin(17)=10d0 ; PI%parmax(17)=100d0

    !
    ! NOTE number order not consistent
    !

    ! minimum temperature for development (oC)
    PI%parmin(26)=(-1d0+273.15d0) ; PI%parmax(26)=(10d0+273.15d0)  ! -10,10
    ! maximum temperature for development (oC)
    PI%parmin(27)=(10d0+273.15d0) ; PI%parmax(27)=(36d0+273.15d0)   ! 20,42
    ! optimum temperature for development (oC)
    PI%parmin(28)=(10d0+273.15d0) ; PI%parmax(28)=(30d0+273.15d0)   ! 10,35

    ! minimum temperature for vernalisation (oC)
    PI%parmin(29)=(-5.3d0+273.15d0) ; PI%parmax(29)=(-0.3d0+273.15d0)   ! -15,10
    ! maximum temperature for vernalisation (oC)
    PI%parmin(30)=(12.7d0+273.15d0) ; PI%parmax(30)=(18.7d0+273.15d0)    ! 5,30
    ! optimum temperature for vernalisation (oC)
    PI%parmin(31)=(2.9d0+273.15d0) ; PI%parmax(31)=(6.9d0+273.15d0)   ! -5,15

    ! critical photoperiod for development (hrs)
    PI%parmin(32)=6d0 ; PI%parmax(32)=13d0
    ! photoperiod sensitivity
    PI%parmin(33)=0.10d0 ; PI%parmax(33)=0.35d0

    ! turnover rate of labile
    PI%parmin(34)=1d-5*24d0  ; PI%parmax(34)=0.00625d0*24d0
    ! turnover rate of autotrophic pool
    PI%parmin(35)=0.001d0*24d0 ; PI%parmax(35)=0.07d0*24d0

    ! BUCKET - root biomass needed to reach 50 % of max depth
    PI%parmin(36)=10d0
    PI%parmax(36)=500d0

    ! BUCKET - maximum rooting depth
    PI%parmin(37)=0.35d0
    PI%parmax(37)=5d0 !20.0

    !
    ! INITIAL VALUES (gC.m-2) DECLARED HERE
    !

    ! C labile
    PI%parmin(18)=1d0 ; PI%parmax(18)=10d0
    ! C foliar
    PI%parmin(19)=1d0 ; PI%parmax(19)=5d0
    ! C roots
    PI%parmin(20)=1d0 ; PI%parmax(20)=5d0
    ! C_wood
    PI%parmin(21)=1d0 ; PI%parmax(21)=5d0
    ! C litter
    PI%parmin(22)=1d0 ; PI%parmax(22)=10d0
    ! C_som
    PI%parmin(23)=100d0 ; PI%parmax(23)=200000d0
    ! C autotrophic pool
    PI%parmin(24)=0.1d0 ; PI%parmax(24)=5d0
    ! C storage organ
    PI%parmin(25)=0.1d0 ; PI%parmax(25)=1d0

  end subroutine crop_parameters
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                        ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)

    ! subroutine reads in the fixed crop development files which are linked the
    ! the development state of the crops. The development model varies between
    ! which species. e.g. winter wheat and barley, spring wheat and barley

    implicit none

    ! declare inputs
    ! crop specific variables
    double precision,intent(inout) :: stock_seed_labile
    double precision, allocatable, dimension(:),intent(inout)  :: DS_shoot, & !
                                                                   DS_root, & !
                                                                  fol_frac, & !
                                                                 stem_frac, & !
                                                                 root_frac, & !
                                                                   DS_LRLV, & !
                                                                      LRLV, & !
                                                                   DS_LRRT, & !
                                                                      LRRT

    ! local variables..
    integer                 :: columns, i, rows, input_crops_unit, ios
    character(100) :: variables,filename

    ! for the moment hard code the file name
    filename="winter_wheat_development.csv"
    input_crops_unit = 20 ; ios = 0

    ! crop development file
    open(unit = input_crops_unit, file=trim(filename),iostat=ios, status='old', action='read')

    ! ensure we are definitely at the beginning
    rewind(input_crops_unit)

    ! read in the amount of carbon available (as labile) in each seed..
    read(unit=input_crops_unit,fmt=*)variables,stock_seed_labile,variables,variables

    ! read in C partitioning/fraction data and corresponding developmental
    ! stages (DS)
    ! shoot
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_shoot(rows) , fol_frac(rows) , stem_frac(rows)  )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_shoot(i), fol_frac(i), stem_frac(i)
    enddo

    ! root
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_root(rows) , root_frac(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_root(i), root_frac(i)
    enddo

    ! loss rates of leaves and roots
    ! leaves
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRLV(rows) , LRLV(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRLV(i), LRLV(i)
    enddo

    ! roots
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRRT(rows) , LRRT(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRRT(i), LRRT(i)
    enddo

    ! rewind and close
    rewind(input_crops_unit) ; close(input_crops_unit)

  end subroutine crop_development_parameters
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
