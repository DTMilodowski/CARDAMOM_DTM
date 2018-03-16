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

    ! contains 6 fields with min max log for par and par

    ! 
    ! declare parameters
    ! 

    ! Fraction of litter turnover to atmosphere
    PI%parmin(1)=0.25
    PI%parmax(1)=0.75

    ! Fraction of GPP respired
    PI%parmin(2)=0.3
    PI%parmax(2)=0.7

    ! GSI sensitivity for leaf growth
    PI%parmin(3)=1.00
    PI%parmax(3)=1.025 !1.05

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4)=0.01
    PI%parmax(4)=1.0

    ! GSI max leaf turnover
    PI%parmin(5)=0.000001
    PI%parmax(5)=0.2

    ! TOR wood* - 1% loss per year value
    PI%parmin(6)=0.00001
    PI%parmax(6)=0.001

    ! TOR roots 
    PI%parmin(7)=0.0001
    PI%parmax(7)=0.01

    ! TOR litter
    PI%parmin(8)=0.0001
    PI%parmax(8)=0.01

    ! TOR SOM
    PI%parmin(9)=0.0000001
    PI%parmax(9)=0.001

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10)=0.018
    PI%parmax(10)=0.08

    ! log10 avg foliar N (gN.m-2) ! Canopy Efficiency
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=-0.50
    PI%parmax(11)=1.1

    ! GSI max labile turnover
    PI%parmin(12)=0.000001
    PI%parmax(12)=0.2

    ! Fraction to Clab*/
    PI%parmin(13)=0.01
    PI%parmax(13)=0.5

    ! GSI min temperature threshold (oC)
    PI%parmin(14)=225.0
    PI%parmax(14)=330.0

    ! GSI max temperature threshold (oC)
    PI%parmin(15)=225.0
    PI%parmax(15)=330.0

    ! GSI min photoperiod threshold (sec)
    PI%parmin(16)=3600.0 !21600d0 ! 6 hours
    PI%parmax(16)=3600.0*10.0 ! 64800d0 ! 18 hours

    ! LMA
    ! Kattge et al. 2011, 
    PI%parmin(17)=10.0 
    PI%parmax(17)=200.0 

    ! GSI max photoperiod threshold (sec)
    PI%parmin(24)=3600.0 !21600d0 ! 6 hours 
    PI%parmax(24)=64800.0 ! 18 hours

    ! GSI min VPD threshold (Pa)
    PI%parmin(25)=1.0
    PI%parmax(25)=5500.0

    ! GSI max VPD threshold (Pa)
    PI%parmin(26)=1.0
    PI%parmax(26)=5500.0

    ! critical GPP for LAI increase (fraction)
    PI%parmin(27)=1e-10
    PI%parmax(27)=0.20

    ! fraction of Cwood which is branch
    PI%parmin(28)=0.05
    PI%parmax(28)=0.65

    ! fraction of Cwood which is coarse root
    PI%parmin(29)=0.15
    PI%parmax(29)=0.45

    ! GSI senstivity for leaf senescence 
    PI%parmin(34)=0.96
    PI%parmax(34)=1.00

    ! GSI - have I just left a growing state (>1)
    PI%parmin(35)=0.50
    PI%parmax(35)=1.50

    ! GSI - initial GSI value 
    PI%parmin(36)=1.0
    PI%parmax(36)=2.0

    ! maximum rate of root labile turnover
    PI%parmin(37)=0.000001
    PI%parmax(37)=0.2

    ! maximum rate of wood labile turnover
    PI%parmin(38)=0.000001
    PI%parmax(38)=0.2

    ! CWD turnover 
    PI%parmin(43)=0.0001
    PI%parmax(43)=0.02

    ! total labile storage fraction in wood
    PI%parmin(45)=0.025
    PI%parmax(45)=0.125

    ! root C:N
    PI%parmin(46)=1.0
    PI%parmax(46)=100.0

    ! wood C:N
    PI%parmin(47)=200.0
    PI%parmax(47)=700.0

    ! soil C:N
    PI%parmin(48)=5.0
    PI%parmax(48)=50.0

    ! DON leaching fraction
    PI%parmin(49)=0.0015*0.5
    PI%parmax(49)=0.0015*1.5
   
    ! DIN leached per day
    PI%parmin(50)=0.00001*0.5
    PI%parmax(50)=0.00001*1.5
 
    ! DIN depostion per day
    PI%parmin(51)=0.5/365.25
    PI%parmax(51)=1.5/365.25

    ! Nitrogen use efficiency
    PI%parmin(52)=5.0**(-2.999299929993)
    PI%parmax(52)=15.0**(-2.999299929993)

    ! foliar N retranslocation fraction
    PI%parmin(53)=0.5
    PI%parmax(53)=0.9

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C foliar labile
    PI%parmin(18)=1.0
    PI%parmax(18)=1000.0

    ! C root labile
    PI%parmin(39)=1.0
    PI%parmax(39)=1000.0

    ! C wood labile
    PI%parmin(40)=1.0
    PI%parmax(40)=1000.0

    ! C foliar
    PI%parmin(19)=1.0
    PI%parmax(19)=1000.0

    ! C roots
    PI%parmin(20)=1.0
    PI%parmax(20)=1000.0
    
    ! C_wood
    PI%parmin(21)=1.0
    PI%parmax(21)=20000.0

    ! C litter
    PI%parmin(22)=1.0
    PI%parmax(22)=10000.0

    ! C cwd
    PI%parmin(44)=1.0
    PI%parmax(44)=10000.0
 
    ! C_som
    PI%parmin(23)=100.0
    PI%parmax(23)=200000.0

    ! N foliar labile 
    PI%parmin(54)=0.01
    PI%parmax(54)=500.0

    ! N_litter
    PI%parmin(55)=0.01
    PI%parmax(55)=1000.0

    ! Disolved inorganic N
    PI%parmin(56)=0.001
    PI%parmax(56)=20.0

    !
    ! Replanting pools values 
    !

    ! C foliar labile
    PI%parmin(30)=1.0
    PI%parmax(30)=1000.0

    ! C root labile
    PI%parmin(41)=1.0
    PI%parmax(41)=1000.0

    ! C wood labile
    PI%parmin(42)=1.0
    PI%parmax(42)=1000.0

    ! C foliar
    PI%parmin(31)=1.0
    PI%parmax(31)=1000.0

    ! C roots 
    PI%parmin(32)=1.0
    PI%parmax(32)=1000.0

    ! C_wood
    PI%parmin(33)=1.0
    PI%parmax(33)=1000.0

    ! N foliar labile
    PI%parmin(57)=0.01
    PI%parmax(57)=500.0
 
  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
