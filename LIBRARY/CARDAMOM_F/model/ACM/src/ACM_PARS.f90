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

!    PI%npars=12 ! dont forget to change in cardamom_io.f90

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    !
    ! Metabolic photosynthesis
    !

    ! Nitrogen use efficiency (gC/gN per m2 at optimum temperature)
    ! Derived from Vcmax reported in Wullschleger (1993), Journal of
    ! Experimental Botany, Vol 44, No. 262, pp. 907-920.
    PI%parmin(1)=03d0 
    PI%parmax(1)=40d0 

    ! max temperature for photosynthesis (oC)
    PI%parmin(2)=45d0
    PI%parmax(2)=70d0
    ! optimum temperature for photosynthesis (oC)
    PI%parmin(3)=20d0
    PI%parmax(3)=40d0
    ! kurtosis of photosynthesis temperature response
    PI%parmin(4)=0.10d0
    PI%parmax(4)=0.30d0

    ! light limited photosynthesis

    ! maximum canopy quantum yield (gC/MJ)
    PI%parmin(5)=1d0  !7.19298-(0.9*7.19298)
    PI%parmax(5)=7d0  !7.19298+(0.9*7.19298)

    !
    ! Canopy longwave radiation absorption
    !

    ! maximum absorbed radiation fraction
    PI%parmin(6)=0.90d0
    PI%parmax(6)=0.99d0

    ! LAI at which radiation absorption is half saturation (m2/m2)
    PI%parmin(7)=0.5d0
    PI%parmax(7)=2.5d0

    !
    ! Canopy NIR shortwave radiation absorption
    !

    ! maximum absorbed radiation (fraction)
    PI%parmin(8)=0.50d0
    PI%parmax(8)=0.90d0
    ! LAI at which radiation absorption is at half saturation (m2/m2)
    ! SPA's radiative transfer scheme absorbs 50 % of shortwave at lai == ~1.65
    PI%parmin(9)=0.5d0
    PI%parmax(9)=2.5d0

    !
    ! canopy conductance (gc) drivers
    !

    ! leafWP-soilWP (MPa); actual SPA parameter = 2+gplant
    PI%parmin(10)=-2.01d0 !-4.0
    PI%parmax(10)=-1.99d0 !-1.0

    !
    ! Canopy PAR shortwave radiation absorption
    !

    ! maximum absorbed
    PI%parmin(11)=0.60d0
    PI%parmax(11)=0.90d0
    ! LAI at which radiation absorption is at half saturation
    PI%parmin(12)=0.5d0
    PI%parmax(12)=2.5d0

    !
    ! Longwave reflectance back to sky
    !

    ! lai at which longwave returned to sky by canopy is half saturation
    PI%parmin(13)=0.5d0
    PI%parmax(13)=1.0d0

    !
    ! GPP / transpiration optimisation
    !

    ! iWUE (gC/m2leaf/day/mmolH2Ogs)
    PI%parmin(14)=1d-3
    PI%parmax(14)=0.1d0

    !
    ! Soil shortwave radiation absorption
    !

    ! soil sw radiation absorption (fraction)
    PI%parmin(15)=0.90d0!0.50d0 !0.10d0 !0.50
    PI%parmax(15)=0.99d0

    ! max sw radiation returned to sky by canopy (fraction)
    PI%parmin(16)=0.05d0
    PI%parmax(16)=0.50d0

    ! max fraction of longwave release from canopy (0.20:0.30)
    PI%parmin(17)=0.20d0
    PI%parmax(17)=0.30d0

    ! lai adjustment for long wave release from canopy (1:3)
    PI%parmin(18)=0.5d0
    PI%parmax(18)=2.5d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
