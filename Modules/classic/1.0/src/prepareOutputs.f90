!> \file
!> Central module that handles all preparation and writing to output files
module prepareOutputs

  ! J. Melton Mar 30 2015

  implicit none

  ! subroutines contained in this module:

  public  :: class_hh_w           ! Prepares and writes the CLASS (physics) half hourly file
  public  :: class_daily_aw       ! Accumulates and writes the CLASS (physics) daily outputs
  public  :: class_monthly_aw     ! Accumulates and writes the CLASS (physics) monthly outputs
  public  :: class_annual_aw      ! Accumulates and writes the CLASS (physics) annual outputs
  public  :: convertUnitsCTEM     ! Converts units prior to output for CTEM (biogeochemistry) variables.
  public  :: ctem_daily_aw        ! Accumulates and writes the CTEM (biogeochemistry) daily outputs
  public  :: ctem_monthly_aw      ! Accumulates and writes the CTEM (biogeochemistry) monthly outputs
  public  :: ctem_annual_aw       ! Accumulates and writes the CTEM (biogeochemistry) annual outputs


contains

  !> \ingroup prepareoutputs_class_halfhourly_aw
  !> @{
  !> Prepares and writes the CLASS (physics) half hourly file
  !! @author J. Melton
  subroutine class_hh_w (lonLocalIndex, latLocalIndex, nltest, nmtest, ncount, nday, iday, realyr)

    use classStateVars, only : class_rot, class_gat, initRowVars
    use ctemStateVars, only : c_switch, vrot
    use classicParams, only : ignd, icc, SBC, TFREZ
    use outputManager, only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: ncount
    integer, intent(in) :: nday
    integer, intent(in) :: iday
    integer, intent(in) :: realyr

    ! local variables
    real, dimension(1) :: timeStamp
    real :: ALTOT           !< Broadband albedo [-] (temporary variable)
    real :: EVAPSUM         !< Total evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$ (temporary variable)
    real :: FLSTAR
    real :: FSSTAR
    real :: TCN, TPN, TSN, TSURF, ZSN
    real :: an_grd, rml_grd, totvegarea
    real, dimension(:), allocatable :: anveggrd, rmlveggrd, fcanctot
    integer :: i, j, m

    ! pointers
    logical, pointer :: ctem_on              !< True if this run includes the biogeochemistry parameterizations (CTEM)
    logical, pointer :: dopertileoutput      !< Switch for making extra output files that are at the per tile level
    logical, pointer :: doperpftoutput       !< Switch for making extra output files that are at the per pft level
    real, pointer, dimension(:,:) :: TSFSGAT !< Ground surface temperature over subarea [K]
    real, pointer, dimension(:) :: FG        !< Subarea fractional coverage of modelled area - bare ground [ ]
    real, pointer, dimension(:) :: FC        !< Subarea fractional coverage of modelled area - ground under canopy [ ]
    real, pointer, dimension(:) :: FCS       !< Subarea fractional coverage of modelled area - snow-covered ground under canopy  [ ]
    real, pointer, dimension(:) :: FGS       !< Subarea fractional coverage of modelled area - snow-covered bare ground [ ]

    real, pointer, dimension(:,:,:) :: ailcgrow     !< Green LAI for CTEM's pfts
    real, pointer, dimension(:,:,:) :: anvegrow     !< Net photosynthesis rate for each pft
    real, pointer, dimension(:,:,:) :: rmlvegrow    !< Leaf maintenance respiration rate for each pft
    real, pointer, dimension(:,:,:) :: ancsvegrow   !< Net photosynthetic rate for CTEM's pfts for canopy over snow subarea
    real, pointer, dimension(:,:,:) :: ancgvegrow   !< Net photosynthetic rate for CTEM's pfts for canopy over ground subarea
    real, pointer, dimension(:,:,:) :: rmlcsvegrow  !< Leaf respiration rate for CTEM' pfts forcanopy over snow subarea
    real, pointer, dimension(:,:,:) :: rmlcgvegrow  !< Leaf respiration rate for CTEM' pfts forcanopy over ground subarea

    real, pointer, dimension(:,:) :: FAREROT !< Fractional coverage of mosaic tile on modelled area
    real, pointer, dimension(:,:) :: CDHROT  !< Surface drag coefficient for heat [ ]
    real, pointer, dimension(:,:) :: CDMROT  !< Surface drag coefficient for momentum [ ]
    real, pointer, dimension(:,:) :: HFSROT  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: TFXROT  !< Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: QEVPROT !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: QFSROT  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: QFXROT  !< Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: PETROT  !< Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: GAROT   !< Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: EFROT   !< Evaporation efficiency at ground surface [ ]
    real, pointer, dimension(:,:) :: GTROT   !< Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:,:) :: QGROT   !< Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:,:) :: ALIRROT !< Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer, dimension(:,:) :: ALVSROT !< Diagnosed total visible albedo of land surface [ ]
    real, pointer, dimension(:,:) :: SFCQROT !< Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:,:) :: SFCTROT !< Diagnosed screen-level air temperature [K]
    real, pointer, dimension(:,:) :: SFCUROT !< Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: SFCVROT !< Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: SFRHROT !<
    real, pointer, dimension(:,:) :: FSNOROT !< Diagnosed fractional snow coverage [ ]
    real, pointer, dimension(:,:) :: FLGGROT !< Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: FLGSROT !< Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: FLGVROT !< Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: FSGGROT !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: FSGSROT !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: FSGVROT !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HEVCROT !< Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HEVGROT !< Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HEVSROT !< Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HFSGROT !< Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HFSCROT !< Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HFSSROT !< Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HMFCROT !< Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HMFNROT !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HTCCROT !< Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HTCSROT !< Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: PCFCROT !< Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: PCLCROT !< Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: PCPGROT !< Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: PCPNROT !< Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: QFCFROT !< Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: QFCLROT !< Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: QFGROT  !< Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: QFNROT  !< Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROFROT  !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROFBROT !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROFCROT !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROFNROT !< Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROFOROT !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROFSROT !< Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROVGROT !< Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: WTRCROT !< Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: WTRGROT !< Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: WTRSROT !< Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: DRROT   !< Surface drag coefficient under neutral stability [ ]
    real, pointer, dimension(:,:) :: ILMOROT !< Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    real, pointer, dimension(:,:) :: UEROT   !< Friction velocity of air \f$[m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: HBLROT  !< Height of the atmospheric boundary layer [m]
    real, pointer, dimension(:,:,:) :: HMFGROT !< Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:,:) :: GFLXROT !< Heat conduction between soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:,:) :: HTCROT  !< Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:,:) :: QFCROT  !< Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:,:) :: TBARROT !< Temperature of soil layers [K]
    real, pointer, dimension(:,:,:) :: THICROT !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:,:) :: THLQROT !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: RHOSROT   !< Density of snow \f$[kg m^{-3}]\f$
    real, pointer, dimension(:,:) :: SCANROT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: RCANROT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: SNOROT  !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, pointer, dimension(:,:) :: WSNOROT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: TCANROT !< Vegetation canopy temperature [K]
    real, pointer, dimension(:,:) :: TSNOROT !< Snowpack temperature [K]
    real, pointer, dimension(:,:) :: TPNDROT !< Temperature of ponded water [K]
    real, pointer, dimension(:,:) :: ZPNDROT !< Depth of ponded water [m]

    real, dimension(:), pointer :: PREROW    !< Surface precipitation rate \f$[kg m^{-2}  s^{-1} ]\f$
    real, dimension(:), pointer :: UVROW     !< Wind speed at reference height \f$[m s^{-1} ]\f$
    real, dimension(:), pointer :: TAROW     !< Air temperature at reference height [K]
    real, pointer, dimension(:) :: QAROW     !< Specific humidity at reference height \f$[kg kg^{-1}]\f$
    real, pointer, dimension(:) :: PRESROW   !< Surface air pressure [Pa]

    ! FLAG some of the row variables here can probably be moved to local vars. JM Nov 2017.
    real, pointer, dimension(:) :: CDHROW  !< Surface drag coefficient for heat [ ]
    real, pointer, dimension(:) :: CDMROW  !< Surface drag coefficient for momentum [ ]
    real, pointer, dimension(:) :: HFSROW  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: TFXROW  !< Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    real, pointer, dimension(:) :: QEVPROW !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: QFSROW  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFXROW  !< Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: PETROW  !< Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: GAROW   !< Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: EFROW   !< Evaporation efficiency at ground surface [ ]
    real, pointer, dimension(:) :: GTROW   !< Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:) :: QGROW   !< Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: ALIRROW !< Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer, dimension(:) :: ALVSROW !< Diagnosed total visible albedo of land surface [ ]
    real, pointer, dimension(:) :: SFCQROW !< Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: SFCTROW !< Diagnosed screen-level air temperature [K]
    real, pointer, dimension(:) :: SFCUROW !< Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: SFCVROW !< Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: SFRHROW !< Diagnosed screen-level relative humidity [%]
    real, pointer, dimension(:) :: FSNOROW !< Diagnosed fractional snow coverage [ ]
    real, pointer, dimension(:) :: FLGGROW !< Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FLGSROW !< Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FLGVROW !< Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGGROW !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGSROW !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGVROW !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HEVCROW !< Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HEVGROW !< Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HEVSROW !< Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSGROW !< Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSCROW !< Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSSROW !< Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HMFCROW !< Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HMFNROW !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HTCCROW !< Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HTCSROW !< Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: PCFCROW !< Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCLCROW !< Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCPGROW !< Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCPNROW !< Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFCFROW !< Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFCLROW !< Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFGROW  !< Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFNROW  !< Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFROW  !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFBROW !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFCROW !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFNROW !< Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFOROW !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFSROW !< Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROVGROW !< Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: WTRCROW !< Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: WTRGROW !< Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: WTRSROW !< Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: DRROW   !< Surface drag coefficient under neutral stability [ ]
    real, pointer, dimension(:) :: ILMOROW !< Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    real, pointer, dimension(:) :: UEROW   !< Friction velocity of air \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: HBLROW  !< Height of the atmospheric boundary layer [m]
    real, pointer, dimension(:,:) :: HMFGROW !< Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: GFLXROW !< Heat conduction between soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HTCROW  !< Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: QFCROW  !< Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: TBARROW !< Temperature of soil layers [K]
    real, pointer, dimension(:,:) :: THALROW !< Total volumetric water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THICROW !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THLQROW !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:) :: FSSROW    !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FDLROW    !< Downwelling longwave sky radiation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: RHOSROW   !< Density of snow \f$[kg m^{-3}]\f$
    real, pointer, dimension(:) :: SNOROW    !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, pointer, dimension(:) :: TCANROW !< Vegetation canopy temperature [K]
    real, pointer, dimension(:) :: SCANROW !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: RCANROW !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: TSNOROW !< Snowpack temperature [K]
    real, pointer, dimension(:) :: TPNDROW !< Temperature of ponded water [K]
    real, pointer, dimension(:) :: ZPNDROW !< Depth of ponded water [m]
    real, pointer, dimension(:) :: WSNOROW !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:,:) :: fcancmxrow

    ctem_on           => c_switch%ctem_on
    dopertileoutput   => c_switch%dopertileoutput
    doperpftoutput    => c_switch%doperpftoutput
    TSFSGAT => class_gat%TSFSGAT
    FG => class_gat%FG
    FC => class_gat%FC
    FGS => class_gat%FGS
    FCS => class_gat%FCS
    ailcgrow => vrot%ailcg
    anvegrow => vrot%anveg
    rmlvegrow=> vrot%rmlveg
    ancsvegrow  => vrot%ancsveg
    ancgvegrow  => vrot%ancgveg
    rmlcsvegrow => vrot%rmlcsveg
    rmlcgvegrow => vrot%rmlcgveg

    FAREROT=> class_rot%FAREROT
    CDHROT => class_rot%CDHROT
    CDMROT => class_rot%CDMROT
    HFSROT => class_rot%HFSROT
    TFXROT => class_rot%TFXROT
    QEVPROT => class_rot%QEVPROT
    QFSROT => class_rot%QFSROT
    QFXROT => class_rot%QFXROT
    PETROT => class_rot%PETROT
    GAROT => class_rot%GAROT
    EFROT => class_rot%EFROT
    GTROT => class_rot%GTROT
    QGROT => class_rot%QGROT
    ALIRROT => class_rot%ALIRROT
    ALVSROT => class_rot%ALVSROT
    SFCQROT => class_rot%SFCQROT
    SFCTROT => class_rot%SFCTROT
    SFCUROT => class_rot%SFCUROT
    SFCVROT => class_rot%SFCVROT
    SFRHROT => class_rot%SFRHROT
    FSNOROT => class_rot%FSNOROT
    FLGGROT => class_rot%FLGGROT
    FLGSROT => class_rot%FLGSROT
    FLGVROT => class_rot%FLGVROT
    FSGGROT => class_rot%FSGGROT
    FSGSROT => class_rot%FSGSROT
    FSGVROT => class_rot%FSGVROT
    HEVCROT => class_rot%HEVCROT
    HEVGROT => class_rot%HEVGROT
    HEVSROT => class_rot%HEVSROT
    HFSCROT => class_rot%HFSCROT
    HFSGROT => class_rot%HFSGROT
    HFSSROT => class_rot%HFSSROT
    HMFCROT => class_rot%HMFCROT
    HMFNROT => class_rot%HMFNROT
    HTCCROT => class_rot%HTCCROT
    HTCSROT => class_rot%HTCSROT
    PCFCROT => class_rot%PCFCROT
    PCLCROT => class_rot%PCLCROT
    PCPGROT => class_rot%PCPGROT
    PCPNROT => class_rot%PCPNROT
    QFCFROT => class_rot%QFCFROT
    QFCLROT => class_rot%QFCLROT
    QFGROT => class_rot%QFGROT
    QFNROT => class_rot%QFNROT
    ROFROT => class_rot%ROFROT
    ROFBROT => class_rot%ROFBROT
    ROFCROT => class_rot%ROFCROT
    ROFNROT => class_rot%ROFNROT
    ROFOROT => class_rot%ROFOROT
    ROFSROT => class_rot%ROFSROT
    ROVGROT => class_rot%ROVGROT
    WTRCROT => class_rot%WTRCROT
    WTRGROT => class_rot%WTRGROT
    WTRSROT => class_rot%WTRSROT
    DRROT => class_rot%DRROT
    ILMOROT => class_rot%ILMOROT
    UEROT => class_rot%UEROT
    HBLROT => class_rot%HBLROT
    HMFGROT => class_rot%HMFGROT
    GFLXROT => class_rot%GFLXROT
    HTCROT => class_rot%HTCROT
    QFCROT => class_rot%QFCROT
    TBARROT=> class_rot%TBARROT
    THICROT=> class_rot%THICROT
    THLQROT=> class_rot%THLQROT
    RHOSROT => class_rot%RHOSROT
    SCANROT => class_rot%SCANROT
    RCANROT => class_rot%RCANROT
    SNOROT => class_rot%SNOROT
    WSNOROT => class_rot%WSNOROT
    TCANROT=> class_rot%TCANROT
    TSNOROT=> class_rot%TSNOROT
    TPNDROT=> class_rot%TPNDROT
    ZPNDROT=> class_rot%ZPNDROT
    PREROW  => class_rot%PREROW
    UVROW => class_rot%UVROW
    TAROW  => class_rot%TAROW
    QAROW => class_rot%QAROW
    PRESROW => class_rot%PRESROW
    FSSROW => class_rot%FSSROW
    FDLROW => class_rot%FDLROW
    RHOSROW => class_rot%RHOSROW
    SNOROW => class_rot%SNOROW
    CDHROW => class_rot%CDHROW
    CDMROW => class_rot%CDMROW
    HFSROW => class_rot%HFSROW
    TFXROW => class_rot%TFXROW
    QEVPROW => class_rot%QEVPROW
    QFSROW => class_rot%QFSROW
    QFXROW => class_rot%QFXROW
    PETROW => class_rot%PETROW
    GAROW => class_rot%GAROW
    EFROW => class_rot%EFROW
    GTROW => class_rot%GTROW
    QGROW => class_rot%QGROW
    ALIRROW => class_rot%ALIRROW
    ALVSROW => class_rot%ALVSROW
    SFCQROW => class_rot%SFCQROW
    SFCTROW => class_rot%SFCTROW
    SFCUROW => class_rot%SFCUROW
    SFCVROW => class_rot%SFCVROW
    SFRHROW => class_rot%SFRHROW
    FSNOROW => class_rot%FSNOROW
    FLGGROW => class_rot%FLGGROW
    FLGSROW => class_rot%FLGSROW
    FLGVROW => class_rot%FLGVROW
    FSGGROW => class_rot%FSGGROW
    FSGSROW => class_rot%FSGSROW
    FSGVROW => class_rot%FSGVROW
    HEVCROW => class_rot%HEVCROW
    HEVGROW => class_rot%HEVGROW
    HEVSROW => class_rot%HEVSROW
    HFSCROW => class_rot%HFSCROW
    HFSGROW => class_rot%HFSGROW
    HFSSROW => class_rot%HFSSROW
    HMFCROW => class_rot%HMFCROW
    HMFNROW => class_rot%HMFNROW
    HTCCROW => class_rot%HTCCROW
    HTCSROW => class_rot%HTCSROW
    PCFCROW => class_rot%PCFCROW
    PCLCROW => class_rot%PCLCROW
    PCPGROW => class_rot%PCPGROW
    PCPNROW => class_rot%PCPNROW
    QFCFROW => class_rot%QFCFROW
    QFCLROW => class_rot%QFCLROW
    QFGROW => class_rot%QFGROW
    QFNROW => class_rot%QFNROW
    ROFROW => class_rot%ROFROW
    ROFBROW => class_rot%ROFBROW
    ROFCROW => class_rot%ROFCROW
    ROFNROW => class_rot%ROFNROW
    ROFOROW => class_rot%ROFOROW
    ROFSROW => class_rot%ROFSROW
    ROVGROW => class_rot%ROVGROW
    WTRCROW => class_rot%WTRCROW
    WTRGROW => class_rot%WTRGROW
    WTRSROW => class_rot%WTRSROW
    DRROW => class_rot%DRROW
    ILMOROW => class_rot%ILMOROW
    UEROW => class_rot%UEROW
    HBLROW => class_rot%HBLROW
    HMFGROW => class_rot%HMFGROW
    GFLXROW => class_rot%GFLXROW
    HTCROW => class_rot%HTCROW
    QFCROW => class_rot%QFCROW
    TBARROW => class_rot%TBARROW
    THALROW => class_rot%THALROW
    THLQROW => class_rot%THLQROW
    THICROW => class_rot%THICROW
    TCANROW => class_rot%TCANROW
    SCANROW => class_rot%SCANROW
    RCANROW => class_rot%RCANROW
    TSNOROW => class_rot%TSNOROW
    WSNOROW => class_rot%WSNOROW
    TPNDROW => class_rot%TPNDROW
    ZPNDROW => class_rot%ZPNDROW
    fcancmxrow => vrot%fcancmx

    ! Calculate grid cell average diagnostic fields.

    ! First set all to zero
    call initRowVars(nltest)

    do I = 1,NLTEST
      do M = 1,NMTEST
        CDHROW(I) = CDHROW(I) + CDHROT(I,M) * FAREROT(I,M)
        CDMROW(I) = CDMROW(I) + CDMROT(I,M) * FAREROT(I,M)
        HFSROW(I) = HFSROW(I) + HFSROT(I,M) * FAREROT(I,M)
        TFXROW(I) = TFXROW(I) + TFXROT(I,M) * FAREROT(I,M)
        QEVPROW(I) = QEVPROW(I) + QEVPROT(I,M) * FAREROT(I,M)
        QFSROW(I) = QFSROW(I) + QFSROT(I,M) * FAREROT(I,M)
        QFXROW(I) = QFXROW(I) + QFXROT(I,M) * FAREROT(I,M)
        PETROW(I) = PETROW(I) + PETROT(I,M) * FAREROT(I,M)
        RHOSROW(I) = RHOSROW(I) + RHOSROT(I,M) * FAREROT(I,M)
        GAROW(I) = GAROW(I) + GAROT(I,M) * FAREROT(I,M)
        EFROW(I) = EFROW(I) + EFROT(I,M) * FAREROT(I,M)
        GTROW(I) = GTROW(I) + GTROT(I,M) * FAREROT(I,M)
        QGROW(I) = QGROW(I) + QGROT(I,M) * FAREROT(I,M)
        ALVSROW(I) = ALVSROW(I) + ALVSROT(I,M) * FAREROT(I,M)
        ALIRROW(I) = ALIRROW(I) + ALIRROT(I,M) * FAREROT(I,M)
        SFCTROW(I) = SFCTROW(I) + SFCTROT(I,M) * FAREROT(I,M)
        SFCUROW(I) = SFCUROW(I) + SFCUROT(I,M) * FAREROT(I,M)
        SFCVROW(I) = SFCVROW(I) + SFCVROT(I,M) * FAREROT(I,M)
        SFCQROW(I) = SFCQROW(I) + SFCQROT(I,M) * FAREROT(I,M)
        SFRHROW(I) = SFRHROW(I) + SFRHROT(I,M) * FAREROT(I,M)
        FSNOROW(I) = FSNOROW(I) + FSNOROT(I,M) * FAREROT(I,M)
        FSGVROW(I) = FSGVROW(I) + FSGVROT(I,M) * FAREROT(I,M)
        FSGSROW(I) = FSGSROW(I) + FSGSROT(I,M) * FAREROT(I,M)
        FSGGROW(I) = FSGGROW(I) + FSGGROT(I,M) * FAREROT(I,M)
        FLGVROW(I) = FLGVROW(I) + FLGVROT(I,M) * FAREROT(I,M)
        FLGSROW(I) = FLGSROW(I) + FLGSROT(I,M) * FAREROT(I,M)
        FLGGROW(I) = FLGGROW(I) + FLGGROT(I,M) * FAREROT(I,M)
        HFSCROW(I) = HFSCROW(I) + HFSCROT(I,M) * FAREROT(I,M)
        HFSSROW(I) = HFSSROW(I) + HFSSROT(I,M) * FAREROT(I,M)
        HFSGROW(I) = HFSGROW(I) + HFSGROT(I,M) * FAREROT(I,M)
        HEVCROW(I) = HEVCROW(I) + HEVCROT(I,M) * FAREROT(I,M)
        HEVSROW(I) = HEVSROW(I) + HEVSROT(I,M) * FAREROT(I,M)
        HEVGROW(I) = HEVGROW(I) + HEVGROT(I,M) * FAREROT(I,M)
        HMFCROW(I) = HMFCROW(I) + HMFCROT(I,M) * FAREROT(I,M)
        HMFNROW(I) = HMFNROW(I) + HMFNROT(I,M) * FAREROT(I,M)
        HTCCROW(I) = HTCCROW(I) + HTCCROT(I,M) * FAREROT(I,M)
        HTCSROW(I) = HTCSROW(I) + HTCSROT(I,M) * FAREROT(I,M)
        PCFCROW(I) = PCFCROW(I) + PCFCROT(I,M) * FAREROT(I,M)
        PCLCROW(I) = PCLCROW(I) + PCLCROT(I,M) * FAREROT(I,M)
        PCPNROW(I) = PCPNROW(I) + PCPNROT(I,M) * FAREROT(I,M)
        PCPGROW(I) = PCPGROW(I) + PCPGROT(I,M) * FAREROT(I,M)
        QFGROW(I) = QFGROW(I) + QFGROT(I,M) * FAREROT(I,M)
        QFNROW(I) = QFNROW(I) + QFNROT(I,M) * FAREROT(I,M)
        QFCLROW(I) = QFCLROW(I) + QFCLROT(I,M) * FAREROT(I,M)
        QFCFROW(I) = QFCFROW(I) + QFCFROT(I,M) * FAREROT(I,M)
        ROFROW(I) = ROFROW(I) + ROFROT(I,M) * FAREROT(I,M)
        ROFOROW(I) = ROFOROW(I) + ROFOROT(I,M) * FAREROT(I,M)
        ROFSROW(I) = ROFSROW(I) + ROFSROT(I,M) * FAREROT(I,M)
        ROFBROW(I) = ROFBROW(I) + ROFBROT(I,M) * FAREROT(I,M)
        ROFCROW(I) = ROFCROW(I) + ROFCROT(I,M) * FAREROT(I,M)
        ROFNROW(I) = ROFNROW(I) + ROFNROT(I,M) * FAREROT(I,M)
        ROVGROW(I) = ROVGROW(I) + ROVGROT(I,M) * FAREROT(I,M)
        WTRCROW(I) = WTRCROW(I) + WTRCROT(I,M) * FAREROT(I,M)
        WTRSROW(I) = WTRSROW(I) + WTRSROT(I,M) * FAREROT(I,M)
        WTRGROW(I) = WTRGROW(I) + WTRGROT(I,M) * FAREROT(I,M)
        DRROW(I) = DRROW(I) + DRROT(I,M) * FAREROT(I,M)
        ! wtableROW(I)=wtableROW(I)+wtableROT(I,M)*FAREROT(I,M)  ! FLAG
        ILMOROW(I) = ILMOROW(I) + ILMOROT(I,M) * FAREROT(I,M)
        UEROW(I) = UEROW(I) + UEROT(I,M) * FAREROT(I,M)
        HBLROW(I) = HBLROW(I) + HBLROT(I,M) * FAREROT(I,M)
        do J = 1,IGND
          HMFGROW(I,J) = HMFGROW(I,J) + HMFGROT(I,M,J) * FAREROT(I,M)
          HTCROW(I,J) = HTCROW(I,J) + HTCROT(I,M,J) * FAREROT(I,M)
          QFCROW(I,J) = QFCROW(I,J) + QFCROT(I,M,J) * FAREROT(I,M)
          GFLXROW(I,J) = GFLXROW(I,J) + GFLXROT(I,M,J) * FAREROT(I,M)
          TBARROW(i,j) = TBARROW(i,j) + TBARROT(i,m,j) * FAREROT(i,m)
          THLQROW(i,j) = THLQROW(i,j) + THLQROT(i,m,j) * FAREROT(i,m)
          THICROW(i,j) = THICROW(i,j) + THICROT(i,m,j) * FAREROT(i,m)
        end do ! loop 550
      end do ! loop 575
    end do ! loop 600


    ! Prepare the timestamp for this timestep.
    timeStamp = consecDays + real(iday - 1) + ((real(ncount) - 1.) / real(nday))

    ! Now prepare and write out the grid averaged physics variables to output files
    do I = 1,NLTEST

      ! First write out the model meteorological forcing so they can be compared to inputs
      ! or to check on those generated by the dissagregation module.
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fss_hh' ,timeStamp,'rsds', [FSSROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'fdl_hh' ,timeStamp,'rlds', [FDLROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'pre_hh' ,timeStamp,'pr',   [PREROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'uv_hh' ,timeStamp,'uvas',  [UVROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ta_hh' ,timeStamp,'tas',   [TAROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qa_hh' ,timeStamp,'huss',  [QAROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'pres_hh' ,timeStamp,'ps',  [PRESROW(I)])

      ALTOT = 0.0
      if (FSSROW(I) > 0.0) ALTOT = (FSSROW(I) - (FSGVROW(I) + FSGSROW(I) + FSGGROW(I)))/FSSROW(I)
      FSSTAR = FSSROW(I) * (1.0 - ALTOT)

      call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_hh' ,timeStamp,'rss', [FSSTAR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'altot_hh' ,timeStamp,'albs', [ALTOT])

      FLSTAR = FDLROW(I) - SBC * GTROW(I) ** 4
      call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_hh' ,timeStamp,'rls', [FLSTAR])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_hh'     ,timeStamp,'hfss', [HFSROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_hh'     ,timeStamp,'hfls', [QEVPROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'snm_hh' ,timeStamp,'snm', [HMFNROW(I)])

      ZSN = 0.0
      if (RHOSROW(I) > 0.0) ZSN = SNOROW(I)/RHOSROW(I)
      call writeOutput1D(lonLocalIndex,latLocalIndex,'sisnthick_hh'   ,timeStamp,'sisnthick', [ZSN])

      TCN = 0.0
      if (TCANROW(I) > 0.01) TCN = TCANROW(I) - TFREZ
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tcs_hh'   ,timeStamp,'tcs', [TCN])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'scan_hh'   ,timeStamp,'scanopy', [SCANROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rcan_hh'   ,timeStamp,'rcanopy', [RCANROW(I)])

      ! Find the surface temperature for the grid average over all subareas (1-4)
      TSURF = FCS(I) * TSFSGAT(I,1) + FGS(I) * TSFSGAT(I,2) + FC(I) * TSFSGAT(I,3) + FG(I) * TSFSGAT(I,4)
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsurf_hh'  ,timeStamp,'ts', [TSURF])

      TSN = 0.0
      if (TSNOROW(I) > 0.01) TSN = TSNOROW(I) - TFREZ
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tsno_hh'   ,timeStamp,'tsn', [TSN])

      TPN = 0.0
      if (TPNDROW(I) > 0.01) TPN = TPNDROW(I) - TFREZ
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tpond_hh'   ,timeStamp,'tpond', [TPN])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'zpond_hh'   ,timeStamp,'zpond', [ZPNDROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gt_hh'   ,timeStamp,'tsblack', [GTROW(I) - TFREZ])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'sno_hh' ,timeStamp,'snw', [SNOROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'wsno_hh',timeStamp,'wsnw', [WSNOROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'snodens_hh',timeStamp,'snwdens', [RHOSROW(I)])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'rof_hh',timeStamp,'mrro', [ROFROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofo_hh',timeStamp,'mrros', [ROFOROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofs_hh',timeStamp,'mrroi', [ROFSROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofb_hh',timeStamp,'mrrob', [ROFBROW(I)])

      call writeOutput1D(lonLocalIndex,latLocalIndex,'cdh_hh',timeStamp,'cdh', [CDHROW(I)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'cdm_hh',timeStamp,'cdm', [CDMROW(I)])
      !             call writeOutput1D(lonLocalIndex,latLocalIndex,'windu_hh',timeStamp,'windu', [SFCUROW(I)])  ! name
      !             call writeOutput1D(lonLocalIndex,latLocalIndex,'windv_hh',timeStamp,'windv', [SFCVROW(I)])  ! name
      !
      call writeOutput1D(lonLocalIndex,latLocalIndex,'tbar_hh',timeStamp,'tsl', [TBARROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'thlq_hh',timeStamp,'mrsll', [THLQROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'thic_hh',timeStamp,'mrsfl', [THICROW(I,:)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gflx_hh',timeStamp,'gflx', [GFLXROW(I,:)])  ! name

      ! If requested also prepare and write out the per tile physics variables
      if (dopertileoutput) then
        do m = 1,nmtest

          ALTOT = 0.0
          if (FSSROW(I) > 0.0) then ! there is no ROT for fss,as it will always be the same for all tiles.
            ALTOT = (FSSROW(I) - (FSGVROT(I,M) + FSGSROT(I,M) + FSGGROT(I,M))) / FSSROW(I)
          end if
          FSSTAR = FSSROW(I) * (1.0 - ALTOT)
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_hh_t' ,timeStamp,'rss', [FSSTAR])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'altot_hh_t' ,timeStamp,'albs', [ALTOT])

          FLSTAR = FDLROW(I) - SBC * GTROT(I,M) ** 4 ! there is no ROT for fdl,as it will always be the same for all tiles.
          call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_hh_t' ,timeStamp,'rls', [FLSTAR])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_hh_t'     ,timeStamp,'hfss', [HFSROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_hh_t'     ,timeStamp,'hfls', [QEVPROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'snm_hh_t'   ,timeStamp,'snm', [HMFNROT(I,M)])

          ZSN = 0.0
          if (RHOSROT(I,M) > 0.0) ZSN = SNOROT(I,M) / RHOSROT(I,M)
          call writeOutput1D(lonLocalIndex,latLocalIndex,'sisnthick_hh_t'   ,timeStamp,'sisnthick', [ZSN])

          TCN = 0.0
          if (TCANROT(I,M) > 0.01) TCN = TCANROT(I,M) - TFREZ
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tcs_hh_t'   ,timeStamp,'tcs', [TCN])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'scan_hh_t'   ,timeStamp,'scanopy', [SCANROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rcan_hh_t'   ,timeStamp,'rcanopy', [RCANROT(I,M)])

          ! TACGAT(I)-TFREZ  ! temp of air within canopy.
          TSN = 0.0
          if (TSNOROT(I,M) > 0.01) TSN = TSNOROT(I,M) - TFREZ
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tsno_hh_t'   ,timeStamp,'tsn', [TSN])

          TPN = 0.0
          if (TPNDROT(I,M) > 0.01) TPN = TPNDROT(I,M) - TFREZ
          call writeOutput1D(lonLocalIndex,latLocalIndex,'tpond_hh_t'   ,timeStamp,'tpond', [TPN])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'zpond_hh_t'   ,timeStamp,'zpond', [ZPNDROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gt_hh_t'   ,timeStamp,'tsblack', [GTROT(I,M) - TFREZ])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'sno_hh_t' ,timeStamp,'snw', [SNOROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_hh_t',timeStamp,'wsnw', [WSNOROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'snodens_hh_t',timeStamp,'snwdens', [RHOSROT(I,M)])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'rof_hh_t',timeStamp,'mrro', [ROFROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rofo_hh_t',timeStamp,'mrros', [ROFOROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rofs_hh_t',timeStamp,'mrroi', [ROFSROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rofb_hh_t',timeStamp,'mrrob', [ROFBROT(I,M)])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'cdh_hh_t',timeStamp,'cdh', [CDHROT(I,M)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'cdm_hh_t',timeStamp,'cdm', [CDMROT(I,M)])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'windu_hh_t',timeStamp,'windu', [SFCUROT(I,M)])  ! name
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'windv_hh_t',timeStamp,'windv', [SFCVROT(I,M)])  ! name

          EVAPSUM = QFCFROT(I,M) + QFCLROT(I,M) + QFNROT(I,M) &
                    + QFGROT(I,M) + QFCROT(I,M,1) + QFCROT(I,M,2) + QFCROT(I,M,3) ! FLAG this only considers the top 3 layers !!
          call writeOutput1D(lonLocalIndex,latLocalIndex,'evspsbl_hh_t'   ,timeStamp,'evspsbl', [EVAPSUM])

          call writeOutput1D(lonLocalIndex,latLocalIndex,'tbar_hh_t',timeStamp,'tsl', [TBARROT(I,M,:)])  ! name
          call writeOutput1D(lonLocalIndex,latLocalIndex,'thlq_hh_t',timeStamp,'mrsll', [THLQROT(I,M,:)])  ! name
          call writeOutput1D(lonLocalIndex,latLocalIndex,'thic_hh_t',timeStamp,'mrsfl', [THICROT(I,M,:)])  ! name
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gflx_hh_t',timeStamp,'gflx', [GFLXROT(I,M,:)])  ! name

        end do
      end if

      ! Other ones to make outputs with later.
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap',timeStamp,'evspsblsoi', [GROUNDEVAP(I)])
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap',timeStamp,'evspsblveg', [CANOPYEVAP(I)])
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_mo',timeStamp,'evspsbl', [EVAPACC_MO(I)])
      !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_mo',timeStamp,'tran', [TRANSPACC_MO(I)])
      ! TROFROT(I,M) !< Temperature of total runoff [K]
      ! TROOROT(I,M) !< Temperature of overland flow from top of soil column [K]
      ! TROSROT(I,M) !< Temperature of interflow from sides of soil column [K]
      ! TROBROT(I,M) !< Temperature of base flow from bottom of soil column [K]

      !                                 !
      !                                     &                   FCS(M),FGS(M),FC(M),FG(M),'
      !                                 WRITE(68,6800) IHOUR,IMIN,IDAY,realyr, &
      !                                     &                   FSGVROT(I,M),FSGSROT(I,M),FSGGROT(I,M), &
      !                                     &                   FLGVROT(I,M),FLGSROT(I,M),FLGGROT(I,M), &
      !                                     &                   HFSCROT(I,M),HFSSROT(I,M),HFSGROT(I,M), &
      !                                     &                   HEVCROT(I,M),HEVSROT(I,M),HEVGROT(I,M), &
      !                                     &                   HMFCROT(I,M),HMFNROT(I,M), &
      !                                     &                   (HMFGROT(I,M,J),J=1,3), &
      !                                     &                   HTCCROT(I,M),HTCSROT(I,M), &
      !                                     &                   (HTCROT(I,M,J),J=1,3),' TILE ',M
      !                                 WRITE(69,6900) IHOUR,IMIN,IDAY,realyr, &
      !                                     &                   PCFCROT(I,M),PCLCROT(I,M),PCPNROT(I,M), &
      !                                     &                   PCPGROT(I,M),QFCFROT(I,M),QFCLROT(I,M), &
      !                                     &                   QFNROT(I,M),QFGROT(I,M),(QFCROT(I,M,J),J=1,3), &
      !                                     &                   ROFCROT(I,M),ROFNROT(I,M),
      !                                     &                   ROFROT(I,M),WTRCROT(I,M),WTRSROT(I,M), &
      !                                     &                   WTRGROT(I,M),' TILE ',M

      ! Write half-hourly CTEM results to file
      !
      ! Net photosynthetic rates (GPP) and leaf maintenance respiration for each pft. however, if ctem_on then physyn subroutine
      ! is using storage lai while actual lai is zero. if actual lai is zero then we make anveg and rmlveg zero as well because these
      ! are imaginary just like storage lai. note that anveg and rmlveg are not passed to ctem. rather ancsveg, ancgveg, rmlcsveg, and
      ! rmlcgveg are passed.
      !
      allocate(anveggrd(icc),rmlveggrd(icc),fcanctot(icc))
      if (ctem_on) then
        do m = 1,nmtest
          do j = 1,icc
            if (ailcgrow(i,m,j) <= 0.0) then
              anvegrow(i,m,j) = 0.0
              rmlvegrow(i,m,j) = 0.0
            else
              anvegrow(i,m,j) = ancsvegrow(i,m,j) * FSNOROT(i,m) + ancgvegrow(i,m,j) * (1. - FSNOROT(i,m))
              rmlvegrow(i,m,j) = rmlcsvegrow(i,m,j) * FSNOROT(i,m) + rmlcgvegrow(i,m,j) * (1. - FSNOROT(i,m))
            end if
          end do
          if (dopertileoutput) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_hh_t',timeStamp,'gpp', [anvegrow(I,M,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'rml_hh_t',timeStamp,'rmLeaf', [rmlvegrow(I,M,:)])
          end if
          anveggrd = 0.0
          rmlveggrd = 0.0
          fcanctot = 0.0
          do j = 1,icc
            anveggrd(j) = anveggrd(j) + anvegrow(i,m,j) * FAREROT(i,m)
            rmlveggrd(j) = rmlveggrd(j) + rmlvegrow(i,m,j) * FAREROT(i,m)
            fcanctot(j) = fcanctot(j) + fcancmxrow(i,m,j) * FAREROT(i,m)
          end do
        end do ! m loop
        if (doperpftoutput) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_hh',timeStamp,'gpp', [anveggrd(:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rml_hh',timeStamp,'rmLeaf', [rmlveggrd(:)])
        end if
        ! Grid cell level outputs
        an_grd = 0.
        rml_grd = 0.
        totvegarea = sum(fcanctot)
        do j = 1,icc
          an_grd  = an_grd + anveggrd(j)  * fcanctot(j) / totvegarea
          rml_grd = rml_grd + rmlveggrd(j) * fcanctot(j) / totvegarea
        end do
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_hh_g',timeStamp,'gpp', [an_grd])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rml_hh_g',timeStamp,'rmLeaf', [rml_grd])
      end if ! ctem_on
    end do ! nltest loop

  end subroutine class_hh_w
  !< @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_class_daily_aw
  !> @{
  !> Accumlates and writes the daily physics variables. These are kept in pointer structures as
  !! this subroutine is called each physics timestep and we increment the timestep values to produce a daily value.
  !! The pointer to the daily data structures (in classStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine class_daily_aw (lonLocalIndex, latLocalIndex, iday, nltest, nmtest, ncount, nday, lastDOY, realyr)

    use classStateVars, only : class_rot, resetAccVars
    use classicParams,  only : ignd, sbc, delt, TFREZ
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: iday
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: ncount
    integer, intent(in) :: nday
    integer, intent(in) :: lastDOY
    integer, intent(in) :: realyr

    ! local variables
    integer :: i, m, j, k
    real, dimension(1) :: timeStamp
    real :: FSSTAR, FLSTAR
    real, allocatable, dimension(:) :: ALIRACC !< Diagnosed total near-infrared albedo of land surface [ ]
    real, allocatable, dimension(:) :: ALVSACC !< Diagnosed total visible albedo of land surface [ ]
    real, allocatable, dimension(:) :: EVAPACC !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLINACC !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLUTACC !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSINACC !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GROACC  !< Vegetation growth index [ ]
    real, allocatable, dimension(:) :: GTACC   !< Diagnosed effective surface black-body temperature [K]
    real, allocatable, dimension(:) :: HFSACC  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HMFNACC !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: OVRACC  !< Overland flow from top of soil column \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: PREACC  !< Surface precipitation rate \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: PRESACC !< Surface air pressure [Pa]
    real, allocatable, dimension(:) :: QAACC   !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: QEVPACC !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: RCANACC !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: RHOSACC !< Density of snow \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: ROFACC  !< Total runoff from soil \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: SCANACC !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: SNOACC  !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: TAACC   !< Air temperature at reference height [K]
    real, allocatable, dimension(:) :: TCANACC !< Vegetation canopy temperature [K]
    real, allocatable, dimension(:) :: TSNOACC !< Snowpack temperature [K]
    real, allocatable, dimension(:) :: UVACC   !< Wind speed \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: WSNOACC !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: WTBLACC !< Depth of water table in soil [m]
    real, allocatable, dimension(:) :: ALTOTACC !< Broadband albedo [-]
    real, allocatable, dimension(:) :: ALSNOACC !< Snow albedo [-]
    real, allocatable, dimension(:,:) :: TBARACC  !< Temperature of soil layers [K] (accumulated)
    real, allocatable, dimension(:,:) :: THLQACC  !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated)
    real, allocatable, dimension(:,:) :: THICACC  !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated)

    real, allocatable, dimension(:,:) :: UVACC_M
    real, allocatable, dimension(:,:) :: PRESACC_M
    real, allocatable, dimension(:,:) :: QAACC_M

    ! pointers
    integer, pointer, dimension(:) :: altotcntr_d   !< Used to count the number of time steps with the sun above the horizon
    real, pointer, dimension(:,:) :: FAREROT !< Fractional coverage of mosaic tile on modelled area
    real, pointer, dimension(:,:) :: FSGVROT        !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: FSGGROT        !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: FSGSROT        !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: GTROT          !< Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:,:) :: QEVPROT        !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: QFSROT         !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: HFSROT         !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HMFNROT        !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: ROFROT         !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:) :: ROFOROT        !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:,:,:) :: TBARROT      !< Temperature of soil layers [K]
    real, pointer, dimension(:,:,:) :: THICROT      !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:,:) :: THLQROT      !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: ALIRROT        !< Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer, dimension(:,:) :: ALVSROT        !< Diagnosed total visible albedo of land surface [ ]
    real, pointer, dimension(:,:) :: ALBSROT        !< Snow albedo [ ]
    real, pointer, dimension(:,:) :: WSNOROT        !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: SNOROT         !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, pointer, dimension(:,:) :: RHOSROT        !< Density of snow \f$[kg m^{-3}]\f$
    real, pointer, dimension(:,:) :: TSNOROT        !< Snowpack temperature [K]
    real, pointer, dimension(:,:) :: TCANROT        !< Vegetation canopy temperature [K]
    real, pointer, dimension(:,:) :: RCANROT        !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: SCANROT        !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: GROROT         !< Vegetation growth index [ ]

    real, pointer, dimension(:) :: FSSROW           !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FDLROW           !< Downwelling longwave sky radiation \f$[W m^{-2} ]\f$
    real, dimension(:), pointer :: PREROW           !< Surface precipitation rate \f$[kg m^{-2}  s^{-1} ]\f$
    real, dimension(:), pointer :: FSVHROW          !< Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, dimension(:), pointer :: FSIHROW          !< Near infrared shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    real, dimension(:), pointer :: TAROW            !< Air temperature at reference height [K]

    real, pointer, dimension(:,:) :: ALTOTACC_M     !< Broadband albedo [-] (accumulated)
    real, pointer, dimension(:,:) :: ALSNOACC_M     !< Snow albedo [-] (accumulated)
    real, pointer, dimension(:,:) :: PREACC_M       !< Surface precipitation rate \f$[kg m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: GTACC_M        !< Diagnosed effective surface black-body temperature [K] (accumulated)
    real, pointer, dimension(:,:) :: QEVPACC_M      !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: EVAPACC_M      !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: HFSACC_M       !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: HMFNACC_M      !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: ROFACC_M       !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: OVRACC_M       !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: WTBLACC_M
    real, pointer, dimension(:,:,:) :: TBARACC_M    !< Temperature of soil layers [K] (accumulated)
    real, pointer, dimension(:,:,:) :: THLQACC_M    !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated)
    real, pointer, dimension(:,:,:) :: THICACC_M    !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: ALVSACC_M      !< Diagnosed total visible albedo of land surface [ ] (accumulated)
    real, pointer, dimension(:,:) :: ALIRACC_M      !< Diagnosed total near-infrared albedo of land surface [ ] (accumulated)
    real, pointer, dimension(:,:) :: WSNOACC_M      !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: TCANACC_M      !< Vegetation canopy temperature [K] (accumulated)
    real, pointer, dimension(:,:) :: SNOACC_M       !< Mass of snow pack \f$[kg m^{-2}]\f$ (accumulated)
    real, pointer, dimension(:,:) :: RHOSACC_M      !< Density of snow \f$[kg m^{-3}]\f$ (accumulated)
    real, pointer, dimension(:,:) :: TSNOACC_M      !< Snowpack temperature [K] (accumulated)
    real, pointer, dimension(:,:) :: RCANACC_M      !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: SCANACC_M      !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: GROACC_M       !< Vegetation growth index [ ] (accumulated)
    real, pointer, dimension(:,:) :: FSINACC_M      !< Shortwave radiation \f$[W m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: FLINACC_M      !< Downwelling longwave sky radiation \f$[W m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: FLUTACC_M      !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$ (accumulated)
    real, pointer, dimension(:,:) :: TAACC_M        !< Air temperature at reference height [K] (accumulated)


    FSSROW => class_rot%FSSROW
    PREROW  => class_rot%PREROW
    FSVHROW => class_rot%FSVHROW
    FSIHROW => class_rot%FSIHROW
    FDLROW => class_rot%FDLROW
    TAROW  => class_rot%TAROW
    TBARROT=> class_rot%TBARROT
    THICROT=> class_rot%THICROT
    THLQROT=> class_rot%THLQROT
    SNOROT => class_rot%SNOROT
    HFSROT => class_rot%HFSROT
    QEVPROT => class_rot%QEVPROT
    QFSROT => class_rot%QFSROT
    GTROT => class_rot%GTROT
    ALIRROT => class_rot%ALIRROT
    ALVSROT => class_rot%ALVSROT
    ALBSROT => class_rot%ALBSROT
    FSGGROT => class_rot%FSGGROT
    FSGSROT => class_rot%FSGSROT
    FSGVROT => class_rot%FSGVROT
    HMFNROT => class_rot%HMFNROT
    ROFROT => class_rot%ROFROT
    GROROT => class_rot%GROROT
    ROFOROT => class_rot%ROFOROT
    RHOSROT => class_rot%RHOSROT
    TSNOROT=> class_rot%TSNOROT
    WSNOROT => class_rot%WSNOROT
    TCANROT=> class_rot%TCANROT
    RCANROT => class_rot%RCANROT
    SCANROT => class_rot%SCANROT
    FAREROT=> class_rot%FAREROT
    PREACC_M  => class_rot%PREACC_M
    GTACC_M   => class_rot%GTACC_M
    QEVPACC_M => class_rot%QEVPACC_M
    HFSACC_M  => class_rot%HFSACC_M
    HMFNACC_M => class_rot%HMFNACC_M
    ROFACC_M  => class_rot%ROFACC_M
    SNOACC_M  => class_rot%SNOACC_M
    OVRACC_M  => class_rot%OVRACC_M
    WTBLACC_M => class_rot%WTBLACC_M
    TBARACC_M => class_rot%TBARACC_M
    THLQACC_M => class_rot%THLQACC_M
    THICACC_M => class_rot%THICACC_M
    ALVSACC_M => class_rot%ALVSACC_M
    ALIRACC_M => class_rot%ALIRACC_M
    RHOSACC_M => class_rot%RHOSACC_M
    TSNOACC_M => class_rot%TSNOACC_M
    WSNOACC_M => class_rot%WSNOACC_M
    TCANACC_M => class_rot%TCANACC_M
    RCANACC_M => class_rot%RCANACC_M
    SCANACC_M => class_rot%SCANACC_M
    GROACC_M  => class_rot%GROACC_M
    FSINACC_M => class_rot%FSINACC_M
    FLINACC_M => class_rot%FLINACC_M
    TAACC_M   => class_rot%TAACC_M
    !         UVACC_M  => class_rot%UVACC_M
    !         PRESACC_M => class_rot%PRESACC_M
    !         QAACC_M  => class_rot%QAACC_M
    ALTOTACC_M => class_rot%ALTOTACC_M
    ALSNOACC_M => class_rot%ALSNOACC_M
    EVAPACC_M   => class_rot%EVAPACC_M
    FLUTACC_M   => class_rot%FLUTACC_M
    altotcntr_d => class_rot%altotcntr_d

    ! Accumulate output data for diurnally averaged fields. Both grid mean and mosaic mean
    do I = 1,NLTEST
      do M = 1,NMTEST
        if (FSSROW(I) > 0.) then
          ALTOTACC_M(I,M) = ALTOTACC_M(I,M) + (FSSROW(I) - (FSGVROT(I,M) &
                            + FSGSROT(I,M) + FSGGROT(I,M)))/FSSROW(I)
          ALSNOACC_M(I,M) = ALSNOACC_M(I,M) + ALBSROT(I,M)
          if (i == 1) altotcntr_d(i) = altotcntr_d(i) + 1 ! only count once per gridcell, not per tile
        end if

        PREACC_M(I,M) = PREACC_M(I,M) + PREROW(I)
        GTACC_M(I,M) = GTACC_M(I,M) + GTROT(I,M)
        QEVPACC_M(I,M) = QEVPACC_M(I,M) + QEVPROT(I,M)
        EVAPACC_M(I,M) = EVAPACC_M(I,M) + QFSROT(I,M)
        HFSACC_M(I,M) = HFSACC_M(I,M) + HFSROT(I,M)
        HMFNACC_M(I,M) = HMFNACC_M(I,M) + HMFNROT(I,M)
        ROFACC_M(I,M) = ROFACC_M(I,M) + ROFROT(I,M)
        OVRACC_M(I,M) = OVRACC_M(I,M) + ROFOROT(I,M)
        ! WTBLACC_M(I,M)=WTBLACC_M(I,M)+wtableROT(I,M)  ! FLAG fix !
        do J = 1,IGND
          TBARACC_M(I,M,J) = TBARACC_M(I,M,J) + TBARROT(I,M,J)
          THLQACC_M(I,M,J) = THLQACC_M(I,M,J) + THLQROT(I,M,J)
          THICACC_M(I,M,J) = THICACC_M(I,M,J) + THICROT(I,M,J)
        end do
        ALVSACC_M(I,M) = ALVSACC_M(I,M) + ALVSROT(I,M) * FSVHROW(I)
        ALIRACC_M(I,M) = ALIRACC_M(I,M) + ALIRROT(I,M) * FSIHROW(I)
        if (SNOROT(I,M) > 0.0) then
          RHOSACC_M(I,M) = RHOSACC_M(I,M) + RHOSROT(I,M)
          TSNOACC_M(I,M) = TSNOACC_M(I,M) + TSNOROT(I,M)
          WSNOACC_M(I,M) = WSNOACC_M(I,M) + WSNOROT(I,M)
        end if
        if (TCANROT(I,M) > 0.5) then
          TCANACC_M(I,M) = TCANACC_M(I,M) + TCANROT(I,M)
        end if
        SNOACC_M(I,M) = SNOACC_M(I,M) + SNOROT(I,M)
        RCANACC_M(I,M) = RCANACC_M(I,M) + RCANROT(I,M)
        SCANACC_M(I,M) = SCANACC_M(I,M) + SCANROT(I,M)
        GROACC_M(I,M) = GROACC_M(I,M) + GROROT(I,M)
        FSINACC_M(I,M) = FSINACC_M(I,M) + FSSROW(I)  ! not per tile
        FLINACC_M(I,M) = FLINACC_M(I,M) + FDLROW(I)  ! not per tile
        FLUTACC_M(I,M) = FLUTACC_M(I,M) + SBC * GTROT(I,M) ** 4
        TAACC_M(I,M) = TAACC_M(I,M) + TAROW(I)  ! not per tile
        !                UVACC_M(I,M)=UVACC_M(I,M)+UVROW(I)  ! not per tile
        !                 PRESACC_M(I,M)=PRESACC_M(I,M)+PRESROW(I)
        !                 QAACC_M(I,M)=QAACC_M(I,M)+QAROW(I)
      end do ! loop 50
    end do ! loop 75

    if (NCOUNT == NDAY) then

      allocate(ALIRACC(nltest),ALVSACC(nltest),EVAPACC(nltest),FLINACC(nltest), &
        FLUTACC(nltest),FSINACC(nltest),GROACC(nltest),GTACC(nltest), &
        HFSACC(nltest),HMFNACC(nltest),OVRACC(nltest),PREACC(nltest), &
        PRESACC(nltest),QAACC(nltest),QEVPACC(nltest),RCANACC(nltest), &
        RHOSACC(nltest),ROFACC(nltest),SCANACC(nltest),SNOACC(nltest), &
        TAACC(nltest),TCANACC(nltest),TSNOACC(nltest),UVACC(nltest), &
        WSNOACC(nltest),WTBLACC(nltest),ALTOTACC(nltest),TBARACC(nltest,ignd), &
        THLQACC(nltest,ignd),THICACC(nltest,ignd),ALSNOACC(nltest))
      ! UVACC_M,PRESACC_M,QAACC_M

      ALIRACC(:) = 0.0 ; ALVSACC(:) = 0.0 ; EVAPACC(:) = 0.0 ;   FLINACC(:) = 0.0
      FLUTACC(:) = 0.0 ; FSINACC(:) = 0.0 ; GROACC(:) = 0.0 ;    GTACC(:) = 0.0
      HFSACC(:) = 0.0 ;  HMFNACC(:) = 0.0 ; OVRACC(:) = 0.0 ;    PREACC(:) = 0.0
      PRESACC(:) = 0.0 ;   QAACC(:) = 0.0 ; QEVPACC(:) = 0.0 ;   RCANACC(:) = 0.0
      RHOSACC(:) = 0.0 ;  ROFACC(:) = 0.0 ; SCANACC(:) = 0.0 ;   SNOACC(:) = 0.0
      TAACC(:) = 0.0   ; TCANACC(:) = 0.0 ; TSNOACC(:) = 0.0 ;   UVACC(:) = 0.0
      WSNOACC(:) = 0.0 ; WTBLACC(:) = 0.0 ; ALTOTACC(:) = 0.0 ; THLQACC(:,:) = 0.0
      THICACC(:,:) = 0.0 ; TBARACC(:,:) = 0.0 ; ALSNOACC(:) = 0.0


      do I = 1,NLTEST
        do M = 1,NMTEST
          PREACC(I) = PREACC(I) + PREACC_M(I,M) * FAREROT(I,M)
          GTACC(I) = GTACC(I) + GTACC_M(I,M) * FAREROT(I,M)
          QEVPACC(I) = QEVPACC(I) + QEVPACC_M(I,M) * FAREROT(I,M)
          EVAPACC(I) = EVAPACC(I) + EVAPACC_M(I,M) * FAREROT(I,M)
          HFSACC(I) = HFSACC(I) + HFSACC_M(I,M) * FAREROT(I,M)
          HMFNACC(I) = HMFNACC(I) + HMFNACC_M(I,M) * FAREROT(I,M)
          ROFACC(I) = ROFACC(I) + ROFACC_M(I,M) * FAREROT(I,M)
          OVRACC(I) = OVRACC(I) + OVRACC_M(I,M) * FAREROT(I,M)
          WTBLACC(I) = WTBLACC(I) + WTBLACC_M(I,M) * FAREROT(I,M)
          ALTOTACC(I) = ALTOTACC(I) + ALTOTACC_M(I,M) * FAREROT(I,M)
          ALSNOACC(I) = ALSNOACC(I) + ALSNOACC_M(I,M) * FAREROT(I,M)
          do J = 1,IGND
            TBARACC(I,J) = TBARACC(I,J) + TBARACC_M(I,M,J) * FAREROT(I,M)
            THLQACC(I,J) = THLQACC(I,J) + THLQACC_M(I,M,J) * FAREROT(I,M)
            THICACC(I,J) = THICACC(I,J) + THICACC_M(I,M,J) * FAREROT(I,M)
          end do
          ALVSACC(I) = ALVSACC(I) + ALVSACC_M(I,M) * FAREROT(I,M)
          ALIRACC(I) = ALIRACC(I) + ALIRACC_M(I,M) * FAREROT(I,M)
          RHOSACC(I) = RHOSACC(I) + RHOSACC_M(I,M) * FAREROT(I,M)
          TSNOACC(I) = TSNOACC(I) + TSNOACC_M(I,M) * FAREROT(I,M)
          WSNOACC(I) = WSNOACC(I) + WSNOACC_M(I,M) * FAREROT(I,M)
          TCANACC(I) = TCANACC(I) + TCANACC_M(I,M) * FAREROT(I,M)
          SNOACC(I) = SNOACC(I) + SNOACC_M(I,M) * FAREROT(I,M)
          RCANACC(I) = RCANACC(I) + RCANACC_M(I,M) * FAREROT(I,M)
          SCANACC(I) = SCANACC(I) + SCANACC_M(I,M) * FAREROT(I,M)
          GROACC(I) = GROACC(I) + GROACC_M(I,M) * FAREROT(I,M)
          FSINACC(I) = FSINACC(I) + FSINACC_M(I,M) * FAREROT(I,M)
          FLINACC(I) = FLINACC(I) + FLINACC_M(I,M) * FAREROT(I,M)
          FLUTACC(I) = FLUTACC(I) + FLUTACC_M(I,M) * FAREROT(I,M)
          TAACC(I) = TAACC(I) + TAACC_M(I,M) * FAREROT(I,M)
          ! UVACC(I)=UVACC(I)+UVACC_M(I)*FAREROT(I,M)
          ! PRESACC(I)=PRESACC(I)+PRESROW(I)*FAREROT(I,M)
          ! QAACC(I)=QAACC(I)+QAROW(I)*FAREROT(I,M)
        end do
      end do

      ! Now write to file the grid average values

      ! Prepare the timestamp for this month. Take one day off since it referenced to 01-01 of the refyr.
      timeStamp = consecDays + iday - 1

      do i = 1,nltest
        if (altotcntr_d(i) > 0) then
          ALTOTACC(I) = ALTOTACC(I)/real(altotcntr_d(i))
          ALSNOACC(I) = ALSNOACC(I)/real(altotcntr_d(i))
        else
          ALTOTACC(I) = 0.
          ALSNOACC(I) = 0.
        end if
        FSSTAR = FSINACC(I)/real(nday) * (1. - ALTOTACC(I))
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_d' ,timeStamp,'rss', [FSSTAR])

        FLSTAR = (FLINACC(I) - FLUTACC(I))/real(NDAY)
        call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_d' ,timeStamp,'rls', [FLSTAR])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_d'     ,timeStamp,'hfss', [HFSACC(I)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_d'     ,timeStamp,'hfls', [QEVPACC(I)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snm_d'    ,timeStamp,'snm', [HMFNACC(I)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tbaracc_d',timeStamp,'tsl', [(TBARACC(I,:)/real(NDAY))])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thlqacc_d',timeStamp,'mrsll', [THLQACC(I,:)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thicacc_d',timeStamp,'mrsfl', [THICACC(I,:)/real(NDAY)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snm_d' ,timeStamp,'snm', [HMFNACC(I)])
        !                 !                 BEG=FSSTAR+FLSTAR-QH-QE

        ! Below are obviously not daily ones but I copy in as a reminder of what likely will be put out daily.
        ! after these writeOutput1D statements is the listed daily file variables. It would be good to do those again
        !                 call writeOutput1D(lonLocalIndex,latLocalIndex,'gflx_hh',timeStamp,'gflx', [GFLXROW(I,:)])  ! name
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'snoacc_mo' ,timeStamp,'snw', [HMFNACC(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_mo',timeStamp,'wsnw', [WSNOACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_mo'  ,timeStamp,'tas', [TAACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap',timeStamp,'evspsblsoi', [GROUNDEVAP(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap',timeStamp,'evspsblveg', [CANOPYEVAP(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_mo' ,timeStamp,'mrro', [ROFACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_mo' ,timeStamp,'pr', [PREACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_mo',timeStamp,'evspsbl', [EVAPACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_mo',timeStamp,'tran', [TRANSPACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'alsnoacc',timeStamp,'albsn', [ALSNOACC(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_mo',timeStamp,'albs', [ALTOTACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_mo',timeStamp,'actlyr', [ACTLYR_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_max_mo',timeStamp,'actlyrmax', [ACTLYR_MAX_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_min_mo',timeStamp,'actlyrmin', [ACTLYR_MIN_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_mo',timeStamp,'ftable', [FTABLE_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_max_mo',timeStamp,'ftablemax', [FTABLE_MAX_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_min_mo',timeStamp,'ftablemin', [FTABLE_MIN_MO(I)])
        !
        ! !                             &                       BEG,GTOUT,SNOACC(I),RHOSACC(I), &
        ! !                             &                       WSNOACC(I),ALTOTACC(I),ROFACC(I),CUMSNO
        ! !                             WRITE(62,6200) IDAY,realyr,(TBARACC(I,J)-TFREZ, &
        ! !                                 &                       THLQACC(I,J),THICACC(I,J),J=1,3), &
        ! !                                 &                       TCN,RCANACC(I),SCANACC(I),TSN,ZSN, &
        ! !                                 &                       ACTLYR_G(I),FTABLE_g(I)
        ! !                         WRITE(63,6300) IDAY,realyr,FSINACC(I),FLINACC(I), &
        ! !                             &                       TAACC(I)-TFREZ,UVACC(I),PRESACC(I), &
        ! !                             &                       QAACC(I),PREACC(I),EVAPACC(I)
        ! !                                 &                    BEG,GTOUT,SNOACC_M(I,M),RHOSACC_M(I,M), &
        ! !                                 &                    WSNOACC_M(I,M),ALTOTACC_M(I,M),ROFACC_M(I,M), &
        ! !                                 &                    CUMSNO,' TILE ',M
        ! !                                     &                  TCN,RCANACC_M(I,M),SCANACC_M(I,M),TSN,ZSN, &
        ! !                                     &                  ' TILE ',M
        ! !                             END IF
        ! !                             WRITE(631,6300) IDAY,realyr,FSINACC_M(I,M),FLINACC_M(I,M), &
        ! !                                 &                  TAACC_M(I,M)-TFREZ,UVACC_M(I,M),PRESACC_M(I,M), &
        ! !                                 &                  QAACC_M(I,M),PREACC_M(I,M),EVAPACC_M(I,M), &
        ! !                                 &                  ' TILE ',M
        !
        ! !
      end do
      ! !             do I=1,NLTEST
      ! !                 PREACC(I)=PREACC(I)
      ! !                 GTACC(I)=GTACC(I)/REAL(NDAY)
      ! !                 EVAPACC(I)=EVAPACC(I)
      ! !                 HMFNACC(I)=HMFNACC(I)/REAL(NDAY)
      ! !                 ROFACC(I)=ROFACC(I)
      ! !                 OVRACC(I)=OVRACC(I)
      ! !                 WTBLACC(I)=WTBLACC(I)/REAL(NDAY)
      ! !                 do J=1,IGND
      ! !                     TBARACC(I,J)=TBARACC(I,J)/REAL(NDAY)
      ! !                     THLQACC(I,J)=THLQACC(I,J)/REAL(NDAY)
      ! !                     THICACC(I,J)=THICACC(I,J)/REAL(NDAY)
      ! !                     THALACC(I,J)=THALACC(I,J)/REAL(NDAY)
      ! ! 725                     CONTINUE
      ! !                 IF (FSINACC(I)>0.0) THEN
      ! !                     ALVSACC(I)=ALVSACC(I)/(FSINACC(I)*0.5)
      ! !                     ALIRACC(I)=ALIRACC(I)/(FSINACC(I)*0.5)
      ! !                 ELSE
      ! !                     ALVSACC(I)=0.0
      ! !                     ALIRACC(I)=0.0
      ! !                 END IF
      ! !                 IF (SNOARE(I)>0.0) THEN
      ! !                     RHOSACC(I)=RHOSACC(I)/SNOARE(I)
      ! !                     TSNOACC(I)=TSNOACC(I)/SNOARE(I)
      ! !                     WSNOACC(I)=WSNOACC(I)/SNOARE(I)
      ! !                 END IF
      ! !                 IF (CANARE(I)>0.0) THEN
      ! !                     TCANACC(I)=TCANACC(I)/CANARE(I)
      ! !                 END IF
      ! !                 SNOACC(I)=SNOACC(I)/REAL(NDAY)
      ! !                 RCANACC(I)=RCANACC(I)/REAL(NDAY)
      ! !                 SCANACC(I)=SCANACC(I)/REAL(NDAY)
      ! !                 GROACC(I)=GROACC(I)/REAL(NDAY)
      ! !                 FSINACC(I)=FSINACC(I)/REAL(NDAY)
      ! !                 FLINACC(I)=FLINACC(I)/REAL(NDAY)
      ! !                 FLUTACC(I)=FLUTACC(I)/REAL(NDAY)
      ! !                 TAACC(I)=TAACC(I)/REAL(NDAY)
      ! !                 UVACC(I)=UVACC(I)/REAL(NDAY)
      ! !                 PRESACC(I)=PRESACC(I)/REAL(NDAY)
      ! !                 QAACC(I)=QAACC(I)/REAL(NDAY)
      ! !                 IF (RHOSACC(I)>0.0) THEN
      ! !                     ZSN=SNOACC(I)/RHOSACC(I)
      ! !                 ELSE
      ! !                     ZSN=0.0
      ! !                 END IF
      ! !                 IF (TCANACC(I)>0.01) THEN
      ! !                     TCN=TCANACC(I)-TFREZ
      ! !                 ELSE
      ! !                     TCN=0.0
      ! !                 END IF
      ! !                 IF (TSNOACC(I)>0.01) THEN
      ! !                     TSN=TSNOACC(I)-TFREZ
      ! !                 ELSE
      ! !                     TSN=0.0
      ! !                 END IF
      ! !                 GTOUT=GTACC(I)-TFREZ
      ! !!
      ! !
      ! !                     end if
      ! !                 END IF
      ! !
      ! !             do I=1,NLTEST
      ! !                 do M=1,NMTEST
      ! !                     PREACC_M(I,M)=PREACC_M(I,M)     ! became [kg m-2 day-1] instead of [kg m-2 s-1]
      ! !                     GTACC_M(I,M)=GTACC_M(I,M)/REAL(NDAY)
      ! !                     QEVPACC_M(I,M)=QEVPACC_M(I,M)/REAL(NDAY)
      ! !                     EVAPACC_M(I,M)=EVAPACC_M(I,M)   ! became [kg m-2 day-1] instead of [kg m-2 s-1]
      ! !                     HFSACC_M(I,M)=HFSACC_M(I,M)/REAL(NDAY)
      ! !                     HMFNACC_M(I,M)=HMFNACC_M(I,M)/REAL(NDAY)
      ! !                     ROFACC_M(I,M)=ROFACC_M(I,M)   ! became [kg m-2 day-1] instead of [kg m-2 s-1
      ! !                     OVRACC_M(I,M)=OVRACC_M(I,M)   ! became [kg m-2 day-1] instead of [kg m-2 s-1]
      ! !                     WTBLACC_M(I,M)=WTBLACC_M(I,M)/REAL(NDAY)
      ! !                     do J=1,IGND
      ! !                         TBARACC_M(I,M,J)=TBARACC_M(I,M,J)/REAL(NDAY)
      ! !                         THLQACC_M(I,M,J)=THLQACC_M(I,M,J)/REAL(NDAY)
      ! !                         THICACC_M(I,M,J)=THICACC_M(I,M,J)/REAL(NDAY)
      ! !                         THALACC_M(I,M,J)=THALACC_M(I,M,J)/REAL(NDAY)
      ! ! 726                         CONTINUE
      ! !
      ! !                     IF (FSINACC_M(I,M)>0.0) THEN
      ! !                         ALVSACC_M(I,M)=ALVSACC_M(I,M)/(FSINACC_M(I,M)*0.5)
      ! !                         ALIRACC_M(I,M)=ALIRACC_M(I,M)/(FSINACC_M(I,M)*0.5)
      ! !                     ELSE
      ! !                         ALVSACC_M(I,M)=0.0
      ! !                         ALIRACC_M(I,M)=0.0
      ! !                     END IF
      ! !
      ! !                     SNOACC_M(I,M)=SNOACC_M(I,M)/REAL(NDAY)
      ! !                     if (SNOARE_M(I,M) > 0.) THEN
      ! !                         RHOSACC_M(I,M)=RHOSACC_M(I,M)/SNOARE_M(I,M)
      ! !                         TSNOACC_M(I,M)=TSNOACC_M(I,M)/SNOARE_M(I,M)
      ! !                         WSNOACC_M(I,M)=WSNOACC_M(I,M)/SNOARE_M(I,M)
      ! !                     END IF
      ! !                     TCANACC_M(I,M)=TCANACC_M(I,M)/REAL(NDAY)
      ! !                     RCANACC_M(I,M)=RCANACC_M(I,M)/REAL(NDAY)
      ! !                     SCANACC_M(I,M)=SCANACC_M(I,M)/REAL(NDAY)
      ! !                     GROACC_M(I,M)=GROACC_M(I,M)/REAL(NDAY)
      ! !                     FSINACC_M(I,M)=FSINACC_M(I,M)/REAL(NDAY)
      ! !                     FLINACC_M(I,M)=FLINACC_M(I,M)/REAL(NDAY)
      ! !                     FLUTACC_M(I,M)=FLUTACC_M(I,M)/REAL(NDAY)
      ! !                     TAACC_M(I,M)=TAACC_M(I,M)/REAL(NDAY)
      ! !                     UVACC_M(I,M)=UVACC_M(I,M)/REAL(NDAY)
      ! !                     PRESACC_M(I,M)=PRESACC_M(I,M)/REAL(NDAY)
      ! !                     QAACC_M(I,M)=QAACC_M(I,M)/REAL(NDAY)
      ! !                     if (altotcntr_d(i) > 0) then ! altotcntr_d(i) could be 0
      ! !                         ALTOTACC_M(I,M)=ALTOTACC_M(I,M)/REAL(altotcntr_d(i))
      ! !                     else
      ! !                         ALTOTACC_M(I,M)=0.
      ! !                     end if
      ! !                     FSSTAR=FSINACC_M(I,M)*(1.-ALTOTACC_M(I,M))
      ! !                     FLSTAR=FLINACC_M(I,M)-FLUTACC_M(I,M)
      ! !                     QH=HFSACC_M(I,M)
      ! !                     QE=QEVPACC_M(I,M)
      ! !                     QEVPACC_M_SAVE(I,M)=QEVPACC_M(I,M)   ! FLAG ! What is the point of this? JM Apr 12015
      ! !                     BEG=FSSTAR+FLSTAR-QH-QE
      ! !                     SNOMLT=HMFNACC_M(I,M)
      ! !
      ! !                     IF (RHOSACC_M(I,M)>0.0) THEN
      ! !                         ZSN=SNOACC_M(I,M)/RHOSACC_M(I,M)
      ! !                     ELSE
      ! !                         ZSN=0.0
      ! !                     END IF
      ! !
      ! !                     IF (TCANACC_M(I,M)>0.01) THEN
      ! !                         TCN=TCANACC_M(I,M)-TFREZ
      ! !                     ELSE
      ! !                         TCN=0.0
      ! !                     END IF
      ! !
      ! !                     IF (TSNOACC_M(I,M)>0.01) THEN
      ! !                         TSN=TSNOACC_M(I,M)-TFREZ
      ! !                     ELSE
      ! !                         TSN=0.0
      ! !                     END IF
      ! !
      ! !                     GTOUT=GTACC_M(I,M)-TFREZ
      ! !
      ! !                             !         WRITE TO OUTPUT FILES
      ! !                             !
      ! !                         end if
      ! !
      ! ! 809                     CONTINUE
      ! ! 808                 CONTINUE

      deallocate(ALIRACC, ALVSACC, EVAPACC, FLINACC, FLUTACC, FSINACC, GROACC, GTACC, &
      HFSACC, HMFNACC, OVRACC, PREACC, PRESACC, QAACC, QEVPACC, RCANACC, &
      RHOSACC, ROFACC, SCANACC, SNOACC, TAACC, TCANACC, TSNOACC, UVACC, &
      WSNOACC, WTBLACC, ALTOTACC)
      ! UVACC_M, PRESACC_M, QAACC_M

      ! RESET ACCUMULATOR ARRAYS (*_M).
      call resetAccVars(nltest, nmtest)

    end if ! IF (NCOUNT==NDAY)

  end subroutine class_daily_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_class_monthly_aw
  !> @{
  !> Accumulate and write out the monthly physics outputs. These are kept in pointer structures as
  !! this subroutine is called each physics timestep and we increment the timestep values to produce a monthly value.
  !! The pointer to the monthly data structures (in classStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine class_monthly_aw (lonLocalIndex, latLocalIndex, IDAY, realyr, NCOUNT, NDAY, nltest, nmtest, lastDOY)

    use classStateVars, only : class_out, resetClassMon, class_rot
    use classicParams,  only : nmon, monthend, nmos, ignd, SBC, DELT, TFREZ
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: IDAY
    integer, intent(in) :: realyr
    integer, intent(in) :: NCOUNT
    integer, intent(in) :: NDAY
    integer, intent(in) :: lastDOY
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    ! pointers
    real, dimension(:,:,:), pointer :: TBARROT
    real, dimension(:,:,:), pointer :: THLQROT
    real, dimension(:,:,:), pointer :: THICROT
    real, dimension(:,:,:), pointer :: QFCROT
    real, dimension(:), pointer :: FSSROW
    real, dimension(:), pointer :: FDLROW
    real, dimension(:), pointer :: FSVHROW
    real, dimension(:), pointer :: FSIHROW
    real, dimension(:), pointer :: TAROW
    real, dimension(:), pointer :: PREROW
    real, dimension(:,:), pointer :: ALVSROT
    real, dimension(:,:), pointer :: ALBSROT
    real, dimension(:,:), pointer :: FAREROT
    real, dimension(:,:), pointer :: ALIRROT
    real, dimension(:,:), pointer :: GTROT
    real, dimension(:,:), pointer :: HFSROT
    real, dimension(:,:), pointer :: QEVPROT
    real, dimension(:,:), pointer :: groundHeatFluxROT !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    real, dimension(:,:), pointer :: SNOROT
    real, dimension(:,:), pointer :: WSNOROT
    real, dimension(:,:), pointer :: ROFROT
    real, dimension(:,:), pointer :: QFSROT
    real, dimension(:,:), pointer :: QFGROT
    real, dimension(:,:), pointer :: QFNROT
    real, dimension(:,:), pointer :: QFCLROT
    real, dimension(:,:), pointer :: QFCFROT
    real, dimension(:,:), pointer :: FSGVROT           !< Diagnosed net shortwave radiation on vegetation canopy
    real, dimension(:,:), pointer :: FSGSROT           !< Diagnosed net shortwave radiation on ground snow surface
    real, dimension(:,:), pointer :: FSGGROT           !< Diagnosed net shortwave radiation on ground surface
    real, pointer, dimension(:,:)  :: ftable      !< Depth to frozen water table (m)
    real, pointer, dimension(:,:)  :: actlyr      !< Active layer depth (m)
    real, pointer, dimension(:,:,:) :: dlzwrot    !< Permeable thickness of soil layer [m]
    real, pointer, dimension(:) :: ALVSACC_MO
    real, pointer, dimension(:) :: ALIRACC_MO
    real, pointer, dimension(:) :: FLUTACC_MO
    real, pointer, dimension(:) :: FSINACC_MO     !< Surface Downwelling Shortwave Radiative flux in air [$W m^{-2}$]
    real, pointer, dimension(:) :: FLINACC_MO
    real, pointer, dimension(:) :: HFSACC_MO
    real, pointer, dimension(:) :: QEVPACC_MO
    real, pointer, dimension(:) :: groundHeatFlux_MO  !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: SNOACC_MO
    real, pointer, dimension(:) :: WSNOACC_MO
    real, pointer, dimension(:) :: ROFACC_MO
    real, pointer, dimension(:) :: PREACC_MO
    real, pointer, dimension(:) :: EVAPACC_MO    !< Diagnosed total surface evaporation water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: TRANSPACC_MO
    real, pointer, dimension(:) :: TAACC_MO
    real, pointer, dimension(:) :: ACTLYR_MO
    real, pointer, dimension(:) :: FTABLE_MO
    real, pointer, dimension(:) :: ACTLYR_MIN_MO
    real, pointer, dimension(:) :: FTABLE_MIN_MO
    real, pointer, dimension(:) :: ACTLYR_MAX_MO
    real, pointer, dimension(:) :: FTABLE_MAX_MO
    real, pointer, dimension(:) :: ALTOTACC_MO
    real, pointer, dimension(:) :: ALSNOACC_MO
    real, pointer, dimension(:) :: GROUNDEVAP
    real, pointer, dimension(:) :: CANOPYEVAP
    real, pointer, dimension(:,:) :: TBARACC_MO
    real, pointer, dimension(:,:) :: THLQACC_MO
    real, pointer, dimension(:,:) :: THICACC_MO
    real, pointer, dimension(:) :: MRSO_MO
    real, pointer, dimension(:,:) :: MRSOL_MO
    integer, pointer, dimension(:) :: altotcntr_m

    ! local

    real :: ALTOT_MO
    integer :: NT
    integer :: NDMONTH
    integer :: i, m, j
    integer :: IMONTH
    real :: tovere
    real :: ACTLYR_tmp
    real :: FTABLE_tmp
    real :: FSSTAR_MO
    real :: FLSTAR_MO
    real :: QH_MO
    real :: QE_MO
    real, dimension(1) :: timeStamp

    ! point pointers
    TBARROT         => class_rot%TBARROT
    THLQROT         => class_rot%THLQROT
    THICROT         => class_rot%THICROT
    QFCROT          => class_rot%QFCROT
    ALVSROT         => class_rot%ALVSROT
    ALBSROT         => class_rot%ALBSROT
    FAREROT         => class_rot%FAREROT
    ALIRROT         => class_rot%ALIRROT
    GTROT           => class_rot%GTROT
    HFSROT          => class_rot%HFSROT
    QEVPROT         => class_rot%QEVPROT
    groundHeatFluxROT => class_rot%groundHeatFluxROT
    SNOROT          => class_rot%SNOROT
    WSNOROT         => class_rot%WSNOROT
    ROFROT          => class_rot%ROFROT
    QFSROT          => class_rot%QFSROT
    QFGROT          => class_rot%QFGROT
    QFNROT          => class_rot%QFNROT
    QFCLROT         => class_rot%QFCLROT
    QFCFROT         => class_rot%QFCFROT
    FSGVROT         => class_rot%FSGVROT
    FSGSROT         => class_rot%FSGSROT
    FSGGROT         => class_rot%FSGGROT
    FSSROW          => class_rot%FSSROW
    FDLROW          => class_rot%FDLROW
    FSVHROW         => class_rot%FSVHROW
    FSIHROW         => class_rot%FSIHROW
    TAROW           => class_rot%TAROW
    PREROW          => class_rot%PREROW
    ftable          => class_rot%ftable
    actlyr          => class_rot%actlyr
    dlzwrot         => class_rot%dlzwrot
    ALVSACC_MO        => class_out%ALVSACC_MO
    ALIRACC_MO        => class_out%ALIRACC_MO
    FLUTACC_MO        => class_out%FLUTACC_MO
    FSINACC_MO        => class_out%FSINACC_MO
    FLINACC_MO        => class_out%FLINACC_MO
    HFSACC_MO         => class_out%HFSACC_MO
    QEVPACC_MO        => class_out%QEVPACC_MO
    groundHeatFlux_MO => class_out%groundHeatFlux_MO
    SNOACC_MO         => class_out%SNOACC_MO
    WSNOACC_MO        => class_out%WSNOACC_MO
    ROFACC_MO         => class_out%ROFACC_MO
    PREACC_MO         => class_out%PREACC_MO
    EVAPACC_MO        => class_out%EVAPACC_MO
    TRANSPACC_MO      => class_out%TRANSPACC_MO
    TAACC_MO          => class_out%TAACC_MO
    TBARACC_MO        => class_out%TBARACC_MO
    THLQACC_MO        => class_out%THLQACC_MO
    THICACC_MO        => class_out%THICACC_MO
    ACTLYR_MO         => class_out%ACTLYR_MO
    FTABLE_MO         => class_out%FTABLE_MO
    ACTLYR_MIN_MO     => class_out%ACTLYR_MIN_MO
    FTABLE_MIN_MO     => class_out%FTABLE_MIN_MO
    ACTLYR_MAX_MO     => class_out%ACTLYR_MAX_MO
    FTABLE_MAX_MO     => class_out%FTABLE_MAX_MO
    GROUNDEVAP        => class_out%GROUNDEVAP
    CANOPYEVAP        => class_out%CANOPYEVAP
    ALTOTACC_MO       => class_out%ALTOTACC_MO
    ALSNOACC_MO       => class_out%ALSNOACC_MO
    altotcntr_m       => class_out%altotcntr_m
    MRSO_MO           => class_out%MRSO_MO
    MRSOL_MO           => class_out%MRSOL_MO

    ! ------------

    !> Accumulate output data for monthly averaged fields for class grid-mean.
    !> for both parallel mode and stand alone mode

    FSSTAR_MO   = 0.0
    FLSTAR_MO   = 0.0
    QH_MO       = 0.0
    QE_MO       = 0.0
    ACTLYR_tmp  = 0.0
    FTABLE_tmp  = 0.0

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do M = 1,NMTEST

      ! These are presently not being outputted but the code is kept in place if the need arises.
      !     ALVSACC_MO(I)=ALVSACC_MO(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
      !     ALIRACC_MO(I)=ALIRACC_MO(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)
      FLUTACC_MO(I) = FLUTACC_MO(I) + SBC * GTROT(I,M) ** 4 * FAREROT(I,M)
      FSINACC_MO(I) = FSINACC_MO(I) + FSSROW(I) * FAREROT(I,M)
      FLINACC_MO(I) = FLINACC_MO(I) + FDLROW(I) * FAREROT(I,M)
      HFSACC_MO(I) = HFSACC_MO(I) + HFSROT(I,M) * FAREROT(I,M)
      QEVPACC_MO(I) = QEVPACC_MO(I) + QEVPROT(I,M) * FAREROT(I,M)
      groundHeatFlux_MO(I) = groundHeatFlux_MO(I) + groundHeatFluxROT(I,M) * FAREROT(I,M) !*()*()*()*()*()*()
      SNOACC_MO(I) = SNOACC_MO(I) + SNOROT(I,M) * FAREROT(I,M)
      TAACC_MO(I) = TAACC_MO(I) + TAROW(I) * FAREROT(I,M)
      ACTLYR_MO(I) = ACTLYR_MO(I) + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_MO(I) = FTABLE_MO(I) + FTABLE(I,M) * FAREROT(I,M)
      ACTLYR_tmp = ACTLYR_tmp + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_tmp = FTABLE_tmp + FTABLE(I,M) * FAREROT(I,M)
      GROUNDEVAP(I) = GROUNDEVAP(I) + (QFGROT(I,M) + QFNROT(I,M)) * FAREROT(I,M) ! ground evap includes both evap and sublimation from snow
      CANOPYEVAP(I) = CANOPYEVAP(I) + (QFCLROT(I,M) + QFCFROT(I,M)) * FAREROT(I,M) ! canopy evap includes both evap and sublimation

      if (SNOROT(I,M) > 0.0) then
        WSNOACC_MO(I) = WSNOACC_MO(I) + WSNOROT(I,M) * FAREROT(I,M)
      end if

      ROFACC_MO(I) = ROFACC_MO(I) + ROFROT(I,M) * FAREROT(I,M)
      PREACC_MO(I) = PREACC_MO(I) + PREROW(I) * FAREROT(I,M)
      EVAPACC_MO(I) = EVAPACC_MO(I) + QFSROT(I,M) * FAREROT(I,M) ! Only evaporation

      if (FSSROW(I) > 0.0) then
        ALTOTACC_MO(I) = ALTOTACC_MO(I) + ( (FSSROW(I) - (FSGVROT(I,M) + FSGSROT(I,M) + FSGGROT(I,M))) &
                         /FSSROW(I) ) * FAREROT(I,M)
        ALSNOACC_MO(I) = ALSNOACC_MO(I) + ALBSROT(I,M)
        altotcntr_m(i) = altotcntr_m(i) + 1
      end if

      do J = 1,IGND
        TBARACC_MO(I,J) = TBARACC_MO(I,J) + TBARROT(I,M,J) * FAREROT(I,M)
        ! Convert from m3/m3 to kg/m2
        THLQACC_MO(I,J) = THLQACC_MO(I,J) + THLQROT(I,M,J) * FAREROT(I,M) * 1000. * DLZWROT(I,M,J)
        THICACC_MO(I,J) = THICACC_MO(I,J) + THICROT(I,M,J) * FAREROT(I,M) * 1000. * DLZWROT(I,M,J)
        ! Find the total soil moisture content
        ! Add up each soil layers moisture and convert from m3/m3 to kg/m2
        MRSO_MO(I) = MRSO_MO(I) + (THLQROT(I,M,J) * FAREROT(I,M) + THICROT(I,M,J) * FAREROT(I,M)) * 1000. * DLZWROT(I,M,J)
        MRSOL_MO(I,J) = MRSOL_MO(I,J) + (THLQROT(I,M,J) * FAREROT(I,M) + THICROT(I,M,J) * FAREROT(I,M)) * 1000. * DLZWROT(I,M,J)
        TRANSPACC_MO(I) = TRANSPACC_MO(I) + QFCROT(I,M,J) * FAREROT(I,M)

      end do ! loop 823

    end do ! loop 821

    ! Check if the active layer has become more shallow or deepened.
    ACTLYR_MAX_MO(I) = max(ACTLYR_MAX_MO(I),ACTLYR_tmp)
    ACTLYR_MIN_MO(I) = min(ACTLYR_MIN_MO(I),ACTLYR_tmp)
    FTABLE_MAX_MO(I) = max(FTABLE_MAX_MO(I),FTABLE_tmp)
    FTABLE_MIN_MO(I) = min(FTABLE_MIN_MO(I),FTABLE_tmp)

    do NT = 1,NMON
      if (IDAY == monthend(NT + 1) .and. NCOUNT == NDAY) then

        IMONTH = NT
        NDMONTH = (monthend(NT + 1) - monthend(NT)) * NDAY

        ! These are presently not being outputted but the code is kept in place if the need arises.
        !             IF (FSINACC_MO(I)>0.0) THEN
        !                 ALVSACC_MO(I)=ALVSACC_MO(I)/(FSINACC_MO(I)*0.5)
        !                 ALIRACC_MO(I)=ALIRACC_MO(I)/(FSINACC_MO(I)*0.5)
        !             ELSE
        !                 ALVSACC_MO(I)=0.0
        !                 ALIRACC_MO(I)=0.0
        !             END IF

        ! Albedo is only counted when sun is above horizon so it uses its own counter.\

        if (altotcntr_m(i) > 0) then
          ALTOTACC_MO(I) = ALTOTACC_MO(I)/real(altotcntr_m(i))
          ALSNOACC_MO(I) = ALSNOACC_MO(I)/real(altotcntr_m(i))
        else
          ALTOTACC_MO(I) = 0.
          ALSNOACC_MO(I) = 0.
        end if

        FLUTACC_MO(I) = FLUTACC_MO(I)/real(NDMONTH)
        FSINACC_MO(I) = FSINACC_MO(I)/real(NDMONTH)
        FLINACC_MO(I) = FLINACC_MO(I)/real(NDMONTH)
        HFSACC_MO(I) = HFSACC_MO(I)/real(NDMONTH)
        QEVPACC_MO(I) = QEVPACC_MO(I)/real(NDMONTH)
        groundHeatFlux_MO(I) = groundHeatFlux_MO(I)/real(NDMONTH)
        SNOACC_MO(I) = SNOACC_MO(I)/real(NDMONTH)
        WSNOACC_MO(I) = WSNOACC_MO(I)/real(NDMONTH)
        TAACC_MO(I) = TAACC_MO(I)/real(NDMONTH)
        ACTLYR_MO(I) = ACTLYR_MO(I)/real(NDMONTH)
        FTABLE_MO(I) = FTABLE_MO(I)/real(NDMONTH)
        MRSO_MO(I) = MRSO_MO(I)/real(NDMONTH)
        ROFACC_MO(I) = ROFACC_MO(I) /real(NDMONTH)
        PREACC_MO(I) = PREACC_MO(I) /real(NDMONTH)
        EVAPACC_MO(I) = EVAPACC_MO(I)/real(NDMONTH)
        TRANSPACC_MO(I) = TRANSPACC_MO(I)/real(NDMONTH)
        GROUNDEVAP(I) = GROUNDEVAP(I) /real(NDMONTH)
        CANOPYEVAP(I) = CANOPYEVAP(I) /real(NDMONTH)

        do J = 1,IGND
          TBARACC_MO(I,J) = TBARACC_MO(I,J)/real(NDMONTH)
          THLQACC_MO(I,J) = THLQACC_MO(I,J)/real(NDMONTH)
          THICACC_MO(I,J) = THICACC_MO(I,J)/real(NDMONTH)
          MRSOL_MO(I,J) = MRSOL_MO(I,J)/real(NDMONTH)
        end do

        FSSTAR_MO = FSINACC_MO(I) * (1. - ALTOTACC_MO(I))
        FLSTAR_MO = FLINACC_MO(I) - FLUTACC_MO(I)
        QH_MO = HFSACC_MO(I)
        QE_MO = QEVPACC_MO(I)

        tovere = 0.
        if (EVAPACC_MO(I) > 0.) tovere = TRANSPACC_MO(I)/EVAPACC_MO(I)

        ! Prepare the timestamp for this month. Take one day off so it is the last day of the month
        ! rather than the first day of the next month.
        timeStamp = consecDays + monthend(imonth + 1) - 1

        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsinacc_mo' ,timeStamp,'rsds', [FSINACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_mo' ,timeStamp,'rss', [FSSTAR_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_mo' ,timeStamp,'rls', [FLSTAR_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_mo'     ,timeStamp,'hfss', [QH_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_mo'     ,timeStamp,'hfls', [QE_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'hfg_mo'    ,timeStamp,'hfg', [groundHeatFlux_MO])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'snoacc_mo' ,timeStamp,'snw', [SNOACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'wsnoacc_mo',timeStamp,'wsnw', [WSNOACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'taacc_mo'  ,timeStamp,'tas', [TAACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'groundevap',timeStamp,'evspsblsoi', [GROUNDEVAP(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'canopyevap',timeStamp,'evspsblveg', [CANOPYEVAP(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'evspsbl_mo',timeStamp,'evspsbl', [CANOPYEVAP(I) + GROUNDEVAP(I) + TRANSPACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_mo' ,timeStamp,'mrro', [ROFACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_mo' ,timeStamp,'pr', [PREACC_MO(I)])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_mo',timeStamp,'evap', [EVAPACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_mo',timeStamp,'tran', [TRANSPACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'alsnoacc_mo',timeStamp,'albsn', [ALSNOACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_mo',timeStamp,'albs', [ALTOTACC_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'tbaracc_mo',timeStamp,'tsl', [TBARACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thlqacc_mo',timeStamp,'mrsll', [THLQACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'thicacc_mo',timeStamp,'mrsfl', [THICACC_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'mrso_mo',timeStamp,'mrso', [MRSO_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'mrsol_mo',timeStamp,'mrsol', [MRSOL_MO(I,:)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_mo',timeStamp,'actlyr', [ACTLYR_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_max_mo',timeStamp,'actlyrmax', [ACTLYR_MAX_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_min_mo',timeStamp,'actlyrmin', [ACTLYR_MIN_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_mo',timeStamp,'ftable', [FTABLE_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_max_mo',timeStamp,'ftablemax', [FTABLE_MAX_MO(I)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_min_mo',timeStamp,'ftablemin', [FTABLE_MIN_MO(I)])
        ! tovere

        call resetClassMon(nltest)

      end if ! if (IDAY==monthend(NT+1).AND.NCOUNT==NDAY)
    end do ! NMON

  end subroutine class_monthly_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_class_annual_aw
  !> @{
  !> Accumulate and write out the annual physics outputs. These are kept in pointer structures as
  !! this subroutine is called each physics timestep and we increment the timestep values to produce annuals values.
  !! The pointer to the annual data structures (in classStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine class_annual_aw (lonLocalIndex, latLocalIndex, IDAY, realyr, NCOUNT, NDAY, &
                              nltest, nmtest, lastDOY)

    use classStateVars, only : class_out, resetClassYr, class_rot
    use classicParams,  only : nmon, monthend, nmos, ignd, SBC, DELT
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: IDAY
    integer, intent(in) :: realyr
    integer, intent(in) :: NCOUNT
    integer, intent(in) :: NDAY
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: lastDOY

    ! pointers
    real, dimension(:), pointer :: FSSROW
    real, dimension(:), pointer :: FDLROW
    real, dimension(:), pointer :: FSVHROW
    real, dimension(:), pointer :: FSIHROW
    real, dimension(:), pointer :: TAROW
    real, dimension(:), pointer :: PREROW
    real, dimension(:,:), pointer :: ALVSROT
    real, dimension(:,:), pointer :: FAREROT
    real, dimension(:,:), pointer :: ALIRROT
    real, dimension(:,:), pointer :: GTROT
    real, dimension(:,:), pointer :: HFSROT
    real, dimension(:,:), pointer :: QEVPROT
    real, dimension(:,:), pointer :: ROFROT
    real, dimension(:,:), pointer :: QFSROT
    real, dimension(:,:,:), pointer :: QFCROT
    real, dimension(:,:), pointer :: FSGVROT           !< Diagnosed net shortwave radiation on vegetation canopy
    real, dimension(:,:), pointer :: FSGSROT           !< Diagnosed net shortwave radiation on ground snow surface
    real, dimension(:,:), pointer :: FSGGROT           !< Diagnosed net shortwave radiation on ground surface
    real, pointer, dimension(:,:)  :: ftable      !< Depth to frozen water table (m)
    real, pointer, dimension(:,:)  :: actlyr      !< Active layer depth (m)

    integer, pointer, dimension(:) :: altotcntr_yr
    real, pointer, dimension(:) :: ALVSACC_YR
    real, pointer, dimension(:) :: ALIRACC_YR
    real, pointer, dimension(:) :: FLUTACC_YR
    real, pointer, dimension(:) :: FSINACC_YR
    real, pointer, dimension(:) :: FLINACC_YR
    real, pointer, dimension(:) :: HFSACC_YR
    real, pointer, dimension(:) :: QEVPACC_YR
    real, pointer, dimension(:) :: ROFACC_YR
    real, pointer, dimension(:) :: PREACC_YR
    real, pointer, dimension(:) :: EVAPACC_YR
    real, pointer, dimension(:) :: TRANSPACC_YR
    real, pointer, dimension(:) :: TAACC_YR
    real, pointer, dimension(:) :: ACTLYR_YR
    real, pointer, dimension(:) :: ACTLYR_MIN_YR
    real, pointer, dimension(:) :: ACTLYR_MAX_YR
    real, pointer, dimension(:) :: FTABLE_YR
    real, pointer, dimension(:) :: FTABLE_MIN_YR
    real, pointer, dimension(:) :: FTABLE_MAX_YR
    real, pointer, dimension(:) :: ALTOTACC_YR

    ! local
    integer :: i, m, j
    real :: tovere
    real :: FSSTAR_YR
    real :: FLSTAR_YR
    real :: QH_YR
    real :: QE_YR
    real :: ACTLYR_tmp
    real :: FTABLE_tmp
    real, dimension(1) :: timeStamp

    ! point pointers
    ALVSROT         => class_rot%ALVSROT
    FAREROT         => class_rot%FAREROT
    ALIRROT         => class_rot%ALIRROT
    GTROT           => class_rot%GTROT
    HFSROT          => class_rot%HFSROT
    QEVPROT         => class_rot%QEVPROT
    ROFROT          => class_rot%ROFROT
    QFSROT          => class_rot%QFSROT
    QFCROT          => class_rot%QFCROT
    FSGVROT         => class_rot%FSGVROT
    FSGSROT         => class_rot%FSGSROT
    FSGGROT         => class_rot%FSGGROT
    FSSROW          => class_rot%FSSROW
    FDLROW          => class_rot%FDLROW
    FSVHROW         => class_rot%FSVHROW
    FSIHROW         => class_rot%FSIHROW
    TAROW           => class_rot%TAROW
    PREROW          => class_rot%PREROW
    ftable          => class_rot%ftable
    actlyr          => class_rot%actlyr
    ALVSACC_YR        => class_out%ALVSACC_YR
    ALIRACC_YR        => class_out%ALIRACC_YR
    FLUTACC_YR        => class_out%FLUTACC_YR
    FSINACC_YR        => class_out%FSINACC_YR
    FLINACC_YR        => class_out%FLINACC_YR
    HFSACC_YR         => class_out%HFSACC_YR
    QEVPACC_YR        => class_out%QEVPACC_YR
    ROFACC_YR         => class_out%ROFACC_YR
    PREACC_YR         => class_out%PREACC_YR
    EVAPACC_YR        => class_out%EVAPACC_YR
    TRANSPACC_YR      => class_out%TRANSPACC_YR
    TAACC_YR          => class_out%TAACC_YR
    ACTLYR_YR         => class_out%ACTLYR_YR
    ACTLYR_MIN_YR     => class_out%ACTLYR_MIN_YR
    ACTLYR_MAX_YR     => class_out%ACTLYR_MAX_YR
    FTABLE_YR         => class_out%FTABLE_YR
    FTABLE_MIN_YR     => class_out%FTABLE_MIN_YR
    FTABLE_MAX_YR     => class_out%FTABLE_MAX_YR
    ALTOTACC_YR       => class_out%ALTOTACC_YR
    altotcntr_yr      => class_out%altotcntr_yr

    !> Accumulate output data for yearly averaged fields for class grid-mean.
    !> for both parallel mode and stand alone mode
    FSSTAR_YR   = 0.0
    FLSTAR_YR   = 0.0
    QH_YR       = 0.0
    QE_YR       = 0.0
    ACTLYR_tmp  = 0.0
    FTABLE_tmp  = 0.0

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do M = 1,NMTEST

      ! These are presently not being outputted but the code is kept in place if the need arises.
      !         ALVSACC_YR(I)=ALVSACC_YR(I)+ALVSROT(I,M)*FAREROT(I,M)*FSVHROW(I)
      !         ALIRACC_YR(I)=ALIRACC_YR(I)+ALIRROT(I,M)*FAREROT(I,M)*FSIHROW(I)

      FLUTACC_YR(I) = FLUTACC_YR(I) + SBC * GTROT(I,M) ** 4 * FAREROT(I,M)
      FSINACC_YR(I) = FSINACC_YR(I) + FSSROW(I) * FAREROT(I,M)
      FLINACC_YR(I) = FLINACC_YR(I) + FDLROW(I) * FAREROT(I,M)
      HFSACC_YR(I) = HFSACC_YR(I) + HFSROT(I,M) * FAREROT(I,M)
      QEVPACC_YR(I) = QEVPACC_YR(I) + QEVPROT(I,M) * FAREROT(I,M)
      TAACC_YR(I) = TAACC_YR(I) + TAROW(I) * FAREROT(I,M)
      ROFACC_YR(I) = ROFACC_YR(I) + ROFROT(I,M) * FAREROT(I,M)
      PREACC_YR(I) = PREACC_YR(I) + PREROW(I) * FAREROT(I,M)
      EVAPACC_YR(I) = EVAPACC_YR(I) + QFSROT(I,M) * FAREROT(I,M)
      ACTLYR_YR(I) = ACTLYR_YR(I) + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_YR(I) = FTABLE_YR(I) + FTABLE(I,M) * FAREROT(I,M)
      ACTLYR_TMP = ACTLYR_TMP + ACTLYR(I,M) * FAREROT(I,M)
      FTABLE_TMP = FTABLE_TMP + FTABLE(I,M) * FAREROT(I,M)

      do J = 1,IGND
        TRANSPACC_YR(I) = TRANSPACC_YR(I) + QFCROT(I,M,J) * FAREROT(I,M)
      end do

      if (FSSROW(I) > 0.0) then
        ALTOTACC_YR(I) = ALTOTACC_YR(I) + ((FSSROW(I) - (FSGVROT(I,M) + FSGSROT(I,M) + FSGGROT(I,M))) &
                         /FSSROW(I) ) * FAREROT(I,M)
        altotcntr_yr(i) = altotcntr_yr(i) + 1
      end if

      ! Check if the active layer has become more shallow or deepened.
      ACTLYR_MAX_YR(I) = max(ACTLYR_MAX_YR(I),ACTLYR_tmp)
      ACTLYR_MIN_YR(I) = min(ACTLYR_MIN_YR(I),ACTLYR_tmp)
      FTABLE_MAX_YR(I) = max(FTABLE_MAX_YR(I),FTABLE_tmp)
      FTABLE_MIN_YR(I) = min(FTABLE_MIN_YR(I),FTABLE_tmp)

    end do ! loop 828

    if (IDAY == lastDOY .and. NCOUNT == NDAY) then

      ! These are presently not being outputted but the code is kept in place if the need arises.
      !             IF (FSINACC_YR(I)>0.0) THEN
      !                 ALVSACC_YR(I)=ALVSACC_YR(I)/(FSINACC_YR(I)*0.5)
      !                 ALIRACC_YR(I)=ALIRACC_YR(I)/(FSINACC_YR(I)*0.5)
      !             ELSE
      !                 ALVSACC_YR(I)=0.0
      !                 ALIRACC_YR(I)=0.0
      !             END IF

      FLUTACC_YR(I) = FLUTACC_YR(I)/(real(NDAY) * real(lastDOY))
      FSINACC_YR(I) = FSINACC_YR(I)/(real(NDAY) * real(lastDOY))
      FLINACC_YR(I) = FLINACC_YR(I)/(real(NDAY) * real(lastDOY))
      HFSACC_YR(I) = HFSACC_YR(I)/(real(NDAY) * real(lastDOY))
      QEVPACC_YR(I) = QEVPACC_YR(I)/(real(NDAY) * real(lastDOY))
      ROFACC_YR(I) = ROFACC_YR(I)/(real(NDAY) * real(lastDOY))
      PREACC_YR(I) = PREACC_YR(I)/(real(NDAY) * real(lastDOY))
      EVAPACC_YR(I) = EVAPACC_YR(I)/(real(NDAY) * real(lastDOY))
      TRANSPACC_YR(I) = TRANSPACC_YR(I)/(real(NDAY) * real(lastDOY))
      TAACC_YR(I) = TAACC_YR(I)/(real(NDAY) * real(lastDOY))
      ACTLYR_YR(I) = ACTLYR_YR(I)/(real(NDAY) * real(lastDOY))
      FTABLE_YR(I) = FTABLE_YR(I)/(real(NDAY) * real(lastDOY))

      ! Albedo is only counted when sun is above horizon so it uses its own counter.
      ALTOTACC_YR(I) = 0.
      if (altotcntr_yr(i) > 0) ALTOTACC_YR(I) = ALTOTACC_YR(I)/(real(altotcntr_yr(i)))

      FSSTAR_YR = FSINACC_YR(I) * (1. - ALTOTACC_YR(I))
      FLSTAR_YR = FLINACC_YR(I) - FLUTACC_YR(I)
      QH_YR = HFSACC_YR(I)
      QE_YR = QEVPACC_YR(I)

      tovere = 0.
      if (EVAPACC_YR(I) > 0.) tovere = TRANSPACC_YR(I)/EVAPACC_YR(I)

      ! Prepare the timestamp for this year
      timeStamp = consecDays

      call writeOutput1D(lonLocalIndex,latLocalIndex,'fsstar_yr' ,timeStamp,'rss', [FSSTAR_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'flstar_yr' ,timeStamp,'rls', [FLSTAR_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qh_yr'     ,timeStamp,'hfss', [QH_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'qe_yr'     ,timeStamp,'hfls', [QE_YR])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rofacc_yr' ,timeStamp,'mrro', [ROFACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'preacc_yr' ,timeStamp,'pr', [PREACC_YR(i)])
      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'evapacc_yr' ,timeStamp,'evap', [EVAPACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'transpacc_yr' ,timeStamp,'tran', [TRANSPACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'altotacc_yr' ,timeStamp,'albs', [ALTOTACC_YR(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_yr' ,timeStamp,'actlyr', [actlyr_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_yr' ,timeStamp,'ftable', [ftable_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_max_yr' ,timeStamp,'actlyrmax', [actlyr_max_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_max_yr' ,timeStamp,'ftablemax', [ftable_max_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'actlyr_min_yr' ,timeStamp,'actlyrmin', [actlyr_min_yr(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ftable_min_yr' ,timeStamp,'ftablemin', [ftable_min_yr(i)])

      !> ADD INITIALIZTION FOR YEARLY ACCUMULATED ARRAYS
      call resetClassYr(nltest)

    end if !> IDAY==365/366 .AND. NDAY

  end subroutine class_annual_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_convertUnitsCTEM
  !! @{
  !> Do some unit conversions for CTEM (biogeochemical processes) so they are ready to be written out
  subroutine convertUnitsCTEM (nltest, nmtest)

    use classicParams, only : icc, ignd, iccp1, wtCH4, convertkgC, iccp2
    use ctemStateVars, only : vrot

    implicit none

    integer, intent(in) :: nltest   ! number of grid cells (offline = 1)
    integer, intent(in) :: nmtest   ! number of tiles per cell

    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:) :: nppmossrow
    real, pointer, dimension(:,:) :: armossrow
    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    ! real, pointer, dimension(:,:,:,:) :: litresvegrow
    ! real, pointer, dimension(:,:,:,:) :: soilcresvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:) :: npprow
    real, pointer, dimension(:,:) :: neprow
    real, pointer, dimension(:,:) :: nbprow
    real, pointer, dimension(:,:) :: gpprow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    real, pointer, dimension(:,:) :: hetroresrow
    real, pointer, dimension(:,:) :: autoresrow
    real, pointer, dimension(:,:) :: litresrow
    real, pointer, dimension(:,:) :: socresrow
    real, pointer, dimension(:,:) :: ch4WetSpecrow
    real, pointer, dimension(:,:) :: ch4WetDynrow
    real, pointer, dimension(:,:) :: ch4soillsrow
    real, pointer, dimension(:,:,:) :: emit_co2row

    integer :: i, m, j, k

    fcancmxrow        => vrot%fcancmx
    gppvegrow         => vrot%gppveg
    nepvegrow         => vrot%nepveg
    nbpvegrow         => vrot%nbpveg
    nppvegrow         => vrot%nppveg
    hetroresvegrow    => vrot%hetroresveg
    autoresvegrow     => vrot%autoresveg
    rmrvegrow         => vrot%rmrveg
    rgvegrow          => vrot%rgveg
    litresvegrow      => vrot%litresveg
    soilcresvegrow    => vrot%soilcresveg
    npprow            => vrot%npp
    neprow            => vrot%nep
    nbprow            => vrot%nbp
    gpprow            => vrot%gpp
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    hetroresrow       => vrot%hetrores
    autoresrow        => vrot%autores
    litresrow         => vrot%litres
    socresrow         => vrot%socres
    ch4WetSpecrow        => vrot%ch4WetSpec
    ch4WetDynrow        => vrot%ch4WetDyn
    ch4soillsrow      => vrot%ch4_soills
    nppmossrow         => vrot%nppmoss
    armossrow          => vrot%armoss
    emit_co2row       => vrot%emit_co2

    !> Some unit conversions:
    !!
    !! We want to go from umol CO2/m2/s to kg C/m2/s so:
    !! umolCO2/m2/s * mol/10^6umol * mol C/ molCO2 * 12.01 g C / mol C * 1 kg/ 1000g = kgC/m2/s
    !! umolCO2/m2/s * 1.201E-8 = kgC/m2/s
    !! convertkgC = 1.201E-8

    do i = 1,nltest
      do m = 1,nmtest

        nppmossrow(i,m) = nppmossrow(i,m) * convertkgC
        armossrow(i,m) = armossrow(i,m) * convertkgC

        do j = 1,icc
          if (fcancmxrow(i,m,j) > 0.0) then

            gppvegrow(i,m,j) = gppvegrow(i,m,j) * convertkgC
            nppvegrow(i,m,j) = nppvegrow(i,m,j) * convertkgC
            nepvegrow(i,m,j) = nepvegrow(i,m,j) * convertkgC
            nbpvegrow(i,m,j) = nbpvegrow(i,m,j) * convertkgC
            hetroresvegrow(i,m,j) = hetroresvegrow(i,m,j) * convertkgC
            autoresvegrow(i,m,j) = autoresvegrow(i,m,j) * convertkgC
            rmrvegrow(i,m,j) = rmrvegrow(i,m,j) * convertkgC
            rgvegrow(i,m,j) = rgvegrow(i,m,j) * convertkgC
            
            ! COMBAK PERLAY
            litresvegrow(i,m,j) = litresvegrow(i,m,j) * convertkgC
            soilcresvegrow(i,m,j) = soilcresvegrow(i,m,j) * convertkgC
            ! do k = 1,ignd
            !   litresvegrow(i,m,j,k)=litresvegrow(i,m,j,k)*convertkgC
            !   soilcresvegrow(i,m,j,k)=soilcresvegrow(i,m,j,k)*convertkgC
            ! end do
            ! COMBAK PERLAY

            ! emit_co2,like all fire gas fluxes is in kg {species} / m2 / s. So we need
            ! to convert from kg CO2/m2/s to kg C/m2/s. Since 1 g C = 0.083 mole CO2 = 3.664 g CO2
            ! kg CO2/m2/s * 1 g C/ 3.664 g CO2 = kg C /m2/s
            emit_co2row(i,m,j) = emit_co2row(i,m,j) / 3.664

          end if
        end do ! loop 30 ! icc

        !> Now for the bare fraction of the grid cell and the LUC products pool.
        hetroresvegrow(i,m,iccp1) = hetroresvegrow(i,m,iccp1) * convertkgC
        nepvegrow(i,m,iccp1) = nepvegrow(i,m,iccp1) * convertkgC
        nbpvegrow(i,m,iccp1) = nbpvegrow(i,m,iccp1) * convertkgC
        ! COMBAK PERLAY
        litresvegrow(i,m,iccp1:iccp2) = litresvegrow(i,m,iccp1:iccp2) * convertkgC
        soilcresvegrow(i,m,iccp1:iccp2) = soilcresvegrow(i,m,iccp1:iccp2) * convertkgC
        ! do k = 1,ignd
        !   litresvegrow(i,m,iccp1:iccp2,k)=litresvegrow(i,m,iccp1:iccp2,k)*convertkgC
        !   soilcresvegrow(i,m,iccp1:iccp2,k)=soilcresvegrow(i,m,iccp1:iccp2,k)*convertkgC
        ! end do
        ! COMBAK PERLAY
        npprow(i,m)     = npprow(i,m) * convertkgC
        gpprow(i,m)     = gpprow(i,m) * convertkgC
        neprow(i,m)     = neprow(i,m) * convertkgC
        nbprow(i,m)     = nbprow(i,m) * convertkgC
        lucemcomrow(i,m) = lucemcomrow(i,m) * convertkgC
        lucltrinrow(i,m) = lucltrinrow(i,m) * convertkgC
        lucsocinrow(i,m) = lucsocinrow(i,m) * convertkgC
        hetroresrow(i,m) = hetroresrow(i,m) * convertkgC
        autoresrow(i,m) = autoresrow(i,m) * convertkgC
        litresrow(i,m)  = litresrow(i,m) * convertkgC
        socresrow(i,m)  = socresrow(i,m) * convertkgC

        ch4WetSpecrow(i,m) = ch4WetSpecrow(i,m) * convertkgC * wtCH4 / 12.01 ! convert from umolch4/m2/s to kg CH4/ m2 /s
        ch4WetDynrow(i,m) = ch4WetDynrow(i,m) * convertkgC * wtCH4 / 12.01 ! convert from umolch4/m2/s to kg CH4/ m2 /s
        ch4soillsrow(i,m) = ch4soillsrow(i,m) * convertkgC * wtCH4 / 12.01 ! convert from umolch4/m2/s to kg CH4/ m2 /s

      end do ! loop 20
    end do ! loop 10

  end subroutine convertUnitsCTEM
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_ctem_daily_aw
  !> @{
  !> Accumulate and write the daily biogeochemical outputs
  !! @author J. Melton

  subroutine ctem_daily_aw (lonLocalIndex, latLocalIndex, nltest, nmtest, iday, ncount, nday, realyr, grclarea, ipeatlandrow)

    ! J. Melton Feb 2016.

    use classStateVars, only : class_rot
    use ctemStateVars,  only : vrot, c_switch !, resetdaily, ctem_grd ctem_tile,
    use classicParams,  only : icc, ignd, nmos, iccp1, wtCH4, seed, convertkgC
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    integer, intent(in) :: ncount
    integer, intent(in) :: nday
    integer, intent(in) :: realyr
    real, intent(in), dimension(:) :: grclarea ! use pointer FLAG
    integer, intent(in), dimension(:,:) :: ipeatlandrow ! use pointer FLAG

    ! pointers

    logical, pointer :: dofire
    logical, pointer :: lnduseon
    logical, pointer :: PFTCompetition
    logical, pointer :: doperpftoutput
    logical, pointer :: dopertileoutput
    logical, pointer :: transientOBSWETF
    integer, pointer :: fixedYearOBSWETF

    real, pointer, dimension(:,:) :: FAREROT !< Fractional coverage of mosaic tile on modelled area

    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    ! real, pointer, dimension(:,:,:,:) :: litresvegrow
    ! real, pointer, dimension(:,:,:,:) :: soilcresvegrow
    ! real, pointer, dimension(:,:,:,:) :: litrmassrow
    ! real, pointer, dimension(:,:,:,:) :: soilcmasrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: rmlvegaccrow
    real, pointer, dimension(:,:,:) :: rmsvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    real, pointer, dimension(:,:,:) :: ailcgrow
    real, pointer, dimension(:,:,:) :: emit_co2row
    real, pointer, dimension(:,:,:) :: emit_corow
    real, pointer, dimension(:,:,:) :: emit_ch4row
    real, pointer, dimension(:,:,:) :: emit_nmhcrow
    real, pointer, dimension(:,:,:) :: emit_h2row
    real, pointer, dimension(:,:,:) :: emit_noxrow
    real, pointer, dimension(:,:,:) :: emit_n2orow
    real, pointer, dimension(:,:,:) :: emit_pm25row
    real, pointer, dimension(:,:,:) :: emit_tpmrow
    real, pointer, dimension(:,:,:) :: emit_tcrow
    real, pointer, dimension(:,:,:) :: emit_ocrow
    real, pointer, dimension(:,:,:) :: emit_bcrow
    real, pointer, dimension(:,:) :: burnfracrow
    real, pointer, dimension(:,:,:) :: burnvegfrow
    real, pointer, dimension(:,:,:) :: smfuncvegrow
    real, pointer, dimension(:,:,:) :: btermrow
    real, pointer, dimension(:,:) :: ltermrow
    real, pointer, dimension(:,:,:) :: mtermrow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    real, pointer, dimension(:,:) :: ch4WetSpecrow
    real, pointer, dimension(:,:) :: wetfdynrow
    real, pointer, dimension(:,:) :: ch4WetDynrow
    real, pointer, dimension(:,:) :: ch4soillsrow
    real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:,:) :: gleafmasrow
    real, pointer, dimension(:,:,:) :: bleafmasrow
    real, pointer, dimension(:,:) :: gavglairow
    real, pointer, dimension(:,:,:) :: slairow
    real, pointer, dimension(:,:,:) :: ailcbrow
    real, pointer, dimension(:,:,:) :: flhrlossrow
    real, pointer, dimension(:,:) :: dstcemls3row
    integer, pointer, dimension(:,:,:) :: lfstatusrow
    real, pointer, dimension(:,:) :: vgbiomasrow
    real, pointer, dimension(:,:) :: gavgltmsrow
    real, pointer, dimension(:,:) :: gavgscmsrow

    real, pointer, dimension(:,:) :: npprow
    real, pointer, dimension(:,:) :: neprow
    real, pointer, dimension(:,:) :: nbprow
    real, pointer, dimension(:,:) :: gpprow
    real, pointer, dimension(:,:) :: hetroresrow
    real, pointer, dimension(:,:) :: autoresrow
    real, pointer, dimension(:,:) :: soilcresprow
    real, pointer, dimension(:,:) :: rgrow
    real, pointer, dimension(:,:) :: litresrow
    real, pointer, dimension(:,:) :: socresrow

    real, pointer, dimension(:,:) :: nppmossrow
    real, pointer, dimension(:,:) :: armossrow
    real, pointer, dimension(:,:,:) :: bmasvegrow
    real, pointer, dimension(:,:,:) :: cmasvegcrow
    real, pointer, dimension(:,:,:) :: veghghtrow
    real, pointer, dimension(:,:,:) :: rootdpthrow
    real, pointer, dimension(:,:) :: rmlrow
    real, pointer, dimension(:,:) :: rmsrow
    real, pointer, dimension(:,:,:) :: tltrleafrow
    real, pointer, dimension(:,:,:) :: tltrstemrow
    real, pointer, dimension(:,:,:) :: tltrrootrow
    real, pointer, dimension(:,:,:) :: leaflitrrow
    real, pointer, dimension(:,:,:) :: roottemprow
    real, pointer, dimension(:,:,:) :: afrleafrow
    real, pointer, dimension(:,:,:) :: afrstemrow
    real, pointer, dimension(:,:,:) :: afrrootrow
    real, pointer, dimension(:,:,:) :: wtstatusrow
    real, pointer, dimension(:,:,:) :: ltstatusrow
    real, pointer, dimension(:,:) :: rmrrow

    real, pointer, dimension(:,:,:,:) :: rmatctemrow
    real, pointer, dimension(:,:) :: dstcemlsrow
    real, pointer, dimension(:,:) :: litrfallrow
    real, pointer, dimension(:,:) :: humiftrsrow

    ! local vars
    real, allocatable, dimension(:) :: gpp_g      !<
    real, allocatable, dimension(:) :: npp_g      !<
    real, allocatable, dimension(:) :: nbp_g      !<
    real, allocatable, dimension(:) :: socres_g   !<
    real, allocatable, dimension(:) :: autores_g  !<
    real, allocatable, dimension(:) :: litres_g   !<
    real, allocatable, dimension(:) :: dstcemls3_g !<
    real, allocatable, dimension(:) :: litrfall_g !<
    real, allocatable, dimension(:) :: rml_g      !<
    real, allocatable, dimension(:) :: rms_g      !<
    real, allocatable, dimension(:) :: rg_g       !<
    real, allocatable, dimension(:) :: leaflitr_g !<
    real, allocatable, dimension(:) :: tltrstem_g !<
    real, allocatable, dimension(:) :: tltrroot_g !<
    real, allocatable, dimension(:) :: nep_g      !<
    real, allocatable, dimension(:) :: hetrores_g !<
    real, allocatable, dimension(:) :: dstcemls_g !<
    real, allocatable, dimension(:) :: humiftrs_g !<
    real, allocatable, dimension(:) :: rmr_g      !<
    real, allocatable, dimension(:) :: tltrleaf_g !<
    real, allocatable, dimension(:) :: gavgltms_g !<
    real, allocatable, dimension(:) :: vgbiomas_g !<
    real, allocatable, dimension(:) :: gavglai_g  !<
    real, allocatable, dimension(:) :: gavgscms_g !<
    real, allocatable, dimension(:) :: gleafmas_g !<
    real, allocatable, dimension(:) :: bleafmas_g !<
    real, allocatable, dimension(:) :: stemmass_g !<
    real, allocatable, dimension(:) :: rootmass_g !<
    real, allocatable, dimension(:) :: litrmass_g !<
    real, allocatable, dimension(:) :: soilcmas_g !<
    real, allocatable, dimension(:) :: slai_g     !<
    real, allocatable, dimension(:) :: ailcg_g    !<
    real, allocatable, dimension(:) :: ailcb_g    !<
    real, allocatable, dimension(:) :: veghght_g  !<
    real, allocatable, dimension(:) :: rootdpth_g !<
    real, allocatable, dimension(:) :: roottemp_g !<
    real, allocatable, dimension(:) :: totcmass_g !<
    real, allocatable, dimension(:) :: tcanoacc_out_g !<
    real, allocatable, dimension(:) :: burnfrac_g !<
    real, allocatable, dimension(:) :: smfuncveg_g !<
    real, allocatable, dimension(:) :: lucemcom_g !<
    real, allocatable, dimension(:) :: lucltrin_g !<
    real, allocatable, dimension(:) :: lucsocin_g !<
    real, allocatable, dimension(:) :: emit_co2_g !<
    real, allocatable, dimension(:) :: emit_ch4_g !<
    real, allocatable, dimension(:) :: ch4WetSpec_g  !<
    real, allocatable, dimension(:) :: wetfdyn_g  !<
    real, allocatable, dimension(:) :: ch4WetDyn_g  !<
    real, allocatable, dimension(:) :: ch4soills_g   !<
    real, allocatable, dimension(:,:) :: afrleaf_g  !<
    real, allocatable, dimension(:,:) :: afrstem_g  !<
    real, allocatable, dimension(:,:) :: afrroot_g  !<
    real, allocatable, dimension(:,:) :: lfstatus_g !<
    real, allocatable, dimension(:,:) :: rmlvegrow_g !<
    real, allocatable, dimension(:,:) :: anvegrow_g !<
    real, allocatable, dimension(:,:) :: rmatctem_g !<

    real, allocatable, dimension(:,:) :: leaflitr_t !<
    real, allocatable, dimension(:,:) :: tltrleaf_t !<
    real, allocatable, dimension(:,:) :: tltrstem_t !<
    real, allocatable, dimension(:,:) :: tltrroot_t !<
    real, allocatable, dimension(:,:) :: ailcg_t    !<
    real, allocatable, dimension(:,:) :: ailcb_t    !<
    real, allocatable, dimension(:,:,:) :: rmatctem_t !< nlat,nmos,ignd
    real, allocatable, dimension(:,:) :: veghght_t  !<
    real, allocatable, dimension(:,:) :: rootdpth_t !<
    real, allocatable, dimension(:,:) :: roottemp_t !<
    real, allocatable, dimension(:,:) :: slai_t     !<
    real, allocatable, dimension(:,:) :: afrroot_t  !<
    real, allocatable, dimension(:,:) :: afrleaf_t  !<
    real, allocatable, dimension(:,:) :: afrstem_t  !<
    real, allocatable, dimension(:,:) :: laimaxg_t  !<
    real, allocatable, dimension(:,:) :: stemmass_t !<
    real, allocatable, dimension(:,:) :: rootmass_t !<
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: litrmass_t !<
    real, allocatable, dimension(:,:) :: soilcmas_t !<
    ! real, allocatable, dimension(:,:,:) :: litrmass_t !<
    ! real, allocatable, dimension(:,:,:) :: soilcmas_t !<
    ! COMBAK PERLAY
    real, allocatable, dimension(:,:) :: gleafmas_t !<
    real, allocatable, dimension(:,:) :: bleafmas_t !<
    real, allocatable, dimension(:,:) :: emit_co2_t !<
    real, allocatable, dimension(:,:) :: emit_ch4_t !<
    real, allocatable, dimension(:,:) :: smfuncveg_t !<

    integer :: i, m, j, nt, k
    real :: barefrac
    real :: sumfare
    real, dimension(1) :: timeStamp

    ! point pointers

    dofire                => c_switch%dofire
    lnduseon              => c_switch%lnduseon
    PFTCompetition        => c_switch%PFTCompetition
    doperpftoutput        => c_switch%doperpftoutput
    dopertileoutput       => c_switch%dopertileoutput
    transientOBSWETF      => c_switch%transientOBSWETF
    fixedYearOBSWETF      => c_switch%fixedYearOBSWETF

    FAREROT => class_rot%FAREROT

    fcancmxrow        => vrot%fcancmx
    gppvegrow         => vrot%gppveg
    nepvegrow         => vrot%nepveg
    nbpvegrow         => vrot%nbpveg
    nppvegrow         => vrot%nppveg
    hetroresvegrow    => vrot%hetroresveg
    autoresvegrow     => vrot%autoresveg
    litresvegrow      => vrot%litresveg
    soilcresvegrow    => vrot%soilcresveg
    rmlvegaccrow      => vrot%rmlvegacc
    rmsvegrow         => vrot%rmsveg
    rmrvegrow         => vrot%rmrveg
    rgvegrow          => vrot%rgveg
    ailcgrow          => vrot%ailcg
    emit_co2row       => vrot%emit_co2
    emit_corow        => vrot%emit_co
    emit_ch4row       => vrot%emit_ch4
    emit_nmhcrow      => vrot%emit_nmhc
    emit_h2row        => vrot%emit_h2
    emit_noxrow       => vrot%emit_nox
    emit_n2orow       => vrot%emit_n2o
    emit_pm25row      => vrot%emit_pm25
    emit_tpmrow       => vrot%emit_tpm
    emit_tcrow        => vrot%emit_tc
    emit_ocrow        => vrot%emit_oc
    emit_bcrow        => vrot%emit_bc
    burnfracrow       => vrot%burnfrac
    burnvegfrow       => vrot%burnvegf
    smfuncvegrow      => vrot%smfuncveg
    btermrow          => vrot%bterm
    ltermrow          => vrot%lterm
    mtermrow          => vrot%mterm
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    ch4WetSpecrow        => vrot%ch4WetSpec
    wetfdynrow        => vrot%wetfdyn
    ch4WetDynrow        => vrot%ch4WetDyn
    ch4soillsrow      => vrot%ch4_soills
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    vgbiomas_vegrow   => vrot%vgbiomas_veg
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    flhrlossrow       => vrot%flhrloss
    dstcemls3row      => vrot%dstcemls3
    lfstatusrow       => vrot%lfstatus

    npprow            => vrot%npp
    neprow            => vrot%nep
    nbprow            => vrot%nbp
    gpprow            => vrot%gpp
    hetroresrow       => vrot%hetrores
    autoresrow        => vrot%autores
    soilcresprow      => vrot%soilcresp
    rgrow             => vrot%rg
    litresrow         => vrot%litres
    socresrow         => vrot%socres
    vgbiomasrow       => vrot%vgbiomas
    gavgltmsrow       => vrot%gavgltms
    gavgscmsrow       => vrot%gavgscms

    nppmossrow         => vrot%nppmoss
    armossrow          => vrot%armoss

    bmasvegrow        => vrot%bmasveg
    cmasvegcrow       => vrot%cmasvegc
    veghghtrow        => vrot%veghght
    rootdpthrow       => vrot%rootdpth
    rmlrow            => vrot%rml
    rmsrow            => vrot%rms
    tltrleafrow       => vrot%tltrleaf
    tltrstemrow       => vrot%tltrstem
    tltrrootrow       => vrot%tltrroot
    leaflitrrow       => vrot%leaflitr
    roottemprow       => vrot%roottemp
    afrleafrow        => vrot%afrleaf
    afrstemrow        => vrot%afrstem
    afrrootrow        => vrot%afrroot
    wtstatusrow       => vrot%wtstatus
    ltstatusrow       => vrot%ltstatus
    rmrrow            => vrot%rmr
    gleafmasrow       => vrot%gleafmas
    bleafmasrow       => vrot%bleafmas
    gavglairow        => vrot%gavglai
    slairow           => vrot%slai
    ailcbrow          => vrot%ailcb
    flhrlossrow       => vrot%flhrloss
    rmatctemrow       => vrot%rmatctem
    dstcemlsrow       => vrot%dstcemls
    litrfallrow       => vrot%litrfall
    humiftrsrow       => vrot%humiftrs

    ! Allocate the grid average and tile average structures
    allocate(gpp_g (nltest),npp_g (nltest),nbp_g (nltest),socres_g (nltest), &
             autores_g (nltest),litres_g (nltest),dstcemls3_g (nltest),litrfall_g (nltest), &
             rml_g (nltest),rms_g (nltest),rg_g (nltest),leaflitr_g (nltest),tltrstem_g (nltest), &
             tltrroot_g (nltest),nep_g (nltest),hetrores_g (nltest),dstcemls_g (nltest), &
             humiftrs_g (nltest),rmr_g (nltest),tltrleaf_g (nltest),gavgltms_g (nltest), &
             vgbiomas_g (nltest),gavglai_g (nltest),gavgscms_g (nltest),gleafmas_g (nltest), &
             bleafmas_g (nltest),stemmass_g (nltest),rootmass_g (nltest),litrmass_g (nltest), &
             soilcmas_g (nltest),slai_g (nltest),ailcg_g (nltest),ailcb_g (nltest),veghght_g (nltest), &
             rootdpth_g (nltest),roottemp_g (nltest),totcmass_g (nltest),tcanoacc_out_g (nltest), &
             burnfrac_g (nltest),smfuncveg_g (nltest),lucemcom_g (nltest),lucltrin_g (nltest),lucsocin_g (nltest), &
             emit_co2_g (nltest),emit_ch4_g (nltest),ch4WetSpec_g (nltest),wetfdyn_g (nltest),ch4WetDyn_g (nltest), &
             ch4soills_g (nltest),afrleaf_g (nltest,icc),afrstem_g (nltest,icc), &
             afrroot_g (nltest,icc),lfstatus_g (nltest,icc),rmlvegrow_g (nltest,icc),anvegrow_g(nltest,icc), &
             rmatctem_g (nltest,ignd) )

    allocate(leaflitr_t (nltest,nmtest),tltrleaf_t (nltest,nmtest),tltrstem_t (nltest,nmtest),tltrroot_t (nltest,nmtest), &
             ailcg_t (nltest,nmtest),ailcb_t (nltest,nmtest),rmatctem_t (nltest,nmtest,ignd),veghght_t (nltest,nmtest), &
             rootdpth_t (nltest,nmtest),roottemp_t (nltest,nmtest),slai_t (nltest,nmtest),afrroot_t (nltest,nmtest), &
             afrleaf_t (nltest,nmtest),afrstem_t (nltest,nmtest),laimaxg_t (nltest,nmtest),stemmass_t (nltest,nmtest), &
             rootmass_t (nltest,nmtest),gleafmas_t (nltest,nmtest),bleafmas_t(nltest,nmtest), &
             emit_co2_t (nltest,nmtest),emit_ch4_t (nltest,nmtest),smfuncveg_t (nltest,nmtest), &
             ! COMBAK PERLAY
             soilcmas_t (nltest,nmtest),litrmass_t (nltest,nmtest) )
    ! soilcmas_t (nltest,nmtest,ignd),litrmass_t (nltest,nmtest,ignd) )
    ! COMBAK PERLAY


    ! First set the local variables to 0.
    gpp_g(:) = 0.0 ; npp_g(:) = 0.0 ; nep_g(:) = 0.0 ; nbp_g(:) = 0.0 ; autores_g(:) = 0.0
    hetrores_g(:) = 0.0 ; litres_g(:) = 0.0 ; socres_g(:) = 0.0 ; dstcemls_g(:) = 0.0
    dstcemls3_g(:) = 0.0 ; litrfall_g(:) = 0.0 ; humiftrs_g(:) = 0.0 ; rml_g(:) = 0.0
    rms_g(:) = 0.0 ; rmr_g(:) = 0.0 ; rg_g(:) = 0.0 ; vgbiomas_g(:) = 0.0 ; totcmass_g(:) = 0.0
    gavglai_g(:) = 0.0 ; gavgltms_g(:) = 0.0 ; gavgscms_g(:) = 0.0 ; ailcg_g(:) = 0.0
    ailcb_g(:) = 0.0 ; tcanoacc_out_g(:) = 0.0 ; burnfrac_g(:) = 0.0 ; smfuncveg_g(:) = 0.0
    lucemcom_g(:) = 0.0 ; lucltrin_g(:) = 0.0 ; lucsocin_g(:) = 0.0 ; emit_co2_g(:) = 0.0
    leaflitr_g(:) = 0.0 ; tltrleaf_g(:) = 0.0 ; tltrstem_g(:) = 0.0 ; tltrroot_g(:) = 0.0
    gleafmas_g(:) = 0.0 ; bleafmas_g(:) = 0.0 ; stemmass_g(:) = 0.0 ; rootmass_g(:) = 0.0
    litrmass_g(:) = 0.0 ; soilcmas_g(:) = 0.0 ; veghght_g(:) = 0.0 ; rootdpth_g(:) = 0.0
    roottemp_g(:) = 0.0 ; slai_g(:) = 0.0 ; ch4WetSpec_G(:) = 0.0
    WETFDYN_G(:) = 0.0 ; ch4WetDyn_G(:) = 0.0 ; ch4soills_g(:) = 0.0
    rmatctem_g(:,:) = 0.0
    afrleaf_g(:,:) = 0.0 ; afrstem_g(:,:) = 0.0 ; afrroot_g(:,:) = 0.0
    leaflitr_t(:,:) = 0.0 ; tltrleaf_t(:,:) = 0.0 ; tltrstem_t(:,:) = 0.0
    tltrroot_t(:,:) = 0.0 ; ailcg_t(:,:) = 0.0 ; ailcb_t(:,:) = 0.0
    afrleaf_t(:,:) = 0.0 ; afrstem_t(:,:) = 0.0 ; afrroot_t(:,:) = 0.0
    veghght_t(:,:) = 0.0 ; rootdpth_t(:,:) = 0.0 ; roottemp_t(:,:) = 0.0
    slai_t(:,:) = 0.0 ; gleafmas_t(:,:) = 0.0 ; bleafmas_t(:,:) = 0.0
    stemmass_t(:,:) = 0.0 ; rootmass_t(:,:) = 0.0
    emit_co2_t(:,:) = 0.0; rmatctem_t(:,:,:) = 0.0
    ! COMBAK PERLAY
    soilcmas_t(:,:) = 0.0 ; litrmass_t(:,:) = 0.0
    ! soilcmas_t(:,:,:) = 0.0 ; litrmass_t(:,:,:) = 0.0
    ! COMBAK PERLAY

    !> Aggregate to the tile avg vars:
    do i = 1,nltest
      do m = 1,nmtest
        barefrac = 1.0
        do j = 1,icc
          barefrac = barefrac - fcancmxrow(i,m,j)
          leaflitr_t(i,m) = leaflitr_t(i,m) + leaflitrrow(i,m,j) * fcancmxrow(i,m,j)
          tltrleaf_t(i,m) = tltrleaf_t(i,m) + tltrleafrow(i,m,j) * fcancmxrow(i,m,j)
          tltrstem_t(i,m) = tltrstem_t(i,m) + tltrstemrow(i,m,j) * fcancmxrow(i,m,j)
          tltrroot_t(i,m) = tltrroot_t(i,m) + tltrrootrow(i,m,j) * fcancmxrow(i,m,j)
          veghght_t(i,m) = veghght_t(i,m) + veghghtrow(i,m,j) * fcancmxrow(i,m,j)
          rootdpth_t(i,m) = rootdpth_t(i,m) + rootdpthrow(i,m,j) * fcancmxrow(i,m,j)
          roottemp_t(i,m) = roottemp_t(i,m) + roottemprow(i,m,j) * fcancmxrow(i,m,j)
          slai_t(i,m) = slai_t(i,m) + slairow(i,m,j) * fcancmxrow(i,m,j)
          afrleaf_t(i,m) = afrleaf_t(i,m) + afrleafrow(i,m,j) * fcancmxrow(i,m,j)
          afrstem_t(i,m) = afrstem_t(i,m) + afrstemrow(i,m,j) * fcancmxrow(i,m,j)
          afrroot_t(i,m) = afrroot_t(i,m) + afrrootrow(i,m,j) * fcancmxrow(i,m,j)
          ailcg_t(i,m) = ailcg_t(i,m) + ailcgrow(i,m,j) * fcancmxrow(i,m,j)
          ailcb_t(i,m) = ailcb_t(i,m) + ailcbrow(i,m,j) * fcancmxrow(i,m,j)
          gleafmas_t(i,m) = gleafmas_t(i,m) + gleafmasrow(i,m,j) * fcancmxrow(i,m,j)
          bleafmas_t(i,m) = bleafmas_t(i,m) + bleafmasrow(i,m,j) * fcancmxrow(i,m,j)
          stemmass_t(i,m) = stemmass_t(i,m) + stemmassrow(i,m,j) * fcancmxrow(i,m,j)
          rootmass_t(i,m) = rootmass_t(i,m) + rootmassrow(i,m,j) * fcancmxrow(i,m,j)
          emit_co2_t(i,m) = emit_co2_t(i,m) + emit_co2row(i,m,j) * fcancmxrow(i,m,j)
          emit_ch4_t(i,m) = emit_ch4_t(i,m) + emit_ch4row(i,m,j) * fcancmxrow(i,m,j)
          smfuncveg_t(i,m) = smfuncveg_t(i,m) + smfuncvegrow(i,m,j) * fcancmxrow(i,m,j)

          do k = 1,ignd
            rmatctem_t(i,m,k) = rmatctem_t(i,m,k) + rmatctemrow(i,m,j,k) * fcancmxrow(i,m,j)
            ! COMBAK PERLAY
            litrmass_t(i,m) = litrmass_t(i,m) + litrmassrow(i,m,j) * fcancmxrow(i,m,j)
            soilcmas_t(i,m) = soilcmas_t(i,m) + soilcmasrow(i,m,j) * fcancmxrow(i,m,j)
            ! litrmass_t(i,m,k) = litrmass_t(i,m,k) + litrmassrow(i,m,j,k)*fcancmxrow(i,m,j)
            ! soilcmas_t(i,m,k) = soilcmas_t(i,m,k) + soilcmasrow(i,m,j,k)*fcancmxrow(i,m,j)
            ! COMBAK PERLAY
          end do
        end do ! icc

        !> Do the bare ground also:
        ! COMBAK PERLAY
        litrmass_t(i,m) = litrmass_t(i,m) + litrmassrow(i,m,iccp1) * barefrac
        soilcmas_t(i,m) = soilcmas_t(i,m) + soilcmasrow(i,m,iccp1) * barefrac
        ! do k = 1,ignd
        !   litrmass_t(i,m,k) = litrmass_t(i,m,k) + litrmassrow(i,m,iccp1,k)*barefrac
        !   soilcmas_t(i,m,k) = soilcmas_t(i,m,k) + soilcmasrow(i,m,iccp1,k)*barefrac
        ! end do
        ! COMBAK PERLAY
        !> Calculation of grid averaged variables

        gpp_g(i) = gpp_g(i) + gpprow(i,m) * FAREROT(i,m)
        npp_g(i) = npp_g(i) + npprow(i,m) * FAREROT(i,m)
        nep_g(i) = nep_g(i) + neprow(i,m) * FAREROT(i,m)
        ! NOTE: This NBP will include LUC product pool contributions
        ! since it is the nbprow varialbe. nbpveg variable does not
        ! include LUC contributions since they exist at the tile level
        ! not at the PFT level.
        nbp_g(i) = nbp_g(i) + nbprow(i,m) * FAREROT(i,m)
        autores_g(i) = autores_g(i) + autoresrow(i,m) * FAREROT(i,m)
        hetrores_g(i) = hetrores_g(i) + hetroresrow(i,m) * FAREROT(i,m)
        litres_g(i) = litres_g(i) + litresrow(i,m) * FAREROT(i,m)
        socres_g(i) = socres_g(i) + socresrow(i,m) * FAREROT(i,m)
        dstcemls_g(i) = dstcemls_g(i) + dstcemlsrow(i,m) * FAREROT(i,m)
        dstcemls3_g(i) = dstcemls3_g(i) + dstcemls3row(i,m) * FAREROT(i,m)
        litrfall_g(i) = litrfall_g(i) + litrfallrow(i,m) * FAREROT(i,m)
        humiftrs_g(i) = humiftrs_g(i) + humiftrsrow(i,m) * FAREROT(i,m)
        rml_g(i) = rml_g(i) + rmlrow(i,m) * FAREROT(i,m)
        rms_g(i) = rms_g(i) + rmsrow(i,m) * FAREROT(i,m)
        rmr_g(i) = rmr_g(i) + rmrrow(i,m) * FAREROT(i,m)
        rg_g(i) = rg_g(i) + rgrow(i,m) * FAREROT(i,m)
        leaflitr_g(i) = leaflitr_g(i) + leaflitr_t(i,m) * FAREROT(i,m)
        tltrleaf_g(i) = tltrleaf_g(i) + tltrleaf_t(i,m) * FAREROT(i,m)
        tltrstem_g(i) = tltrstem_g(i) + tltrstem_t(i,m) * FAREROT(i,m)
        tltrroot_g(i) = tltrroot_g(i) + tltrroot_t(i,m) * FAREROT(i,m)
        slai_g(i) = slai_g(i) + slai_t(i,m) * FAREROT(i,m)
        ailcg_g(i) = ailcg_g(i) + ailcg_t(i,m) * FAREROT(i,m)
        ailcb_g(i) = ailcb_g(i) + ailcb_t(i,m) * FAREROT(i,m)
        vgbiomas_g(i) = vgbiomas_g(i) + vgbiomasrow(i,m) * FAREROT(i,m)
        veghght_g(i) = veghght_g(i) + veghght_t(i,m) * FAREROT(i,m)
        gavglai_g(i) = gavglai_g(i) + gavglairow(i,m) * FAREROT(i,m)
        gavgltms_g(i) = gavgltms_g(i) + gavgltmsrow(i,m) * FAREROT(i,m)
        gavgscms_g(i) = gavgscms_g(i) + gavgscmsrow(i,m) * FAREROT(i,m)
        totcmass_g(i) = vgbiomas_g(i) + gavgltms_g(i) + gavgscms_g(i)
        gleafmas_g(i) = gleafmas_g(i) + gleafmas_t(i,m) * FAREROT(i,m)
        bleafmas_g(i) = bleafmas_g(i) + bleafmas_t(i,m) * FAREROT(i,m)
        stemmass_g(i) = stemmass_g(i) + stemmass_t(i,m) * FAREROT(i,m)
        rootmass_g(i) = rootmass_g(i) + rootmass_t(i,m) * FAREROT(i,m)
        rootdpth_g(i) = rootdpth_g(i) + rootdpth_t(i,m) * FAREROT(i,m)
        roottemp_g(i) = roottemp_g(i) + roottemp_t(i,m) * FAREROT(i,m)
        burnfrac_g(i) = burnfrac_g(i) + burnfracrow(i,m) * FAREROT(i,m)
        smfuncveg_g(i) = smfuncveg_g(i) + smfuncveg_t(i,m) * FAREROT(i,m)
        lucemcom_g(i) = lucemcom_g(i) + lucemcomrow(i,m) * FAREROT(i,m)
        lucltrin_g(i) = lucltrin_g(i) + lucltrinrow(i,m) * FAREROT(i,m)
        lucsocin_g(i) = lucsocin_g(i) + lucsocinrow(i,m) * FAREROT(i,m)
        ch4WetSpec_g(i) = ch4WetSpec_g(i) + ch4WetSpecrow(i,m) * farerot(i,m)
        wetfdyn_g(i) = wetfdyn_g(i) + wetfdynrow(i,m) * farerot(i,m)
        ch4WetDyn_g(i) = ch4WetDyn_g(i) + ch4WetDynrow(i,m) * farerot(i,m)
        ch4soills_g(i) = ch4soills_g(i) + ch4soillsrow(i,m) * farerot(i,m)
        emit_co2_g(i) = emit_co2_g(i) + emit_co2_t(i,m) * FAREROT(i,m)
        emit_ch4_g(i) = emit_ch4_g(i) + emit_ch4_t(i,m) * FAREROT(i,m)
        ! nppmoss_g(i)  = nppmoss_g(i) +nppmossrow(i,m)*FAREROT(i,m)
        ! armoss_g(i)   = armoss_g(i) + armossrow(i,m)*FAREROT(i,m)

        do k = 1,ignd
          rmatctem_g(i,k) = rmatctem_g(i,k) + rmatctem_t(i,m,k) * FAREROT(i,m)
          ! FLAG not putting as per layer yet since it is not presently written out.
          ! COMBAK PERLAY
          litrmass_g(i) = litrmass_g(i) + litrmass_t(i,m) * FAREROT(i,m)
          soilcmas_g(i) = soilcmas_g(i) + soilcmas_t(i,m) * FAREROT(i,m)
          ! litrmass_g(i) = litrmass_g(i) + litrmass_t(i,m,k) * FAREROT(i,m)
          ! soilcmas_g(i) = soilcmas_g(i) + soilcmas_t(i,m,k) * FAREROT(i,m)
          ! COMBAK PERLAY
        end do
      end do ! loop 70 ! nmtest
    end do ! loop 60 ! nltest

    i = 1 ! offline nltest is always 1.

    ! Prepare the timestamp for this timestep.
    timeStamp = consecDays + real(iday) - 1.

    !> Write grid average values

    call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_d_g' ,timeStamp,'gpp', [gpp_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_d_g' ,timeStamp,'npp', [npp_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_d_g' ,timeStamp,'nep', [nep_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_d_g' ,timeStamp,'nbp', [nbp_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_d_g' ,timeStamp,'ra', [autores_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_d_g' ,timeStamp,'rh', [hetrores_g(i)])
    if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_d_g' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_g(i)])
    end if
    call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_d_g' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_d_g' ,timeStamp,'soilCH4cons',[ch4soills_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_d_g' ,timeStamp,'fFireCH4',[emit_ch4_g(i)])
    call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_d_g' ,timeStamp,'wetlandFrac',[wetfdyn_g(i)])
    ! litres_g(i),
    ! socres_g(i), &
    !         (dstcemls_g(i)+dstcemls3_g(i)), &
    !         litrfall_g(i),
    ! humiftrs_g(i),' GRDAV'

    if (doperpftoutput) then

      !                     !> File: .CT01D
      !                     write(72,8200)iday,realyr,gppvegrow(i,m,j),nppvegrow(i,m,j), &
      !                     nepvegrow(i,m,j),nbpvegrow(i,m,j),autoresvegrow(i,m,j), &
      !                     hetroresvegrow(i,m,j),litresvegrow(i,m,j),soilcresvegrow(i,m,j), &
      !                     (dstcemlsrow(i,m)+dstcemls3row(i,m)), &   ! FLAG at present dstcemls are only per tile values
      !                     litrfallrow(i,m),humiftrsrow(i,m), & ! same with litrfall and humiftrs.
      !                     ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
      !                 !> File: .CT01D
      !                 write(72,8200)iday,realyr,0.0,0.0, &
      !                 nepvegrow(i,m,iccp1),nbpvegrow(i,m,iccp1),0.0, &
      !                 hetroresvegrow(i,m,iccp1),litresvegrow(i,m,iccp1),soilcresvegrow(i,m,iccp1), &
      !                 (dstcemlsrow(i,m)+dstcemls3row(i,m)), &   ! FLAG at present dstcemls are only per tile values
      !                 litrfallrow(i,m),humiftrsrow(i,m), & ! same with litrfall and humiftrs.
      !                 ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac

    end if

    if (dopertileoutput) then
      !                 !> File: .CT01D
      !                 write(72,8200)iday,realyr,gpprow(i,m),npprow(i,m), &
      !                 neprow(i,m),nbprow(i,m),autoresrow(i,m), &
      !                 hetroresrow(i,m),litresrow(i,m),socresrow(i,m), &
      !                 (dstcemlsrow(i,m)+dstcemls3row(i,m)), &
      !                 litrfallrow(i,m),humiftrsrow(i,m), &
      !                 ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)

    end if


    !     do i=1,nltest
    !         do m=1,nmtest
    !
    !             barefrac = 1.0
    !
    !             !> First the per PFT values to file .CT01D
    !             do j=1,icc
    !
    !                 ! call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo_g' ,timeStamp,'lai', [laimaxg_mo_g(i)])
    !                 if (fcancmxrow(i,m,j) > seed) then
    !
    !                     barefrac = barefrac - fcancmxrow(i,m,j)
    !
    !
    !                     !> File .CT02D
    !                     write(73,8300)iday,realyr,rmlvegaccrow(i,m,j), &
    !                     rmsvegrow(i,m,j),rmrvegrow(i,m,j),rgvegrow(i,m,j), &
    !                     leaflitrrow(i,m,j),tltrleafrow(i,m,j), &
    !                     tltrstemrow(i,m,j),tltrrootrow(i,m,j), &
    !                     ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
    !
    !                     !> File *.CT03D
    !                     write(74,8401)iday,realyr,vgbiomas_vegrow(i,m,j), &
    !                     ailcgrow(i,m,j),gleafmasrow(i,m,j), &
    !                     bleafmasrow(i,m,j),stemmassrow(i,m,j), &
    !                     rootmassrow(i,m,j),litrmassrow(i,m,j),  &
    !                     soilcmasrow(i,m,j), &
    !                     ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
    !
    !                     !> File .CT04D
    !                     write(75,8500)iday,realyr,ailcgrow(i,m,j),  &
    !                     ailcbrow(i,m,j),(rmatctemrow(i,m,j,k),k=1,3), &
    !                     veghghtrow(i,m,j),rootdpthrow(i,m,j), &
    !                     roottemprow(i,m,j),slairow(i,m,j), &
    !                     ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
    !
    !                     ! File .CT05D
    !                     ! write(76,8600)iday,realyr,afrleafrow(i,m,j),  &
    !                     ! afrstemrow(i,m,j),afrrootrow(i,m,j),  &
    !                     !, lfstatusrow(i,m,j), &
    !                     !' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
    !
    !                     !> File *.CT06D
    !                     if (dofire .or. lnduseon) then
    !                         write(77,8800)iday,realyr, &
    !                         emit_co2row(i,m,j),emit_corow(i,m,j),emit_ch4row(i,m,j), &
    !                         emit_nmhcrow(i,m,j),emit_h2row(i,m,j),emit_noxrow(i,m,j), &
    !                         emit_n2orow(i,m,j),emit_pm25row(i,m,j), &
    !                         emit_tpmrow(i,m,j),emit_tcrow(i,m,j),emit_ocrow(i,m,j), &
    !                         emit_bcrow(i,m,j),burnvegfrow(i,m,j)*100., &
    !                         smfuncvegrow(i,m,j),lucemcom_g(i), &  ! FLAG only per grid values for these last ones.
    !                         lucltrin_g(i),lucsocin_g(i), &
    !                         grclarea(i),btermrow(i,m,j),lterm_g(i),mtermrow(i,m,j), &
    !                         ' TILE ',m,' PFT ',j,' FRAC ',fcancmxrow(i,m,j)
    !                     end if
    !
    !                 end if ! fcancmx
    !
    !             end do ! icc
    !
    !             !> Now write out the bare fraction values ( only needed if you have vars that are affected by barefrac values)
    !             if (barefrac > seed) then
    !
    !
    !                 !> File *.CT03D
    !                 write(74,8401)iday,realyr,0.0, &
    !                 0.0,0.0, &
    !                 0.0,0.0, &
    !                 0.0,litrmassrow(i,m,iccp1),  &
    !                 soilcmasrow(i,m,iccp1), &
    !                 ' TILE ',m,' PFT ',iccp1,' FRAC ',barefrac
    !
    !             end if
    !
    !             !> Now write out the tile average values for each tile if the tile number
    !             !> is greater than 1 (nmtest > 1).
    !             if (nmtest > 1) then
    !
    !
    !                 !> File .CT02D
    !                 write(73,8300)iday,realyr,rmlrow(i,m),rmsrow(i,m), &
    !                 rmrrow(i,m),rgrow(i,m),leaflitr_t(i,m),tltrleaf_t(i,m), &
    !                 tltrstem_t(i,m),tltrroot_t(i,m), &
    !                 ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)
    !
    !                 !> File .CT03D
    !                 write(74,8401)iday,realyr,vgbiomasrow(i,m), &
    !                 ailcg_t(i,m),gleafmas_t(i,m), &
    !                 bleafmas_t(i,m),stemmass_t(i,m), &
    !                 rootmass_t(i,m),litrmass_t(i,m), &
    !                 soilcmas_t(i,m), &
    !                 ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)
    !
    !                 !> File .CT04D
    !                 write(75,8500)iday,realyr,ailcg_t(i,m), &
    !                 ailcb_t(i,m),(rmatctem_t(i,m,k),k=1,3), &
    !                 veghght_t(i,m),rootdpth_t(i,m), &
    !                 roottemp_t(i,m),slai_t(i,m), &
    !                 ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)
    !
    !                 ! File .CT05D
    !                 ! write(76,8601)iday,realyr,afrleaf_t(i,m), &
    !                 !    afrstem_t(i,m),afrroot_t(i,m),  &
    !                 !    -999,   & ! lfstatus is kinda meaningless grid avg so set to -999
    !                 !    ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)
    !
    !                 if (dofire .or. lnduseon) then
    !                     write(77,8800)iday,realyr,  &
    !                     emit_co2_t(i,m),emit_co_t(i,m),emit_ch4_t(i,m), &
    !                     emit_nmhc_t(i,m),emit_h2_t(i,m),emit_nox_t(i,m), &
    !                     emit_n2o_t(i,m),emit_pm25_t(i,m),emit_tpm_t(i,m), &
    !                     emit_tc_t(i,m),emit_oc_t(i,m),emit_bc_t(i,m), &
    !                     burnfrac_g(i)*100., smfuncveg_t(i,m),lucemcom_g(i), & ! FLAG only per grid values for these last ones.
    !                     lucltrin_g(i),lucsocin_g(i), &
    !                     grclarea(i),bterm_t(i,m),lterm_g(i),mterm_t(i,m), &
    !                     ' TILE ',m,' OF ',nmtest,' TFRAC ',FAREROT(i,m)
    !                 end if
    !
    !             end if ! nmtest>1
    !
    ! 90              continue ! nmtest
    !
    !         !> Finally do the grid avg values:
    !
    !
    !         !> File .CT02D
    !         write(73,8300)iday,realyr,rml_g(i),rms_g(i), &
    !         rmr_g(i),rg_g(i),leaflitr_g(i),tltrleaf_g(i), &
    !         tltrstem_g(i),tltrroot_g(i),' GRDAV'
    !
    !         !> File .CT03D
    !         write(74,8401)iday,realyr,vgbiomas_g(i), &
    !         gavglai_g(i), &
    !         gleafmas_g(i),bleafmas_g(i),stemmass_g(i), &
    !         rootmass_g(i),litrmass_g(i),soilcmas_g(i),' GRDAV'
    !
    !         !> File .CT04D
    !         write(75,8500)iday,realyr,ailcg_g(i),  &
    !         ailcb_g(i),(rmatctem_g(i,k),k=1,3), &
    !         veghght_g(i),rootdpth_g(i),roottemp_g(i), &
    !         slai_g(i),' GRDAV'
    !
    !         ! File .CT05D
    !         ! write(76,8601)iday,realyr,afrleaf_t(i,m), &
    !         !    afrstem_t(i,m),afrroot_t(i,m),  &
    !         !   , -999,   & ! lfstatus is kinda meaningless grid avg so set to -999
    !         !    ' GRDAV'
    !
    !         !> File *.CT06D
    !         if (dofire .or. lnduseon) then
    !             write(77,8800)iday,realyr,  &
    !             emit_co2_g(i),emit_co_g(i),emit_ch4_g(i), &
    !             emit_nmhc_g(i),emit_h2_g(i),emit_nox_g(i), &
    !             emit_n2o_g(i),emit_pm25_g(i),emit_tpm_g(i), &
    !             emit_tc_g(i),emit_oc_g(i),emit_bc_g(i), &
    !             burnfrac_g(i)*100., smfuncveg_g(i),lucemcom_g(i), &
    !             lucltrin_g(i),lucsocin_g(i), &
    !             grclarea(i),bterm_g(i),lterm_g(i),mterm_g(i),' GRDAV'
    !         end if
    !
    !
    !         !> File .CT08D
    !         if (.or. obswetf) then
    !             write(79,8810)iday,realyr,ch4WetSpec_g(i),  &
    !              wetfdyn_g(i),  &
    !             ch4WetDyn_g(i),  &
    !             ch4soills_g(i),' GRDAV'
    !         end if
    !
    !         ! FLAG FLAG FLAG
    !         !   ---------------peatland outputs-----------------------------------\
    !         !   - Note that YW's original code used a mixture of gat/row variables.
    !         !     The gat variables are replaced with row versions except for
    !         !     gppmosac_g which has no row version equivalent. Needs to be scattered out?
    !         !   - Also note that in YW's original code the arrays were hard-coded
    !         !     for just the 1st tile of the 1st grid point.
    !         !   - Leave gppmosac_g hard-coded to 1 for now,but should be changed !
    !         !   EC - Feb 2106.
    !         ! I think I will just remove these unless they are needed. JM Nov 2016.
    !         do m=1,nmtest
    !
    !             !   CT11D_G   convert moss gpp from umol/m2/s to g/m2/day
    !             write (93,6993) iday,realyr, &
    !             nppmossrow(i,m),armossrow(i,m), &!,gppmosac_g(1)*1.0377, &
    !             (fcancmxrow(i,m,j)*gppvegrow(i,m,j),j=1,icc),      &
    !             (fcancmxrow(i,m,j)*nppvegrow(i,m,j),j=1,icc),      &
    !             (fcancmxrow(i,m,j)*autoresvegrow(i,m,j),j=1,icc),  &
    !             (fcancmxrow(i,m,j)*hetroresvegrow(i,m,j),j=1,icc), &
    !             (fcancmxrow(i,m,j),j=1,icc)
    !             !   CT12D_G
    !             write (94,6993) iday,realyr,(veghghtrow(i,m,j),j=1,icc), &
    !             (rootdpthrow(i,m,j),j=1,icc),(ailcgrow(i,m,j),j=1,icc), &
    !             (fcancmxrow(i,m,j)*stemmassrow(i,m,j),j=1,icc), &
    !             (fcancmxrow(i,m,j)*rootmassrow(i,m,j),j=1,icc), &
    !             (fcancmxrow(i,m,j)*litrmassrow(i,m,j),j=1,icc), &
    !             (fcancmxrow(i,m,j)*gleafmasrow(i,m,j),j=1,icc), &
    !             (fcancmxrow(i,m,j)*bleafmasrow(i,m,j),j=1,icc)
    !
    !         end do
    !
    ! 6993            format(2i5,100f12.6)
    !
    !         !   ----------------YW March 27,2015 -------------------------------/
    !
    ! !                 if (PFTCompetition .or. lnduseon) then
    ! !
    ! !                     sumfare=0.0
    ! !                     if (onetile_perPFT) then
    ! !                         do m=1,nmos
    ! !                             sumfare=sumfare+FAREROT(i,m)
    ! !                         end do
    ! !                         write(78,8200)iday,realyr,(FAREROT(i,m)*100.,m=1,nmos),sumfare
    ! !                     else ! composite
    ! !                         do m=1,nmos
    ! !                             sumfare=0.0
    ! !                             do j=1,icc  ! m = 1
    ! !                                 sumfare=sumfare+fcancmxrow(i,m,j)
    ! !                             end do ! j
    ! !                             write(78,8200)iday,realyr,(fcancmxrow(i,m,j)*100.,j=1,icc),(1.0-sumfare)*100.,sumfare,' TILE ',m
    ! !                         end do ! m
    ! !                     end if ! mosaic/composite
    ! !                 end if ! PFTCompetition/lnduseon

    ! 80          continue ! nltest

    !         end if ! if write daily
    !     end if ! if write daily

    ! Deallocate the grid average and tile average structures
    deallocate(gpp_g,npp_g,nbp_g,socres_g,autores_g,litres_g,dstcemls3_g,litrfall_g , &
               rml_g,rms_g,rg_g,leaflitr_g,tltrstem_g,tltrroot_g,nep_g,hetrores_g,dstcemls_g , &
               humiftrs_g,rmr_g,tltrleaf_g,gavgltms_g,vgbiomas_g,gavglai_g,gavgscms_g,gleafmas_g , &
               bleafmas_g,stemmass_g,rootmass_g,litrmass_g,soilcmas_g,slai_g,ailcg_g,ailcb_g,veghght_g , &
               rootdpth_g,roottemp_g,totcmass_g,tcanoacc_out_g,burnfrac_g,smfuncveg_g , &
               lucemcom_g,lucltrin_g,lucsocin_g,emit_co2_g,emit_ch4_g,ch4WetSpec_g,wetfdyn_g,ch4WetDyn_g , &
               ch4soills_g,afrleaf_g,afrstem_g,afrroot_g,lfstatus_g,rmlvegrow_g,anvegrow_g, &
               rmatctem_g)

    deallocate(leaflitr_t,tltrleaf_t,tltrstem_t,tltrroot_t,ailcg_t,ailcb_t,rmatctem_t,veghght_t , &
               rootdpth_t,roottemp_t,slai_t,afrroot_t,afrleaf_t,afrstem_t,laimaxg_t,stemmass_t , &
               rootmass_t,litrmass_t,gleafmas_t,bleafmas_t,soilcmas_t,emit_co2_t,emit_ch4_t,smfuncveg_t)

  end subroutine ctem_daily_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_ctem_monthly_aw
  !> @{
  !> Accumulate and write out the monthly CTEM outputs. These are kept in pointer structures as
  !! this subroutine is called daily and we increment the daily values to produce a monthly value.
  !! The pointer to the monthly data structures (in ctemStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine ctem_monthly_aw (lonLocalIndex, latLocalIndex, nltest, nmtest, iday, realyr, nday, lastDOY)

    ! J. Melton Feb 2016.

    use classStateVars, only : class_rot
    use ctemStateVars,  only : ctem_tile_mo, vrot, ctem_grd_mo, c_switch, &
                               resetMonthEnd, ctem_mo, tracer
    use classicParams,  only : icc, iccp1, nmon, mmday, monthend, monthdays, seed, iccp2, ignd
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    integer, intent(in) :: realyr
    integer, intent(in) :: nday
    integer, intent(in) :: lastDOY

    ! pointers

    logical, pointer :: dofire
    logical, pointer :: lnduseon
    logical, pointer :: PFTCompetition
    logical, pointer :: doperpftoutput
    logical, pointer :: dopertileoutput
    logical, pointer :: transientOBSWETF
    integer, pointer :: fixedYearOBSWETF
    integer, pointer :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.

    real, pointer, dimension(:,:) :: FAREROT !< Fractional coverage of mosaic tile on modelled area

    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: laimaxg_mo
    real, pointer, dimension(:,:,:) :: stemmass_mo
    real, pointer, dimension(:,:,:) :: leafmass_mo
    real, pointer, dimension(:,:,:) :: rootmass_mo
    real, pointer, dimension(:,:,:) :: litrfallveg_mo
    real, pointer, dimension(:,:,:) :: humiftrsveg_mo
    real, pointer, dimension(:,:,:) :: npp_mo
    real, pointer, dimension(:,:,:) :: gpp_mo
    real, pointer, dimension(:,:,:) :: vgbiomas_mo
    real, pointer, dimension(:,:,:) :: autores_mo
    real, pointer, dimension(:,:,:) :: soilres_mo
    real, pointer, dimension(:,:,:) :: totcmass_mo
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litrmass_mo
    real, pointer, dimension(:,:,:) :: soilcmas_mo
    real, pointer, dimension(:,:,:) :: litres_mo
    real, pointer, dimension(:,:,:) :: soilcres_mo
    ! real, pointer, dimension(:,:,:,:) :: litrmass_mo
    ! real, pointer, dimension(:,:,:,:) :: soilcmas_mo
    ! real, pointer, dimension(:,:,:,:) :: litres_mo
    ! real, pointer, dimension(:,:,:,:) :: soilcres_mo
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: nep_mo
    real, pointer, dimension(:,:,:) :: hetrores_mo
    real, pointer, dimension(:,:,:) :: nbp_mo
    real, pointer, dimension(:,:,:) :: emit_co2_mo
    real, pointer, dimension(:,:,:) :: emit_co_mo
    real, pointer, dimension(:,:,:) :: emit_ch4_mo
    real, pointer, dimension(:,:,:) :: emit_nmhc_mo
    real, pointer, dimension(:,:,:) :: emit_h2_mo
    real, pointer, dimension(:,:,:) :: emit_nox_mo
    real, pointer, dimension(:,:,:) :: emit_n2o_mo
    real, pointer, dimension(:,:,:) :: emit_pm25_mo
    real, pointer, dimension(:,:,:) :: emit_tpm_mo
    real, pointer, dimension(:,:,:) :: emit_tc_mo
    real, pointer, dimension(:,:,:) :: emit_oc_mo
    real, pointer, dimension(:,:,:) :: emit_bc_mo
    real, pointer, dimension(:,:,:) :: bterm_mo
    real, pointer, dimension(:,:,:) :: mterm_mo
    real, pointer, dimension(:,:,:) :: burnfrac_mo
    real, pointer, dimension(:,:,:) :: smfuncveg_mo

    real, pointer, dimension(:,:) :: laimaxg_mo_t
    real, pointer, dimension(:,:) :: leafmass_mo_t
    real, pointer, dimension(:,:) :: stemmass_mo_t
    real, pointer, dimension(:,:) :: rootmass_mo_t
    real, pointer, dimension(:,:) :: litrfall_mo_t
    real, pointer, dimension(:,:) :: humiftrs_mo_t
    real, pointer, dimension(:,:) :: npp_mo_t
    real, pointer, dimension(:,:) :: gpp_mo_t
    real, pointer, dimension(:,:) :: vgbiomas_mo_t
    real, pointer, dimension(:,:) :: autores_mo_t
    real, pointer, dimension(:,:) :: soilres_mo_t
    real, pointer, dimension(:,:) :: totcmass_mo_t
    ! COMBAK PERLAY
    real, pointer, dimension(:,:) :: litrmass_mo_t
    real, pointer, dimension(:,:) :: soilcmas_mo_t
    real, pointer, dimension(:,:) :: litres_mo_t
    real, pointer, dimension(:,:) :: soilcres_mo_t
    ! real, pointer, dimension(:,:,:) :: litrmass_mo_t
    ! real, pointer, dimension(:,:,:) :: soilcmas_mo_t
    ! real, pointer, dimension(:,:,:) :: litres_mo_t
    ! real, pointer, dimension(:,:,:) :: soilcres_mo_t
    ! COMBAK PERLAY
    real, pointer, dimension(:,:) :: nep_mo_t
    real, pointer, dimension(:,:) :: hetrores_mo_t
    real, pointer, dimension(:,:) :: nbp_mo_t
    real, pointer, dimension(:,:) :: emit_co2_mo_t
    real, pointer, dimension(:,:) :: emit_co_mo_t
    real, pointer, dimension(:,:) :: emit_ch4_mo_t
    real, pointer, dimension(:,:) :: emit_nmhc_mo_t
    real, pointer, dimension(:,:) :: emit_h2_mo_t
    real, pointer, dimension(:,:) :: emit_nox_mo_t
    real, pointer, dimension(:,:) :: emit_n2o_mo_t
    real, pointer, dimension(:,:) :: emit_pm25_mo_t
    real, pointer, dimension(:,:) :: emit_tpm_mo_t
    real, pointer, dimension(:,:) :: emit_tc_mo_t
    real, pointer, dimension(:,:) :: emit_oc_mo_t
    real, pointer, dimension(:,:) :: emit_bc_mo_t
    real, pointer, dimension(:,:) :: burnfrac_mo_t
    real, pointer, dimension(:,:) :: smfuncveg_mo_t
    real, pointer, dimension(:,:) :: bterm_mo_t
    real, pointer, dimension(:,:) :: luc_emc_mo_t
    real, pointer, dimension(:,:) :: lterm_mo_t
    real, pointer, dimension(:,:) :: lucsocin_mo_t
    real, pointer, dimension(:,:) :: mterm_mo_t
    real, pointer, dimension(:,:) :: lucltrin_mo_t
    real, pointer, dimension(:,:) :: ch4WetSpec_mo_t
    real, pointer, dimension(:,:) :: wetfdyn_mo_t
    real, pointer, dimension(:,:) :: wetfpres_mo_t
    real, pointer, dimension(:,:) :: ch4WetDyn_mo_t
    real, pointer, dimension(:,:) :: ch4soills_mo_t
    real, pointer, dimension(:,:) :: wind_mo_t
    real, pointer, dimension(:,:) :: fProductDecomp_mo_t    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$


    logical, pointer, dimension(:,:,:) :: pftexistrow
    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:)   :: nbprow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    real, pointer, dimension(:,:,:) :: humiftrsvegrow
    ! real, pointer, dimension(:,:,:,:) :: litresvegrow
    ! real, pointer, dimension(:,:,:,:) :: soilcresvegrow
    ! real, pointer, dimension(:,:,:,:) :: litrmassrow
    ! real, pointer, dimension(:,:,:,:) :: soilcmasrow
    ! real, pointer, dimension(:,:,:,:) :: humiftrsvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: rmlvegaccrow
    real, pointer, dimension(:,:,:) :: rmsvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    real, pointer, dimension(:,:,:) :: afrrootrow
    real, pointer, dimension(:,:,:) :: ailcgrow
    real, pointer, dimension(:,:,:) :: emit_co2row
    real, pointer, dimension(:,:,:) :: emit_corow
    real, pointer, dimension(:,:,:) :: emit_ch4row
    real, pointer, dimension(:,:,:) :: emit_nmhcrow
    real, pointer, dimension(:,:,:) :: emit_h2row
    real, pointer, dimension(:,:,:) :: emit_noxrow
    real, pointer, dimension(:,:,:) :: emit_n2orow
    real, pointer, dimension(:,:,:) :: emit_pm25row
    real, pointer, dimension(:,:,:) :: emit_tpmrow
    real, pointer, dimension(:,:,:) :: emit_tcrow
    real, pointer, dimension(:,:,:) :: emit_ocrow
    real, pointer, dimension(:,:,:) :: emit_bcrow
    real, pointer, dimension(:,:) :: burnfracrow
    real, pointer, dimension(:,:,:) :: burnvegfrow
    real, pointer, dimension(:,:,:) :: smfuncvegrow
    real, pointer, dimension(:,:,:) :: btermrow
    real, pointer, dimension(:,:) :: ltermrow
    real, pointer, dimension(:,:,:) :: mtermrow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    real, pointer, dimension(:,:) :: ch4WetSpecrow
    real, pointer, dimension(:,:) :: wetfdynrow
    real, pointer, dimension(:,:) :: wetfrac_presrow
    real, pointer, dimension(:,:) :: ch4WetDynrow
    real, pointer, dimension(:,:) :: ch4soillsrow
    real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
    real, pointer, dimension(:,:,:) :: gleafmassrow
    real, pointer, dimension(:,:,:) :: bleafmassrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:,:) :: litrfallvegrow
    real, pointer, dimension(:,:) :: uvaccrow_m
    real, pointer, dimension(:,:) :: vvaccrow_m
    real, pointer, dimension(:) :: wetfrac_presgat

    real, pointer, dimension(:) :: laimaxg_mo_g
    real, pointer, dimension(:) :: stemmass_mo_g
    real, pointer, dimension(:) :: leafmass_mo_g
    real, pointer, dimension(:) :: rootmass_mo_g
    ! COMBAK PERLAY
    real, pointer, dimension(:) :: litres_mo_g
    real, pointer, dimension(:) :: soilcres_mo_g
    real, pointer, dimension(:) :: litrmass_mo_g
    real, pointer, dimension(:) :: soilcmas_mo_g
    ! real, pointer, dimension(:,:) :: litres_mo_g
    ! real, pointer, dimension(:,:) :: soilcres_mo_g
    ! real, pointer, dimension(:,:) :: litrmass_mo_g
    ! real, pointer, dimension(:,:) :: soilcmas_mo_g
    ! COMBAK PERLAY
    real, pointer, dimension(:) :: litrfall_mo_g
    real, pointer, dimension(:) :: humiftrs_mo_g
    real, pointer, dimension(:) :: npp_mo_g
    real, pointer, dimension(:) :: gpp_mo_g
    real, pointer, dimension(:) :: nep_mo_g
    real, pointer, dimension(:) :: nbp_mo_g
    real, pointer, dimension(:) :: hetrores_mo_g
    real, pointer, dimension(:) :: autores_mo_g
    real, pointer, dimension(:) :: soilres_mo_g
    real, pointer, dimension(:) :: vgbiomas_mo_g
    real, pointer, dimension(:) :: totcmass_mo_g
    real, pointer, dimension(:) :: emit_co2_mo_g
    real, pointer, dimension(:) :: emit_co_mo_g
    real, pointer, dimension(:) :: emit_ch4_mo_g
    real, pointer, dimension(:) :: emit_nmhc_mo_g
    real, pointer, dimension(:) :: emit_h2_mo_g
    real, pointer, dimension(:) :: emit_nox_mo_g
    real, pointer, dimension(:) :: emit_n2o_mo_g
    real, pointer, dimension(:) :: emit_pm25_mo_g
    real, pointer, dimension(:) :: emit_tpm_mo_g
    real, pointer, dimension(:) :: emit_tc_mo_g
    real, pointer, dimension(:) :: emit_oc_mo_g
    real, pointer, dimension(:) :: emit_bc_mo_g
    real, pointer, dimension(:) :: smfuncveg_mo_g
    real, pointer, dimension(:) :: luc_emc_mo_g
    real, pointer, dimension(:) :: lucltrin_mo_g
    real, pointer, dimension(:) :: lucsocin_mo_g
    real, pointer, dimension(:) :: burnfrac_mo_g
    real, pointer, dimension(:) :: bterm_mo_g
    real, pointer, dimension(:) :: lterm_mo_g
    real, pointer, dimension(:) :: mterm_mo_g
    real, pointer, dimension(:) :: ch4WetSpec_mo_g
    real, pointer, dimension(:) :: wetfdyn_mo_g
    real, pointer, dimension(:) :: wetfpres_mo_g
    real, pointer, dimension(:) :: ch4WetDyn_mo_g
    real, pointer, dimension(:) :: ch4soills_mo_g
    real, pointer, dimension(:) :: cProduct_mo_g          !< Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    real, pointer, dimension(:) :: fProductDecomp_mo_g    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-

    real, pointer, dimension(:,:) :: tracermossCMassrot      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer, dimension(:,:) :: tracermossLitrMassrot   !< Tracer mass in moss litter, \f$kg C/m^2\f$


    real, pointer, dimension(:,:,:) :: tracergLeafMassrot      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:) :: tracerbLeafMassrot      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:) :: tracerstemMassrot       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:) :: tracerrootMassrot       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    ! allocated with nlat,nmos,iccp2,ignd:
    real, pointer, dimension(:,:,:,:) :: tracerlitrMassrot       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:,:) :: tracersoilCMassrot      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    ! local
    integer :: i, m, j, nt, k
    real :: barefrac
    real :: sumfare
    integer :: NDMONTH
    integer :: imonth
    real, dimension(1) :: timeStamp
    real, dimension(icc) :: pftExist
    real :: oneOverDPM
    real, dimension(:), allocatable :: bulkLitterCarbon !< Temporary variable used to produce the bulk soil litter carbon quantity for output \f$[kg C m^{-2}]\f$
    real, dimension(:), allocatable :: bulkSoilCarbon   !< Temporary variable used to produce the bulk soil carbon quantity for output \f$[kg C m^{-2}]\f$
    real, dimension(:), allocatable :: bulkLitterResp   !< Temporary variable used to produce the bulk litter respiration quantity for output \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(:), allocatable :: bulkSoilResp     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output \f$[kg C m^{-2} s^{-1}]\f$

    ! point pointers

    dofire                => c_switch%dofire
    lnduseon              => c_switch%lnduseon
    PFTCompetition        => c_switch%PFTCompetition
    doperpftoutput        => c_switch%doperpftoutput
    dopertileoutput       => c_switch%dopertileoutput
    transientOBSWETF      => c_switch%transientOBSWETF
    fixedYearOBSWETF      => c_switch%fixedYearOBSWETF
    useTracer             => c_switch%useTracer

    FAREROT => class_rot%FAREROT

    pftexistrow           => vrot%pftexist
    fcancmxrow            => vrot%fcancmx
    laimaxg_mo            =>ctem_mo%laimaxg_mo
    stemmass_mo           =>ctem_mo%stemmass_mo
    leafmass_mo           =>ctem_mo%leafmass_mo
    rootmass_mo           =>ctem_mo%rootmass_mo
    litrfallveg_mo        =>ctem_mo%litrfallveg_mo
    humiftrsveg_mo        =>ctem_mo%humiftrsveg_mo
    npp_mo                =>ctem_mo%npp_mo
    gpp_mo                =>ctem_mo%gpp_mo
    vgbiomas_mo           =>ctem_mo%vgbiomas_mo
    autores_mo            =>ctem_mo%autores_mo
    soilres_mo            =>ctem_mo%soilres_mo
    totcmass_mo           =>ctem_mo%totcmass_mo
    litrmass_mo           =>ctem_mo%litrmass_mo
    soilcmas_mo           =>ctem_mo%soilcmas_mo
    nep_mo                =>ctem_mo%nep_mo
    litres_mo             =>ctem_mo%litres_mo
    soilcres_mo           =>ctem_mo%soilcres_mo
    hetrores_mo           =>ctem_mo%hetrores_mo
    nbp_mo                =>ctem_mo%nbp_mo
    emit_co2_mo           =>ctem_mo%emit_co2_mo
    emit_co_mo            =>ctem_mo%emit_co_mo
    emit_ch4_mo           =>ctem_mo%emit_ch4_mo
    emit_nmhc_mo          =>ctem_mo%emit_nmhc_mo
    emit_h2_mo            =>ctem_mo%emit_h2_mo
    emit_nox_mo           =>ctem_mo%emit_nox_mo
    emit_n2o_mo           =>ctem_mo%emit_n2o_mo
    emit_pm25_mo          =>ctem_mo%emit_pm25_mo
    emit_tpm_mo           =>ctem_mo%emit_tpm_mo
    emit_tc_mo            =>ctem_mo%emit_tc_mo
    emit_oc_mo            =>ctem_mo%emit_oc_mo
    emit_bc_mo            =>ctem_mo%emit_bc_mo
    bterm_mo              =>ctem_mo%bterm_mo
    mterm_mo              =>ctem_mo%mterm_mo
    burnfrac_mo           =>ctem_mo%burnfrac_mo
    smfuncveg_mo          =>ctem_mo%smfuncveg_mo

    laimaxg_mo_t          =>ctem_tile_mo%laimaxg_mo_t
    stemmass_mo_t         =>ctem_tile_mo%stemmass_mo_t
    leafmass_mo_t         =>ctem_tile_mo%leafmass_mo_t
    rootmass_mo_t         =>ctem_tile_mo%rootmass_mo_t
    litrfall_mo_t         =>ctem_tile_mo%litrfall_mo_t
    humiftrs_mo_t         =>ctem_tile_mo%humiftrs_mo_t
    npp_mo_t              =>ctem_tile_mo%npp_mo_t
    gpp_mo_t              =>ctem_tile_mo%gpp_mo_t
    vgbiomas_mo_t         =>ctem_tile_mo%vgbiomas_mo_t
    autores_mo_t          =>ctem_tile_mo%autores_mo_t
    soilres_mo_t          =>ctem_tile_mo%soilres_mo_t
    totcmass_mo_t         =>ctem_tile_mo%totcmass_mo_t
    litrmass_mo_t         =>ctem_tile_mo%litrmass_mo_t
    soilcmas_mo_t         =>ctem_tile_mo%soilcmas_mo_t
    nep_mo_t              =>ctem_tile_mo%nep_mo_t
    litres_mo_t           =>ctem_tile_mo%litres_mo_t
    soilcres_mo_t         =>ctem_tile_mo%soilcres_mo_t
    hetrores_mo_t         =>ctem_tile_mo%hetrores_mo_t
    nbp_mo_t              =>ctem_tile_mo%nbp_mo_t
    emit_co2_mo_t         =>ctem_tile_mo%emit_co2_mo_t
    emit_co_mo_t          =>ctem_tile_mo%emit_co_mo_t
    emit_ch4_mo_t         =>ctem_tile_mo%emit_ch4_mo_t
    emit_nmhc_mo_t        =>ctem_tile_mo%emit_nmhc_mo_t
    emit_h2_mo_t          =>ctem_tile_mo%emit_h2_mo_t
    emit_nox_mo_t         =>ctem_tile_mo%emit_nox_mo_t
    emit_n2o_mo_t         =>ctem_tile_mo%emit_n2o_mo_t
    emit_pm25_mo_t        =>ctem_tile_mo%emit_pm25_mo_t
    emit_tpm_mo_t         =>ctem_tile_mo%emit_tpm_mo_t
    emit_tc_mo_t          =>ctem_tile_mo%emit_tc_mo_t
    emit_oc_mo_t          =>ctem_tile_mo%emit_oc_mo_t
    emit_bc_mo_t          =>ctem_tile_mo%emit_bc_mo_t
    burnfrac_mo_t         =>ctem_tile_mo%burnfrac_mo_t
    smfuncveg_mo_t        =>ctem_tile_mo%smfuncveg_mo_t
    bterm_mo_t            =>ctem_tile_mo%bterm_mo_t
    luc_emc_mo_t          =>ctem_tile_mo%luc_emc_mo_t
    lterm_mo_t            =>ctem_tile_mo%lterm_mo_t
    lucsocin_mo_t         =>ctem_tile_mo%lucsocin_mo_t
    mterm_mo_t            =>ctem_tile_mo%mterm_mo_t
    lucltrin_mo_t         =>ctem_tile_mo%lucltrin_mo_t
    ch4WetSpec_mo_t       =>ctem_tile_mo%ch4WetSpec_mo_t
    wetfdyn_mo_t          =>ctem_tile_mo%wetfdyn_mo_t
    wetfpres_mo_t         =>ctem_tile_mo%wetfpres_mo_t
    ch4WetDyn_mo_t        =>ctem_tile_mo%ch4WetDyn_mo_t
    ch4soills_mo_t        =>ctem_tile_mo%ch4soills_mo_t
    wind_mo_t             =>ctem_tile_mo%wind_mo_t
    fProductDecomp_mo_t   =>ctem_tile_mo%fProductDecomp_mo_t

    gppvegrow         => vrot%gppveg
    nepvegrow         => vrot%nepveg
    nbpvegrow         => vrot%nbpveg
    nbprow            => vrot%nbp
    nppvegrow         => vrot%nppveg
    hetroresvegrow    => vrot%hetroresveg
    autoresvegrow     => vrot%autoresveg
    litresvegrow      => vrot%litresveg
    soilcresvegrow    => vrot%soilcresveg
    rmlvegaccrow      => vrot%rmlvegacc
    rmsvegrow         => vrot%rmsveg
    rmrvegrow         => vrot%rmrveg
    rgvegrow          => vrot%rgveg
    afrrootrow        => vrot%afrroot
    ailcgrow          => vrot%ailcg
    emit_co2row       => vrot%emit_co2
    emit_corow        => vrot%emit_co
    emit_ch4row       => vrot%emit_ch4
    emit_nmhcrow      => vrot%emit_nmhc
    emit_h2row        => vrot%emit_h2
    emit_noxrow       => vrot%emit_nox
    emit_n2orow       => vrot%emit_n2o
    emit_pm25row      => vrot%emit_pm25
    emit_tpmrow       => vrot%emit_tpm
    emit_tcrow        => vrot%emit_tc
    emit_ocrow        => vrot%emit_oc
    emit_bcrow        => vrot%emit_bc
    burnfracrow       => vrot%burnfrac
    burnvegfrow       => vrot%burnvegf
    smfuncvegrow      => vrot%smfuncveg
    btermrow          => vrot%bterm
    ltermrow          => vrot%lterm
    mtermrow          => vrot%mterm
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    ch4WetSpecrow     => vrot%ch4WetSpec
    wetfdynrow        => vrot%wetfdyn
    wetfrac_presrow   => vrot%wetfrac_pres

    ch4WetDynrow      => vrot%ch4WetDyn
    ch4soillsrow      => vrot%ch4_soills
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    vgbiomas_vegrow   => vrot%vgbiomas_veg
    gleafmassrow      => vrot%gleafmas
    bleafmassrow      => vrot%bleafmas
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    uvaccrow_m        => vrot%uvaccrow_m
    vvaccrow_m        => vrot%vvaccrow_m
    litrfallvegrow    => vrot%litrfallveg
    humiftrsvegrow    => vrot%humiftrsveg

    tracerGLeafMassrot   => tracer%gLeafMassrot
    tracerBLeafMassrot   => tracer%bLeafMassrot
    tracerStemMassrot    => tracer%stemMassrot
    tracerRootMassrot    => tracer%rootMassrot
    tracerLitrMassrot    => tracer%litrMassrot
    tracerSoilCMassrot   => tracer%soilCMassrot
    tracerMossCMassrot   => tracer%mossCMassrot
    tracerMossLitrMassrot => tracer%mossLitrMassrot

    laimaxg_mo_g        =>ctem_grd_mo%laimaxg_mo_g
    stemmass_mo_g       =>ctem_grd_mo%stemmass_mo_g
    leafmass_mo_g       =>ctem_grd_mo%leafmass_mo_g
    rootmass_mo_g       =>ctem_grd_mo%rootmass_mo_g
    litrmass_mo_g       =>ctem_grd_mo%litrmass_mo_g
    soilcmas_mo_g       =>ctem_grd_mo%soilcmas_mo_g
    litrfall_mo_g       =>ctem_grd_mo%litrfall_mo_g
    humiftrs_mo_g       =>ctem_grd_mo%humiftrs_mo_g
    npp_mo_g            =>ctem_grd_mo%npp_mo_g
    gpp_mo_g            =>ctem_grd_mo%gpp_mo_g
    nep_mo_g            =>ctem_grd_mo%nep_mo_g
    nbp_mo_g            =>ctem_grd_mo%nbp_mo_g
    hetrores_mo_g       =>ctem_grd_mo%hetrores_mo_g
    autores_mo_g        =>ctem_grd_mo%autores_mo_g
    soilres_mo_g        =>ctem_grd_mo%soilres_mo_g
    litres_mo_g         =>ctem_grd_mo%litres_mo_g
    soilcres_mo_g       =>ctem_grd_mo%soilcres_mo_g
    vgbiomas_mo_g       =>ctem_grd_mo%vgbiomas_mo_g
    totcmass_mo_g       =>ctem_grd_mo%totcmass_mo_g
    emit_co2_mo_g       =>ctem_grd_mo%emit_co2_mo_g
    emit_co_mo_g        =>ctem_grd_mo%emit_co_mo_g
    emit_ch4_mo_g       =>ctem_grd_mo%emit_ch4_mo_g
    emit_nmhc_mo_g      =>ctem_grd_mo%emit_nmhc_mo_g
    emit_h2_mo_g        =>ctem_grd_mo%emit_h2_mo_g
    emit_nox_mo_g       =>ctem_grd_mo%emit_nox_mo_g
    emit_n2o_mo_g       =>ctem_grd_mo%emit_n2o_mo_g
    emit_pm25_mo_g      =>ctem_grd_mo%emit_pm25_mo_g
    emit_tpm_mo_g       =>ctem_grd_mo%emit_tpm_mo_g
    emit_tc_mo_g        =>ctem_grd_mo%emit_tc_mo_g
    emit_oc_mo_g        =>ctem_grd_mo%emit_oc_mo_g
    emit_bc_mo_g        =>ctem_grd_mo%emit_bc_mo_g
    smfuncveg_mo_g      =>ctem_grd_mo%smfuncveg_mo_g
    luc_emc_mo_g        =>ctem_grd_mo%luc_emc_mo_g
    lucltrin_mo_g       =>ctem_grd_mo%lucltrin_mo_g
    lucsocin_mo_g       =>ctem_grd_mo%lucsocin_mo_g
    burnfrac_mo_g       =>ctem_grd_mo%burnfrac_mo_g
    bterm_mo_g          =>ctem_grd_mo%bterm_mo_g
    lterm_mo_g          =>ctem_grd_mo%lterm_mo_g
    mterm_mo_g          =>ctem_grd_mo%mterm_mo_g
    ch4WetSpec_mo_g     =>ctem_grd_mo%ch4WetSpec_mo_g
    wetfdyn_mo_g        =>ctem_grd_mo%wetfdyn_mo_g
    wetfpres_mo_g       =>ctem_grd_mo%wetfpres_mo_g
    ch4WetDyn_mo_g      =>ctem_grd_mo%ch4WetDyn_mo_g
    ch4soills_mo_g      =>ctem_grd_mo%ch4soills_mo_g
    cProduct_mo_g       =>ctem_grd_mo%cProduct_mo_g
    fProductDecomp_mo_g =>ctem_grd_mo%fProductDecomp_mo_g

    !> ------------

    !> Accumulate monthly outputs

    ! Find which month you are in so we know the number of days per month.
    do nt = 1,nmon
      if (iday <= monthend(nt + 1)) then
        oneOverDPM = 1./real(monthdays(nt))
        exit
      else
        cycle
      end if
    end do

    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do m = 1,nmtest
      do j = 1,icc

        !> Accumulate monthly outputs at the per PFT level.
        if (ailcgrow(i,m,j) > laimaxg_mo(i,m,j)) then
          laimaxg_mo(i,m,j) = ailcgrow(i,m,j)
        end if

        npp_mo(i,m,j) = npp_mo(i,m,j) + nppvegrow(i,m,j) * oneOverDPM
        gpp_mo(i,m,j) = gpp_mo(i,m,j) + gppvegrow(i,m,j) * oneOverDPM
        nep_mo(i,m,j) = nep_mo(i,m,j) + nepvegrow(i,m,j) * oneOverDPM
        ! NOTE: This NBP does not include LUC product pool contributions since they are
        ! not per PFT but rather per tile
        nbp_mo(i,m,j) = nbp_mo(i,m,j) + nbpvegrow(i,m,j) * oneOverDPM
        hetrores_mo(i,m,j) = hetrores_mo(i,m,j) + hetroresvegrow(i,m,j) * oneOverDPM
        autores_mo(i,m,j) = autores_mo(i,m,j) + autoresvegrow(i,m,j) * oneOverDPM
        ! Calculate the total soil respiration (Rh + root contributions of Rm + Rg)
        soilres_mo(i,m,j) = soilres_mo(i,m,j) + (hetroresvegrow(i,m,j) + rmrvegrow(i,m,j) &
                                              + rgvegrow(i,m,j) * afrrootrow(i,m,j)) * oneOverDPM
        ! COMBAK PERLAY
        litres_mo(i,m,j) = litres_mo(i,m,j) + litresvegrow(i,m,j) * oneOverDPM
        soilcres_mo(i,m,j) = soilcres_mo(i,m,j) + soilcresvegrow(i,m,j) * oneOverDPM
        ! Let accumulate,not fluxes nor meant to be mean values.
        humiftrsveg_mo(i,m,j) = humiftrsveg_mo(i,m,j) + humiftrsvegrow(i,m,j)
        ! do k = 1,ignd
        !   litres_mo(i,m,j,k) = litres_mo(i,m,j,k) + litresvegrow(i,m,j,k)*oneOverDPM
        !   soilcres_mo(i,m,j,k) = soilcres_mo(i,m,j,k) + soilcresvegrow(i,m,j,k)*oneOverDPM
        !   ! Let accumulate,not fluxes nor meant to be mean values.
        !   humiftrsveg_mo(i,m,j) = humiftrsveg_mo(i,m,j) + humiftrsvegrow(i,m,j,k)
        ! end do
        ! COMBAK PERLAY
        emit_co2_mo(i,m,j) = emit_co2_mo(i,m,j) + emit_co2row(i,m,j) * oneOverDPM
        emit_co_mo(i,m,j) = emit_co_mo(i,m,j) + emit_corow(i,m,j) * oneOverDPM
        emit_ch4_mo(i,m,j) = emit_ch4_mo(i,m,j) + emit_ch4row(i,m,j) * oneOverDPM
        emit_nmhc_mo(i,m,j) = emit_nmhc_mo(i,m,j) + emit_nmhcrow(i,m,j) * oneOverDPM
        emit_h2_mo(i,m,j) = emit_h2_mo(i,m,j) + emit_h2row(i,m,j) * oneOverDPM
        emit_nox_mo(i,m,j) = emit_nox_mo(i,m,j) + emit_noxrow(i,m,j) * oneOverDPM
        emit_n2o_mo(i,m,j) = emit_n2o_mo(i,m,j) + emit_n2orow(i,m,j) * oneOverDPM
        emit_pm25_mo(i,m,j) = emit_pm25_mo(i,m,j) + emit_pm25row(i,m,j) * oneOverDPM
        emit_tpm_mo(i,m,j) = emit_tpm_mo(i,m,j) + emit_tpmrow(i,m,j) * oneOverDPM
        emit_tc_mo(i,m,j) = emit_tc_mo(i,m,j) + emit_tcrow(i,m,j) * oneOverDPM
        emit_oc_mo(i,m,j) = emit_oc_mo(i,m,j) + emit_ocrow(i,m,j) * oneOverDPM
        emit_bc_mo(i,m,j) = emit_bc_mo(i,m,j) + emit_bcrow(i,m,j) * oneOverDPM
        bterm_mo(i,m,j) = bterm_mo(i,m,j) + btermrow(i,m,j) * oneOverDPM
        mterm_mo(i,m,j) = mterm_mo(i,m,j) + mtermrow(i,m,j) * oneOverDPM
        smfuncveg_mo(i,m,j) = smfuncveg_mo(i,m,j) + smfuncvegrow(i,m,j) * oneOverDPM

        ! Let accumulate,not fluxes nor meant to be mean values.
        burnfrac_mo(i,m,j) = burnfrac_mo(i,m,j) + burnvegfrow(i,m,j)
        litrfallveg_mo(i,m,j) = litrfallveg_mo(i,m,j) + litrfallvegrow(i,m,j)

      end do ! j

      !> Also do the bare ground
      nep_mo(i,m,iccp1) = nep_mo(i,m,iccp1) + nepvegrow(i,m,iccp1) * oneOverDPM
      nbp_mo(i,m,iccp1) = nbp_mo(i,m,iccp1) + nbpvegrow(i,m,iccp1) * oneOverDPM
      hetrores_mo(i,m,iccp1) = hetrores_mo(i,m,iccp1) + hetroresvegrow(i,m,iccp1) * oneOverDPM
      ! COMBAK PERLAY
      litres_mo(i,m,iccp1) = litres_mo(i,m,iccp1) + litresvegrow(i,m,iccp1) * oneOverDPM
      soilcres_mo(i,m,iccp1) = soilcres_mo(i,m,iccp1) &
                               + soilcresvegrow(i,m,iccp1) * oneOverDPM
      ! do k = 1,ignd
      !   litres_mo(i,m,iccp1,k) = litres_mo(i,m,iccp1,k) + litresvegrow(i,m,iccp1,k) * oneOverDPM
      !   soilcres_mo(i,m,iccp1,k) = soilcres_mo(i,m,iccp1,k) &
      !                            + soilcresvegrow(i,m,iccp1,k) * oneOverDPM
      ! end do
      ! COMBAK PERLAY
      !> Accumulate monthly outputs at the per tile level.
      luc_emc_mo_t(i,m) = luc_emc_mo_t(i,m) + lucemcomrow(i,m) * oneOverDPM
      lucsocin_mo_t(i,m) = lucsocin_mo_t(i,m) + lucsocinrow(i,m) * oneOverDPM
      lucltrin_mo_t(i,m) = lucltrin_mo_t(i,m) + lucltrinrow(i,m) * oneOverDPM
      ch4WetSpec_mo_t(i,m) = ch4WetSpec_mo_t(i,m) + ch4WetSpecrow(i,m) * oneOverDPM
      wetfdyn_mo_t(i,m) = wetfdyn_mo_t(i,m) + wetfdynrow(i,m) * oneOverDPM
      wetfpres_mo_t(i,m) = wetfpres_mo_t(i,m) + wetfrac_presrow(i,m) * oneOverDPM
      ch4WetDyn_mo_t(i,m) = ch4WetDyn_mo_t(i,m) + ch4WetDynrow(i,m) * oneOverDPM
      ch4soills_mo_t(i,m) = ch4soills_mo_t(i,m) + ch4soillsrow(i,m) * oneOverDPM
      lterm_mo_t(i,m) = lterm_mo_t(i,m) + ltermrow(i,m) * oneOverDPM
      ! wind_mo_t(i,m) = wind_mo_t(i,m) + (sqrt(uvaccrow_m(i,m)**2.0 + vvaccrow_m(i,m)**2.0))*3.6 !> take mean wind speed and convert to km/h

      ! COMBAK PERLAY
      ! NOTE: LUC product pools are kept in layer 1.
      fProductDecomp_mo_t(i,m) = fProductDecomp_mo_t(i,m) &
                                 + (soilcresvegrow(i,m,iccp2) + litresvegrow(i,m,iccp2)) * oneOverDPM
      ! fProductDecomp_mo_t(i,m) = fProductDecomp_mo_t(i,m) &
      !                           + (soilcresvegrow(i,m,iccp2,1) + litresvegrow(i,m,iccp2,1)) * oneOverDPM
      ! COMBAK PERLAY

    end do ! loop 863 ! m

    do nt = 1,nmon

      if (iday == mmday(nt)) then

        !> Do the mid-month variables (these are not accumulated,we just keep the mid month value for printing in the monthly file)

        do m = 1,nmtest
          do j = 1,icc

            vgbiomas_mo(i,m,j) = vgbiomas_vegrow(i,m,j)
            stemmass_mo(i,m,j) = stemmassrow(i,m,j)
            leafmass_mo(i,m,j) = gleafmassrow(i,m,j) + bleafmassrow(i,m,j)
            rootmass_mo(i,m,j) = rootmassrow(i,m,j)
            totcmass_mo(i,m,j) = vgbiomas_vegrow(i,m,j)
            ! COMBAK PERLAY
            litrmass_mo(i,m,j) = litrmassrow(i,m,j)
            soilcmas_mo(i,m,j) = soilcmasrow(i,m,j)
            totcmass_mo(i,m,j) = totcmass_mo(i,m,j) + litrmassrow(i,m,j) + soilcmasrow(i,m,j)
            ! do k = 1,ignd
            !   litrmass_mo(i,m,j,k) = litrmassrow(i,m,j,k)
            !   soilcmas_mo(i,m,j,k) = soilcmasrow(i,m,j,k)
            !   totcmass_mo(i,m,j) = totcmass_mo(i,m,j) + litrmassrow(i,m,j,k) + soilcmasrow(i,m,j,k)
            ! end do
            ! COMBAK PERLAY
          end do ! loop 867

          !> Do the bare fraction too
          ! COMBAK PERLAY
          litrmass_mo(i,m,iccp1) = litrmassrow(i,m,iccp1)
          soilcmas_mo(i,m,iccp1) = soilcmasrow(i,m,iccp1)
          totcmass_mo(i,m,iccp1) = totcmass_mo(i,m,iccp1) &
                                   + soilcmasrow(i,m,iccp1) + litrmassrow(i,m,iccp1)
          ! do k = 1,ignd
          !   litrmass_mo(i,m,iccp1,k) = litrmassrow(i,m,iccp1,k)
          !   soilcmas_mo(i,m,iccp1,k) = soilcmasrow(i,m,iccp1,k)
          !   totcmass_mo(i,m,iccp1) = totcmass_mo(i,m,iccp1) &
          !                           + soilcmasrow(i,m,iccp1,k) + litrmassrow(i,m,iccp1,k)
          ! end do
          ! COMBAK PERLAY
          barefrac = 1.0

          !> Now find the per tile values:
          do j = 1,icc
            vgbiomas_mo_t(i,m) = vgbiomas_mo_t(i,m) + vgbiomas_mo(i,m,j) * fcancmxrow(i,m,j)
            ! COMBAK PERLAY
            litrmass_mo_t(i,m) = litrmass_mo_t(i,m) + litrmass_mo(i,m,j) * fcancmxrow(i,m,j)
            soilcmas_mo_t(i,m) = soilcmas_mo_t(i,m) + soilcmas_mo(i,m,j) * fcancmxrow(i,m,j)
            ! do k = 1,ignd
            !   litrmass_mo_t(i,m,k) = litrmass_mo_t(i,m,k) + litrmass_mo(i,m,j,k) * fcancmxrow(i,m,j)
            !   soilcmas_mo_t(i,m,k) = soilcmas_mo_t(i,m,k) + soilcmas_mo(i,m,j,k) * fcancmxrow(i,m,j)
            ! end do
            ! COMBAK PERLAY
            leafmass_mo_t(i,m) = leafmass_mo_t(i,m) + leafmass_mo(i,m,j) * fcancmxrow(i,m,j)
            stemmass_mo_t(i,m) = stemmass_mo_t(i,m) + stemmass_mo(i,m,j) * fcancmxrow(i,m,j)
            rootmass_mo_t(i,m) = rootmass_mo_t(i,m) + rootmass_mo(i,m,j) * fcancmxrow(i,m,j)
            totcmass_mo_t(i,m) = totcmass_mo_t(i,m) + totcmass_mo(i,m,j) * fcancmxrow(i,m,j)
            barefrac = barefrac - fcancmxrow(i,m,j)
          end do

          !> Also add in the bare fraction contributions.
          ! COMBAK PERLAY
          litrmass_mo_t(i,m) = litrmass_mo_t(i,m) + litrmass_mo(i,m,iccp1) * barefrac
          soilcmas_mo_t(i,m) = soilcmas_mo_t(i,m) + soilcmas_mo(i,m,iccp1) * barefrac
          totcmass_mo_t(i,m) = totcmass_mo_t(i,m) &
                               + (litrmass_mo(i,m,iccp1) + soilcmas_mo(i,m,iccp1)) * barefrac
          ! do k = 1,ignd
          !   litrmass_mo_t(i,m,k) = litrmass_mo_t(i,m,k) + litrmass_mo(i,m,iccp1,k) * barefrac
          !   soilcmas_mo_t(i,m,k) = soilcmas_mo_t(i,m,k) + soilcmas_mo(i,m,iccp1,k) * barefrac
          !   totcmass_mo_t(i,m) = totcmass_mo_t(i,m) &
          !                       + (litrmass_mo(i,m,iccp1,k) + soilcmas_mo(i,m,iccp1,k)) * barefrac
          ! end do
          ! COMBAK PERLAY

          !> Now find the gridcell level values:
          vgbiomas_mo_g(i) = vgbiomas_mo_g(i) + vgbiomas_mo_t(i,m) * FAREROT(i,m)
          ! COMBAK PERLAY
          litrmass_mo_g(i) = litrmass_mo_g(i) + litrmass_mo_t(i,m) * FAREROT(i,m)
          soilcmas_mo_g(i) = soilcmas_mo_g(i) + soilcmas_mo_t(i,m) * FAREROT(i,m)
          ! do k = 1,ignd
          !   litrmass_mo_g(i,k) = litrmass_mo_g(i,k) + litrmass_mo_t(i,m,k) * FAREROT(i,m)
          !   soilcmas_mo_g(i,k) = soilcmas_mo_g(i,k) + soilcmas_mo_t(i,m,k) * FAREROT(i,m)
          ! end do
          ! COMBAK PERLAY
          leafmass_mo_g(i) = leafmass_mo_g(i) + leafmass_mo_t(i,m) * FAREROT(i,m)
          stemmass_mo_g(i) = stemmass_mo_g(i) + stemmass_mo_t(i,m) * FAREROT(i,m)
          rootmass_mo_g(i) = rootmass_mo_g(i) + rootmass_mo_t(i,m) * FAREROT(i,m)
          totcmass_mo_g(i) = totcmass_mo_g(i) + totcmass_mo_t(i,m) * FAREROT(i,m)

          ! Including the LUC product pools. They are per tile values and
          ! are assumed to occupy the whole tile. Only kept in layer 1.
          ! COMBAK PERLAY
          cProduct_mo_g(i) = cProduct_mo_g(i) &
                             + (litrmassrow(i,m,iccp2) &
                             + soilcmasrow(i,m,iccp2)) * FAREROT(i,m)
          ! cProduct_mo_g(i) = cProduct_mo_g(i) &
          !                    + (litrmassrow(i,m,iccp2,1) &
          !                    + soilcmasrow(i,m,iccp2,1)) * FAREROT(i,m)
          ! COMBAK PERLAY

        end do ! loop 866  ! nmtest loop.

      end if ! mmday (mid-month instantaneous value)

      if (iday == monthend(nt + 1)) then

        !> Do the end of month variables
        ndmonth = (monthend(nt + 1) - monthend(nt)) * nday

        do m = 1,nmtest

          barefrac = 1.0

          do j = 1,icc

            !> Find the monthly outputs at the per tile level from the outputs at the per PFT level
            npp_mo_t(i,m) = npp_mo_t(i,m) + npp_mo(i,m,j) * fcancmxrow(i,m,j)
            gpp_mo_t(i,m) = gpp_mo_t(i,m) + gpp_mo(i,m,j) * fcancmxrow(i,m,j)
            nep_mo_t(i,m) = nep_mo_t(i,m) + nep_mo(i,m,j) * fcancmxrow(i,m,j)
            ! nbp_mo_t(i,m)=nbp_mo_t(i,m)+nbp_mo(i,m,j)*fcancmxrow(i,m,j)
            hetrores_mo_t(i,m) = hetrores_mo_t(i,m) + hetrores_mo(i,m,j) * fcancmxrow(i,m,j)
            autores_mo_t(i,m) = autores_mo_t(i,m) + autores_mo(i,m,j) * fcancmxrow(i,m,j)
            soilres_mo_t(i,m) = soilres_mo_t(i,m) + soilres_mo(i,m,j) * fcancmxrow(i,m,j)

            ! COMBAK PERLAY
            litres_mo_t(i,m) = litres_mo_t(i,m) + litres_mo(i,m,j) * fcancmxrow(i,m,j)
            soilcres_mo_t(i,m) = soilcres_mo_t(i,m) + soilcres_mo(i,m,j) * fcancmxrow(i,m,j)
            ! do k = 1,ignd
            !   litres_mo_t(i,m,k) = litres_mo_t(i,m,k) + litres_mo(i,m,j,k) * fcancmxrow(i,m,j)
            !   soilcres_mo_t(i,m,k) = soilcres_mo_t(i,m,k) + soilcres_mo(i,m,j,k) * fcancmxrow(i,m,j)
            ! end do
            ! COMBAK PERLAY
            emit_co2_mo_t(i,m) = emit_co2_mo_t(i,m) + emit_co2_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_co_mo_t(i,m) = emit_co_mo_t(i,m) + emit_co_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_ch4_mo_t(i,m) = emit_ch4_mo_t(i,m) + emit_ch4_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_nmhc_mo_t(i,m) = emit_nmhc_mo_t(i,m) + emit_nmhc_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_h2_mo_t(i,m) = emit_h2_mo_t(i,m) + emit_h2_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_nox_mo_t(i,m) = emit_nox_mo_t(i,m) + emit_nox_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_n2o_mo_t(i,m) = emit_n2o_mo_t(i,m) + emit_n2o_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_pm25_mo_t(i,m) = emit_pm25_mo_t(i,m) + emit_pm25_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_tpm_mo_t(i,m) = emit_tpm_mo_t(i,m) + emit_tpm_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_tc_mo_t(i,m) = emit_tc_mo_t(i,m) + emit_tc_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_oc_mo_t(i,m) = emit_oc_mo_t(i,m) + emit_oc_mo(i,m,j) * fcancmxrow(i,m,j)
            emit_bc_mo_t(i,m) = emit_bc_mo_t(i,m) + emit_bc_mo(i,m,j) * fcancmxrow(i,m,j)
            bterm_mo_t(i,m) = bterm_mo_t(i,m) + bterm_mo(i,m,j) * fcancmxrow(i,m,j)
            mterm_mo_t(i,m) = mterm_mo_t(i,m) + mterm_mo(i,m,j) * fcancmxrow(i,m,j)
            smfuncveg_mo_t(i,m) = smfuncveg_mo_t(i,m) + smfuncveg_mo(i,m,j) * fcancmxrow(i,m,j)
            burnfrac_mo_t(i,m) = burnfrac_mo_t(i,m) + burnfrac_mo(i,m,j) * fcancmxrow(i,m,j)
            laimaxg_mo_t(i,m) = laimaxg_mo_t(i,m) + laimaxg_mo(i,m,j) * fcancmxrow(i,m,j)
            litrfall_mo_t(i,m) = litrfall_mo_t(i,m) + litrfallveg_mo(i,m,j) * fcancmxrow(i,m,j)
            humiftrs_mo_t(i,m) = humiftrs_mo_t(i,m) + humiftrsveg_mo(i,m,j) * fcancmxrow(i,m,j)
            barefrac = barefrac - fcancmxrow(i,m,j)

          end do ! j

          nep_mo_t(i,m) = nep_mo_t(i,m) + nep_mo(i,m,iccp1) * barefrac
          ! nbp_mo_t(i,m)=nbp_mo_t(i,m)+nbp_mo(i,m,iccp1)*barefrac
          hetrores_mo_t(i,m) = hetrores_mo_t(i,m) + hetrores_mo(i,m,iccp1) * barefrac
          humiftrs_mo_t(i,m) = humiftrs_mo_t(i,m) + humiftrsveg_mo(i,m,iccp1) * barefrac
          ! COMBAK PERLAY
          litres_mo_t(i,m) = litres_mo_t(i,m) + litres_mo(i,m,iccp1) * barefrac
          soilcres_mo_t(i,m) = soilcres_mo_t(i,m) + soilcres_mo(i,m,iccp1) * barefrac
          ! do k = 1,ignd
          !   litres_mo_t(i,m,k) = litres_mo_t(i,m,k) + litres_mo(i,m,iccp1,k) * barefrac
          !   soilcres_mo_t(i,m,k) = soilcres_mo_t(i,m,k) + soilcres_mo(i,m,iccp1,k) * barefrac
          ! end do
          ! COMBAK PERLAY

          ! NOTE: NBP is a special case here. The LUC product pool contributions are not
          ! per PFT, they exist uniformly across a tile, so they are not inclued in the
          ! nbp_mo calculation. Instead we need to use the nbp, not nbpveg variable
          ! for per tile and per gridcell outputting.
          nbp_mo_t(i,m) = nbprow(i,m)

          !> Find the monthly outputs at the per grid cell level from the outputs at the per tile level
          npp_mo_g(i) = npp_mo_g(i) + npp_mo_t(i,m) * FAREROT(i,m)
          gpp_mo_g(i) = gpp_mo_g(i) + gpp_mo_t(i,m) * FAREROT(i,m)
          nep_mo_g(i) = nep_mo_g(i) + nep_mo_t(i,m) * FAREROT(i,m)
          nbp_mo_g(i) = nbp_mo_g(i) + nbp_mo_t(i,m) * FAREROT(i,m)
          hetrores_mo_g(i) = hetrores_mo_g(i) + hetrores_mo_t(i,m) * FAREROT(i,m)
          autores_mo_g(i) = autores_mo_g(i) + autores_mo_t(i,m) * FAREROT(i,m)
          soilres_mo_g(i) = soilres_mo_g(i) + soilres_mo_t(i,m) * FAREROT(i,m)
          ! COMBAK PERLAY
          litres_mo_g(i) = litres_mo_g(i) + litres_mo_t(i,m) * FAREROT(i,m)
          soilcres_mo_g(i) = soilcres_mo_g(i) + soilcres_mo_t(i,m) * FAREROT(i,m)
          ! do k = 1,ignd
          !   litres_mo_g(i,k) = litres_mo_g(i,k) + litres_mo_t(i,m,k) * FAREROT(i,m)
          !   soilcres_mo_g(i,k) = soilcres_mo_g(i,k) + soilcres_mo_t(i,m,k) * FAREROT(i,m)
          ! end do
          ! COMBAK PERLAY
          laimaxg_mo_g(i) = laimaxg_mo_g(i) + laimaxg_mo_t(i,m) * FAREROT(i,m)
          emit_co2_mo_g(i) = emit_co2_mo_g(i) + emit_co2_mo_t(i,m) * FAREROT(i,m)
          emit_co_mo_g(i) = emit_co_mo_g(i) + emit_co_mo_t(i,m) * FAREROT(i,m)
          emit_ch4_mo_g(i) = emit_ch4_mo_g(i) + emit_ch4_mo_t(i,m) * FAREROT(i,m)
          emit_nmhc_mo_g(i) = emit_nmhc_mo_g(i) + emit_nmhc_mo_t(i,m) * FAREROT(i,m)
          emit_h2_mo_g(i) = emit_h2_mo_g(i) + emit_h2_mo_t(i,m) * FAREROT(i,m)
          emit_nox_mo_g(i) = emit_nox_mo_g(i) + emit_nox_mo_t(i,m) * FAREROT(i,m)
          emit_n2o_mo_g(i) = emit_n2o_mo_g(i) + emit_n2o_mo_t(i,m) * FAREROT(i,m)
          emit_pm25_mo_g(i) = emit_pm25_mo_g(i) + emit_pm25_mo_t(i,m) * FAREROT(i,m)
          emit_tpm_mo_g(i) = emit_tpm_mo_g(i) + emit_tpm_mo_t(i,m) * FAREROT(i,m)
          emit_tc_mo_g(i) = emit_tc_mo_g(i) + emit_tc_mo_t(i,m) * FAREROT(i,m)
          emit_oc_mo_g(i) = emit_oc_mo_g(i) + emit_oc_mo_t(i,m) * FAREROT(i,m)
          emit_bc_mo_g(i) = emit_bc_mo_g(i) + emit_bc_mo_t(i,m) * FAREROT(i,m)
          burnfrac_mo_g(i) = burnfrac_mo_g(i) + burnfrac_mo_t(i,m) * FAREROT(i,m)
          luc_emc_mo_g(i) = luc_emc_mo_g(i) + luc_emc_mo_t(i,m) * FAREROT(i,m)
          lucsocin_mo_g(i) = lucsocin_mo_g(i) + lucsocin_mo_t(i,m) * FAREROT(i,m)
          lucltrin_mo_g(i) = lucltrin_mo_g(i) + lucltrin_mo_t(i,m) * FAREROT(i,m)
          ch4WetSpec_mo_g(i) = ch4WetSpec_mo_g(i) + ch4WetSpec_mo_t(i,m) * FAREROT(i,m)
          wetfdyn_mo_g(i) = wetfdyn_mo_g(i) + wetfdyn_mo_t(i,m) * FAREROT(i,m)
          wetfpres_mo_g(i) = wetfpres_mo_g(i) + wetfpres_mo_t(i,m) * FAREROT(i,m)
          ch4WetDyn_mo_g(i) = ch4WetDyn_mo_g(i) + ch4WetDyn_mo_t(i,m) * FAREROT(i,m)
          ch4soills_mo_g(i) = ch4soills_mo_g(i) + ch4soills_mo_t(i,m) * FAREROT(i,m)
          smfuncveg_mo_g(i) = smfuncveg_mo_g(i) + smfuncveg_mo_t(i,m) * FAREROT(i,m)
          bterm_mo_g(i) = bterm_mo_g(i) + bterm_mo_t(i,m) * FAREROT(i,m)
          lterm_mo_g(i) = lterm_mo_g(i) + lterm_mo_t(i,m) * FAREROT(i,m)
          mterm_mo_g(i) = mterm_mo_g(i) + mterm_mo_t(i,m) * FAREROT(i,m)
          litrfall_mo_g(i) = litrfall_mo_g(i) + litrfall_mo_t(i,m) * FAREROT(i,m)
          humiftrs_mo_g(i) = humiftrs_mo_g(i) + humiftrs_mo_t(i,m) * FAREROT(i,m)

          fProductDecomp_mo_g(i) = fProductDecomp_mo_g(i) + fProductDecomp_mo_t(i,m) * FAREROT(i,m)

        end do ! loop 900

        imonth = nt

        ! Prepare the timestamp for this month

        ! Take one day off so it is the last day of the month rather than the first day of the next month.
        timeStamp(1) = consecDays + monthend(imonth + 1) - 1

        call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo_g' ,timeStamp,'lai', [laimaxg_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo_g',timeStamp,'cVeg',[vgbiomas_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo_g'     ,timeStamp,'npp',[npp_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo_g'     ,timeStamp,'gpp',[gpp_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo_g'     ,timeStamp,'nep',[nep_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo_g'     ,timeStamp,'nbp',[nbp_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo_g',timeStamp,'rh',[hetrores_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo_g' ,timeStamp,'ra',[autores_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilres_mo_g' ,timeStamp,'rSoil',[soilres_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_mo_g' ,timeStamp,'fVegLitter',[litrfall_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrs_mo_g' ,timeStamp,'fLitterSoil',[humiftrs_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_mo_g' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_mo_g' ,timeStamp,'wetlandFrac',[wetfdyn_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_mo_g' ,timeStamp,'soilCH4cons',[ch4soills_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'leafmass_mo_g',timeStamp,'cLeaf',[leafmass_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_mo_g',timeStamp,'cStem',[stemmass_mo_g(i)])

        ! COMBAK PERLAY
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_g',timeStamp,'cLitter',[litrmass_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_g',timeStamp,'cSoil',[soilcmas_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_g'  ,timeStamp,'rhLitter',[litres_mo_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_g',timeStamp,'rhSoil',[soilcres_mo_g(i)])
        ! Make the bulk litter and soil C pool and respiration temporary variables:
        ! allocate(bulkLitterCarbon(1))
        ! allocate(bulkSoilCarbon(1))
        ! allocate(bulkLitterResp(1))
        ! allocate(bulkSoilResp(1))
        ! bulkLitterCarbon(1) = sum(litrmass_mo_g(i,:))
        ! bulkSoilCarbon(1) = sum(soilcmas_mo_g(i,:))
        ! bulkLitterResp(1) = sum(litres_mo_g(i,:))
        ! bulkSoilResp(1) = sum(soilcres_mo_g(i,:))
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_g',timeStamp,'cLitter',[bulkLitterCarbon])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_g',timeStamp,'cSoil',[bulkSoilCarbon])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_g'  ,timeStamp,'rhLitter',[bulkLitterResp])
        ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_g',timeStamp,'rhSoil',[bulkSoilResp])
        ! deallocate(bulkLitterCarbon)
        ! deallocate(bulkSoilCarbon)
        ! deallocate(bulkLitterResp)
        ! deallocate(bulkSoilResp)
        ! COMBAK PERLAY

        if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_mo_g' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfpres_mo_g' ,timeStamp,'wetlandFracPresc',[wetfpres_mo_g(i)])
        end if

        do m = 1,nmtest
          sumfare = 0.0
          do j = 1,icc
            sumfare = sumfare + fcancmxrow(i,m,j)
          end do ! j
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_mo_g' ,timeStamp,'landCoverFrac',[fcancmxrow(i,m,1:icc),1. - sumfare])
        end do ! m

        if (dofire) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_mo_g' ,timeStamp,'fFireCH4',[emit_ch4_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co_mo_g' ,timeStamp,'fFireCO',[emit_co_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo_g' ,timeStamp,'fFire',[emit_co2_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_bc_mo_g' ,timeStamp,'fFireBC',[emit_bc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_oc_mo_g' ,timeStamp,'fFireOC',[emit_oc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_nox_mo_g' ,timeStamp,'fFireNOX',[emit_nox_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_nmhc_mo_g' ,timeStamp,'fFireNMHC',[emit_nmhc_mo_g(i)])          
          call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo_g' ,timeStamp,'burntFractionAll',[burnfrac_mo_g(i)])
        end if
        if (lnduseon) then
          call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_mo_g' ,timeStamp,'fDeforestToAtmos',[luc_emc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_mo_g' ,timeStamp,'fDeforestToLitter',[lucltrin_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_mo_g' ,timeStamp,'fDeforestToSoil',[lucsocin_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_mo_g' ,timeStamp,'fDeforestTotal', &
                             [lucsocin_mo_g(i) + lucltrin_mo_g(i) + luc_emc_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fProductDecomp_mo_g' ,timeStamp,'fProductDecomp',[fProductDecomp_mo_g(i)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'cProduct_mo_g' ,timeStamp,'cProduct',[cProduct_mo_g(i)])
        end if

        if (PFTCompetition) then
          pftExist = 0.0
          do j = 1,icc
            if (pftexistrow(i,1,j)) pftExist(j) = 1.0
          end do
          call writeOutput1D(lonLocalIndex,latLocalIndex,'pftexistrow_mo_g' ,timeStamp,'landCoverExist',[pftExist]) ! flag only set up for one tile !
        end if

        if (doperpftoutput) then
          if (nmtest > 1) then
            print * ,'Per PFT and per tile outputs together not implemented yet'
          else
            m = 1
            call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo' ,timeStamp,'lai', [laimaxg_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'leafmass_mo'  ,timeStamp,'cLeaf',[leafmass_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo',timeStamp,'cVeg',[vgbiomas_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo'     ,timeStamp,'npp',[npp_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo'     ,timeStamp,'gpp',[gpp_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo'     ,timeStamp,'nep',[nep_mo(i,m,:)])
            ! NOTE: This NBP does not include LUC product pool contributions since they are
            ! not per PFT but rather per tile
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo'     ,timeStamp,'nbp',[nbp_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo',timeStamp,'rh',[hetrores_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo' ,timeStamp,'ra',[autores_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilres_mo' ,timeStamp,'rSoil',[soilres_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfallveg_mo' ,timeStamp,'fVegLitter',[litrfallveg_mo(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrsveg_mo' ,timeStamp,'fLitterSoil',[humiftrsveg_mo(i,m,1:iccp1)])

            ! COMBAK PERLAY
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo',timeStamp,'cLitter',[litrmass_mo(i,m,1:iccp1)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo',timeStamp,'cSoil',[soilcmas_mo(i,m,1:iccp1)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo'  ,timeStamp,'rhLitter',[litres_mo(i,m,1:iccp1)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo',timeStamp,'rhSoil',[soilcres_mo(i,m,1:iccp1)])
            ! Make the bulk litter and soil C pool and respiration temporary variables:
            ! allocate(bulkLitterCarbon(iccp1))
            ! allocate(bulkSoilCarbon(iccp1))
            ! allocate(bulkLitterResp(iccp1))
            ! allocate(bulkSoilResp(iccp1))
            ! do k = 1,iccp1
            !   bulkLitterCarbon(k) = sum(litrmass_mo(i,m,k,:))
            !   bulkSoilCarbon(k) = sum(soilcmas_mo(i,m,k,:))
            !   bulkLitterResp(k) = sum(litres_mo(i,m,k,:))
            !   bulkSoilResp(k) = sum(soilcres_mo(i,m,k,:))
            ! end do
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo',timeStamp,'cLitter',[bulkLitterCarbon])
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo',timeStamp,'cSoil',[bulkSoilCarbon])
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo'  ,timeStamp,'rhLitter',[bulkLitterResp])
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo',timeStamp,'rhSoil',[bulkSoilResp])
            ! deallocate(bulkLitterCarbon)
            ! deallocate(bulkSoilCarbon)
            ! deallocate(bulkLitterResp)
            ! deallocate(bulkSoilResp)
            ! COMBAK PERLAY

            if (dofire .or. lnduseon) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_mo' ,timeStamp,'fFireCH4',[emit_ch4_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo' ,timeStamp,'fFire',[emit_co2_mo(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo' ,timeStamp,'burntFractionAll',[burnfrac_mo(i,m,:)])
            end if

            if (useTracer > 0) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerGLeafMass' ,timeStamp,'cLeafTracer',[tracerGLeafMassrot(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerStemMass' ,timeStamp,'cStemTracer',[tracerStemMassrot(i,m,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerRootMass' ,timeStamp,'cRootTracer',[tracerRootMassrot(i,m,:)])
              ! call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerLitrMass' ,timeStamp,'cLitterTracer',[tracerLitrMassrot(i,m,:)])
              ! call writeOutput1D(lonLocalIndex,latLocalIndex,'tracerSoilCMass' ,timeStamp,'cSoilTracer',[tracerSoilCMassrot(i,m,:)])
            end if

          end if
        end if

        if (dopertileoutput) then
          if (nmtest > 1) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo_t' ,timeStamp,'lai', [laimaxg_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_mo_t',timeStamp,'cVeg',[vgbiomas_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_mo_t'     ,timeStamp,'npp',[npp_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_mo_t'     ,timeStamp,'gpp',[gpp_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_mo_t'     ,timeStamp,'nep',[nep_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_mo_t'     ,timeStamp,'nbp',[nbp_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_mo_t',timeStamp,'rh',[hetrores_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_mo_t' ,timeStamp,'ra',[autores_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilres_mo_t' ,timeStamp,'rSoil',[soilres_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrfall_mo_t' ,timeStamp,'fVegLitter',[litrfall_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'humiftrs_mo_t' ,timeStamp,'fLitterSoil',[humiftrs_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'leafmass_mo_t',timeStamp,'cLeaf',[leafmass_mo_t(i,:)])

            ! COMBAK PERLAY

            call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_t',timeStamp,'cLitter',[litrmass_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_t',timeStamp,'cSoil',[soilcmas_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_t'  ,timeStamp,'rhLitter',[litres_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_t',timeStamp,'rhSoil',[soilcres_mo_t(i,:)])
            ! Make the bulk litter and soil C pool and respiration temporary variables:
            ! allocate(bulkLitterCarbon(nmtest))
            ! allocate(bulkSoilCarbon(nmtest))
            ! allocate(bulkLitterResp(nmtest))
            ! allocate(bulkSoilResp(nmtest))
            !
            ! do m = 1,nmtest
            !   bulkLitterCarbon(m) = sum(litrmass_mo_t(i,m,:))
            !   bulkSoilCarbon(m) = sum(soilcmas_mo_t(i,m,:))
            !   bulkLitterResp(m) = sum(litres_mo_t(i,m,:))
            !   bulkSoilResp(m) = sum(soilcres_mo_t(i,m,:))
            ! end do
            !
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_mo_t',timeStamp,'cLitter',[bulkLitterCarbon])
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_mo_t',timeStamp,'cSoil',[bulkSoilCarbon])
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_mo_t'  ,timeStamp,'rhLitter',[bulkLitterResp])
            ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_mo_t',timeStamp,'rhSoil',[bulkSoilResp])
            !
            ! deallocate(bulkLitterCarbon)
            ! deallocate(bulkSoilCarbon)
            ! deallocate(bulkLitterResp)
            ! deallocate(bulkSoilResp)
            ! COMBAK PERLAY

            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_mo_t' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_mo_t' ,timeStamp,'wetlandFrac',[wetfdyn_mo_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_mo_t' ,timeStamp,'soilCH4cons',[ch4soills_mo_t(i,:)])

            if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfpres_mo_t' ,timeStamp,'wetlandFracPresc',[wetfpres_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_mo_t' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_mo_t(i,:)])
            end if

            if (dofire) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_mo_t' ,timeStamp,'fFireCH4',[emit_ch4_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_mo_t' ,timeStamp,'fFire',[emit_co2_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_mo_t' ,timeStamp,'burntFractionAll',[burnfrac_mo_t(i,:)])
            end if

            if (lnduseon) then
              call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_mo_t' ,timeStamp,'fDeforestToAtmos',[luc_emc_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_mo_t' ,timeStamp,'fDeforestToLitter',[lucltrin_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_mo_t' ,timeStamp,'fDeforestToSoil',[lucsocin_mo_t(i,:)])
              call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_mo_t' ,timeStamp,'fDeforestTotal', &
                                 [lucsocin_mo_t(i,:) + lucltrin_mo_t(i,:) + luc_emc_mo_t(i,:)])
            end if
          end if
        end if

        !> Reset all end of month accumulated arrays
        call resetMonthEnd(nltest ,nmtest)

      end if ! end of month
    end do ! loop 865 ! nmon

  end subroutine ctem_monthly_aw
  !! @}

  !==============================================================================================================

  !> \ingroup prepareoutputs_ctem_annual_aw
  !> @{
  !> Accumulate and write out the annual biogeochemical (CTEM) outputs. These are kept in pointer structures as
  !! this subroutine is called daily and we increment the daily values to produce annual values.
  !! The pointer to the annual data structures (in ctemStateVars) keeps the data between calls.
  !! @author J. Melton

  subroutine ctem_annual_aw (lonLocalIndex, latLocalIndex, iday, realyr, nltest, nmtest, lastDOY)

    use classStateVars, only : class_rot
    use ctemStateVars,  only : ctem_tile_yr, vrot, ctem_grd_yr, c_switch, ctem_yr, &
                               resetYearEnd, vgat
    use classicParams,  only : icc, iccp1, seed, iccp2, ignd, nlat, nmos
    use outputManager,  only : writeOutput1D, consecDays

    implicit none

    ! arguments
    integer, intent(in) :: lonLocalIndex, latLocalIndex
    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest
    integer, intent(in) :: iday
    integer, intent(in) :: realyr
    integer, intent(in) :: lastDOY

    ! pointers

    logical, pointer :: dofire
    logical, pointer :: lnduseon
    logical, pointer :: PFTCompetition
    logical, pointer :: doperpftoutput
    logical, pointer :: dopertileoutput
    logical, pointer :: transientOBSWETF
    integer, pointer :: fixedYearOBSWETF

    real, pointer, dimension(:,:) :: FAREROT !< Fractional coverage of mosaic tile on modelled area
    real, pointer, dimension(:,:,:) :: laimaxg_yr
    real, pointer, dimension(:,:,:) :: stemmass_yr
    real, pointer, dimension(:,:,:) :: rootmass_yr
    real, pointer, dimension(:,:,:) :: npp_yr
    real, pointer, dimension(:,:,:) :: gpp_yr
    real, pointer, dimension(:,:,:) :: vgbiomas_yr
    real, pointer, dimension(:,:,:) :: autores_yr
    real, pointer, dimension(:,:,:) :: totcmass_yr
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litrmass_yr
    real, pointer, dimension(:,:,:) :: soilcmas_yr
    real, pointer, dimension(:,:,:) :: litres_yr
    real, pointer, dimension(:,:,:) :: soilcres_yr
    ! real, pointer, dimension(:,:,:,:) :: litrmass_yr
    ! real, pointer, dimension(:,:,:,:) :: soilcmas_yr
    ! real, pointer, dimension(:,:,:,:) :: litres_yr
    ! real, pointer, dimension(:,:,:,:) :: soilcres_yr
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: nep_yr
    real, pointer, dimension(:,:,:) :: hetrores_yr
    real, pointer, dimension(:,:,:) :: nbp_yr
    real, pointer, dimension(:,:,:) :: emit_co2_yr
    real, pointer, dimension(:,:,:) :: emit_co_yr
    real, pointer, dimension(:,:,:) :: emit_ch4_yr
    real, pointer, dimension(:,:,:) :: emit_nmhc_yr
    real, pointer, dimension(:,:,:) :: emit_h2_yr
    real, pointer, dimension(:,:,:) :: emit_nox_yr
    real, pointer, dimension(:,:,:) :: emit_n2o_yr
    real, pointer, dimension(:,:,:) :: emit_pm25_yr
    real, pointer, dimension(:,:,:) :: emit_tpm_yr
    real, pointer, dimension(:,:,:) :: emit_tc_yr
    real, pointer, dimension(:,:,:) :: emit_oc_yr
    real, pointer, dimension(:,:,:) :: emit_bc_yr
    real, pointer, dimension(:,:,:) :: bterm_yr
    real, pointer, dimension(:,:,:) :: mterm_yr
    real, pointer, dimension(:,:,:) :: smfuncveg_yr
    real, pointer, dimension(:,:,:) :: burnfrac_yr
    real, pointer, dimension(:,:,:) :: veghght_yr
    real, pointer, dimension(:,:) :: laimaxg_yr_t
    real, pointer, dimension(:,:) :: stemmass_yr_t
    real, pointer, dimension(:,:) :: rootmass_yr_t
    real, pointer, dimension(:,:) :: npp_yr_t
    real, pointer, dimension(:,:) :: gpp_yr_t
    real, pointer, dimension(:,:) :: vgbiomas_yr_t
    real, pointer, dimension(:,:) :: autores_yr_t
    real, pointer, dimension(:,:) :: totcmass_yr_t
    ! COMBAK PERLAY
    real, pointer, dimension(:,:) :: litrmass_yr_t
    real, pointer, dimension(:,:) :: soilcmas_yr_t
    real, pointer, dimension(:,:) :: litres_yr_t
    real, pointer, dimension(:,:) :: soilcres_yr_t
    ! real, pointer, dimension(:,:,:) :: litrmass_yr_t
    ! real, pointer, dimension(:,:,:) :: soilcmas_yr_t
    ! real, pointer, dimension(:,:,:) :: litres_yr_t
    ! real, pointer, dimension(:,:,:) :: soilcres_yr_t
    ! COMBAK PERLAY
    real, pointer, dimension(:,:) :: nep_yr_t
    real, pointer, dimension(:,:) :: hetrores_yr_t
    real, pointer, dimension(:,:) :: nbp_yr_t
    real, pointer, dimension(:,:) :: emit_co2_yr_t
    real, pointer, dimension(:,:) :: emit_co_yr_t
    real, pointer, dimension(:,:) :: emit_ch4_yr_t
    real, pointer, dimension(:,:) :: emit_nmhc_yr_t
    real, pointer, dimension(:,:) :: emit_h2_yr_t
    real, pointer, dimension(:,:) :: emit_nox_yr_t
    real, pointer, dimension(:,:) :: emit_n2o_yr_t
    real, pointer, dimension(:,:) :: emit_pm25_yr_t
    real, pointer, dimension(:,:) :: emit_tpm_yr_t
    real, pointer, dimension(:,:) :: emit_tc_yr_t
    real, pointer, dimension(:,:) :: emit_oc_yr_t
    real, pointer, dimension(:,:) :: emit_bc_yr_t
    real, pointer, dimension(:,:) :: burnfrac_yr_t
    real, pointer, dimension(:,:) :: smfuncveg_yr_t
    real, pointer, dimension(:,:) :: bterm_yr_t
    real, pointer, dimension(:,:) :: luc_emc_yr_t
    real, pointer, dimension(:,:) :: lterm_yr_t
    real, pointer, dimension(:,:) :: lucsocin_yr_t
    real, pointer, dimension(:,:) :: mterm_yr_t
    real, pointer, dimension(:,:) :: lucltrin_yr_t
    real, pointer, dimension(:,:) :: ch4WetSpec_yr_t
    real, pointer, dimension(:,:) :: wetfdyn_yr_t
    real, pointer, dimension(:,:) :: ch4WetDyn_yr_t
    real, pointer, dimension(:,:) :: ch4soills_yr_t
    real, pointer, dimension(:,:) :: veghght_yr_t
    real, pointer, dimension(:,:) :: peatdep_yr_t
    real, pointer, dimension(:,:) :: fProductDecomp_yr_t    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$
    logical, pointer, dimension(:,:,:) :: pftexistrow
    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:)   :: nbprow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    ! real, pointer, dimension(:,:,:,:) :: litrmassrow
    ! real, pointer, dimension(:,:,:,:) :: soilcmasrow
    ! real, pointer, dimension(:,:,:,:) :: litresvegrow
    ! real, pointer, dimension(:,:,:,:) :: soilcresvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: rmlvegaccrow
    real, pointer, dimension(:,:,:) :: rmsvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    real, pointer, dimension(:,:,:) :: ailcgrow
    real, pointer, dimension(:,:,:) :: emit_co2row
    real, pointer, dimension(:,:,:) :: emit_corow
    real, pointer, dimension(:,:,:) :: emit_ch4row
    real, pointer, dimension(:,:,:) :: emit_nmhcrow
    real, pointer, dimension(:,:,:) :: emit_h2row
    real, pointer, dimension(:,:,:) :: emit_noxrow
    real, pointer, dimension(:,:,:) :: emit_n2orow
    real, pointer, dimension(:,:,:) :: emit_pm25row
    real, pointer, dimension(:,:,:) :: emit_tpmrow
    real, pointer, dimension(:,:,:) :: emit_tcrow
    real, pointer, dimension(:,:,:) :: emit_ocrow
    real, pointer, dimension(:,:,:) :: emit_bcrow
    real, pointer, dimension(:,:) :: burnfracrow
    real, pointer, dimension(:,:,:) :: burnvegfrow
    real, pointer, dimension(:,:,:) :: smfuncvegrow
    real, pointer, dimension(:,:,:) :: btermrow
    real, pointer, dimension(:,:) :: ltermrow
    real, pointer, dimension(:,:,:) :: mtermrow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow
    real, pointer, dimension(:,:) :: ch4WetSpecrow
    real, pointer, dimension(:,:) :: wetfdynrow
    real, pointer, dimension(:,:) :: ch4WetDynrow
    real, pointer, dimension(:,:) :: ch4soillsrow
    real, pointer, dimension(:,:,:) :: vgbiomas_vegrow
    real, pointer, dimension(:,:,:) :: stemmassrow
    real, pointer, dimension(:,:,:) :: rootmassrow
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:,:) :: veghghtrow
    real, pointer, dimension(:,:) :: peatdeprow
    real, pointer, dimension(:) :: laimaxg_yr_g
    real, pointer, dimension(:) :: stemmass_yr_g
    real, pointer, dimension(:) :: rootmass_yr_g
    ! COMBAK PERLAY
    real, pointer, dimension(:) :: litrmass_yr_g
    real, pointer, dimension(:) :: soilcmas_yr_g
    real, pointer, dimension(:) :: litres_yr_g
    real, pointer, dimension(:) :: soilcres_yr_g
    ! real, pointer, dimension(:,:) :: litrmass_yr_g
    ! real, pointer, dimension(:,:) :: soilcmas_yr_g
    ! real, pointer, dimension(:,:) :: litres_yr_g
    ! real, pointer, dimension(:,:) :: soilcres_yr_g
    ! COMBAK PERLAY
    real, pointer, dimension(:) :: npp_yr_g
    real, pointer, dimension(:) :: gpp_yr_g
    real, pointer, dimension(:) :: nep_yr_g
    real, pointer, dimension(:) :: nbp_yr_g
    real, pointer, dimension(:) :: hetrores_yr_g
    real, pointer, dimension(:) :: autores_yr_g
    real, pointer, dimension(:) :: vgbiomas_yr_g
    real, pointer, dimension(:) :: totcmass_yr_g
    real, pointer, dimension(:) :: emit_co2_yr_g
    real, pointer, dimension(:) :: emit_co_yr_g
    real, pointer, dimension(:) :: emit_ch4_yr_g
    real, pointer, dimension(:) :: emit_nmhc_yr_g
    real, pointer, dimension(:) :: emit_h2_yr_g
    real, pointer, dimension(:) :: emit_nox_yr_g
    real, pointer, dimension(:) :: emit_n2o_yr_g
    real, pointer, dimension(:) :: emit_pm25_yr_g
    real, pointer, dimension(:) :: emit_tpm_yr_g
    real, pointer, dimension(:) :: emit_tc_yr_g
    real, pointer, dimension(:) :: emit_oc_yr_g
    real, pointer, dimension(:) :: emit_bc_yr_g
    real, pointer, dimension(:) :: smfuncveg_yr_g
    real, pointer, dimension(:) :: luc_emc_yr_g
    real, pointer, dimension(:) :: lucltrin_yr_g
    real, pointer, dimension(:) :: lucsocin_yr_g
    real, pointer, dimension(:) :: burnfrac_yr_g
    real, pointer, dimension(:) :: bterm_yr_g
    real, pointer, dimension(:) :: lterm_yr_g
    real, pointer, dimension(:) :: mterm_yr_g
    real, pointer, dimension(:) :: ch4WetSpec_yr_g
    real, pointer, dimension(:) :: wetfdyn_yr_g
    real, pointer, dimension(:) :: ch4WetDyn_yr_g
    real, pointer, dimension(:) :: ch4soills_yr_g
    real, pointer, dimension(:) :: veghght_yr_g
    real, pointer, dimension(:) :: peatdep_yr_g
    real, pointer, dimension(:) :: cProduct_yr_g          !< Carbon in the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2}]\f$
    real, pointer, dimension(:) :: fProductDecomp_yr_g    !< Respiration of carbon from the LUC product pools (litter and soil C iccp2 position) \f$[kg C m^{-2} s^{-1}]\f$

    ! local
    integer :: i, m, j, nt, k
    real :: barefrac
    real :: sumfare
    real, dimension(1) :: timeStamp
    real, dimension(icc) :: pftExist
    real, dimension(icc) :: fcancmxNoSeed
    real, dimension(:), allocatable :: bulkLitterCarbon !< Temporary variable used to produce the bulk soil litter carbon quantity for output \f$[kg C m^{-2}]\f$
    real, dimension(:), allocatable :: bulkSoilCarbon   !< Temporary variable used to produce the bulk soil carbon quantity for output \f$[kg C m^{-2}]\f$
    real, dimension(:), allocatable :: bulkLitterResp   !< Temporary variable used to produce the bulk litter respiration quantity for output \f$[kg C m^{-2} s^{-1}]\f$
    real, dimension(:), allocatable :: bulkSoilResp     !< Temporary variable used to produce the bulk soil carbon respiration quantity for output \f$[kg C m^{-2} s^{-1}]\f$

    real :: oneOverDPY
    
    real :: coverLow, coverHigh

    ! point pointers

    dofire                => c_switch%dofire
    lnduseon              => c_switch%lnduseon
    PFTCompetition        => c_switch%PFTCompetition
    doperpftoutput        => c_switch%doperpftoutput
    dopertileoutput       => c_switch%dopertileoutput
    transientOBSWETF      => c_switch%transientOBSWETF
    fixedYearOBSWETF      => c_switch%fixedYearOBSWETF

    FAREROT => class_rot%FAREROT

    laimaxg_yr          =>ctem_yr%laimaxg_yr
    stemmass_yr         =>ctem_yr%stemmass_yr
    rootmass_yr         =>ctem_yr%rootmass_yr
    npp_yr              =>ctem_yr%npp_yr
    gpp_yr              =>ctem_yr%gpp_yr
    vgbiomas_yr         =>ctem_yr%vgbiomas_yr
    autores_yr          =>ctem_yr%autores_yr
    totcmass_yr         =>ctem_yr%totcmass_yr
    litrmass_yr         =>ctem_yr%litrmass_yr
    soilcmas_yr         =>ctem_yr%soilcmas_yr
    nep_yr              =>ctem_yr%nep_yr
    litres_yr           =>ctem_yr%litres_yr
    soilcres_yr         =>ctem_yr%soilcres_yr
    hetrores_yr         =>ctem_yr%hetrores_yr
    nbp_yr              =>ctem_yr%nbp_yr
    emit_co2_yr         =>ctem_yr%emit_co2_yr
    emit_co_yr          =>ctem_yr%emit_co_yr
    emit_ch4_yr         =>ctem_yr%emit_ch4_yr
    emit_nmhc_yr        =>ctem_yr%emit_nmhc_yr
    emit_h2_yr          =>ctem_yr%emit_h2_yr
    emit_nox_yr         =>ctem_yr%emit_nox_yr
    emit_n2o_yr         =>ctem_yr%emit_n2o_yr
    emit_pm25_yr        =>ctem_yr%emit_pm25_yr
    emit_tpm_yr         =>ctem_yr%emit_tpm_yr
    emit_tc_yr          =>ctem_yr%emit_tc_yr
    emit_oc_yr          =>ctem_yr%emit_oc_yr
    emit_bc_yr          =>ctem_yr%emit_bc_yr
    bterm_yr            =>ctem_yr%bterm_yr
    mterm_yr            =>ctem_yr%mterm_yr
    burnfrac_yr         =>ctem_yr%burnfrac_yr
    smfuncveg_yr        =>ctem_yr%smfuncveg_yr
    veghght_yr          =>ctem_yr%veghght_yr

    laimaxg_yr_t          =>ctem_tile_yr%laimaxg_yr_t
    stemmass_yr_t         =>ctem_tile_yr%stemmass_yr_t
    rootmass_yr_t         =>ctem_tile_yr%rootmass_yr_t
    npp_yr_t              =>ctem_tile_yr%npp_yr_t
    gpp_yr_t              =>ctem_tile_yr%gpp_yr_t
    vgbiomas_yr_t         =>ctem_tile_yr%vgbiomas_yr_t
    autores_yr_t          =>ctem_tile_yr%autores_yr_t
    totcmass_yr_t         =>ctem_tile_yr%totcmass_yr_t
    litrmass_yr_t         =>ctem_tile_yr%litrmass_yr_t
    soilcmas_yr_t         =>ctem_tile_yr%soilcmas_yr_t
    nep_yr_t              =>ctem_tile_yr%nep_yr_t
    litres_yr_t           =>ctem_tile_yr%litres_yr_t
    soilcres_yr_t         =>ctem_tile_yr%soilcres_yr_t
    hetrores_yr_t         =>ctem_tile_yr%hetrores_yr_t
    nbp_yr_t              =>ctem_tile_yr%nbp_yr_t
    emit_co2_yr_t         =>ctem_tile_yr%emit_co2_yr_t
    emit_co_yr_t          =>ctem_tile_yr%emit_co_yr_t
    emit_ch4_yr_t         =>ctem_tile_yr%emit_ch4_yr_t
    emit_nmhc_yr_t        =>ctem_tile_yr%emit_nmhc_yr_t
    emit_h2_yr_t          =>ctem_tile_yr%emit_h2_yr_t
    emit_nox_yr_t         =>ctem_tile_yr%emit_nox_yr_t
    emit_n2o_yr_t         =>ctem_tile_yr%emit_n2o_yr_t
    emit_pm25_yr_t        =>ctem_tile_yr%emit_pm25_yr_t
    emit_tpm_yr_t         =>ctem_tile_yr%emit_tpm_yr_t
    emit_tc_yr_t          =>ctem_tile_yr%emit_tc_yr_t
    emit_oc_yr_t          =>ctem_tile_yr%emit_oc_yr_t
    emit_bc_yr_t          =>ctem_tile_yr%emit_bc_yr_t
    burnfrac_yr_t         =>ctem_tile_yr%burnfrac_yr_t
    smfuncveg_yr_t        =>ctem_tile_yr%smfuncveg_yr_t
    bterm_yr_t            =>ctem_tile_yr%bterm_yr_t
    luc_emc_yr_t          =>ctem_tile_yr%luc_emc_yr_t
    lterm_yr_t            =>ctem_tile_yr%lterm_yr_t
    lucsocin_yr_t         =>ctem_tile_yr%lucsocin_yr_t
    mterm_yr_t            =>ctem_tile_yr%mterm_yr_t
    lucltrin_yr_t         =>ctem_tile_yr%lucltrin_yr_t
    ch4WetSpec_yr_t       =>ctem_tile_yr%ch4WetSpec_yr_t
    wetfdyn_yr_t          =>ctem_tile_yr%wetfdyn_yr_t
    ch4WetDyn_yr_t        =>ctem_tile_yr%ch4WetDyn_yr_t
    ch4soills_yr_t        =>ctem_tile_yr%ch4soills_yr_t
    veghght_yr_t          =>ctem_tile_yr%veghght_yr_t
    peatdep_yr_t          =>ctem_tile_yr%peatdep_yr_t
    fProductDecomp_yr_t   =>ctem_tile_yr%fProductDecomp_yr_t

    pftexistrow       => vrot%pftexist
    gppvegrow         => vrot%gppveg
    nepvegrow         => vrot%nepveg
    nbpvegrow         => vrot%nbpveg
    nbprow            => vrot%nbp
    nppvegrow         => vrot%nppveg
    hetroresvegrow    => vrot%hetroresveg
    autoresvegrow     => vrot%autoresveg
    litresvegrow      => vrot%litresveg
    soilcresvegrow    => vrot%soilcresveg
    rmlvegaccrow      => vrot%rmlvegacc
    rmsvegrow         => vrot%rmsveg
    rmrvegrow         => vrot%rmrveg
    rgvegrow          => vrot%rgveg
    ailcgrow          => vrot%ailcg
    emit_co2row       => vrot%emit_co2
    emit_corow        => vrot%emit_co
    emit_ch4row       => vrot%emit_ch4
    emit_nmhcrow      => vrot%emit_nmhc
    emit_h2row        => vrot%emit_h2
    emit_noxrow       => vrot%emit_nox
    emit_n2orow       => vrot%emit_n2o
    emit_pm25row      => vrot%emit_pm25
    emit_tpmrow       => vrot%emit_tpm
    emit_tcrow        => vrot%emit_tc
    emit_ocrow        => vrot%emit_oc
    emit_bcrow        => vrot%emit_bc
    burnfracrow       => vrot%burnfrac
    burnvegfrow       => vrot%burnvegf
    smfuncvegrow      => vrot%smfuncveg
    btermrow          => vrot%bterm
    ltermrow          => vrot%lterm
    mtermrow          => vrot%mterm
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin
    ch4WetSpecrow     => vrot%ch4WetSpec
    wetfdynrow        => vrot%wetfdyn
    ch4WetDynrow      => vrot%ch4WetDyn
    ch4soillsrow      => vrot%ch4_soills
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    vgbiomas_vegrow   => vrot%vgbiomas_veg
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    fcancmxrow        => vrot%fcancmx
    veghghtrow        => vrot%veghght

    peatdeprow            => vrot%peatdep

    laimaxg_yr_g          =>ctem_grd_yr%laimaxg_yr_g
    stemmass_yr_g         =>ctem_grd_yr%stemmass_yr_g
    rootmass_yr_g         =>ctem_grd_yr%rootmass_yr_g
    litrmass_yr_g         =>ctem_grd_yr%litrmass_yr_g
    soilcmas_yr_g         =>ctem_grd_yr%soilcmas_yr_g
    npp_yr_g              =>ctem_grd_yr%npp_yr_g
    gpp_yr_g              =>ctem_grd_yr%gpp_yr_g
    nep_yr_g              =>ctem_grd_yr%nep_yr_g
    nbp_yr_g              =>ctem_grd_yr%nbp_yr_g
    hetrores_yr_g         =>ctem_grd_yr%hetrores_yr_g
    autores_yr_g          =>ctem_grd_yr%autores_yr_g
    litres_yr_g           =>ctem_grd_yr%litres_yr_g
    soilcres_yr_g         =>ctem_grd_yr%soilcres_yr_g
    vgbiomas_yr_g         =>ctem_grd_yr%vgbiomas_yr_g
    totcmass_yr_g         =>ctem_grd_yr%totcmass_yr_g
    emit_co2_yr_g         =>ctem_grd_yr%emit_co2_yr_g
    emit_co_yr_g          =>ctem_grd_yr%emit_co_yr_g
    emit_ch4_yr_g         =>ctem_grd_yr%emit_ch4_yr_g
    emit_nmhc_yr_g        =>ctem_grd_yr%emit_nmhc_yr_g
    emit_h2_yr_g          =>ctem_grd_yr%emit_h2_yr_g
    emit_nox_yr_g         =>ctem_grd_yr%emit_nox_yr_g
    emit_n2o_yr_g         =>ctem_grd_yr%emit_n2o_yr_g
    emit_pm25_yr_g        =>ctem_grd_yr%emit_pm25_yr_g
    emit_tpm_yr_g         =>ctem_grd_yr%emit_tpm_yr_g
    emit_tc_yr_g          =>ctem_grd_yr%emit_tc_yr_g
    emit_oc_yr_g          =>ctem_grd_yr%emit_oc_yr_g
    emit_bc_yr_g          =>ctem_grd_yr%emit_bc_yr_g
    smfuncveg_yr_g        =>ctem_grd_yr%smfuncveg_yr_g
    luc_emc_yr_g          =>ctem_grd_yr%luc_emc_yr_g
    lucltrin_yr_g         =>ctem_grd_yr%lucltrin_yr_g
    lucsocin_yr_g         =>ctem_grd_yr%lucsocin_yr_g
    burnfrac_yr_g         =>ctem_grd_yr%burnfrac_yr_g
    bterm_yr_g            =>ctem_grd_yr%bterm_yr_g
    lterm_yr_g            =>ctem_grd_yr%lterm_yr_g
    mterm_yr_g            =>ctem_grd_yr%mterm_yr_g
    ch4WetSpec_yr_g       =>ctem_grd_yr%ch4WetSpec_yr_g
    wetfdyn_yr_g          =>ctem_grd_yr%wetfdyn_yr_g
    ch4WetDyn_yr_g        =>ctem_grd_yr%ch4WetDyn_yr_g
    ch4soills_yr_g        =>ctem_grd_yr%ch4soills_yr_g
    veghght_yr_g          =>ctem_grd_yr%veghght_yr_g
    peatdep_yr_g          =>ctem_grd_yr%peatdep_yr_g
    cProduct_yr_g         =>ctem_grd_yr%cProduct_yr_g
    fProductDecomp_yr_g   =>ctem_grd_yr%fProductDecomp_yr_g

    !------------
    !> Accumulate yearly outputs
    oneOverDPY = 1./real(lastDOY)
    i = 1 ! offline nlat is always 1 so this array position is always 1.
    do m = 1,nmtest
      do j = 1,icc

        !> Accumulate the variables at the per PFT level

        if (ailcgrow(i,m,j) > laimaxg_yr(i,m,j)) then
          laimaxg_yr(i,m,j) = ailcgrow(i,m,j)
        end if

        npp_yr(i,m,j) = npp_yr(i,m,j) + nppvegrow(i,m,j) * oneOverDPY
        gpp_yr(i,m,j) = gpp_yr(i,m,j) + gppvegrow(i,m,j) * oneOverDPY
        nep_yr(i,m,j) = nep_yr(i,m,j) + nepvegrow(i,m,j) * oneOverDPY
        ! NOTE: This NBP does not include LUC product pool contributions since they are
        ! not per PFT but rather per tile
        nbp_yr(i,m,j) = nbp_yr(i,m,j) + nbpvegrow(i,m,j) * oneOverDPY
        emit_co2_yr(i,m,j) = emit_co2_yr(i,m,j) + emit_co2row(i,m,j) * oneOverDPY
        emit_co_yr(i,m,j) = emit_co_yr(i,m,j) + emit_corow(i,m,j) * oneOverDPY
        emit_ch4_yr(i,m,j) = emit_ch4_yr(i,m,j) + emit_ch4row(i,m,j) * oneOverDPY
        emit_nmhc_yr(i,m,j) = emit_nmhc_yr(i,m,j) + emit_nmhcrow(i,m,j) * oneOverDPY
        emit_h2_yr(i,m,j) = emit_h2_yr(i,m,j) + emit_h2row(i,m,j) * oneOverDPY
        emit_nox_yr(i,m,j) = emit_nox_yr(i,m,j) + emit_noxrow(i,m,j) * oneOverDPY
        emit_n2o_yr(i,m,j) = emit_n2o_yr(i,m,j) + emit_n2orow(i,m,j) * oneOverDPY
        emit_pm25_yr(i,m,j) = emit_pm25_yr(i,m,j) + emit_pm25row(i,m,j) * oneOverDPY
        emit_tpm_yr(i,m,j) = emit_tpm_yr(i,m,j) + emit_tpmrow(i,m,j) * oneOverDPY
        emit_tc_yr(i,m,j) = emit_tc_yr(i,m,j) + emit_tcrow(i,m,j) * oneOverDPY
        emit_oc_yr(i,m,j) = emit_oc_yr(i,m,j) + emit_ocrow(i,m,j) * oneOverDPY
        emit_bc_yr(i,m,j) = emit_bc_yr(i,m,j) + emit_bcrow(i,m,j) * oneOverDPY

        bterm_yr(i,m,j) = bterm_yr(i,m,j) + btermrow(i,m,j) * oneOverDPY
        mterm_yr(i,m,j) = mterm_yr(i,m,j) + mtermrow(i,m,j) * oneOverDPY
        smfuncveg_yr(i,m,j) = smfuncveg_yr(i,m,j) + smfuncvegrow(i,m,j) * oneOverDPY
        hetrores_yr(i,m,j) = hetrores_yr(i,m,j) + hetroresvegrow(i,m,j) * oneOverDPY
        autores_yr(i,m,j) = autores_yr(i,m,j) + autoresvegrow(i,m,j) * oneOverDPY
        ! COMBAK PERLAY
        litres_yr(i,m,j) = litres_yr(i,m,j) + litresvegrow(i,m,j) * oneOverDPY
        soilcres_yr(i,m,j) = soilcres_yr(i,m,j) + soilcresvegrow(i,m,j) * oneOverDPY
        ! do k = 1,ignd
        !   litres_yr(i,m,j,k) = litres_yr(i,m,j,k) + litresvegrow(i,m,j,k) * oneOverDPY
        !   soilcres_yr(i,m,j,k) = soilcres_yr(i,m,j,k) + soilcresvegrow(i,m,j,k) * oneOverDPY
        ! end do
        ! COMBAK PERLAY
        ! Let accumulate, not a flux or a mean value.
        burnfrac_yr(i,m,j) = burnfrac_yr(i,m,j) + burnvegfrow(i,m,j)

      end do ! loop 884

      !>   Also do the bare fraction amounts
      hetrores_yr(i,m,iccp1) = hetrores_yr(i,m,iccp1) + hetroresvegrow(i,m,iccp1) * oneOverDPY
      ! COMBAK PERLAY
      litres_yr(i,m,iccp1) = litres_yr(i,m,iccp1) + litresvegrow(i,m,iccp1) * oneOverDPY
      soilcres_yr(i,m,iccp1) = soilcres_yr(i,m,iccp1) &
                               + soilcresvegrow(i,m,iccp1) * oneOverDPY
      ! do k = 1,ignd
      !   litres_yr(i,m,iccp1,k) = litres_yr(i,m,iccp1,k) + litresvegrow(i,m,iccp1,k) * oneOverDPY
      !   soilcres_yr(i,m,iccp1,k) = soilcres_yr(i,m,iccp1,k) &
      !                              + soilcresvegrow(i,m,iccp1,k) * oneOverDPY
      ! end do
      ! COMBAK PERLAY
      nep_yr(i,m,iccp1) = nep_yr(i,m,iccp1) + nepvegrow(i,m,iccp1) * oneOverDPY
      nbp_yr(i,m,iccp1) = nbp_yr(i,m,iccp1) + nbpvegrow(i,m,iccp1) * oneOverDPY

      peatdep_yr_t(i,m) = peatdeprow(i,m)      ! YW September 04,2015

      !> Accumulate the variables at the per tile level
      lterm_yr_t(i,m) = lterm_yr_t(i,m) + ltermrow(i,m) * oneOverDPY
      wetfdyn_yr_t(i,m) = wetfdyn_yr_t(i,m) + wetfdynrow(i,m) * oneOverDPY
      luc_emc_yr_t(i,m) = luc_emc_yr_t(i,m) + lucemcomrow(i,m) * oneOverDPY
      lucsocin_yr_t(i,m) = lucsocin_yr_t(i,m) + lucsocinrow(i,m) * oneOverDPY
      lucltrin_yr_t(i,m) = lucltrin_yr_t(i,m) + lucltrinrow(i,m) * oneOverDPY
      ch4WetSpec_yr_t(i,m) = ch4WetSpec_yr_t(i,m) + ch4WetSpecrow(i,m) * oneOverDPY
      ch4WetDyn_yr_t(i,m) = ch4WetDyn_yr_t(i,m) + ch4WetDynrow(i,m) * oneOverDPY
      ch4soills_yr_t(i,m) = ch4soills_yr_t(i,m) + ch4soillsrow(i,m) * oneOverDPY

      ! COMBAK PERLAY
      fProductDecomp_yr_t(i,m) = fProductDecomp_yr_t(i,m) &
                                 + (soilcresvegrow(i,m,iccp2) &
                                 + litresvegrow(i,m,iccp2)) * oneOverDPY
      ! NOTE: LUC product pools are only in layer 1.
      ! fProductDecomp_yr_t(i,m) = fProductDecomp_yr_t(i,m) &
      !                           + (soilcresvegrow(i,m,iccp2,1) &
      !                           + litresvegrow(i,m,iccp2,1)) * oneOverDPY
      ! COMBAK PERLAY
      ! NOTE: NBP is a special case here. The LUC product pool contributions are not
      ! per PFT, they exist uniformly across a tile, so they are not included in the
      ! nbp_yr calculation. Instead we need to use the nbp, not nbpveg variable
      ! for per tile and per gridcell outputting.
      nbp_yr_t(i,m) = nbp_yr_t(i,m) + nbprow(i,m) * oneOverDPY

    end do ! loop 883 ! m

    if (iday == lastDOY) then
      do m = 1,nmtest
        do j = 1,icc

          !> The pools are looked at just at the end of the year.
          stemmass_yr(i,m,j) = stemmassrow(i,m,j)
          rootmass_yr(i,m,j) = rootmassrow(i,m,j)
          veghght_yr(i,m,j) = veghghtrow(i,m,j)
          vgbiomas_yr(i,m,j) = vgbiomas_vegrow(i,m,j)
          totcmass_yr(i,m,j) = vgbiomas_yr(i,m,j)
          ! COMBAK PERLAY
          litrmass_yr(i,m,j) = litrmassrow(i,m,j)
          soilcmas_yr(i,m,j) = soilcmasrow(i,m,j)
          totcmass_yr(i,m,j) = totcmass_yr(i,m,j) + litrmass_yr(i,m,j) &
                               + soilcmas_yr(i,m,j)
          ! do k = 1,ignd
          !   litrmass_yr(i,m,j,k) = litrmassrow(i,m,j,k)
          !   soilcmas_yr(i,m,j,k) = soilcmasrow(i,m,j,k)
          !   totcmass_yr(i,m,j) = totcmass_yr(i,m,j) + litrmass_yr(i,m,j,k) &
          !                                           + soilcmas_yr(i,m,j,k)
          ! end do
          ! COMBAK PERLAY
        end do ! loop 925

        peatdep_yr_g(i) = peatdep_yr_g(i) + peatdep_yr_t(i,m) * farerot(i,m)    ! YW September 04,2015

        ! COMBAK PERLAY
        do k = 1,ignd
          litrmass_yr(i,m,iccp1) = litrmassrow(i,m,iccp1)
          soilcmas_yr(i,m,iccp1) = soilcmasrow(i,m,iccp1)
          totcmass_yr(i,m,iccp1) = totcmass_yr(i,m,iccp1) + litrmassrow(i,m,iccp1) + soilcmasrow(i,m,iccp1)
        end do
        ! do k = 1,ignd
        !   litrmass_yr(i,m,iccp1,k) = litrmassrow(i,m,iccp1,k)
        !   soilcmas_yr(i,m,iccp1,k) = soilcmasrow(i,m,iccp1,k)
        !   totcmass_yr(i,m,iccp1) = totcmass_yr(i,m,iccp1) + litrmassrow(i,m,iccp1,k) + soilcmasrow(i,m,iccp1,k)
        ! end do
        ! COMBAK PERLAY
        barefrac = 1.0

        !> Add values to the per tile vars
        ! NOTE: This implictly assumes that the fcancmx is only changing annually !
        do j = 1,icc

          laimaxg_yr_t(i,m) = laimaxg_yr_t(i,m) + laimaxg_yr(i,m,j) * fcancmxrow(i,m,j)
          stemmass_yr_t(i,m) = stemmass_yr_t(i,m) + stemmass_yr(i,m,j) * fcancmxrow(i,m,j)
          rootmass_yr_t(i,m) = rootmass_yr_t(i,m) + rootmass_yr(i,m,j) * fcancmxrow(i,m,j)
          vgbiomas_yr_t(i,m) = vgbiomas_yr_t(i,m) + vgbiomas_yr(i,m,j) * fcancmxrow(i,m,j)
          totcmass_yr_t(i,m) = totcmass_yr_t(i,m) + totcmass_yr(i,m,j) * fcancmxrow(i,m,j)
          npp_yr_t(i,m) = npp_yr_t(i,m) + npp_yr(i,m,j) * fcancmxrow(i,m,j)
          gpp_yr_t(i,m) = gpp_yr_t(i,m) + gpp_yr(i,m,j) * fcancmxrow(i,m,j)
          nep_yr_t(i,m) = nep_yr_t(i,m) + nep_yr(i,m,j) * fcancmxrow(i,m,j)
          ! nbp_yr_t(i,m)=nbp_yr_t(i,m)+nbp_yr(i,m,j)*fcancmxrow(i,m,j)
          emit_co2_yr_t(i,m) = emit_co2_yr_t(i,m) + emit_co2_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_co_yr_t(i,m) = emit_co_yr_t(i,m) + emit_co_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_ch4_yr_t(i,m) = emit_ch4_yr_t(i,m) + emit_ch4_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_nmhc_yr_t(i,m) = emit_nmhc_yr_t(i,m) + emit_nmhc_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_h2_yr_t(i,m) = emit_h2_yr_t(i,m) + emit_h2_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_nox_yr_t(i,m) = emit_nox_yr_t(i,m) + emit_nox_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_n2o_yr_t(i,m) = emit_n2o_yr_t(i,m) + emit_n2o_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_pm25_yr_t(i,m) = emit_pm25_yr_t(i,m) + emit_pm25_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_tpm_yr_t(i,m) = emit_tpm_yr_t(i,m) + emit_tpm_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_tc_yr_t(i,m) = emit_tc_yr_t(i,m) + emit_tc_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_oc_yr_t(i,m) = emit_oc_yr_t(i,m) + emit_oc_yr(i,m,j) * fcancmxrow(i,m,j)
          emit_bc_yr_t(i,m) = emit_bc_yr_t(i,m) + emit_bc_yr(i,m,j) * fcancmxrow(i,m,j)
          bterm_yr_t(i,m) = bterm_yr_t(i,m) + bterm_yr(i,m,j) * fcancmxrow(i,m,j)
          mterm_yr_t(i,m) = mterm_yr_t(i,m) + mterm_yr(i,m,j) * fcancmxrow(i,m,j)
          smfuncveg_yr_t(i,m) = smfuncveg_yr_t(i,m) + smfuncveg_yr(i,m,j) * fcancmxrow(i,m,j)
          hetrores_yr_t(i,m) = hetrores_yr_t(i,m) + hetrores_yr(i,m,j) * fcancmxrow(i,m,j)
          autores_yr_t(i,m) = autores_yr_t(i,m) + autores_yr(i,m,j) * fcancmxrow(i,m,j)
          burnfrac_yr_t(i,m) = burnfrac_yr_t(i,m) + burnfrac_yr(i,m,j) * fcancmxrow(i,m,j)
          veghght_yr_t(i,m) = veghght_yr_t(i,m) + veghght_yr(i,m,j) * fcancmxrow(i,m,j)

          barefrac = barefrac - fcancmxrow(i,m,j)
          ! COMBAK PERLAY
          litrmass_yr_t(i,m) = litrmass_yr_t(i,m) + litrmass_yr(i,m,j) * fcancmxrow(i,m,j)
          soilcmas_yr_t(i,m) = soilcmas_yr_t(i,m) + soilcmas_yr(i,m,j) * fcancmxrow(i,m,j)
          litres_yr_t(i,m) = litres_yr_t(i,m) + litres_yr(i,m,j) * fcancmxrow(i,m,j)
          soilcres_yr_t(i,m) = soilcres_yr_t(i,m) + soilcres_yr(i,m,j) * fcancmxrow(i,m,j)
          ! do k = 1,ignd
          !   litrmass_yr_t(i,m,k) = litrmass_yr_t(i,m,k) + litrmass_yr(i,m,j,k) * fcancmxrow(i,m,j)
          !   soilcmas_yr_t(i,m,k) = soilcmas_yr_t(i,m,k) + soilcmas_yr(i,m,j,k) * fcancmxrow(i,m,j)
          !   litres_yr_t(i,m,k) = litres_yr_t(i,m,k) + litres_yr(i,m,j,k) * fcancmxrow(i,m,j)
          !   soilcres_yr_t(i,m,k) = soilcres_yr_t(i,m,k) + soilcres_yr(i,m,j,k) * fcancmxrow(i,m,j)
          ! end do
          ! COMBAK PERLAY

        end do ! j
        ! COMBAK PERLAY
        litrmass_yr_t(i,m) = litrmass_yr_t(i,m) + litrmass_yr(i,m,iccp1) * barefrac
        soilcmas_yr_t(i,m) = soilcmas_yr_t(i,m) + soilcmas_yr(i,m,iccp1) * barefrac
        litres_yr_t(i,m) = litres_yr_t(i,m) + litres_yr(i,m,iccp1) * barefrac
        soilcres_yr_t(i,m) = soilcres_yr_t(i,m) + soilcres_yr(i,m,iccp1) * barefrac
        totcmass_yr_t(i,m) = totcmass_yr_t(i,m) &
                             + (litrmass_yr(i,m,iccp1) + soilcmas_yr(i,m,iccp1)) * barefrac
        ! do k = 1,ignd
        !   litrmass_yr_t(i,m,k) = litrmass_yr_t(i,m,k) + litrmass_yr(i,m,iccp1,k) * barefrac
        !   soilcmas_yr_t(i,m,k) = soilcmas_yr_t(i,m,k) + soilcmas_yr(i,m,iccp1,k) * barefrac
        !   litres_yr_t(i,m,k) = litres_yr_t(i,m,k) + litres_yr(i,m,iccp1,k) * barefrac
        !   soilcres_yr_t(i,m,k) = soilcres_yr_t(i,m,k) + soilcres_yr(i,m,iccp1,k) * barefrac
        !   totcmass_yr_t(i,m) = totcmass_yr_t(i,m) &
        !                     + (litrmass_yr(i,m,iccp1,k) + soilcmas_yr(i,m,iccp1,k)) * barefrac
        ! end do
        ! COMBAK PERLAY
        hetrores_yr_t(i,m) = hetrores_yr_t(i,m) + hetrores_yr(i,m,iccp1) * barefrac
        nep_yr_t(i,m) = nep_yr_t(i,m) + nep_yr(i,m,iccp1) * barefrac
        ! nbp_yr_t(i,m)=nbp_yr_t(i,m)+nbp_yr(i,m,iccp1)*barefrac

        !> Add values to the per gridcell vars
        laimaxg_yr_g(i) = laimaxg_yr_g(i) + laimaxg_yr_t(i,m) * FAREROT(i,m)
        stemmass_yr_g(i) = stemmass_yr_g(i) + stemmass_yr_t(i,m) * FAREROT(i,m)
        rootmass_yr_g(i) = rootmass_yr_g(i) + rootmass_yr_t(i,m) * FAREROT(i,m)
        ! COMBAK PERLAY
        litrmass_yr_g(i) = litrmass_yr_g(i) + litrmass_yr_t(i,m) * FAREROT(i,m)
        soilcmas_yr_g(i) = soilcmas_yr_g(i) + soilcmas_yr_t(i,m) * FAREROT(i,m)
        litres_yr_g(i) = litres_yr_g(i) + litres_yr_t(i,m) * FAREROT(i,m)
        soilcres_yr_g(i) = soilcres_yr_g(i) + soilcres_yr_t(i,m) * FAREROT(i,m)
        ! do k = 1,ignd
        !   litrmass_yr_g(i,k) = litrmass_yr_g(i,k) + litrmass_yr_t(i,m,k) * FAREROT(i,m)
        !   soilcmas_yr_g(i,k) = soilcmas_yr_g(i,k) + soilcmas_yr_t(i,m,k) * FAREROT(i,m)
        !   litres_yr_g(i,k) = litres_yr_g(i,k) + litres_yr_t(i,m,k) * FAREROT(i,m)
        !   soilcres_yr_g(i,k) = soilcres_yr_g(i,k) + soilcres_yr_t(i,m,k) * FAREROT(i,m)
        ! end do
        ! COMBAK PERLAY
        vgbiomas_yr_g(i) = vgbiomas_yr_g(i) + vgbiomas_yr_t(i,m) * FAREROT(i,m)
        totcmass_yr_g(i) = totcmass_yr_g(i) + totcmass_yr_t(i,m) * FAREROT(i,m)
        npp_yr_g(i) = npp_yr_g(i) + npp_yr_t(i,m) * FAREROT(i,m)
        gpp_yr_g(i) = gpp_yr_g(i) + gpp_yr_t(i,m) * FAREROT(i,m)
        nep_yr_g(i) = nep_yr_g(i) + nep_yr_t(i,m) * FAREROT(i,m)
        nbp_yr_g(i) = nbp_yr_g(i) + nbp_yr_t(i,m) * FAREROT(i,m)
        emit_co2_yr_g(i) = emit_co2_yr_g(i) + emit_co2_yr_t(i,m) * FAREROT(i,m)
        emit_co_yr_g(i) = emit_co_yr_g(i) + emit_co_yr_t(i,m) * FAREROT(i,m)
        emit_ch4_yr_g(i) = emit_ch4_yr_g(i) + emit_ch4_yr_t(i,m) * FAREROT(i,m)
        emit_nmhc_yr_g(i) = emit_nmhc_yr_g(i) + emit_nmhc_yr_t(i,m) * FAREROT(i,m)
        emit_h2_yr_g(i) = emit_h2_yr_g(i) + emit_h2_yr_t(i,m) * FAREROT(i,m)
        emit_nox_yr_g(i) = emit_nox_yr_g(i) + emit_nox_yr_t(i,m) * FAREROT(i,m)
        emit_n2o_yr_g(i) = emit_n2o_yr_g(i) + emit_n2o_yr_t(i,m) * FAREROT(i,m)
        emit_pm25_yr_g(i) = emit_pm25_yr_g(i) + emit_pm25_yr_t(i,m) * FAREROT(i,m)
        emit_tpm_yr_g(i) = emit_tpm_yr_g(i) + emit_tpm_yr_t(i,m) * FAREROT(i,m)
        emit_tc_yr_g(i) = emit_tc_yr_g(i) + emit_tc_yr_t(i,m) * FAREROT(i,m)
        emit_oc_yr_g(i) = emit_oc_yr_g(i) + emit_oc_yr_t(i,m) * FAREROT(i,m)
        emit_bc_yr_g(i) = emit_bc_yr_g(i) + emit_bc_yr_t(i,m) * FAREROT(i,m)
        hetrores_yr_g(i) = hetrores_yr_g(i) + hetrores_yr_t(i,m) * FAREROT(i,m)
        autores_yr_g(i) = autores_yr_g(i) + autores_yr_t(i,m) * FAREROT(i,m)
        burnfrac_yr_g(i) = burnfrac_yr_g(i) + burnfrac_yr_t(i,m) * FAREROT(i,m)
        smfuncveg_yr_g(i) = smfuncveg_yr_g(i) + smfuncveg_yr_t(i,m) * FAREROT(i,m)
        bterm_yr_g(i) = bterm_yr_g(i) + bterm_yr_t(i,m) * FAREROT(i,m)
        lterm_yr_g(i) = lterm_yr_g(i) + lterm_yr_t(i,m) * FAREROT(i,m)
        mterm_yr_g(i) = mterm_yr_g(i) + mterm_yr_t(i,m) * FAREROT(i,m)
        luc_emc_yr_g(i) = luc_emc_yr_g(i) + luc_emc_yr_t(i,m) * FAREROT(i,m)
        lucsocin_yr_g(i) = lucsocin_yr_g(i) + lucsocin_yr_t(i,m) * FAREROT(i,m)
        lucltrin_yr_g(i) = lucltrin_yr_g(i) + lucltrin_yr_t(i,m) * FAREROT(i,m)
        ch4WetSpec_yr_g(i) = ch4WetSpec_yr_g(i) + ch4WetSpec_yr_t(i,m) * FAREROT(i,m)
        wetfdyn_yr_g(i) = wetfdyn_yr_g(i) + wetfdyn_yr_t(i,m) * FAREROT(i,m)
        ch4WetDyn_yr_g(i) = ch4WetDyn_yr_g(i) + ch4WetDyn_yr_t(i,m) * FAREROT(i,m)
        ch4soills_yr_g(i) = ch4soills_yr_g(i) + ch4soills_yr_t(i,m) * FAREROT(i,m)
        veghght_yr_g(i) = veghght_yr_g(i) + veghght_yr_t(i,m) * FAREROT(i,m)

        ! Including the LUC product pools. They are per tile values and
        ! are assumed to occupy the whole tile.
        fProductDecomp_yr_g(i) = fProductDecomp_yr_g(i) &
                                 + fProductDecomp_yr_t(i,m) * FAREROT(i,m)
        ! COMBAK PERLAY
        cProduct_yr_g(i) = cProduct_yr_g(i) &
                           + (litrmassrow(i,m,iccp2) &
                           + soilcmasrow(i,m,iccp2)) * FAREROT(i,m)
        ! cProduct_yr_g(i) = cProduct_yr_g(i) &
        !                    + (litrmassrow(i,m,iccp2,1) &
        !                    + soilcmasrow(i,m,iccp2,1)) * FAREROT(i,m)
        ! COMBAK PERLAY

      end do ! loop 900 ! m

      !> Write to annual output files:

      ! Prepare the timestamp for this year
      timeStamp = consecDays

      !> First write out the per gridcell values
      call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr_g' ,timeStamp,'lai', [laimaxg_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr_g',timeStamp,'cVeg',[vgbiomas_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr_g',timeStamp,'cStem',[stemmass_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr_g',timeStamp,'cRoot',[rootmass_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr_g',timeStamp,'cLand',[totcmass_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr_g'     ,timeStamp,'npp',[npp_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr_g'     ,timeStamp,'gpp',[gpp_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr_g'     ,timeStamp,'nep',[nep_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr_g'     ,timeStamp,'nbp',[nbp_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr_g',timeStamp,'rh',[hetrores_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr_g' ,timeStamp,'ra',[autores_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr_g' ,timeStamp,'vegHeight',[veghght_yr_g(i)])
      
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_yr_g' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_yr_g' ,timeStamp,'wetlandFrac',[wetfdyn_yr_g(i)])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_yr_g' ,timeStamp,'soilCH4cons',[ch4soills_yr_g(i)])

      ! COMBAK PERLAY
      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litrMassPerLay_yr_g',timeStamp,'cLitterperlay',[litrmass_yr_g(i,:)])
      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilCMasPerLay_yr_g',timeStamp,'cSoilperlay',[soilcmas_yr_g(i,:)])
      !
      ! ! Make the bulk litter and soil C pool and respiration temporary variables:
      ! allocate(bulkLitterCarbon(1))
      ! allocate(bulkSoilCarbon(1))
      ! allocate(bulkLitterResp(1))
      ! allocate(bulkSoilResp(1))
      ! bulkLitterCarbon(1) = sum(litrmass_yr_g(i,:))
      ! bulkSoilCarbon(1) = sum(soilcmas_yr_g(i,:))
      ! bulkLitterResp(1) = sum(litres_yr_g(i,:))
      ! bulkSoilResp(1) = sum(soilcres_yr_g(i,:))
      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_g',timeStamp,'cLitter',[bulkLitterCarbon])
      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_g',timeStamp,'cSoil',[bulkSoilCarbon])
      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_g'  ,timeStamp,'rhLitter',[bulkLitterResp])
      ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_g',timeStamp,'rhSoil',[bulkSoilResp])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_g',timeStamp,'cLitter',[litrmass_yr_g])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_g',timeStamp,'cSoil',[soilcmas_yr_g])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_g'  ,timeStamp,'rhLitter',[litres_yr_g])
      call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_g',timeStamp,'rhSoil',[soilcres_yr_g])
      ! deallocate(bulkLitterCarbon)
      ! deallocate(bulkSoilCarbon)
      ! deallocate(bulkLitterResp)
      ! deallocate(bulkSoilResp)
      ! COMBAK PERLAY
      if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
        call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_yr_g' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_yr_g(i)])
      end if

      if (PFTCompetition) then 
        do m = 1,nmtest
          sumfare = 0.0
          fcancmxNoSeed = 0.0
          pftExist = 0.0
          do j = 1,icc
            if (pftexistrow(i,1,j)) then
              pftExist(j) = 1.0
            end if
              sumfare = sumfare + fcancmxrow(i,m,j)
              fcancmxNoSeed(j) = fcancmxrow(i,m,j)
          end do
          call writeOutput1D(lonLocalIndex,latLocalIndex,'pftexistrow_yr_g' ,timeStamp,'landCoverExist',[pftExist])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_yr_g' ,timeStamp,'landCoverFrac',[fcancmxNoSeed(1:icc),1 - sumfare])
        end do
      else
        do m = 1,nmtest
          sumfare = 0.0
          do j = 1,icc
            sumfare = sumfare + fcancmxrow(i,m,j)
          end do ! j
          call writeOutput1D(lonLocalIndex,latLocalIndex,'fcancmxrow_yr_g' ,timeStamp,'landCoverFrac',[fcancmxrow(i,m,1:icc),1 - sumfare])
        end do ! m
      end if

      if (dofire .or. lnduseon) then
        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_yr_g' ,timeStamp,'fFireCH4',[emit_ch4_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr_g' ,timeStamp,'fFire',[emit_co2_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr_g' ,timeStamp,'burntFractionAll',[burnfrac_yr_g(i)])
      end if
      if (lnduseon) then
        call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_yr_g' ,timeStamp,'fDeforestToAtmos',[luc_emc_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_yr_g' ,timeStamp,'fDeforestToLitter',[lucltrin_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_yr_g' ,timeStamp,'fDeforestToSoil',[lucsocin_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_yr_g' ,timeStamp,'fDeforestTotal', &
                           [lucsocin_yr_g(i) + lucltrin_yr_g(i) + luc_emc_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'fProductDecomp_yr_g' ,timeStamp,'fProductDecomp',[fProductDecomp_yr_g(i)])
        call writeOutput1D(lonLocalIndex,latLocalIndex,'cProduct_yr_g' ,timeStamp,'cProduct',[cProduct_yr_g(i)])
      end if

      if (doperpftoutput) then
        if (nmtest > 1) then
          print * ,'Per PFT and per tile outputs not implemented yet'
        else

          m = 1 ! FLAG only implemented for composite mode,tiles == 1 !! !

          call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr' ,timeStamp,'lai', [laimaxg_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr',timeStamp,'cVeg',[vgbiomas_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr',timeStamp,'cStem',[stemmass_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr',timeStamp,'cRoot',[rootmass_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr',timeStamp,'cLand',[totcmass_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr'     ,timeStamp,'npp',[npp_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr'     ,timeStamp,'gpp',[gpp_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr'     ,timeStamp,'nep',[nep_yr(i,m,:)])
          ! NOTE: This NBP does not include LUC product pool contributions since they are
          ! not per PFT but rather per tile
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr'     ,timeStamp,'nbp',[nbp_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr',timeStamp,'rh',[hetrores_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr' ,timeStamp,'ra',[autores_yr(i,m,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr' ,timeStamp,'vegHeight',[veghght_yr(i,m,:)])
          
          ! For FireMIP calc the low vs high cover.
          coverLow = 0.
          coverHigh = 0.
          do j = 1, icc
            if (fcancmxrow(i,m,j) > 0.) then 
              if (veghght_yr(i,m,j) < 5.) then
                coverLow = coverLow + fcancmxrow(i,m,j)
              else 
                coverHigh = coverHigh + fcancmxrow(i,m,j)
              end if 
            end if 
          end do
          call writeOutput1D(lonLocalIndex,latLocalIndex,'coverLow' ,timeStamp,'coverLow',[coverLow])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'coverHigh' ,timeStamp,'coverHigh',[coverHigh])


          ! COMBAK PERLAY
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr',timeStamp,'cLitter',[litrmass_yr(i,m,1:iccp1)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr',timeStamp,'cSoil',[soilcmas_yr(i,m,1:iccp1)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr'  ,timeStamp,'rhLitter',[litres_yr(i,m,1:iccp1)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr',timeStamp,'rhSoil',[soilcres_yr(i,m,1:iccp1)])
          ! ! Make the bulk litter and soil C pool and respiration temporary variables:
          ! allocate(bulkLitterCarbon(iccp1))
          ! allocate(bulkSoilCarbon(iccp1))
          ! allocate(bulkLitterResp(iccp1))
          ! allocate(bulkSoilResp(iccp1))
          !
          ! do k = 1,iccp1
          !   bulkLitterCarbon(k) = sum(litrmass_yr(i,m,k,:))
          !   bulkSoilCarbon(k) = sum(soilcmas_yr(i,m,k,:))
          !   bulkLitterResp(k) = sum(litres_yr(i,m,k,:))
          !   bulkSoilResp(k) = sum(soilcres_yr(i,m,k,:))
          ! end do
          !
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr',timeStamp,'cLitter',[bulkLitterCarbon])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr',timeStamp,'cSoil',[bulkSoilCarbon])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr'  ,timeStamp,'rhLitter',[bulkLitterResp])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr',timeStamp,'rhSoil',[bulkSoilResp])
          !
          ! deallocate(bulkLitterCarbon)
          ! deallocate(bulkSoilCarbon)
          ! deallocate(bulkLitterResp)
          ! deallocate(bulkSoilResp)
          ! COMBAK PERLAY

          if (dofire .or. lnduseon) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_yr' ,timeStamp,'fFireCH4',[emit_ch4_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr' ,timeStamp,'fFire',[emit_co2_yr(i,m,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr' ,timeStamp,'burntFractionAll',[burnfrac_yr(i,m,:)])
          end if

        end if
      end if

      if (dopertileoutput) then

        if (nmtest == 1) then
          print * ,'Switch selected for per tile output but number of tiles is only one.'
        else

          !> Write out the per tile values
          call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_yr_t' ,timeStamp,'lai', [laimaxg_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'vgbiomas_yr_t',timeStamp,'cVeg',[vgbiomas_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'stemmass_yr_t',timeStamp,'cStem',[stemmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'rootmass_yr_t',timeStamp,'cRoot',[rootmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'totcmass_yr_t',timeStamp,'cLand',[totcmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'npp_yr_t'     ,timeStamp,'npp',[npp_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'gpp_yr_t'     ,timeStamp,'gpp',[gpp_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nep_yr_t'     ,timeStamp,'nep',[nep_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'nbp_yr_t'     ,timeStamp,'nbp',[nbp_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'hetrores_yr_t',timeStamp,'rh',[hetrores_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'autores_yr_t' ,timeStamp,'ra',[autores_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'veghght_yr_t' ,timeStamp,'vegHeight',[veghght_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetDyn_yr_t' ,timeStamp,'wetlandCH4dyn',[ch4WetDyn_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'wetfdyn_yr_t' ,timeStamp,'wetlandFrac',[wetfdyn_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4soills_yr_t' ,timeStamp,'soilCH4cons',[ch4soills_yr_t(i,:)])

          ! COMBAK PERLAY
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_t',timeStamp,'cLitter',[litrmass_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_t',timeStamp,'cSoil',[soilcmas_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_t'  ,timeStamp,'rhLitter',[litres_yr_t(i,:)])
          call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_t',timeStamp,'rhSoil',[soilcres_yr_t(i,:)])
          ! ! Make the bulk litter and soil C pool and respiration temporary variables:
          ! allocate(bulkLitterCarbon(nmtest))
          ! allocate(bulkSoilCarbon(nmtest))
          ! allocate(bulkLitterResp(nmtest))
          ! allocate(bulkSoilResp(nmtest))
          !
          ! do m = 1,nmtest
          !   bulkLitterCarbon(m) = sum(litrmass_yr_t(i,m,:))
          !   bulkSoilCarbon(m) = sum(soilcmas_yr_t(i,m,:))
          !   bulkLitterResp(m) = sum(litres_yr_t(i,m,:))
          !   bulkSoilResp(m) = sum(soilcres_yr_t(i,m,:))
          ! end do
          !
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litrmass_yr_t',timeStamp,'cLitter',[bulkLitterCarbon])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcmas_yr_t',timeStamp,'cSoil',[bulkSoilCarbon])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'litres_yr_t'  ,timeStamp,'rhLitter',[bulkLitterResp])
          ! call writeOutput1D(lonLocalIndex,latLocalIndex,'soilcres_yr_t',timeStamp,'rhSoil',[bulkSoilResp])
          !
          ! deallocate(bulkLitterCarbon)
          ! deallocate(bulkSoilCarbon)
          ! deallocate(bulkLitterResp)
          ! deallocate(bulkSoilResp)
          ! COMBAK PERLAY

          if (transientOBSWETF .or. fixedYearOBSWETF /= - 9999) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'ch4WetSpec_yr_t' ,timeStamp,'wetlandCH4spec',[ch4WetSpec_yr_t(i,:)])
          end if
          if (dofire .or. lnduseon) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_ch4_yr_t' ,timeStamp,'fFireCH4',[emit_ch4_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'emit_co2_yr_t' ,timeStamp,'fFire',[emit_co2_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'burnfrac_yr_t' ,timeStamp,'burntFractionAll',[burnfrac_yr_t(i,:)])
          end if
          if (lnduseon) then
            call writeOutput1D(lonLocalIndex,latLocalIndex,'luc_emc_yr_t' ,timeStamp,'fDeforestToAtmos',[luc_emc_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'lucltrin_yr_t' ,timeStamp,'fDeforestToLitter',[lucltrin_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'lucsocin_yr_t' ,timeStamp,'fDeforestToSoil',[lucsocin_yr_t(i,:)])
            call writeOutput1D(lonLocalIndex,latLocalIndex,'luctot_yr_t' ,timeStamp,'fDeforestTotal', &
                               [lucsocin_yr_t(i,:) + lucltrin_yr_t(i,:) + luc_emc_yr_t(i,:)])
          end if
        end if
      end if

      !> Reset all annual vars in preparation for next year
      call resetYearEnd(nltest, nmtest)

    end if ! if iday=365/366

  end subroutine ctem_annual_aw
  !! @}
  !> \namespace prepareoutputs
  !> Central module that handles all CTEM preparation and writing of output files

end module prepareOutputs
