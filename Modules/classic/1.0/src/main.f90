!> \file
!> Main model driver for CLASSIC in stand-alone mode using specified boundary
!! conditions and atmospheric forcing.
!! @author D. Verseghy, V. Arora, J. Melton, E. Chan
!! This driver program initializes the run, reads in CLASSIC input files,
!! manages the run and the coupling between CLASS and CTEM, calls subroutines
!! that aggregate and write outputs, and closes the run for this grid cell.

module main

  implicit none

  public :: main_driver

contains

  !> \ingroup main_main_driver
  !> @{
  !>
  !! ## Dimension statements.
  !!
  !!     ### First set of definitions:
  !!     Background variables, and prognostic and diagnostic
  !!     variables normally provided by and/or used by the GCM.
  !!     The suffix "rot" refers to variables existing on the
  !!     mosaic grid on the current latitude circle.  The suffix
  !!     "gat" refers to the same variables after they have undergone
  !!     a "gather" operation in which the two mosaic dimensions
  !!     are collapsed into one.  The suffix "row" refers both to
  !!     grid-constant input variables. and to grid-averaged
  !!     diagnostic variables.
  !!
  !!     The first dimension element of the "rot" variables
  !!     refers to the number of grid cells on the current
  !!     latitude circle.  In this stand-alone version, this
  !!     number is set to 1, the second dimension
  !!     element of the "rot" variables refers to the maximum
  !!     number of tiles in the mosaic.  The first
  !!     dimension element in the "gat" variables is given by
  !!     the product of the first two dimension elements in the
  !!     "rot" variables.
  !!
  !!     The majority of CTEM parameters are stored in classicParams.f90.
  !!     Also the CLASS and CTEM variables are stored in modules that we point to
  !!     in this driver. We access the variables and parameters
  !!     through use statements for modules:

  subroutine main_driver (longitude, latitude, lonIndex, latIndex, lonLocalIndex, latLocalIndex)

    use classicParams,      only : nlat, nmos, ilg, nmon, ican, ignd, icc, monthend, &
                                   modelpft, l2max, deltat, NBS, readInParams, nol2pfts, &
                                   DELT, TFREZ, zbldJobOpt, zrfhJobOpt, zrfmJobOpt, &
                                   REFF0_LAND, ZSNMIN, ZSNMAX2
    use landuseChange,      only : initializeLandCover
    use ctemStateVars,      only : vrot, vgat, c_switch, initRowVars, &
                                   resetMonthEnd, resetYearEnd, &
                                   ctem_tile, resetMosaicAccum, tracer
    use classStateVars,     only : class_gat, class_rot, resetAccVars, &
                                   resetClassMon, resetClassYr, initDiagnosticVars
    use prepareOutputs,     only : class_monthly_aw, ctem_annual_aw, ctem_monthly_aw, &
                                   ctem_daily_aw, class_annual_aw, class_hh_w, class_daily_aw, &
                                   convertUnitsCTEM
    use modelStateDrivers,  only : read_initialstate, write_restart, getInput, updateInput, deallocInput, getMet, updateMet
    use generalUtils,       only : findDaylength, findLeapYears, run_model, findCloudiness, &
                                   findPermafrostVars, initRandomSeed, checksumCalc
    use classGatherScatter, only : classGather, classGatherPrep, classScatter
    use ctemGatherScatter,  only : ctems2, ctemg1, ctemg2
    use ctemUtilities,      only : dayEndCTEMPreparation, accumulateForCTEM, ctemInit
    use metDisaggModule,    only : disaggMet
    use outputManager,      only : consecDays
    use ctemDriver,         only : ctem
    use tracerModule,       only : decay14C
    use applyAllometry,     only : allometry
    use fourBandAlbedo,     only : fourBandDriver

    implicit none

    ! Arguments
    real, intent(in) :: longitude, latitude                 !< Longitude/latitude of grid cell (degrees)
    integer, intent(in) :: lonIndex, latIndex               !< Index of grid cell being run on the input files grid
    integer, intent(in) :: lonLocalIndex, latLocalIndex     !< Index of grid cell being run on the output files grid

    ! Local variables
    integer :: lastDOY             !< Initialized to 365 days, can be overwritten later is leap = true and it is a leap year.
    integer :: metTimeIndex        !< Counter used to move through the meteorological input arrays
    logical :: metDone             !< Logical switch when the end of the stored meteorological array is reached.
    integer :: runyr               !< Year of the model run (counts up starting with readMetStartYear continously, even if metLoop > 1)

    integer :: NLTEST  !< Number of grid cells being modelled for this run
    integer :: NMTEST  !< Number of mosaic tiles per grid cell being modelled for this run
    integer :: NCOUNT  !< Counter for daily averaging
    integer :: NDAY    !< Number of short (physics) timesteps in one day. e.g., if physics timestep is 15 min this is 48.
    integer :: IMONTH !< Month of the year simulation is in.
    integer :: DOM  !< Day of month counter
    integer :: NT      !<
    integer :: IHOUR   !< Hour of day
    integer :: IMIN    !< Minutes elapsed in current hour
    integer :: IDAY    !< Julian day of the year
    integer :: IYEAR   !< Year of run
    integer :: NML     !< Counter representing number of mosaic tiles on modelled domain that are land
    integer :: NMW     !< Counter representing number of mosaic tiles on modelled domain that are lakes
    integer :: JLAT    !< Integer index corresponding to latitude of grid cell
    integer :: NLANDCS !< Number of modelled areas that contain subareas of canopy over snow
    integer :: NLANDGS !< Number of modelled areas that contain subareas of snow over bare ground
    integer :: NLANDC  !< Number of modelled areas that contain subareas of canopy over bare ground
    integer :: NLANDG  !< Number of modelled areas that contain subareas of bare ground
    integer :: NLANDI  !< Number of modelled areas that are ice sheets
    integer :: I, J, K, L, M, N
    integer :: NTLD    !<

    ! Pointers
    integer, pointer :: readMetStartYear    !< First year of meteorological forcing to read in from the met file

    ! The following are stored in the data structure: class_gat
    ! they are allocatted in allocClassVars in the classStateVars
    ! module and pointed to here.

    ! These will be allocated the dimension: 'ignd'

    real, pointer, dimension(:) :: DELZ    !<
    real, pointer, dimension(:) :: ZBOT    !<

    ! These will be allocated the dimension: 'ilg'

    integer, pointer, dimension(:) :: ILMOS !< Index of gridcell corresponding to current element of gathered vector of land surface variables [ ]
    integer, pointer, dimension(:) :: JLMOS !< Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ]
    integer, pointer, dimension(:) :: IWMOS !< Index of gridcell corresponding to current element of gathered vector of inland water body variables [ ]
    integer, pointer, dimension(:) :: JWMOS !< Index of mosaic tile corresponding to current element of gathered vector of inland water body variables [ ]
    integer, pointer, dimension(:) :: IGDRGAT   !< Index of soil layer in which bedrock is encountered

    real, pointer, dimension(:) :: GCGAT   !< Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
    real, pointer, dimension(:) :: TZSGAT  !< Vertical temperature gradient in a snow pack
    real, pointer, dimension(:) :: PCSNGAT !< Snow fall flux \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ALBSGAT !< Snow albedo [ ]
    real, pointer, dimension(:) :: CMAIGAT !< Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: GROGAT  !< Vegetation growth index [ ]
    real, pointer, dimension(:) :: QACGAT  !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: RCANGAT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: RHOSGAT !< Density of snow \f$[kg m^{-3} ]\f$
    real, pointer, dimension(:) :: SCANGAT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: SNOGAT  !< Mass of snow pack [kg m^{-2} ]\f$
    real, pointer, dimension(:) :: TACGAT  !< Temperature of air within vegetation canopy [K]
    real, pointer, dimension(:) :: TBASGAT !< Temperature of bedrock in third soil layer [K]
    real, pointer, dimension(:) :: TCANGAT !< Vegetation canopy temperature [K]
    real, pointer, dimension(:) :: TPNDGAT !< Temperature of ponded water [K]
    real, pointer, dimension(:) :: TSNOGAT !< Snowpack temperature [K]
    real, pointer, dimension(:) :: WSNOGAT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:) :: ZPNDGAT !< Depth of ponded water on surface [m]
    real, pointer, dimension(:) :: REFGAT  !< Snow grain size (for ISNOALB=1 option)  [m]
    real, pointer, dimension(:) :: BCSNGAT !<
    real, pointer, dimension(:) :: AGIDGAT !< Optional user-specified value of ground near-infrared albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: AGVDGAT !< Optional user-specified value of ground visible albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: ALGDGAT !< Reference albedo for dry soil [ ]
    real, pointer, dimension(:) :: ALGWGAT !< Reference albedo for saturated soil [ ]
    real, pointer, dimension(:) :: ASIDGAT !< Optional user-specified value of snow near-infrared albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: ASVDGAT !< Optional user-specified value of snow visible albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:) :: DRNGAT  !< Drainage index at bottom of soil profile [ ]
    real, pointer, dimension(:) :: GRKFGAT !< WATROF parameter used when running MESH code [ ]
    real, pointer, dimension(:) :: WFCIGAT !< WATROF parameter used when running MESH code [ ]
    real, pointer, dimension(:) :: WFSFGAT !< WATROF parameter used when running MESH code [ ]
    real, pointer, dimension(:) :: XSLPGAT !< Surface slope (used when running MESH code) [degrees]
    real, pointer, dimension(:) :: ZPLGGAT !< Maximum water ponding depth for snow-free subareas (user-specified when running MESH code) [m]
    real, pointer, dimension(:) :: ZPLSGAT !< Maximum water ponding depth for snow-covered subareas (user-specified when running MESH code) [m]
    real, pointer, dimension(:) :: ZSNLGAT !< Limiting snow depth below which coverage is < 100% [m]
    real, pointer, dimension(:) :: ALGWVGAT !<
    real, pointer, dimension(:) :: ALGWNGAT !<
    real, pointer, dimension(:) :: ALGDVGAT !<
    real, pointer, dimension(:) :: ALGDNGAT !<
    real, pointer, dimension(:) :: EMISGAT  !<
    real, pointer, dimension(:) :: CSZGAT  !< Cosine of solar zenith angle [ ]
    real, pointer, dimension(:) :: DLONGAT !< Longitude of grid cell (east of Greenwich) [degrees]
    real, pointer, dimension(:) :: DLATGAT !< Latitude of grid cell [degrees]
    real, pointer, dimension(:) :: FCLOGAT !< Fractional cloud cover [ ]
    real, pointer, dimension(:) :: FDLGAT  !< Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSIHGAT !< Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSVHGAT !< Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: GGEOGAT !< Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: PADRGAT !< Partial pressure of dry air [Pa]
    real, pointer, dimension(:) :: PREGAT  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PRESGAT !< Surface air pressure [Pa]
    real, pointer, dimension(:) :: QAGAT   !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: RADJGAT !< Latitude of grid cell (positive north of equator) [rad]
    real, pointer, dimension(:) :: RHOAGAT !< Density of air \f$[kg m^{-3} ]\f$
    real, pointer, dimension(:) :: RHSIGAT !< Density of fresh snow \f$[kg m^{-3} ]\f$
    real, pointer, dimension(:) :: RPCPGAT !< Rainfall rate over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: SPCPGAT !< Snowfall rate over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: TAGAT   !< Air temperature at reference height [K]
    real, pointer, dimension(:) :: TADPGAT !< Dew point temperature of air [K]
    real, pointer, dimension(:) :: TRPCGAT !< Rainfall temperature [K]
    real, pointer, dimension(:) :: TSPCGAT !< Snowfall temperature [K]
    real, pointer, dimension(:) :: ULGAT   !< Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: VLGAT   !< Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: VMODGAT !< Wind speed at reference height \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: VPDGAT  !< Vapour pressure deficit [mb]
    real, pointer, dimension(:) :: Z0ORGAT !< Orographic roughness length [m]
    real, pointer, dimension(:) :: ZBLDGAT !< Atmospheric blending height for surface roughness length averaging [m]
    real, pointer, dimension(:) :: ZDHGAT  !< User-specified height associated with diagnosed screen-level variables [m]
    real, pointer, dimension(:) :: ZDMGAT  !< User-specified height associated with diagnosed anemometer-level wind speed [m]
    real, pointer, dimension(:) :: ZRFHGAT !< Reference height associated with forcing air temperature and humidity [m]
    real, pointer, dimension(:) :: ZRFMGAT !< Reference height associated with forcing wind speed [m]
    real, pointer, dimension(:) :: FSGGAT  !<
    real, pointer, dimension(:) :: FLGGAT  !<
    real, pointer, dimension(:) :: GUSTGAT !<
    real, pointer, dimension(:) :: DEPBGAT !<
    real, pointer, dimension(:) :: GTBS    !<
    real, pointer, dimension(:) :: SFCUBS  !<
    real, pointer, dimension(:) :: SFCVBS  !<
    real, pointer, dimension(:) :: USTARBS !<
    real, pointer, dimension(:) :: TCSNOW  !< Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
    real, pointer, dimension(:) :: GSNOW   !< Diagnostic heat flux at snow surface for use in CCCma black carbon deposition scheme \f$[W m^{-2}]\f$
    real, pointer, dimension(:) :: ALIRGAT !< Diagnosed total near-infrared albedo of land surface [ ]
    real, pointer, dimension(:) :: ALVSGAT !< Diagnosed total visible albedo of land surface [ ]
    real, pointer, dimension(:) :: CDHGAT  !< Surface drag coefficient for heat [ ]
    real, pointer, dimension(:) :: CDMGAT  !< Surface drag coefficient for momentum [ ]
    real, pointer, dimension(:) :: DRGAT   !< Surface drag coefficient under neutral stability [ ]
    real, pointer, dimension(:) :: EFGAT   !< Evaporation efficiency at ground surface [ ]
    real, pointer, dimension(:) :: FLGGGAT !< Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FLGSGAT !< Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FLGVGAT !< Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGGGAT !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGSGAT !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSGVGAT !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: FSNOGAT !< Diagnosed fractional snow coverage [ ]
    real, pointer, dimension(:) :: GAGAT   !< Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: GTGAT   !< Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:) :: HBLGAT  !< Height of the atmospheric boundary layer [m]
    real, pointer, dimension(:) :: HEVCGAT !< Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HEVGGAT !< Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HEVSGAT !< Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSGAT  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSCGAT !< Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSGGAT !< Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HFSSGAT !< Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HMFCGAT !< Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HMFNGAT !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HTCCGAT !< Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: HTCSGAT !< Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: ILMOGAT !< Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    real, pointer, dimension(:) :: PCFCGAT !< Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCLCGAT !< Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCPGGAT !< Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PCPNGAT !< Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PETGAT  !< Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QEVPGAT !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: QFCFGAT !< Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFCLGAT !< Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFGGAT  !< Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFNGAT  !< Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFSGAT  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: QFXGAT  !< Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: QGGAT   !< Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: ROFGAT  !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFBGAT !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFCGAT !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFNGAT !< Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFOGAT !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROFSGAT !< Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: ROVGGAT !< Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: SFCQGAT !< Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    real, pointer, dimension(:) :: SFCTGAT !< Diagnosed screen-level air temperature [K]
    real, pointer, dimension(:) :: SFCUGAT !< Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: SFCVGAT !< Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: TFXGAT  !< Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    real, pointer, dimension(:) :: TROBGAT !< Temperature of base flow from bottom of soil column [K]
    real, pointer, dimension(:) :: TROFGAT !< Temperature of total runoff [K]
    real, pointer, dimension(:) :: TROOGAT !< Temperature of overland flow from top of soil column [K]
    real, pointer, dimension(:) :: TROSGAT !< Temperature of interflow from sides of soil column [K]
    real, pointer, dimension(:) :: UEGAT   !< Friction velocity of air \f$[m s^{-1} ]\f$
    real, pointer, dimension(:) :: WTABGAT !< Depth of water table in soil [m]
    real, pointer, dimension(:) :: WTRCGAT !< Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: WTRGGAT !< Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: WTRSGAT !< Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: wtableGAT !< Depth of water table in soil [m]
    real, pointer, dimension(:) :: maxAnnualActLyrGAT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, pointer, dimension(:) :: QLWOGAT !<
    real, pointer, dimension(:) :: SFRHGAT !<
    real, pointer, dimension(:) :: FTEMP   !<
    real, pointer, dimension(:) :: FVAP    !<
    real, pointer, dimension(:) :: RIB     !<
    real, pointer, dimension(:) :: FC      !<
    real, pointer, dimension(:) :: FG      !<
    real, pointer, dimension(:) :: FCS     !<
    real, pointer, dimension(:) :: FGS     !<
    real, pointer, dimension(:) :: RBCOEF  !<
    real, pointer, dimension(:) :: ZSNOW   !< Depth of snow pack \f$[m] (z_s)\f$
    real, pointer, dimension(:) :: FSVF    !<
    real, pointer, dimension(:) :: FSVFS   !<
    real, pointer, dimension(:) :: ALVSCN  !<
    real, pointer, dimension(:) :: ALIRCN  !<
    real, pointer, dimension(:) :: ALVSG   !<
    real, pointer, dimension(:) :: ALIRG   !<
    real, pointer, dimension(:) :: ALVSCS  !<
    real, pointer, dimension(:) :: ALIRCS  !<
    real, pointer, dimension(:) :: ALVSSN  !<
    real, pointer, dimension(:) :: ALIRSN  !<
    real, pointer, dimension(:) :: ALVSGC  !<
    real, pointer, dimension(:) :: ALIRGC  !<
    real, pointer, dimension(:) :: ALVSSC  !<
    real, pointer, dimension(:) :: ALIRSC  !<
    real, pointer, dimension(:) :: TRVSCN  !<
    real, pointer, dimension(:) :: TRIRCN  !<
    real, pointer, dimension(:) :: TRVSCS  !<
    real, pointer, dimension(:) :: TRIRCS  !<
    real, pointer, dimension(:) :: RC      !<
    real, pointer, dimension(:) :: RCS     !<
    real, pointer, dimension(:) :: FRAINC  !<
    real, pointer, dimension(:) :: FSNOWC  !<
    real, pointer, dimension(:) :: FRAICS  !<
    real, pointer, dimension(:) :: FSNOCS  !<
    real, pointer, dimension(:) :: CMASSC  !<
    real, pointer, dimension(:) :: CMASCS  !<
    real, pointer, dimension(:) :: DISP    !<
    real, pointer, dimension(:) :: DISPS   !<
    real, pointer, dimension(:) :: ZOMLNC  !<
    real, pointer, dimension(:) :: ZOELNC  !<
    real, pointer, dimension(:) :: ZOMLNG  !<
    real, pointer, dimension(:) :: ZOELNG  !<
    real, pointer, dimension(:) :: ZOMLCS  !<
    real, pointer, dimension(:) :: ZOELCS  !<
    real, pointer, dimension(:) :: ZOMLNS  !<
    real, pointer, dimension(:) :: ZOELNS  !<
    real, pointer, dimension(:) :: TRSNOWC !<
    real, pointer, dimension(:) :: CHCAP   !< Heat capacity of vegetation canopy \f$[J m^{-2} K^{-1} ] (C_c)\f$
    real, pointer, dimension(:) :: CHCAPS  !<
    real, pointer, dimension(:) :: GZEROC  !< Vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: GZEROG  !< Bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: GZROCS  !< Snow-covered vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: GZROGS  !< Snow-covered bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: groundHeatFlux !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: G12C    !<
    real, pointer, dimension(:) :: G12G    !<
    real, pointer, dimension(:) :: G12CS   !<
    real, pointer, dimension(:) :: G12GS   !<
    real, pointer, dimension(:) :: G23C    !<
    real, pointer, dimension(:) :: G23G    !<
    real, pointer, dimension(:) :: G23CS   !<
    real, pointer, dimension(:) :: G23GS   !<
    real, pointer, dimension(:) :: QFREZC  !<
    real, pointer, dimension(:) :: QFREZG  !<
    real, pointer, dimension(:) :: QMELTC  !<
    real, pointer, dimension(:) :: QMELTG  !<
    real, pointer, dimension(:) :: EVAPC   !<
    real, pointer, dimension(:) :: EVAPCG  !<
    real, pointer, dimension(:) :: EVAPG   !<
    real, pointer, dimension(:) :: EVAPCS  !<
    real, pointer, dimension(:) :: EVPCSG  !<
    real, pointer, dimension(:) :: EVAPGS  !<
    real, pointer, dimension(:) :: TCANO   !<
    real, pointer, dimension(:) :: TCANS   !<
    real, pointer, dimension(:) :: RAICAN  !<
    real, pointer, dimension(:) :: SNOCAN  !<
    real, pointer, dimension(:) :: RAICNS  !<
    real, pointer, dimension(:) :: SNOCNS  !<
    real, pointer, dimension(:) :: CWLCAP  !<
    real, pointer, dimension(:) :: CWFCAP  !<
    real, pointer, dimension(:) :: CWLCPS  !<
    real, pointer, dimension(:) :: CWFCPS  !<
    real, pointer, dimension(:) :: TSNOCS  !<
    real, pointer, dimension(:) :: TSNOGS  !<
    real, pointer, dimension(:) :: RHOSCS  !<
    real, pointer, dimension(:) :: RHOSGS  !<
    real, pointer, dimension(:) :: WSNOCS  !<
    real, pointer, dimension(:) :: WSNOGS  !<
    real, pointer, dimension(:) :: TPONDC  !<
    real, pointer, dimension(:) :: TPONDG  !<
    real, pointer, dimension(:) :: TPNDCS  !<
    real, pointer, dimension(:) :: TPNDGS  !<
    real, pointer, dimension(:) :: ZPLMCS  !<
    real, pointer, dimension(:) :: ZPLMGS  !<
    real, pointer, dimension(:) :: ZPLIMC  !<
    real, pointer, dimension(:) :: ZPLIMG  !<
    !
    !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
    !     * BALANCES.
    !
    real, pointer, dimension(:) :: CTVSTP !<
    real, pointer, dimension(:) :: CTSSTP !<
    real, pointer, dimension(:) :: CT1STP !<
    real, pointer, dimension(:) :: CT2STP !<
    real, pointer, dimension(:) :: CT3STP !<
    real, pointer, dimension(:) :: WTVSTP !<
    real, pointer, dimension(:) :: WTSSTP !<
    real, pointer, dimension(:) :: WTGSTP !<

    ! These will be allocated the dimension: 'ilg,ignd'
    integer, pointer, dimension(:,:) :: ISNDGAT !< Integer identifier associated with sand content
    real, pointer, dimension(:,:) :: TBARGAT !< Temperature of soil layers [K]
    real, pointer, dimension(:,:) :: THICGAT !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THLQGAT !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: BIGAT   !< Clapp and Hornberger empirical “b” parameter [ ]
    real, pointer, dimension(:,:) :: DLZWGAT !< Permeable thickness of soil layer [m]
    real, pointer, dimension(:,:) :: GRKSGAT !< Saturated hydraulic conductivity of soil layers \f$[m s^{-1} ]\f$
    real, pointer, dimension(:,:) :: HCPSGAT !< Volumetric heat capacity of soil particles \f$[J m^{-3} ]\f$
    real, pointer, dimension(:,:) :: PSISGAT !< Soil moisture suction at saturation [m]
    real, pointer, dimension(:,:) :: PSIWGAT !< Soil moisture suction at wilting point [m]
    real, pointer, dimension(:,:) :: TCSGAT  !< Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$\
    real, pointer, dimension(:,:) :: THFCGAT !< Field capacity \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THMGAT  !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THPGAT  !< Pore volume in soil layer \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THRGAT  !< Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THRAGAT !< Fractional saturation of soil behind the wetting front [ ]
    real, pointer, dimension(:,:) :: ZBTWGAT !< Depth to permeable bottom of soil layer [m]
    real, pointer, dimension(:,:) :: THLWGAT !< Soil water content at wilting point, \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: GFLXGAT !< Heat conduction between soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HMFGGAT !< Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HTCGAT  !< Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: QFCGAT  !< Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: TBARC  !< Temperature of soil layers for ground under canopy subarea [K]
    real, pointer, dimension(:,:) :: TBARG  !< Temperature of soil layers for bareground subarea [K]
    real, pointer, dimension(:,:) :: TBARCS !< Temperature of soil layers for snow-covered ground under canopy subarea [K]
    real, pointer, dimension(:,:) :: TBARGS !< Temperature of soil layers for bareground subarea [K]
    real, pointer, dimension(:,:) :: THLIQC !< Volumetric liquid water content of soil layers for ground under canopy subarea \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THLIQG !< Volumetric liquid water content of soil layers for bareground subarea \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THICEC !< Volumetric ice content of soil layers for ground under canopy subarea \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: THICEG !< Volumetric ice content of soil layers for bareground subarea \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:) :: FROOT  !<
    real, pointer, dimension(:,:) :: HCPC   !<
    real, pointer, dimension(:,:) :: HCPG   !<
    real, pointer, dimension(:,:) :: FROOTS !<
    real, pointer, dimension(:,:) :: TCTOPC !<
    real, pointer, dimension(:,:) :: TCBOTC !<
    real, pointer, dimension(:,:) :: TCTOPG !<
    real, pointer, dimension(:,:) :: TCBOTG !<

    ! These will be allocated the dimension: 'ilg,ican'
    real, pointer, dimension(:,:) :: ACIDGAT !< Optional user-specified value of canopy near-infrared albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:,:) :: ACVDGAT !< Optional user-specified value of canopy visible albedo to override CLASS-calculated value [ ]
    real, pointer, dimension(:,:) :: CMASGAT !< Maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$
    real, pointer, dimension(:,:) :: HGTDGAT !< Optional user-specified values of height of vegetation categories to override CLASS-calculated values [m]
    real, pointer, dimension(:,:) :: PAIDGAT !< Optional user-specified value of plant area indices of vegetation categories to override CLASS-calculated values [ ]
    real, pointer, dimension(:,:) :: PAMNGAT !< Minimum plant area index of vegetation category [ ]
    real, pointer, dimension(:,:) :: PAMXGAT !< Minimum plant area index of vegetation category [ ]
    real, pointer, dimension(:,:) :: PSGAGAT !< Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, pointer, dimension(:,:) :: PSGBGAT !< Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, pointer, dimension(:,:) :: QA50GAT !< Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$
    real, pointer, dimension(:,:) :: ROOTGAT !< Maximum rooting depth of vegetation category [m]
    real, pointer, dimension(:,:) :: RSMNGAT !< Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$
    real, pointer, dimension(:,:) :: VPDAGAT !< Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, pointer, dimension(:,:) :: VPDBGAT !< Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]


    ! These will be allocated the dimension: 'ilg,icp1'
    real, pointer, dimension(:,:) :: ALICGAT !< Background average near-infrared albedo of vegetation category [ ]
    real, pointer, dimension(:,:) :: ALVCGAT !< Background average visible albedo of vegetation category [ ]
    real, pointer, dimension(:,:) :: FCANGAT !< Maximum fractional coverage of modelled area by vegetation category [ ]
    real, pointer, dimension(:,:) :: LNZ0GAT !< Natural logarithm of maximum roughness length of vegetation category [ ]

    ! These will be allocated the dimension: 'ilg,nbs'
    real, pointer, dimension(:,:) :: FSDBGAT !<
    real, pointer, dimension(:,:) :: FSFBGAT !<
    real, pointer, dimension(:,:) :: FSSBGAT !<
    real, pointer, dimension(:,:) :: SALBGAT !<
    real, pointer, dimension(:,:) :: CSALGAT !<
    real, pointer, dimension(:,:) :: ALTG    !<
    real, pointer, dimension(:,:) :: ALSNO   !<
    real, pointer, dimension(:,:) :: TRSNOWG !<

    ! These will be allocated the dimension: 'ilg,4'
    real, pointer, dimension(:,:) :: TSFSGAT !< Ground surface temperature over subarea [K]

    ! These will be allocated the dimension: 'ilg,6,50'
    integer, pointer, dimension(:,:,:) :: ITCTGAT !< Counter of number of iterations required to solve surface energy balance for the elements of the four subareas

    ! The following are stored in the data structure: class_rot
    ! they are allocatted in allocClassVars in the classStateVars
    ! module and pointed to here.

    ! These will be allocated the dimension: 'nlat'

    real, pointer, dimension(:) :: CSZROW  !< Cosine of solar zenith angle [ ]
    real, pointer, dimension(:) :: DLONROW !<
    real, pointer, dimension(:) :: DLATROW !<
    real, pointer, dimension(:) :: FCLOROW !< Fractional cloud cover [ ]
    real, pointer, dimension(:) :: FDLROW  !<
    real, pointer, dimension(:) :: FSIHROW !<
    real, pointer, dimension(:) :: FSVHROW !<
    real, pointer, dimension(:) :: GCROW   !< Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
    real, pointer, dimension(:) :: GGEOROW !< The geothermal heat flux
    real, pointer, dimension(:) :: PADRROW !<
    real, pointer, dimension(:) :: PREROW  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: PRESROW !<
    real, pointer, dimension(:) :: QAROW   !<
    real, pointer, dimension(:) :: RADJROW !< Latitude of grid cell (positive north of equator) [rad]
    real, pointer, dimension(:) :: RHOAROW !<
    real, pointer, dimension(:) :: RHSIROW !<
    real, pointer, dimension(:) :: RPCPROW !<
    real, pointer, dimension(:) :: RPREROW !< Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: SPCPROW !<
    real, pointer, dimension(:) :: SPREROW !< Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, pointer, dimension(:) :: TAROW   !<
    real, pointer, dimension(:) :: TADPROW !<
    real, pointer, dimension(:) :: TRPCROW !<
    real, pointer, dimension(:) :: TSPCROW !<
    real, pointer, dimension(:) :: ULROW   !<
    real, pointer, dimension(:) :: VLROW   !<
    real, pointer, dimension(:) :: VMODROW !<
    real, pointer, dimension(:) :: VPDROW  !<
    real, pointer, dimension(:) :: ZBLDROW !<
    real, pointer, dimension(:) :: ZDHROW  !<
    real, pointer, dimension(:) :: ZDMROW  !<
    real, pointer, dimension(:) :: ZRFHROW !<
    real, pointer, dimension(:) :: ZRFMROW !<
    real, pointer, dimension(:) :: UVROW   !<
    real, pointer, dimension(:) :: XDIFFUS !< Fraction of diffused radiation
    real, pointer, dimension(:) :: Z0ORROW !< The orographic roughness length
    real, pointer, dimension(:) :: FSSROW  !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, pointer, dimension(:) :: PRENROW !<
    real, pointer, dimension(:) :: CLDTROW !<
    real, pointer, dimension(:) :: FSGROL  !<
    real, pointer, dimension(:) :: FLGROL  !<
    real, pointer, dimension(:) :: GUSTROL !<
    real, pointer, dimension(:) :: DEPBROW !<
    real, pointer, dimension(:) :: WTABROW !<  FLAG make so only one of WTABROW and wtableROW ! JM.
    real, pointer, dimension(:) :: wtableROW !< Depth of water table in soil [m]

    ! These will be allocated the dimension: 'nlat,nmos'

    integer, pointer, dimension(:,:) :: IGDRROT !<
    integer, pointer, dimension(:,:) :: MIDROT  !< Mosaic tile type identifier (1 for land surface, 0 for inland lake)
    real, pointer, dimension(:,:) :: ALBSROT !<
    real, pointer, dimension(:,:) :: CMAIROT !<
    real, pointer, dimension(:,:) :: GROROT  !<
    real, pointer, dimension(:,:) :: QACROT  !<
    real, pointer, dimension(:,:) :: RCANROT !<
    real, pointer, dimension(:,:) :: RHOSROT !<
    real, pointer, dimension(:,:) :: SCANROT !<
    real, pointer, dimension(:,:) :: SNOROT  !<
    real, pointer, dimension(:,:) :: TACROT  !<
    real, pointer, dimension(:,:) :: TBASROT !<
    real, pointer, dimension(:,:) :: TCANROT !<
    real, pointer, dimension(:,:) :: TPNDROT !<
    real, pointer, dimension(:,:) :: TSNOROT !<
    real, pointer, dimension(:,:) :: WSNOROT !<
    real, pointer, dimension(:,:) :: ZPNDROT !<
    real, pointer, dimension(:,:) :: REFROT  !< Snow grain size (for ISNOALB=1 option)  [m]
    real, pointer, dimension(:,:) :: BCSNROT !<
    real, pointer, dimension(:,:) :: AGIDROT !<
    real, pointer, dimension(:,:) :: AGVDROT !<
    real, pointer, dimension(:,:) :: ALGDROT !<
    real, pointer, dimension(:,:) :: ALGWROT !<
    real, pointer, dimension(:,:) :: ASIDROT !<
    real, pointer, dimension(:,:) :: ASVDROT !<
    real, pointer, dimension(:,:) :: DRNROT  !<
    real, pointer, dimension(:,:) :: FAREROT !< Fractional coverage of mosaic tile on modelled area
    real, pointer, dimension(:,:) :: GRKFROT !<
    real, pointer, dimension(:,:) :: WFCIROT !<
    real, pointer, dimension(:,:) :: WFSFROT !<
    real, pointer, dimension(:,:) :: XSLPROT !<
    real, pointer, dimension(:,:) :: ZPLGROT !<
    real, pointer, dimension(:,:) :: ZPLSROT !<
    real, pointer, dimension(:,:) :: ZSNLROT !< Limiting snow depth (m)
    real, pointer, dimension(:,:) :: ZSNOROT  !<
    real, pointer, dimension(:,:) :: ALGWVROT !<
    real, pointer, dimension(:,:) :: ALGWNROT !<
    real, pointer, dimension(:,:) :: ALGDVROT !<
    real, pointer, dimension(:,:) :: ALGDNROT !<
    real, pointer, dimension(:,:) :: EMISROT  !<
    real, pointer, dimension(:,:) :: ALIRROT !<
    real, pointer, dimension(:,:) :: ALVSROT !<
    real, pointer, dimension(:,:) :: CDHROT  !<
    real, pointer, dimension(:,:) :: CDMROT  !<
    real, pointer, dimension(:,:) :: DRROT   !<
    real, pointer, dimension(:,:) :: EFROT   !<
    real, pointer, dimension(:,:) :: FLGGROT !<
    real, pointer, dimension(:,:) :: FLGSROT !<
    real, pointer, dimension(:,:) :: FLGVROT !<
    real, pointer, dimension(:,:) :: FSGGROT !<
    real, pointer, dimension(:,:) :: FSGSROT !<
    real, pointer, dimension(:,:) :: FSGVROT !<
    real, pointer, dimension(:,:) :: FSNOROT !<
    real, pointer, dimension(:,:) :: GAROT   !<
    real, pointer, dimension(:,:) :: GTROT   !< Diagnosed effective surface black-body temperature [K]
    real, pointer, dimension(:,:) :: HBLROT  !<
    real, pointer, dimension(:,:) :: HEVCROT !<
    real, pointer, dimension(:,:) :: HEVGROT !<
    real, pointer, dimension(:,:) :: HEVSROT !<
    real, pointer, dimension(:,:) :: HFSROT  !<
    real, pointer, dimension(:,:) :: HFSCROT !<
    real, pointer, dimension(:,:) :: HFSGROT !<
    real, pointer, dimension(:,:) :: HFSSROT !<
    real, pointer, dimension(:,:) :: HMFCROT !<
    real, pointer, dimension(:,:) :: HMFNROT !<
    real, pointer, dimension(:,:) :: HTCCROT !<
    real, pointer, dimension(:,:) :: SDEPROT !< Depth to bedrock in the soil profile
    real, pointer, dimension(:,:) :: SOCIROT  !<
    real, pointer, dimension(:,:) :: HTCSROT !<
    real, pointer, dimension(:,:) :: ILMOROT !<
    real, pointer, dimension(:,:) :: PCFCROT !<
    real, pointer, dimension(:,:) :: PCLCROT !<
    real, pointer, dimension(:,:) :: PCPGROT !<
    real, pointer, dimension(:,:) :: PCPNROT !<
    real, pointer, dimension(:,:) :: PETROT  !<
    real, pointer, dimension(:,:) :: QEVPROT !<
    real, pointer, dimension(:,:) :: QFCFROT !<
    real, pointer, dimension(:,:) :: QFCLROT !<
    real, pointer, dimension(:,:) :: QFGROT  !<
    real, pointer, dimension(:,:) :: QFNROT  !<
    real, pointer, dimension(:,:) :: QFSROT  !<
    real, pointer, dimension(:,:) :: QFXROT  !<
    real, pointer, dimension(:,:) :: QGROT   !<
    real, pointer, dimension(:,:) :: ROFROT  !<
    real, pointer, dimension(:,:) :: ROFBROT !<
    real, pointer, dimension(:,:) :: ROFCROT !<
    real, pointer, dimension(:,:) :: ROFNROT !<
    real, pointer, dimension(:,:) :: ROFOROT !<
    real, pointer, dimension(:,:) :: ROFSROT !<
    real, pointer, dimension(:,:) :: ROVGROT !<
    real, pointer, dimension(:,:) :: SFCQROT !<
    real, pointer, dimension(:,:) :: SFCTROT !<
    real, pointer, dimension(:,:) :: SFCUROT !<
    real, pointer, dimension(:,:) :: SFCVROT !<
    real, pointer, dimension(:,:) :: TFXROT  !<
    real, pointer, dimension(:,:) :: TROBROT !<
    real, pointer, dimension(:,:) :: TROFROT !<
    real, pointer, dimension(:,:) :: TROOROT !<
    real, pointer, dimension(:,:) :: TROSROT !<
    real, pointer, dimension(:,:) :: UEROT   !<
    real, pointer, dimension(:,:) :: WTABROT !<
    real, pointer, dimension(:,:) :: WTRCROT !<
    real, pointer, dimension(:,:) :: WTRGROT !<
    real, pointer, dimension(:,:) :: WTRSROT !<
    real, pointer, dimension(:,:) :: SFRHROT !<
    real, pointer, dimension(:,:) :: wtableROT !< Depth of water table in soil [m]
    real, pointer, dimension(:,:) :: maxAnnualActLyrROT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, pointer, dimension(:,:) :: groundHeatFluxROT !< Heat flux at soil surface \f$[W m^{-2} ]\f$

    ! There will be allocated the dimension: 'nlat,nmos,ignd'
    integer, pointer, dimension(:,:,:) :: ISNDROT !<
    real, pointer, dimension(:,:,:) :: TBARROT !<
    real, pointer, dimension(:,:,:) :: THICROT !<
    real, pointer, dimension(:,:,:) :: THLQROT !<
    real, pointer, dimension(:,:,:) :: BIROT   !<
    real, pointer, dimension(:,:,:) :: DLZWROT !< Permeable thickness of soil layer [m]
    real, pointer, dimension(:,:,:) :: GRKSROT !<
    real, pointer, dimension(:,:,:) :: HCPSROT !<
    real, pointer, dimension(:,:,:) :: SANDROT !< Percentage sand content of soil
    real, pointer, dimension(:,:,:) :: CLAYROT !< Percentage clay content of soil
    real, pointer, dimension(:,:,:) :: ORGMROT !< Percentage organic matter content of soil
    real, pointer, dimension(:,:,:) :: PSISROT !<
    real, pointer, dimension(:,:,:) :: PSIWROT !<
    real, pointer, dimension(:,:,:) :: TCSROT  !<
    real, pointer, dimension(:,:,:) :: THFCROT !<
    real, pointer, dimension(:,:,:) :: THMROT  !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, pointer, dimension(:,:,:) :: THPROT  !<
    real, pointer, dimension(:,:,:) :: THRROT  !<
    real, pointer, dimension(:,:,:) :: THRAROT !<
    real, pointer, dimension(:,:,:) :: ZBTWROT !<
    real, pointer, dimension(:,:,:) :: THLWROT !<
    real, pointer, dimension(:,:,:) :: GFLXROT !<
    real, pointer, dimension(:,:,:) :: HMFGROT !<
    real, pointer, dimension(:,:,:) :: HTCROT  !<
    real, pointer, dimension(:,:,:) :: QFCROT  !<

    ! These will be allocated the dimension: 'nlat,nmos,ican'
    real, pointer, dimension(:,:,:) :: ACIDROT !<
    real, pointer, dimension(:,:,:) :: ACVDROT !<
    real, pointer, dimension(:,:,:) :: CMASROT !<
    real, pointer, dimension(:,:,:) :: HGTDROT !<
    real, pointer, dimension(:,:,:) :: PAIDROT !<
    real, pointer, dimension(:,:,:) :: PAMNROT !<
    real, pointer, dimension(:,:,:) :: PAMXROT !<
    real, pointer, dimension(:,:,:) :: PSGAROT !<
    real, pointer, dimension(:,:,:) :: PSGBROT !<
    real, pointer, dimension(:,:,:) :: QA50ROT !<
    real, pointer, dimension(:,:,:) :: ROOTROT !<
    real, pointer, dimension(:,:,:) :: RSMNROT !<
    real, pointer, dimension(:,:,:) :: VPDAROT !<
    real, pointer, dimension(:,:,:) :: VPDBROT !<

    ! These will be allocated the dimension: 'nlat,nmos,icp1'
    real, pointer, dimension(:,:,:) :: ALICROT !<
    real, pointer, dimension(:,:,:) :: ALVCROT !<
    real, pointer, dimension(:,:,:) :: FCANROT !<
    real, pointer, dimension(:,:,:) :: LNZ0ROT !<

    ! These will be allocated the dimension: 'nlat,nmos,nbs'
    real, pointer, dimension(:,:,:)  :: SALBROT  !<
    real, pointer, dimension(:,:,:)  :: CSALROT  !<

    ! These will be allocated the dimension: 'nlat,nbs'
    real, pointer, dimension(:,:) :: FSDBROL  !<
    real, pointer, dimension(:,:) :: FSFBROL  !<
    real, pointer, dimension(:,:) :: FSSBROL  !<

    ! These will be allocated the dimension: 'nlat,nmos,6,50'
    integer, pointer, dimension(:,:,:,:) :: ITCTROT !<

    ! These will be allocated the dimension: 'nlat,nmos,4'
    real, pointer, dimension(:,:,:)  :: TSFSROT !<
    !
    real :: CUMSNO

    !================= CTEM array declaration ===============================\
    !
    !     Local variables for coupling CLASS and CTEM
    !
    integer   :: lopcount ! xday, month1, month2,

    integer, pointer :: spinfast !< set this to a higher number up to 10 to spin up
    !! soil carbon pool faster
    integer, pointer :: useTracer !< Switch for use of a model tracer. If useTracer is 0 then the tracer code is not used.
    !! useTracer = 1 turns on a simple tracer that tracks pools and fluxes. The simple tracer then requires that the tracer values in
    !!               the init_file and the tracerCO2file are set to meaningful values for the experiment being run.
    !! useTracer = 2 means the tracer is 14C and will then call a 14C decay scheme.
    !! useTracer = 3 means the tracer is 13C and will then call a 13C fractionation scheme.
    integer, pointer :: metLoop !< no. of times the .met file is to be read. this
    !! option is useful to see how ctem's c pools
    !! equilibrate when driven with same climate data
    !! over and over again.
    integer, pointer :: jhhstd  !< day of the year to start writing the half-hourly output
    integer, pointer :: jhhendd !< day of the year to stop writing the half-hourly output
    integer, pointer :: jdstd   !< day of the year to start writing the daily output
    integer, pointer :: jdendd  !< day of the year to stop writing the daily output
    integer, pointer :: jhhsty  !< simulation year (runyr) to start writing the half-hourly output
    integer, pointer :: jhhendy !< simulation year (runyr) to stop writing the half-hourly output
    integer, pointer :: jdsty   !< simulation year (runyr) to start writing the daily output
    integer, pointer :: jdendy  !< simulation year (runyr) to stop writing the daily output
    integer, pointer :: jmosty    !< Year to start writing out the monthly output files. If you want to write monthly outputs right
    integer, pointer :: fixedYearLUC  !< Set the year to use for land cover if lnduseon is false. If set to -9999,
    !! we use the PFT distribution found in the initialization file. Any other year
    !1 we search for that year in the LUCFile
    integer, pointer :: fixedYearOBSWETF !< set the year to use for observed wetland fraction if transientOBSWETF is false.

    integer, pointer :: idisp    !< if idisp=0, vegetation displacement heights are ignored,
    !! because the atmospheric model considers these to be part
    !! of the "terrain".
    !! if idisp=1, vegetation displacement heights are calculated.
    integer, pointer :: izref    !< if izref=1, the bottom of the atmospheric model is taken
    !! to lie at the ground surface.
    !! if izref=2, the bottom of the atmospheric model is taken
    !! to lie at the local roughness height.
    integer, pointer :: islfd    !< if islfd=0, drcoef is called for surface stability corrections
    !! and the original gcm set of screen-level diagnostic calculations
    !! is done.
    !! if islfd=1, drcoef is called for surface stability corrections
    !! and sldiag is called for screen-level diagnostic calculations.
    !1 if islfd=2, flxsurfz is called for surface stability corrections
    !! and diasurf is called for screen-level diagnostic calculations.
    integer, pointer :: ipcp     !< if ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 c.
    !! if ipcp=2, a linear partitioning of precipitation betweeen
    !! rainfall and snowfall is done between 0 c and 2 c.
    !! if ipcp=3, rainfall and snowfall are partitioned according to
    !! a polynomial curve between 0 c and 6 c.
    integer, pointer :: iwf     !< if iwf=0, only overland flow and baseflow are modelled, and
    !! the ground surface slope is not modelled.
    !! if iwf=n (0<n<4), the watflood calculations of overland flow
    !! and interflow are performed; interflow is drawn from the top
    !! n soil layers.
    integer, pointer :: ITC !< itc, itcg and itg are switches to choose the iteration scheme to
    !! be used in calculating the canopy or ground surface temperature
    !! respectively.  if the switch is set to 1, a bisection method is
    !! used; if to 2, the newton-raphson method is used.
    integer, pointer :: ITCG !< itc, itcg and itg are switches to choose the iteration scheme to
    !! be used in calculating the canopy or ground surface temperature
    !! respectively.  if the switch is set to 1, a bisection method is
    !! used; if to 2, the newton-raphson method is used.
    integer, pointer :: ITG !< itc, itcg and itg are switches to choose the iteration scheme to
    !! be used in calculating the canopy or ground surface temperature
    !! respectively.  if the switch is set to 1, a bisection method is
    !! used; if to 2, the newton-raphson method is used.
    integer, pointer :: IPAI !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !! plant area index, vegetation height, canopy albedo, snow albedo
    !! and soil albedo respectively calculated by class are used.
    !! if any of these switches is set to 1, the value of the
    !! corresponding parameter calculated by class is overridden by
    !! a user-supplied input value.
    integer, pointer :: IHGT !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !! plant area index, vegetation height, canopy albedo, snow albedo
    !! and soil albedo respectively calculated by class are used.
    !! if any of these switches is set to 1, the value of the
    !! corresponding parameter calculated by class is overridden by
    !! a user-supplied input value.
    integer, pointer :: IALC !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
     !! plant area index, vegetation height, canopy albedo, snow albedo
     !! and soil albedo respectively calculated by class are used.
    !! if any of these switches is set to 1, the value of the
    !! corresponding parameter calculated by class is overridden by
    !! a user-supplied input value.
    integer, pointer :: IALS !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
    !! plant area index, vegetation height, canopy albedo, snow albedo
    !! and soil albedo respectively calculated by class are used.
    !! if any of these switches is set to 1, the value of the
    !! corresponding parameter calculated by class is overridden by
    !! a user-supplied input value.
    integer, pointer :: IALG !< if ipai, ihgt, ialc, ials and ialg are zero, the values of
                            !! plant area index, vegetation height, canopy albedo, snow albedo
                            !! and soil albedo respectively calculated by class are used.
                            !! if any of these switches is set to 1, the value of the
                            !! corresponding parameter calculated by class is overridden by
                            !! a user-supplied input value.
    integer, pointer :: isnoalb !< if isnoalb is set to 0, the original two-band snow albedo algorithms are used.
                                !! if it is set to 1, the new four-band routines are used.
    integer, pointer, dimension(:) :: altotcount_ctm ! nlat
    real, pointer, dimension(:,:)  :: todfrac  !(ilg,icc)
    real, pointer, dimension(:)    :: fsinacc_gat !(ilg)
    real, pointer, dimension(:)    :: flutacc_gat !(ilg)
    real, pointer, dimension(:)    :: flinacc_gat !(ilg)
    real, pointer, dimension(:)    :: altotacc_gat !(ilg)
    real, pointer, dimension(:)    :: netrad_gat !(ilg)
    real, pointer, dimension(:)    :: preacc_gat !(ilg)
    real, pointer, dimension(:)    :: sdepgat !(ilg)
    real, pointer, dimension(:,:)  :: sandgat !(ilg,ignd)
    real, pointer, dimension(:,:)  :: claygat !(ilg,ignd)
    real, pointer, dimension(:,:)  :: orgmgat !(ilg,ignd)
    real, pointer, dimension(:)    :: xdiffusgat !(ilg) ! the corresponding ROW is CLASS's XDIFFUS
    real, pointer, dimension(:)    :: faregat !(ilg)   ! the ROT is FAREROT
    !
    ! Model switches:
    logical, pointer :: ctem_on
    logical, pointer :: dofire
    logical, pointer :: PFTCompetition
    logical, pointer :: start_bare
    logical, pointer :: lnduseon
    logical, pointer :: transientCO2
    logical, pointer :: doMethane
    logical, pointer :: transientCH4
    logical, pointer :: transientPOPD
    logical, pointer :: transientLGHT
    logical, pointer :: transientOBSWETF
    logical, pointer :: inibioclim
    logical, pointer :: leap
    logical, pointer :: doAnnualOutput
    logical, pointer :: doMonthOutput
    logical, pointer :: doDayOutput
    logical, pointer :: doHhOutput
    logical, pointer :: projectedGrid    !< True if you have a projected lon lat grid, false if not. Projected grids can only have
                                          !! regions referenced by the indexes, not coordinates, when running a sub-region

    ! ROW vars:
    logical, pointer, dimension(:,:,:) :: pftexistrow
    integer, pointer, dimension(:,:,:) :: colddaysrow
    integer, pointer, dimension(:,:,:) :: lfstatusrow
    integer, pointer, dimension(:,:,:) :: pandaysrow
    real, pointer, dimension(:,:) :: tcanrs
    real, pointer, dimension(:,:) :: tsnors
    real, pointer, dimension(:,:) :: tpndrs
    real, pointer, dimension(:,:,:) :: csum
    real, pointer, dimension(:,:,:) :: tbaraccrow_m

    real, pointer, dimension(:,:,:) :: gleafmasrow        !
    real, pointer, dimension(:,:,:) :: bleafmasrow        !
    real, pointer, dimension(:,:,:) :: stemmassrow        !
    real, pointer, dimension(:,:,:) :: rootmassrow        !
    real, pointer, dimension(:,:,:) :: pstemmassrow       !
    real, pointer, dimension(:,:,:) :: pgleafmassrow      !
    real, pointer, dimension(:,:,:) :: fcancmxrow
    real, pointer, dimension(:,:) :: gavglairow
    real, pointer, dimension(:,:,:) :: zolncrow
    real, pointer, dimension(:,:,:) :: ailcrow
    real, pointer, dimension(:,:,:) :: ailcgrow
    real, pointer, dimension(:,:,:) :: ailcgsrow
    real, pointer, dimension(:,:,:) :: fcancsrow
    real, pointer, dimension(:,:,:) :: fcancrow
    real, pointer, dimension(:,:) :: co2concrow
    real, pointer, dimension(:,:) :: ch4concrow
    real, pointer, dimension(:,:,:) :: co2i1cgrow
    real, pointer, dimension(:,:,:) :: co2i1csrow
    real, pointer, dimension(:,:,:) :: co2i2cgrow
    real, pointer, dimension(:,:,:) :: co2i2csrow
    real, pointer, dimension(:,:,:) :: ancsvegrow
    real, pointer, dimension(:,:,:) :: ancgvegrow
    real, pointer, dimension(:,:,:) :: rmlcsvegrow
    real, pointer, dimension(:,:,:) :: rmlcgvegrow
    real, pointer, dimension(:,:,:) :: slairow
    real, pointer, dimension(:,:,:) :: ailcbrow
    real, pointer, dimension(:,:) :: canresrow
    real, pointer, dimension(:,:,:) :: flhrlossrow

    real, pointer, dimension(:,:,:) :: grwtheffrow
    real, pointer, dimension(:,:,:) :: lystmmasrow
    real, pointer, dimension(:,:,:) :: lyrotmasrow
    real, pointer, dimension(:,:,:) :: tymaxlairow
    real, pointer, dimension(:,:) :: vgbiomasrow
    real, pointer, dimension(:,:) :: gavgltmsrow
    real, pointer, dimension(:,:) :: gavgscmsrow
    real, pointer, dimension(:,:,:) :: stmhrlosrow
    real, pointer, dimension(:,:,:,:) :: rmatcrow
    real, pointer, dimension(:,:,:,:) :: rmatctemrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: litrmassrow
    real, pointer, dimension(:,:,:) :: soilcmasrow
    real, pointer, dimension(:,:,:) :: litresvegrow
    real, pointer, dimension(:,:,:) :: soilcresvegrow
    real, pointer, dimension(:,:,:) :: humiftrsvegrow
    ! real, pointer, dimension(:,:,:,:) :: litrmassrow
    ! real, pointer, dimension(:,:,:,:) :: soilcmasrow
    ! real, pointer, dimension(:,:,:,:) :: litresvegrow
    ! real, pointer, dimension(:,:,:,:) :: soilcresvegrow
    ! real, pointer, dimension(:,:,:,:) :: humiftrsvegrow
    ! COMBAK PERLAY
    real, pointer, dimension(:,:,:) :: vgbiomas_vegrow

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
    real, pointer, dimension(:,:) :: popdinrow
    real, pointer, dimension(:,:,:) :: btermrow
    real, pointer, dimension(:,:) :: ltermrow
    real, pointer, dimension(:,:,:) :: mtermrow

    real, pointer, dimension(:,:) :: extnprobrow
    real, pointer, dimension(:,:) :: prbfrhucrow
    real, pointer, dimension(:) :: dayl_maxrow
    real, pointer, dimension(:) :: daylrow
    real, pointer, dimension(:) :: grclarea

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

    real, pointer, dimension(:,:,:) :: slopefracrow
    real, pointer, dimension(:,:) :: wetfrac_presrow
    real, pointer, dimension(:,:) :: ch4WetSpecrow
    real, pointer, dimension(:,:) :: wetfdynrow
    real, pointer, dimension(:,:) :: ch4WetDynrow
    real, pointer, dimension(:,:) :: ch4soillsrow

    real, pointer, dimension(:,:) :: peatdeprow
    real, pointer, dimension(:,:) :: lucemcomrow
    real, pointer, dimension(:,:) :: lucltrinrow
    real, pointer, dimension(:,:) :: lucsocinrow

    real, pointer, dimension(:,:) :: npprow
    real, pointer, dimension(:,:) :: neprow
    real, pointer, dimension(:,:) :: nbprow
    real, pointer, dimension(:,:) :: gpprow
    real, pointer, dimension(:,:) :: hetroresrow
    real, pointer, dimension(:,:) :: autoresrow
    real, pointer, dimension(:,:) :: soilcresprow
    real, pointer, dimension(:,:) :: rmrow
    real, pointer, dimension(:,:) :: rgrow
    real, pointer, dimension(:,:) :: litresrow
    real, pointer, dimension(:,:) :: socresrow
    real, pointer, dimension(:,:) :: dstcemlsrow
    real, pointer, dimension(:,:) :: litrfallrow
    real, pointer, dimension(:,:) :: humiftrsrow

    real, pointer, dimension(:,:,:) :: gppvegrow
    real, pointer, dimension(:,:,:) :: nepvegrow
    real, pointer, dimension(:,:,:) :: nbpvegrow
    real, pointer, dimension(:,:,:) :: nppvegrow
    real, pointer, dimension(:,:,:) :: hetroresvegrow
    real, pointer, dimension(:,:,:) :: autoresvegrow
    real, pointer, dimension(:,:,:) :: rmlvegaccrow
    real, pointer, dimension(:,:,:) :: rmsvegrow
    real, pointer, dimension(:,:,:) :: rmrvegrow
    real, pointer, dimension(:,:,:) :: rgvegrow
    real, pointer, dimension(:,:,:) :: litrfallvegrow
    real, pointer, dimension(:,:,:) :: rothrlosrow
    real, pointer, dimension(:,:,:) :: pfcancmxrow
    real, pointer, dimension(:,:,:) :: nfcancmxrow
    real, pointer, dimension(:,:,:) :: alvsctmrow
    real, pointer, dimension(:,:,:) :: paicrow
    real, pointer, dimension(:,:,:) :: slaicrow
    real, pointer, dimension(:,:,:) :: alirctmrow
    real, pointer, dimension(:,:) :: cfluxcgrow
    real, pointer, dimension(:,:) :: cfluxcsrow
    real, pointer, dimension(:,:) :: dstcemls3row
    real, pointer, dimension(:,:,:) :: anvegrow
    real, pointer, dimension(:,:,:) :: rmlvegrow

    real, pointer, dimension(:,:) :: twarmmrow
    real, pointer, dimension(:,:) :: tcoldmrow
    real, pointer, dimension(:,:) :: gdd5row
    real, pointer, dimension(:,:) :: aridityrow
    real, pointer, dimension(:,:) :: srplsmonrow
    real, pointer, dimension(:,:) :: defctmonrow
    real, pointer, dimension(:,:) :: anndefctrow
    real, pointer, dimension(:,:) :: annsrplsrow
    real, pointer, dimension(:,:) :: annpcprow
    real, pointer, dimension(:,:) :: dry_season_lengthrow


    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
    ! GAT version:

    logical, pointer, dimension(:,:) :: pftexistgat
    integer, pointer, dimension(:,:) :: colddaysgat
    integer, pointer, dimension(:,:) :: lfstatusgat
    integer, pointer, dimension(:,:) :: pandaysgat
    real, pointer, dimension(:) :: lightng

    real, pointer, dimension(:,:) :: gleafmasgat        !
    real, pointer, dimension(:,:) :: bleafmasgat        !
    real, pointer, dimension(:,:) :: stemmassgat        !
    real, pointer, dimension(:,:) :: rootmassgat        !
    real, pointer, dimension(:,:) :: pstemmassgat       !
    real, pointer, dimension(:,:) :: pgleafmassgat      !
    real, pointer, dimension(:,:) :: fcancmxgat
    real, pointer, dimension(:) :: gavglaigat
    real, pointer, dimension(:,:) :: zolncgat
    real, pointer, dimension(:,:) :: ailcgat
    real, pointer, dimension(:,:) :: ailcggat
    real, pointer, dimension(:,:) :: ailcgsgat
    real, pointer, dimension(:,:) :: fcancsgat
    real, pointer, dimension(:,:) :: fcancgat
    real, pointer, dimension(:) :: co2concgat
    real, pointer, dimension(:) :: ch4concgat
    real, pointer, dimension(:,:) :: co2i1cggat
    real, pointer, dimension(:,:) :: co2i1csgat
    real, pointer, dimension(:,:) :: co2i2cggat
    real, pointer, dimension(:,:) :: co2i2csgat
    real, pointer, dimension(:,:) :: ancsveggat
    real, pointer, dimension(:,:) :: ancgveggat
    real, pointer, dimension(:,:) :: rmlcsveggat
    real, pointer, dimension(:,:) :: rmlcgveggat
    real, pointer, dimension(:,:) :: slaigat
    real, pointer, dimension(:,:) :: ailcbgat
    real, pointer, dimension(:) :: canresgat
    real, pointer, dimension(:,:) :: flhrlossgat

    real, pointer, dimension(:,:) :: grwtheffgat
    real, pointer, dimension(:,:) :: lystmmasgat
    real, pointer, dimension(:,:) :: lyrotmasgat
    real, pointer, dimension(:,:) :: tymaxlaigat
    real, pointer, dimension(:) :: vgbiomasgat
    real, pointer, dimension(:) :: gavgltmsgat
    real, pointer, dimension(:) :: gavgscmsgat
    real, pointer, dimension(:,:) :: stmhrlosgat
    real, pointer, dimension(:,:,:) :: rmatcgat
    real, pointer, dimension(:,:,:) :: rmatctemgat
    ! COMBAK PERLAY
    real, pointer, dimension(:,:) :: litrmassgat
    real, pointer, dimension(:,:) :: soilcmasgat
    real, pointer, dimension(:,:) :: litresveggat
    real, pointer, dimension(:,:) :: soilcresveggat
    real, pointer, dimension(:,:) :: humiftrsveggat
    ! real, pointer, dimension(:,:,:) :: litrmassgat
    ! real, pointer, dimension(:,:,:) :: soilcmasgat
    ! real, pointer, dimension(:,:,:) :: litresveggat
    ! real, pointer, dimension(:,:,:) :: soilcresveggat
    ! real, pointer, dimension(:,:,:) :: humiftrsveggat
    ! COMBAK PERLAY
    real, pointer, dimension(:,:) :: vgbiomas_veggat

    real, pointer, dimension(:,:) :: emit_co2gat
    real, pointer, dimension(:,:) :: emit_cogat
    real, pointer, dimension(:,:) :: emit_ch4gat
    real, pointer, dimension(:,:) :: emit_nmhcgat
    real, pointer, dimension(:,:) :: emit_h2gat
    real, pointer, dimension(:,:) :: emit_noxgat
    real, pointer, dimension(:,:) :: emit_n2ogat
    real, pointer, dimension(:,:) :: emit_pm25gat
    real, pointer, dimension(:,:) :: emit_tpmgat
    real, pointer, dimension(:,:) :: emit_tcgat
    real, pointer, dimension(:,:) :: emit_ocgat
    real, pointer, dimension(:,:) :: emit_bcgat
    real, pointer, dimension(:) :: burnfracgat
    real, pointer, dimension(:,:) :: burnvegfgat
    real, pointer, dimension(:,:) :: smfuncveggat
    real, pointer, dimension(:) :: popdingat
    real, pointer, dimension(:,:) :: btermgat
    real, pointer, dimension(:) :: ltermgat
    real, pointer, dimension(:,:) :: mtermgat
    real, pointer, dimension(:,:) :: glcaemls  !< green leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: blcaemls  !< brown leaf carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: rtcaemls  !< root carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: stcaemls  !< stem carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: ltrcemls  !< litter carbon emission disturbance losses, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: blfltrdt  !< brown leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, pointer, dimension(:,:) :: glfltrdt  !< green leaf litter generated due to disturbance \f$(kg c/m^2)\f$
    real, pointer, dimension(:,:) :: ntchlveg  !< fluxes for each pft: Net change in leaf biomass, u-mol CO2/m2.sec
    real, pointer, dimension(:,:) :: ntchsveg  !< fluxes for each pft: Net change in stem biomass, u-mol CO2/m2.sec
    real, pointer, dimension(:,:) :: ntchrveg  !< fluxes for each pft: Net change in root biomass,
                                              !! the net change is the difference between allocation and
                                              !! autotrophic respiratory fluxes, u-mol CO2/m2.sec
    real, pointer, dimension(:) :: extnprobgat
    real, pointer, dimension(:) :: prbfrhucgat
    real, pointer, dimension(:) :: dayl_maxgat
    real, pointer, dimension(:) :: daylgat

    real, pointer, dimension(:,:) :: bmasveggat
    real, pointer, dimension(:,:) :: cmasvegcgat
    real, pointer, dimension(:,:) :: veghghtgat
    real, pointer, dimension(:,:) :: rootdpthgat
    real, pointer, dimension(:) :: rmlgat
    real, pointer, dimension(:) :: rmsgat
    real, pointer, dimension(:,:) :: tltrleafgat
    real, pointer, dimension(:,:) :: tltrstemgat
    real, pointer, dimension(:,:) :: tltrrootgat
    real, pointer, dimension(:,:) :: leaflitrgat
    real, pointer, dimension(:,:) :: roottempgat
    real, pointer, dimension(:,:) :: afrleafgat
    real, pointer, dimension(:,:) :: afrstemgat
    real, pointer, dimension(:,:) :: afrrootgat
    real, pointer, dimension(:,:) :: wtstatusgat
    real, pointer, dimension(:,:) :: ltstatusgat
    real, pointer, dimension(:) :: rmrgat

    real, pointer, dimension(:,:) :: slopefracgat
    real, pointer, dimension(:) :: wetfrac_presgat
    real, pointer, dimension(:) :: ch4WetSpecgat
    real, pointer, dimension(:) :: wetfdyngat
    real, pointer, dimension(:) :: ch4WetDyngat
    real, pointer, dimension(:) :: ch4soillsgat

    real, pointer, dimension(:) :: lucemcomgat
    real, pointer, dimension(:) :: lucltringat
    real, pointer, dimension(:) :: lucsocingat

    real, pointer, dimension(:) :: nppgat
    real, pointer, dimension(:) :: nepgat
    real, pointer, dimension(:) :: nbpgat
    real, pointer, dimension(:) :: gppgat
    real, pointer, dimension(:) :: hetroresgat
    real, pointer, dimension(:) :: autoresgat
    real, pointer, dimension(:) :: soilcrespgat
    real, pointer, dimension(:) :: rmgat
    real, pointer, dimension(:) :: rggat
    real, pointer, dimension(:) :: litresgat
    real, pointer, dimension(:) :: socresgat
    real, pointer, dimension(:) :: dstcemlsgat
    real, pointer, dimension(:) :: litrfallgat
    real, pointer, dimension(:) :: humiftrsgat

    real, pointer, dimension(:,:) :: gppveggat
    real, pointer, dimension(:,:) :: nepveggat
    real, pointer, dimension(:,:) :: nbpveggat
    real, pointer, dimension(:,:) :: nppveggat
    real, pointer, dimension(:,:) :: hetroresveggat
    real, pointer, dimension(:,:) :: autoresveggat
    real, pointer, dimension(:,:) :: rmlvegaccgat
    real, pointer, dimension(:,:) :: rmsveggat
    real, pointer, dimension(:,:) :: rmrveggat
    real, pointer, dimension(:,:) :: rgveggat
    real, pointer, dimension(:,:) :: litrfallveggat
    real, pointer, dimension(:,:) :: reprocost   !< Cost of making reproductive tissues, only non-zero when NPP is positive (\f$\mu mol CO_2 m^{-2} s^{-1}\f$)

    real, pointer, dimension(:,:) :: rothrlosgat
    real, pointer, dimension(:,:) :: pfcancmxgat
    real, pointer, dimension(:,:) :: nfcancmxgat
    real, pointer, dimension(:,:) :: alvsctmgat
    real, pointer, dimension(:,:) :: paicgat
    real, pointer, dimension(:,:) :: slaicgat
    real, pointer, dimension(:,:) :: alirctmgat
    real, pointer, dimension(:) :: cfluxcggat
    real, pointer, dimension(:) :: cfluxcsgat
    real, pointer, dimension(:) :: dstcemls3gat
    real, pointer, dimension(:,:) :: anveggat
    real, pointer, dimension(:,:) :: rmlveggat

    real, pointer, dimension(:) :: twarmmgat
    real, pointer, dimension(:) :: tcoldmgat
    real, pointer, dimension(:) :: gdd5gat
    real, pointer, dimension(:) :: ariditygat
    real, pointer, dimension(:) :: srplsmongat
    real, pointer, dimension(:) :: defctmongat
    real, pointer, dimension(:) :: anndefctgat
    real, pointer, dimension(:) :: annsrplsgat
    real, pointer, dimension(:) :: annpcpgat
    real, pointer, dimension(:) :: dry_season_lengthgat

    real, pointer, dimension(:) :: tcurm
    real, pointer, dimension(:) :: srpcuryr
    real, pointer, dimension(:) :: dftcuryr
    real, pointer, dimension(:,:) :: tmonth
    real, pointer, dimension(:) :: anpcpcur
    real, pointer, dimension(:) :: anpecur
    real, pointer, dimension(:) :: gdd5cur
    real, pointer, dimension(:) :: surmncur
    real, pointer, dimension(:) :: defmncur
    real, pointer, dimension(:) :: srplscur
    real, pointer, dimension(:) :: defctcur

    real, pointer, dimension(:,:) :: geremortgat
    real, pointer, dimension(:,:) :: intrmortgat
    real, pointer, dimension(:,:) :: lambdagat
    real, pointer, dimension(:,:) :: ccgat
    real, pointer, dimension(:,:) :: mmgat

    !      Outputs
    real, pointer, dimension(:,:) :: qevpacc_m_save

    !     -----------------------
    !      Tile-level variables (denoted by an ending of "_t")
    real, pointer, dimension(:) :: fsnowacc_t
    real, pointer, dimension(:) :: taaccgat_t
    real, pointer, dimension(:) :: uvaccgat_t
    real, pointer, dimension(:) :: vvaccgat_t
    real, pointer, dimension(:,:) :: tbaraccgat_t
    real, pointer, dimension(:,:) :: thliqacc_t
    real, pointer, dimension(:,:) :: thiceacc_t
    real, pointer, dimension(:,:) :: ancgvgac_t
    real, pointer, dimension(:,:) :: rmlcgvga_t

    !============= CTEM array declaration done =============================/

    ! leap year flag (if the switch 'leap' is true, this will be used, otherwise it remains false)
    logical :: leapnow

    integer, pointer, dimension(:,:) :: ipeatlandrow ! This is first set in read_from_ctm.
    integer, pointer, dimension(:) :: ipeatlandgat
    real, pointer, dimension(:) :: peatdepgat
    real, pointer, dimension(:,:) :: anmossrow
    real, pointer, dimension(:) :: anmossgat
    real, pointer, dimension(:,:) :: rmlmossrow
    real, pointer, dimension(:) :: rmlmossgat
    real, pointer, dimension(:,:) :: gppmossrow
    real, pointer, dimension(:) :: gppmossgat
    real, pointer, dimension(:,:) :: nppmossrow
    real, pointer, dimension(:) :: nppmossgat
    real, pointer, dimension(:,:) :: armossrow
    real, pointer, dimension(:) :: armossgat
    real, pointer, dimension(:,:) :: litrmsmossrow
    real, pointer, dimension(:) :: litrmsmossgat
    real, pointer, dimension(:,:) :: Cmossmasrow
    real, pointer, dimension(:) :: Cmossmasgat
    real, pointer, dimension(:,:) :: dmossrow
    real, pointer, dimension(:) :: dmossgat
    real, pointer, dimension(:,:) :: pddrow
    real, pointer, dimension(:) :: pddgat
    real, pointer, dimension(:) :: ancsmoss
    real, pointer, dimension(:) :: angsmoss
    real, pointer, dimension(:) :: ancmoss
    real, pointer, dimension(:) :: angmoss
    real, pointer, dimension(:) :: rmlcsmoss
    real, pointer, dimension(:) :: rmlgsmoss
    real, pointer, dimension(:) :: rmlcmoss
    real, pointer, dimension(:) :: rmlgmoss

    real, pointer, dimension(:) :: anmossac_t
    real, pointer, dimension(:) :: rmlmossac_t
    real, pointer, dimension(:) :: gppmossac_t

    real, pointer, dimension(:,:) :: tracermossCMassrot      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer, dimension(:,:) :: tracermossLitrMassrot   !< Tracer mass in moss litter, \f$kg C/m^2\f$

    real, pointer, dimension(:,:,:) :: tracergLeafMassrot      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:) :: tracerbLeafMassrot      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:) :: tracerstemMassrot       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:) :: tracerrootMassrot       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    ! allocated with nlat, nmos, iccp2, ignd:
    real, pointer, dimension(:,:,:,:) :: tracerlitrMassrot       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:,:) :: tracersoilCMassrot      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: tracerCO2rot     !< Atmopspheric tracer CO2 concentration (units vary)

    ! allocated with ilg, ...:
    real, pointer, dimension(:) :: tracermossCMassgat      !< Tracer mass in moss biomass, \f$kg C/m^2\f$
    real, pointer, dimension(:) :: tracermossLitrMassgat   !< Tracer mass in moss litter, \f$kg C/m^2\f$
    real, pointer, dimension(:) :: tracerCO2gat     !< Atmopspheric tracer CO2 concentration (units vary)

    ! allocated with nlat, nmos, icc:
    real, pointer, dimension(:,:) :: tracergLeafMassgat      !< Tracer mass in the green leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: tracerbLeafMassgat      !< Tracer mass in the brown leaf pool for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: tracerstemMassgat       !< Tracer mass in the stem for each of the CTEM pfts, \f$kg c/m^2\f$
    real, pointer, dimension(:,:) :: tracerrootMassgat       !< Tracer mass in the roots for each of the CTEM pfts, \f$kg c/m^2\f$
    ! allocated with nlat, nmos, iccp2, ignd:
    real, pointer, dimension(:,:,:) :: tracerlitrMassgat       !< Tracer mass in the litter pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$
    real, pointer, dimension(:,:,:) :: tracersoilCMassgat      !< Tracer mass in the soil carbon pool for each of the CTEM pfts + bareground and LUC products, \f$kg c/m^2\f$

    ! Point the CLASS pointers

    ILMOS   => class_gat%ILMOS
    JLMOS   => class_gat%JLMOS
    IWMOS   => class_gat%IWMOS
    JWMOS   => class_gat%JWMOS
    IGDRGAT => class_gat%IGDRGAT
    DELZ    => class_gat%DELZ
    ZBOT    => class_gat%ZBOT
    GCGAT   => class_gat%GCGAT
    TZSGAT  => class_gat%TZSGAT
    PCSNGAT => class_gat%PCSNGAT
    ALBSGAT => class_gat%ALBSGAT
    CMAIGAT => class_gat%CMAIGAT
    GROGAT  => class_gat%GROGAT
    QACGAT  => class_gat%QACGAT
    RCANGAT => class_gat%RCANGAT
    RHOSGAT=> class_gat%RHOSGAT
    SCANGAT=> class_gat%SCANGAT
    SNOGAT=> class_gat%SNOGAT
    TACGAT=> class_gat%TACGAT
    TBASGAT=> class_gat%TBASGAT
    TCANGAT=> class_gat%TCANGAT
    TPNDGAT=> class_gat%TPNDGAT
    TSNOGAT=> class_gat%TSNOGAT
    WSNOGAT=> class_gat%WSNOGAT
    ZPNDGAT=> class_gat%ZPNDGAT
    REFGAT=> class_gat%REFGAT
    BCSNGAT => class_gat%BCSNGAT
    AGIDGAT => class_gat%AGIDGAT
    AGVDGAT => class_gat%AGVDGAT
    ALGDGAT => class_gat%ALGDGAT
    ALGWGAT => class_gat%ALGWGAT
    ASIDGAT => class_gat%ASIDGAT
    ASVDGAT => class_gat%ASVDGAT
    DRNGAT => class_gat%DRNGAT
    GRKFGAT => class_gat%GRKFGAT
    WFCIGAT => class_gat%WFCIGAT
    WFSFGAT => class_gat%WFSFGAT
    XSLPGAT => class_gat%XSLPGAT
    ZPLGGAT => class_gat%ZPLGGAT
    ZPLSGAT => class_gat%ZPLSGAT
    ZSNLGAT => class_gat%ZSNLGAT
    ALGWVGAT => class_gat%ALGWVGAT
    ALGWNGAT => class_gat%ALGWNGAT
    ALGDVGAT => class_gat%ALGDVGAT
    ALGDNGAT => class_gat%ALGDNGAT
    EMISGAT => class_gat%EMISGAT
    CSZGAT => class_gat%CSZGAT
    DLONGAT => class_gat%DLONGAT
    DLATGAT => class_gat%DLATGAT
    FCLOGAT => class_gat%FCLOGAT
    FDLGAT => class_gat%FDLGAT
    FSIHGAT => class_gat%FSIHGAT
    FSVHGAT => class_gat%FSVHGAT
    GGEOGAT => class_gat%GGEOGAT
    PADRGAT => class_gat%PADRGAT
    PREGAT => class_gat%PREGAT
    PRESGAT => class_gat%PRESGAT
    QAGAT => class_gat%QAGAT
    RADJGAT => class_gat%RADJGAT
    RHOAGAT => class_gat%RHOAGAT
    RHSIGAT => class_gat%RHSIGAT
    RPCPGAT => class_gat%RPCPGAT
    SPCPGAT => class_gat%SPCPGAT
    TAGAT => class_gat%TAGAT
    TADPGAT => class_gat%TADPGAT
    TRPCGAT => class_gat%TRPCGAT
    TSPCGAT => class_gat%TSPCGAT
    ULGAT => class_gat%ULGAT
    VLGAT => class_gat%VLGAT
    VMODGAT => class_gat%VMODGAT
    VPDGAT => class_gat%VPDGAT
    Z0ORGAT => class_gat%Z0ORGAT
    ZBLDGAT => class_gat%ZBLDGAT
    ZDHGAT => class_gat%ZDHGAT
    ZDMGAT => class_gat%ZDMGAT
    ZRFHGAT => class_gat%ZRFHGAT
    ZRFMGAT => class_gat%ZRFMGAT
    FSGGAT => class_gat%FSGGAT
    FLGGAT => class_gat%FLGGAT
    GUSTGAT => class_gat%GUSTGAT
    DEPBGAT => class_gat%DEPBGAT
    GTBS => class_gat%GTBS
    SFCUBS => class_gat%SFCUBS
    SFCVBS => class_gat%SFCVBS
    USTARBS => class_gat%USTARBS
    TCSNOW => class_gat%TCSNOW
    GSNOW => class_gat%GSNOW
    ALIRGAT => class_gat%ALIRGAT
    ALVSGAT => class_gat%ALVSGAT
    CDHGAT => class_gat%CDHGAT
    CDMGAT => class_gat%CDMGAT
    DRGAT => class_gat%DRGAT
    EFGAT => class_gat%EFGAT
    FLGGGAT => class_gat%FLGGGAT
    FLGSGAT => class_gat%FLGSGAT
    FLGVGAT => class_gat%FLGVGAT
    FSGGGAT => class_gat%FSGGGAT
    FSGSGAT => class_gat%FSGSGAT
    FSGVGAT => class_gat%FSGVGAT
    FSNOGAT => class_gat%FSNOGAT
    GAGAT => class_gat%GAGAT
    GTGAT => class_gat%GTGAT
    HBLGAT => class_gat%HBLGAT
    HEVCGAT => class_gat%HEVCGAT
    HEVGGAT => class_gat%HEVGGAT
    HEVSGAT => class_gat%HEVSGAT
    HFSGAT => class_gat%HFSGAT
    HFSCGAT => class_gat%HFSCGAT
    HFSGGAT => class_gat%HFSGGAT
    HFSSGAT => class_gat%HFSSGAT
    HMFCGAT => class_gat%HMFCGAT
    HMFNGAT => class_gat%HMFNGAT
    HTCCGAT => class_gat%HTCCGAT
    HTCSGAT => class_gat%HTCSGAT
    ILMOGAT => class_gat%ILMOGAT
    PCFCGAT => class_gat%PCFCGAT
    PCLCGAT => class_gat%PCLCGAT
    PCPGGAT => class_gat%PCPGGAT
    PCPNGAT => class_gat%PCPNGAT
    PETGAT => class_gat%PETGAT
    QEVPGAT => class_gat%QEVPGAT
    QFCFGAT => class_gat%QFCFGAT
    QFCLGAT => class_gat%QFCLGAT
    QFGGAT => class_gat%QFGGAT
    QFNGAT => class_gat%QFNGAT
    QFSGAT => class_gat%QFSGAT
    QFXGAT => class_gat%QFXGAT
    QGGAT => class_gat%QGGAT
    ROFGAT => class_gat%ROFGAT
    ROFBGAT => class_gat%ROFBGAT
    ROFCGAT => class_gat%ROFCGAT
    ROFNGAT => class_gat%ROFNGAT
    ROFOGAT => class_gat%ROFOGAT
    ROFSGAT => class_gat%ROFSGAT
    ROVGGAT => class_gat%ROVGGAT
    SFCQGAT => class_gat%SFCQGAT
    SFCTGAT => class_gat%SFCTGAT
    SFCUGAT => class_gat%SFCUGAT
    SFCVGAT => class_gat%SFCVGAT
    TFXGAT => class_gat%TFXGAT
    TROBGAT => class_gat%TROBGAT
    TROFGAT => class_gat%TROFGAT
    TROOGAT => class_gat%TROOGAT
    TROSGAT => class_gat%TROSGAT
    UEGAT => class_gat%UEGAT
    WTABGAT => class_gat%WTABGAT
    WTRCGAT => class_gat%WTRCGAT
    WTRGGAT => class_gat%WTRGGAT
    WTRSGAT => class_gat%WTRSGAT
    wtableGAT => class_gat%wtableGAT
    maxAnnualActLyrGAT => class_gat%maxAnnualActLyrGAT
    QLWOGAT => class_gat%QLWOGAT
    SFRHGAT => class_gat%SFRHGAT
    FTEMP => class_gat%FTEMP
    FVAP => class_gat%FVAP
    RIB => class_gat%RIB
    FC => class_gat%FC
    FG => class_gat%FG
    FCS => class_gat%FCS
    FGS => class_gat%FGS
    RBCOEF => class_gat%RBCOEF
    ZSNOW => class_gat%ZSNOW
    FSVF => class_gat%FSVF
    FSVFS => class_gat%FSVFS
    ALVSCN => class_gat%ALVSCN
    ALIRCN => class_gat%ALIRCN
    ALVSG => class_gat%ALVSG
    ALIRG => class_gat%ALIRG
    ALVSCS => class_gat%ALVSCS
    ALIRCS => class_gat%ALIRCS
    ALVSSN => class_gat%ALVSSN
    ALIRSN => class_gat%ALIRSN
    ALVSGC => class_gat%ALVSGC
    ALIRGC => class_gat%ALIRGC
    ALVSSC => class_gat%ALVSSC
    ALIRSC => class_gat%ALIRSC
    TRVSCN => class_gat%TRVSCN
    TRIRCN => class_gat%TRIRCN
    TRVSCS => class_gat%TRVSCS
    TRIRCS => class_gat%TRIRCS
    RC => class_gat%RC
    RCS => class_gat%RCS
    FRAINC => class_gat%FRAINC
    FSNOWC => class_gat%FSNOWC
    FRAICS => class_gat%FRAICS
    FSNOCS => class_gat%FSNOCS
    CMASSC => class_gat%CMASSC
    CMASCS => class_gat%CMASCS
    DISP => class_gat%DISP
    DISPS => class_gat%DISPS
    ZOMLNC => class_gat%ZOMLNC
    ZOELNC => class_gat%ZOELNC
    ZOMLNG => class_gat%ZOMLNG
    ZOELNG => class_gat%ZOELNG
    ZOMLCS => class_gat%ZOMLCS
    ZOELCS => class_gat%ZOELCS
    ZOMLNS => class_gat%ZOMLNS
    ZOELNS => class_gat%ZOELNS
    TRSNOWC => class_gat%TRSNOWC
    CHCAP => class_gat%CHCAP
    CHCAPS => class_gat%CHCAPS
    GZEROC => class_gat%GZEROC
    GZEROG => class_gat%GZEROG
    GZROCS => class_gat%GZROCS
    GZROGS => class_gat%GZROGS
    groundHeatFlux => class_gat%groundHeatFlux
    G12C => class_gat%G12C
    G12G => class_gat%G12G
    G12CS => class_gat%G12CS
    G12GS => class_gat%G12GS
    G23C => class_gat%G23C
    G23G => class_gat%G23G
    G23CS => class_gat%G23CS
    G23GS => class_gat%G23GS
    QFREZC => class_gat%QFREZC
    QFREZG => class_gat%QFREZG
    QMELTC => class_gat%QMELTC
    QMELTG => class_gat%QMELTG
    EVAPC => class_gat%EVAPC
    EVAPCG => class_gat%EVAPCG
    EVAPG => class_gat%EVAPG
    EVAPCS => class_gat%EVAPCS
    EVPCSG => class_gat%EVPCSG
    EVAPGS => class_gat%EVAPGS
    TCANO => class_gat%TCANO
    TCANS => class_gat%TCANS
    RAICAN => class_gat%RAICAN
    SNOCAN => class_gat%SNOCAN
    RAICNS => class_gat%RAICNS
    SNOCNS => class_gat%SNOCNS
    CWLCAP => class_gat%CWLCAP
    CWFCAP => class_gat%CWFCAP
    CWLCPS => class_gat%CWLCPS
    CWFCPS => class_gat%CWFCPS
    TSNOCS => class_gat%TSNOCS
    TSNOGS => class_gat%TSNOGS
    RHOSCS => class_gat%RHOSCS
    RHOSGS => class_gat%RHOSGS
    WSNOCS => class_gat%WSNOCS
    WSNOGS => class_gat%WSNOGS
    TPONDC => class_gat%TPONDC
    TPONDG => class_gat%TPONDG
    TPNDCS => class_gat%TPNDCS
    TPNDGS => class_gat%TPNDGS
    ZPLMCS => class_gat%ZPLMCS
    ZPLMGS => class_gat%ZPLMGS
    ZPLIMC => class_gat%ZPLIMC
    ZPLIMG => class_gat%ZPLIMG
    CTVSTP => class_gat%CTVSTP
    CTSSTP => class_gat%CTSSTP
    CT1STP => class_gat%CT1STP
    CT2STP => class_gat%CT2STP
    CT3STP => class_gat%CT3STP
    WTVSTP => class_gat%WTVSTP
    WTSSTP => class_gat%WTSSTP
    WTGSTP => class_gat%WTGSTP
    ISNDGAT => class_gat%ISNDGAT
    TBARGAT => class_gat%TBARGAT
    THICGAT => class_gat%THICGAT
    THLQGAT => class_gat%THLQGAT
    BIGAT => class_gat%BIGAT
    DLZWGAT => class_gat%DLZWGAT
    GRKSGAT => class_gat%GRKSGAT
    HCPSGAT => class_gat%HCPSGAT
    PSISGAT => class_gat%PSISGAT
    PSIWGAT => class_gat%PSIWGAT
    TCSGAT => class_gat%TCSGAT
    THFCGAT => class_gat%THFCGAT
    THMGAT => class_gat%THMGAT
    THPGAT => class_gat%THPGAT
    THRGAT => class_gat%THRGAT
    THRAGAT => class_gat%THRAGAT
    ZBTWGAT => class_gat%ZBTWGAT
    THLWGAT => class_gat%THLWGAT
    GFLXGAT => class_gat%GFLXGAT
    HMFGGAT => class_gat%HMFGGAT
    HTCGAT => class_gat%HTCGAT
    QFCGAT => class_gat%QFCGAT
    TBARC => class_gat%TBARC
    TBARG => class_gat%TBARG
    TBARCS => class_gat%TBARCS
    TBARGS => class_gat%TBARGS
    THLIQC => class_gat%THLIQC
    THLIQG => class_gat%THLIQG
    THICEC => class_gat%THICEC
    THICEG => class_gat%THICEG
    FROOT => class_gat%FROOT
    HCPC => class_gat%HCPC
    HCPG => class_gat%HCPG
    FROOTS => class_gat%FROOTS
    TCTOPC => class_gat%TCTOPC
    TCBOTC => class_gat%TCBOTC
    TCTOPG => class_gat%TCTOPG
    TCBOTG => class_gat%TCBOTG
    ACIDGAT => class_gat%ACIDGAT
    ACVDGAT => class_gat%ACVDGAT
    CMASGAT => class_gat%CMASGAT
    HGTDGAT => class_gat%HGTDGAT
    PAIDGAT => class_gat%PAIDGAT
    PAMNGAT => class_gat%PAMNGAT
    PAMXGAT => class_gat%PAMXGAT
    PSGAGAT => class_gat%PSGAGAT
    PSGBGAT => class_gat%PSGBGAT
    QA50GAT => class_gat%QA50GAT
    ROOTGAT => class_gat%ROOTGAT
    RSMNGAT => class_gat%RSMNGAT
    VPDAGAT => class_gat%VPDAGAT
    VPDBGAT => class_gat%VPDBGAT
    ALICGAT => class_gat%ALICGAT
    ALVCGAT => class_gat%ALVCGAT
    FCANGAT => class_gat%FCANGAT
    LNZ0GAT => class_gat%LNZ0GAT
    FSDBGAT => class_gat%FSDBGAT
    FSFBGAT => class_gat%FSFBGAT
    FSSBGAT => class_gat%FSSBGAT
    SALBGAT => class_gat%SALBGAT
    CSALGAT => class_gat%CSALGAT
    ALTG => class_gat%ALTG
    ALSNO => class_gat%ALSNO
    TRSNOWG => class_gat%TRSNOWG
    TSFSGAT => class_gat%TSFSGAT
    ITCTGAT => class_gat%ITCTGAT

    CSZROW => class_rot%CSZROW
    DLONROW => class_rot%DLONROW
    DLATROW => class_rot%DLATROW
    FCLOROW => class_rot%FCLOROW
    FDLROW => class_rot%FDLROW
    FSIHROW => class_rot%FSIHROW
    FSVHROW => class_rot%FSVHROW
    GCROW => class_rot%GCROW
    GGEOROW => class_rot%GGEOROW
    PADRROW => class_rot%PADRROW
    PREROW => class_rot%PREROW
    PRESROW => class_rot%PRESROW
    QAROW => class_rot%QAROW
    RADJROW => class_rot%RADJROW
    RHOAROW => class_rot%RHOAROW
    RHSIROW => class_rot%RHSIROW
    RPCPROW => class_rot%RPCPROW
    RPREROW => class_rot%RPREROW
    SPCPROW => class_rot%SPCPROW
    SPREROW => class_rot%SPREROW
    TAROW => class_rot%TAROW
    TADPROW => class_rot%TADPROW
    TRPCROW => class_rot%TRPCROW
    TSPCROW => class_rot%TSPCROW
    ULROW => class_rot%ULROW
    VLROW => class_rot%VLROW
    VMODROW => class_rot%VMODROW
    VPDROW => class_rot%VPDROW
    ZBLDROW => class_rot%ZBLDROW
    ZDHROW => class_rot%ZDHROW
    ZDMROW => class_rot%ZDMROW
    ZRFHROW => class_rot%ZRFHROW
    ZRFMROW => class_rot%ZRFMROW
    UVROW => class_rot%UVROW
    XDIFFUS => class_rot%XDIFFUS
    Z0ORROW => class_rot%Z0ORROW
    FSSROW => class_rot%FSSROW
    PRENROW => class_rot%PRENROW
    CLDTROW => class_rot%CLDTROW
    FSGROL => class_rot%FSGROL
    FLGROL => class_rot%FLGROL
    GUSTROL => class_rot%GUSTROL
    DEPBROW => class_rot%DEPBROW
    WTABROW => class_rot%WTABROW  ! FLAG
    wtableROW => class_rot%wtableROW  ! FLAG
    maxAnnualActLyrROT => class_rot%maxAnnualActLyrROT
    groundHeatFluxROT => class_rot%groundHeatFluxROT
    IGDRROT => class_rot%IGDRROT
    MIDROT => class_rot%MIDROT
    ALBSROT => class_rot%ALBSROT
    CMAIROT => class_rot%CMAIROT
    GROROT => class_rot%GROROT
    QACROT => class_rot%QACROT
    RCANROT => class_rot%RCANROT
    RHOSROT => class_rot%RHOSROT
    SCANROT => class_rot%SCANROT
    SNOROT => class_rot%SNOROT
    TACROT => class_rot%TACROT
    TBASROT => class_rot%TBASROT
    TCANROT => class_rot%TCANROT
    TPNDROT => class_rot%TPNDROT
    TSNOROT => class_rot%TSNOROT
    WSNOROT => class_rot%WSNOROT
    ZPNDROT => class_rot%ZPNDROT
    REFROT => class_rot%REFROT
    BCSNROT => class_rot%BCSNROT
    AGIDROT => class_rot%AGIDROT
    AGVDROT => class_rot%AGVDROT
    ALGDROT => class_rot%ALGDROT
    ALGWROT => class_rot%ALGWROT
    ASIDROT => class_rot%ASIDROT
    ASVDROT => class_rot%ASVDROT
    DRNROT => class_rot%DRNROT
    FAREROT => class_rot%FAREROT
    GRKFROT => class_rot%GRKFROT
    WFCIROT => class_rot%WFCIROT
    WFSFROT => class_rot%WFSFROT
    XSLPROT => class_rot%XSLPROT
    ZPLGROT => class_rot%ZPLGROT
    ZPLSROT => class_rot%ZPLSROT
    ZSNLROT => class_rot%ZSNLROT
    ZSNOROT => class_rot%ZSNOROT
    ALGWVROT => class_rot%ALGWVROT
    ALGWNROT => class_rot%ALGWNROT
    ALGDVROT => class_rot%ALGDVROT
    ALGDNROT => class_rot%ALGDNROT
    EMISROT => class_rot%EMISROT
    ALIRROT => class_rot%ALIRROT
    ALVSROT => class_rot%ALVSROT
    CDHROT => class_rot%CDHROT
    CDMROT => class_rot%CDMROT
    DRROT => class_rot%DRROT
    EFROT => class_rot%EFROT
    FLGGROT => class_rot%FLGGROT
    FLGSROT => class_rot%FLGSROT
    FLGVROT => class_rot%FLGVROT
    FSGGROT => class_rot%FSGGROT
    FSGSROT => class_rot%FSGSROT
    FSGVROT => class_rot%FSGVROT
    FSNOROT => class_rot%FSNOROT
    GAROT => class_rot%GAROT
    GTROT => class_rot%GTROT
    HBLROT => class_rot%HBLROT
    HEVCROT => class_rot%HEVCROT
    HEVGROT => class_rot%HEVGROT
    HEVSROT => class_rot%HEVSROT
    HFSROT => class_rot%HFSROT
    HFSCROT => class_rot%HFSCROT
    HFSGROT => class_rot%HFSGROT
    HFSSROT => class_rot%HFSSROT
    HMFCROT => class_rot%HMFCROT
    HMFNROT => class_rot%HMFNROT
    HTCCROT => class_rot%HTCCROT
    SDEPROT => class_rot%SDEPROT
    SOCIROT => class_rot%SOCIROT
    HTCSROT => class_rot%HTCSROT
    ILMOROT => class_rot%ILMOROT
    PCFCROT => class_rot%PCFCROT
    PCLCROT => class_rot%PCLCROT
    PCPGROT => class_rot%PCPGROT
    PCPNROT => class_rot%PCPNROT
    PETROT => class_rot%PETROT
    QEVPROT => class_rot%QEVPROT
    QFCFROT => class_rot%QFCFROT
    QFCLROT => class_rot%QFCLROT
    QFGROT => class_rot%QFGROT
    QFNROT => class_rot%QFNROT
    QFSROT => class_rot%QFSROT
    QFXROT => class_rot%QFXROT
    QGROT => class_rot%QGROT
    ROFROT => class_rot%ROFROT
    ROFBROT => class_rot%ROFBROT
    ROFCROT => class_rot%ROFCROT
    ROFNROT => class_rot%ROFNROT
    ROFOROT => class_rot%ROFOROT
    ROFSROT => class_rot%ROFSROT
    ROVGROT => class_rot%ROVGROT
    SFCQROT => class_rot%SFCQROT
    SFCTROT => class_rot%SFCTROT
    SFCUROT => class_rot%SFCUROT
    SFCVROT => class_rot%SFCVROT
    TFXROT => class_rot%TFXROT
    TROBROT => class_rot%TROBROT
    TROFROT => class_rot%TROFROT
    TROOROT => class_rot%TROOROT
    TROSROT => class_rot%TROSROT
    UEROT => class_rot%UEROT
    WTABROT => class_rot%WTABROT
    WTRCROT => class_rot%WTRCROT
    WTRGROT => class_rot%WTRGROT
    WTRSROT => class_rot%WTRSROT
    SFRHROT => class_rot%SFRHROT
    wtableROT => class_rot%wtableROT
    ISNDROT => class_rot%ISNDROT
    TBARROT => class_rot%TBARROT
    THICROT => class_rot%THICROT
    THLQROT => class_rot%THLQROT
    BIROT => class_rot%BIROT
    DLZWROT => class_rot%DLZWROT
    GRKSROT => class_rot%GRKSROT
    HCPSROT => class_rot%HCPSROT
    SANDROT => class_rot%SANDROT
    CLAYROT => class_rot%CLAYROT
    ORGMROT => class_rot%ORGMROT
    PSISROT => class_rot%PSISROT
    PSIWROT => class_rot%PSIWROT
    TCSROT => class_rot%TCSROT
    THFCROT => class_rot%THFCROT
    THMROT => class_rot%THMROT
    THPROT => class_rot%THPROT
    THRROT => class_rot%THRROT
    THRAROT => class_rot%THRAROT
    ZBTWROT => class_rot%ZBTWROT
    THLWROT => class_rot%THLWROT
    GFLXROT => class_rot%GFLXROT
    HMFGROT => class_rot%HMFGROT
    HTCROT => class_rot%HTCROT
    QFCROT => class_rot%QFCROT
    ACIDROT => class_rot%ACIDROT
    ACVDROT => class_rot%ACVDROT
    CMASROT => class_rot%CMASROT
    HGTDROT => class_rot%HGTDROT
    PAIDROT => class_rot%PAIDROT
    PAMNROT => class_rot%PAMNROT
    PAMXROT => class_rot%PAMXROT
    PSGAROT => class_rot%PSGAROT
    PSGBROT => class_rot%PSGBROT
    QA50ROT => class_rot%QA50ROT
    ROOTROT => class_rot%ROOTROT
    RSMNROT => class_rot%RSMNROT
    VPDAROT => class_rot%VPDAROT
    VPDBROT => class_rot%VPDBROT
    ALICROT => class_rot%ALICROT
    ALVCROT => class_rot%ALVCROT
    FCANROT => class_rot%FCANROT
    LNZ0ROT => class_rot%LNZ0ROT
    SALBROT => class_rot%SALBROT
    CSALROT => class_rot%CSALROT
    FSDBROL => class_rot%FSDBROL
    FSFBROL => class_rot%FSFBROL
    FSSBROL => class_rot%FSSBROL
    ITCTROT => class_rot%ITCTROT
    TSFSROT => class_rot%TSFSROT

    ! Point CTEM pointers

    ctem_on           => c_switch%ctem_on
    dofire            => c_switch%dofire
    PFTCompetition    => c_switch%PFTCompetition
    start_bare        => c_switch%start_bare
    lnduseon          => c_switch%lnduseon
    transientCO2      => c_switch%transientCO2
    doMethane         => c_switch%doMethane
    transientCH4      => c_switch%transientCH4
    transientPOPD     => c_switch%transientPOPD
    transientLGHT     => c_switch%transientLGHT
    fixedYearLUC      => c_switch%fixedYearLUC
    transientOBSWETF  => c_switch%transientOBSWETF
    fixedYearOBSWETF  => c_switch%fixedYearOBSWETF
    inibioclim        => c_switch%inibioclim
    leap              => c_switch%leap
    jhhstd            => c_switch%jhhstd
    jhhendd           => c_switch%jhhendd
    jdstd             => c_switch%jdstd
    jdendd            => c_switch%jdendd
    jhhsty            => c_switch%jhhsty
    jhhendy           => c_switch%jhhendy
    jdsty             => c_switch%jdsty
    jdendy            => c_switch%jdendy
    jmosty            => c_switch%jmosty
    metLoop           => c_switch%metLoop
    spinfast          => c_switch%spinfast
    useTracer         => c_switch%useTracer
    IDISP             => c_switch%IDISP
    IZREF             => c_switch%IZREF
    ISLFD             => c_switch%ISLFD
    IPCP              => c_switch%IPCP
    ITC               => c_switch%ITC
    ITCG              => c_switch%ITCG
    ITG               => c_switch%ITG
    IWF               => c_switch%IWF
    IPAI              => c_switch%IPAI
    IHGT              => c_switch%IHGT
    IALC              => c_switch%IALC
    IALS              => c_switch%IALS
    IALG              => c_switch%IALG
    isnoalb           => c_switch%isnoalb
    doAnnualOutput    => c_switch%doAnnualOutput
    doMonthOutput     => c_switch%doMonthOutput
    doDayOutput       => c_switch%doDayOutput
    doHhOutput        => c_switch%doHhOutput
    readMetStartYear  => c_switch%readMetStartYear
    projectedGrid     => c_switch%projectedGrid

    tcanrs            => vrot%tcanrs
    tsnors            => vrot%tsnors
    tpndrs            => vrot%tpndrs
    csum              => vrot%csum
    tbaraccrow_m      => class_rot%tbaraccrow_m

    ! ROW:
    gleafmasrow       => vrot%gleafmas
    bleafmasrow       => vrot%bleafmas
    stemmassrow       => vrot%stemmass
    rootmassrow       => vrot%rootmass
    pstemmassrow      => vrot%pstemmass
    pgleafmassrow     => vrot%pgleafmass
    fcancmxrow        => vrot%fcancmx
    gavglairow        => vrot%gavglai
    zolncrow          => vrot%zolnc
    ailcrow           => vrot%ailc
    ailcgrow          => vrot%ailcg
    ailcgsrow         => vrot%ailcgs
    fcancsrow         => vrot%fcancs
    fcancrow          => vrot%fcanc
    co2concrow        => vrot%co2conc
    ch4concrow        => vrot%ch4conc
    co2i1cgrow        => vrot%co2i1cg
    co2i1csrow        => vrot%co2i1cs
    co2i2cgrow        => vrot%co2i2cg
    co2i2csrow        => vrot%co2i2cs
    ancsvegrow        => vrot%ancsveg
    ancgvegrow        => vrot%ancgveg
    rmlcsvegrow       => vrot%rmlcsveg
    rmlcgvegrow       => vrot%rmlcgveg
    slairow           => vrot%slai
    ailcbrow          => vrot%ailcb
    canresrow         => vrot%canres
    flhrlossrow       => vrot%flhrloss

    qevpacc_m_save    => vrot%qevpacc_m_save

    grwtheffrow       => vrot%grwtheff
    lystmmasrow       => vrot%lystmmas
    lyrotmasrow       => vrot%lyrotmas
    tymaxlairow       => vrot%tymaxlai
    vgbiomasrow       => vrot%vgbiomas
    gavgltmsrow       => vrot%gavgltms
    gavgscmsrow       => vrot%gavgscms
    stmhrlosrow       => vrot%stmhrlos
    rmatcrow          => vrot%rmatc
    rmatctemrow       => vrot%rmatctem
    litrmassrow       => vrot%litrmass
    soilcmasrow       => vrot%soilcmas
    vgbiomas_vegrow   => vrot%vgbiomas_veg

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
    popdinrow         => vrot%popdin
    btermrow          => vrot%bterm
    ltermrow          => vrot%lterm
    mtermrow          => vrot%mterm

    extnprobrow       => vrot%extnprob
    prbfrhucrow       => vrot%prbfrhuc
    daylrow           => vrot%dayl
    dayl_maxrow       => vrot%dayl_max

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

    slopefracrow      => vrot%slopefrac
    wetfrac_presrow   => vrot%wetfrac_pres
    ch4WetSpecrow     => vrot%ch4WetSpec
    wetfdynrow        => vrot%wetfdyn
    ch4WetDynrow      => vrot%ch4WetDyn
    ch4soillsrow      => vrot%ch4_soills

    peatdeprow        => vrot%peatdep
    lucemcomrow       => vrot%lucemcom
    lucltrinrow       => vrot%lucltrin
    lucsocinrow       => vrot%lucsocin

    npprow            => vrot%npp
    neprow            => vrot%nep
    nbprow            => vrot%nbp
    gpprow            => vrot%gpp
    hetroresrow       => vrot%hetrores
    autoresrow        => vrot%autores
    soilcresprow      => vrot%soilcresp
    rmrow             => vrot%rm
    rgrow             => vrot%rg
    litresrow         => vrot%litres
    socresrow         => vrot%socres
    dstcemlsrow       => vrot%dstcemls
    litrfallrow       => vrot%litrfall
    humiftrsrow       => vrot%humiftrs

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
    litrfallvegrow    => vrot%litrfallveg
    humiftrsvegrow    => vrot%humiftrsveg

    rothrlosrow       => vrot%rothrlos
    pfcancmxrow       => vrot%pfcancmx
    nfcancmxrow       => vrot%nfcancmx
    alvsctmrow        => vrot%alvsctm
    paicrow           => vrot%paic
    slaicrow          => vrot%slaic
    alirctmrow        => vrot%alirctm
    cfluxcgrow        => vrot%cfluxcg
    cfluxcsrow        => vrot%cfluxcs
    dstcemls3row      => vrot%dstcemls3
    anvegrow          => vrot%anveg
    rmlvegrow         => vrot%rmlveg

    pftexistrow       => vrot%pftexist
    colddaysrow       => vrot%colddays
    lfstatusrow       => vrot%lfstatus
    pandaysrow        => vrot%pandays

    twarmmrow            => vrot%twarmm
    tcoldmrow            => vrot%tcoldm
    gdd5row              => vrot%gdd5
    aridityrow           => vrot%aridity
    srplsmonrow          => vrot%srplsmon
    defctmonrow          => vrot%defctmon
    anndefctrow          => vrot%anndefct
    annsrplsrow          => vrot%annsrpls
    annpcprow            => vrot%annpcp
    dry_season_lengthrow => vrot%dry_season_length

    ipeatlandrow     => vrot%ipeatland
    anmossrow        => vrot%anmoss
    rmlmossrow       => vrot%rmlmoss
    gppmossrow       => vrot%gppmoss
    nppmossrow       => vrot%nppmoss
    armossrow        => vrot%armoss
    litrmsmossrow    => vrot%litrmsmoss
    Cmossmasrow      => vrot%Cmossmas
    dmossrow         => vrot%dmoss
    pddrow           => vrot%pdd

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>
    ! GAT:

    grclarea          => vgat%grclarea
    lightng           => vgat%lightng

    gleafmasgat       => vgat%gleafmas
    bleafmasgat       => vgat%bleafmas
    stemmassgat       => vgat%stemmass
    rootmassgat       => vgat%rootmass
    pstemmassgat      => vgat%pstemmass
    pgleafmassgat     => vgat%pgleafmass
    fcancmxgat        => vgat%fcancmx
    gavglaigat        => vgat%gavglai
    zolncgat          => vgat%zolnc
    ailcgat           => vgat%ailc
    ailcggat          => vgat%ailcg
    ailcgsgat         => vgat%ailcgs
    fcancsgat         => vgat%fcancs
    fcancgat          => vgat%fcanc
    co2concgat        => vgat%co2conc
    ch4concgat        => vgat%ch4conc
    co2i1cggat        => vgat%co2i1cg
    co2i1csgat        => vgat%co2i1cs
    co2i2cggat        => vgat%co2i2cg
    co2i2csgat        => vgat%co2i2cs
    ancsveggat        => vgat%ancsveg
    ancgveggat        => vgat%ancgveg
    rmlcsveggat       => vgat%rmlcsveg
    rmlcgveggat       => vgat%rmlcgveg
    slaigat           => vgat%slai
    ailcbgat          => vgat%ailcb
    canresgat         => vgat%canres
    flhrlossgat       => vgat%flhrloss

    grwtheffgat       => vgat%grwtheff
    lystmmasgat       => vgat%lystmmas
    lyrotmasgat       => vgat%lyrotmas
    tymaxlaigat       => vgat%tymaxlai
    vgbiomasgat       => vgat%vgbiomas
    gavgltmsgat       => vgat%gavgltms
    gavgscmsgat       => vgat%gavgscms
    stmhrlosgat       => vgat%stmhrlos
    rmatcgat          => vgat%rmatc
    rmatctemgat       => vgat%rmatctem
    litrmassgat       => vgat%litrmass
    soilcmasgat       => vgat%soilcmas
    vgbiomas_veggat   => vgat%vgbiomas_veg
    litrfallveggat    => vgat%litrfallveg
    humiftrsveggat    => vgat%humiftrsveg
    reprocost         => vgat%reprocost

    altotcount_ctm    => vgat%altotcount_ctm
    todfrac           => vgat%todfrac
    fsinacc_gat       => vgat%fsinacc_gat
    flutacc_gat       => vgat%flutacc_gat
    flinacc_gat       => vgat%flinacc_gat
    altotacc_gat      => vgat%altotacc_gat
    netrad_gat        => vgat%netrad_gat
    preacc_gat        => vgat%preacc_gat
    sdepgat           => vgat%sdepgat
    sandgat           => vgat%sandgat
    claygat           => vgat%claygat
    orgmgat           => vgat%orgmgat
    xdiffusgat        => vgat%xdiffusgat
    faregat           => vgat%faregat

    emit_co2gat       => vgat%emit_co2
    emit_cogat        => vgat%emit_co
    emit_ch4gat       => vgat%emit_ch4
    emit_nmhcgat      => vgat%emit_nmhc
    emit_h2gat        => vgat%emit_h2
    emit_noxgat       => vgat%emit_nox
    emit_n2ogat       => vgat%emit_n2o
    emit_pm25gat      => vgat%emit_pm25
    emit_tpmgat       => vgat%emit_tpm
    emit_tcgat        => vgat%emit_tc
    emit_ocgat        => vgat%emit_oc
    emit_bcgat        => vgat%emit_bc
    burnfracgat       => vgat%burnfrac
    burnvegfgat       => vgat%burnvegf
    popdingat         => vgat%popdin
    smfuncveggat      => vgat%smfuncveg
    btermgat          => vgat%bterm
    ltermgat          => vgat%lterm
    mtermgat          => vgat%mterm
    blfltrdt          => vgat%blfltrdt
    glfltrdt          => vgat%glfltrdt
    glcaemls          => vgat%glcaemls
    blcaemls          => vgat%blcaemls
    rtcaemls          => vgat%rtcaemls
    stcaemls          => vgat%stcaemls
    ltrcemls          => vgat%ltrcemls
    ntchlveg          => vgat%ntchlveg
    ntchsveg          => vgat%ntchsveg
    ntchrveg          => vgat%ntchrveg
    extnprobgat       => vgat%extnprob
    prbfrhucgat       => vgat%prbfrhuc
    daylgat           => vgat%dayl
    dayl_maxgat       => vgat%dayl_max

    bmasveggat        => vgat%bmasveg
    cmasvegcgat       => vgat%cmasvegc
    veghghtgat        => vgat%veghght
    rootdpthgat       => vgat%rootdpth
    rmlgat            => vgat%rml
    rmsgat            => vgat%rms
    tltrleafgat       => vgat%tltrleaf
    tltrstemgat       => vgat%tltrstem
    tltrrootgat       => vgat%tltrroot
    leaflitrgat       => vgat%leaflitr
    roottempgat       => vgat%roottemp
    afrleafgat        => vgat%afrleaf
    afrstemgat        => vgat%afrstem
    afrrootgat        => vgat%afrroot
    wtstatusgat       => vgat%wtstatus
    ltstatusgat       => vgat%ltstatus
    rmrgat            => vgat%rmr

    slopefracgat      => vgat%slopefrac
    ch4WetSpecgat     => vgat%ch4WetSpec
    wetfdyngat        => vgat%wetfdyn
    wetfrac_presgat   => vgat%wetfrac_pres
    ch4WetDyngat      => vgat%ch4WetDyn
    ch4soillsgat      => vgat%ch4_soills

    lucemcomgat       => vgat%lucemcom
    lucltringat       => vgat%lucltrin
    lucsocingat       => vgat%lucsocin

    nppgat            => vgat%npp
    nepgat            => vgat%nep
    nbpgat            => vgat%nbp
    gppgat            => vgat%gpp
    hetroresgat       => vgat%hetrores
    autoresgat        => vgat%autores
    soilcrespgat      => vgat%soilcresp
    rmgat             => vgat%rm
    rggat             => vgat%rg
    litresgat         => vgat%litres
    socresgat         => vgat%socres
    dstcemlsgat       => vgat%dstcemls
    litrfallgat       => vgat%litrfall
    humiftrsgat       => vgat%humiftrs

    gppveggat         => vgat%gppveg
    nepveggat         => vgat%nepveg
    nbpveggat         => vgat%nbpveg
    nppveggat         => vgat%nppveg
    hetroresveggat    => vgat%hetroresveg
    autoresveggat     => vgat%autoresveg
    litresveggat      => vgat%litresveg
    soilcresveggat    => vgat%soilcresveg
    rmlvegaccgat      => vgat%rmlvegacc
    rmsveggat         => vgat%rmsveg
    rmrveggat         => vgat%rmrveg
    rgveggat          => vgat%rgveg

    rothrlosgat       => vgat%rothrlos
    pfcancmxgat       => vgat%pfcancmx
    nfcancmxgat       => vgat%nfcancmx
    alvsctmgat        => vgat%alvsctm
    paicgat           => vgat%paic
    slaicgat          => vgat%slaic
    alirctmgat        => vgat%alirctm
    cfluxcggat        => vgat%cfluxcg
    cfluxcsgat        => vgat%cfluxcs
    dstcemls3gat      => vgat%dstcemls3
    anveggat          => vgat%anveg
    rmlveggat         => vgat%rmlveg

    twarmmgat            => vgat%twarmm
    tcoldmgat            => vgat%tcoldm
    gdd5gat              => vgat%gdd5
    ariditygat           => vgat%aridity
    srplsmongat          => vgat%srplsmon
    defctmongat          => vgat%defctmon
    anndefctgat          => vgat%anndefct
    annsrplsgat          => vgat%annsrpls
    annpcpgat            => vgat%annpcp
    dry_season_lengthgat => vgat%dry_season_length

    tcurm             => vgat%tcurm
    srpcuryr          => vgat%srpcuryr
    dftcuryr          => vgat%dftcuryr
    tmonth            => vgat%tmonth
    anpcpcur          => vgat%anpcpcur
    anpecur           => vgat%anpecur
    gdd5cur           => vgat%gdd5cur
    surmncur          => vgat%surmncur
    defmncur          => vgat%defmncur
    srplscur          => vgat%srplscur
    defctcur          => vgat%defctcur

    geremortgat       => vgat%geremort
    intrmortgat       => vgat%intrmort
    lambdagat         => vgat%lambda
    ccgat             => vgat%cc
    mmgat             => vgat%mm

    pftexistgat       => vgat%pftexist
    colddaysgat       => vgat%colddays
    lfstatusgat       => vgat%lfstatus
    pandaysgat        => vgat%pandays

    ipeatlandgat     => vgat%ipeatland
    peatdepgat       => vgat%peatdep
    anmossgat        => vgat%anmoss
    rmlmossgat       => vgat%rmlmoss
    gppmossgat       => vgat%gppmoss
    nppmossgat       => vgat%nppmoss
    armossgat        => vgat%armoss
    litrmsmossgat    => vgat%litrmsmoss
    Cmossmasgat      => vgat%Cmossmas
    dmossgat         => vgat%dmoss
    pddgat           => vgat%pdd
    ancsmoss         => vgat%ancsmoss
    angsmoss         => vgat%angsmoss
    ancmoss          => vgat%ancmoss
    angmoss          => vgat%angmoss
    rmlcsmoss        => vgat%rmlcsmoss
    rmlgsmoss        => vgat%rmlgsmoss
    rmlcmoss         => vgat%rmlcmoss
    rmlgmoss         => vgat%rmlgmoss

    ! mosaic level variables (CLASS):

    fsnowacc_t        => ctem_tile%fsnowacc_t
    taaccgat_t        => ctem_tile%taaccgat_t
    uvaccgat_t        => ctem_tile%uvaccgat_t
    vvaccgat_t        => ctem_tile%vvaccgat_t
    tbaraccgat_t      => ctem_tile%tbaraccgat_t
    thliqacc_t        => ctem_tile%thliqacc_t
    thiceacc_t        => ctem_tile%thiceacc_t  ! Added in place of YW's thicaccgat_m. EC Dec 23 2016.
    ancgvgac_t        => ctem_tile%ancgvgac_t
    rmlcgvga_t        => ctem_tile%rmlcgvga_t
    anmossac_t        => ctem_tile%anmossac_t
    rmlmossac_t       => ctem_tile%rmlmossac_t
    gppmossac_t       => ctem_tile%gppmossac_t

    tracerGLeafMassrot   => tracer%gLeafMassrot
    tracerBLeafMassrot   => tracer%bLeafMassrot
    tracerStemMassrot    => tracer%stemMassrot
    tracerRootMassrot    => tracer%rootMassrot
    tracerLitrMassrot    => tracer%litrMassrot
    tracerSoilCMassrot   => tracer%soilCMassrot
    tracerMossCMassrot   => tracer%mossCMassrot
    tracerMossLitrMassrot => tracer%mossLitrMassrot
    tracerCO2rot         => tracer%tracerCO2rot

    tracerGLeafMassgat   => tracer%gLeafMassgat
    tracerBLeafMassgat   => tracer%bLeafMassgat
    tracerStemMassgat    => tracer%stemMassgat
    tracerRootMassgat    => tracer%rootMassgat
    tracerLitrMassgat    => tracer%litrMassgat
    tracerSoilCMassgat   => tracer%soilCMassgat
    tracerMossCMassgat   => tracer%mossCMassgat
    tracerMossLitrMassgat => tracer%mossLitrMassgat
    tracerCO2gat          => tracer%tracerCO2gat

    !    =================================================================================

    !> NLTEST and NMTEST are the number of grid cells and the number of mosaic tiles per grid cell for this test run, respectively.
    !! This driver is set up to handle one grid cell with any number of mosaic tiles. These are given the values then
    !! of nlat and nmos.
    nltest = nlat
    nmtest = nmos
    NTLD = NMOS

    !> The parameter JLAT is calculated from DLATROW as the nearest integer :: value,
    DLATROW = latitude
    JLAT = NINT(DLATROW(1))
    DLONROW = longitude

    !> The timestep counter N for the run is initialized to 0, the daily
    !! averaging counter NCOUNT is set to 1, and the total number of
    !! timesteps in the day NDAY is calculated as the number of seconds
    !! in a day (86400) divided by the timestep length DELT.
    N = 0
    NCOUNT = 1
    NDAY = 86400/NINT(DELT)
    metTimeIndex = 1    !< Counter used to move through the meteorological input arrays
    metDone = .false.   !< Logical switch when the end of the stored meteorological array is reached.
    run_model = .true.  !< Simple logical switch to either keep run going or finish
    IMONTH = 0          !< Month of the year simulation is in.
    DOM = 1             !< Day of month counter
    CUMSNO = 0.0
    lopcount = 1
    leapnow = .false.
    lastDOY = 365
    wetfrac_presgat = - 9999. !< If transientOBSWETF or fixedYearOBSWETF != -9999 this variable will be overwritten with
    !! real wetland fractions. Otherwise the negative is used as a switch so the dynamic
    !! wetland extent is used instead of the prescribed.

    !> The grid-average height for the momentum diagnostic variables, ZDMROW, and for the
    !! energy diagnostic variables, ZDHROW, are hard-coded to the standard anemometer
    !! height of 10 m and to the screen height of 2 m respectively.
    ZDMROW(:) = 10.0
    ZDHROW(:) = 2.0

    !> ZRFMROW and ZRFHROW, the reference heights at which the momentum variables (wind speed) and energy variables
    !! (temperature and specific humidity) are provided.  In a run using atmospheric model forcing data, these heights
    !! would vary by time step, but since this version of the driver is set up to use field data, ZRFMROW and ZRFHROW
    !! refer to the measurement height of these variables, which is fixed. The value is read in from the job options file.
    ZRFMROW = zrfmJobOpt
    ZRFHROW = zrfhJobOpt

    !> ZBLDROW, the atmospheric blending height.  Technically this variable depends on the length scale of the
    !! patches of roughness elements on the land surface, but this is difficult to ascertain.  Usually it is assigned a value of 50 m.
    !!  The value is read in from the job options file.
    ZBLDROW = zbldJobOpt

    !> Initialize variables in preparation for the run
    call initRowVars
    call resetAccVars(nlat, nmos)
    call resetMosaicAccum

    !> Read in the model initial state
    call read_initialstate(lonIndex, latIndex)
    
    if (ctem_on) then
      !> Read in the inputs for a run with biogeochemical component turned on
      call getInput('CO2') ! CO2 atmospheric concentration
      if (doMethane) call getInput('CH4') ! CH4 atmospheric concentration
      if (useTracer > 0) call getInput('tracerCO2',longitude,latitude) ! tracer atmospheric values
      if (.not. projectedGrid) then
        ! regular lon/lat grid
        if (dofire) call getInput('POPD',longitude,latitude) ! Population density
        if (dofire) call getInput('LGHT',longitude,latitude) ! Cloud-to-ground lightning frequency
        if (doMethane .and. transientOBSWETF .or. fixedYearOBSWETF /= - 9999) call getInput('OBSWETF',longitude,latitude) ! Observed wetland distribution
        if (lnduseon .or. (fixedYearLUC /= - 9999)) call getInput('LUC',longitude,latitude) ! Land use change
      else
        ! Projected grids use the lon and lat indexes, not the actual coordinates
        if (dofire) call getInput('POPD',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Population density
        if (dofire) call getInput('LGHT',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Cloud-to-ground lightning frequency
        if (doMethane .and. transientOBSWETF .or. fixedYearOBSWETF /= - 9999) &
            call getInput('OBSWETF',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Observed wetland distribution
        if (lnduseon .or. (fixedYearLUC /= - 9999)) &
            call getInput('LUC',longitude,latitude,projLonInd = lonIndex,projLatInd = latIndex) ! Land use change
      end if
      !> Regardless of whether lnduseon or not, we need to check the land cover that was read in
      !! and assign the CLASS PFTs as they are not read in when ctem_on.
      call initializeLandCover
    end if

    !> Read in the meteorological forcing data to a suite of arrays
    if (.not. projectedGrid) then
      ! regular lon lat grid
      call getMet(longitude, latitude, nday)
    else
      ! Projected grids use the lon and lat indexes, not the actual coordinates
      call getMet(longitude, latitude, nday, projLonInd = lonIndex, projLatInd = latIndex)
    end if

    !> In preparation for the use of the random number generator by disaggMet,
    !! we need to provide a seed to allow repeatable results.
    call initRandomSeed()

    !> Now disaggregate the meteorological forcing to the right timestep
    !! for this model run (if needed; this is checked for in the subroutine)
    call disaggMet(longitude, latitude)

    ! Initialize accumulated array for monthly & yearly outputs
    call resetClassMon(nltest)
    call resetClassYr(nltest)
    if (ctem_on) then
      call resetMonthEnd(nltest, nmtest)
      call resetYearEnd(nltest, nmtest)
    end if

    !> As the last step in the initialization sequence, the subroutine soilProperties is
    !! called, to assign soil thermal and hydraulic properties on the basis of the
    !! textural information read in for each of the soil layers.

    call soilProperties(THPROT, THRROT, THMROT, BIROT, PSISROT, GRKSROT, & ! Formerly CLASSB
                        THRAROT, HCPSROT, TCSROT, THFCROT, THLWROT, PSIWROT, &
                        DLZWROT, ZBTWROT, &
                        ALGWVROT, ALGWNROT, ALGDVROT, ALGDNROT, &
                        SANDROT, CLAYROT, ORGMROT, SOCIROT, DELZ, ZBOT, &
                        SDEPROT, ISNDROT, IGDRROT, &
                        NLAT, NMOS, 1, NLTEST, NMTEST, IGND, ipeatlandrow)

    ! ctem initializations.
    if (ctem_on) then
      call ctemInit(nltest, nmtest)

      call classGatherPrep(ILMOS, JLMOS, IWMOS, JWMOS, & ! Formerly GATPREP
                           NML, NMW, GCROW, FAREROT, MIDROT, &
                           NLAT, NMOS, ILG, 1, NLTEST, NMTEST)

      ! ctemg1 converts variables from the 'row' format (nlat, nmos, ...)
      ! to the 'gat' format (ilg, ...) which is what the model calculations
      ! are performed on. The ctemg1 subroutine is used to transform the
      ! read in state variables (which come in with the 'row' format from the
      ! various input files).

      call ctemg1(gleafmasgat, bleafmasgat, stemmassgat, &  ! Out
                         rootmassgat, fcancmxgat, zbtwgat, & ! Out
                         dlzwgat, sdepgat, ailcggat, & ! Out
                         ailcbgat, ailcgat, zolncgat, rmatcgat, & ! Out
                         rmatctemgat, slaigat, bmasveggat, cmasvegcgat, & ! Out
                         veghghtgat, rootdpthgat, alvsctmgat, alirctmgat, & ! Out
                         paicgat, slaicgat, faregat, & ! Out
                         ipeatlandgat, maxAnnualActLyrGAT, & ! Out
                         tracergLeafMassgat, tracerBLeafMassgat, tracerStemMassgat, & ! Out
                         tracerRootMassgat, tracerLitrMassgat, tracerSoilCMassgat, & ! Out
                         tracerMossCMassgat, tracerMossLitrMassgat, & ! Out                     
                         twarmmgat, tcoldmgat, gdd5gat, & ! Out
                         ariditygat, srplsmongat, defctmongat, anndefctgat, & ! Out
                         annsrplsgat, annpcpgat, dry_season_lengthgat, & ! Out
                         litrmsmossgat, Cmossmasgat, dmossgat, & ! Out                     
                         pandaysgat, lfstatusgat, slopefracgat, pstemmassgat, & ! Out
                         pgleafmassgat,litrmassgat, soilcmasgat, grwtheffgat, & ! Out
                         ilmos, jlmos, iwmos, jwmos, nml, &! In
                         gleafmasrow, bleafmasrow, stemmassrow, rootmassrow, &! In
                         fcancmxrow, zbtwrot, dlzwrot, sdeprot, &! In
                         ailcgrow, ailcbrow, ailcrow, zolncrow, &! In
                         rmatcrow, rmatctemrow, slairow, bmasvegrow, &! In
                         cmasvegcrow, veghghtrow, rootdpthrow, alvsctmrow, &! In
                         alirctmrow, paicrow, slaicrow, FAREROT, &! In
                         ipeatlandrow, maxAnnualActLyrROT, &! In
                         tracergLeafMassrot, tracerBLeafMassrot, tracerStemMassrot, &! In
                         tracerRootMassrot, tracerLitrMassrot, tracerSoilCMassrot, &! In
                         tracerMossCMassrot, tracerMossLitrMassrot, & ! In
                         twarmmrow, tcoldmrow, gdd5row, & ! In
                         aridityrow, srplsmonrow, defctmonrow, anndefctrow, & ! In
                         annsrplsrow, annpcprow, dry_season_lengthrow, & ! In
                         litrmsmossrow, Cmossmasrow, dmossrow, & ! In
                         pandaysrow, lfstatusrow, slopefracrow, pstemmassrow, & ! In
                         pgleafmassrow, litrmassrow, soilcmasrow, grwtheffrow)! In

      !> Find mosaic tile (grid) average vegetation biomass, litter mass, and soil c mass.
      !! Set growth efficiency to some large number so that no growth related mortality
      !! occurs in first year. For peatlands determine the peatdepth and the peat soil 
      !! carbon amounts. Lastly find the maximum daylength for this location.
      
      !call ctemInit

      call allometry(gleafmasgat, bleafmasgat, stemmassgat, rootmassgat, & ! In
                     1, nml, ilg, zbtwgat, & ! In
                     sdepgat, fcancmxgat, & ! In
                     ipeatlandgat, maxAnnualActLyrGAT, & ! In
                     ailcggat, ailcbgat, ailcgat, zolncgat, & ! Out
                     rmatcgat, rmatctemgat, slaigat, bmasveggat, & ! Out
                     cmasvegcgat, veghghtgat, rootdpthgat, alvsctmgat, & ! Out
                     alirctmgat, paicgat, slaicgat) ! Out

    end if   ! if (ctem_on)

    !     ctem initial preparation done

    !     **** LAUNCH RUN. ****

    !> The do while loop marks the beginning of the time stepping loop
    !! for the actual run.  N is incremented by 1, and the atmospheric forcing
    !! data for the current time step are updated for each grid cell or modelled
    !! area (see the manual section on “Data Requirements”).

    runyr = readMetStartYear  ! Initialize the runyr as the first year of met forcing.

    ! start up the main model loop

    mainModelLoop: do while (run_model)

      !
      !> Update the meteorological forcing data for current time step
      !
      call updateMet(metTimeIndex, iyear, iday, ihour, imin, metDone)

      ! print*, 'year=', iyear, 'day=', iday, ' hour=', ihour, ' min=', imin

      N = N + 1

      !> Generally only the total incoming shortwave radiation FSDOWN
      !! is available; so it is partitioned 50:50 between the incoming visible (FSVHROW)
      !! and near-infrared (FSIHROW) radiation.  The first two elements of the
      !! generalized incoming radiation array, FSSBROL (used for both the ISNOALB=0
      !! and ISNOALB=1 options) are set to FSVHROW and FSIHROW respectively.
      !! The air temperature TAROW is converted from degrees C to K.  The zonal
      !! (ULROW) and meridional (VLROW) components of the wind speed are generally not
      !! used; only the overall wind speed UVROW is
      !! measured.  However, CLASS does not require wind direction for its calculations,
      !! so ULROW is arbitrarily assigned the value of UVROW and VLROW is set to zero for
      !! this run.  The input wind speed VMODROW is assigned the value of UVROW.

      do I = 1,NLTEST

        FSVHROW(I) = 0.5 * FSSROW(I)
        FSIHROW(I) = 0.5 * FSSROW(I)
        TAROW(I) = TAROW(I) + TFREZ
        ULROW(I) = UVROW(I)
        VLROW(I) = 0.0
        VMODROW(I) = UVROW(I)
        FSSBROL(I,1) = FSVHROW(I)
        FSSBROL(I,2) = FSIHROW(I)

      end do ! loop 250

      ! Check if we are on the first timestep of the day
      if (ihour == 0 .and. imin == 0) then

        ! Find the daylength of this day
        daylrow = findDaylength(real(iday),radjrow(1)) ! following rest of code, radjrow is always given index of 1 offline.
        
        ! Update the lightning if fire is on and transientLGHT is true
        if (dofire .and. ctem_on) call updateInput('LGHT',runyr,imonth = imonth,iday = iday,dom = DOM)

        ! Update the wetland fractions if we are using read-in wetland fractions
        if (ctem_on .and. doMethane .and. (transientOBSWETF .or. fixedYearOBSWETF /= - 9999)) then
          call updateInput('OBSWETF',runyr,imonth = imonth,iday = iday,dom = DOM)
        end if

        ! Check if this is the first day of the year
        if (iday == 1) then

          ! Check if this year is a leap year, and if so adjust the monthdays, monthend and mmday values.
          if (leap) call findLeapYears(iyear,leapnow,lastDOY)

          ! If needed, update values that were read in from the accessory input files (popd, wetlands, lightning...)
          if (ctem_on) then

            if (transientCO2) call updateInput('CO2',runyr)
            if (useTracer > 0 .and. transientCO2) call updateInput('tracerCO2',runyr)
            if (doMethane .and. transientCH4) call updateInput('CH4',runyr)
            if (dofire .and. transientPOPD) call updateInput('POPD',runyr)
            if (lnduseon) then
              call updateInput('LUC',runyr)
            else ! If landuse change is not on, then set the next years landcover to be
              ! the same as this years.
              nfcancmxrow = pfcancmxrow
            end if

            ! Reset the peatland degree days counter (JM:FLAG is this ok for S. Hemi to be Jan 1?)
            pddrow = 0

          end if
        end if ! first day
      end if   ! first timestep

      !> The cosine of the solar zenith angle COSZ is calculated from the day of
      !> the year, the hour, the minute and the latitude using basic radiation geometry,
      !> and (avoiding vanishingly small numbers) is assigned to CSZROW.  The fractional
      !> cloud cover FCLOROW is commonly not available so a rough estimate is
      !> obtained by setting it to 1 when precipitation is occurring, and to the fraction
      !> of incoming diffuse radiation XDIFFUS otherwise (assumed to be 1 when the sun
      !> is at the horizon, and 0.10 when it is at the zenith). These calculations are
      !> done in findCloudiness

      call findCloudiness(nltest, imin, ihour, iday, lastDOY)

      !> atmosphericVarsCalc evaluates a series of derived atmospheric variables

      call atmosphericVarsCalc(VPDROW, TADPROW, PADRROW, RHOAROW, RHSIROW, & ! Formerly CLASSI
                               RPCPROW, TRPCROW, SPCPROW, TSPCROW, TAROW, QAROW, &
                               PREROW, RPREROW, SPREROW, PRESROW, &
                               IPCP, NLAT, 1, NLTEST)

      CUMSNO = CUMSNO + SPCPROW(1) * RHSIROW(1) * DELT

      !> classGatherPrep assigns values to vectors governing the gather-scatter operations

      call classGatherPrep(ILMOS, JLMOS, IWMOS, JWMOS, & ! Formerly GATPREP
                           NML, NMW, GCROW, FAREROT, MIDROT, &
                           NLAT, NMOS, ILG, 1, NLTEST, NMTEST)

      !> classGather performs the gather operation, gathering variables from their
      !> positions as mosaic tiles within the modelled areas to long vectors of mosaic tiles
      
      call classGather(TBARGAT, THLQGAT, THICGAT, TPNDGAT, ZPNDGAT, & ! Formerly CLASSG
                       TBASGAT, ALBSGAT, TSNOGAT, RHOSGAT, SNOGAT, &
                       TCANGAT, RCANGAT, SCANGAT, GROGAT, CMAIGAT, &
                       FCANGAT, LNZ0GAT, ALVCGAT, ALICGAT, PAMXGAT, &
                       PAMNGAT, CMASGAT, ROOTGAT, RSMNGAT, QA50GAT, &
                       VPDAGAT, VPDBGAT, PSGAGAT, PSGBGAT, PAIDGAT, &
                       HGTDGAT, ACVDGAT, ACIDGAT, TSFSGAT, WSNOGAT, &
                       THPGAT, THRGAT, THMGAT, BIGAT, PSISGAT, &
                       GRKSGAT, THRAGAT, HCPSGAT, TCSGAT, IGDRGAT, &
                       THFCGAT, THLWGAT, PSIWGAT, DLZWGAT, ZBTWGAT, &
                       VMODGAT, ZSNLGAT, ZPLGGAT, ZPLSGAT, TACGAT, &
                       QACGAT, DRNGAT, XSLPGAT, GRKFGAT, WFSFGAT, &
                       WFCIGAT, ALGWVGAT, ALGWNGAT, ALGDVGAT, &
                       ALGDNGAT, ASVDGAT, ASIDGAT, AGVDGAT, &
                       AGIDGAT, ISNDGAT, RADJGAT, ZBLDGAT, Z0ORGAT, &
                       ZRFMGAT, ZRFHGAT, ZDMGAT, ZDHGAT, FSVHGAT, &
                       FSIHGAT, FSDBGAT, FSFBGAT, FSSBGAT, CSZGAT, &
                       FSGGAT, FLGGAT, FDLGAT, ULGAT, VLGAT, &
                       TAGAT, QAGAT, PRESGAT, PREGAT, PADRGAT, &
                       VPDGAT, TADPGAT, RHOAGAT, RPCPGAT, TRPCGAT, &
                       SPCPGAT, TSPCGAT, RHSIGAT, FCLOGAT, DLONGAT, &
                       GGEOGAT, GUSTGAT, REFGAT, BCSNGAT, DEPBGAT, &
                       DLATGAT, maxAnnualActLyrGAT, ILMOS, JLMOS, &
                       NML, NLAT, NTLD, NMOS, ILG, IGND, ICAN, ICAN + 1, NBS, &
                       TBARROT, THLQROT, THICROT, TPNDROT, ZPNDROT, &
                       TBASROT, ALBSROT, TSNOROT, RHOSROT, SNOROT, &
                       TCANROT, RCANROT, SCANROT, GROROT, CMAIROT, &
                       FCANROT, LNZ0ROT, ALVCROT, ALICROT, PAMXROT, &
                       PAMNROT, CMASROT, ROOTROT, RSMNROT, QA50ROT, &
                       VPDAROT, VPDBROT, PSGAROT, PSGBROT, PAIDROT, &
                       HGTDROT, ACVDROT, ACIDROT, TSFSROT, WSNOROT, &
                       THPROT, THRROT, THMROT, BIROT, PSISROT, &
                       GRKSROT, THRAROT, HCPSROT, TCSROT, IGDRROT, &
                       THFCROT, THLWROT, PSIWROT, DLZWROT, ZBTWROT, &
                       VMODROW, ZSNLROT, ZPLGROT, ZPLSROT, TACROT, &
                       QACROT, DRNROT, XSLPROT, GRKFROT, WFSFROT, &
                       WFCIROT, ALGWVROT, ALGWNROT, ALGDVROT, &
                       ALGDNROT, ASVDROT, ASIDROT, AGVDROT, &
                       AGIDROT, ISNDROT, RADJROW, ZBLDROW, Z0ORROW, &
                       ZRFMROW, ZRFHROW, ZDMROW, ZDHROW, FSVHROW, &
                       FSIHROW, FSDBROL, FSFBROL, FSSBROL, CSZROW, &
                       FSGROL, FLGROL, FDLROW, ULROW, VLROW, &
                       TAROW, QAROW, PRESROW, PREROW, PADRROW, &
                       VPDROW, TADPROW, RHOAROW, RPCPROW, TRPCROW, &
                       SPCPROW, TSPCROW, RHSIROW, FCLOROW, DLONROW, &
                       GGEOROW, GUSTROL, REFROT, BCSNROT, DEPBROW, &
                       DLATROW, maxAnnualActLyrROT)

      ! FLAG New for Jason's Albedo
      ! Calculate snow fall flux in kg/m2/s
      if (isnoalb == 1) then
        do K = 1,NML
          PCSNGAT(K) = SPCPGAT(K) * RHSIGAT(K)  ! snowfall rate * density of fresh snow
        end do
      end if
      ! End flag

      !    * INITIALIZATION OF DIAGNOSTIC VARIABLES SPLIT OUT OF classGather
      !    * FOR CONSISTENCY WITH GCM APPLICATIONS.
      call initDiagnosticVars(nml, ilg)

      !========================================================================
      
      !> energyWaterBalanceCheck does the initial calculations for the energy and water balance checks

      call energyWaterBalanceCheck(0, CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP, & ! Formerly CLASSZ
                                   WTVSTP, WTSSTP, WTGSTP, &
                                   FSGVGAT, FLGVGAT, HFSCGAT, HEVCGAT, HMFCGAT, HTCCGAT, &
                                   FSGSGAT, FLGSGAT, HFSSGAT, HEVSGAT, HMFNGAT, HTCSGAT, &
                                   FSGGGAT, FLGGGAT, HFSGGAT, HEVGGAT, HMFGGAT, HTCGAT, &
                                   PCFCGAT, PCLCGAT, QFCFGAT, QFCLGAT, ROFCGAT, WTRCGAT, &
                                   PCPNGAT, QFNGAT, ROFNGAT, WTRSGAT, PCPGGAT, QFGGAT, &
                                   QFCGAT, ROFGAT, WTRGGAT, CMAIGAT, RCANGAT, SCANGAT, &
                                   TCANGAT, SNOGAT, WSNOGAT, TSNOGAT, THLQGAT, THICGAT, &
                                   HCPSGAT, THPGAT, DLZWGAT, TBARGAT, ZPNDGAT, TPNDGAT, &
                                   DELZ, FCS, FGS, FC, FG, &
                                   1, NML, ILG, IGND, N)

      ! ctemg2 takes variables in the 'row' format (nlat, nmos, ...)
      ! and converts them to the 'gat' format (ilg, ...). At present
      ! ctemg2 is bloated with many variables that do not require
      ! gathering. This subroutine should ideally be only used for
      ! state variables that are updated from external files as
      ! the run progresses. Since the model calculations operate
      ! on the 'gat' form, any other variables need not be gathered
      ! as they will already be in the correct format from the previous
      ! model timestep.
      call ctemg2(fcancmxgat, rmatcgat, zolncgat, paicgat, &
                  ailcgat, ailcggat, cmasvegcgat, slaicgat, &
                  ailcgsgat, fcancsgat, fcancgat, rmatctemgat, &
                  co2concgat, co2i1cggat, co2i1csgat, co2i2cggat, &
                  co2i2csgat, xdiffusgat, slaigat, cfluxcggat, &
                  cfluxcsgat, ancsveggat, ancgveggat, rmlcsveggat, &
                  rmlcgveggat, canresgat, sdepgat, ch4concgat, &
                  sandgat, claygat, orgmgat, &
                  anveggat, rmlveggat, tbaraccgat_t, prbfrhucgat, &
                  extnprobgat, pfcancmxgat, nfcancmxgat, &
                  stemmassgat, rootmassgat, litrmassgat, gleafmasgat, &
                  bleafmasgat, soilcmasgat, ailcbgat, flhrlossgat, &
                  pandaysgat, lfstatusgat, grwtheffgat, lystmmasgat, &
                  lyrotmasgat, tymaxlaigat, vgbiomasgat, gavgltmsgat, &
                  stmhrlosgat, bmasveggat, colddaysgat, rothrlosgat, &
                  alvsctmgat, alirctmgat, gavglaigat, nppgat, &
                  nepgat, hetroresgat, autoresgat, soilcrespgat, &
                  rmgat, rggat, nbpgat, litresgat, &
                  socresgat, gppgat, dstcemlsgat, litrfallgat, &
                  humiftrsgat, veghghtgat, rootdpthgat, rmlgat, &
                  rmsgat, rmrgat, tltrleafgat, tltrstemgat, &
                  tltrrootgat, leaflitrgat, roottempgat, afrleafgat, &
                  afrstemgat, afrrootgat, wtstatusgat, ltstatusgat, &
                  burnfracgat, smfuncveggat, lucemcomgat, lucltringat, &
                  lucsocingat, dstcemls3gat, popdingat, &
                  faregat, gavgscmsgat, rmlvegaccgat, pftexistgat, &
                  rmsveggat, rmrveggat, rgveggat, vgbiomas_veggat, &
                  gppveggat, nepveggat, &
                  emit_co2gat, emit_cogat, emit_ch4gat, emit_nmhcgat, &
                  emit_h2gat, emit_noxgat, emit_n2ogat, emit_pm25gat, &
                  emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &
                  btermgat, ltermgat, mtermgat, daylgat, dayl_maxgat, &
                  nbpveggat, hetroresveggat, autoresveggat, litresveggat, &
                  soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat, &
                  ch4WetSpecgat, slopefracgat, &
                  wetfdyngat, ch4WetDyngat, ch4soillsgat, &
                  twarmmgat, tcoldmgat, gdd5gat, &
                  ariditygat, srplsmongat, defctmongat, anndefctgat, &
                  annsrplsgat, annpcpgat, dry_season_lengthgat, &
                  anmossgat, rmlmossgat, gppmossgat, armossgat, nppmossgat, &
                  litrmsmossgat, peatdepgat, Cmossmasgat, dmossgat, & ! thlqaccgat_m, &
                  ipeatlandgat, pddgat, tracerCO2gat, &
                  !          thicaccgat_m, ipeatlandgat, pddgat, & this line commented out.
                  ilmos, jlmos, iwmos, jwmos, &
                  nml, fcancmxrow, rmatcrow, zolncrow, paicrow, &
                  ailcrow, ailcgrow, cmasvegcrow, slaicrow, &
                  ailcgsrow, fcancsrow, fcancrow, rmatctemrow, &
                  co2concrow, co2i1cgrow, co2i1csrow, co2i2cgrow, &
                  co2i2csrow, xdiffus, slairow, cfluxcgrow, &
                  cfluxcsrow, ancsvegrow, ancgvegrow, rmlcsvegrow, &
                  rmlcgvegrow, canresrow, SDEPROT, ch4concrow, &
                  SANDROT, CLAYROT, ORGMROT, &
                  anvegrow, rmlvegrow, tbaraccrow_m, prbfrhucrow, &
                  extnprobrow, pfcancmxrow, nfcancmxrow, &
                  stemmassrow, rootmassrow, litrmassrow, gleafmasrow, &
                  bleafmasrow, soilcmasrow, ailcbrow, flhrlossrow, &
                  pandaysrow, lfstatusrow, grwtheffrow, lystmmasrow, &
                  lyrotmasrow, tymaxlairow, vgbiomasrow, gavgltmsrow, &
                  stmhrlosrow, bmasvegrow, colddaysrow, rothrlosrow, &
                  alvsctmrow, alirctmrow, gavglairow, npprow, &
                  neprow, hetroresrow, autoresrow, soilcresprow, &
                  rmrow, rgrow, nbprow, litresrow, &
                  socresrow, gpprow, dstcemlsrow, litrfallrow, &
                  humiftrsrow, veghghtrow, rootdpthrow, rmlrow, &
                  rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                  tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                  afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                  burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                  lucsocinrow, dstcemls3row, popdinrow, &
                  FAREROT, gavgscmsrow, rmlvegaccrow, pftexistrow, &
                  rmsvegrow, rmrvegrow, rgvegrow, vgbiomas_vegrow, &
                  gppvegrow, nepvegrow, &
                  emit_co2row, emit_corow, emit_ch4row, emit_nmhcrow, &
                  emit_h2row, emit_noxrow, emit_n2orow, emit_pm25row, &
                  emit_tpmrow, emit_tcrow, emit_ocrow, emit_bcrow, &
                  btermrow, ltermrow, mtermrow, daylrow, dayl_maxrow, &
                  nbpvegrow, hetroresvegrow, autoresvegrow, litresvegrow, &
                  soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow, &
                  ch4WetSpecrow, slopefracrow, &
                  wetfdynrow, ch4WetDynrow, ch4soillsrow, &
                  twarmmrow, tcoldmrow, gdd5row, &
                  aridityrow, srplsmonrow, defctmonrow, anndefctrow, &
                  annsrplsrow, annpcprow, dry_season_lengthrow, &
                  anmossrow, rmlmossrow, gppmossrow, armossrow, nppmossrow, &
                  litrmsmossrow, peatdeprow, Cmossmasrow, dmossrow, &
                  ipeatlandrow, pddrow, tracerCO2rot)
      !    5      thlqaccrow_m, thicaccrow_m, ipeatlandrow, pddrow) this line commented out.

      !-----------------------------------------------------------------------
      !* ALBEDO AND TRANSMISSIVITY CALCULATIONS; GENERAL VEGETATION
      !* CHARACTERISTICS.

      !     * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM by including: zolnc,
      !     * cmasvegc, alvsctm, alirctm, ipeatlandgat in the arguments.

      !> radiationDriver manages the calculation of albedos and other surface parameters

      call radiationDriver(FC, FG, FCS, FGS, ALVSCN, ALIRCN, & ! Formerly CLASSA
                           ALVSG, ALIRG, ALVSCS, ALIRCS, ALVSSN, ALIRSN, &
                           ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN, &
                           TRVSCS, TRIRCS, FSVF, FSVFS, &
                           RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC, &
                           FRAICS, FSNOCS, DISP, DISPS, ZOMLNC, ZOMLCS, &
                           ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, &
                           CHCAP, CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP, &
                           CWLCPS, CWFCPS, RC, RCS, RBCOEF, FROOT, &
                           FROOTS, ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZSNOW, &
                           WSNOGAT, ALVSGAT, ALIRGAT, HTCCGAT, HTCSGAT, HTCGAT, &
                           ALTG, ALSNO, TRSNOWC, TRSNOWG, &
                           WTRCGAT, WTRSGAT, WTRGGAT, CMAIGAT, FSNOGAT, &
                           FCANGAT, LNZ0GAT, ALVCGAT, ALICGAT, PAMXGAT, PAMNGAT, &
                           CMASGAT, ROOTGAT, RSMNGAT, QA50GAT, VPDAGAT, VPDBGAT, &
                           PSGAGAT, PSGBGAT, PAIDGAT, HGTDGAT, ACVDGAT, ACIDGAT, &
                           ASVDGAT, ASIDGAT, AGVDGAT, AGIDGAT, &
                           ALGWVGAT, ALGWNGAT, ALGDVGAT, ALGDNGAT, &
                           THLQGAT, THICGAT, TBARGAT, RCANGAT, SCANGAT, TCANGAT, &
                           GROGAT, SNOGAT, TSNOGAT, RHOSGAT, ALBSGAT, ZBLDGAT, &
                           Z0ORGAT, ZSNLGAT, ZPLGGAT, ZPLSGAT, &
                           FCLOGAT, TAGAT, VPDGAT, RHOAGAT, CSZGAT, &
                           FSDBGAT, FSFBGAT, REFGAT, BCSNGAT, &
                           FSVHGAT, RADJGAT, DLONGAT, RHSIGAT, DELZ, DLZWGAT, &
                           ZBTWGAT, THPGAT, THMGAT, PSISGAT, BIGAT, PSIWGAT, &
                           HCPSGAT, ISNDGAT, &
                           FCANCMXGAT, ICC, ctem_on, RMATCGAT, ZOLNCGAT, &
                           CMASVEGCGAT, AILCGAT, PAICGAT, NOL2PFTS, &
                           SLAICGAT, AILCGGAT, AILCGSGAT, FCANCGAT, FCANCSGAT, &
                           IDAY, ILG, 1, NML, NBS, &
                           JLAT, N, ICAN, ICAN + 1, IGND, IDISP, IZREF, &
                           IWF, IPAI, IHGT, IALC, IALS, IALG, &
                           ISNOALB, alvsctmgat, alirctmgat, ipeatlandgat)

      !-----------------------------------------------------------------------
      !          * SURFACE TEMPERATURE AND FLUX CALCULATIONS.

      !          * ADAPTED TO COUPLING OF CLASS3.6 AND CTEM
      !          * by including in the arguments: lfstatus

      !> energyBudgetDriver calls the subroutines associated with the surface energy balance calculations

      call energyBudgetDriver(TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, & ! Formerly CLASST
                              THICEC, THICEG, HCPC, HCPG, TCTOPC, TCBOTC, TCTOPG, TCBOTG, &
                              GZEROC, GZEROG, GZROCS, GZROGS, G12C, G12G, G12CS, G12GS, &
                              G23C, G23G, G23CS, G23GS, QFREZC, QFREZG, QMELTC, QMELTG, &
                              EVAPC, EVAPCG, EVAPG, EVAPCS, EVPCSG, EVAPGS, TCANO, TCANS, &
                              RAICAN, SNOCAN, RAICNS, SNOCNS, CHCAP, CHCAPS, TPONDC, TPONDG, &
                              TPNDCS, TPNDGS, TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                              ITCTGAT, CDHGAT, CDMGAT, HFSGAT, TFXGAT, QEVPGAT, QFSGAT, &
                              PETGAT, GAGAT, EFGAT, GTGAT, QGGAT, &
                              SFCTGAT, SFCUGAT, SFCVGAT, SFCQGAT, SFRHGAT, &
                              GTBS, SFCUBS, SFCVBS, USTARBS, &
                              FSGVGAT, FSGSGAT, FSGGGAT, FLGVGAT, FLGSGAT, FLGGGAT, &
                              HFSCGAT, HFSSGAT, HFSGGAT, HEVCGAT, HEVSGAT, HEVGGAT, HMFCGAT, HMFNGAT, &
                              HTCCGAT, HTCSGAT, HTCGAT, QFCFGAT, QFCLGAT, DRGAT, wtableGAT, ILMOGAT, &
                              UEGAT, HBLGAT, TACGAT, QACGAT, ZRFMGAT, ZRFHGAT, ZDMGAT, ZDHGAT, &
                              VPDGAT, TADPGAT, RHOAGAT, FSVHGAT, FSIHGAT, FDLGAT, ULGAT, VLGAT, &
                              TAGAT, QAGAT, PADRGAT, FC, FG, FCS, FGS, RBCOEF, &
                              FSVF, FSVFS, PRESGAT, VMODGAT, ALVSCN, ALIRCN, ALVSG, ALIRG, &
                              ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, ALIRSC, &
                              TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC, RCS, WTRGGAT, groundHeatFlux, QLWOGAT, &
                              FRAINC, FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP, DISPS, &
                              ZOMLNC, ZOELNC, ZOMLNG, ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS, &
                              TBARGAT, THLQGAT, THICGAT, TPNDGAT, ZPNDGAT, TBASGAT, TCANGAT, TSNOGAT, &
                              ZSNOW, RHOSGAT, WSNOGAT, THPGAT, THRGAT, THMGAT, THFCGAT, THLWGAT, &
                              TRSNOWC, TRSNOWG, ALSNO, FSSBGAT, FROOT, FROOTS, &
                              RADJGAT, PREGAT, HCPSGAT, TCSGAT, TSFSGAT, DELZ, DLZWGAT, ZBTWGAT, &
                              FTEMP, FVAP, RIB, ISNDGAT, &
                              AILCGGAT, AILCGSGAT, FCANCGAT, FCANCSGAT, CO2CONCGAT, CO2I1CGGAT, &
                              CO2I1CSGAT, CO2I2CGGAT, CO2I2CSGAT, CSZGAT, XDIFFUSGAT, SLAIGAT, ICC, &
                              ctem_on, RMATCTEMGAT, FCANCMXGAT, L2MAX, NOL2PFTS, CFLUXCGGAT, &
                              CFLUXCSGAT, ANCSVEGGAT, ANCGVEGGAT, RMLCSVEGGAT, RMLCGVEGGAT, &
                              TCSNOW, GSNOW, ITC, ITCG, ITG, ILG, 1, NML, JLAT, N, ICAN, &
                              IGND, IZREF, ISLFD, NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI, &
                              NBS, ISNOALB, daylgat, dayl_maxgat, &
                              ipeatlandgat, ancsmoss, angsmoss, ancmoss, angmoss, &
                              rmlcsmoss, rmlgsmoss, rmlcmoss, rmlgmoss, &
                              Cmossmasgat, dmossgat, iday, pddgat)

      !-----------------------------------------------------------------------
      !          * WATER BUDGET CALCULATIONS.

      !> waterBudgetDriver calls the subroutines associated with the surface water balance calculations

      call waterBudgetDriver(THLQGAT, THICGAT, TBARGAT, TCANGAT, RCANGAT, SCANGAT, & ! Formerly CLASSW
                             ROFGAT, TROFGAT, SNOGAT, TSNOGAT, RHOSGAT, ALBSGAT, &
                             WSNOGAT, ZPNDGAT, TPNDGAT, GROGAT, TBASGAT, GFLXGAT, &
                             PCFCGAT, PCLCGAT, PCPNGAT, PCPGGAT, QFCFGAT, QFCLGAT, &
                             QFNGAT, QFGGAT, QFCGAT, HMFCGAT, HMFGGAT, HMFNGAT, &
                             HTCCGAT, HTCSGAT, HTCGAT, ROFCGAT, ROFNGAT, ROVGGAT, &
                             WTRSGAT, WTRGGAT, ROFOGAT, ROFSGAT, ROFBGAT, &
                             TROOGAT, TROSGAT, TROBGAT, QFSGAT, QFXGAT, RHOAGAT, &
                             TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, &
                             THICEC, THICEG, HCPC, HCPG, RPCPGAT, TRPCGAT, &
                             SPCPGAT, TSPCGAT, PREGAT, TAGAT, RHSIGAT, GGEOGAT, &
                             FC, FG, FCS, FGS, TPONDC, TPONDG, &
                             TPNDCS, TPNDGS, EVAPC, EVAPCG, EVAPG, EVAPCS, &
                             EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG, &
                             RAICAN, SNOCAN, RAICNS, SNOCNS, FSVF, FSVFS, &
                             CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO, &
                             TCANS, CHCAP, CHCAPS, CMASSC, CMASCS, ZSNOW, &
                             GZEROC, GZEROG, GZROCS, GZROGS, G12C, G12G, &
                             G12CS, G12GS, G23C, G23G, G23CS, G23GS, &
                             TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                             ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TSFSGAT, &
                             TCTOPC, TCBOTC, TCTOPG, TCBOTG, FROOT, FROOTS, &
                             THPGAT, THRGAT, THMGAT, BIGAT, PSISGAT, GRKSGAT, &
                             THRAGAT, THFCGAT, DRNGAT, HCPSGAT, DELZ, &
                             DLZWGAT, ZBTWGAT, XSLPGAT, GRKFGAT, WFSFGAT, WFCIGAT, &
                             ISNDGAT, IGDRGAT, &
                             IWF, ILG, 1, NML, N, &
                             JLAT, ICAN, IGND, IGND + 1, IGND + 2, &
                             NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI)

      !========================================================================

      !> energyWaterBalanceCheck completes the energy and water balance checks for the current time step

      call energyWaterBalanceCheck(1, CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP, & ! Formerly CLASSZ
                                   WTVSTP, WTSSTP, WTGSTP, &
                                   FSGVGAT, FLGVGAT, HFSCGAT, HEVCGAT, HMFCGAT, HTCCGAT, &
                                   FSGSGAT, FLGSGAT, HFSSGAT, HEVSGAT, HMFNGAT, HTCSGAT, &
                                   FSGGGAT, FLGGGAT, HFSGGAT, HEVGGAT, HMFGGAT, HTCGAT, &
                                   PCFCGAT, PCLCGAT, QFCFGAT, QFCLGAT, ROFCGAT, WTRCGAT, &
                                   PCPNGAT, QFNGAT, ROFNGAT, WTRSGAT, PCPGGAT, QFGGAT, &
                                   QFCGAT, ROFGAT, WTRGGAT, CMAIGAT, RCANGAT, SCANGAT, &
                                   TCANGAT, SNOGAT, WSNOGAT, TSNOGAT, THLQGAT, THICGAT, &
                                   HCPSGAT, THPGAT, DLZWGAT, TBARGAT, ZPNDGAT, TPNDGAT, &
                                   DELZ, FCS, FGS, FC, FG, &
                                   1, NML, ILG, IGND, N)

      if (ctem_on) then

        !> Accumulate variables not already accumulated but which are required by CTEM.
        call accumulateForCTEM(nml,ILMOS)

        if (ncount == nday) then

          ! Find daily averages of accumulated variables for CTEM
          call dayEndCTEMPreparation(nml, nday, ILMOS)

          ! Call Canadian Terrestrial Ecosystem Model which operates at a daily time step,
          ! and uses daily accumulated values of variables simulated by CLASS.
          call ctem(fsnowacc_t, sandgat, ILMOS, & ! In
                    ilg, 1, nml, iday, radjgat, &! In
                    taaccgat_t, dlzwgat, ancgvgac_t, rmlcgvga_t, & ! In
                    zbtwgat, doMethane, & ! In
                    uvaccgat_t, vvaccgat_t, lightng, tbaraccgat_t, &! In
                    sdepgat, spinfast, todfrac, & ! In
                    netrad_gat, preacc_gat, PSISGAT, &! In
                    grclarea, popdingat, isndgat, &! In
                    wetfrac_presgat, slopefracgat, BIGAT, &! In
                    THPGAT, DLATGAT, ch4concgat, &! In
                    THFCGAT, THLWGAT, thliqacc_t, thiceacc_t, &! In
                    ipeatlandgat, anmossac_t, rmlmossac_t, gppmossac_t, &! In
                    wtablegat, maxAnnualActLyrGAT, & ! In
                    PFTCompetition, dofire, lnduseon, inibioclim, & ! In
                    leapnow, useTracer, tracerCO2gat, &! In
                    pfcancmxgat, nfcancmxgat, & ! In
                    stemmassgat, rootmassgat, litrmassgat, gleafmasgat, & ! In/Out
                    bleafmasgat, soilcmasgat, ailcggat, ailcgat, & ! In/Out
                    zolncgat, rmatctemgat, rmatcgat, ailcbgat, & ! In/Out
                    flhrlossgat, pandaysgat, lfstatusgat, grwtheffgat, & ! In/Out
                    lystmmasgat, lyrotmasgat, tymaxlaigat, vgbiomasgat, & ! In/Out
                    gavgltmsgat, gavgscmsgat, stmhrlosgat, slaigat, & ! In/Out
                    bmasveggat, cmasvegcgat, colddaysgat, rothrlosgat, & ! In/Out
                    fcangat, alvsctmgat, alirctmgat, gavglaigat, &! In/Out
                    Cmossmasgat, litrmsmossgat, peatdepgat, fcancmxgat, &! In/Out
                    geremortgat, intrmortgat, pstemmassgat, pgleafmassgat, &! In/Out
                    tcurm, srpcuryr, dftcuryr, lambdagat, &! In/Out
                    tmonth, anpcpcur, anpecur, gdd5cur, &! In/Out
                    surmncur, defmncur, srplscur, defctcur, &! In/Out
                    ariditygat, srplsmongat, defctmongat, anndefctgat, &! In/Out
                    annsrplsgat, annpcpgat, dry_season_lengthgat, &! In/Out
                    pftexistgat, twarmmgat, tcoldmgat, gdd5gat, nppveggat, &! In/Out
                    tracerStemMassgat, tracerRootMassgat, tracerGLeafMassgat, tracerBLeafMassgat, & ! In/Out
                    tracerSoilCMassgat, tracerLitrMassgat, tracerMossCMassgat, tracerMossLitrMassgat, & ! In/Out
                    nppgat, nepgat, hetroresgat, autoresgat, &! Out (Primary)
                    soilcrespgat, rmgat, rggat, nbpgat, &! Out (Primary)
                    litresgat, socresgat, gppgat, dstcemlsgat, &! Out (Primary)
                    litrfallgat, humiftrsgat, veghghtgat, rootdpthgat, &! Out (Primary)
                    rmlgat, rmsgat, rmrgat, tltrleafgat, &! Out (Primary)
                    tltrstemgat, tltrrootgat, leaflitrgat, roottempgat, &! Out (Primary)
                    burnfracgat, lucemcomgat, lucltringat, &! Out (Primary)
                    lucsocingat, dstcemls3gat, &! Out (Primary)
                    ch4WetSpecgat, ch4WetDyngat, wetfdyngat, ch4soillsgat, &! Out (Primary)
                    paicgat, slaicgat, &! Out (Primary)
                    emit_co2gat, emit_ch4gat, reprocost, blfltrdt, glfltrdt, &! Out (Primary)
                    glcaemls, blcaemls, rtcaemls, stcaemls, ltrcemls, &  ! Out (Primary)
                    ntchlveg, ntchsveg, ntchrveg, &  ! Out (Primary)
                    emit_cogat, emit_nmhcgat, smfuncveggat, &! Out (Secondary)
                    emit_h2gat, emit_noxgat, emit_n2ogat, emit_pm25gat, &! Out (Secondary)
                    emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &! Out (Secondary)
                    btermgat, ltermgat, mtermgat, burnvegfgat, &! Out (Secondary)
                    litrfallveggat, humiftrsveggat, ltstatusgat,  &! Out (Secondary)
                    afrleafgat, afrstemgat, afrrootgat, wtstatusgat, &! Out (Secondary)
                    rmlvegaccgat, rmsveggat, rmrveggat, rgveggat, &! Out (Secondary)
                    vgbiomas_veggat, gppveggat, nepveggat, nbpveggat, &! Out (Secondary)
                    hetroresveggat, autoresveggat, litresveggat, soilcresveggat, &! Out (Secondary)
                    nppmossgat, armossgat, &! Out (Secondary)
                    ccgat, mmgat)

          !     reset mosaic accumulator arrays. These are scattered in ctems2 so we need
          !     to reset here,prior to ctems2.
          do i = 1,nml
            vvaccgat_t(i) = 0.0  !
            uvaccgat_t(i) = 0.0  !
            do j = 1,ignd
              tbaraccgat_t(i,j) = 0.0 !
            end do
          end do

          ! Once a year,calculate the 14C lost to decay if using the 14C tracer.
          if (useTracer == 2 .and. &
              iday == lastdoy .and. ncount == nday) call decay14C(1,nml)

        end if  ! if (ncount==nday)
      end if  ! if (ctem_on)

      if (isnoalb == 1) then
        call fourBandDriver(NML, SNOGAT, GCGAT, TSNOGAT, & ! in
                            ZSNOW, TCSNOW, GSNOW, PCSNGAT, WSNOGAT, & ! in
                            ALBSGAT, RHOSGAT, & ! in/out
                            TZSGAT, REFGAT) ! out
      end if

      !> classScatter performs the scatter operation, scattering the variables from
      !> the long vectors of mosaic tiles back onto the configuration of mosaic tiles within grid cells.

      call classScatter(TBARROT, THLQROT, THICROT, TSFSROT, TPNDROT, & ! Formerly CLASSS
                        ZPNDROT, TBASROT, ALBSROT, TSNOROT, RHOSROT, &
                        SNOROT, GTROT, TCANROT, RCANROT, SCANROT, &
                        GROROT, CMAIROT, TACROT, QACROT, WSNOROT, &
                        REFROT, BCSNROT, EMISROT, SALBROT, CSALROT, &
                        groundHeatFluxROT, &
                        ILMOS, JLMOS, NML, NLAT, NTLD, NMOS, &
                        ILG, IGND, ICAN, ICAN + 1, NBS, &
                        TBARGAT, THLQGAT, THICGAT, TSFSGAT, TPNDGAT, &
                        ZPNDGAT, TBASGAT, ALBSGAT, TSNOGAT, RHOSGAT, &
                        SNOGAT, GTGAT, TCANGAT, RCANGAT, SCANGAT, &
                        GROGAT, CMAIGAT, TACGAT, QACGAT, WSNOGAT, &
                        REFGAT, BCSNGAT, EMISGAT, SALBGAT, CSALGAT,&
                        groundHeatFlux)

      !
      !    * SCATTER OPERATION ON DIAGNOSTIC VARIABLES SPLIT OUT OF
      !    * classScatter FOR CONSISTENCY WITH GCM APPLICATIONS.
      !
      do K = 1,NML
        CDHROT (ILMOS(K),JLMOS(K)) = CDHGAT (K)
        CDMROT (ILMOS(K),JLMOS(K)) = CDMGAT (K)
        HFSROT (ILMOS(K),JLMOS(K)) = HFSGAT (K)
        TFXROT (ILMOS(K),JLMOS(K)) = TFXGAT (K)
        QEVPROT(ILMOS(K),JLMOS(K)) = QEVPGAT(K)
        QFSROT (ILMOS(K),JLMOS(K)) = QFSGAT (K)
        QFXROT (ILMOS(K),JLMOS(K)) = QFXGAT (K)
        PETROT (ILMOS(K),JLMOS(K)) = PETGAT (K)
        GAROT  (ILMOS(K),JLMOS(K)) = GAGAT  (K)
        EFROT  (ILMOS(K),JLMOS(K)) = EFGAT  (K)
        QGROT  (ILMOS(K),JLMOS(K)) = QGGAT  (K)
        ALVSROT(ILMOS(K),JLMOS(K)) = ALVSGAT(K)
        ALIRROT(ILMOS(K),JLMOS(K)) = ALIRGAT(K)
        SFCTROT(ILMOS(K),JLMOS(K)) = SFCTGAT(K)
        SFCUROT(ILMOS(K),JLMOS(K)) = SFCUGAT(K)
        SFCVROT(ILMOS(K),JLMOS(K)) = SFCVGAT(K)
        SFCQROT(ILMOS(K),JLMOS(K)) = SFCQGAT(K)
        FSNOROT(ILMOS(K),JLMOS(K)) = FSNOGAT(K)
        FSGVROT(ILMOS(K),JLMOS(K)) = FSGVGAT(K)
        FSGSROT(ILMOS(K),JLMOS(K)) = FSGSGAT(K)
        FSGGROT(ILMOS(K),JLMOS(K)) = FSGGGAT(K)
        FLGVROT(ILMOS(K),JLMOS(K)) = FLGVGAT(K)
        FLGSROT(ILMOS(K),JLMOS(K)) = FLGSGAT(K)
        FLGGROT(ILMOS(K),JLMOS(K)) = FLGGGAT(K)
        HFSCROT(ILMOS(K),JLMOS(K)) = HFSCGAT(K)
        HFSSROT(ILMOS(K),JLMOS(K)) = HFSSGAT(K)
        HFSGROT(ILMOS(K),JLMOS(K)) = HFSGGAT(K)
        HEVCROT(ILMOS(K),JLMOS(K)) = HEVCGAT(K)
        HEVSROT(ILMOS(K),JLMOS(K)) = HEVSGAT(K)
        HEVGROT(ILMOS(K),JLMOS(K)) = HEVGGAT(K)
        HMFCROT(ILMOS(K),JLMOS(K)) = HMFCGAT(K)
        HMFNROT(ILMOS(K),JLMOS(K)) = HMFNGAT(K)
        HTCCROT(ILMOS(K),JLMOS(K)) = HTCCGAT(K)
        HTCSROT(ILMOS(K),JLMOS(K)) = HTCSGAT(K)
        PCFCROT(ILMOS(K),JLMOS(K)) = PCFCGAT(K)
        PCLCROT(ILMOS(K),JLMOS(K)) = PCLCGAT(K)
        PCPNROT(ILMOS(K),JLMOS(K)) = PCPNGAT(K)
        PCPGROT(ILMOS(K),JLMOS(K)) = PCPGGAT(K)
        QFGROT (ILMOS(K),JLMOS(K)) = QFGGAT (K)
        QFNROT (ILMOS(K),JLMOS(K)) = QFNGAT (K)
        QFCLROT(ILMOS(K),JLMOS(K)) = QFCLGAT(K)
        QFCFROT(ILMOS(K),JLMOS(K)) = QFCFGAT(K)
        ROFROT (ILMOS(K),JLMOS(K)) = ROFGAT (K)
        ROFOROT(ILMOS(K),JLMOS(K)) = ROFOGAT(K)
        ROFSROT(ILMOS(K),JLMOS(K)) = ROFSGAT(K)
        ROFBROT(ILMOS(K),JLMOS(K)) = ROFBGAT(K)
        TROFROT(ILMOS(K),JLMOS(K)) = TROFGAT(K)
        TROOROT(ILMOS(K),JLMOS(K)) = TROOGAT(K)
        TROSROT(ILMOS(K),JLMOS(K)) = TROSGAT(K)
        TROBROT(ILMOS(K),JLMOS(K)) = TROBGAT(K)
        ROFCROT(ILMOS(K),JLMOS(K)) = ROFCGAT(K)
        ROFNROT(ILMOS(K),JLMOS(K)) = ROFNGAT(K)
        ROVGROT(ILMOS(K),JLMOS(K)) = ROVGGAT(K)
        WTRCROT(ILMOS(K),JLMOS(K)) = WTRCGAT(K)
        WTRSROT(ILMOS(K),JLMOS(K)) = WTRSGAT(K)
        WTRGROT(ILMOS(K),JLMOS(K)) = WTRGGAT(K)
        DRROT  (ILMOS(K),JLMOS(K)) = DRGAT  (K)
        wtableROT(ILMOS(K),JLMOS(K)) = wtableGAT(K)
        ILMOROT(ILMOS(K),JLMOS(K)) = ILMOGAT(K)
        UEROT  (ILMOS(K),JLMOS(K)) = UEGAT(K)
        HBLROT (ILMOS(K),JLMOS(K)) = HBLGAT(K)
      end do ! loop 380

      do L = 1,IGND
        do K = 1,NML
          HMFGROT(ILMOS(K),JLMOS(K),L) = HMFGGAT(K,L)
          HTCROT (ILMOS(K),JLMOS(K),L) = HTCGAT (K,L)
          QFCROT (ILMOS(K),JLMOS(K),L) = QFCGAT (K,L)
          GFLXROT(ILMOS(K),JLMOS(K),L) = GFLXGAT(K,L)
        end do
      end do ! loop 390

      do M = 1,50
        do L = 1,6
          do K = 1,NML
            ITCTROT(ILMOS(K),JLMOS(K),L,M) = ITCTGAT(K,L,M)
          end do ! loop 410
        end do ! loop 420
      end do ! loop 430

      ! ctems2 converts variables from the 'gat' format to the
      ! 'row' format, which is suitable for writing to output/restart
      ! files. If a variable is not written to either of those files,
      ! there is no need to scatter the variable as it will be in the
      ! correct format for model calclations ('gat').
      call ctems2(fcancmxrow, rmatcrow, zolncrow, paicrow, &
                  ailcrow, ailcgrow, cmasvegcrow, slaicrow, &
                  ailcgsrow, fcancsrow, fcancrow, rmatctemrow, &
                  co2concrow, co2i1cgrow, co2i1csrow, co2i2cgrow, &
                  co2i2csrow, xdiffus, slairow, cfluxcgrow, &
                  cfluxcsrow, ancsvegrow, ancgvegrow, rmlcsvegrow, &
                  rmlcgvegrow, canresrow, SDEPROT, ch4concrow, &
                  SANDROT, CLAYROT, ORGMROT, &
                  anvegrow, rmlvegrow, tbaraccrow_m, prbfrhucrow, &
                  extnprobrow, pfcancmxrow, nfcancmxrow, &
                  stemmassrow, rootmassrow, litrmassrow, gleafmasrow, &
                  bleafmasrow, soilcmasrow, ailcbrow, flhrlossrow, &
                  pandaysrow, lfstatusrow, grwtheffrow, lystmmasrow, &
                  lyrotmasrow, tymaxlairow, vgbiomasrow, gavgltmsrow, &
                  stmhrlosrow, bmasvegrow, colddaysrow, rothrlosrow, &
                  alvsctmrow, alirctmrow, gavglairow, npprow, &
                  neprow, hetroresrow, autoresrow, soilcresprow, &
                  rmrow, rgrow, nbprow, litresrow, &
                  socresrow, gpprow, dstcemlsrow, litrfallrow, &
                  humiftrsrow, veghghtrow, rootdpthrow, rmlrow, &
                  litrfallvegrow, humiftrsvegrow, &
                  rmsrow, rmrrow, tltrleafrow, tltrstemrow, &
                  tltrrootrow, leaflitrrow, roottemprow, afrleafrow, &
                  afrstemrow, afrrootrow, wtstatusrow, ltstatusrow, &
                  burnfracrow, smfuncvegrow, lucemcomrow, lucltrinrow, &
                  lucsocinrow, nppvegrow, dstcemls3row, &
                  FAREROT, gavgscmsrow, &
                  rmlvegaccrow, rmsvegrow, rmrvegrow, rgvegrow, &
                  vgbiomas_vegrow, gppvegrow, nepvegrow, &
                  FCANROT, pftexistrow, &
                  emit_co2row, emit_corow, emit_ch4row, emit_nmhcrow, &
                  emit_h2row, emit_noxrow, emit_n2orow, emit_pm25row, &
                  emit_tpmrow, emit_tcrow, emit_ocrow, emit_bcrow, &
                  btermrow, ltermrow, mtermrow, &
                  nbpvegrow, hetroresvegrow, autoresvegrow, litresvegrow, &
                  soilcresvegrow, burnvegfrow, pstemmassrow, pgleafmassrow, &
                  ch4WetSpecrow, wetfdynrow, ch4WetDynrow, ch4soillsrow, &
                  twarmmrow, tcoldmrow, gdd5row, &
                  aridityrow, srplsmonrow, defctmonrow, anndefctrow, &
                  annsrplsrow, annpcprow, dry_season_lengthrow, &
                  anmossrow, rmlmossrow, gppmossrow, armossrow, nppmossrow, &
                  peatdeprow, litrmsmossrow, Cmossmasrow, dmossrow, &
                  ipeatlandrow, pddrow, wetfrac_presrow, &
                  tracergLeafMassrot, tracerBLeafMassrot, tracerStemMassrot, &
                  tracerRootMassrot, tracerLitrMassrot, tracerSoilCMassrot, &
                  tracerMossCMassrot, tracerMossLitrMassrot, &
                  !    ----
                  ilmos, jlmos, iwmos, jwmos, &
                  nml, fcancmxgat, rmatcgat, zolncgat, paicgat, &
                  ailcgat, ailcggat, cmasvegcgat, slaicgat, &
                  ailcgsgat, fcancsgat, fcancgat, rmatctemgat, &
                  co2concgat, co2i1cggat, co2i1csgat, co2i2cggat, &
                  co2i2csgat, xdiffusgat, slaigat, cfluxcggat, &
                  cfluxcsgat, ancsveggat, ancgveggat, rmlcsveggat, &
                  rmlcgveggat, canresgat, sdepgat, ch4concgat, &
                  sandgat, claygat, orgmgat, &
                  anveggat, rmlveggat, tbaraccgat_t, prbfrhucgat, &
                  extnprobgat, pfcancmxgat, nfcancmxgat, &
                  stemmassgat, rootmassgat, litrmassgat, gleafmasgat, &
                  bleafmasgat, soilcmasgat, ailcbgat, flhrlossgat, &
                  pandaysgat, lfstatusgat, grwtheffgat, lystmmasgat, &
                  lyrotmasgat, tymaxlaigat, vgbiomasgat, gavgltmsgat, &
                  stmhrlosgat, bmasveggat, colddaysgat, rothrlosgat, &
                  alvsctmgat, alirctmgat, gavglaigat, nppgat, &
                  nepgat, hetroresgat, autoresgat, soilcrespgat, &
                  rmgat, rggat, nbpgat, litresgat, &
                  socresgat, gppgat, dstcemlsgat, litrfallgat, &
                  humiftrsgat, veghghtgat, rootdpthgat, rmlgat, &
                  litrfallveggat, humiftrsveggat, &
                  rmsgat, rmrgat, tltrleafgat, tltrstemgat, &
                  tltrrootgat, leaflitrgat, roottempgat, afrleafgat, &
                  afrstemgat, afrrootgat, wtstatusgat, ltstatusgat, &
                  burnfracgat, smfuncveggat, lucemcomgat, lucltringat, &
                  lucsocingat, nppveggat, dstcemls3gat, &
                  faregat, gavgscmsgat, &
                  rmlvegaccgat, rmsveggat, rmrveggat, rgveggat, &
                  vgbiomas_veggat, gppveggat, nepveggat, &
                  fcangat, pftexistgat, &
                  emit_co2gat, emit_cogat, emit_ch4gat, emit_nmhcgat, &
                  emit_h2gat, emit_noxgat, emit_n2ogat, emit_pm25gat, &
                  emit_tpmgat, emit_tcgat, emit_ocgat, emit_bcgat, &
                  btermgat, ltermgat, mtermgat, &
                  nbpveggat, hetroresveggat, autoresveggat, litresveggat, &
                  soilcresveggat, burnvegfgat, pstemmassgat, pgleafmassgat, &
                  ch4WetSpecgat, wetfdyngat, ch4WetDyngat, ch4soillsgat, &
                  twarmmgat, tcoldmgat, gdd5gat, &
                  ariditygat, srplsmongat, defctmongat, anndefctgat, &
                  annsrplsgat, annpcpgat, dry_season_lengthgat, &
                  anmossgat, rmlmossgat, gppmossgat, armossgat, nppmossgat, &
                  peatdepgat, litrmsmossgat, Cmossmasgat, dmossgat, &
                  ipeatlandgat, pddgat, wetfrac_presgat, &
                  tracergLeafMassgat, tracerBLeafMassgat, tracerStemMassgat, &
                  tracerRootMassgat, tracerLitrMassgat, tracerSoilCMassgat, &
                  tracerMossCMassgat, tracerMossLitrMassgat)

      if (ncount == nday) then

        DOM = DOM + 1 ! increment the day of month counter

        ! Reset mosaic accumulator arrays.
        if (ctem_on) call resetMosaicAccum

        ! Determine the active layer depth and depth to the frozen water table.
        ! This only occurs once per day since they don't change rapidly.
        call findPermafrostVars(nmtest, nltest, iday)

      end if ! ncount eq nday

      !=======================================================================

      ! Half-hourly physics outputs
      if  (doHhOutput .and. &
          (runyr >= jhhsty) .and. &
          (runyr <= jhhendy) .and. &
          (iday >= jhhstd) .and. &
          (iday <= jhhendd) ) call class_hh_w(lonLocalIndex, latLocalIndex, nltest, &
                                              nmtest, ncount, nday, iday, runyr)

      ! Daily physics outputs
      if (doDayOutput .and. &
          (runyr >= jdsty) .and. &
          (runyr <= jdendy) .and. &
          (iday  >= jdstd) .and. &
          (iday  <= jdendd))  call class_daily_aw(lonLocalIndex, latLocalIndex, iday, nltest, nmtest, &
                                                  ncount, nday, lastDOY, runyr)

      do NT = 1,NMON
        if ((IDAY == monthend(NT + 1)) .and. (NCOUNT == NDAY)) then
          IMONTH = NT
          DOM = 1 ! reset the day of month counter
        end if
      end do

      ! Monthly physics outputs
      if (doMonthOutput .and. (runyr >= jmosty)) call class_monthly_aw(lonLocalIndex, &
                                                        latLocalIndex, IDAY, runyr, NCOUNT, &
                                                        NDAY, nltest, nmtest, lastDOY)

      ! Annual physics outputs
      if (doAnnualOutput) call class_annual_aw(lonLocalIndex, latLocalIndex, IDAY, runyr, NCOUNT, NDAY, &
                                              nltest, nmtest, lastDOY)

      if (ctem_on .and. (ncount == nday)) then

        ! Convert units in preparation for output:
        call convertUnitsCTEM(nltest,nmtest)

        ! Daily outputs from biogeochem (CTEM)
        if (doDayOutput .and. &
            (runyr >= jdsty) .and. &
            (runyr <= jdendy) .and. &
            (iday   >= jdstd) .and. &
            (iday   <= jdendd)) call ctem_daily_aw(lonLocalIndex, latLocalIndex, nltest, &
                                                  nmtest, iday, ncount, nday, &
                                                  runyr, grclarea, ipeatlandrow)

        ! Monthly biogeochem outputs
        if (doMonthOutput .and. &
            (runyr >= jmosty)) call ctem_monthly_aw(lonLocalIndex, latLocalIndex, nltest, &
                                                    nmtest, iday, runyr, nday, lastDOY)

        ! Annual biogeochem outputs
        if (doAnnualOutput) call ctem_annual_aw(lonLocalIndex, latLocalIndex, iday, runyr, nltest, nmtest, lastDOY)
        
      end if

      ! Check if it is the last timestep of the last day of the year
      if ((IDAY == lastDOY) .and. (NCOUNT == NDAY)) then

        write( * , * )'IYEAR = ',IYEAR,'runyr = ',runyr,'Loop count = ',lopcount,'/',metLoop

        ! Write to the restart file
        call write_restart(lonIndex,latIndex)

        ! Increment the timestamp year (it is the number of consecutive days since the refyr. refyr is set
        ! to the first year of the run in outputManager so consecDays is 0 initially then increments up.)
        consecDays = consecDays + lastDOY

        ! Increment the runyr
        runyr = runyr + 1

        ! Reset the month
        imonth = 0

      end if ! last day of year check

      ! Increment the counter for timestep of the day, or reset it to 1.
      NCOUNT = NCOUNT + 1
      if (NCOUNT > NDAY) then
        NCOUNT = 1
      end if

      !> Now check if the met file is done, needs to loop more, or just continues to the next timestep
      if (.not. metDone) then
        !> Increment the metTimeIndex to advance to the next timestep on the next time around
        metTimeIndex = metTimeIndex + 1
      else if (metDone) then
        !> End of met array read-in reached, decide what to do
        if (lopcount == metLoop) then
          !> The lopcount is reached so the run must be over
          run_model = .false.
        else
          !> Loop again so reset the metTimeIndex
          lopcount = lopcount + 1
          metTimeIndex = 1
        end if
      end if


    end do mainModelLoop ! MAIN MODEL LOOP

    ! If we've enabled checksums, now is the time to calculate them.
    if (c_switch%doChecksums) call checksumCalc(lonIndex, latIndex)

    ! deallocate arrays used for input files
    call deallocInput

    ! This is stored in outputManager so that means it retains its values between
    ! grid cells run. It must be reset here to ensure the value doesn't carry over to
    ! the next grid cell !
    consecDays = 0.

    return

  end subroutine main_driver
  !! @}
  !> \namespace main
  !> Main model driver for CLASSIC in stand-alone mode using specified boundary
  !! conditions and atmospheric forcing.
  !!
  !! This driver program initializes the run, reads in CLASSIC input files,
  !! manages the run and the coupling between CLASS and CTEM, calls subroutines
  !! that aggregate and write outputs, and closes the run for this grid cell.

end module main
