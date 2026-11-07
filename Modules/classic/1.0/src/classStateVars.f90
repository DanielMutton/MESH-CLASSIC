!> \file
!! Contains the physics variable type structures.
!! @author J. Melton
!!
!! 1. class_rot - CLASS's 'rot' and 'row' vars
!! 2. class_gat - CLASS's 'gat' vars
!! 3. class_out - CLASS's monthly outputs

module classStateVars

  ! J. Melton Nov 2016

  use classicParams,       only : ican, icp1, NBS
  use model_variables,     only : vs

  implicit none

  public :: allocClassVars
  public :: resetClassMon
  public :: resetClassYr
  public :: resetAccVars
  public :: initDiagnosticVars
  public :: initRowVars


 ! Everything between here and class_gather was added by Daniel Mutton for MESH-CLASSIC Coupling. Retrieved from RUNCLASS36_variables
 
!     character(4) TITLE1, TITLE2, TITLE3, TITLE4, TITLE5, &
!        TITLE6, NAME1, NAME2, NAME3, NAME4, NAME5, NAME6, &
!        PLACE1, PLACE2, PLACE3, PLACE4, PLACE5, PLACE6
 
!	INTEGER ICTEMMOD
    !> PDMROF.
	!todo: isolate PDMROF from CLASS code
!     real, dimension(:), allocatable :: &
!        CMINPDM, CMAXPDM, BPDM, K1PDM, K2PDM, &
!        ZPNDPRECS, ZPONDPREC, ZPONDPREG, ZPNDPREGS, &
!        UM1CS, UM1C, UM1G, UM1GS, &
!        QM1CS, QM1C, QM1G, QM1GS, &
!        QM2CS, QM2C, QM2G, QM2GS, UMQ, &
!        FSTRCS, FSTRC, FSTRG, FSTRGS!, &
!		FREZTH, SWELIM, SNDENLIM, ICE, TICE, TSF
	
	
	!Retrieved from RUNCLASS36_variables - Daniel Mutton
!	type CLASS_forcing_input
!        real, dimension(:), allocatable :: &
!            FDL, FSIH, FSVH, PRE, PRES, QA, TA, UL, VL, VMOD
!            FCLO, ZBLD, ZRFH, ZRFM
!    end type

!    type CLASS_prognostic_variables
!        real, dimension(:), allocatable :: &
!            ALBS, CMAI, GRO, QAC, RCAN, RHOS, SNCAN, SNO, TAC, TBAS, &
!            TCAN, TPND, TSNO, WSNO, ZPND
!        real, dimension(:, :), allocatable :: &
!            TBAR, THIC, THLQ, TSFS
!    end type
	
!	integer, dimension(:, :), allocatable :: IORG, ISND   !DANIEL MUTTON ADDED FOR INITIALIZATION
	!integer, dimension(:), allocatable :: INFILTYPE		  !DANIEL MUTTON ADDED FOR INITIALIZATION
	
	!> Dimension: NML, IGND
!    integer, dimension(:), allocatable :: IGDR
!    real, dimension(:, :), allocatable :: BI, CLAY, DELZW, GRKS, HCPS, ORGM, PSIS, PSIW, SAND, TCS, THFC, THM, THP, THR, THRA, ZBTW

!    type CLASS_surface_variables

!WATROF GRKF
!WATROF WFCI
!WATROF WFSF
!WATROF XSLP

        !> Dimension: NML
!        real, dimension(:), allocatable :: &
!            AGID, AGVD, ALGD, ALGW, ASID, ASVD, DRN, FARE, GRKF, &
!            SDEP, WFCI, WFSF, XSLP, ZPLG, ZPLS, ZSNL
!        integer, dimension(:), allocatable :: MID, IWF

        !> Dimension: NML, ICAN
!        real, dimension(:, :), allocatable :: &
!            ACID, ACVD, CMAS, HGTD, PAID, PAMN, PAMX, PSGA, PSGB, &
!            QA50, ROOT, RSMN, VPDA, VPDB

        !> Dimension: NML, ICP1
!        real, dimension(:, :), allocatable :: &
!            ALIC, ALVC, FCAN, LNZ0
!    end type

!    type CLASS_glacier_variables
!        real, dimension(:), allocatable :: FREZTH, SWELIM, SNDENLIM
!    end type

!    type CLASS_atmospheric_variables
!        real, dimension(:), allocatable :: &
!            CSZ, DLON, FCLO, GC, GGEO, PADR, RADJ, RHOA, RHSI, RPCP, &
!            RPRE, SPCP, SPRE, TADP, TRPC, TSPC, VPD, Z0OR, ZBLD, ZDH, &
!            ZDM, ZRFH, ZRFM
!            FDL, FSIH, FSVH, PRE, PRES, QA, TA, UL, UV, VL
!    end type

!    type CLASS_diagnostic_variables
!        integer ISUM
!        integer, dimension(:, :, :), allocatable :: ITCT
!        real, dimension(:), allocatable :: &
!            ALIR, ALVS, CDH, CDM, DR, EF, FCS, FGS, FC, FG, FLGG, &
!            FLGS, FLGV, FSGG, FSGS, FSGV, FSNO, GA, GTE, HBL, HEVC, &
!            HEVG, HEVS, HFS, HFSC, HFSG, HFSS, HMFC, HMFN, HTCC, HTCS, &
!            ILMO, PCFC, PCLC, PCPG, PCPN, PET, QEVP, QFCF, QFCL, QFG, &
!            QFN, QFS, QFX, QG, ROF, ROFB, ROFC, ROFN, ROFO, &
!            ROVG, SFCQ, SFCT, SFCU, SFCV, TFX, TROB, TROF, TROO, &
!            TSF, UE, WTAB, WTRC, WTRG, WTRS, ICE, TICE
!        real, dimension(:, :), allocatable :: &
!            ROFS, TROS, GFLX, HMFG, HTC, QFC
!    end type
	
	
!	integer, dimension(:), allocatable :: MID, IWF                   ! I DONT KNOW WHAT THESE ARE FOR OR HOW THEY'RE DIFFERENT FROM THE INT IWF FLAG OR MIDROT - DANIEL MUTTON
	
	
!    type(CLASS_forcing_input), save :: cfi
!    type(CLASS_prognostic_variables), save :: cpv
!    type(CLASS_surface_variables), save :: csfv
!    type(CLASS_glacier_variables), save :: cglv
!    type(CLASS_atmospheric_variables), save :: catv
!    type(CLASS_diagnostic_variables), save :: cdv

  !=================================================================================
  !> Physics variables in the 'gather' structure
  type class_gather

    ! These will be allocated the dimension: 'ilg'

    integer, allocatable, dimension(:) :: ILMOS     !< Index of grid cell corresponding to current element of gathered vector of land surface variables [ ]
    integer, allocatable, dimension(:) :: JLMOS     !< Index of mosaic tile corresponding to current element of gathered vector of land surface variables [ ]
    integer, allocatable, dimension(:) :: IWMOS     !< Index of grid cell corresponding to current element of gathered vector of inland water body variables [ ]
    integer, allocatable, dimension(:) :: JWMOS     !< Index of mosaic tile corresponding to current element of gathered vector of inland water body variables [ ]
    integer, allocatable, dimension(:) :: IGDRGAT   !< Index of soil layer in which bedrock is encountered

	integer, allocatable, dimension(:) :: TEST    !< Test to ensure bridge exists
    real, allocatable, dimension(:) :: GCGAT   !< Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
    real, allocatable, dimension(:) :: TZSGAT  !< Vertical temperature gradient in a snow pack
    real, allocatable, dimension(:) :: PCSNGAT !< Snow fall flux \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ALBSGAT !< Snow albedo [ ]
    real, allocatable, dimension(:) :: CMAIGAT !< Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: GROGAT  !< Vegetation growth index [ ]
    real, allocatable, dimension(:) :: QACGAT  !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: RCANGAT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: RHOSGAT !< Density of snow \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: SCANGAT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: SNOGAT  !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: TACGAT  !< Temperature of air within vegetation canopy [K]
    real, allocatable, dimension(:) :: TBASGAT !< Temperature of bedrock in third soil layer [K]
    real, allocatable, dimension(:) :: TCANGAT !< Vegetation canopy temperature [K]
    real, allocatable, dimension(:) :: TPNDGAT !< Temperature of ponded water [K]
    real, allocatable, dimension(:) :: TSNOGAT !< Snowpack temperature [K]
    real, allocatable, dimension(:) :: WSNOGAT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: maxAnnualActLyrGAT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, allocatable, dimension(:) :: ZPNDGAT !< Depth of ponded water on surface [m]
    real, allocatable, dimension(:) :: REFGAT  !<
    real, allocatable, dimension(:) :: BCSNGAT !<
    real, allocatable, dimension(:) :: AGIDGAT !< Optional user-specified value of ground near-infrared albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: AGVDGAT !< Optional user-specified value of ground visible albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: ALGDGAT !< Reference albedo for dry soil [ ]
    real, allocatable, dimension(:) :: ALGWGAT !< Reference albedo for saturated soil [ ]
    real, allocatable, dimension(:) :: ASIDGAT !< Optional user-specified value of snow near-infrared albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: ASVDGAT !< Optional user-specified value of snow visible albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:) :: DRNGAT  !< Drainage index at bottom of soil profile [ ]
    real, allocatable, dimension(:) :: GRKFGAT !< WATROF parameter used when running MESH code [ ]
    real, allocatable, dimension(:) :: WFCIGAT !< WATROF parameter used when running MESH code [ ]
    real, allocatable, dimension(:) :: WFSFGAT !< WATROF parameter used when running MESH code [ ]
    real, allocatable, dimension(:) :: XSLPGAT !< Surface slope (used when running MESH code) [degrees]
    real, allocatable, dimension(:) :: ZPLGGAT !< Maximum water ponding depth for snow-free subareas (user-specified when running MESH code) [m]
    real, allocatable, dimension(:) :: ZPLSGAT !< Maximum water ponding depth for snow-covered subareas (user-specified when running MESH code) [m]
    real, allocatable, dimension(:) :: ZSNLGAT !< Limiting snow depth below which coverage is < 100% [m]
    real, allocatable, dimension(:) :: ALGWVGAT !<
    real, allocatable, dimension(:) :: ALGWNGAT !<
    real, allocatable, dimension(:) :: ALGDVGAT !<
    real, allocatable, dimension(:) :: ALGDNGAT !<
    real, allocatable, dimension(:) :: EMISGAT  !<
    real, allocatable, dimension(:) :: CSZGAT  !< Cosine of solar zenith angle [ ]
    real, allocatable, dimension(:) :: DLONGAT !< Longitude of grid cell (east of Greenwich) [degrees]
    real, allocatable, dimension(:) :: DLATGAT !< Latitude of grid cell [degrees]
    real, allocatable, dimension(:) :: FCLOGAT !< Fractional cloud cover [ ]
    real, allocatable, dimension(:) :: FDLGAT  !< Downwelling longwave radiation at bottom of atmosphere (i.e. incident on modelled land surface elements \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSSGAT  !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSIHGAT !< Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSVHGAT !< Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GGEOGAT !< Geothermal heat flux at bottom of soil profile \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: PADRGAT !< Partial pressure of dry air [Pa]
    real, allocatable, dimension(:) :: PREGAT  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PRESGAT !< Surface air pressure [Pa]
    real, allocatable, dimension(:) :: QAGAT   !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: RADJGAT !< Latitude of grid cell (positive north of equator) [rad]
    real, allocatable, dimension(:) :: RHOAGAT !< Density of air \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: RHSIGAT !< Density of fresh snow \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:) :: RPCPGAT !< Rainfall rate over modelled area \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: SPCPGAT !< Snowfall rate over modelled area \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: TAGAT   !< Air temperature at reference height [K]
    real, allocatable, dimension(:) :: TADPGAT !< Dew point temperature of air [K]
    real, allocatable, dimension(:) :: TRPCGAT !< Rainfall temperature [K]
    real, allocatable, dimension(:) :: TSPCGAT !< Snowfall temperature [K]
    real, allocatable, dimension(:) :: ULGAT   !< Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VLGAT   !< Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VMODGAT !< Wind speed at reference height \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VPDGAT  !< Vapour pressure deficit [mb]
    real, allocatable, dimension(:) :: Z0ORGAT !< Orographic roughness length [m]
    real, allocatable, dimension(:) :: ZBLDGAT !< Atmospheric blending height for surface roughness length averaging [m]
    real, allocatable, dimension(:) :: ZDHGAT  !< User-specified height associated with diagnosed screen-level variables [m]
    real, allocatable, dimension(:) :: ZDMGAT  !< User-specified height associated with diagnosed anemometer-level wind speed [m]
    real, allocatable, dimension(:) :: ZRFHGAT !< Reference height associated with forcing air temperature and humidity [m]
    real, allocatable, dimension(:) :: ZRFMGAT !< Reference height associated with forcing wind speed [m]
    real, allocatable, dimension(:) :: FSGGAT  !<
    real, allocatable, dimension(:) :: FLGGAT  !<
    real, allocatable, dimension(:) :: GUSTGAT !<
    real, allocatable, dimension(:) :: DEPBGAT !<
    real, allocatable, dimension(:) :: GTBS    !<
    real, allocatable, dimension(:) :: SFCUBS  !<
    real, allocatable, dimension(:) :: SFCVBS  !<
    real, allocatable, dimension(:) :: USTARBS !<
    real, allocatable, dimension(:) :: TCSNOW  !<
    real, allocatable, dimension(:) :: GSNOW   !<
    real, allocatable, dimension(:) :: ALIRGAT !< Diagnosed total near-infrared albedo of land surface [ ]
    real, allocatable, dimension(:) :: ALVSGAT !< Diagnosed total visible albedo of land surface [ ]
    real, allocatable, dimension(:) :: CDHGAT  !< Surface drag coefficient for heat [ ]
    real, allocatable, dimension(:) :: CDMGAT  !< Surface drag coefficient for momentum [ ]
    real, allocatable, dimension(:) :: DRGAT   !< Surface drag coefficient under neutral stability [ ]
    real, allocatable, dimension(:) :: EFGAT   !< Evaporation efficiency at ground surface [ ]
    real, allocatable, dimension(:) :: FLGGGAT !< Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLGSGAT !< Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLGVGAT !< Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSGGGAT !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSGSGAT !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSGVGAT !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSNOGAT !< Diagnosed fractional snow coverage [ ]
    real, allocatable, dimension(:) :: GAGAT   !< Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: GTGAT   !< Diagnosed effective surface black-body temperature [K]
    real, allocatable, dimension(:) :: HBLGAT  !< Height of the atmospheric boundary layer [m]
    real, allocatable, dimension(:) :: HEVCGAT !< Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HEVGGAT !< Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HEVSGAT !< Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSGAT  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSCGAT !< Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSGGAT !< Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSSGAT !< Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HMFCGAT !< Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HMFNGAT !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HTCCGAT !< Diagnosed internal energy change of vegetation canopy due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HTCSGAT !< Diagnosed internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: ILMOGAT !< Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
    real, allocatable, dimension(:) :: PCFCGAT !< Diagnosed frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PCLCGAT !< Diagnosed liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PCPGGAT !< Diagnosed precipitation incident on ground \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PCPNGAT !< Diagnosed precipitation incident on snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PETGAT  !< Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QEVPGAT !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: QFCFGAT !< Diagnosed vapour flux from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFCLGAT !< Diagnosed vapour flux from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFGGAT  !< Diagnosed water vapour flux from ground \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFNGAT  !< Diagnosed water vapour flux from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFSGAT  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: QFXGAT  !< Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: QGGAT   !< Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFGAT  !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFBGAT !< Base flow from bottom of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFCGAT !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFNGAT !< Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROFOGAT !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:,:) :: ROFSGAT !< Interflow from sides of soil column \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: ROVGGAT !< Diagnosed liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: SFCQGAT !< Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:) :: SFCTGAT !< Diagnosed screen-level air temperature [K]
    real, allocatable, dimension(:) :: SFCUGAT !< Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: SFCVGAT !< Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: TFXGAT  !< Product of surface drag coefficient, wind speed and surface-air temperature difference \f$[K m s^{-1} ]\f$
    real, allocatable, dimension(:) :: TROBGAT !< Temperature of base flow from bottom of soil column [K]
    real, allocatable, dimension(:) :: TROFGAT !< Temperature of total runoff [K]
    real, allocatable, dimension(:) :: TROOGAT !< Temperature of overland flow from top of soil column [K]
    real, allocatable, dimension(:,:) :: TROSGAT !< Temperature of interflow from sides of soil column [K]
    real, allocatable, dimension(:) :: UEGAT   !< Friction velocity of air \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: WTABGAT !< Depth of water table in soil [m]
    real, allocatable, dimension(:) :: WTRCGAT !< Diagnosed residual water transferred off the vegetation canopy \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: WTRGGAT !< Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: WTRSGAT !< Diagnosed residual water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: wtableGAT !< Depth of water table in soil [m]
    real, allocatable, dimension(:) :: FTABLE !< Depth to frozen water table (m)
    real, allocatable, dimension(:) :: ACTLYR !< Active layer depth (m)
    real, allocatable, dimension(:) :: actLyrThisYrGAT !< Annual active layer depth maximum starting from summer solstice for the present year (m)
    real, allocatable, dimension(:) :: QLWOGAT !<
    real, allocatable, dimension(:) :: SFRHGAT !<
    real, allocatable, dimension(:) :: FTEMP   !<
    real, allocatable, dimension(:) :: FVAP    !<
    real, allocatable, dimension(:) :: RIB     !<
    real, allocatable, dimension(:) :: FC      !< Subarea fractional coverage of modelled area - ground under canopy [ ]
    real, allocatable, dimension(:) :: FG      !< Subarea fractional coverage of modelled area - bare ground [ ]
    real, allocatable, dimension(:) :: FCS     !< Subarea fractional coverage of modelled area - snow-covered ground under canopy  [ ]
    real, allocatable, dimension(:) :: FGS     !< Subarea fractional coverage of modelled area - snow-covered bare ground [ ]
    real, allocatable, dimension(:) :: RBCOEF  !<
    real, allocatable, dimension(:) :: ZSNOW   !<
    real, allocatable, dimension(:) :: FSVF    !<
    real, allocatable, dimension(:) :: FSVFS   !<
    real, allocatable, dimension(:) :: ALVSCN  !<
    real, allocatable, dimension(:) :: ALIRCN  !<
    real, allocatable, dimension(:) :: ALVSG   !<
    real, allocatable, dimension(:) :: ALIRG   !<
    real, allocatable, dimension(:) :: ALVSCS  !<
    real, allocatable, dimension(:) :: ALIRCS  !<
    real, allocatable, dimension(:) :: ALVSSN  !<
    real, allocatable, dimension(:) :: ALIRSN  !<
    real, allocatable, dimension(:) :: ALVSGC  !<
    real, allocatable, dimension(:) :: ALIRGC  !<
    real, allocatable, dimension(:) :: ALVSSC  !<
    real, allocatable, dimension(:) :: ALIRSC  !<
    real, allocatable, dimension(:) :: TRVSCN  !<
    real, allocatable, dimension(:) :: TRIRCN  !<
    real, allocatable, dimension(:) :: TRVSCS  !<
    real, allocatable, dimension(:) :: TRIRCS  !<
    real, allocatable, dimension(:) :: RC      !<
    real, allocatable, dimension(:) :: RCS     !<
    real, allocatable, dimension(:) :: FRAINC  !<
    real, allocatable, dimension(:) :: FSNOWC  !<
    real, allocatable, dimension(:) :: FRAICS  !<
    real, allocatable, dimension(:) :: FSNOCS  !<
    real, allocatable, dimension(:) :: CMASSC  !<
    real, allocatable, dimension(:) :: CMASCS  !<
    real, allocatable, dimension(:) :: DISP    !<
    real, allocatable, dimension(:) :: DISPS   !<
    real, allocatable, dimension(:) :: ZOMLNC  !<
    real, allocatable, dimension(:) :: ZOELNC  !<
    real, allocatable, dimension(:) :: ZOMLNG  !<
    real, allocatable, dimension(:) :: ZOELNG  !<
    real, allocatable, dimension(:) :: ZOMLCS  !<
    real, allocatable, dimension(:) :: ZOELCS  !<
    real, allocatable, dimension(:) :: ZOMLNS  !<
    real, allocatable, dimension(:) :: ZOELNS  !<
    real, allocatable, dimension(:) :: TRSNOWC !<
    real, allocatable, dimension(:) :: CHCAP   !<
    real, allocatable, dimension(:) :: CHCAPS  !<
    real, allocatable, dimension(:) :: GZEROC  !< Vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GZEROG  !< Bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GZROCS  !< Snow-covered vegetated subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GZROGS  !<  Snow-covered bare ground subarea heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: groundHeatFlux  !<  !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: G12C    !<
    real, allocatable, dimension(:) :: G12G    !<
    real, allocatable, dimension(:) :: G12CS   !<
    real, allocatable, dimension(:) :: G12GS   !<
    real, allocatable, dimension(:) :: G23C    !<
    real, allocatable, dimension(:) :: G23G    !<
    real, allocatable, dimension(:) :: G23CS   !<
    real, allocatable, dimension(:) :: G23GS   !<
    real, allocatable, dimension(:) :: QFREZC  !<
    real, allocatable, dimension(:) :: QFREZG  !<
    real, allocatable, dimension(:) :: QMELTC  !<
    real, allocatable, dimension(:) :: QMELTG  !<
    real, allocatable, dimension(:) :: EVAPC   !<
    real, allocatable, dimension(:) :: EVAPCG  !<
    real, allocatable, dimension(:) :: EVAPG   !<
    real, allocatable, dimension(:) :: EVAPCS  !<
    real, allocatable, dimension(:) :: EVPCSG  !<
    real, allocatable, dimension(:) :: EVAPGS  !<
    real, allocatable, dimension(:) :: TCANO   !< Temperature of canopy over ground [K]
    real, allocatable, dimension(:) :: TCANS   !< Temperature of canopy over snow [K]
    real, allocatable, dimension(:) :: RAICAN  !<
    real, allocatable, dimension(:) :: SNOCAN  !<
    real, allocatable, dimension(:) :: RAICNS  !<
    real, allocatable, dimension(:) :: SNOCNS  !<
    real, allocatable, dimension(:) :: CWLCAP  !<
    real, allocatable, dimension(:) :: CWFCAP  !<
    real, allocatable, dimension(:) :: CWLCPS  !<
    real, allocatable, dimension(:) :: CWFCPS  !<
    real, allocatable, dimension(:) :: TSNOCS  !<
    real, allocatable, dimension(:) :: TSNOGS  !<
    real, allocatable, dimension(:) :: RHOSCS  !<
    real, allocatable, dimension(:) :: RHOSGS  !<
    real, allocatable, dimension(:) :: WSNOCS  !<
    real, allocatable, dimension(:) :: WSNOGS  !<
    real, allocatable, dimension(:) :: TPONDC  !<
    real, allocatable, dimension(:) :: TPONDG  !<
    real, allocatable, dimension(:) :: TPNDCS  !<
    real, allocatable, dimension(:) :: TPNDGS  !<
    real, allocatable, dimension(:) :: ZPLMCS  !<
    real, allocatable, dimension(:) :: ZPLMGS  !<
    real, allocatable, dimension(:) :: ZPLIMC  !<
    real, allocatable, dimension(:) :: ZPLIMG  !<
    !
    !     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
    !     * BALANCES.
    !
    real, allocatable, dimension(:) :: CTVSTP !<
    real, allocatable, dimension(:) :: CTSSTP !<
    real, allocatable, dimension(:) :: CT1STP !<
    real, allocatable, dimension(:) :: CT2STP !<
    real, allocatable, dimension(:) :: CT3STP !<
    real, allocatable, dimension(:) :: WTVSTP !<
    real, allocatable, dimension(:) :: WTSSTP !<
    real, allocatable, dimension(:) :: WTGSTP !<

    ! These will be allocated the dimension: 'ignd'
    real, allocatable, dimension(:) :: DELZ    !< Overall thickness of soil layer [m]
    real, allocatable, dimension(:) :: ZBOT    !< Depth of to the bottom of soil layer [m]

    ! These will be allocated the dimension: 'ilg,ignd'
    integer, allocatable, dimension(:,:) :: ISNDGAT !< Integer identifier associated with sand content
    real, allocatable, dimension(:,:) :: TBARGAT !< Temperature of soil layers [K]
    real, allocatable, dimension(:,:) :: THICGAT !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THLQGAT !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: BIGAT   !< Clapp and Hornberger empirical “b” parameter [ ]
    real, allocatable, dimension(:,:) :: DLZWGAT !< Permeable thickness of soil layer [m]
    real, allocatable, dimension(:,:) :: GRKSGAT !< Saturated hydraulic conductivity of soil layers \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:,:) :: HCPSGAT !< Volumetric heat capacity of soil particles \f$[J m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: PSISGAT !< Soil moisture suction at saturation [m]
    real, allocatable, dimension(:,:) :: PSIWGAT !< Soil moisture suction at wilting point [m]
    real, allocatable, dimension(:,:) :: TCSGAT  !< Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$\
    real, allocatable, dimension(:,:) :: THFCGAT !< Field capacity \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THMGAT  !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THPGAT  !< Pore volume in soil layer \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THRGAT  !< Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THRAGAT !< Fractional saturation of soil behind the wetting front [ ]
    real, allocatable, dimension(:,:) :: ZBTWGAT !< Depth to permeable bottom of soil layer [m]
    real, allocatable, dimension(:,:) :: THLWGAT !< Soil water content at wilting point, \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: GFLXGAT !< Heat conduction between soil layers \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HMFGGAT !< Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HTCGAT  !< Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: QFCGAT  !< Diagnosed vapour flux from transpiration over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: TBARC  !<
    real, allocatable, dimension(:,:) :: TBARG  !<
    real, allocatable, dimension(:,:) :: TBARCS !<
    real, allocatable, dimension(:,:) :: TBARGS !<
    real, allocatable, dimension(:,:) :: THLIQC !<
    real, allocatable, dimension(:,:) :: THLIQG !<
    real, allocatable, dimension(:,:) :: THICEC !<
    real, allocatable, dimension(:,:) :: THICEG !<
    real, allocatable, dimension(:,:) :: FROOT  !<
    real, allocatable, dimension(:,:) :: HCPC   !<
    real, allocatable, dimension(:,:) :: HCPG   !<
    real, allocatable, dimension(:,:) :: FROOTS !<
    real, allocatable, dimension(:,:) :: TCTOPC !<
    real, allocatable, dimension(:,:) :: TCBOTC !<
    real, allocatable, dimension(:,:) :: TCTOPG !<
    real, allocatable, dimension(:,:) :: TCBOTG !<

    ! These will be allocated the dimension: 'ilg,ican'
    real, allocatable, dimension(:,:) :: ACIDGAT !< Optional user-specified value of canopy near-infrared albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:,:) :: ACVDGAT !< Optional user-specified value of canopy visible albedo to override CLASS-calculated value [ ]
    real, allocatable, dimension(:,:) :: CMASGAT !< Maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HGTDGAT !< Optional user-specified values of height of vegetation categories to override CLASS-calculated values [m]
    real, allocatable, dimension(:,:) :: PAIDGAT !< Optional user-specified value of plant area indices of vegetation categories to override CLASS-calculated values [ ]
    real, allocatable, dimension(:,:) :: PAMNGAT !< Minimum plant area index of vegetation category [ ]
    real, allocatable, dimension(:,:) :: PAMXGAT !< Minimum plant area index of vegetation category [ ]
    real, allocatable, dimension(:,:) :: PSGAGAT !< Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, allocatable, dimension(:,:) :: PSGBGAT !< Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, allocatable, dimension(:,:) :: QA50GAT !< Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: ROOTGAT !< Maximum rooting depth of vegetation category [m]
    real, allocatable, dimension(:,:) :: RSMNGAT !< Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$
    real, allocatable, dimension(:,:) :: VPDAGAT !< Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]
    real, allocatable, dimension(:,:) :: VPDBGAT !< Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) [ ]


    ! These will be allocated the dimension: 'ilg,icp1'
    real, allocatable, dimension(:,:) :: ALICGAT !< Background average near-infrared albedo of vegetation category [ ]
    real, allocatable, dimension(:,:) :: ALVCGAT !< Background average visible albedo of vegetation category [ ]
    real, allocatable, dimension(:,:) :: FCANGAT !< Maximum fractional coverage of modelled area by vegetation category [ ]
    real, allocatable, dimension(:,:) :: LNZ0GAT !< Natural logarithm of maximum roughness length of vegetation category [ ]

    ! These will be allocated the dimension: 'ilg,nbs'
    real, allocatable, dimension(:,:) :: FSDBGAT !<
    real, allocatable, dimension(:,:) :: FSFBGAT !<
    real, allocatable, dimension(:,:) :: FSSBGAT !<
    real, allocatable, dimension(:,:) :: SALBGAT !<
    real, allocatable, dimension(:,:) :: CSALGAT !<
    real, allocatable, dimension(:,:) :: ALTG    !<
    real, allocatable, dimension(:,:) :: ALSNO   !<
    real, allocatable, dimension(:,:) :: TRSNOWG !<

    ! These will be allocated the dimension: 'ilg,4'
    real, allocatable, dimension(:,:) :: TSFSGAT !< Ground surface temperature over subarea [K]

    ! These will be allocated the dimension: 'ilg,6,50'
    integer, allocatable, dimension(:,:,:) :: ITCTGAT !< Counter of number of iterations required to solve surface energy balance for the elements of the four subareas

  end type class_gather

  type (class_gather), save, target :: class_gat

  ! ================================================================================
  !> Physics variables in the 'rotated' (rot) structure
  type class_rotated

    ! These will be allocated the dimension: 'nlat'

    real, allocatable, dimension(:) :: WTBLACC !< Depth of water table in soil [m]
    real, allocatable, dimension(:) :: CSZROW  !<
    real, allocatable, dimension(:) :: DLONROW !<
    real, allocatable, dimension(:) :: DLATROW !<
    real, allocatable, dimension(:) :: FCLOROW !< Fractional cloud cover [ ]
    real, allocatable, dimension(:) :: RHOSROW !< Density of snow \f$[kg m^{-3}]\f$
    real, allocatable, dimension(:) :: FDLROW  !< Downwelling longwave sky radiation \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSIHROW !< Near infrared shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSVHROW !< Visible shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: GCROW   !< Type identifier for grid cell (1 = sea ice, 0 = ocean, -1 = land)
    real, allocatable, dimension(:) :: GGEOROW !<
    real, allocatable, dimension(:) :: PADRROW !<
    real, allocatable, dimension(:) :: PREROW  !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PRESROW !< Surface air pressure \f$[P_a]\f$
    real, allocatable, dimension(:) :: QAROW   !< Specific humidity at reference height \f$[kg kg^{-1}]\f$
    real, allocatable, dimension(:) :: RADJROW !<
    real, allocatable, dimension(:) :: RHOAROW !<
    real, allocatable, dimension(:) :: RHSIROW !<
    real, allocatable, dimension(:) :: RPCPROW !<
    real, allocatable, dimension(:) :: RPREROW !< Rainfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: SPCPROW !<
    real, allocatable, dimension(:) :: SPREROW !< Snowfall rate over modelled area \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: TAROW   !< Air temperature at reference height [K]
    real, allocatable, dimension(:) :: TSNOROW !< Snowpack temperature [K]
    real, allocatable, dimension(:) :: TCANROW !< Vegetation canopy temperature [K]
    real, allocatable, dimension(:) :: TPNDROW !< Temperature of ponded water [K]
    real, allocatable, dimension(:) :: ZPNDROW !< Depth of ponded water [m]
    real, allocatable, dimension(:) :: SCANROW !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: RCANROW !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: TADPROW !<
    real, allocatable, dimension(:) :: TRPCROW !<
    real, allocatable, dimension(:) :: TSPCROW !<
    real, allocatable, dimension(:) :: ULROW   !< Zonal component of wind velocity \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VLROW   !< Meridional component of wind velocity \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VMODROW !< Wind speed at reference height \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: VPDROW  !<
    real, allocatable, dimension(:) :: ZBLDROW !< Atmospheric blending height for surface roughness length averaging [m]
    real, allocatable, dimension(:) :: ZDHROW  !<
    real, allocatable, dimension(:) :: ZDMROW  !<
    real, allocatable, dimension(:) :: ZRFHROW !< Reference height associated with forcing air temperature and humidity [m]
    real, allocatable, dimension(:) :: ZRFMROW !< Reference height associated with forcing wind speed [m]
    real, allocatable, dimension(:) :: UVROW   !< Wind speed at reference height \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:) :: XDIFFUS !<
    real, allocatable, dimension(:) :: Z0ORROW !<
    real, allocatable, dimension(:) :: FSSROW  !< Shortwave radiation \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: PRENROW !<
    real, allocatable, dimension(:) :: CLDTROW !<
    real, allocatable, dimension(:) :: FSGROL  !<
    real, allocatable, dimension(:) :: FLGROL  !<
    real, allocatable, dimension(:) :: GUSTROL !<
    real, allocatable, dimension(:) :: DEPBROW !<
    real, allocatable, dimension(:) :: ALIRROW !<
    real, allocatable, dimension(:) :: ALVSROW !<
    real, allocatable, dimension(:) :: CDHROW  !<
    real, allocatable, dimension(:) :: CDMROW  !<
    real, allocatable, dimension(:) :: DRROW   !<
    real, allocatable, dimension(:) :: EFROW   !<
    real, allocatable, dimension(:) :: FLGGROW !<
    real, allocatable, dimension(:) :: FLGSROW !<
    real, allocatable, dimension(:) :: FLGVROW !<
    real, allocatable, dimension(:) :: FSGGROW !<
    real, allocatable, dimension(:) :: FSGSROW !<
    real, allocatable, dimension(:) :: FSGVROW !<
    real, allocatable, dimension(:) :: FSNOROW !<
    real, allocatable, dimension(:) :: GAROW   !<
    real, allocatable, dimension(:) :: GTROW   !<
    real, allocatable, dimension(:) :: HBLROW  !<
    real, allocatable, dimension(:) :: HEVCROW !<
    real, allocatable, dimension(:) :: HEVGROW !<
    real, allocatable, dimension(:) :: HEVSROW !<
    real, allocatable, dimension(:) :: HFSROW  !<
    real, allocatable, dimension(:) :: HFSCROW !<
    real, allocatable, dimension(:) :: HFSGROW !<
    real, allocatable, dimension(:) :: HFSSROW !<
    real, allocatable, dimension(:) :: HMFCROW !<
    real, allocatable, dimension(:) :: HMFNROW !<
    real, allocatable, dimension(:) :: HTCCROW !<
    real, allocatable, dimension(:) :: HTCSROW !<
    real, allocatable, dimension(:) :: ILMOROW !<
    real, allocatable, dimension(:) :: PCFCROW !<
    real, allocatable, dimension(:) :: PCLCROW !<
    real, allocatable, dimension(:) :: PCPGROW !<
    real, allocatable, dimension(:) :: PCPNROW !<
    real, allocatable, dimension(:) :: PETROW  !<
    real, allocatable, dimension(:) :: QEVPROW !<
    real, allocatable, dimension(:) :: QFCFROW !<
    real, allocatable, dimension(:) :: QFCLROW !<
    real, allocatable, dimension(:) :: QFGROW  !<
    real, allocatable, dimension(:) :: QFNROW  !<
    real, allocatable, dimension(:) :: QFSROW  !<
    real, allocatable, dimension(:) :: QFXROW  !<
    real, allocatable, dimension(:) :: QGROW   !<
    real, allocatable, dimension(:) :: ROFROW  !<
    real, allocatable, dimension(:) :: ROFBROW !<
    real, allocatable, dimension(:) :: ROFCROW !<
    real, allocatable, dimension(:) :: ROFNROW !<
    real, allocatable, dimension(:) :: ROFOROW !<
    real, allocatable, dimension(:) :: ROFSROW !<
    real, allocatable, dimension(:) :: ROVGROW !<
    real, allocatable, dimension(:) :: SFCQROW !<
    real, allocatable, dimension(:) :: SFCTROW !<
    real, allocatable, dimension(:) :: SFCUROW !<
    real, allocatable, dimension(:) :: SFCVROW !<
    real, allocatable, dimension(:) :: TFXROW  !<
    real, allocatable, dimension(:) :: UEROW   !<
    real, allocatable, dimension(:) :: WTABROW !<
    real, allocatable, dimension(:) :: WTRCROW !<
    real, allocatable, dimension(:) :: WTRGROW !<
    real, allocatable, dimension(:) :: WTRSROW !<
    real, allocatable, dimension(:) :: SFRHROW !<
    real, allocatable, dimension(:) :: SNOROW  !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, allocatable, dimension(:) :: WSNOROW !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: wtableROW !<
    integer, allocatable, dimension(:) :: altotcntr_d     !< Used to count the number of time steps with the sun above the horizon

    ! These will be allocated the dimension: 'nlat,nmos'

    integer, allocatable, dimension(:,:) :: IGDRROT !<
    integer, allocatable, dimension(:,:) :: MIDROT  !< Mosaic tile type identifier (1 for land surface, 0 for inland lake)
    real, allocatable, dimension(:,:) :: ALBSROT !< Snow albedo [ ]
    real, allocatable, dimension(:,:) :: CMAIROT !<
    real, allocatable, dimension(:,:) :: GROROT  !< Vegetation growth index [ ]
    real, allocatable, dimension(:,:) :: QACROT  !<
    real, allocatable, dimension(:,:) :: RCANROT !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: RHOSROT !< Density of snow \f$[kg m^{-3}]\f$
    real, allocatable, dimension(:,:) :: SCANROT !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: SNOROT  !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, allocatable, dimension(:,:) :: TACROT  !<
    real, allocatable, dimension(:,:) :: TBASROT !<
    real, allocatable, dimension(:,:) :: TCANROT !< Vegetation canopy temperature [K]
    real, allocatable, dimension(:,:) :: TPNDROT !< Temperature of ponded water [K]
    real, allocatable, dimension(:,:) :: TSNOROT !< Snowpack temperature [K]
    real, allocatable, dimension(:,:) :: WSNOROT !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: ZPNDROT !< Depth of ponded water [m]
    real, allocatable, dimension(:,:) :: REFROT  !<
    real, allocatable, dimension(:,:) :: BCSNROT !<
    real, allocatable, dimension(:,:) :: AGIDROT !<
    real, allocatable, dimension(:,:) :: AGVDROT !<
    real, allocatable, dimension(:,:) :: ALGDROT !<
    real, allocatable, dimension(:,:) :: ALGWROT !<
    real, allocatable, dimension(:,:) :: ASIDROT !<
    real, allocatable, dimension(:,:) :: ASVDROT !<
    real, allocatable, dimension(:,:) :: DRNROT  !<
    real, allocatable, dimension(:,:) :: FAREROT !< Fractional coverage of mosaic tile on modelled area
    real, allocatable, dimension(:,:) :: GRKFROT !<
    real, allocatable, dimension(:,:) :: WFCIROT !<
    real, allocatable, dimension(:,:) :: WFSFROT !<
    real, allocatable, dimension(:,:) :: XSLPROT !<
    real, allocatable, dimension(:,:) :: ZPLGROT !<
    real, allocatable, dimension(:,:) :: ZPLSROT !<
    real, allocatable, dimension(:,:) :: ZSNLROT !<
    real, allocatable, dimension(:,:) :: ZSNOROT  !<
    real, allocatable, dimension(:,:) :: ALGWVROT !<
    real, allocatable, dimension(:,:) :: ALGWNROT !<
    real, allocatable, dimension(:,:) :: ALGDVROT !<
    real, allocatable, dimension(:,:) :: ALGDNROT !<
    real, allocatable, dimension(:,:) :: EMISROT  !<
    real, allocatable, dimension(:,:) :: ALIRROT !<
    real, allocatable, dimension(:,:) :: ALVSROT !<
    real, allocatable, dimension(:,:) :: CDHROT  !<
    real, allocatable, dimension(:,:) :: CDMROT  !<
    real, allocatable, dimension(:,:) :: DRROT   !<
    real, allocatable, dimension(:,:) :: EFROT   !<
    real, allocatable, dimension(:,:) :: FLGGROT !<
    real, allocatable, dimension(:,:) :: FLGSROT !<
    real, allocatable, dimension(:,:) :: FLGVROT !<
    real, allocatable, dimension(:,:) :: FSGGROT !<
    real, allocatable, dimension(:,:) :: FSGSROT !<
    real, allocatable, dimension(:,:) :: FSGVROT !<
    real, allocatable, dimension(:,:) :: FSNOROT !<
    real, allocatable, dimension(:,:) :: GAROT   !<
    real, allocatable, dimension(:,:) :: GTROT   !<
    real, allocatable, dimension(:,:) :: HBLROT  !<
    real, allocatable, dimension(:,:) :: HEVCROT !<
    real, allocatable, dimension(:,:) :: HEVGROT !<
    real, allocatable, dimension(:,:) :: HEVSROT !<
    real, allocatable, dimension(:,:) :: HFSROT  !<
    real, allocatable, dimension(:,:) :: HFSCROT !<
    real, allocatable, dimension(:,:) :: HFSGROT !<
    real, allocatable, dimension(:,:) :: HFSSROT !<
    real, allocatable, dimension(:,:) :: HMFCROT !<
    real, allocatable, dimension(:,:) :: HMFNROT !<
    real, allocatable, dimension(:,:) :: HTCCROT !<
    real, allocatable, dimension(:,:) :: SDEPROT !< Depth to bedrock in the soil profile
    real, allocatable, dimension(:,:) :: SOCIROT !<
    real, allocatable, dimension(:,:) :: HTCSROT !<
    real, allocatable, dimension(:,:) :: ILMOROT !<
    real, allocatable, dimension(:,:) :: PCFCROT !<
    real, allocatable, dimension(:,:) :: PCLCROT !<
    real, allocatable, dimension(:,:) :: PCPGROT !<
    real, allocatable, dimension(:,:) :: PCPNROT !<
    real, allocatable, dimension(:,:) :: PETROT  !<
    real, allocatable, dimension(:,:) :: QEVPROT !<
    real, allocatable, dimension(:,:) :: QFCFROT !<
    real, allocatable, dimension(:,:) :: QFCLROT !<
    real, allocatable, dimension(:,:) :: QFGROT  !<
    real, allocatable, dimension(:,:) :: QFNROT  !<
    real, allocatable, dimension(:,:) :: QFSROT  !<
    real, allocatable, dimension(:,:) :: QFXROT  !<
    real, allocatable, dimension(:,:) :: QGROT   !<
    real, allocatable, dimension(:,:) :: ROFROT  !<
    real, allocatable, dimension(:,:) :: ROFBROT !<
    real, allocatable, dimension(:,:) :: ROFCROT !<
    real, allocatable, dimension(:,:) :: ROFNROT !<
    real, allocatable, dimension(:,:) :: ROFOROT !<
    real, allocatable, dimension(:,:) :: ROFSROT !<
    real, allocatable, dimension(:,:) :: ROVGROT !<
    real, allocatable, dimension(:,:) :: SFCQROT !<
    real, allocatable, dimension(:,:) :: SFCTROT !<
    real, allocatable, dimension(:,:) :: SFCUROT !<
    real, allocatable, dimension(:,:) :: SFCVROT !<
    real, allocatable, dimension(:,:) :: TFXROT  !<
    real, allocatable, dimension(:,:) :: TROBROT !<
    real, allocatable, dimension(:,:) :: TROFROT !<
    real, allocatable, dimension(:,:) :: TROOROT !<
    real, allocatable, dimension(:,:) :: TROSROT !<
    real, allocatable, dimension(:,:) :: UEROT   !<
    real, allocatable, dimension(:,:) :: WTABROT !<
    real, allocatable, dimension(:,:) :: WTRCROT !<
    real, allocatable, dimension(:,:) :: WTRGROT !<
    real, allocatable, dimension(:,:) :: WTRSROT !<
    real, allocatable, dimension(:,:) :: SFRHROT !<
    real, allocatable, dimension(:,:) :: wtableROT !<
    real, allocatable, dimension(:,:) :: FTABLE !< Depth to frozen water table (m)
    real, allocatable, dimension(:,:) :: ACTLYR !< Active layer depth (m)
    real, allocatable, dimension(:,:) :: maxAnnualActLyrROT  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    real, allocatable, dimension(:,:) :: actLyrThisYrROT !< Annual active layer depth maximum starting from summer solstice for the present year (m)
    real, allocatable, dimension(:,:) :: groundHeatFluxROT !< Heat flux at soil surface \f$[W m^{-2} ]\f$

    ! There will be allocated the dimension: 'nlat,nmos,ignd'
    integer, allocatable, dimension(:,:,:) :: ISNDROT !< Sand content flag, used to delineate non-soils.
    real, allocatable, dimension(:,:,:) :: TBARROT !< Temperature of soil layers [K]
    real, allocatable, dimension(:,:,:) :: THICROT !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:,:) :: THLQROT !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:,:) :: BIROT   !<
    real, allocatable, dimension(:,:,:) :: DLZWROT !< Permeable thickness of soil layer [m]
    real, allocatable, dimension(:,:,:) :: GRKSROT !<
    real, allocatable, dimension(:,:,:) :: HCPSROT !<
    real, allocatable, dimension(:,:,:) :: SANDROT !< Percentage sand content of soil
    real, allocatable, dimension(:,:,:) :: CLAYROT !< Percentage clay content of soil
    real, allocatable, dimension(:,:,:) :: ORGMROT !< Percentage organic matter content of soil
    real, allocatable, dimension(:,:,:) :: PSISROT !<
    real, allocatable, dimension(:,:,:) :: PSIWROT !<
    real, allocatable, dimension(:,:,:) :: TCSROT  !<
    real, allocatable, dimension(:,:,:) :: THFCROT !<
    real, allocatable, dimension(:,:,:) :: THMROT  !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:,:) :: THPROT  !<
    real, allocatable, dimension(:,:,:) :: THRROT  !<
    real, allocatable, dimension(:,:,:) :: THRAROT !<
    real, allocatable, dimension(:,:,:) :: ZBTWROT !<
    real, allocatable, dimension(:,:,:) :: THLWROT !< Soil water content at wilting point, \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:,:) :: GFLXROT !<
    real, allocatable, dimension(:,:,:) :: HMFGROT !<
    real, allocatable, dimension(:,:,:) :: HTCROT  !<
    real, allocatable, dimension(:,:,:) :: QFCROT  !< Water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$

    ! allocated with nlat,nmos,ignd:
    real, allocatable, dimension(:,:,:) :: TBARACC_M        !< Temperature of soil layers [K] (accumulated for means)
    real, allocatable, dimension(:,:,:) :: THLQACC_M        !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated for means)
    real, allocatable, dimension(:,:,:) :: THICACC_M        !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated for means)
    real, allocatable, dimension(:,:,:) :: tbaraccrow_m     !< Temperature of soil layers [K] (accumulated for CTEM)


    ! These will be allocated the dimension: 'nlat,nmos,ican'
    real, allocatable, dimension(:,:,:) :: ACIDROT !<
    real, allocatable, dimension(:,:,:) :: ACVDROT !<
    real, allocatable, dimension(:,:,:) :: CMASROT !<
    real, allocatable, dimension(:,:,:) :: HGTDROT !<
    real, allocatable, dimension(:,:,:) :: PAIDROT !<
    real, allocatable, dimension(:,:,:) :: PAMNROT !<
    real, allocatable, dimension(:,:,:) :: PAMXROT !<
    real, allocatable, dimension(:,:,:) :: PSGAROT !<
    real, allocatable, dimension(:,:,:) :: PSGBROT !<
    real, allocatable, dimension(:,:,:) :: QA50ROT !<
    real, allocatable, dimension(:,:,:) :: ROOTROT !<
    real, allocatable, dimension(:,:,:) :: RSMNROT !<
    real, allocatable, dimension(:,:,:) :: VPDAROT !<
    real, allocatable, dimension(:,:,:) :: VPDBROT !<

    ! These will be allocated the dimension: 'nlat,nmos,icp1'
    real, allocatable, dimension(:,:,:) :: ALICROT !<
    real, allocatable, dimension(:,:,:) :: ALVCROT !<
    real, allocatable, dimension(:,:,:) :: FCANROT !<
    real, allocatable, dimension(:,:,:) :: LNZ0ROT !<

    ! These will be allocated the dimension: 'nlat,nmos,nbs'
    real, allocatable, dimension(:,:,:)  :: SALBROT  !<
    real, allocatable, dimension(:,:,:)  :: CSALROT  !<

    ! These will be allocated the dimension: 'nlat,nbs'
    real, allocatable, dimension(:,:) :: FSDBROL  !<
    real, allocatable, dimension(:,:) :: FSFBROL  !<
    real, allocatable, dimension(:,:) :: FSSBROL  !<

    ! These will be allocated the dimension: 'nlat,ignd'

    real, allocatable, dimension(:,:) :: TBARROW !< Temperature of soil layers [K]
    real, allocatable, dimension(:,:) :: THALROW !< Total volumetric water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THICROW !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: THLQROW !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: GFLXROW !<
    real, allocatable, dimension(:,:) :: HMFGROW !<
    real, allocatable, dimension(:,:) :: HTCROW  !<
    real, allocatable, dimension(:,:) :: QFCROW  !< Water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$

    ! These will be allocated the dimension: 'nlat,nmos,6,50'
    integer, allocatable, dimension(:,:,:,:) :: ITCTROT !<

    ! These will be allocated the dimension: 'nlat,nmos,4'
    real, allocatable, dimension(:,:,:)  :: TSFSROT !<

    ! allocated with nlat,nmos:
    real, allocatable, dimension(:,:) :: PREACC_M              !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:,:) :: GTACC_M               !< Diagnosed effective surface black-body temperature [K]
    real, allocatable, dimension(:,:) :: QEVPACC_M             !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HFSACC_M              !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: HMFNACC_M             !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: ROFACC_M              !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:,:) :: SNOACC_M              !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: OVRACC_M              !< Overland flow from top of soil column \f$[kg m^{-2} s^{-1}]\f$
    real, allocatable, dimension(:,:) :: WTBLACC_M             !< Depth of water table in soil [m]
    real, allocatable, dimension(:,:) :: ALVSACC_M             !< Diagnosed total visible albedo of land surface [ ]
    real, allocatable, dimension(:,:) :: ALIRACC_M             !< Diagnosed total near-infrared albedo of land surface [ ]
    real, allocatable, dimension(:,:) :: RHOSACC_M             !< Density of snow \f$[kg m^{-3} ]\f$
    real, allocatable, dimension(:,:) :: TSNOACC_M             !< Snowpack temperature [K]
    real, allocatable, dimension(:,:) :: WSNOACC_M             !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: TCANACC_M             !< Vegetation canopy temperature [K]
    real, allocatable, dimension(:,:) :: RCANACC_M             !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: SCANACC_M             !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: ALTOTACC_M            !< Daily broadband albedo
    real, allocatable, dimension(:,:) :: ALSNOACC_M            !< Daily snow albedo [ ]
    real, allocatable, dimension(:,:) :: GROACC_M              !< Vegetation growth index [ ]
    real, allocatable, dimension(:,:) :: FSINACC_M             !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: FLINACC_M             !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:,:) :: TAACC_M               !< Air temperature at reference height [K]
    real, allocatable, dimension(:,:) :: UVACC_M               !< Wind speed \f$[m s^{-1} ]\f$
    real, allocatable, dimension(:,:) :: PRESACC_M             !< Surface air pressure [Pa]
    real, allocatable, dimension(:,:) :: QAACC_M               !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
    real, allocatable, dimension(:,:) :: EVAPACC_M             !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
    real, allocatable, dimension(:,:) :: FLUTACC_M             !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$

  end type class_rotated

  type (class_rotated), save, target :: class_rot

  !=================================================================================
  !> CLASS's monthly and annual outputs
  type class_moyr_output

    !   MONTHLY OUTPUT FOR CLASS GRID-MEAN

    ! allocated with nlat:
    real, allocatable, dimension(:) :: ALVSACC_MO   !< Diagnosed total visible albedo of land surface [ ]
    real, allocatable, dimension(:) :: ALIRACC_MO   !< Diagnosed total near-infrared albedo of land surface [ ]
    real, allocatable, dimension(:) :: FLUTACC_MO   !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSINACC_MO   !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLINACC_MO   !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSACC_MO    !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: QEVPACC_MO   !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: groundHeatFlux_MO  !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: SNOACC_MO    !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: WSNOACC_MO   !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, allocatable, dimension(:) :: ROFACC_MO    !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PREACC_MO    !< Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$
    real, allocatable, dimension(:) :: EVAPACC_MO   !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
    real, allocatable, dimension(:) :: TRANSPACC_MO !<
    real, allocatable, dimension(:) :: TAACC_MO     !< Air temperature at reference height [K]
    real, allocatable, dimension(:) :: ACTLYR_MO
    real, allocatable, dimension(:) :: FTABLE_MO
    real, allocatable, dimension(:) :: ACTLYR_MIN_MO
    real, allocatable, dimension(:) :: ACTLYR_MAX_MO
    real, allocatable, dimension(:) :: FTABLE_MIN_MO
    real, allocatable, dimension(:) :: FTABLE_MAX_MO
    real, allocatable, dimension(:) :: MRSO_MO      !< Total Soil Moisture Content [kg $m^{-2}$]
    real, allocatable, dimension(:,:) :: MRSOL_MO     !< Total water content of soil layer [kg $m^{-2}$]
    real, allocatable, dimension(:) :: ALTOTACC_MO  !< Broadband albedo [ ]
    real, allocatable, dimension(:) :: ALSNOACC_MO  !< Snow albedo [ ]
    real, allocatable, dimension(:) :: GROUNDEVAP   !< evaporation and sublimation from the ground surface (formed from QFG and QFN), kg /m/mon
    real, allocatable, dimension(:) :: CANOPYEVAP   !< evaporation and sublimation from the canopy (formed from QFCL and QFCF), kg /m/mon
    integer, allocatable, dimension(:) :: altotcntr_m !< Used to count the number of time steps with the sun above the horizon

    ! allocated with nlat,ignd:
    real, allocatable, dimension(:,:) :: TBARACC_MO !< Temperature of soil layers [K] (accumulated for means)
    real, allocatable, dimension(:,:) :: THLQACC_MO !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated for means)
    real, allocatable, dimension(:,:) :: THICACC_MO !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$ (accumulated for means)

    !   YEARLY OUTPUT FOR CLASS GRID-MEAN

    real, allocatable, dimension(:) :: ALVSACC_YR  !< Diagnosed total visible albedo of land surface [ ]
    real, allocatable, dimension(:) :: ALIRACC_YR  !< Diagnosed total near-infrared albedo of land surface [ ]
    real, allocatable, dimension(:) :: FLUTACC_YR  !< Upwelling longwave radiation from surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FSINACC_YR  !< Downwelling shortwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: FLINACC_YR  !< Downwelling longwave radiation above surface \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: HFSACC_YR   !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: QEVPACC_YR  !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
    real, allocatable, dimension(:) :: ROFACC_YR   !< Total runoff from soil \f$[kg m^{-2} s^{-1} ]\f$
    real, allocatable, dimension(:) :: PREACC_YR   !< Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$
    real, allocatable, dimension(:) :: EVAPACC_YR  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
    real, allocatable, dimension(:) :: TRANSPACC_YR !<
    real, allocatable, dimension(:) :: TAACC_YR    !< Air temperature at reference height [K]
    real, allocatable, dimension(:) :: ACTLYR_YR
    real, allocatable, dimension(:) :: ACTLYR_MIN_YR
    real, allocatable, dimension(:) :: ACTLYR_MAX_YR
    real, allocatable, dimension(:) :: FTABLE_YR
    real, allocatable, dimension(:) :: FTABLE_MIN_YR
    real, allocatable, dimension(:) :: FTABLE_MAX_YR
    real, allocatable, dimension(:) :: ALTOTACC_YR !< Broadband albedo
    integer, allocatable, dimension(:) :: altotcntr_yr !< Used to count the number of time steps with the sun above the horizon

  end type class_moyr_output

  type (class_moyr_output), save, target :: class_out

  !=================================================================================

contains

  !> \ingroup classstatevars_alloc
  !! @{
  !> Allocates the CLASS (physics) variables in preparation for the simulation
  subroutine allocClassVars

    use classicParams,     only : ican, nbs, icp1, nlat, nmos, ignd, ilg

    implicit none

    ! ------------------

    ! ilg
	!nlat = vs%active_cell_count
!	nlat = vs%cell_count
!	nmos = vs%gru_count
!	ilg = vs%landtile_count
	!ilg = nlat*nmos
!    ignd = vs%sol_count
	
    allocate(class_gat%TEST (4), &
			 class_gat%ILMOS   (ilg), &
             class_gat%JLMOS   (ilg), &
             class_gat%IWMOS   (ilg), &
             class_gat%JWMOS   (ilg), &
             class_gat%IGDRGAT  (ilg), &
             class_gat%ALBSGAT (ilg), &
             class_gat%GCGAT (ilg), &
             class_gat%TZSGAT (ilg), &
             class_gat%PCSNGAT (ilg), &
             class_gat%CMAIGAT (ilg), &
             class_gat%GROGAT  (ilg), &
             class_gat%QACGAT  (ilg), &
             class_gat%RCANGAT(ilg), &
             class_gat%RHOSGAT(ilg), &
             class_gat%SCANGAT(ilg), &
             class_gat%SNOGAT(ilg), &
             class_gat%TACGAT(ilg), &
             class_gat%TBASGAT(ilg), &
             class_gat%TCANGAT(ilg), &
             class_gat%TPNDGAT(ilg), &
             class_gat%TSNOGAT(ilg), &
             class_gat%WSNOGAT(ilg), &
             class_gat%maxAnnualActLyrGAT(ilg), &
             class_gat%ZPNDGAT(ilg), &
             class_gat%REFGAT(ilg), &
             class_gat%BCSNGAT(ilg), &
             class_gat%AGIDGAT(ilg), &
             class_gat%AGVDGAT(ilg), &
             class_gat%ALGDGAT(ilg), &
             class_gat%ALGWGAT(ilg), &
             class_gat%ASIDGAT(ilg), &
             class_gat%ASVDGAT(ilg), &
             class_gat%DRNGAT(ilg), &
             class_gat%GRKFGAT(ilg), &
             class_gat%WFCIGAT(ilg), &
             class_gat%WFSFGAT(ilg), &
             class_gat%XSLPGAT(ilg), &
             class_gat%ZPLGGAT(ilg), &
             class_gat%ZPLSGAT(ilg), &
             class_gat%ZSNLGAT(ilg), &
             class_gat%ALGWVGAT(ilg), &
             class_gat%ALGWNGAT(ilg), &
             class_gat%ALGDVGAT(ilg), &
             class_gat%ALGDNGAT(ilg), &
             class_gat%EMISGAT(ilg), &
             class_gat%CSZGAT(ilg), &
             class_gat%DLONGAT(ilg), &
             class_gat%DLATGAT(ilg), &
             class_gat%FCLOGAT(ilg), &
             class_gat%FDLGAT(ilg), &
             class_gat%FSSGAT(ilg), &
             class_gat%FSIHGAT(ilg), &
             class_gat%FSVHGAT(ilg), &
             class_gat%GGEOGAT(ilg), &
             class_gat%PADRGAT(ilg), &
             class_gat%PREGAT(ilg), &
             class_gat%PRESGAT(ilg), &
             class_gat%QAGAT(ilg), &
             class_gat%RADJGAT(ilg), &
             class_gat%RHOAGAT(ilg), &
             class_gat%RHSIGAT(ilg), &
             class_gat%RPCPGAT(ilg), &
             class_gat%SPCPGAT(ilg), &
             class_gat%TAGAT(ilg), &
             class_gat%TADPGAT(ilg), &
             class_gat%TRPCGAT(ilg), &
             class_gat%TSPCGAT(ilg), &
             class_gat%ULGAT(ilg), &
             class_gat%VLGAT(ilg), &
             class_gat%VMODGAT(ilg), &
             class_gat%VPDGAT(ilg), &
             class_gat%Z0ORGAT(ilg), &
             class_gat%ZBLDGAT(ilg), &
             class_gat%ZDHGAT(ilg), &
             class_gat%ZDMGAT(ilg), &
             class_gat%ZRFHGAT(ilg), &
             class_gat%ZRFMGAT(ilg), &
             class_gat%FSGGAT(ilg), &
             class_gat%FLGGAT(ilg), &
             class_gat%GUSTGAT(ilg), &
             class_gat%DEPBGAT(ilg), &
             class_gat%GTBS(ilg), &
             class_gat%SFCUBS(ilg), &
             class_gat%SFCVBS(ilg), &
             class_gat%USTARBS(ilg), &
             class_gat%TCSNOW(ilg), &
             class_gat%GSNOW(ilg), &
             class_gat%ALIRGAT(ilg), &
             class_gat%ALVSGAT(ilg), &
             class_gat%CDHGAT(ilg), &
             class_gat%CDMGAT(ilg), &
             class_gat%DRGAT(ilg), &
             class_gat%EFGAT(ilg), &
             class_gat%FLGGGAT(ilg), &
             class_gat%FLGSGAT(ilg), &
             class_gat%FLGVGAT(ilg), &
             class_gat%FSGGGAT(ilg), &
             class_gat% FSGSGAT (ilg), &
             class_gat% FSGVGAT (ilg), &
             class_gat% FSNOGAT (ilg), &
             class_gat% GAGAT   (ilg), &
             class_gat% GTGAT   (ilg), &
             class_gat% HBLGAT  (ilg), &
             class_gat% HEVCGAT (ilg), &
             class_gat% HEVGGAT (ilg), &
             class_gat% HEVSGAT (ilg), &
             class_gat% HFSGAT  (ilg), &
             class_gat% HFSCGAT (ilg), &
             class_gat% HFSGGAT (ilg), &
             class_gat% HFSSGAT (ilg), &
             class_gat% HMFCGAT (ilg), &
             class_gat% HMFNGAT (ilg), &
             class_gat% HTCCGAT (ilg), &
             class_gat% HTCSGAT (ilg), &
             class_gat% ILMOGAT (ilg), &
             class_gat% PCFCGAT (ilg), &
             class_gat% PCLCGAT (ilg), &
             class_gat% PCPGGAT (ilg), &
             class_gat% PCPNGAT (ilg), &
             class_gat% PETGAT  (ilg), &
             class_gat% QEVPGAT (ilg), &
             class_gat% QFCFGAT (ilg), &
             class_gat% QFCLGAT (ilg), &
             class_gat% QFGGAT  (ilg), &
             class_gat% QFNGAT  (ilg), &
             class_gat% QFSGAT  (ilg), &
             class_gat% QFXGAT  (ilg), &
             class_gat% QGGAT   (ilg), &
             class_gat% ROFGAT  (ilg), &
             class_gat% ROFBGAT (ilg), &
             class_gat% ROFCGAT (ilg), &
             class_gat% ROFNGAT (ilg), &
             class_gat% ROFOGAT (ilg), &
             class_gat% ROFSGAT (ilg,ignd), &
             class_gat% ROVGGAT (ilg), &
             class_gat% SFCQGAT (ilg), &
             class_gat% SFCTGAT (ilg), &
             class_gat% SFCUGAT (ilg), &
             class_gat% SFCVGAT (ilg), &
             class_gat% TFXGAT  (ilg), &
             class_gat% TROBGAT (ilg), &
             class_gat% TROFGAT (ilg), &
             class_gat% TROOGAT (ilg), &
             class_gat% TROSGAT (ilg,ignd), &
             class_gat% UEGAT   (ilg), &
             class_gat% WTABGAT (ilg), &
             class_gat% WTRCGAT (ilg), &
             class_gat% WTRGGAT (ilg), &
             class_gat% WTRSGAT (ilg), &
             class_gat% wtableGAT (ilg), &
             class_gat%FTABLE(ilg), &
             class_gat%ACTLYR(ilg), &
             class_gat%actLyrThisYrGAT(ilg), &
             class_gat% QLWOGAT (ilg), &
             class_gat% SFRHGAT (ilg), &
             class_gat% FTEMP   (ilg), &
             class_gat% FVAP    (ilg), &
             class_gat% RIB     (ilg), &
             class_gat% FC      (ilg), &
             class_gat% FG      (ilg), &
             class_gat% FCS     (ilg), &
             class_gat% FGS     (ilg), &
             class_gat% RBCOEF  (ilg), &
             class_gat% ZSNOW   (ilg), &
             class_gat% FSVF    (ilg), &
             class_gat% FSVFS   (ilg), &
             class_gat% ALVSCN  (ilg), &
             class_gat% ALIRCN  (ilg), &
             class_gat% ALVSG   (ilg), &
             class_gat% ALIRG   (ilg), &
             class_gat% ALVSCS  (ilg), &
             class_gat% ALIRCS  (ilg), &
             class_gat% ALVSSN  (ilg), &
             class_gat% ALIRSN  (ilg), &
             class_gat% ALVSGC  (ilg), &
             class_gat% ALIRGC  (ilg), &
             class_gat% ALVSSC  (ilg), &
             class_gat% ALIRSC  (ilg), &
             class_gat% TRVSCN  (ilg), &
             class_gat% TRIRCN  (ilg), &
             class_gat%TRVSCS  (ilg), &
             class_gat%TRIRCS  (ilg), &
             class_gat%RC      (ilg), &
             class_gat%RCS     (ilg), &
             class_gat%FRAINC  (ilg), &
             class_gat%FSNOWC  (ilg), &
             class_gat%FRAICS  (ilg), &
             class_gat%FSNOCS  (ilg), &
             class_gat%CMASSC  (ilg), &
             class_gat%CMASCS  (ilg), &
             class_gat%DISP    (ilg), &
             class_gat%DISPS   (ilg), &
             class_gat%ZOMLNC  (ilg), &
             class_gat%ZOELNC  (ilg), &
             class_gat%ZOMLNG  (ilg), &
             class_gat%ZOELNG  (ilg), &
             class_gat%ZOMLCS  (ilg), &
             class_gat%ZOELCS  (ilg), &
             class_gat%ZOMLNS  (ilg), &
             class_gat%ZOELNS  (ilg), &
             class_gat%TRSNOWC (ilg), &
             class_gat%CHCAP   (ilg), &
             class_gat%CHCAPS  (ilg), &
             class_gat%GZEROC  (ilg), &
             class_gat%GZEROG  (ilg), &
             class_gat%GZROCS  (ilg), &
             class_gat%GZROGS  (ilg), &
             class_gat%groundHeatFlux  (ilg), &
             class_gat%G12C    (ilg), &
             class_gat%G12G    (ilg), &
             class_gat%G12CS   (ilg), &
             class_gat%G12GS   (ilg), &
             class_gat%G23C    (ilg), &
             class_gat%G23G    (ilg), &
             class_gat%G23CS   (ilg), &
             class_gat%G23GS   (ilg), &
             class_gat%QFREZC  (ilg), &
             class_gat%QFREZG  (ilg), &
             class_gat%QMELTC  (ilg), &
             class_gat%QMELTG  (ilg), &
             class_gat%EVAPC   (ilg), &
             class_gat%EVAPCG  (ilg), &
             class_gat%EVAPG   (ilg), &
             class_gat%EVAPCS  (ilg), &
             class_gat%EVPCSG  (ilg), &
             class_gat%EVAPGS  (ilg), &
             class_gat%TCANO   (ilg), &
             class_gat%TCANS   (ilg), &
             class_gat%RAICAN  (ilg), &
             class_gat%SNOCAN  (ilg), &
             class_gat%RAICNS  (ilg), &
             class_gat%SNOCNS  (ilg), &
             class_gat%CWLCAP  (ilg), &
             class_gat%CWFCAP  (ilg), &
             class_gat%CWLCPS  (ilg), &
             class_gat%CWFCPS  (ilg), &
             class_gat%TSNOCS  (ilg), &
             class_gat%TSNOGS  (ilg), &
             class_gat%RHOSCS  (ilg), &
             class_gat%RHOSGS  (ilg), &
             class_gat%WSNOCS  (ilg), &
             class_gat%WSNOGS  (ilg), &
             class_gat%TPONDC  (ilg), &
             class_gat%TPONDG  (ilg), &
             class_gat%TPNDCS  (ilg), &
             class_gat%TPNDGS  (ilg), &
             class_gat%ZPLMCS  (ilg), &
             class_gat%ZPLMGS  (ilg), &
             class_gat%ZPLIMC  (ilg), &
             class_gat%ZPLIMG  (ilg), &
             class_gat%CTVSTP (ilg), &
             class_gat%CTSSTP (ilg), &
             class_gat%CT1STP (ilg), &
             class_gat%CT2STP (ilg), &
             class_gat%CT3STP (ilg), &
             class_gat%WTVSTP (ilg), &
             class_gat%WTSSTP (ilg), &
             class_gat%WTGSTP (ilg))

    ! These will be allocated the dimension: 'ignd'
    allocate(class_gat% DELZ   (ignd), &
             class_gat% ZBOT   (ignd))

    ! These will be allocated the dimension: 'ilg,ignd'
    allocate(class_gat% ISNDGAT (ilg,ignd), &
             class_gat% TBARGAT (ilg,ignd), &
             class_gat% THICGAT (ilg,ignd), &
             class_gat% THLQGAT (ilg,ignd), &
             class_gat% BIGAT   (ilg,ignd), &
             class_gat% DLZWGAT (ilg,ignd), &
             class_gat% GRKSGAT (ilg,ignd), &
             class_gat% HCPSGAT (ilg,ignd), &
             class_gat% PSISGAT (ilg,ignd), &
             class_gat% PSIWGAT (ilg,ignd), &
             class_gat% TCSGAT  (ilg,ignd), &
             class_gat% THFCGAT (ilg,ignd), &
             class_gat% THMGAT  (ilg,ignd), &
             class_gat% THPGAT  (ilg,ignd), &
             class_gat% THRGAT  (ilg,ignd), &
             class_gat% THRAGAT (ilg,ignd), &
             class_gat% ZBTWGAT (ilg,ignd), &
             class_gat% THLWGAT (ilg,ignd), &
             class_gat% GFLXGAT (ilg,ignd), &
             class_gat% HMFGGAT (ilg,ignd), &
             class_gat% HTCGAT  (ilg,ignd), &
             class_gat% QFCGAT  (ilg,ignd), &
             class_gat% TBARC  (ilg,ignd), &
             class_gat% TBARG  (ilg,ignd), &
             class_gat% TBARCS (ilg,ignd), &
             class_gat% TBARGS (ilg,ignd), &
             class_gat% THLIQC (ilg,ignd), &
             class_gat% THLIQG (ilg,ignd), &
             class_gat% THICEC (ilg,ignd), &
             class_gat% THICEG (ilg,ignd), &
             class_gat% FROOT  (ilg,ignd), &
             class_gat% HCPC   (ilg,ignd), &
             class_gat% HCPG   (ilg,ignd), &
             class_gat% FROOTS (ilg,ignd), &
             class_gat% TCTOPC (ilg,ignd), &
             class_gat% TCBOTC (ilg,ignd), &
             class_gat% TCTOPG (ilg,ignd), &
             class_gat% TCBOTG (ilg,ignd))

    ! These will be allocated the dimension: 'ilg,ican'
    allocate(class_gat% ACIDGAT (ilg,ican), &
             class_gat% ACVDGAT (ilg,ican), &
             class_gat% CMASGAT (ilg,ican), &
             class_gat% HGTDGAT (ilg,ican), &
             class_gat% PAIDGAT (ilg,ican), &
             class_gat% PAMNGAT (ilg,ican), &
             class_gat% PAMXGAT (ilg,ican), &
             class_gat% PSGAGAT (ilg,ican), &
             class_gat% PSGBGAT (ilg,ican), &
             class_gat% QA50GAT (ilg,ican), &
             class_gat% ROOTGAT (ilg,ican), &
             class_gat% RSMNGAT (ilg,ican), &
             class_gat% VPDAGAT (ilg,ican), &
             class_gat% VPDBGAT (ilg,ican))


    ! These will be allocated the dimension: 'ilg,icp1'
    allocate(class_gat% ALICGAT (ilg,icp1), &
             class_gat% ALVCGAT (ilg,icp1), &
             class_gat% FCANGAT (ilg,icp1), &
             class_gat% LNZ0GAT (ilg,icp1))

    ! These will be allocated the dimension: 'ilg,nbs'
    allocate(class_gat% FSDBGAT (ilg,nbs), &
             class_gat% FSFBGAT (ilg,nbs), &
             class_gat% FSSBGAT (ilg,nbs), &
             class_gat% SALBGAT (ilg,nbs), &
             class_gat% CSALGAT (ilg,nbs), &
             class_gat% ALTG    (ilg,nbs), &
             class_gat% ALSNO   (ilg,nbs), &
             class_gat% TRSNOWG (ilg,nbs))

    ! These will be allocated the dimension: 'ilg,4'
    allocate(class_gat% TSFSGAT (ilg,4))

    ! These will be allocated the dimension: 'ilg,6,50'
    allocate(class_gat% ITCTGAT (ilg,6,50))

    ! -----------------------------------------------------------
    ! Now allocate the class_rot structure:

    ! These will be allocated the dimension: 'nlat'

    allocate(class_rot% CSZROW  (nlat), &
             class_rot% DLONROW (nlat), &
             class_rot% DLATROW (nlat), &
             class_rot% FCLOROW (nlat), &
             class_rot% RHOSROW (nlat), &
             class_rot% FDLROW  (nlat), &
             class_rot% FSIHROW (nlat), &
             class_rot% FSVHROW (nlat), &
             class_rot% GCROW   (nlat), &
             class_rot% GGEOROW (nlat), &
             class_rot% PADRROW (nlat), &
             class_rot% PREROW  (nlat), &
             class_rot% PRESROW (nlat), &
             class_rot% QAROW   (nlat), &
             class_rot% RADJROW (nlat), &
             class_rot% RHOAROW (nlat), &
             class_rot% RHSIROW (nlat), &
             class_rot% RPCPROW (nlat), &
             class_rot% RPREROW (nlat), &
             class_rot% SPCPROW (nlat), &
             class_rot% SPREROW (nlat), &
             class_rot% TAROW   (nlat), &
             class_rot% TSNOROW (nlat), &
             class_rot% WSNOROW (nlat), &
             class_rot% TCANROW (nlat), &
             class_rot% SCANROW (nlat), &
             class_rot% RCANROW (nlat), &
             class_rot% TPNDROW (nlat), &
             class_rot% ZPNDROW (nlat), &
             class_rot% TADPROW (nlat), &
             class_rot% TRPCROW (nlat), &
             class_rot% TSPCROW (nlat), &
             class_rot% ULROW   (nlat), &
             class_rot% VLROW   (nlat), &
             class_rot% VMODROW (nlat), &
             class_rot% VPDROW  (nlat), &
             class_rot% ZBLDROW (nlat), &
             class_rot% ZDHROW  (nlat), &
             class_rot% ZDMROW  (nlat), &
             class_rot% ZRFHROW (nlat), &
             class_rot% ZRFMROW (nlat), &
             class_rot% UVROW   (nlat), &
             class_rot% XDIFFUS (nlat), &
             class_rot% Z0ORROW (nlat), &
             class_rot% FSSROW  (nlat), &
             class_rot% PRENROW (nlat), &
             class_rot% CLDTROW (nlat), &
             class_rot% FSGROL  (nlat), &
             class_rot% FLGROL  (nlat), &
             class_rot% GUSTROL (nlat), &
             class_rot% DEPBROW (nlat), &
             class_rot% ALIRROW (nlat), &
             class_rot% ALVSROW (nlat), &
             class_rot% CDHROW  (nlat), &
             class_rot% CDMROW  (nlat), &
             class_rot% DRROW   (nlat), &
             class_rot% EFROW   (nlat), &
             class_rot% FLGGROW (nlat), &
             class_rot% FLGSROW (nlat), &
             class_rot% FLGVROW (nlat), &
             class_rot% FSGGROW (nlat), &
             class_rot% FSGSROW (nlat), &
             class_rot% FSGVROW (nlat), &
             class_rot% FSNOROW (nlat), &
             class_rot% GAROW   (nlat), &
             class_rot% GTROW   (nlat), &
             class_rot% HBLROW  (nlat), &
             class_rot% HEVCROW (nlat), &
             class_rot% HEVGROW (nlat), &
             class_rot% HEVSROW (nlat), &
             class_rot% HFSROW  (nlat), &
             class_rot% HFSCROW (nlat), &
             class_rot% HFSGROW (nlat), &
             class_rot% HFSSROW (nlat), &
             class_rot% HMFCROW (nlat), &
             class_rot% HMFNROW (nlat), &
             class_rot% HTCCROW (nlat), &
             class_rot% HTCSROW (nlat), &
             class_rot% ILMOROW (nlat), &
             class_rot% PCFCROW (nlat), &
             class_rot% PCLCROW (nlat), &
             class_rot% PCPGROW (nlat), &
             class_rot% PCPNROW (nlat), &
             class_rot% PETROW  (nlat), &
             class_rot% QEVPROW (nlat), &
             class_rot% QFCFROW (nlat), &
             class_rot% QFCLROW (nlat), &
             class_rot% QFGROW  (nlat), &
             class_rot% QFNROW  (nlat), &
             class_rot% QFSROW  (nlat), &
             class_rot% QFXROW  (nlat), &
             class_rot% QGROW   (nlat), &
             class_rot% ROFROW  (nlat), &
             class_rot% ROFBROW (nlat), &
             class_rot% ROFCROW (nlat), &
             class_rot% ROFNROW (nlat), &
             class_rot% ROFOROW (nlat), &
             class_rot% ROFSROW (nlat), &
             class_rot% ROVGROW (nlat), &
             class_rot% SFCQROW (nlat), &
             class_rot% SFCTROW (nlat), &
             class_rot% SFCUROW (nlat), &
             class_rot% SFCVROW (nlat), &
             class_rot% TFXROW  (nlat), &
             class_rot% UEROW   (nlat), &
             class_rot% WTABROW (nlat), &
             class_rot% WTRCROW (nlat), &
             class_rot% WTRGROW (nlat), &
             class_rot% WTRSROW (nlat), &
             class_rot% SFRHROW (nlat), &
             class_rot% SNOROW (nlat), &
             class_rot% wtableROW (nlat), &
             class_rot%altotcntr_d(nlat), &

    ! allocated with nlat:
             class_out%ALVSACC_MO(nlat), &
             class_out%ALIRACC_MO (nlat), &
             class_out%FLUTACC_MO (nlat), &
             class_out%FSINACC_MO (nlat), &
             class_out%FLINACC_MO (nlat), &
             class_out%HFSACC_MO (nlat), &
             class_out%QEVPACC_MO (nlat), &
             class_out%groundHeatFlux_MO (nlat), &
             class_out%SNOACC_MO (nlat), &
             class_out%WSNOACC_MO (nlat), &
             class_out%ROFACC_MO (nlat), &
             class_out%PREACC_MO (nlat), &
             class_out%EVAPACC_MO (nlat), &
             class_out%TRANSPACC_MO (nlat), &
             class_out%TAACC_MO (nlat), &
             class_out%ACTLYR_MO (nlat), &
             class_out%FTABLE_MO (nlat), &
             class_out%ACTLYR_MIN_MO (nlat), &
             class_out%ACTLYR_MAX_MO (nlat), &
             class_out%FTABLE_MIN_MO (nlat), &
             class_out%FTABLE_MAX_MO (nlat), &
             class_out%ALTOTACC_MO (nlat), &
             class_out%ALSNOACC_MO (nlat), &
             class_out%MRSO_MO (nlat), &
             class_out%GROUNDEVAP (nlat), &
             class_out%CANOPYEVAP (nlat), &
             class_out%altotcntr_m (nlat), &

    ! allocated with nlat,ignd:
             class_out%TBARACC_MO (nlat,ignd), &
             class_out%THLQACC_MO (nlat,ignd), &
             class_out%THICACC_MO (nlat,ignd), &
             class_out%MRSOL_MO (nlat,ignd), &

             class_out%ALVSACC_YR (nlat), &
             class_out%ALIRACC_YR (nlat), &
             class_out%FLUTACC_YR (nlat), &
             class_out%FSINACC_YR (nlat), &
             class_out%FLINACC_YR (nlat), &
             class_out%HFSACC_YR (nlat), &
             class_out%QEVPACC_YR (nlat), &
             class_out%ROFACC_YR (nlat), &
             class_out%PREACC_YR (nlat), &
             class_out%EVAPACC_YR (nlat), &
             class_out%TRANSPACC_YR (nlat), &
             class_out%TAACC_YR (nlat), &
             class_out%ACTLYR_YR (nlat), &
             class_out%ACTLYR_MIN_YR (nlat), &
             class_out%ACTLYR_MAX_YR (nlat), &
             class_out%FTABLE_YR (nlat), &
             class_out%FTABLE_MIN_YR (nlat), &
             class_out%FTABLE_MAX_YR (nlat), &
             class_out%ALTOTACC_YR (nlat), &
             class_out%altotcntr_yr (nlat))


    ! These will be allocated the dimension: 'nlat,nmos'

    allocate(class_rot% IGDRROT (nlat,nmos), &
             class_rot% MIDROT  (nlat,nmos), &
             class_rot% ALBSROT (nlat,nmos), &
             class_rot% CMAIROT (nlat,nmos), &
             class_rot% GROROT  (nlat,nmos), &
             class_rot% QACROT  (nlat,nmos), &
             class_rot% RCANROT (nlat,nmos), &
             class_rot% RHOSROT (nlat,nmos), &
             class_rot% SCANROT (nlat,nmos), &
             class_rot% SNOROT  (nlat,nmos), &
             class_rot% TACROT  (nlat,nmos), &
             class_rot% TBASROT (nlat,nmos), &
             class_rot% TCANROT (nlat,nmos), &
             class_rot% TPNDROT (nlat,nmos), &
             class_rot% TSNOROT (nlat,nmos), &
             class_rot% WSNOROT (nlat,nmos), &
             class_rot% ZPNDROT (nlat,nmos), &
             class_rot% REFROT  (nlat,nmos), &
             class_rot% BCSNROT (nlat,nmos), &
             class_rot% AGIDROT (nlat,nmos), &
             class_rot% AGVDROT (nlat,nmos), &
             class_rot% ALGDROT (nlat,nmos), &
             class_rot% ALGWROT (nlat,nmos), &
             class_rot% ASIDROT (nlat,nmos), &
             class_rot% ASVDROT (nlat,nmos), &
             class_rot% DRNROT  (nlat,nmos), &
             class_rot% FAREROT (nlat,nmos), &
             class_rot% GRKFROT (nlat,nmos), &
             class_rot% WFCIROT (nlat,nmos), &
             class_rot% WFSFROT (nlat,nmos), &
             class_rot% XSLPROT (nlat,nmos), &
             class_rot% ZPLGROT (nlat,nmos), &
             class_rot% ZPLSROT (nlat,nmos), &
             class_rot% ZSNLROT (nlat,nmos), &
             class_rot% ZSNOROT  (nlat,nmos), &
             class_rot% ALGWVROT (nlat,nmos), &
             class_rot% ALGWNROT (nlat,nmos), &
             class_rot% ALGDVROT (nlat,nmos), &
             class_rot% ALGDNROT (nlat,nmos), &
             class_rot% EMISROT  (nlat,nmos), &
             class_rot% ALIRROT (nlat,nmos), &
             class_rot% ALVSROT (nlat,nmos), &
             class_rot% CDHROT  (nlat,nmos), &
             class_rot% CDMROT  (nlat,nmos), &
             class_rot% DRROT   (nlat,nmos), &
             class_rot% EFROT   (nlat,nmos), &
             class_rot% FLGGROT (nlat,nmos), &
             class_rot% FLGSROT (nlat,nmos), &
             class_rot% FLGVROT (nlat,nmos), &
             class_rot% FSGGROT (nlat,nmos), &
             class_rot% FSGSROT (nlat,nmos), &
             class_rot% FSGVROT (nlat,nmos), &
             class_rot% FSNOROT (nlat,nmos), &
             class_rot% GAROT   (nlat,nmos), &
             class_rot% GTROT   (nlat,nmos), &
             class_rot% HBLROT  (nlat,nmos), &
             class_rot% HEVCROT (nlat,nmos), &
             class_rot% HEVGROT (nlat,nmos), &
             class_rot% HEVSROT (nlat,nmos), &
             class_rot% HFSROT  (nlat,nmos), &
             class_rot% HFSCROT (nlat,nmos), &
             class_rot% HFSGROT (nlat,nmos), &
             class_rot% HFSSROT (nlat,nmos), &
             class_rot% HMFCROT (nlat,nmos), &
             class_rot% HMFNROT (nlat,nmos), &
             class_rot% HTCCROT (nlat,nmos), &
             class_rot% SDEPROT (nlat,nmos), &
             class_rot% SOCIROT (nlat,nmos), &
             class_rot% HTCSROT (nlat,nmos), &
             class_rot% ILMOROT (nlat,nmos), &
             class_rot% PCFCROT (nlat,nmos), &
             class_rot% PCLCROT (nlat,nmos), &
             class_rot% PCPGROT (nlat,nmos), &
             class_rot% PCPNROT (nlat,nmos), &
             class_rot% PETROT  (nlat,nmos), &
             class_rot% QEVPROT (nlat,nmos), &
             class_rot% QFCFROT (nlat,nmos), &
             class_rot% QFCLROT (nlat,nmos), &
             class_rot% QFGROT  (nlat,nmos), &
             class_rot% QFNROT  (nlat,nmos), &
             class_rot% QFSROT  (nlat,nmos), &
             class_rot% QFXROT  (nlat,nmos), &
             class_rot% QGROT   (nlat,nmos), &
             class_rot% ROFROT  (nlat,nmos), &
             class_rot% ROFBROT (nlat,nmos), &
             class_rot% ROFCROT (nlat,nmos), &
             class_rot% ROFNROT (nlat,nmos), &
             class_rot% ROFOROT (nlat,nmos), &
             class_rot% ROFSROT (nlat,nmos), &
             class_rot% ROVGROT (nlat,nmos), &
             class_rot% SFCQROT (nlat,nmos), &
             class_rot% SFCTROT (nlat,nmos), &
             class_rot% SFCUROT (nlat,nmos), &
             class_rot% SFCVROT (nlat,nmos), &
             class_rot% TFXROT  (nlat,nmos), &
             class_rot% TROBROT (nlat,nmos), &
             class_rot% TROFROT (nlat,nmos), &
             class_rot% TROOROT (nlat,nmos), &
             class_rot% TROSROT (nlat,nmos), &
             class_rot% UEROT   (nlat,nmos), &
             class_rot% WTABROT (nlat,nmos), &
             class_rot% WTRCROT (nlat,nmos), &
             class_rot% WTRGROT (nlat,nmos), &
             class_rot% WTRSROT (nlat,nmos), &
             class_rot% SFRHROT (nlat,nmos), &
             class_rot% wtableROT(nlat,nmos), &
             class_rot% ACTLYR(nlat,nmos), &
             class_rot% maxAnnualActLyrROT(nlat,nmos), &
             class_rot% actLyrThisYrROT(nlat,nmos), &
             class_rot% groundHeatFluxROT(nlat,nmos), &
             class_rot% FTABLE(nlat,nmos), &
             class_rot%PREACC_M(nlat,nmos), &
             class_rot%GTACC_M (nlat,nmos), &
             class_rot%QEVPACC_M (nlat,nmos), &
             class_rot%HFSACC_M(nlat,nmos), &
             class_rot%HMFNACC_M (nlat,nmos), &
             class_rot%ROFACC_M(nlat,nmos), &
             class_rot%SNOACC_M(nlat,nmos), &
             class_rot%OVRACC_M(nlat,nmos), &
             class_rot%WTBLACC_M(nlat,nmos), &
             class_rot%ALVSACC_M(nlat,nmos), &
             class_rot%ALIRACC_M(nlat,nmos), &
             class_rot%RHOSACC_M(nlat,nmos), &
             class_rot%TSNOACC_M(nlat,nmos), &
             class_rot%WSNOACC_M(nlat,nmos), &
             class_rot%TCANACC_M(nlat,nmos), &
             class_rot%RCANACC_M(nlat,nmos), &
             class_rot%SCANACC_M(nlat,nmos), &
             class_rot%ALTOTACC_M(nlat,nmos), &
             class_rot%ALSNOACC_M(nlat,nmos), &
             class_rot%GROACC_M(nlat,nmos), &
             class_rot%FSINACC_M (nlat,nmos), &
             class_rot%FLINACC_M(nlat,nmos), &
             class_rot%TAACC_M (nlat,nmos), &
             class_rot%UVACC_M (nlat,nmos), &
             class_rot%PRESACC_M (nlat,nmos), &
             class_rot%QAACC_M (nlat,nmos), &
             class_rot%EVAPACC_M (nlat,nmos), &
             class_rot%FLUTACC_M(nlat,nmos))

    ! There will be allocated the dimension: 'nlat,nmos,ignd'
    allocate(class_rot% ISNDROT (nlat,nmos,ignd), &
             class_rot% TBARROT (nlat,nmos,ignd), &
             class_rot% THICROT (nlat,nmos,ignd), &
             class_rot% THLQROT (nlat,nmos,ignd), &
             class_rot% BIROT   (nlat,nmos,ignd), &
             class_rot% DLZWROT (nlat,nmos,ignd), &
             class_rot% GRKSROT (nlat,nmos,ignd), &
             class_rot% HCPSROT (nlat,nmos,ignd), &
             class_rot% SANDROT (nlat,nmos,ignd), &
             class_rot% CLAYROT (nlat,nmos,ignd), &
             class_rot% ORGMROT (nlat,nmos,ignd), &
             class_rot% PSISROT (nlat,nmos,ignd), &
             class_rot% PSIWROT (nlat,nmos,ignd), &
             class_rot% TCSROT  (nlat,nmos,ignd), &
             class_rot% THFCROT (nlat,nmos,ignd), &
             class_rot% THMROT  (nlat,nmos,ignd), &
             class_rot% THPROT  (nlat,nmos,ignd), &
             class_rot% THRROT  (nlat,nmos,ignd), &
             class_rot% THRAROT (nlat,nmos,ignd), &
             class_rot% ZBTWROT (nlat,nmos,ignd), &
             class_rot% THLWROT (nlat,nmos,ignd), &
             class_rot% GFLXROT (nlat,nmos,ignd), &
             class_rot% HMFGROT (nlat,nmos,ignd), &
             class_rot% HTCROT  (nlat,nmos,ignd), &
             class_rot% QFCROT  (nlat,nmos,ignd), &
             class_rot%TBARACC_M(nlat,nmos,ignd), &
             class_rot%THLQACC_M(nlat,nmos,ignd), &
             class_rot%THICACC_M(nlat,nmos,ignd), &
             class_rot%tbaraccrow_m(nlat,nmos,ignd))

    ! These will be allocated the dimension: 'nlat,nmos,ican'
    allocate(class_rot% ACIDROT (nlat,nmos,ican), &
             class_rot% ACVDROT (nlat,nmos,ican), &
             class_rot% CMASROT (nlat,nmos,ican), &
             class_rot% HGTDROT (nlat,nmos,ican), &
             class_rot% PAIDROT (nlat,nmos,ican), &
             class_rot% PAMNROT (nlat,nmos,ican), &
             class_rot% PAMXROT (nlat,nmos,ican), &
             class_rot% PSGAROT (nlat,nmos,ican), &
             class_rot% PSGBROT (nlat,nmos,ican), &
             class_rot% QA50ROT (nlat,nmos,ican), &
             class_rot% ROOTROT (nlat,nmos,ican), &
             class_rot% RSMNROT (nlat,nmos,ican), &
             class_rot% VPDAROT (nlat,nmos,ican), &
             class_rot% VPDBROT (nlat,nmos,ican))

    ! These will be allocated the dimension: 'nlat,nmos,icp1'
    allocate(class_rot% ALICROT (nlat,nmos,icp1), &
             class_rot% ALVCROT (nlat,nmos,icp1), &
             class_rot% FCANROT (nlat,nmos,icp1), &
             class_rot% LNZ0ROT (nlat,nmos,icp1))

    ! These will be allocated the dimension: 'nlat,nmos,nbs'
    allocate(class_rot% SALBROT  (nlat,nmos,nbs), &
            class_rot% CSALROT  (nlat,nmos,nbs))

    ! These will be allocated the dimension: 'nlat,nbs'
    allocate(class_rot% FSDBROL  (nlat,nbs), &
             class_rot% FSFBROL  (nlat,nbs), &
             class_rot% FSSBROL  (nlat,nbs))

    ! These will be allocated the dimension: 'nlat,ignd'

    allocate(class_rot% TBARROW (nlat,ignd), &
             class_rot% THALROW (nlat,ignd), &
             class_rot% THICROW (nlat,ignd), &
             class_rot% THLQROW (nlat,ignd), &
             class_rot% GFLXROW (nlat,ignd), &
             class_rot% HMFGROW (nlat,ignd), &
             class_rot% HTCROW  (nlat,ignd), &
             class_rot% QFCROW  (nlat,ignd))

    ! These will be allocated the dimension: 'nlat,nmos,6,50'
    allocate(class_rot% ITCTROT (nlat,nmos,6,50))

    ! These will be allocated the dimension: 'nlat,nmos,4'
    allocate(class_rot% TSFSROT (nlat,nmos,4))

  end subroutine allocClassVars
  !! @}
  !==================================================

  !> \ingroup classstatevars_resetClassMon
  !! @{
  !> Resets the CLASS (physics) monthly variables in preparation for the next month
  subroutine resetClassMon (nltest)

    use classicParams,    only : ignd

    implicit none

    integer, intent(in) :: nltest

    integer :: i, j

    do i = 1,nltest
      class_out%ALVSACC_MO(I) = 0.
      class_out%ALIRACC_MO(I) = 0.
      class_out%FLUTACC_MO(I) = 0.
      class_out%FSINACC_MO(I) = 0.
      class_out%FLINACC_MO(I) = 0.
      class_out%HFSACC_MO(I) = 0.
      class_out%QEVPACC_MO(I) = 0.
      class_out%groundHeatFlux_MO(I) = 0.
      class_out%TRANSPACC_MO(I) = 0.
      class_out%SNOACC_MO(I) = 0.
      class_out%WSNOACC_MO(I) = 0.
      class_out%ROFACC_MO(I) = 0.
      class_out%PREACC_MO(I) = 0.
      class_out%EVAPACC_MO(I) = 0.
      class_out%TAACC_MO(I) = 0.
      class_out%ACTLYR_MO(I) = 0.
      class_out%FTABLE_MO(I) = 0.
      class_out%ACTLYR_MIN_MO(I) = 100000.
      class_out%FTABLE_MIN_MO(I) = 100000.
      class_out%ACTLYR_MAX_MO(I) = 0.
      class_out%FTABLE_MAX_MO(I) = 0.
      class_out%CANOPYEVAP(I) = 0.
      class_out%GROUNDEVAP(I) = 0.
      class_out%ALTOTACC_MO(I) = 0.
      class_out%ALSNOACC_MO(I) = 0.
      class_out%altotcntr_m(i) = 0
      class_out%MRSO_MO(i) = 0.

      do J = 1,IGND
        class_out%TBARACC_MO(I,J) = 0.
        class_out%THLQACC_MO(I,J) = 0.
        class_out%THICACC_MO(I,J) = 0.
        class_out%MRSOL_MO(i,j) = 0.
      end do
    end do

  end subroutine resetClassMon
  !! @}
  !==================================================

  !> \ingroup classstatevars_resetClassYr
  !! @{
  !> Resets the CLASS (physics) annual variables in preparation for the next year
  subroutine resetClassYr (nltest)

    implicit none

    integer, intent(in) :: nltest

    integer :: i

    do i = 1,nltest
      class_out%ALVSACC_YR(I) = 0.
      class_out%ALIRACC_YR(I) = 0.
      class_out%FLUTACC_YR(I) = 0.
      class_out%FSINACC_YR(I) = 0.
      class_out%FLINACC_YR(I) = 0.
      class_out%HFSACC_YR(I) = 0.
      class_out%QEVPACC_YR(I) = 0.
      class_out%ROFACC_YR(I) = 0.
      class_out%PREACC_YR(I) = 0.
      class_out%EVAPACC_YR(I) = 0.
      class_out%TRANSPACC_YR(I) = 0.
      class_out%TAACC_YR(I) = 0.
      class_out%ALTOTACC_YR(I) = 0.
      class_out%altotcntr_yr(i) = 0
      class_out%ACTLYR_yr(I) = 0.
      class_out%FTABLE_yr(I) = 0.
      class_out%ACTLYR_MIN_yr(I) = 100000.
      class_out%FTABLE_MIN_yr(I) = 100000.
      class_out%ACTLYR_MAX_yr(I) = 0.
      class_out%FTABLE_MAX_yr(I) = 0.

    end do

  end subroutine resetClassYr
  !! @}
  !==================================================

  !> \ingroup classstatevars_resetAccVars
  !! @{
  !> Resets the CLASS (physics) aggregation variables in preparation for the next period
  subroutine resetAccVars (nltest,nmtest)

    use classicParams,         only : ignd

    implicit none

    integer, intent(in) :: nltest
    integer, intent(in) :: nmtest

    integer :: i, m, j

    do I = 1,NLTEST

      class_rot%altotcntr_d(i) = 0

      do M = 1,NMTEST

        class_rot%PREACC_M(i,m) = 0.
        class_rot%GTACC_M(i,m) = 0.
        class_rot%QEVPACC_M(i,m) = 0.
        class_rot%HFSACC_M(i,m) = 0.
        class_rot%HMFNACC_M(i,m) = 0.
        class_rot%ROFACC_M(i,m) = 0.
        class_rot%SNOACC_M(i,m) = 0.
        class_rot%OVRACC_M(i,m) = 0.
        class_rot%WTBLACC_M(i,m) = 0.
        class_rot%ALVSACC_M(i,m) = 0.
        class_rot%ALIRACC_M(i,m) = 0.
        class_rot%RHOSACC_M(i,m) = 0.
        class_rot%TSNOACC_M(i,m) = 0.
        class_rot%WSNOACC_M(i,m) = 0.
        class_rot%TCANACC_M(i,m) = 0.
        class_rot%RCANACC_M(i,m) = 0.
        class_rot%SCANACC_M(i,m) = 0.
        class_rot%GROACC_M(i,m) = 0.
        class_rot%FSINACC_M(i,m) = 0.
        class_rot%FLINACC_M(i,m) = 0.
        class_rot%TAACC_M(i,m) = 0.
        class_rot%UVACC_M(i,m) = 0.
        class_rot%PRESACC_M(i,m) = 0.
        class_rot%QAACC_M(i,m) = 0.
        class_rot%ALTOTACC_M(i,m) = 0.
        class_rot%ALSNOACC_M(i,m) = 0.
        class_rot%EVAPACC_M(i,m) = 0.
        class_rot%FLUTACC_M(i,m) = 0.

        do J = 1,IGND
          class_rot%TBARACC_M(I,M,J) = 0.
          class_rot%THLQACC_M(I,M,J) = 0.
          class_rot%THICACC_M(I,M,J) = 0.
          class_rot%tbaraccrow_m(i,m,j)  = 0.0

        end do
      end do
    end do

  end subroutine resetAccVars
  !! @}
  !==================================================
  !> \ingroup classstatevars_initDiagnosticVars
  !! @{
  !> Initialization of diagnostic variables split out of classGather for consistency with gcm applications.
  subroutine initDiagnosticVars (nml, ilg)

    use classicParams,      only : ignd

    implicit none

    integer, intent(in) :: nml, ilg

    integer :: k, m, l

    !    * INITIALIZATION OF DIAGNOSTIC VARIABLES SPLIT OUT OF classGather
    !    * FOR CONSISTENCY WITH GCM APPLICATIONS.

    do K = 1,ILG ! loop 330
      class_gat%CDHGAT (K) = 0.0
      class_gat%CDMGAT (K) = 0.0
      class_gat%HFSGAT (K) = 0.0
      class_gat%TFXGAT (K) = 0.0
      class_gat%QEVPGAT(K) = 0.0
      class_gat%QFSGAT (K) = 0.0
      class_gat%QFXGAT (K) = 0.0
      class_gat%PETGAT (K) = 0.0
      class_gat%GAGAT  (K) = 0.0
      class_gat%EFGAT  (K) = 0.0
      class_gat%GTGAT  (K) = 0.0
      class_gat%QGGAT  (K) = 0.0
      class_gat%ALVSGAT(K) = 0.0
      class_gat%ALIRGAT(K) = 0.0
      class_gat%SFCTGAT(K) = 0.0
      class_gat%SFCUGAT(K) = 0.0
      class_gat%SFCVGAT(K) = 0.0
      class_gat%SFCQGAT(K) = 0.0
      class_gat%FSNOGAT(K) = 0.0
      class_gat%FSGVGAT(K) = 0.0
      class_gat%FSGSGAT(K) = 0.0
      class_gat%FSGGGAT(K) = 0.0
      class_gat%FLGVGAT(K) = 0.0
      class_gat%FLGSGAT(K) = 0.0
      class_gat%FLGGGAT(K) = 0.0
      class_gat%HFSCGAT(K) = 0.0
      class_gat%HFSSGAT(K) = 0.0
      class_gat%HFSGGAT(K) = 0.0
      class_gat%HEVCGAT(K) = 0.0
      class_gat%HEVSGAT(K) = 0.0
      class_gat%HEVGGAT(K) = 0.0
      class_gat%HMFCGAT(K) = 0.0
      class_gat%HMFNGAT(K) = 0.0
      class_gat%HTCCGAT(K) = 0.0
      class_gat%HTCSGAT(K) = 0.0
      class_gat%PCFCGAT(K) = 0.0
      class_gat%PCLCGAT(K) = 0.0
      class_gat%PCPNGAT(K) = 0.0
      class_gat%PCPGGAT(K) = 0.0
      class_gat%QFGGAT (K) = 0.0
      class_gat%QFNGAT (K) = 0.0
      class_gat%QFCFGAT(K) = 0.0
      class_gat%QFCLGAT(K) = 0.0
      class_gat%ROFGAT (K) = 0.0
      class_gat%ROFOGAT(K) = 0.0
      class_gat%ROFBGAT(K) = 0.0
      class_gat%TROFGAT(K) = 0.0
      class_gat%TROOGAT(K) = 0.0
      class_gat%TROBGAT(K) = 0.0
      class_gat%ROFCGAT(K) = 0.0
      class_gat%ROFNGAT(K) = 0.0
      class_gat%ROVGGAT(K) = 0.0
      class_gat%WTRCGAT(K) = 0.0
      class_gat%WTRSGAT(K) = 0.0
      class_gat%WTRGGAT(K) = 0.0
      class_gat%DRGAT  (K) = 0.0
      class_gat%wtableGAT(K) = 0.0
    end do ! loop 330

    do L = 1,IGND ! loop 334
      do K = 1,ILG ! loop 332
        class_gat%HMFGGAT(K,L) = 0.0
        class_gat%HTCGAT (K,L) = 0.0
        class_gat%QFCGAT (K,L) = 0.0
        class_gat%GFLXGAT(K,L) = 0.0
        class_gat%ROFSGAT(K,L) = 0.0
        class_gat%TROSGAT(K,L) = 0.0
      end do ! loop 332
    end do ! loop 334

    do M = 1,50 ! loop 340
      do L = 1,6 ! loop 338
        do K = 1,NML ! loop 336
          class_gat%ITCTGAT(K,L,M) = 0
        end do ! loop 336
      end do ! loop 338
    end do ! loop 340

  end subroutine initDiagnosticVars
  !! @}

  !> \ingroup classstatevars_initRowVars
  !! @{
  !> Initialization of diagnostic variables split out of classGather for consistency with gcm applications.
  subroutine initRowVars (nml)

    use classicParams,      only : ignd

    implicit none

    integer, intent(in) :: nml

    integer :: j, i

    do I = 1,nml ! loop 525
      class_rot%CDHROW(I) = 0.
      class_rot%CDMROW(I) = 0.
      class_rot%HFSROW(I) = 0.
      class_rot%TFXROW(I) = 0.
      class_rot%QEVPROW(I) = 0.
      class_rot%QFSROW(I) = 0.
      class_rot%QFXROW(I) = 0.
      class_rot%PETROW(I) = 0.
      class_rot%GAROW(I) = 0.
      class_rot%EFROW(I) = 0.
      class_rot%GTROW(I) = 0.
      class_rot%QGROW(I) = 0.
      class_rot%ALVSROW(I) = 0.
      class_rot%ALIRROW(I) = 0.
      class_rot%SFCTROW(I) = 0.
      class_rot%SFCUROW(I) = 0.
      class_rot%SFCVROW(I) = 0.
      class_rot%SFCQROW(I) = 0.
      class_rot%SFRHROW(I) = 0.
      class_rot%SNOROW(I) = 0.
      class_rot%FSNOROW(I) = 0.
      class_rot%FSGVROW(I) = 0.
      class_rot%FSGSROW(I) = 0.
      class_rot%FSGGROW(I) = 0.
      class_rot%FLGVROW(I) = 0.
      class_rot%FLGSROW(I) = 0.
      class_rot%FLGGROW(I) = 0.
      class_rot%HFSCROW(I) = 0.
      class_rot%HFSSROW(I) = 0.
      class_rot%HFSGROW(I) = 0.
      class_rot%HEVCROW(I) = 0.
      class_rot%HEVSROW(I) = 0.
      class_rot%HEVGROW(I) = 0.
      class_rot%HMFCROW(I) = 0.
      class_rot%HMFNROW(I) = 0.
      class_rot%HTCCROW(I) = 0.
      class_rot%HTCSROW(I) = 0.
      class_rot%PCFCROW(I) = 0.
      class_rot%PCLCROW(I) = 0.
      class_rot%PCPNROW(I) = 0.
      class_rot%PCPGROW(I) = 0.
      class_rot%QFGROW(I) = 0.
      class_rot%QFNROW(I) = 0.
      class_rot%QFCLROW(I) = 0.
      class_rot%QFCFROW(I) = 0.
      class_rot%ROFROW(I) = 0.
      class_rot%ROFOROW(I) = 0.
      class_rot%ROFSROW(I) = 0.
      class_rot%ROFBROW(I) = 0.
      class_rot%ROFCROW(I) = 0.
      class_rot%ROFNROW(I) = 0.
      class_rot%ROVGROW(I) = 0.
      class_rot%RHOSROW(I) = 0.
      class_rot%WTRCROW(I) = 0.
      class_rot%WTRSROW(I) = 0.
      class_rot%WTRGROW(I) = 0.
      class_rot%DRROW(I) = 0.
      class_rot%TCANROW(I) = 0.
      class_rot%SCANROW(I) = 0.
      class_rot%RCANROW(I) = 0.
      class_rot%TSNOROW(I) = 0.
      class_rot%WSNOROW(I) = 0.
      class_rot%TPNDROW(I) = 0.
      class_rot%ZPNDROW(I) = 0.

      class_rot%wtableROW(I) = 0. ! FLAG

      class_rot%ILMOROW(I) = 0.
      class_rot%UEROW(I) = 0.
      class_rot%HBLROW(I) = 0.

      ! G12GRD(I)= 0.       ! YW March 27, 2015
      ! G23GRD(I)= 0.       ! YW March 27, 2015
      do J = 1,IGND ! loop 500
        class_rot%HMFGROW(I,J) = 0.
        class_rot%HTCROW(I,J) = 0.
        class_rot%QFCROW(I,J) = 0.
        class_rot%GFLXROW(I,J) = 0.
        class_rot%TBARROW(I,J) = 0.
        class_rot%THALROW(I,J) = 0.
        class_rot%THICROW(I,J) = 0.
        class_rot%THLQROW(I,J) = 0.
      end do ! loop 500
    end do ! loop 525

    ! Initialize to 0 for the start of a run.
    class_rot%actLyrThisYrROT(:,:) = 0.

  end subroutine initRowVars
  !! @}

  !> \namespace classstatevars
  !!
  !! Contains the physics variable type structures.
end module classStateVars
