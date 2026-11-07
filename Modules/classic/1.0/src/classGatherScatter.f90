!> \file
!> Transfers information between the 'gathered' and 'scattered' form of the CLASS data arrays.
!!
module classGatherScatter

  implicit none

  public :: classGather
  public :: classGatherPrep
  public :: classScatter

contains

  !> \ingroup classgatherscatter_classGather
  !! @{
  !> Gathers variables from two-dimensional arrays (latitude
  !! circle x mosaic tiles) onto long vectors for optimum processing
  !! efficiency on vector supercomputers.
  !> @author D. Verseghy, M. Lazare
  subroutine classGather (TBARGAT, THLQGAT, THICGAT, TPNDGAT, ZPNDGAT, & ! Formerly CLASSG
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
                          NML, NL, NT, NM, ILG, IG, IC, ICP1, NBS, &
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
                          VMODL, ZSNLROT, ZPLGROT, ZPLSROT, TACROT, &
                          QACROT, DRNROT, XSLPROT, GRKFROT, WFSFROT, &
                          WFCIROT, ALGWVROT, ALGWNROT, ALGDVROT, &
                          ALGDNROT, ASVDROT, ASIDROT, AGVDROT, &
                          AGIDROT, ISNDROT, RADJ, ZBLDROW, Z0ORROW, &
                          ZRFMROW, ZRFHROW, ZDMROW, ZDHROW, FSVHROW, &
                          FSIHROW, FSDBROL, FSFBROL, FSSBROL, CSZROW, &
                          FSGROL, FLGROL, FDLROL, ULROW, VLROW, &
                          TAROW, QAROW, PRESROW, PREROW, PADRROW, &
                          VPDROW, TADPROW, RHOAROW, RPCPROW, TRPCROW, &
                          SPCPROW, TSPCROW, RHSIROW, FCLOROW, DLONROW, &
                          GGEOROW, GUSTROL, REFROT, BCSNROT, DEPBROW, &
                          DLATROW, maxAnnualActLyrROT)

    !
    !     * DEC 23/16 - M.LAZARE.  PROMOTE DIMENSIONS OF WSNOROT, ASVDROT,
    !     *                        ASIDROT TO NLAT, NMOS (FOR LAKE MODEL)
    !     * Jan 16, 2015 - M.Lazare. New version called by "sfcproc3":
    !     *                          - Add THLW.
    !     *                          - {ALGWV, ALGWN, ALGDV, ALGDN} replace
    !     *                            {ALGW, ALGD}.
    !     *                          - FSG, FLG, GUST added.
    !     *                          - FDLROW changed to FDLROL (cosmetic).
    !     *                          - Adds GTGAT/GTROT.
    !     *                          - Adds NT (NTLD in sfcproc2
    !     *                            call) to dimension land-only
    !     *                            ROT fields, consistent with
    !     *                            new comrow12.
    !     *                          - Unused IWMOS, JWMOS removed.
    !     * Jun 13, 2013 - M.Lazare. CLASS gather routine called by
    !     *                          "sfcproc" in new version gcm17.
    !     * NOTE: This contains the following changes compared to the
    !     *       working temporary version used in conjunction with
    !     *       updates to gcm16 (ie not official):
    !     *         1) {DEPB, REF, BCSN} added for Maryam's new code.
    !     *         2) {FSDB, FSFB, FSSB} added for Jason's new code.
    !     * OCT 18/11 - M.LAZARE.  ADD IGDR.
    !     * OCT 07/11 - M.LAZARE.  ADD VMODL->VMODGAT.
    !     * OCT 05/11 - M.LAZARE.  PUT BACK IN PRESGROW->PRESGAT
    !     *                        REQUIRED FOR ADDED SURFACE RH
    !     *                        CALCULATION.
    !     * OCT 03/11 - M.LAZARE.  REMOVE ALL INITIALIZATION TO
    !     *                        ZERO OF GAT ARRAYS (NOW DONE
    !     *                        IN CLASS DRIVER).
    !     * SEP 16/11 - M.LAZARE.  - ROW->ROT AND GRD->ROW.
    !     *                        - REMOVE INITIALIZATION OF
    !     *                          {ALVS, ALIR} TO ZERO.
    !     *                        - REMOVE PRESGROW->PRESGAT
    !     *                          (OCEAN-ONLY NOW).
    !     *                        - RADJROW (64-BIT) NOW RADJ
    !     *                          (32-BIT).
    !     * MAR 23/06 - D.VERSEGHY. ADD WSNO, FSNO, GGEO.
    !     * MAR 18/05 - D.VERSEGHY. ADDITIONAL VARIABLES.
    !     * FEB 18/05 - D.VERSEGHY. ADD "TSFS" VARIABLES.
    !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
    !     * AUG 15/02 - D.VERSEGHY. GATHER OPERATION ON CLASS
    !     *                         VARIABLES.
    !
    implicit none
    !
    !     (Suffix GAT refers to variables on gathered long vectors; suffix
    !     ROT refers to variables on original two-dimensional arrays.)
    !
    !     * INTEGER CONSTANTS.
    !
    integer, intent(in)  :: NML, NL, NM, NT, ILG, IG, IC, ICP1, NBS
    integer  :: K, L, M
    !
    !     * LAND SURFACE PROGNOSTIC VARIABLES.
    !
    real, intent(in) :: TBARROT(NL,NT,IG)    !< Temperature of soil layers [K]
    real, intent(in) :: THLQROT(NL,NT,IG)    !< Volumetric liquid water content of soil
    !! layers \f$[m^3 m^{-3}]\f$
    real, intent(in) :: THICROT(NL,NT,IG)    !< Frozen water content of soil layers
    !! under vegetation \f$[m^3 m^{-3}]\f$
    real, intent(in) :: TPNDROT(NL,NT)   !< Temperature of ponded water [K]
    real, intent(in) :: ZPNDROT(NL,NT)   !< Depth of ponded water on surface [m]
    real, intent(in) :: TBASROT(NL,NT)   !< Temperature of bedrock in third soil layer [K]
    real, intent(in) :: ALBSROT(NL,NM)   !< Snow albedo [ ]
    real, intent(in) :: TSNOROT(NL,NM)   !< Snowpack temperature [K]
    real, intent(in) :: RHOSROT(NL,NM)   !< Density of snow \f$[kg m^{-3}]\f$
    real, intent(in) :: SNOROT (NL,NM)   !< Mass of snow pack \f$[kg m^{-2}]\f$
    real, intent(in) :: TCANROT(NL,NT)   !< Vegetation canopy temperature [K]
    real, intent(in) :: RCANROT(NL,NT)   !< Intercepted liquid water stored on canopy \f$[kg m^{-2}]\f$
    real, intent(in) :: SCANROT(NL,NT)   !< Intercepted frozen water stored on canopy \f$[kg m^{-2}]\f$
    real, intent(in) :: GROROT (NL,NT)   !< Vegetation growth index [ ]
    real, intent(in) :: CMAIROT(NL,NT)   !< Aggregated mass of vegetation canopy \f$[kg m^{-2}]\f$
    real, intent(in) :: TSFSROT(NL,NT,4) !< Ground surface temperature over subarea [K]
    real, intent(in) :: TACROT (NL,NT)   !< Temperature of air within vegetation canopy [K]
    real, intent(in) :: QACROT (NL,NT)   !< Specific humidity of air within vegetation
    !! canopy \f$[kg kg^{-1}]\f$
    real, intent(in) :: WSNOROT(NL,NM)   !< Liquid water content of snow pack \f$[kg m^{-2}]\f$
    real, intent(in) :: REFROT(NL,NM)    !< Snow grain size (for ISNOALB=1 option)
    real, intent(in) :: BCSNROT(NL,NM)   !< Black carbon mixing ratio (for ISNOALB=1 option)  [kg m-3]
    real, intent(in) :: maxAnnualActLyrROT(NL,NT)  !< Active layer depth maximum over the e-folding period specified by parameter eftime (m).
    !
    real, intent(out)   :: TBARGAT(ILG,IG), THLQGAT(ILG,IG), THICGAT(ILG,IG), &
                           TPNDGAT(ILG), ZPNDGAT(ILG), TBASGAT(ILG), &
                           ALBSGAT(ILG), TSNOGAT(ILG), RHOSGAT(ILG), &
                           SNOGAT(ILG), TCANGAT(ILG), RCANGAT(ILG), &
                           SCANGAT(ILG), GROGAT(ILG), CMAIGAT(ILG), &
                           TSFSGAT(ILG,4), TACGAT(ILG), QACGAT(ILG), &
                           WSNOGAT(ILG), REFGAT(ILG), BCSNGAT(ILG), &
                           maxAnnualActLyrGAT(ILG)
    !
    !     * GATHER-SCATTER INDEX ARRAYS.
    !
    integer, intent(in)  :: ILMOS (ILG)  !< Index of latitude grid cell corresponding
    !! to current element of gathered vector of
    !! land surface variables [ ]
    integer, intent(in)  :: JLMOS (ILG)  !< Index of mosaic tile corresponding to
    !! current element of gathered vector of land
    !! surface variables [ ]

    !
    !     * CANOPY AND SOIL INFORMATION ARRAYS.
    !     * (THE LENGTH OF THESE ARRAYS IS DETERMINED BY THE NUMBER
    !     * OF SOIL LAYERS (3) AND THE NUMBER OF BROAD VEGETATION
    !     * CATEGORIES (4, OR 5 INCLUDING URBAN AREAS).)
    !
    real, intent(in) :: FCANROT(NL,NT,ICP1)  !< Maximum fractional coverage of modelled
    !! area by vegetation category [ ]
    real, intent(in) :: LNZ0ROT(NL,NT,ICP1)  !< Natural logarithm of maximum roughness
    !! length of vegetation category [ ]
    real, intent(in) :: ALVCROT(NL,NT,ICP1)  !< Background average visible albedo of
    !! vegetation category [ ]
    real, intent(in) :: ALICROT(NL,NT,ICP1)  !< Background average near-infrared albedo
    !! of vegetation category [ ]
    real, intent(in) :: PAMXROT(NL,NT,IC)    !< Maximum plant area index of vegetation
    !! category [ ]
    real, intent(in) :: PAMNROT(NL,NT,IC)    !< Minimum plant area index of vegetation
    !! category [ ]
    real, intent(in) :: CMASROT(NL,NT,IC)    !< Maximum canopy mass for vegetation
    !! category \f$[kg m^{-2}]\f$
    real, intent(in) :: ROOTROT(NL,NT,IC)    !< Maximum rooting depth of vegetation
    !! category [m]
    real, intent(in) :: RSMNROT(NL,NT,IC)    !< Minimum stomatal resistance of
    !! vegetation category \f$[s m^{-1}]\f$
    real, intent(in) :: QA50ROT(NL,NT,IC)    !< Reference value of incoming shortwave
    !! radiation for vegetation category (used
    !! in stomatal resistance calculation) \f$[W m^{-2}]\f$
    real, intent(in) :: VPDAROT(NL,NT,IC)    !< Vapour pressure deficit coefficient for
    !! vegetation category (used in stomatal
    !! resistance calculation) [ ]
    real, intent(in) :: VPDBROT(NL,NT,IC)    !< Vapour pressure deficit coefficient for
    !! vegetation category (used in stomatal
    !! resistance calculation) [ ]
    real, intent(in) :: PSGAROT(NL,NT,IC)    !< Soil moisture suction coefficient for
    !! vegetation category (used in stomatal
    !! resistance calculation) [ ]
    real, intent(in) :: PSGBROT(NL,NT,IC)    !< Soil moisture suction coefficient for
    !! vegetation category (used in stomatal
    !! resistance calculation) [ ]
    real, intent(in) :: PAIDROT(NL,NT,IC)    !< Optional user-specified value of plant
    !! area indices of vegetation categories
    !! to override CLASS-calculated values [ ]
    real, intent(in) :: HGTDROT(NL,NT,IC)    !< Optional user-specified values of
    !! height of vegetation categories to
    !! override CLASS-calculated values [m]
    real, intent(in) :: ACVDROT(NL,NT,IC)    !< Optional user-specified value of canopy
    !! visible albedo to override CLASS-
    !! calculated value [ ]
    real, intent(in) :: ACIDROT(NL,NT,IC)    !< Optional user-specified value of canopy
    !! near-infrared albedo to override CLASS-
    !! calculated value [ ]
    !
    real, intent(out) :: FCANGAT(ILG,ICP1), LNZ0GAT(ILG,ICP1), &
                         ALVCGAT(ILG,ICP1), ALICGAT(ILG,ICP1), &
                         PAMXGAT(ILG,IC), PAMNGAT(ILG,IC), &
                         CMASGAT(ILG,IC), ROOTGAT(ILG,IC), &
                         RSMNGAT(ILG,IC), QA50GAT(ILG,IC), &
                         VPDAGAT(ILG,IC), VPDBGAT(ILG,IC), &
                         PSGAGAT(ILG,IC), PSGBGAT(ILG,IC), &
                         PAIDGAT(ILG,IC), HGTDGAT(ILG,IC), &
                         ACVDGAT(ILG,IC), ACIDGAT(ILG,IC)
    !
    real, intent(in) :: THPROT (NL,NT,IG)    !< Pore volume in soil layer \f$[m^3 m^{-3}]\f$
    real, intent(in) :: THRROT (NL,NT,IG)    !< Liquid water retention capacity for
    !! organic soil \f$[m^3 m^{-3}]\f$
    real, intent(in) :: THMROT (NL,NT,IG)    !< Residual soil liquid water content
    !! remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
    real, intent(in) :: BIROT  (NL,NT,IG)    !< Clapp and Hornberger empirical "b"
    !! parameter [ ]
    real, intent(in) :: PSISROT(NL,NT,IG)    !< Soil moisture suction at saturation [m]
    real, intent(in) :: GRKSROT(NL,NT,IG)    !< Saturated hydraulic conductivity of
    !! soil layers \f$[m s^{-1}]\f$
    real, intent(in) :: THRAROT(NL,NT,IG)    !< Fractional saturation of soil behind
    !! the wetting front [ ]
    real, intent(in) :: HCPSROT(NL,NT,IG)    !< Volumetric heat capacity of soil
    !! particles \f$[J m^{-3}]\f$
    real, intent(in) :: TCSROT (NL,NT,IG)    !< Thermal conductivity of soil particles
    !! \f$[W m^{-1} K^{-1}]\f$
    real, intent(in) :: THFCROT(NL,NT,IG)    !< Field capacity \f$[m^3 m^{-3}]\f$
    real, intent(in) :: THLWROT(NL,NT,IG)    !< Liquid water content at wilting point \f$[m^3 m^{-3}]\f$
    real, intent(in) :: PSIWROT(NL,NT,IG)    !< Soil moisture suction at wilting point
    !! [m]
    real, intent(in) :: DLZWROT(NL,NT,IG)    !< Permeable thickness of soil layer [m]
    real, intent(in) :: ZBTWROT(NL,NT,IG)    !< Depth to permeable bottom of soil layer [m]
    real, intent(in) :: DRNROT (NL,NT)       !< Drainage index at bottom of soil profile [ ]
    real, intent(in) :: XSLPROT(NL,NT)       !< Surface slope (used when running MESH code) [degrees]
    real, intent(in) :: GRKFROT(NL,NT)       !< WATROF parameter used when running MESH code [ ]
    real, intent(in) :: WFSFROT(NL,NT)       !< WATROF parameter used when running MESH code [ ]
    real, intent(in) :: WFCIROT(NL,NT)       !< WATROF parameter used when running MESH code [ ]
    real, intent(in) :: ALGWVROT(NL,NT)      !< Reference visible albedo for saturated soil  [  ]
    real, intent(in) :: ALGWNROT(NL,NT)      !< Reference near-infrared albedo for saturated soil  [  ]
    real, intent(in) :: ALGDVROT(NL,NT)      !< Reference visible albedo for dry soil  [  ]
    real, intent(in) :: ALGDNROT(NL,NT)      !< Reference near-infrared albedo for dry soil  [  ]
    real, intent(in) :: ASVDROT(NL,NT)       !< Optional user-specified value of snow
    !! visible albedo to override CLASS-
    !! calculated value [ ]
    real, intent(in) :: ASIDROT(NL,NM)       !< Optional user-specified value of snow
    !! near-infrared albedo to override CLASS-
    !! calculated value [ ]
    real, intent(in) :: AGVDROT(NL,NT)       !< Optional user-specified value of ground
    !! visible albedo to override CLASS-
    !! calculated value [ ]
    real, intent(in) :: AGIDROT(NL,NT)       !< Optional user-specified value of ground
    !! near-infrared albedo to override CLASS-
    !! calculated value [ ]
    real, intent(in) :: ZSNLROT(NL,NT)       !< Limiting snow depth below which
    !! coverage is < 100% [m]
    real, intent(in) :: ZPLGROT(NL,NT)       !< Maximum water ponding depth for snow-
    !! free subareas (user-specified when
    !! running MESH code) [m]
    real, intent(in) :: ZPLSROT(NL,NT)       !< Maximum water ponding depth for snow-
    !! covered subareas (user-specified when
    !! running MESH code) [m]
    !

    real, intent(out) :: THPGAT (ILG,IG), THRGAT (ILG,IG), THMGAT (ILG,IG), &
                         BIGAT  (ILG,IG), PSISGAT(ILG,IG), GRKSGAT(ILG,IG), &
                         THRAGAT(ILG,IG), HCPSGAT(ILG,IG), &
                         TCSGAT (ILG,IG), THFCGAT(ILG,IG), THLWGAT(ILG,IG), &
                         PSIWGAT(ILG,IG), DLZWGAT(ILG,IG), ZBTWGAT(ILG,IG), &
                         DRNGAT (ILG), XSLPGAT(ILG), GRKFGAT(ILG), &
                         WFSFGAT(ILG), WFCIGAT(ILG), ALGWVGAT(ILG), &
                         ALGWNGAT(ILG), ALGDVGAT(ILG), ALGDNGAT(ILG), &
                         ASVDGAT(ILG), ASIDGAT(ILG), &
                         AGVDGAT(ILG), AGIDGAT(ILG), ZSNLGAT(ILG), &
                         ZPLGGAT(ILG), ZPLSGAT(ILG)
    !
    integer, intent(in)  :: ISNDROT(NL,NT,IG)
    integer, intent(out) :: ISNDGAT(ILG,IG) !< Sand content flag
    integer, intent(in)  :: IGDRROT(NL,NT)
    integer, intent(out) :: IGDRGAT(ILG) !< Index of soil layer in
    !! which bedrock is
    !! encountered
    !
    !     * ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES.
    !
    real, intent(in) :: ZRFMROW(NL) !< Reference height associated with forcing wind
    !! speed [m]
    real, intent(in) :: ZRFHROW(NL) !< Reference height associated with forcing air
    !! temperature and humidity [m]
    real, intent(in) :: ZDMROW (NL) !< User-specified height associated with diagnosed
    !! anemometer-level wind speed [m]
    real, intent(in) :: ZDHROW (NL) !< User-specified height associated with diagnosed
    !! screen-level variables [m]
    real, intent(in) :: FSVHROW(NL) !< Visible radiation incident on horizontal
    !! surface \f$[W m^{-2}]\f$
    real, intent(in) :: FSIHROW(NL) !< Near-infrared radiation incident on horizontal
    !! surface \f$[W m^{-2}]\f$
    real, intent(in) :: CSZROW (NL) !< Cosine of solar zenith angle [ ]
    real, intent(in) :: FSGROL (NL) !< Total shortwave radiation absorbed by land surface  [W m-2]
    real, intent(in) :: FLGROL (NL) !< Total longwave radiation absorbed by land surface  [W m-2]
    real, intent(in) :: FDLROL (NL) !< Downwelling longwave radiation at bottom of
    !! atmosphere \f$[W m^{-2}]\f$

    real, intent(in) :: ULROW  (NL) !< Zonal component of wind speed \f$[m s^{-1}]\f$
    real, intent(in) :: VLROW  (NL) !< Meridional component of wind speed \f$[m s^{-1}]\f$
    real, intent(in) :: TAROW  (NL) !< Air temperature at reference height [K]
    real, intent(in) :: QAROW  (NL) !< Specific humidity at reference height \f$[kg kg^{-1}]\f$
    real, intent(in) :: PRESROW(NL) !< Surface air pressure [Pa]
    real, intent(in) :: PREROW (NL) !< Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$
    real, intent(in) :: PADRROW(NL) !< Partial pressure of dry air [Pa]
    real, intent(in) :: VPDROW (NL) !< Vapour pressure deficit [mb]
    real, intent(in) :: TADPROW(NL) !< Dew point temperature of air [K]
    real, intent(in) :: RHOAROW(NL) !< Density of air \f$[kg m^{-3}]\f$
    real, intent(in) :: ZBLDROW(NL) !< Atmospheric blending height for surface
    !! roughness length averaging [m]
    real, intent(in) :: Z0ORROW(NL) !< Orographic roughness length [m]
    real, intent(in) :: RPCPROW(NL) !< Rainfall rate over modelled area \f$[m s^{-1}]\f$
    real, intent(in) :: TRPCROW(NL) !< Rainfall temperature [K]
    real, intent(in) :: SPCPROW(NL) !< Snowfall rate over modelled area \f$[m s^{-1}]\f$
    real, intent(in) :: TSPCROW(NL) !< Snowfall temperature [K]
    real, intent(in) :: RHSIROW(NL) !< Density of fresh snow \f$[kg m^{-3}]\f$
    real, intent(in) :: FCLOROW(NL) !< Fractional cloud cover [ ]
    real, intent(in) :: DLONROW(NL) !< Longitude of grid cell (east of Greenwich) [degrees]
    real, intent(in) :: DLATROW(NL) !< Latitude of grid cell [degrees]
    real, intent(in) :: GGEOROW(NL) !< Geothermal heat flux at bottom of soil profile
    !! \f$[W m^{-2}]\f$
    real, intent(in) :: GUSTROL(NL) !< Wind gustiness factor  [  ]
    real, intent(in) :: RADJ   (NL) !< Latitude of grid cell (positive north of equator) [rad]
    real, intent(in) :: VMODL  (NL) !< Wind speed at reference height \f$[m s^{-1}]\f$
    real, intent(in) :: DEPBROW(NL) ! Black carbon deposition rate  [kg m-2 s-1]

    real, intent(in) ::  FSDBROL(NL,NBS) !< Direct solar radiation in each modelled wavelength band  [W m-2]
    real, intent(in) ::  FSFBROL(NL,NBS) !< Diffuse solar radiation in each modelled wavelength band  [W m-2]
    real, intent(in) ::  FSSBROL(NL,NBS) !< Total solar radiation in each modelled wavelength band  [W m-2]

    !
    real, intent(out) :: ZRFMGAT(ILG), ZRFHGAT(ILG), ZDMGAT (ILG), ZDHGAT (ILG), &
                         FSVHGAT(ILG), FSIHGAT(ILG), CSZGAT (ILG), DLATGAT(ILG), &
                         FSGGAT (ILG), FLGGAT (ILG), FDLGAT (ILG), &
                         ULGAT  (ILG), VLGAT  (ILG), TAGAT  (ILG), QAGAT  (ILG), &
                         PRESGAT(ILG), PREGAT (ILG), PADRGAT(ILG), VPDGAT (ILG), &
                         TADPGAT(ILG), RHOAGAT(ILG), ZBLDGAT(ILG), Z0ORGAT(ILG), &
                         RPCPGAT(ILG), TRPCGAT(ILG), SPCPGAT(ILG), TSPCGAT(ILG), &
                         RHSIGAT(ILG), FCLOGAT(ILG), DLONGAT(ILG), GGEOGAT(ILG), &
                         GUSTGAT(ILG), RADJGAT(ILG), VMODGAT(ILG), DEPBGAT(ILG)

    real, intent(out) :: FSDBGAT(ILG,NBS), FSFBGAT(ILG,NBS), FSSBGAT(ILG,NBS)
    !----------------------------------------------------------------------
    !
    ! The prognostic, background and input variables are gathered into
    ! long arrays (collapsing the latitude and mosaic dimensions into
    ! one, but retaining the soil level and canopy category dimensions)
    ! using the pointer vectors generated in classGatherPrep.
    !
    do K = 1,NML ! loop 100
      TPNDGAT(K) = TPNDROT(ILMOS(K),JLMOS(K))
      ZPNDGAT(K) = ZPNDROT(ILMOS(K),JLMOS(K))
      TBASGAT(K) = TBASROT(ILMOS(K),JLMOS(K))
      ALBSGAT(K) = ALBSROT(ILMOS(K),JLMOS(K))
      TSNOGAT(K) = TSNOROT(ILMOS(K),JLMOS(K))
      RHOSGAT(K) = RHOSROT(ILMOS(K),JLMOS(K))
      SNOGAT (K) = SNOROT (ILMOS(K),JLMOS(K))
      REFGAT (K) = REFROT (ILMOS(K),JLMOS(K))
      BCSNGAT(K) = BCSNROT(ILMOS(K),JLMOS(K))
      maxAnnualActLyrGAT(k) = maxAnnualActLyrROT(ILMOS(K),JLMOS(K))
      WSNOGAT(K) = WSNOROT(ILMOS(K),JLMOS(K))
      TCANGAT(K) = TCANROT(ILMOS(K),JLMOS(K))
      RCANGAT(K) = RCANROT(ILMOS(K),JLMOS(K))
      SCANGAT(K) = SCANROT(ILMOS(K),JLMOS(K))
      GROGAT (K) = GROROT (ILMOS(K),JLMOS(K))
      CMAIGAT(K) = CMAIROT(ILMOS(K),JLMOS(K))
      DRNGAT (K) = DRNROT (ILMOS(K),JLMOS(K))
      XSLPGAT(K) = XSLPROT(ILMOS(K),JLMOS(K))
      GRKFGAT(K) = GRKFROT(ILMOS(K),JLMOS(K))
      WFSFGAT(K) = WFSFROT(ILMOS(K),JLMOS(K))
      WFCIGAT(K) = WFCIROT(ILMOS(K),JLMOS(K))
      ALGWVGAT(K) = ALGWVROT(ILMOS(K),JLMOS(K))
      ALGWNGAT(K) = ALGWNROT(ILMOS(K),JLMOS(K))
      ALGDVGAT(K) = ALGDVROT(ILMOS(K),JLMOS(K))
      ALGDNGAT(K) = ALGDNROT(ILMOS(K),JLMOS(K))
      ASVDGAT(K) = ASVDROT(ILMOS(K),JLMOS(K))
      ASIDGAT(K) = ASIDROT(ILMOS(K),JLMOS(K))
      AGVDGAT(K) = AGVDROT(ILMOS(K),JLMOS(K))
      AGIDGAT(K) = AGIDROT(ILMOS(K),JLMOS(K))
      ZSNLGAT(K) = ZSNLROT(ILMOS(K),JLMOS(K))
      ZPLGGAT(K) = ZPLGROT(ILMOS(K),JLMOS(K))
      ZPLSGAT(K) = ZPLSROT(ILMOS(K),JLMOS(K))
      TACGAT (K) = TACROT (ILMOS(K),JLMOS(K))
      QACGAT (K) = QACROT (ILMOS(K),JLMOS(K))
      IGDRGAT(K) = IGDRROT(ILMOS(K),JLMOS(K))
      ZBLDGAT(K) = ZBLDROW(ILMOS(K))
      Z0ORGAT(K) = Z0ORROW(ILMOS(K))
      ZRFMGAT(K) = ZRFMROW(ILMOS(K))
      ZRFHGAT(K) = ZRFHROW(ILMOS(K))
      ZDMGAT (K) = ZDMROW(ILMOS(K))
      ZDHGAT (K) = ZDHROW(ILMOS(K))
      FSVHGAT(K) = FSVHROW(ILMOS(K))
      FSIHGAT(K) = FSIHROW(ILMOS(K))
      CSZGAT (K) = CSZROW (ILMOS(K))
      FSGGAT (K) = FSGROL (ILMOS(K))
      FLGGAT (K) = FLGROL (ILMOS(K))
      FDLGAT (K) = FDLROL (ILMOS(K))
      ULGAT  (K) = ULROW  (ILMOS(K))
      VLGAT  (K) = VLROW  (ILMOS(K))
      TAGAT  (K) = TAROW  (ILMOS(K))
      QAGAT  (K) = QAROW  (ILMOS(K))
      PRESGAT(K) = PRESROW(ILMOS(K))
      PREGAT (K) = PREROW (ILMOS(K))
      PADRGAT(K) = PADRROW(ILMOS(K))
      VPDGAT (K) = VPDROW (ILMOS(K))
      TADPGAT(K) = TADPROW(ILMOS(K))
      RHOAGAT(K) = RHOAROW(ILMOS(K))
      RPCPGAT(K) = RPCPROW(ILMOS(K))
      TRPCGAT(K) = TRPCROW(ILMOS(K))
      SPCPGAT(K) = SPCPROW(ILMOS(K))
      TSPCGAT(K) = TSPCROW(ILMOS(K))
      RHSIGAT(K) = RHSIROW(ILMOS(K))
      FCLOGAT(K) = FCLOROW(ILMOS(K))
      DLONGAT(K) = DLONROW(ILMOS(K))
      DLATGAT(K) = DLATROW(ILMOS(K))
      GGEOGAT(K) = GGEOROW(ILMOS(K))
      GUSTGAT(K) = GUSTROL(ILMOS(K))
      RADJGAT(K) = RADJ   (ILMOS(K))
      VMODGAT(K) = VMODL  (ILMOS(K))
      DEPBGAT(K) = DEPBROW(ILMOS(K))
    end do ! loop 100
    !
    do L = 1,IG ! loop 250
      do K = 1,NML ! loop 200
        TBARGAT(K,L) = TBARROT(ILMOS(K),JLMOS(K),L)
        THLQGAT(K,L) = THLQROT(ILMOS(K),JLMOS(K),L)
        THICGAT(K,L) = THICROT(ILMOS(K),JLMOS(K),L)
        THPGAT (K,L) = THPROT (ILMOS(K),JLMOS(K),L)
        THRGAT (K,L) = THRROT (ILMOS(K),JLMOS(K),L)
        THMGAT (K,L) = THMROT (ILMOS(K),JLMOS(K),L)
        BIGAT  (K,L) = BIROT  (ILMOS(K),JLMOS(K),L)
        PSISGAT(K,L) = PSISROT(ILMOS(K),JLMOS(K),L)
        GRKSGAT(K,L) = GRKSROT(ILMOS(K),JLMOS(K),L)
        THRAGAT(K,L) = THRAROT(ILMOS(K),JLMOS(K),L)
        HCPSGAT(K,L) = HCPSROT(ILMOS(K),JLMOS(K),L)
        TCSGAT (K,L) = TCSROT (ILMOS(K),JLMOS(K),L)
        THFCGAT(K,L) = THFCROT(ILMOS(K),JLMOS(K),L)
        THLWGAT(K,L) = THLWROT(ILMOS(K),JLMOS(K),L)
        PSIWGAT(K,L) = PSIWROT(ILMOS(K),JLMOS(K),L)
        DLZWGAT(K,L) = DLZWROT(ILMOS(K),JLMOS(K),L)
        ZBTWGAT(K,L) = ZBTWROT(ILMOS(K),JLMOS(K),L)
        ISNDGAT(K,L) = ISNDROT(ILMOS(K),JLMOS(K),L)
      end do ! loop 200
    end do ! loop 250
    !
    do L = 1,ICP1 ! loop 300
      do K = 1,NML
        FCANGAT(K,L) = FCANROT(ILMOS(K),JLMOS(K),L)
        LNZ0GAT(K,L) = LNZ0ROT(ILMOS(K),JLMOS(K),L)
        ALVCGAT(K,L) = ALVCROT(ILMOS(K),JLMOS(K),L)
        ALICGAT(K,L) = ALICROT(ILMOS(K),JLMOS(K),L)
      end do
    end do ! loop 300
    !
    do L = 1,IC ! loop 400
      do K = 1,NML
        PAMXGAT(K,L) = PAMXROT(ILMOS(K),JLMOS(K),L)
        PAMNGAT(K,L) = PAMNROT(ILMOS(K),JLMOS(K),L)
        CMASGAT(K,L) = CMASROT(ILMOS(K),JLMOS(K),L)
        ROOTGAT(K,L) = ROOTROT(ILMOS(K),JLMOS(K),L)
        RSMNGAT(K,L) = RSMNROT(ILMOS(K),JLMOS(K),L)
        QA50GAT(K,L) = QA50ROT(ILMOS(K),JLMOS(K),L)
        VPDAGAT(K,L) = VPDAROT(ILMOS(K),JLMOS(K),L)
        VPDBGAT(K,L) = VPDBROT(ILMOS(K),JLMOS(K),L)
        PSGAGAT(K,L) = PSGAROT(ILMOS(K),JLMOS(K),L)
        PSGBGAT(K,L) = PSGBROT(ILMOS(K),JLMOS(K),L)
        PAIDGAT(K,L) = PAIDROT(ILMOS(K),JLMOS(K),L)
        HGTDGAT(K,L) = HGTDROT(ILMOS(K),JLMOS(K),L)
        ACVDGAT(K,L) = ACVDROT(ILMOS(K),JLMOS(K),L)
        ACIDGAT(K,L) = ACIDROT(ILMOS(K),JLMOS(K),L)
      end do
    end do ! loop 400

    do L = 1,4 ! over the four subareas ! FLAG HACK - legit bugfix.
      do K = 1,NML
        TSFSGAT(K,L) = TSFSROT(ILMOS(K),JLMOS(K),L)
      end do
    end do
    !
    do L = 1,NBS
      do K = 1,NML
        FSDBGAT(K,L) = FSDBROL(ILMOS(K),L)
        FSFBGAT(K,L) = FSFBROL(ILMOS(K),L)
        FSSBGAT(K,L) = FSSBROL(ILMOS(K),L)
      end do ! K
    end do ! L
    return
  end subroutine classGather
  !! @}

  !> \ingroup classgatherscatter_classGatherPrep
  !! @{
  !> Assigns values to pointer vectors relating the location of elements on the "gathered" variable
  !! vectors to elements on the original two-dimensional arrays (latitude circle x mosaic tiles) for land grid
  !! cells.
  !> @author D. Verseghy, M. Lazare, E. Chan
  subroutine classGatherPrep (ILMOS, JLMOS, IWMOS, JWMOS, & ! Formerly GATPREP
                              NML, NMW, GCROW, FAREA, MOSID, &
                              NL, NM, ILG, IL1, IL2, IM)
    !
    !     * JAN 12/17 - D.VERSEGHY. NOTE: THIS VERSION OF classGatherPrep
    !     *                         IS DESIGNED SPECIFICALLY FOR LAND.
    !     *                         THE VERSION CURRENTLY USED IN THE AGCM
    !     *                         (SINCE JAN. 2014) HAS BEEN GENERALIZED
    !     *                         FOR LAND, LAKES AND WATER/ICE.
    !     * DEC 28/11 - D.VERSEGHY. CHANGE ILGM BACK TO ILG AND
    !     *                         ILG TO NL FOR CONSISTENCY WITH
    !     *                         BOTH STAND-ALONE AND GCM
    !     *                         CONVENTIONS.
    !     * OCT 22/11 - M.LAZARE. REMOVE OCEAN/ICE CODE (NOW DONE
    !     *                       IN COISS).
    !     * OCT 21/11 - M.LAZARE. COSMETIC: ILG->ILGM AND NLAT->ILG,
    !     *                       TO BE CONSISTENT WITH MODEL
    !     *                       CONVENTION. ALSO GCGRD->GCROW.
    !     * JUN 12/06 - E.CHAN.  DIMENSION IWAT AND IICE BY ILG.
    !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
    !     * AUG 09/02 - D.VERSEGHY/M.LAZARE. DETERMINE INDICES FOR
    !     *                        GATHER-SCATTER OPERATIONS ON
    !     *                        CURRENT LATITUDE LOOP.
    !
    implicit none
    !
    !     * INTEGER CONSTANTS.
    !
    integer, intent(inout)  :: NML  !< Total number of mosaic tiles in land surface gather vectors
    integer, intent(inout)  :: NMW  !< Total number of mosaic tiles in inland water gather vectors
    integer, intent(in)  :: NL   !< Number of latitude grid cells in land surface scatter vectors
    integer, intent(in)  :: NM   !< Total number of CLASS mosaic tiles in land surface gather vectors
    integer, intent(in)  :: ILG  !< Total number of mosaic tiles per latitude grid cell in land surface scatter vector
    integer, intent(in)  :: IL1,IL2  !<
    integer, intent(in)  :: IM   !< Maximum number of mosaic tiles within the grid cells in the array under consideration
    integer  :: I,J
    !
    !     * OUTPUT FIELDS.
    !
    integer, intent(out)  :: ILMOS  (ILG) !< Index of grid cell corresponding to current element
    !< of gathered vector of land surface variables [ ]
    integer, intent(out)  :: JLMOS  (ILG) !< Index of mosaic tile corresponding to current element
    !< of gathered vector of land surface variables [ ]
    integer, intent(out)  :: IWMOS  (ILG) !< Index of grid cell corresponding to current element of gathered vector
    !< of inland water body variables [ ]
    integer, intent(out)  :: JWMOS  (ILG) !< Index of mosaic tile corresponding to current element of gathered vector
    !< of inland water body variables [ ]
    !
    !     * INPUT FIELDS.
    !
    real     :: GCROW (NL)    !< Real number identifier indicating whether the grid cell
    !< is land (-1.0), sea ice (+1.0), or ocean (0.0)
    real     :: FAREA (NL,NM) !< Fractional coverage of mosaic tile on grid cell [ ]
    !
    integer  :: MOSID (NL,NM) !< Mosaic tile type identifier (1 for land, 0 for inland water)
    !---------------------------------------------------------------------
    NML = 0
    NMW = 0
    !>
    !! A looping operation is performed over the latitude circle, or array of grid cells, under consideration. If
    !! the grid cell is a land one (GCROW = -1.0), an additional internal loop is performed over all the mosaic
    !! tiles present. For each mosaic tile, if its fractional coverage is greater than zero, then if the mosaic type
    !! identifier MOSID is equal to 1 (indicating land), the counter of total mosaic tiles in the land surface gather
    !! vectors, NML, is incremented by one, and the elements of the vectors ILMOS and JLMOS corresponding
    !! to NML are set to the indices of the current grid cell and mosaic tile respectively. If MOSID is equal to
    !! zero (indicating inland water), the counter of total mosaic tiles in the inland water gather vectors, NMW,
    !! is incremented by one, and the elements of the vectors IWMOS and JWMOS corresponding to NMW are
    !! set to the indices of the current grid cell and mosaic tile respectively.
    !!
    do I = IL1,IL2 ! loop 200
      if (GCROW(I) <= - 0.5) then
        do J = 1,IM ! loop 100
          if (FAREA(I,J) > 0.0) then
            if (MOSID(I,J) > 0) then
              NML = NML + 1
              ILMOS(NML) = I
              JLMOS(NML) = J
            else
              NMW = NMW + 1
              IWMOS(NMW) = I
              JWMOS(NMW) = J
            end if
          end if
        end do ! loop 100
      end if
    end do ! loop 200

    return
  end subroutine classGatherPrep
  !! @}

  !> \ingroup classgatherscatter_classScatter
  !! @{
  !> Scatters variables from long, gathered vectors back onto original two-dimensional arrays (latitude
  !! circle x mosaic tiles). The suffix ROT refers to variables on original two-dimensional arrays.
  !! The suffix GAT refers to variables on gathered long vectors.
  !> @author D. Verseghy, M. Lazare
  subroutine classScatter (TBARROT, THLQROT, THICROT, TSFSROT, TPNDROT, & ! Formerly CLASSS
                           ZPNDROT, TBASROT, ALBSROT, TSNOROT, RHOSROT, &
                           SNOROT, GTROT, TCANROT, RCANROT, SCANROT, &
                           GROROT, CMAIROT, TACROT, QACROT, WSNOROT, &
                           REFROT, BCSNROT, EMISROT, SALBROT, CSALROT, &
                           groundHeatFluxROT, &
                           ILMOS, JLMOS, &
                           NML, NL, NT, NM, ILG, IG, IC, ICP1, NBS, &
                           TBARGAT, THLQGAT, THICGAT, TSFSGAT, TPNDGAT, &
                           ZPNDGAT, TBASGAT, ALBSGAT, TSNOGAT, RHOSGAT, &
                           SNOGAT, GTGAT, TCANGAT, RCANGAT, SCANGAT, &
                           GROGAT, CMAIGAT, TACGAT, QACGAT, WSNOGAT, &
                           REFGAT, BCSNGAT, EMISGAT, SALBGAT, CSALGAT,&
                           groundHeatFlux)
    !
    !     * DEC 23/16 - M.LAZARE.  PROMOTE DIMENSIONS OF WSNOROT TO
    !     *                        NLAT, NMOS (FOR LAKE MODEL)
    !     * Jun 20, 2014 - M.Lazare. New version for gcm18, called
    !     *                          by new "sfcproc2":
    !     *                          - Adds SALBGAT/SALBROT and
    !     *                            CSALGAT, CSALROT (need to pass
    !     *                            NBS as well).
    !     *                          - Adds EMISGAT/EMISROT.
    !     *                          - Adds GTGAT/GTROT.
    !     *                          - Adds NT (NTLD in sfcproc2
    !     *                            call) to dimension land-only
    !     *                            ROT fields, consistent with
    !     *                            new comrow12.
    !     *                          - Unused IWMOS, JWMOS removed.
    !     * Jun 12, 2013 - M.Lazare. Previous version for gcm17,
    !     *                          called by "sfcproc".
    !     *                          CLASS scatter routine called by
    !     *                          "sfcproc" in new version gcm17.
    !     * NOTE: This contains the following changes compared to the
    !     *       working temporary version used in conjunction with
    !     *       updates to gcm16 (ie not official):
    !     *         1) {REF, BCSN} added for Maryam's new code.
    !     *         2) GFLX removed.
    !
    !     * OCT 25/11 - M.LAZARE.   REMOVE OPERATIONS ON INTERNAL
    !     *                         ROT ARRAYS (NOW DONE DIRECTLY
    !     *                         GAT->ROW IN SFCPROC).
    !     * OCT 07/11 - M.LAZARE.   REMOVE TSF.
    !     * OCT 05/11 - M.LAZARE.   ADD SFCH.
    !     * OCT 04/11 - M.LAZARE.   REMOVE ITCT.
    !     * MAR 23/06 - D.VERSEGHY. ADD WSNO, FSNO.
    !     * MAR 18/05 - D.VERSEGHY. ADDITIONAL VARIABLES.
    !     * FEB 18/05 - D.VERSEGHY. ADD "TSFS" VARIABLES.
    !     * AUG 05/04 - D.VERSEGHY. ADD NEW DIAGNOSTIC VARIABLES
    !     *                         ILMO, UE AND HBL.
    !     * AUG 15/02 - D.VERSEGHY. SCATTER OPERATION ON CLASS
    !     *                         VARIABLES.
    !
    implicit none
    !
    !     * INTEGER CONSTANTS.
    !
    integer, intent(in) :: NML, NL, NT, NM, ILG, IG, IC, ICP1, NBS
    integer             :: K, L
    !
    !     * LAND SURFACE PROGNOSTIC VARIABLES.
    !

    !! Suffix ROT refers to variables on original two-dimensional arrays.
    real, intent(out)    :: SALBROT(NL,NM,NBS) !< All-sky albedo  [  ]
    real, intent(out)    :: CSALROT(NL,NM,NBS) !< Clear-sky albedo  [  ]
    real, intent(out)    :: TBARROT(NL,NT,IG) !< Temperature of soil layers [K]
    real, intent(out)    :: THLQROT(NL,NT,IG) !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, intent(out)    :: THICROT(NL,NT,IG) !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, intent(out)    :: TSFSROT(NL,NT,4) !< Ground surface temperature over subarea [K]
    real, intent(out)    :: TPNDROT(NL,NT) !< Temperature of ponded water [K]
    real, intent(out)    :: ZPNDROT(NL,NT) !< Depth of ponded water on surface [m]
    real, intent(out)    :: TBASROT(NL,NT) !< Temperature of bedrock in third soil layer [K]
    real, intent(out)    :: ALBSROT(NL,NM) !< Snow albedo [  ]
    real, intent(out)    :: TSNOROT(NL,NM) !< Snowpack temperature [K]
    real, intent(out)    :: RHOSROT(NL,NM) !< Density of snow \f$[kg m^{-3} ]\f$
    real, intent(out)    :: SNOROT (NL,NM) !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, intent(out)    :: GTROT  (NL,NM) !< Effective surface black-body temperature  [K]
    real, intent(out)    :: TCANROT(NL,NT) !< Vegetation canopy temperature [K]
    real, intent(out)    :: RCANROT(NL,NT) !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(out)    :: SCANROT(NL,NT) !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(out)    :: GROROT (NL,NT) !< Vegetation growth index [  ]
    real, intent(out)    :: TACROT (NL,NT) !< Temperature of air within vegetation canopy [K]
    real, intent(out)    :: QACROT (NL,NT) !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
    real, intent(out)    :: WSNOROT(NL,NT) !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, intent(out)    :: CMAIROT(NL,NT) !< Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
    real, intent(out)    :: REFROT (NL,NM) !< Snow grain size  [m]
    real, intent(out)    :: BCSNROT(NL,NM) !< Black carbon mixing ratio \f$[kg m^{-3} ]\f$
    real, intent(out)    :: EMISROT(NL,NM) !< Surface emissivity  [  ]
    real, intent(out)    :: groundHeatFluxROT(NL, NM) !< Heat flux at soil surface \f$[W m^{-2} ]\f$
    !
    real, intent(in)     :: SALBGAT(ILG,NBS) !< All-sky albedo  [  ]
    real, intent(in)     :: CSALGAT(ILG,NBS) !< Clear-sky albedo  [  ]
    real, intent(in)     :: TBARGAT(ILG,IG) !< Temperature of soil layers [K]
    real, intent(in)     :: THLQGAT(ILG,IG) !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, intent(in)     :: THICGAT(ILG,IG) !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
    real, intent(in)     :: TSFSGAT(ILG,4) !< Ground surface temperature over subarea [K]
    real, intent(in)     :: TPNDGAT(ILG) !< Temperature of ponded water [K]
    real, intent(in)     :: ZPNDGAT(ILG) !< Depth of ponded water on surface [m]
    real, intent(in)     :: TBASGAT(ILG) !< Temperature of bedrock in third soil layer [K]
    real, intent(in)     :: ALBSGAT(ILG) !< Snow albedo [  ]
    real, intent(in)     :: TSNOGAT(ILG) !< Snowpack temperature [K]
    real, intent(in)     :: RHOSGAT(ILG) !< Density of snow \f$[kg m^{-3} ]\f$
    real, intent(in)     :: SNOGAT (ILG) !< Mass of snow pack \f$[kg m^{-2} ]\f$
    real, intent(in)     :: GTGAT  (ILG) !< Effective surface black-body temperature  [K]
    real, intent(in)     :: TCANGAT(ILG) !< Vegetation canopy temperature [K]
    real, intent(in)     :: RCANGAT(ILG) !< Intercepted liquid water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in)     :: SCANGAT(ILG) !< Intercepted frozen water stored on canopy \f$[kg m^{-2} ]\f$
    real, intent(in)     :: GROGAT (ILG) !< Vegetation growth index [  ]
    real, intent(in)     :: TACGAT (ILG) !< Temperature of air within vegetation canopy [K]
    real, intent(in)     :: QACGAT (ILG) !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
    real, intent(in)     :: WSNOGAT(ILG) !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
    real, intent(in)     :: CMAIGAT(ILG) !< Aggregated mass of vegetation canopy \f$[kg m^{-2} ]\f$
    real, intent(in)     :: REFGAT (ILG) !< Snow grain size  [m]
    real, intent(in)     :: BCSNGAT(ILG) !< Black carbon mixing ratio \f$[kg m^{-3} ]\f$
    real, intent(in)     :: EMISGAT(ILG) !< Surface emissivity  [  ]
    real, intent(in)     :: groundHeatFlux(ILG) !< Heat flux at soil surface \f$[W m^{-2} ]\f$

    !
    !     * GATHER-SCATTER INDEX ARRAYS.
    !
    integer, intent(in) :: ILMOS(ILG), JLMOS(ILG)
    !----------------------------------------------------------------------
    do K = 1,NML ! loop 100
      TPNDROT(ILMOS(K),JLMOS(K)) = TPNDGAT(K)
      ZPNDROT(ILMOS(K),JLMOS(K)) = ZPNDGAT(K)
      TBASROT(ILMOS(K),JLMOS(K)) = TBASGAT(K)
      ALBSROT(ILMOS(K),JLMOS(K)) = ALBSGAT(K)
      TSNOROT(ILMOS(K),JLMOS(K)) = TSNOGAT(K)
      RHOSROT(ILMOS(K),JLMOS(K)) = RHOSGAT(K)
      SNOROT (ILMOS(K),JLMOS(K)) = SNOGAT (K)
      GTROT  (ILMOS(K),JLMOS(K)) = GTGAT  (K)
      WSNOROT(ILMOS(K),JLMOS(K)) = WSNOGAT(K)
      TCANROT(ILMOS(K),JLMOS(K)) = TCANGAT(K)
      RCANROT(ILMOS(K),JLMOS(K)) = RCANGAT(K)
      SCANROT(ILMOS(K),JLMOS(K)) = SCANGAT(K)
      GROROT (ILMOS(K),JLMOS(K)) = GROGAT (K)
      TACROT (ILMOS(K),JLMOS(K)) = TACGAT (K)
      QACROT (ILMOS(K),JLMOS(K)) = QACGAT (K)
      CMAIROT(ILMOS(K),JLMOS(K)) = CMAIGAT(K)
      REFROT (ILMOS(K),JLMOS(K)) = REFGAT (K)
      BCSNROT(ILMOS(K),JLMOS(K)) = BCSNGAT(K)
      EMISROT(ILMOS(K),JLMOS(K)) = EMISGAT(K)
      groundHeatFluxROT(ILMOS(K),JLMOS(K)) = groundHeatFlux(K)

      !>
      !! The prognostic variables are scattered from the long, gathered arrays (collapsing the latitude and mosaic
      !! dimensions into one) back onto the original arrays using the pointer vectors generated in classGatherPrep.
      !!

    end do ! loop 100

    do L = 1,NBS ! loop 200
      do K = 1,NML
        SALBROT(ILMOS(K),JLMOS(K),L) = SALBGAT(K,L)
        CSALROT(ILMOS(K),JLMOS(K),L) = CSALGAT(K,L)
      end do
    end do ! loop 200

    do L = 1,IG ! loop 300
      do K = 1,NML
        TBARROT(ILMOS(K),JLMOS(K),L) = TBARGAT(K,L)
        THLQROT(ILMOS(K),JLMOS(K),L) = THLQGAT(K,L)
        THICROT(ILMOS(K),JLMOS(K),L) = THICGAT(K,L)
      end do
    end do ! loop 300

    do L = 1,4 ! loop 400
      do K = 1,NML
        TSFSROT(ILMOS(K),JLMOS(K),L) = TSFSGAT(K,L)
      end do
    end do ! loop 400

    return
  end subroutine classScatter
  !! @}

  !> \namespace classgatherscatter
  !! Transfers information between the 'gathered' and 'scattered' form of the CLASS data arrays.


end module classGatherScatter
