!> \file
!> Organizes calculation of radiation-related and other
!! surface parameters.
!! @author D. Verseghy, M. Lazare, J. Cole, Y. Wu, J. Melton
!!
subroutine radiationDriver (FC, FG, FCS, FGS, ALVSCN, ALIRCN, & ! Formerly CLASSA
                            ALVSG, ALIRG, ALVSCS, ALIRCS, ALVSSN, ALIRSN, &
                            ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN, &
                            TRVSCS, TRIRCS, FSVF, FSVFS, &
                            RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC, &
                            FRAICS, FSNOCS, DISP, DISPS, ZOMLNC, ZOMLCS, &
                            ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, &
                            CHCAP, CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP, &
                            CWLCPS, CWFCPS, RC, RCS, RBCOEF, FROOT, &
                            FROOTS, ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZSNOW, &
                            WSNOW, ALVS, ALIR, HTCC, HTCS, HTC, &
                            ALTG, ALSNO, TRSNOWC, TRSNOWG, &
                            WTRC, WTRS, WTRG, CMAI, FSNOW, &
                            FCANMX, ZOLN, ALVSC, ALIRC, PAIMAX, PAIMIN, &
                            CWGTMX, ZRTMAX, RSMIN, QA50, VPDA, VPDB, &
                            PSIGA, PSIGB, PAIDAT, HGTDAT, ACVDAT, ACIDAT, &
                            ASVDAT, ASIDAT, AGVDAT, AGIDAT, &
                            ALGWV, ALGWN, ALGDV, ALGDN, &
                            THLIQ, THICE, TBAR, RCAN, SNCAN, TCAN, &
                            GROWTH, SNO, TSNOW, RHOSNO, ALBSNO, ZBLEND, &
                            Z0ORO, SNOLIM, ZPLMG0, ZPLMS0, &
                            FCLOUD, TA, VPD, RHOAIR, COSZS, &
                            FSDB, FSFB, REFSNO, BCSNO, &
                            QSWINV, RADJ, DLON, RHOSNI, DELZ, DELZW, &
                            ZBOTW, THPOR, THLMIN, PSISAT, BI, PSIWLT, &
                            HCPS, ISAND, &
                            FCANCMX, ICTEM, ctem_on, RMATC, ZOLNC, CMASVEGC, &
                            AILC, PAIC, NOL2PFTS, SLAIC, &
                            AILCG, AILCGS, FCANC, FCANCS, &
                            IDAY, ILG, IL1, IL2, NBS, &
                            JL, N, IC, ICP1, IG, IDISP, IZREF, &
                            IWF, IPAI, IHGT, IALC, IALS, IALG, &
                            ISNOALB, ALVSCTM, ALIRCTM, ipeatland)

  !     * OCT  3/16 - J.Melton    Implementing Yuanqiao Wu's peatland code, added
  !                               ipeatland
  !     * AUG 30/16 - J.Melton    Replace ICTEMMOD with ctem_on (logical switch).
  !     * AUG 04/15 - M.LAZARE.   SPLIT FROOT INTO TWO ARRAYS, FOR CANOPY
  !     *                         AREAS WITH AND WITHOUT SNOW.
  !     * FEB 09/15 - D.VERSEGHY. New version for gcm18 and class 3.6:
  !     *                         - {ALGWV, ALGWN, ALGDV, ALGDN} are passed
  !     *                           in (originating in soilProperties) and then
  !     *                           passed into groundAlbedo, instead of
  !     *                           {ALGWET, ALGDRY}.
  !     * NOV 16/13 - J.COLE.     FINAL VERSION FOR GCM17:
  !     *                         - PASS "RHOSNO"IN TO snowAlbedoTransmiss TO
  !     *                           CALCULATE THE PROPER BC MIXING RATIO
  !     *                           IN SNOW.
  !     *                         - NEW "ISNOALBA" OPTION, BASED ON
  !     *                           4-BAND SOLAR.
  !     * NOV 14/11 - M.LAZARE.   IMPLEMENT CTEM SUPPORT, PRIMARILY
  !     *                         INVOLVING ADDITIONAL FIELDS TO PASS
  !     *                         IN/OUT OF NEW calcLandSurfParams ROUTINE. THIS
  !     *                         INCLUDES NEW INPUT ARRAY "PAIC".
  !     * NOV 30/06 - D.VERSEGHY. CONVERT RADJ TO REGULAR PRECISION.
  !     * APR 13/06 - D.VERSEGHY. SEPARATE GROUND AND SNOW ALBEDOS FOR
  !     *                         OPEN AND CANOPY-COVERED AREAS; KEEP
  !     *                         FSNOW AS OUTPUT ARRAY.
  !     * APR 06/06 - D.VERSEGHY. INTRODUCE MODELLING OF WSNOW.
  !     * MAR 14/05 - D.VERSEGHY. RENAME SCAN TO SNCAN (RESERVED NAME
  !     *                         IN F90); CHANGE SNOLIM FROM CONSTANT
  !     *                         TO VARIABLE.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * DEC 05/02 - D.VERSEGHY. NEW PARAMETERS FOR calcLandSurfParams.
  !     * JUL 31/02 - D.VERSEGHY. MODIFICATIONS ASSOCIATED WITH NEW
  !     *                         CALCULATION OF STOMATAL RESISTANCE.
  !     *                         SHORTENED CLASS3 COMMON BLOCK.
  !     * JUL 23/02 - D.VERSEGHY. MODIFICATIONS TO MOVE ADDITION OF AIR
  !     *                         TO CANOPY MASS INTO calcLandSurfParams; SHORTENED
  !     *                         CLASS4 COMMON BLOCK.
  !     * MAR 18/02 - D.VERSEGHY. NEW CALLS TO ALL SUBROUTINES TO ENABLE
  !     *                         ASSIGNMENT OF USER-SPECIFIED VALUES TO
  !     *                         ALBEDOS AND VEGETATION PROPERTIES; NEW
  !     *                         "CLASS8" COMMON BLOCK; MOVE CALCULATION
  !     *                         OF "FCLOUD" INTO CLASS DRIVER.
  !     * SEP 19/00 - D.VERSEGHY. PASS ADDITIONAL ARRAYS TO calcLandSurfParams IN COMMON
  !     *                         BLOCK CLASS7, FOR CALCULATION OF NEW
  !     *                         STOMATAL RESISTANCE COEFFICIENTS USED
  !     *                         IN energyBudgetPrep.
  !     * APR 12/00 - D.VERSEGHY. RCMIN NOW VARIES WITH VEGETATION TYPE:
  !     *                         PASS IN BACKGROUND ARRAY "RCMINX".
  !     * DEC 16/99 - D.VERSEGHY. ADD "XLEAF" ARRAY TO CLASS7 COMMON BLOCK
  !     *                         AND CALCULATION OF LEAF DIMENSION PARAMETER
  !     *                         "DLEAF" IN calcLandSurfParams.
  !     * NOV 16/98 - M.LAZARE.   "DLON" NOW PASSED IN AND USED DIRECTLY
  !     *                         (INSTEAD OF INFERRING FROM "LONSL" AND
  !     *                         "ILSL" WHICH USED TO BE PASSED) TO PASS
  !     *                         TO calcLandSurfParams TO CALCULATE GROWTH INDEX. THIS
  !     *                         IS DONE TO MAKE THE PHYSICS PLUG COMPATIBLE
  !     *                         FOR USE WITH THE RCM WHICH DOES NOT HAVE
  !     *                         EQUALLY-SPACED LONGITUDES.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         FIX BUG TO CALCULATE GROUND ALBEDO
  !     *                         UNDER CANOPIES AS WELL AS OVER BARE
  !     *                         SOIL.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     *                         ALSO, PASS IDISP TO SUBROUTINE calcLandSurfParams.
  !     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         VARIABLE SURFACE DETENTION CAPACITY
  !     *                         IMPLEMENTED.
  !     * AUG 16/95 - D.VERSEGHY. THREE NEW ARRAYS TO COMPLETE WATER
  !     *                         BALANCE DIAGNOSTICS.
  !     * OCT 14/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         REVISE CALCULATION OF FCLOUD TO
  !     *                         HANDLE CASES WHERE INCOMING SOLAR
  !     *                         RADIATION IS ZERO AT LOW SUN ANGLES.
  !     * NOV 24/92 - M.LAZARE.   CLASS - VERSION 2.1.
  !     *                         MODIFIED FOR MULTIPLE LATITUDES.
  !     * OCT 13/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. VISIBLE AND NEAR-IR ALBEDOS AND
  !     *                         TRANSMISSIVITIES FOR COMPONENTS OF
  !     *                         LAND SURFACE.
  !
  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IDAY, ILG, IL1, IL2, JL, IC, ICP1, IG, IDISP, IZREF, IWF, &
                         IPAI, IHGT, IALC, IALS, IALG, N, NBS
  integer :: I, J

  integer, intent(in) :: ISNOALB !< Switch to model snow albedo in two or more wavelength bands
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out) :: FC    (ILG)  !< Subarea fractional coverage of modelled area [ ]
  real, intent(out) :: FG    (ILG)  !< Subarea fractional coverage of modelled area [ ]
  real, intent(out) :: FCS   (ILG)  !< Subarea fractional coverage of modelled area [ ]
  real, intent(out) :: FGS   (ILG)  !< Subarea fractional coverage of modelled area [ ]
  real, intent(out) :: ALVSCN(ILG)  !< Visible albedo of vegetation over bare ground [ ]
  real, intent(out) :: ALIRCN(ILG)  !< Near-IR albedo of vegetation over bare ground [ ]
  real, intent(out) :: ALVSG (ILG)  !< Visible albedo of open bare ground [ ]
  real, intent(out) :: ALIRG (ILG)  !< Near-IR albedo of open bare ground [ ]
  real, intent(out) :: ALVSCS(ILG)  !< Visible albedo of vegetation over snow [ ]
  real, intent(out) :: ALIRCS(ILG)  !< Near-IR albedo of vegetation over snow [ ]
  real, intent(out) :: ALVSSN(ILG)  !< Visible albedo of open snow cover [ ]
  real, intent(out) :: ALIRSN(ILG)  !< Near-IR albedo of open snow cover [ ]
  real, intent(out) :: ALVSGC(ILG)  !< Visible albedo of bare ground under vegetation [ ]
  real, intent(out) :: ALIRGC(ILG)  !< Near-IR albedo of bare ground under vegetation [ ]
  real, intent(out) :: ALVSSC(ILG)  !< Visible albedo of snow under vegetation [ ]
  real, intent(out) :: ALIRSC(ILG)  !< Near-IR albedo of snow under vegetation [ ]
  real, intent(out) :: TRVSCN(ILG)  !< Visible transmissivity of vegetation over bare ground [ ]
  real, intent(out) :: TRIRCN(ILG)  !< Near-IR transmissivity of vegetation over bare ground [ ]
  real, intent(out) :: TRVSCS(ILG)  !< Visible transmissivity of vegetation over snow [ ]
  real, intent(out) :: TRIRCS(ILG)  !< Near-IR transmissivity of vegetation over snow [ ]
  real, intent(out) :: FSVF  (ILG)  !< Sky view factor for bare ground under canopy [ ]
  real, intent(out) :: FSVFS (ILG)  !< Sky view factor for snow under canopy [ ]
  real, intent(out) :: RAICAN(ILG)  !< Intercepted liquid water stored on canopy over
  !! bare ground \f$[kg m^{-2}]\f$
  real, intent(out) :: RAICNS(ILG)  !< Intercepted liquid water stored on canopy over
  !! snow \f$[kg m^{-2}]\f$
  real, intent(out) :: SNOCAN(ILG)  !< Intercepted frozen water stored on canopy over
  !! bare soil \f$[kg m^{-2}]\f$
  real, intent(out) :: SNOCNS(ILG)  !< Intercepted frozen water stored on canopy over
  !! snow \f$[kg m^{-2}]\f$
  real, intent(out) :: FRAINC(ILG)  !< Fractional coverage of canopy by liquid water
  !! over snow-free subarea [ ]
  real, intent(out) :: FSNOWC(ILG)  !< Fractional coverage of canopy by frozen water
  !! over snow-free subarea [ ]
  real, intent(out) :: FRAICS(ILG)  !< Fractional coverage of canopy by liquid water
  !! over snow-covered subarea [ ]
  real, intent(out) :: FSNOCS(ILG)  !< Fractional coverage of canopy by frozen water
  !! over snow-covered subarea [ ]
  real, intent(out) :: DISP  (ILG)  !< Displacement height of vegetation over bare ground [m]
  real, intent(out) :: DISPS (ILG)  !< Displacement height of vegetation over snow [m]
  real, intent(out) :: ZOMLNC(ILG)  !< Logarithm of roughness length for momentum of
  !! vegetation over bare ground [ ]
  real, intent(out) :: ZOMLCS(ILG)  !< Logarithm of roughness length for momentum of
  !! vegetation over snow [ ]
  real, intent(out) :: ZOELNC(ILG)  !< Logarithm of roughness length for heat of
  !! vegetation over bare ground [ ]
  real, intent(out) :: ZOELCS(ILG)  !< Logarithm of roughness length for heat of
  !! vegetation over snow [ ]
  real, intent(out) :: ZOMLNG(ILG)  !< Logarithm of roughness length for momentum of bare ground [ ]
  real, intent(out) :: ZOMLNS(ILG)  !< Logarithm of roughness length for momentum of snow [ ]
  real, intent(out) :: ZOELNG(ILG)  !< Logarithm of roughness length for heat of bare ground [ ]
  real, intent(out) :: ZOELNS(ILG)  !< Logarithm of roughness length for heat of snow [ ]
  real, intent(out) :: CHCAP (ILG)  !< Heat capacity of canopy over bare ground \f$ [J m^{-2} K^{-1}]\f$
  real, intent(out) :: CHCAPS(ILG)  !< Heat capacity of canopy over snow \f$[J m^{-2} K^{-1}]\f$
  real, intent(out) :: CMASSC(ILG)  !< Mass of canopy over bare ground \f$[kg m^{-2}]\f$
  real, intent(out) :: CMASCS(ILG)  !< Mass of canopy over snow \f$[kg m^{-2}]\f$
  real, intent(out) :: CWLCAP(ILG)  !< Storage capacity of canopy over bare ground for
  !! liquid water \f$[kg m^{-2}]\f$
  real, intent(out) :: CWFCAP(ILG)  !< Storage capacity of canopy over bare ground for
  !! frozen water \f$[kg m^{-2}]\f$
  real, intent(out) :: CWLCPS(ILG)  !< Storage capacity of canopy over snow for liquid
  !! water \f$[kg m^{-2}]\f$
  real, intent(out) :: CWFCPS(ILG)  !< Storage capacity of canopy over snow for frozen
  !! water \f$[kg m^{-2}]\f$
  real, intent(out) :: RC    (ILG)  !< Stomatal resistance of vegetation over bare
  !! ground \f$[s m^{-1}]\f$
  real, intent(out) :: RCS   (ILG)  !< Stomatal resistance of vegetation over snow \f$[s m^{-1}]\f$
  real, intent(out) :: ZPLIMC(ILG)  !< Maximum water ponding depth for ground under canopy [m]
  real, intent(out) :: ZPLIMG(ILG)  !< Maximum water ponding depth for bare ground [m]
  real, intent(out) :: ZPLMCS(ILG)  !< Maximum water ponding depth for ground under
  !! snow under canopy [m]
  real, intent(out) :: ZPLMGS(ILG)  !< Maximum water ponding depth for ground under snow [m]
  real, intent(out) :: RBCOEF(ILG)  !< Parameter for calculation of leaf boundary resistance
  real, intent(out) :: TRSNOWC(ILG) !< Short-wave transmissivity of snow pack under vegetation  [  ]
  real, intent(out) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_s)\f$
  real, intent(out) :: WSNOW (ILG)  !< Liquid water content of snow pack \f$[kg m^{-2}]\f$
  real, intent(out) :: ALVS  (ILG)  !< Diagnosed total visible albedo of land surface [ ]
  real, intent(out) :: ALIR  (ILG)  !< Diagnosed total near-infrared albedo of land surface [ ]
  real, intent(out) :: HTCC  (ILG)  !< Diagnosed internal energy change of vegetation
  !! canopy due to conduction and/or change in mass \f$[W m^{-2}]\f$
  real, intent(out) :: HTCS  (ILG)  !< Diagnosed internal energy change of snow pack
  !! due to conduction and/or change in mass \f$[W m^{-2}]\f$
  real, intent(out) :: WTRC  (ILG)  !< Diagnosed residual water transferred off the
  !! vegetation canopy \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: WTRS  (ILG)  !< Diagnosed residual water transferred into or
  !! out of the snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: WTRG  (ILG)  !< Diagnosed residual water transferred into or
  !! out of the soil \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: CMAI  (ILG)  !< Aggregated mass of vegetation canopy \f$[kg m^{-2}]\f$
  real, intent(out) :: FSNOW (ILG)  !< Diagnosed fractional snow coverage [ ]
  !
  real, intent(out) :: FROOT (ILG,IG)   !< Fraction of total transpiration contributed by soil layer over snow-free subarea  [  ]
  real, intent(out) :: FROOTS(ILG,IG)   !< Fraction of total transpiration contributed
  !! by snow-covered subarea [ ]
  real, intent(out) :: HTC   (ILG,IG)   !< Diagnosed internal energy change of soil
  !! layer due to conduction and/or change in mass \f$[W m^{-2}]\f$

  real, intent(out) :: TRSNOWG(ILG,NBS) !< Short-wave transmissivity of snow pack in bare areas  [  ]
  real, intent(out) :: ALTG(ILG,NBS)    !< Total albedo in each modelled wavelength band  [  ]
  real, intent(out) :: ALSNO(ILG,NBS)   !< Albedo of snow in each modelled wavelength band  [  ]
  !
  !     * INPUT ARRAYS DEPENDENT ON LONGITUDE.
  !
  real, intent(in) :: FCANMX(ILG,ICP1) !< Maximum fractional coverage of modelled
  !! area by vegetation category [ ]
  real, intent(in) :: ZOLN  (ILG,ICP1) !< Natural logarithm of maximum roughness
  !! length of vegetation category [ ]
  real, intent(inout) :: ALVSC (ILG,ICP1) !< Background average visible albedo of
  !! vegetation category [ ]
  real, intent(inout) :: ALIRC (ILG,ICP1) !< Background average near-infrared albedo of
  !! vegetation category [ ]
  real, intent(in) :: PAIMAX(ILG,IC)   !< Maximum plant area index of vegetation
  !! category [ ]
  real, intent(in) :: PAIMIN(ILG,IC)   !< Minimum plant area index of vegetation
  !! category [ ]
  real, intent(in) :: CWGTMX(ILG,IC)   !< Maximum canopy mass for vegetation category \f$[kg m^{-2}]\f$
  real, intent(in) :: ZRTMAX(ILG,IC)   !< Maximum rooting depth of vegetation
  !! category [m]
  real, intent(in) :: RSMIN (ILG,IC)   !< Minimum stomatal resistance of vegetation
  !! category \f$[s m^{-1}]\f$
  real, intent(in) :: QA50  (ILG,IC)   !< Reference value of incoming shortwave
  !! radiation for vegetation category (used in
  !! stomatal resistance calculation) \f$[W m^{-2}]\f$
  real, intent(in) :: VPDA  (ILG,IC)   !< Vapour pressure deficit coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) [ ]
  real, intent(in) :: VPDB  (ILG,IC)   !< Vapour pressure deficit coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) [ ]
  real, intent(in) :: PSIGA (ILG,IC)   !< Soil moisture suction coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) [ ]
  real, intent(in) :: PSIGB (ILG,IC)   !< Soil moisture suction coefficient for
  !! vegetation category (used in stomatal
  !! resistance calculation) [ ]
  real, intent(in) :: PAIDAT(ILG,IC)   !< Optional user-specified value of plant area
  !! indices of vegetation categories to
  !! override CLASS-calculated values [ ]
  real, intent(in) :: HGTDAT(ILG,IC)   !< Optional user-specified values of height of
  !! vegetation categories to override CLASS-
  !! calculated values [m]
  real, intent(in) :: ACVDAT(ILG,IC)   !< Optional user-specified value of canopy
  !! visible albedo to override CLASS-calculated
  !! value [ ]
  real, intent(in) :: ACIDAT(ILG,IC)   !< Optional user-specified value of canopy
  !! near-infrared albedo to override CLASS-
  !! calculated value [ ]
  real, intent(in) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
  !! layers \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THICE (ILG,IG)   !< Volumetric frozen water content of soil
  !! layers \f$[m^3 m^{-3}]\f$
  real, intent(in) :: TBAR  (ILG,IG)   !< Temperature of soil layers [K]
  !
  real, intent(in) :: ASVDAT(ILG)  !< Optional user-specified value of snow visible
  !! albedo to override CLASS-calculated value [ ]
  real, intent(in) :: ASIDAT(ILG)  !< Optional user-specified value of snow near-
  !! infrared albedo to override CLASS-calculated value [ ]
  real, intent(in) :: AGVDAT(ILG)  !< Optional user-specified value of ground visible
  !! albedo to override CLASS-calculated value [ ]
  real, intent(in) :: AGIDAT(ILG)  !< Optional user-specified value of ground near-
  !! infrared albedo to override CLASS-calculated value [ ]
  real, intent(in) :: ALGWV(ILG)   !< Reference albedo for saturated soil (visible) [ ]
  real, intent(in) :: ALGWN(ILG)   !< Reference albedo for saturated soil (NIR) [ ]
  real, intent(in) :: ALGDV(ILG)   !< Reference albedo for dry soil (visible) [ ]
  real, intent(in) :: ALGDN(ILG)   !< Reference albedo for dry soil (NIR) [ ]
  real, intent(in) :: RHOSNI(ILG)  !< Density of fresh snow \f$[kg m^{-3}]\f$
  real, intent(in) :: Z0ORO (ILG)  !< Orographic roughness length [m]
  real, intent(in) :: RCAN  (ILG)  !< Intercepted liquid water stored on canopy \f$[kg m^{-2}]\f$
  real, intent(in) :: SNCAN (ILG)  !< Intercepted frozen water stored on canopy \f$[kg m^{-2}]\f$
  real, intent(in) :: TCAN  (ILG)  !< Vegetation canopy temperature [K]
  real, intent(in) :: GROWTH(ILG)  !< Vegetation growth index [ ]
  real, intent(in) :: SNO   (ILG)  !< Mass of snow pack \f$[kg m^{-2}] (W_s)\f$
  real, intent(in) :: TSNOW (ILG)  !< Snowpack temperature [K]
  real, intent(in) :: RHOSNO(ILG)  !< Density of snow \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(in) :: ALBSNO(ILG)  !< Snow albedo [ ]
  real, intent(in) :: FCLOUD(ILG)  !< Fractional cloud cover [ ]
  real, intent(in) :: TA    (ILG)  !< Air temperature at reference height [K]
  real, intent(in) :: VPD   (ILG)  !< Vapour pressure deficit of air [mb]
  real, intent(in) :: RHOAIR(ILG)  !< Density of air \f$[kg m^{-3}]\f$
  real, intent(in) :: COSZS (ILG)  !< Cosine of solar zenith angle [ ]
  real, intent(in) :: QSWINV(ILG)  !< Visible radiation incident on horizontal surface \f$[W m^{-2}]\f$
  real, intent(in) :: DLON  (ILG)  !< Longitude of grid cell (east of Greenwich) [degrees]
  real, intent(in) :: ZBLEND(ILG)  !< Atmospheric blending height for surface
  !! roughness length averaging [m]
  real, intent(in) :: SNOLIM(ILG)  !< Limiting snow depth below which coverage is
  !! < 100% \f$[m] (z_{s,lim})\f$
  real, intent(in) :: ZPLMG0(ILG)  !< Maximum water ponding depth for snow-free
  !! subareas (user-specified when running MESH code) [m]
  real, intent(in) :: ZPLMS0(ILG)  !< Maximum water ponding depth for snow-covered
  !! subareas (user-specified when running MESH code) [m]
  real, intent(in) :: RADJ  (ILG)  !< Latitude of grid cell (positive north of equator) [rad]
  real, intent(in) :: REFSNO(ILG)  !< Snow grain size (for ISNOALB=1 option)  [m]
  real, intent(in) :: BCSNO(ILG)   !< Black carbon mixing ratio (for ISNOALB=1 option)  \f$[kg m^{-3}]\f$
  real, intent(in) :: FSDB(ILG,NBS) !< Direct solar radiation in each modelled wavelength band  \f$[W m^{-2}]\f$
  real, intent(in) :: FSFB(ILG,NBS) !< Diffuse solar radiation in each modelled wavelength band \f$[W m^{-2}]\f$
  !
  !    * SOIL PROPERTY ARRAYS.
  !
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable thickness of soil layer [m]
  real, intent(in) :: DELZ  (IG)       !< Soil layer thickness [m]
  real, intent(in) :: ZBOTW (ILG,IG)   !< Depth to permeable bottom of soil layer [m]
  real, intent(in) :: THPOR (ILG,IG)   !< Pore volume in soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !! remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: PSISAT(ILG,IG)   !< Soil moisture suction at saturation [m]
  real, intent(in) :: BI    (ILG,IG)   !< Clapp and Hornberger empirical "b" parameter [ ]
  real, intent(in) :: PSIWLT(ILG,IG)   !< Soil moisture suction at wilting point [m]
  real, intent(in) :: HCPS  (ILG,IG)   !< Volumetric heat capacity of soil particles \f$[J m^{-3}]\f$
  !
  integer, intent(in)   :: ISAND (ILG,IG)  !< Sand content flag
  !
  !     * CTEM-RELATED FIELDS.

  !     * AILCG  - GREEN LAI FOR USE IN PHOTOSYNTHESIS
  !     * AILCGS - GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
  !     * AILCMIN- MIN. LAI FOR CTEM PFTs
  !     * AILCMAX- MAX. LAI FOR CTEM PFTs
  !     * L2MAX  - MAXIMUM NUMBER OF LEVEL 2 CTEM PFTs
  !     * NOL2PFTS - NUMBER OF LEVEL 2 CTEM PFTs
  !     * FCANC  - FRACTION OF CANOPY OVER GROUND FOR CTEM's 9 PFTs
  !     * FCANCS - FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
  !     * SEE BIO2STR SUBROUTINE FOR DEFINITION OF OTHER CTEM VARIABLES
  !
  real, intent(in) :: FCANCMX(ILG,ICTEM), RMATC(ILG,IC,IG), &
                      AILC(ILG,IC), PAIC(ILG,IC), &
                      AILCG(ILG,ICTEM), AILCGS(ILG,ICTEM), &
                      FCANC(ILG,ICTEM), FCANCS(ILG,ICTEM), &
                      ZOLNC(ILG,IC), CMASVEGC(ILG,IC), &
                      SLAIC(ILG,IC), ALVSCTM(ILG,IC), &
                      ALIRCTM(ILG,IC)

  integer, intent(in) :: ICTEM, NOL2PFTS(IC)

  logical, intent(in) :: ctem_on
  integer, intent(in)  :: ipeatland(ilg) !< Peatland flag: 0 = not a peatland, 1 = bog, 2 = fen
  !
  !     * INTERNAL WORK ARRAYS FOR THIS AND ASSOCIATED SUBROUTINES.
  !
  real :: RMAT(ILG,IC,IG), H(ILG,IC),      HS(ILG,IC), &
          PAI(ILG,IC),     PAIS(ILG,IC),   FCAN(ILG,IC), &
          FCANS(ILG,IC),   CXTEFF(ILG,IC), AIL(ILG,IC), &
          RCACC(ILG,IC),   RCG(ILG,IC),    RCV(ILG,IC)
  !
  real :: PSIGND(ILG), CWCPAV(ILG), FRTOTS(ILG), &
          GROWA (ILG), GROWN (ILG), GROWB (ILG), &
          RRESID(ILG), SRESID(ILG), FRTOT (ILG), &
          TRVS  (ILG), TRIR  (ILG), RCT   (ILG), &
          GC    (ILG), PAICAN(ILG), PAICNS(ILG)
  !
  !------------------------------------------------------------------
  !>
  !! In the first loop, the depth of snow \f$z_s\f$ is calculated from the
  !! snow mass \f$W_s\f$ and density \f$\varrho_s\f$ as:
  !! \f$ z_s = W_s / \varrho_s. \f$
  !!
  !! If the calculated value of \f$z_s\f$ is less than the limiting snow
  !! depth \f$z_{s, lim}\f$, the snow cover is deemed to be discontinuous. The
  !! fractional snow coverage \f$X_s\f$ of the modelled area is evaluated
  !! as
  !! \f$ X_s = z_s / z_{s, lim} \f$
  !! and the snow depth is reset to \f$z_{s, lim}\f$. The water content of the
  !! snow pack is corrected according to the new snow fractional area.
  !!
  !! The subarea albedo and transmissivity arrays (for canopy, bare
  !! ground, canopy over snow and snow over bare ground) are next
  !! initialized to zero, and the four radiationDriver subsidiary subroutines
  !! are called in turn: calcLandSurfParams to evaluate various model parameters
  !! for the four subareas, groundAlbedo to calculate the ground surface albedo,
  !! snowAlbedoTransmiss to calculate the snow albedo and transmissivity, and canopyAlbedoTransmiss
  !! to calculate the canopy albedo, transmissivity and stomatal resistance.
  !! Finally, the overall visible, near-infrared and total albedos for the
  !! modelled area are determined as weighted averages over the four subareas.
  !!
  !     * CALCULATION OF SNOW DEPTH ZSNOW AND FRACTIONAL SNOW COVER
  !     * FSNOW; INITIALIZATION OF COMPUTATIONAL ARRAYS.
  !
  do I = IL1,IL2
    if (SNO(I) > 0.0) then
      ZSNOW(I) = SNO(I) / RHOSNO(I)
      if (ZSNOW(I) >= (SNOLIM(I) - 0.00001)) then
        FSNOW(I) = 1.0
      else
        FSNOW(I) = ZSNOW(I) / SNOLIM(I)
        ZSNOW(I) = SNOLIM(I)
        WSNOW(I) = WSNOW(I) / FSNOW(I)
      end if
    else
      ZSNOW(I) = 0.0
      FSNOW(I) = 0.0
    end if
    !
    ALVSCN(I) = 0.0
    ALIRCN(I) = 0.0
    ALVSCS(I) = 0.0
    ALIRCS(I) = 0.0
    TRVSCN(I) = 0.0
    TRIRCN(I) = 0.0
    TRVSCS(I) = 0.0
    TRIRCS(I) = 0.0
    ALVSSN(I) = 0.0
    ALIRSN(I) = 0.0
    ALVSG (I) = 0.0
    ALIRG (I) = 0.0
    ALVSGC(I) = 0.0
    ALIRGC(I) = 0.0
    ALVSSC(I) = 0.0
    ALIRSC(I) = 0.0
    TRSNOWC(I) = 0.0

    TRSNOWG(I,1:NBS) = 0.0
    ALTG(I,1:NBS)    = 0.0
    ALSNO(I,1:NBS)   = 0.0

  end do ! loop 100
  !
  ! ===================== CTEM =====================================\
  !     IF USING DYNAMIC VEGETATION COMPONENT OF CTEM, REPLACE ALBEDOS
  !     THAT ARE BASED ON CTEM.

  if (ctem_on) then
    do J = 1,IC
      do I = IL1,IL2
        ALVSC(I,J) = ALVSCTM(I,J)
        ALIRC(I,J) = ALIRCTM(I,J)
      end do
    end do
  end if
  !===================== CTEM =====================================/
  !
  !     * PREPARATION.
  !
  call calcLandSurfParams(FC, FG, FCS, FGS, PAICAN, PAICNS, FSVF, FSVFS, & ! Formerly called APREP
                          FRAINC, FSNOWC, FRAICS, FSNOCS, RAICAN, RAICNS, SNOCAN, &
                          SNOCNS, DISP, DISPS, ZOMLNC, ZOMLCS, ZOELNC, ZOELCS, &
                          ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, CHCAP, CHCAPS, CMASSC, &
                          CMASCS, CWLCAP, CWFCAP, CWLCPS, CWFCPS, RBCOEF, &
                          ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, HTCC, HTCS, HTC, &
                          FROOT, FROOTS, &
                          WTRC, WTRS, WTRG, CMAI, PAI, PAIS, AIL, FCAN, FCANS, PSIGND, &
                          FCANMX, ZOLN, PAIMAX, PAIMIN, CWGTMX, ZRTMAX, &
                          PAIDAT, HGTDAT, THLIQ, THICE, TBAR, RCAN, SNCAN, &
                          TCAN, GROWTH, ZSNOW, TSNOW, FSNOW, RHOSNO, SNO, Z0ORO, &
                          ZBLEND, ZPLMG0, ZPLMS0, &
                          TA, RHOAIR, RADJ, DLON, RHOSNI, DELZ, DELZW, ZBOTW, &
                          THPOR, THLMIN, PSISAT, BI, PSIWLT, HCPS, ISAND, &
                          ILG, IL1, IL2, JL, IC, ICP1, IG, IDAY, IDISP, IZREF, IWF, &
                          IPAI, IHGT, RMAT, H, HS, CWCPAV, GROWA, GROWN, GROWB, &
                          RRESID, SRESID, FRTOT, FRTOTS, &
                          FCANCMX, ICTEM, ctem_on, RMATC, &
                          AILC, PAIC, AILCG, NOL2PFTS, &
                          AILCGS, FCANCS, FCANC, ZOLNC, CMASVEGC, SLAIC, &
                          ipeatland)
  !
  !     * BARE SOIL ALBEDOS.
  !
  call groundAlbedo(ALVSG, ALIRG, ALVSGC, ALIRGC, & ! Formerly GRALB
                    ALGWV, ALGWN, ALGDV, ALGDN, &
                    THLIQ, FSNOW, ALVSC(1,ICP1), ALIRC(1,ICP1), &
                    FCANMX(1, ICP1), AGVDAT, AGIDAT, FG, ISAND, &
                    ILG, IG, IL1, IL2, JL, IALG)


  !     * SNOW ALBEDOS AND TRANSMISSIVITY.
  !
  call snowAlbedoTransmiss(ALVSSN, ALIRSN, ALVSSC, ALIRSC, ALBSNO, & ! Formerly SNOALBA
                           TRSNOWC, ALSNO, TRSNOWG, FSDB, FSFB, RHOSNO, &
                           REFSNO, BCSNO, SNO, COSZS, ZSNOW, FSNOW, ASVDAT, ASIDAT, &
                           ALVSG, ALIRG, &
                           ILG, IG, IL1, IL2, JL, IALS, NBS, ISNOALB)

  !     * CANOPY ALBEDOS AND TRANSMISSIVITIES, AND VEGETATION
  !     * STOMATAL RESISTANCE.
  !
  call canopyAlbedoTransmiss(ALVSCN, ALIRCN, ALVSCS, ALIRCS, TRVSCN, TRIRCN, & ! Formerly CANALB
                             TRVSCS, TRIRCS, RC, RCS, &
                             ALVSC, ALIRC, RSMIN, QA50, VPDA, VPDB, PSIGA, PSIGB, &
                             FC, FCS, FSNOW, FSNOWC, FSNOCS, FCAN, FCANS, PAI, PAIS, &
                             AIL, PSIGND, FCLOUD, COSZS, QSWINV, VPD, TA, &
                             ACVDAT, ACIDAT, ALVSGC, ALIRGC, ALVSSC, ALIRSC, &
                             ILG, IL1, IL2, JL, IC, ICP1, IG, IALC, &
                             CXTEFF, TRVS, TRIR, RCACC, RCG, RCV, RCT, GC)
  !
  !     * EFFECTIVE WHOLE-SURFACE VISIBLE AND NEAR-IR ALBEDOS.
  !
  do I = IL1,IL2
    ALVS(I) = FC(I) * ALVSCN(I) + FG(I) * ALVSG(I) + FCS(I) * ALVSCS(I) + &
              FGS(I) * ALVSSN(I)
    ALIR(I) = FC(I) * ALIRCN(I) + FG(I) * ALIRG(I) + FCS(I) * ALIRCS(I) + &
              FGS(I) * ALIRSN(I)
  end do ! loop 500
  !
  if (ISNOALB == 0) then
    do I = IL1,IL2
      ALTG(I,1) = ALVS(I)
      ALTG(I,2) = ALIR(I)
      ALTG(I,3) = ALIR(I)
      ALTG(I,4) = ALIR(I)
    end do ! I
  else if (ISNOALB == 1) then
    do I = IL1,IL2
      ALTG(I,1) = FC(I) * ALVSCN(I) + FG(I) * ALVSG(I) + FCS(I) * ALVSCS(I) + &
                  FGS(I) * ALSNO(I,1)
      ALTG(I,2) = FC(I) * ALIRCN(I) + FG(I) * ALIRG(I) + FCS(I) * ALIRCS(I) + &
                  FGS(I) * ALSNO(I,2)
      ALTG(I,3) = FC(I) * ALIRCN(I) + FG(I) * ALIRG(I) + FCS(I) * ALIRCS(I) + &
                  FGS(I) * ALSNO(I,3)
      ALTG(I,4) = FC(I) * ALIRCN(I) + FG(I) * ALIRG(I) + FCS(I) * ALIRCS(I) + &
                  FGS(I) * ALSNO(I,4)
    end do ! I
  end if ! ISNOALB
  return
end subroutine radiationDriver
