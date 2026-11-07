!> \file
!! Calculates various land surface parameters.
!! @author D. Verseghy, M. Lazare, V. Fortin, V. Arora, E. Chan, P. Bartlett, Y. Wu, J. Melton, A. Wu, Y. Delage
!!
!! This subroutine is adaptable to any number of vegetation categories recognized by CLASS
!! (e.g. needleleaf trees, broadleaf trees, crops and grass), if an unknown PFT is present, a call to abort  number of
!! is performed.

subroutine calcLandSurfParams (FC, FG, FCS, FGS, PAICAN, PAICNS, FSVF, FSVFS, & ! Formerly APREP
                               FRAINC, FSNOWC, FRAICS, FSNOCS, RAICAN, RAICNS, SNOCAN, &
                               SNOCNS, DISP, DISPS, ZOMLNC, ZOMLCS, ZOELNC, ZOELCS, &
                               ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, CHCAP, CHCAPS, CMASSC, &
                               CMASCS, CWLCAP, CWFCAP, CWLCPS, CWFCPS, RBCOEF, ZPLIMC, &
                               ZPLIMG, ZPLMCS, ZPLMGS, HTCC, HTCS, HTC, FROOT, FROOTS, &
                               WTRC, WTRS, WTRG, CMAI, PAI, PAIS, AIL, FCAN, FCANS, &
                               PSIGND, FCANMX, ZOLN, PAIMAX, PAIMIN, CWGTMX, ZRTMAX, &
                               PAIDAT, HGTDAT, THLIQ, THICE, TBAR, RCAN, SNCAN, TCAN, &
                               GROWTH, ZSNOW, TSNOW, FSNOW, RHOSNO, SNO, Z0ORO, ZBLEND, &
                               ZPLMG0, ZPLMS0, TA, RHOAIR, RADJ, DLON, RHOSNI, DELZ, &
                               DELZW, ZBOTW, THPOR, THLMIN, PSISAT, BI, PSIWLT, HCPS, &
                               ISAND, ILG, IL1, IL2, JL, IC, ICP1, IG, IDAY, IDISP, &
                               IZREF, IWF, IPAI, IHGT, RMAT, H, HS, CWCPAV, GROWA, &
                               GROWN, GROWB, RRESID, SRESID, FRTOT, FRTOTS, FCANCMX, &
                               ICTEM, ctem_on, RMATC, AILC, PAIC, AILCG, &
                               NOL2PFTS, AILCGS, FCANCS, FCANC, ZOLNC, CMASVEGC, &
                               SLAIC, ipeatland)

  !     * JAN 2019 - J. Melton    Remove common block parameters, use classicParams instead.
  !     * Nov 2018 - J. Melton/S.Sun  Allow >4 original PFTs. Revert change of JAN 05/16 as the XLEAF
  !                               bug makes this unneccesary.
  !     * SEP  3/16 - J.Melton/Yuanqiao Wu - Bring in peatlands code
  !     * AUG 30/16 - J.Melton    Replace ICTEMMOD with ctem_on (logical switch).
  !
  !     * JAN 14/16 - J.MELTON    IN LOOP 450, MODIFIED SO IT COULD HANDLE >3 SOIL LAYERS
  !                               ALSO REMOVED SOME HARDCODED CTEM CODE FOR MORE FLEXIBLE FORMS
  !     * JAN 05/16 - J.MELTON.   TREE PFTS NOW HAVE A MINIMUM PAI OF 1 (LIKE
  !     *                         CROPS AND GRASSES) TO PREVENT WILD CANOPY TEMPERATURE
  !     *                         VALUES WHEN THE CANOPY IS SMALL.
  !     * AUG 04/15 - D.VERSEGHY. SPLIT FROOT INTO TWO ARRAYS, FOR CANOPY
  !     *                         AREAS WITH AND WITHOUT SNOW.
  !     * SEP 05/12 - J.MELTON.   CHANGED IDAY
  !                               CONVERSION FROM FLOAT TO real, REINTEGRATED
  !                               CTEM
  !     * NOV 15/11 - M.LAZARE.   CTEM ADDED. CALCULATIONS ARE DIFFERENT
  !     *                         IN SEVERAL AREAS, UNDER CONTROL OF
  !     *                         "ICTEMMOD" SWITCH (ICTEMMOD=0 REVERTS
  !     *                         BACK TO APREP4 FORMULATION). THIS
  !     *                         INCLUDES NEW INPUT "PAIC".
  !     * OCT 07/11 - V.FORTIN/D.VERSEGHY. MAKE THE LIMITING PONDING DEPTH
  !     *                         CALCULATION OVER ORGANIC SOILS THE SAME
  !     *                         AS OVER MINERAL SOILS (LOOP 175).
  !     * DEC 23/09 - D.VERSEGHY. IN LIMITING PONDING DEPTH CALCULATIONS,
  !     *                         IDENTIFY PEATLANDS WHERE ISAND(I, 2)=-2
  !     * JAN 06/09 - D.VERSEGHY. REINTRODUCE CHECKS ON FRACTIONAL AREAS.
  !     * MAR 25/08 - D.VERSEGHY. DISTINGUISH BETWEEN LEAF AREA INDEX
  !     *                         AND PLANT AREA INDEX.
  !     * JAN 17/08 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS; REMOVE
  !     *                         SUPERFLUOUS CHECKS ON FRACTIONAL AREAS.
  !     * NOV 30/06 - E.CHAN/M.LAZARE/D.VERSEGHY. CHANGE RADJ TO REAL
  !     *                         ENSURE CONSISTENCY IN CALCULATION
  !     *                         OF FRACTIONAL CANOPY AREAS.
  !     * SEP 13/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 100,
  !     *                         450 LOOPS.
  !     * MAR 14/05 - D.VERSEGHY. RENAME SCAN TO SNCAN (RESERVED NAME
  !     *                         IN F90); TREAT SOIL FROZEN WATER AS ICE
  !     *                         VOLUME RATHER THAN AS EQUIVALENT WATER.
  !     * MAR 03/05 - Y.DELAGE.   ADD CONTRIBUTION OF SUBGRID-SCALE
  !     *                         OROGRAPHY TO ROUGHNESS LENGTH.
  !     * JAN 12/05 - P.BARTLETT/D.VERSEGHY. DETERMINE SEPARATE CANOPY
  !     *                         WATER INTERCEPTION CAPACITIES FOR
  !     *                         RAIN AND SNOW, AND NEW FRACTIONAL
  !     *                         CANOPY COVERAGE OF INTERCEPTED RAIN
  !     *                         AND SNOW; DEFINE NEW PARAMETER RBCOEF
  !     *                         FOR RBINV CALCULATION IN energBalVegSolve.
  !     * NOV 03/04 - D.VERSEGHY. CHANGE RADJ AND DLON TO GATHERED
  !     *                         VARIABLES AND REMOVE ILAND ARRAY
  !     *                         ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 05/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS
  !     *                         AGAINST ROUNDOFF ERRORS.
  !     * JUL 02/03 - D.VERSEGHY. RATIONALIZE ASSIGNMENT OF RESIDUAL
  !     *                         CANOPY MOISTURE TO SOIL LAYERS.
  !     * DEC 05/02 - Y.DELAGE/D.VERSEGHY. ADD PARTS OF CANOPY AIR MASS TO
  !     *                         CANOPY MASS ONLY IF IDISP=0 OR IZREF=2.
  !     *                         ALSO, REPLACE LOGARITHMIC AVERAGING OF
  !     *                         ROUGHNESS HEIGHTS WITH BLENDING HEIGHT
  !     *                         AVERAGING.
  !     * JUL 31/02 - D.VERSEGHY. MOVE CALCULATION OF PSIGND AND FULL
  !     *                         CALCULATION OF FROOT INTO THIS ROUTINE
  !     *                         FROM energyBudgetPrep; REMOVE CALCULATION OF RCMIN.
  !     *                         SHORTENED CLASS3 COMMON BLOCK.
  !     * JUL 23/02 - D.VERSEGHY. MOVE ADDITION OF AIR TO CANOPY MASS
  !     *                         INTO THIS ROUTINE; SHORTENED CLASS4
  !     *                         COMMON BLOCK.
  !     * MAR 18/02 - D.VERSEGHY. MOVE CALCULATION OF SOIL PROPERTIES INTO
  !     *                         ROUTINE "soilProperties"; ALLOW FOR ASSIGNMENT
  !     *                         OF SPECIFIED TIME-VARYING VEGETATION
  !     *                         HEIGHT AND LEAF AREA INDEX.
  !     * SEP 19/00 - D.VERSEGHY. ADD CALCULATION OF VEGETATION-DEPENDENT
  !     *                         COEFFICIENTS FOR DETERMINATION OF STOMATAL
  !     *                         RESISTANCE.
  !     * APR 12/00 - D.VERSEGHY. RCMIN NOW VARIES WITH VEGETATION TYPE:
  !     *                         PASS IN BACKGROUND ARRAY "RCMINX".
  !     * DEC 16/99 - A.WU/D.VERSEGHY. ADD CALCULATION OF NEW LEAF DIMENSION
  !     *                              PARAMETER FOR REVISED CANOPY TURBULENT
  !     *                              TRANSFER FORMULATION.
  !     * NOV 16/98 - M.LAZARE.   "DLON" NOW PASSED IN AND USED DIRECTLY
  !     *                         (INSTEAD OF INFERRING FROM "LONSL" AND
  !     *                         "ILSL" WHICH USED TO BE PASSED) TO CALCULATE
  !     *                         GROWTH INDEX. THIS IS DONE TO MAKE THE PHYSICS
  !     *                         PLUG COMPATIBLE FOR USE WITH THE RCM WHICH
  !     *                         DOES NOT HAVE EQUALLY-SPACED LONGITUDES.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * OCT 11/96 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         BUG FIX: TO AVOID ROUND-OFF ERRORS,
  !     *                         SET CANOPY COVER EQUAL TO 1 IF THE
  !     *                         CALCULATED SUM OF FC AND FCS IS
  !     *                         VERY CLOSE TO 1.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     *                         ALSO CORRECT BUG IN CALCULATION OF
  !     *                         DEGLON, AND USE IDISP TO DETERMINE
  !     *                         METHOD OF CALCULATING DISP AND DISPS.
  !     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         VARIABLE SURFACE DETENTION CAPACITY
  !     *                         IMPLEMENTED.
  !     * AUG 16/95 - D.VERSEGHY. THREE NEW ARRAYS TO COMPLETE WATER
  !     *                         BALANCE DIAGNOSTICS.
  !     * NOV 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         RATIONALIZE CALCULATION OF RCMIN.
  !     * NOV 12/94 - D.VERSEGHY. FIX BUGS IN SENESCING LIMB OF CROP
  !     *                         GROWTH INDEX AND IN CANOPY MASS
  !     *                         CALCULATION.
  !     * MAY 06/93 - M.LAZARE/D.VERSEGHY. CLASS - VERSION 2.1.
  !     *                                  USE NEW "CANEXT" CANOPY
  !     *                                  EXTINCTION ARRAY TO DEFINE
  !     *                                  SKY-VIEW FACTORS. ALSO, CORRECT
  !     *                                  MINOR BUG WHERE HAD "IF (IN<=9)..."
  !     *                                  INSTEAD OF "IF (IN>9)...".
  !     * DEC 12/92 - M.LAZARE.   MODIFIED FOR MULTIPLE LATITUDES.
  !     * OCT 24/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CALCULATION OF LAND SURFACE CANOPY
  !     *                         PARAMETERS.
  !

  use classicParams,    only : zolnmoss, DELT, HCPW, HCPICE, &
                               HCPSND, SPHW, SPHICE, SPHVEG, SPHAIR, RHOW, RHOICE, &
                               PI, ZOLNG, ZOLNS, ZOLNI, ZORATG, GROWYR, ZORAT, &
                               CANEXT, XLEAF, classpfts, CL4CTEM, ctempfts

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL, IC, ICP1, IG, IDAY, IDISP, IZREF, IWF, IPAI, IHGT
  !
  integer :: I, J, K, IN, NL
  !
  !     * OUTPUT ARRAYS USED ELSEWHERE IN CLASS.
  !
  real, intent(out) :: FC    (ILG) !< Subarea fractional coverage of modelled area (X) [ ]
  real, intent(out) :: FG    (ILG) !< Subarea fractional coverage of modelled area (X) [ ]
  real, intent(out) :: FCS   (ILG) !< Subarea fractional coverage of modelled area (X) [ ]
  real, intent(out) :: FGS   (ILG) !< Subarea fractional coverage of modelled area (X) [ ]

  real, intent(out) :: PAICAN(ILG) !< Plant area index of canopy over bare ground (\f$\Lambda_p\f$) [ ]
  real, intent(out) :: PAICNS(ILG) !< Plant area index of canopy over snow (\f$\Lambda_p\f$) [ ]

  real, intent(out) :: FSVF  (ILG) !< Sky view factor for bare ground under canopy (\f$\chi\f$) [ ]
  real, intent(out) :: FSVFS (ILG) !< Sky view factor for snow under canopy (\f$\chi\f$) [ ]
  real, intent(out) :: FRAINC(ILG) !< Fractional coverage of canopy by liquid water over snow-free subarea [ ]
  real, intent(out) :: FSNOWC(ILG) !< Fractional coverage of canopy by frozen water over snow-free subarea [ ]
  real, intent(out) :: FRAICS(ILG) !< Fractional coverage of canopy by liquid water over snow-covered subarea [ ]
  real, intent(out) :: FSNOCS(ILG) !< Fractional coverage of canopy by frozen water over snow-covered subarea [ ]

  real, intent(out) :: RAICAN(ILG) !< Intercepted liquid water stored on canopy over bare ground (\f$W_l\f$) [\f$kg m^{-2}\f$]
  real, intent(out) :: RAICNS(ILG) !< Intercepted liquid water stored on canopy over snow (\f$W_l\f$) [\f$kg m^{-2}\f$]
  real, intent(out) :: SNOCAN(ILG) !< Intercepted frozen water stored on canopy over bare soil (\f$W_f\f$) [\f$kg m^{-2}\f$]
  real, intent(out) :: SNOCNS(ILG) !< Intercepted frozen water stored on canopy over snow (\f$W_f\f$) [\f$kg m^{-2}\f$]

  real, intent(out) :: DISP  (ILG) !< Displacement height of vegetation over bare ground (d) [m]
  real, intent(out) :: DISPS (ILG) !< Displacement height of vegetation over snow (d) [m]

  real, intent(out) :: ZOMLNC(ILG) !< Logarithm of roughness length for momentum of vegetation over bare ground [ ]
  real, intent(out) :: ZOMLCS(ILG) !< Logarithm of roughness length for momentum of vegetation over snow [ ]
  real, intent(out) :: ZOELNC(ILG) !< Logarithm of roughness length for heat of vegetation over bare ground [ ]
  real, intent(out) :: ZOELCS(ILG) !< Logarithm of roughness length for heat of vegetation over snow [ ]
  real, intent(out) :: ZOMLNG(ILG) !< Logarithm of roughness length for momentum of bare ground [ ]
  real, intent(out) :: ZOMLNS(ILG) !< Logarithm of roughness length for momentum of snow [ ]
  real, intent(out) :: ZOELNG(ILG) !< Logarithm of roughness length for heat of bare ground [ ]
  real, intent(out) :: ZOELNS(ILG) !< Logarithm of roughness length for heat of snow [ ]
  real, intent(out) :: RBCOEF(ILG) !< Parameter for calculation of leaf boundary resistance (\f$C_{rb}\f$)
  real, intent(out) :: CHCAP (ILG) !< Heat capacity of canopy over bare ground [\f$J m^{-2} K^{-1}\f$]
  real, intent(out) :: CHCAPS(ILG) !< Heat capacity of canopy over snow [\f$J m^{-2} K^{-1}\f$]
  real, intent(out) :: CMASSC(ILG) !< Mass of canopy over bare ground [\f$kg m^{-2}\f$]
  real, intent(out) :: CMASCS(ILG) !< Mass of canopy over snow [\f$kg m^{-2}\f$]
  real, intent(out) :: CWLCAP(ILG) !< Storage capacity of canopy over bare ground for liquid water (\f$W_{l,max}\f$) [\f$kg m^{-2}\f$]
  real, intent(out) :: CWFCAP(ILG) !< Storage capacity of canopy over bare ground for frozen water (\f$W_{f,max}\f$) [\f$kg m^{-2}\f$]
  real, intent(out) :: CWLCPS(ILG) !< Storage capacity of canopy over snow for liquid water (\f$W_{l,max}\f$) [\f$kg m^{-2}\f$]
  real, intent(out) :: CWFCPS(ILG) !< Storage capacity of canopy over snow for frozen water (\f$W_{f,max}\f$) [\f$kg m^{-2}\f$]

  real, intent(out) :: ZPLIMC(ILG) !< Maximum water ponding depth for ground under canopy [m]
  real, intent(out) :: ZPLIMG(ILG) !< Maximum water ponding depth for bare ground [m]
  real, intent(out) :: ZPLMCS(ILG) !< Maximum water ponding depth for ground under snow under canopy [m]
  real, intent(out) :: ZPLMGS(ILG) !< Maximum water ponding depth for ground under snow [m]

  real, intent(out) :: HTCC  (ILG) !< Diagnosed internal energy change of vegetation canopy
  !! due to conduction and/or change in mass [\f$W m^{-2}\f$]
  real, intent(out) :: HTCS  (ILG) !< Diagnosed internal energy change of snow pack
  !! due to conduction and/or change in mass [\f$W m^{-2}\f$]
  real, intent(out) :: WTRC  (ILG) !< Diagnosed residual water transferred off the vegetation canopy [\f$kg m^{-2} s^{-1}\f$]
  real, intent(out) :: WTRS  (ILG) !< Diagnosed residual water transferred into or out of the snow pack [\f$kg m^{-2} s^{-1}\f$]
  real, intent(out) :: WTRG  (ILG) !< Diagnosed residual water transferred into or out of the soil [\f$kg m^{-2} s^{-1}\f$]
  real, intent(out) :: CMAI  (ILG) !< Aggregated mass of vegetation canopy [\f$kg m^{-2}\f$]

  real, intent(out) :: FROOT (ILG,IG) !< Fraction of total transpiration contributed by soil layer over snow-free subarea  [  ]
  real, intent(out) :: FROOTS(ILG,IG) !< Fraction of total transpiration contributed by soil layer over snow-covered subarea  [  ]
  real, intent(out) :: HTC   (ILG,IG) !< Diagnosed internal energy change of soil layer
  !! due to conduction and/or change in mass [\f$W m^{-2}\f$]

  !
  !     * OUTPUT ARRAYS ONLY USED ELSEWHERE IN radiationDriver.
  !
  real, intent(inout) :: PAI   (ILG,IC) !< Plant area index of vegetation category over bare ground [ ]
  real, intent(inout) :: PAIS  (ILG,IC) !< Plant area index of vegetation category over snow [ ]
  real, intent(inout) :: AIL   (ILG,IC) !< Leaf area index of vegetation category over bare ground [ ]
  real, intent(inout) :: FCAN  (ILG,IC) !< Fractional coverage of vegetation category over bare ground (\f$X_i\f$) [ ]
  real, intent(inout) :: FCANS (ILG,IC) !< Fractional coverage of vegetation category over snow (\f$X_i\f$) [ ]
  real, intent(inout) :: PSIGND(ILG)    !< Minimum liquid moisture suction in soil layers [m]

  !
  !     * INPUT ARRAYS.
  !
  integer, intent(in) :: ipeatland (ilg) !< Peatland flag: 0 = not a peatland,1= bog,2 = fen
  real, intent(in) :: FCANMX(ILG,ICP1) !< Maximum fractional coverage of modelled area by vegetation category [ ]
  real, intent(inout) :: ZOLN  (ILG,ICP1) !< Natural logarithm of maximum roughness length of vegetation category [ ]
  real, intent(in) :: PAIMAX(ILG,IC)   !< Maximum plant area index of vegetation category [ ]
  real, intent(in) :: PAIMIN(ILG,IC)   !< Minimum plant area index of vegetation category [ ]
  real, intent(in) :: CWGTMX(ILG,IC)   !< Maximum canopy mass for vegetation category [\f$kg m^{-2}\f$]
  real, intent(in) :: ZRTMAX(ILG,IC)   !< Maximum rooting depth of vegetation category [m]
  real, intent(in) :: PAIDAT(ILG,IC)   !< Optional user-specified value of plant area indices of
  !! vegetation categories to override CLASS-calculated values [ ]
  real, intent(in) :: HGTDAT(ILG,IC)   !< Optional user-specified values of height of
  !! vegetation categories to override CLASS-calculated values [m]
  real, intent(inout) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil layers (\f$\theta\f$ l) [\f$m^3 m^{-3}\f$]
  real, intent(inout) :: THICE (ILG,IG)   !< Frozen water content of soil layers under vegetation [\f$m^3 m^{-3}\f$]
  real, intent(inout) :: TBAR  (ILG,IG)   !< Temperature of soil layers [K]

  real, intent(inout) :: RCAN  (ILG) !< Intercepted liquid water stored on canopy (\f$W_l\f$) [\f$kg m^{-2}\f$]
  real, intent(inout) :: SNCAN (ILG) !< Intercepted frozen water stored on canopy (\f$W_f\f$) [\f$kg m^{-2}\f$]
  real, intent(inout) :: TCAN  (ILG) !< Vegetation canopy temperature [K]
  real, intent(in) :: GROWTH(ILG) !< Vegetation growth index [ ]
  real, intent(inout) :: ZSNOW (ILG) !< Depth of snow pack (\f$z_s\f$) [m]
  real, intent(inout) :: TSNOW (ILG) !< Snowpack temperature [K]
  real, intent(inout) :: FSNOW (ILG) !< Diagnosed fractional snow coverage [ ]
  real, intent(in) :: RHOSNO(ILG) !< Density of snow (\f$ s\f$) [\f$kg m^{-3}\f$]
  real, intent(inout) :: SNO   (ILG) !< Mass of snow pack (\f$W_s\f$) [\f$kg m^{-2}\f$]
  real, intent(in) :: TA    (ILG) !< Air temperature at reference height [K]
  real, intent(in) :: RHOAIR(ILG) !< Density of air [\f$kg m^{-3}\f$]
  real, intent(in) :: DLON  (ILG) !< Longitude of grid cell (east of Greenwich) [degrees]
  real, intent(in) :: Z0ORO (ILG) !< Orographic roughness length [m]
  real, intent(in) :: ZBLEND(ILG) !< Atmospheric blending height for surface roughness length averaging (\f$z_b\f$) [m]
  real, intent(in) :: RHOSNI(ILG) !< Density of fresh snow (\f$\rho\f$ s,f) [\f$kg m^{-3}\f$]
  real, intent(in) :: ZPLMG0(ILG) !< Maximum water ponding depth for snow-free subareas
  !! (user-specified when running MESH code) [m]
  real, intent(in) :: ZPLMS0(ILG) !< Maximum water ponding depth for snow-covered subareas
  !! (user-specified when running MESH code) [m]
  real, intent(in) :: RADJ  (ILG) !< Latitude of grid cell (positive north of equator) [rad]

  !
  !     * SOIL PROPERTY ARRAYS.
  !
  real, intent(in) :: DELZW (ILG,IG) !< Permeable thickness of soil layer [m]
  real, intent(in) :: ZBOTW (ILG,IG) !< Depth to permeable bottom of soil layer [m]
  real, intent(in) :: THPOR (ILG,IG) !< Pore volume in soil layer (\f$\theta\f$ p) [\f$m^3 m^{-3}\f$]
  real, intent(in) :: THLMIN(ILG,IG) !< Residual soil liquid water content remaining after freezing or evaporation [\f$m^3 m^{-3}\f$]
  real, intent(in) :: PSISAT(ILG,IG) !< Soil moisture suction at saturation (\f$\Psi\f$ sat) [m]
  real, intent(in) :: BI    (ILG,IG) !< Clapp and Hornberger empirical "b" parameter [ ]
  real, intent(in) :: PSIWLT(ILG,IG) !< Soil moisture suction at wilting point (\f$\Psi\f$ w) [m]
  real, intent(in) :: HCPS  (ILG,IG) !< Volumetric heat capacity of soil particles [\f$J m^{-3}\f$]
  real, intent(in) :: DELZ  (IG)     !< Soil layer thickness [m]
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  !     * WORK ARRAYS NOT USED ELSEWHERE IN radiationDriver.
  !
  real, intent(inout) :: RMAT (ILG,IC,IG), H(ILG,IC), HS(ILG,IC), &
                         CWCPAV(ILG), GROWA(ILG), GROWN(ILG), &
                         GROWB(ILG), RRESID(ILG), SRESID(ILG), &
                         FRTOT(ILG), FRTOTS(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: DAY, GROWG, FSUM, SNOI, ZSNADD, THSUM, THICEI, THLIQI, ZROOT, &
          ZROOTG, FCOEFF, PSII, LZ0ORO, THR_LAI, PSIRAT, buriedcrp, &
          buriedgras, sumfcanmx, temp1
  !
  !     * CTEM-RELATED FIELDS.
  !
  real, intent(in)   :: AILC(ILG,IC)       !<
  real, intent(in)   :: PAIC(ILG,IC)       !<
  real, intent(in)   :: AILCG(ILG,ICTEM)   !< GREEN LAI FOR USE WITH PHTSYN SUBROUTINE
  real, intent(out)  :: AILCGS(ILG,ICTEM)  !< GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
  real, intent(in)   :: RMATC(ILG,IC,IG)   !<
  real, intent(in)   :: FCANCMX(ILG,ICTEM) !<
  real, intent(out)  :: FCANC(ILG,ICTEM)   !< FRACTION OF CANOPY OVER GROUND FOR CTEM's 9 PFTs
  real, intent(out)  :: FCANCS(ILG,ICTEM)  !< FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
  real, intent(in)   :: ZOLNC(ILG,IC)      !<
  real, intent(in)   :: CMASVEGC(ILG,IC)   !<
  real, intent(in)   :: SLAIC(ILG,IC)      !<
  !
  !     * NOL2PFTS - NUMBER OF LEVEL 2 CTEM PFTs
  !     * SEE BIO2STR SUBROUTINE FOR EXPLANATION OF OTHER CTEM VARIABLES

  !     * INTERNAL WORK FIELD.
  !
  integer, intent(in) :: ICTEM, NOL2PFTS(IC)   !<
  real  :: SFCANCMX(ILG,IC)
  !
  logical, intent(in) :: ctem_on

  integer :: M, N, K1, K2
  !
  !-----------------------------------------------------------------------
  !
  !     * INITIALIZE DIAGNOSTIC AND OTHER ARRAYS.
  !
  do I = IL1,IL2 ! loop 100
    HTCC(I) = 0.0
    HTCS(I) = 0.0
    do J = 1,IG ! loop 50
      HTC(I,J) = 0.0
    end do ! loop 50
    WTRC(I) = 0.0
    WTRS(I) = 0.0
    WTRG(I) = 0.0
    FRTOT(I) = 0.0
    FRTOTS(I) = 0.0
    DISP  (I) = 0.
    ZOMLNC(I) = 0.
    ZOELNC(I) = 1.
    DISPS (I) = 0.
    ZOMLCS(I) = 0.
    ZOELCS(I) = 1.
    ZOMLNG(I) = 0.
    ZOELNG(I) = 0.
    ZOMLNS(I) = 0.
    ZOELNS(I) = 0.
    CMASSC(I) = 0.
    CMASCS(I) = 0.
    PSIGND(I) = 1.0E+5
  end do ! loop 100
  !
  !     * DETERMINE GROWTH INDEX FOR CROPS (VEGETATION TYPE 3).
  !     * MUST USE UN-GATHERED LONGITUDES TO COMPUTE ACTUAL LONGITUDE/
  !     * LATITUDE VALUES.
  !
  DAY = real(IDAY)
  !
  !     * FOR CTEM, CROP GROWTH IS BUILT IN, SO GROWA=1.
  !
  if (.not. ctem_on) then
    !>
    !! In the 120 loop, the growth index for crops, GROWA, is calculated (if CLASS is not being run coupled to
    !! CTEM). This is done by referring to the three-dimensional array GROWYR, which contains values
    !! corresponding to the four Julian days of the year on which crops are planted, on which they reach
    !! maturity, on which harvesting begins, and on which the harvest is complete, for each ten-degree latitude
    !! half-circle in each hemisphere. These are generic, average dates, approximated using information gleaned
    !! from annual UN FAO (Food and Agriculture Organization) reports. (In the tropics, it is assumed that
    !! areas classified as agricultural are constantly under cultivation, so all four values are set to zero.)
    !!
    !! First, the latitude of the modelled area is converted from a value in radians to a value from 1 to 18, IN,
    !! corresponding to the index of the latitude circle (1 for latitudes \f$80-90^o S\f$, 18 for latitudes \f$80-90^o N\f$). Then
    !! the hemisphere index, NL, is set to 1 for the Eastern Hemisphere, and 2 for the Western Hemisphere. If
    !! the planting date for the modelled area is zero (indicating a location in the tropics), GROWA is set to 1.
    !! Otherwise, GROWA is set to 1 if the day of the year lies between the maturity date and the start of the
    !! harvest, and to zero if the day of the year lies between the end of the harvest and the planting date. For
    !! dates in between, the value of GROWA is interpolated between 0 and 1. Checks are performed at the
    !! end to ensure that GROWA is not less than 0 or greater than 1. If the calculated value of GROWA is
    !! vanishingly small, it is set to zero.
    !!
    do I = IL1,IL2 ! loop 120
      IN = INT( (RADJ(I) + PI / 2.0) * 18.0 / PI) + 1
      if (DLON(I) > 190. .and. DLON(I) < 330.) then
        NL = 2
      else
        NL = 1
      end if
      if (GROWYR(IN,1,NL) < 0.1) then
        GROWA(I) = 1.0
      else
        if (IN > 9) then
          if (DAY >= GROWYR(IN,2,NL) .and. DAY < GROWYR(IN,3,NL)) &
              GROWA(I) = 1.0
          if (DAY >= GROWYR(IN,4,NL) .or. DAY < GROWYR(IN,1,NL)) &
              GROWA(I) = 0.0
        else
          if (DAY >= GROWYR(IN,2,NL) .or. DAY < GROWYR(IN,3,NL)) &
              GROWA(I) = 1.0
          if (DAY >= GROWYR(IN,4,NL) .and. DAY < GROWYR(IN,1,NL)) &
              GROWA(I) = 0.0
        end if
        if (DAY >= GROWYR(IN,1,NL) .and. DAY < GROWYR(IN,2,NL)) &
            GROWA(I) = (DAY - GROWYR(IN,1,NL)) / (GROWYR(IN,2,NL) - &
            GROWYR(IN,1,NL))
        if (DAY >= GROWYR(IN,3,NL) .and. DAY < GROWYR(IN,4,NL)) &
            GROWA(I) = (GROWYR(IN,4,NL) - DAY) / (GROWYR(IN,4,NL) - &
            GROWYR(IN,3,NL))
        GROWA(I) = MAX(0.0,MIN(GROWA(I),1.0))
        if (GROWA(I) < 1.0E-5) GROWA(I) = 0.0
      end if
    end do ! loop 120
  else
    do I = IL1,IL2
      GROWA(I) = 1.
    end do
  end if
  !
  !     * DETERMINE GROWTH INDICES FOR NEEDLELEAF TREES, BROADLEAF
  !     * TREES AND GRASS (VEGETATION TYPES 1, 2 AND 4); CALCULATE
  !     * VEGETATION HEIGHT, CORRECTED FOR GROWTH STAGE FOR CROPS
  !     * AND FOR SNOW COVER FOR CROPS AND GRASS; CALCULATE CURRENT
  !     * LEAF AREA INDEX FOR FOUR VEGETATION TYPES.
  !

  !>
  !! In the 150 loop the other three growth indices are evaluated, as well as the vegetation heights and plant
  !! area indices for the four vegetation categories over snow-covered and snow-free ground. The
  !! background growth index for trees, GROWTH, is evaluated separately in subroutine classGrowthIndex. It varies
  !! from a value of 0 for dormant or leafless periods to 1 for fully-leafed periods, with a sixty-day transition
  !! between the two. When senescence begins, it is set instantaneously to -1 and thereafter increases over a
  !! sixty-day period back to 0. (The onset of spring budburst and fall senescence are triggered by near-zero
  !! values of the air temperature and the first soil layer temperature.) For needleleaf trees, the growth index
  !! GROWN is simply set to the absolute value of GROWTH. For broadleaf trees, the transition period is
  !! assumed to last thirty days instead of sixty, and so the growth index GROWB is set to the absolute value
  !! of double the value of GROWTH, with upper and lower limits of 1 and 0. Finally, the growth index of
  !! grasses is set to 1 all year round.
  !!
  do I = IL1,IL2 ! loop 150

    GROWN(I) = ABS(GROWTH(I))
    if (GROWTH(I) > 0.0) then
      GROWB(I) = MIN(1.0,GROWTH(I) * 2.0)
    else
      GROWB(I) = MAX(0.0,(ABS(GROWTH(I)) * 2.0 - 1.0))
    end if
    GROWG = 1.0

    !    IF USING CTEM's STRUCTURAL ATTRIBUTES OVERWRITE ZOLN
    if (ctem_on) then
      do J = 1,IC
        ZOLN(I,J) = ZOLNC(I,J)
      end do
    end if

    !>
    !! A branch in the code occurs next, depending on the value of the flag IHGT. If IHGT=0, the values of
    !! vegetation height calculated by CLASS are to be used. For trees and grass, the vegetation height under
    !! snow-free conditions is assumed to be constant year-round, and is calculated as 10 times the exponential
    !! of ZOLN, the logarithm of the maximum vegetation roughness length. For crops, this maximum height
    !! is multiplied by GROWA. If IHGT=1, vegetation heights specified by the user are utilized instead. This
    !! height H for each of the four vegetation categories is used to calculate the height HS over snow-covered
    !! areas. For needleleaf and broadleaf trees, HS is set to H. For crops and grass, HS is calculated by
    !! subtracting the snow depth ZSNOW from H, to account for the burying of short vegetation by snow.
    !!
    if (IHGT == 0) then ! Vegetation height from CLASS
      do J = 1,IC
        select case (classpfts(J))
        case ('NdlTr' , 'BdlTr', 'Grass', 'BdlSh')  ! <tree, grass, and shrub
          H(I,J) = 10.0 * EXP(ZOLN(I,J))
        case ('Crops') ! <Crops
          H(I,J) = 10.0 * EXP(ZOLN(I,J)) * GROWA(I)
        case default
          print * ,'Unknown PFT in calcLandSurfParams ',classpfts(J)
          call errorHandler('calcLandSurfParams', - 1)
        end select
      end do
    else ! Vegetation height read in from external file
      do J = 1,IC
        H(I,J) = HGTDAT(I,J) ! <Vegetation heights are user specified.
      end do
    end if

    !> Now account for burying under snow:
    do J = 1,IC
      select case (classpfts(J))
      case ('NdlTr' , 'BdlTr', 'BdlSh')   ! tree, shrub (shrub with no bending or burying here)
        HS(I,J) = H(I,J)
      case ('Crops','Grass') ! Crops,Grass
        HS(I,J) = MAX(H(I,J) - ZSNOW(I),1.0E-3)
      case default
        print * ,'Unknown PFT in calcLandSurfParams ',classpfts(J)
        call errorHandler('calcLandSurfParams', - 2)
      end select
    end do
    !>
    !! If CLASS is being run uncoupled to CTEM, a second branch now occurs, depending on the value of the
    !! flag IPAI. If IPAI=0, the values of plant area index calculated by CLASS are to be used. For all four
    !! vegetation categories, the plant area index over snow-free ground, PAI, is determined by interpolating
    !! between the annual maximum and minimum plant area indices using the growth index. If IPAI=1, plant
    !! area index values specified by the user are utilized instead. For trees, the plant area index over snow-
    !! covered ground, PAIS, is set to PAI. For crops and grass, if H>0, PAIS is set to PAI scaled by the ratio
    !! of HS/H; otherwise, it is set to zero. Lastly, the leaf area indices for the four vegetation categories over
    !! snow-free ground, AIL, are determined from the PAI values. For needleleaf trees, AIL is estimated as
    !! 0.90 PAI; for broadleaf trees it is estimated as the excess PAI over the annual minimum value. For crops
    !! and grass AIL is assumed to be equal to PAI. (If CLASS is being run coupled to CTEM, the CTEM-
    !! generated values of PAI and AIL are used instead.)
    !!
    if (IPAI == 0) then
      !             USE CTEM GENERATED PAI OR CLASS' OWN SPECIFIED PAI
      if (ctem_on) then
        do J = 1,IC
          PAI(I,J) = PAIC(I,J)
        end do
      else ! CLASS (physics) only
        do J = 1,IC
          select case (classpfts(J))
          case ('BdlTr','BdlSh')   ! <Broadleaf tree and shrub
            PAI(I,J) = PAIMIN(I,J) + GROWB(I) * (PAIMAX(I,J) - &
                       PAIMIN(I,J))
          case ('NdlTr') ! <Needleleaf tree
            PAI(I,J) = PAIMIN(I,J) + GROWN(I) * &
                       (PAIMAX(I,J) - PAIMIN(I,J))
          case ('Crops') ! <Crops
            PAI(I,J) = PAIMIN(I,J) + GROWA(I) * (PAIMAX(I,J) - &
                       PAIMIN(I,J))
          case ('Grass ') ! <Grass
            PAI(I,J) = PAIMIN(I,J) + GROWG   * (PAIMAX(I,J) - &
                       PAIMIN(I,J))
          case default
            ! Error checking,if no known CLASS PFT given,bail.
            print * ,'calcLandSurfParams says: Unknown CLASS pft=>',classpfts(j)
            call errorHandler('calcLandSurfParams', - 3)
          end select
        end do
      end if ! CTEM on/off
    else ! IPAI = 1,so
      ! Use a read-in time series of PAI.
      do J = 1,IC
        PAI(I,J) = PAIDAT(I,J)
      end do
    end if

    !        Now account for burying due to snow:
    do J = 1,IC
      select case (classpfts(J))
      case ('NdlTr' , 'BdlTr', 'BdlSh')   ! <trees and shrub (NO SNOW BURIAL)
        PAIS(I,J) = PAI(I,J)
      case ('Crops','Grass') ! <Crops,Grass
        if (H(I,J) > 0.0) then
          PAIS(I,J) = PAI(I,J) * HS(I,J) / H(I,J)
        else
          PAIS(I,J) = 0.0
        end if
      case default
        ! Error checking,if no known CLASS PFT given,bail.
        print * ,'calcLandSurfParams says: Unknown CLASS pft => ',classpfts(j)
        call errorHandler('calcLandSurfParams', - 4)
      end select
    end do

    if (ctem_on) then
      do J = 1,IC
        AIL(I,J) = MAX(AILC(I,J),SLAIC(I,J))
      end do
    else ! CLASS (physics) only
      do J = 1,IC
        select case (classpfts(J))
        case ('Crops','Grass')   ! <Crops and Grass
          AIL(I,J) = PAI(I,J)
        case ('NdlTr') ! <Needleleaf
          AIL(I,J) = PAI(I,J) * 0.90    
        case ('BdlTr','BdlSh') ! <Broadleaf
          AIL(I,J) = MAX((PAI(I,J) - PAIMIN(I,J)),0.0)
        case default
          ! Error checking,if no known CLASS PFT given,bail.
          print * ,'calcLandSurfParams says: Unknown CLASS pft => ',classpfts(j)
          call errorHandler('calcLandSurfParams', - 5)
        end select
      end do
    end if ! CTEM on/off

    !         ESTIMATE GREEN LAI FOR CANOPY OVER SNOW FRACTION FOR CTEM's
    !         9 PFTs, JUST LIKE CLASS DOES.
    !
    if (ctem_on) then
      do J = 1,ICTEM
        select case (ctempfts(J))
        case ('NdlEvgTr','NdlDcdTr','BdlEvgTr','BdlDCoTr', &
                'BdlDDrTr','BdlDCoSh','BdlEvgSh') ! <Tree and shrub
          AILCGS(I,J) = AILCG(I,J)
        case ('CropC3  ','CropC4  ','GrassC3 ', &
                'GrassC4 ','Sedge   ')
          if (H(I,CL4CTEM(J)) > 0.0) then
            AILCGS(I,J) = AILCG(I,J) * HS(I,CL4CTEM(J)) / &
                          H(I,CL4CTEM(J))
          else
            AILCGS(I,J) = 0.0
          end if
        case default
          ! Error checking, if no known CTEM PFT given, bail.
          print * ,'calcLandSurfParams says: Unknown CTEM pft => ',ctempfts(j)
          call errorHandler('calcLandSurfParams', - 6)
        end select
      end do
    end if

  end do ! loop 150
  !
  !     * ADJUST FRACTIONAL COVERAGE OF GRID CELL FOR CROPS AND
  !     * GRASS IF LAI FALLS BELOW A SET THRESHOLD VALUE DUE TO
  !     * GROWTH STAGE OR SNOW COVER; RESET LAI TO THE THRESHOLD
  !     * VALUE; CALCULATE RESULTANT GRID CELL COVERAGE BY CANOPY,
  !     * BARE GROUND,cCANOPY OVER SNOW AND SNOW OVER BARE GROUND.
  !     *
  !     * ALSO CALCULATE SURFACE DETENTION CAPACITY FOR FOUR
  !     * GRID CELL SUBAREAS BASED ON VALUES SUPPLIED BY
  !     * U. OF WATERLOO:
  !     *        IMPERMEABLE SURFACES: 0.001 M.
  !     *        BARE SOIL:            0.002 M.
  !     *        LOW VEGETATION:       0.003 M.
  !     *        FOREST:               0.01  M.
  !
  THR_LAI = 1.0  
  !>
  !! In the 175 loop, the fractional coverage of the modelled area by each of the four vegetation categories is
  !! calculated, for snow-free (FCAN) and snow-covered ground (FCANS). For needleleaf and broadleaf
  !! trees, FCAN is set to the maximum coverage FCANMX of the vegetation category, scaled by the snow-
  !! free fraction of the modelled area, 1-FSNOW. For crops and grass, this calculation is modified for cases
  !! where the plant area index has been calculated as falling below a threshold value owing to growth stage or
  !! burying by snow. (If CLASS is being run coupled to CTEM, this threshold value is set to 0.05; otherwise
  !! it is set to 1.) In such cases the vegetation coverage is assumed to become discontinuous, and so an
  !! additional multiplication by PAI is performed to produce a reduced value of FCAN, and PAI is reset to
  !! the threshold value. An identical procedure is followed to determine the FCANS values.
  !!
  !! The areal :: fractions of each of the four CLASS subareas, vegetation over bare soil (FC), bare soil (FG),
  !! vegetation over snow (FCS) and snow (FGS) are then calculated, using the FCAN and FCANS values and
  !! FSNOW. Checks are carried out, and adjustments performed if necessary, to ensure that none of the
  !! four subareas is vanishingly small. The values of FSNOW and of the four FCANs and FCANSs are
  !! recalculated accordingly. Finally, checks are carried out to ensure that each of the four subareas is greater
  !! than zero, and that their sum is unity. If this is not the case, a call to abort is performed.
  !! In the last part of the 175 loop, the limiting ponding depth for surface water is determined for each of the
  !! four subareas. If the flag IWF is zero, indicating that lateral flow of water within the soil is to be
  !! neglected, these values are assigned as follows. If the index ISAND of the first soil layer is -3 or -4,
  !! indicating a rock surface or an ice sheet, the bare soil ponding limit, ZPLIMG, is set to 1 mm; otherwise
  !! ZPLIMG is set to 2 mm. If the fractional area of snow on bare soil is greater than zero, the subarea
  !! ponding limit ZPLMGS is set to the weighted average of ZPLIMG over the areas where snow has not
  !! buried vegetation and where it has buried crops, and to 3 mm over areas where it has buried grass
  !! otherwise to zero. If the fractional area of canopy over bare soil is greater than zero, the subarea ponding
  !! depth ZPLIMC is set to the weighted average of 1 cm under trees and 3 mm under crops and grass
  !! otherwise to zero. If the fractional area of canopy over snow is greater than zero, the subarea ponding
  !! depth ZPLMCS is also currently set to the weighted average of 1 cm under trees and 3 mm under crops
  !! and grass; otherwise to zero. Finally, if the flag IWF is greater than zero, indicating that lateral flow of soil
  !! water is being modelled, externally derived user-specified values of the ponding limit for the four subareas
  !! are assigned.
  !!

  do I = IL1,IL2 ! loop 175
    do J = 1,IC
      FCAN(I,J) = FCANMX(I,J) * (1.0 - FSNOW(I))
      if (FCAN(I,J) < 1.0E-5) FCAN(I,J) = 0.0
      select case (classpfts(J))
      case ('Crops','Grass')
        if (PAI(I,J) < THR_LAI) then
          FCAN(I,J) = FCANMX(I,J) * (1.0 - FSNOW(I)) * PAI(I,J)
          PAI (I,J) = THR_LAI
          if (FCAN(I,J) < 1.0E-5) FCAN(I,J) = 0.0
        end if
      case ('NdlTr' , 'BdlTr', 'BdlSh')
        if (PAI(I,J) < THR_LAI) then ! FLAG HACK
          FCAN(I,J) = FCANMX(I,J) * (1.0 - FSNOW(I)) * PAI(I,J)
          PAI (I,J) = THR_LAI
          if (FCAN(I,J) < 1.0E-5) FCAN(I,J) = 0.0
        end if
        ! Do nothing.
      end select
    end do
    do J = 1,IC
      FCANS(I,J) = FCANMX(I,J) * FSNOW(I)
      if (FCANS(I,J) < 1.0E-5) FCANS(I,J) = 0.0
      select case (classpfts(J))
      case ('Crops','Grass')
        if (PAIS(I,J) < THR_LAI) then
          FCANS(I,J) = FCANMX(I,J) * FSNOW(I) * PAIS(I,J)
          PAIS (I,J) = THR_LAI
          if (FCANS(I,J) < 1.0E-5) FCANS(I,J) = 0.0
        end if
      case ('NdlTr' , 'BdlTr', 'BdlSh')
        if (PAIS(I,J) < THR_LAI) then ! FLAG HACK
          FCANS(I,J) = FCANMX(I,J) * FSNOW(I) * PAIS(I,J)
          PAIS (I,J) = THR_LAI
          if (FCANS(I,J) < 1.0E-5) FCANS(I,J) = 0.0
        end if
        ! Do nothing
      end select
    end do

    FC (I) = SUM(FCAN(I,:))
    FCS(I) = SUM(FCANS(I,:))
    FG (I) = 1.0 - FSNOW(I) - FC(I)
    FGS(I) = FSNOW(I) - FCS(I)
    if (ABS(1.0 - FCS(I) - FC(I)) < 8.0E-5) then
      if (FCS(I) < 1.0E-5) then
        FSNOW (I) = 0.0
      else if (FC(I) < 1.0E-5) then
        FSNOW(I) = 1.0
      end if
      if (FCS(I) > 0.) then
        do J = 1,IC
          FCANS(I,J) = FCANS(I,J) * FSNOW(I) / FCS(I)
        end do
      end if
      if (FC(I) > 0.) then
        do J = 1,IC
          FCAN(I,J) = FCAN(I,J) * (1.0 - FSNOW(I)) / FC(I)
        end do
      end if
      FCS(I) = MIN(FSNOW(I),1.0)
      FC(I) = 1.0 - FCS(I)
      FGS(I) = 0.0
      FG(I) = 0.0
    end if
    FC (I) = MAX(FC (I),0.0)
    FG (I) = MAX(FG (I),0.0)
    FCS(I) = MAX(FCS(I),0.0)
    FGS(I) = MAX(FGS(I),0.0)
    FSUM = (FCS(I) + FGS(I) + FC(I) + FG(I))
    FC (I) = FC (I) / FSUM
    FG (I) = FG (I) / FSUM
    FCS(I) = FCS(I) / FSUM
    FGS(I) = FGS(I) / FSUM

    if (ABS(1.0 - FCS(I) - FGS(I) - FC(I) - FG(I)) > 1.0E-5) &
        call errorHandler('calcLandSurfParams',1)
    !
    if (IWF == 0) then
      if (ISAND(I,1) == - 4) then
        ZPLIMG(I) = 0.001  
      else if (ISAND(I,1) == - 3) then
        ZPLIMG(I) = 0.001
      else
        ZPLIMG(I) = 0.002
      end if

      ZPLMGS(I) = 0.0
      if (FGS(I) > 0.0) then
        sumfcanmx = 0.
        buriedgras = 0.
        buriedcrp = 0.
        do J = 1,IC
          sumfcanmx = sumfcanmx + FCANMX(I,J)
          select case (classpfts(J))
          case ('Crops')
            buriedcrp = FSNOW(I) * FCANMX(I,J) - FCANS(I,J)
          case ('Grass')
            buriedgras = FSNOW(I) * FCANMX(I,J) - FCANS(I,J)
          case ('NdlTr' , 'BdlTr', 'BdlSh')
            ! Assume not buried so do nothing.
          case default
            ! Error checking, if no known CLASS PFT given, bail.
            print * ,'calcLandSurfParams says: Unknown CLASS pft => ' &
                         ,classpfts(j)
            call errorHandler('calcLandSurfParams', - 7)
          end select
        end do
        ZPLMGS(I) = (ZPLIMG(I) * FSNOW(I) * (1.0 - sumfcanmx) &
                    + ZPLIMG(I) * buriedcrp + 0.003 &
                    * buriedgras) &
                    / FGS(I)
      end if
      !
      ZPLIMC(I) = 0.0
      if (FC(I) > 0.0) then
        do J = 1,IC
          select case (classpfts(J))
          case ('NdlTr' , 'BdlTr', 'BdlSh') ! assume trees and shrubs don't differ here.
            ZPLIMC(I) = ZPLIMC(I) + 0.01 * FCAN(I,J)
          case ('Crops','Grass')
            ZPLIMC(I) = ZPLIMC(I) + 0.003 * FCAN(I,J)
          case default
            ! Error checking, if no known CLASS PFT given, bail.
            print * ,'calcLandSurfParams says: Unknown CLASS pft => ' &
                        ,classpfts(j)
            call errorHandler('calcLandSurfParams', - 8)
          end select
        end do
        ZPLIMC(I) = ZPLIMC(I) / FC(I)
      end if
      !
      ZPLMCS(I) = 0.0
      if (FCS(I) > 0.0) then
        do J = 1,IC
          select case (classpfts(J))
          case ('NdlTr' , 'BdlTr', 'BdlSh') ! assume trees and shrubs don't differ here.
            ZPLMCS(I) = ZPLMCS(I) + 0.01 * FCANS(I,J)
          case ('Crops','Grass')
            ZPLMCS(I) = ZPLMCS(I) + 0.003 * FCANS(I,J)
          case default
            ! Error checking, if no known CLASS PFT given, bail.
            print * ,'calcLandSurfParams says: Unknown CLASS pft => ' &
                    ,classpfts(j)
            call errorHandler('calcLandSurfParams', - 9)
          end select
        end do
        ZPLMCS(I) = ZPLMCS(I) / FCS(I)
      end if
    else
      ZPLMCS(I) = ZPLMS0(I)
      ZPLMGS(I) = ZPLMS0(I)
      ZPLIMC(I) = ZPLMG0(I)
      ZPLIMG(I) = ZPLMG0(I)
    end if
  end do ! loop 175
  !
  !     * PARTITION INTERCEPTED LIQUID AND FROZEN MOISTURE BETWEEN
  !     * CANOPY OVERLYING BARE GROUND AND CANOPY OVERLYING SNOW,
  !     * USING DIFFERENT EFFECTIVE LEAF AREAS FOR EACH.  ADD
  !     * RESIDUAL TO SOIL MOISTURE OR SNOW (IF PRESENT); CALCULATE
  !     * RELATIVE FRACTIONS OF LIQUID AND FROZEN INTERCEPTED
  !     * MOISTURE ON CANOPY.
  !

  !>
  !! In loop 200, calculations are done related to the interception of water on vegetation. First, the plant area
  !! indices of the composite vegetation canopy over the snow-free and snow-covered subareas are calculated
  !! as weighted averages of the plant area indices of the four vegetation categories over each subarea. The
  !! liquid water interception capacity \f$W_{l, max}\f$ on each of the two subareas is calculated as
  !! \f$W_{l, max} = 0.20 \Lambda_p\f$
  !! where \f$\Lambda_p\f$ is the plant area index of the composite canopy. This simple relation has been found to work
  !! well for a wide range of canopy types and precipitation events (see Verseghy et al, 1993). If either the
  !! average amount of liquid water on the canopy, RCAN, or the total cancpy coverage, FC+FCS, is less than
  !! a small threshold value, the value of RCAN is stored in a residual water array RRESID, and RCAN is set
  !! to zero. Next the intercepted liquid water is partitioned between FC and FCS. First RCAN is re-
  !! evaluated as an average value over the canopy-covered area only, rather than over the whole modelled
  !! area. Then the intercepted liquid water amounts on vegetation over snow-free (RAICAN) and snow-
  !! covered areas (RAICNS) are calculated by making use of the relations
  !!
  !! \f$W_{L, 0} / \Lambda_{p, 0} = W_{L, s} / \Lambda_{p, s}\f$ and
  !!
  !! \f$W_l (X_0 + X_s) = W_{l, 0} X_0 + W_{l, s} X_s\f$
  !!
  !! where \f$W_l\f$ is the liquid water on the canopy, \f$X\f$ is the fractional area, and the subscripts \f$0\f$ and \f$s\f$ refer to
  !! snow-free and snow-covered areas respectively.
  !!
  !! For snow interception on the canopy, a modified calculation of the plant area indices \f$\Lambda_{p, 0}\f$ and \f$\Lambda_{p, s}\f$ is
  !! performed, assigning a weight of 0.7 to the plant area index of needleleaf trees, to account for the effect
  !! of needle clumping. The interception capacity for snow, \f$W_{f, max}\f$, is calculated following Bartlett et al.
  !! (2006) \cite Bartlett2006-xp, using a relation developed by Schmidt and Gluns (1991):
  !! \f$W_{f, max} = 6.0 \Lambda_p [0.27 + 46.0 \rho_{s, f} ]\f$
  !! where \f$\rho_{s, f}\f$ is the density of fresh snow. As was done for the intercepted liquid water, if either the average
  !! amount of snow on the canopy, SNCAN, or the total cancpy coverage is less than a small threshold value,
  !! the value of SNCAN is stored in a residual water array SRESID, and SNCAN is set to zero. Next the
  !! intercepted snow is partitioned between FC and FCS. First SNCAN is recalculated as an average over the
  !! canopy-covered area only, rather than over the whole modelled area. Then the intercepted snow amounts
  !! on vegetation over snow-free (SNOCAN) and snow-covered areas (SNOCNS) are calculated in the same
  !! way as for RAICAN and RAICNS.
  !!
  !! The fractional canopy coverages of snow and liquid water are calculated as the ratio of the intercepted
  !! snow or liquid water to their respective interception capacities. Separate values are determined for the
  !! snow-covered (FSNOCS, FRAICS) and snow-free (FSNOWC, FRAINC) subareas. If intercepted snow
  !! and liquid water are both present on the canopy, part of the liquid water cover is assumed to underlie the
  !! snow cover, so the fractional liquid water coverage is decreased by the fractional snow coverage, to yield
  !! the fractional coverage of liquid water that is exposed to the air.
  !!
  !! Next, tests are performed to ascertain whether the liquid water and snow on the canopy exceed their
  !! respective interception capacities. If so, the excess is assigned to RRESID for liquid water and SRESID
  !! for snow, and the intercepted liquid water or snow is reset to the respective interception capacity. The
  !! sum of RRESID and SRESID is added to WTRC, the diagnosed residual water transferred off the
  !! canopy, and the diagnosed change in internal energy of the canopy, HTCC, is updated. If the fractional
  !! coverage of snow on the modelled area is greater than zero, SRESID is added to the snow pack; the snow
  !! depth, mass, and temperature are recalculated, the diagnosed change in internal energy of the snow pack
  !! HTCS is updated, and SRESID is added to WTRS, the diagnosed residual water transferred to or from
  !! the snow pack. The remaining amounts of RRESID and SRESID are added to the soil. For each layer in
  !! turn whose permeable depth is greater than zero, if the sum of RRESID, SRESID and the ambient soil
  !! liquid and frozen moisture contents is less than the pore volume THPOR, RRESID is added to the liquid
  !! water content and SRESID is added to the frozen moisture content. The layer temperature is
  !! recalculated, the diagnosed change in internal energy HTC is updated, and RRESID and SRESID are
  !! added to WTRG, the diagnosed residual water transferred into or out of the soil, and are then set to zero.
  !!
  do I = IL1,IL2 ! loop 200
    PAICAN(I) = 0.
    if (FC(I) > 0.) then
      do J = 1,IC
        PAICAN(I) = PAICAN(I) + FCAN(I,J) * PAI(I,J)
      end do
      PAICAN(I) = PAICAN(I) / FC(I)
    end if

    PAICNS(I) = 0.0
    if (FCS(I) > 0.) then
      do J = 1,IC
        PAICNS(I) = PAICNS(I) + FCANS(I,J) * PAIS(I,J)
      end do
      PAICNS(I) = PAICNS(I) / FCS(I)
    end if
    !
    CWLCAP(I) = 0.20 * PAICAN(I)  
    CWLCPS(I) = 0.20 * PAICNS(I)  
    !
    RRESID(I) = 0.0
    if (RCAN(I) < 1.0E-5 .or. (FC(I) + FCS(I)) < 1.0E-5) then
      RRESID(I) = RRESID(I) + RCAN(I)
      RCAN(I) = 0.0
    end if
    !
    if (RCAN(I) > 0. .and. (FC(I) + FCS(I)) > 0.) then
      RCAN(I) = RCAN(I) / (FC(I) + FCS(I))
      if (PAICAN(I) > 0.0) then
        RAICAN(I) = RCAN(I) * (FC(I) + FCS(I)) / (FC(I) + FCS(I) * &
                    PAICNS(I) / PAICAN(I))
      else
        RAICAN(I) = 0.0
      end if
      if (PAICNS(I) > 0.0) then
        RAICNS(I) = RCAN(I) * (FC(I) + FCS(I)) / (FCS(I) + FC(I) * &
                    PAICAN(I) / PAICNS(I))
      else
        RAICNS(I) = 0.0
      end if
    else
      RAICAN(I) = 0.0
      RAICNS(I) = 0.0
    end if

    ! For snow interception on the canopy, a modified calculation of the plant area indices \f$\Lambda_{p, 0}\f$ and \f$\Lambda_{p, s}\f$ is performed, assigning a weight of 0.7 to the plant area index of needleleaf trees, to account for the effect
    ! of needle clumping.
    PAICAN(I) = 0.
    if (FC(I) > 0.) then
      do J = 1,IC
        select case (classpfts(J))
        case ('BdlTr', 'Crops', 'Grass', 'BdlSh')
          PAICAN(I) = PAICAN(I) + FCAN(I,J) * PAI(I,J)
        case ('NdlTr')
          PAICAN(I) = PAICAN(I) + 0.7 * FCAN(I,J) * PAI(I,J)
        case default
          ! Error checking, if no known CLASS PFT given, bail.
          print * ,'calcLandSurfParams says: Unknown CLASS pft => ' &
                    ,classpfts(j)
          call errorHandler('calcLandSurfParams', - 10)
        end select
      end do
      PAICAN(I) = PAICAN(I) / FC(I)
    end if

    PAICNS(I) = 0.0
    if (FCS(I) > 0.) then
      do J = 1,IC
        select case (classpfts(J))
        case ('BdlTr', 'Crops', 'Grass', 'BdlSh')
          PAICNS(I) = PAICNS(I) + FCANS(I,J) * PAIS(I,J)
        case ('NdlTr')
          PAICNS(I) = PAICNS(I) + 0.7 * FCANS(I,J) * PAIS(I,J)
        case default
          ! Error checking, if no known CLASS PFT given, bail.
          print * ,'calcLandSurfParams says: Unknown CLASS pft => ' &
                    ,classpfts(j)
          call errorHandler('calcLandSurfParams', - 11)
        end select
      end do
      PAICNS(I) = PAICNS(I) / FCS(I)
    end if
    !
    ! Calculate the canopy storage capacity for precip:

    CWFCAP(I) = 6.0 * PAICAN(I) * (0.27 + 46.0 / RHOSNI(I))   
    CWFCPS(I) = 6.0 * PAICNS(I) * (0.27 + 46.0 / RHOSNI(I))   
    !
    SRESID(I) = 0.0
    if (SNCAN(I) < 1.0E-5 .or. (FC(I) + FCS(I)) < 1.0E-5) then
      SRESID(I) = SRESID(I) + SNCAN(I)
      SNCAN(I) = 0.0
    end if
    !
    if (SNCAN(I) > 0. .and. (FC(I) + FCS(I)) > 0.) then
      SNCAN(I) = SNCAN(I) / (FC(I) + FCS(I))
      if (PAICAN(I) > 0.0) then
        SNOCAN(I) = SNCAN(I) * (FC(I) + FCS(I)) / (FC(I) + FCS(I) * &
                    PAICNS(I) / PAICAN(I))
      else
        SNOCAN(I) = 0.0
      end if
      if (PAICNS(I) > 0.0) then
        SNOCNS(I) = SNCAN(I) * (FC(I) + FCS(I)) / (FCS(I) + FC(I) * &
                    PAICAN(I) / PAICNS(I))
      else
        SNOCNS(I) = 0.0
      end if
    else
      SNOCAN(I) = 0.0
      SNOCNS(I) = 0.0
    end if
    !
    if (CWFCAP(I) > 0.0) then
      FSNOWC(I) = MIN(SNOCAN(I) / CWFCAP(I),1.0)
    else
      FSNOWC(I) = 0.0
    end if
    if (CWFCPS(I) > 0.0) then
      FSNOCS(I) = MIN(SNOCNS(I) / CWFCPS(I),1.0)
    else
      FSNOCS(I) = 0.0
    end if
    !
    if (CWLCAP(I) > 0.0) then
      FRAINC(I) = MIN(RAICAN(I) / CWLCAP(I),1.0)
    else
      FRAINC(I) = 0.0
    end if
    if (CWLCPS(I) > 0.0) then
      FRAICS(I) = MIN(RAICNS(I) / CWLCPS(I),1.0)
    else
      FRAICS(I) = 0.0
    end if
    FRAINC(I) = MAX(0.0,MIN(FRAINC(I) - FSNOWC(I),1.0))
    FRAICS(I) = MAX(0.0,MIN(FRAICS(I) - FSNOCS(I),1.0))
    !
    if (RAICAN(I) > CWLCAP(I)) then
      RRESID(I) = RRESID(I) + FC(I) * (RAICAN(I) - CWLCAP(I))
      RAICAN(I) = CWLCAP(I)
    end if
    if (SNOCAN(I) > CWFCAP(I)) then
      SRESID(I) = SRESID(I) + FC(I) * (SNOCAN(I) - CWFCAP(I))
      SNOCAN(I) = CWFCAP(I)
    end if
    !
    if (RAICNS(I) > CWLCPS(I)) then
      RRESID(I) = RRESID(I) + FCS(I) * (RAICNS(I) - CWLCPS(I))
      RAICNS(I) = CWLCPS(I)
    end if
    if (SNOCNS(I) > CWFCPS(I)) then
      SRESID(I) = SRESID(I) + FCS(I) * (SNOCNS(I) - CWFCPS(I))
      SNOCNS(I) = CWFCPS(I)
    end if
    !
    WTRC (I) = WTRC(I) - (RRESID(I) + SRESID(I)) / DELT
    HTCC (I) = HTCC(I) - TCAN(I) * (SPHW * RRESID(I) + SPHICE * SRESID(I)) / &
               DELT
    if (FSNOW(I) > 0.0) then
      SNOI = SNO(I)
      ZSNADD = SRESID(I) / (RHOSNO(I) * FSNOW(I))
      ZSNOW(I) = ZSNOW(I) + ZSNADD
      SNO(I) = ZSNOW(I) * FSNOW(I) * RHOSNO(I)
      TSNOW(I) = (TCAN(I) * SPHICE * SRESID(I) + TSNOW(I) * HCPICE * &
                 SNOI / RHOICE) / (HCPICE * SNO(I) / RHOICE)
      HTCS (I) = HTCS(I) + TCAN(I) * SPHICE * SRESID(I) / DELT
      WTRS (I) = WTRS(I) + SRESID(I) / DELT
      SRESID(I) = 0.0
    end if
    !
    do J = 1,IG ! loop 190
      if (DELZW(I,J) > 0.0 .and. (RRESID(I) > 0.0 &
          .or. SRESID(I) > 0.0)) then
        THSUM = THLIQ(I,J) + THICE(I,J) + &
                (RRESID(I) + SRESID(I)) / (RHOW * DELZW(I,J))
        if (THSUM < THPOR(I,J)) then
          THICEI = THICE(I,J)
          THLIQI = THLIQ(I,J)
          THICE(I,J) = THICE(I,J) + SRESID(I) / &
                       (RHOICE * DELZW(I,J))
          THLIQ(I,J) = THLIQ(I,J) + RRESID(I) / &
                       (RHOW * DELZW(I,J))
          TBAR(I,J) = (TBAR(I,J) * ((DELZ(J) - DELZW(I,J)) * &
                      HCPSND + DELZW(I,J) * (THLIQI * HCPW + THICEI * &
                      HCPICE + (1.0 - THPOR(I,J)) * HCPS(I,J))) + TCAN(I) * &
                      (RRESID(I) * HCPW / RHOW + SRESID(I) * HCPICE / RHOICE)) &
                      / ((DELZ(J) - DELZW(I,J)) * HCPSND + DELZW(I,J) * &
                      (HCPW * THLIQ(I,J) + HCPICE * THICE(I,J) + HCPS(I,J) * &
                      (1.0 - THPOR(I,J))))
          HTC(I,J) = HTC(I,J) + TCAN(I) * (RRESID(I) * HCPW / RHOW + &
                     SRESID(I) * HCPICE / RHOICE) / DELT
          WTRG (I) = WTRG(I) + (RRESID(I) + SRESID(I)) / DELT
          RRESID(I) = 0.0
          SRESID(I) = 0.0
        end if
      end if
    end do ! loop 190
    !
  end do ! loop 200
  !
  !     * CALCULATION OF ROUGHNESS LENGTHS FOR HEAT AND MOMENTUM AND
  !     * ZERO-PLANE DISPLACEMENT FOR CANOPY OVERLYING BARE SOIL AND
  !     * CANOPY OVERLYING SNOW.
  !

  !>
  !! In loops 250 and 275, calculations of the displacement height and the logarithms of the roughness lengths
  !! for heat and momentum are performed for the vegetated subareas. The displacement height \f$d_i\f$ and the
  !! roughness length \f$z_{0, i}\f$ for the separate vegetation categories are obtained as simple ratios of the canopy
  !! height H:
  !! \f$d_i = 0.70 H\f$
  !! \f$z_{0, i} = 0.10 H\f$
  !!
  !! The averaged displacement height d over the vegetated subareas is only calculated if the flag IDISP has
  !! been set to 1. If IDISP = 0, this indicates that the atmospheric model is using a terrain-following
  !! coordinate system, and thus the displacement height is treated as part of the "terrain". If DISP = 1, d is
  !! calculated as a logarithmic average over the vegetation categories:
  !! \f$X ln(d) = \Sigma [X_i ln(d_i)]\f$
  !! where X is the fractional coverage of the subarea. The averaged roughness length for momentum \f$z_{0m}\f$
  !! over the subarea is determined based on the assumption that averaging should be performed on the basis
  !! of the drag coefficient formulation. Thus, following Delage et al. (1999) \cite Delage1999-vj, and after Mason (1988):
  !! \f$X/ln^2 (z_b /z_{0m}) = \Sigma [X_i /ln^2 (z_b /z_{0i})]\f$
  !!
  !! The averaged roughness length for heat \f$z_{0e}\f$ over the subarea is calculated as a geometric mean over the
  !! vegetation categories:
  !! \f$z_{0e} z_{0mX} = \Pi ( z_{0i}^{2Xi} )\f$
  !!

  do J = 1,IC ! loop 250
    do I = IL1,IL2
      if (FC(I) > 0. .and. H(I,J) > 0.) then
        if (IDISP == 1)   DISP(I) = DISP(I) + FCAN (I,J) * &
            LOG(0.7 * H(I,J))   
        ZOMLNC(I) = ZOMLNC(I) + FCAN (I,J) / &
                    ((LOG(ZBLEND(I) / (0.1 * H(I,J)))) ** 2)   
        ZOELNC(I) = ZOELNC(I) * &
                    (0.01 * H(I,J) * H(I,J) / ZORAT(IC)) ** FCAN(I,J)
      end if
      if (FCS(I) > 0. .and. HS(I,J) > 0.) then
        if (IDISP == 1)   DISPS(I) = DISPS (I) + FCANS(I,J) * &
            LOG(0.7 * HS(I,J))
        ZOMLCS(I) = ZOMLCS(I) + FCANS(I,J) / &
                    ((LOG(ZBLEND(I) / (0.1 * HS(I,J)))) ** 2)
        ZOELCS(I) = ZOELCS(I) * &
                    (0.01 * HS(I,J) * HS(I,J) / ZORAT(IC)) ** FCANS(I,J)
      end if
    end do
  end do ! loop 250
  !
  do I = IL1,IL2 ! loop 275
    if (FC(I) > 0.) then
      if (IDISP == 1)   DISP(I) = EXP(DISP(I) / FC(I))
      ZOMLNC(I) = ZBLEND(I) / EXP(SQRT(1.0 / (ZOMLNC(I) / FC(I))))
      ZOELNC(I) = LOG(ZOELNC(I) ** (1.0 / FC(I)) / ZOMLNC(I))
      ZOMLNC(I) = LOG(ZOMLNC(I))
    end if
    if (FCS(I) > 0.) then
      if (IDISP == 1)   DISPS(I) = EXP(DISPS(I) / FCS(I))
      ZOMLCS(I) = ZBLEND(I) / EXP(SQRT(1.0 / (ZOMLCS(I) / FCS(I))))
      ZOELCS(I) = LOG(ZOELCS(I) ** (1.0 / FCS(I)) / ZOMLCS(I))
      ZOMLCS(I) = LOG(ZOMLCS(I))
    end if
  end do ! loop 275
  !
  !     * ADJUST ROUGHNESS LENGTHS OF BARE SOIL AND SNOW-COVERED BARE
  !     * SOIL FOR URBAN ROUGHNESS IF PRESENT.
  !
  !>
  !! In loop 300, calculations of the logarithms of the roughness lengths for heat and momentum are
  !! performed for the bare ground and snow-covered subareas. Background values of \f$ln(z_{om})\f$ for soil, snow
  !! cover, ice sheets and urban areas are passed into the subroutine through common blocks. In CLASS,
  !! urban areas are treated very simply, as areas of bare soil with a high roughness length. The subarea values
  !! of \f$ln(z_{om})\f$ for bare soil and snow are therefore adjusted for the fractional coverage of urban area. Values
  !! for the ratio between the roughness lengths for momentum and heat for bare soil and snow are also
  !! passed in via common blocks. These are used to derive subarea values of \f$ln(z_{oe})\f$ from \f$ln(z_{om})\f$.
  !!
  !! In the same loop the roughness length for peatlands is also calculated assuming
  !! a natural log of the roughness length of the moss surface is -6.57 (parameter stored in classicParams.f90)
  !!
  do I = IL1,IL2 ! loop 300
    if (FG(I) > 0.) then
      if (ISAND(I,1) /= - 4) then
        ZOMLNG(I) = ((FG(I) - FCANMX(I,ICP1) * (1.0 - FSNOW(I))) * ZOLNG + &
                    FCANMX(I,ICP1) * (1.0 - FSNOW(I)) * ZOLN(I,ICP1)) &
                    / FG(I)
        if (ipeatland(i) > 0) then ! roughness length of moss surface in peatlands.
          ZOMLNG(I) = ((FG(I) - FCANMX(I,ICP1) * (1.0 - FSNOW(I))) &
                      * zolnmoss + FCANMX(I,ICP1) * (1.0 - FSNOW(I)) * &
                      ZOLN(I,ICP1)) / FG(I)
        end if
      else
        ZOMLNG(I) = ZOLNI
      end if
      ZOELNG(I) = ZOMLNG(I) - LOG(ZORATG)
    end if
    if (FGS(I) > 0.) then
      ZOMLNS(I) = ((FGS(I) - FCANMX(I,ICP1) * FSNOW(I)) * ZOLNS + &
                  FCANMX(I,ICP1) * FSNOW(I) * ZOLN(I,ICP1)) / FGS(I)
      ZOELNS(I) = ZOMLNS(I) - LOG(ZORATG)
    end if
  end do ! loop 300
  !
  !     * ADD CONTRIBUTION OF OROGRAPHY TO MOMENTUM ROUGHNESS LENGTH
  !

  !>
  !! In loop 325, an adjustment is applied to \f$ln(z_{om})\f$ if the effect of terrain roughness needs to be taken into
  !! account. If the surface orographic roughness length is not vanishingly small, its logarithm, LZ0ORO, is
  !! calculated. If it is greater than the calculated logarithm of the roughness length for momentum of any of
  !! the subareas, these are reset to LZ0ORO.
  !!
  do I = IL1,IL2 ! loop 325
    if (Z0ORO(I) > 1.0E-4) then
      LZ0ORO = LOG(Z0ORO(I))
    else
      LZ0ORO = - 10.0
    end if
    ZOMLNC(I) = MAX(ZOMLNC(I),LZ0ORO)
    ZOMLCS(I) = MAX(ZOMLCS(I),LZ0ORO)
    ZOMLNG(I) = MAX(ZOMLNG(I),LZ0ORO)
    ZOMLNS(I) = MAX(ZOMLNS(I),LZ0ORO)
  end do ! loop 325
  !
  !     * CALCULATE HEAT CAPACITY FOR CANOPY OVERLYING BARE SOIL AND
  !     * CANOPY OVERLYING SNOW.
  !     * ALSO CALCULATE INSTANTANEOUS GRID-CELL AVERAGED CANOPY MASS.
  !

  !>
  !! In loop 350, the canopy mass is calculated as a weighted average over the vegetation categories, for
  !! canopy over bare soil (CMASSC) and over snow (CMASCS). (For crops over bare soil, the mass is
  !! adjusted according to the growth index; for crops and grass over snow, the mass is additionally adjusted
  !! to take into account burying by snow.) If IDISP = 0, indicating that the vegetation displacement height is
  !! part of the "terrain", the mass of air within the displacement height is normalized by the vegetation heat
  !! capacity and added to the canopy mass. If IZREF = 2, indicating that the bottom of the atmosphere is
  !! taken to lie at the local roughness length rather than at the ground surface, the mass of air within the
  !! roughness length is likewise normalized by the vegetation heat capacity and added to the canopy mass.
  !! The canopy heat capacities over bare soil (CHCAP) and over snow (CHCAPS) are evaluated from the
  !! respective values of canopy mass and of intercepted liquid water and snow. The aggregated canopy mass
  !! CMAI is recalculated, and is used to determine the change in internal energy of the canopy, HTCC, owing
  !! to growth or disappearance of the vegetation.
  !!
  do I = IL1,IL2 ! loop 350
    if (FC(I) > 0.) then ! For non-snow covered ground
      if (ctem_on) then
        CMASSC(I) = 0.0
        do J = 1,IC
          CMASSC(I) = CMASSC(I) + FCAN(I,J) * CMASVEGC(I,J)
        end do
        CMASSC(I) = CMASSC(I) / FC(I)
      else ! CLASS (physics) only
        CMASSC(I) = 0.0
        do J = 1,IC
          select case (classpfts(J))
          case ('NdlTr','BdlTr','Grass', 'BdlSh')
            CMASSC(I) = CMASSC(I) + FCAN(I,J) * CWGTMX(I,J)
          case ('Crops')
            CMASSC(I) = CMASSC(I) + FCAN(I,J) * CWGTMX(I,J) * GROWA(I)
          case default
            ! Error checking,if no known CLASS PFT given,bail.
            print * ,'calcLandSurfParams says: Unknown CLASS pft => ' &
                     ,classpfts(j)
            call errorHandler('calcLandSurfParams', - 12)
          end select
        end do
        CMASSC(I) = CMASSC(I) / FC(I)
      end if ! CTEM on/off

      !> if idisp=0, vegetation displacement heights are ignored, because the atmospheric model considers these to be part of the "terrain". if idisp=1, vegetation displacement heights are calculated.
      if (IDISP == 0) then
        temp1 = 0.0
        do J = 1,IC
          temp1 = temp1 + FCAN(I,J) * H(I,J)
        end do
        CMASSC(I) = CMASSC(I) + RHOAIR(I) * (SPHAIR / SPHVEG) * 0.7 * & 
                    temp1 / FC(I)
      end if
      !> if izref=1, the bottom of the atmospheric model is taken lie at the ground surface.
      !> if izref=2, the bottom of the atmospheric model is taken to lie at the local roughness height.
      if (IZREF == 2) then
        temp1 = 0.0
        do J = 1,IC
          temp1 = temp1 + FCAN(I,J) * H(I,J)
        end do
        CMASSC(I) = CMASSC(I) + RHOAIR(I) * (SPHAIR / SPHVEG) * 0.1 * &
                    temp1 / FC(I)
      end if
    end if

    if (FCS(I) > 0.) then ! For snow covered ground

      if (ctem_on) then
        CMASCS(I) = 0.0
        do J = 1,IC
          select case (classpfts(J))
          case ('NdlTr','BdlTr','BdlSh')
            CMASCS(I) = CMASCS(I) + FCANS(I,J) * CMASVEGC(I,J)
          case ('Crops','Grass')
            CMASCS(I) = CMASCS(I) + FCANS(I,J) * CMASVEGC(I,J) * &
                        HS(I,J) / MAX(H(I,J),HS(I,J))
          case default
            ! Error checking, if no known CLASS PFT given, bail.
            print * ,'calcLandSurfParams says: Unknown CLASS pft => ', &
                       classpfts(j)
            call errorHandler('calcLandSurfParams', - 13)
          end select
        end do
        CMASCS(I) = CMASCS(I) / FCS(I)
      else ! CLASS (physics) only
        CMASCS(I) = 0.0
        do J = 1,IC
          select case (classpfts(J))
          case ('NdlTr','BdlTr','BdlSh')
            CMASCS(I) = CMASCS(I) + FCANS(I,J) * CWGTMX(I,J)
          case ('Crops')
            CMASCS(I) = CMASCS(I) + FCANS(I,J) * CWGTMX(I,J) * &
                        GROWA(I) * HS(I,J) / MAX(H(I,J),HS(I,J))
          case ('Grass')
            CMASCS(I) = CMASCS(I) + FCANS(I,J) * CWGTMX(I,J) * &
                        HS(I,J) / MAX(H(I,J),HS(I,J))
          case default
            ! Error checking, if no known CLASS PFT given, bail.
            print * ,'calcLandSurfParams says: Unknown CLASS pft => ', &
                      classpfts(j)
            call errorHandler('calcLandSurfParams', - 14)
          end select
        end do
        CMASCS(I) = CMASCS(I) / FCS(I)
      end if   ! CTEM on/off
      !
      !> if idisp=0, vegetation displacement heights are ignored, because the atmospheric model considers these to be part of the "terrain". if idisp=1, vegetation displacement heights are calculated.
      if (IDISP == 0) then
        temp1 = 0.0
        do J = 1,IC
          temp1 = temp1 + FCANS(I,J) * HS(I,J)
        end do
        CMASCS(I) = CMASCS(I) + RHOAIR(I) * (SPHAIR / SPHVEG) * 0.7 * &
                    temp1 / FCS(I)
      end if

      !> if izref=1, the bottom of the atmospheric model is taken lie at the ground surface.
      !> if izref=2, the bottom of the atmospheric model is taken to lie at the local roughness height.
      if (IZREF == 2) then
        temp1 = 0.0
        do J = 1,IC
          temp1 = temp1 + FCANS(I,J) * HS(I,J)
        end do
        CMASCS(I) = CMASCS(I) + RHOAIR(I) * (SPHAIR / SPHVEG) * 0.1 * &
                    temp1 / FCS(I)
      end if
    end if

    CHCAP (I) = SPHVEG * CMASSC(I) + SPHW * RAICAN(I) + SPHICE * SNOCAN(I)
    CHCAPS(I) = SPHVEG * CMASCS(I) + SPHW * RAICNS(I) + SPHICE * SNOCNS(I)
    HTCC  (I) = HTCC(I) - SPHVEG * CMAI(I) * TCAN(I) / DELT

    if (ctem_on) then
      CMAI  (I) = FC(I) * CMASSC(I) + FCS(I) * CMASCS(I)
      if (CMAI(I) < 1.0E-5 .and. (CMASSC(I) > 0.0 .or. &
          CMASCS(I) > 0.0)) TCAN(I) = TA(I)
    else ! CLASS (physics) only
      if (CMAI(I) < 1.0E-5 .and. (CMASSC(I) > 0.0 .or. &
          CMASCS(I) > 0.0)) TCAN(I) = TA(I)
      CMAI  (I) = FC(I) * CMASSC(I) + FCS(I) * CMASCS(I)
    end if

    HTCC  (I) = HTCC(I) + SPHVEG * CMAI(I) * TCAN(I) / DELT
    RBCOEF(I) = 0.0
  end do ! loop 350
  !
  !     * CALCULATE VEGETATION ROOTING DEPTH AND FRACTION OF ROOTS
  !     * IN EACH SOIL LAYER (SAME FOR SNOW/BARE SOIL CASES).
  !     * ALSO CALCULATE LEAF BOUNDARY RESISTANCE PARAMETER RBCOEF.
  !

  !>
  !! In the 450 and 500 loops, the fraction of plant roots in each soil layer is calculated. If CLASS is being run
  !! coupled to CTEM, the CTEM-derived values are assigned. Otherwise, for each vegetation category the
  !! rooting depth ZROOT is set to the background maximum value, except in the case of crops, for which it
  !! is set to the maximum scaled by GROWA. If the soil permeable depth is less than ZROOT, ZROOT is
  !! set to this depth instead. Values are then assigned in the matrix RMAT, which stores the fraction of roots
  !! in each vegetation category for each soil layer. According to Feddes et al. (1974) \cite Feddes1974-ff, the fractional root
  !! volume R(z) below a depth z is well represented for many varieties of plants by the following exponential
  !! function:
  !! \f$R(z) = a_1 exp(-3.0z) + a_2.\f$
  !!
  !! Making use of the boundary conditions \f$R(0) = 1\f$ and \f$R(z_r) = 0\f$, where \f$z_r\f$ is the rooting depth ZROOT, it
  !! can be seen that the fraction of roots within a soil depth interval \f$\Delta z\f$ can be obtained as the difference
  !! between R(z) evaluated at the top \f$(z_T)\f$ and bottom \f$(z_B)\f$ of the interval:
  !! \f$R(\Delta z) = [exp(-3.0z_T) - exp(-3.0z_B)]/ [1 - exp(-3.0z_r)]\f$
  !!
  !! The total fraction of roots in each soil layer, FROOT for snow-free areas and FROOTS for snow-covered areas,
  !! can then be determined as weighted averages over the four vegetation categories.
  !!
  !! In loop 450, a leaf boundary resistance parameter \f$C_{rb}\f$ , incorporating the plant area indices of the four
  !! vegetation subareas, is also calculated for later use in subroutine energBalVegSolve:
  !! \f$C_{rb} = C_l \Lambda_{p, i}^{0.5} /0.75 \bullet [1 - exp(-0.75 \Lambda_{p, i}^{0.5} )]\f$
  !! where \f$C_l\f$ is a parameter that varies with the vegetation category. The aggregated value of \f$C_{rb}\f$ is obtained
  !! as a weighted average over the four vegetation categories over bare ground and snow cover.
  !!
  do J = 1,IC ! loop 450
    do I = IL1,IL2
      if (ctem_on) then
        do K = 1,IG
          RMAT(I,J,K) = RMATC(I,J,K)
        end do
      else
        ZROOT = ZRTMAX(I,J)
        if (J == 3) ZROOT = ZRTMAX(I,J) * GROWA(I)
        ZROOTG = 0.0
        do K = 1,IG ! loop 375
          ZROOTG = ZROOTG + DELZW(I,K)
        end do ! loop 375
        ZROOT = MIN(ZROOT,ZROOTG)
        do K = 1,IG ! loop 400
          if (ZROOT <= (ZBOTW(I,K) - DELZW(I,K) + 0.0001)) then
            RMAT(I,J,K) = 0.0
          else if (ZROOT <= ZBOTW(I,K)) then
            RMAT(I,J,K) = (EXP( - 3.0 * (ZBOTW(I,K) - DELZW(I,K))) - & 
                          EXP( - 3.0 * ZROOT)) / (1.0 - EXP( - 3.0 * ZROOT))
          else
            RMAT(I,J,K) = (EXP( - 3.0 * (ZBOTW(I,K) - DELZW(I,K))) - &
                          EXP( - 3.0 * ZBOTW(I,K))) / (1.0 - EXP( - 3.0 * ZROOT))
          end if
        end do ! loop 400
      end if
      if ((FC(I) + FCS(I)) > 0.) then
        RBCOEF(I) = RBCOEF(I) + &
                    (FCAN(I,J) * XLEAF(J) * (SQRT(PAI(I,J)) / 0.75) * &
                    (1.0 - EXP( - 0.75 * SQRT(PAI(I,J)))) + &
                    FCANS(I,J) * XLEAF(J) * (SQRT(PAIS(I,J)) / 0.75) * &
                    (1.0 - EXP( - 0.75 * SQRT(PAIS(I,J))))) / &
                    (FC(I) + FCS(I))
      end if
    end do
  end do ! loop 450
  !
  do K = 1,IG ! loop 500
    do I = IL1,IL2
      FROOT(I,K) = 0.0
      FROOTS(I,K) = 0.0
      if (FC(I) > 0.) then
        do J = 1,IC
          FROOT(I,K) = FROOT(I,K) + FCAN(I,J) * RMAT(I,J,K)
        end do
        FROOT(I,K) = FROOT(I,K) / FC(I)
      end if
      if (FCS(I) > 0.) then
        do J = 1,IC
          FROOTS(I,K) = FROOTS(I,K) + FCANS(I,J) * RMAT(I,J,K)
        end do
        FROOTS(I,K) = FROOTS(I,K) / FCS(I)
      end if
    end do
  end do ! loop 500
  !
  !     * CALCULATE SKY-VIEW FACTORS FOR BARE GROUND AND SNOW
  !     * UNDERLYING CANOPY.
  !

  !>
  !! In loop 600, the sky view factor \f$\chi\f$ of the ground underlying the canopy is calculated for the vegetated
  !! subareas. The standard approach is to determine \f$\chi\f$ as an exponential function of the plant area index \f$\Lambda_p\f$ :
  !! \f$\chi = exp[-c \Lambda_p ]\f$
  !! where c is a constant depending on the vegetation category. The subarea values of \f$\chi\f$ are obtained as
  !! weighted averages over the four vegetation categories.
  !!

  do I = IL1,IL2 ! loop 600
    FSVF (I) = 0.0
    if (FC(I) > 0.) then
      do J = 1,IC
        FSVF(I) = FSVF(I) + FCAN(I,J) * EXP(CANEXT(J) * PAI(I,J))
      end do
      FSVF(I) = FSVF(I) / FC(I)
    end if
    FSVFS(I) = 0.0
    if (FCS(I) > 0.) then
      do J = 1,IC
        FSVFS(I) = FSVFS(I) + FCANS(I,J) * EXP(CANEXT(J) * PAIS(I,J))
      end do
      FSVFS(I) = FSVFS(I) / FCS(I)
    end if
  end do ! loop 600
  !
  !     * CALCULATE BULK SOIL MOISTURE SUCTION FOR STOMATAL RESISTANCE.
  !     * CALCULATE FRACTIONAL TRANSPIRATION EXTRACTED FROM SOIL LAYERS.
  !

  !>
  !! In the 650 loop, the fraction of the total transpiration of water by plants that is extracted from each soil
  !! layer is determined. This is done by weighting the values of FROOT and FROOTS calculated above by the relative soil
  !! moisture suction in each layer: \f$( \Psi_w - \Psi_i)/( \Psi_w - \Psi_{sat} )\f$
  !! where \f$\Psi_i\f$ , the soil moisture suction in the layer, is obtained as
  !! \f$\Psi_i = \Psi_{sat} ( \theta_{l, i} / \theta_p)^{-b}\f$
  !! In these equations \f$\Psi_w\f$ is the soil moisture suction at the wilting point, \f$\Psi_{sat}\f$ is the suction at
  !! saturation, \f$\theta_{l, i}\f$ is the volumetric liquid water content of the soil layer, \f$\theta_p\f$ is the pore
  !! volume, and b is an empirical parameter developed by Clapp and Hornberger (1978) \cite Clapp1978-898. The layer values of FROOT and FROOTS are then
  !! re-normalized so that their sum adds up to unity. In this loop, the representative soil moisture suction PSIGND is also
  !! calculated for later use in the vegetation stomatal resistance formulation, as the minimum value of \f$\Psi_i\f$ and
  !! \f$\Psi_w\f$ over all the soil layers.
  !!
  do J = 1,IG ! loop 650
    do I = IL1,IL2
      if (FCS(I) > 0.0 .or. FC(I) > 0.0) then
        if (THLIQ(I,J) > (THLMIN(I,J) + 0.01)) then
          PSII = PSISAT(I,J) * (THLIQ(I,J) / THPOR(I,J)) ** ( - BI(I,J))
          PSII = MIN(PSII,PSIWLT(I,J))
          if (FROOT(I,J) > 0.0) PSIGND(I) = MIN(PSIGND(I),PSII)
          PSIRAT = (PSIWLT(I,J) - PSII) / (PSIWLT(I,J) - PSISAT(I,J))
          FROOT(I,J) = FROOT(I,J) * PSIRAT
          FROOTS(I,J) = FROOTS(I,J) * PSIRAT
          FRTOT(I) = FRTOT(I) + FROOT(I,J)
          FRTOTS(I) = FRTOTS(I) + FROOTS(I,J)
        else
          FROOT(I,J) = 0.0
          FROOTS(I,J) = 0.0
        end if
      end if
    end do
  end do ! loop 650

  do J = 1,IG ! loop 700
    do I = IL1,IL2
      if (FRTOT(I) > 0.) then
        FROOT(I,J) = FROOT(I,J) / FRTOT(I)
      end if
      if (FRTOTS(I) > 0.) then
        FROOTS(I,J) = FROOTS(I,J) / FRTOTS(I)
      end if
    end do
  end do ! loop 700
  !
  !     * CALCULATE EFFECTIVE LEAF AREA INDICES FOR TRANSPIRATION.
  !

  !>
  !! Finally, in loop 800 the aggregated canopy plant area indices PAICAN and PAICNS are set back to their
  !! original values, from the modified values used for the snow interception calculations above; and if CLASS
  !! is being run coupled with CTEM, a set of CTEM-related calculations is performed.
  !!

  do I = IL1,IL2 ! loop 800
    PAICAN(I) = 0.0
    if (FC(I) > 0.) then
      do J = 1,IC
        PAICAN(I) = PAICAN(I) + FCAN(I,J) * PAI(I,J)
      end do
      PAICAN(I) = PAICAN(I) / FC(I)
    end if

    PAICNS(I) = 0.0
    if (FCS(I) > 0.) then
      do J = 1,IC
        PAICNS(I) = PAICNS(I) + FCANS(I,J) * PAIS(I,J)
      end do
      PAICNS(I) = PAICNS(I) / FCS(I)
    end if
  end do ! loop 800
  !
  if (ctem_on) then
    !
    !       * ESTIMATE FCANC AND FCANCS FOR USE BY PHTSYN SUBROUTINE BASED ON
    !       * FCAN AND FCANS FOR CTEM PFTS.
    !
    SFCANCMX(:,:) = 0.0  ! SUM OF FCANCMXS
    K1 = 0
    do J = 1,IC ! loop 830
      if (J == 1) then
        K1 = K1 + 1
      else
        K1 = K1 + NOL2PFTS(J - 1)
      end if
      K2 = K1 + NOL2PFTS(J) - 1
      do M = K1,K2 ! loop 820
        do I = IL1,IL2
          SFCANCMX(I,J) = SFCANCMX(I,J) + FCANCMX(I,M)
        end do
      end do ! loop 820
    end do ! loop 830
    K1 = 0
    do J = 1,IC ! loop 860
      if (J == 1) then
        K1 = K1 + 1
      else
        K1 = K1 + NOL2PFTS(J - 1)
      end if
      K2 = K1 + NOL2PFTS(J) - 1
      do M = K1,K2 ! loop 850
        do I = IL1,IL2
          if (SFCANCMX(I,J) > 1.E-20) then
            FCANC(I,M)  = FCAN(I,J) * (FCANCMX(I,M) / SFCANCMX(I,J))
            FCANCS(I,M) = FCANS(I,J) * (FCANCMX(I,M) / SFCANCMX(I,J))
          else
            FCANC(I,M)  = 0.0
            FCANCS(I,M) = 0.0
          end if
        end do
      end do ! loop 850
    end do ! loop 860
  end if

  return
end subroutine calcLandSurfParams
