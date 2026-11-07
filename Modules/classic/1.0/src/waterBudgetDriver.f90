!> \file
!! Calls subroutines to perform surface water budget calculations
!! @author D. Verseghy, P. Bartlett, M. Lazare
!
subroutine waterBudgetDriver (THLIQ, THICE, TBAR, TCAN, RCAN, SNCAN, & ! Formerly CLASSW
                              RUNOFF, TRUNOF, SNO, TSNOW, RHOSNO, ALBSNO, &
                              WSNOW, ZPOND, TPOND, GROWTH, TBASE, GFLUX, &
                              PCFC, PCLC, PCPN, PCPG, QFCF, QFCL, &
                              QFN, QFG, QFC, HMFC, HMFG, HMFN, &
                              HTCC, HTCS, HTC, ROFC, ROFN, ROVG, &
                              WTRS, WTRG, OVRFLW, SUBFLW, BASFLW, &
                              TOVRFL, TSUBFL, TBASFL, EVAP, QFLUX, RHOAIR, &
                              TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, &
                              THICEC, THICEG, HCPC, HCPG, RPCP, TRPCP, &
                              SPCP, TSPCP, PCPR, TA, RHOSNI, GGEO, &
                              FC, FG, FCS, FGS, TPONDC, TPONDG, &
                              TPNDCS, TPNDGS, EVAPC, EVAPCG, EVAPG, EVAPCS, &
                              EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG, &
                              RAICAN, SNOCAN, RAICNS, SNOCNS, FSVF, FSVFS, &
                              CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO, &
                              TCANS, CHCAP, CHCAPS, CMASSC, CMASCS, ZSNOW, &
                              GZEROC, GZEROG, GZROCS, GZROGS, G12C, G12G, &
                              G12CS, G12GS, G23C, G23G, G23CS, G23GS, &
                              TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                              ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TSFSAV, &
                              TCTOPC, TCBOTC, TCTOPG, TCBOTG, FROOT, FROOTS, &
                              THPOR, THLRET, THLMIN, BI, PSISAT, GRKSAT, &
                              THLRAT, THFC, XDRAIN, HCPS, DELZ, &
                              DELZW, ZBOTW, XSLOPE, XDRAINH, WFSURF, KSAT, &
                              MANNING_N, DD, &
                              ISAND, IGDR, &
                              IWF, ILG, IL1, IL2, N, &
                              JL, IC, IG, IGP1, IGP2, &
                              NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI)
  !
  !     * AUG 04/15 - M.LAZARE.   SPLIT FROOT INTO TWO ARRAYS, FOR CANOPY
  !     *                         AREAS WITH AND WITHOUT SNOW.
  !     * OCT 03/14 - D.VERSEGHY. CHANGE LIMITING VALUE OF SNOW PACK
  !     *                         FROM 100 KG/M2 TO 10 M.
  !     * AUG 19/13 - M.LAZARE.   ADD CALCULATION OF "QFLUX" (NOW PASSED
  !     *                         IN ALONG WITH "RHOAIR") PREVIOUSLY DONE
  !     *                         IN energyBudgetDriver.
  !     * JUN 21/13 - M.LAZARE.   SET ZSNOW=0. IF THERE IS NO
  !     *                         SNOW IN ANY OF THE 4 SUBCLASSES,
  !     *                         SIMILAR TO WHAT IS DONE FOR THE
  !     *                         OTHER SNOW-RELATED FIELDS.
  !     * OCT 18/11 - M.LAZARE.   PASS IN IGDR THROUGH CALLS TO
  !     *                         waterFlowNonInfiltrate/waterFlowInfiltrate (ORIGINATES NOW
  !     *                         IN soilProperties - ONE CONSISTENT
  !     *                         CALCULATION).
  !     * APR 04/11 - D.VERSEGHY. ADD DELZ TO waterFlowInfiltrate CALL.
  !     * DEC 07/09 - D.VERSEGHY. ADD RADD AND SADD TO waterCalcPrep CALL.
  !     * JAN 06/09 - D.VERSEGHY. INCREASE LIMITING SNOW AMOUNT.
  !     * FEB 25/08 - D.VERSEGHY. MODIFICATIONS REFLECTING CHANGES
  !     *                         ELSEWHERE IN CODE.
  !     * MAR 23/06 - D.VERSEGHY. CHANGES TO ADD MODELLING OF WSNOW
  !     *                         PASS IN GEOTHERMAL HEAT FLUX.
  !     * MAR 21/06 - P.BARTLETT. PASS ADDITIONAL VARIABLES TO waterCalcPrep.
  !     * DEC 07/05 - D.VERSEGHY. REVISIONS TO CALCULATION OF TBASE.
  !     * OCT 05/05 - D.VERSEGHY. MODIFICATIONS TO ALLOW OPTION OF SUB-
  !     *                         DIVIDING THIRD SOIL LAYER.
  !     * MAR 23/05 - D.VERSEGHY. ADD VARIABLES TO SUBROUTINE CALLS.
  !     * MAR 14/05 - D.VERSEGHY. RENAME SCAN TO SNCAN (RESERVED NAME
  !     *                         IN F90).
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 08/04 - D.VERSEGHY. NEW LOWER LIMITS FOR RCAN, SCAN, ZPOND
  !     *                         AND SNOW.
  !     * DEC 09/02 - D.VERSEGHY. SWITCH CALLING ORDER OF pondedWaterFreeze AND
  !     *                         snowSublimation FOR CONSISTENCY WITH DIAGNOSTICS.
  !     * SEP 26.02 - D.VERSEGHY. CHANGED CALL TO waterUnderCanopy.
  !     * AUG 01/02 - D.VERSEGHY. ADD CALL TO WATROF, NEW SUBROUTINE
  !     *                         CONTAINING WATERLOO OVERLAND FLOW
  !     *                         AND INTERFLOW CALCULATIONS.
  !     *                         SHORTENED CLASS3 COMMON BLOCK.
  !     * JUL 03/02 - D.VERSEGHY. STREAMLINE SUBROUTINE CALLS; MOVE
  !     *                         CALCULATION OF BACKGROUND SOIL
  !     *                         PROPERTIES INTO "soilProperties"; CHANGE
  !     *                         RHOSNI FROM CONSTANT TO VARIABLE.
  !     * OCT 04/01 - M.LAZARE.   REMOVE SEVERAL OLD DIAGNOSTIC FIELDS
  !     *                         AND ADD NEW FIELD "ROVG".
  !     * MAY 14/01 - M.LAZARE.   ADD CALLS TO SUBROUTINE "snowSublimation" FOR
  !     *                         FC AND FG SUBAREAS OF GRID CELL.
  !     * OCT 20/00 - D.VERSEGHY. ADD WORK ARRAY "RHOMAX" FOR snowAging.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         CHANGES RELATED TO VARIABLE SOIL DEPTH
  !     *                         (MOISTURE HOLDING CAPACITY) AND DEPTH-
  !     *                         VARYING SOIL PROPERTIES.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS; INTRODUCE CALCULATION OF
  !     *                         OVERLAND FLOW.
  !     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         VARIABLE SURFACE DETENTION CAPACITY
  !     *                         IMPLEMENTED.
  !     * AUG 24/95 - D.VERSEGHY. UPDATE ARRAY "EVAP" TO TAKE INTO
  !     *                         ACCOUNT "WLOST"; RATIONALIZE
  !     *                         CALCULATION OF THE LATTER.
  !     *                         COMPLETION OF WATER BUDGET DIAGNOSTICS.
  !     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
  !     *                         ORGANIC MATTER CONTENT.
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         CHANGES TO SUBROUTINE CALLS ASSOCIATED
  !     *                         WITH REVISIONS TO DIAGNOSTICS.
  !     *                         ALLOW SPECIFICATION OF LIMITING POND
  !     *                         DEPTH "PNDLIM" (PARALLEL CHANGES MADE
  !     *                         SIMULTANEOUSLY IN waterUpdates).
  !     * DEC 16/94 - D.VERSEGHY. TWO NEW DIAGNOSTIC FIELDS.
  !     * NOV 18/93 - D.VERSEGHY. LOCAL VERSION WITH INTERNAL WORK ARRAYS
  !     *                         HARD-CODED FOR USE ON PCS.
  !     * NOV 01/93 - D.VERSEGHY. CLASS - VERSION 2.2.
  !     *                         REVISIONS ASSOCIATED WITH NEW VERSION
  !     *                         OF waterUpdates.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. NUMEROUS NEW DIAGNOSTIC FIELDS.
  !     * MAY 06/93 - D.VERSEGHY/M.LAZARE. CORRECT BUG IN CALL TO waterUpdates
  !     *                                  FOR CANOPY-SNOW CASE, WHERE
  !     *                                  SHOULD BE PASSING "HCPCS"
  !     *                                  INSTEAD OF "HCPGS".
  !     * MAY 15/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !                               CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. LAND SURFACE WATER BUDGET CALCULATIONS.
  !
  use classicParams, only : DELT, TFREZ, HCPW, HCPSND, SPHW, SPHICE, RHOW

  implicit none

  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IWF(ILG)   !< Flag governing lateral soil water flow calculations
  integer, intent(in) :: ILG, IL1, IL2, JL, IC, IG, IGP1, IGP2, N
  integer, intent(in) :: NLANDCS   !< Number of modelled areas that contain subareas of canopy over snow
  integer, intent(in) :: NLANDGS   !< Number of modelled areas that contain subareas of snow
  integer, intent(in) :: NLANDC    !< Number of modelled areas that contain subareas of canopy over bare ground
  integer, intent(in) :: NLANDG    !< Number of modelled areas that contain subareas of bare ground
  integer, intent(in) :: NLANDI    !< Number of modelled areas that are ice sheets [ ]
  integer :: IPTBAD, JPTBAD, KPTBAD, LPTBAD, I, J
  !
  !     * MAIN OUTPUT FIELDS.
  !
  real, intent(inout) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3}]\f$
  real, intent(out) :: THICE (ILG,IG)   !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3}]\f$
  real, intent(inout) :: TBAR  (ILG,IG)   !< Temperature of soil layers [K]
  real, intent(out) :: GFLUX (ILG,IG)   !< Heat flux at interfaces between soil layers \f$[W m^{-2}]\f$
  !
  real, intent(inout) :: TCAN  (ILG)  !< Vegetation canopy temperature [K]
  real, intent(inout) :: RCAN  (ILG)  !< Intercepted liquid water stored on canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: SNCAN (ILG)  !< Intercepted frozen water stored on canopy \f$[kg m^{-2}]\f$
  real, intent(inout) :: RUNOFF(ILG)  !< Total runoff from soil \f$[m or kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: SNO   (ILG)  !< Mass of snow pack \f$[kg m^{-2}]\f$
  real, intent(inout) :: TSNOW (ILG)  !< Snowpack temperature [K]
  real, intent(inout) :: RHOSNO(ILG)  !< Density of snow \f$[kg m^{-3}]\f$
  real, intent(inout) :: ALBSNO(ILG)  !< Snow albedo [ ]
  real, intent(inout) :: ZPOND (ILG)  !< Depth of ponded water on surface [m]
  real, intent(inout) :: TPOND (ILG)  !< Temperature of ponded water [K]
  real, intent(inout) :: GROWTH(ILG)  !< Vegetation growth index [ ]
  real, intent(inout) :: TBASE (ILG)  !< Temperature of bedrock in third soil layer [K]
  real, intent(inout) :: TRUNOF(ILG)  !< Temperature of total runoff [K]
  real, intent(inout) :: WSNOW (ILG)  !< Liquid water content of snow pack \f$[kg m^{-2}]\f$
  !
  !     * DIAGNOSTIC ARRAYS.
  !
  real, intent(in) :: PCFC  (ILG)  !< Frozen precipitation intercepted by vegetation
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: PCLC  (ILG)  !< Liquid precipitation intercepted by vegetation
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: PCPN  (ILG)  !< Precipitation incident on snow pack
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCPG  (ILG)  !< Precipitation incident on ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: QFCF  (ILG)  !< Sublimation from frozen water on vegetation
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: QFCL  (ILG)  !< Evaporation from liquid water on vegetation
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: QFN   (ILG)  !< Sublimation from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: QFG   (ILG)  !< Evaporation from ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: HMFC  (ILG)  !< Diagnosed energy associated with phase change
  !< of water on vegetation \f$[W m^{-2}]\f$
  real, intent(in) :: HMFN  (ILG)  !< Diagnosed energy associated with phase change
  !< of water in snow pack \f$[W m^{-2}]\f$
  real, intent(inout) :: HTCC  (ILG)  !< Diagnosed internal energy change of vegetation
  !< canopy due to conduction and/or change in mass \f$[W m^{-2}]\f$
  real, intent(inout) :: HTCS  (ILG)  !< Diagnosed internal energy change of snow pack
  !< due to conduction and/or change in mass \f$[W m^{-2}]\f$
  real, intent(inout) :: ROFC  (ILG)  !< Liquid/frozen water runoff from vegetation
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: ROFN  (ILG)  !< Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: ROVG  (ILG)  !< Liquid/frozen water runoff from vegetation to
  !< ground surface \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: WTRS  (ILG)  !< Diagnosed residual water transferred into or
  !< out of the snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: WTRG  (ILG)  !< Diagnosed residual water transferred into or
  !< out of the soil \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: OVRFLW(ILG)  !< Overland flow from top of soil column
  !\f$[m or kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: SUBFLW(ILG,IG)  !< Interflow from sides of soil column
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: BASFLW(ILG)  !< Base flow from bottom of soil column
  !< \f$[m or kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: TOVRFL(ILG)  !< Temperature of overland flow from top of soil column [K]
  real, intent(in) :: TSUBFL(ILG,IG)  !< Temperature of interflow from sides of soil column [K]
  real, intent(in) :: TBASFL(ILG)  !< Temperature of base flow from bottom of soil column [K]
  real, intent(inout) :: EVAP  (ILG)  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFLUX (ILG)  !< Product of surface drag coefficient, wind speed and surface-air specific humidity difference \f$[m s^{-1}]\f$
  real, intent(in) :: RHOAIR (ILG) !< Density of air \f$[kg m^{-3}]\f$
  !
  real, intent(in) :: QFC  (ILG,IG) !< Water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: HMFG (ILG,IG) !< Diagnosed energy associated with phase change of water in soil layers \f$[W m^{-2}]\f$
  real, intent(inout) :: HTC  (ILG,IG) !< Diagnosed internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2}]\f$
  !
  !     * I/O FIELDS PASSED THROUGH CLASS.
  !
  !(In composite definitions, suffix C or CO = vegetation over
  ! ground; G or GO = bare ground; CS = vegetation over snow cover
  ! GS = bare snow cover.)
  !
  real, intent(in) :: RPCP  (ILG)  !< Rainfall rate over modelled area \f$[m s^{-1}]\f$
  real, intent(in) :: TRPCP (ILG)  !< Rainfall temperature over modelled area [C]
  real, intent(in) :: SPCP  (ILG)  !< Snowfall rate over modelled area \f$[m s^{-1}]\f$
  real, intent(in) :: TSPCP (ILG)  !< Snowfall temperature over modelled area [C]
  real, intent(in) :: PCPR  (ILG)  !< Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: TA    (ILG)  !< Air temperature at reference height [K]
  !
  real, intent(in) :: TBARC(ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(in) :: TBARG(ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(in) :: TBARCS(ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(in) :: TBARGS(ILG,IG) !< Subarea temperatures of soil layers [C]
  !
  real, intent(in) :: THLIQC(ILG,IG)   !< Liquid water content of soil layers under vegetation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLIQG(ILG,IG)   !< Liquid water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THICEC(ILG,IG)   !< Frozen water content of soil layers under vegetation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THICEG(ILG,IG)   !< Frozen water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
  real, intent(in) :: HCPC  (ILG,IG)   !< Heat capacity of soil layers under vegetation \f$[J m^{-3} K^{-1}]\f$
  real, intent(in) :: HCPG  (ILG,IG)   !< Heat capacity of soil layers in bare areas \f$[J m^{-3} K^{-1}]\f$
  real, intent(in) :: TCTOPC(ILG,IG)   !< Thermal conductivity of soil at top of layer (vegetation over ground) \f$[W m^{-1} K^{-1}]\f$
  real, intent(in) :: TCBOTC(ILG,IG)   !< Thermal conductivity of soil at bottom of layer (vegetation over ground) \f$[W m^{-1} K^{-1}]\f$
  real, intent(in) :: TCTOPG(ILG,IG)   !< Thermal conductivity of soil at top of layer (bare ground) \f$[W m^{-1} K^{-1}]\f$
  real, intent(in) :: TCBOTG(ILG,IG)   !< Thermal conductivity of soil at bottom of layer (bare ground) \f$[W m^{-1} K^{-1}]\f$
  real, intent(in) :: FROOT (ILG,IG)   !< Fraction of total transpiration contributed by soil layer over snow-free subarea [ ]
  real, intent(in) :: FROOTS (ILG,IG)  !< Fraction of total transpiration contributed by soil layer over snow-covered subarea [ ]
  real, intent(in) :: TSFSAV(ILG,4)    !< Ground surface temperature over subarea [K]
  !
  real, intent(in) :: FC    (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FG    (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FCS   (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FGS   (ILG) !< Subarea fractional coverage of modelled area [ ]
  !
  real, intent(in) :: TPONDC(ILG) !< Subarea temperature of surface ponded water [C]
  real, intent(in) :: TPONDG(ILG) !< Subarea temperature of surface ponded water [C]
  real, intent(in) :: TPNDCS(ILG) !< Subarea temperature of surface ponded water [C]
  real, intent(in) :: TPNDGS(ILG) !< Subarea temperature of surface ponded water [C]
  !
  real, intent(in) :: EVAPC (ILG)  !< Evaporation from vegetation over ground \f$[m s^{-1}]\f$
  real, intent(in) :: EVAPCG(ILG)  !< Evaporation from ground under vegetation \f$[m s^{-1}]\f$
  real, intent(in) :: EVAPG (ILG)  !< Evaporation from bare ground \f$[m s^{-1}]\f$
  real, intent(in) :: EVAPCS(ILG)  !< Evaporation from vegetation over snow \f$[m s^{-1}]\f$
  real, intent(in) :: EVPCSG(ILG)  !< Evaporation from snow under vegetation \f$[m s^{-1}]\f$
  real, intent(in) :: EVAPGS(ILG)  !< Evaporation from snow on bare ground \f$[m s^{-1}]\f$
  real, intent(in) :: QFREZC(ILG)  !< Heat sink to be used for freezing water on ground under canopy \f$[W m^{-2}]\f$
  real, intent(in) :: QFREZG(ILG)  !< Heat sink to be used for freezing water on bare ground \f$[W m^{-2}]\f$
  real, intent(in) :: QMELTC(ILG)  !< Heat to be used for melting snow under canopy \f$[W m^{-2}]\f$
  real, intent(in) :: QMELTG(ILG)  !< Heat to be used for melting snow on bare ground \f$[W m^{-2}]\f$
  real, intent(in) :: RAICAN(ILG)  !< Intercepted liquid water stored on canopy over ground \f$[kg m^{-2}]\f$
  real, intent(in) :: SNOCAN(ILG)  !< Intercepted frozen water stored on canopy over ground \f$[kg m^{-2}]\f$
  real, intent(in) :: RAICNS(ILG)  !< Intercepted liquid water stored on canopy over snow \f$[kg m^{-2}]\f$
  real, intent(in) :: SNOCNS(ILG)  !< Intercepted frozen water stored on canopy over snow \f$[kg m^{-2}]\f$
  real, intent(in) :: FSVF  (ILG)  !< Sky view factor of ground under vegetation canopy [ ]
  real, intent(in) :: FSVFS (ILG)  !< Sky view factor of snow under vegetation canopy [ ]
  real, intent(in) :: CWLCAP(ILG)  !< Storage capacity of canopy over bare ground for liquid water \f$[kg m^{-2}]\f$
  real, intent(in) :: CWFCAP(ILG)  !< Storage capacity of canopy over bare ground for frozen water \f$[kg m^{-2}]\f$
  real, intent(in) :: CWLCPS(ILG)  !< Storage capacity of canopy over snow for liquid water \f$[kg m^{-2}]\f$
  real, intent(in) :: CWFCPS(ILG)  !< Storage capacity of canopy over snow for frozen water \f$[kg m^{-2}]\f$
  real, intent(in) :: TCANO (ILG)  !< Temperature of canopy over ground [K]
  real, intent(in) :: TCANS (ILG)  !< Temperature of canopy over snow [K]
  real, intent(in) :: CHCAP (ILG)  !< Heat capacity of canopy over bare ground \f$[J m^{-2} K^{-1}] \f$
  real, intent(in) :: CHCAPS(ILG)  !< Heat capacity of canopy over snow \f$[J m^{-2} K^{-1}] \f$
  real, intent(in) :: CMASSC(ILG)  !< Mass of canopy over bare ground \f$[kg m^{-2}]\f$
  real, intent(in) :: CMASCS(ILG)  !< Mass of canopy over snow \f$[kg m^{-2}]\f$
  real, intent(inout) :: ZSNOW (ILG)  !< Depth of snow pack [m]
  real, intent(in) :: RHOSNI(ILG)  !< Density of fresh snow \f$[kg m^{-3}]\f$
  real, intent(in) :: GZEROC(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(in) :: GZEROG(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(in) :: GZROCS(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(in) :: GZROGS(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  !
  real, intent(in) :: G12C  (ILG)  !< Subarea heat flux between first and second soil layers \f$[W m^{-2}]\f$
  real, intent(in) :: G12G  (ILG)  !< Subarea heat flux between first and second soil layers \f$[W m^{-2}]\f$
  real, intent(in) :: G12CS (ILG)  !< Subarea heat flux between first and second soil layers \f$[W m^{-2}]\f$
  real, intent(in) :: G12GS (ILG)  !< Subarea heat flux between first and second soil layers \f$[W m^{-2}]\f$
  !
  real, intent(in) :: G23C  (ILG)  !< Subarea heat flux between second and third soil layers \f$[W m^{-2}]\f$
  real, intent(in) :: G23G  (ILG)  !< Subarea heat flux between second and third soil layers \f$[W m^{-2}]\f$
  real, intent(in) :: G23CS (ILG)  !< Subarea heat flux between second and third soil layers \f$[W m^{-2}]\f$
  real, intent(in) :: G23GS (ILG)  !< Subarea heat flux between second and third soil layers \f$[W m^{-2}]\f$
  !
  real, intent(in) :: TSNOCS(ILG)  !< Temperature of snow pack under vegetation [K]
  real, intent(in) :: TSNOGS(ILG)  !< Temperature of snow pack in bare areas [K]
  real, intent(in) :: WSNOCS(ILG)  !< Liquid water content of snow pack under vegetation \f$[kg m^{-2}]\f$
  real, intent(in) :: WSNOGS(ILG)  !< Liquid water content of snow pack in bare areas \f$[kg m^{-2}]\f$
  real, intent(in) :: RHOSCS(ILG)  !< Density of snow under vegetation \f$[kg m^{-3}]\f$
  real, intent(in) :: RHOSGS(ILG)  !< Density of snow in bare areas \f$[kg m^{-3}]\f$
  real, intent(in) :: ZPLIMC(ILG)  !< Subarea maximum ponding depth [m]
  real, intent(in) :: ZPLIMG(ILG)  !< Subarea maximum ponding depth [m]
  real, intent(in) :: ZPLMCS(ILG)  !< Subarea maximum ponding depth [m]
  real, intent(in) :: ZPLMGS(ILG)  !< Subarea maximum ponding depth [m]
  real, intent(in) :: GGEO  (ILG)  !< Geothermal heat flux at bottom of soil profile \f$[W m^{-2}]\f$
  !
  !     * SOIL PROPERTY ARRAYS.
  !
  real, intent(in) :: THPOR (ILG,IG)   !< Pore volume in soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLRET(ILG,IG)   !< Liquid water retention capacity for organic soil [m3 m-3 ]
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !< remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: BI    (ILG,IG)   !< Clapp and Hornberger empirical "b" parameter [ ]
  real, intent(in) :: GRKSAT(ILG,IG)   !< Saturated hydraulic conductivity of soil layer \f$[m s^{-1}]\f$
  real, intent(in) :: PSISAT(ILG,IG)   !< Soil moisture suction at saturation [m]
  real, intent(in) :: THLRAT(ILG,IG)   !< Fractional saturation of soil behind the wetting front [ ]
  real, intent(in) :: THFC  (ILG,IG)   !< Field capacity \f$[m^3 m^{-3}]\f$
  real, intent(in) :: HCPS  (ILG,IG)   !< Heat capacity of soil material \f$[J m^{-3} K^{-1}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Overall thickness of soil layer [m]
  real :: DELZZ (ILG,IG)   !< Permeable thickness of soil layer [m]
  real, intent(in) :: ZBOTW (ILG,IG)   !< Depth to permeable bottom of soil layer [m]
  real, intent(in) :: XDRAIN(ILG)      !< Drainage index at bottom of soil profile [ ]
  real, intent(in) :: XSLOPE(ILG)      !< Surface slope (used when running MESH code) [degrees]
  real, intent(in) :: XDRAINH(ILG)     !< WATROF parameter used when running MESH code [ ]
  real, intent(in) :: WFSURF(ILG)      !< WATROF parameter used when running MESH code [ ] (not used)
  real, intent(in) :: KSAT(ILG)        !< WATROF parameter used when running MESH code [ ]
  real, intent(in) :: DD(ILG)          !< WATROF parameter used when running MESH code [ ]
  real, intent(in) :: MANNING_N(ILG)   !< WATROF parameter used when running MESH code [ ]
  real, intent(in) :: DELZ  (IG)       !< Overall thickness of soil layer [m]
  !
  integer, intent(in) :: ISAND(ILG,IG) !< Sand content flag
  integer, intent(in) :: IGDR  (ILG)   !< Index of soil layer in which bedrock is encountered
  !
  !     * INTERNAL WORK ARRAYS USED THROUGHOUT waterBudgetDriver.
  !
  real :: TBARWC(ILG,IG), TBARWG(ILG,IG), TBRWCS(ILG,IG), TBRWGS(ILG,IG), &
          THLQCO(ILG,IG), THLQGO(ILG,IG), THLQCS(ILG,IG), THLQGS(ILG,IG), &
          THICCO(ILG,IG), THICGO(ILG,IG), THICCS(ILG,IG), THICGS(ILG,IG), &
          HCPCO (ILG,IG), HCPGO (ILG,IG), HCPCS (ILG,IG), HCPGS (ILG,IG), &
          GRKSC (ILG,IG), GRKSG (ILG,IG), GRKSCS(ILG,IG), GRKSGS(ILG,IG), &
          GFLXC (ILG,IG), GFLXG (ILG,IG), GFLXCS(ILG,IG), GFLXGS(ILG,IG)
  !
  real :: SPCC  (ILG), SPCG  (ILG), SPCCS (ILG), SPCGS (ILG), &
          TSPCC (ILG), TSPCG (ILG), TSPCCS(ILG), TSPCGS(ILG), &
          RPCC  (ILG), RPCG  (ILG), RPCCS (ILG), RPCGS (ILG), &
          TRPCC (ILG), TRPCG (ILG), TRPCCS(ILG), TRPCGS(ILG), &
          EVPIC (ILG), EVPIG (ILG), EVPICS(ILG), EVPIGS(ILG), &
          ZPONDC(ILG), ZPONDG(ILG), ZPNDCS(ILG), ZPNDGS(ILG), &
          XSNOWC(ILG), XSNOWG(ILG), XSNOCS(ILG), XSNOGS(ILG), &
          ZSNOWC(ILG), ZSNOWG(ILG), ZSNOCS(ILG), ZSNOGS(ILG), &
          ALBSC (ILG), ALBSG (ILG), ALBSCS(ILG), ALBSGS(ILG), &
          RHOSC (ILG), RHOSG (ILG), &
          HCPSC (ILG), HCPSG (ILG), HCPSCS(ILG), HCPSGS(ILG), &
          RUNFC (ILG), RUNFG (ILG), RUNFCS(ILG), RUNFGS(ILG), &
          TRUNFC(ILG), TRUNFG(ILG), TRNFCS(ILG), TRNFGS(ILG), &
          TBASC (ILG), TBASG (ILG), TBASCS(ILG), TBASGS(ILG)
  !
  real :: SUBLC (ILG), SUBLCS(ILG), WLOSTC(ILG), WLOSTG(ILG), &
          WLSTCS(ILG), WLSTGS(ILG), RAC   (ILG), RACS  (ILG), &
          SNC   (ILG), SNCS  (ILG), TSNOWC(ILG), TSNOWG(ILG), &
          DT    (ILG), ZERO  (ILG), RALB  (ILG), ZFAV  (ILG), &
          THLINV(ILG)
  !
  integer :: LZFAV (ILG)
  !
  !     * INTERNAL WORK ARRAYS FOR waterCalcPrep AND canopyInterception.
  !
  real :: RADD(ILG), SADD(ILG)
  !
  !     * INTERNAL WORK FIELDS FOR waterFlowInfiltrate/waterFlowNonInfiltrate (AND THEIR CALLED
  !     * ROUTINES (I.E. waterInfiltrateUnsat, waterInfiltrateSat, waterBaseflow) AND iceSheetBalance.
  !
  real :: ZMAT(ILG,IGP2,IGP1)
  !
  real :: WMOVE(ILG,IGP2),  TMOVE (ILG,IGP2)
  !
  real :: THLIQX(ILG,IGP1), THICEX(ILG,IGP1), TBARWX(ILG,IGP1), &
          DELZX (ILG,IGP1), ZBOTX (ILG,IGP1), FDT   (ILG,IGP1), &
          TFDT  (ILG,IGP1), PSIF  (ILG,IGP1), THLINF(ILG,IGP1), &
          GRKINF(ILG,IGP1), FDUMMY(ILG,IGP1), TDUMMY(ILG,IGP1), &
          ZRMDR (ILG,IGP1)
  !
  real :: THLMAX(ILG,IG), THTEST(ILG,IG), THLDUM(ILG,IG), &
          THIDUM(ILG,IG), TDUMW (ILG,IG)
  !
  real :: TRMDR (ILG), ZF    (ILG), FMAX  (ILG), TUSED (ILG), &
          RDUMMY(ILG), WEXCES(ILG), FDTBND(ILG), WADD  (ILG), &
          TADD  (ILG), WADJ  (ILG), TIMPND(ILG), DZF   (ILG), &
          DTFLOW(ILG), THLNLZ(ILG), THLQLZ(ILG), DZDISP(ILG), &
          WDISP (ILG), WABS  (ILG), ZMOVE (ILG), TBOT  (ILG)
  !
  integer :: IGRN  (ILG), IGRD  (ILG), IZERO (ILG), &
             IFILL (ILG), LZF   (ILG), NINF  (ILG), &
             IFIND (ILG), ITER  (ILG), NEND  (ILG), &
             ISIMP (ILG), ICONT (ILG)
  !
  !     * INTERNAL WORK ARRAYS FOR canopyWaterUpdate AND snowAging.
  !
  real :: EVLOST(ILG), RLOST (ILG), RHOMAX(ILG)
  !
  integer              :: IROOT (ILG)
  !
  !     * INTERNAL WORK ARRAYS FOR WATROF.
  !
  real :: THCRIT(ILG,IG), DODRN(ILG), DOVER(ILG), DIDRN(ILG,IG), &
          DIDRNMX(ILG,IG), BULK_FC(ILG,IG)
  !
  !     * INTERNAL WORK ARRAYS FOR checkWaterBudget.
  !
  real :: BAL(ILG)
  !
  !     * INTERNAL SCALARS.
  !
  real :: SNOROF, WSNROF
  !
  !-----------------------------------------------------------------------
  !>
  !! First, subroutine waterCalcPrep is called to initialize various arrays
  !! and produce parameters for the four subareas of canopy over snow
  !! (CS), snow on ground (GS), canopy over ground (C) and bare ground
  !! (G). Then, for each of the four subareas, if the number of
  !! modelled areas containing that subarea is greater than zero, a
  !! series of subroutines is called. The subroutines associated with
  !! each subarea are listed in the table below.
  !!
  !! \f[
  !! \begin{array} { | l | l | c | }
  !! \hline
  !! \text{canopyWaterUpdate}  & \text{Evaporation/sublimation of water from vegetation canopy}        &   \text{CS,C}    \\ \hline
  !! \text{canopyInterception}  & \text{Addition of rainfall/snowfall to canopy; throughfall and drip}  &   \text{CS,C}    \\ \hline
  !! \text{canopyPhaseChange}  & \text{Freezing/thawing of liquid/frozen water on canopy}              &   \text{CS,C}    \\ \hline
  !! \text{waterUnderCanopy}  & \text{Precipitaiton and condensation under canopy}                    &   \text{CS,C}    \\ \hline
  !! \text{soilWaterPhaseChg}  & \text{Freezing/thawing of liquid/frozen water in soil}                & \text{CS,GS,C,G} \\ \hline
  !! \text{snowSublimation}  & \text{Sublimaiton from snow pack}                                     & \text{CS,GS,C,G} \\ \hline
  !! \text{pondedWaterFreeze}  & \text{Freezing of ponded water on soil}                               & \text{CS,GS,C,G} \\ \hline
  !! \text{snowMelt}   & \text{Melting of snow pack}                                           &   \text{CS,GS}   \\ \hline
  !! \text{snowAddNew}  & \text{Accumulation of snow on ground}                                 & \text{CS,GS,C,G} \\ \hline
  !! \text{snowInfiltrateRipen}  & \text{Infiltration of rain into snow pack}                            &   \text{CS,GS}   \\ \hline
  !! \text{iceSheetBalance}  & \text{Energy and water budget of ice sheets}                          &   \text{GS,G}    \\ \hline
  !! \text{waterFlowInfiltrate}  & \text{Infiltraiton of water into soil}                                & \text{CS,GS,C,G} \\ \hline
  !! \text{waterFlowNonInfiltrate}  & \text{Soil water movement in response to gravity and suction forces}  & \text{CS,GS,C,G} \\ \hline
  !! \text{waterUpdates}  & \text{Step ahead soil layer temperatures, check for freezing/thawing} & \text{CS,GS,C,G} \\ \hline
  !! \text{checkWaterBudget}  & \text{Check subarea moisture balances for closure}                    & \text{CS,GS,C,G} \\ \hline
  !! \text{snowAging} & \text{Temporal variation of snow albedo and density}                  &   \text{CS,GS}   \\ \hline
  !! \end{array}
  !! \f]
  !     * PREPARATION.
  !

  call waterCalcPrep(THLQCO, THLQGO, THLQCS, THLQGS, THICCO, THICGO, & ! Formerly WPREP
                     THICCS, THICGS, HCPCO, HCPGO, HCPCS, HCPGS, &
                     GRKSC, GRKSG, GRKSCS, GRKSGS, &
                     SPCC, SPCG, SPCCS, SPCGS, TSPCC, TSPCG, &
                     TSPCCS, TSPCGS, RPCC, RPCG, RPCCS, RPCGS, &
                     TRPCC, TRPCG, TRPCCS, TRPCGS, EVPIC, EVPIG, &
                     EVPICS, EVPIGS, ZPONDC, ZPONDG, ZPNDCS, ZPNDGS, &
                     XSNOWC, XSNOWG, XSNOCS, XSNOGS, ZSNOWC, ZSNOWG, &
                     ZSNOCS, ZSNOGS, ALBSC, ALBSG, ALBSCS, ALBSGS, &
                     RHOSC, RHOSG, HCPSC, HCPSG, HCPSCS, HCPSGS, &
                     RUNFC, RUNFG, RUNFCS, RUNFGS, &
                     TRUNFC, TRUNFG, TRNFCS, TRNFGS, TBASC, TBASG, &
                     TBASCS, TBASGS, GFLXC, GFLXG, GFLXCS, GFLXGS, &
                     SUBLC, SUBLCS, WLOSTC, WLOSTG, WLSTCS, WLSTGS, &
                     RAC, RACS, SNC, SNCS, TSNOWC, TSNOWG, &
                     OVRFLW, SUBFLW, BASFLW, TOVRFL, TSUBFL, TBASFL, &
                     PCFC, PCLC, PCPN, PCPG, QFCF, QFCL, &
                     QFN, QFG, QFC, HMFG, &
                     ROVG, ROFC, ROFN, TRUNOF, &
                     THLIQX, THICEX, THLDUM, THIDUM, &
                     DT, RDUMMY, ZERO, IZERO, DELZZ, &
                     FC, FG, FCS, FGS, &
                     THLIQC, THLIQG, THICEC, THICEG, HCPC, HCPG, &
                     TBARC, TBARG, TBARCS, TBARGS, TBASE, TSFSAV, &
                     FSVF, FSVFS, RAICAN, SNOCAN, RAICNS, SNOCNS, &
                     EVAPC, EVAPCG, EVAPG, EVAPCS, EVPCSG, EVAPGS, &
                     RPCP, TRPCP, SPCP, TSPCP, RHOSNI, ZPOND, &
                     ZSNOW, ALBSNO, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                     THPOR, HCPS, GRKSAT, ISAND, DELZW, DELZ, &
                     ILG, IL1, IL2, JL, IG, IGP1, &
                     NLANDCS, NLANDGS, NLANDC, NLANDG, RADD, SADD, &
                     BI, PSISAT, DD, XSLOPE, BULK_FC)
  !
  !
  !     * CALCULATIONS FOR CANOPY OVER SNOW.
  !
  if (NLANDCS > 0) then
    call canopyWaterUpdate(EVAPCS, SUBLCS, RAICNS, SNOCNS, TCANS, THLQCS, & ! Formerly CANVAP
                           TBARCS, ZSNOCS, WLSTCS, CHCAPS, QFCF, QFCL, QFN, QFC, &
                           HTCC, HTCS, HTC, FCS, CMASCS, TSNOCS, HCPSCS, RHOSCS, &
                           FROOTS, THPOR, THLMIN, DELZW, EVLOST, RLOST, IROOT, &
                           IG, ILG, IL1, IL2, JL, N)
    call canopyInterception(2, RPCCS, TRPCCS, SPCCS, TSPCCS, RAICNS, SNOCNS, & ! Formerly CANADD
                            TCANS, CHCAPS, HTCC, ROFC, ROVG, PCPN, PCPG, &
                            FCS, FSVFS, CWLCPS, CWFCPS, CMASCS, RHOSNI, &
                            TSFSAV(1, 1), RADD, SADD, ILG, IL1, IL2, JL)
    call canopyPhaseChange(TCANS, RAICNS, SNOCNS, RDUMMY, RDUMMY, CHCAPS, & ! Formerly CWCALC
                           HMFC, HTCC, FCS, CMASCS, ILG, IL1, IL2, JL)
    call waterUnderCanopy(2, RPCCS, TRPCCS, SPCCS, TSPCCS, RHOSNI, EVPCSG, & ! Formerly SUBCAN
                          QFN, QFG, PCPN, PCPG, FCS, ILG, IL1, IL2, JL)
    call soilWaterPhaseChg(TBARCS, THLQCS, THICCS, HCPCS, TBRWCS, HMFG, HTC, & ! Formerly TWCALC
                           FCS, ZERO, THPOR, THLMIN, HCPS, DELZW, DELZZ, ISAND, &
                           IG, ILG, IL1, IL2, JL)
    call snowSublimation(RHOSCS, ZSNOCS, HCPSCS, TSNOCS, EVPCSG, QFN, QFG, & ! Formerly SNOVAP
                         HTCS, WLSTCS, TRNFCS, RUNFCS, TOVRFL, OVRFLW, &
                         FCS, RPCCS, SPCCS, RHOSNI, WSNOCS, ILG, IL1, IL2, JL)
    call pondedWaterFreeze(ZPNDCS, TPNDCS, ZSNOCS, TSNOCS, ALBSCS, & ! Formerly TFREEZ
                           RHOSCS, HCPSCS, GZROCS, HMFG, HTCS, HTC, &
                           WTRS, WTRG, FCS, ZERO, WSNOCS, TA, TBARCS, &
                           ISAND, IG, ILG, IL1, IL2, JL)
    call snowMelt(ZSNOCS, TSNOCS, QMELTC, RPCCS, TRPCCS, & ! Formerly TMELT
                  GZROCS, RALB, HMFN, HTCS, HTC, FCS, HCPSCS, &
                  RHOSCS, WSNOCS, ISAND, IG, ILG, IL1, IL2, JL)
    call snowAddNew(ALBSCS, TSNOCS, RHOSCS, ZSNOCS, & ! Formerly SNOADD
                    HCPSCS, HTCS, FCS, SPCCS, TSPCCS, RHOSNI, WSNOCS, &
                    ILG, IL1, IL2, JL)
    call snowInfiltrateRipen(RPCCS, TRPCCS, ZSNOCS, TSNOCS, RHOSCS, HCPSCS, & ! Formerly SNINFL
                             WSNOCS, HTCS, HMFN, PCPG, ROFN, FCS, ILG, IL1, IL2, JL)
    call waterFlowInfiltrate(1, THLQCS, THICCS, TBRWCS, BASFLW, TBASFL, RUNFCS, & ! Formerly GRINFL
                             TRNFCS, ZFAV, LZFAV, THLINV, QFG, WLSTCS, &
                             FCS, EVPCSG, RPCCS, TRPCCS, TPNDCS, ZPNDCS, &
                             DT, ZMAT, WMOVE, TMOVE, THLIQX, THICEX, TBARWX, &
                             DELZX, ZBOTX, FDT, TFDT, PSIF, THLINF, GRKINF, &
                             THLMAX, THTEST, ZRMDR, FDUMMY, TDUMMY, THLDUM, &
                             THIDUM, TDUMW, TRMDR, ZF, FMAX, TUSED, RDUMMY, &
                             ZERO, WEXCES, FDTBND, WADD, TADD, WADJ, TIMPND, &
                             DZF, DTFLOW, THLNLZ, THLQLZ, DZDISP, WDISP, WABS, &
                             THPOR, THLRET, THLMIN, BI, PSISAT, GRKSCS, &
                             THLRAT, THFC, DELZW, ZBOTW, XDRAIN, DELZ, ISAND, &
                             IGRN, IGRD, IFILL, IZERO, LZF, NINF, IFIND, ITER, &
                             NEND, ISIMP, IGDR, &
                             IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call waterFlowNonInfiltrate(1, THLQCS, THICCS, TBRWCS, FDUMMY, TDUMMY, & ! Formerly GRDRAN
                                BASFLW, TBASFL, RUNFCS, TRNFCS, &
                                QFG, WLSTCS, FCS, EVPCSG, RPCCS, ZPNDCS, &
                                DT, WEXCES, THLMAX, THTEST, THPOR, THLRET, THLMIN, &
                                BI, PSISAT, GRKSCS, THFC, DELZW, XDRAIN, ISAND, &
                                IZERO, IGRN, IGRD, IGDR, &
                                IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call WATROF(THLQCS, THICCS, ZPNDCS, TPNDCS, OVRFLW, TOVRFL, &
                SUBFLW, TSUBFL, RUNFCS, TRNFCS, FCS, ZPLMCS, &
                XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBRWCS, &
                DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN, &
                ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
    call waterUpdates(TBARCS, THLQCS, THICCS, HCPCS, TPNDCS, ZPNDCS, & ! Formerly TMCALC
                      TSNOCS, ZSNOCS, ALBSCS, RHOSCS, HCPSCS, TBASCS, &
                      OVRFLW, TOVRFL, RUNFCS, TRNFCS, HMFG, HTC, HTCS, &
                      WTRS, WTRG, FCS, TBRWCS, GZROCS, G12CS, &
                      G23CS, GGEO, TA, WSNOCS, TCTOPC, TCBOTC, GFLXCS, &
                      ZPLMCS, THPOR, THLMIN, HCPS, DELZW, DELZZ, DELZ, &
                      ISAND, IWF, IG, ILG, IL1, IL2, JL, N)
    call checkWaterBudget(1, PCPR, EVPICS, RUNFCS, WLSTCS, RAICNS, SNOCNS, & ! Formerly CHKWAT
                          RACS, SNCS, ZPNDCS, ZPOND, THLQCS, THICCS, &
                          THLIQC, THICEC, ZSNOCS, RHOSCS, XSNOCS, SNO, &
                          WSNOCS, WSNOW, FCS, FGS, FCS, BAL, THPOR, THLMIN, &
                          DELZW, ISAND, IG, ILG, IL1, IL2, JL, N)
    call snowAging(ALBSCS, RHOSCS, ZSNOCS, HCPSCS, & ! Formerly SNOALBW
                   TSNOCS, FCS, SPCCS, RALB, WSNOCS, RHOMAX, &
                   ISAND, ILG, IG, IL1, IL2, JL)
  end if
  !
  !     * CALCULATIONS FOR SNOW-COVERED GROUND.
  !
  if (NLANDGS > 0) then
    call soilWaterPhaseChg(TBARGS, THLQGS, THICGS, HCPGS, TBRWGS, HMFG, HTC, & ! Formerly TWCALC
                           FGS, ZERO, THPOR, THLMIN, HCPS, DELZW, DELZZ, ISAND, &
                           IG, ILG, IL1, IL2, JL)
    call snowSublimation(RHOSGS, ZSNOGS, HCPSGS, TSNOGS, EVAPGS, QFN, QFG, & ! Formerly SNOVAP
                         HTCS, WLSTGS, TRNFGS, RUNFGS, TOVRFL, OVRFLW, &
                         FGS, RPCGS, SPCGS, RHOSNI, WSNOGS, ILG, IL1, IL2, JL)
    call pondedWaterFreeze(ZPNDGS, TPNDGS, ZSNOGS, TSNOGS, ALBSGS, & ! Formerly TFREEZ
                           RHOSGS, HCPSGS, GZROGS, HMFG, HTCS, HTC, &
                           WTRS, WTRG, FGS, ZERO, WSNOGS, TA, TBARGS, &
                           ISAND, IG, ILG, IL1, IL2, JL)
    call snowMelt(ZSNOGS, TSNOGS, QMELTG, RPCGS, TRPCGS, & ! Formerly TMELT
                  GZROGS, RALB, HMFN, HTCS, HTC, FGS, HCPSGS, &
                  RHOSGS, WSNOGS, ISAND, IG, ILG, IL1, IL2, JL)
    call snowAddNew(ALBSGS, TSNOGS, RHOSGS, ZSNOGS, & ! Formerly SNOADD
                    HCPSGS, HTCS, FGS, SPCGS, TSPCGS, RHOSNI, WSNOGS, &
                    ILG, IL1, IL2, JL)
    call snowInfiltrateRipen(RPCGS, TRPCGS, ZSNOGS, TSNOGS, RHOSGS, HCPSGS, & ! Formerly SNINFL
                             WSNOGS, HTCS, HMFN, PCPG, ROFN, FGS, ILG, IL1, IL2, JL)
    if (NLANDI /= 0) then
      call iceSheetBalance(TBARGS, TPNDGS, ZPNDGS, TSNOGS, RHOSGS, ZSNOGS, & ! Formerly ICEBAL
                           HCPSGS, ALBSGS, HMFG, HTCS, HTC, WTRS, WTRG, GFLXGS, &
                           RUNFGS, TRNFGS, OVRFLW, TOVRFL, ZPLMGS, GGEO, &
                           FGS, EVAPGS, RPCGS, TRPCGS, GZROGS, G12GS, G23GS, &
                           HCPGS, QMELTG, WSNOGS, ZMAT, TMOVE, WMOVE, ZRMDR, &
                           TADD, ZMOVE, TBOT, DELZ, ISAND, ICONT, &
                           IWF, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    end if
    call waterFlowInfiltrate(2, THLQGS, THICGS, TBRWGS, BASFLW, TBASFL, RUNFGS, & ! Formerly GRINFL
                             TRNFGS, ZFAV, LZFAV, THLINV, QFG, WLSTGS, &
                             FGS, EVAPGS, RPCGS, TRPCGS, TPNDGS, ZPNDGS, &
                             DT, ZMAT, WMOVE, TMOVE, THLIQX, THICEX, TBARWX, &
                             DELZX, ZBOTX, FDT, TFDT, PSIF, THLINF, GRKINF, &
                             THLMAX, THTEST, ZRMDR, FDUMMY, TDUMMY, THLDUM, &
                             THIDUM, TDUMW, TRMDR, ZF, FMAX, TUSED, RDUMMY, &
                             ZERO, WEXCES, FDTBND, WADD, TADD, WADJ, TIMPND, &
                             DZF, DTFLOW, THLNLZ, THLQLZ, DZDISP, WDISP, WABS, &
                             THPOR, THLRET, THLMIN, BI, PSISAT, GRKSGS, &
                             THLRAT, THFC, DELZW, ZBOTW, XDRAIN, DELZ, ISAND, &
                             IGRN, IGRD, IFILL, IZERO, LZF, NINF, IFIND, ITER, &
                             NEND, ISIMP, IGDR, &
                             IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call waterFlowNonInfiltrate(2, THLQGS, THICGS, TBRWGS, FDUMMY, TDUMMY, & ! Formerly GRDRAN
                                BASFLW, TBASFL, RUNFGS, TRNFGS, &
                                QFG, WLSTGS, FGS, EVAPGS, RPCGS, ZPNDGS, &
                                DT, WEXCES, THLMAX, THTEST, THPOR, THLRET, THLMIN, &
                                BI, PSISAT, GRKSGS, THFC, DELZW, XDRAIN, ISAND, &
                                IZERO, IGRN, IGRD, IGDR, &
                                IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call WATROF(THLQGS, THICGS, ZPNDGS, TPNDGS, OVRFLW, TOVRFL, &
                SUBFLW, TSUBFL, RUNFGS, TRNFGS, FGS, ZPLMGS, &
                XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBRWGS, &
                DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN, &
                ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
    call waterUpdates(TBARGS, THLQGS, THICGS, HCPGS, TPNDGS, ZPNDGS, & ! Formerly TMCALC
                      TSNOGS, ZSNOGS, ALBSGS, RHOSGS, HCPSGS, TBASGS, &
                      OVRFLW, TOVRFL, RUNFGS, TRNFGS, HMFG, HTC, HTCS, &
                      WTRS, WTRG, FGS, TBRWGS, GZROGS, G12GS, &
                      G23GS, GGEO, TA, WSNOGS, TCTOPG, TCBOTG, GFLXGS, &
                      ZPLMGS, THPOR, THLMIN, HCPS, DELZW, DELZZ, DELZ, &
                      ISAND, IWF, IG, ILG, IL1, IL2, JL, N)
    call checkWaterBudget(2, PCPR, EVPIGS, RUNFGS, WLSTGS, RAICNS, SNOCNS, &
                          RACS, SNCS, ZPNDGS, ZPOND, THLQGS, THICGS, &
                          THLIQG, THICEG, ZSNOGS, RHOSGS, XSNOGS, SNO, &
                          WSNOGS, WSNOW, FCS, FGS, FGS, BAL, THPOR, THLMIN, &
                          DELZW, ISAND, IG, ILG, IL1, IL2, JL, N)
    call snowAging(ALBSGS, RHOSGS, ZSNOGS, HCPSGS, & ! Formerly SNOALBW
                   TSNOGS, FGS, SPCGS, RALB, WSNOGS, RHOMAX, &
                   ISAND, ILG, IG, IL1, IL2, JL)
  end if
  !
  !     * CALCULATIONS FOR CANOPY OVER BARE GROUND.
  !
  if (NLANDC > 0) then
    call canopyWaterUpdate(EVAPC, SUBLC, RAICAN, SNOCAN, TCANO, THLQCO, & ! Formerly CANVAP
                           TBARC, ZSNOWC, WLOSTC, CHCAP, QFCF, QFCL, QFN, QFC, &
                           HTCC, HTCS, HTC, FC, CMASSC, TSNOWC, HCPSC, RHOSC, &
                           FROOT, THPOR, THLMIN, DELZW, EVLOST, RLOST, IROOT, &
                           IG, ILG, IL1, IL2, JL, N)
    call canopyInterception(1, RPCC, TRPCC, SPCC, TSPCC, RAICAN, SNOCAN, & ! Formerly CANADD
                            TCANO, CHCAP, HTCC, ROFC, ROVG, PCPN, PCPG, &
                            FC, FSVF, CWLCAP, CWFCAP, CMASSC, RHOSNI, &
                            TSFSAV(1,3), RADD, SADD, ILG, IL1, IL2, JL)
    call canopyPhaseChange(TCANO, RAICAN, SNOCAN, RDUMMY, RDUMMY, CHCAP, & ! Formerly CWCALC
                           HMFC, HTCC, FC, CMASSC, ILG, IL1, IL2, JL)
    call waterUnderCanopy(1, RPCC, TRPCC, SPCC, TSPCC, RHOSNI, EVAPCG, & ! Formerly SUBCAN
                          QFN, QFG, PCPN, PCPG, FC, ILG, IL1, IL2, JL)
    call soilWaterPhaseChg(TBARC, THLQCO, THICCO, HCPCO, TBARWC, HMFG, HTC, & ! Formerly TWCALC
                           FC, EVAPCG, THPOR, THLMIN, HCPS, DELZW, DELZZ, &
                           ISAND, IG, ILG, IL1, IL2, JL)
    call snowSublimation(RHOSC, ZSNOWC, HCPSC, TSNOWC, EVAPCG, QFN, QFG, & ! Formerly SNOVAP
                         HTCS, WLOSTC, TRUNFC, RUNFC, TOVRFL, OVRFLW, &
                         FC, RPCC, SPCC, RHOSNI, ZERO, ILG, IL1, IL2, JL)
    call pondedWaterFreeze(ZPONDC, TPONDC, ZSNOWC, TSNOWC, ALBSC, & ! Formerly TFREEZ
                           RHOSC, HCPSC, GZEROC, HMFG, HTCS, HTC, &
                           WTRS, WTRG, FC, QFREZC, ZERO, TA, TBARC, &
                           ISAND, IG, ILG, IL1, IL2, JL)
    call snowAddNew(ALBSC, TSNOWC, RHOSC, ZSNOWC, & ! Formerly SNOADD
                    HCPSC, HTCS, FC, SPCC, TSPCC, RHOSNI, ZERO, &
                    ILG, IL1, IL2, JL)
    call waterFlowInfiltrate(3, THLQCO, THICCO, TBARWC, BASFLW, TBASFL, RUNFC, & ! Formerly GRINFL
                             TRUNFC, ZFAV, LZFAV, THLINV, QFG, WLOSTC, &
                             FC, EVAPCG, RPCC, TRPCC, TPONDC, ZPONDC, &
                             DT, ZMAT, WMOVE, TMOVE, THLIQX, THICEX, TBARWX, &
                             DELZX, ZBOTX, FDT, TFDT, PSIF, THLINF, GRKINF, &
                             THLMAX, THTEST, ZRMDR, FDUMMY, TDUMMY, THLDUM, &
                             THIDUM, TDUMW, TRMDR, ZF, FMAX, TUSED, RDUMMY, &
                             ZERO, WEXCES, FDTBND, WADD, TADD, WADJ, TIMPND, &
                             DZF, DTFLOW, THLNLZ, THLQLZ, DZDISP, WDISP, WABS, &
                             THPOR, THLRET, THLMIN, BI, PSISAT, GRKSC, &
                             THLRAT, THFC, DELZW, ZBOTW, XDRAIN, DELZ, ISAND, &
                             IGRN, IGRD, IFILL, IZERO, LZF, NINF, IFIND, ITER, &
                             NEND, ISIMP, IGDR, &
                             IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call waterFlowNonInfiltrate(3, THLQCO, THICCO, TBARWC, FDUMMY, TDUMMY, & ! Formerly GRDRAN
                                BASFLW, TBASFL, RUNFC, TRUNFC, &
                                QFG, WLOSTC, FC, EVAPCG, RPCC, ZPONDC, &
                                DT, WEXCES, THLMAX, THTEST, THPOR, THLRET, THLMIN, &
                                BI, PSISAT, GRKSC, THFC, DELZW, XDRAIN, ISAND, &
                                IZERO, IGRN, IGRD, IGDR, &
                                IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call WATROF(THLQCO, THICCO, ZPONDC, TPONDC, OVRFLW, TOVRFL, &
                SUBFLW, TSUBFL, RUNFC, TRUNFC, FC, ZPLIMC, &
                XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBARWC, &
                DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN, &
                ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
    call waterUpdates(TBARC, THLQCO, THICCO, HCPCO, TPONDC, ZPONDC, & ! Formerly TMCALC
                      TSNOWC, ZSNOWC, ALBSC, RHOSC, HCPSC, TBASC, &
                      OVRFLW, TOVRFL, RUNFC, TRUNFC, HMFG, HTC, HTCS, &
                      WTRS, WTRG, FC, TBARWC, GZEROC, G12C, &
                      G23C, GGEO, TA, ZERO, TCTOPC, TCBOTC, GFLXC, &
                      ZPLIMC, THPOR, THLMIN, HCPS, DELZW, DELZZ, DELZ, &
                      ISAND, IWF, IG, ILG, IL1, IL2, JL, N)
    call checkWaterBudget(3, PCPR, EVPIC, RUNFC, WLOSTC, RAICAN, SNOCAN, & ! Formerly CHKWAT
                          RAC, SNC, ZPONDC, ZPOND, THLQCO, THICCO, &
                          THLIQC, THICEC, ZSNOWC, RHOSC, XSNOWC, SNO, &
                          ZERO, ZERO, FCS, FGS, FC, BAL, THPOR, THLMIN, &
                          DELZW, ISAND, IG, ILG, IL1, IL2, JL, N)
    !
  end if
  !
  !     * CALCULATIONS FOR BARE GROUND.
  !
  if (NLANDG > 0) then
    call soilWaterPhaseChg(TBARG, THLQGO, THICGO, HCPGO, TBARWG, HMFG, HTC, & ! Formerly TWCALC
                           FG, EVAPG, THPOR, THLMIN, HCPS, DELZW, DELZZ, &
                           ISAND, IG, ILG, IL1, IL2, JL)
    call snowSublimation(RHOSG, ZSNOWG, HCPSG, TSNOWG, EVAPG, QFN, QFG, & ! Formerly SNOVAP
                         HTCS, WLOSTG, TRUNFG, RUNFG, TOVRFL, OVRFLW, &
                         FG, RPCG, SPCG, RHOSNI, ZERO, ILG, IL1, IL2, JL)
    call pondedWaterFreeze(ZPONDG, TPONDG, ZSNOWG, TSNOWG, ALBSG, & ! Formerly TFREEZ
                           RHOSG, HCPSG, GZEROG, HMFG, HTCS, HTC, &
                           WTRS, WTRG, FG, QFREZG, ZERO, TA, TBARG, &
                           ISAND, IG, ILG, IL1, IL2, JL)
    call snowAddNew(ALBSG, TSNOWG, RHOSG, ZSNOWG, & ! Formerly SNOADD
                    HCPSG, HTCS, FG, SPCG, TSPCG, RHOSNI, ZERO, &
                    ILG, IL1, IL2, JL)
    if (NLANDI /= 0) then
      call iceSheetBalance(TBARG, TPONDG, ZPONDG, TSNOWG, RHOSG, ZSNOWG, & ! Formerly ICEBAL
                           HCPSG, ALBSG, HMFG, HTCS, HTC, WTRS, WTRG, GFLXG, &
                           RUNFG, TRUNFG, OVRFLW, TOVRFL, ZPLIMG, GGEO, &
                           FG, EVAPG, RPCG, TRPCG, GZEROG, G12G, G23G, &
                           HCPGO, QFREZG, ZERO, ZMAT, TMOVE, WMOVE, ZRMDR, &
                           TADD, ZMOVE, TBOT, DELZ, ISAND, ICONT, &
                           IWF, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    end if
    call waterFlowInfiltrate(4, THLQGO, THICGO, TBARWG, BASFLW, TBASFL, RUNFG, & ! Formerly GRINFL
                             TRUNFG, ZFAV, LZFAV, THLINV, QFG, WLOSTG, &
                             FG, EVAPG, RPCG, TRPCG, TPONDG, ZPONDG, &
                             DT, ZMAT, WMOVE, TMOVE, THLIQX, THICEX, TBARWX, &
                             DELZX, ZBOTX, FDT, TFDT, PSIF, THLINF, GRKINF, &
                             THLMAX, THTEST, ZRMDR, FDUMMY, TDUMMY, THLDUM, &
                             THIDUM, TDUMW, TRMDR, ZF, FMAX, TUSED, RDUMMY, &
                             ZERO, WEXCES, FDTBND, WADD, TADD, WADJ, TIMPND, &
                             DZF, DTFLOW, THLNLZ, THLQLZ, DZDISP, WDISP, WABS, &
                             THPOR, THLRET, THLMIN, BI, PSISAT, GRKSG, &
                             THLRAT, THFC, DELZW, ZBOTW, XDRAIN, DELZ, ISAND, &
                             IGRN, IGRD, IFILL, IZERO, LZF, NINF, IFIND, ITER, &
                             NEND, ISIMP, IGDR, &
                             IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call waterFlowNonInfiltrate(4, THLQGO, THICGO, TBARWG, FDUMMY, TDUMMY, & ! Formerly GRDRAN
                                BASFLW, TBASFL, RUNFG, TRUNFG, &
                                QFG, WLOSTG, FG, EVAPG, RPCG, ZPONDG, &
                                DT, WEXCES, THLMAX, THTEST, THPOR, THLRET, THLMIN, &
                                BI, PSISAT, GRKSG, THFC, DELZW, XDRAIN, ISAND, &
                                IZERO, IGRN, IGRD, IGDR, &
                                IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
    call WATROF(THLQGO, THICGO, ZPONDG, TPONDG, OVRFLW, TOVRFL, &
                SUBFLW, TSUBFL, RUNFG, TRUNFG, FG, ZPLIMG, &
                XSLOPE, XDRAINH, MANNING_N, DD, KSAT, TBARWG, &
                DELZW, THPOR, THLMIN, BI, DODRN, DOVER, DIDRN, &
                ISAND, IWF, IG, ILG, IL1, IL2, BULK_FC)
    call waterUpdates(TBARG, THLQGO, THICGO, HCPGO, TPONDG, ZPONDG, & ! Formerly TMCALC
                      TSNOWG, ZSNOWG, ALBSG, RHOSG, HCPSG, TBASG, &
                      OVRFLW, TOVRFL, RUNFG, TRUNFG, HMFG, HTC, HTCS, &
                      WTRS, WTRG, FG, TBARWG, GZEROG, G12G, &
                      G23G, GGEO, TA, ZERO, TCTOPG, TCBOTG, GFLXG, &
                      ZPLIMG, THPOR, THLMIN, HCPS, DELZW, DELZZ, DELZ, &
                      ISAND, IWF, IG, ILG, IL1, IL2, JL, N)
    call checkWaterBudget(4, PCPR, EVPIG, RUNFG, WLOSTG, RAICAN, SNOCAN, & ! Formerly CHKWAT
                          RAC, SNC, ZPONDG, ZPOND, THLQGO, THICGO, &
                          THLIQG, THICEG, ZSNOWG, RHOSG, XSNOWG, SNO, &
                          ZERO, ZERO, FCS, FGS, FG, BAL, THPOR, THLMIN, &
                          DELZW, ISAND, IG, ILG, IL1, IL2, JL, N)
    !
  end if
  !>
  !! After these calls have been done, average values of the main
  !! prognostic variables over the modelled area are determined by
  !! performing weighted averages over the four subareas, and checks
  !! are carried out to identify and remove vanishingly small values.
  !! First the bedrock temperature in the third soil layer, the total
  !! runoff and the runoff temperature are calculated. Then the total
  !! runoff and the overland flow, interflow and baseflow are
  !! converted from units of m to \f$kg m^{-2} s^{-1}\f$. The total surface water
  !! vapour flux over the modelled area is updated to account for the
  !! residual amounts of evaporative demand over the four subareas
  !! that could not be supplied by surface stores (WLSTCS, WLSTGS,
  !! WLOSTC and WLOSTG, variables that are defined internally in this
  !! subroutine), and the diagnostic variable QFLUX is evaluated.
  !!
  !! The temperature of the vegetation canopy TCAN and the amount of
  !! intercepted liquid water RCAN are calculated as weighted averages
  !! over the two canopy subareas. A flag is set to trigger a call to
  !! abort if TCAN is less than -100 C or greater than 100 C. If RCAN
  !! is vanishingly small, it is added to the overland flow and to the
  !! total runoff, and their respective temperatures are recalculated.
  !! The diagnostic arrays ROFC, ROVG, PCPG and HTCC are updated, and
  !! RCAN is set to zero. The amount of intercepted snow SNCAN is
  !! likewise calculated as a weighted average over the two canopy
  !! subareas. If SNCAN is vanishingly small, it is added to the
  !! overland flow and to the total runoff, and their respective
  !! temperatures are recalculated. The diagnostic arrays ROFC, ROVG,
  !! PCPG and HTCC are updated, and SNCAN is set to zero. If there is
  !! no canopy present, TCAN is set to zero.
  !!
  !! At the end of the 600 loop, the depth of ponded water ZPOND and
  !! its temperature TPOND over the modelled area are calculated as
  !! weighted averages over the four subareas. If ZPOND is vanishingly
  !! small, then as in the case of intercepted water, it is added to
  !! the overland flow and to the total runoff, and their respective
  !! temperatures are recalculated. The diagnostic array HTC is
  !! updated, ZPOND is set to zero and TPOND is set to 0 \f$^o\f$C.
  !!
  !
  !     * AVERAGE RUNOFF AND PROGNOSTIC VARIABLES OVER FOUR GRID CELL
  !     * SUBAREAS.
  !
  JPTBAD = 0
  KPTBAD = 0
  LPTBAD = 0
  do I = IL1,IL2
    TBASE (I) = FCS(I) * (TBASCS(I) + TFREZ) + &
                FGS(I) * (TBASGS(I) + TFREZ) + &
                FC (I) * (TBASC (I) + TFREZ) + &
                FG (I) * (TBASG (I) + TFREZ)
    RUNOFF(I) = FCS(I) * RUNFCS(I) + FGS(I) * RUNFGS(I) + &
                FC (I) * RUNFC (I) + FG (I) * RUNFG (I)
    if (RUNOFF(I) > 0.0) &
        TRUNOF(I) = (FCS(I) * RUNFCS(I) * TRNFCS(I) + &
        FGS(I) * RUNFGS(I) * TRNFGS(I) + &
        FC (I) * RUNFC (I) * TRUNFC(I) + &
        FG (I) * RUNFG (I) * TRUNFG(I)) / RUNOFF(I)
    RUNOFF(I) = RUNOFF(I) * RHOW / DELT
    OVRFLW(I) = OVRFLW(I) * RHOW / DELT
    BASFLW(I) = BASFLW(I) * RHOW / DELT
    EVAP  (I) = EVAP(I) - (FCS(I) * WLSTCS(I) + FGS(I) * WLSTGS(I) + &
                FC(I) * WLOSTC(I) + FG(I) * WLOSTG(I)) / DELT
    QFLUX(I) = - EVAP(I) / RHOAIR(I)
    if ((FC(I) + FCS(I)) > 0.) then
      TCAN(I) = (FCS(I) * TCANS(I) * CHCAPS(I) + FC(I) * TCANO(I) * &
                CHCAP(I)) / (FCS(I) * CHCAPS(I) + FC(I) * CHCAP(I))
      RCAN(I) = FCS(I) * RAICNS(I) + FC (I) * RAICAN(I)
      if (TCAN(I) < 173.16 .or. TCAN(I) > 373.16) JPTBAD = I
      if (RCAN(I) < 0.0) RCAN(I) = 0.0
      if (RCAN(I) < 1.0E-5 .and. RCAN(I) > 0.0) then
        TOVRFL(I) = (TOVRFL(I) * OVRFLW(I) + TCAN(I) * RCAN(I) / &
                    DELT) / (OVRFLW(I) + RCAN(I) / DELT)
        OVRFLW(I) = OVRFLW(I) + RCAN(I) / DELT
        TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + TCAN(I) * RCAN(I) / &
                    DELT) / (RUNOFF(I) + RCAN(I) / DELT)
        RUNOFF(I) = RUNOFF(I) + RCAN(I) / DELT
        ROFC(I) = ROFC(I) + RCAN(I) / DELT
        ROVG(I) = ROVG(I) + RCAN(I) / DELT
        PCPG(I) = PCPG(I) + RCAN(I) / DELT
        HTCC(I) = HTCC(I) - TCAN(I) * SPHW * RCAN(I) / DELT
        RCAN(I) = 0.0
      end if
      SNCAN  (I) = FCS(I) * SNOCNS(I) + FC (I) * SNOCAN(I)
      if (SNCAN(I) < 0.0) SNCAN(I) = 0.0
      if (SNCAN(I) < 1.0E-5 .and. SNCAN(I) > 0.0) then
        TOVRFL(I) = (TOVRFL(I) * OVRFLW(I) + TCAN(I) * SNCAN(I) / &
                    DELT) / (OVRFLW(I) + SNCAN(I) / DELT)
        OVRFLW(I) = OVRFLW(I) + SNCAN(I) / DELT
        TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + TCAN(I) * SNCAN(I) / &
                    DELT) / (RUNOFF(I) + SNCAN(I) / DELT)
        RUNOFF(I) = RUNOFF(I) + SNCAN(I) / DELT
        ROFC(I) = ROFC(I) + SNCAN(I) / DELT
        ROVG(I) = ROVG(I) + SNCAN(I) / DELT
        PCPG(I) = PCPG(I) + SNCAN(I) / DELT
        HTCC(I) = HTCC(I) - TCAN(I) * SPHICE * SNCAN(I) / DELT
        SNCAN(I) = 0.0
      end if
    else
      TCAN(I) = 0.0
    end if
    if (ZPNDCS(I) > 0. .or. ZPNDGS(I) > 0. .or. &
        ZPONDC(I) > 0. .or. ZPONDG(I) > 0.) then
      ZPOND(I) = (FCS(I) * ZPNDCS(I) + FGS(I) * ZPNDGS(I) + &
                 FC (I) * ZPONDC(I) + FG (I) * ZPONDG(I))
      TPOND(I) = (FCS(I) * (TPNDCS(I) + TFREZ) * ZPNDCS(I) + &
                 FGS(I) * (TPNDGS(I) + TFREZ) * ZPNDGS(I) + &
                 FC (I) * (TPONDC(I) + TFREZ) * ZPONDC(I) + &
                 FG (I) * (TPONDG(I) + TFREZ) * ZPONDG(I)) / &
                 ZPOND(I)
      if (ZPOND(I) < 0.0) ZPOND(I) = 0.0
      if (ZPOND(I) < 1.0E-8 .and. ZPOND(I) > 0.0) then
        TOVRFL(I) = (TOVRFL(I) * OVRFLW(I) + TPOND(I) * RHOW * &
                    ZPOND(I) / DELT) / (OVRFLW(I) + RHOW * ZPOND(I) / DELT)
        OVRFLW(I) = OVRFLW(I) + RHOW * ZPOND(I) / DELT
        TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + TPOND(I) * RHOW * &
                    ZPOND(I) / DELT) / (RUNOFF(I) + RHOW * ZPOND(I) / DELT)
        RUNOFF(I) = RUNOFF(I) + RHOW * ZPOND(I) / DELT
        HTC(I,1) = HTC(I,1) - TPOND(I) * HCPW * ZPOND(I) / DELT
        TPOND(I) = TFREZ
        ZPOND(I) = 0.0
      end if
    else
      ZPOND(I) = 0.0
      TPOND(I) = TFREZ
    end if
  end do ! loop 600
  !
  !>
  !! In the 650 loop, values of the snow prognostic variables are
  !! calculated as weighted averages over the four subareas. The
  !! weightings for the subareas include the four internally-defined
  !! waterBudgetDriver variables XSNOCS, XSNOGS, XSNOWC and XSNOWG, which are set
  !! in subroutine checkWaterBudget to 1 if the subarea snow depth is greater
  !! than zero, and to zero otherwise. If the snow depth over the CS
  !! and GS subareas is greater than zero (meaning that there was a
  !! pre-existing snow cover at the beginning of the time step), the
  !! average snow albedo ALBSNO is preferentially set to the average
  !! over these two subareas. Otherwise ALBSNO is set to the average
  !! over the C and G subareas (where snow has just been added in the
  !! current time step). The snow temperature TSNOW and density RHOSNO
  !! are set to weighted averages over the four subareas, using the
  !! internally-defined subarea volumetric heat capacities
  !! HCPSCS/GS/C/G and RHOSCS/GS/C/G. Finally the snow depth ZSNOW is
  !! calculated from the subarea depths; the liquid water content of
  !! the snow pack WSNOW is obtained as a weighted average over the CS
  !! and GS subareas (assuming that freshly fallen snow does not yet
  !! contain liquid water); and the snow mass is determined from ZSNOW
  !! and RHOSNO. Upper and lower limits apply to the snow pack.  If the
  !! snow depth exceeds 10 m, the excess snow and its associated liquid
  !! water content are added to the overland flow and the total runoff,
  !! the respective temperatures of the latter are recalculated and the
  !! diagnostic arrays ROFN, PCPG and HTCS are updated. As in the case
  !! of intercepted and ponded water, if
  !! the snow mass is vanishingly small it and its liquid water
  !! content are added to the overland flow and to the total runoff,
  !! and their respective temperatures are recalculated, the
  !! diagnostic arrays ROFN, PCPG and HTCS are updated, and TSNOW,
  !! RHOSNO, SNO and WSNOW are set to zero. Flags are set to trigger
  !! calls to abort if TSNOW is less than 0 K or greater than 0.001 C.
  !! Finally, the three abort flags set thus far are checked, and
  !! calls to abort are performed if they are greater than zero.
  !!
  do I = IL1,IL2
    if (ZSNOCS(I) > 0. .or. ZSNOGS(I) > 0. .or. &
        ZSNOWC(I) > 0. .or. ZSNOWG(I) > 0.) then
      if (ZSNOCS(I) > 0. .or. ZSNOGS(I) > 0.) then
        ALBSNO(I) = (FCS(I) * ALBSCS(I) * XSNOCS(I) + &
                    FGS(I) * ALBSGS(I) * XSNOGS(I)) / &
                    (FCS(I) * XSNOCS(I) + FGS(I) * XSNOGS(I))
      else
        ALBSNO(I) = (FC (I) * ALBSC(I) * XSNOWC(I) + &
                    FG (I) * ALBSG(I) * XSNOWG(I)) / &
                    (FC (I) * XSNOWC(I) + FG (I) * XSNOWG(I))
      end if
      TSNOW(I) = (FCS(I) * (TSNOCS(I) + TFREZ) * HCPSCS(I) * &
                 ZSNOCS(I) * XSNOCS(I) + &
                 FGS(I) * (TSNOGS(I) + TFREZ) * HCPSGS(I) * &
                 ZSNOGS(I) * XSNOGS(I) + &
                 FC (I) * (TSNOWC(I) + TFREZ) * HCPSC(I) * &
                 ZSNOWC(I) * XSNOWC(I) + &
                 FG (I) * (TSNOWG(I) + TFREZ) * HCPSG(I) * &
                 ZSNOWG(I) * XSNOWG(I)) / &
                 (FCS(I) * HCPSCS(I) * ZSNOCS(I) * XSNOCS(I) + &
                 FGS(I) * HCPSGS(I) * ZSNOGS(I) * XSNOGS(I) + &
                 FC (I) * HCPSC(I) * ZSNOWC(I) * XSNOWC(I) + &
                 FG (I) * HCPSG(I) * ZSNOWG(I) * XSNOWG(I))
      RHOSNO(I) = (FCS(I) * RHOSCS(I) * ZSNOCS(I) * XSNOCS(I) + &
                  FGS(I) * RHOSGS(I) * ZSNOGS(I) * XSNOGS(I) + &
                  FC (I) * RHOSC(I) * ZSNOWC(I) * XSNOWC(I) + &
                  FG (I) * RHOSG(I) * ZSNOWG(I) * XSNOWG(I)) / &
                  (FCS(I) * ZSNOCS(I) * XSNOCS(I) + &
                  FGS(I) * ZSNOGS(I) * XSNOGS(I) + &
                  FC (I) * ZSNOWC(I) * XSNOWC(I) + &
                  FG (I) * ZSNOWG(I) * XSNOWG(I))
      ZSNOW(I) = FCS(I) * ZSNOCS(I) + FGS(I) * ZSNOGS(I) + &
                 FC (I) * ZSNOWC(I) + FG (I) * ZSNOWG(I)
      WSNOW(I) = FCS(I) * WSNOCS(I) + FGS(I) * WSNOGS(I)
      SNO(I) = ZSNOW(I) * RHOSNO(I)
      if (SNO(I) < 0.0) SNO(I) = 0.0
      !
      !           * LIMIT SNOW MASS TO A MAXIMUM OF 10 METRES TO AVOID
      !           * SNOW PILING UP AT EDGES OF GLACIERS AT HIGH ELEVATIONS.
      !           * THIS IS TARGETTED AS OVERLAND RUNOFF.
      !
      if (ZSNOW(I) > 10.0) then
        SNOROF = (ZSNOW(I) - 10.0) * RHOSNO(I)
        WSNROF = WSNOW(I) * SNOROF / SNO(I)
        TOVRFL(I) = (TOVRFL(I) * OVRFLW(I) + TSNOW(I) * (SNOROF + &
                    WSNROF) / DELT) / (OVRFLW(I) + (SNOROF + WSNROF) / &
                    DELT)
        OVRFLW(I) = OVRFLW(I) + (SNOROF + WSNROF) / DELT
        TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + TSNOW(I) * (SNOROF + &
                    WSNROF) / DELT) / (RUNOFF(I) + (SNOROF + WSNROF) / &
                    DELT)
        RUNOFF(I) = RUNOFF(I) + (SNOROF + WSNROF) / DELT
        ROFN(I) = ROFN(I) + (SNOROF + WSNROF) / DELT
        PCPG(I) = PCPG(I) + (SNOROF + WSNROF) / DELT
        HTCS(I) = HTCS(I) - TSNOW(I) * (SPHICE * SNOROF + SPHW * &
                  WSNROF) / DELT
        SNO(I) = SNO(I) - SNOROF
        WSNOW(I) = WSNOW(I) - WSNROF
        ZSNOW(I) = 10.0
      end if
      if (SNO(I) < 1.0E-2 .and. SNO(I) > 0.0) then
        TOVRFL(I) = (TOVRFL(I) * OVRFLW(I) + TSNOW(I) * (SNO(I) + &
                    WSNOW(I)) / DELT) / (OVRFLW(I) + (SNO(I) + WSNOW(I)) / &
                    DELT)
        OVRFLW(I) = OVRFLW(I) + (SNO(I) + WSNOW(I)) / DELT
        TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + TSNOW(I) * (SNO(I) + &
                    WSNOW(I)) / DELT) / (RUNOFF(I) + (SNO(I) + WSNOW(I)) / &
                    DELT)
        RUNOFF(I) = RUNOFF(I) + (SNO(I) + WSNOW(I)) / DELT
        ROFN(I) = ROFN(I) + (SNO(I) + WSNOW(I)) / DELT
        PCPG(I) = PCPG(I) + (SNO(I) + WSNOW(I)) / DELT
        HTCS(I) = HTCS(I) - TSNOW(I) * (SPHICE * SNO(I) + SPHW * &
                  WSNOW(I)) / DELT
        TSNOW(I) = 0.0
        RHOSNO(I) = 0.0
        SNO(I) = 0.0
        WSNOW(I) = 0.0
      end if
    else
      TSNOW(I) = 0.0
      RHOSNO(I) = 0.0
      SNO(I) = 0.0
      WSNOW(I) = 0.0
      ZSNOW(I) = 0.0
    end if
    !
    if (TSNOW(I) < 0.0) KPTBAD = I
    if ((TSNOW(I) - TFREZ) > 1.0E-3) LPTBAD = I
  end do ! loop 650
  !
  if (JPTBAD /= 0) then
    write(6,6625) JPTBAD,JL,TCAN(JPTBAD)
6625 format('0AT (I,J) = (',I3,',',I3,'),TCAN = ',F10.5)
    call errorHandler('waterBudgetDriver2', - 2)
  end if
  !
  if (KPTBAD /= 0) then
    write(6,6626) KPTBAD,JL,TSNOW(KPTBAD)
6626 format('0AT (I,J) = (',I3,',',I3,'),TSNOW = ',F10.5)
    call errorHandler('waterBudgetDriver2', - 3)
  end if
  !
  if (LPTBAD /= 0) then
    write(6,6626) LPTBAD,JL,TSNOW(LPTBAD)
    call errorHandler('waterBudgetDriver2', - 4)
  end if
  !
  !>
  !! In the 700 loop, the temperature of each soil layer is calculated
  !! as a weighted average over the four subareas. In the case of the
  !! third soil layer., if the standard three-layer configuration is
  !! being modelled (with a very thick third soil layer of 3.75 m),
  !! the subarea layer temperatures TBARCS/GS/C/G and the layer heat
  !! capacities HCPCS/GS/C/G apply to the permeable depth DELZW of the
  !! layer, and the bedrock temperature TBASE and the rock heat
  !! capacity HCPSND to the remainder, DELZ-DELZW. The averaging is
  !! carried out accordingly. In all other soil layers, the layer
  !! temperature applies to the whole thickness, whose heat capacity
  !! is a weighted average of HCPCS/GS/C/G over DELZW and HCPSND over
  !! DELZ-DELZW. The volumetric liquid water content THLIQ, the
  !! volumetric frozen water content THICE, and the heat flux at the
  !! soil layer interfaces GFLUX are calculated as simple weighted
  !! averages over the subareas. A flag is set to trigger a call to
  !! abort if the soil layer temperature is less than -100 C or
  !! greater than 100 C, and after the end of the loop, a call to
  !! abort is performed if the flag is greater than zero.
  !!
  IPTBAD = 0
  do J = 1,IG ! loop 700
    do I = IL1,IL2
      if (IG == 3 .and. J == IG .and. ISAND(I,1) > - 4) then
        TBAR(I,J) = ((FCS(I) * (TBARCS(I,J) + TFREZ) * HCPCS(I,J) + &
                    FGS(I) * (TBARGS(I,J) + TFREZ) * HCPGS(I,J) + &
                    FC (I) * (TBARC (I,J) + TFREZ) * HCPCO(I,J) + &
                    FG (I) * (TBARG (I,J) + TFREZ) * HCPGO(I,J)) * &
                    DELZW(I,J) + TBASE(I) * HCPSND * &
                    (DELZ(J) - DELZW(I,J))) / &
                    ((FCS(I) * HCPCS(I,J) + FGS(I) * HCPGS(I,J) + &
                    FC (I) * HCPCO(I,J) + FG (I) * HCPGO(I,J)) * &
                    DELZW(I,J) + HCPSND * (DELZ(J) - DELZW(I,J)))
      else
        TBAR(I,J) = (FCS(I) * (TBARCS(I,J) + TFREZ) * (DELZW(I,J) * &
                    HCPCS(I,J) + (DELZ(J) - DELZW(I,J)) * HCPSND) + &
                    FGS(I) * (TBARGS(I,J) + TFREZ) * (DELZW(I,J) * &
                    HCPGS(I,J) + (DELZ(J) - DELZW(I,J)) * HCPSND) + &
                    FC (I) * (TBARC (I,J) + TFREZ) * (DELZW(I,J) * &
                    HCPCO(I,J) + (DELZ(J) - DELZW(I,J)) * HCPSND) + &
                    FG (I) * (TBARG (I,J) + TFREZ) * (DELZW(I,J) * &
                    HCPGO(I,J) + (DELZ(J) - DELZW(I,J)) * HCPSND)) / &
                    (FCS(I) * (DELZW(I,J) * HCPCS(I,J) + &
                    (DELZ(J) - DELZW(I,J)) * HCPSND) + &
                    FGS(I) * (DELZW(I,J) * HCPGS(I,J) + &
                    (DELZ(J) - DELZW(I,J)) * HCPSND) + &
                    FC (I) * (DELZW(I,J) * HCPCO(I,J) + &
                    (DELZ(J) - DELZW(I,J)) * HCPSND) + &
                    FG (I) * (DELZW(I,J) * HCPGO(I,J) + &
                    (DELZ(J) - DELZW(I,J)) * HCPSND))
      end if
      THLIQ(I,J) = FCS(I) * THLQCS(I,J) + FGS(I) * THLQGS(I,J) + &
                   FC (I) * THLQCO(I,J) + FG (I) * THLQGO(I,J)
      THICE(I,J) = FCS(I) * THICCS(I,J) + FGS(I) * THICGS(I,J) + &
                   FC (I) * THICCO(I,J) + FG (I) * THICGO(I,J)
      GFLUX(I,J) = FCS(I) * GFLXCS(I,J) + FGS(I) * GFLXGS(I,J) + &
                   FC (I) * GFLXC (I,J) + FG (I) * GFLXG (I,J)
      SUBFLW(I,J) = SUBFLW(I,J) * RHOW / DELT
      !     ipy test
      !          IF (THLIQ(I,J)>THFC(I,J))                               THEN
      !              BASFLW(I)=BASFLW(I)+(THLIQ(I,J)-THFC(I,J))*DELZW(I,J)*
      !     1            RHOW/DELT
      !              RUNOFF(I)=RUNOFF(I)+(THLIQ(I,J)-THFC(I,J))*DELZW(I,J)*
      !     1            RHOW/DELT
      !              HTC(I,J)=HTC(I,J)-TBAR(I,J)*(THLIQ(I,J)-THFC(I,J))*
      !     1            HCPW*DELZW(I,J)/DELT
      !              THLIQ(I,J)=THFC(I,J)
      !          END IF
      if (TBAR(I,1) < 173.16 .or. TBAR(I,1) > 373.16) IPTBAD = I
    end do
  end do ! loop 700
  !
  if (IPTBAD /= 0) then
    write(6,6600) IPTBAD,JL,TBAR(IPTBAD,1)
6600 format('0AT (I,J) = (',I3,',',I3,'),TBAR(1) = ',F10.5)
    call errorHandler('waterBudgetDriver2', - 1)
  end if
  !
  !> Finally, subroutine classGrowthIndex is called to update the vegetation
  !! growth index.
  !!
  call classGrowthIndex(GROWTH, TBAR, TA, FC, FCS, ILG, IG, IL1, IL2, JL) ! Formerly CGROW
  !
  return
end subroutine waterBudgetDriver
