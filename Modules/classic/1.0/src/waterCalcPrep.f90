!> \file
!! Initializes subarea variables for surface water budget
!! calculations, and performs preliminary calculations for diagnostic
!! variables.
!! @author D. Verseghy, M. Lazare, Y. Delage, R. Harvey, R. Soulis, P. Bartlett
!
subroutine waterCalcPrep (THLQCO, THLQGO, THLQCS, THLQGS, THICCO, THICGO, & ! Formerly WPREP
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
                          TBARC, TBARG, TBARCS, TBARGS, TBASE, TSURX, &
                          FSVF, FSVFS, RAICAN, SNOCAN, RAICNS, SNOCNS, &
                          EVAPC, EVAPCG, EVAPG, EVAPCS, EVPCSG, EVAPGS, &
                          RPCP, TRPCP, SPCP, TSPCP, RHOSNI, ZPOND, &
                          ZSNOW, ALBSNO, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                          THPOR, HCPS, GRKSAT, ISAND, DELZW, DELZ, &
                          ILG, IL1, IL2, JL, IG, IGP1, &
                          NLANDCS, NLANDGS, NLANDC, NLANDG, RADD, SADD, &
                          BI, PSISAT, DD, XSLOPE, BULK_FC)
  !
  !     * AUG 25/11 - D.VERSEGHY. REFINE CALCULATION OF TEMPERATURE OF
  !     *                         LUMPED PRECIPITATION.
  !     * NOV 24/09 - D.VERSEGHY. RESTORE EVAPOTRANSPIRATION WHEN
  !     *                         PRECIPITATION IS OCCURRING.
  !     * MAR 27/08 - D.VERSEGHY. MOVE MODIFICATION OF GRKSAT IN PRESENCE
  !     *                         OF ICE TO waterFlowInfiltrate AND waterFlowNonInfiltrate.
  !     * FEB 19/07 - D.VERSEGHY. MODIFICATIONS TO REFLECT SHIFT OF CANOPY
  !     *                         WATER DEPOSITION CALCULATIONS TO energBalVegSolve,
  !     *                         AND SUPPRESSION OF ALL EVAPOTRANSPIRATION
  !     *                         WHEN PRECIPITATION IS OCCURRING.
  !     * MAR 23/06 - D.VERSEGHY. MODIFY CALCULATIONS OF SNOW THERMAL
  !     *                         PROPERTIES TO ACCOUNT FOR WATER CONTENT.
  !     * MAR 21/06 - P.BARTLETT. INITIALIZE ADDITIONAL VARIABLES TO ZERO.
  !     * DEC 07/05 - D.VERSEGHY. ADD INITIALIZATION OF TBASE SUBAREAS.
  !     * OCT 05/05 - D.VERSEGHY. MODIFY DELZZ CALCULATION FOR IG>3.
  !     * APR 15/05 - D.VERSEGHY. SUBLIMATION OF INTERCEPTED SNOW TAKES
  !     *                         PLACE BEFORE EVAPORATION OF INTERCEPTED
  !     *                         RAIN.
  !     * MAR 30/05 - D.VERSEGHY/R.SOULIS. ADD RUNOFF TEMPERATURE
  !     *                         INITIALIZATIONS AND MODIFICATION OF
  !     *                         GRKSAT FOR TEMPERATURE AND PRESENCE
  !     *                         OF ICE.
  !     * SEP 13/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 21/04 - Y.DELAGE, R.HARVEY, D.VERSEGHY. NEW LOWER LIMIT
  !     *                         ON RADD AND SADD.
  !     * SEP 26/02 - D.VERSEGHY. MODIFICATIONS ASSOCIATED WITH BUGFIX
  !     *                         IN waterUnderCanopy.
  !     * AUG 06/02 - D.VERSEGHY. SHORTENED CLASS3 COMMON BLOCK.
  !     * JUN 18/02 - D.VERSEGHY. MOVE PARTITIONING OF PRECIPITATION
  !     *                         BETWEEN RAINFALL AND SNOWFALL INTO
  !     *                         "atmosphericVarsCalc"; TIDY UP SUBROUTINE CALL
  !     *                         CHANGE RHOSNI FROM CONSTANT TO
  !     *                         VARIABLE.
  !     * OCT 04/01 - M.LAZARE.   NEW DIAGNOSTIC FIELD "ROVG".
  !     * NOV 09/00 - D.VERSEGHY. MOVE DIAGNOSTIC CALCULATIONS FROM
  !     *                         waterUnderCanopy INTO THIS ROUTINE.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         CHANGES RELATED TO VARIABLE SOIL DEPTH
  !     *                         (MOISTURE HOLDING CAPACITY) AND DEPTH-
  !     *                         VARYING SOIL PROPERTIES.
  !     * JAN 02/95 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS; INTRODUCE CALCULATION OF
  !     *                         OVERLAND FLOW.
  !     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         RATIONALIZE USE OF "WLOST":
  !     *                         COMPLETION OF WATER BUDGET DIAGNOSTICS.
  !     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
  !     *                         ORGANIC MATTER CONTENT.
  !     * DEC 16/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         INITIALIZE TWO NEW DIAGNOSTIC FIELDS.
  !     * AUG 20/93 - D.VERSEGHY. CLASS - VERSION 2.2.
  !     *                         REVISED CALCULATION OF CANOPY
  !     *                         SUBLIMATION RATE.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. NEW DIAGNOSTIC FIELDS.
  !     * APR 15/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. PREPARATION AND INITIALIZATION FOR
  !     *                         LAND SURFACE WATER BUDGET CALCULATIONS.
  use classicParams, only : DELT, TFREZ, HCPW, HCPICE, RHOW, RHOICE

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL, IG, IGP1, NLANDCS, NLANDGS, NLANDC, NLANDG
  integer             :: I, J
  !
  !     * OUTPUT ARRAYS.
  !
  !
  !(Suffix CS = vegetation over snow cover; GS = bare snow cover; C
  ! or CO = vegetation over ground; G or GO = bare ground.)
  !
  real, intent(out) :: THLQCO(ILG,IG) !< Subarea volumetric liquid water content of soil layers \f$[m^3 m^{-3}]\f$
  real, intent(out) :: THLQGO(ILG,IG) !< Subarea volumetric liquid water content of soil layers \f$[m^3 m^{-3}]\f$
  real, intent(out) :: THLQCS(ILG,IG) !< Subarea volumetric liquid water content of soil layers \f$[m^3 m^{-3}]\f$
  real, intent(out) :: THLQGS(ILG,IG) !< Subarea volumetric liquid water content of soil layers \f$[m^3 m^{-3}]\f$
  !
  real, intent(out) :: THICCO(ILG,IG) !< Subarea volumetric frozen water content of soil layers
  !< \f$[m^3 m^{-3}] (\theta_i)\f$
  real, intent(out) :: THICGO(ILG,IG) !< Subarea volumetric frozen water content of soil layers
  !< \f$[m^3 m^{-3}] (\theta_i)\f$
  real, intent(out) :: THICCS(ILG,IG) !< Subarea volumetric frozen water content of soil layers
  !< \f$[m^3 m^{-3}] (\theta_i)\f$
  real, intent(out) :: THICGS(ILG,IG) !< Subarea volumetric frozen water content of soil layers
  !< \f$[m^3 m^{-3}] (\theta_i)\f$
  !
  real, intent(out) :: HCPCO (ILG,IG) !< Subarea heat capacity of soil layers \f$[J m^{-3} K^{-1}]\f$
  real, intent(out) :: HCPGO (ILG,IG) !< Subarea heat capacity of soil layers \f$[J m^{-3} K^{-1}]\f$
  real, intent(out) :: HCPCS (ILG,IG) !< Subarea heat capacity of soil layers \f$[J m^{-3} K^{-1}]\f$
  real, intent(out) :: HCPGS (ILG,IG) !< Subarea heat capacity of soil layers \f$[J m^{-3} K^{-1}]\f$
  !
  real, intent(out) :: GRKSC (ILG,IG) !< Subarea saturated hydraulic conductivity \f$[m s^{-1}] (K_{sat})\f$
  real, intent(out) :: GRKSG (ILG,IG) !< Subarea saturated hydraulic conductivity \f$[m s^{-1}] (K_{sat})\f$
  real, intent(out) :: GRKSCS(ILG,IG) !< Subarea saturated hydraulic conductivity \f$[m s^{-1}] (K_{sat})\f$
  real, intent(out) :: GRKSGS(ILG,IG) !< Subarea saturated hydraulic conductivity \f$[m s^{-1}] (K_{sat})\f$
  !
  real, intent(out) :: GFLXC (ILG,IG) !< Subarea heat flux between soil layers \f$[W m^{-2}]\f$
  real, intent(out) :: GFLXG (ILG,IG) !< Subarea heat flux between soil layers \f$[W m^{-2}]\f$
  real, intent(out) :: GFLXCS(ILG,IG) !< Subarea heat flux between soil layers \f$[W m^{-2}]\f$
  real, intent(out) :: GFLXGS(ILG,IG) !< Subarea heat flux between soil layers \f$[W m^{-2}]\f$
  !
  real, intent(out) :: THLDUM(ILG,IG) !< Internal waterBudgetDriver dummy variable for soil frozen water \f$[m^3 m^{-3}]\f$
  real, intent(out) :: THIDUM(ILG,IG) !< Internal waterBudgetDriver dummy variable for soil frozen water \f$[m^3 m^{-3}]\f$
  real, intent(out) :: QFC   (ILG,IG) !< Water removed from soil layers by transpiration \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: HMFG  (ILG,IG) !< Energy associated with phase change of water in soil layers \f$[W m^{-2}]\f$
  !
  real, intent(out) :: THLIQX(ILG,IGP1) !< Internal waterBudgetDriver work array for soil frozen water \f$[m^3 m^{-3}]\f$
  real, intent(out) :: THICEX(ILG,IGP1) !< Internal waterBudgetDriver work array for soil frozen water \f$[m^3 m^{-3}]\f$
  !
  real, intent(out) :: SPCC  (ILG) !< Subarea snowfall rate \f$[m s^{-1}]\f$
  real, intent(out) :: SPCG  (ILG) !< Subarea snowfall rate \f$[m s^{-1}]\f$
  real, intent(out) :: SPCCS (ILG) !< Subarea snowfall rate \f$[m s^{-1}]\f$
  real, intent(out) :: SPCGS (ILG) !< Subarea snowfall rate \f$[m s^{-1}]\f$
  !
  real, intent(out) :: TSPCC (ILG) !< Subarea snowfall temperature [K/C]
  real, intent(out) :: TSPCG (ILG) !< Subarea snowfall temperature [K/C]
  real, intent(out) :: TSPCCS(ILG) !< Subarea snowfall temperature [K/C]
  real, intent(out) :: TSPCGS(ILG) !< Subarea snowfall temperature [K/C]
  !
  real, intent(out) :: RPCC  (ILG) !< Subarea rainfall rate \f$[m s^{-1}]\f$
  real, intent(out) :: RPCG  (ILG) !< Subarea rainfall rate \f$[m s^{-1}]\f$
  real, intent(out) :: RPCCS (ILG) !< Subarea rainfall rate \f$[m s^{-1}]\f$
  real, intent(out) :: RPCGS (ILG) !< Subarea rainfall rate \f$[m s^{-1}]\f$
  !
  real, intent(out) :: TRPCC (ILG) !< Subarea rainfall temperature [K/C]
  real, intent(out) :: TRPCG (ILG) !< Subarea rainfall temperature [K/C]
  real, intent(out) :: TRPCCS(ILG) !< Subarea rainfall temperature [K/C]
  real, intent(out) :: TRPCGS(ILG) !< Subarea rainfall temperature [K/C]
  !
  real, intent(out) :: EVPIC (ILG) !< Subarea evapotranspiration rate going into waterBudgetDriver \f$[m s^{-1}]\f$
  real, intent(out) :: EVPIG (ILG) !< Subarea evapotranspiration rate going into waterBudgetDriver \f$[m s^{-1}]\f$
  real, intent(out) :: EVPICS(ILG) !< Subarea evapotranspiration rate going into waterBudgetDriver \f$[m s^{-1}]\f$
  real, intent(out) :: EVPIGS(ILG) !< Subarea evapotranspiration rate going into waterBudgetDriver \f$[m s^{-1}]\f$
  !
  real, intent(out) :: ZPONDC(ILG) !< Subarea depth of surface ponded water [m]
  real, intent(out) :: ZPONDG(ILG) !< Subarea depth of surface ponded water [m]
  real, intent(out) :: ZPNDCS(ILG) !< Subarea depth of surface ponded water [m]
  real, intent(out) :: ZPNDGS(ILG) !< Subarea depth of surface ponded water [m]
  !
  real, intent(out) :: XSNOWC(ILG) !< Subarea fractional snow coverage [ ]
  real, intent(out) :: XSNOWG(ILG) !< Subarea fractional snow coverage [ ]
  real, intent(out) :: XSNOCS(ILG) !< Subarea fractional snow coverage [ ]
  real, intent(out) :: XSNOGS(ILG) !< Subarea fractional snow coverage [ ]
  !
  real, intent(out) :: ZSNOWC(ILG) !< Subarea depth of snow pack
  real, intent(out) :: ZSNOWG(ILG) !< Subarea depth of snow pack
  real, intent(out) :: ZSNOCS(ILG) !< Subarea depth of snow pack
  real, intent(out) :: ZSNOGS(ILG) !< Subarea depth of snow pack
  !
  real, intent(out) :: ALBSC (ILG) !< Subarea snow albedo [ ]
  real, intent(out) :: ALBSG (ILG) !< Subarea snow albedo [ ]
  real, intent(out) :: ALBSCS(ILG) !< Subarea snow albedo [ ]
  real, intent(out) :: ALBSGS(ILG) !< Subarea snow albedo [ ]
  !
  real, intent(out) :: RHOSC (ILG) !< Subarea snow density \f$[kg m^{-3}]\f$
  real, intent(out) :: RHOSG (ILG) !< Subarea snow density \f$[kg m^{-3}]\f$
  !
  real, intent(out) :: HCPSC (ILG) !< Subarea heat capacity of snow \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(out) :: HCPSG (ILG) !< Subarea heat capacity of snow \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(out) :: HCPSCS(ILG) !< Subarea heat capacity of snow \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(out) :: HCPSGS(ILG) !< Subarea heat capacity of snow \f$[J m^{-3} K^{-1}] (C_s)\f$
  !
  real, intent(out) :: RUNFC (ILG) !< Subarea total runoff [m]
  real, intent(out) :: RUNFG (ILG) !< Subarea total runoff [m]
  real, intent(out) :: RUNFCS(ILG) !< Subarea total runoff [m]
  real, intent(out) :: RUNFGS(ILG) !< Subarea total runoff [m]
  !
  real, intent(out) :: TRUNFC(ILG) !< Subarea total runoff temperature [K]
  real, intent(out) :: TRUNFG(ILG) !< Subarea total runoff temperature [K]
  real, intent(out) :: TRNFCS(ILG) !< Subarea total runoff temperature [K]
  real, intent(out) :: TRNFGS(ILG) !< Subarea total runoff temperature [K]
  !
  real, intent(out) :: TBASC (ILG) !< Subarea temperature of bedrock in third soil layer [C]
  real, intent(out) :: TBASG (ILG) !< Subarea temperature of bedrock in third soil layer [C]
  real, intent(out) :: TBASCS(ILG) !< Subarea temperature of bedrock in third soil layer [C]
  real, intent(out) :: TBASGS(ILG) !< Subarea temperature of bedrock in third soil layer [C]
  !
  !
  real, intent(inout) :: SUBLC (ILG) !< Subarea sublimation rate from vegetation \f$[m s^{-1}]\f$
  real, intent(inout) :: SUBLCS(ILG) !< Subarea sublimation rate from vegetation \f$[m s^{-1}]\f$

  real, intent(out) :: WLOSTC(ILG) !< Subarea residual water not met by surface stores \f$[kg m^{-2}]\f$
  real, intent(out) :: WLOSTG(ILG) !< Subarea residual water not met by surface stores \f$[kg m^{-2}]\f$
  real, intent(out) :: WLSTCS(ILG) !< Subarea residual water not met by surface stores \f$[kg m^{-2}]\f$
  real, intent(out) :: WLSTGS(ILG) !< Subarea residual water not met by surface stores \f$[kg m^{-2}]\f$
  !
  real, intent(out) :: RAC   (ILG) !< Subarea liquid water on canopy going into waterBudgetDriver \f$[kg m^{-2}]\f$
  real, intent(out) :: RACS  (ILG) !< Subarea liquid water on canopy going into waterBudgetDriver \f$[kg m^{-2}]\f$
  real, intent(out) :: SNC   (ILG) !< Subarea frozen water on canopy going into waterBudgetDriver \f$[kg m^{-2}]\f$
  real, intent(out) :: SNCS  (ILG) !< Subarea frozen water on canopy going into waterBudgetDriver \f$[kg m^{-2}]\f$
  real, intent(out) :: TSNOWC(ILG) !< Subarea snowpack temperature [K]
  real, intent(out) :: TSNOWG(ILG) !< Subarea snowpack temperature [K]

  real, intent(out) :: OVRFLW(ILG)  !< Overland flow from top of soil column [m]
  real, intent(out) :: SUBFLW(ILG,IG)  !< Interflow from sides of soil column [m]
  real, intent(out) :: BASFLW(ILG)  !< Base flow from bottom of soil column [m]
  real, intent(out) :: TOVRFL(ILG)  !< Temperature of overland flow from top of soil column [K]
  real, intent(out) :: TSUBFL(ILG,IG)  !< Temperature of interflow from sides of soil column [K]
  real, intent(out) :: TBASFL(ILG)  !< Temperature of base flow from bottom of soil column [K]
  real, intent(inout) :: PCFC  (ILG)  !< Frozen precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCLC  (ILG)  !< Liquid precipitation intercepted by vegetation \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCPN  (ILG)  !< Precipitation incident on snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCPG  (ILG)  !< Precipitation incident on ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFCF  (ILG)  !< Sublimation from frozen water on vegetation \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFCL  (ILG)  !< Evaporation from liquid water on vegetation \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFN   (ILG)  !< Sublimation from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFG   (ILG)  !< Evaporation from ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: ROVG  (ILG)  !< Liquid/frozen water runoff from vegetation to ground surface \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: ROFC  (ILG)  !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: ROFN  (ILG)  !< Liquid water runoff from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: TRUNOF(ILG)  !< Temperature of total runoff [K]
  real, intent(out) :: DT    (ILG)  !< Time stepping variable used in waterFlowInfiltrate/waterFlowNonInfiltrate [s]
  real, intent(out) :: RDUMMY(ILG)  !< Dummy variable
  real, intent(out) :: ZERO  (ILG)  !< Zero vector used in several subroutines [ ]
  !
  integer, intent(out) :: IZERO (ILG) !< Zero integer :: flag used in waterFlowInfiltrate
  !
  !     * INPUT ARRAYS.
  !
  !(Suffix CS = vegetation over snow cover; GS = bare snow cover; C
  ! or CO = vegetation over ground; G or GO = bare ground.)
  !
  real, intent(in) :: FC    (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FG    (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FCS   (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FGS   (ILG) !< Subarea fractional coverage of modelled area [ ]
  !
  real, intent(in) :: FSVF  (ILG)  !< Sky view factor of ground under vegetation canopy [ ]
  real, intent(in) :: FSVFS (ILG)  !< Sky view factor of snow under vegetation canopy [ ]
  real, intent(in) :: RAICAN(ILG)  !< Intercepted liquid water stored on canopy over ground \f$[kg m^{-2}]\f$
  real, intent(in) :: SNOCAN(ILG)  !< Intercepted frozen water stored on canopy over ground \f$[kg m^{-2}]\f$
  real, intent(in) :: RAICNS(ILG)  !< Intercepted liquid water stored on canopy over snow \f$[kg m^{-2}]\f$
  real, intent(in) :: SNOCNS(ILG)  !< Intercepted frozen water stored on canopy over snow \f$[kg m^{-2}]\f$
  real, intent(inout) :: EVAPC (ILG)  !< Evaporation from vegetation over ground \f$[m s^{-1}]\f$
  real, intent(in) :: EVAPCG(ILG)  !< Evaporation from ground under vegetation \f$[m s^{-1}]\f$
  real, intent(inout) :: EVAPG (ILG)  !< Evaporation from bare ground \f$[m s^{-1}]\f$
  real, intent(inout) :: EVAPCS(ILG)  !< Evaporation from vegetation over snow \f$[m s^{-1}]\f$
  real, intent(in) :: EVPCSG(ILG)  !< Evaporation from snow under vegetation \f$[m s^{-1}]\f$
  real, intent(inout) :: EVAPGS(ILG)  !< Evaporation from snow on bare ground \f$[m s^{-1}]\f$
  real, intent(in) :: RPCP  (ILG)  !< Rainfall rate over modelled area \f$[m s^{-1}]\f$
  real, intent(in) :: TRPCP (ILG)  !< Rainfall temperature over modelled area [C]
  real, intent(in) :: SPCP  (ILG)  !< Snowfall rate over modelled area \f$[m s^{-1}]\f$
  real, intent(in) :: TSPCP (ILG)  !< Snowfall temperature over modelled area [C]
  real, intent(in) :: RHOSNI(ILG)  !< Density of fresh snow \f$[kg m^{-3}]\f$
  real, intent(in) :: ZPOND (ILG)  !< Depth of ponded water on surface [m]
  real, intent(in) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_s)\f$
  real, intent(in) :: ALBSNO(ILG)  !< Albedo of snow [ ]
  real, intent(in) :: WSNOCS(ILG)  !< Liquid water content of snow pack under vegetation \f$[kg m^{-2}] (w_s)\f$
  real, intent(in) :: WSNOGS(ILG)  !< Liquid water content of snow pack in bare areas \f$[kg m^{-2}] (w_s)\f$
  real, intent(in) :: RHOSCS(ILG)  !< Density of snow under vegetation \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(in) :: RHOSGS(ILG)  !< Density of snow in bare areas \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(in) :: TBASE (ILG)  !< Temperature of bedrock in third soil layer [K]
  real, intent(in) :: TSURX(ILG,4) !< Ground surface temperature over subarea [K]
  !
  real, intent(in) :: THLIQC(ILG,IG) !< Liquid water content of soil layers under vegetation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLIQG(ILG,IG) !< Liquid water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THICEC(ILG,IG) !< Frozen water content of soil layers under vegetation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THICEG(ILG,IG) !< Frozen water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
  real, intent(in) :: TBARC (ILG,IG) !< Subarea temperatures of soil layers \f$[C] (T_g)\f$
  real, intent(in) :: TBARG (ILG,IG) !< Subarea temperatures of soil layers \f$[C] (T_g)\f$
  real, intent(in) :: TBARCS(ILG,IG) !< Subarea temperatures of soil layers \f$[C] (T_g)\f$
  real, intent(in) :: TBARGS(ILG,IG) !< Subarea temperatures of soil layers \f$[C] (T_g)\f$
  real, intent(in) :: HCPC  (ILG,IG) !< Heat capacity of soil layers under vegetation \f$[J m^{-3} K^{-1}]\f$
  real, intent(in) :: HCPG  (ILG,IG) !<
  !
  !     * SOIL INFORMATION ARRAYS.
  !
  real, intent(in)  :: THPOR (ILG,IG) !< Pore volume in soil layer \f$[m^3 m^{-3}] (\theta_p)\f$
  real, intent(in)  :: HCPS  (ILG,IG) !< Heat capacity of soil material \f$[J m^{-3} K^{-1}]\f$
  real, intent(in)  :: GRKSAT(ILG,IG) !< Saturated hydraulic conductivity of soil layers \f$[m s^{-1}]\f$
  real, intent(out) :: DELZZ (ILG,IG) !< Soil layer depth variable used in soilWaterPhaseChg/waterUpdates [m]
  real, intent(in)  :: DELZW (ILG,IG) !< Permeable thickness of soil layer [m]
  real, intent(in)  :: DELZ  (IG)     !< Overall thickness of soil layer [m]
  !
  real BI(ILG,IG)
  real PSISAT(ILG,IG)
  real DD(ILG)
  real XSLOPE(ILG)
  real BULK_FC(ILG,IG)
  !
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  !     * INTERNAL WORK ARRAYS.
  !
  real, intent(inout) :: RADD(ILG), SADD(ILG)
  !
  !-----------------------------------------------------------------------
  !>
  !! In the first three loops, various subarea arrays and internal
  !! waterBudgetDriver variables are initialized.
  !!
  !     * INITIALIZE 2-D ARRAYS.
  !
  do J = 1,IG ! loop 50
    do I = IL1,IL2
      THLQCO(I,J) = 0.0
      THLQGO(I,J) = 0.0
      THLQCS(I,J) = 0.0
      THLQGS(I,J) = 0.0
      THICCO(I,J) = 0.0
      THICGO(I,J) = 0.0
      THICCS(I,J) = 0.0
      THICGS(I,J) = 0.0
      HCPCO (I,J) = 0.0
      HCPGO (I,J) = 0.0
      HCPCS (I,J) = 0.0
      HCPGS (I,J) = 0.0
      GRKSC (I,J) = 0.0
      GRKSG (I,J) = 0.0
      GRKSCS(I,J) = 0.0
      GRKSGS(I,J) = 0.0
      GFLXC (I,J) = 0.0
      GFLXG (I,J) = 0.0
      GFLXCS(I,J) = 0.0
      GFLXGS(I,J) = 0.0
      QFC   (I,J) = 0.0
      HMFG  (I,J) = 0.0
      THLDUM(I,J) = 0.0
      THIDUM(I,J) = 0.0
      SUBFLW(I,J) = 0.0
      TSUBFL(I,J) = 0.0
      if (J == 3 .and. IG == 3) then
        DELZZ (I,J) = DELZW(I,J)
      else
        DELZZ (I,J) = DELZ(J)
      end if
!c         calculate BULK_FC, make sure we don't divide by 0 (for slope = 0 or bedrock i.e. %sand = -2)
          if (XSLOPE(I) > 0.0 .and. THPOR(I,J) > 0.0) then
              BULK_FC(I,J) = (THPOR(I,J)/(BI(I,J)-1))* &
               ((2*DD(I)*PSISAT(I,J)*BI(I,J)/XSLOPE(I))**(1/BI(I,J)))* &
               ((3*BI(I,J)+2)**((BI(I,J)-1)/BI(I,J))- &
                (2*BI(I,J)+2)**((BI(I,J)-1)/BI(I,J)))
          else
!c         make sure that BULK_FC is high to shut off interflow
              BULK_FC(I,J) = 2.0
          end if
!d          print *, "I = ", I, "J = ", J, "BULK_FC(I,J) = ", BULK_FC(I,J)
    end do
  end do ! loop 50
  !
  do J = 1,IGP1 ! loop 75
    do I = IL1,IL2
      THLIQX(I,J) = 0.0
      THICEX(I,J) = 0.0
    end do
  end do ! loop 75
  !
  !     * INITIALIZE OTHER DIAGNOSTIC AND WORK ARRAYS.
  !
  do I = IL1,IL2 ! loop 100
    EVPICS(I) = EVAPCS(I) + EVPCSG(I)
    EVPIGS(I) = EVAPGS(I)
    EVPIC (I) = EVAPC (I) + EVAPCG(I)
    EVPIG (I) = EVAPG (I)
    TSNOWC(I) = 0.0
    TSNOWG(I) = 0.0
    WLOSTC(I) = 0.0
    WLOSTG(I) = 0.0
    WLSTCS(I) = 0.0
    WLSTGS(I) = 0.0
    RAC   (I) = RAICAN(I)
    RACS  (I) = RAICNS(I)
    SNC   (I) = SNOCAN(I)
    SNCS  (I) = SNOCNS(I)
    PCFC  (I) = 0.0
    PCLC  (I) = 0.0
    PCPN  (I) = 0.0
    PCPG  (I) = 0.0
    QFN   (I) = 0.0
    QFG   (I) = 0.0
    ROVG  (I) = 0.0
    ROFC  (I) = 0.0
    ROFN  (I) = 0.0
    OVRFLW(I) = 0.0
    BASFLW(I) = 0.0
    TOVRFL(I) = 0.0
    TBASFL(I) = 0.0
    ZPONDC(I) = 0.0
    ZPONDG(I) = 0.0
    ZPNDCS(I) = 0.0
    ZPNDGS(I) = 0.0
    XSNOWC(I) = 0.0
    XSNOWG(I) = 0.0
    XSNOCS(I) = 0.0
    XSNOGS(I) = 0.0
    ZSNOWC(I) = 0.0
    ZSNOWG(I) = 0.0
    ZSNOCS(I) = 0.0
    ZSNOGS(I) = 0.0
    ALBSC (I) = 0.0
    ALBSG (I) = 0.0
    ALBSCS(I) = 0.0
    ALBSGS(I) = 0.0
    RHOSC (I) = 0.0
    RHOSG (I) = 0.0
    HCPSC (I) = 0.0
    HCPSG (I) = 0.0
    HCPSCS(I) = 0.0
    HCPSGS(I) = 0.0
    RUNFC (I) = 0.0
    RUNFG (I) = 0.0
    RUNFCS(I) = 0.0
    RUNFGS(I) = 0.0
    TRUNFC(I) = 0.0
    TRUNFG(I) = 0.0
    TRNFCS(I) = 0.0
    TRNFGS(I) = 0.0
    TRUNOF(I) = 0.0
    TBASC (I) = TBASE(I) - TFREZ
    TBASG (I) = TBASE(I) - TFREZ
    TBASCS(I) = TBASE(I) - TFREZ
    TBASGS(I) = TBASE(I) - TFREZ
    DT    (I) = DELT
    RDUMMY(I) = 0.
    ZERO  (I) = 0.
    IZERO (I) = 0
    !
    !     * PRECIPITATION DIAGNOSTICS.
    !
    !>
    !! At the end of the 100 loop, a preliminary calculation of the
    !! precipitation diagnostics over the four subareas is carried
    !! out, as follows:
    !!
    !! The rainfall incident on the vegetation, PCLC, is summed over
    !! the vegetated subareas FC and FCS, minus the fraction that
    !! falls through the gaps in the canopy (denoted by the sky view
    !! factors FSVF and FSVFS respectively). The rainfall incident
    !! on the snowpack PCPN is the sum of the rainfall on the snow-
    !! covered bare area FGS, and that on the snow under the gaps in
    !! the canopy in subarea FCS. The rainfall incident on bare
    !! ground, PCPG, is the sum of the rainfall on the snow-free
    !! bare area FG, and that on the ground under the gaps in the
    !! canopy in subarea FC. The snowfall incident on the
    !! vegetation, PCFC, as in the case of rainfall, is summed over
    !! the vegetation subareas FC and FCS, minus the fraction that
    !! falls through the gaps in the canopy. The remaining amount is
    !! assigned to snowfall incident on the snow pack, PCPN.
    !!
    if (RPCP(I) > 0.) then
      PCLC(I) = (FCS(I) * (1.0 - FSVFS(I)) + FC(I) * (1.0 - FSVF(I))) * &
                RPCP(I) * RHOW
      PCPN(I) = (FCS(I) * FSVFS(I) + FGS(I)) * RPCP(I) * RHOW
      PCPG(I) = (FC(I) * FSVF(I) + FG(I)) * RPCP(I) * RHOW
    end if
    !
    if (SPCP(I) > 0.) then
      PCFC(I) = (FCS(I) * (1.0 - FSVFS(I)) + FC(I) * (1.0 - FSVF(I))) * &
                SPCP(I) * RHOSNI(I)
      PCPN(I) = PCPN(I) + (FCS(I) * FSVFS(I) + FGS(I) + &
                FC(I) * FSVF(I) + FG(I)) * SPCP(I) * RHOSNI(I)
    end if
  end do ! loop 100
  !>
  !! In loops 200 to 550, each of the four subareas is addressed in
  !! turn. Additional variables are initialized, and an empirical
  !! correction is applied to the saturated hydraulic conductivity in
  !! each of the soil layers over the subarea, to account for the
  !! viscosity of water at the layer temperature (Dingman, 2002) \cite Lawrence_Dingman2002-sq :
  !!
  !! \f$K_{sat}' = (1.7915 * 10^{-3}) K_{sat} / [(2.0319 * 10^{-4})+(1.5883 * 10^{-3}) exp(-T_g^{0.9}/22.0)]\f$
  !!
  !
  !     * RAINFALL/SNOWFALL RATES AND OTHER INITIALIZATION PROCEDURES
  !     * OVER GRID CELL SUBAREAS. DOWNWARD WATER FLUXES ARE LUMPED
  !     * TOGETHER WITH PRECIPITATION, AND UPWARD AND DOWNWARD WATER
  !     * FLUXES CANCEL OUT. CORRECTION MADE TO SOIL SATURATED HYDRAULIC
  !     * CONDUCTIVITY FOR WATER VISCOSITY EFFECTS.
  !
  !     * CALCULATIONS FOR CANOPY OVER SNOW.
  !
  if (NLANDCS > 0) then
    !
    do J = 1,IG ! loop 200
      do I = IL1,IL2
        if (FCS(I) > 0.) then
          THLQCS(I,J) = THLIQC(I,J)
          THICCS(I,J) = THICEC(I,J)
          HCPCS (I,J) = HCPC  (I,J)
          if (THPOR(I,J) > 0.0001) then
            GRKSCS(I,J) = GRKSAT(I,J) * (1.7915E-03 / &
                          (2.0319E-04 + 1.5883E-03 * EXP( - ((MAX(0.0, &
                          MIN(100.,TBARCS(I,J))) ** 0.9) / 22.))))
          else
            GRKSCS(I,J) = GRKSAT(I,J)
          end if
        end if
      end do
    end do ! loop 200
    !
    !>
    !! Next, a preliminary calculation of the water vapour flux
    !! diagnostics for each subarea is carried out. Over the
    !! vegetated subareas, first the canopy evaporative flux is
    !! assigned to sublimation if there is snow present on the
    !! canopy (since it is assumed that any water present exists
    !! within or underneath the snow). Then, if the sublimation rate
    !! is upward, it is assigned to the diagnostic variable QFCF; if
    !! it is downward, the portion of the flux corresponding to the
    !! canopy-covered area is assigned to QFCF and that
    !! corresponding to the gap area to QFN. Similarly, if the
    !! evaporation rate is upward it is assigned to the diagnostic
    !! variable QFCL; if it is downward, the flux corresponding to
    !! the canopy-covered area is assigned to QFCL and that
    !! corresponding to the gap area to QFN over subarea FCS and to
    !! QFG over subarea FC. Over the non-vegetated subareas, the
    !! evaporation rates are assigned to QFN for subarea GS, and to
    !! QFG for subarea G.
    !!
    !! For the purposes of the subsequent water balance calculations
    !! done in the other waterBudgetDriver subroutines, the subarea snowfall is
    !! lumped together with any simultaneously occurring
    !! sublimation, and the subarea rainfall with any simultaneously
    !! occurring evaporation. Depending on whether the sum of the
    !! snowfall and the sublimation, and the sum of the rainfall and
    !! the evaporation, are positive (downward) or negative
    !! (upward), corrections are applied to the appropriate
    !! diagnostic variables, and the snowfall/rainfall are set to
    !! the net flux if downward and the sublimation/evaporation are
    !! set to the net flux if upward. The smaller of the two fluxes
    !! in the sums are set to zero.
    !!
    !! Finally, ponded water and snow pack physical characteristics
    !! are set, including the snow heat capacity, which is
    !! calculated from the heat capacities of ice and water \f$C_i\f$
    !! and \f$C_w\f$, the snow, ice and water densities \f$\rho_s\f$ \f$\rho_i\f$, and
    !! \f$\rho_w\f$, and the water content and depth of the snow pack \f$w_s\f$
    !! and \f$z_s\f$, as follows:
    !!
    !! \f$C_s = C_i [\rho_s /\rho_i] + C_w w_s /[\rho_w z_s]\f$
    !!
    do I = IL1,IL2 ! loop 250
      if (FCS(I) > 0.) then
        if (SNOCNS(I) > 0.) then
          SUBLCS(I) = EVAPCS(I)
          EVAPCS(I) = 0.0
        else
          SUBLCS(I) = 0.0
        end if
        if (SUBLCS(I) > 0.0) then
          QFCF(I) = QFCF(I) + FCS(I) * SUBLCS(I) * RHOW
        else
          QFCF(I) = QFCF(I) + FCS(I) * (1.0 - FSVFS(I)) * SUBLCS(I) * &
                    RHOW
          QFN(I) = QFN(I) + FCS(I) * FSVFS(I) * SUBLCS(I) * RHOW
        end if
        if (EVAPCS(I) > 0.0) then
          QFCL(I) = QFCL(I) + FCS(I) * EVAPCS(I) * RHOW
        else
          QFCL(I) = QFCL(I) + FCS(I) * (1.0 - FSVFS(I)) * EVAPCS(I) * &
                    RHOW
          QFN(I) = QFN(I) + FCS(I) * FSVFS(I) * EVAPCS(I) * RHOW
        end if
        !
        if (SPCP(I) > 0. .or. SUBLCS(I) < 0.) then
          SADD(I) = SPCP(I) - SUBLCS(I) * RHOW / RHOSNI(I)
          if (ABS(SADD(I)) < 1.0E-12) SADD(I) = 0.0
          if (SADD(I) > 0.0) then
            if (SUBLCS(I) > 0.) then
              QFCF(I) = QFCF(I) - FCS(I) * FSVFS(I) * &
                        SUBLCS(I) * RHOW
              QFN(I) = QFN(I) + FCS(I) * FSVFS(I) * &
                       SUBLCS(I) * RHOW
            end if
            SPCCS (I) = SADD(I)
            if (SPCP(I) > 0.0) then
              TSPCCS(I) = TSPCP(I) + TFREZ
            else
              TSPCCS(I) = MIN(TSURX(I,1),TFREZ)
            end if
            SUBLCS(I) = 0.0
          else
            PCPN(I) = PCPN(I) - FCS(I) * FSVFS(I) * SPCP(I) * &
                      RHOSNI(I)
            PCFC(I) = PCFC(I) + FCS(I) * FSVFS(I) * SPCP(I) * &
                      RHOSNI(I)
            SUBLCS(I) = - SADD(I) * RHOSNI(I) / RHOW
            SPCCS (I) = 0.0
            TSPCCS(I) = 0.0
          end if
        else
          SPCCS(I) = 0.0
          TSPCCS(I) = 0.0
        end if
        !
        if (RPCP(I) > 0. .or. EVAPCS(I) < 0.) then
          RADD(I) = RPCP(I) - EVAPCS(I)
          if (ABS(RADD(I)) < 1.0E-12) RADD(I) = 0.0
          if (RADD(I) > 0.) then
            if (EVAPCS(I) > 0.) then
              QFCL(I) = QFCL(I) - FCS(I) * FSVFS(I) * &
                        EVAPCS(I) * RHOW
              QFN(I) = QFN(I) + FCS(I) * FSVFS(I) * &
                       EVAPCS(I) * RHOW
            end if
            RPCCS (I) = RADD(I)
            if (RPCP(I) > 0.0) then
              TRPCCS(I) = TRPCP(I) + TFREZ
            else
              TRPCCS(I) = MAX(TSURX(I,1),TFREZ)
            end if
            EVAPCS(I) = 0.0
          else
            PCPN(I) = PCPN(I) - FCS(I) * FSVFS(I) * RPCP(I) * RHOW
            PCLC(I) = PCLC(I) + FCS(I) * FSVFS(I) * RPCP(I) * RHOW
            EVAPCS(I) = - RADD(I)
            RPCCS (I) = 0.0
            TRPCCS(I) = 0.0
          end if
        else
          RPCCS(I) = 0.0
          TRPCCS(I) = 0.0
        end if
        ZPNDCS(I) = ZPOND (I)
        ZSNOCS(I) = ZSNOW (I)
        ALBSCS(I) = ALBSNO(I)
        HCPSCS(I) = HCPICE * RHOSCS(I) / RHOICE + HCPW * WSNOCS(I) / &
                    (RHOW * ZSNOCS(I))
        QFN   (I) = QFN(I) + FCS(I) * EVPCSG(I) * RHOW
      end if
    end do ! loop 250
  end if
  !
  !     * CALCULATIONS FOR SNOW-COVERED GROUND.
  !
  if (NLANDGS > 0) then
    !
    do J = 1,IG ! loop 300
      do I = IL1,IL2
        if (FGS(I) > 0.) then
          THLQGS(I,J) = THLIQG(I,J)
          THICGS(I,J) = THICEG(I,J)
          HCPGS (I,J) = HCPG  (I,J)
          if (THPOR(I,J) > 0.0001) then
            GRKSGS(I,J) = GRKSAT(I,J) * (1.7915E-03 / &
                          (2.0319E-04 + 1.5883E-03 * EXP( - ((MAX(0.0, &
                          MIN(100.,TBARGS(I,J))) ** 0.9) / 22.))))
          else
            GRKSGS(I,J) = GRKSAT(I,J)
          end if
        end if
      end do
    end do ! loop 300
    !
    do I = IL1,IL2 ! loop 350
      if (FGS(I) > 0.) then
        QFN(I) = QFN(I) + FGS(I) * EVAPGS(I) * RHOW
        if (SPCP(I) > 0. .or. EVAPGS(I) < 0.) then
          SADD(I) = SPCP(I) - EVAPGS(I) * RHOW / RHOSNI(I)
          if (ABS(SADD(I)) < 1.0E-12) SADD(I) = 0.0
          if (SADD(I) > 0.0) then
            SPCGS (I) = SADD(I)
            if (SPCP(I) > 0.0) then
              TSPCGS(I) = TSPCP(I)
            else
              TSPCGS(I) = MIN((TSURX(I,2) - TFREZ),0.0)
            end if
            EVAPGS(I) = 0.0
          else
            EVAPGS(I) = - SADD(I) * RHOSNI(I) / RHOW
            SPCGS (I) = 0.0
            TSPCGS(I) = 0.0
          end if
        else
          SPCGS (I) = 0.0
          TSPCGS(I) = 0.0
        end if
        !
        if (RPCP(I) > 0.) then
          RADD(I) = RPCP(I) - EVAPGS(I)
          if (ABS(RADD(I)) < 1.0E-12) RADD(I) = 0.0
          if (RADD(I) > 0.) then
            RPCGS (I) = RADD(I)
            TRPCGS(I) = TRPCP(I)
            EVAPGS(I) = 0.0
          else
            EVAPGS(I) = - RADD(I)
            RPCGS (I) = 0.0
            TRPCGS(I) = 0.0
          end if
        else
          RPCGS (I) = 0.0
          TRPCGS(I) = 0.0
        end if
        ZPNDGS(I) = ZPOND (I)
        ZSNOGS(I) = ZSNOW (I)
        ALBSGS(I) = ALBSNO(I)
        HCPSGS(I) = HCPICE * RHOSGS(I) / RHOICE + HCPW * WSNOGS(I) / &
                    (RHOW * ZSNOGS(I))
      end if
    end do ! loop 350
  end if
  !
  !     * CALCULATIONS FOR CANOPY OVER BARE GROUND.
  !
  if (NLANDC > 0) then
    !
    do J = 1,IG ! loop 400
      do I = IL1,IL2
        if (FC(I) > 0.) then
          THLQCO(I,J) = THLIQC(I,J)
          THICCO(I,J) = THICEC(I,J)
          HCPCO (I,J) = HCPC  (I,J)
          if (THPOR(I,J) > 0.0001) then
            GRKSC (I,J) = GRKSAT(I,J) * (1.7915E-03 / &
                          (2.0319E-04 + 1.5883E-03 * EXP( - ((MAX(0.0, &
                          MIN(100.,TBARC(I,J))) ** 0.9) / 22.))))
          else
            GRKSC (I,J) = GRKSAT(I,J)
          end if
        end if
      end do
    end do ! loop 400
    !
    do I = IL1,IL2 ! loop 450
      if (FC(I) > 0.) then
        if (SNOCAN(I) > 0.) then
          SUBLC(I) = EVAPC(I)
          EVAPC(I) = 0.0
        else
          SUBLC(I) = 0.0
        end if
        if (SUBLC(I) > 0.0) then
          QFCF(I) = QFCF(I) + FC(I) * SUBLC(I) * RHOW
        else
          QFCF(I) = QFCF(I) + FC(I) * (1.0 - FSVF(I)) * SUBLC(I) * &
                    RHOW
          QFN(I) = QFN(I) + FC(I) * FSVF(I) * SUBLC(I) * RHOW
        end if
        if (EVAPC(I) > 0.0) then
          QFCL(I) = QFCL(I) + FC(I) * EVAPC(I) * RHOW
        else
          QFCL(I) = QFCL(I) + FC(I) * (1.0 - FSVF(I)) * EVAPC(I) * &
                    RHOW
          QFG(I) = QFG(I) + FC(I) * FSVF(I) * EVAPC(I) * RHOW
        end if
        !
        if (SPCP(I) > 0. .or. SUBLC(I) < 0.) then
          SADD(I) = SPCP(I) - SUBLC(I) * RHOW / RHOSNI(I)
          if (ABS(SADD(I)) < 1.0E-12) SADD(I) = 0.0
          if (SADD(I) > 0.0) then
            if (SUBLC(I) > 0.) then
              QFCF(I) = QFCF(I) - FC(I) * FSVF(I) * SUBLC(I) * &
                        RHOW
              QFN(I) = QFN(I) + FC(I) * FSVF(I) * SUBLC(I) * &
                       RHOW
            end if
            SPCC  (I) = SADD(I)
            if (SPCP(I) > 0.0) then
              TSPCC (I) = TSPCP(I) + TFREZ
            else
              TSPCC(I) = MIN(TSURX(I,3),TFREZ)
            end if
            SUBLC (I) = 0.0
          else
            PCPN(I) = PCPN(I) - FC(I) * FSVF(I) * SPCP(I) * &
                      RHOSNI(I)
            PCFC(I) = PCFC(I) + FC(I) * FSVF(I) * SPCP(I) * &
                      RHOSNI(I)
            SUBLC (I) = - SADD(I) * RHOSNI(I) / RHOW
            SPCC  (I) = 0.0
            TSPCC (I) = 0.0
          end if
        else
          SPCC  (I) = 0.0
          TSPCC (I) = 0.0
        end if
        !
        if (RPCP(I) > 0. .or. EVAPC(I) < 0.) then
          RADD(I) = RPCP(I) - EVAPC(I)
          if (ABS(RADD(I)) < 1.0E-12) RADD(I) = 0.0
          if (RADD(I) > 0.) then
            if (EVAPC(I) > 0.) then
              QFCL(I) = QFCL(I) - FC(I) * FSVF(I) * EVAPC(I) * &
                        RHOW
              QFG(I) = QFG(I) + FC(I) * FSVF(I) * EVAPC(I) * &
                       RHOW
            end if
            RPCC  (I) = RADD(I)
            if (RPCP(I) > 0.0) then
              TRPCC (I) = TRPCP(I) + TFREZ
            else
              TRPCC (I) = MAX(TSURX(I,3),TFREZ)
            end if
            EVAPC (I) = 0.0
          else
            PCPG(I) = PCPG(I) - FC(I) * FSVF(I) * RPCP(I) * RHOW
            PCLC(I) = PCLC(I) + FC(I) * FSVF(I) * RPCP(I) * RHOW
            EVAPC (I) = - RADD(I)
            RPCC  (I) = 0.0
            TRPCC (I) = 0.0
          end if
        else
          RPCC  (I) = 0.0
          TRPCC (I) = 0.0
        end if
        ZPONDC(I) = ZPOND (I)
        ZSNOWC(I) = 0.
        RHOSC (I) = 0.
        HCPSC (I) = 0.
        QFG   (I) = QFG(I) + FC(I) * EVAPCG(I) * RHOW
      end if
    end do ! loop 450
  end if
  !
  !     * CALCULATIONS FOR BARE GROUND.
  !
  if (NLANDG > 0) then
    !
    do J = 1,IG ! loop 500
      do I = IL1,IL2
        if (FG(I) > 0.) then
          THLQGO(I,J) = THLIQG(I,J)
          THICGO(I,J) = THICEG(I,J)
          HCPGO (I,J) = HCPG  (I,J)
          if (THPOR(I,J) > 0.0001) then
            GRKSG (I,J) = GRKSAT(I,J) * (1.7915E-03 / &
                          (2.0319E-04 + 1.5883E-03 * EXP( - ((MAX(0.0, &
                          MIN(100.,TBARG(I,J))) ** 0.9) / 22.))))
          else
            GRKSG (I,J) = GRKSAT(I,J)
          end if
        end if
      end do
    end do ! loop 500
    !
    do I = IL1,IL2 ! loop 550
      if (FG(I) > 0.) then
        QFG(I) = QFG(I) + FG(I) * EVAPG(I) * RHOW
        if (SPCP(I) > 0.) then
          SADD(I) = SPCP(I) - EVAPG(I) * RHOW / RHOSNI(I)
          if (ABS(SADD(I)) < 1.0E-12) SADD(I) = 0.0
          if (SADD(I) > 0.0) then
            QFN(I) = QFN(I) + FG(I) * EVAPG(I) * RHOW
            QFG(I) = QFG(I) - FG(I) * EVAPG(I) * RHOW
            SPCG  (I) = SADD(I)
            TSPCG (I) = TSPCP(I)
            EVAPG (I) = 0.0
          else
            PCPN(I) = PCPN(I) - FG(I) * SPCP(I) * RHOSNI(I)
            PCPG(I) = PCPG(I) + FG(I) * SPCP(I) * RHOSNI(I)
            EVAPG (I) = - SADD(I) * RHOSNI(I) / RHOW
            SPCG  (I) = 0.0
            TSPCG (I) = 0.0
          end if
        else
          SPCG  (I) = 0.0
          TSPCG (I) = 0.0
        end if
        !
        if (RPCP(I) > 0. .or. EVAPG(I) < 0.) then
          RADD(I) = RPCP(I) - EVAPG(I)
          if (ABS(RADD(I)) < 1.0E-12) RADD(I) = 0.0
          if (RADD(I) > 0.) then
            RPCG  (I) = RADD(I)
            if (RPCP(I) > 0.0) then
              TRPCG (I) = TRPCP(I)
            else
              TRPCG (I) = MAX((TSURX(I,4) - TFREZ),0.0)
            end if
            EVAPG (I) = 0.0
          else
            EVAPG (I) = - RADD(I)
            RPCG  (I) = 0.0
            TRPCG (I) = 0.0
          end if
        else
          RPCG (I) = 0.0
          TRPCG(I) = 0.0
        end if
        ZPONDG(I) = ZPOND (I)
        ZSNOWG(I) = 0.
        RHOSG (I) = 0.
        HCPSG (I) = 0.
      end if
    end do ! loop 550
  end if

  return
end subroutine waterCalcPrep
