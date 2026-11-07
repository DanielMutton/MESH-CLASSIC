!> \file
!! Initializes subarea variables and calculate various
!! parameters for surface energy budget calculations.
!! @author D. Verseghy, M. Lazare, A. Wu, Y. Delage, J. P. Paquin,
!! R. Harvey, J. Melton, G. Meyer, E. Humphreys
!
subroutine energyBudgetPrep (THLIQC, THLIQG, THICEC, THICEG, TBARC, TBARG, & ! Formerly TPREP
                             TBARCS, TBARGS, HCPC, HCPG, TCTOPC, TCBOTC, &
                             TCTOPG, TCBOTG, HCPSCS, HCPSGS, TCSNOW, TSNOCS, &
                             TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, TCANO, &
                             TCANS, CEVAP, IEVAP, TBAR1P, WTABLE, ZERO, &
                             EVAPC, EVAPCG, EVAPG, EVAPCS, EVPCSG, EVAPGS, &
                             GSNOWC, GSNOWG, GZEROC, GZEROG, GZROCS, GZROGS, &
                             QMELTC, QMELTG, EVAP, GSNOW, &
                             TPONDC, TPONDG, TPNDCS, TPNDGS, QSENSC, QSENSG, &
                             QEVAPC, QEVAPG, TACCO, QACCO, TACCS, QACCS, &
                             ILMOX, UEX, HBLX, &
                             ILMO, UE, HBL, &
                             ST, SU, SV, SQ, SRH, &
                             CDH, CDM, QSENS, QEVAP, QLWAVG, &
                             FSGV, FSGS, FSGG, FLGV, FLGS, FLGG, &
                             HFSC, HFSS, HFSG, HEVC, HEVS, HEVG, &
                             HMFC, HMFN, QFCF, QFCL, EVPPOT, ACOND, &
                             DRAG, THLIQ, THICE, TBAR, ZPOND, TPOND, &
                             THPOR, THLMIN, THLRET, THFC, HCPS, TCS, &
                             TA, RHOSNO, TSNOW, ZSNOW, WSNOW, TCAN, &
                             FC, FCS, DELZ, DELZW, ZBOTW, &
                             ISAND, ILG, IL1, IL2, JL, IG, &
                             FVEG, TCSATU, TCSATF, FTEMP, FTEMPX, FVAP, &
                             FVAPX, RIB, RIBX)
  !
  !     * Aug 21/19 - G. Meyer, E. Humphreys, J. Melton.  Change calculation of
  !     *                         CEVAP to Merlin et al. (2011)
  !     * JUN 21/13 - M.LAZARE.   PASS IN AND INITIALIZE TO ZERO "GSNOW".
  !     * NOV 24/11 - R.HARVEY.   NEW SNOW THERMAL CONDUCTIVITY FROM
  !     *                         STURM ET AL. (1997).
  !     * OCT 12/11 - M.LAZARE.   REMOVED TSURF.
  !     * AUG   /08 - J.P.PAQUIN. ADD CALCULATION FOR FTEMP, FVAP AND
  !     *                         RIB FOR OUTPUT IN GEM (IMPLEMENTED BY
  !     *                         L. DUARTE ON OCT. 28/08).
  !     * MAR 20/08 - D.VERSEGHY. REMOVE TBAR3, TCTOP3, TCBOT3.
  !     * DEC 12/07 - D.VERSEGHY. MAJOR REVISIONS TO CALCULATION OF
  !     *                         SOIL THERMAL CONDUCTIVITY.
  !     * MAY 18/06 - D.VERSEGHY. ADJUST CALCULATION OF TBAR1P FOR ROCK
  !     *                         SOILS; LIMIT CALCULATION OF TBAR3(I, 3)
  !     *                         TO UPPER 4.1 M OF SOIL; CORRECT WTABLE
  !     *                         TO ACCOUNT FOR PRESENCE OF ICE.
  !     * MAR 23/06 - D.VERSEGHY. MODIFY CALCULATION OF HCPSNO TO ACCOUNT
  !     *                         FOR PRESENCE OF WATER IN SNOWPACK.
  !     * MAR 21/06 - P.BARTLETT. INITIALIZE ADDITIONAL VARIABLES TO ZERO.
  !     * OCT 04/05 - D.VERSEGHY. NEW VARIABLES TBAR3, TCTOP3, TCBOT3.
  !     * APR 08/05 - Y.DELAGE. TCTOP VARIES GRADUALLY WITH ZPOND TO TCW.
  !     * MAR 16/05 - D.VERSEGHY. TREAT FROZEN SOIL WATER AS ICE
  !     *                         VOLUME RATHER THAN AS EQUIVALENT
  !     *                         LIQUID WATER VOLUME; REVERSE ORDER
  !     *                         IN LOOP 500.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 05/04 - Y.DELAGE/D.VERSEGHY. INITIALIZE NEW DIAGNOSTIC
  !     *                         VARIABLES ILMO, UE AND HBL.
  !     * JUL 30/02 - D.VERSEGHY. MOVE CALCULATION OF VEGETATION
  !     *                         STOMATAL RESISTANCE INTO calcLandSurfParams
  !     *                         AND canopyAlbedoTransmiss; SHORTENED CLASS3
  !     *                         COMMON BLOCK.
  !     * JUN 17/02 - D.VERSEGHY. NEW THERMAL ARRAYS FOR SURFACE
  !     *                         TEMPERATURE ITERATION, WITH PONDED
  !     *                         WATER ROLLED INTO SOIL UPPER LAYER
  !     *                         SHORTENED CLASS4 COMMON BLOCK.
  !     * MAR 20/02 - D.VERSEGHY. MOVE CALCULATION OF BACKGROUND SOIL
  !     *                         PROPERTIES INTO "soilProperties"; UPDATES
  !     *                         TO MAKE ZPOND A PROGNOSTIC VARIABLE.
  !     * FEB 27/02 - D.VERSEGHY. RECALCULATE WILTING POINT BASED ON
  !     *                         FIELD CAPACITY.
  !     * JAN 18/02 - D.VERSEGHY. INTRODUCTION OF CALCULATION OF FIELD
  !     *                         CAPACITY AND NEW BARE SOIL EVAPORATION
  !     *                         PARAMETERS.
  !     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
  !     * NOV 01/00 - A.WU/D.VERSEGHY. EXTEND MINERAL SOIL CALCULATION
  !     *                              OF SOIL EVAPORATION "BETA" TO
  !     *                              ORGANIC SOILS.
  !     * SEP 19/00 - A.WU/D.VERSEGHY. CHANGE CALCULATION OF THERMAL
  !     *                              CONDUCTIVITY FOR ORGANIC SOILS,
  !     *                              USING METHOD OF FAROUKI (1981).
  !     *                              ALSO, CALCULATE STOMATAL RESISTANCE
  !     *                              USING VEGETATION-VARYING
  !     *                              COEFFICIENTS FOR ENVIRONMENTAL
  !     *                              VARIABLES.
  !     * FEB 14/00 - D.VERSEGHY. INSERT CALCULATION OF WATER TABLE DEPTH
  !     *                         FOR ORGANIC SOILS.
  !     * DEC 07/99 - A.WU/D.VERSEGHY.  INCORPORATE CALCULATION OF "BETA"
  !     *                               PARAMETER FOR NEW SOIL EVAPORATION
  !     *                               FORMULATION.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         CHANGES RELATED TO VARIABLE SOIL DEPTH
  !     *                         (MOISTURE HOLDING CAPACITY) AND DEPTH-
  !     *                         VARYING SOIL PROPERTIES.
  !     * JAN 24/97 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         SET RC AND RCS TO ZERO FOR GRID CELLS
  !     *                         WITH NO VEGETATION.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REMOVE SUBTRACTION OF RESIDUAL SOIL
  !     *                         MOISTURE CONTENT IN CALCULATIONS OF
  !     *                         "PSIZRO" AND "PSII".
  !     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
  !     *                         ORGANIC MATTER CONTENT.
  !     * DEC 16/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         INITIALIZE THREE NEW DIAGNOSTIC FIELDS.
  !     * NOV 12/94 - D.VERSEGHY. SET INITIAL TEMPERATURE OF EMERGING
  !     *                         CANOPY TO TA INSTEAD OF TO ZERO.
  !     * JAN 31/94 - D.VERSEGHY. CLASS - VERSION 2.2.
  !     *                         INTRODUCE LIMITING VALUES INTO
  !     *                         CALCULATION OF "PSIZRO" TO AVOID
  !     *                         OVERFLOWS.
  !     * JUL 27/93 - D.VERSEGHY/M.LAZARE. INITIALIZE NEW DIAGNOSTIC
  !     *                                  FIELDS FSGV, FSGG, FLGV, FLGG,
  !     *                                  HFSC, HFSG, HMFC.
  !     * MAY 06/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  MODIFICATIONS TO CANOPY
  !     *                                  RESISTANCE TO ADD "RCS"
  !     *                                  FIELD FOR SNOW-COVERED
  !     *                                  CANOPY.
  !     * JUL 04/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. PREPARATION AND INITIALIZATION FOR
  !     *                         LAND SURFACE ENERGY BUDGET
  !     *                         CALCULATIONS.
  !
  use classicParams,       only : TCW, TCICE, TCSAND, HCPW, HCPICE, &
                                  HCPSND, RHOW, RHOICE, TCGLAC

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL, IG
  integer :: I, J
  !
  !     * OUTPUT ARRAYS.
  !
  !     (Suffix CS = vegetation over snow cover; GS = bare snow cover; C
  !     or CO = vegetation over ground; G or GO = bare ground.)
  !
  real, intent(out)   :: TBARC (ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(out)   :: TBARG (ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(out)   :: TBARCS(ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(out)   :: TBARGS(ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(inout) :: THLIQC(ILG,IG) !< Liquid water content of soil layers under vegetation \f$[m^3 m^{-3}]\f$
  real, intent(inout) :: THLIQG(ILG,IG) !< Liquid water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
  real, intent(inout) :: THICEC(ILG,IG) !< Frozen water content of soil layers under vegetation \f$[m^3 m^{-3}]\f$
  real, intent(inout) :: THICEG(ILG,IG) !< Frozen water content of soil layers in bare areas \f$[m^3 m^{-3}]\f$
  real, intent(out)   :: HCPC  (ILG,IG) !< Heat capacity of soil layers under vegetation \f$[J m^{-3} K^{-1}] (C_g)\f$
  real, intent(inout) :: HCPG  (ILG,IG) !< Heat capacity of soil layers in bare areas \f$[J m^{-3} K^{1}] (Cg)\f$
  real, intent(inout) :: TCTOPC(ILG,IG) !< Thermal conductivity of soil at top of layer in canopy-covered subareas
  !< \f$[W m^{-1} K^{-1}] (\lambda)\f$
  real, intent(out)   :: TCBOTC(ILG,IG) !< Thermal conductivity of soil at bottom of layer in canopy-covered subareas
  !< \f$[W m^{-1} K^{-1}] (\lambda)\f$
  real, intent(out)   :: TCTOPG(ILG,IG) !< Thermal conductivity of soil at top of
  !< layer in bare ground subareas \f$[W m^{-1} K^{-1}] (\lambda)\f$
  real, intent(out)   :: TCBOTG(ILG,IG) !< Thermal conductivity of soil at bottom of
  !< layer in bare ground subareas \f$[W m^{-1} K^{-1}] (\lambda)\f$
  !
  real, intent(inout) :: HCPSCS(ILG)    !< Heat capacity of snow pack under vegetation canopy \f$[J m^{-3} K^1] (C_s)\f$
  real, intent(out)   :: HCPSGS(ILG)    !< Heat capacity of snow pack in bare areas \f$[J m^{-3} K^{1}] (C_s)\f$
  real, intent(out)   :: TCSNOW(ILG)    !< Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
  real, intent(out)   :: TSNOGS(ILG)    !< Temperature of snow pack in bare areas [K]
  real, intent(out)   :: TSNOCS(ILG)    !< Temperature of snow pack under vegetation canopy [K]
  real, intent(out)   :: WSNOCS(ILG)    !< Liquid water content of snow pack under vegetation \f$[kg m^{-2}]\f$
  real, intent(out)   :: WSNOGS(ILG)    !< Liquid water content of snow pack in bare areas \f$[kg m^{-2}]\f$
  real, intent(out)   :: RHOSCS(ILG)    !< Density of snow pack under vegetation canopy \f$[kg m^{-3}]\f$
  real, intent(out)   :: RHOSGS(ILG)    !< Density of snow pack in bare areas \f$[kg m^{-3}]\f$
  real, intent(out)   :: TCANO (ILG)    !< Temperature of canopy over ground [K]
  real, intent(out)   :: TCANS (ILG)    !< Temperature of canopy over snow [K]
  real, intent(out)   :: CEVAP (ILG)    !< Soil evaporation efficiency coefficient \f$[ ] (\beta)\f$
  real, intent(out)   :: TBAR1P(ILG)    !< Lumped temperature of ponded water and first soil layer [K]
  real, intent(out)   :: WTABLE(ILG)    !< Depth of water table in soil \f$[m] (z_{wt})\f$
  real, intent(out)   :: ZERO  (ILG)    !< Dummy vector containing all zeros
  real, intent(out)   :: TPONDC(ILG)    !< Subarea temperature of surface ponded water [C]
  real, intent(out)   :: TPONDG(ILG)    !< Subarea temperature of surface ponded water [C]
  real, intent(out)   :: TPNDCS(ILG)    !< Subarea temperature of surface ponded water [C]
  real, intent(out)   :: TPNDGS(ILG)    !< Subarea temperature of surface ponded water [C]
  !
  integer, intent(out) :: IEVAP (ILG) !< Flag indicating whether soil evaporation is occurring or not
  !
  !     * OUTPUT ARRAYS WHICH ARE INTERNAL WORK ARRAYS FOR energyBudgetDriver
  !     * AND ARE INITIALIZED TO ZERO HERE.
  !
  real, intent(out) :: EVAPC (ILG)  !< Evaporation from vegetation over ground \f$[m s^{-1}]\f$
  real, intent(out) :: EVAPCG(ILG)  !< Evaporation from ground under vegetation \f$[m s^{-1}]\f$
  real, intent(out) :: EVAPG (ILG)  !< Evaporation from bare ground \f$[m s^{-1}]\f$
  real, intent(out) :: EVAPCS(ILG)  !< Evaporation from vegetation over snow \f$[m s^{-1}]\f$
  real, intent(out) :: EVPCSG(ILG)  !< Evaporation from snow under vegetation \f$[m s^{-1}]\f$
  real, intent(out) :: EVAPGS(ILG)  !< Evaporation from snow on bare ground \f$[m s^{-1}]\f$
  real, intent(out) :: GSNOWC(ILG)  !< Heat flux at top of snow pack under canopy \f$[W m^{-2}]\f$
  real, intent(out) :: GSNOWG(ILG)  !< Heat flux at top of snow pack over bare ground \f$[W m^{-2}]\f$
  real, intent(out) :: GZEROC(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: GZEROG(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: GZROCS(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: GZROGS(ILG)  !< Subarea heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: QMELTC(ILG)  !< Heat to be used for melting snow under canopy \f$[W m^{-2}]\f$
  real, intent(out) :: QMELTG(ILG)  !< Heat to be used for melting snow on bare ground \f$[W m^{-2}]\f$
  real, intent(out) :: EVAP  (ILG)  !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: QSENSC(ILG)  !< Sensible heat flux from vegetation canopy over subarea \f$[W m^{-2}]\f$
  real, intent(out) :: QSENSG(ILG)  !< Sensible heat flux from ground over subarea \f$[W m^{-2}]\f$
  real, intent(out) :: QEVAPC(ILG)  !< Latent heat flux from vegetation canopy over subarea \f$[W m^{-2}]\f$
  real, intent(out) :: QEVAPG(ILG)  !< Latent heat flux from ground over subarea \f$[W m^{-2}]\f$
  real, intent(out) :: TACCO (ILG)  !< Temperature of air within vegetation canopy space over bare ground [K]
  real, intent(out) :: QACCO (ILG)  !< Specific humidity of air within vegetation canopy space over bare ground \f$[kg kg^{-1}]\f$
  real, intent(out) :: TACCS (ILG)  !< Temperature of air within vegetation canopy space over snow [K]
  real, intent(out) :: QACCS (ILG)  !< Specific humidity of air within vegetation canopy space over snow \f$[kg kg^{-1}]\f$
  real, intent(out) :: GSNOW (ILG)  !< Diagnostic heat flux at snow surface for use in CCCma black carbon deposition scheme \f$[W m^{-2}]\f$

  !
  !     * DIAGNOSTIC ARRAYS.
  !
  real, intent(out) :: ST    (ILG)  !< Diagnosed screen-level air temperature [K]
  real, intent(out) :: SU    (ILG)  !< Diagnosed anemometer-level zonal wind \f$[m s^{-1}]\f$
  real, intent(out) :: SV    (ILG)  !< Diagnosed anemometer-level meridional wind \f$[m s^{-1}]\f$
  real, intent(out) :: SQ    (ILG)  !< Diagnosed screen-level specific humidity \f$[kg kg^{-1}]\f$
  real, intent(out) :: SRH   (ILG)  !< Diagnosed screen-level relative humidity  [%]
  real, intent(out) :: CDH   (ILG)  !< Surface drag coefficient for heat [ ]
  real, intent(out) :: CDM   (ILG)  !< Surface drag coefficient for momentum [ ]
  real, intent(out) :: QSENS (ILG)  !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2}]\f$
  real, intent(out) :: QEVAP (ILG)  !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2}]\f$
  real, intent(out) :: QLWAVG(ILG)  !< Upwelling longwave radiation over modelled area \f$[W m^{-2}]\f$
  real, intent(out) :: FSGV  (ILG)  !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2}]\f$
  real, intent(out) :: FSGS  (ILG)  !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2}]\f$
  real, intent(out) :: FSGG  (ILG)  !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: FLGV  (ILG)  !< Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2}]\f$
  real, intent(out) :: FLGS  (ILG)  !< Diagnosed net longwave radiation at snow surface \f$[W m^{-2}]\f$
  real, intent(out) :: FLGG  (ILG)  !< Diagnosed net longwave radiation at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: HFSC  (ILG)  !< Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2}]\f$
  real, intent(out) :: HFSS  (ILG)  !< Diagnosed sensible heat flux at snow surface \f$[W m^{-2}]\f$
  real, intent(out) :: HFSG  (ILG)  !< Diagnosed sensible heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: HEVC  (ILG)  !< Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2}]\f$
  real, intent(out) :: HEVS  (ILG)  !< Diagnosed latent heat flux at snow surface \f$[W m^{-2}]\f$
  real, intent(out) :: HEVG  (ILG)  !< Diagnosed latent heat flux at soil surface \f$[W m^{-2}]\f$
  real, intent(out) :: HMFC  (ILG)  !< Diagnosed energy associated with phase change
  !< of water on vegetation \f$[W m^{-2}]\f$
  real, intent(out) :: HMFN  (ILG)  !< Diagnosed energy associated with phase change
  !< of water in snow pack \f$[W m^{-2}]\f$
  real, intent(out) :: EVPPOT(ILG)  !< Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1}]\f$
  real, intent(out) :: ACOND (ILG)  !< Diagnosed product of drag coefficient and wind
  !< speed over modelled area \f$[m s^{-1}]\f$
  real, intent(out) :: DRAG  (ILG)  !< Surface drag coefficient under neutral stability [ ]
  real, intent(out) :: ILMO  (ILG)  !< Surface drag coefficient under neutral stability [ ]
  real, intent(out) :: UE    (ILG)  !< Friction velocity of air \f$[m s^{-1}]\f$
  real, intent(out) :: HBL   (ILG)  !< Height of the atmospheric boundary layer [m]
  real, intent(out) :: ILMOX (ILG)  !< Inverse of Monin-Obukhov roughness length over each subarea \f$[m^{-1}]\f$
  real, intent(out) :: UEX   (ILG)  !< Friction velocity of air over each subarea \f$[m s^{-1}]\f$
  real, intent(out) :: HBLX  (ILG)  !< Height of the atmospheric boundary layer over each subarea [m]

  real, intent(out) :: QFCF(ILG), QFCL(ILG), FTEMP(ILG), FTEMPX(ILG), &
                       FVAP(ILG), FVAPX(ILG), RIB(ILG), RIBX(ILG)

  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
  !< layers \f$[m^3 m^{-3}] (\theta_l)\f$
  real, intent(in) :: THICE (ILG,IG)   !< Volumetric frozen water content of soil
  !< layers \f$[m^3 m^{-3}] (\theta_i)\f$
  real, intent(in) :: TBAR  (ILG,IG)   !< Temperature of soil layers [K]
  real, intent(in) :: ZPOND (ILG)      !< Depth of ponded water on surface [m]
  real, intent(in) :: TPOND (ILG)      !< Temperature of ponded water [K]
  !
  real, intent(in) :: TA    (ILG)      !< Air temperature at reference height [K]
  real, intent(in) :: RHOSNO(ILG)      !< Density of snow \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(in) :: TSNOW (ILG)      !< Snowpack temperature [K]
  real, intent(in) :: ZSNOW (ILG)      !< Depth of snow pack \f$[m] (z_s)\f$
  real, intent(in) :: WSNOW (ILG)      !< Liquid water content of snow pack \f$[kg m^{-2}] (w_s)\f$
  real, intent(in) :: TCAN  (ILG)      !< Vegetation canopy temperature [K]
  real, intent(in) :: FC    (ILG)      !< Fractional coverage of canopy over bare ground for modelled area [ ]
  real, intent(in) :: FCS   (ILG)      !< Fractional coverage of canopy over snow for modelled area [ ]
  !
  !     * SOIL PROPERTY ARRAYS.
  !
  real, intent(in) :: THPOR(ILG,IG)    !< Pore volume in soil layer \f$[m^3 m^{-3}] (\theta_p)\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !< remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLRET(ILG,IG)   !< Liquid water retention capacity for organic soil \f$[m^3 m^{-3}] (\theta_{ret})\f$
  real, intent(in) :: THFC  (ILG,IG)   !< Field capacity \f$[m^3 m^{-3}] (theta_{fc})\f$
  real, intent(in) :: HCPS  (ILG,IG)   !< Heat capacity of soil material \f$[J m^{-3} K^{-1}] (C_m)\f$
  real, intent(in) :: TCS   (ILG,IG)   !< Thermal conductivity of soil particles \f$[W m^{-1} K^{-1}] (\theta_s)\f$
  real, intent(in) :: DELZW(ILG,IG)    !< Permeable thickness of soil layer \f$[m] (\Delta z_w)\f$
  real, intent(in) :: ZBOTW(ILG,IG)    !< Depth to permeable bottom of soil layer \f$[m] (z_{b,w})\f$
  real, intent(in) :: DELZ(IG)         !< Overall thickness of soil layer [m]
  !
  integer, intent(in) :: ISAND (ILG,IG)!< Sand content flag
  !
  !     * INTERNAL WORK FIELDS FOR THIS ROUTINE.
  !
  real, intent(inout) :: FVEG(ILG), TCSATU(ILG), TCSATF(ILG)
  !
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: SATRAT, THLSAT, THISAT, TCDRY, TCKAPU, TCKAPF, TCRATU, TCRATF, &
          TCSOLU, TCSOLF, TCSOIL, TSUM1, TSUM2, ZSUM
  integer :: IWTABL(ILG)
  !
  !----------------------------------------------------------------------
  !
  !>
  !! In the first two loops, various subarea arrays and internal
  !! energyBudgetDriver variables are initialized. The initial temperatures of the
  !! vegetation canopy above snow and above bare ground (TCANS and
  !! TCANO) are set to the temperature of the vegetation over the
  !! whole modelled area (TCAN) if TCAN is not effectively 0 K (the
  !! value it is assigned if vegetation is not present). Otherwise,
  !! the canopy temperatures are initialized to the air temperature
  !! TA.
  !!
  !     * INITIALIZE 2-D AND 3-D ARRAYS.
  !
  do J = 1,IG ! loop 50
    do I = IL1,IL2
      THLIQG(I,J) = THLIQ(I,J)
      THICEG(I,J) = THICE(I,J)
      THLIQC(I,J) = THLIQ(I,J)
      THICEC(I,J) = THICE(I,J)
      TBARCS(I,J) = 0.0
      TBARGS(I,J) = 0.0
      TBARC (I,J) = 0.0
      TBARG (I,J) = 0.0
      TCTOPC(I,J) = 0.0
      TCBOTC(I,J) = 0.0
      TCTOPG(I,J) = 0.0
      TCBOTG(I,J) = 0.0
    end do
  end do ! loop 50
  !
  !     * INITIALIZE 1-D INTERNAL WORK FIELDS AND DIAGNOSTIC ARRAYS.
  !
  do I = IL1,IL2 ! loop 100
    FVEG  (I) = FC(I) + FCS(I)
    if (TCAN(I) > 5.0) then
      TCANS (I) = TCAN(I)
      TCANO (I) = TCAN(I)
    else
      TCANS (I) = TA(I)
      TCANO (I) = TA(I)
    end if
    EVAPC (I) = 0.
    EVAPCG(I) = 0.
    EVAPG (I) = 0.
    EVAPCS(I) = 0.
    EVPCSG(I) = 0.
    EVAPGS(I) = 0.
    GSNOWC(I) = 0.
    GSNOWG(I) = 0.
    GZEROC(I) = 0.
    GZEROG(I) = 0.
    GZROCS(I) = 0.
    GZROGS(I) = 0.
    QMELTC(I) = 0.
    QMELTG(I) = 0.
    QSENSC(I) = 0.
    QSENSG(I) = 0.
    QEVAPC(I) = 0.
    QEVAPG(I) = 0.
    TPONDC(I) = 0.
    TPONDG(I) = 0.
    TPNDCS(I) = 0.
    TPNDGS(I) = 0.
    TACCS (I) = 0.
    QACCS (I) = 0.
    TACCO (I) = 0.
    QACCO (I) = 0.
    ST    (I) = 0.
    SU    (I) = 0.
    SV    (I) = 0.
    SQ    (I) = 0.
    SRH   (I) = 0.
    CDH   (I) = 0.
    CDM   (I) = 0.
    QSENS (I) = 0.
    QEVAP (I) = 0.
    EVAP  (I) = 0.
    QLWAVG(I) = 0.
    FSGV  (I) = 0.
    FSGS  (I) = 0.
    FSGG  (I) = 0.
    FLGV  (I) = 0.
    FLGS  (I) = 0.
    FLGG  (I) = 0.
    HFSC  (I) = 0.
    HFSS  (I) = 0.
    HFSG  (I) = 0.
    HEVC  (I) = 0.
    HEVS  (I) = 0.
    HEVG  (I) = 0.
    HMFC  (I) = 0.
    HMFN  (I) = 0.
    QFCF  (I) = 0.
    QFCL  (I) = 0.
    EVPPOT(I) = 0.
    ACOND (I) = 0.
    DRAG  (I) = 0.
    ILMO  (I) = 0.
    UE    (I) = 0.
    HBL   (I) = 0.
    ILMOX (I) = 0.
    UEX   (I) = 0.
    HBLX  (I) = 0.
    ZERO  (I) = 0.
    FTEMP (I) = 0.
    FVAP  (I) = 0.
    RIB   (I) = 0.
    GSNOW (I) = 0.
    FTEMPX(I) = 0.
    FVAPX (I) = 0.
    RIBX  (I) = 0.
    WTABLE(I) = 9999.
    IWTABL(I) = 0
  end do ! loop 100
  !
  !     * SURFACE EVAPORATION EFFICIENCY FOR BARE SOIL ENERGY BALANCE
  !     * CALCULATIONS.
  !
  !>
  !! In loop 200 the soil surface evaporation flag IEVAP and the
  !! evaporation efficiency coefficient CEVAP are assigned. If the
  !! liquid water content of the first soil layer is effectively equal
  !! to the minimum water content THLMIN, IEVAP and CEVAP are set to
  !! zero. If the liquid water content of the first soil layer is
  !! greater than the soil porosity (THPOR), IEVAP and CEVAP are set to
  !! unity. Otherwise, IEVAP is set to 1 and CEVAP (or \f$\beta\f$ as it is
  !! typically symbolized in the literature) is calculated using a
  !! relation presented by Merlin et al. (2011) \cite Merlin2011-xy:
  !!
  !! \f$\beta = 0.25 [1 – cos(\theta_l \pi / \theta_{p})]^2\f$
  !!
  !! where \f$\theta_l\f$ is the liquid water content of the first soil layer
  !! and \f$\theta_{p}\f$ is its porosity. We follow Merlin et al. in using a
  !! 'P' value (their terminology for exponent term) of 2.
  !!
  do I = IL1,IL2 ! loop 200
    if (THLIQG(I,1) < (THLMIN(I,1) + 0.001)) then
      IEVAP(I) = 0
      CEVAP(I) = 0.0
    else if (THLIQG(I,1) > THPOR(I,1)) then
      IEVAP(I) = 1
      CEVAP(I) = 1.0
    else
      IEVAP(I) = 1
      CEVAP(I) = 0.25 * (1.0 - COS(3.14159 * THLIQG(I,1) / THPOR(I,1))) ** 2
    end if
  end do ! loop 200
  !
  !     * VOLUMETRIC HEAT CAPACITIES OF SOIL LAYERS.
  !
  !>
  !! In loop 300 the volumetric heat capacities Cg of the soil layers
  !! under a bare surface (HCPG) and under vegetation (HCPC) are
  !! calculated, from their respective liquid and frozen water
  !! contents \f$\theta_l\f$ and \f$\theta_i\f$:
  !!
  !! \f$C_g = C_w \theta_l + C_i \theta_i + C_m(1 - \theta_p)\f$
  !!
  !! where \f$C_m\f$ is the heat capacity of the soil matter and \f$\theta_p\f$ is
  !! the pore volume. (The heat capacity of air is neglected.)
  !!
  do J = 1,IG ! loop 300
    do I = IL1,IL2
      if (ISAND(I,1) > - 4) then
        HCPG(I,J) = HCPW * THLIQG(I,J) + HCPICE * THICEG(I,J) + &
                    HCPS(I,J) * (1.0 - THPOR(I,J))
        HCPC(I,J) = HCPW * THLIQC(I,J) + HCPICE * THICEC(I,J) + &
                    HCPS(I,J) * (1.0 - THPOR(I,J))
      else
        HCPC(I,J) = HCPICE
        HCPG(I,J) = HCPICE
      end if
    end do
  end do ! loop 300
  !
  !     * THERMAL PROPERTIES OF SNOW.
  !
  !>
  !! In loop 400, the thermal properties of the snow pack under the
  !! vegetation canopy and over bare soil are assigned on the basis of
  !! the properties of the snow pack over the whole modelled area. The
  !! heat capacity of the snow pack Cs is calculated from the volume
  !! fractions of snow particles and liquid water in the pack. The
  !! former is obtained from the ratio of the densities of the snow
  !! pack and ice, and the latter from the ratio of the liquid water
  !! content, normalized by the snow depth, and the density of water:
  !!
  !! \f$C_s = C_i [\rho_s /\rho_i ] + C_w w_s /[\rho_w z_s]\f$
  !!
  !! The thermal conductivity of snow \f$\lambda_s\f$ is obtained from the
  !! snow density using an empirical relationship derived by Sturm et
  !! al. (1997):
  !!
  !! \f$\lambda_s = 3.233 x 10^{-6} \rho_s^2 – 1.01 x 10^{-3} \rho_s + 0.138      \rho_s \geq 156.0 \f$
  !!
  !! \f$\lambda_s = 0.234 x 10^{-3} \rho_s + 0.023                                \rho_s < 156.0 \f$
  !!
  do I = IL1,IL2 ! loop 400
    if (ZSNOW(I) > 0.) then
      HCPSCS(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                  (RHOW * ZSNOW(I))
      HCPSGS(I) = HCPSCS(I)
      !             TCSNOW(I)=2.576E-6*RHOSNO(I)*RHOSNO(I)+0.074
      if (RHOSNO(I) < 156.0) then
        TCSNOW(I) = 0.234E-3 * RHOSNO(I) + 0.023
      else
        TCSNOW(I) = 3.233E-6 * RHOSNO(I) * RHOSNO(I) - 1.01E-3 * &
                    RHOSNO(I) + 0.138
      end if
      if (FVEG(I) < 1.) then
        TSNOGS(I) = TSNOW(I)
        WSNOGS(I) = WSNOW(I)
        RHOSGS(I) = RHOSNO(I)
      else
        TSNOGS(I) = 0.0
        WSNOGS(I) = 0.0
        RHOSGS(I) = 0.0
      end if
      if (FVEG(I) > 0.) then
        TSNOCS(I) = TSNOW(I)
        WSNOCS(I) = WSNOW(I)
        RHOSCS(I) = RHOSNO(I)
      else
        TSNOCS(I) = 0.0
        WSNOCS(I) = 0.0
        RHOSCS(I) = 0.0
      end if
    else
      TSNOGS(I) = 0.0
      WSNOGS(I) = 0.0
      RHOSGS(I) = 0.0
      TSNOCS(I) = 0.0
      WSNOCS(I) = 0.0
      RHOSCS(I) = 0.0
      TCSNOW(I) = 0.0
    end if
  end do ! loop 400
  !
  !     * THERMAL CONDUCTIVITIES OF SOIL LAYERS AND DEPTH OF WATER
  !     * TABLE IN ORGANIC SOILS.
  !
  !>
  !! In loop 500, the thermal conductivities of the soil layers are
  !! assigned. If the ISAND flag for the first soil layer is -4
  !! (indicating glacier or ice sheet), or if the ISAND flag is -3
  !! (indicating rock), then literature values for glacier ice or sand
  !! particles respectively are assigned. If the ISAND flag is equal
  !! to -2, indicating organic soil, the depth of the water table
  !! \f$z_{wt}\f$ is first calculated. This is taken to lie within the first
  !! layer, counting from the bottom of the soil profile, in which the
  !! soil water content is larger than the retention capacity \f$\theta_{ret}\f$.
  !! The water table depth is deduced by assuming that the soil is
  !! saturated below the water table, and that the water content is at
  !! the retention capacity above it. Thus, if \f$\theta_l + \theta_i = \theta_p\f$
  !! for the soil layer, the water table is located at the top of the
  !! soil layer; if \f$\theta_l + \theta_i = \theta_{ret}\f$, it is located at the
  !! permeable bottom of the soil layer; and if \f$\theta_l + \theta_i\f$ is
  !! between these two values, its location is given by:
  !!
  !! \f$z_{wt} = z_{b, w} - \Delta z_w [(\theta_l + \theta_i - \theta_{ret})/(\theta_p - \theta_{ret})]\f$
  !!
  !! where \f$\Delta z_w\f$ is the permeable thickness of the soil layer.
  !!
  !! The thermal conductivities of organic and mineral soils are
  !! calculated following the analysis of Côté and Konrad (2005).
  !! They model the soil thermal conductivity \f$\lambda\f$ using the concept of a
  !! relative thermal conductivity \f$\lambda_r\f$ which has a value of 0 for
  !! dry soils and 1 at saturation:
  !!
  !! \f$\lambda = [ \lambda_{sat} – \lambda_{dry} ] \lambda_r + \lambda_{dry}\f$
  !!
  !! The relative thermal conductivity is obtained from the degree of
  !! saturation (the water content divided by the pore volume) \f$S_r\f$,
  !! using the following generalized relationship:
  !!
  !! \f$\lambda_r = \kappa S_r/[1 + (\kappa-1) S_r ]\f$
  !!
  !! The empirical coefficient kappa takes the following values:
  !!
  !! Unfrozen coarse mineral soils:   \f$\kappa = 4.0\f$  \n
  !! Frozen coarse mineral soils:     \f$\kappa = 1.2\f$  \n
  !! Unfrozen fine mineral soils:     \f$\kappa = 1.9\f$  \n
  !! Frozen fine mineral soils:       \f$\kappa = 0.85\f$ \n
  !! Unfrozen organic soils:          \f$\kappa = 0.6\f$  \n
  !! Frozen organic soils:            \f$\kappa = 0.25\f$
  !!
  !! The dry thermal conductivity \f$lambda_{dry}\f$ is calculated using an
  !! empirical relationship based on the pore volume \f$\theta_p\f$, with
  !! different coefficients for mineral and organic soils:
  !!
  !! \f$\lambda_{dry} = 0.75 exp(-2.76 \theta_p)   (mineral)\f$
  !! \f$\lambda_{dry} = 0.30 exp(-2.0 \theta_p)    (organic)\f$
  !!
  !! The saturated thermal conductivity \f$\lambda_{sat}\f$ is calculated by
  !! Cote and Konrad (2005) \cite Cote2005-ew as a geometric mean of the conductivities of the
  !! soil components. However, other researchers (e.g. Zhang et al.,
  !! 2008) have found the linear averaging used by de Vries (1963) \cite Vries1963-ti to
  !! be more generally accurate:
  !!
  !! \f$lambda_{sat} = lambda_w \theta_p + \lambda_s (1 - \theta_p)   (unfrozen)\f$
  !! \f$lambda_{sat} = lambda_i \theta_p + \lambda_s (1 - \theta_p)   (frozen)\f$
  !!
  !! where \f$\lambda_w\f$ is the thermal conductivity of water, \f$\lambda_i\f$ is
  !! that of ice and \f$\lambda_s\f$ is that of the soil particles.
  !!
  !! In the 500 loop, thermal conductivities are calculated for the
  !! top and bottom of each soil layer. The degree of saturation
  !! SATRAT is calculated as the sum of liquid and frozen water
  !! contents, \f$\theta_w\f$ and \f$\theta_i\f$, divided by the pore volume. In
  !! organic soils, if the liquid water content of the soil layer is
  !! above the retention capacity \f$\theta_{ret}\f$, \f$\theta_w\f$ at the top of the
  !! soil layer is assumed to be equal to \f$\theta_{re}\f$ and \f$S_r\f$ at the
  !! bottom of the layer is assumed to be 1. The relative liquid and
  !! frozen water contents, THLSAT and THISAT, are calculated from
  !! \f$\theta_w\f$ and \f$\theta_i\f$ normalized by \f$\theta_w + \theta_i\f$. The dry thermal
  !! conductivity, and the saturated thermal conductivity for unfrozen
  !! and frozen conditions, are evaluated using the equations above.
  !! The unfrozen and frozen relative thermal conductivity, TCRATU and
  !! TCRATF, are obtained from SATRAT and the appropriate values of
  !! the empirical coefficient \f$\kappa\f$. For mineral soils, \f$\kappa\f$ is
  !! obtained as a weighted average over the percent sand content
  !! (ISAND converted to a real :: value) and the percentage of fine
  !! material (assumed to be 100-ISAND). The unfrozen and frozen soil
  !! thermal conductivities, TCSOLU and TCSOLF, are then calculated
  !! from TCRATU, TCRATF, and the dry and saturated thermal
  !! conductivities; and the actual thermal conductivity of the soil,
  !! TCSOIL, is determined as the average of TCSOLU and TCSOLF,
  !! weighted according to the relative liquid and frozen water
  !! contents THLSAT and THISAT. If the permeable thickness of the
  !! layer, DELZW, is greater than zero, the thermal conductivity at
  !! the top of the layer is set to TCSOIL; otherwise it is set to the
  !! rock value, TCSAND. If DELZW is less than the thermal thickness
  !! of the layer DELZ, the thermal conductivity at the bottom of the
  !! layer is set to TCSAND; otherwise it is set to TCSOIL. (In the
  !! case of organic soils in the latter branch, if \f$\theta_w\f$ was
  !! greater than \f$\theta_{ret}\f$, the thermal conductivity at the bottom of
  !! the layer is set to the average of the saturated unfrozen and
  !! frozen values, weighted by THLSAT and THISAT.) Finally, if there
  !! is ponded water present on the soil surface, the thermal
  !! conductivity at the top of the first soil layer is treated as
  !! varying linearly from the calculated soil thermal conductivity
  !! if the pond depth ZPOND is zero, to the thermal conductivity of
  !! water if ZPOND \f$\geq 10^{-2} m\f$.
  !!
  do J = IG,1, - 1 ! loop 500
    do I = IL1,IL2
      if    (ISAND(I,1) == - 4) then
        TCTOPG(I,J) = TCGLAC
        TCBOTG(I,J) = TCGLAC
      else if (ISAND(I,J) == - 3) then
        TCTOPC(I,J) = TCSAND
        TCTOPG(I,J) = TCSAND
        TCBOTC(I,J) = TCSAND
        TCBOTG(I,J) = TCSAND
      else if (ISAND(I,J) == - 2) then
        !             FLAG - Needs to be reviewed.   EC Feb 032017.
        !                    For some peatland cases, thliqg+thiceg > thpor-0.01 at bottom layer
        !                    and wtable remains 9999, causing wrong values of socres_peat
        !                    to be computed in hetresPeat, leading to crash in ctem loop 1020
        !                    during calculation of peatdep.
        !                    Testing shows that commenting out the following IF condition
        !                    produces the same WTABLE as YW's original code, except that the
        !                    crash is avoided.
        !                    Removing the check for WTABLE<9000 also doesn't change the result.
        !                    This ensures that the water table is at least in the last layer.
        !             IF (J==IG .AND. (THLIQG(I,J)+THICEG(I,J))<
        !    1            (THPOR(I,J)-0.01)) IWTABL(I)=1
        if (IWTABL(I) == 0) then ! YW
          if ((THLIQG(I,J) + THICEG(I,J)) > (THPOR(I,J) - 0.01)) then
            WTABLE(I) = ZBOTW(I,J) - DELZW(I,J)
            ! ELSEIF (WTABLE(I)<9000.0)             THEN
          else
            WTABLE(I) = ZBOTW(I,J) - DELZW(I,J) * MIN(1.0, &
                        (THLIQG(I,J) + THICEG(I,J) - THLRET(I,J)) / &
                        (THPOR(I,J) - THLRET(I,J)))
            IWTABL(I) = 1
          end if
        end if
        if (THLIQG(I,J) > (THLRET(I,J) + 0.0001)) then
          SATRAT = MIN((THLRET(I,J) + THICEG(I,J)) / &
                   THPOR(I,J),1.0)
          THLSAT = THLIQG(I,J) / (THLIQG(I,J) + THICEG(I,J))
          THISAT = THICEG(I,J) / (THLIQG(I,J) + THICEG(I,J))
          TCDRY = 0.30 * EXP( - 2.0 * THPOR(I,J))
          TCSATU(I) = TCW * THPOR(I,J) + TCS(I,J) * (1.0 - THPOR(I,J))
          TCSATF(I) = TCICE * THPOR(I,J) + TCS(I,J) * (1.0 - THPOR(I,J))
          TCRATU = 0.6 * SATRAT / (1.0 - 0.4 * SATRAT)
          TCRATF = 0.25 * SATRAT / (1.0 - 0.75 * SATRAT)
          TCSOLU = (TCSATU(I) - TCDRY) * TCRATU + TCDRY
          TCSOLF = (TCSATF(I) - TCDRY) * TCRATF + TCDRY
          TCSOIL = TCSOLU * THLSAT + TCSOLF * THISAT
          if (DELZW(I,J) > 0.0) then
            TCTOPC(I,J) = TCSOIL
            TCTOPG(I,J) = TCSOIL
          else
            TCTOPC(I,J) = TCSAND
            TCTOPG(I,J) = TCSAND
          end if
          if (DELZW(I,J) < (DELZ(J) - 0.01)) then
            TCBOTC(I,J) = TCSAND
            TCBOTG(I,J) = TCSAND
          else
            TCBOTC(I,J) = TCSATU(I) * THLSAT + TCSATF(I) * THISAT
            TCBOTG(I,J) = TCSATU(I) * THLSAT + TCSATF(I) * THISAT
          end if
          if (J == 1) then
            TCTOPC(I,J) = TCTOPC(I,J) + (TCW - TCTOPC(I,J)) * &
                          MIN(ZPOND(I),1.0E-2) * 100.0
            TCTOPG(I,J) = TCTOPC(I,J)
          end if
        else
          SATRAT = MIN((THLIQG(I,J) + THICEG(I,J)) / &
                   THPOR(I,J),1.0)
          THLSAT = THLIQG(I,J) / (THLIQG(I,J) + THICEG(I,J))
          THISAT = THICEG(I,J) / (THLIQG(I,J) + THICEG(I,J))
          TCDRY = 0.30 * EXP( - 2.0 * THPOR(I,J))
          TCSATU(I) = TCW * THPOR(I,J) + TCS(I,J) * (1.0 - THPOR(I,J))
          TCSATF(I) = TCICE * THPOR(I,J) + TCS(I,J) * (1.0 - THPOR(I,J))
          TCRATU = 0.6 * SATRAT / (1.0 - 0.4 * SATRAT)
          TCRATF = 0.25 * SATRAT / (1.0 - 0.75 * SATRAT)
          TCSOLU = (TCSATU(I) - TCDRY) * TCRATU + TCDRY
          TCSOLF = (TCSATF(I) - TCDRY) * TCRATF + TCDRY
          TCSOIL = TCSOLU * THLSAT + TCSOLF * THISAT
          if (DELZW(I,J) > 0.0) then
            TCTOPC(I,J) = TCSOIL
            TCTOPG(I,J) = TCSOIL
          else
            TCTOPC(I,J) = TCSAND
            TCTOPG(I,J) = TCSAND
          end if
          if (DELZW(I,J) < (DELZ(J) - 0.01)) then
            TCBOTC(I,J) = TCSAND
            TCBOTG(I,J) = TCSAND
          else
            TCBOTC(I,J) = TCSOIL
            TCBOTG(I,J) = TCSOIL
          end if
          if (J == 1) then
            TCTOPC(I,J) = TCTOPC(I,J) + (TCW - TCTOPC(I,J)) * &
                          MIN(ZPOND(I),1.0E-2) * 100.0
            TCTOPG(I,J) = TCTOPC(I,J)
          end if
        end if
      else

        SATRAT = MIN((THLIQG(I,J) + THICEG(I,J)) / &
                 THPOR(I,J),1.0)
        THLSAT = THLIQG(I,J) / (THLIQG(I,J) + THICEG(I,J))
        THISAT = THICEG(I,J) / (THLIQG(I,J) + THICEG(I,J))
        TCDRY = 0.75 * EXP( - 2.76 * THPOR(I,J))
        TCSATU(I) = TCW * THPOR(I,J) + TCS(I,J) * (1.0 - THPOR(I,J))
        TCSATF(I) = TCICE * THPOR(I,J) + TCS(I,J) * (1.0 - THPOR(I,J))
        TCKAPU = (4.0 * real(ISAND(I,J)) + 1.9 * real(100 - ISAND(I,J))) / &
                 100.0
        TCKAPF = (1.2 * real(ISAND(I,J)) + 0.85 * real(100 - ISAND(I,J))) / &
                 100.0
        TCRATU = TCKAPU * SATRAT / (1.0 + (TCKAPU - 1.0) * SATRAT)
        TCRATF = TCKAPF * SATRAT / (1.0 + (TCKAPF - 1.0) * SATRAT)
        TCSOLU = (TCSATU(I) - TCDRY) * TCRATU + TCDRY
        TCSOLF = (TCSATF(I) - TCDRY) * TCRATF + TCDRY
        TCSOIL = TCSOLU * THLSAT + TCSOLF * THISAT
        if (DELZW(I,J) > 0.0) then
          TCTOPC(I,J) = TCSOIL
          TCTOPG(I,J) = TCSOIL
          !                  IF (J==1) TCTOPC(I,J)=TCTOPC(I,J)*0.1
        else
          TCTOPC(I,J) = TCSAND
          TCTOPG(I,J) = TCSAND
        end if
        if (DELZW(I,J) < (DELZ(J) - 0.01)) then
          TCBOTC(I,J) = TCSAND
          TCBOTG(I,J) = TCSAND
        else
          TCBOTC(I,J) = TCSOIL
          TCBOTG(I,J) = TCSOIL
        end if
        if (J == 1) then
          TCTOPC(I,J) = TCTOPC(I,J) + (TCW - TCTOPC(I,J)) * &
                        MIN(ZPOND(I),1.0E-2) * 100.0
          TCTOPG(I,J) = TCTOPC(I,J)
        end if
      end if
    end do
  end do ! loop 500
  !
  !     * ADD PONDED WATER TEMPERATURE TO FIRST SOIL LAYER FOR USE
  !     * IN GROUND HEAT FLUX CALCULATIONS.
  !
  !>
  !! Finally, in loop 600, a variable TBAR1P is evaluated,
  !! representing the weighted average value of the first layer soil
  !! temperature and the ponded water, if any. (The heat capacity of
  !! the soil is determined as the weighted average of HCPG over the
  !! permeable thickness DELZW, and the heat capacity of rock, HCPSND,
  !! over the impermeable thickness, DELZ-DELZW.)
  !!
  do I = IL1,IL2
    if (ZPOND(I) > 0.) then
      TBAR1P(I) = (TPOND(I) * HCPW * ZPOND(I) + &
                  TBAR(I,1) * (HCPG(I,1) * DELZW(I,1) + &
                  HCPSND * (DELZ(1) - DELZW(I,1)))) / &
                  (HCPW * ZPOND(I) + HCPG(I,1) * DELZW(I,1) + &
                  HCPSND * (DELZ(1) - DELZW(I,1)))
    else
      TBAR1P(I) = TBAR(I,1)
    end if
  end do ! loop 600
  !
  ! 6990 format(I3,F6.2)
  return
end subroutine energyBudgetPrep
