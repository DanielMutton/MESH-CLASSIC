!> \file
!> Calls subroutines to perform surface energy budget calculations.
!! @author D. Verseghy, M. Lazare, K. Abdella, P. Bartlett, M. Namazi, M. MacDonald,
!! E. Chan, J. Melton, A. Wu, JP Paquin, L Duarte, Y. Delage
!
subroutine energyBudgetDriver (TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, & ! Formerly CLASST
                               THICEC, THICEG, HCPC, HCPG, TCTOPC, TCBOTC, TCTOPG, TCBOTG, &
                               GZEROC, GZEROG, GZROCS, GZROGS, G12C, G12G, G12CS, G12GS, &
                               G23C, G23G, G23CS, G23GS, QFREZC, QFREZG, QMELTC, QMELTG, &
                               EVAPC, EVAPCG, EVAPG, EVAPCS, EVPCSG, EVAPGS, TCANO, TCANS, &
                               RAICAN, SNOCAN, RAICNS, SNOCNS, CHCAP, CHCAPS, TPONDC, TPONDG, &
                               TPNDCS, TPNDGS, TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                               ITERCT, CDH, CDM, QSENS, TFLUX, QEVAP, EVAP, &
                               EVPPOT, ACOND, EVAPB, GT, QG, &
                               ST, SU, SV, SQ, SRH, &
                               GTBS, SFCUBS, SFCVBS, USTARBS, &
                               FSGV, FSGS, FSGG, FLGV, FLGS, FLGG, &
                               HFSC, HFSS, HFSG, HEVC, HEVS, HEVG, HMFC, HMFN, &
                               HTCC, HTCS, HTC, QFCF, QFCL, DRAG, WTABLE, ILMO, &
                               UE, HBL, TAC, QAC, ZREFM, ZREFH, ZDIAGM, ZDIAGH, &
                               VPD, TADP, RHOAIR, QSWINV, QSWINI, QLWIN, UWIND, VWIND, &
                               TA, QA, PADRY, FC, FG, FCS, FGS, RBCOEF, &
                               FSVF, FSVFS, PRESSG, VMOD, ALVSCN, ALIRCN, ALVSG, ALIRG, &
                               ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, ALIRSC, &
                               TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC, RCS, WTRG, groundHeatFlux, QLWAVG, &
                               FRAINC, FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP, DISPS, &
                               ZOMLNC, ZOELNC, ZOMLNG, ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS, &
                               TBAR, THLIQ, THICE, TPOND, ZPOND, TBASE, TCAN, TSNOW, &
                               ZSNOW, RHOSNO, WSNOW, THPOR, THLRET, THLMIN, THFC, THLW, &
                               TRSNOWC, TRSNOWG, ALSNO, FSSB, FROOT, FROOTS, &
                               RADJ, PCPR, HCPS, TCS, TSFSAV, DELZ, DELZW, ZBOTW, &
                               FTEMP, FVAP, RIB, &
                               ISAND, &
                               AILCG, AILCGS, FCANC, FCANCS, &
                               CO2CONC, CO2I1CG, CO2I1CS, CO2I2CG, &
                               CO2I2CS, COSZS, XDIFFUS, SLAI, &
                               ICTEM, ctem_on, RMATCTEM, FCANCMX, &
                               L2MAX, NOL2PFTS, CFLUXCG, CFLUXCS, &
                               ANCSVEG, ANCGVEG, RMLCSVEG, RMLCGVEG, &
                               TCSNOW, GSNOW, &
                               ITC, ITCG, ITG, ILG, IL1, IL2, JL, N, IC, &
                               IG, IZREF, ISLFD, NLANDCS, NLANDGS, NLANDC, NLANDG, NLANDI, &
                               NBS, ISNOALB, DAYL, DAYL_MAX, &
                               ipeatland, ancsmoss, angsmoss, ancmoss, angmoss, &
                               rmlcsmoss, rmlgsmoss, rmlcmoss, rmlgmoss, Cmossmas, dmoss, &
                               iday, pdd)

  !
  !     * OCT 26/16 - D.VERSEGHY. ADD ZPOND TO energBalNoVegSolve CALLS.
  !     * AUG 30/16 - J.Melton    Replace ICTEMMOD with ctem_on (logical switch).
  !     * AUG 04/15 - M.LAZARE.   SPLIT FROOT INTO TWO ARRAYS, FOR CANOPY
  !     *                         AREAS WITH AND WITHOUT SNOW.
  !     * JUL 22/15 - D.VERSEGHY. CHANGES TO energBalVegSolve AND energBalNoVegSolve CALLS.
  !     * FEB 27/15 - J. MELTON - WILTSM AND FIELDSM ARE RENAMED THLW AND THFC, RESPECTIVELY.
  !     * FEB 09/15 - D.VERSEGHY. New version for gcm18 and class 3.6:
  !     *                         - Revised calls to revised energyBudgetPrep for
  !     *                           initialization of SRH and SLDIAG.
  !     *                         - Input {THFC, THLW} (from soilProperties) replace
  !     *                           work arrays {FIELDSM, WILTSM}.
  !     *                         - Calculation of new bare-soil fields
  !     *                           {GTBS, SFCUBS, SFCVBS, USTARBS}.
  !     * SEP 09/14 - D.VERSEGHY/M.LAZARE. CORRECTIONS TO SCREEN LEVEL
  !     *                         DIAGNOSTIC CALCULATIONS.
  !     * AUG 19/13 - M.LAZARE.   REMOVE CALCULATION AND REFERENCES TO
  !     *                         "QFLUX" (NOW DONE IN waterBudgetDriver).
  !     * JUN 21/13 - M.LAZARE.   REVISED CALL TO energyBudgetPrep TO SUPPORT ADDING
  !     *                         INITIALIZATION OF "GSNOW".
  !     * JUN 10/13 - M.LAZARE/   ADD SUPPORT FOR "ISNOALB" FORMULATION.
  !     *             M.NAMAZI.
  !
  !     * NOV 11/11 - M.LAZARE.   IMPLEMENT CTEM (INITIALIZATION OF FIELDS
  !     *                         NEAR BEGINNING AND TWO REVISED CALLS TO
  !     *                         energBalVegSolve).
  !     * OCT 12/11 - M.LAZARE.   REMOVED "TSURF" (REQUIRED CHANGE
  !     *                         TO energyBudgetPrep INITIALIZATION AS WELL).
  !     * OCT 07/11 - M.LAZARE.   - CHANGE QLWAVG FROM AN INTERNAL WORK
  !     *                           ARRAY TO ONE PASSED OUT TO THE CLASS
  !     *                           DRIVER, TO ACCOMODATE RPN.
  !     *                         - WIND SPEED NOW PASSED IN (POSSIBLY
  !     *                           CONTAINING GUSTINESS FACTOR) AS "VMOD",
  !     *                           INSTEAD OF CALCULATING IT LOCALLY.
  !     * OCT 05/11 - M.LAZARE.   ADD CALCULATION OF SRH (REQUIRES PASSING
  !     *                         IN OF PRESSG PLUS ADDITIONAL INTERNAL
  !     *                         WORK ARRAYS).
  !     * APR 28/10 - D.VERSEGHY. REVISE CALCULATION OF QG.
  !     * APR 28/10 - M.MACDONALD/D.VERSEGHY. CORRECT CALCULATIONS OF
  !     *                         CRIB, DRAG AND VAC FOR ISLFD=1.
  !     * APR 28/10 - E.CHAN/D/VERSEGHY. CORRECT CALCULATIONS OF ST AND
  !     *                         SQ FOR ISLFD=0.
  !     * DEC 21/09 - D.VERSEGHY. CORRECT BUG IN CALL TO energBalVegSolve IN CS
  !     *                         SUBAREA (CALL WITH FSNOCS AND RAICNS).
  !     * DEC 07/09 - D.VERSEGHY. ADD EVAP TO energBalVegSolve CALL.
  !     * JAN 06/09 - D.VERSEGHY. INSERT UPDATES TO HTC AND WTRG
  !     *                         BRACKETTING LOOP 60; CORRECT TPOTA,
  !     *                         ZRSLDM AND ZRSLDH CALCULATIONS; USE
  !     *                         TPOTA IN SLDIAG CALL; ASSUME THAT TA IS
  !     *                         ADIABATIC EXTRAPOLATE TO SURFACE FOR
  !     *                         ATMOSPHERIC MODELS.
  !     * NOV 03/08 - L.DUARTE    CORRECTED CALL TO energBalVegSolve.
  !     * AUG    08 - JP PAQUIN   OUTPUT FTEMP, FVAP AND RIB FOR GEM
  !                               (IMPLEMENTED BY L.DUARTE ON OCT 28/08)
  !     * FEB 25/08 - D.VERSEGHY. MODIFICATIONS REFLECTING CHANGES
  !     *                         ELSEWHERE IN CODE.
  !     * NOV 24/06 - D.VERSEGHY. REMOVE CALL TO TZTHRM; MAKE RADJ REAL.
  !     * AUG 16/06 - D.VERSEGHY. NEW CALLS TO snowHeatCond AND snowTempUpdate.
  !     * APR 13/06 - D.VERSEGHY. SEPARATE GROUND AND SNOW ALBEDOS FOR
  !     *                         OPEN AND CANOPY-COVERED AREAS.
  !     * MAR 23/06 - D.VERSEGHY. CHANGES TO ADD MODELLING OF WSNOW.
  !     * MAR 21/06 - P.BARTLETT. PASS ADDITIONAL VARIABLES TO energyBudgetPrep.
  !     * OCT 04/05 - D.VERSEGHY. MODIFICATIONS TO ALLOW OPTION OF SUB-
  !     *                         DIVIDING THIRD SOIL LAYER.
  !     * APR 12/05 - D.VERSEGHY. VARIOUS NEW FIELDS; ADD CALL TO NEW
  !     *                         SUBROUTINE TZTHRM; MOVE CALCULATION
  !     *                         OF CPHCHC INTO energBalVegSolve.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 05/04 - Y.DELAGE/D.VERSEGHY. NEW DIAGNOSTIC VARIABLES
  !     *                         ILMO, UE AND HBL.
  !     * NOV 07/02 - Y.DELAGE/D.VERSEGHY. CALLS TO NEW DIAGNOSTIC
  !     *                         SUBROUTINES "SLDIAG" AND "DIASURF"
  !     *                         MODIFICATIONS TO ACCOMMODATE DIFFERENT
  !     *                         SURFACE REFERENCE HEIGHT CONVENTIONS.
  !     * JUL 31/02 - D.VERSEGHY. MOVE CALCULATION OF VEGETATION STOMATAL
  !     *                         RESISTANCE FROM energyBudgetPrep INTO calcLandSurfParams AND
  !     *                         canopyAlbedoTransmiss; SHORTENED CLASS3 COMMON BLOCK.
  !     * JUL 23/02 - D.VERSEGHY. MOVE ADDITION OF AIR TO CANOPY MASS
  !     *                         INTO radiationDriver; SHORTENED CLASS4
  !     *                         COMMON BLOCK.
  !     * MAR 28/02 - D.VERSEGHY. STREAMLINED SUBROUTINE CALLS.
  !     * MAR 22/02 - D.VERSEGHY. MOVE CALCULATION OF BACKGROUND SOIL
  !     *                         PROPERTIES INTO "soilProperties"; ADD NEW
  !     *                         DIAGNOSTIC VARIABLES "EVPPOT", "ACOND"
  !     *                         AND "TSURF"; MODIFY CALCULATIONS OF VAC,
  !     *                         EVAPB AND QG.
  !     * JAN 18/02 - D.VERSEGHY. CHANGES TO INCORPORATE NEW BARE SOIL
  !     *                         EVAPORATION FORMULATION.
  !     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
  !     * SEP 19/00 - D.VERSEGHY. PASS VEGETATION-VARYING COEFFICIENTS
  !     *                         TO energyBudgetPrep FOR CALCULATION OF STOMATAL
  !     *                         RESISTANCE.
  !     * DEC 16/99 - A.WU/D.VERSEGHY. CHANGES MADE TO INCORPORATE NEW SOIL
  !     *                              EVAPORATION ALGORITHMS AND NEW CANOPY
  !     *                              TURBULENT FLUX FORMULATION.  MODIFY
  !     *                              CALCULATION OF BULK RICHARDSON NUMBER
  !     *                              AND CANOPY MASS.
  !     * APR 15/99 - M.LAZARE.   CORRECT SCREEN-LEVEL CALCULATION FOR WINDS
  !     *                         TO HOLD AT ANEMOMETER LEVEL (10M) INSTEAD
  !     *                         OF SCREEN LEVEL (2M).
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         CHANGES RELATED TO VARIABLE SOIL DEPTH
  !     *                         (MOISTURE HOLDING CAPACITY) AND DEPTH-
  !     *                         VARYING SOIL PROPERTIES.
  !     *                         ALSO, APPLY UPPER BOUND ON "RATFC1").
  !     * OCT 11/96 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         REVISE CALCULATION OF SLTHKEF AND
  !     *                         DEFINITION OF ZREF FOR INTERNAL
  !     *                         CONSISTENCY.
  !     * SEP 27/96 - D.VERSEGHY. FIX BUG IN CALCULATION OF FLUXES
  !     *                         BETWEEN SOIL LAYERS (PRESENT SINCE
  !     *                         RELEASE OF CLASS VERSION 2.5).
  !     * MAY 21/96 - K.ABDELLA.  CORRECT EXPRESSION FOR ZOSCLH (4 PLACES).
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS; ALSO, PASS IN ZREF AND
  !     *                         ILW THROUGH SUBROUTINE CALL.
  !     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
  !     *                         ORGANIC MATTER CONTENT.
  !     * DEC 16/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         ADD THREE NEW DIAGNOSTIC FIELDS
  !     *                         REVISE CALCULATION OF HTCS, HTC.
  !     * DEC 06/94 - M.LAZARE. - PASS "CFLUX" TO energBalNoVegSolve INSTEAD OF
  !     *                         "CLIMIT" IN CONJUNCTION WITH CHANGES
  !     *                         TO THAT ROUTINE.
  !     *                       - REVISE CALCULATION OF "ZREF" TO INCLUDE
  !     *                         VIRTUAL TEMPERATURE EFFECTS.
  !     *                       - REVISE CALCULATION OF "SLTHKEF".
  !     * NOV 28/94 - M.LAZARE.   FORM DRAG "CDOM" MODIFICATION REMOVED.
  !     * NOV 18/93 - D.VERSEGHY. CLASS - VERSION 2.2.
  !     *                         LOCAL VERSION WITH INTERNAL WORK ARRAYS
  !     *                         HARD-CODED FOR USE ON PCS.
  !     * NOV 05/93 - M.LAZARE.   ADD NEW DIAGNOSTIC OUTPUT FIELD: DRAG.
  !     * JUL 27/93 - D.VERSEGHY/M.LAZARE. PREVIOUS VERSION energyBudgetDriverO.
  !
  use classicParams,      only : DELT, TFREZ, GRAV, SBC, VKC, VMIN, HCPICE, &
                                 HCPW, SPHAIR, RHOW, RHOICE
  use generalutilsGAT,  only : calcEsat
  
  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(inout) :: NLANDCS !< Number of modelled areas that contain subareas of canopy over
  !< bare ground / bare ground / canopy over snow /snow
  integer, intent(inout) :: NLANDGS !< Number of modelled areas that contain subareas of canopy over
  !< bare ground / bare ground / canopy over snow /snow
  integer, intent(inout) :: NLANDC !< Number of modelled areas that contain subareas of canopy over
  !< bare ground / bare ground / canopy over snow /snow
  integer, intent(inout) :: NLANDG !< Number of modelled areas that contain subareas of canopy over
  !< bare ground / bare ground / canopy over snow /snow
  integer, intent(inout) :: NLANDI !< Number of modelled areas that are ice sheets [ ]
  integer, intent(in) :: N      !<
  integer, intent(in) :: NBS    !<
  integer, intent(in) :: ISNOALB !< Switch to model snow albedo in two or more wavelength bands
  !
  integer, intent(in) :: ITC   !< Flag to select iteration scheme for canopy temperature
  integer, intent(in) :: ITCG  !< Flag to select iteration scheme for surface under canopy
  integer, intent(in) :: ITG   !< Flag to select iteration scheme for ground or snow surface
  integer, intent(in) :: ILG   !<
  integer, intent(in) :: IL1   !<
  integer, intent(in) :: IL2   !<
  integer, intent(in) :: JL    !<
  integer, intent(in) :: IC    !<
  integer, intent(in) :: IG    !<
  integer, intent(in) :: IZREF !< Flag governing treatment of surface roughness length
  integer, intent(in) :: ISLFD !< Flag governing options for surface stability functions and diagnostic calculations
  integer :: I, J, ISNOW
  !
  !     * OUTPUT FIELDS.
  !
  real, intent(out) :: ACOND (ILG) !< Diagnosed product of drag coefficient and wind speed over modelled area \f$[m s^{-1} ]\f$
  real, intent(out) :: CDH   (ILG) !< Surface drag coefficient for heat \f$[ ] (C_{DH} )\f$
  real, intent(out) :: CDM   (ILG) !< Surface drag coefficient for momentum \f$[ ] (C_{DM} )\f$
  real, intent(out) :: CHCAP (ILG) !< Surface drag coefficient for momentum \f$[ ] (C_{DM} )\f$
  real, intent(out) :: CHCAPS(ILG) !< Heat capacity of canopy over bare ground \f$[J m^{-2} K^{-1} ]\f$
  real, intent(out) :: DRAG  (ILG) !< Surface drag coefficient under neutral stability [ ]

  real, intent(out) :: EVAP  (ILG) !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(out) :: EVAPB (ILG) !< Evaporation efficiency at ground surface [ ]
  real, intent(out) :: EVAPC (ILG) !< Evaporation from vegetation over ground \f$[m s^{-1} ]\f$
  real, intent(out) :: EVAPCG(ILG) !< Evaporation from ground under vegetation \f$[m s^{-1} ]\f$
  real, intent(out) :: EVAPCS(ILG) !< Evaporation from vegetation over snow \f$[m s^{-1} ]\f$
  real, intent(out) :: EVAPG (ILG) !< Evaporation from bare ground \f$[m s^{-1} ]\f$
  real, intent(out) :: EVAPGS(ILG) !< Evaporation from snow on bare ground \f$[m s^{-1} ]\f$
  real, intent(out) :: EVPCSG(ILG) !< Evaporation from snow under vegetation \f$[m s^{-1} ]\f$
  real, intent(out) :: EVPPOT(ILG) !< Diagnosed potential evapotranspiration \f$[kg m^{-2} s^{-1} ] (E_p)\f$
  real, intent(out) :: FLGG  (ILG) !< Diagnosed net longwave radiation at soil surface \f$[W m^{-2} ]\f$
  real, intent(out) :: FLGS  (ILG) !< Diagnosed net longwave radiation at snow surface \f$[W m^{-2} ]\f$
  real, intent(out) :: FLGV  (ILG) !< Diagnosed net longwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
  real, intent(out) :: FSGG  (ILG) !< Diagnosed net shortwave radiation at soil surface \f$[W m^{-2} ]\f$
  real, intent(out) :: FSGS  (ILG) !< Diagnosed net shortwave radiation at snow surface \f$[W m^{-2} ]\f$
  real, intent(out) :: FSGV  (ILG) !< Diagnosed net shortwave radiation on vegetation canopy \f$[W m^{-2} ]\f$
  real, intent(out) :: FTEMP (ILG) !<
  real, intent(out) :: FVAP  (ILG) !<

  real, intent(out) :: G12C  (ILG) !< Subarea heat flux between first and second soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: G12CS (ILG) !< Subarea heat flux between first and second soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: G12G  (ILG) !< Subarea heat flux between first and second soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: G12GS (ILG) !< Subarea heat flux between first and second soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: G23C  (ILG) !< Subarea heat flux between second and third soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: G23CS (ILG) !< Subarea heat flux between second and third soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: G23G  (ILG) !< Subarea heat flux between second and third soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: G23GS (ILG) !< Subarea heat flux between second and third soil layers \f$[W m^{-2} ]\f$
  real, intent(out) :: GT    (ILG) !< Diagnosed effective surface black-body temperature \f$[K] (T_{0, eff} )\f$
  real, intent(out) :: GTBS  (ILG) !< Surface temperature for CCCma black carbon scheme  [K]
  real, intent(out) :: GSNOW (ILG) !< Heat conduction into surface of snow pack  \f$[W m^{-2} ]\f$
  real, intent(out) :: GZEROC(ILG) !< Subarea heat flux at soil surface \f$[W m^{-2} ]\f$
  real, intent(out) :: GZEROG(ILG) !< Subarea heat flux at soil surface \f$[W m^{-2} ]\f$
  real, intent(out) :: GZROCS(ILG) !< Subarea heat flux at soil surface \f$[W m^{-2} ]\f$
  real, intent(out) :: GZROGS(ILG) !< Subarea heat flux at soil surface \f$[W m^{-2} ]\f$

  real, intent(out) :: HBL   (ILG)    !< Height of the atmospheric boundary layer [m]
  real, intent(out) :: HCPC  (ILG,IG) !< Heat capacity of soil layers under vegetation \f$[J m^{-3} K^{-1} ]\f$
  real, intent(out) :: HCPG  (ILG,IG) !< Heat capacity of soil layers in bare areas \f$[J m^{-3} K^{-1} ]\f$
  real, intent(out) :: HEVC  (ILG)    !< Diagnosed latent heat flux on vegetation canopy \f$[W m^{-2} ]\f$
  real, intent(out) :: HEVG  (ILG)    !< Diagnosed latent heat flux at soil surface \f$[W m^{-2} ]\f$
  real, intent(out) :: HEVS  (ILG)    !< Diagnosed latent heat flux at snow surface \f$[W m^{-2} ]\f$
  real, intent(out) :: HFSC  (ILG)    !< Diagnosed sensible heat flux on vegetation canopy \f$[W m^{-2} ]\f$
  real, intent(out) :: HFSG  (ILG)    !< Diagnosed sensible heat flux at soil surface \f$[W m^{-2} ]\f$
  real, intent(out) :: HFSS  (ILG)    !< Diagnosed sensible heat flux at snow surface \f$[W m^{-2} ]\f$
  real, intent(out) :: HMFC  (ILG)    !< Diagnosed energy associated with phase change of water on vegetation \f$[W m^{-2} ]\f$
  real, intent(out) :: HMFN  (ILG)    !< Diagnosed energy associated with phase change of water in snow pack \f$[W m^{-2} ]\f$
  real, intent(out) :: HTC   (ILG,IG) !< Diagnosed internal energy change of soil layer due to conduction
  !< and/or change in mass \f$[W m^{-2} ]\f$
  real, intent(out) :: HTCC  (ILG)    !< Diagnosed internal energy change of vegetation canopy due to conduction
  !< and/or change in mass \f$[W m^{-2} ]\f$
  real, intent(out) :: HTCS  (ILG)    !< Diagnosed internal energy change of snow pack due to conduction
  !< and/or change in mass \f$[W m^{-2} ]\f$

  real, intent(out) :: ILMO  (ILG)         !< Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
  integer :: ITERCT(ILG,6,50) !< Counter of number of iterations required to solve surface energy balance
  !< for the elements of the four subareas

  real, intent(out) :: QAC   (ILG) !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ]\f$
  real, intent(out) :: QEVAP (ILG) !< Diagnosed total surface latent heat flux over modelled area \f$[W m^{-2} ]\f$
  real, intent(out) :: QFCF  (ILG) !< Sublimation from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(out) :: QFCL  (ILG) !< Evaporation from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$

  real, intent(out) :: QFREZC(ILG) !< Heat sink to be used for freezing water on ground under canopy \f$[W m^{-2} ]\f$
  real, intent(out) :: QFREZG(ILG) !< Heat sink to be used for freezing water on bare ground \f$[W m^{-2} ]\f$
  real, intent(out) :: QG    (ILG) !< Diagnosed surface specific humidity \f$[kg kg^{-1} ]\f$
  real, intent(out) :: QLWAVG(ILG) !< Upwelling longwave radiation from land surface \f$[W m^{-2} ]\f$
  real, intent(out) :: QMELTC(ILG) !< Heat to be used for melting snow under canopy \f$[W m^{-2} ]\f$
  real, intent(out) :: QMELTG(ILG) !< Heat to be used for melting snow on bare ground \f$[W m^{-2} ]\f$
  real, intent(out) :: QSENS (ILG) !< Diagnosed total surface sensible heat flux over modelled area \f$[W m^{-2} ]\f$
  real, intent(out) :: RAICAN(ILG) !< Intercepted liquid water stored on canopy over ground \f$[kg m^{-2} ]\f$
  real, intent(out) :: RAICNS(ILG) !< Intercepted liquid water stored on canopy over snow \f$[kg m^{-2} ]\f$
  real, intent(out) :: RHOSCS(ILG) !< Density of snow under vegetation \f$[kg m^{-3} ]\f$
  real, intent(out) :: RHOSGS(ILG) !< Density of snow in bare areas \f$[kg m^{-3} ]\f$
  real, intent(out) :: RIB   (ILG) !<
  real, intent(out) :: SFCUBS(ILG) !< Zonal surface wind velocity for CCCma black carbon scheme  \f$[m s^{-1} ]\f$
  real, intent(out) :: SFCVBS(ILG) !< Meridional surface wind velocity for CCCma black carbon scheme  \f$[m s^{-1} ]\f$
  real, intent(out) :: SNOCAN(ILG) !< Intercepted frozen water stored on canopy over ground \f$[kg m^{-2} ]\f$
  real, intent(out) :: SNOCNS(ILG) !< Intercepted frozen water stored on canopy over snow \f$[kg m^{-2} ]\f$
  real, intent(out) :: SQ    (ILG) !< Diagnosed screen-level specific humidity \f$[kg kg^{-1} ]\f$
  real, intent(out) :: SRH   (ILG) !< Diagnosed screen-level relative humidity [%]
  real, intent(out) :: ST    (ILG) !< Diagnosed screen-level air temperature [K]
  real, intent(out) :: SU    (ILG) !< Diagnosed anemometer-level zonal wind \f$[m s^{-1} ]\f$
  real, intent(out) :: SV    (ILG) !< Diagnosed anemometer-level meridional wind \f$[m s^{-1} ]\f$

  real, intent(out) :: TAC   (ILG)    !< Temperature of air within vegetation canopy [K]
  real, intent(out) :: TBARC (ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(out) :: TBARCS(ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(out) :: TBARG (ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(out) :: TBARGS(ILG,IG) !< Subarea temperatures of soil layers [C]
  real, intent(out) :: TCANO (ILG)    !< Temperature of canopy over ground [K]
  real, intent(out) :: TCANS (ILG)    !< Temperature of canopy over snow [K]
  real, intent(out) :: TCBOTC(ILG,IG) !< Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$
  real, intent(out) :: TCBOTG(ILG,IG) !< Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1} ]\f$
  real, intent(out) :: TCTOPC(ILG,IG) !< Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$
  real, intent(out) :: TCTOPG(ILG,IG) !< Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1} ]\f$
  real, intent(out) :: TCSNOW(ILG)    !< Thermal conductivity of snow  \f$[W m^{-1} K^{-1} ]\f$
  real, intent(out) :: TFLUX (ILG)    !< Product of surface drag coefficient, wind speed and surface-air
  !< temperature difference \f$[K m s^{-1} ]\f$
  real, intent(out) :: THICEC(ILG,IG) !< Frozen water content of soil layers under vegetation \f$[m^3 m^{-3} ]\f$
  real, intent(out) :: THICEG(ILG,IG) !< Frozen water content of soil layers in bare areas \f$[m^3 m^{-3} ]\f$
  real, intent(out) :: THLIQC(ILG,IG) !< Liquid water content of soil layers under vegetation \f$[m^3 m^{-3} ]\f$
  real, intent(out) :: THLIQG(ILG,IG) !< Liquid water content of soil layers in bare areas \f$[m^3 m^{-3} ]\f$
  real, intent(out) :: TPONDC(ILG)    !< Subarea temperature of surface ponded water [C]
  real, intent(out) :: TPONDG(ILG)    !< Subarea temperature of surface ponded water [C]
  real, intent(out) :: TPNDCS(ILG)    !< Subarea temperature of surface ponded water [C]
  real, intent(out) :: TPNDGS(ILG)    !< Subarea temperature of surface ponded water [C]
  real, intent(out) :: TSFSAV(ILG,4)  !< Ground surface temperature over subarea [K]
  real, intent(out) :: TSNOCS(ILG)    !< Temperature of snow pack under vegetation [K]
  real, intent(out) :: TSNOGS(ILG)    !< Temperature of snow pack in bare areas [K]

  real, intent(out) :: UE    (ILG)  !< Friction velocity of air \f$[m s^{-1} ]\f$
  real, intent(out) :: USTARBS(ILG) !< Friction velocity for CCCma black carbon scheme  \f$[m s^{-1} ]\f$
  real, intent(out) :: WSNOCS(ILG)  !< Liquid water content of snow pack under vegetation \f$[kg m^{-2} ]\f$
  real, intent(out) :: WSNOGS(ILG)  !< Liquid water content of snow pack in bare areas \f$[kg m^{-2} ]\f$
  real, intent(out) :: WTABLE(ILG)  !< Depth of water table in soil [m]
  real, intent(out) :: WTRG  (ILG)  !< Diagnosed residual water transferred into or out of the soil \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(out) :: groundHeatFlux(ILG) !< Heat flux at soil surface \f$[W m^{-2} ]\f$
  !
  !
  !     * INPUT FIELDS.
  !
  real, intent(in) :: ZREFM (ILG) !< Reference height associated with forcing wind speed [m]
  real, intent(in) :: ZREFH (ILG) !< Reference height associated with forcing air temperature and humidity [m]
  real, intent(in) :: ZDIAGM(ILG) !< User-specified height associated with diagnosed anemometer-level wind speed [m]
  real, intent(in) :: ZDIAGH(ILG) !< User-specified height associated with diagnosed screen-level variables [m]
  real, intent(in) :: VPD   (ILG) !< Vapour pressure deficit [mb]
  real, intent(in) :: TADP  (ILG) !< Dew point temperature of air [K]
  real, intent(in) :: RHOAIR(ILG) !< Density of air \f$[kg m^{-3} ] (\rho_a)\f$
  real, intent(in) :: QSWINV(ILG) !< Visible radiation incident on horizontal surface \f$[W m^{-2} ]\f$
  real, intent(in) :: QSWINI(ILG) !< Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ]\f$
  real, intent(in) :: QLWIN (ILG) !< Downwelling longwave radiation at bottom of atmosphere \f$[W m^{-2} ]\f$
  real, intent(in) :: UWIND (ILG) !< Zonal component of wind speed \f$[m s^{-1} ] (U_a)\f$
  real, intent(in) :: VWIND (ILG) !< Meridional component of wind speed \f$[m s^{-1} ] (V_a)\f$
  real, intent(in) :: TA    (ILG) !< Air temperature at reference height \f$[K] (T_a)\f$
  real, intent(in) :: QA    (ILG) !< Specific humidity at reference height \f$[kg kg^{-1} ]\f$
  real, intent(in) :: PADRY (ILG) !< Partial pressure of dry air \f$[Pa] (p_{dry} )\f$
  real, intent(in) :: FC    (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FG    (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FCS   (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: FGS   (ILG) !< Subarea fractional coverage of modelled area [ ]
  real, intent(in) :: RBCOEF(ILG) !< Parameter for calculation of leaf boundary resistance
  real, intent(in) :: FSVF  (ILG) !< Sky view factor for bare ground under canopy [ ]
  real, intent(in) :: FSVFS (ILG) !< Sky view factor for snow under canopy [ ]
  real, intent(in) :: PRESSG(ILG) !< Surface atmospheric pressure [Pa]
  real, intent(in) :: VMOD  (ILG) !< Wind speed at reference height \f$[m s^{-1} ]\f$
  real, intent(in) :: ALSNO(ILG,NBS) !< Albedo of snow in each modelled wavelength band  [  ]
  real, intent(in) :: ALVSCN(ILG) !< Visible/near-IR albedo of vegetation over bare ground [ ]
  real, intent(in) :: ALIRCN(ILG) !< Visible/near-IR albedo of vegetation over bare ground [ ]
  real, intent(in) :: ALVSG (ILG) !< Visible/near-IR albedo of open bare ground [ ]
  real, intent(in) :: ALIRG (ILG) !< Visible/near-IR albedo of open bare ground [ ]
  real, intent(in) :: ALVSCS(ILG) !< Visible/near-IR albedo of vegetation over snow [ ]
  real, intent(in) :: ALIRCS(ILG) !< Visible/near-IR albedo of vegetation over snow [ ]
  real, intent(in) :: ALVSSN(ILG) !< Visible/near-IR albedo of open snow cover [ ]
  real, intent(in) :: ALIRSN(ILG) !< Visible/near-IR albedo of open snow cover [ ]
  real, intent(in) :: ALVSGC(ILG) !< Visible/near-IR albedo of bare ground under vegetation [ ]
  real, intent(in) :: ALIRGC(ILG) !< Visible/near-IR albedo of bare ground under vegetation [ ]
  real, intent(in) :: ALVSSC(ILG) !< Visible/near-IR albedo of snow under vegetation [ ]
  real, intent(in) :: ALIRSC(ILG) !< Visible/near-IR albedo of snow under vegetation [ ]
  real, intent(in) :: TRVSCN(ILG) !< Visible/near-IR transmissivity of vegetation over bare ground [ ]
  real, intent(in) :: TRIRCN(ILG) !< Visible/near-IR transmissivity of vegetation over bare ground [ ]
  real, intent(in) :: TRVSCS(ILG) !< Visible/near-IR transmissivity of vegetation over snow [ ]
  real, intent(in) :: TRIRCS(ILG) !< Visible/near-IR transmissivity of vegetation over snow [ ]
  real, intent(in) :: RC    (ILG) !< Stomatal resistance of vegetation over bare ground \f$[s m^{-1} ]\f$
  real, intent(in) :: RCS   (ILG) !< Stomatal resistance of vegetation over snow \f$[s m^{-1} ]\f$
  real, intent(in) :: FRAINC(ILG) !< Fractional coverage of canopy by liquid water over snow-free subarea [ ]
  real, intent(in) :: FSNOWC(ILG) !< Fractional coverage of canopy by frozen water over snow-free subarea [ ]
  real, intent(in) :: FRAICS(ILG) !< Fractional coverage of canopy by liquid water over snow-covered subarea [ ]
  real, intent(in) :: FSNOCS(ILG) !< Fractional coverage of canopy by frozen water over snow-covered subarea [ ]
  real, intent(in) :: CMASSC(ILG) !< Mass of canopy over bare ground \f$[kg m^{-2} ]\f$
  real, intent(in) :: CMASCS(ILG) !< Mass of canopy over snow \f$[kg m^{-2} ]\f$
  real, intent(in) :: DISP  (ILG) !< Displacement height of vegetation over bare ground [m] (d)
  real, intent(in) :: DISPS (ILG) !< Displacement height of vegetation over snow [m] (d)
  real, intent(in) :: ZOMLNC(ILG) !< Logarithm of roughness length for momentum of vegetation over bare ground [ ]
  real, intent(in) :: ZOELNC(ILG) !< Logarithm of roughness length for heat of vegetation over bare ground [ ]
  real, intent(in) :: ZOMLNG(ILG) !< Logarithm of roughness length for momentum of bare ground [ ]
  real, intent(in) :: ZOELNG(ILG) !< Logarithm of roughness length for heat of bare ground [ ]
  real, intent(in) :: ZOMLCS(ILG) !< Logarithm of roughness length for momentum of vegetation over snow [ ]
  real, intent(in) :: ZOELCS(ILG) !< Logarithm of roughness length for heat of vegetation over snow [ ]
  real, intent(in) :: ZOMLNS(ILG) !< Logarithm of roughness length for momentum of snow [ ]
  real, intent(in) :: ZOELNS(ILG) !< Logarithm of roughness length for heat of snow [ ]
  real, intent(in) :: TPOND (ILG) !< Temperature of ponded water [K]
  real, intent(inout) :: ZPOND (ILG) !< Depth of ponded water on surface [m]
  real, intent(in) :: TBASE (ILG) !< Temperature of bedrock in third soil layer [K]
  real, intent(in) :: TCAN  (ILG) !< Vegetation canopy temperature [K]
  real, intent(in) :: TSNOW (ILG) !< Snowpack temperature [K]
  real, intent(in) :: ZSNOW (ILG) !< Depth of snow pack [m]
  real, intent(in) :: TRSNOWC(ILG)!< Transmissivity of snow under vegetation to shortwave radiation  [  ]
  real, intent(in) :: TRSNOWG(ILG,NBS) !< Transmissivity of snow in bare areas to shortwave radiation  [  ]
  real, intent(in) :: RHOSNO(ILG) !< Density of snow \f$[kg m^{-3} ]\f$
  real, intent(in) :: WSNOW (ILG) !< Liquid water content of snow pack \f$[kg m^{-2} ]\f$
  real, intent(in) :: RADJ  (ILG) !< Latitude of grid cell (positive north of equator) [rad] \f$(\varphi)\f$
  real, intent(in) :: PCPR  (ILG) !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(in) :: FSSB(ILG,NBS)  !< Total solar radiation in each modelled wavelength band  \f$[W m^{-2} ]\f$
  real, intent(in) :: TBAR  (ILG,IG) !< Temperature of soil layers [K]
  real, intent(inout) :: THLIQ (ILG,IG) !< Volumetric liquid water content of soil layers \f$[m^3 m^{-3} ]\f$
  real, intent(inout) :: THICE (ILG,IG) !< Volumetric frozen water content of soil layers \f$[m^3 m^{-3} ]\f$
  !
  !     * SOIL PROPERTY ARRAYS.
  !
  real, intent(in) :: THPOR (ILG,IG) !< Pore volume in soil layer \f$[m^3 m^{-3} ]\f$
  real, intent(in) :: THLRET(ILG,IG) !< Liquid water retention capacity for organic soil \f$[m^3 m^{-3} ]\f$
  real, intent(in) :: THLMIN(ILG,IG) !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^3 m^{-3} ]\f$
  real, intent(in) :: THFC  (ILG,IG) !< Field capacity \f$[m^3 m^{-3} ]\f$
  real, intent(in) :: THLW  (ILG,IG) !< Soil water content at wilting point, \f$[m^3 m^{-3} ]\f$
  real, intent(in) :: HCPS  (ILG,IG) !< Heat capacity of soil material \f$[J m^{-3} K^{-1} ]\f$
  real, intent(in) :: TCS   (ILG,IG) !< Thermal conductivity of soil particles \f$[W m^{-1} K^{-1} ]\f$
  real, intent(in) :: DELZ  (IG)     !< Overall thickness of soil layer [m]
  real, intent(in) :: DELZW (ILG,IG) !< Permeable thickness of soil layer [m]
  real, intent(in) :: ZBOTW (ILG,IG) !< Depth to permeable bottom of soil layer [m]
  real, intent(in) :: FROOT (ILG,IG) !< Fraction of total transpiration contributed by soil layer over snow-free subarea  [  ]
  real, intent(in) :: FROOTS(ILG,IG) !< Fraction of total transpiration contributed by soil layer over snow-covered subarea  [  ]
  !
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  !     * CTEM-RELATED I/O FIELDS.
  !

  real, intent(in) :: AILCG(ILG,ICTEM)   !< GREEN LAI FOR USE WITH PHOTOSYNTHESIS SUBTROUTINE FOR CANOPY OVER GROUND SUBAREA
  real, intent(in) :: AILCGS(ILG,ICTEM)  !< GREEN LAI FOR USE WITH PHOTOSYNTHESIS SUBTROUTINE FOR CANOPY OVER SNOW SUBAREA
  real, intent(in) :: FCANC(ILG,ICTEM)   !< FRACTIONAL COVERAGE OF 8 CARBON PFTs, CANOPY OVER GROUND
  real, intent(in) :: FCANCS(ILG,ICTEM)  !< FRACTIONAL COVERAGE OF 8 CARBON PFTs, CANOPY OVER SNOW
  real, intent(in) :: CO2CONC(ILG)       !< ATMOS. CO2 CONC. IN PPM
  real, intent(in) :: CO2I1CG(ILG,ICTEM) !< INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER GROUND SUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
  real, intent(in) :: CO2I1CS(ILG,ICTEM) !< SAME AS ABOVE BUT FOR SHADED LEAF
  real, intent(in) :: CO2I2CG(ILG,ICTEM) !< INTERCELLULAR CO2 CONC FOR 8 PFTs FOR CANOPY OVER SNOWSUBAREA (Pa) - FOR SINGLE/SUNLIT LEAF
  real, intent(in) :: CO2I2CS(ILG,ICTEM) !< SAME AS ABOVE BUT FOR SHADED LEAF
  real, intent(in) :: COSZS(ILG)         !< COSINE OF SUN'S ZENITH ANGLE
  real, intent(in) :: XDIFFUS(ILG)       !< FRACTION OF DIFFUSED RADIATION
  real, intent(in) :: SLAI(ILG,ICTEM)    !< STORAGE LAI. SEE PHTSYN SUBROUTINE FOR MORE DETAILS.
  real, intent(in) :: RMATCTEM(ILG,ICTEM,IG) !< FRACTION OF ROOTS IN EACH SOIL LAYER FOR EACH OF CTEM's 8 PFTs
  real, intent(in) :: FCANCMX(ILG,ICTEM) !< MAX. FRACTIONAL COVERAGE OF CTEM PFTs
  real, intent(inout) :: ANCSVEG(ILG,ICTEM) !< NET PHOTOSYNTHETIC RATE FOR CTEM's 8 PFTs FOR CANOPY OVER SNOW SUBAREA
  real, intent(inout) :: ANCGVEG(ILG,ICTEM) !< NET PHOTOSYNTHETIC RATE FOR CTEM's 8 PFTs FOR CANOPY OVER GROUND SUBAREA
  real, intent(inout) :: RMLCSVEG(ILG,ICTEM)!< LEAF RESPIRATION RATE FOR CTEM's 8 PFTs FOR CANOPY OVER SNOW SUBAREA
  real, intent(inout) :: RMLCGVEG(ILG,ICTEM)!< LEAF RESPIRATION RATE FOR CTEM's 8 PFTs FOR CANOPY OVER GROUND SUBAREA

  real, intent(in) :: CFLUXCG(ILG)
  real, intent(in) :: CFLUXCS(ILG)
  !
  integer, intent(in) :: ICTEM               !< 8 (CTEM's PLANT FUNCTIONAL TYPES)
  logical, intent(in) :: ctem_on             !< TRUE GIVES COUPLING TO CTEM
  real, intent(in) :: DAYL_MAX(ILG)          !< MAXIMUM DAYLENGTH FOR THAT LOCATION
  real, intent(in) :: DAYL(ILG)              !< DAYLENGTH FOR THAT LOCATION

  integer, intent(in) :: L2MAX, NOL2PFTS(IC)
  !
  !     * INTERNAL WORK ARRAYS FOR THIS ROUTINE.
  !
  real :: VA    (ILG), ZRSLDM(ILG), ZRSLDH(ILG), ZRSLFM(ILG), &
          ZRSLFH(ILG), ZDSLM (ILG), ZDSLH (ILG), TPOTA (ILG), &
          TVIRTA(ILG), CRIB  (ILG), CPHCHC(ILG), CPHCHG(ILG), &
          HCPSCS(ILG), HCPSGS(ILG), CEVAP (ILG), &
          TBAR1P(ILG), GSNOWC(ILG), GSNOWG(ILG), &
          GDENOM(ILG), GCOEFF(ILG), GCONST(ILG), &
          TSNBOT(ILG), GCOEFFS(ILG), GCONSTS(ILG), &
          A1    (ILG), A2    (ILG), B1    (ILG), B2    (ILG), &
          C2    (ILG), ZOM   (ILG), ZOH   (ILG), ZOSCLM(ILG), &
          ZOSCLH(ILG), VAC   (ILG), FCOR  (ILG), &
          CFLUX (ILG), CDHX  (ILG), CDMX  (ILG), &
          QSWX  (ILG), QSWNC (ILG), QSWNG (ILG), QLWX  (ILG), &
          QLWOC (ILG), QLWOG (ILG), QTRANS(ILG), &
          QSENSX(ILG), QSENSC(ILG), QSENSG(ILG), QEVAPX(ILG), &
          QEVAPC(ILG), QEVAPG(ILG), QPHCHC(ILG), QCANX (ILG), &
          TSURX (ILG), QSURX (ILG), &
          TACCS (ILG), QACCS (ILG), TACCO (ILG), QACCO (ILG), &
          ILMOX (ILG), UEX   (ILG), HBLX  (ILG), ZERO  (ILG), &
          STT   (ILG), SQT   (ILG), SUT   (ILG), SVT   (ILG), &
          SHT   (ILG)
  !
  integer             :: IEVAP (ILG), IWATER(ILG)
  !
  !     * INTERNAL WORK ARRAYS FOR energyBudgetPrep.
  !
  real :: FVEG  (ILG), TCSATU(ILG), TCSATF(ILG)
  !
  !     * INTERNAL WORK ARRAYS FOR energBalVegSolve/energBalNoVegSolve.
  !
  real :: TSTEP (ILG), TVIRTC(ILG), TVIRTG(ILG), TVIRTS(ILG), &
          EVBETA(ILG), XEVAP (ILG), EVPWET(ILG), Q0SAT (ILG), &
          RA    (ILG), RB    (ILG), RAGINV(ILG), RBINV (ILG), &
          RBTINV(ILG), RBCINV(ILG), &
          TVRTAC(ILG), TPOTG (ILG), RESID (ILG), &
          TCANP (ILG), TRTOP (ILG), TRTOPG(ILG, NBS), QSTOR (ILG), &
          AC    (ILG), BC    (ILG), &
          LZZ0  (ILG), LZZ0T (ILG), FM    (ILG), FH    (ILG), &
          DCFLXM(ILG), CFLUXM(ILG), WZERO (ILG), XEVAPM(ILG), &
          WC    (ILG), DRAGIN(ILG), CFSENS(ILG), CFEVAP(ILG), &
          QSGADD(ILG), CFLX  (ILG), &
          FTEMPX(ILG), FVAPX (ILG), RIBX  (ILG)
  !
  integer  :: ITER  (ILG), NITER (ILG), JEVAP (ILG), &
              KF    (ILG), KF1   (ILG), KF2   (ILG), &
              IEVAPC(ILG)
  !
  !    Peatland variables
  integer, intent(in) :: ipeatland(ilg), iday
  real, intent(in) :: Cmossmas(ilg), dmoss(ilg), pdd(ilg)
  real, intent(inout) :: ancsmoss(ilg), angsmoss(ilg), &
                         ancmoss(ilg), angmoss(ilg), &
                         rmlcsmoss(ilg), rmlgsmoss(ilg), &
                         rmlcmoss(ilg), rmlgmoss(ilg)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: THTOT, CA, CB, WACSAT, QACSAT, RATIOM, RATIOH, FACTM, FACTH
  !
  !----------------------------------------------------------------------
  !
  !     * CALCULATION OF ATMOSPHERIC INPUT FIELDS REQUIRED BY CLASS FROM
  !     * VARIABLES SUPPLIED BY GCM.
  !
  !> First, two parameters are calculated for later use in the energyBudgetDriver subroutines: the corrected wind speed \f$v_a\f$ ,
  !! and the Coriolis parameter \f$f_{cor}\f$ (describing the effect of the earth’s rotation on the movement of air
  !! according to the reference frame of the earth’s surface). The wind speed correction is applied because it
  !! is assumed that air is never completely still, so \f$v_a\f$ is specified as the maximum of VMOD and a limiting
  !! lower value of \f$0.1 m s^{-1}\f$ . The Coriolis parameter is calculated from the angular velocity \f$\Omega\f$ of the earth’s
  !! rotation (7.29 x 10 -5 radians/s), and the latitude \f$\varphi\f$:
  !! \f$f_{cor} = 2 \Omega sin \varphi\f$
  !!
  !! The packing and unpacking of binary files may cause small shifts in the values of variables at certain
  !! points in the model run, so checks are performed on the depth of ponded water and on the soil liquid and
  !! frozen moisture contents to ensure that unphysical values have not occurred. If the ponded water depth
  !! is vanishingly small or less than zero, it is set to zero. If the soil liquid water content falls below the set
  !! minimum value, it is set to the minimum value. If the soil frozen water content is less than zero, it is set
  !! back to zero. If the sum of the liquid water content and frozen water content converted to an equivalent
  !! liquid amount is greater than the pore volume, both are re-normalized by the pore volume. (This
  !! treatment of frozen water is employed in recognition of the fact that water expands upon freezing and can
  !! therefore occupy a greater volume than the nominal pore volume of the soil.) The small changes in
  !! internal energy and water content of the soil layers resulting from these operations are accounted for by
  !! updating the diagnostic variables HTC and WTRG.
  !!
  !! If CLASS is being run in coupled mode with CTEM, several CTEM variables are initialized to zero.
  !! Subroutine energyBudgetPrep is then called, to carry out the initialization of a number of variables and to do
  !! preparatory calculations of various parameters associated with the energy and water budget calculations.
  !!
  !! The energy budget calculations that follow are performed over the four subareas, of canopy over snow
  !! (CS), snow over ground (GS), canopy over bare ground (C) and bare ground (G). First, a counter
  !! NLANDx is defined for each subarea, containing the number of modelled areas in which the subarea
  !! occurs. The subarea calculations are only done if the relevant counter is greater than zero. A counter
  !! NLANDI is also set to the number of modelled areas that are ice sheets (indicated by ISAND = -4).
  !! (This is used later in subroutine waterBudgetDriver to toggle the call to subroutine iceSheetBalance.)
  !!
  !! For each subarea, the roughness lengths for momentum and heat, ZOM and ZOH, are obtained from
  !! their logarithms, and are then used to determine various height intervals in the atmosphere for subsequent
  !! calculations. These heights are derived from the input variables ZREFM and ZREFH, the reference
  !! heights corresponding to the input values of wind speed and of air temperature and humidity respectively,
  !! and ZDIAGM and ZDIAGH, the heights at which the diagnostic values of anemometer wind speed and
  !! screen level temperature and humidity are to be determined. The form of the calculations depends on the
  !! value of the flag IZREF. If IZREF = 1, the zero plane is taken to lie at the physical ground surface (as
  !! with field measurements); if IZREF = 2, the zero plane is taken to lie at the local roughness height for
  !! momentum(as with atmospheric models). The variables ZRSLDM and ZRSLDH are the height
  !! differences used in subroutine DRCOEF to express the interval between the conceptual bottom of the
  !! atmosphere and the reference heights for wind speed and for temperature and humidity respectively; the
  !! variables ZRSLFM and ZRSLFH are the corresponding height differences used in subroutine
  !! FLXSURFZ. If IZREF = 1, ZRSLDM and ZRSLDH are set to ZREFM and ZREFH minus the
  !! displacement height DISP or DISPS, and ZRSLFM and ZRSLFH are set to ZREFM and ZREFH minus
  !! the roughness length ZOM and DISP or DISPS. If IZREF = 2, ZRSLDM and ZRSLDH are set to
  !! ZREFM and ZREFH plus the roughness height ZOM, and ZRSLFM and ZRSLFH are set to ZREFM
  !! and ZREFH minus DISP or DISPS. (In the absence of a vegetation canopy, the displacement height is
  !! zero.) The variables ZDSLM and ZDSLH are the heights above the bottom of the modelled atmosphere
  !! at which the diagnostic values are to be calculated. If IZREF = 1, they are set to ZDIAGM and
  !! ZDIAGH minus ZOM; if IZREF = 2 they are simply set to ZDIAGM and ZDIAGH. At the end of the
  !! branch in the code, the ratios ZOSCLM and ZOSCLH, which are used in subroutine DRCOEF, are
  !! calculated as ZOM/ZRSLDM and ZOH/ZRSLDH respectively.
  !!
  !! Several other local parameters are also calculated. The potential temperature is the temperature that air
  !! would have if brought adiabatically (without addition or removal of heat) to a given height. The potential
  !! temperature of the air at the reference height, \f$T_{a, pot}\f$ , is calculated relative to the height at which the
  !! horizontal wind speed goes to zero, using the dry adiabatic lapse rate, \f$dT/dz = -g/c_p\f$ , where g is the
  !! acceleration due to gravity and \f$c_p\f$ is the specific heat at constant pressure:. Thus,
  !! \f$T_{a, pot} = T_a + z_{ref, h} g/c_p\f$
  !! where \f$T_a\f$ is the temperature of the air at the reference height and \f$z_{ref, h}\f$ is the height interval, equivalent to
  !! ZRSLFH defined above. If CLASS is being run coupled to an atmospheric model, i.e. if IZREF=2, the
  !! air temperature at the reference height has already been adiabatically extrapolated before being passed to
  !! CLASS. Otherwise, the correction is performed using the above equation.
  !!
  !! The virtual potential temperature of the air at the reference height, \f$T_{a, v}\f$ , is the potential temperature
  !! adjusted for the reduction in air density due to the presence of water vapour. This is applied in order to
  !! enable the use of the equation of state for dry air. \f$T_{a, v}\f$ can be approximated as:
  !! \f$T_{a, v} = T_{a, pot} [1 + 0.61 q_a ]\f$
  !! where \f$q_a\f$ is the specific humidity of the air at the reference height.
  !!
  !! The bulk Richardson number \f$Ri_B\f$ , used in the calculations of the atmospheric stability functions in
  !! subroutine DRCOEF, is formulated as:
  !! \f$Ri_B = [T_0 – T_{a, v} ] (-g z_{ref} )/(T_{a, v} v_a^2)\f$
  !! where \f$T_0\f$ is the surface temperature. For ease of calculation later on, the factor multiplying \f$[T_0 – T_{a, v} ]\f$ on
  !! the right-hand side of the above equation is evaluated and assigned to a coefficient CRIB, using ZRSLDM
  !! for \f$z_{ref}\f$ . The drag coefficient under neutral stability, \f$C_{DN}\f$ , is expressed using basic flux-gradient analysis as:
  !! \f$C_{DN} = k^2 /[ln(z_{ref} ) – ln(z_0)]^2\f$
  !! where k is the von Karman constant and \f$z_0\f$ is the roughness length. ZRSLDM is used for \f$z_{ref}\f$ and the
  !! logarithm of the local roughness length for \f$ln(z_0)\f$, and the neutral drag coefficient DRAG over the
  !! modelled area is obtained as a weighted average over the four subareas.
  !!
  !! For the two subareas with canopy cover, the wind speed of the air at the canopy top, \f$v_{a, c}\f$ , is obtained by
  !! applying the classic logarithmic wind law for the wind speed v(z) at a height z:
  !! \f$kv(z)/v_* = ln[(z – d)/z_0 ]\f$
  !! where \f$v_*\f$ is the friction velocity and d is the displacement height. Thus, \f$v_{a, c}\f$ at the canopy height H can be
  !! related to v a at the reference height \f$z_{ref}\f$ as:
  !! \f$v_{a, c} = v_a [ln(H – d) – ln(z_0)]/[ln(z_{ref} ) – ln(z_0)]\f$
  !!
  !! The vegetation height is calculated as \f$10z_0\f$ . Local values of the temperature of the canopy air TAC and
  !! the humidity of the canopy air QAC are assigned to variables TACCS/TACCO and QACCS/QACCO
  !! respectively.
  !!
  !! At this point calls are made to a series of subroutines addressing the calculation of the energy balance
  !! components of the subarea in question. The calls are summarized in the table below.
  !!
  !! \f[
  !! \begin{array} { | l | l | c | c | c | c | }
  !! \hline
  !! & & \text{CS} & \text{GS} & \text{C} & \text{G} \\
  !! \hline
  !! \text{canopyPhaseChange} & \text{Freezing/thawing of liquid/frozen water on canopy} & \text{YES} & & \text{YES} & \\ \hline
  !! \text{soilHeatFluxPrep} & \text{Set coefficients for temperature calculations in soil} & \text{YES} & \text{YES} & \text{YES} & \text{YES} \\ \hline
  !! \text{snowHeatCond} & \text{Set coefficients for temperature calculations of snow} & \text{YES} & \text{YES} & & \\ \hline
  !! \text{energBalVegSolve} & \text{Calculate components of canopy energy balance} & \text{YES} & & \text{YES} & \\ \hline
  !! \text{energBalNoVegSolve} & \text{Calculate components of ground or snow energy balance} & & \text{YES} & & \text{YES} \\ \hline
  !! \text{snowTempUpdate} & \text{Heat conduction in snow pack} & \text{YES} & \text{YES} & & \\ \hline
  !! \text{soilHeatFluxCleanup} & \text{Heat conduction in soil} & \text{YES} & \text{YES} & \text{YES} & \text{YES} \\ \hline
  !!
  !! \end{array}
  !! \f]
  !!
  !! After these calls, various diagnostic calculations are performed. First the screen-level temperature and
  !! humidity, and the anemometer-level zonal and meridional wind speeds, are calculated. Three options are
  !! provided for doing this, indicated by the flag ISLFD.  If ISLFD = 0, a simple similarity-based approach
  !! is used.  (This is currently the standard option in the AGCM.)  The ratio (RATIOM) of the square root
  !! of the surface drag coefficient for momentum \f$C_{DM}\f$ to that of the neutral drag coefficient \f$C_{DN}\f$ for the
  !! anemometer height is calculated, and an analogous calculation is performed for the ratio (RATIOH) of
  !! the square root of the surface drag coefficient for heat \f$C_{DH}\f$ to that of the neutral drag coefficient
  !! for the screen height \f$z_s\f$, to give a measure of the degree of atmospheric instability.  If the bulk
  !! Richardson number RIB  is negative (indicating unstable conditions), the ratio used for the temperature
  !! and humidity is the minimum of RATIOH and (\f$z_s\f$ / \f$z_{ref}\f$)^(1/3), a measure of the depth of convection.  These
  !! ratios are applied to the calculation of the screen and anemometer level variables.  If the ratios are
  !! large, indicating strong coupling with the atmosphere, the screen level variables tend toward the values
  !! at the reference height; if the ratio is small, they tend to the values at the surface.  At the end of the
  !! block of code, the CCCma subroutine screenRelativeHumidity is called to evaluate the screen-level relative humidity, and
  !! the grid-cell average values of the screen and anemometer level diagnostic variables are incremented.  For
  !! the bare soil subarea, the zonal and meridional anemometer-level wind speeds, the friction velocity and the
  !! surface temperature are saved separately as inputs to the CCCma black carbon deposition scheme
  !! (Namazi et al., 2015) \cite Namazi2015-wz.
  !!
  !! If ISLFD= 1 or 2, the more complex :: calculations in subroutines SLDIAG and DIASURFZ are followed
  !! for the screen-level and anemometer-level diagnostics.
  !! The calculations done in SLDIAG are consistent with the modelling approach used in subroutine
  !! DRCOEF to determine the atmospheric stability functions, so when ISLFD = 1, DRCOEF and
  !! SLDIAG are called. The calculations done in DIASURFZ are consistent with the modelling approach
  !! used in subroutine FLXSURFZ for the atmospheric stability functions, so when ISLFD = 2, FLXSURFZ
  !! and DIASURFZ are called.
  !!
  !! A number of additional diagnostic variables are calculated as weighted averages over the four subareas.
  !! For the most part, these calculations are straightforward; only the calculation of the porential
  !! evapotranspiration \f$E_p\f$ (EVPPOT) involves some complexity. \f$E_p\f$ is defined as the evapotranspiration that
  !! would occur under ambient atmospheric conditions if the soil were completely saturated and the
  !! vegetation canopy were completely water-covered, i.e. if there were no surface resistances to evaporation:
  !! \f$E_p = \rho_a C_{DH} v_a [q_{0, sat} – q_a ]\f$
  !! where \f$\rho_a\f$ is the density of air and \f$q_{0, sat}\f$ is the saturated specific humidity at the surface. For the ground or
  !! snow surface \f$q_{0, sat}\f$ was calculated in subroutine energBalNoVegSolve. For the canopy, the saturated specific humidity
  !! at the canopy air temperature, \f$q_{ac, sat}\f$ , is used. This is obtained from the mixing ratio at saturation, \f$w_{ac, sat}\f$ :
  !! \f$q_{ac, sat} = w_{ac, sat} /[1 + w_{ac, sat} ]\f$
  !!
  !! The mixing ratio is a function of the saturation vapour pressure \f$e_{ac, sat}\f$ at the canopy air temperature:
  !! \f$w_{ac, sat} = 0.622 e_{ac, sat} /(p_{dry} )\f$
  !!
  !! For the saturated vapour pressure, following Emanuel (1994) \cite Emanuel1994-dt
  !! \f$e_{sat}\f$ is from the temperature \f$T_a\f$ and the freezing
  !! point \f$T_f\f$:
  !!
  !! \f$e_{sat} = exp[53.67957 - 6743.769 / T - 4.8451 * ln(T)]       T \geq T_f\f$
  !!
  !! \f$e_{sat} = exp[23.33086 - 6111.72784 / T + 0.15215 * log(T)]    T < T_f\f$
  !!
  !! At the end of the blocks of code dealing with the four subareas, several more diagnostic variables are evaluated.
  !! Again, these calculations are generally straightforward. The effective black-body surface temperature \f$T_{0, eff}\f$
  !! is obtained by inverting the Stefan-Boltzmann equation:
  !! \f$L\uparrow = \sigma T_{0, eff}^4\f$
  !! where \f$L\uparrow\f$ is the outgoing longwave radiation and \f$\sigma\f$ is the Stefan-Boltzmann constant. The evaporation
  !! efficiency parameter EVAPB is calculated as the ratio of the actual evapotranspiration to the potential
  !! evapotranspiration.
  !!
  do I = IL1,IL2 ! loop 50
    VA(I) = MAX(VMIN,VMOD(I))
    FCOR(I) = 2.0 * 7.29E-5 * SIN(RADJ(I))
    !
    !     * CHECK DEPTH OF PONDED WATER FOR UNPHYSICAL VALUES.
    !
    if (ZPOND(I) < 1.0E-8) ZPOND(I) = 0.0
    QG(I) = 0.0
  end do ! loop 50
  !
  !     * CHECK LIQUID AND FROZEN SOIL MOISTURE CONTENTS FOR SMALL
  !     * ABERRATIONS CAUSED BY PACKING/UNPACKING.
  !
  do J = 1,IG ! loop 60
    do I = IL1,IL2
      if (ISAND(I,1) > - 4) then
        HTC(I,J) = HTC(I,J) - TBAR(I,J) * (HCPW * THLIQ(I,J) + &
                   HCPICE * THICE(I,J)) * DELZW(I,J) / DELT
        WTRG(I) = WTRG(I) - (RHOW * THLIQ(I,J) + RHOICE * THICE(I,J)) * &
                  DELZW(I,J) / DELT
        if (THLIQ(I,J) < THLMIN(I,J)) &
            THLIQ(I,J) = THLMIN(I,J)
        if (THICE(I,J) < 0.0) THICE(I,J) = 0.0
        THTOT = THLIQ(I,J) + THICE(I,J) * RHOICE / RHOW
        if (THTOT > THPOR(I,J)) then
          THLIQ(I,J) = MAX(THLIQ(I,J) * THPOR(I,J) / &
                       THTOT,THLMIN(I,J))
          THICE(I,J) = (THPOR(I,J) - THLIQ(I,J)) * &
                       RHOW / RHOICE
          if (THICE(I,J) < 0.0) THICE(I,J) = 0.0
        end if
        HTC(I,J) = HTC(I,J) + TBAR(I,J) * (HCPW * THLIQ(I,J) + &
                   HCPICE * THICE(I,J)) * DELZW(I,J) / DELT
        WTRG(I) = WTRG(I) + (RHOW * THLIQ(I,J) + RHOICE * THICE(I,J)) * &
                  DELZW(I,J) / DELT
      end if
    end do
  end do ! loop 60
  !
  if (ctem_on) then
    !
    !       * INITIALIZE VARIABLES ESTIMATED BY THE PHOTOSYNTHESIS SUBROUTINE
    !       * CALLED FROM WITHIN energBalVegSolve. Also those of moss for peatlands.
    !
    do I = IL1,IL2 ! loop 65
      ancsmoss(i) = 0.0
      angsmoss(i) = 0.0
      ancmoss(i)  = 0.0
      angmoss(i)  = 0.0
      rmlcsmoss(i) = 0.0
      rmlgsmoss(i) = 0.0
      rmlcmoss(i)  = 0.0
      rmlgmoss(i)  = 0.0

      do J = 1,ICTEM
        ANCSVEG(I,J) = 0.0
        ANCGVEG(I,J) = 0.0
        RMLCSVEG(I,J) = 0.0
        RMLCGVEG(I,J) = 0.0
      end do
    end do ! loop 65
  end if
  !
  !     * PREPARATION.
  !
  call  energyBudgetPrep(THLIQC, THLIQG, THICEC, THICEG, TBARC, TBARG, & ! Formerly TPREP
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
  !     * DEFINE NUMBER OF PIXELS OF EACH LAND SURFACE SUBAREA
  !     * (CANOPY-COVERED, CANOPY-AND-SNOW-COVERED, BARE SOIL, AND
  !     * SNOW OVER BARE SOIL) AND NUMBER OF LAND ICE PIXELS FOR
  !     * CALCULATIONS IN energyBudgetDriver/waterBudgetDriver.

  NLANDC = 0
  NLANDCS = 0
  NLANDG = 0
  NLANDGS = 0
  NLANDI = 0

  do I = IL1,IL2 ! loop 70
    if (FC (I) > 0.)            NLANDC = NLANDC + 1
    if (FCS(I) > 0.)            NLANDCS = NLANDCS + 1
    if (FG (I) > 0.)            NLANDG = NLANDG + 1
    if (FGS(I) > 0.)            NLANDGS = NLANDGS + 1
    if (ISAND(I,1) == - 4)        NLANDI = NLANDI + 1
  end do ! loop 70
  !
  !     * CALCULATIONS FOR CANOPY OVER SNOW.
  !
  if (NLANDCS > 0) then
    do I = IL1,IL2 ! loop 100
      if (FCS(I) > 0.) then
        ZOM(I) = EXP(ZOMLCS(I))
        ZOH(I) = EXP(ZOELCS(I))
        if (IZREF == 1) then
          ZRSLDM(I) = ZREFM(I) - DISPS(I)
          ZRSLDH(I) = ZREFH(I) - DISPS(I)
          ZRSLFM(I) = ZREFM(I) - ZOM(I) - DISPS(I)
          ZRSLFH(I) = ZREFH(I) - ZOM(I) - DISPS(I)
          ZDSLM(I) = ZDIAGM(I) - ZOM(I)
          ZDSLH(I) = ZDIAGH(I) - ZOM(I)
          TPOTA(I) = TA(I) + ZRSLFH(I) * GRAV / SPHAIR
        else
          ZRSLDM(I) = ZREFM(I) + ZOM(I)
          ZRSLDH(I) = ZREFH(I) + ZOM(I)
          ZRSLFM(I) = ZREFM(I) - DISPS(I)
          ZRSLFH(I) = ZREFH(I) - DISPS(I)
          ZDSLM(I) = ZDIAGM(I)
          ZDSLH(I) = ZDIAGH(I)
          TPOTA(I) = TA(I)
        end if
        ZOSCLM(I) = ZOM(I) / ZRSLDM(I)
        ZOSCLH(I) = ZOH(I) / ZRSLDH(I)
        TVIRTA(I) = TPOTA(I) * (1.0 + 0.61 * QA(I))
        CRIB(I) = - GRAV * ZRSLDM(I) / (TVIRTA(I) * &
                  VA(I) ** 2)
        DRAG(I) = DRAG(I) + FCS(I) * (VKC / (LOG(ZRSLDM(I)) - &
                  ZOMLCS(I))) ** 2
        VAC(I) = VA(I) * (LOG(10.0 * ZOM(I) - DISPS(I)) - ZOMLCS(I)) / &
                 (LOG(ZRSLDM(I)) - ZOMLCS(I))
        TACCS(I) = TAC(I)
        QACCS(I) = QAC(I)
      end if
    end do ! loop 100
    !
    call canopyPhaseChange(TCANS, RAICNS, SNOCNS, FRAICS, FSNOCS, CHCAPS, & ! Formerly CWCALC
                           HMFC, HTCC, FCS, CMASCS, ILG, IL1, IL2, JL)
    call soilHeatFluxPrep(A1, A2, B1, B2, C2, GDENOM, GCOEFF, & ! Formerly TNPREP
                          GCONST, CPHCHG, IWATER, &
                          TBAR, TCTOPC, TCBOTC, &
                          FCS, ZPOND, TBAR1P, DELZ, TCSNOW, ZSNOW, &
                          ISAND, ILG, IL1, IL2, JL, IG)
    call snowHeatCond(GCOEFFS, GCONSTS, CPHCHG, IWATER, & ! Formerly TSPREP
                      FCS, ZSNOW, TSNOW, TCSNOW, &
                      ILG, IL1, IL2, JL)
    ISNOW = 1
    call energBalVegSolve(ISNOW, FCS, & ! Formerly TSOLVC
                          QSWX, QSWNC, QSWNG, QLWX, QLWOC, QLWOG, QTRANS, &
                          QSENSX, QSENSC, QSENSG, QEVAPX, QEVAPC, QEVAPG, EVAPCS, &
                          EVPCSG, EVAP, TCANS, QCANX, TSURX, QSURX, GSNOWC, QPHCHC, &
                          QMELTC, RAICNS, SNOCNS, CDHX, CDMX, RIBX, TACCS, QACCS, &
                          CFLUX, FTEMPX, FVAPX, ILMOX, UEX, HBLX, QFCF, QFCL, HTCC, &
                          QSWINV, QSWINI, QLWIN, TPOTA, TA, QA, VA, VAC, PADRY, &
                          RHOAIR, ALVSCS, ALIRCS, ALVSSC, ALIRSC, TRVSCS, TRIRCS, &
                          FSVFS, CRIB, CPHCHC, CPHCHG, CEVAP, TADP, TVIRTA, RCS, &
                          RBCOEF, ZOSCLH, ZOSCLM, ZRSLFH, ZRSLFM, ZOH, ZOM, &
                          FCOR, GCONSTS, GCOEFFS, TSFSAV(1, 1), TRSNOWC, FSNOCS, &
                          FRAICS, CHCAPS, CMASCS, PCPR, FROOTS, THLMIN, DELZW, &
                          RHOSCS, ZSNOW, IWATER, IEVAP, ITERCT, &
                          ISLFD, ITC, ITCG, ILG, IL1, IL2, JL, N, &
                          TSTEP, TVIRTC, TVIRTG, EVBETA, XEVAP, EVPWET, Q0SAT, &
                          RA, RB, RAGINV, RBINV, RBTINV, RBCINV, TVRTAC, TPOTG, &
                          RESID, TCANP, &
                          WZERO, XEVAPM, DCFLXM, WC, DRAGIN, CFLUXM, CFLX, IEVAPC, &
                          TRTOP, QSTOR, CFSENS, CFEVAP, QSGADD, AC, BC, &
                          LZZ0, LZZ0T, FM, FH, ITER, NITER, KF1, KF2, &
                          AILCGS, FCANCS, CO2CONC, RMATCTEM, &
                          THLIQC, THFC, THLW, ISAND, IG, COSZS, PRESSG, &
                          XDIFFUS, ICTEM, IC, CO2I1CS, CO2I2CS, &
                          ctem_on, SLAI, FCANCMX, L2MAX, &
                          NOL2PFTS, CFLUXCS, ANCSVEG, RMLCSVEG, &
                          DAYL, DAYL_MAX, ipeatland, Cmossmas, dmoss, &
                          ancsmoss, rmlcsmoss, iday, pdd)

    call snowTempUpdate(GSNOWC, TSNOCS, WSNOCS, RHOSCS, QMELTC, & ! Formerly TSPOST
                        GZROCS, TSNBOT, HTCS, HMFN, &
                        GCONSTS, GCOEFFS, GCONST, GCOEFF, TBAR, &
                        TSURX, ZSNOW, TCSNOW, HCPSCS, QTRANS, &
                        FCS, DELZ, ILG, IL1, IL2, JL, IG)
    call soilHeatFluxCleanup(TBARCS, G12CS, G23CS, TPNDCS, GZROCS, ZERO, GCONST, & ! Formerly TNPOST
                             GCOEFF, TBAR, TCTOPC, TCBOTC, HCPC, ZPOND, TSNBOT, &
                             TBASE, TBAR1P, A1, A2, B1, B2, C2, FCS, IWATER, &
                             ISAND, DELZ, DELZW, ILG, IL1, IL2, JL, IG)
    !
    !     * DIAGNOSTICS.
    !
    if (ISLFD == 0) then
      do I = IL1,IL2 ! loop 150
        if (FCS(I) > 0.) then
          FACTM = ZDSLM(I) + ZOM(I)
          FACTH = ZDSLH(I) + ZOM(I)
          RATIOM = SQRT(CDMX(I)) * LOG(FACTM / ZOM(I)) / VKC
          RATIOM = MIN(RATIOM,1.)
          RATIOH = SQRT(CDMX(I)) * LOG(FACTH / ZOH(I)) / VKC
          RATIOH = MIN(RATIOH,1.)
          if (RIBX(I) < 0.) then
            RATIOH = RATIOH * CDHX(I) / CDMX(I)
            RATIOH = MIN(RATIOH,(FACTH / ZRSLDH(I)) ** (1. / 3.))
          end if
          STT(I) = TACCS(I) - (MIN(RATIOH,1.)) * (TACCS(I) - TA(I))
          SQT(I) = QACCS(I) - (MIN(RATIOH,1.)) * (QACCS(I) - QA(I))
          SUT(I) = RATIOM * UWIND(I)
          SVT(I) = RATIOM * VWIND(I)
        end if
      end do ! loop 150
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FCS, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FCS(I) > 0.) then
          ST (I) = ST (I) + FCS(I) * STT(I)
          SQ (I) = SQ (I) + FCS(I) * SQT(I)
          SU (I) = SU (I) + FCS(I) * SUT(I)
          SV (I) = SV (I) + FCS(I) * SVT(I)
          SRH(I) = SRH(I) + FCS(I) * SHT(I)
        end if
      end do
      !
    else if (ISLFD == 1) then
      call sldiag1(SUT, SVT, STT, SQT, &
                  CDMX, CDHX, UWIND, VWIND, TPOTA, QA, &
                  TACCS, QACCS, ZOM, ZOH, FCS, ZRSLDM, &
                  ZDSLM, ZDSLH, ILG, IL1, IL2, JL)
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FCS, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FCS(I) > 0.) then
          ST (I) = ST (I) + FCS(I) * STT(I)
          SQ (I) = SQ (I) + FCS(I) * SQT(I)
          SU (I) = SU (I) + FCS(I) * SUT(I)
          SV (I) = SV (I) + FCS(I) * SVT(I)
          SRH(I) = SRH(I) + FCS(I) * SHT(I)
        end if
      end do
    else if (ISLFD == 2) then
      call DIASURFZ1(SU, SV, ST, SQ, ILG, UWIND, VWIND, TACCS, QACCS, &
                    ZOM, ZOH, ILMOX, ZRSLFM, HBLX, UEX, FTEMPX, FVAPX, &
                    ZDSLM, ZDSLH, RADJ, FCS, IL1, IL2, JL)
    end if
    !
    do I = IL1,IL2 ! loop 175
      if (FCS(I) > 0.) then
        WACSAT = 0.622 * calcEsat(TACCS(I)) / PADRY(I)
        QACSAT = WACSAT / (1.0 + WACSAT)
        EVPPOT(I) = EVPPOT(I) + FCS(I) * RHOAIR(I) * CFLUX(I) * &
                    (QACSAT - QA(I))
        ACOND(I) = ACOND(I) + FCS(I) * CFLUX(I)
        ILMO(I) = ILMO(I) + FCS(I) * ILMOX(I)
        UE(I)   = UE(I) + FCS(I) * UEX(I)
        HBL(I)  = HBL(I) + FCS(I) * HBLX(I)
        CDH (I) = CDH(I) + FCS(I) * CDHX(I)
        CDM (I) = CDM(I) + FCS(I) * CDMX(I)
        TSFSAV(I,1) = TSURX(I)
        QG(I) = QG(I) + FCS(I) * QACCS(I)
        QSENS(I) = QSENS(I) + FCS(I) * QSENSX(I)
        QEVAP(I) = QEVAP(I) + FCS(I) * QEVAPX(I)
        QLWAVG(I) = QLWAVG(I) + FCS(I) * QLWX(I)
        FSGV(I) = FSGV(I) + FCS(I) * QSWNC(I)
        FSGS(I) = FSGS(I) + FCS(I) * QSWNG(I)
        FSGG(I) = FSGG(I) + FCS(I) * QTRANS(I)
        FLGV(I) = FLGV(I) + FCS(I) * (QLWIN(I) + QLWOG(I) - 2.0 * &
                  QLWOC(I)) * (1.0 - FSVFS(I))
        FLGS(I) = FLGS(I) + FCS(I) * (QLWOC(I) * (1.0 - FSVFS(I)) + &
                  QLWIN(I) * FSVFS(I) - QLWOG(I))
        if (ITC == 1) then
          HFSC(I) = HFSC(I) + FCS(I) * QSENSC(I)
        else
          HFSC(I) = HFSC(I) + FCS(I) * (QSENSC(I) - QSENSG(I))
        end if
        HFSS(I) = HFSS(I) + FCS(I) * QSENSG(I)
        HEVC(I) = HEVC(I) + FCS(I) * QEVAPC(I)
        HEVS(I) = HEVS(I) + FCS(I) * QEVAPG(I)
        HMFC(I) = HMFC(I) + FCS(I) * QPHCHC(I)
        HTCS(I) = HTCS(I) + FCS(I) * ( - GZROCS(I) + &
                  QTRANS(I))
        HTC(I,1) = HTC(I,1) + FCS(I) * (GZROCS(I) - QTRANS(I) - &
                   G12CS(I))
        HTC(I,2) = HTC(I,2) + FCS(I) * (G12CS(I) - G23CS(I))
        HTC(I,3) = HTC(I,3) + FCS(I) * G23CS(I)
        FTEMP(I) = FTEMP(I) + FCS(I) * FTEMPX(I)
        FVAP (I) = FVAP (I) + FCS(I) * FVAPX (I)
        RIB  (I) = RIB  (I) + FCS(I) * RIBX  (I)
        GSNOW(I) = GSNOW(I) + FCS(I) / (FCS(I) + FGS(I)) * GSNOWC(I)
      end if
    end do ! loop 175
  end if
  !
  !     * CALCULATIONS FOR SNOW-COVERED GROUND.
  !
  if (NLANDGS > 0) then
    do I = IL1,IL2 ! loop 200
      if (FGS(I) > 0.) then
        ZOM(I) = EXP(ZOMLNS(I))
        ZOH(I) = EXP(ZOELNS(I))
        if (IZREF == 1) then
          ZRSLDM(I) = ZREFM(I)
          ZRSLDH(I) = ZREFH(I)
          ZRSLFM(I) = ZREFM(I) - ZOM(I)
          ZRSLFH(I) = ZREFH(I) - ZOM(I)
          ZDSLM(I) = ZDIAGM(I) - ZOM(I)
          ZDSLH(I) = ZDIAGH(I) - ZOM(I)
          TPOTA(I) = TA(I) + ZRSLFH(I) * GRAV / SPHAIR
        else
          ZRSLDM(I) = ZREFM(I) + ZOM(I)
          ZRSLDH(I) = ZREFH(I) + ZOM(I)
          ZRSLFM(I) = ZREFM(I)
          ZRSLFH(I) = ZREFH(I)
          ZDSLM(I) = ZDIAGM(I)
          ZDSLH(I) = ZDIAGH(I)
          TPOTA(I) = TA(I)
        end if
        ZOSCLM(I) = ZOM(I) / ZRSLDM(I)
        ZOSCLH(I) = ZOH(I) / ZRSLDH(I)
        TVIRTA(I) = TPOTA(I) * (1.0 + 0.61 * QA(I))
        CRIB(I) = - GRAV * ZRSLDM(I) / (TVIRTA(I) * VA(I) ** 2)
        DRAG(I) = DRAG(I) + FGS(I) * (VKC / (LOG(ZRSLDM(I)) - &
                  ZOMLNS(I))) ** 2
      end if
    end do ! loop 200
    !
    call soilHeatFluxPrep(A1, A2, B1, B2, C2, GDENOM, GCOEFF, & ! Formerly TNPREP
                          GCONST, CPHCHG, IWATER, &
                          TBAR, TCTOPG, TCBOTG, &
                          FGS, ZPOND, TBAR1P, DELZ, TCSNOW, ZSNOW, &
                          ISAND, ILG, IL1, IL2, JL, IG)
    call snowHeatCond(GCOEFFS, GCONSTS, CPHCHG, IWATER, & ! Formerly TSPREP
                      FGS, ZSNOW, TSNOW, TCSNOW, &
                      ILG, IL1, IL2, JL)
    ISNOW = 1
    call energBalNoVegSolve(ISNOW, FGS, & ! Formerly TSOLVE
                            QSWX, QLWX, QTRANS, QSENSX, QEVAPX, EVAPGS, &
                            TSURX, QSURX, GSNOWG, QMELTG, CDHX, CDMX, RIBX, CFLUX, &
                            FTEMPX, FVAPX, ILMOX, UEX, HBLX, &
                            QLWIN, TPOTA, QA, VA, PADRY, RHOAIR, &
                            ALVSSN, ALIRSN, CRIB, CPHCHG, CEVAP, TVIRTA, &
                            ZOSCLH, ZOSCLM, ZRSLFH, ZRSLFM, ZOH, ZOM, FCOR, &
                            GCONSTS, GCOEFFS, TSFSAV(1, 2), PCPR, &
                            TRSNOWG, FSSB, ALSNO, &
                            THLIQG, THLMIN, DELZW, RHOSGS, ZSNOW, ZPOND, &
                            IWATER, IEVAP, ITERCT, ISAND, &
                            ISLFD, ITG, ILG, IG, IL1, IL2, JL, NBS, ISNOALB, &
                            TSTEP, TVIRTS, EVBETA, Q0SAT, RESID, &
                            DCFLXM, CFLUXM, WZERO, TRTOPG, AC, BC, &
                            LZZ0, LZZ0T, FM, FH, ITER, NITER, JEVAP, KF, &
                            ipeatland, co2conc, pressg, coszs, Cmossmas, dmoss, &
                            angsmoss, rmlgsmoss, iday, DAYL, pdd)
    !
    call snowTempUpdate(GSNOWG, TSNOGS, WSNOGS, RHOSGS, QMELTG, & ! Formerly TSPOST
                        GZROGS, TSNBOT, HTCS, HMFN, &
                        GCONSTS, GCOEFFS, GCONST, GCOEFF, TBAR, &
                        TSURX, ZSNOW, TCSNOW, HCPSGS, QTRANS, &
                        FGS, DELZ, ILG, IL1, IL2, JL, IG)
    call soilHeatFluxCleanup(TBARGS, G12GS, G23GS, TPNDGS, GZROGS, ZERO, GCONST, & ! Formerly TNPOST
                             GCOEFF, TBAR, TCTOPG, TCBOTG, HCPG, ZPOND, TSNBOT, &
                             TBASE, TBAR1P, A1, A2, B1, B2, C2, FGS, IWATER, &
                             ISAND, DELZ, DELZW, ILG, IL1, IL2, JL, IG)
    !
    !     * DIAGNOSTICS.
    !
    if (ISLFD == 0) then
      do I = IL1,IL2 ! loop 250
        if (FGS(I) > 0.) then
          FACTM = ZDSLM(I) + ZOM(I)
          FACTH = ZDSLH(I) + ZOM(I)
          RATIOM = SQRT(CDMX(I)) * LOG(FACTM / ZOM(I)) / VKC
          RATIOM = MIN(RATIOM,1.)
          RATIOH = SQRT(CDMX(I)) * LOG(FACTH / ZOH(I)) / VKC
          RATIOH = MIN(RATIOH,1.)
          if (RIBX(I) < 0.) then
            RATIOH = RATIOH * CDHX(I) / CDMX(I)
            RATIOH = MIN(RATIOH,(FACTH / ZRSLDH(I)) ** (1. / 3.))
          end if
          STT(I) = TSURX(I) - (MIN(RATIOH,1.)) * (TSURX(I) - TA(I))
          SQT(I) = QSURX(I) - (MIN(RATIOH,1.)) * (QSURX(I) - QA(I))
          SUT(I) = RATIOM * UWIND(I)
          SVT(I) = RATIOM * VWIND(I)
        end if
      end do ! loop 250
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FGS, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FGS(I) > 0.) then
          ST (I) = ST (I) + FGS(I) * STT(I)
          SQ (I) = SQ (I) + FGS(I) * SQT(I)
          SU (I) = SU (I) + FGS(I) * SUT(I)
          SV (I) = SV (I) + FGS(I) * SVT(I)
          SRH(I) = SRH(I) + FGS(I) * SHT(I)
        end if
      end do
      !
    else if (ISLFD == 1) then
      call sldiag1(SUT, SVT, STT, SQT, &
                  CDMX, CDHX, UWIND, VWIND, TPOTA, QA, &
                  TSURX, QSURX, ZOM, ZOH, FGS, ZRSLDM, &
                  ZDSLM, ZDSLH, ILG, IL1, IL2, JL)
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FGS, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FGS(I) > 0.) then
          ST (I) = ST (I) + FGS(I) * STT(I)
          SQ (I) = SQ (I) + FGS(I) * SQT(I)
          SU (I) = SU (I) + FGS(I) * SUT(I)
          SV (I) = SV (I) + FGS(I) * SVT(I)
          SRH(I) = SRH(I) + FGS(I) * SHT(I)
        end if
      end do
    else if (ISLFD == 2) then
      call DIASURFZ1(SU, SV, ST, SQ, ILG, UWIND, VWIND, TSURX, QSURX, &
                    ZOM, ZOH, ILMOX, ZRSLFM, HBLX, UEX, FTEMPX, FVAPX, &
                    ZDSLM, ZDSLH, RADJ, FGS, IL1, IL2, JL)
    end if
    !
    do I = IL1,IL2 ! loop 275
      if (FGS(I) > 0.) then
        EVPPOT(I) = EVPPOT(I) + FGS(I) * RHOAIR(I) * CFLUX(I) * &
                    (Q0SAT(I) - QA(I))
        ACOND(I) = ACOND(I) + FGS(I) * CFLUX(I)
        ILMO(I) = ILMO(I) + FGS(I) * ILMOX(I)
        UE(I)   = UE(I) + FGS(I) * UEX(I)
        HBL(I)  = HBL(I) + FGS(I) * HBLX(I)
        CDH (I) = CDH(I) + FGS(I) * CDHX(I)
        CDM (I) = CDM(I) + FGS(I) * CDMX(I)
        TSFSAV(I,2) = TSURX(I)
        QG(I) = QG(I) + FGS(I) * QSURX(I)
        QSENS(I) = QSENS(I) + FGS(I) * QSENSX(I)
        QEVAP(I) = QEVAP(I) + FGS(I) * QEVAPX(I)
        QLWAVG(I) = QLWAVG(I) + FGS(I) * QLWX(I)
        FSGS(I) = FSGS(I) + FGS(I) * (QSWX(I) - QTRANS(I))
        FSGG(I) = FSGG(I) + FGS(I) * QTRANS(I)
        FLGS(I) = FLGS(I) + FGS(I) * (QLWIN(I) - QLWX(I))
        HFSS(I) = HFSS(I) + FGS(I) * QSENSX(I)
        HEVS(I) = HEVS(I) + FGS(I) * QEVAPX(I)
        HTCS(I) = HTCS(I) + FGS(I) * ( - GZROGS(I) + &
                  QTRANS(I))
        HTC(I,1) = HTC(I,1) + FGS(I) * (GZROGS(I) - QTRANS(I) - &
                   G12GS(I))
        HTC(I,2) = HTC(I,2) + FGS(I) * (G12GS(I) - G23GS(I))
        HTC(I,3) = HTC(I,3) + FGS(I) * G23GS(I)
        FTEMP(I) = FTEMP(I) + FGS(I) * FTEMPX(I)
        FVAP (I) = FVAP (I) + FGS(I) * FVAPX (I)
        RIB  (I) = RIB  (I) + FGS(I) * RIBX  (I)
        GSNOW(I) = GSNOW(I) + FGS(I) / (FCS(I) + FGS(I)) * GSNOWG(I)
      end if
    end do ! loop 275
  end if
  !
  !     * CALCULATIONS FOR CANOPY OVER BARE GROUND.
  !
  if (NLANDC > 0) then
    do I = IL1,IL2 ! loop 300
      if (FC(I) > 0.) then
        ZOM(I) = EXP(ZOMLNC(I))
        ZOH(I) = EXP(ZOELNC(I))
        if (IZREF == 1) then
          ZRSLDM(I) = ZREFM(I) - DISP(I)
          ZRSLDH(I) = ZREFH(I) - DISP(I)
          ZRSLFM(I) = ZREFM(I) - ZOM(I) - DISP(I)
          ZRSLFH(I) = ZREFH(I) - ZOM(I) - DISP(I)
          ZDSLM(I) = ZDIAGM(I) - ZOM(I)
          ZDSLH(I) = ZDIAGH(I) - ZOM(I)
          TPOTA(I) = TA(I) + ZRSLFH(I) * GRAV / SPHAIR
        else
          ZRSLDM(I) = ZREFM(I) + ZOM(I)
          ZRSLDH(I) = ZREFH(I) + ZOM(I)
          ZRSLFM(I) = ZREFM(I) - DISP(I)
          ZRSLFH(I) = ZREFH(I) - DISP(I)
          ZDSLM(I) = ZDIAGM(I)
          ZDSLH(I) = ZDIAGH(I)
          TPOTA(I) = TA(I)
        end if
        ZOSCLM(I) = ZOM(I) / ZRSLDM(I)
        ZOSCLH(I) = ZOH(I) / ZRSLDH(I)
        TVIRTA(I) = TPOTA(I) * (1.0 + 0.61 * QA(I))
        CRIB(I) = - GRAV * ZRSLDM(I) / (TVIRTA(I) * VA(I) ** 2)
        DRAG(I) = DRAG(I) + FC(I) * (VKC / (LOG(ZRSLDM(I)) - &
                  ZOMLNC(I))) ** 2
        VAC(I) = VA(I) * (LOG(10.0 * ZOM(I) - DISP(I)) - ZOMLNC(I)) / &
                 (LOG(ZRSLDM(I)) - ZOMLNC(I))
        TACCO(I) = TAC(I)
        QACCO(I) = QAC(I)
      end if
    end do ! loop 300
    !
    call canopyPhaseChange(TCANO, RAICAN, SNOCAN, FRAINC, FSNOWC, CHCAP, & ! Formerly CWCALC
                           HMFC, HTCC, FC, CMASSC, ILG, IL1, IL2, JL)
    call soilHeatFluxPrep(A1, A2, B1, B2, C2, GDENOM, GCOEFF, & ! Formerly TNPREP
                          GCONST, CPHCHG, IWATER, &
                          TBAR, TCTOPC, TCBOTC, &
                          FC, ZPOND, TBAR1P, DELZ, TCSNOW, ZSNOW, &
                          ISAND, ILG, IL1, IL2, JL, IG)
    ISNOW = 0

    call energBalVegSolve(ISNOW, FC, & ! Formerly TSOLVC
                          QSWX, QSWNC, QSWNG, QLWX, QLWOC, QLWOG, QTRANS, &
                          QSENSX, QSENSC, QSENSG, QEVAPX, QEVAPC, QEVAPG, EVAPC, &
                          EVAPCG, EVAP, TCANO, QCANX, TSURX, QSURX, GZEROC, QPHCHC, &
                          QFREZC, RAICAN, SNOCAN, CDHX, CDMX, RIBX, TACCO, QACCO, &
                          CFLUX, FTEMPX, FVAPX, ILMOX, UEX, HBLX, QFCF, QFCL, HTCC, &
                          QSWINV, QSWINI, QLWIN, TPOTA, TA, QA, VA, VAC, PADRY, &
                          RHOAIR, ALVSCN, ALIRCN, ALVSGC, ALIRGC, TRVSCN, TRIRCN, &
                          FSVF, CRIB, CPHCHC, CPHCHG, CEVAP, TADP, TVIRTA, RC, &
                          RBCOEF, ZOSCLH, ZOSCLM, ZRSLFH, ZRSLFM, ZOH, ZOM, &
                          FCOR, GCONST, GCOEFF, TSFSAV(1, 3), TRSNOWC, FSNOWC, &
                          FRAINC, CHCAP, CMASSC, PCPR, FROOT, THLMIN, DELZW, &
                          ZERO, ZERO, IWATER, IEVAP, ITERCT, &
                          ISLFD, ITC, ITCG, ILG, IL1, IL2, JL, N, &
                          TSTEP, TVIRTC, TVIRTG, EVBETA, XEVAP, EVPWET, Q0SAT, &
                          RA, RB, RAGINV, RBINV, RBTINV, RBCINV, TVRTAC, TPOTG, &
                          RESID, TCANP, &
                          WZERO, XEVAPM, DCFLXM, WC, DRAGIN, CFLUXM, CFLX, IEVAPC, &
                          TRTOP, QSTOR, CFSENS, CFEVAP, QSGADD, AC, BC, &
                          LZZ0, LZZ0T, FM, FH, ITER, NITER, KF1, KF2, &
                          AILCG, FCANC, CO2CONC, RMATCTEM, &
                          THLIQC, THFC, THLW, ISAND, IG, COSZS, PRESSG, &
                          XDIFFUS, ICTEM, IC, CO2I1CG, CO2I2CG, &
                          ctem_on, SLAI, FCANCMX, L2MAX, &
                          NOL2PFTS, CFLUXCG, ANCGVEG, RMLCGVEG, &
                          DAYL, DAYL_MAX, ipeatland, Cmossmas, dmoss, &
                          ancmoss, rmlcmoss, iday, pdd)


    call soilHeatFluxCleanup(TBARC, G12C, G23C, TPONDC, GZEROC, QFREZC, GCONST, & ! Formerly TNPOST
                             GCOEFF, TBAR, TCTOPC, TCBOTC, HCPC, ZPOND, TSURX, &
                             TBASE, TBAR1P, A1, A2, B1, B2, C2, FC, IWATER, &
                             ISAND, DELZ, DELZW, ILG, IL1, IL2, JL, IG)
    !
    !     * DIAGNOSTICS.
    !
    if (ISLFD == 0) then
      do I = IL1,IL2 ! loop 350
        if (FC(I) > 0.) then
          FACTM = ZDSLM(I) + ZOM(I)
          FACTH = ZDSLH(I) + ZOM(I)
          RATIOM = SQRT(CDMX(I)) * LOG(FACTM / ZOM(I)) / VKC
          RATIOM = MIN(RATIOM,1.)
          RATIOH = SQRT(CDMX(I)) * LOG(FACTH / ZOH(I)) / VKC
          RATIOH = MIN(RATIOH,1.)
          if (RIBX(I) < 0.) then
            RATIOH = RATIOH * CDHX(I) / CDMX(I)
            RATIOH = MIN(RATIOH,(FACTH / ZRSLDH(I)) ** (1. / 3.))
          end if
          STT(I) = TACCO(I) - (MIN(RATIOH,1.)) * (TACCO(I) - TA(I))
          SQT(I) = QACCO(I) - (MIN(RATIOH,1.)) * (QACCO(I) - QA(I))
          SUT(I) = RATIOM * UWIND(I)
          SVT(I) = RATIOM * VWIND(I)
        end if
      end do ! loop 350
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FC, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FC(I) > 0.) then
          ST (I) = ST (I) + FC(I) * STT(I)
          SQ (I) = SQ (I) + FC(I) * SQT(I)
          SU (I) = SU (I) + FC(I) * SUT(I)
          SV (I) = SV (I) + FC(I) * SVT(I)
          SRH(I) = SRH(I) + FC(I) * SHT(I)
        end if
      end do
    else if (ISLFD == 1) then
      call sldiag1(SUT, SVT, STT, SQT, &
                  CDMX, CDHX, UWIND, VWIND, TPOTA, QA, &
                  TACCO, QACCO, ZOM, ZOH, FC, ZRSLDM, &
                  ZDSLM, ZDSLH, ILG, IL1, IL2, JL)
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FC, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FC(I) > 0.) then
          ST (I) = ST (I) + FC(I) * STT(I)
          SQ (I) = SQ (I) + FC(I) * SQT(I)
          SU (I) = SU (I) + FC(I) * SUT(I)
          SV (I) = SV (I) + FC(I) * SVT(I)
          SRH(I) = SRH(I) + FC(I) * SHT(I)
        end if
      end do
    else if (ISLFD == 2) then
      call DIASURFZ1(SU, SV, ST, SQ, ILG, UWIND, VWIND, TACCO, QACCO, &
                    ZOM, ZOH, ILMOX, ZRSLFM, HBLX, UEX, FTEMPX, FVAPX, &
                    ZDSLM, ZDSLH, RADJ, FC, IL1, IL2, JL)
    end if
    !
    do I = IL1,IL2 ! loop 375
      if (FC(I) > 0.) then
        WACSAT = 0.622 * calcEsat(TACCO(I)) / PADRY(I)
        QACSAT = WACSAT / (1.0 + WACSAT)
        EVPPOT(I) = EVPPOT(I) + FC(I) * RHOAIR(I) * CFLUX(I) * &
                    (QACSAT - QA(I))
        ACOND(I) = ACOND(I) + FC(I) * CFLUX(I)
        ILMO(I) = ILMO(I) + FC(I) * ILMOX(I)
        UE(I)   = UE(I) + FC(I) * UEX(I)
        HBL(I)  = HBL(I) + FC(I) * HBLX(I)
        CDH (I) = CDH(I) + FC(I) * CDHX(I)
        CDM (I) = CDM(I) + FC(I) * CDMX(I)
        TSFSAV(I,3) = TSURX(I)
        QG(I) = QG(I) + FC(I) * QACCO(I)
        QSENS(I) = QSENS(I) + FC(I) * QSENSX(I)
        QEVAP(I) = QEVAP(I) + FC(I) * QEVAPX(I)
        QLWAVG(I) = QLWAVG(I) + FC(I) * QLWX(I)
        FSGV(I) = FSGV(I) + FC(I) * QSWNC(I)
        FSGG(I) = FSGG(I) + FC(I) * QSWNG(I)
        FLGV(I) = FLGV(I) + FC(I) * (QLWIN(I) + QLWOG(I) - 2.0 * &
                  QLWOC(I)) * (1.0 - FSVF(I))
        FLGG(I) = FLGG(I) + FC(I) * (FSVF(I) * QLWIN(I) + &
                  (1.0 - FSVF(I)) * QLWOC(I) - QLWOG(I))
        if (ITC == 1) then
          HFSC(I) = HFSC(I) + FC(I) * QSENSC(I)
        else
          HFSC(I) = HFSC(I) + FC(I) * (QSENSC(I) - QSENSG(I))
        end if
        HFSG(I) = HFSG(I) + FC(I) * QSENSG(I)
        HEVC(I) = HEVC(I) + FC(I) * QEVAPC(I)
        HEVG(I) = HEVG(I) + FC(I) * QEVAPG(I)
        HMFC(I) = HMFC(I) + FC(I) * QPHCHC(I)
        HTC(I,1) = HTC(I,1) + FC(I) * ( - G12C(I))
        HTC(I,2) = HTC(I,2) + FC(I) * (G12C(I) - G23C(I))
        HTC(I,3) = HTC(I,3) + FC(I) * G23C(I)
        FTEMP(I) = FTEMP(I) + FC(I) * FTEMPX(I)
        FVAP (I) = FVAP (I) + FC(I) * FVAPX (I)
        RIB  (I) = RIB  (I) + FC(I) * RIBX   (I)
      end if
    end do ! loop 375
  end if
  !
  !     * CALCULATIONS FOR BARE GROUND.
  !
  if (NLANDG > 0) then
    do I = IL1,IL2 ! loop 400
      if (FG(I) > 0.) then
        ZOM(I) = EXP(ZOMLNG(I))
        ZOH(I) = EXP(ZOELNG(I))
        if (IZREF == 1) then
          ZRSLDM(I) = ZREFM(I)
          ZRSLDH(I) = ZREFH(I)
          ZRSLFM(I) = ZREFM(I) - ZOM(I)
          ZRSLFH(I) = ZREFH(I) - ZOM(I)
          ZDSLM(I) = ZDIAGM(I) - ZOM(I)
          ZDSLH(I) = ZDIAGH(I) - ZOM(I)
          TPOTA(I) = TA(I) + ZRSLFH(I) * GRAV / SPHAIR
        else
          ZRSLDM(I) = ZREFM(I) + ZOM(I)
          ZRSLDH(I) = ZREFH(I) + ZOM(I)
          ZRSLFM(I) = ZREFM(I)
          ZRSLFH(I) = ZREFH(I)
          ZDSLM(I) = ZDIAGM(I)
          ZDSLH(I) = ZDIAGH(I)
          TPOTA(I) = TA(I)
        end if
        ZOSCLM(I) = ZOM(I) / ZRSLDM(I)
        ZOSCLH(I) = ZOH(I) / ZRSLDH(I)
        TVIRTA(I) = TPOTA(I) * (1.0 + 0.61 * QA(I))
        CRIB(I) = - GRAV * ZRSLDM(I) / (TVIRTA(I) * VA(I) ** 2)
        DRAG(I) = DRAG(I) + FG(I) * (VKC / (LOG(ZRSLDM(I)) - &
                  ZOMLNG(I))) ** 2
      end if
    end do ! loop 400
    !
    call soilHeatFluxPrep(A1, A2, B1, B2, C2, GDENOM, GCOEFF, & ! Formerly TNPREP
                          GCONST, CPHCHG, IWATER, &
                          TBAR, TCTOPG, TCBOTG, &
                          FG, ZPOND, TBAR1P, DELZ, TCSNOW, ZSNOW, &
                          ISAND, ILG, IL1, IL2, JL, IG)
    ISNOW = 0
    call energBalNoVegSolve(ISNOW, FG, & ! Formerly TSOLVE
                            QSWX, QLWX, QTRANS, QSENSX, QEVAPX, EVAPG, &
                            TSURX, QSURX, GZEROG, QFREZG, CDHX, CDMX, RIBX, CFLUX, &
                            FTEMPX, FVAPX, ILMOX, UEX, HBLX, &
                            QLWIN, TPOTA, QA, VA, PADRY, RHOAIR, &
                            ALVSG, ALIRG, CRIB, CPHCHG, CEVAP, TVIRTA, &
                            ZOSCLH, ZOSCLM, ZRSLFH, ZRSLFM, ZOH, ZOM, FCOR, &
                            GCONST, GCOEFF, TSFSAV(1, 4), PCPR, &
                            TRSNOWG, FSSB, ALSNO, &
                            THLIQG, THLMIN, DELZW, ZERO, ZERO, ZPOND, &
                            IWATER, IEVAP, ITERCT, ISAND, &
                            ISLFD, ITG, ILG, IG, IL1, IL2, JL, NBS, ISNOALB, &
                            TSTEP, TVIRTS, EVBETA, Q0SAT, RESID, &
                            DCFLXM, CFLUXM, WZERO, TRTOPG, AC, BC, &
                            LZZ0, LZZ0T, FM, FH, ITER, NITER, JEVAP, KF, &
                            ipeatland, co2conc, pressg, coszs, Cmossmas, dmoss, &
                            angmoss, rmlgmoss, iday, DAYL, pdd)
    !
    call soilHeatFluxCleanup(TBARG, G12G, G23G, TPONDG, GZEROG, QFREZG, GCONST, & ! Formerly TNPOST
                             GCOEFF, TBAR, TCTOPG, TCBOTG, HCPG, ZPOND, TSURX, &
                             TBASE, TBAR1P, A1, A2, B1, B2, C2, FG, IWATER, &
                             ISAND, DELZ, DELZW, ILG, IL1, IL2, JL, IG)
    !
    !     * DIAGNOSTICS.
    !
    if (ISLFD == 0) then
      do I = IL1,IL2 ! loop 450
        if (FG(I) > 0.) then
          FACTM = ZDSLM(I) + ZOM(I)
          FACTH = ZDSLH(I) + ZOM(I)
          RATIOM = SQRT(CDMX(I)) * LOG(FACTM / ZOM(I)) / VKC
          RATIOM = MIN(RATIOM,1.)
          RATIOH = SQRT(CDMX(I)) * LOG(FACTH / ZOH(I)) / VKC
          RATIOH = MIN(RATIOH,1.)
          if (RIBX(I) < 0.) then
            RATIOH = RATIOH * CDHX(I) / CDMX(I)
            RATIOH = MIN(RATIOH,(FACTH / ZRSLDH(I)) ** (1. / 3.))
          end if
          STT(I) = TSURX(I) - (MIN(RATIOH,1.)) * (TSURX(I) - TA(I))
          SQT(I) = QSURX(I) - (MIN(RATIOH,1.)) * (QSURX(I) - QA(I))
          SUT(I) = RATIOM * UWIND(I)
          SVT(I) = RATIOM * VWIND(I)
        end if
      end do ! loop 450
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FG, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FG(I) > 0.) then
          ST (I) = ST (I) + FG(I) * STT(I)
          SQ (I) = SQ (I) + FG(I) * SQT(I)
          SU (I) = SU (I) + FG(I) * SUT(I)
          SV (I) = SV (I) + FG(I) * SVT(I)
          SRH(I) = SRH(I) + FG(I) * SHT(I)
          SFCUBS (I) = SUT(I)
          SFCVBS (I) = SVT(I)
          USTARBS(I) = VA(I) * SQRT(CDMX(I))
        end if
      end do
    else if (ISLFD == 1) then
      call sldiag1(SUT, SVT, STT, SQT, &
                  CDMX, CDHX, UWIND, VWIND, TPOTA, QA, &
                  TSURX, QSURX, ZOM, ZOH, FG, ZRSLDM, &
                  ZDSLM, ZDSLH, ILG, IL1, IL2, JL)
      !
      call screenRelativeHumidity(SHT, STT, SQT, PRESSG, FG, ILG, IL1, IL2) ! Formerly SCREENRH
      !
      do I = IL1,IL2
        if (FG(I) > 0.) then
          ST (I) = ST (I) + FG(I) * STT(I)
          SQ (I) = SQ (I) + FG(I) * SQT(I)
          SU (I) = SU (I) + FG(I) * SUT(I)
          SV (I) = SV (I) + FG(I) * SVT(I)
          SRH(I) = SRH(I) + FG(I) * SHT(I)
          SFCUBS (I) = SUT(I)
          SFCVBS (I) = SVT(I)
          USTARBS(I) = VA(I) * SQRT(CDMX(I))
        end if
      end do
    else if (ISLFD == 2) then
      call DIASURFZ1(SU, SV, ST, SQ, ILG, UWIND, VWIND, TSURX, QSURX, &
                    ZOM, ZOH, ILMOX, ZRSLFM, HBLX, UEX, FTEMPX, FVAPX, &
                    ZDSLM, ZDSLH, RADJ, FG, IL1, IL2, JL)
    end if
    !
    do I = IL1,IL2 ! loop 475
      if (FG(I) > 0.) then
        EVPPOT(I) = EVPPOT(I) + FG(I) * RHOAIR(I) * CFLUX(I) * &
                    (Q0SAT(I) - QA(I))
        ACOND(I) = ACOND(I) + FG(I) * CFLUX(I)
        ILMO(I) = ILMO(I) + FG(I) * ILMOX(I)
        UE(I)   = UE(I) + FG(I) * UEX(I)
        HBL(I)  = HBL(I) + FG(I) * HBLX(I)
        CDH (I) = CDH(I) + FG(I) * CDHX(I)
        CDM (I) = CDM(I) + FG(I) * CDMX(I)
        TSFSAV(I,4) = TSURX(I)
        GTBS(I) = TSURX(I)
        QG(I) = QG(I) + FG(I) * QSURX(I)
        QSENS(I) = QSENS(I) + FG(I) * QSENSX(I)
        QEVAP(I) = QEVAP(I) + FG(I) * QEVAPX(I)
        QLWAVG(I) = QLWAVG(I) + FG(I) * QLWX(I)
        FSGG(I) = FSGG(I) + FG(I) * QSWX(I)
        FLGG(I) = FLGG(I) + FG(I) * (QLWIN(I) - QLWX(I))
        HFSG(I) = HFSG(I) + FG(I) * QSENSX(I)
        HEVG(I) = HEVG(I) + FG(I) * QEVAPX(I)
        HTC(I,1) = HTC(I,1) + FG(I) * ( - G12G(I))
        HTC(I,2) = HTC(I,2) + FG(I) * (G12G(I) - G23G(I))
        HTC(I,3) = HTC(I,3) + FG(I) * G23G(I)
        FTEMP(I) = FTEMP(I) + FG(I) * FTEMPX(I)
        FVAP (I) = FVAP (I) + FG(I) * FVAPX (I)
        RIB  (I) = RIB  (I) + FG(I) * RIBX  (I)
      end if
    end do ! loop 475
  end if
  !
  !     * ADDITIONAL DIAGNOSTIC VARIABLES.
  !
  do I = IL1,IL2 ! loop 500
    GT(I) = (QLWAVG(I) / SBC) ** 0.25
    TFLUX(I) = - QSENS(I) / (RHOAIR(I) * SPHAIR)
    EVAP(I) = EVAP(I) + RHOW * &
              (FCS(I) * (EVAPCS(I) + EVPCSG(I)) + FGS(I) * EVAPGS(I) + &
              FC (I) * (EVAPC (I) + EVAPCG(I)) + FG (I) * EVAPG(I))
    if (EVPPOT(I) /= 0.0) then
      EVAPB(I) = EVAP(I) / EVPPOT(I)
    else
      EVAPB(I) = 0.0
    end if
    if ((FCS(I) + FC(I)) > 1.0E-5) then
      TAC(I) = (FCS(I) * TACCS(I) + FC(I) * TACCO(I)) / (FCS(I) + FC(I))
      QAC(I) = (FCS(I) * QACCS(I) + FC(I) * QACCO(I)) / (FCS(I) + FC(I))
    else
      TAC(I) = TA(I)
      QAC(I) = QA(I)
    end if
    ! Calculate the ground heat flux across the sub-regions.
    groundHeatFlux(I) = FC(I) * GZEROC(I) & !vegetated 
                       + FG(I) * GZEROG(I) & !bare ground 
                       + FCS(I) * GZROCS(I) & ! snow covered vegetated
                       + FGS(I) * GZROGS(I) ! bare ground snow covered.
  end do ! loop 500

  return
end subroutine energyBudgetDriver
