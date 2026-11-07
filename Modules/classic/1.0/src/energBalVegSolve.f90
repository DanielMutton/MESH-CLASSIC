!> \file
!> Solves surface energy balance for vegetated subareas.
!! @author D. Verseghy, M. Lazare, A. Wu, P. Bartlett, Y. Delage, V. Arora, E. Chan, J. Melton, Y. Wu
!
subroutine energBalVegSolve (ISNOW, FI, & ! Formerly TSOLVC
                             QSWNET, QSWNC, QSWNG, QLWOUT, QLWOC, QLWOG, QTRANS, &
                             QSENS, QSENSC, QSENSG, QEVAP, QEVAPC, QEVAPG, EVAPC, &
                             EVAPG, EVAP, TCAN, QCAN, TZERO, QZERO, GZERO, QMELTC, &
                             QMELTG, RAICAN, SNOCAN, CDH, CDM, RIB, TAC, QAC, &
                             CFLUX, FTEMP, FVAP, ILMO, UE, H, QFCF, QFCL, HTCC, &
                             QSWINV, QSWINI, QLWIN, TPOTA, TA, QA, VA, VAC, PADRY, &
                             RHOAIR, ALVISC, ALNIRC, ALVISG, ALNIRG, TRVISC, TRNIRC, &
                             FSVF, CRIB, CPHCHC, CPHCHG, CEVAP, TADP, TVIRTA, RC, &
                             RBCOEF, ZOSCLH, ZOSCLM, ZRSLFH, ZRSLFM, ZOH, ZOM, &
                             FCOR, GCONST, GCOEFF, TGND, TRSNOW, FSNOWC, FRAINC, &
                             CHCAP, CMASS, PCPR, FROOT, THLMIN, DELZW, RHOSNO, ZSNOW, &
                             IWATER, IEVAP, ITERCT, &
                             ISLFD, ITC, ITCG, ILG, IL1, IL2, JL, N, &
                             TSTEP, TVIRTC, TVIRTG, EVBETA, XEVAP, EVPWET, Q0SAT, &
                             RA, RB, RAGINV, RBINV, RBTINV, RBCINV, TVRTAC, &
                             TPOTG, RESID, &
                             TCANO, WZERO, XEVAPM, DCFLXM, WC, DRAGIN, CFLUXM, CFLX, &
                             IEVAPC, TRTOP, QSTOR, CFSENS, CFEVAP, QSGADD, A, B, &
                             LZZ0, LZZ0T, FM, FH, ITER, NITER, KF1, KF2, &
                             AILCG, FCANC, CO2CONC, RMATCTEM, &
                             THLIQ, THFC, THLW, ISAND, IG, COSZS, PRESSG, &
                             XDIFFUS, ICTEM, IC, CO2I1, CO2I2, &
                             ctem_on, SLAI, FCANCMX, L2MAX, &
                             NOL2PFTS, CFLUXV, ANVEG, RMLVEG, &
                             DAYL, DAYL_MAX, ipeatland, Cmossmas, dmoss, &
                             anmoss, rmlmoss, iday, pdd)
  !
  !     * OCT 30/16 - J.Melton    Finish implementation of peatlands by Yuanqiao Wu
  !     * AUG 30/16 - J.Melton    Replace ICTEMMOD with ctem_on (logical switch).
  !     * JUL 22/15 - D.VERSEGHY. LIMIT CALCULATED EVAPORATION RATES
  !     *                         ACCORDING TO WATER AVAILABILITY.
  !     * FEB 27/15 - J. MELTON - WILTSM AND FIELDSM ARE RENAMED THLW AND THFC, RESPECTIVELY.
  !     * JUN 27/14 - D.VERSEGHY. CHANGE ITERATION LIMIT BACK TO 50 FOR
  !     *                         BISECTION SCHEME; BUGFIX IN CALCULATION
  !     *                         OF EVPWET.
  !     * OCT 30/12 - V. ARORA  - CFLUXV WAS BEING INITIALIZED TO ZERO INAPPROPRIATELY
  !     *                         FOR MOSAIC RUNS. NOT A PROBLEM WITH COMPOSITE
  !     *                         RUNS. CREATED A TEMPORARY STORAGE VAR TO ALLOW
  !     *                         AN APPROPRIATE VALUE FOR THE INITIALIZATION
  !     * SEP 05/12 - J.MELTON. - MADE AN IMPLICIT INT TO REAL CONVERSION
  !     *                         EXPLICIT. ALSO PROBLEM WITH TAC, SEE NOTE
  !     *                         IN THE CODE BEFORE CALL TO PHTSYN.
  !     * NOV 11/11 - M.LAZARE. - INCORPORATES CTEM. THIS INVOLVES
  !     *                         SEVERAL CHANGES AND NEW OUTPUT ROUTINES.
  !     *                         QSWNVC IS PROMOTED TO A WORK ARRAY
  !     *                         SINCE PASSED AS INPUT TO THE NEW CALLED
  !     *                         PHOTOSYNTHESIS ROUTINE "photosynCanopyConduct". THE
  !     *                         CTEM CANOPY RESISTANCE COMING OUT OF
  !     *                         THIS ROUTINE, "RCPHTSYN" IS STORED INTO
  !     *                         THE USUAL "RC" ARRAY AS LONG AS THE
  !     *                         BONE-DRY SOIL FLAG IS NOT SET (RC=1.E20).
  !     *                         WE ALSO HAVE TO PASS "TA" THROUGH FROM
  !     *                         energyBudgetDriver. FINALLY, "ISAND", "FIELDSM" AND
  !     *                         "WILTSM" ARE PASSED THROUGH TO photosynCanopyConduct
  !     *                         FOR CTEM.
  !     * OCT 14/11 - D.VERSEGHY. FOR POST-ITERATION CLEANUP WITH N-R SCHEME,
  !     *                         REMOVE CONDITION INVOLVING LAST ITERATION
  !     *                         TEMPERATURE.
  !     * DEC 07/09 - D.VERSEGHY. RESTORE EVAPOTRANSPIRATION WHEN
  !     *                         PRECIPITATION IS OCCURRING; ADD EVAPC
  !     *                         TO EVAP WHEN DEPOSITION OF WATER ON
  !     *                         CANOPY IS OCCURRING.
  !     * MAR 13/09 - D.VERSEGHY. REPLACE COMMON BLOCK SURFCON WITH CLASSD2
  !     *                         REVISED CALL TO FLXSURFZ.
  !     * JAN 20/09 - D.VERSEGHY. CORRECT CALCULATION OF TPOTG.
  !     * JAN 06/09 - E.CHAN/D.VERSEGHY. SET UPPER LIMIT FOR TSTEP IN
  !     *                         N-R ITERATION SCHEME.
  !     * FEB 26/08 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS; REMOVE
  !     *                         "ILW" SWITCH; SUPPRESS WATER VAPOUR FLUX
  !     *                         IF PRECIPITATION IS OCCURRING.
  !     * FEB 19/07 - D.VERSEGHY. UPDATE CANOPY WATER STORES IN THIS
  !     *                         ROUTINE INSTEAD OF canopyInterception FOR CASES
  !     *                         OF WATER DEPOSITION.
  !     * MAY 17/06 - D.VERSEGHY. ADD IL1 AND IL2 TO CALL TO FLXSURFZ
  !     *                         REMOVE JL FROM CALL TO DRCOEF.
  !     * APR 15/05 - D.VERSEGHY. SUBLIMATION OF INTERCEPTED SNOW TAKES
  !     *                         PLACE BEFORE EVAPORATION OF INTERCEPTED
  !     *                         RAIN.
  !     * APR 14/05 - Y.DELAGE.   REFINEMENTS TO N-R ITERATION SCHEME.
  !     * FEB 23/05 - D.VERSEGHY. INCORPORATE A SWITCH TO USE EITHER THE
  !     *                         BISECTION ITERATION SCHEME WITH CANOPY
  !     *                         AIR PARAMETRIZATION, OR THE NEWTON-
  !     *                         RAPHSON ITERATION SCHEME WITH MODIFIED
  !     *                         ZOH.
  !     * JAN 31/05 - Y.DELAGE.   USE THE CANOPY AIR RESISTANCE TO CALCULATE A
  !     *                         ROUGHNESS LENGTH FOR TEMPERATURE AND HUMIDITY.
  !     *                         REPLACE SECANT METHOD BY NEWTON-RAPHSON SCHEME
  !     *                         FOR BOTH ITERATION LOOPS.
  !     *                         LIMIT NUMBER OF ITERATIONS (ITERMX) TO 5 AND
  !     *                         APPLY CORRECTIONS IF RESIDUE REMAINS.
  !     * JAN 12/05 - P.BARTLETT/D.VERSEGHY. MODIFICATION TO CALCULATION
  !     *                         OF RBINV; ALLOW SUBLIMATION OF FROZEN
  !     *                         WATER ONLY ONTO SNOW-COVERED PORTION
  !     *                         OF CANOPY.
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 06/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS
  !     *                         FROM ROUNDOFF ERRORS.
  !     * NOV 07/02 - Y.DELAGE/D.VERSEGHY. NEW CALL TO FLXSURFZ; VIRTUAL
  !     *                         AND POTENTIAL TEMPERATURE CORRECTIONS.
  !     * NOV 01/02 - P.BARTLETT. MODIFICATIONS TO CALCULATIONS OF QAC
  !     *                         AND RB.
  !     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
  !     * MAR 28/02 - D.VERSEGHY. STREAMLINED SUBROUTINE CALL.
  !     * MAR 10/02 - M.LAZARE.   VECTORIZE LOOP 650 BY SPLITTING INTO TWO.
  !     * JAN 18/02 - P.BARTLETT/D.VERSEGHY. NEW "BETA" FORMULATION FOR
  !     *                         BARE SOIL EVAPORATION BASED ON LEE AND
  !     *                         PIELKE.
  !     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
  !     * OCT 06/00 - D.VERSEGHY. CONDITIONAL "IF" IN ITERATION SEQUENCE
  !     *                         TO AVOID DIVIDE BY ZERO.
  !     * DEC 16/99 - A.WU/D.VERSEGHY. REVISED CANOPY TURBULENT FLUX
  !     *                              FORMULATION: ADD PARAMETRIZATION
  !     *                              OF CANOPY AIR TEMPERATURE.
  !     * DEC 07/99 - A.WU/D.VERSEGHY. NEW SOIL EVAPORATION FORMULATION.
  !     * JUL 24/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         REPLACE BISECTION METHOD IN SURFACE
  !     *                         TEMPERATURE ITERATION SCHEME WITH
  !     *                         SECANT METHOD FOR FIRST TEN ITERATIONS.
  !     *                         PASS QZERO, QA, ZOMS, ZOHS TO REVISED
  !     *                         DRCOEF (ZOMS AND ZOHS ALSO NEW WORK ARRAYS
  !     *                         PASSED TO THIS ROUTINE).
  !     * JUN 20/97 - D.VERSEGHY. PASS IN NEW "CLASS4" COMMON BLOCK.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.  ALSO, PASS SWITCH "ILW"
  !     *                         THROUGH SUBROUTINE CALL, SPECIFYING
  !     *                         WHETHER QLWIN REPRESENTS INCOMING
  !     *                         (ILW=1) OR NET (ILW=2) LONGWAVE
  !     *                         RADIATION ABOVE THE GROUND.
  !     * NOV 30/94 - M.LAZARE.   CLASS - VERSION 2.3.
  !     *                         NEW DRAG COEFFICIENT AND RELATED FIELDS,
  !     *                         NOW DETERMINED IN ROUTINE "DRCOEF".
  !     * OCT 04/94 - D.VERSEGHY. CHANGE "CALL ABORT" TO "CALL errorHandler" TO
  !     *                         ENABLE RUNNING ON PCS.
  !     * JAN 24/94 - M.LAZARE.   UNFORMATTED I/O COMMENTED OUT IN LOOPS
  !     *                         200 AND 600.
  !     * JUL 29/93 - D.VERSEGHY. CLASS - VERSION 2.2.
  !     *                         ADD TRANSMISSION THROUGH SNOWPACK TO
  !     *                         "QSWNET" FOR DIAGNOSTIC PURPOSES.
  !     * OCT 15/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. ITERATIVE TEMPERATURE CALCULATIONS
  !     *                         FOR VEGETATION CANOPY AND UNDERLYING
  !     *                         SURFACE.
  !
  use peatlandsMod,        only : mossPht
  use classicParams,       only : DELT, TFREZ, GRAV, SBC, SPHW, SPHICE, &
                                  SPHVEG, SPHAIR, RHOW, CLHMLT, &
                                  CLHVAP, DELTA, BETA
  use generalutilsGAT,        only : calcEsat
  
  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  ! input variable
  integer, intent(in) :: ISNOW !< Flag indicating presence or absence of snow
  integer, intent(in) :: ISLFD, ITC, ITCG, ILG, IL1, IL2, JL, N
  integer :: I, J
  !
  integer :: NUMIT, IBAD, NIT, ITERMX
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out)   :: QSWNET(ILG) !< Total net shortwave radiation of canopy and underlying surface \f$[W m^{-2} ]\f$
  real, intent(inout) :: QSWNC (ILG) !< Net shortwave radiation on vegetation canopy \f$[W m^{-2} ] (K_{*c} )\f$
  real, intent(inout) :: QSWNG (ILG) !< Net shortwave radiation at underlying surface \f$[W m^{-2} ] (K_{*g} )\f$
  real, intent(out)   :: QLWOUT(ILG) !< Upwelling longwave radiation from canopy and underlying surface \f$[W m^{-2} ]\f$
  real, intent(inout) :: QLWOC (ILG) !< Upwelling longwave radiation from vegetation canopy \f$[W m^{-2} ] (L \uparrow_c, L \downarrow_c)\f$
  real, intent(inout) :: QLWOG (ILG) !< Upwelling longwave radiation from underlying surface \f$[W m^{-2} ] (L \uparrow_g)\f$
  real, intent(inout) :: QTRANS(ILG) !< Shortwave radiation transmitted into surface \f$[W m^{-2} ] \f$
  real, intent(out)   :: QSENS (ILG) !< Sensible heat flux from canopy and underlying surface \f$[W m^{-2} ] (Q_H)\f$
  real, intent(inout) :: QSENSC(ILG) !< Sensible heat flux from vegetation canopy \f$[W m^{-2} ] (Q_{H,c} )\f$
  real, intent(inout) :: QSENSG(ILG) !< Sensible heat flux from underlying surface \f$[W m^{-2} ] (Q_{H,g} )\f$
  real, intent(inout) :: QEVAP (ILG) !< Latent heat flux from canopy and underlying surface \f$[W m^{-2} ] (Q_E)\f$
  real, intent(inout) :: QEVAPC(ILG) !< Latent heat flux from vegetation canopy \f$[W m^{-2} ] (Q_{E,c} )\f$
  real, intent(inout) :: QEVAPG(ILG) !< Latent heat flux from underlying surface \f$[W m^{-2} ] (Q_{E,g} )\f$
  real, intent(inout) :: EVAPC (ILG) !< Evaporation rate from vegetation \f$[kg m^{-2} s^{-1} ] (E_c)\f$
  real, intent(inout) :: EVAPG (ILG) !< Evaporation rate from underlying surface \f$[kg m^{-2} s^{-1} ] (E(0))\f$
  real, intent(inout) :: TCAN  (ILG) !< Vegetation canopy temperature \f$[K] (T_c)\f$
  real, intent(inout) :: QCAN  (ILG) !< Saturated specific humidity at canopy temperature \f$[kg kg^{-1} ] (q_c)\f$
  real, intent(inout) :: TZERO (ILG) !< Temperature at surface \f$[K] (T(0))\f$
  real, intent(inout) :: QZERO (ILG) !< Specific humidity at surface \f$[kg kg^{-1} ] (q(0))\f$
  real, intent(inout) :: GZERO (ILG) !< Heat flux into surface \f$[W m^{-2} ] (G(0))\f$
  real, intent(inout) :: QMELTC(ILG) !< Heat available for melting snow or freezing water on the vegetation \f$[W m^{-2} ]\f$
  real, intent(out)   :: QMELTG(ILG) !< Heat available for melting snow or freezing water on the underlying surface \f$[W m^{-2} ]\f$
  real, intent(inout) :: RAICAN(ILG) !< Intercepted liquid water stored on vegetation canopy \f$[kg m^{-2} ]\f$
  real, intent(inout) :: SNOCAN(ILG) !< Intercepted frozen water stored on vegetation canopy \f$[kg m^{-2} ]\f$
  real, intent(in)    :: CDH   (ILG) !< Surface drag coefficient for heat [ ]
  real, intent(in)    :: CDM   (ILG) !< Surface drag coefficient for momentum [ ]
  real, intent(in)    :: RIB   (ILG) !< Bulk Richardson number at surface [ ]
  real, intent(inout) :: TAC   (ILG) !< Temperature of air within vegetation canopy \f$[K] (T_{ac} )\f$
  real, intent(inout) :: QAC   (ILG) !< Specific humidity of air within vegetation canopy space \f$[kg kg^{-1} ] (q_{ac} )\f$
  real, intent(in)    :: CFLUX (ILG) !< Product of surface drag coefficient and wind speed \f$[m s^{-1} ]\f$
  real, intent(in)    :: FTEMP (ILG) !< Product of surface-air temperature gradient,
  !< drag coefficient and wind speed \f$[K m s^{-1} ]\f$
  real, intent(in)    :: FVAP  (ILG) !< Product of surface-air humidity gradient, drag coefficient
  !< and wind speed \f$[kg kg^{-1} m s^{-1} ]\f$
  real, intent(in)    :: ILMO  (ILG) !< Inverse of Monin-Obukhov roughness length \f$(m^{-1} ]\f$
  real, intent(in)    :: UE    (ILG) !< Friction velocity of air \f$[m s^{-1} ]\f$
  real, intent(in)    :: H     (ILG) !< Height of the atmospheric boundary layer [m]
  real, intent(inout) :: QFCF  (ILG) !< Sublimation from frozen water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(inout) :: QFCL  (ILG) !< Evaporation from liquid water on vegetation \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(inout) :: HTCC  (ILG) !< Internal energy change of canopy due to changes in temperature and/or mass \f$[W m^{-2} ]\f$
  real, intent(inout) :: EVAP  (ILG) !< Diagnosed total surface water vapour flux over modelled area \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(inout)    :: anmoss(ilg) !<
  real, intent(inout)    :: rmlmoss(ilg)!<
  real   :: cevapmoss(ilg)!<
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI    (ILG) !< Fractional coverage of subarea in question on modelled area [ ]
  real, intent(in) :: QSWINV(ILG) !< Visible radiation incident on horizontal surface \f$[W m^{-2} ] (K \downarrow)\f$
  real, intent(in) :: QSWINI(ILG) !< Near-infrared radiation incident on horizontal surface \f$[W m^{-2} ] (K \downarrow)\f$
  real, intent(in) :: QLWIN (ILG) !< Downwelling longwave radiation at bottom of atmosphere \f$[W m^{-2} ] (L \downarrow)\f$
  real, intent(in) :: TPOTA (ILG) !< Potential temperature of air at reference height \f$[K] (T_{a,pot} )\f$
  real, intent(in) :: TA    (ILG) !< Air temperature at reference height [K]
  real, intent(in) :: QA    (ILG) !< Specific humidity at reference height \f$[kg kg^{-1} ] (q_a)\f$
  real, intent(in) :: VA    (ILG) !< Wind speed at reference height \f$[m s^{-1} ]\f$
  real, intent(in) :: VAC   (ILG) !< Wind speed within vegetation canopy space \f$[m s^{-1} ] (v_{ac} )\f$
  real, intent(in) :: PADRY (ILG) !< Partial pressure of dry air \f$[Pa] (p_{dry} )\f$
  real, intent(in) :: RHOAIR(ILG) !< Density of air \f$[kg m^{-3} ] ( \rho_a)\f$
  real, intent(in) :: ALVISC(ILG) !< Visible albedo of vegetation canopy \f$[ ] ( \alpha_c)\f$
  real, intent(in) :: ALNIRC(ILG) !< Near-IR albedo of vegetation canopy \f$[ ] ( \alpha_c)\f$
  real, intent(in) :: ALVISG(ILG) !< Visible albedo of underlying surface \f$[ ] ( \alpha_g)\f$
  real, intent(in) :: ALNIRG(ILG) !< Near-IR albedo of underlying surface \f$[ ] ( \alpha_g)\f$
  real, intent(in) :: TRVISC(ILG) !< Visible transmissivity of vegetation canopy \f$[ ] ( \tau_c)\f$
  real, intent(in) :: TRNIRC(ILG) !< Near-IR transmissivity of vegetation canopy \f$[ ] ( \tau_c)\f$
  real, intent(in) :: FSVF  (ILG) !< Sky view factor of ground underlying canopy \f$[ ] ( \chi)\f$
  real, intent(in) :: CRIB  (ILG) !< Richardson number coefficient \f$[K^{-1} ]\f$
  real, intent(inout) :: CPHCHC(ILG) !< Latent heat of vaporization on vegetation canopy \f$[J kg^{-1} ] (L_v)\f$
  real, intent(in) :: CPHCHG(ILG) !< Latent heat of vaporization on underlying ground \f$[J kg^{-1} ] (L_v)\f$
  real, intent(in) :: CEVAP (ILG) !< Soil evaporation efficiency coefficient [ ]
  real, intent(in) :: TADP  (ILG) !< Dew point of air at reference height [K]
  real, intent(in) :: TVIRTA(ILG) !< Virtual potential temperature of air at reference height \f$[K] (T_{a,v} )\f$
  real, intent(inout) :: RC    (ILG) !< Stomatal resistance of vegetation \f$[s m^{-1}] (r_c)\f$
  real, intent(in) :: RBCOEF(ILG) !< Parameter for calculation of leaf boundary resistance
  real, intent(in) :: ZOSCLH(ILG) !< Ratio of roughness length for heat to reference height
  !< for temperature and humidity [ ]
  real, intent(in) :: ZOSCLM(ILG) !< Ratio of roughness length for momentum to reference height for wind speed [ ]
  real, intent(in) :: ZRSLFH(ILG) !< Difference between reference height for temperature and humidity and
  !< height at which extrapolated wind speed goes to zero [m]
  real, intent(in) :: ZRSLFM(ILG) !< Difference between reference height for wind speed and height at which
  !< extrapolated wind speed goes to zero [m]
  real, intent(in) :: ZOH   (ILG) !< Surface roughness length for heat [m]
  real, intent(in) :: ZOM   (ILG) !< Surface roughness length for momentum \f$[m] (z_{0,m} )\f$
  real, intent(in) :: FCOR  (ILG) !< Coriolis parameter \f$[s_{-1} ]\f$
  real, intent(in) :: GCONST(ILG) !< Intercept used in equation relating surface heat flux
  !< to surface temperature \f$[W m^{-2} ]\f$
  real, intent(in) :: GCOEFF(ILG) !< Multiplier used in equation relating surface heat flux
  !< to surface temperature \f$[W m^{-2} K^{-1} ]\f$
  real, intent(in) :: TGND  (ILG) !< Starting point for surface temperature iteration [K]
  real, intent(in) :: TRSNOW(ILG) !< Short-wave transmissivity of snow pack [ ]
  real, intent(inout) :: FSNOWC(ILG) !< Fractional coverage of canopy by frozen water \f$[ ] (F_s)\f$
  real, intent(inout) :: FRAINC(ILG) !< Fractional coverage of canopy by liquid water \f$[ ] (F_l)\f$
  real, intent(inout) :: CHCAP (ILG) !< Heat capacity of vegetation canopy \f$[J m^{-2} K^{-1} ] (C_c)\f$

  real, intent(in) :: CMASS (ILG) !< Mass of vegetation canopy \f$[kg m^{-2} ]\f$
  real, intent(in) :: PCPR  (ILG) !< Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
  real, intent(in) :: RHOSNO(ILG) !< Density of snow  \f$[kg m^{-3} ]\f$
  real, intent(in) :: ZSNOW (ILG) !< Depth of snow pack  [m] \f$(z_s)\f$
  !
  real, intent(inout) :: FROOT (ILG,IG) !< Fraction of total transpiration contributed by soil layer  [  ]
  real, intent(in) :: THLMIN(ILG,IG) !< Residual soil liquid water content remaining after freezing or evaporation \f$[m^{3} m^{-3} ]\f$
  real, intent(in) :: THLIQ(ILG,IG)  !< Volumetric liquid water content of soil layers  \f$[m^{3} m^{-3} ]\f$
  real, intent(in) :: DELZW(ILG,IG) !< Permeable thickness of soil layer  [m]
  !
  integer, intent(in) :: IWATER(ILG) !< Flag indicating condition of surface (dry, water-covered or snow-covered)
  integer, intent(inout) :: IEVAP (ILG) !< Flag indicating whether surface evaporation is occurring or not output variable
  integer, intent(inout) :: ITERCT(ILG,6,50) !< Counter of number of iterations required to solve energy balance for four subareas
  integer, intent(in)  :: ipeatland(ilg)!<
  integer :: ievapmoss(ilg)!<
  integer, intent(in) :: iday          !<
  real :: thmin(ilg,ig)    !<
  real, intent(in) :: Cmossmas(ilg)    !<
  real, intent(in) :: dmoss(ilg)       !<
  real, intent(inout) :: pdd(ilg)         !<

  !
  !     * ARRAYS FOR CTEM.
  !
  !     * AILCG    - GREEN LAI FOR CARBON PURPOSES
  !     * FCANC    - FRACTIONAL COVERAGE OF 8 CARBON PFTs
  !     * CO2CONC  - ATMOS. CO2 CONC. IN PPM
  !     * RMATCTEM - FRACTION OF ROOTS IN EACH SOIL LAYER FOR EACH OF THE 8 PFTs
  !                  FOR CARBON RELATED PURPOSES.
  !     * RCPHTSYN - STOMATAL RESISTANCE ESTIMATED BY THE PHTSYN SUBROUTINE, S/M
  !     * COSZS    - COS OF SUN'S ZENITH ANGLE
  !     * XDIFFUS  - FRACTION OF DIFFUSED RADIATION
  !     * CO2I1    - INTERCELLULAR CO2 CONC.(PA) FOR THE SINGLE/SUNLIT LEAF
  !     * CO2I2    - INTERCELLULAR CO2 CONC.(PA) FOR THE SHADED LEAF
  !     * CTEM1    - LOGICAL BOOLEAN FOR USING CTEM's STOMATAL RESISTANCE
  !                  OR NOT
  !     * CTEM2    - LOGICAL BOOLEAN FOR USING CTEM's STRUCTURAL ATTRIBUTES
  !                  OR NOT
  !     * SLAI     - STORAGE LAI. SEE PHTSYN SUBROUTINE FOR MORE DETAILS.
  !     * FCANCMX  - MAX. FRACTIONAL COVERAGE OF CTEM PFTs
  !     * L2MAX    - MAX. NUMBER OF LEVEL 2 CTEM PFTs
  !     * NOL2PFTS - NUMBER OF LEVEL 2 CTEM PFTs
  !     * ANVEG    - NET PHTOSYNTHETIC RATE, u-MOL/M^2/S, FOR CTEM's 8 PFTs
  !     * RMLVEG   - LEAF MAINTENANCE RESP. RATE, u-MOL/M^2/S, FOR CTEM's 8 PFTs
  !
  real, intent(in) :: AILCG(ILG,ICTEM), FCANC(ILG,ICTEM), CO2CONC(ILG), &
                      CO2I1(ILG,ICTEM), CO2I2(ILG,ICTEM), COSZS(ILG), &
                      SLAI(ILG,ICTEM), PRESSG(ILG), XDIFFUS(ILG), &
                      RMATCTEM(ILG,ICTEM,IG), FCANCMX(ILG,ICTEM), &
                      ANVEG(ILG,ICTEM), RMLVEG(ILG,ICTEM), &
                      THFC(ILG,IG), THLW(ILG,IG)
  real, intent(inout) :: CFLUXV(ILG)
  real :: CFLUXV_IN(ILG)

  real, intent(in) :: DAYL_MAX(ILG)      ! MAXIMUM DAYLENGTH FOR THAT LOCATION
  real, intent(in) :: DAYL(ILG)          ! DAYLENGTH FOR THAT LOCATION


  integer, intent(in) :: ISAND(ILG,IG)
  !
  logical, intent(in) :: ctem_on

  integer, intent(in) :: ICTEM, L2MAX, NOL2PFTS(IC), IC, IG
  !
  !     * LOCAL WORK ARRAYS FOR CTEM.
  !
  real :: RCPHTSYN(ILG), QSWNVC(ILG)
  !
  !     * GENERAL INTERNAL WORK ARRAYS.
  !
  real, intent(inout) :: TSTEP (ILG), TVIRTC(ILG), TVIRTG(ILG), &
                         EVBETA(ILG), XEVAP (ILG), EVPWET(ILG), Q0SAT (ILG), &
                         RA    (ILG), RB    (ILG), RAGINV(ILG), RBINV (ILG), &
                         RBTINV(ILG), RBCINV(ILG), TVRTAC(ILG), &
                         TPOTG (ILG), RESID (ILG), TCANO (ILG), &
                         TRTOP (ILG), QSTOR (ILG), A     (ILG), B     (ILG), &
                         LZZ0  (ILG), LZZ0T (ILG), &
                         FM    (ILG), FH    (ILG), WZERO (ILG), XEVAPM(ILG), &
                         DCFLXM(ILG), WC    (ILG), DRAGIN(ILG), CFLUXM(ILG), &
                         CFSENS(ILG), CFEVAP(ILG), QSGADD(ILG), CFLX  (ILG)

  !
  real :: WAVAIL(ILG,IG), WROOT(ILG,IG), WTRTOT(ILG), EVPMAX(ILG)
  !
  integer, intent(inout) :: ITER(ILG), NITER(ILG), IEVAPC(ILG), &
                            KF1(ILG), KF2(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: qswnvg(ilg), &
          QSWNIG, QSWNIC, HFREZ, HCONV, &
          RCONV, HCOOL, HMELT, SCONV, HWARM, WCAN, DQ0DT, &
          DRDT0, QEVAPT, BOWEN, DCFLUX, DXEVAP, TCANT, QEVAPCT, &
          TZEROT, YEVAP, RAGCO, EZERO, WTRANSP, WTEST
  !
  !-----------------------------------------------------------------------
  !     * INITIALIZATION AND PRE-ITERATION SEQUENCE.
  !===================== CTEM =====================================\
  !
  do I = 1,ILG
    QSWNVC(I) = 0.0
    !    initialize QSWNVG to be used in mossPht.f90 ---
    if (ipeatland(i) > 0) then
      qswnvg(i) = 0.0
    end if
  end do

  !===================== CTEM =====================================/
  !
  if (ITCG < 2) then
    ITERMX = 50
  else
    ITERMX = 12      ! was 5 YW March 27, 2015
  end if
  !      IF (ISNOW==0) THEN
  !          EZERO=0.0
  !      ELSE
  !          EZERO=2.0
  !      END IF
  EZERO = 0.0
  RAGCO = 1.9E-3
  !
  !>
  !! For the subcanopy surface temperature iteration, two alternative schemes are offered: the bisection
  !! method (selected if the flag ITCG = 1) and the Newton-Raphson method (selected if ITCG = 2). In the
  !! first case, the maximum number of iterations ITERMX is set to 12, and in the second case it is set to 5.
  !! An optional windless transfer coefficient EZERO is made available, which can be used, following the
  !! recommendations of Brown et al. (2006) \cite Brown2006-ec, to prevent the sensible heat flux over snow packs from
  !! becoming vanishingly small under highly stable conditions. If the snow cover flag ISNOW is zero
  !! (indicating bare ground), EZERO is set to zero; if ISNOW=1, EZERO is set to \f$2.0 W m^{-2} K^{-1}\f$ . The
  !! surface transfer coefficient under conditions of free convection, RAGCO, is set to \f$1.9 x 10^{-3}\f$ (see the
  !! calculation of RAGINV below).
  !!
  !! In the 50 loop, some preliminary calculations are done. The shortwave transmissivity at the surface,
  !! TRTOP, is set to zero in the absence of a snow pack, and to the transmissivity of snow, TRSNOW,
  !! otherwise. The net shortwave radiation at the surface, QSWNG, is calculated as the sum of the net visible
  !! and net near-infrared shortwave radiation. Each of these is obtained as:
  !! \f$K_{*g} = K \downarrow \tau_c [1 - \alpha_g ]\f$
  !! where \f$K_{*g} is the net radiation at the surface, \f$K \downarrow\f$ is the incoming shortwave radiation above the canopy, \f$\tau_c\f$
  !! is the canopy transmissivity and \f$\alpha_g\f$ is the surface albedo. This average value is corrected for the amount
  !! of radiation transmitted into the surface, QTRANS, obtained using TRTOP. The net shortwave radiation
  !! for the vegetation canopy, QSWNC, is calculated as the sum of the net visible and net near-infrared
  !! shortwave radiation. Each of these is determined as:
  !! \f$K_{*c} = K \downarrow [1 - \alpha_c ] - K_{*g}\f$
  !! where \f$K_{*c}\f$ is the net radiation on the canopy and \f$\alpha_c\f$ is the canopy albedo. If the canopy temperature is
  !! essentially 0 K, indicating that a canopy was not present in the previous time step but has now appeared,
  !! the canopy temperature is initialized to the potential temperature of the air, TPOTA. The outgoing
  !! longwave radiation emitted upward \f$(L \uparrow_c)\f$ or downward \f$(L \downarrow_c)\f$ from the canopy is calculated using the
  !! standard Stefan-Boltzmann equation:
  !! \f$L \uparrow_c = L \downarrow_c = \sigma T_c^4\f$
  !! where \f$\sigma\f$ is the Stefan-Boltzmann constant and \f$T_c\f$ is the canopy temperature.
  !!
  !! Virtual temperature is defined as temperature adjusted for the reduction in air density due to the presence
  !! of water vapour. This is applied in order to enable the use of the equation of state for dry air. The virtual
  !! temperature can be approximated by multiplying the actual temperature by a factor [1 + 0.61 q], where q
  !! is the specific humidity of the air in question. Thus, the virtual temperature of air at the vegetation
  !! canopy temperature, \f$T_{c, v}\f$ , is obtained as:
  !! \f$T_{c, v} = T_c [1 + 0.61 q_c ]\f$
  !! where \f$q_c\f$ is the saturated specific humidity at the canopy temperature. This is determined from the
  !! saturation mixing ratio at the canopy temperature, \f$w_c\f$ :
  !! \f$q_c = w_c /[1 + w_c ]\f$
  !! The saturation mixing ratio is a function of the saturation vapour pressure \f$e_c\f$ at the canopy temperature:
  !! \f$w_c = 0.622 e_c /(p_{dry} )\f$
  !! where \f$p_{dry}\f$ is the partial pressure of dry air. A standard empirical equation for the saturation vapour
  !! pressure dependence on the temperature T is used following Emanuel (1994) \cite Emanuel1994-dt
  !! \f$e_{sat} = exp[53.67957 - 6743.769 / T - 4.8451 * ln(T)]       T \geq T_f\f$
  !!
  !! \f$e_{sat} = exp[23.33086 - 6111.72784 / T + 0.15215 * log(T)]    T < T_f\f$
  !!
  !! where \f$T_f\f$ is the freezing point. The virtual temperature of the air in the canopy space, \f$T_{ac, v}\f$ , is likewise
  !! calculated from the canopy air temperature \f$T_{ac}\f$ and the specific humidity in the canopy air space, \f$q_{ac}\f$ , as
  !!
  !! \f$T_{ac, v} = T_{ac} [1 + 0.61 q_{ac} ]\f$
  !!
  !! If the Newton-Raphson method is being used for the canopy temperature iteration (pre-selected by
  !! setting the flag ITC to 2), the temperature of the air in the canopy space is approximated as the canopy
  !! temperature, and the specific humidity in the canopy air space is approximated as the specific humidity of
  !! the air above the canopy.
  !!
  !! If there is intercepted snow on the vegetation, the latent heat associated with water flux from the canopy,
  !! CPHCHC, is set to the latent heat of sublimation (the sum of the heat of vaporization and the heat of
  !! melting). Otherwise the latent heat is set to that of vaporization alone. The leaf boundary layer resistance
  !! RB, and its inverse RBINV, are calculated using the wind speed in the canopy air space, VAC, and a
  !! coefficient RBCOEF evaluated in subroutine calcLandSurfParams. This coefficient is formulated after Bartlett (2004),
  !! who developed an expression for the inverse of the leaf boundary resistance, \f$1/r_b\f$ , drawing on the analysis
  !! of Bonan (1996) \cite Bonan1996-as and McNaughton and van den Hurk (1995), of the form:
  !!
  !! \f$1/r_b = v_{ac}^{1/2} \sigma f_i \gamma_i \Lambda_i^{1/2} /0.75 [1 - exp(-0.75 \Lambda_i^{1/2})]\f$
  !!
  !! where \f$v_{ac}\f$ is the wind speed in the canopy air space, \f$f_i\f$ is the fractional coverage of each vegetation type i
  !! over the subarea in question, \f$\Lambda_i\f$ is its leaf area index, and \f$\gamma_i\f$ is a vegetation-dependent parameter which
  !! incorporates the effects of leaf dimension and sheltering. The initial value of the surface temperature
  !! TZERO is set to TGND, which contains the value of TZERO from the previous time step, and the
  !! initial temperature of the canopy, TCANO, is set to the current canopy temperature. The first step in the
  !! iteration sequence, TSTEP, is set to 1.0 K. The flag ITER is set to 1 for each element of the set of
  !! modelled areas, indicating that its surface temperature has not yet been found. The iteration counter
  !! NITER is initialized to 1 for each element. Initial values are assigned to other variables.In particular,
  !! the maximum evaporation from the ground surface that can be sustained for the current time step is
  !! specified as the total snow mass if snow is present, otherwise as the total available water in the first
  !! soil layer.
  !!
  !! (At this point in the code, if CLASS is being run in coupled mode with CTEM, the CTEM subroutine
  !! photosynCanopyConduct is called. These lines are commented out for uncoupled runs.)
  !!
  do I = IL1,IL2 ! loop 50
    if (FI(I) > 0.) then
      if (ISNOW == 0) then
        TRTOP(I) = 0.
      else
        TRTOP(I) = TRSNOW(I)
      end if

      !    calculate visible short wave radiation QSWNVG on the ground for moss
      !    photosynthesis ---------------------------------------------------\
      qswnvg(i) = QSWINV(I) * TRVISC(I) * (1.0 - ALVISG(I))
      QSWNIG = QSWINI(I) * TRNIRC(I) * (1.0 - ALNIRG(I))
      QSWNG(i) = qswnvg(i) + QSWNIG
      QTRANS(I) = QSWNG(I) * TRTOP(I)
      QSWNG(I) = QSWNG(I) - QTRANS(I)
      QSWNVC(I) = QSWINV(I) * (1.0 - ALVISC(I)) - qswnvg(i)
      !    QSWNVG is changed to qswnvg(i) YW March 20, 2015 -----------------/
      QSWNIC = QSWINI(I) * (1.0 - ALNIRC(I)) - QSWNIG
      QSWNC(I) = QSWNVC(I) + QSWNIC
      if (ABS(TCAN(I)) < 1.0E-3)        TCAN(I) = TPOTA(I)
      QLWOC(I) = SBC * TCAN(I) * TCAN(I) * TCAN(I) * TCAN(I)
      WCAN = 0.622 * calcEsat(TCAN(I)) / PADRY(I)
      QCAN(I) = WCAN / (1.0 + WCAN)
      TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
      if (ITC == 2) then
        TAC(I) = TCAN(I)
        QAC(I) = QA(I)
      end if
      TVRTAC(I) = TAC(I) * (1.0 + 0.61 * QAC(I))
      !
      if (SNOCAN(I) > 0.) then
        CPHCHC(I) = CLHVAP + CLHMLT
      else
        CPHCHC(I) = CLHVAP
      end if
      RBINV(I) = RBCOEF(I) * SQRT(VAC(I))
      if (RBINV(I) <= 0.) call errorHandler('energBalVegSolve',0)
      ! If RBINV is <= 0, it is possible your INI file is missing
      ! the PAIMIN and PAIMAX values.
      RB(I) = 1.0 / RBINV(I)
      TZERO(I) = TGND(I)
      TCANO(I) = TCAN(I)
      TSTEP(I) = 1.0
      ITER(I) = 1
      NITER(I) = 1
      QMELTC(I) = 0.0
      QMELTG(I) = 0.0
      if (ISNOW == 1) then
        KF1(I) = 1
        KF2(I) = 2
        EVPMAX(I) = RHOSNO(I) * ZSNOW(I) / DELT
      else
        KF1(I) = 4
        KF2(I) = 5
        EVPMAX(I) = RHOW * (THLIQ(I,1) - THLMIN(I,1)) * DELZW(I,1) / &
                    DELT
        EVPMAX(I) = MAX(EVPMAX(I),0.)
      end if
    end if
  end do ! loop 50
  !
  !     * CALL PHOTOSYNTHESIS SUBROUTINE HERE TO GET A NEW ESTIMATE OF
  !     * RC BASED ON PHOTOSYNTHESIS.
  !
  if (ctem_on) then

    !
    !       STORE CFLUXV NUMBERS IN A TEMPORARY ARRAY
    do I = IL1,IL2
      CFLUXV_IN(I) = CFLUXV(I)
    end do
    !
    call photosynCanopyConduct(AILCG, FCANC, TCAN, CO2CONC, & ! In ! Formerly PHTSYN3
                               PRESSG, FI, CFLUXV, QA, QSWNVC, IC, THLIQ, & ! In
                               ISAND, TA, RMATCTEM, COSZS, XDIFFUS, ILG, & ! In
                               IL1, IL2, IG, ICTEM, ISNOW, SLAI, & ! In
                               THFC, THLW, FCANCMX, L2MAX, NOL2PFTS, & ! In
                               CO2I1, CO2I2, & ! In/Out
                               RCPHTSYN, ANVEG, RMLVEG, DAYL, DAYL_MAX) ! Out
    !
    !       * KEEP CLASS RC FOR BONEDRY POINTS (DIANA'S FLAG OF 1.E20) SUCH
    !       * THAT WE GET (BALT-BEG) CONSERVATION.
    !
    do I = IL1,IL2
      RC(I) = MIN(RCPHTSYN(I),4999.999)
    end do ! loop 70

  end if

  !    Do moss photosynthesis:
  call  mossPht(il1, il2, iday, qswnvg, thliq, co2conc, tgnd, zsnow, &
                pressg, Cmossmas, dmoss, anmoss, rmlmoss, &
                cevapmoss, ievapmoss, ipeatland, DAYL, pdd)

  !
  !     * ITERATION FOR SURFACE TEMPERATURE OF GROUND UNDER CANOPY.
  !     * LOOP IS REPEATED UNTIL SOLUTIONS HAVE BEEN FOUND FOR ALL POINTS
  !     * ON THE CURRENT LATITUDE CIRCLE(S).
  !
  !>
  !! The 100 continuation line marks the beginning of the surface temperature iteration sequence. First the
  !! flag NUMIT (indicating whether there are still locations at the end of the current iteration step for which
  !! the surface temperature has not yet been found) is set to zero. Loop 125 is then performed over the
  !! vector of modelled areas. If ITER=1, the saturated specific humidity at the ground surface temperature,
  !! \f$q(0)_{sat}\f$ , is calculated using the same method as that outlined above for the saturated specific humidity at
  !! the canopy temperature. If there is a snow cover or ponded water present on the surface (IWATER > 0),
  !! the surface evaporation efficiency EVBETA is set to 1 and the surface specific humidity q(0) is set to
  !! \f$q(0)_{sat}\f$ . Otherwise EVBETA is set to CEVAP, the value obtained in subroutine energyBudgetPrep on the basis of
  !! ambient conditions, and q(0) is calculated from \f$q(0)_{sat}\f$ by making use of the definition of the surface
  !! evaporation efficiency \f$\beta\f$:
  !!
  !! \f$\beta = [q(0) - q_{ac} ]/[q(0)_{sat} - q_{ac} ]\f$
  !!
  !! which is inverted to obtain an expression for q(0). If \f$q(0) > q_{ac}\f$ and the evaporation flag IEVAP has been
  !! set to zero, EVBETA is reset to zero and q(0) is reset to \f$q_{ac}\f$.
  !!
  !! Next, the potential temperature of air at the ground surface temperature, \f$T(0)_{pot}\f$ , is calculated relative to
  !! the sum of the displacement height d and the roughness length for momentum \f$z_{0, m}\f$ , the height at which
  !! the canopy temperature is assumed to apply. (This is also the height corresponding to the potential
  !! temperature of the air at the reference height above the canopy). \f$T(0)_{pot}\f$ is found using the dry adiabatic
  !! lapse rate, \f$dT/dz = -g/c_p\f$ , where g is the acceleration due to gravity and \f$c_p\f$ is the specific heat at constant
  !! pressure. Since the displacement height d is assumed to lie at 0.7 of the canopy height, and the roughness
  !! length for momentum is assumed to lie at 0.1 of the canopy height, their sum can be obtained as \f$8.0 z_{0, m}\f$ .
  !! Thus,
  !! \f$T(0)_{pot} = T(0) - 8.0 z_{0, m} g/c_p\f$
  !!
  !! The virtual potential temperature at the surface is obtained using the same expression as above:
  !! \f$T(0)_v = T(0)_{pot} [1 + 0.61 q(0)]\f$
  !!
  !! Since wind speeds under the canopy are expected to be relatively small, it is assumed that under stable or
  !! neutral conditions, turbulent fluxes are negligible. If the virtual potential temperature of the surface is
  !! more than 1 K greater than that of the canopy air, free convection conditions are assumed. Townsend's
  !! (1964) equation for the surface-air transfer coefficient, or the inverse of the surface resistance \f$r_{a, , g}\f$ , is used
  !! in a form derived from the analysis of Deardorff (1972) \cite Deardorff1972-ay :
  !! \f$1/r_{a, , g} = 1.9 x 10^{-3} [T(0)_v - T_{ac, v} ]^{1/3}\f$
  !!
  !! The first derivative of the transfer coefficient with respect to the surface temperature is calculated for use
  !! with the Newton-Raphson iteration scheme:
  !! \f$d(1/r_{a, , g} )/dT = 1.9 x 10^{-3} [T(0)_v - T_{ac, v} ]^{-2/3} /3\f$
  !!
  !! If the virtual potential temperature of the surface is between 0.001 and 1 K greater than that of the
  !! canopy air, a simple diffusion relation is assumed:
  !! \f$1/r_{a, , g} = 1.9 x 10^{-3} [T(0)_v - T_{ac, v} ]\f$
  !! Thus,
  !! \f$d(1/r_{a, , g} )/dT = 1.9 x 10^{-3}\f$
  !!
  !! The remaining terms of the surface energy balance are now evaluated. The energy balance equation is
  !! expressed as:
  !! \f$K_{*g} + L_{*g} - Q_{H, g} - Q_{E, g} - G(0) = 0\f$
  !! where \f$K_{*g}\f$ is the net surface shortwave radiation, \f$L_{*g}\f$ is the net longwave radiation, \f$Q_{H, g}\f$ is the sensible heat
  !! flux, \f$Q_{E, g}\f$ is the latent heat flux, and G(0) is the conduction into the surface. \f$K_{*g}\f$ was evaluated earlier in
  !! loop 50. \f$L_{*g}\f$ is calculated as the difference between the downwelling radiation \f$L \downarrow_g\f$ at the surface and the
  !! upwelling radiation \f$L \uparrow_g\f$ . The downwelling radiation incident on the surface is determined from the
  !! downwelling sky radiation above the canopy \f$L \downarrow\f$ and the downwelling radiation from the canopy itself,
  !! \f$L \downarrow_c\f$ , weighted according to the sky view factor \f$\chi\f$:
  !! \f$L \downarrow_g = \chi L \downarrow + [1 - \chi] L \downarrow_c\f$
  !! The upwelling radiation is obtained using the Stefan-Boltzmann equation:
  !! \f$L \uparrow_g = \sigma T(0)^4\f$
  !!
  !! The sensible heat flux is given by
  !! \f$Q_{H, g} = \rho_a c_p [T(0)_{pot} - T_{a, c} ]/r_{a, , g}\f$
  !! where \f$\rho_a\f$ is the density of the air and \f$c_p\f$ is its specific heat. (The windless transfer coefficient defined at
  !! the beginning of the subroutine is currently not used under vegetation canopies.) The evaporation rate at
  !! the surface, E(0), is calculated as
  !! \f$E(0) = \rho_a [q(0) - q_{a, c} ]/r_{a, , g}\f$
  !!
  !! The evaporation rate is constrained to be less than or equal to the maximum rate evaluated earlier in the subroutine.
  !! \f$Q_E\f$ is obtained by multiplying E(0) by the latent heat of vaporization at the surface. The heat flux into the
  !! surface G(0) is determined as a linear function of T(0) (see documentation for subroutines soilHeatFluxPrep and
  !! snowHeatCond). It can be seen that each of the terms of the surface energy balance is a function of a single
  !! unknown, T(0) or TZERO. The residual RESID of the energy balance is now evaluated on the basis of
  !! the current estimation for TZERO. If the absolute value of RESID is less than \f$5.0 W m^{-2}\f$ , or if the
  !! absolute value of the iteration step TSTEP most recently used is less than 0.01 K, the surface temperature
  !! is deemed to have been found and ITER is set to 0. If the iteration counter NITER is equal to the
  !! maximum number and ITER is still 1, ITER is set to -1.
  !!
  !! In the following section, the iteration sequence is moved ahead a step. If ITCG = 1, the calculations for
  !! the bisection method of solution in loop 150 are performed, over each of the modelled areas for which
  !! ITER = 1. If NITER = 1 (indicating that this is the first step in the iteration), then if RESID > 0
  !! (indicating that the current value for TZERO had undershot the correct value), TZERO is incremented
  !! by 1 K; otherwise it is decremented by 1 K. If this is not the first step in the iteration, then if RESID >0
  !! and TSTEP < 0 (indicating that TZERO has undershot the correct value and the last temperature
  !! increment had been a negative one) or if RESID < 0 and TSTEP > 0 (indicating that TZERO has
  !! overshot the correct value and the last temperature increment had been a positive one), TSTEP is divided
  !! in half and its sign changed. TSTEP is then added to TZERO. The iteration counter NITER and the
  !! flag NUMIT are each incremented by one. Finally, if NUMIT > 0, the iteration cycle is repeated from
  !! line 100 on.
  !!
  !! If ITCG = 2, the calculations for the Newton-Raphson method of iteration in loop 175 are performed,
  !! over each of the modelled areas for which ITER = 1. In this approach, the value \f$x_{n+1}\f$ used at each
  !! iteration step is obtained from the value \f$x_n\f$ at the previous step as follows:
  !! \f$x_{n+1} = x_n - f(x_n)/f'(x_n)\f$
  !!
  !! Identifying x n with TZERO and \f$f(x_n)\f$ with the surface energy balance equation, it can be seen that the
  !! second term on the right-hand side corresponds to TSTEP; the numerator is equal to RESID and the
  !! denominator to the first derivative of the energy balance equation evaluated at TZERO, which in turn is
  !! equal to the sum of the derivatives of the individual terms:
  !! \f[
  !! d(L \uparrow_g)/dT = -4 \sigma T(0)^3 \f]
  !! \f[
  !! d(Q_{H, g} )/dT = \rho_a c_p {1/r_{a, , g} + [T(0)_{pot} - T_{ac} ] d(1/r_{a, , g} )/dT}
  !! \f]
  !! \f[
  !! d(Q_{E, g} )/dT = L_v \rho_a {1/r_{a, , g} dq(0)/dT+ [q(0) - q_{ac} ] d(1/r_{a, , g} )/dT}
  !! \f]
  !! and dG(0)/dT is equal to the coefficient multiplying TZERO in the equation for G(0). (\f$L_v\f$ is the latent
  !! heat of vaporization at the surface.) At the end of the calculations the iteration counter NITER and the
  !! flag NUMIT are each incremented by one, and upon exiting the loop, if NUMIT > 0, the iteration cycle is
  !! repeated from line 100 on.
  !!
  !ignoreLint(1)
100 continue ! anchor for a GO TO statement that needs to be removed
  !
  NUMIT = 0
  do I = IL1,IL2 ! loop 125
    if (FI(I) > 0. .and. ITER(I) == 1) then
      WZERO(I) = 0.622 * calcEsat(TZERO(I)) / PADRY(I)
      Q0SAT(I) = WZERO(I) / (1.0 + WZERO(I))
      if (IWATER(I) > 0) then
        EVBETA(I) = 1.0
        QZERO(I) = Q0SAT(I)
      else
        !    evaporation coefficient is moss-controlled for peatland--
        if (ipeatland(i) == 0) then
          EVBETA(I) = CEVAP(I)
        else
          ievap(i) = ievapmoss(i)
          evbeta(i) = cevapmoss(i)
        end if
        QZERO(I) = EVBETA(I) * Q0SAT(I) + (1.0 - EVBETA(I)) * QAC(I)
        if (QZERO(I) > QAC(I) .and. IEVAP(I) == 0) then
          EVBETA(I) = 0.0
          QZERO(I) = QAC(I)
        end if
      end if
      !
      TPOTG(I) = TZERO(I) - 8.0 * ZOM(I) * GRAV / SPHAIR
      TVIRTG(I) = TPOTG(I) * (1.0 + 0.61 * QZERO(I))
      if (TVIRTG(I) > TVRTAC(I) + 1.) then
        RAGINV(I) = RAGCO * (TVIRTG(I) - TVRTAC(I)) ** 0.333333
        DRAGIN(I) = 0.333 * RAGCO * (TVIRTG(I) - TVRTAC(I)) ** ( - .667)
      else if (TVIRTG(I) > (TVRTAC(I) + 0.001)) then
        RAGINV(I) = RAGCO * (TVIRTG(I) - TVRTAC(I))
        DRAGIN(I) = RAGCO
      else
        RAGINV(I) = 0.0
        DRAGIN(I) = 0.0
      end if
      !
      QLWOG(I) = SBC * TZERO(I) * TZERO(I) * TZERO(I) * TZERO(I)
      QSENSG(I) = RHOAIR(I) * SPHAIR * RAGINV(I) * &
                  (TPOTG(I) - TAC(I))
      EVAPG (I) = RHOAIR(I) * (QZERO(I) - QAC(I)) * RAGINV(I)
      if (EVAPG(I) > EVPMAX(I)) EVAPG(I) = EVPMAX(I)
      QEVAPG(I) = CPHCHG(I) * EVAPG(I)
      GZERO(I) = GCOEFF(I) * TZERO(I) + GCONST(I)
      RESID(I) = QSWNG(I) + FSVF(I) * QLWIN(I) + (1.0 - FSVF(I)) * &
                 QLWOC(I) - QLWOG(I) - QSENSG(I) - QEVAPG(I) - GZERO(I)
      if (ABS(RESID(I)) < 5.0)                     ITER(I) = 0
      if (ABS(TSTEP(I)) < 1.0E-2)                  ITER(I) = 0
      if (NITER(I) == ITERMX .and. ITER(I) == 1)    ITER(I) = - 1
    end if
  end do ! loop 125
  !
  if (ITCG < 2) then
    !
    !     * OPTION #1: BISECTION ITERATION METHOD.
    !
    do I = IL1,IL2 ! loop 150
      if (FI(I) > 0. .and. ITER(I) == 1) then
        if (NITER(I) == 1) then
          if (RESID(I) > 0.0) then
            TZERO(I) = TZERO(I) + TSTEP(I)
          else
            TZERO(I) = TZERO(I) - TSTEP(I)
          end if
        else
          if ((RESID(I) > 0. .and. TSTEP(I) < 0.) .or. &
              (RESID(I) < 0. .and. TSTEP(I) > 0.)) then
            TSTEP(I) = - TSTEP(I) / 2.0
          end if
          TZERO(I) = TZERO(I) + TSTEP(I)
        end if
        NITER(I) = NITER(I) + 1
        NUMIT = NUMIT + 1
      end if
    end do ! loop 150
    !
  else
    !
    !     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.
    !
    do I = IL1,IL2 ! loop 175
      if (FI(I) > 0. .and. ITER(I) == 1) then
        DQ0DT = - WZERO(I) * A(I) * (B(I) - TFREZ) / ((TZERO(I) - B(I)) * &
                (1.0 + WZERO(I))) ** 2 * EVBETA(I)
        DRDT0 = - 4.0 * SBC * TZERO(I) ** 3 &
                - GCOEFF(I) - RHOAIR(I) * SPHAIR * &
                (RAGINV(I) + (TPOTG(I) - TAC(I)) * DRAGIN(I)) - &
                CPHCHG(I) * RHOAIR(I) * (DQ0DT * RAGINV(I) &
                + (QZERO(I) - QAC(I)) * DRAGIN(I))
        TSTEP(I) = - RESID(I) / DRDT0
        if (ABS(TSTEP(I)) > 20.0) TSTEP(I) = SIGN(10.0,TSTEP(I))
        TZERO(I) = TZERO(I) + TSTEP(I)
        NITER(I) = NITER(I) + 1
        NUMIT = NUMIT + 1
      end if
    end do ! loop 175
    !
  end if
  !
  if (NUMIT > 0) GO TO 100
  !
  !     * IF CONVERGENCE HAS NOT BEEN REACHED FOR ITERATION METHOD #2,
  !     * CALCULATE TEMPERATURE AND FLUXES ASSUMING NEUTRAL STABILITY
  !     * AND USING BOWEN RATIO APPROACH.
  !
  if (ITCG == 2) then
    !
    !>
    !! After the iteration has been completed, if the Newton-Raphson method has been used, a check is carried
    !! out in loop 200 to ascertain whether convergence has not been reached (i.e. whether ITER = -1) for any
    !! location. In such cases it is assumed that conditions of near-neutral stability at the surface are the cause of
    !! the difficulty in finding a solution. A trial value of TZERO is calculated using the virtual potential
    !! temperature of the canopy. If the absolute value of RESID is \f$> 15 W m^{-2}\f$ , TZERO is set to this trial
    !! value. The values of q(0) and the components of the surface energy balance are recalculated as above,
    !! except that \f$Q_{H, g}\f$ and \f$Q_{E, g}\f$ are assumed as a first approximation to be zero. RESID is determined on this
    !! basis, and is then partitioned between \f$Q_{H, g}\f$ and \f$Q_{E, g}\f$ on the basis of the Bowen ratio B, the ratio of \f$Q_{H, g}\f$
    !! over \f$Q_{E, g}\f$ . Setting the residual R equal to \f$Q_{H, g} + Q_{E, g}\f$ , and substituting for \f$Q_{H, g}\f$ using \f$B = Q_{H, g} /Q_{E, g}\f$ ,
    !! results in:
    !! \f$Q_{E, g} = R/(1 + B)\f$
    !! \f$Q_{H, g}\f$ is then obtained as \f$R - Q_{E, g}\f$ , the residual is reset to zero, and E(0) is recalculated from \f$Q_{E, g}\f$ .
    !!
    !! At this point a check is performed for unphysical values of the surface temperature, i.e. for values greater
    !! than 100 C or less than -100 C. If such values are encountered, an error message is printed and a call to
    !! abort is carried out.
    !!
    do I = IL1,IL2
      if (ITER(I) == - 1) then
        TZEROT = TVIRTC(I) / (1.0 + 0.61 * QZERO(I))
        if (ABS(RESID(I)) > 15.) then
          TZERO(I) = TZEROT
          WZERO(I) = 0.622 * calcEsat(TZERO(I)) / PADRY(I)
          Q0SAT(I) = WZERO(I) / (1.0 + WZERO(I))
          QZERO(I) = EVBETA(I) * Q0SAT(I) + (1.0 - EVBETA(I)) * QAC(I)
          QLWOG(I) = SBC * TZERO(I) * TZERO(I) * TZERO(I) * TZERO(I)
          GZERO(I) = GCOEFF(I) * TZERO(I) + GCONST(I)
          RESID(I) = QSWNG(I) + FSVF(I) * QLWIN(I) + (1.0 - FSVF(I)) * &
                     QLWOC(I) - QLWOG(I) - GZERO(I)
          QEVAPT = CPHCHG(I) * (QZERO(I) - QAC(I))
          BOWEN = SPHAIR * (TZERO(I) - TAC(I)) / &
                  SIGN(MAX(ABS(QEVAPT),1.E-6),QEVAPT)
          QEVAPG(I) = RESID(I) / SIGN(MAX(ABS(1. + BOWEN),0.1),1. + BOWEN)
          QSENSG(I) = RESID(I) - QEVAPG(I)
          RESID(I) = 0.
          EVAPG(I) = QEVAPG(I) / CPHCHG(I)
        end if
      end if
    end do ! loop 200
    !
  end if
  !
  IBAD = 0
  !
  do I = IL1,IL2
    !          IF (FI(I)>0. .AND. ITER(I)==-1)                     THEN
    !              WRITE(6,6250) I,JL,NITER(I),RESID(I),TZERO(I),RIB(I)
    ! 6250          FORMAT('0waterUnderCanopy ITERATION LIMIT',3X,3I3,3(F8.2,E12.4))
    !          END IF
    if (FI(I) > 0.) then
      if (TZERO(I) < 173.16 .or. TZERO(I) > 373.16) then
        IBAD = I
      end if
    end if
  end do ! loop 225
  !
  if (IBAD /= 0) then
    write(6,6370) IBAD,N,TZERO(IBAD),NITER(IBAD),ISNOW
6370 format('0BAD GROUND ITERATION TEMPERATURE',3X,2I8,F16.2,2I4)
    write(6,6380) QSWNG(IBAD),FSVF(IBAD),QLWIN(IBAD),QLWOC(IBAD), &
         QLWOG(IBAD),QSENSG(IBAD),QEVAPG(IBAD),GZERO(IBAD)
    write(6,6380) TCAN(IBAD)
    call errorHandler('energBalVegSolve', - 1)
  end if
  !
  !     * POST-ITERATION CLEAN-UP.
  !
  !>
  !! Finally, clean-up calculations are performed in loop 250. A check is carried out to ensure that TZERO is
  !! not less than 0 C if ponded water is present on the surface (IWATER = 1) or greater than 0 C if snow is
  !! present on the surface (IWATER = 2). If either is the case, TZERO is reset to the freezing point, and
  !! q(0), \f$T(0)_{pot}\f$ and \f$T(0)_v\f$ are re-evaluated. The components of the surface energy balance are recalculated
  !! using the above equations. The residual of the energy balance equation is assigned to the energy
  !! associated with phase change of water at the surface, QMELTG, and RESID is set to zero.
  !!
  !! In the last part of the loop, some final adjustments are made to a few other variables. If the evaporation
  !! flux is vanishingly small, it is added to RESID and reset to zero. Any remaining RESID is added to
  !! \f$Q_{H,g}\f$ . Lastly, the iteration counter ITERCT is updated for the level corresponding to the subarea type and
  !! the value of NITER.
  !!
  do I = IL1,IL2
    if (FI(I) > 0.) then
      if ((IWATER(I) == 1 .and. TZERO(I) < TFREZ) .or. &
          (IWATER(I) == 2 .and. TZERO(I) > TFREZ)) then
        TZERO(I) = TFREZ
        WZERO(I) = 0.622 * 611.0 / PADRY(I)
        QZERO(I) = WZERO(I) / (1.0 + WZERO(I))
        TPOTG(I) = TZERO(I) - 8.0 * ZOM(I) * GRAV / SPHAIR
        TVIRTG(I) = TPOTG(I) * (1.0 + 0.61 * QZERO(I))
        !
        QLWOG(I) = SBC * TZERO(I) * TZERO(I) * TZERO(I) * TZERO(I)
        GZERO(I) = GCOEFF(I) * TZERO(I) + GCONST(I)
        if (TVIRTG(I) > (TVRTAC(I) + 0.001)) then
          RAGINV(I) = RAGCO * (TVIRTG(I) - TVRTAC(I)) ** 0.333333
          QSENSG(I) = RHOAIR(I) * SPHAIR * RAGINV(I) * &
                      (TPOTG(I) - TAC(I))
          EVAPG (I) = RHOAIR(I) * (QZERO(I) - QAC(I)) * RAGINV(I)
        else
          RAGINV(I) = 0.0
          QSENSG(I) = 0.0
          EVAPG (I) = 0.0
        end if
        if (EVAPG(I) > EVPMAX(I)) EVAPG(I) = EVPMAX(I)
        QEVAPG(I) = CPHCHG(I) * EVAPG(I)
        QMELTG(I) = QSWNG(I) + FSVF(I) * QLWIN(I) + (1.0 - FSVF(I)) * &
                    QLWOC(I) - QLWOG(I) - QSENSG(I) - QEVAPG(I) - GZERO(I)
        RESID(I) = 0.0
      end if
      !
      if (ABS(EVAPG(I)) < 1.0E-8) then
        RESID(I) = RESID(I) + QEVAPG(I)
        EVAPG(I) = 0.0
        QEVAPG(I) = 0.0
      end if
      !              IF (RESID(I)>15. .AND. QEVAPG(I)>10. .AND. PCPR(I)
      !     1                   <1.0E-8)                 THEN
      !                  QEVAPG(I)=QEVAPG(I)+RESID(I)
      !              ELSE
      QSENSG(I) = QSENSG(I) + RESID(I)
      !              END IF
      ITERCT(I,KF2(I),NITER(I)) = ITERCT(I,KF2(I),NITER(I)) + 1
    end if
  end do ! loop 250
  !
  !     * PRE-ITERATION SEQUENCE FOR VEGETATION CANOPY.
  !
  !>
  !! In the 300 loop, preliminary calculations are done in preparation for the canopy temperature iteration.
  !! The sensible heat flux from the surface is treated differently if ITC = 1 (bisection method of solution) and
  !! ITC = 2 (Newton-Raphson method of solution). In the first instance the surface sensible heat flux is
  !! applied to heating the air in the canopy space; in the second, the canopy and the air space within it are
  !! treated as one aggregated mass, and the sensible heat flux from below is assumed to be added to its energy
  !! balance. Thus, if ITC = 1, the sensible heat flux that is added to the canopy, QSGADD, is set to 0. If
  !! ITC = 2, it is set to \f$Q_{H, g}\f$ ; and \f$T_{ac}\f$ is set to \f$T_c\f$ , \f$q_{ac}\f$ to \f$q_c\f$ , and \f$T_{ac, v}\f$ to \f$T_{c, v}\f$ .
  !! The flag ITER is set to 1 for each
  !! element of the set of modelled areas, indicating that its surface temperature has not yet been found. The
  !! iteration counter NITER is initialized to 1 for each element. The first step in the iteration sequence,
  !! TSTEP, is set to 1.0 K. Initial values are assigned to other variables. In the following 350 loop, the water
  !! available for transpiration, WAVAIL, is evaluated for each soil layer.  Lastly, after exiting the loop, the
  !! maximum number of iterations ITERMX is set to 50 if ITC = 1, and to 5 if ITC = 2.
  !!
  do I = IL1,IL2
    if (FI(I) > 0.) then
      QSGADD(I) = 0.0
      if (ITC == 2) then
        QSGADD(I) = QSENSG(I)
        TAC(I) = TCAN(I)
        QAC(I) = QCAN(I)
        TVRTAC(I) = TVIRTC(I)
      end if
      ITER(I) = 1
      NITER(I) = 1
      TSTEP(I) = 1.0
      CFLUXM(I) = 0.0
      DCFLXM(I) = 0.0
      WTRTOT(I) = 0.0
    end if
  end do ! loop 300
  !
  do J = 1,IG ! loop 350
    do I = IL1,IL2
      if (FI(I) > 0.) then
        WAVAIL(I,J) = RHOW * (THLIQ(I,J) - THLMIN(I,J)) * DELZW(I,J)
        if (J == 1 .and. EVAPG(I) > 0.0) &
            WAVAIL(I,J) = WAVAIL(I,J) - EVAPG(I) * DELT
        WAVAIL(I,J) = MAX(WAVAIL(I,J),0.)
        WROOT(I,J) = 0.0
      end if
    end do
  end do ! loop 350
  !
  if (ITC < 2) then
    ITERMX = 50
  else
    ITERMX = 5
  end if
  !
  !     * ITERATION FOR CANOPY TEMPERATURE.
  !     * LOOP IS REPEATED UNTIL SOLUTIONS HAVE BEEN FOUND FOR ALL POINTS
  !     * ON THE CURRENT LATITUDE CIRCLE(S).
  !
  !>
  !! The 400 continuation line marks the beginning of the canopy temperature iteration sequence. First the
  !! flags NIT (indicating whether there are still locations at the beginning of the current iteration step for
  !! which the surface temperature has not yet been found) and NUMIT (indicating whether there are still
  !! locations at the end of the current iteration step for which the surface temperature has not yet been
  !! found) are set to zero. Loop 450 is then performed over the set of modelled areas. If ITER=1, NIT is
  !! incremented by one; if ITC = 1, the vegetation virtual temperature \f$T_{c, v}\f$ is recalculated.
  !!
  !! If NIT > 0, a subroutine is called to evaluate the stability-corrected surface drag coefficients for heat and
  !! momentum. The subroutine selection is made on the basis of the flag ISLFD. If ISLFD=1, indicating
  !! that the calculations are to be consistent with CCCma conventions, subroutine DRCOEF is called; if
  !! ISLFD=2, indicating that the calculations are to be consistent with RPN conventions, subroutine
  !! FLXSURFZ is called.
  !!
  !! Next, canopy parameters and turbulent transfer coefficients are calculated prior to evaluating the terms of
  !! the surface energy balance equation. If ITC = 1, the analysis of Garratt (1992) \cite Garratt1992-dt is followed. The sensible
  !! and latent heat fluxes from the canopy air to the overlying atmosphere, \f$Q_H\f$ and \f$Q_E\f$ , are obtained as the
  !! sums of the sensible and latent heat fluxes from the canopy to the canopy air, \f$Q_{H, c}\f$ and \f$Q_{E, c}\f$ , and from the
  !! underlying surface to the canopy air, \f$Q_{H, g}\f$ and \f$Q_{E, g}\f$ :
  !!
  !! \f$Q_H = Q_{H, c} + Q_{H, g}\f$
  !!
  !! \f$Q_E = Q_{E, c} + Q_{E, g}\f$
  !!
  !! where
  !!
  !! \f$Q_H = \rho_a c_p [T_{a, c} - T_{a, pot, } ]/r_a\f$
  !!
  !! \f$Q_{H, c} = \rho_a c_p [T_c - T_{a, c, } ]/r_b\f$
  !!
  !! and
  !!
  !! \f$Q_E = L_v \rho_a [q_{a, c} - q_{a, } ]/r_a\f$
  !!
  !! \f$Q_{E, c} = L_v \rho_a [q_{, c} - q_{a, c, } ]/(r_b + r_c)\f$
  !!
  !! The equations for the sensible and latent heat fluxes from the surface were presented above. In these
  !! expressions, \f$T_{a, pot}\f$ and \f$q_a\f$ are the potential temperature and specific humidity of the air overlying the
  !! canopy, \f$r_a\f$ is the aerodynamic resistance to turbulent transfer between the canopy air and the overlying air,
  !! and \f$r_c\f$ is the stomatal resistance to transpiration. (The term CFLUX that is generated by the subroutines
  !! DRCOEF and FLXSURFZ is equivalent to \f$1/r_a\f$ .) Thus, \f$T_{a, c}\f$ and \f$q_{a, c}\f$ can be evaluated as
  !!
  !! \f$T_{a, c} = [T_{a, pot} /r_a + T_c /r_b + T(0)_{pot} /r_{a, , g} ]/[1/r_a + 1/r_b + 1/r_{a, , g} ]\f$
  !!
  !! \f$q_{a, c} = [q_a /r_a + q_c /(r_b + r_c) + q(0)/r_{a, , g} ]/[1/r_a + 1/(r_b + r_c) + 1/r_{a, , g} ]\f$
  !!
  !! If the water vapour flux is towards the canopy leaves, or if the canopy is snow-covered or rain-covered, \f$r_c\f$
  !! is zero. If the water vapour flux is away from the canopy leaves and the fractional snow or water
  !! coverage on the canopy, \f$F_s\f$ or \f$F_l\f$ , is less than 1, the term \f$1/(r_b + r_c)\f$ adjusted for the presence of
  !! intercepted snow or water, \f$X_E\f$ (XEVAP in the code) is calculated as a weighted average, as follows. If \f$F_s > 0\f$,
  !! the canopy must be at a temperature of 0 C or less, and so it is deduced that no transpiration can be
  !! occurring. Thus,
  !! \f$X_E = (F_s + F_l)/r_b\f$ .
  !!
  !! Otherwise, \f$X_E\f$ is calculated on the assumption that the resistance is equal to \f$r_b\f$ over the water-covered
  !! area and \f$r_b + r_c\f$ over the rest:
  !! \f$X_E = F_l /r_b + [1 - F_l ]/[r_b + r_c ]\f$
  !!
  !ignoreLint(1)
400 continue ! anchor for a GO TO statement that needs to be removed

  NUMIT = 0
  NIT = 0
  !>
  !! In the 450 loop, XEVAP is first set as a trial value to RBINV (neglecting stomatal resistance), and QAC is
  !! calculated using this value. If \f$q_{a, c} < q_c\f$ , indicating that the vapour flux is away from the canopy leaves,
  !! XEVAP is recalculated as above and QAC is re-evaluated. Otherwise, if \f$F_s > 0\f$, XEVAP is scaled by \f$F_s\f$ ,
  !! since it is assumed that snow on the canopy will be at a lower temperature than the canopy itself, and
  !! deposition via sublimation will occur preferentially onto it. \f$T_{a, c}\f$ and \f$T_{ac, v}\f$ are calculated as above, and the
  !! canopy-air turbulent transfer coefficients for sensible and latent heat flux, CFSENS and CFEVAP, are set
  !! to RBINV and XEVAP respectively.
  !!
  !! If ITC = 2, the canopy parameters and turbulent transfer coefficients are calculated, as noted above, on
  !! the basis of the assumption that the vegetation canopy and the air space within it can be treated as one
  !! aggregated mass, with a single representative temperature. Thus, the resistances \f$r_a\f$ and \f$r_b\f$ are considered as
  !! acting in series upon the sensible and latent heat fluxes between the canopy and the overlying air:
  !!
  !! \f$Q_{H, c} = \rho_a c_p [T_{, c} - T_{a, pot, } ]/(r_a + r_b)\f$
  !!
  !! \f$Q_{E, c} = L_v \rho_a [q_{, c} - q_{a, } ]/(r_a + r_b + r_c)\f$
  !!
  do I = IL1,IL2
    if (FI(I) > 0. .and. ITER(I) == 1) then
      NIT = NIT + 1
      if (ITC == 1) then
        WCAN = 0.622 * calcEsat(TCAN(I)) / PADRY(I)
        QCAN(I) = WCAN / (1.0 + WCAN)
        TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
      end if
    end if
  end do ! loop 450
  !
  if (NIT > 0) then
    !
    !     * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT)
    !     * AND OTHER RELATED QUANTITIES BETWEEN CANOPY AIR SPACE AND
    !     * ATMOSPHERE.
    !
    if (ISLFD < 2) then
      call drcoef1(CDM, CDH, RIB, CFLUX, QAC, QA, ZOSCLM, ZOSCLH, &
                  CRIB, TVRTAC, TVIRTA, VA, FI, ITER, &
                  ILG, IL1, IL2)
    else
      call flxsurfz1(CDM, CDH, CFLUX, RIB, FTEMP, FVAP, ILMO, &
                    UE, FCOR, TPOTA, QA, ZRSLFM, ZRSLFH, VA, &
                    TAC, QAC, H, ZOM, ZOH, &
                    LZZ0, LZZ0T, FM, FH, ILG, IL1, IL2, FI, ITER, JL)
    end if
    !
    !     * CALCULATE CANOPY AIR TEMPERATURE AND SPECIFIC HUMIDITY OF
    !     * CANOPY AIR (FIRST WITHOUT RC TO CHECK FOR CONDENSATION
    !     * IF NO CONDENSATION EXISTS, RECALCULATE).
    !
    if (ITC == 1) then
      !
      do I = IL1,IL2
        if (FI(I) > 0. .and. ITER(I) == 1) then
          XEVAP(I) = RBINV(I)
          QAC(I) = (QCAN(I) * XEVAP(I) + QZERO(I) * RAGINV(I) + &
                   QA(I) * CFLUX(I)) / (XEVAP(I) + RAGINV(I) + CFLUX(I))
          if (QAC(I) < QCAN(I)) then
            if (FSNOWC(I) > 0.0) then
              XEVAP(I) = (FRAINC(I) + FSNOWC(I)) / RB(I)
            else
              XEVAP(I) = FRAINC(I) / RB(I) + (1.0 - FRAINC(I)) / &
                         (RB(I) + RC(I))
              QAC(I) = (QCAN(I) * XEVAP(I) + QZERO(I) * RAGINV(I) + &
                       QA(I) * CFLUX(I)) / (XEVAP(I) + RAGINV(I) + &
                       CFLUX(I))
            end if
          else
            if (FSNOWC(I) > 1.0E-5) then
              XEVAP(I) = FSNOWC(I) / RB(I)
            else
              XEVAP(I) = 1.0 / RB(I)
            end if
          end if
          TAC(I) = (TCAN(I) * RBINV(I) + TPOTG(I) * RAGINV(I) + &
                   TPOTA(I) * CFLUX(I)) / (RBINV(I) + RAGINV(I) + CFLUX(I))
          TVRTAC(I) = TAC(I) * (1.0 + 0.61 * QAC(I))
          CFSENS(I) = RBINV(I)
          CFEVAP(I) = XEVAP(I)
        end if
      end do ! loop 475
      !
    else
      !
      !>
      !! In loop 500 the term \f$1/(r_a + r_b)\f$ is calculated from the variables CFLUX (the inverse of \f$r_a\f$ ) and RBINV
      !! (the inverse of \f$r_b\f$ ), and assigned to the temporary variable CFLX. If the incoming visible shortwave
      !! radiation QSWINV is greater than or equal to \f$25 W m^{-2}\f$ , the calculated value of CFLX is retained; if
      !! QSWINV is zero, it is reset to CFLUX; and between these two limits it varies linearly between the two.
      !!
      !! Thus the effect of the calculated leaf boundary resistance is suppressed during conditions of zero or low
      !! solar heating. (This is done to avoid unrealistically low calculated turbulent fluxes at night, which can lead
      !! to anomalously low canopy temperatures.) The overall aerodynamic resistance \f$r_A = r_a + r_b\f$ is obtained as
      !! 1/CFLX and assigned to the variable RA. As with ITC = 1, if \f$q_a < q_c\f$ , XEVAP is recalculated as above
      !! (except that \f$r_A\f$ is substituted for \f$r_b\f$ ). If \f$F_s > 0\f$, XEVAP is again scaled by \f$F_s\f$ .
      !!
      !! In the approach used here, the specific humidity of the aggregated canopy \f$q_{0, c}\f$ is not assumed to be equal
      !! to the saturated specific humidity at the canopy temperature, but is rather determined using \f$X_E\f$ . If the
      !! two methods of calculating \f$Q_{E, c}\f$ are assumed to be analogous:
      !!
      !! \f$Q_{E, c} = L_v \rho_a X_E [q_{, c} - q_{a, } ]\f$
      !!
      !! \f$Q_{E, c} = L_v \rho_a [q_{, 0, c} - q_{a, } ]/r_A\f$
      !!
      !! then solving for \f$q_{0, c}\f$ leads to the expression
      !!
      !! \f$q_{0, c} = r_A X_E q_c + [1 - r_A X_E ]q_a\f$
      !!
      !! In the second part of the 500 loop the saturated specific humidity of the canopy is calculated as before
      !! and adjusted using the above equation to obtain the specific humidity of the aggregated canopy. This is
      !! then used to calculate \f$T_{c, v}\f$ . The canopy-air turbulent transfer coefficients for sensible and latent heat flux,
      !! CFSENS and CFEVAP, are both set to CFLX, and for calculation purposes in the following loop \f$T_{a, c}\f$ is
      !! set to the potential temperature of the air above the canopy, and \f$q_{a, c}\f$ to the specific humidity of the air
      !! above the canopy.
      !!
      do I = IL1,IL2
        if (FI(I) > 0. .and. ITER(I) == 1) then
          CFLX(I) = RBINV(I) * CFLUX(I) / (RBINV(I) + CFLUX(I))
          CFLX(I) = CFLUX(I) + (CFLX(I) - CFLUX(I)) * &
                    MIN(1.0,QSWINV(I) * 0.04)
          RA(I) = 1.0 / CFLX(I)
          if (QA(I) < QCAN(I)) then
            if (FSNOWC(I) > 0.0) then
              XEVAP(I) = (FRAINC(I) + FSNOWC(I)) / RA(I)
            else
              XEVAP(I) = FRAINC(I) / RA(I) + (1.0 - FRAINC(I)) / &
                         (RA(I) + RC(I))
            end if
          else
            if (FSNOWC(I) > 1.0E-5) then
              XEVAP(I) = FSNOWC(I) / RA(I)
            else
              XEVAP(I) = 1.0 / RA(I)
            end if
          end if
          WCAN = 0.622 * calcEsat(TCAN(I)) / PADRY(I)
          WC(I) = WCAN
          QCAN(I) = WCAN / (1.0 + WCAN)
          QCAN(I) = RA(I) * XEVAP(I) * QCAN(I) + (1.0 - RA(I) * XEVAP(I)) * &
                    QA(I)
          TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
          CFSENS(I) = CFLX(I)
          CFEVAP(I) = CFLX(I)
          TAC(I) = TPOTA(I)
          QAC(I) = QA(I)
        end if
      end do ! loop 500
      !
    end if
    !
    !     * CALCULATE THE TERMS IN THE ENERGY BALANCE AND SOLVE.
    !
    !>
    !! In the 525 loop, the terms of the canopy energy balance are evaluated. The energy balance is expressed
    !! as:
    !!
    !! \f$K_{*c} + L_{*c} - Q_{H, c} + Q_{H, g+} - Q_{E, c} - \Delta Q_{S, c} = 0\f$
    !!
    !! (The term \f$Q_{H, g+}\f$ corresponds to the variable QSGADD discussed above; the term \f$\Delta Q_{S, c}\f$ represents the
    !! change of energy storage in the canopy.) The net shortwave radiation \f$K_{*c}\f$ was evaluated earlier in the 50
    !! loop. The net longwave radiation is obtained as:
    !!
    !! \f$L_{*c} = (1 - \chi) [L \downarrow + L \uparrow_g - 2 L \downarrow_c ]\f$
    !!
    !! \f$Q_{H, c}\f$ is calculated as above. If there is intercepted liquid or frozen water on the canopy, or if the vapour
    !! flux is downward, or if the stomatal resistance is less than a limiting high value of \f$5000 s m^{-1}\f$ , the canopy
    !! vapour flux \f$E_c\f$ (that is, \f$Q_{E, c} /L_v\f$ ) is calculated as above and the flag IEVAPC is set to 1; otherwise \f$E_c\f$ and
    !! IEVAPC are both set to zero and \f$q_c\f$ is set to \f$q_a\f$ . If the water vapour flux is towards the canopy and the
    !! canopy temperature is greater than the air dew point temperature, the flux is set to zero. In the following IF block,
    !! checks are done to ensure that the calculated canopy evapotranspiration flux will not exceed the available water.
    !! If intercepted snow is present on the canopy, the maximum evapotranspiration is set to this amount.  Otherwise,
    !! if the calculated evapotranspiration rate exceeds the liquid water stored on the canopy, the rate is tested to
    !! ensure that the additional transpiration does not exhaust the available soil water weighted according to the root
    !! distribution. \f$Q_{E, c}\f$ is calculated from \f$E_c\f$ and \f$\Delta Q_{S, c}\f$ is
    !! obtained as
    !!
    !! \f$\Delta Q_{S, c} = C_c [T_c - T_{c, o} ]/ \Delta t\f$
    !!
    !! where \f$C_c\f$ is the canopy heat capacity, \f$T_{c, o}\f$ is the canopy temperature from the previous time step and \f$\Delta t\f$ is
    !! the length of the time step. The residual RESID of the energy balance is evaluated on the basis of the
    !! current estimation for the canopy temperature TCAN. If the absolute value of RESID is less than \f$5.0 W
    !! m^{-2}\f$ , or if the absolute value of the iteration step TSTEP most recently used is less than 0.01 K, the surface
    !! temperature is deemed to have been found and ITER is set to 0. If the iteration counter NITER is equal
    !! to the maximum number and ITER is still 1, ITER is set to -1.
    !!
    !! In the following section, the iteration sequence is moved ahead a step. If ITC = 1, the calculations for the
    !! bisection method of solution in loop 550 are performed, over each of the modelled areas for which ITER
    !! = 1. If NITER = 1 (indicating that this is the first step in the iteration), then if RESID > 0 (indicating
    !! that the current value for TCAN had undershot the correct value), TCAN is incremented by 1 K
    !! otherwise it is decremented by 1 K. If this is not the first step in the iteration, then if RESID >0 and
    !! TSTEP < 0 (indicating that TCAN has undershot the correct value and the last temperature increment
    !! had been a negative one) or if RESID < 0 and TSTEP > 0 (indicating that TCAN has overshot the
    !! correct value and the last temperature increment had been a positive one), TSTEP is divided in half and
    !! its sign changed. TSTEP is then added to TCAN. If TCAN is vanishingly close to 0 C, it is reset to that
    !! value. The iteration counter NITER and the flag NUMIT are each incremented by one. Finally, if
    !! NUMIT > 0, the iteration cycle is repeated from line 400 on.
    !!
    !! If ITC = 2, the calculations for the Newton-Raphson method of iteration in loop 575 are performed, over
    !! each of the modelled areas for which ITER = 1. As outlined above, in this approach the value \f$x_{n+1}\f$ used
    !! at each iteration step is obtained from the value \f$x_n\f$ at the previous step as follows:
    !!
    !! \f$x_{n+1} = x_n - f(x_n)/f'(x_n)\f$
    !!
    !! Identifying \f$x_n\f$ with TCAN and \f$f(x_n)\f$ with the surface energy balance equation, it can be seen that the
    !! second term on the right-hand side corresponds to TSTEP; the numerator is equal to RESID and the
    !! denominator to the first derivative of the energy balance equation evaluated at TCAN, which in turn is
    !! equal to the sum of the derivatives of the individual terms:
    !!
    !! \f$d(L_{*c} )/dT = -8 \sigma T_c^3 (1 - \chi)\f$
    !!
    !! \f$d(Q_{H, c} )/dT = \rho_a c_p {1/r_A + [T_c - T_{a, pot} ] d(1/r_A)/dT}\f$
    !!
    !! \f$d(Q_{E, c} )/dT = L_v \rho_a {X_E dq_c /dT + [q_c - q_a ] dX_E /dT}\f$
    !!
    !! \f$d \Delta Q_{S, c} /dT = C_c / \Delta t\f$
    !!
    !! The term \f$d(1/r_A)/dT\f$ is represented by the variable DCFLXM, which is approximated as the difference
    !! between CFLX and its value for the previous iteration, CFLUXM, divided by TSTEP. The term \f$dX_E /dT\f$
    !! is represented by the variable DXEVAP, which is approximated as the difference between XEVAP and
    !! its value for the previous iteration, XEVAPM, divided by TSTEP. The calculated value of TSTEP
    !! obtained from the above calculations is constrained to be between -10 and 5 K to dampen any spurious
    !! oscillations, and is then added to TCAN. If the resulting value of TCAN is vanishingly close to 0 C, it is
    !! reset to that value. At the end of the calculations the iteration counter NITER and the flag NUMIT are
    !! each incremented by one. The values of \f$T_{a, c}\f$ , \f$q_{a, c}\f$ and \f$T_{ac, v}\f$ are reset to \f$T_c\f$ , \f$q_c\f$ and \f$T_{c, v}\f$ respectively. Upon
    !! exiting the loop, if NUMIT > 0, the iteration cycle is repeated from line 400 on.
    !!
    do I = IL1,IL2
      if (FI(I) > 0. .and. ITER(I) == 1) then
        QLWOC(I) = SBC * TCAN(I) * TCAN(I) * TCAN(I) * TCAN(I)
        QSENSC(I) = RHOAIR(I) * SPHAIR * CFSENS(I) * (TCAN(I) - TAC(I))
        if (FRAINC(I) > 0. .or. FSNOWC(I) > 0. .or. &
            RC(I) <= 5000. .or. QAC(I) > QCAN(I)) then
          EVAPC(I) = RHOAIR(I) * CFEVAP(I) * (QCAN(I) - QAC(I))
          IEVAPC(I) = 1
        else
          EVAPC(I) = 0.0
          IEVAPC(I) = 0
          QCAN(I) = QA(I)
        end if
        if (EVAPC(I) < 0. .and. TCAN(I) > TADP(I)) EVAPC(I) = 0.0
        if (SNOCAN(I) > 0.) then
          EVPWET(I) = SNOCAN(I) / DELT
          if (EVAPC(I) > EVPWET(I)) EVAPC(I) = EVPWET(I)
        else
          EVPWET(I) = RAICAN(I) / DELT
          if (EVAPC(I) > EVPWET(I)) then
            WTRANSP = (EVAPC(I) - EVPWET(I)) * DELT
            EVPMAX(I) = EVPWET(I)
            WTRTOT(I) = 0.0
            do J = 1,IG
              WTEST = WTRANSP * FROOT(I,J)
              WROOT(I,J) = MIN(WTEST,WAVAIL(I,J))
              WTRTOT(I) = WTRTOT(I) + WROOT(I,J)
              EVPMAX(I) = EVPMAX(I) + WROOT(I,J) / DELT
            end do
            if (EVAPC(I) > EVPMAX(I)) EVAPC(I) = EVPMAX(I)
          end if
        end if
        QEVAPC(I) = CPHCHC(I) * EVAPC(I)
        QSTOR (I) = CHCAP(I) * (TCAN(I) - TCANO(I)) / DELT
        RESID(I) = QSWNC(I) + (QLWIN(I) + QLWOG(I) - 2.0 * QLWOC(I)) * &
                   (1.0 - FSVF(I)) + QSGADD(I) - QSENSC(I) - QEVAPC(I) - &
                   QSTOR(I) - QMELTC(I)
        if (ABS(RESID(I)) < 5.0)                       ITER(I) = 0
        if (ABS(TSTEP(I)) < 1.0E-2)                   ITER(I) = 0
        if (NITER(I) == ITERMX .and. ITER(I) == 1)      ITER(I) = - 1
      end if
    end do ! loop 525
    !
    if (ITC < 2) then
      !
      !     * OPTION #1: SECANT/BISECTION ITERATION METHOD.
      !
      do I = IL1,IL2
        if (FI(I) > 0. .and. ITER(I) == 1) then
          if (NITER(I) == 1) then
            if (RESID(I) > 0.0) then
              TCAN(I) = TCAN(I) + TSTEP(I)
            else
              TCAN(I) = TCAN(I) - TSTEP(I)
            end if
          else
            if ((RESID(I) > 0. .and. TSTEP(I) < 0.) .or. &
                (RESID(I) < 0. .and. TSTEP(I) > 0.)) then
              TSTEP(I) = - TSTEP(I) / 2.0
            end if
            TCAN(I) = TCAN(I) + TSTEP(I)
          end if
          if (ABS(TCAN(I) - TFREZ) < 1.0E-6)             TCAN(I) = TFREZ
          NITER(I) = NITER(I) + 1
          NUMIT = NUMIT + 1
        end if
      end do ! loop 550
      !
    else
      !
      !     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.
      !
      do I = IL1,IL2
        if (FI(I) > 0. .and. ITER(I) == 1) then
          if (NITER(I) > 1) then
            DCFLUX = (CFLX(I) - CFLUXM(I)) / &
                     SIGN(MAX(.001,ABS(TSTEP(I))),TSTEP(I))
            if (ABS(TVIRTA(I) - TVIRTC(I)) < 0.4) &
                DCFLUX = MAX(DCFLUX,0.8 * DCFLXM(I))
            DXEVAP = (XEVAP(I) - XEVAPM(I)) / &
                     SIGN(MAX(.001,ABS(TSTEP(I))),TSTEP(I))
          else
            DCFLUX = 0.
            DXEVAP = 0.
          end if
          XEVAPM(I) = XEVAP(I)
          CFLUXM(I) = CFLX(I)
          DCFLXM(I) = DCFLUX
          DRDT0 = - 4.0 * SBC * TCAN(I) * TCAN(I) * TCAN(I) * &
                  (1.0 - FSVF(I)) * 2.0 - RHOAIR(I) * SPHAIR * (CFLX(I) + MAX(0., &
                  TCAN(I) - TPOTA(I)) * DCFLUX) + real(IEVAPC(I)) * CPHCHC(I) * &
                  RHOAIR(I) * (XEVAP(I) * WC(I) * A(I) * (B(I) - TFREZ) / &
                  ((TCAN(I) - B(I)) * (1.0 + WC(I))) ** 2 - (QCAN(I) - QA(I)) * &
                  DXEVAP) - CHCAP(I) / DELT
          TSTEP(I) = - RESID(I) / DRDT0
          TSTEP(I) = MAX( - 10.,MIN(5.,TSTEP(I)))
          TCAN(I) = TCAN(I) + TSTEP(I)
          if (ABS(TCAN(I) - TFREZ) < 1.0E-3)             TCAN(I) = TFREZ
          NITER(I) = NITER(I) + 1
          NUMIT = NUMIT + 1
          TAC(I) = TCAN(I)
          QAC(I) = QCAN(I)
          TVRTAC(I) = TVIRTC(I)
        end if
      end do ! loop 575
      !
    end if
    !
  end if
  if (NUMIT > 0)                                    GO TO 400
  !
  !     * IF CONVERGENCE HAS NOT BEEN REACHED FOR ITERATION METHOD #2,
  !     * CALCULATE TEMPERATURE AND FLUXES ASSUMING NEUTRAL STABILITY.
  !
  if (ITC == 2) then
    !
    NUMIT = 0
    !>
    !! After the iteration has been completed, if the Newton-Raphson method has been used, a check is carried
    !! out in loop 600 to ascertain whether convergence has not been reached (i.e. whether ITER = -1) for any
    !! location. In such cases it is assumed that conditions of near-neutral stability at the surface are the cause of
    !! the difficulty in finding a solution. The flags NUMIT and IEVAPC are set to zero, and a trial value of
    !! TCAN is calculated using the virtual potential temperature of the air and the canopy specific humidity. If
    !! the absolute value of RESID is \f$> 100 W m^{-2}\f$ , TCAN is set to this trial value. The values of \f$q_{0, c}\f$ and the
    !! components of the surface energy balance are recalculated as above, except that \f$Q_{H, c}\f$ and \f$Q_{E, c}\f$ are assumed
    !! as a first approximation to be zero. RESID is determined on this basis, and is then assigned to \f$Q_{H, c}\f$ and
    !! \f$Q_{E, c}\f$ . If RESID > 0, \f$Q_{E, c}\f$ is set to this value; otherwise it is equally divided between \f$Q_{H, c}\f$ and \f$Q_{E, c}\f$ . The
    !! residual is then reset to zero, and E(0) and \f$T_{v, c}\f$ are recalculated. NUMIT is incremented by 1, and the flag
    !! IEVAPC for the current location is set to 1.
    !!
    do I = IL1,IL2
      IEVAPC(I) = 0
      if (ITER(I) == - 1) then
        TCANT = TVIRTA(I) / (1.0 + 0.61 * QCAN(I))
        if (ABS(RESID(I)) > 100.) then
          TCAN(I) = TCANT
          WCAN = 0.622 * calcEsat(TCAN(I)) / PADRY(I)
          QCAN(I) = WCAN / (1.0 + WCAN)
          if (FSNOWC(I) > 0.0) then
            YEVAP = FRAINC(I) + FSNOWC(I)
          else
            YEVAP = FRAINC(I) + (1.0 - FRAINC(I)) * 10. / (10. + RC(I))
          end if
          QCAN(I) = YEVAP * QCAN(I) + (1.0 - YEVAP) * QA(I)
          QSTOR(I) = CHCAP(I) * (TCAN(I) - TCANO(I)) / DELT
          QLWOC(I) = SBC * TCAN(I) * TCAN(I) * TCAN(I) * TCAN(I)
          RESID(I) = QSWNC(I) + (QLWIN(I) + QLWOG(I) - 2.0 * QLWOC(I)) * &
                     (1.0 - FSVF(I)) + QSENSG(I) - QSTOR(I)
          if (RESID(I) > 0.) then
            QEVAPC(I) = RESID(I)
          else
            QEVAPC(I) = RESID(I) * 0.5
          end if
          QSENSC(I) = RESID(I) - QEVAPC(I)
          RESID(I) = 0.
          EVAPC(I) = QEVAPC(I) / CPHCHC(I)
          TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
          NUMIT = NUMIT + 1
          IEVAPC(I) = 1
        end if
      end if
    end do ! loop 600
    !>
    !! After loop 600, calls to DRCOEF or FLXSURFZ are performed to re-evaluate the surface turbulent
    !! transfer coefficients for any locations for which the fluxes were modified in the previous loop, i.e. for any
    !! locations for which IEVAPC was set to 1. After this a check is performed for unphysical values of the
    !! canopy temperature, i.e. for values greater than 100 C or less than -100 C. If such values are encountered,
    !! an error message is printed and a call to abort is carried out.
    !!
    !! Next a check is carried out to determine whether freezing or melting of intercepted water has occurred
    !! over the current time step. If this is the case, adjustments are required to the canopy temperature, the
    !! intercepted liquid and frozen water amounts and fractional coverages, and to \f$C_c\f$ and \f$\Delta Q_{S, c}\f$ . In loop 650, if
    !! there is liquid water stored on the canopy and the canopy temperature is less than 0 C, the first half of the
    !! adjustment to \f$\Delta Q_{S, c}\f$ is performed, and the flags ITER and NIT are set to 1. The available energy sink
    !! HFREZ is calculated from CHCAP and the difference between TCAN and 0 C, and compared with
    !! HCONV, calculated as the energy sink required to freeze all of the liquid water on the canopy. If
    !! HFREZ \f$\leq\f$ HCONV, the amount of water that can be frozen is calculated using the latent heat of melting.
    !! The fractional coverages of frozen and liquid water FSNOWC and FRAINC and their masses SNOCAN
    !! and RAICAN are adjusted accordingly, TCAN is set to 0 C, and the amount of energy involved is stored
    !! in the diagnostic variable QMELTC. Otherwise all of the intercepted liquid water is converted to frozen
    !! water, and the energy available for cooling the canopy is calculated as HCOOL = HFREZ - HCONV.
    !! This available energy is applied to decreasing the temperature of the canopy, using the specific heat of the
    !! canopy elements, and the amount of energy that was involved in the phase change is stored in the
    !! diagnostic variable QMELTC. In both cases \f$q_c\f$ and \f$T_{c, v}\f$ are recalculated, and at the end of the loop \f$C_c\f$ and
    !! \f$\Delta Q_{S, c}\f$ are re-evaluated.
    !!
    !
    if (NUMIT > 0) then
      if (ISLFD < 2) then
        call drcoef1(CDM, CDH, RIB, CFLUX, QA, QA, ZOSCLM, ZOSCLH, &
                    CRIB, TVIRTC, TVIRTA, VA, FI, IEVAPC, &
                    ILG, IL1, IL2)
      else
        call flxsurfz1(CDM, CDH, CFLUX, RIB, FTEMP, FVAP, ILMO, &
                      UE, FCOR, TPOTA, QA, ZRSLFM, ZRSLFH, VA, &
                      TCAN, QCAN, H, ZOM, ZOH, &
                      LZZ0, LZZ0T, FM, FH, ILG, IL1, IL2, FI, IEVAPC, JL)
      end if
    end if
    !
  end if
  !
  IBAD = 0
  !
  do I = IL1,IL2
    !         IF (FI(I)>0. .AND. ITER(I)==-1)                      THEN
    !             WRITE(6,6350) I,JL,NITER(I),RESID(I),TCAN(I),RIB(I)
    ! 6350         FORMAT('0CANOPY ITERATION LIMIT',3X,3I3,3(F8.2,E12.4))
    !         END IF
    if (FI(I) > 0. .and. (TCAN(I) < 173.16 .or. &
        TCAN(I) > 373.16)) then
      IBAD = I
    end if
  end do ! loop 625
  !
  if (IBAD /= 0) then
    write(6,6375) IBAD,JL,TCAN(IBAD),NITER(IBAD),ISNOW
6375 format('0BAD CANOPY ITERATION TEMPERATURE',3X,2I3,F16.2,2I4)
    write(6,6380) QSWNC(IBAD),QLWIN(IBAD),QLWOG(IBAD), &
                   QLWOC(IBAD),QSENSG(IBAD),QSENSC(IBAD), &
                   QEVAPC(IBAD),QSTOR(IBAD),QMELTC(IBAD)
    write(6,6380) TCAN(IBAD),TPOTA(IBAD),TZERO(IBAD)
6380 format(2X,9F10.2)
    call errorHandler('energBalVegSolve', - 2)
  end if
  !
  !     * POST-ITERATION CLEAN-UP.
  !
  NIT = 0
  do I = IL1,IL2
    if (FI(I) > 0.) then
      if (RAICAN(I) > 0. .and. TCAN(I) < TFREZ) then
        QSTOR(I) = - CHCAP(I) * TCANO(I) / DELT
        ITER(I) = 1
        NIT = NIT + 1
        HFREZ = CHCAP(I) * (TFREZ - TCAN(I))
        HCONV = RAICAN(I) * CLHMLT
        if (HFREZ <= HCONV) then
          RCONV = HFREZ / CLHMLT
          FSNOWC(I) = FSNOWC(I) + FRAINC(I) * RCONV / RAICAN(I)
          FRAINC(I) = FRAINC(I) - FRAINC(I) * RCONV / RAICAN(I)
          SNOCAN(I) = SNOCAN(I) + RCONV
          RAICAN(I) = RAICAN(I) - RCONV
          TCAN  (I) = TFREZ
          QMELTC(I) = - CLHMLT * RCONV / DELT
          WCAN = 0.622 * 611.0 / PADRY(I)
          QCAN(I) = WCAN / (1.0 + WCAN)
          TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
        else
          HCOOL = HFREZ - HCONV
          SNOCAN(I) = SNOCAN(I) + RAICAN(I)
          FSNOWC(I) = FSNOWC(I) + FRAINC(I)
          FRAINC(I) = 0.0
          TCAN  (I) = - HCOOL / (SPHVEG * CMASS(I) + SPHICE * &
                      SNOCAN(I)) + TFREZ
          QMELTC(I) = - CLHMLT * RAICAN(I) / DELT
          RAICAN(I) = 0.0
          WCAN = 0.622 * calcEsat(TCAN(I)) / PADRY(I)
          QCAN(I) = WCAN / (1.0 + WCAN)
          TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
        end if
        CHCAP(I) = SPHVEG * CMASS(I) + SPHICE * SNOCAN(I) + &
                   SPHW * RAICAN(I)
        QSTOR(I) = QSTOR(I) + CHCAP(I) * TCAN(I) / DELT
      else
        ITER(I) = 0
      end if
      if (ITC == 2) then
        TAC(I) = TCAN(I)
        QAC(I) = QCAN(I)
        TVRTAC(I) = TVIRTC(I)
      end if
    end if
  end do ! loop 650
  !
  !>
  !! In loop 675, if there is frozen water stored on the canopy and the canopy temperature is greater than 0 C,
  !! the first half of the adjustment to \f$\Delta Q_{S, c}\f$ is performed, and the flags ITER and NIT are set to 1. The
  !! available energy for melting, HMELT, is calculated from CHCAP and the difference between TCAN and
  !! 0 C, and compared with HCONV, calculated as the energy required to melt all of the frozen water on the
  !! canopy. If HMELT \f$\leq\f$ HCONV, the amount of frozen water that can be melted is calculated using the
  !! latent heat of melting. The fractional coverages of frozen and liquid water FSNOWC and FRAINC and
  !! their masses SNOCAN and RAICAN are adjusted accordingly, TCAN is set to 0 C, and the amount of
  !! energy involved is stored in the diagnostic variable QMELTC. Otherwise, all of the intercepted frozen
  !! water is converted to liquid water, and the energy available for warming the canopy is calculated as
  !! HWARM = HMELT - HCONV. This available energy is applied to increasing the temperature of the
  !! canopy, using the specific heats of the canopy elements, and the amount of energy that was involved in
  !! the phase change is stored in the diagnostic variable QMELTC. In both cases \f$q_c\f$ and \f$T_{c, v}\f$ are recalculated,
  !! and at the end of the loop \f$C_c\f$ and \f$\Delta Q_{S, c}\f$ are re-evaluated.
  !!
  do I = IL1,IL2
    if (FI(I) > 0.) then
      if (SNOCAN(I) > 0. .and. TCAN(I) > TFREZ) then
        QSTOR(I) = - CHCAP(I) * TCANO(I) / DELT
        ITER(I) = 1
        NIT = NIT + 1
        HMELT = CHCAP(I) * (TCAN(I) - TFREZ)
        HCONV = SNOCAN(I) * CLHMLT
        if (HMELT <= HCONV) then
          SCONV = HMELT / CLHMLT
          FRAINC(I) = FRAINC(I) + FSNOWC(I) * SCONV / SNOCAN(I)
          FSNOWC(I) = FSNOWC(I) - FSNOWC(I) * SCONV / SNOCAN(I)
          SNOCAN(I) = SNOCAN(I) - SCONV
          RAICAN(I) = RAICAN(I) + SCONV
          TCAN  (I) = TFREZ
          QMELTC(I) = CLHMLT * SCONV / DELT
          WCAN = 0.622 * 611.0 / PADRY(I)
          QCAN(I) = WCAN / (1.0 + WCAN)
          TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
        else
          HWARM = HMELT - HCONV
          RAICAN(I) = RAICAN(I) + SNOCAN(I)
          FRAINC(I) = FRAINC(I) + FSNOWC(I)
          FSNOWC(I) = 0.0
          TCAN  (I) = HWARM / (SPHVEG * CMASS(I) + SPHW * &
                      RAICAN(I)) + TFREZ
          QMELTC(I) = CLHMLT * SNOCAN(I) / DELT
          SNOCAN(I) = 0.0
          WCAN = 0.622 * calcEsat(TCAN(I)) / PADRY(I)
          QCAN(I) = WCAN / (1.0 + WCAN)
          TVIRTC(I) = TCAN(I) * (1.0 + 0.61 * QCAN(I))
        end if
        CHCAP(I) = SPHVEG * CMASS(I) + SPHW * RAICAN(I) + &
                   SPHICE * SNOCAN(I)
        QSTOR(I) = QSTOR(I) + CHCAP(I) * TCAN(I) / DELT
      end if
      if (ITC == 2) then
        TAC(I) = TCAN(I)
        QAC(I) = QCAN(I)
        TVRTAC(I) = TVIRTC(I)
      end if
    end if
  end do ! loop 675
  !
  if (NIT > 0) then
    !
    !     * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT)
    !     * AND OTHER RELATED QUANTITIES BETWEEN CANOPY AIR SPACE AND
    !     * ATMOSPHERE.
    !
    if (ISLFD < 2) then
      call drcoef1(CDM, CDH, RIB, CFLUX, QAC, QA, ZOSCLM, ZOSCLH, &
                  CRIB, TVRTAC, TVIRTA, VA, FI, ITER, &
                  ILG, IL1, IL2)
    else
      call flxsurfz1(CDM, CDH, CFLUX, RIB, FTEMP, FVAP, ILMO, &
                    UE, FCOR, TPOTA, QA, ZRSLFM, ZRSLFH, VA, &
                    TAC, QAC, H, ZOM, ZOH, &
                    LZZ0, LZZ0T, FM, FH, ILG, IL1, IL2, FI, ITER, JL)
    end if
  end if
  !
  !     * REMAINING CALCULATIONS.
  !
  if (ITC == 1) then
    !
    !>
    !! For the locations over which melting or freezing of water on the canopy has occurred, the surface fluxes
    !! must now be recalculated to reflect the modified canopy temperature and humidity. If NIT > 0, first
    !! DRCOEF or FLXSURFZ is called to re-evaluate the surface turbulent transfer coefficients over all
    !! locations where ITER has been set to 1. Loops 700 and 750 repeat the calculations done in loops 475
    !! and 500, to obtain the surface transfer coefficients. Loop 800 repeats the calculations of the surface
    !! fluxes done in loop 525.
    !!
    do I = IL1,IL2
      if (FI(I) > 0. .and. ITER(I) == 1) then
        XEVAP(I) = RBINV(I)
        QAC(I) = (QCAN(I) * XEVAP(I) + QZERO(I) * RAGINV(I) + &
                 QA(I) * CFLUX(I)) / (XEVAP(I) + RAGINV(I) + CFLUX(I))
        if (QAC(I) < QCAN(I)) then
          if (FSNOWC(I) > 0.0) then
            XEVAP(I) = (FRAINC(I) + FSNOWC(I)) / RB(I)
          else
            XEVAP(I) = FRAINC(I) / RB(I) + (1.0 - FRAINC(I)) / &
                       (RB(I) + RC(I))
            QAC(I) = (QCAN(I) * XEVAP(I) + QZERO(I) * RAGINV(I) + &
                     QA(I) * CFLUX(I)) / (XEVAP(I) + RAGINV(I) + &
                     CFLUX(I))
          end if
        else
          if (FSNOWC(I) > 1.0E-5) then
            XEVAP(I) = FSNOWC(I) / RB(I)
          else
            XEVAP(I) = 1.0 / RB(I)
          end if
        end if
        TAC(I) = (TCAN(I) * RBINV(I) + TPOTG(I) * RAGINV(I) + &
                 TPOTA(I) * CFLUX(I)) / (RBINV(I) + RAGINV(I) + CFLUX(I))
        TVRTAC(I) = TAC(I) * (1.0 + 0.61 * QAC(I))
        CFSENS(I) = RBINV(I)
        CFEVAP(I) = XEVAP(I)
      end if
    end do ! loop 700
    !
  else
    !
    do I = IL1,IL2
      if (FI(I) > 0. .and. ITER(I) == 1) then
        CFLX(I) = RBINV(I) * CFLUX(I) / (RBINV(I) + CFLUX(I))
        CFLX(I) = CFLUX(I) + (CFLX(I) - CFLUX(I)) * &
                  MIN(1.0,QSWINV(I) * 0.04)
        RA(I) = 1.0 / CFLX(I)
        if (QA(I) < QCAN(I)) then
          if (FSNOWC(I) > 0.0) then
            XEVAP(I) = (FRAINC(I) + FSNOWC(I)) / RA(I)
          else
            XEVAP(I) = FRAINC(I) / RA(I) + (1.0 - FRAINC(I)) / &
                       (RA(I) + RC(I))
          end if
        else
          if (FSNOWC(I) > 1.0E-5) then
            XEVAP(I) = FSNOWC(I) / RA(I)
          else
            XEVAP(I) = 1.0 / RA(I)
          end if
        end if
        QCAN(I) = RA(I) * XEVAP(I) * QCAN(I) + (1.0 - RA(I) * &
                  XEVAP(I)) * QA(I)
        CFSENS(I) = CFLX(I)
        CFEVAP(I) = CFLX(I)
        TAC(I) = TPOTA(I)
        QAC(I) = QA(I)
      end if
    end do ! loop 750
    !
  end if
  !
  do I = IL1,IL2
    if (FI(I) > 0. .and. ITER(I) == 1) then
      if (SNOCAN(I) > 0.) then
        CPHCHC(I) = CLHVAP + CLHMLT
      else
        CPHCHC(I) = CLHVAP
      end if
      QLWOC(I) = SBC * TCAN(I) * TCAN(I) * TCAN(I) * TCAN(I)
      QSENSC(I) = RHOAIR(I) * SPHAIR * CFSENS(I) * (TCAN(I) - TAC(I))
      if (FRAINC(I) > 0. .or. FSNOWC(I) > 0. .or. &
          RC(I) <= 5000. .or. QAC(I) > QCAN(I)) then
        EVAPC(I) = RHOAIR(I) * CFEVAP(I) * (QCAN(I) - QAC(I))
      else
        EVAPC(I) = 0.0
      end if
      if (EVAPC(I) < 0. .and. TCAN(I) >= TADP(I)) EVAPC(I) = 0.0
      if (SNOCAN(I) > 0.) then
        EVPWET(I) = SNOCAN(I) / DELT
        if (EVAPC(I) > EVPWET(I)) EVAPC(I) = EVPWET(I)
      else
        EVPWET(I) = RAICAN(I) / DELT
        if (EVAPC(I) > EVPWET(I)) then
          WTRANSP = (EVAPC(I) - EVPWET(I)) * DELT
          EVPMAX(I) = EVPWET(I)
          WTRTOT(I) = 0.0
          do J = 1,IG
            WTEST = WTRANSP * FROOT(I,J)
            WROOT(I,J) = MIN(WTEST,WAVAIL(I,J))
            WTRTOT(I) = WTRTOT(I) + WROOT(I,J)
            EVPMAX(I) = EVPMAX(I) + WROOT(I,J) / DELT
          end do
          if (EVAPC(I) > EVPMAX(I)) EVAPC(I) = EVPMAX(I)
        end if
      end if
      QEVAPC(I) = CPHCHC(I) * EVAPC(I)
      RESID(I) = QSWNC(I) + (QLWIN(I) + QLWOG(I) - 2.0 * QLWOC(I)) * &
                 (1.0 - FSVF(I)) + QSGADD(I) - QSENSC(I) - QEVAPC(I) - &
                 QSTOR(I) - QMELTC(I)
    end if
  end do ! loop 800
  !
  !>
  !! At the end of the subroutine, some final diagnostic and clean-up calculations are performed. If the
  !! evaporation rate from the canopy, EVAPC, is very small, it is added to the residual of the canopy energy
  !! balance equation, RESID, and then reset to zero. The overall residual is added to the sensible heat flux
  !! from the canopy. The energy balance of the surface underlying the canopy is then re-evaluated to take
  !! into account the new value of the canopy longwave radiation. If the surface temperature is close to 0 C,
  !! the residual of the equation is assigned to QMELTG, representing the energy associated with melting or
  !! freezing of water at the surface; otherwise it is assigned to the ground heat flux. If the water vapour flux
  !! is towards the canopy, the water sublimated or condensed is assigned to interception storage: to
  !! SNOCAN if SNOCAN > 0 (for consistency with the definition of CPHCHC), and to RAICAN
  !! otherwise. The diagnostic water vapour flux variables QFCF and QFCL, and the diagnostic change of
  !! canopy heat storage HTCC, are updated accordingly, and the canopy heat capacity CHCAP is
  !! recalculated; EVAPC is added to the diagnostic variable EVAP and then reset to zero. The net
  !! shortwave radiation, outgoing longwave radiation, sensible and latent heat fluxes are calculated for the
  !! whole canopy-ground surface ensemble; the evaporation rates are converted to \f$m s^{-1}\f$ ; and the iteration
  !! counter ITERCT is updated for the level corresponding to the subarea type and the value of NITER. Finally,
  !! the fraction of water removed by transpiration from each soil layer is updated according to the fractions
  !! determined above on the basis of water availability.
  !!
  do I = IL1,IL2
    if (FI(I) > 0.) then
      if (ABS(EVAPC(I)) < 1.0E-8) then
        RESID(I) = RESID(I) + QEVAPC(I)
        EVAPC(I) = 0.0
        QEVAPC(I) = 0.0
      end if
      QSENSC(I) = QSENSC(I) + RESID(I)
      if (ABS(TZERO(I) - TFREZ) < 1.0E-3) then
        QMELTG(I) = QSWNG(I) + FSVF(I) * QLWIN(I) + (1.0 - FSVF(I)) * &
                    QLWOC(I) - QLWOG(I) - QSENSG(I) - QEVAPG(I) - GZERO(I)
      else
        GZERO(I) = QSWNG(I) + FSVF(I) * QLWIN(I) + (1.0 - FSVF(I)) * &
                   QLWOC(I) - QLWOG(I) - QSENSG(I) - QEVAPG(I)
      end if
      if (EVAPC(I) < 0.) then
        if (SNOCAN(I) > 0.) then
          SNOCAN(I) = SNOCAN(I) - EVAPC(I) * DELT
          QFCF(I) = QFCF(I) + FI(I) * EVAPC(I)
          HTCC(I) = HTCC(I) - FI(I) * TCAN(I) * SPHICE * EVAPC(I)
        else
          RAICAN(I) = RAICAN(I) - EVAPC(I) * DELT
          QFCL(I) = QFCL(I) + FI(I) * EVAPC(I)
          HTCC(I) = HTCC(I) - FI(I) * TCAN(I) * SPHW * EVAPC(I)
        end if
        EVAP(I) = EVAP(I) + FI(I) * EVAPC(I)
        EVAPC(I) = 0.0
        CHCAP(I) = SPHVEG * CMASS(I) + SPHICE * SNOCAN(I) + &
                   SPHW * RAICAN(I)
      end if
      QSWNET(I) = QSWNG(I) + QSWNC(I) + QTRANS(I)
      QLWOUT(I) = FSVF(I) * QLWOG(I) + (1.0 - FSVF(I)) * QLWOC(I)
      QSENS(I) = QSENSC(I) + QSENSG(I) - QSGADD(I)
      QEVAP(I) = QEVAPC(I) + QEVAPG(I)
      EVAPC(I) = EVAPC(I) / RHOW
      EVAPG(I) = EVAPG(I) / RHOW
      ITERCT(I,KF1(I),NITER(I)) = ITERCT(I,KF1(I),NITER(I)) + 1
    end if
  end do ! loop 850

  if (ctem_on) then
    !
    !       * STORE AERODYNAMIC CONDUCTANCE FOR USE IN NEXT TIME STEP
    !       * OVERWRITE OLDER NUMBERS ONLY WHEN FRACTION OF CANOPY
    !       * OR FRACTION OF CANOPY OVER SNOW (AKA FI) IS > 0.
    !
    do I = IL1,IL2
      if (FI(I) > 0.) then
        CFLUXV(I) = CFLUX(I)
      else
        CFLUXV(I) = CFLUXV_IN(I)
      end if
    end do ! loop 900
  end if
  !
  do J = 1,IG ! loop 950
    do I = IL1,IL2
      if (FI(I) > 0.) then
        if (WTRTOT(I) > 0.0) FROOT(I,J) = WROOT(I,J) / WTRTOT(I)
      end if
    end do
  end do ! loop 950
  !

  return
end subroutine energBalVegSolve
