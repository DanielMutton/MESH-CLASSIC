!> \file
!> Solves surface energy balance for non-vegetated
!! subareas.
!! @author D. Verseghy, M. Lazare, A. Wu, P. Bartlett, Y. Delage, J. Cole, R. Brown, J. Melton, Y. Wu
!
subroutine energBalNoVegSolve (ISNOW, FI, & ! Formerly TSOLVE
                               QSWNET, QLWOUT, QTRANS, QSENS, QEVAP, EVAP, &
                               TZERO, QZERO, GZERO, QMELT, CDH, CDM, RIB, CFLUX, &
                               FTEMP, FVAP, ILMO, UE, H, &
                               QLWIN, TPOTA, QA, VA, PADRY, RHOAIR, &
                               ALVISG, ALNIRG, CRIB, CPHCH, CEVAP, TVIRTA, &
                               ZOSCLH, ZOSCLM, ZRSLFH, ZRSLFM, ZOH, ZOM, FCOR, &
                               GCONST, GCOEFF, TSTART, PCPR, TRSNOWG, FSSB, ALSNO, &
                               THLIQ, THLMIN, DELZW, RHOSNO, ZSNOW, ZPOND, &
                               IWATER, IEVAP, ITERCT, ISAND, &
                               ISLFD, ITG, ILG, IG, IL1, IL2, JL, NBS, ISNOALB, &
                               TSTEP, TVIRTS, EVBETA, Q0SAT, RESID, &
                               DCFLXM, CFLUXM, WZERO, TRTOP, A, B, &
                               LZZ0, LZZ0T, FM, FH, ITER, NITER, JEVAP, KF, &
                               ipeatland, co2conc, pressg, coszs, Cmossmas, &
                               dmoss, anmoss, rmlmoss, iday, daylength, pdd)
  !
  !     * OCT 30/16 - J. MELTON. Finish implementation of peatland code by
  !                              Yuanqiao Wu.
  !     * OCT 26/16 - D.VERSEGHY. ADD ZPOND TO CALCULATION OF EVPMAX.
  !     * JUL 22/15 - D.VERSEGHY. LIMIT CALCULATED EVAPORATION RATE
  !     *                         ACCORDING TO WATER AVAILABILITY.
  !     * JAN 09/15 - D.VERSEGHY. FIX TO SUPPRESS EVAPORATION FROM ROCK.
  !     * JUN 27/14 - D.VERSEGHY. CHANGE ITERATION LIMIT BACK TO 50 FOR
  !     *                         BISECTION SCHEME.
  !     * NOV 16/13 - J.COLE/     FINAL VERSION FOR GCM17:
  !     *             M.LAZARE.   - FIX COMPUTATION OF QSWNI OVER SNOW FREE
  !     *                           BARE SOIL for ISNOW=0 and ISNOALB=1 (NEED
  !     *                           TO SUM OVER THE 3 NEAR-IR BANDS).
  !     * JUN 22/13 - J.COLE/     - ADD "ISNOALB" OPTION (4-BAND SOLAR).
  !     *             M.LAZARE.   - MODIFY ABORT CONDITION FOR TOO COLD
  !     *                           TEMPS FROM 173 TO 123, SO WON'T
  !     *                           BLOW UP OVER ANTARCTICA.
  !     * OCT 14/11 - D.VERSEGHY. FOR POST-ITERATION CLEANUP WITH N-R SCHEME,
  !     *                         REMOVE CONDITION INVOLVING LAST ITERATION
  !     *                         TEMPERATURE.
  !     * DEC 07/09 - D.VERSEGHY. RESTORE EVAPORATION WHEN PRECIPITATION
  !     *                         IS OCCURRING.
  !     * MAR 13/09 - D.VERSEGHY. REPLACE SURFCON COMMON BLOCK WITH CLASSD2
  !     *                         REVISED CALL TO FLXSURFZ.
  !     * JAN 06/09 - D.VERSEGHY/M.LAZARE. SPLIT IF CONDITIONS FRAMING
  !     *                         300 LOOP.
  !     * FEB 25/08 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS; REMOVE
  !     *                         "ILW" SWITCH; SUPPRESS WATER VAPOUR FLUX
  !     *                         IF PRECIPITATION IS OCCURRING.
  !     * MAY 17/06 - D.VERSEGHY. SUPPRESS EVAPORATION WHEN PONDED WATER
  !     *                         IS FREEZING; ADD IL1 AND IL2 TO CALL TO
  !     *                         FLXSURFZ; REMOVE JL FROM CALL TO DRCOEF.
  !     * APR 13/05 - R.BROWN. ADD WINDLESS TRANFER COEFFICIENT TO QSENS
  !     *                         CALCULATION FOR SNOW PACKS.
  !     * DEC 17/04 - Y.DELAGE/D.VERSEGHY. ADD SWITCH TO USE EITHER SECANT/
  !     *                         BISECTION OR NEWTON-RAPHSON ITERATION
  !     *                         SCHEME (WITH NUMBER OF ITERATIONS LIMITED
  !     *                         TO FIVE AND CORRECTION FOR REMAINING
  !     *                         RESIDUAL).
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 06/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS
  !     *                         FROM ROUNDOFF ERRORS.
  !     * NOV 07/02 - Y.DELAGE/D.VERSEGHY. NEW CALL TO FLXSURFZ.
  !     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
  !     * MAR 28/02 - D.VERSEGHY. STREAMLINED SUBROUTINE CALL.
  !     *                         BYPASS EVAPORATION EFFICIENCY PARAMETER
  !     *                         IN CASES OF CONDENSATION.
  !     * JAN 18/02 - P.BARTLETT/D.VERSEGHY. NEW "BETA" FORMULATION FOR
  !     *                         BARE SOIL EVAPORATION BASED ON LEE AND
  !     *                         PIELKE.
  !     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
  !     * OCT 06/00 - D.VERSEGHY. CONDITIONAL "IF" IN ITERATION SEQUENCE
  !     *                         TO AVOID DIVIDE BY ZERO.
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
  !     *                         NOW DETERMINED IN ROUTINE "DRCOEF"
  !     *                         "CFLUX" NOW WORK FIELD INSTEAD OF "CLIMIT".
  !     * OCT 04/94 - D.VERSEGHY. CHANGE "CALL ABORT" TO "CALL errorHandler" TO
  !     *                         ENABLE RUNNING ON PCS.
  !     * JAN 24/94 - M.LAZARE.   UNFORMATTED I/O COMMENTED OUT IN LOOP 200.
  !     * JUL 29/93 - D.VERSEGHY. CLASS - VERSION 2.2.
  !     *                         REMOVE RE-DEFINITION OF QMELT NEAR END
  !     *                         (SINCE DONE ELSEWHERE ALREADY) AND
  !     *                         REDEFINE QSWNET FOR DIAGNOSTIC PURPOSES
  !     *                         TO INCLUDE TRANSMISSION THROUGH
  !     *                         SNOWPACK.
  !     * OCT 15/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. ITERATIVE SURFACE TEMPERATURE
  !     *                         CALCULATIONS FOR SNOW/SOIL.
  !
  use classicParams, only : DELT, TFREZ, SBC, SPHAIR, RHOW, BETA
  use peatlandsMod,  only : mossPht
  use generalutilsGAT,  only : calcEsat

  implicit none

  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ISNOW !< Flag indicating presence or absence of snow
  integer, intent(in) :: ISLFD, ITG, ILG, IG, IL1, IL2, JL, NBS
  integer :: I, IB
  integer, intent(in) :: ISNOALB !< Switch to model snow albedo in two or more wavelength bands
  !
  integer :: NUMIT,NIT,IBAD,ITERMX
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(inout) :: QSWNET(ILG)  !< Net shortwave radiation at surface \f$[W m^{-2}]\f$
  real, intent(inout) :: QLWOUT(ILG)  !< Upwelling longwave radiation at surface \f$[W m^{-2}]\f$ (L)
  real, intent(inout) :: QTRANS(ILG)  !< Shortwave radiation transmitted into surface \f$[W m^{-2}]\f$
  real, intent(inout) :: QSENS (ILG)  !< Sensible heat flux from surface \f$[W m^{-2}] (Q_H)\f$
  real, intent(inout) :: QEVAP (ILG)  !< Latent heat flux from surface \f$[W m^{-2}] (Q_E)\f$
  real, intent(inout) :: EVAP  (ILG)  !< Evaporation rate at surface \f$[kg m^{-2} s^{-1}] (E(0))\f$
  real, intent(inout) :: TZERO (ILG)  !< Temperature at surface \f$[K] (T(0))\f$
  real, intent(inout) :: QZERO (ILG)  !< Specific humidity at surface \f$[kg kg^{-1}] (q(0))\f$
  real, intent(inout) :: GZERO (ILG)  !< Heat flux into surface \f$[W m^{-2}] (G(0))\f$
  real, intent(inout) :: QMELT (ILG)  !< Heat available for melting snow or freezing
  !< water at the surface \f$[W m^{-2}]\f$
  real, intent(in) :: CDH   (ILG)  !< Surface drag coefficient for heat \f$[ ] (C_{DH}) \f$
  real, intent(in) :: CDM   (ILG)  !< Surface drag coefficient for momentum [ ]
  real, intent(in) :: RIB   (ILG)  !< Bulk Richardson number at surface [ ]
  real, intent(inout) :: CFLUX (ILG)  !< Product of surface drag coefficient and wind speed \f$[m s^{-1}]\f$
  real, intent(in) :: FTEMP (ILG)  !< Product of surface-air temperature gradient,
  !< drag coefficient and wind speed \f$[K m s^{-1}]\f$
  real, intent(in) :: FVAP  (ILG)  !< Product of surface-air humidity gradient, drag
  !< coefficient and wind speed \f$[kg kg^{-1} m s^{-1}]\f$
  real, intent(in) :: ILMO  (ILG)  !< Inverse of Monin-Obukhov roughness length \f$(m-1]\f$
  real, intent(in) :: UE    (ILG)  !< Friction velocity of air \f$[m s^{-1}]\f$
  real, intent(in) :: H     (ILG)  !< Height of the atmospheric boundary layer [m]
  real, intent(in) :: ZPOND (ILG)  !< Depth of ponded water on surface [m]
  real, intent(inout) :: anmoss(ilg)  !<
  real, intent(inout) :: rmlmoss(ilg) !<
  real :: cevapmoss(ilg)!<
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on
  !< modelled area [ ]
  real, intent(in) :: QLWIN (ILG)  !< Downwelling longwave radiation at bottom of
  !< atmosphere \f$[W m^{-2}]\f$
  real, intent(in) :: TPOTA (ILG)  !< Potential temperature of air at reference
  !< height \f$[K] (T_{a,pot})\f$
  real, intent(in) :: QA    (ILG)  !< Specific humidity at reference height \f$[kg kg^{-1}] (q_a)\f$
  real, intent(in) :: VA    (ILG)  !< Wind speed at reference height \f$[m s^{-1}] (v_a)\f$
  real, intent(in) :: PADRY (ILG)  !< Partial pressure of dry air \f$[Pa] (p_{dry})\f$
  real, intent(in) :: RHOAIR(ILG)  !< Density of air \f$[kg m^{-3}] (\rho_a)\f$
  real, intent(in) :: ALVISG(ILG)  !< Visible albedo of ground surface [ ]
  real, intent(in) :: ALNIRG(ILG)  !< Near-IR albedo of ground surface [ ]
  real, intent(in) :: CRIB  (ILG)  !< Richardson number coefficient \f$[K^{-1}]\f$
  real, intent(in) :: CPHCH (ILG)  !< Latent heat of vaporization at surface \f$[J kg^{-1}]\f$
  real, intent(in) :: CEVAP (ILG)  !< Soil evaporation efficiency coefficient \f$[ ] (\beta)\f$
  real, intent(in) :: TVIRTA(ILG)  !< Virtual potential temperature of air at
  !< reference height [K]
  real, intent(in) :: ZOSCLH(ILG)  !< Ratio of roughness length for heat to reference
  !< height for temperature and humidity [ ]
  real, intent(in) :: ZOSCLM(ILG)  !< Ratio of roughness length for momentum to
  !< reference height for wind speed [ ]
  real, intent(in) :: ZRSLFH(ILG)  !< Difference between reference height for
  !< temperature and humidity and height at which
  !< extrapolated wind speed goes to zero [m]
  real, intent(in) :: ZRSLFM(ILG)  !< Difference between reference height for wind
  !< speed and height at which extrapolated wind
  !< speed goes to zero [m]
  real, intent(in) :: ZOH   (ILG)  !< Surface roughness length for heat [m]
  real, intent(in) :: ZOM   (ILG)  !< Surface roughness length for momentum [m]
  real, intent(in) :: GCONST(ILG)  !< Intercept used in equation relating surface
  !< heat flux to surface temperature \f$[W m^{-2}]\f$
  real, intent(in) :: GCOEFF(ILG)  !< Multiplier used in equation relating surface
  !< heat flux to surface temperature \f$[W m^{-2} K^{-1}]\f$
  real, intent(in) :: TSTART(ILG)  !< Starting point for surface temperature
  !< iteration [K]
  real, intent(in) :: FCOR  (ILG)  !< Coriolis parameter \f$[s^{-1}]\f$
  real, intent(in) :: PCPR  (ILG)  !< Surface precipitation rate \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: RHOSNO(ILG)  !< Density of snow  \f$[kg m^{-3}]\f$
  real, intent(in) :: ZSNOW(ILG)   !< Depth of snow pack  [m]  \f$(z_s)\f$

  real, intent(in) :: THLIQ(ILG,IG)    !< Volumetric liquid water content of soil layers \f$[m^{3} m^{-3}]\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content remaining after freezing or evaporation  \f$[m^{3} m^{-3}]\f$
  real, intent(in) :: DELZW(ILG,IG)    !< Permeable thickness of soil layer  [m]
  real, intent(in) :: TRSNOWG(ILG,NBS) !< Short-wave transmissivity of snow pack over bare ground  [  ]
  real, intent(in) :: ALSNO(ILG,NBS)   !< Albedo of snow in each modelled wavelength band  [  ]
  real, intent(in) :: FSSB(ILG,NBS)    !< Total solar radiation in each modelled wavelength band  \f$[W m^{-2}]\f$
  real, intent(inout) :: TRTOP(ILG,NBS)   !<

  !
  integer, intent(in) :: IWATER(ILG)  !< Flag indicating condition of surface
  !< (dry, water-covered or snow-covered)
  integer, intent(inout) :: IEVAP (ILG)  !< Flag indicating whether surface
  !< evaporation is occurring or not
  integer, intent(inout) :: ITERCT(ILG,6,50) !< Counter of number of iterations required to
  !< solve energy balance for four subareas
  integer, intent(in) :: ISAND(ILG,IG)    !< Sand content flag

  integer, intent(in)  :: ipeatland(ilg)
  integer, intent(in) :: iday
  integer :: ievapmoss(ilg)
  real, intent(in) :: co2conc(ilg)
  real, intent(in)  :: pressg(ilg)
  real, intent(in) :: coszs(ilg)
  real, intent(in) :: Cmossmas(ilg)
  real, intent(in) :: dmoss(ilg)
  real, intent(in) :: daylength(ilg)
  real, intent(inout) :: pdd(ilg)
  !
  !     * INTERNAL WORK ARRAYS.
  !
  real, intent(inout) :: TSTEP (ILG), TVIRTS(ILG), EVBETA(ILG), Q0SAT (ILG), &
                         RESID (ILG), DCFLXM(ILG), CFLUXM(ILG), &
                         A     (ILG), B     (ILG), &
                         LZZ0  (ILG), LZZ0T (ILG), FM    (ILG), FH    (ILG), &
                         WZERO (ILG)
  !
  integer, intent(inout) :: ITER(ILG), NITER(ILG), JEVAP (ILG), KF(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: QSWNV(ilg), QSWNI, DCFLUX, DRDT0, TZEROT, QEVAPT, BOWEN, EZERO, EVPMAX(ILG)
  !
  !>
  !! For the surface temperature iteration, two alternative schemes
  !! are offered: the bisection method (selected if the flag ITG = 1)
  !! and the Newton-Raphson method (selected if ITG = 2). In the first
  !! case, the maximum number of iterations ITERMX is set to 50, and
  !! in the second case it is set to 5. An optional windless transfer
  !! coefficient EZERO is available, which can be used, following the
  !! recommendations of Brown et al. (2006) \cite Brown2006-ec, to prevent the sensible
  !! heat flux over snow packs from becoming vanishingly small under
  !! highly stable conditions. Currently, EZERO is set to zero.
  !!
  !     * INITIALIZATION AND PRE-ITERATION SEQUENCE.
  !
  if (ITG < 2) then
    ITERMX = 50
  else
    ITERMX = 12      ! was 5 YW March 27, 2015
  end if
  !
  !      IF (ISNOW==0) THEN
  !          EZERO=0.0
  !      ELSE
  !          EZERO=2.0
  !      END IF
  EZERO = 0.0
  !

  do I = IL1,IL2
    QSWNET(I) = 0.0
    QTRANS(I) = 0.0
    if (ipeatland(i) > 0) then
      qswnv(i) = 0.0
    end if
  end do
  !
  !>
  !!
  !! In the next section, calculations of surface net and transmitted
  !! shortwave radiation are done depending on whether a snow pack is present
  !! (ISNOW=1) and whether the incoming shortwave radiation is being supplied
  !! in two bands (one visible and one near-IR, ISNOALB=0) or four bands
  !! (one visible and three near-IR, ISNOALB=1).  For each band the net
  !! shortwave radiation is calculated as the incoming radiation multiplied
  !! by one minus the appropriate albedo.  If there is no snow present the
  !! shortwave transmissivity at the surface is set to zero and the absorbed
  !! visible radiation is calculated from the first radiation band.  If
  !! ISNOALB=0 the absorbed near-IR radiation is calculated from the second
  !! band and if ISNOALB=1 it is calculated from the sum of the second to
  !! fourth bands.  The total net shortwave radiation QSWNET  is then evaluated.
  !! If a snow pack is present, then if ISNOALB=0 the calculations are done
  !! separately for the visible and near-IR bands; if ISNOALB=1, they are done
  !! separately over the four wavelength bands. The transmissivity is calculated
  !! as an average value for ISNOALB=0, and separately for the four wavelength
  !! ranges for ISNOALB=1. The net shortwave radiation at the surface is
  !! obtained from the sum of the net values in each wavelength band,
  !! corrected for the respective loss by transmission into the surface.

  !> In the 50 loop, the initial value of the surface temperature TZERO
  !! is set to TSTART, which contains the value of TZERO from the previous
  !! time step, and the first step in the iteration sequence, TSTEP, is
  !! set to 1.0 K.  The flag ITER is set to 1 for each element of the set
  !! of modelled areas, indicating that its surface temperature has not
  !! yet been found.  The iteration counter NITER is initialized to 1 for
  !! each element.  Initial values are assigned to several other variables.
  !! In particular, the maximum evaporation from the ground surface that
  !! can be sustained for the current time step is specified as the total
  !! snow mass if snow is present, otherwise as the water ponded on the
  !! surface plus the total available water in the first soil layer.
  !!
  if (ISNOW == 0) then ! Use usual snow-free bare soil formulation
    do I = IL1,IL2
      if (FI(I) > 0.) then
        TRTOP(I,1) = 0.
        QSWNV(i) = FSSB(I,1) * (1.0 - ALVISG(I)) ! YW QSWNV to QSWNV(i)
        if (ISNOALB == 0) then
          QSWNI = FSSB(I,2) * (1.0 - ALNIRG(I))
        else if (ISNOALB == 1) then
          QSWNI = 0.0
          do IB = 2,NBS
            QSWNI = QSWNI + FSSB(I,IB) * (1.0 - ALNIRG(I))
          end do ! IB
        end if
        QSWNET(I) = QSWNV(i) + QSWNI
        QTRANS(I) = QSWNET(I) * TRTOP(I,1)
        QSWNET(I) = QSWNET(I) - QTRANS(I)
      end if
    end do ! I

  else

    if (ISNOALB == 0) then ! Use the existing snow albedo and transmission
      do I = IL1,IL2
        if (FI(I) > 0.) then
          TRTOP(I,1) = TRSNOWG(I,1)
          QSWNV(i) = FSSB(I,1) * (1.0 - ALSNO(I,1))  ! YW QSWNV to QSWNV(i)
          QSWNI = FSSB(I,2) * (1.0 - ALSNO(I,2))
          QSWNET(I) = QSWNV(i) + QSWNI
          QTRANS(I) = QSWNET(I) * TRTOP(I,1)
          QSWNET(I) = QSWNET(I) - QTRANS(I)
        end if
      end do ! I
    else if (ISNOALB == 1) then ! Use the band-by-band snow albedo and transmission
      do I = IL1,IL2
        QTRANS(I) = 0.0
        QSWNET(I) = 0.0
      end do ! I
      do IB = 1,NBS
        do I = IL1,IL2
          if (FI(I) > 0.) then
            TRTOP(I,IB) = TRSNOWG(I,IB)
            QSWNV(i) = FSSB(I,IB) * (1.0 - ALSNO(I,IB))
            QSWNET(I) = QSWNET(I) + FSSB(I,IB) * (1.0 - ALSNO(I,IB))
            QTRANS(I) = QTRANS(I) + QSWNV(i) * TRTOP(I,IB)
          end if
        end do ! I
      end do ! IB
      do I = IL1,IL2
        if (FI(I) > 0.) then
          QSWNET(I) = QSWNET(I) - QTRANS(I)
        end if
      end do ! I
    end if ! ISNOALB
  end if ! ISNOW
  !
  do I = IL1,IL2
    if (FI(I) > 0.) then
      TZERO(I) = TSTART(I)
      TSTEP(I) = 1.0
      ITER(I) = 1
      NITER(I) = 1
      !
      QMELT(I) = 0.0
      RESID(I) = 999999.
      DCFLXM(I) = 0.0
      CFLUX(I) = 0.0
      if (ISNOW == 1) then
        KF(I) = 3
        EVPMAX(I) = RHOSNO(I) * ZSNOW(I) / DELT
      else
        KF(I) = 6
        EVPMAX(I) = RHOW * (ZPOND(I) + (THLIQ(I,1) - THLMIN(I,1)) * &
                    DELZW(I,1)) / DELT
      end if
    end if
  end do ! loop 50

  !    Do moss photosynthesis:
  !!      moss subroutine finds ground evaporation rate and photosynthesis--\
  call mossPht(il1, il2, iday, qswnv, thliq, co2conc, tstart, zsnow, & ! EC Jan 302017.
               pressg, Cmossmas, dmoss, anmoss, rmlmoss, &
               cevapmoss, ievapmoss, ipeatland, daylength, pdd)

  !>
  !! The 100 continuation line marks the beginning of the surface
  !! temperature iteration sequence. First the flags NIT (indicating
  !! that there are still locations at the beginning of the current
  !! iteration step for which the surface temperature has not yet been
  !! found) and NUMIT (indicating that there are still locations at
  !! the end of the current iteration step for which the surface
  !! temperature has not yet been found) are set to zero. Loop 150 is
  !! then performed over the set of modelled areas. If ITER=1, NIT is
  !! incremented by one, and the initial value of the surface transfer
  !! coefficient CFLUXM for this iteration pass is set to its value
  !! from the previous pass. The virtual temperature at the surface,
  !! \f$T(0)_v\f$, is obtained using the standard expression (see
  !! documentation for subroutine energyBudgetDriver):
  !!
  !! \f$T(0)_v = T(0) [1 + 0.61 q(0)]\f$
  !!
  !! where T(0) is the surface temperature and q(0) is the specific
  !! humidity at the surface. The surface humidity can be obtained
  !! from the saturated specific humidity \f$q(0)_{sat}\f$ by making use of the
  !! definition of the surface evaporation efficiency \f$\beta\f$:
  !!
  !! \f$\beta = [q(0) – q_a]/[q(0)_{sat} - q_a]\f$
  !!
  !! where \f$q_a\f$ is the specific humidity of the air at the reference
  !! height. This expression is inverted to obtain an expression for
  !! q(0). The saturated specific humidity \f$q(0)_{sat}\f$ is determined from
  !! the mixing ratio at saturation, \f$w(0)_{sat}\f$:
  !!
  !! \f$q(0)_{sat} = w(0)_{sat}/[1 + w(0)_{sat}]\f$
  !!
  !! The saturation mixing ratio is a function of the saturation
  !! vapour pressure \f$e(0)_{sat}\f$ at the surface:
  !!
  !! \f$w(0)_{sat} = 0.622 e(0)_{sat}/(p_{dry})\f$
  !!
  !! where \f$p_{dry}\f$ is the partial pressure of dry air. For the
  !! saturated vapour pressure, following Emanuel (1994) \cite Emanuel1994-dt
  !! \f$e_{sat}\f$ is from the temperature \f$T_a\f$ and the freezing
  !! point \f$T_f\f$:
  !!
  !! \f$e_{sat} = exp[53.67957 - 6743.769 / T - 4.8451 * ln(T)]       T \geq T_f\f$
  !!
  !! \f$e_{sat} = exp[23.33086 - 6111.72784 / T + 0.15215 * log(T)]    T < T_f\f$
  !!
  !! where \f$T_f\f$ is the freezing point. If there is a snow cover or
  !! ponded water present on the surface (IWATER > 0), the surface
  !! evaporation efficiency EVBETA is set to 1 and q(0) is set to
  !! \f$q(0)_{sat}\f$. Otherwise EVBETA is set to CEVAP, the value obtained in
  !! subroutine energyBudgetPrep on the basis of ambient conditions, and q(0) is
  !! calculated as above. If \f$q(0) > q_a\f$ and the evaporation flag IEVAP
  !! has been set to zero, EVBETA is reset to zero and q(0) is reset
  !! to \f$q_a\f$. Finally, \f$T(0)_v\f$ is determined using the equation above.
  !
  !     * ITERATION SECTION.
  !     * LOOP IS REPEATED UNTIL SOLUTIONS HAVE BEEN FOUND FOR ALL POINTS
  !     * ON THE CURRENT LATITUDE CIRCLE(S).
  !
  !ignoreLint(1)
100 continue ! anchor for a GO TO statement that needs to be removed
  !
  NUMIT = 0
  NIT = 0
  do I = IL1,IL2 ! loop 150
    if (FI(I) > 0. .and. ITER(I) == 1) then
      NIT = NIT + 1
      CFLUXM(I) = CFLUX(I)
      WZERO(I) = 0.622 * calcEsat(TZERO(I)) / PADRY(I)
      Q0SAT(I) = WZERO(I) / (1.0 + WZERO(I))
      if (IWATER(I) > 0) then
        EVBETA(I) = 1.0
        QZERO(I) = Q0SAT(I)
      else
        !    evaporation coefficient evbeta is controled by moss in peatland
        if (ipeatland(i) == 0) then
          EVBETA(I) = CEVAP(I)
        else
          evbeta(i) = cevapmoss(i)
          ievap(i) = ievapmoss(i)
        end if

        QZERO(I) = EVBETA(I) * Q0SAT(I) + (1.0 - EVBETA(I)) * QA(I)
        if (QZERO(I) > QA(I) .and. IEVAP(I) == 0) then
          EVBETA(I) = 0.0
          QZERO(I) = QA(I)
        end if
      end if
      TVIRTS(I) = TZERO(I) * (1.0 + 0.61 * QZERO(I))
    end if
  end do ! loop 150
  !
  !>
  !! If NIT > 0, a subroutine is called to evaluate the stability-
  !! corrected surface drag coefficients for heat and momentum. The
  !! subroutine selection is made on the basis of the flag ISLFD. If
  !! ISLFD=1, indicating that the calculations are to be consistent
  !! with CCCma conventions, subroutine DRCOEF is called; if ISLFD=2,
  !! indicating that the calculations are to be consistent with RPN
  !! conventions, subroutine FLXSURFZ is called.
  !!
  if (NIT > 0) then
    !
    !     * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND
    !     * OTHER RELATED QUANTITIES.
    !
    if (ISLFD < 2) then
      call drcoef1(CDM, CDH, RIB, CFLUX, QZERO, QA, ZOSCLM, ZOSCLH, &
                  CRIB, TVIRTS, TVIRTA, VA, FI, ITER, &
                  ILG, IL1, IL2)
    else
      call flxsurfz1(CDM, CDH, CFLUX, RIB, FTEMP, FVAP, ILMO, &
                    UE, FCOR, TPOTA, QA, ZRSLFM, ZRSLFH, VA, &
                    TZERO, QZERO, H, ZOM, ZOH, &
                    LZZ0, LZZ0T, FM, FH, ILG, IL1, IL2, FI, ITER, JL)
    end if
    !
    !     * REMAINING CALCULATIONS.
    !
    !>
    !! In loop 175, the terms of the surface energy balance are
    !! evaluated. The energy balance equation is written as:
    !!
    !! \f$K_* + L_* - Q_H – Q_E – G(0) = 0\f$
    !!
    !! where \f$K_*\f$ is the net shortwave radiation, \f$L_*\f$ is the net longwave
    !! radiation, \f$Q_H\f$ is the sensible heat flux, \f$Q_E\f$ is the latent
    !! heat flux, and G(0) is the conduction into the surface. \f$K_*\f$ was
    !! evaluated earlier in loop 50. \f$L_*\f$ is obtained as the difference
    !! between the downwelling radiation \f$L \downarrow\f$ and the upwelling
    !! radiation \f$L \uparrow\f$, which in turn is determined using the Stefan-
    !! Boltzmann equation:
    !!
    !! \f$L \uparrow = \sigma T(0)^4\f$
    !!
    !! where \f$\sigma\f$ is the Stefan-Boltzmann constant. (It is assumed that
    !! natural surfaces, because of their radiative complexity, act as
    !! effective black bodies, so that their emissivity can be taken
    !! to be 1.) The sensible heat flux is given by
    !!
    !! \f$Q_H = [\rho_a c_p C_{DH} v_a + \epsilon_0] [T(0) – T_{a, pot}]\f$
    !!
    !! where \f$\rho_a\f$ is the density of the air, \f$c_p\f$ is its specific
    !! heat, \f$C_{DH}\f$ is the surface drag coefficient for heat, and \f$v_a\f$ and
    !! \f$T_{a, pot}\f$ are the wind speed and potential temperature respectively
    !! at the reference height. (Note in the code that the variable
    !! CFLUX, evaluated in subroutine DRCOEF or FLXSURFZ, represents
    !! the product of \f$C_{DH}\f$ and \f$v_a\f$.) The windless transfer coefficient
    !! \f$\epsilon_0\f$, evaluated at the beginning of the subroutine, is
    !! used only under stable conditions, i.e. if \f$T(0) < T_{a, pot}\f$. The
    !! evaporation rate at the surface, E(0), is calculated as
    !!
    !! \f$E(0) = \rho_a C_{DH} v_a [Q(0) – q_a]\f$
    !!
    !! The evaporation rate is constrained to be less than or equal to the
    !! maximum rate evaluated earlier in the subroutine.  \f$Q_E\f$ is
    !! obtained by multiplying E(0) by the latent heat of
    !! vaporization at the surface. The ground heat flux G(0) is
    !! determined as a linear function of T(0) (see documentation for
    !! subroutines soilHeatFluxPrep and snowHeatCond). It can be seen that each of the
    !! terms of the surface energy balance is a function of a single
    !! unknown, T(0) or TZERO. The residual RESID of the energy balance is now
    !! evaluated on the basis of the current estimation for TZERO. If
    !! the absolute value of RESID is less than \f$5.0 W m^{-2}\f$, or if the
    !! absolute value of the iteration step TSTEP most recently used
    !! is less than 0.01 K, the surface temperature is deemed to have
    !! been found and ITER is set to 0. If the iteration counter NITER
    !! is equal to the maximum number and ITER is still 1, ITER is set
    !! to -1.
    !!
    !! In the following section, the iteration sequence is moved ahead a step. If ITG = 1, then if NIT > 0 and if
    !! ITER for the array element in question is 1, the calculations for the bisection method of solution are
    !! performed. If NITER = 1 (indicating that this is the first step in the iteration), then if RESID > 0
    !! (indicating that the current value for TZERO had undershot the correct value), TZERO is incremented
    !! by 1 K; otherwise it is decremented by 1 K. If this is not the first step in the iteration, then if RESID >0
    !! and TSTEP < 0 (indicating that TZERO has undershot the correct value and the last temperature
    !! increment had been a negative one) or if RESID < 0 and TSTEP > 0 (indicating that TZERO has
    !! overshot the correct value and the last temperature increment had been a positive one), TSTEP is divided
    !! in half and its sign changed. TSTEP is then added to TZERO. The iteration counter NITER and the
    !! flag NUMIT are each incremented by one. The next loop contains an optional set of print statements
    !! that are executed if ITER=-1, that is if a solution for the surface temperature has not been found within
    !! the prescribed maximum number of iteration steps. Finally, if NUMIT > 0, the iteration cycle is repeated
    !! from line 100 on.
    !!
    do I = IL1,IL2
      if (FI(I) > 0. .and. ITER(I) == 1) then
        QLWOUT(I) = SBC * TZERO(I) * TZERO(I) * TZERO(I) * TZERO(I)
        if (TZERO(I) < TPOTA(I)) then
          QSENS(I) = (RHOAIR(I) * SPHAIR * CFLUX(I) + EZERO) * (TZERO(I) - &
                     TPOTA(I))
        else
          QSENS(I) = RHOAIR(I) * SPHAIR * CFLUX(I) * (TZERO(I) - &
                     TPOTA(I))
        end if
        EVAP(I) = RHOAIR(I) * CFLUX(I) * (QZERO(I) - QA(I))
        if (EVAP(I) > EVPMAX(I)) EVAP(I) = EVPMAX(I)
        QEVAP(I) = CPHCH(I) * EVAP(I)
        GZERO(I) = GCOEFF(I) * TZERO(I) + GCONST(I)
        RESID(I) = QSWNET(I) + QLWIN(I) - QLWOUT(I) - QSENS(I) - QEVAP(I) - &
                   GZERO(I)
        if (ABS(RESID(I)) < 5.0)                       ITER(I) = 0
        if (ABS(TSTEP(I)) < 1.0E-2)                   ITER(I) = 0
        if (NITER(I) == ITERMX .and. ITER(I) == 1)      ITER(I) = - 1
      end if
    end do ! loop 175
  end if
  !
  !>
  !! If ITG = 2, then if NIT > 0 and if ITER for the array element in
  !! question is 1, the calculations for the Newton-Raphson method of
  !! iteration are performed. In this approach, the value \f$x_{n+1}\f$ used at
  !! each iteration step is obtained from the value \f$x_n\f$ at the previous
  !! step as follows:
  !!
  !! \f$x_{n+1} = x_n – f(x_n)/f'(x_n)\f$
  !!
  !! Identifying \f$x_n\f$ with TZERO and \f$f(x_n)\f$ with the surface energy
  !! balance equation, it can be seen that the second term on the
  !! right-hand side corresponds to TSTEP; the numerator is equal to
  !! RESID and the denominator to the first derivative of the energy
  !! balance equation evaluated at TZERO, which in turn is equal to
  !! the sum of the derivatives of the individual terms:
  !!
  !! \f$d( L \uparrow)/dT = -4 \sigma T(0)^3\f$
  !!
  !! \f$d(Q_H)/dT = \rho_a c_p {C_{DH} v_a + [T(0) – T_{a, pot}] d(C_{DH} v_a)/dT}\f$
  !!
  !! \f$d(Q_E)/dT = L_v \rho_a{C_{DH} v_a dq(0)/dT+ [q(0) – q_a] d(C_{DH} v_a)/dT}\f$
  !!
  !! and dG(0)/dT is equal to the coefficient multiplying TZERO in the
  !! equation for G(0). (\f$L_v\f$ is the latent heat of vaporization at the
  !! surface.) At the end of the calculations the iteration counter
  !! NITER and the flag NUMIT are each incremented by one, and upon
  !! exiting the loop, if NUMIT > 0, the iteration cycle is repeated
  !! from line 100 on.
  !!
  if (ITG < 2) then
    !
    !     * OPTION #1: BISECTION ITERATION METHOD.
    !
    if (NIT > 0) then
      do I = IL1,IL2
        if (FI(I) > 0. .and. ITER(I) == 1) then
          if (NITER(I) == 1) then
            if (RESID(I) > 0.0) then
              TZERO(I) = TZERO(I) + 1.0
            else
              TZERO(I) = TZERO(I) - 1.0
            end if
          else
            if ((RESID(I) > 0. .and. TSTEP(I) < 0.) .or. &
                (RESID(I) < 0. .and. TSTEP(I) > 0.)) then
              TSTEP(I) = - TSTEP(I) / 2.0
            end if
            TZERO(I) = TZERO(I) + TSTEP(I)
          end if
        end if
        !
        if (FI(I) > 0. .and. ITER(I) == 1) then
          NITER(I) = NITER(I) + 1
          NUMIT = NUMIT + 1
        end if
      end do ! loop 180
    end if
    !
    !     do I=IL1,IL2
    !         IF (FI(I)>0. .AND. ITER(I)==-1)                      THEN
    !             WRITE(6,6250) I,JL,RESID(I),TZERO(I),RIB(I)
    ! 6250         FORMAT('0GROUND ITERATION LIMIT',3X,2I3,3(F8.2,E12.4))
    !         END IF
    ! 185 CONTINUE
    !
    if (NUMIT > 0) GO TO 100
    !
  else
    !
    !     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.
    !
    if (NIT > 0) then
      do I = IL1,IL2
        if (FI(I) > 0. .and. ITER(I) == 1) then
          if (NITER(I) > 1) then
            DCFLUX = (CFLUX(I) - CFLUXM(I)) / &
                     SIGN(MAX(.001,ABS(TSTEP(I))),TSTEP(I))
            if (ABS(TVIRTA(I) - TVIRTS(I)) < 0.4) &
                DCFLUX = MAX(DCFLUX,0.8 * DCFLXM(I))
            DCFLXM(I) = DCFLUX
          else
            DCFLUX = 0.
          end if
          DRDT0 = - 4.0 * SBC * TZERO(I) ** 3 &
                  - RHOAIR(I) * SPHAIR * (CFLUX(I) + MAX(0.,TZERO(I) - TPOTA(I)) &
                  * DCFLUX) - GCOEFF(I) &
                  + CPHCH(I) * RHOAIR(I) * (CFLUX(I) * WZERO(I) * A(I) &
                  * EVBETA(I) * (B(I) - TFREZ) / ((TZERO(I) - B(I)) * &
                  (1.0 + WZERO(I))) ** 2 - (QZERO(I) - QA(I)) * DCFLUX)
          TSTEP(I) = - RESID(I) / DRDT0
          TSTEP(I) = MAX( - 10.,MIN(5.,TSTEP(I)))
          TZERO(I) = TZERO(I) + TSTEP(I)
          NITER(I) = NITER(I) + 1
          NUMIT = NUMIT + 1
        end if
      end do ! loop 190
    end if
    !
    if (NUMIT > 0) GO TO 100
    !>
    !! After the iteration has been completed, NUMIT is reset to zero
    !! and a check is carried out to ascertain whether convergence has
    !! not been reached (i.e. whether ITER = -1) for any location. In
    !! such cases it is assumed that conditions of near-neutral
    !! stability at the surface are the cause of the difficulty in
    !! finding a solution. A trial value of TZERO is calculated using
    !! the virtual potential temperature of the air. If
    !! RESID > \f$50 W m^{-2}\f$, TZERO is set to this trial value. The values of
    !! q(0) and the components of the surface energy balance are
    !! recalculated as above, except that \f$Q_H\f$ and \f$Q_E\f$ are assumed as
    !! a first approximation to be zero. RESID is determined on this
    !! basis. If RESID is positive, it is assigned to \f$Q_E\f$; otherwise
    !! RESID is divided equally between \f$Q_H\f$ and \f$Q_E\f$, except in the
    !! case of an absolutely dry surface, in which case \f$Q_E\f$ is set to
    !! zero and \f$Q_H\f$ to RESID. RESID is reset to zero, E(0) is obtained
    !! from \f$Q_E\f$, and \f$T(0)_v\f$ is recalculated. The flag JEVAP for the
    !! location is set to 1, and NUMIT is incremented by 1. If NUMIT > 0
    !! at the end of this loop, DRCOEF or FLXSURF are called again, and
    !! their calculations are performed for any location where JEVAP is
    !! 1, to ensure consistency with the new surface temperature and
    !! humidity.
    !!
    !
    !     * IF CONVERGENCE HAS NOT BEEN REACHED, CALCULATE TEMPERATURE AND
    !     * FLUXES ASSUMING NEUTRAL STABILITY.
    !
    do I = IL1,IL2
      NUMIT = 0
      JEVAP(I) = 0
      if (FI(I) > 0. .and. ITER(I) == - 1) then
        TZEROT = TVIRTA(I) / (1.0 + 0.61 * QZERO(I))
        if (ABS(RESID(I)) > 50.) then
          TZERO(I) = TZEROT
          WZERO(I) = 0.622 * calcEsat(TZERO(I)) / PADRY(I)
          Q0SAT(I) = WZERO(I) / (1.0 + WZERO(I))
          QZERO(I) = EVBETA(I) * Q0SAT(I) + (1.0 - EVBETA(I)) * QA(I)
          QLWOUT(I) = SBC * TZERO(I) * TZERO(I) * TZERO(I) * TZERO(I)
          GZERO(I) = GCOEFF(I) * TZERO(I) + GCONST(I)
          RESID(I) = QSWNET(I) + QLWIN(I) - QLWOUT(I) - GZERO(I)
          if (RESID(I) > 0.) then
            QEVAP(I) = RESID(I)
          else
            QEVAP(I) = RESID(I) * 0.5
          end if
          if (IEVAP(I) == 0) QEVAP(I) = 0.0
          QSENS(I) = RESID(I) - QEVAP(I)
          RESID(I) = 0.
          EVAP(I) = QEVAP(I) / CPHCH(I)
          TVIRTS(I) = TZERO(I) * (1.0 + 0.61 * QZERO(I))
          JEVAP(I) = 1
          NUMIT = NUMIT + 1
        end if
      end if
    end do ! loop 195
    !
    if (NUMIT > 0) then
      if (ISLFD < 2) then
        call drcoef1(CDM, CDH, RIB, CFLUX, QZERO, QA, ZOSCLM, ZOSCLH, &
                    CRIB, TVIRTS, TVIRTA, VA, FI, JEVAP, &
                    ILG, IL1, IL2)
      else
        call flxsurfz1(CDM, CDH, CFLUX, RIB, FTEMP, FVAP, ILMO, &
                      UE, FCOR, TPOTA, QA, ZRSLFM, ZRSLFH, VA, &
                      TZERO, QZERO, H, ZOM, ZOH, &
                      LZZ0, LZZ0T, FM, FH, ILG, IL1, IL2, FI, JEVAP, JL)
      end if
    end if
    !
  end if
  !>
  !! At this point a check is performed for unphysical values of the
  !! surface temperature, i.e. for values greater than 100 C or less
  !! than -250 C. If such values are encountered, an error message is
  !! printed and a call to abort is carried out.
  !!
  !
  !     * CHECK FOR BAD ITERATION TEMPERATURES.
  !
  IBAD = 0
  do I = IL1,IL2
    if (FI(I) > 0. .and. (TZERO(I) < 123.16 .or. &
        TZERO(I) > 373.16)) then
      IBAD = I
    end if
  end do ! loop 200
  !
  if (IBAD /= 0) then
    write(6,6275) IBAD,JL,TZERO(IBAD),NITER(IBAD),ISNOW
6275 format('0BAD ITERATION TEMPERATURE',3X,2I3,F16.2,2I4)
    write(6,6280) QSWNET(IBAD),QLWIN(IBAD),QSENS(IBAD), &
         QEVAP(IBAD),GZERO(IBAD),CFLUX(IBAD),RIB(IBAD)
6280 format(2X,7F12.4)
    call errorHandler('energBalNoVegSolve', - 1)
  end if
  !>
  !! Finally, a check is performed to ensure that TZERO is not less
  !! than 0 C if ponded water is present on the surface (IWATER = 1)
  !! or greater than 0 C if snow is present on the surface
  !! (IWATER = 2), or greater than zero if the surface is an ice sheet
  !! (ISAND = -4). If any of these cases is true, TZERO is reset to
  !! the freezing point, and q(0) and \f$T(0)_v\f$ are recalculated. DRCOEF
  !! or FLXSURFZ are called, and their calculations are performed for
  !! all locations meeting these criteria. The components of the
  !! surface energy balance are recalculated; the residual amount is
  !! assigned to the energy associated with phase change of water at
  !! the surface, QMELT, and RESID is set to zero. If QMELT < 0, i.e.
  !! if there is freezing taking place, the evaporation flux QEVAP is
  !! added to it and then reset to zero (since if ponded water is
  !! freezing, it is unavailable for evaporation).
  !!
  !
  !     * POST-ITERATION CLEAN-UP.
  !
  NIT = 0
  do I = IL1,IL2
    if (FI(I) > 0.) then
      if (((IWATER(I) == 1 .and. TZERO(I) < TFREZ) .or. &
          (IWATER(I) == 2 .and. TZERO(I) > TFREZ)) .or. &
          (ISAND(I,1) == - 4 .and. TZERO(I) > TFREZ)) then
        TZERO(I) = TFREZ
        WZERO(I) = 0.622 * 611.0 / PADRY(I)
        QZERO(I) = WZERO(I) / (1.0 + WZERO(I))
        TVIRTS(I) = TZERO(I) * (1.0 + 0.61 * QZERO(I))
        ITER(I) = 1
        NIT = NIT + 1
      else
        ITER(I) = 0
      end if
    end if
  end do ! loop 300
  !
  if (NIT > 0) then
    !
    !       * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND
    !       * OTHER RELATED QUANTITIES.
    !
    if (ISLFD < 2) then
      call drcoef1(CDM, CDH, RIB, CFLUX, QZERO, QA, ZOSCLM, ZOSCLH, &
                  CRIB, TVIRTS, TVIRTA, VA, FI, ITER, &
                  ILG, IL1, IL2)
    else
      call flxsurfz1(CDM, CDH, CFLUX, RIB, FTEMP, FVAP, ILMO, &
                    UE, FCOR, TPOTA, QA, ZRSLFM, ZRSLFH, VA, &
                    TZERO, QZERO, H, ZOM, ZOH, &
                    LZZ0, LZZ0T, FM, FH, ILG, IL1, IL2, FI, ITER, JL)
    end if
  end if
  !>
  !! In the last half of the loop, some final adjustments are made to
  !! a few variables. If the evaporation flux is vanishingly small, it
  !! is added to RESID and reset to zero. If an anomalous case has
  !! arisen in which QMELT < 0 over a snow-covered surface or
  !! QMELT > 0 over a snow-free surface, QMELT is added to the heat
  !! flux into the surface and then reset to zero. Any remaining
  !! residual flux is added to \f$Q_H\f$. The shortwave radiation transmitted
  !! into the surface is added back to the net shortwave radiation for
  !! diagnostic purposes. The surface vapour flux is converted into
  !! units of \f$m s^{-1}\f$. Lastly, the iteration counter ITERCT is updated
  !! for the level corresponding to the subarea type and the value of
  !! NITER.
  !!
  !
  !     * REMAINING CALCULATIONS.
  !
  do I = IL1,IL2
    if (FI(I) > 0. .and. ITER(I) == 1) then
      QLWOUT(I) = SBC * TZERO(I) * TZERO(I) * TZERO(I) * TZERO(I)
      if (TZERO(I) < TPOTA(I)) then
        QSENS(I) = (RHOAIR(I) * SPHAIR * CFLUX(I) + EZERO) * (TZERO(I) - &
                   TPOTA(I))
      else
        QSENS(I) = RHOAIR(I) * SPHAIR * CFLUX(I) * (TZERO(I) - &
                   TPOTA(I))
      end if
      EVAP(I) = RHOAIR(I) * CFLUX(I) * (QZERO(I) - QA(I))
      if (EVAP(I) > EVPMAX(I)) EVAP(I) = EVPMAX(I)
      QEVAP(I) = CPHCH(I) * EVAP(I)
      GZERO(I) = GCOEFF(I) * TZERO(I) + GCONST(I)
      QMELT(I) = QSWNET(I) + QLWIN(I) - QLWOUT(I) - QSENS(I) - QEVAP(I) - &
                 GZERO(I)
      RESID(I) = 0.0
      if (QMELT(I) < 0.0) then
        QMELT(I) = QMELT(I) + QEVAP(I)
        QEVAP(I) = 0.0
        EVAP(I) = 0.0
      end if
    end if
    !
    if (FI(I) > 0.) then
      if (ABS(EVAP(I)) < 1.0E-8) then
        RESID(I) = RESID(I) + QEVAP(I)
        EVAP(I) = 0.0
        QEVAP(I) = 0.0
      end if
      if ((ISNOW == 1 .and. QMELT(I) < 0.0) .or. &
          (ISNOW == 0 .and. QMELT(I) > 0.0)) then
        GZERO(I) = GZERO(I) + QMELT(I)
        QMELT(I) = 0.0
      end if
      !              QSENS(I)=QSENS(I)+0.5*RESID(I)
      !              GZERO(I)=GZERO(I)+0.5*RESID(I)
      QSENS(I) = QSENS(I) + RESID(I)
      QSWNET(I) = QSWNET(I) + QTRANS(I)
      EVAP(I) = EVAP(I) / RHOW
      ITERCT(I,KF(I),NITER(I)) = ITERCT(I,KF(I),NITER(I)) + 1
    end if
  end do ! loop 350
  !
  return
end subroutine energBalNoVegSolve
