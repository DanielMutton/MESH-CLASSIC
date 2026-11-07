!> \file
!! Assesses water flux elements at the ground surface under
!! the vegetation canopy.
!! @author D. Verseghy, M. Lazare
!
subroutine waterUnderCanopy (IWATER, R, TR, S, TS, RHOSNI, EVAPG, QFN, QFG, & ! Formerly SUBCAN
                             PCPN, PCPG, FI, ILG, IL1, IL2, JL)
  !
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 21/04 - D.VERSEGHY. NEW LOWER LIMITS ON RADD AND SADD,
  !     *                         CONSISTENT WITH waterCalcPrep.
  !     * SEP 26/02 - D.VERSEGHY. BUGFIX IN CALCULATIONS OF QFN/QFG.
  !     * JUL 24/02 - D.VERSEGHY. MODIFICAITONS TO ALLOW FOR
  !     *                         SIMULTANEOUS RAINFALL AND SNOWFALL
  !     *                         CHANGE RHOSNI FROM CONSTANT TO
  !     *                         VARIABLE.
  !     * JUN 20/02 - D.VERSEGHY. UPDATE SUBROUTINE CALL.
  !     * NOV 09/00 - D.VERSEGHY. MOVE DIAGNOSTIC CALCULATIONS INTO
  !     *                         waterCalcPrep.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
  !     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         INCORPORATE DIAGNOSTICS.
  !     * APR 21/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. PERFORM "waterCalcPrep" CALCULATIONS UNDER
  !     *                         CANOPY: LUMP DOWNWARD WATER VAPOUR
  !     *                         FLUXES TOGETHER WITH PRECIPITATION
  !     *                         REACHING GROUND.
  !
  use classicParams, only : RHOW

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IWATER, ILG, IL1, IL2, JL
  integer             :: I
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: R     (ILG)  !< Rainfall rate incident on ground \f$[m s^{-1}]\f$
  real, intent(inout) :: S     (ILG)  !< Snowfall rate incident on ground \f$[m s^{-1}]\f$
  real, intent(out)   :: TR    (ILG)  !< Temperature of rainfall [C]
  real, intent(out)   :: TS    (ILG)  !< Temperature of snowfall [C]
  real, intent(inout) :: EVAPG (ILG)  !< Evaporation rate from surface \f$[m s^{-1}]\f$
  real, intent(inout) :: QFN   (ILG)  !< Sublimation from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFG   (ILG)  !< Evaporation from ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCPN  (ILG)  !< Precipitation incident on ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCPG  (ILG)  !< Precipitation incident on ground \f$[kg m^{-2} s^{-1}]\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: RHOSNI(ILG)  !< Density of fresh snow \f$[kg m^{-3}]\f$
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on modelled area [ ]
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: SADD, RADD
  !
  !-----------------------------------------------------------------------
  !>
  !! This subroutine starts with the precipitation rate under the
  !! canopy (a result of throughfall and unloading)
  !! and calculates the resulting overall evaporation or deposition
  !! rates.
  !!
  do I = IL1,IL2 ! loop 100
    !>
    !! For IWATER = 2 (snow on the ground under the canopy), the
    !! water vapour flux EVAPG at the ground
    !! surface is in the first instance assumed to be sublimation.
    !! Thus the first step is to compare it to the
    !! snowfall rate. The sum of the snowfall rate and the
    !! evaporation rate, SADD, is calculated as S – EVAPG,
    !! with EVAPG converted from a liquid water flux (the standard
    !! output from energBalVegSolve) to a snow flux. If
    !! SADD is greater than zero (indicating a downward flux) the
    !! snowfall rate is set to SADD and EVAPG is
    !! set to zero. Otherwise EVAPG is set to -SADD (converted back
    !! to a liquid water flux), and S and TS are
    !! set to zero.
    !!
    if (FI(I) > 0. .and. IWATER == 2) then
      if (S(I) > 0. .or. EVAPG(I) < 0.) then
        SADD = S(I) - EVAPG(I) * RHOW / RHOSNI(I)
        if (ABS(SADD) < 1.0E-12) SADD = 0.0
        if (SADD > 0.) then
          S(I) = SADD
          EVAPG(I) = 0.0
        else
          EVAPG(I) = - SADD * RHOSNI(I) / RHOW
          S(I) = 0.0
          TS(I) = 0.0
        end if
      else
        S(I) = 0.0
        TS(I) = 0.0
      end if
      !
      !>
      !! After this section, any remaining evaporative flux is
      !! compared to the rainfall rate. The sum RADD is
      !! calculated as R – EVAPG. If RADD is greater than zero,
      !! the rainfall rate is set to RADD and EVAPG is
      !! set to zero. Otherwise EVAPG is set to –RADD, and R and
      !! TR are set to zero.
      !!
      if (R(I) > 0. .or. EVAPG(I) < 0.) then
        RADD = R(I) - EVAPG(I)
        if (ABS(RADD) < 1.0E-12) RADD = 0.0
        if (RADD > 0.) then
          R(I) = RADD
          EVAPG(I) = 0.0
        else
          EVAPG(I) = - RADD
          R(I) = 0.0
          TR(I) = 0.0
        end if
      else
        R(I) = 0.0
        TR(I) = 0.0
      end if
    end if
  end do ! loop 100
  !
  !>
  !! Analogous calculations are done for IWATER = 1 (bare ground under
  !! the canopy). In this case EVAPG
  !! is assumed in the first instance to be evaporation or
  !! condensation. Thus the first step is to compare it to
  !! the rainfall rate, and the same steps are followed as in the
  !! paragraph above. Afterwards, any remaining
  !! vapour flux is compared to the snowfall rate. If SADD is positive
  !! (downward), EVAPG, which is now
  !! considered to be absorbed into the snowfall rate, is subtracted
  !! from the ground vapour flux QFG and
  !! added to the snow vapour flux QFN. If SADD is negative (upward),
  !! S, which has now been absorbed
  !! into the evaporative flux, is subtracted from the snow
  !! precipitation flux PCPN and added to the ground
  !! precipitation flux PCPG.
  !!
  do I = IL1,IL2 ! loop 200
    if (FI(I) > 0. .and. IWATER == 1) then
      if (R(I) > 0. .or. EVAPG(I) < 0.) then
        RADD = R(I) - EVAPG(I)
        if (ABS(RADD) < 1.0E-12) RADD = 0.0
        if (RADD > 0.) then
          R(I) = RADD
          EVAPG(I) = 0.0
        else
          EVAPG(I) = - RADD
          R(I) = 0.0
          TR(I) = 0.0
        end if
      else
        R(I) = 0.0
        TR(I) = 0.0
      end if
      !
      if (S(I) > 0. .or. EVAPG(I) < 0.) then
        SADD = S(I) - EVAPG(I) * RHOW / RHOSNI(I)
        if (ABS(SADD) < 1.0E-12) SADD = 0.0
        if (SADD > 0.) then
          S(I) = SADD
          QFN(I) = QFN(I) + FI(I) * EVAPG(I) * RHOW
          QFG(I) = QFG(I) - FI(I) * EVAPG(I) * RHOW
          EVAPG(I) = 0.0
        else
          EVAPG(I) = - SADD * RHOSNI(I) / RHOW
          PCPN(I) = PCPN(I) - FI(I) * S(I) * RHOSNI(I)
          PCPG(I) = PCPG(I) + FI(I) * S(I) * RHOSNI(I)
          S(I) = 0.0
          TS(I) = 0.0
        end if
      else
        S(I) = 0.0
        TS(I) = 0.0
      end if
    end if
  end do ! loop 200
  !
  return
end subroutine waterUnderCanopy
