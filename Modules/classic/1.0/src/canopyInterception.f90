!> \file
!! Calculates canopy interception of rainfall and snowfall,
!! and determines rainfall/snowfall rates at ground surface as a
!! result of throughfall and unloading.
!! @author D. Verseghy, M. Lazare, S. Fassnacht, E. Chan, P. Bartlett
!
subroutine canopyInterception (IWATER, R, TR, S, TS, RAICAN, SNOCAN, TCAN, CHCAP, & ! Formerly CANADD
                               HTCC, ROFC, ROVG, PCPN, PCPG, FI, FSVF, &
                               CWLCAP, CWFCAP, CMASS, RHOSNI, TSURX, RDRIP, SDRIP, &
                               ILG, IL1, IL2, JL)

  !     * NOV 22/06 - E.CHAN/D.VERSEGHY. UNCONDITIONALLY SET TR AND TS.
  !     * JAN 05/05 - P.BARTLETT. CORRECT/REFINE SNOW INTERCEPTION
  !     *                         CALCULATIONS.
  !     * SEP 13/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 29/02 - D.VERSEGHY/S.FASSNACHT. NEW SNOW INTERCEPTION
  !     *                                     ALGORITHM,
  !     * JUL 24/02 - D.VERSEGHY. MOVE DIAGNOSTIC CALCULATIONS FROM
  !     *                         waterBudgetDriver INTO THIS ROUTINE; CHANGE
  !     *                         RHOSNI FROM CONSTANT TO VARIABLE.
  !     * JUN 20/02 - D.VERSEGHY. ADDITIONAL DIAGNOSTIC CALCULATIONS.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CALCULATE CANOPY INTERCEPTION; ADD
  !     *                         THROUGHFALL AND CANOPY DRIP TO
  !     *                         PRECIPITATION REACHING GROUND.
  !     *                         ADJUST CANOPY TEMPERATURE AND HEAT
  !     *                         CAPACITY.
  !
  use classicParams,       only : DELT, TFREZ, SPHW, SPHICE, SPHVEG, RHOW

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IWATER, ILG, IL1, IL2, JL
  integer :: I
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: R(ILG)  !< Rainfall rate over subarea in question \f$[m s^{-1}]\f$
  real, intent(inout) :: TR(ILG)  !< Temperature of rainfall [C]
  real, intent(inout) :: S(ILG)  !< Snowfall rate over subarea in question \f$[m s^{-1}]\f$
  real, intent(inout) :: TS(ILG)  !< Temperature of snowfall [C]
  real, intent(inout) :: RAICAN(ILG)  !< Intercepted liquid water stored on the canopy
  !! \f$[kg m^{-2}] (W_{l,c})\f$
  real, intent(inout) :: SNOCAN(ILG)  !< Intercepted frozen water stored on the canopy
  !! \f$[kg m^{-2}] (W_{f,c})\f$
  real, intent(inout) :: TCAN(ILG)  !< Temperature of vegetation canopy \f$[K] (T_c)\f$
  real, intent(inout) :: CHCAP(ILG)  !< Heat capacity of vegetation canopy \f$[J m^{-2} K^{-1}] (C_c)\f$
  real, intent(inout) :: HTCC(ILG)  !< Internal energy change of canopy due to changes
  !! in temperature and/or mass \f$[W m^{-2}] (I_c)\f$
  real, intent(inout) :: ROFC(ILG)  !< Liquid/frozen water runoff from vegetation \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: ROVG(ILG)  !< Liquid/frozen water runoff from vegetation to
  !! ground surface \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCPN(ILG)  !< Precipitation incident on snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: PCPG(ILG)  !< Precipitation incident on ground \f$[kg m^{-2} s^{-1}]\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI(ILG)  !< Fractional coverage of subarea in question on
  !! modelled area \f$[ ] (X_i)\f$
  real, intent(in) :: FSVF(ILG)  !< Sky view factor of surface under vegetation canopy [ ]
  real, intent(in) :: CWLCAP(ILG)  !< Interception storage capacity of vegetation for
  !! liquid water \f$[kg m^{-2}]\f$
  real, intent(in) :: CWFCAP(ILG)  !< Interception storage capacity of vegetation for
  !! frozen water \f$[kg m^{-2}] (W_{f,max})\f$
  real, intent(in) :: CMASS(ILG)  !< Mass of vegetation canopy \f$[kg m^{-2}]\f$
  real, intent(in) :: RHOSNI(ILG)  !< Density of fresh snow \f$[kg m^{-3}]\f$
  real, intent(in) :: TSURX(ILG)  !< Ground or snow surface temperature of subarea [K]
  !
  !     * INTERNAL WORK ARRAYS.
  !
  real, intent(inout) :: RDRIP(ILG), SDRIP(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: RTHRU, RINT, STHRU, SINT, TRCAN, TSCAN, RWXCES, SLOAD, SWXCES, &
          SNUNLD, CHCAPI, TCANI

  !-----------------------------------------------------------------------
  !>
  !> The calculations in this subroutine are performed if the rainfall
  !> or snowfall rates over the modelled area are greater than zero,
  !> or if the intercepted liquid water RAICAN or frozen water SNOCAN
  !> is greater than zero (to allow for unloading). The throughfall of
  !> rainfall or snowfall incident on the canopy, RTHRU or STHRU, is
  !> calculated from FSVF, the canopy gap fraction or sky view factor.
  !> The remaining rainfall or snowfall is assigned to interception,
  !> as RINT and SINT. The resulting temperature of liquid water on
  !> the canopy, TRCAN, is calculated as a weighted average of RAICAN
  !> at the canopy temperature TCAN, and RINT at the rainfall
  !> temperature TR. The resulting temperature of frozen water on the
  !> canopy, TSCAN, is calculated as a weighted average of SNOCAN at
  !> the canopy temperature TCAN, and SINT at the snowfall temperature
  !> TS.
  !>
  do I = IL1,IL2 ! loop 100
    RDRIP(I) = 0.0
    SDRIP(I) = 0.0
    if (FI(I) > 0. .and. (R(I) > 0. .or. S(I) > 0. .or. &
        RAICAN(I) > 0. .or. SNOCAN(I) > 0.)) then
      RTHRU = R(I) * FSVF(I)
      RINT = (R(I) - RTHRU) * DELT * RHOW
      STHRU = S(I) * FSVF(I)
      SINT = (S(I) - STHRU) * DELT * RHOSNI(I)
      if ((RAICAN(I) + RINT) > 0.) then
        TRCAN = (RAICAN(I) * TCAN(I) + RINT * TR(I)) / (RAICAN(I) + RINT)
      else
        TRCAN = 0.0
      end if
      if ((SNOCAN(I) + SINT) > 0.) then
        TSCAN = (SNOCAN(I) * TCAN(I) + SINT * TS(I)) / (SNOCAN(I) + SINT)
      else
        TSCAN = 0.0
      end if
      !
      !>
      !> Calculations are now done to ascertain whether the total
      !> liquid water on the canopy exceeds the liquid water
      !> interception capacity CWLCAP. If such is the case, this
      !> excess is assigned to RDRIP, water dripping off the
      !> canopy. The rainfall rate reaching the surface under the
      !> canopy is calculated as RDRIP + RTHRU, and the
      !> temperature of this water flux is calculated as a
      !> weighted average of RDRIP at a temperature of TRCAN, and
      !> RTHRU at the original rainfall temperature TR. The
      !> remaining intercepted water becomes CWLCAP. Otherwise,
      !> the rainfall rate reaching the surface under the canopy
      !> is set to RTHRU, and the liquid water on the canopy
      !> RAICAN is augmented by RINT.
      !>
      RWXCES = RINT + RAICAN(I) - CWLCAP(I)
      if (RWXCES > 0.) then
        RDRIP(I) = RWXCES / (DELT * RHOW)
        if ((RDRIP(I) + RTHRU) > 0.) then
          TR(I) = (RDRIP(I) * TRCAN + RTHRU * TR(I)) / &
                  (RDRIP(I) + RTHRU)
        else
          TR(I) = 0.0
        end if
        R(I) = RDRIP(I) + RTHRU
        RAICAN(I) = CWLCAP(I)
      else
        R(I) = RTHRU
        RAICAN(I) = RAICAN(I) + RINT
      end if
      !
      !>
      !> Interception and unloading of snow on the canopy is calculated using a more complex :: method. The
      !> amount of snow intercepted during a snowfall event over a time step, \f$\Delta W_{f, i}\f$ , or SLOAD, is obtained from
      !> the initial intercepted snow amount \f$W_{f, c}\f$ and the interception capacity \f$W_{f, max}\f$ , following Hedstrom and
      !> Pomeroy (1998), as:
      !>
      !> \f$\Delta W_{f, i} = (W_{f, max} – W_{f, c} ) [1 – exp(-S_{int} /W_{f, max} )]\f$
      !>
      !> where \f$S_{int}\f$ is the amount of snow incident on the canopy
      !> during the time step. The amount of snow not stored by
      !> interception, SWXCES, is calculated as SINT – SLOAD.
      !> Between and during precipitation events, snow is unloaded
      !> from the canopy through wind gusts and snow
      !> densification. These effects of these processes are
      !> estimated using an empirical exponential relationship for
      !> the snow unloading rate \f$W_{f, u}\f$ or SNUNLD, again following
      !> Hedstrom and Pomeroy (1998):
      !>
      !> \f$W_{f, u} = {W_{f, c} + \Delta W_{f, i} } exp (-U \Delta t)\f$
      !>
      !> where U is a snow unloading coefficient, assigned a value
      !> of \f$0.1 d^{-1}\f$ or \f$1.157 * 10^{-6} s{-1}\f$. The sum of SWXCES and
      !> SNUNLD is assigned to SDRIP, the snow or frozen water
      !> falling off the canopy. The snowfall rate reaching the
      !> surface under the canopy is calculated as SDRIP + STHRU,
      !> and the temperature of this water flux is calculated as
      !> a weighted average of SDRIP at a temperature of TSCAN,
      !> and STHRU at the original snowfall temperature TS. The
      !> frozen water stored on the canopy is recalculated as
      !> SNOCAN + SINT – SWXCES - SNUNLD. Otherwise, the snowfall
      !> rate reaching the surface under the canopy is set to
      !> STHRU, and SNOCAN is augmented by SINT.
      !>
      SLOAD = (CWFCAP(I) - SNOCAN(I)) * (1.0 - EXP( - SINT / CWFCAP(I)))
      SWXCES = SINT - SLOAD
      SNUNLD = (SLOAD + SNOCAN(I)) * (1.0 - EXP( - 1.157E-6 * DELT))  
      if (SWXCES > 0. .or. SNUNLD > 0.) then
        SDRIP(I) = (MAX(SWXCES,0.0) + SNUNLD) / (DELT * RHOSNI(I))
        if ((SDRIP(I) + STHRU) > 0.) then
          TS(I) = (SDRIP(I) * TSCAN + STHRU * TS(I)) / &
                  (SDRIP(I) + STHRU)
        else
          TS(I) = 0.0
        end if
        S(I) = SDRIP(I) + STHRU
        SNOCAN(I) = SNOCAN(I) + SINT - MAX(SWXCES,0.0) - SNUNLD
      else
        S(I) = STHRU
        SNOCAN(I) = SNOCAN(I) + SINT
      end if
      !
      !>
      !> In the final section of the subroutine, the initial heat
      !> capacity and temperature of the canopy are saved in
      !> temporary variables. The new canopy heat capacity is
      !> calculated as a weighted average over the specific heats
      !> of the liquid and frozen water stored on the canopy and
      !> the canopy mass. The canopy temperature is calculated as
      !> a weighted average over the stored liquid and frozen
      !> water at the updated temperatures TRCAN and TSCAN, and
      !> the vegetation mass at the original temperature TCAN.
      !> Then the change in internal energy \f$I_c\f$ of the vegetation
      !> canopy as a result of the water movement above is
      !> calculated as the difference in \f$I_c\f$ before and after
      !> these processes:
      !> \f[
      !> \Delta I_c = X_i \Delta [C_c T_c ] / \Delta t
      !> \f]
      !> where \f$C_c\f$ represents the canopy heat capacity, \f$T_c\f$ the
      !> canopy temperature, \f$\Delta t\f$ the length of the time step, and
      !> \f$X_i\f$ the fractional coverage of the subarea under
      !> consideration relative to the modelled area.
      !>
      !> Finally, the rainfall and snowfall temperatures are
      !> converted to degrees C. (In the absence of precipitation,
      !> both are set equal to the surface temperature of the
      !> subarea to avoid floating point errors in later
      !> subroutines.) For subareas with a snow cover
      !> (IWATER = 2), the water running off the canopy and the
      !> precipitation incident on the snow pack are updated using
      !> RDRIP and SDRIP. For subareas without snow cover
      !> (IWATER = 1), the water running off the canopy is updated
      !> using RDRIP and SDRIP, the precipitation incident on the
      !> snow pack is augmented by SDRIP, and the precipitation
      !> incident on bare ground is augmented by RDRIP.
      !>
      CHCAPI  = CHCAP(I)
      TCANI   = TCAN(I)
      CHCAP(I) = RAICAN(I) * SPHW + SNOCAN(I) * SPHICE + CMASS(I) * SPHVEG
      TCAN (I) = (RAICAN(I) * SPHW * TRCAN + SNOCAN(I) * SPHICE * TSCAN + &
                 CMASS(I) * SPHVEG * TCAN(I)) / CHCAP(I)
      HTCC (I) = HTCC(I) + FI(I) * (CHCAP(I) * TCAN(I) - CHCAPI * TCANI) / &
                 DELT
      if (R(I) > 0.0) then
        TR(I) = TR(I) - TFREZ
      else
        TR(I) = MAX(TSURX(I) - TFREZ,0.0)
      end if
      if (S(I) > 0.0) then
        TS(I) = TS(I) - TFREZ
      else
        TS(I) = MIN(TSURX(I) - TFREZ,0.0)
      end if
      if (IWATER == 2) then
        ROFC(I) = ROFC(I) + FI(I) * (RDRIP(I) * RHOW + SDRIP(I) * &
                  RHOSNI(I))
        PCPN(I) = PCPN(I) + FI(I) * (RDRIP(I) * RHOW + SDRIP(I) * &
                  RHOSNI(I))
      end if
      if (IWATER == 1) then
        ROFC(I) = ROFC(I) + FI(I) * (RDRIP(I) * RHOW + SDRIP(I) * &
                  RHOSNI(I))
        ROVG(I) = ROVG(I) + FI(I) * (RDRIP(I) * RHOW + SDRIP(I) * &
                  RHOSNI(I))
        PCPN(I) = PCPN(I) + FI(I) * SDRIP(I) * RHOSNI(I)
        PCPG(I) = PCPG(I) + FI(I) * RDRIP(I) * RHOW
      end if
    end if
  end do ! loop 100
  !
  return
end subroutine canopyInterception
