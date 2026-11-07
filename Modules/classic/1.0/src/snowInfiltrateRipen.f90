!> \file
!! Addresses infiltration of rain and meltwater into snow
!! pack, and snow ripening.
!! @author D. Verseghy, M. Lazare
!
subroutine snowInfiltrateRipen (R, TR, ZSNOW, TSNOW, RHOSNO, HCPSNO, WSNOW, & ! Formerly SNINFL
                                HTCS, HMFN, PCPG, ROFN, FI, ILG, IL1, IL2, JL)
  !
  !     * DEC 23/09 - D.VERSEGHY. RESET WSNOW TO ZERO WHEN SNOW
  !     *                         PACK DISAPPEARS.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * DEC 16/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         NEW DIAGNOSTIC FIELD "ROFN" ADDED.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. RAIN INFILTRATION INTO SNOWPACK.
  !
  use classicParams, only : DELT, TFREZ, HCPW, HCPICE, RHOW, &
                            RHOICE, CLHMLT, WSNCAP

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL
  integer :: I
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: R     (ILG)  !< Rainfall rate incident on snow pack \f$[m s^{-1}]\f$
  real, intent(inout) :: TR    (ILG)  !< Temperature of rainfall [C]
  real, intent(inout) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_g)\f$
  real, intent(inout) :: TSNOW (ILG)  !< Temperature of the snow pack \f$[C] (T_s)\f$
  real, intent(inout) :: RHOSNO(ILG)  !< Density of snow pack \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(inout) :: HCPSNO(ILG)  !< Heat capacity of snow pack \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(inout) :: WSNOW (ILG)  !< Liquid water content of snow pack \f$[kg m^{-2}] (w_s)\f$
  real, intent(inout) :: HTCS  (ILG)  !< Internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2}] (I_s)\f$
  real, intent(inout) :: HMFN  (ILG)  !< Energy associated with freezing or thawing of water in the snow pack \f$[W m^{-2}]\f$
  real, intent(inout) :: PCPG  (ILG)  !< Precipitation incident on ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: ROFN  (ILG)  !< Runoff reaching the ground surface from the bottom of the snow pack \f$[kg m^{-2} s^{-1}]\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on modelled area \f$[ ] (X_i)\f$
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: RAIN, HRCOOL, HRFREZ, HSNWRM, HSNMLT, ZMELT, ZFREZ, WAVAIL
  !
  !-----------------------------------------------------------------------
  !>
  !! The rainfall rate, i.e. the liquid water precipitation rate
  !! incident on the snow pack from the atmosphere,
  !! from canopy drip and/or from melting of the top of the snow pack,
  !! may cause warming and/or melting
  !! of the snow pack as a whole. The overall change of internal
  !! energy of the snow pack as a result of the
  !! rainfall added to it, \f$I_s\f$ or HTCS, is calculated as the difference
  !! in Is between the beginning and end of the
  !! subroutine:
  !!
  !! \f$\Delta I_s = X_i \Delta [C_s z_s T_s]/Delta t\f$
  !!
  !! where \f$C_s\f$ represents the volumetric heat capacity of the snow
  !! pack, \f$T_s\f$ its temperature, \f$\Delta\f$ the length of the
  !! time step, and \f$X_i\f$ the fractional coverage of the subarea under
  !! consideration relative to the modelled area.
  !!
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0. .and. R(I) > 0. .and. ZSNOW(I) > 0.) &
        then
      HTCS(I) = HTCS(I) - FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
      RAIN = R(I) * DELT
      !>
      !! Four diagnostic variables are evaluated at the outset.
      !! HRCOOL, the energy sink required to cool the
      !! whole rainfall amount to 0 C, is calculated from the
      !! rainfall rate and temperature, and the heat capacity of
      !! water. HRFREZ, the energy sink required to freeze all the
      !! rainfall, is obtained from the latent heat of
      !! melting, the density of water and the rainfall rate.
      !! HSNWRM, the energy required to warm the whole
      !! snow pack to 0 C, is calculated from the temperature,
      !! heat capacity and depth of the snow pack.
      !! HSNMLT, the energy required to melt all the snow, is
      !! obtained from the latent heat of melting and the
      !! heat capacity and depth of the snow pack.
      !!
      HRCOOL = TR(I) * HCPW * RAIN
      HRFREZ = CLHMLT * RHOW * RAIN
      HSNWRM = (0.0 - TSNOW(I)) * HCPSNO(I) * ZSNOW(I)
      HSNMLT = CLHMLT * RHOSNO(I) * ZSNOW(I)
      !>
      !! If HRCOOL \f$\geq\f$ (HSNWRM + HSNMLT), the energy contributed by
      !! the temperature of the rainfall is
      !! sufficient to warm to 0 C and melt the whole snow pack.
      !! HRCOOL is recalculated as the difference
      !! between HRCOOL and (HSNWRM + HSNMLT), and the snow depth
      !! is converted to a water depth
      !! ZMELT. The energy used to melt the snow is added to the
      !! diagnostic variables HMFN, representing the
      !! energy associated with water phase changes in the snow
      !! pack, and HTCS. The new temperature of the
      !! rainfall is calculated by applying HRCOOL over the new
      !! rainfall rate reaching the soil, which now
      !! includes the original rainfall rate, the melted snow pack
      !! and the liquid water that was contained in the
      !! snow pack. The depth, temperature, density, heat capacity
      !! and liquid water content of the snow pack are
      !! set to zero.
      !!
      if (HRCOOL >= (HSNWRM + HSNMLT)) then
        HRCOOL = HRCOOL - (HSNWRM + HSNMLT)
        ZMELT = ZSNOW(I) * RHOSNO(I) / RHOW
        HMFN(I) = HMFN(I) + FI(I) * CLHMLT * ZMELT * RHOW / DELT
        HTCS(I) = HTCS(I) + FI(I) * CLHMLT * ZMELT * RHOW / DELT
        TR(I) = HRCOOL / (HCPW * (ZMELT + RAIN + WSNOW(I) / RHOW))
        R(I) = R(I) + (ZMELT + WSNOW(I) / RHOW) / DELT
        ZSNOW(I) = 0.0
        TSNOW(I) = 0.0
        RHOSNO(I) = 0.0
        HCPSNO(I) = 0.0
        WSNOW(I) = 0.0
        !>
        !! If HRCOOL \f$\geq\f$ HSNWRM but HRCOOL < (HSNWRM + HSNMLT), the
        !! energy contributed by the
        !! temperature of the rainfall is sufficient to warm the
        !! whole snowpack to 0 C but not to melt all of it.
        !! HSNMLT is therefore recalculated as HRCOOL - HSNWRM, and
        !! used to determine a melted depth of
        !! the snowpack ZMELT, which is subtracted from the snow
        !! depth ZSNOW. The energy used to melt this
        !! depth of snow is added to HMFN and HTCS. The total water
        !! now available for retention in the
        !! snowpack, WAVAIL, is obtained as the sum of the mass of
        !! melted water and the mass of water originally
        !! retained in the snow pack, WSNOW. This amount is compared
        !! to the water retention capacity of the
        !! snow pack, calculated from the maximum retention
        !! percentage by weight, WSNCAP . If WAVAIL is greater than the water retention
        !! capacity, WSNOW is set to the capacity value and
        !! the excess is reassigned to ZMELT. Otherwise WSNOW is set
        !! to WAVAIL and ZMELT is set to zero.
        !! The temperature of the snow and the temperature TR of the
        !! rainfall reaching the ground surface are each
        !! set to 0 C, HCPSNO is recalculated, and ZMELT is added to
        !! the rainfall rate R.
        !!
      else if (HRCOOL >= HSNWRM .and. HRCOOL < (HSNWRM + HSNMLT)) &
                         then
        HSNMLT = HRCOOL - HSNWRM
        ZMELT = HSNMLT / (CLHMLT * RHOSNO(I))
        HMFN(I) = HMFN(I) + FI(I) * CLHMLT * ZMELT * RHOSNO(I) / DELT
        HTCS(I) = HTCS(I) + FI(I) * CLHMLT * ZMELT * RHOSNO(I) / DELT
        ZSNOW(I) = ZSNOW(I) - ZMELT
        WAVAIL = ZMELT * RHOSNO(I) + WSNOW(I)
        if (WAVAIL > (WSNCAP * ZSNOW(I) * RHOSNO(I))) then
          WSNOW(I) = WSNCAP * ZSNOW(I) * RHOSNO(I)
          ZMELT = (WAVAIL - WSNOW(I)) / RHOW
        else
          WSNOW(I) = WAVAIL
          ZMELT = 0.0
        end if
        TSNOW(I) = 0.0
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                    (RHOW * ZSNOW(I))
        TR(I) = 0.0
        R(I) = R(I) + ZMELT / DELT
        !>
        !! If HSNWRM \f$\geq\f$ (HRCOOL + HRFREZ), the energy sink of the
        !! snow pack is sufficient to cool to 0 C
        !! and freeze all of the rainfall. HSNWRM is recalculated
        !! as the difference between HSNWRM and
        !! (HRCOOL + HRFREZ). The energy used in the freezing,
        !! HRFREZ, is added to HMFN and HTCS.
        !! The rainfall is applied to increasing the density of the
        !! snow pack RHOSNO; if the new density is greater
        !! than the density of ice, RHOSNO is reset to the ice
        !! density and ZSNOW is recalculated. HCPSNO is
        !! also recalculated, and the new snow temperature is
        !! obtained from HSNWRM, HCPSNO and ZSNOW.
        !! R and TR are set to zero.
        !!
      else if (HSNWRM >= (HRCOOL + HRFREZ)) then
        HSNWRM = (HRCOOL + HRFREZ) - HSNWRM
        HMFN(I) = HMFN(I) - FI(I) * HRFREZ / DELT
        HTCS(I) = HTCS(I) - FI(I) * HRFREZ / DELT
        RHOSNO(I) = (RHOSNO(I) * ZSNOW(I) + RHOW * RAIN) / ZSNOW(I)
        if (RHOSNO(I) > RHOICE) then
          ZSNOW(I) = RHOSNO(I) * ZSNOW(I) / RHOICE
          RHOSNO(I) = RHOICE
        end if
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                    (RHOW * ZSNOW(I))
        TSNOW(I) = HSNWRM / (HCPSNO(I) * ZSNOW(I))
        TR(I) = 0.0
        R(I) = 0.0
        !>
        !! If HSNWRM > HRCOOL and HSNWRM < (HRCOOL + HRFREZ), the
        !! ! energy sink of the snow pack
        !! is sufficient to cool the rainfall to 0 C, but not to
        !! freeze all of it. HRFREZ is therefore recalculated as
        !! HSNWRM – HRCOOL, and used to determine a depth of rain to
        !! be frozen, ZFREZ. The energy used
        !! in the freezing is added to HMFN and HTCS. The frozen
        !! rainfall is applied to increasing the density of
        !! the snow pack as above; if the calculated density exceeds
        !! that of ice, RHOSNO and ZSNOW are
        !! recalculated. The water available for retention in the
        !! snow pack, WAVAIL, is obtained as the sum of the
        !! unfrozen rainfall and WSNOW, and compared to the water
        !! retention capacity of the snow pack. If
        !! WAVAIL is greater than the water retention capacity,
        !! WSNOW is set to the capacity value and WAVAIL
        !! is recalculated. Otherwise WSNOW is set to WAVAIL and
        !! WAVAIL is set to zero. The heat capacity of
        !! the snow is recalculated, R is calculated from WAVAIL,
        !! and TR and TSNOW are set to zero.
        !!
      else if (HSNWRM >= HRCOOL .and. HSNWRM < (HRCOOL + HRFREZ)) &
                         then
        HRFREZ = HSNWRM - HRCOOL
        ZFREZ = HRFREZ / (CLHMLT * RHOW)
        HMFN(I) = HMFN(I) - FI(I) * CLHMLT * ZFREZ * RHOW / DELT
        HTCS(I) = HTCS(I) - FI(I) * CLHMLT * ZFREZ * RHOW / DELT
        RHOSNO(I) = (RHOSNO(I) * ZSNOW(I) + RHOW * ZFREZ) / ZSNOW(I)
        if (RHOSNO(I) > RHOICE) then
          ZSNOW(I) = RHOSNO(I) * ZSNOW(I) / RHOICE
          RHOSNO(I) = RHOICE
        end if
        WAVAIL = (RAIN - ZFREZ) * RHOW + WSNOW(I)
        if (WAVAIL > (WSNCAP * ZSNOW(I) * RHOSNO(I))) then
          WSNOW(I) = WSNCAP * ZSNOW(I) * RHOSNO(I)
          WAVAIL = WAVAIL - WSNOW(I)
        else
          WSNOW(I) = WAVAIL
          WAVAIL = 0.0
        end if
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                    (RHOW * ZSNOW(I))
        R(I) = WAVAIL / (RHOW * DELT)
        TR(I) = 0.0
        TSNOW(I) = 0.0
      end if
      !>
      !! Finally, the calculation of the change in internal energy
      !! is completed, and the rainfall rate leaving the
      !! bottom of the snow pack and reaching the soil is added to
      !! the diagnostic variables PCPG and ROFN.
      !!
      HTCS(I) = HTCS(I) + FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
      PCPG(I) = PCPG(I) + FI(I) * R(I) * RHOW
      ROFN(I) = ROFN(I) + FI(I) * R(I) * RHOW
    end if
  end do ! loop 100
  !
  return
end subroutine snowInfiltrateRipen
