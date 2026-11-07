!> \file
!> Sublimation calculations for the snow pack on the ground.
!
subroutine snowSublimation (RHOSNO, ZSNOW, HCPSNO, TSNOW, EVAP, QFN, QFG, HTCS, & ! Formerly SNOVAP
                            WLOST, TRUNOF, RUNOFF, TOVRFL, OVRFLW, &
                            FI, R, S, RHOSNI, WSNOW, ILG, IL1, IL2, JL)
  !
  !     * AUG 25/11 - D.VERSEGHY. CORRECT CALCULATION OF TRUNOF
  !     *                         AND TOVRFL.
  !     * FEB 22/07 - D.VERSEGHY. NEW ACCURACY LIMITS FOR R AND S.
  !     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 26/02 - D.VERSEGHY. CHANGE RHOSNI FROM CONSTANT TO
  !     *                         VARIABLE.
  !     * APR 11/01 - M.LAZARE.   CHECK FOR EXISTENCE OF SNOW BEFORE
  !     *                         PERFORMING CALCULATIONS.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * AUG 16/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         INCORPORATE DIAGNOSTIC ARRAY "WLOST".
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         ADDITIONAL DIAGNOSTIC CALCULATION -
  !     *                         UPDATE HTCS.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. SUBLIMATION FROM SNOWPACK.
  !
  use classicParams, only : DELT, TFREZ, HCPW, HCPICE, RHOW, &
                             RHOICE, CLHMLT, CLHVAP

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL
  integer :: I
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: RHOSNO(ILG) !< Density of snow pack \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(inout) :: ZSNOW (ILG) !< Depth of snow pack \f$[m] (z_g)\f$
  real, intent(inout) :: HCPSNO(ILG) !< Heat capacity of snow pack \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(inout) :: TSNOW (ILG) !< Temperature of the snow pack \f$[C] (T_s)\f$
  real, intent(inout) :: EVAP  (ILG) !< Sublimation rate from snow surface at start of subroutine \f$[m s^{-1}]\f$
  real, intent(inout) :: QFN   (ILG) !< Sublimation from snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: QFG   (ILG) !< Evaporation from ground \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: HTCS  (ILG) !< Internal energy change of snow pack due to
  !< conduction and/or change in mass \f$[W m^{-2}] (I_s)\f$
  real, intent(inout) :: WLOST (ILG) !< Residual amount of water that cannot be supplied by surface stores \f$[kg m^{-2}]\f$
  real, intent(inout) :: TRUNOF(ILG) !< Temperature of total runoff [K]
  real, intent(inout) :: RUNOFF(ILG) !< Total runoff \f$[m s^{-1}]\f$
  real, intent(inout) :: TOVRFL(ILG) !< Temperature of overland flow [K]
  real, intent(inout) :: OVRFLW(ILG) !< Overland flow from top of soil column \f$[m s^{-1}]\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI    (ILG) !< Fractional coverage of subarea in question on modelled area \f$[ ] (X_i)\f$
  real, intent(in) :: R     (ILG) !< Rainfall rate incident on snow pack \f$[m s^{-1}]\f$
  real, intent(in) :: S     (ILG) !< Snowfall rate incident on snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: RHOSNI(ILG) !< Density of fresh snow \f$[kg m^{-3}]\f$
  real, intent(inout) :: WSNOW (ILG) !< Liquid water content of snow pack \f$[kg m^{-2}] (w_s)\f$
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: ZADD, ZLOST, ZREM
  !
  !-----------------------------------------------------------------------
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0. .and. (S(I) < 1.0E-11 .or. R(I) < 1.0E-11) &
        .and. ZSNOW(I) > 0.) then
      HTCS(I) = HTCS(I) - FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
      if (EVAP(I) < 0.) then
        ZADD = - EVAP(I) * DELT * RHOW / RHOSNI(I)
        RHOSNO(I) = (ZSNOW(I) * RHOSNO(I) + ZADD * RHOSNI(I)) / &
                    (ZSNOW(I) + ZADD)
        ZSNOW (I) = ZSNOW(I) + ZADD
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                    (RHOW * ZSNOW(I))
        EVAP  (I) = 0.0
      else
        ZLOST = EVAP(I) * DELT * RHOW / RHOSNO(I)
        if (ZLOST <= ZSNOW(I)) then
          ZSNOW(I) = ZSNOW(I) - ZLOST
          EVAP (I) = 0.0
          HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                      (RHOW * ZSNOW(I))
        else
          ZREM = (ZLOST - ZSNOW(I)) * RHOSNO(I) / RHOW
          ZSNOW(I) = 0.0
          HCPSNO(I) = 0.0
          EVAP(I) = ZREM * (CLHMLT + CLHVAP) / (CLHVAP * DELT)
          WLOST(I) = WLOST(I) - ZREM * RHOW * CLHMLT / CLHVAP
          if (RUNOFF(I) > 0. .or. WSNOW(I) > 0.) &
              TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + (TSNOW(I) + TFREZ) * &
              WSNOW(I) / RHOW) / (RUNOFF(I) + WSNOW(I) / RHOW)
          RUNOFF(I) = RUNOFF(I) + WSNOW(I) / RHOW
          if (OVRFLW(I) > 0. .or. WSNOW(I) > 0.) &
              TOVRFL(I) = (TOVRFL(I) * OVRFLW(I) + (TSNOW(I) + TFREZ) * &
              FI(I) * WSNOW(I) / RHOW) / (OVRFLW(I) + FI(I) * &
              WSNOW(I) / RHOW)
          OVRFLW(I) = OVRFLW(I) + FI(I) * WSNOW(I) / RHOW
          TSNOW(I) = 0.0
          WSNOW(I) = 0.0
          QFN(I) = QFN(I) - FI(I) * ZREM * RHOW / DELT
          QFG(I) = QFG(I) + FI(I) * EVAP(I) * RHOW
        end if
      end if
      HTCS(I) = HTCS(I) + FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                ZSNOW(I) / DELT
    end if
  end do ! loop 100
  return
end subroutine snowSublimation
!> \file
!!
!! @author D. Verseghy, M. Lazare
!!
!! These calculations are done if a snowpack is present and there is
!! no rainfall or snowfall occurring. The change of internal energy
!! \f$I_s\f$ of the snow pack as a result of the change in its mass is
!! calculated as the difference in \f$I_s\f$ between the beginning and end
!! of the subroutine:
!!
!! \f$\Delta I_s = X_i \Delta [ C_s z_s T_s ] / \Delta t\f$
!!
!! where \f$C_s\f$ represents the volumetric heat capacity of the snow
!! pack, \f$T_s\f$ its temperature, \f$\Delta t\f$ the length of the time step,
!! and \f$X_i\f$ the fractional coverage of the subarea under consideration
!! relative to the modelled area.
!!
!! If the sublimation rate EVAP over the snow pack is negative
!! (downward), the deposited depth of snow ZADD is calculated from
!! EVAP by converting it from a liquid water flux to a fresh snow
!! depth using RHOSNI, the fresh snow density. The snowpack density
!! is updated as a weighted average of the original snow density
!! RHOSNO and RHOSNI. The new snow depth is calculated as the sum of
!! the old snow depth ZSNOW and ZADD. The new volumetric heat
!! capacity of the snow pack is obtained from the heat capacities of
!! ice and water \f$C_i\f$ and \f$C_w\f$, the snow, ice and water densities
!! \f$\rho_s\f$ \f$\rho_i\f$, and \f$\rho_w\f$, and the water content and depth of the
!! snow pack \f$w_s\f$ and \f$z_s\f$, as:
!!
!! \f$C_s = C_i [\rho_s / \rho_i ] + C_w w_s /[\rho_w z_s]\f$
!!
!! If the sublimation rate is positive, the depth of the snow pack
!! ZLOST that is sublimated over the time step is calculated from
!! EVAP using RHOSNO. If ZLOST \f$\leq\f$ ZSNOW, the snow depth is reduced
!! and HCPSNO is recalculated. Otherwise the deficit amount ZREM is
!! calculated from ZLOST – ZSNOW and converted to a depth of water.
!! This amount is further converted to an evaporation rate by
!! applying a correction factor of \f$(L_m+ L_v)/L_v\f$, where \f$L_m\f$ is the
!! latent heat of melting and \f$L_v\f$ is the latent heat of vaporization
!! (to account for the fact that the energy is now being used to
!! evaporate water instead of sublimate snow). This necessarily
!! leads to a small discrepancy between the overall vapour flux for
!! the subarea that was originally calculated in energyBudgetDriver, and the
!! actual change of water storage in the subarea, and therefore this
!! discrepancy is added to the housekeeping variable WLOST for use
!! in the water balance checks done later in checkWaterBudget. If there was
!! liquid water in the snow pack, WSNOW, it is assigned to overall
!! runoff RUNOFF, and to overland flow OVRFLW. The resulting
!! temperatures of the runoff and overland flow, TRUNOF and TOVRFL,
!! are recalculated as weighted averages using the original runoff
!! amounts and temperatures, and the original snow temperature TSNOW
!! for WSNOW. The snow depth, heat capacity, temperature and water
!! content are all set to zero. Finally, since ZREM now becomes soil
!! evaporation rather than snow sublimation, the diagnostic
!! variables QFN and QFG, representing the vapour flux from snow and
!! soil respectively, are adjusted to reflect this.
