!> \file
!! Adds snow incident on the ground surface to the snow pack.
!
subroutine snowAddNew (ALBSNO, TSNOW, RHOSNO, ZSNOW, HCPSNO, HTCS, & ! Formerly SNOADD
                       FI, S, TS, RHOSNI, WSNOW, ILG, IL1, IL2, JL)
  !
  !     * NOV 17/11 - M.LAZARE.   CHANGE SNOW ALBEDO REFRESHMENT
  !     *                         THRESHOLD (SNOWFALL IN CURRENT
  !     *                         TIMESTEP) FROM 0.005 TO 1.E-4 M.
  !     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
  !     * SEP 10/04 - R.HARVEY/D.VERSEGHY. INCREASE SNOW ALBEDO
  !     *                         REFRESHMENT THRESHOLD; ADD
  !     *                         "IMPLICIT NONE" COMMAND.
  !     * JUL 26/02 - D.VERSEGHY. CHANGE RHOSNI FROM CONSTANT TO
  !     *                         VARIABLE.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. ACCUMULATION OF SNOW ON GROUND.
  !
  use classicParams, only : DELT, TFREZ, HCPW, HCPICE, RHOW, RHOICE

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL
  integer :: I
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(out)   :: ALBSNO(ILG)  !< Albedo of snow [ ]
  real, intent(inout) :: TSNOW (ILG)  !< Temperature of the snow pack \f$[C] (T_s)\f$
  real, intent(inout) :: RHOSNO(ILG)  !< Density of snow pack \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(inout) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_s)\f$
  real, intent(inout) :: HCPSNO(ILG)  !< Heat capacity of snow pack \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(inout) :: HTCS  (ILG)  !< Internal energy change of snow pack due to
  !< conduction and/or change in mass \f$[W m^{-2}] (I_s)\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on modelled area \f$[ ] (X_i)\f$
  real, intent(in) :: S     (ILG)  !< Snowfall rate incident on snow pack \f$[m s^{-1}]\f$
  real, intent(in) :: TS    (ILG)  !< Temperature of snowfall [C]
  real, intent(in) :: RHOSNI(ILG)  !< Density of fresh snow \f$[kg m^{-3}]\f$
  real, intent(in) :: WSNOW (ILG)  !< Liquid water content of snow pack \f$[kg m^{-2}] (w_s)\f$

  !     * TEMPORARY VARIABLES.
  !
  real :: SNOFAL, HCPSNP
  !

  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0. .and. S(I) > 0.) then
      HTCS  (I) = HTCS(I) - FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                  ZSNOW(I) / DELT
      SNOFAL = S(I) * DELT
      if (SNOFAL >= 1.E-4) then
        ALBSNO(I) = 0.84
      else if (.not.(ZSNOW(I) > 0.)) then
        ALBSNO(I) = 0.50
      end if
      HCPSNP = HCPICE * RHOSNI(I) / RHOICE
      TSNOW (I) = ((TSNOW(I) + TFREZ) * ZSNOW(I) * HCPSNO(I) + &
                  (TS   (I) + TFREZ) * SNOFAL  * HCPSNP) / &
                  (ZSNOW(I) * HCPSNO(I) + SNOFAL * HCPSNP) - &
                  TFREZ
      RHOSNO(I) = (ZSNOW(I) * RHOSNO(I) + SNOFAL * RHOSNI(I)) / &
                  (ZSNOW(I) + SNOFAL)
      ZSNOW (I) = ZSNOW(I) + SNOFAL
      HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                  (RHOW * ZSNOW(I))
      HTCS  (I) = HTCS(I) + FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * &
                  ZSNOW(I) / DELT
    end if
  end do ! loop 100
  !
  return
end subroutine snowAddNew
!> \file
!!
!! @author D. Verseghy, M. Lazare, R. Harvey
!!
!! The change of internal energy HTCS of the snow pack as a result of
!! the snowfall added to it is calculated as the difference in Is
!! between the beginning and end of the subroutine:
!!
!! \f$\Delta I_s = X_i \Delta [C_s z_s T_s] / \Delta t\f$
!!
!! where \f$C_s\f$ represents the volumetric heat capacity of the snow
!! pack, \f$T_s\f$ its temperature, \f$\Delta\f$ the length of the time step,
!! and \f$X_i\f$ the fractional coverage of the subarea under consideration
!! relative to the modelled area.
!!
!! The amount of snow incident at the given time step, SNOFAL, is
!! calculated from S and the timestep length DELT. If
!! SNOFAL \f$\geq\f$ 0.1 mm, the snow albedo is set to the fresh snow value
!! of 0.84. Otherwise, if the snow is falling on bare ground, its
!! initial albedo is set to the old snow value of 0.50. The heat
!! capacity of the precipitating snow, HCPSNP, is calculated from
!! the fresh snow density RHOSNI and the heat capacity and density
!! of ice. The new temperature of the snow pack is calculated as a
!! weighted average of its old temperature, weighted by the snow
!! depth ZSNOW and heat capacity HCPSNO, and the snowfall
!! temperature, weighted by SNOFAL and HCPSNP. The new density of
!! snow is calculated as a weighted average of the original density
!! RHOSNO and RHOSNI, and the new snow depth is calculated as
!! ZSNOW + SNOFAL. Finally, the new heat capacity of the snow pack
!! is obtained from the heat capacities of ice and water \f$C_i\f$ and
!! \f$C_w\f$, the snow, ice and water densities \f$\rho_s\f$ \f$rho_i\f$, and \f$\rho_w\f$,
!! and the water content and depth of the snow pack \f$w_s\f$ and \f$z_s\f$,
!! as:
!!
!! \f$C_s = C_i [ \rho_s /\rho_i ] + C_w w_s /[\rho_w z_s]\f$
!!
