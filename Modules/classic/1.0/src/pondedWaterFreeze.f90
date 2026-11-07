!> \file
!> Addresses freezing of water ponded on ground surface.
!! @author D. Verseghy, M. Lazare
!
subroutine pondedWaterFreeze (ZPOND, TPOND, ZSNOW, TSNOW, ALBSNO, RHOSNO, HCPSNO, & ! Formerly TFREEZ
                              GZERO, HMFG, HTCS, HTC, WTRS, WTRG, FI, QFREZ, WSNOW, &
                              TA, TBAR, ISAND, IG, ILG, IL1, IL2, JL)
  !
  !     * JAN 06/09 - D.VERSEGHY. SET QFREZ TO ZERO AFTER CALCULATION
  !     *                         OF HADD.
  !     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUN 20/02 - D.VERSEGHY. COSMETIC CHANGES TO SUBROUTINE CALL
  !     *                         SHORTENED CLASS4 COMMON BLOCK.
  !     * MAY 24/02 - D.VERSEGHY. PASS IN ENTIRE SOIL TEMPERATURE
  !     *                         ARRAY.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS.
  !     * AUG 16/95 - D.VERSEGHY. TWO NEW ARRAYS TO COMPLETE WATER
  !     *                         BALANCE DIAGNOSTICS.
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         REVISE CALCULATION OF HTC.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. FREEZING OF PONDED WATER.
  !
  use classicParams, only : DELT, TFREZ, HCPW, HCPICE, RHOW, &
                            RHOICE, CLHMLT

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IG, ILG, IL1, IL2, JL
  integer             :: I
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: ZPOND (ILG)  !< Depth of ponded water \f$[m] (z_p)\f$
  real, intent(inout) :: TPOND (ILG)  !< Temperature of ponded water \f$[C] (T_p)\f$
  real, intent(inout) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (z_g)\f$
  real, intent(inout) :: TSNOW (ILG)  !< Temperature of the snow pack \f$[C] (T_s)\f$
  real, intent(out)   :: ALBSNO(ILG)  !< Albedo of snow [ ]
  real, intent(inout) :: RHOSNO(ILG)  !< Density of snow pack \f$[kg m^{-3}] (\rho_s)\f$
  real, intent(inout) :: HCPSNO(ILG)  !< Heat capacity of snow pack \f$[J m^{-3} K^{-1}] (C_s)\f$
  real, intent(inout) :: GZERO (ILG)  !< Heat flow into soil surface \f$[W m^{-2}]\f$
  real, intent(inout) :: HTCS  (ILG)  !< Internal energy change of snow pack due to conduction and/or change in mass \f$[W m^{-2}]\f$ (Is)
  real, intent(inout) :: WTRS  (ILG)  !< Water transferred into or out of the snow pack \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: WTRG  (ILG)  !< Water transferred into or out of the soil \f$[kg m^{-2} s^{-1}]\f$
  !
  real, intent(inout) :: HMFG  (ILG,IG) !< Energy associated with phase change of water in soil layers \f$[W m^{-2}]\f$
  real, intent(inout) :: HTC   (ILG,IG) !< Internal energy change of soil layer due to conduction and/or change in mass \f$[W m^{-2}]\f$ (Ig)
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in)    :: FI    (ILG)    !< Fractional coverage of subarea in question on modelled area \f$[ ] (X_i)\f$
  real, intent(inout) :: QFREZ (ILG)    !< Energy sink for freezing of water at the ground surface \f$[W m^{-2}]\f$
  real, intent(in)    :: WSNOW (ILG)    !< Liquid water content of snow pack \f$[kg m^{-2}] (w_s)\f$
  real, intent(in)    :: TA    (ILG)    !< Air temperature [K]
  real, intent(in)    :: TBAR  (ILG,IG) !< Temperature of soil layer \f$[C] (T_g)\f$
  !
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: ZFREZ, HADD, HCOOL, HCONV, TTEST, TLIM, HEXCES
  !
  !-----------------------------------------------------------------------
  !>
  !! Freezing of water ponded on the ground surface occurs if an
  !! energy sink QFREZ is produced as a result
  !! of the solution of the surface energy balance, or if the pond
  !! temperature at the beginning of the
  !! subroutine has been projected to be below 0 C. The change of
  !! internal energy I in the snow and first soil
  !! layer (which for the purposes of diagnostic calculations includes
  !! the ponded water) as a result of these
  !! processes is calculated as the difference in I between the
  !! beginning and end of the subroutine:
  !!
  !! \f$\Delta I_s = X_i \Delta(C_s T_s z_s) / \Delta t\f$
  !!
  !! \f$\Delta I_g = X_i \Delta(C_w T_p z_p)/\Delta t\f$
  !!
  !! where the C terms represent volumetric heat capacities, the T
  !! terms temperatures, and the z terms depths
  !! of the snow pack and the ponded water respectively, \f$\Delta t\f$ is the
  !! length of the time step, and \f$X_i\f$ the
  !! fractional coverage of the subarea under consideration relative
  !! to the modelled area.
  !!
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0. .and. ZPOND(I) > 0. .and. (TPOND(I) < 0. &
        .or. QFREZ(I) < 0.)) then
      HTCS(I) = HTCS(I) - FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * ZSNOW(I) / &
                DELT
      ZFREZ = 0.0
      HADD = - QFREZ(I) * DELT
      QFREZ(I) = 0.0
      if (TPOND(I) < 0.) then
        HADD = HADD - TPOND(I) * HCPW * ZPOND(I)
        TPOND(I) = 0.0
      end if
      !>
      !! The energy sink HADD to be applied to the ponded water is
      !! calculated from QFREZ and the pond
      !! temperature TPOND (if it is below 0 C). Two diagnostic
      !! variables, HCOOL and HCONV, are calculated as the energy
      !! sink required to cool the ponded water to 0 C, and that
      !! required both to cool and to freeze
      !! the ponded water, respectively. If HADD \f$\leq\f$ HCOOL, the
      !! available energy sink is only sufficient to
      !! decrease the temperature of the ponded water. This
      !! decrease is applied, and the energy used is added to
      !! the internal energy HTC of the first soil layer.
      !!
      HCOOL = TPOND(I) * HCPW * ZPOND(I)
      HCONV = HCOOL + CLHMLT * RHOW * ZPOND(I)
      HTC (I,1) = HTC (I,1) - FI(I) * HCPW * (TPOND(I) + TFREZ) * &
                  ZPOND(I) / DELT
      if (HADD <= HCOOL) then
        TPOND(I) = TPOND(I) - HADD / (HCPW * ZPOND(I))
        HTC(I,1) = HTC(I,1) + FI(I) * HADD / DELT
        !>
        !! If HADD > HCOOL but HADD \f$\leq\f$ HCONV, the available energy
        !! sink is sufficient to decrease the
        !! ponded water temperature to 0 C and to freeze part of it.
        !! The energy used in freezing is calculated as
        !! HADD – HCOOL, and is used to calculate the depth of frozen
        !! water ZFREZ, which is then subtracted
        !! from the ponded water depth ZPOND. HCOOL is added to HTC,
        !! and ZFREZ is converted from a
        !! depth of water to a depth of ice and added to the snow
        !! pack. If there is not a pre-existing snow pack, the
        !! snow albedo is set to the background value for old snow,
        !! 0.50. The temperature, density and heat
        !! capacity of the snow pack are recalculated.
        !!
      else if (HADD <= HCONV) then
        HADD = HADD - HCOOL
        ZFREZ = HADD / (CLHMLT * RHOW)
        ZPOND(I) = ZPOND(I) - ZFREZ
        HTC(I,1) = HTC(I,1) + FI(I) * HCOOL / DELT
        ZFREZ = ZFREZ * RHOW / RHOICE
        if (.not.(ZSNOW(I) > 0.0)) ALBSNO(I) = 0.50
        TSNOW(I) = TSNOW(I) * HCPSNO(I) * ZSNOW(I) / (HCPSNO(I) * ZSNOW(I) &
                   + HCPICE * ZFREZ)
        RHOSNO(I) = (RHOSNO(I) * ZSNOW(I) + RHOICE * ZFREZ) / (ZSNOW(I) &
                    + ZFREZ)
        ZSNOW(I) = ZSNOW(I) + ZFREZ
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                    (RHOW * ZSNOW(I))
        TPOND(I) = 0.0
        !>
        !! If HADD > HCONV, the available energy sink is sufficient
        !! to cool and freeze the whole depth of
        !! ponded water, and also to decrease the temperature of the
        !! ice thus formed. The ponded water depth is
        !! converted to a depth of ice ZFREZ, and HCOOL is added to
        !! HTC. In order to avoid unphysical
        !! temperature decreases caused by the length of the time
        !! step and the small ponding depth, a limit is set on
        !! the allowable decrease. The theoretical temperature of
        !! the newly formed ice, TTEST, is calculated from
        !! HADD – HCONV. If there is a pre-existing snow pack, the
        !! limiting temperature TLIM is set to the
        !! minimum of the snow pack temperature TSNOW, the first soil
        !! layer temperature, and 0 C; otherwise it is
        !! set to the minimum of the air temperature, the first soil
        !! layer temperature and 0 C. If TTEST < TLIM,
        !! the new ice is assigned TLIM as its temperature in the
        !! recalculation of TSNOW, the excess heat sink
        !! HEXCES is calculated from HADD and TLIM and assigned to
        !! the ground heat flux GZERO, and
        !! HADD – HEXCES is used to update HTC. Otherwise the new ice
        !! is assigned TTEST as its temperature
        !! in the recalculation of TSNOW, and HADD is used to update
        !! HTC. If there is not a pre-existing snow
        !! pack, the snow albedo is set to the background value. The
        !! density, depth and heat capacity of the snow
        !! pack are recalculated, and the pond depth and temperature
        !! are set to zero.
        !!
      else
        HADD = HADD - HCONV
        ZFREZ = ZPOND(I) * RHOW / RHOICE
        HTC(I,1) = HTC(I,1) + FI(I) * HCOOL / DELT
        TTEST = - HADD / (HCPICE * ZFREZ)
        if (ZSNOW(I) > 0.0) then
          TLIM = MIN(TSNOW(I),TBAR(I,1))
        else
          TLIM = MIN(TA(I) - TFREZ,TBAR(I,1))
        end if
        TLIM = MIN(TLIM,0.0)
        if (TTEST < TLIM) then
          HEXCES = HADD + TLIM * HCPICE * ZFREZ
          GZERO(I) = GZERO(I) - HEXCES / DELT
          HTC(I,1) = HTC(I,1) + FI(I) * (HADD - HEXCES) / DELT
          TSNOW(I) = (TSNOW(I) * HCPSNO(I) * ZSNOW(I) + &
                     TLIM * HCPICE * ZFREZ) &
                     / (HCPSNO(I) * ZSNOW(I) + HCPICE * ZFREZ)
        else
          TSNOW(I) = (TSNOW(I) * HCPSNO(I) * ZSNOW(I) + TTEST * HCPICE * &
                     ZFREZ) / (HCPSNO(I) * ZSNOW(I) + HCPICE * ZFREZ)
          HTC(I,1) = HTC(I,1) + FI(I) * HADD / DELT
        end if
        if (.not.(ZSNOW(I) > 0.0)) ALBSNO(I) = 0.50
        RHOSNO(I) = (RHOSNO(I) * ZSNOW(I) + RHOICE * ZFREZ) / (ZSNOW(I) + &
                    ZFREZ)
        ZSNOW(I) = ZSNOW(I) + ZFREZ
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                    (RHOW * ZSNOW(I))
        ZPOND(I) = 0.0
        TPOND(I) = 0.0
      end if
      !>
      !! At the end of the subroutine, the internal energy
      !! calculations are completed, and ZFREZ is used to
      !! update the diagnostic variable HMFG describing the energy
      !! associated with phase changes of water in soil
      !! layers, and also the diagnostic variables WTRS and WTRG
      !! describing transfers of water into or out of the
      !! snow and soil respectively. Finally, the initial step in
      !! the calculation of the internal energy change for the
      !! ponded water over the following subroutines waterFlowInfiltrate and
      !! waterFlowNonInfiltrate is performed (the calculation is
      !! completed in subroutine waterUpdates).
      !!
      HTC (I,1) = HTC (I,1) + FI(I) * HCPW * (TPOND(I) + TFREZ) * &
                  ZPOND(I) / DELT
      HMFG(I,1) = HMFG(I,1) - FI(I) * CLHMLT * RHOICE * ZFREZ / DELT
      WTRS(I) = WTRS(I) + FI(I) * ZFREZ * RHOICE / DELT
      WTRG(I) = WTRG(I) - FI(I) * ZFREZ * RHOICE / DELT
      HTCS(I) = HTCS(I) + FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * ZSNOW(I) / &
                DELT
    end if
    if (FI(I) > 0. .and. ISAND(I,1) > - 4) then
      HTC (I,1) = HTC (I,1) - FI(I) * HCPW * (TPOND(I) + TFREZ) * &
                  ZPOND(I) / DELT
    end if
  end do ! loop 100
  !
  return
end subroutine pondedWaterFreeze
