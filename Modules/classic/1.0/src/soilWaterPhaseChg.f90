!> \file
!! Check for freezing or thawing of liquid or frozen water in the soil
!! layers, and adjust layer temperatures and water stores accordingly.
!! @author D. Verseghy, M. Lazare, Y. Delage
!
subroutine soilWaterPhaseChg (TBAR, THLIQ, THICE, HCP, TBARW, HMFG, HTC, & ! Formerly TWCALC
                              FI, EVAP, THPOR, THLMIN, HCPS, DELZW, &
                              DELZZ, ISAND, IG, ILG, IL1, IL2, JL)
  !
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * MAY 16/03 - Y.DELAGE/D.VERSEGHY. BUGFIX IN FREEZING/
  !     *                                  THAWING CALCULATIONS
  !     *                                  (PRESENT SINCE V.2.7)
  !     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
  !     * JUN 20/97 - D.VERSEGHY. COSMETIC REARRANGEMENT OF
  !     *                         SUBROUTINE CALL.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
  !     *                         ORGANIC MATTER CONTENT.
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         REVISE CALCULATION OF HTC.
  !     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  NEW DIAGNOSTIC FIELDS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. ADJUST SOIL LAYER TEMPERATURES
  !     *                         AND LIQUID/FROZEN MOISTURE CONTENTS
  !     *                         FOR FREEZING/THAWING.
  !
  use classicParams, only : DELT, TFREZ, HCPW, HCPICE, HCPSND, &
                                    RHOW, RHOICE, CLHMLT

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IG, ILG, IL1, IL2, JL
  integer :: I, J
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: TBAR  (ILG,IG)    !< Temperature of soil layer \f$[C] (T_g)\f$
  real, intent(inout) :: THLIQ (ILG,IG)    !< Volumetric liquid water content of soil layer \f$[m^3 m^{-3}] (\theta_l)\f$
  real, intent(inout) :: THICE (ILG,IG)    !< Volumetric frozen water content of soil layer \f$[m^3 m^{-3}] (\theta_i)\f$
  real, intent(inout) :: HCP   (ILG,IG)    !< Heat capacity of soil layer \f$[J m^{-3} K^{-1}] (C_g)\f$
  real, intent(out)   :: TBARW (ILG,IG)    !< Temperature of water in soil layer [C]
  real, intent(inout) :: HMFG  (ILG,IG)    !< Energy associated with freezing or thawing of water in soil layer \f$[W m^{-2}]\f$
  real, intent(inout) :: HTC   (ILG,IG)    !< Internal energy change of soil layer due to conduction and/or
  !< change in mass \f$[W m^{-2}] (I_g)\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI    (ILG)       !< Fractional coverage of subarea in question on modelled area \f$[ ] (X_i)\f$
  real, intent(in) :: EVAP  (ILG)       !< Calculated evaporation rate from soil surface \f$[m s^{-1}]\f$
  real, intent(in) :: THPOR (ILG,IG)    !< Pore volume in soil layer \f$[mm] (\theta_p)\f$
  real, intent(in) :: THLMIN(ILG,IG)    !< Residual soil liquid water content remaining after
  !< freezing or evaporation \f$[m^3 m^{-3}] (\theta_r)\f$
  real, intent(in) :: HCPS  (ILG,IG)    !< Heat capacity of soil material \f$[J m^{-3} K^{-1}] (C_m)\f$
  real, intent(in) :: DELZW (ILG,IG)    !< Permeable thickness of soil layer \f$[m] (\Delta z_{g,w})\f$
  real, intent(in) :: DELZZ (ILG,IG)    !< Soil layer thicknesses to bottom of permeable depth for standard
  !< three-layer configuration,or to bottom of thermal depth for multiple
  !< layers \f$[m] (\Delta z_{g,z})\f$
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: THFREZ, THEVAP, HADD, THMELT
  !
  !-----------------------------------------------------------------------
  !>
  !!   The adjustments of soil layer temperature and water content in this
  !!   routine are done over the whole soil profile in the case of multiple
  !!   soil layers (see the section on assignment of background data), but
  !!   only to the bottom of the permeable depth in the case of the
  !!   standard three-layer configuration (0.10, 0.25 and 3.75 m). This is
  !!   because if the permeable depth lies within the thick third soil
  !!   layer, it is recognized as desirable to apply the temperature
  !!   changes only to that upper part of the layer in which the phase
  !!   change is occurring, in order to avoid systematic damping of the
  !!   temperature response of the layer. Thus the local array DELZZ (set
  !!   in subroutine waterCalcPrep) is used here instead of DELZ when referring to
  !!   the total thermal thickness of the soil layer, where DELZZ=DELZW for
  !!   the third soil layer when the three-layer configuration is being
  !!   used, but DELZZ=DELZ for all other cases.
  !!
  do J = 1,IG ! loop 100
    do I = IL1,IL2
      if (FI(I) > 0. .and. DELZW(I,J) > 0. .and. ISAND(I,1) > - 4) &
          then
        !>
        !! The heat capacity \f$C_g\f$ of the permeable part \f$\Delta z_{g, w}\f$ of the soil
        !! layer under consideration is calculated here and in various
        !! other places in the subroutine as the weighted average of the
        !! heat capacities of the liquid water content \f$\theta_l\f$, the frozen
        !! water content \f$\theta_i\f$, and the soil material (taken to apply to
        !! the volume fraction not occupied by the pore volume \f$\theta_p\f$).
        !! The heat capacity of air in the pores is neglected:
        !!
        !! \f$C_g = C_w \theta_l + C_i \theta_i + C_m (1 - \theta_p)\f$
        !!
        !! Over the impermeable portion of the layer, the heat capacity
        !! of rock \f$C_r\f$ is assumed to apply. Thus an effective heat
        !! capacity \f$C_{g, e}\f$ (in units of \f$J m^{-2} K^{-1}\f$) over the soil layer in
        !! question, \f$\Delta z_{g, z}\f$, can be calculated as:
        !!
        !! \f$C_{g, e} = C_g \Delta z_{g, w} + C_r(\Delta z_{g, z} - \Delta z_{g, w})\f$
        !!
        !! The change of internal energy I in the soil layers as a
        !! result of freezing and thawing is calculated as the
        !! difference in I between the beginning and end of the
        !! subroutine:
        !!
        !! \f$\Delta I_j = X_i \Delta (C_{g, e} T_g)/\Delta t\f$
        !!
        !! where \f$T_g\f$ is the temperature of the layer, \f$\Delta t\f$ the
        !! length of the time step, and \f$X_i\f$ the fractional coverage of
        !! the subarea under consideration relative to the modelled
        !! area.
        !!
        HCP  (I,J) = HCPW * THLIQ(I,J) + HCPICE * THICE(I,J) + &
                     HCPS(I,J) * (1. - THPOR(I,J))
        HTC  (I,J) = HTC(I,J) - FI(I) * (HCP(I,J) * DELZW(I,J) + &
                     HCPSND * (DELZZ(I,J) - DELZW(I,J))) * &
                     (TBAR(I,J) + TFREZ) / DELT
        !>
        !! If the soil layer temperature is less than 0 C and the
        !! volumetric liquid water content \f$\theta_l\f$ of the layer is
        !! greater than the residual water content \f$\theta_f\f$, the water
        !! content THFREZ that can be frozen by the available energy
        !! sink is calculated from \f$C_e\f$ and \f$T_g\f$. The volumetric water
        !! content THEVAP of the first layer that is required to
        !! satisfy the surface evaporative flux is determined. For
        !! each layer, if THLIQ is found to exceed THLMIN + THEVAP,
        !! THFREZ is compared to the available water. If THFREZ \f$\leq\f$
        !! THLIQ – THLMIN – THEVAP, all of the available energy sink
        !! is used to freeze part of the liquid water content in
        !! the permeable part of the soil layer, the amount of
        !! energy involved is subtracted from HTC and added to HMFG,
        !! \f$C_g\f$ is recalculated and the layer temperature is set to
        !! 0 C. Otherwise, all of the liquid water content of the
        !! layer above THLMIN + THEVAP is converted to frozen water,
        !! and HMFG and HTC are recalculated to reflect this. Then
        !! \f$C_g\f$ is recomputed, and the remaining energy sink is
        !! applied to decreasing the temperature of the soil layer
        !! (both the permeable and impermeable portions) using \f$C_e\f$.
        !!
        if (TBAR(I,J) < 0. .and. THLIQ(I,J) > THLMIN(I,J)) then
          THFREZ = - (HCP(I,J) * DELZW(I,J) + HCPSND * (DELZZ(I,J) - &
                   DELZW(I,J))) * TBAR(I,J) / (CLHMLT * RHOW * &
                   DELZW(I,J))
          if (J == 1) then
            THEVAP = EVAP(I) * DELT / DELZW(I,J)
          else
            THEVAP = 0.0
          end if
          if ((THLIQ(I,J) - THLMIN(I,J) - THEVAP) > 0.0) then
            if (THFREZ <= (THLIQ(I,J) - THLMIN(I,J) - THEVAP)) then
              HMFG(I,J) = HMFG(I,J) - FI(I) * THFREZ * CLHMLT * &
                          RHOW * DELZW(I,J) / DELT
              HTC(I,J) = HTC(I,J) - FI(I) * THFREZ * CLHMLT * &
                         RHOW * DELZW(I,J) / DELT
              THLIQ(I,J) = THLIQ(I,J) - THFREZ
              THICE(I,J) = THICE(I,J) + THFREZ * RHOW / RHOICE
              HCP  (I,J) = HCPW * THLIQ(I,J) + HCPICE * THICE(I,J) + &
                           HCPS(I,J) * (1. - THPOR(I,J))
              TBAR (I,J) = 0.0
            else
              HMFG(I,J) = HMFG(I,J) - FI(I) * (THLIQ(I,J) - &
                          THLMIN(I,J) - THEVAP) * CLHMLT * RHOW * DELZW(I,J) / &
                          DELT
              HTC(I,J) = HTC(I,J) - FI(I) * (THLIQ(I,J) - THLMIN(I,J) - &
                         THEVAP) * CLHMLT * RHOW * DELZW(I,J) / DELT
              HADD = (THFREZ - (THLIQ(I,J) - THLMIN(I,J) - THEVAP)) * &
                     CLHMLT * RHOW * DELZW(I,J)
              THICE(I,J) = THICE(I,J) + (THLIQ(I,J) - THLMIN(I,J) - &
                           THEVAP) * RHOW / RHOICE
              THLIQ(I,J) = THLMIN(I,J) + THEVAP
              HCP  (I,J) = HCPW * THLIQ(I,J) + HCPICE * THICE(I,J) + &
                           HCPS(I,J) * (1. - THPOR(I,J))
              TBAR (I,J) = - HADD / (HCP(I,J) * DELZW(I,J) + HCPSND * &
                           (DELZZ(I,J) - DELZW(I,J)))
            end if
          end if
        end if
        !
        !>
        !! If the soil layer temperature is greater than 0 C and the
        !! volumetric ice content \f$\theta_i\f$ of the layer is greater than
        !! zero, the ice content THMELT that can be melted by the
        !! available energy is calculated from \f$C_e\f$ and \f$T_g\f$. For
        !! each layer, if THMELT \f$\leq\f$ THICE, all of the available
        !! energy is used to melt part of the frozen water content
        !! of the permeable part of the layer, the amount of energy
        !! involved is subtracted from HTC and added to HMFG, \f$C_g\f$ is
        !! recalculated and the layer temperature is set to 0 C.
        !! Otherwise, all of the frozen water content of the layer
        !! is converted to liquid water, and HMFG and HTC are
        !! recalculated to reflect this. Then \f$C_g\f$ is recomputed,
        !! and the remaining energy is applied to increasing the
        !! temperature of the soil layer (both the permeable and
        !! impermeable portions) using \f$C_e\f$.
        !!
        if (TBAR(I,J) > 0. .and. THICE(I,J) > 0.) then
          THMELT = (HCP(I,J) * DELZW(I,J) + HCPSND * (DELZZ(I,J) - &
                   DELZW(I,J))) * TBAR(I,J) / (CLHMLT * RHOICE * &
                   DELZW(I,J))
          if (THMELT <= THICE(I,J)) then
            HMFG(I,J) = HMFG(I,J) + FI(I) * THMELT * CLHMLT * &
                        RHOICE * DELZW(I,J) / DELT
            HTC(I,J) = HTC(I,J) + FI(I) * THMELT * CLHMLT * &
                       RHOICE * DELZW(I,J) / DELT
            THICE(I,J) = THICE(I,J) - THMELT
            THLIQ(I,J) = THLIQ(I,J) + THMELT * RHOICE / RHOW
            HCP  (I,J) = HCPW * THLIQ(I,J) + HCPICE * THICE(I,J) + &
                         HCPS(I,J) * (1. - THPOR(I,J))
            TBAR (I,J) = 0.0
          else
            HMFG(I,J) = HMFG(I,J) + FI(I) * THICE(I,J) * CLHMLT * &
                        RHOICE * DELZW(I,J) / DELT
            HTC(I,J) = HTC(I,J) + FI(I) * THICE(I,J) * CLHMLT * &
                       RHOICE * DELZW(I,J) / DELT
            HADD = (THMELT - THICE(I,J)) * CLHMLT * RHOICE * &
                   DELZW(I,J)
            THLIQ(I,J) = THLIQ(I,J) + THICE(I,J) * RHOICE / RHOW
            THICE(I,J) = 0.0
            HCP  (I,J) = HCPW * THLIQ(I,J) + HCPICE * THICE(I,J) + &
                         HCPS(I,J) * (1. - THPOR(I,J))
            TBAR (I,J) = HADD / (HCP(I,J) * DELZW(I,J) + HCPSND * &
                         (DELZZ(I,J) - DELZW(I,J)))
          end if
        end if
        !>
        !! In the final cleanup, the internal energy calculations
        !! for this subroutine are completed, and the first half of
        !! a new set of internal energy calculations is done to span
        !! the subroutines treating ground water movement, which
        !! will be completed in subroutine waterUpdates. Lastly, TBARW,
        !! the liquid water temperature of each soil layer, is
        !! assigned using TBAR.
        !!
        HTC  (I,J) = HTC(I,J) + FI(I) * (HCP(I,J) * DELZW(I,J) + &
                     HCPSND * (DELZZ(I,J) - DELZW(I,J))) * &
                     (TBAR(I,J) + TFREZ) / DELT
        HTC(I,J) = HTC(I,J) - FI(I) * (TBAR(I,J) + TFREZ) * &
                   (HCPW * THLIQ(I,J) + HCPICE * THICE(I,J)) * &
                   DELZW(I,J) / DELT
      end if
      TBARW(I,J) = TBAR(I,J)
    end do
  end do ! loop 100
  !
  return
end subroutine soilWaterPhaseChg
