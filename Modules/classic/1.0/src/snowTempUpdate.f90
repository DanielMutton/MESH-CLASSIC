!> \file
!! Snow temperature calculations and cleanup after surface
!! energy budget calculations.
!! @author D. Verseghy, M. Lazare
!
subroutine snowTempUpdate (GSNOW, TSNOW, WSNOW, RHOSNO, QMELTG, & ! Formerly TSPOST
                           GZERO, TSNBOT, HTCS, HMFN, &
                           GCONSTS, GCOEFFS, GCONST, GCOEFF, TBAR, &
                           TSURF, ZSNOW, TCSNOW, HCPSNO, QTRANS, &
                           FI, DELZ, ILG, IL1, IL2, JL, IG)
  !
  !     * AUG 16/06 - D.VERSEGHY. MAJOR REVISION TO IMPLEMENT THERMAL
  !     *                         SEPARATION OF SNOW AND SOIL.
  !     * MAR 23/06 - D.VERSEGHY. ADD CALCULATIONS TO ALLOW FOR WATER
  !     *                         FREEZING IN SNOWPACK.
  !     * OCT 04/05 - D.VERSEGHY. MODIFY 300 LOOP FOR CASES WHERE IG>3.
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUN 17/02 - D.VERSEGHY. RESET PONDED WATER TEMPERATURE
  !     *                         USING CALCULATED GROUND HEAT FLUX
  !     *                         SHORTENED CLASS4 COMMON BLOCK.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         INCORPORATE EXPLICITLY CALCULATED
  !     *                         THERMAL CONDUCTIVITIES AT TOPS AND
  !     *                         BOTTOMS OF SOIL LAYERS, AND
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         FIX BUG IN CALCULATION OF FLUXES
  !     *                         BETWEEN SOIL LAYERS (PRESENT SINCE
  !     *                         RELEASE OF CLASS VERSION 2.5).
  !     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         COMPLETION OF ENERGY BALANCE
  !     *                         DIAGNOSTICS.
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         REVISE CALCULATION OF TBARPR(I, 1).
  !     * APR 10/92 - M.LAZARE.   CLASS - VERSION 2.2.
  !     *                         DIVIDE PREVIOUS SUBROUTINE "T4LAYR" INTO
  !     *                         "snowHeatCond" AND "snowTempUpdate" AND VECTORIZE.
  !     * APR 11/89 - D.VERSEGHY. CALCULATE HEAT FLUXES BETWEEN SNOW/SOIL
  !     *                         LAYERS; CONSISTENCY CHECK ON CALCULATED
  !     *                         SURFACE LATENT HEAT OF MELTING/
  !     *                         FREEZING; STEP AHEAD SNOW LAYER
  !     *                         TEMPERATURE AND ASSIGN EXCESS HEAT TO
  !     *                         MELTING IF NECESSARY; DISAGGREGATE
  !     *                         FIRST SOIL LAYER TEMPERATURE INTO
  !     *                         PONDED WATER AND SOIL TEMPERATURES
  !     *                         ADD SHORTWAVE RADIATION TRANSMITTED
  !     *                         THROUGH SNOWPACK TO HEAT FLUX AT TOP
  !     *                         OF FIRST SOIL LAYER; CONVERT LAYER
  !     *                         TEMPERATURES TO DEGREES C.
  !
  use classicParams, only : DELT, TFREZ, HCPW, HCPICE, HCPSOL, &
                                  HCPOM, HCPSND, HCPCLY, RHOW, &
                                  RHOICE, CLHMLT

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL, IG
  integer :: I, J
  !
  !     * OUTPUT ARRAYS.
  !
  real,intent(inout) :: GZERO (ILG)  !< Heat conduction into soil surface \f$[W m^{-2}] (G(\Delta z_s))\f$
  real,intent(inout) :: TSNBOT(ILG)  !< Temperature at bottom of snow pack [K]
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real,intent(inout) :: GSNOW (ILG)  !< Heat conduction into surface of snow pack \f$[W m^{-2}] (G(0))\f$
  real,intent(inout) :: TSNOW (ILG)  !< Snowpack temperature \f$[K/C] (T_s)\f$
  real,intent(inout) :: WSNOW (ILG)  !< Liquid water content of snow pack \f$[kg m^{-2}] (w_s)\f$
  real,intent(inout) :: RHOSNO(ILG)  !< Density of snow \f$[kg m^{-3}] (\rho_s)\f$
  real,intent(inout) :: QMELTG(ILG)  !< Available energy to be applied to melting of snow \f$[W m^{-2}]\f$
  real,intent(inout) :: HTCS  (ILG)  !< Internal energy change of snow pack due to
  !< conduction and/or change in mass \f$[W m^{-2}] (I_s)\f$
  real,intent(inout) :: HMFN  (ILG)  !< Energy associated with phase change of water in
  !< snow pack \f$[W m^{-2}]\f$
  !
  !     * INPUT ARRAYS.
  !
  real,intent(in) :: TSURF (ILG)  !< Snow surface temperature [K]
  real,intent(in) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (\Delta z_s)\f$
  real,intent(in) :: TCSNOW(ILG)  !< Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
  real,intent(inout) :: HCPSNO(ILG)  !< Heat capacity of snow \f$[J m^{-3} K^1] (C_s)\f$
  real,intent(in) :: QTRANS(ILG)  !< Shortwave radiation transmitted through the
  !< snow pack \f$[W m^{-2}]\f$
  real,intent(in) :: GCONST(ILG)  !< Intercept used in equation relating snow
  !< surface heat flux to snow surface temperature \f$[W m^{-2}]\f$
  real,intent(in) :: GCOEFF(ILG)  !< Multiplier used in equation relating snow
  !< surface heat flux to snow surface temperature \f$[W m^{-2} K^{-1}]\f$
  real,intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on
  !< modelled area \f$[ ] (X_i)\f$
  real,intent(in) :: GCONSTS(ILG) !< Intercept used in equation relating snow
  !< surface heat flux to snow surface temperature \f$[W m^{-2}]\f$
  real,intent(in) :: GCOEFFS(ILG) !< Multiplier used in equation relating snow
  !< surface heat flux to snow surface temperature \f$[W m^{-2} K^{-1}]\f$
  real,intent(in) :: TBAR(ILG,IG) !< Temperatures of soil layers, averaged over
  !< modelled area [K]
  real,intent(in) :: DELZ  (IG)   !< Overall thickness of soil layer [m]
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: GSNOLD, HADD, HCONV, WFREZ
  !
  !-----------------------------------------------------------------------
  !
  !>
  !! In the 100 loop, the heat flux into the snow surface (without
  !! adjustments that may have been applied relating to partitioning
  !! of the residual of the surface energy balance among the surface
  !! flux terms) is calculated from the snow surface temperature
  !! TSURF, using the GCOEFFS and GCONSTS terms (see documentation of
  !! subroutine snowHeatCond). The temperature at the bottom of the snow
  !! pack, TSNBOT, is then calculated. Currently TSNBOT is determined
  !! as a simple average of the temperatures of the snow and the first
  !! soil layer, weighted according to their respective depths (and
  !! constrained to be \f$\leq\f$ 0 C), but this is under review. The heat
  !! flux into the soil surface is then evaluated from TSNBOT and the
  !! GCOEFF and GCONST terms (see documentation of subroutine soilHeatFluxPrep).
  !! If the energy to be applied to the melting of snow, QMELTG, is
  !! negative (indicating an energy sink), QMELTG is added to the heat
  !! flux into the ground, GZERO, and reset to zero. The temperature
  !! of the snow pack is then stepped forward using the heat fluxes at
  !! the top and bottom of the snow pack, G(0) and \f$G(\Delta z_s)\f$:
  !!
  !! \f$\Delta T_s = [G(0) - G(\Delta z_s)] \Delta t /(C_s \Delta z_s)\f$
  !!
  !! where \f$C_s\f$ is the snow heat capacity, \f$\Delta t\f$ the time step and
  !! \f$\Delta z_s\f$ the snow depth. If the new snow temperature is greater than
  !! zero, the excess amount of heat is calculated and added to QMELTG
  !! and subtracted from GSNOW, and TSNOW is reset to 0 C. Finally,
  !! the shortwave radiation transmitted through the snow pack,
  !! QTRANS, is added to GZERO.
  !!
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0.) then
      GSNOLD = GCOEFFS(I) * TSURF(I) + GCONSTS(I)
      TSNBOT(I) = (ZSNOW(I) * TSNOW(I) + DELZ(1) * TBAR(I,1)) / &
                  (ZSNOW(I) + DELZ(1))
      !              TSNBOT(I)=0.90*TSNOW(I)+0.10*TBAR(I,1)
      !              TSNBOT(I)=TSURF(I)-GSNOLD*ZSNOW(I)/(2.0*TCSNOW(I))
      TSNBOT(I) = MIN(TSNBOT(I),TFREZ)
      GZERO(I) = GCOEFF(I) * TSNBOT(I) + GCONST(I)
      if (QMELTG(I) < 0.) then
        GSNOW(I) = GSNOW(I) + QMELTG(I)
        QMELTG(I) = 0.
      end if
      TSNOW(I) = TSNOW(I) + (GSNOW(I) - GZERO(I)) * DELT / &
                 (HCPSNO(I) * ZSNOW(I)) - TFREZ
      if (TSNOW(I) > 0.) then
        QMELTG(I) = QMELTG(I) + TSNOW(I) * HCPSNO(I) * ZSNOW(I) / DELT
        GSNOW(I) = GSNOW(I) - TSNOW(I) * HCPSNO(I) * ZSNOW(I) / DELT
        TSNOW(I) = 0.
      end if
      GZERO(I) = GZERO(I) + QTRANS(I)
    end if
  end do ! loop 100
  !
  !>
  !! In the 200 loop, since liquid water is assumed only to exist in
  !! the snow pack if it is at 0 C, a check is carried out to
  !! determine whether the liquid water content WSNOW > 0 at the same
  !! time as the snow temperature TSNOW < 0. If so, the change of
  !! internal energy \f$I_s\f$ of the snow pack as a result of this phase
  !! change is calculated as the difference in \f$I_s\f$ between the
  !! beginning and end of the loop:
  !!
  !! \f$\Delta I_s = X_i \Delta [C_s T_s]/ \Delta t\f$
  !!
  !! where \f$X_i\f$ represents the fractional coverage of the subarea under
  !! consideration relative to the modelled area. The total energy
  !! sink HADD available to freeze liquid water in the snow pack is
  !! calculated from TSNOW, and the amount of energy HCONV required to
  !! freeze all the available water is calculated from WSNOW. If
  !! HADD < HCONV, only part of WSNOW is frozen; this amount WFREZ is
  !! calculated from HADD and subtracted from WSNOW, the snow
  !! temperature is reset to 0 C, the frozen water is used to update
  !! the snow density, and the snow heat capacity is recalculated:
  !!
  !! \f$C_s = C_i [\rho_s /\rho_i] + C_w w_s/[\rho_w \Delta z_s]\f$
  !!
  !! where \f$C_i\f$ and \f$C_w\f$ are the heat capacities of ice and water
  !! respectively, \f$w_s\f$ is the snow water content and \f$\rho_s\f$, \f$\rho_i\f$
  !! and \f$\rho_w\f$ are the densities of snow, ice and water respectively.
  !! If HADD > HCONV, the available energy sink is sufficient to
  !! freeze all of WSNOW. HADD is recalculated as HADD – HCONV, WFREZ
  !! is set to WSNOW and added to the snow density, WSNOW is set to
  !! zero, the snow heat capacity is recalculated and HADD is used to
  !! determine a new value of TSNOW. Finally, WFREZ is used to update
  !! the diagnostic variables HMFN describing phase changes of water
  !! in the snow pack, and the change in internal energy HTCS.
  !!
  do I = IL1,IL2 ! loop 200
    if (FI(I) > 0. .and. TSNOW(I) < 0. .and. WSNOW(I) > 0.) &
        then
      HTCS(I) = HTCS(I) - FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * ZSNOW(I) / &
                DELT
      HADD = - TSNOW(I) * HCPSNO(I) * ZSNOW(I)
      HCONV = CLHMLT * WSNOW(I)
      if (HADD <= HCONV) then
        WFREZ = HADD / CLHMLT
        HADD = 0.0
        WSNOW(I) = MAX(0.0,WSNOW(I) - WFREZ)
        TSNOW(I) = 0.0
        RHOSNO(I) = RHOSNO(I) + WFREZ / ZSNOW(I)
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE + HCPW * WSNOW(I) / &
                    (RHOW * ZSNOW(I))
      else
        HADD = HADD - HCONV
        WFREZ = WSNOW(I)
        WSNOW(I) = 0.0
        RHOSNO(I) = RHOSNO(I) + WFREZ / ZSNOW(I)
        HCPSNO(I) = HCPICE * RHOSNO(I) / RHOICE
        TSNOW(I) = - HADD / (HCPSNO(I) * ZSNOW(I))
      end if
      HMFN(I) = HMFN(I) - FI(I) * CLHMLT * WFREZ / DELT
      HTCS(I) = HTCS(I) - FI(I) * CLHMLT * WFREZ / DELT
      HTCS(I) = HTCS(I) + FI(I) * HCPSNO(I) * (TSNOW(I) + TFREZ) * ZSNOW(I) / &
                DELT
    end if
  end do ! loop 200
  !
  return
end subroutine snowTempUpdate
