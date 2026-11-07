!> \file
!! Soil heat flux calculations and cleanup after surface
!! energy budget calculations.
!! @author D. Verseghy, M. Lazare
!
subroutine soilHeatFluxCleanup (TBARPR, G12, G23, TPOND, GZERO, QFREZG, GCONST, & ! Formerly TNPOST
                                GCOEFF, TBAR, TCTOP, TCBOT, HCP, ZPOND, TSURF, &
                                TBASE, TBAR1P, A1, A2, B1, B2, C2, FI, IWATER, &
                                ISAND, DELZ, DELZW, ILG, IL1, IL2, JL, IG)
  !
  !     * NOV 01/06 - D.VERSEGHY. ALLOW PONDING ON ICE SHEETS.
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
  !     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         FIX BUG IN CALCULATION OF FLUXES
  !     *                         BETWEEN SOIL LAYERS (PRESENT SINCE
  !     *                         FIRST RELEASE OF VERSION 2.5).
  !     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
  !     *                         REVISE CALCULATION OF TBARPR(I, 1).
  !     * APR 10/92 - M.LAZARE.   CLASS - VERSION 2.2.
  !     *                         DIVIDE PREVIOUS SUBROUTINE "T3LAYR" INTO
  !     *                         "soilHeatFluxPrep" AND "soilHeatFluxCleanup" AND VECTORIZE.
  !     * APR 11/89 - D.VERSEGHY. CALCULATE HEAT FLUXES BETWEEN SOIL
  !     *                         LAYERS; DISAGGREGATE FIRST SOIL LAYER
  !     *                         TEMPERATURE INTO PONDED WATER AND
  !     *                         SOIL TEMPERATURES; CONSISTENCY CHECK
  !     *                         ON CALCULATED SURFACE LATENT HEAT OF
  !     *                         MELTING/FREEZING; CONVERT SOIL LAYER
  !     *                         TEMPERATURES TO DEGREES C.
  !
  use classicParams, only : TFREZ, HCPW, HCPSND

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL, IG
  integer :: I, J
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out) :: TBARPR(ILG,IG) !< Temperatures of soil layers for subarea [C]
  !
  real, intent(inout) :: G12   (ILG)    !< Heat conduction between first and second soil layers \f$[W m^{-2}] (G(\Delta z))\f$
  real, intent(out) :: G23   (ILG)    !< Heat conduction between second and third soil layers \f$[W m^{-2}]\f$
  real, intent(inout) :: TPOND (ILG)    !< Temperature of ponded water \f$[C] (T_p)\f$
  !
  !     * INPUT/OUTPUT ARRAYS.
  !
  real, intent(inout) :: GZERO (ILG)    !< Heat conduction into soil surface \f$[W m^{-2}] (G(0))\f$
  real, intent(inout) :: QFREZG(ILG)    !< Energy sink to be applied to freezing of ponded water \f$[W m^{-2}]\f$
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: TBAR  (ILG,IG) !< Temperatures of soil layers, averaged over modelled area [K]
  real, intent(in) :: TCTOP (ILG,IG) !< Thermal conductivity of soil at top of layer \f$[W m^{-1} K^{-1}] (\lambda)\f$
  real, intent(in) :: TCBOT (ILG,IG) !< Thermal conductivity of soil at bottom of layer \f$[W m^{-1} K^{-1}] (\lambda)\f$
  real, intent(in) :: HCP   (ILG,IG) !< Heat capacity of soil layer \f$[J m^{-3} K^{-1}]\f$
  real, intent(in) :: DELZW (ILG,IG) !< Permeable thickness of soil layer [m]
  !
  real, intent(in) :: ZPOND (ILG)    !< Depth of ponded water on surface \f$[m] (\Delta z_p)\f$
  real, intent(in) :: TSURF (ILG)    !< Ground surface temperature [K]
  real, intent(in) :: TBASE (ILG)    !< Temperature of bedrock in third soil layer [K]
  real, intent(in) :: TBAR1P(ILG)    !< Lumped temperature of ponded water and first soil layer \f$[K] (T_{1p})\f$
  real, intent(in) :: A1    (ILG)    !< Work array used in calculation of GCONST and GCOEFF
  real, intent(in) :: A2    (ILG)    !< Work array used in calculation of GCONST and GCOEFF
  real, intent(in) :: B1    (ILG)    !< Work array used in calculation of GCONST and GCOEFF
  real, intent(in) :: B2    (ILG)    !< Work array used in calculation of GCONST and GCOEFF
  real, intent(in) :: C2    (ILG)    !< Work array used in calculation of GCONST and GCOEFF
  real, intent(in) :: FI    (ILG)    !< Fractional coverage of subarea in question on modelled area [ ]
  real, intent(in) :: GCONST(ILG)    !< Intercept used in equation relating ground
  !< surface heat flux to surface temperature \f$[W m^{-2}]\f$
  real, intent(in) :: GCOEFF(ILG)    !< Multiplier used in equation relating ground
  !< surface heat flux to surface temperature \f$[W m^{-2} K^{-1}]\f$
  !
  integer, intent(in) :: IWATER(ILG) !< Flag indicating condition of surface (dry, water-covered or snow-covered)
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  real, intent(in) :: DELZ  (IG)     !< Overall thickness of soil layer \f$[m] (\Delta z)\f$
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: GZROLD, DELZ1
  !
  !-----------------------------------------------------------------------
  !
  !>
  !! In the 100 loop, the heat flux into the ground (without
  !! adjustments that may have been applied owing to phase changes of
  !! water at the surface or partitioning the residual of the surface
  !! energy balance among the surface flux terms) is calculated from
  !! the ground surface temperature TSURF, using the GCOEFF and GCONST
  !! terms (see documentation of subroutine soilHeatFluxPrep). This flux is then
  !! used to back-calculate the heat conduction between the first and
  !! second, and the second and third soil layers.
  !!
  !! If the depth of water ponded on the surface is greater than zero,
  !! the total thickness of the lumped layer consisting of the first
  !! soil layer and the ponded water, \f$\Delta z_{1p}\f$, is calculated as the
  !! sum of the two depths. The temperature of the ponded water over
  !! the subarea in question is disaggregated from the temperature of
  !! the lumped layer, \f$T_{1p}\f$, by making use of the calculated heat
  !! fluxes at the top and bottom of the layer, and the assumption
  !! (discussed in the documentation of subroutine soilHeatFluxPrep) that that
  !! the variation of temperature T with depth z within a soil layer
  !! can be modelled using a quadratic equation:
  !!
  !! \f$T(z) = 1/2 T''(0)z^2 + T'(0)z +T(0)\f$
  !!
  !! Integrating this equation separately over the ponded water depth
  !! \f$\Delta z_p\f$ and over the thickness of the lumped layer produces
  !! expressions for the ponded water temperature \f$T_p\f$ and the
  !! temperature of the lumped layer. Making use of the fact that
  !!
  !! \f$T''(0) = [T'(\Delta z) - T'(0)]/\Delta z\f$
  !!
  !! where \f$\Delta z\f$ is a depth interval, and
  !!
  !! \f$G(z) = - \lambda(z)T'(z)\f$
  !! where \f$\lambda\f$ represents the thermal conductivity, an expression
  !! for \f$T_p\f$ can be derived:
  !!
  !! \f$T_p = [G(0)/ \lambda(0) – G(\Delta z_{1p})/ \lambda(\Delta z_{1p})]
  !! [\Delta z_p^2 - \Delta z_{1p}^2]/ 6 \Delta z_{1p} –
  !! G(0) [\Delta z_p - \Delta z_{1p}]/ 2 \lambda(0) – T_{1P}\f$
  !!
  !! The temperature TBARPR of the first soil layer over the subarea
  !! in question can then be obtained by disaggregating the ponded
  !! water temperature from the lumped layer temperature using the
  !! respective heat capacities of the ponded water and the soil
  !! layer. (The heat capacity of the soil is determined as the
  !! weighted average of HCP over the permeable thickness DELZW, and
  !! the heat capacity of rock, HCPSND, over the impermeable
  !! thickness, DELZ-DELZW.) Both the ponded water temperature and the
  !! soil layer temperature are converted to C. Lastly, if there is
  !! ponded water on the surface (IWATER = 1) and QFREZG, the surface
  !! energy available for phase change of water, is positive
  !! (indicating an energy source), or if there is snow on the ground
  !! (IWATER = 2) and QFREZG is negative (indicating an energy sink),
  !! or if there is no liquid or frozen water on the surface
  !! (IWATER = 0), QFREZG is added to the heat flux into the ground,
  !! GZERO, and then reset to zero.
  !!
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0.) then
      GZROLD = GCOEFF(I) * TSURF(I) + GCONST(I)
      G12(I) = (TSURF(I) - TBAR1P(I) - A1(I) * GZROLD) / B1(I)
      G23(I) = (TSURF(I) - TBAR(I,2) - A2(I) * GZROLD - B2(I) * G12(I)) / &
               C2(I)
      if (ZPOND(I) > 0.) then
        DELZ1 = DELZ(1) + ZPOND(I)
        TPOND(I) = (GZROLD / TCTOP(I,1) - G12(I) / TCBOT(I,1)) * &
                   (ZPOND(I) * ZPOND(I) - DELZ1 * DELZ1) / (6.0 * DELZ1) - &
                   GZROLD * (ZPOND(I) - DELZ1) / (2.0 * TCTOP(I,1)) + &
                   TBAR1P(I) - TFREZ
        TBARPR(I,1) = ((HCP(I,1) * DELZW(I,1) + HCPSND * (DELZ(1) - &
                      DELZW(I,1)) + HCPW * ZPOND(I)) * TBAR1P(I) - &
                      HCPW * ZPOND(I) * (TPOND(I) + TFREZ)) / &
                      (HCP(I,1) * DELZW(I,1) + HCPSND * (DELZ(1) - &
                      DELZW(I,1))) - TFREZ
      else
        TPOND(I) = 0.
        TBARPR(I,1) = TBAR(I,1) - TFREZ
      end if
      !
      if ((IWATER(I) == 1 .and. QFREZG(I) > 0.) .or. &
          (IWATER(I) == 2 .and. QFREZG(I) < 0.) .or. &
          IWATER(I) == 0) then
        GZERO(I) = GZERO(I) + QFREZG(I)
        QFREZG(I) = 0.
      end if
    end if
  end do ! loop 100
  !
  !>
  !! In loop 200, the subarea soil layer temperatures TBARPR are set
  !! for the remaining soil layers. In all cases the temperature of
  !! the layer is set to that for the modelled area, TBAR, converted
  !! to C, except in the case of the third soil layer if the standard
  !! three-layer configuration is being modelled (with a very thick
  !! third soil layer of 3.75 m). In this case TBARPR and the layer
  !! heat capacity HCP are considered to apply to the permeable depth
  !! DELZW of the layer, and the bedrock temperature TBASE and the
  !! rock heat capacity HCPSND apply to the remainder, DELZ-DELZW. The
  !! disaggregation of TBARPR from TBAR and TBASE is carried out on
  !! this basis. TBARPR for this layer is also converted to C.
  !!
  do I = IL1,IL2 ! loop 200
    if (FI(I) > 0.) then
      TBARPR(I,2) = TBAR(I,2) - TFREZ
      if (DELZW(I,3) > 0.0 .and. DELZW(I,3) < DELZ(3) &
          .and. IG == 3) then
        TBARPR(I,3) = (TBAR(I,3) * (HCP(I,3) * DELZW(I,3) + &
                      HCPSND * (DELZ(3) - DELZW(I,3))) - TBASE(I) * &
                      HCPSND * (DELZ(3) - DELZW(I,3))) / (HCP(I,3) * &
                      DELZW(I,3)) - TFREZ
      else
        do J = 3,IG ! loop 150
          TBARPR(I,J) = TBAR(I,J) - TFREZ
        end do ! loop 150
      end if
    end if
  end do ! loop 200
  !
  return
end subroutine soilHeatFluxCleanup
