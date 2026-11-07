!> \file
!! Calculate coefficients for solution of heat conduction
!! into soil.
!
subroutine soilHeatFluxPrep (A1, A2, B1, B2, C2, GDENOM, GCOEFF, & ! Formerly TNPREP
                             GCONST, CPHCHG, IWATER, &
                             TBAR, TCTOP, TCBOT, &
                             FI, ZPOND, TBAR1P, DELZ, TCSNOW, ZSNOW, &
                             ISAND, ILG, IL1, IL2, JL, IG)
  !
  !     * MAR 03/08 - D.VERSEGHY. ASSIGN TCTOP3 AND TCBOT3 ON THE BASIS
  !     *                         OF SUBAREA VALUES FROM energyBudgetPrep; REPLACE
  !     *                         THREE-LEVEL TEMPERATURE AND THERMAL
  !     *                         CONDUCTIVITY VECTORS WITH STANDARD
  !     *                         VALUES.
  !     * AUG 16/06 - D.VERSEGHY. REMOVE TSTART.
  !     * MAY 24/05 - D.VERSEGHY. LIMIT DELZ3 TO <= 4.1 M.
  !     * OCT 04/05 - D.VERSEGHY. USE THREE-LAYER TBAR, TCTOP, TCBOT.
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 06/02 - D.VERSEGHY. SHORTENED CLASS3 COMMON BLOCK,
  !     * JUN 17/02 - D.VERSEGHY. USE NEW LUMPED SOIL AND PONDED
  !     *                         WATER TEMPERATURE FOR FIRST LAYER
  !     *                         SHORTENED COMMON BLOCK.
  !     * MAR 28/02 - D.VERSEGHY. CHANGE POND THRESHOLD VALUE FOR
  !     *                         CALCULATION OF "IWATER".
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         INCORPORATE EXPLICITLY CALCULATED
  !     *                         THERMAL CONDUCTIVITIES AT TOPS AND
  !     *                         BOTTOMS OF SOIL LAYERS.
  !     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         SURFACE TREATED AS WATER ONLY IF
  !     *                         ZPOND > 1 MM.
  !     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
  !     *                         ORGANIC MATTER CONTENT.
  !     * NOV 28/94 - M. LAZARE.  CLASS - VERSION 2.3.
  !     *                         TCSATW, TCSATI DECLARED REAL(16).
  !     * APR 10/92 - M. LAZARE.  CLASS - VERSION 2.1.
  !     *                         DIVIDE PREVIOUS SUBROUTINE "T3LAYR"
  !     *                         INTO "soilHeatFluxPrep" AND "soilHeatFluxCleanup" AND
  !     *                         VECTORIZE.
  !     * APR 11/89 - D.VERSEGHY. CALCULATE COEFFICIENTS FOR GROUND HEAT
  !     *                         FLUX, EXPRESSED AS A LINEAR FUNCTION
  !     *                         OF SURFACE TEMPERATURE.  COEFFICIENTS
  !     *                         ARE CALCULATED FROM LAYER TEMPERATURES,
  !     *                         THICKNESSES AND THERMAL CONDUCTIVITIES,
  !     *                         ASSUMING A QUADRATIC VARIATION OF
  !     *                         TEMPERATURE WITH DEPTH WITHIN EACH
  !     *                         SOIL LAYER. SET THE SURFACE LATENT
  !     *                         HEAT OF VAPORIZATION OF WATER AND
  !     *                         THE STARTING TEMPERATURE FOR THE
  !     *                         ITERATION IN "energBalVegSolve"/"energBalNoVegSolve".

  use classicParams, only : CLHMLT, CLHVAP

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL, IG
  integer :: I, J
  !
  !     * OUTPUT ARRAYS.
  !
  ! Work arrays used in calculation of GCONST and GCOEFF
  real, intent(inout) :: A1(ILG), A2(ILG), B1(ILG), &
                         B2(ILG), C2(ILG)
  !
  real, intent(inout) :: GDENOM(ILG)    !< Work array used in calculation of GCONST and GCOEFF
  real, intent(out)   :: GCOEFF(ILG)    !< Multiplier used in equation relating ground
  !< surface heat flux to surface temperature \f$[W m^{-2} K^{-1}]\f$
  real, intent(out)   :: GCONST(ILG)    !< Intercept used in equation relating ground
  !< surface heat flux to surface temperature \f$[W m^{-2}]\f$
  real, intent(out)   :: CPHCHG(ILG)    !< Latent heat of sublimation \f$[J kg^{-1}]\f$
  !
  integer, intent(inout) :: IWATER(ILG) !< Flag indicating condition of surface (dry, water-covered or snow-covered)
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: TBAR  (ILG,IG) !< Temperatures of soil layers, averaged over modelled area [K]
  real, intent(in) :: TCTOP (ILG,IG) !< Thermal conductivity of soil at top of
  !< layer \f$[W m^{-1} K^{-1}] (\lambda_t)\f$
  real, intent(in) :: TCBOT (ILG,IG) !< Thermal conductivity of soil at bottom of
  !< layer \f$[W m^{-1} K^{-1}] (\lambda_b)\f$
  !
  real, intent(in) :: FI    (ILG)    !< Fractional coverage of subarea in question on modelled area [ ]
  real, intent(in) :: ZPOND (ILG)    !< Depth of ponded water on surface [m]
  real, intent(in) :: TBAR1P(ILG)    !< Lumped temperature of ponded water and first soil layer [K]
  real, intent(in) :: TCSNOW(ILG)    !< Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
  real, intent(in) :: ZSNOW (ILG)    !< Depth of snow pack [m]
  !
  integer, intent(in) :: ISAND (ILG,IG) !< Sand content flag
  !
  real, intent(in) :: DELZ  (IG)     !< Overall thickness of soil layer \f$[m] (\Delta_z)\f$
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: DELZ1, A3, B3, C3, TCZERO
  !
  !-----------------------------------------------------------------------
  !     * INITIALIZATION OF ARRAYS.
  !

  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0.) then
      DELZ1 = DELZ(1) + ZPOND(I)
      if (ZPOND(I) > 0.5E-3) then
        IWATER(I) = 1
      else
        if (ISAND(I,1) > - 4) then
          IWATER(I) = 0
        else
          IWATER(I) = 2
        end if
      end if
      !
      if (IWATER(I) == 2) then
        CPHCHG(I) = CLHVAP + CLHMLT
      else
        CPHCHG(I) = CLHVAP
      end if
      !
      if (ZSNOW(I) > 0.0) then
        TCZERO = 1.0 / (0.5 / TCSNOW(I) + 0.5 / TCTOP(I,1))
      else
        TCZERO = TCTOP(I,1)
      end if
      A1(I) = DELZ1 / (3.0 * TCZERO)
      A2(I) = DELZ1 / (2.0 * TCZERO)
      A3 = A2(I)
      B1(I) = DELZ1 / (6.0 * TCBOT(I,1))
      B2(I) = DELZ1 / (2.0 * TCBOT(I,1)) + DELZ(2) / (3.0 * TCTOP(I,2))
      B3 = DELZ1 / (2.0 * TCBOT(I,1)) + DELZ(2) / (2.0 * TCTOP(I,2))
      C2(I) = DELZ(2) / (6.0 * TCBOT(I,2))
      C3 = DELZ(2) / (2.0 * TCBOT(I,2)) + DELZ(3) / (3.0 * TCTOP(I,3))
      GDENOM(I) = A1(I) * (B2(I) * C3 - B3 * C2(I)) - B1(I) * (A2(I) * C3 - &
                  A3 * C2(I))
      GCOEFF(I) = (B2(I) * C3 - B3 * C2(I) - B1(I) * (C3 - C2(I))) / GDENOM(I)
      GCONST(I) = ( - TBAR1P(I) * (B2(I) * C3 - B3 * C2(I)) + &
                  TBAR(I,2) * B1(I) * C3 - &
                  TBAR(I,3) * B1(I) * C2(I)) / GDENOM(I)
    end if
  end do ! loop 100
  !
  return
end subroutine soilHeatFluxPrep
!> \file
!!
!! @author D. Verseghy, M. Lazare
!!
!! In this subroutine, coefficients are derived for an equation
!! relating the heat flux at the ground surface to the ground
!! surface temperature, using the average temperatures and the
!! thermal conductivities of the underlying first three soil layers.
!! It is assumed that the variation of temperature T with depth z
!! within each soil layer can be modelled by using a quadratic
!! equation:
!!
!! \f$T(z) = (1/2) a z^2 + b z + c\f$
!!
!! By substituting 0 for z in the above equation and in the
!! expressions for its first and second derivatives, it can be shown
!! that \f$a = T''(0)\f$, \f$b = T'(0)\f$, and \f$c = T(0)\f$. The term \f$T''(0)\f$ can be
!! evaluated from the expression for the first derivative evaluated
!! at the bottom of the soil layer, \f$T(\Delta z)\f$:
!!
!! \f$T''(0) = [T'(\Delta z) - T'(0)]/ \Delta z\f$
!!
!! The temperature gradient \f$T'(0)\f$ at the top of each layer is
!! related to the heat flux G(0) through the thermal conductivity
!! \f$\lambda_t\f$; and the temperature gradient and heat flux at the bottom
!! of the layer, \f$G(\Delta z)\f$ and \f$T(\Delta z)\f$, are similarly related through
!! the bottom thermal conductivity \f$\lambda_b\f$:
!!
!! \f$G(0) = - \lambda_t T'(0)\f$
!!
!! \f$G(\Delta z) = - \lambda_b T'(\Delta z)\f$
!!
!! The average soil layer temperature, \f$T_{av}(\Delta z)\f$, can be obtained by
!! integrating the resulting equation for T(z) between 0 and \f$\Delta z\f$.
!! Making use of all of the above expressions, recognizing that the
!! heat fluxes and temperatures at the bottoms of layers 1 and 2
!! must equal the heat fluxes and temperatures at the tops of layers
!! 2 and 3 respectively, and neglecting as a first approximation the
!! heat flux at the bottom of the third layer, a linear equation can
!! be derived relating G(0) to T(0) at the soil surface, where the
!! slope and intercept of the equation are functions only of the
!! average temperatures, thicknesses, and top and bottom thermal
!! conductivities of the three soil layers.
!!
!! In the subroutine loop, first the depth corresponding to TBAR1P
!! (the lumped temperature of the first soil layer and the ponded
!! water) is calculated, as the sum of the first soil layer
!! thickness and the ponded water depth. If the ponded water depth
!! is not vanishingly small, the surface water flag IWATER is set to
!! 1; otherwise it is set to 0 for soils and 2 for ice sheets
!! (indicated by ISAND = -4). If IWATER = 2, indicating a frozen
!! water surface, the latent heat of vaporization, CPHCHG, is set to
!! the value for sublimation (by adding the latent heat of melting
!! to the latent heat of vaporization). If there is a snow pack
!! present, the thermal conductivity at the top of the ground
!! surface is calculated as the harmonic mean of the thermal
!! conductivity at the top of the first soil layer and that of the
!! snow pack. Finally, a series of work arrays is evaluated and is
!! used to calculate the slope and intercept, GCOEFF and GCONST, of
!! the equation relating G(0) to T(0) at the ground surface.
!!
