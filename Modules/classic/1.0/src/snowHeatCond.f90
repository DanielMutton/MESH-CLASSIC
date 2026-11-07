!> \file
!! Calculate coefficients for solution of snow pack heat conduction.
!
subroutine snowHeatCond (GCOEFFS, GCONSTS, CPHCHG, IWATER, & ! Formerly TSPREP
                         FI, ZSNOW, TSNOW, TCSNOW, &
                         ILG, IL1, IL2, JL)
  !
  !     * AUG 16/06 - D.VERSEGHY. MAJOR REVISION TO IMPLEMENT THERMAL
  !     *                         SEPARATION OF SNOW AND SOIL.
  !     * MAY 24/06 - D.VERSEGHY. LIMIT DELZ3 TO <= 4.1 M.
  !     * OCT 04/05 - D.VERSEGHY. USE THREE-LAYER TBAR, TCTOP, TCBOT.
  !     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 06/02 - D.VERSEGHY. SHORTENED CLASS3 COMMON BLOCK.
  !     * JUN 17/02 - D.VERSEGHY. USE NEW LUMPED SOIL AND PONDED WATER
  !     *                         TEMPERATURE FOR FIRST LAYER; SHORTENED
  !     *                         CLASS4 COMMON BLOCK.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         INCORPORATE EXPLICITLY CALCULATED
  !     *                         THERMAL CONDUCTIVITIES AT TOPS AND
  !     *                         BOTTOMS OF SOIL LAYERS.
  !     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
  !     *                         ORGANIC MATTER CONTENT.
  !     * NOV 28/94 - M. LAZARE.  CLASS - VERSION 2.3.
  !     *                         TCSATW, TCSATI DECLARED REAL(16).
  !     * APR 10/92 - M. LAZARE.  CLASS - VERSION 2.1.
  !     *                         DIVIDE PREVIOUS SUBROUTINE "T4LAYR"
  !     *                         INTO "snowHeatCond" AND "snowTempUpdate" AND
  !     *                         VECTORIZE.
  !     * APR 11/89 - D.VERSEGHY. CALCULATE COEFFICIENTS FOR GROUND HEAT
  !     *                         FLUX, EXPRESSED AS A LINEAR FUNCTION OF
  !     *                         SURFACE TEMPERATURE. COEFFICIENTS ARE
  !     *                         CALCULATED FROM LAYER TEMPERATURES,
  !     *                         THICKNESSES AND THERMAL CONDUCTIVITIES,
  !     *                         ASSUMING A QUADRATIC VARIATION OF
  !     *                         TEMPERATURE WITH DEPTH WITHIN EACH
  !     *                         SOIL/SNOW LAYER. SET THE SURFACE
  !     *                         LATENT HEAT OF VAPORIZATION OF WATER
  !     *                         AND THE STARTING TEMPERATURE FOR THE
  !     *                         ITERATION IN "energBalVegSolve"/"energBalNoVegSolve".
  !
  use classicParams, only : CLHMLT, CLHVAP

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IL1, IL2, JL
  integer :: I, J
  !
  !     * OUTPUT ARRAYS.
  !
  real, intent(out) :: GCOEFFS(ILG)     !< Multiplier used in equation relating snow surface heat flux to snow surface
  !< temperature \f$[W m^{-2} K^{-1}]\f$
  real, intent(out) :: GCONSTS(ILG)     !< Intercept used in equation relating snow surface heat flux to snow surface
  !< temperature \f$[W m^{-2} ]\f$
  real, intent(out) :: CPHCHG(ILG)      !< Latent heat of sublimation \f$[J kg^{-1}]\f$

  !
  integer, intent(out) :: IWATER(ILG)   !< Flag indicating condition of surface (dry, water-covered or snow-covered)
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on modelled area [ ]
  real, intent(in) :: ZSNOW (ILG)  !< Depth of snow pack \f$[m] (\Delta z_s)\f$
  real, intent(in) :: TSNOW (ILG)  !< Snowpack temperature \f$[K] (T_s(\Delta z_s))\f$
  real, intent(in) :: TCSNOW(ILG)  !< Thermal conductivity of snow \f$[W m^{-1} K^{-1}]\f$
  !
  !-----------------------------------------------------------------------
  !
  !     * CALCULATE COEFFICIENTS.
  !
  do I = IL1,IL2 ! loop 100
    if (FI(I) > 0.) then
      GCOEFFS(I) = 3.0 * TCSNOW(I) / ZSNOW(I)
      GCONSTS(I) = - 3.0 * TCSNOW(I) * TSNOW(I) / ZSNOW(I)
      IWATER(I) = 2
      CPHCHG(I) = CLHVAP + CLHMLT
    end if
  end do ! loop 100
  !
  return
end subroutine snowHeatCond
!> \file
!!
!! @author D. Verseghy, M. Lazare
!!
!!   In this subroutine, coefficients are derived for an equation
!!   relating the heat flux at the snow surface to the snow surface
!!   temperature. It is assumed that the variation of temperature T with
!!   depth z in the snow pack can be modelled by using a quadratic
!!   equation:
!!
!!   \f$T(z) = (1/2) a z^2 + b z + c\f$
!!
!!   By substituting 0 for z in the above equation and in the expressions
!!   for its first and second derivatives, it can be shown that a =
!!   T"(0), b = T'(0), and c = T(0). The term T"(0) can be evaluated from
!!   the expression for the first derivative evaluated at the bottom of
!!   the snow pack, \f$T(\Delta z_s)\f$:
!!
!!   \f$T"(0) = [T'(\Delta z_s) - T'(0)]/ \Delta z_s\f$
!!
!!   The temperature gradient T'(0) at the snow surface is related to the
!!   surface heat flux G(0) by the snow
!!   thermal conductivity \f$\lambda_s\f$:
!!
!!   \f$G(0) = -\lambda_s T'(0)\f$
!!
!!   The average snow temperature, \f$T_s(\Delta z_s)\f$, can be obtained by integrating
!!   the resulting equation for T(z) between 0 and \f$\Delta z_s\f$. Making use of
!!   all of the above expressions, and assuming as a first approximation
!!   that the heat flux at the bottom of the snow pack is zero, a linear
!!   equation can be derived relating G(0) to T(0):
!!
!!   \f$G(0) = 3 \lambda_s / \Delta z_s [T(0) - T_s(\Delta z_s)]\f$
!!
!!   Just four calculations are performed in this subroutine. The slope
!!   and intercept of the G(0) vs. T(0) relation, GCOEFF and GCONST, are
!!   evaluated as \f$3 \lambda_s / \Delta z_s\f$ and \f$-3 \lambda_s T_s(\Delta z_s)/ \Delta z_s\f$
!!   respectively; the flag IWATER is set to 2, indicating a snow
!!   surface; and the latent heat of vaporization at the surface,
!!   CPHCHG, is set to the value for sublimation (by adding the latent
!!   heat of melting to the latent heat of vaporization).
