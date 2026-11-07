!> \file
!! Quantifies movement of liquid water between soil layers
!! under conditions of infiltration.
!! @author D. Verseghy, M. Lazare, Y. Delage, R. Soulis, J. Melton
!
subroutine waterFlowInfiltrate (IVEG, THLIQ, THICE, TBARW, BASFLW, TBASFL, & ! Formerly GRINFL
                                RUNOFF, TRUNOF, ZFAV, LZFAV, THLINV, QFG, &
                                WLOST, FI, EVAP, R, TR, TPOND, ZPOND, DT, &
                                ZMAT, WMOVE, TMOVE, THLIQX, THICEX, TBARWX, &
                                DELZX, ZBOTX, FDT, TFDT, PSIF, THLINF, GRKINF, &
                                THLMAX, THTEST, ZRMDR, FDUMMY, TDUMMY, THLDUM, &
                                THIDUM, TDUMW, TRMDR, ZF, FMAX, TUSED, RDUMMY, &
                                ZERO, WEXCES, FDTBND, WADD, TADD, WADJ, TIMPND, &
                                DZF, DTFLOW, THLNLZ, THLQLZ, DZDISP, WDISP, WABS, &
                                THPOR, THLRET, THLMIN, BI, PSISAT, GRKSAT, &
                                THLRAT, THFC, DELZW, ZBOTW, XDRAIN, DELZ, ISAND, &
                                IGRN, IGRD, IFILL, IZERO, LZF, NINF, IFIND, ITER, &
                                NEND, ISIMP, IGDR, &
                                IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !
  !     * JAN 10/17 - J. Melton   Fix below no longer needed. Removed.
  !     * JUL 06/12 - D.VERSEGHY. FIX FOR EVAPORATION OVER ROCK.
  !     * OCT 18/11 - M.LAZARE.   PASS IN "IGDR" AS AN INPUT FIELD
  !     *                         (ORIGINATING IN soilProperties) TO
  !     *                         waterFlowNonInfiltrate AND waterBaseflow.
  !     * APR 04/11 - D.VERSEGHY. MODIFY TEST IN 150 LOOP TO USE DELZ
  !     *                         INSTEAD OF XDRAIN.
  !     * JAN 06/09 - D.VERSEGHY. MODIFY DELZX AND ZBOTX OF BOTTOM LAYER
  !     *                         ADDITIONAL THLIQ CHECK IN 350 LOOP
  !     *                         PASS ADDITIONAL VARIABLES TO waterBaseflow.
  !     * MAR 27/08 - D.VERSEGHY. MOVE VISCOSITY ADJUSTMENT TO waterCalcPrep.
  !     * OCT 31/06 - R.SOULIS.   ADJUST GRKSAT FOR VISCOSITY OF WATER
  !     *                         AND PRESENCE OF ICE; ADJUST THPOR FOR
  !     *                         PRESENCE OF ICE.
  !     * MAR 22/06 - D.VERSEGHY. UNCONDITIONALLY DEFINE VARIABLES FOR
  !     *                         ALL "IF" STATEMENTS.
  !     * SEP 28/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 400 LOOP.
  !     * MAR 23/05 - D.VERSEGHY.R.SOULIS. PASS ADDITIONAL VARIABLES
  !     *                         TO waterFlowNonInfiltrate AND waterBaseflow; PASS OUT ZFAV,
  !     *                         LZFAV, THLINV; CALCULATE GRKTLD
  !     *                         INTERNALLY; REVISE CALCULATION OF
  !     *                         THLINF.
  !     * MAR 16/04 - D.VERSEGHY. TREAT FROZEN SOIL WATER AS ICE
  !     *                         VOLUME RATHER THAN AS EQUIVALENT
  !     *                         LIQUID WATER VOLUME.
  !     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 27/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE
  !     *                         CALCULATIONS AGAINST ROUNDOFF ERRORS.
  !     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
  !     * DEC 12/01 - D.VERSEGHY. PASS NEW VARIABLE IN FOR CALCULATION
  !     *                         OF BASEFLOW.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * APR 17/96 - D.VERSEGHY. CLASS - VERSION 2.5.
  !     *                         BUG FIX: INITIALIZE FDT AND TFDT
  !     *                         TO ZERO.
  !     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. UPDATE SOIL LAYER TEMPERATURES AND
  !     *                         LIQUID MOISTURE CONTENTS FOR
  !     *                         INFILTRATING CONDITIONS (I.E.
  !     *                         PONDED WATER OR RAINFALL OCCURRING
  !     *                         WITHIN CURRENT TIMESTEP).
  !
  use classicParams, only : DELT

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IVEG  !< Subarea type flag
  integer, intent(in) :: IG, IGP1, IGP2, ILG, IL1, IL2, JL, N
  integer :: I, J
  !
  !     * INPUT/OUTPUT FIELDS.
  !
  real, intent(inout) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
  !< layer \f$[m^3 m^{-3}] (\theta_l)\f$
  real, intent(inout) :: THICE (ILG,IG)   !< Volumetric frozen water content of soil
  !< layer \f$[m^3 m^{-3}] (\theta_i)\f$
  real, intent(inout) :: TBARW (ILG,IG)   !< Temperature of water in soil layer [C]
  !
  real, intent(in) :: BASFLW(ILG)  !< Base flow from bottom of soil column [m]
  real, intent(in) :: TBASFL(ILG)  !< Temperature of base flow from bottom of
  !< soil column [K]
  real, intent(in) :: RUNOFF(ILG)  !< Total runoff from soil column [m]
  real, intent(in) :: TRUNOF(ILG)  !< Temperature of total runoff from soil column [K]
  real, intent(in) :: QFG   (ILG)  !< Evaporation from soil surface (diagnostic) \f$[kg m^{-2} s^{-1}]\f$
  real, intent(in) :: WLOST (ILG)  !< Residual amount of water that cannot be
  !< supplied by surface stores \f$[kg m^{-2}]\f$
  real, intent(out) :: ZFAV  (ILG)  !< Average depth of wetting front over current time step [m]
  real, intent(out) :: THLINV(ILG)  !< Liquid water content behind the wetting front \f$[m^3 m^{-3}]\f$
  !
  integer, intent(out) :: LZFAV (ILG)!< Soil layer index in which the average position of the wetting front occurs
  !
  !     * INPUT FIELDS.
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on modelled area [ ]
  real, intent(in) :: EVAP  (ILG)  !< Evaporation rate from ground surface \f$[m s^{-1}]\f$
  real, intent(in) :: R     (ILG)  !< Rainfall rate at ground surface \f$[m s^{-1}]\f$
  real, intent(in) :: TR    (ILG)  !< Temperature of rainfall [C]
  real, intent(in) :: TPOND (ILG)  !< Temperature of ponded water [C]
  real, intent(in) :: ZPOND (ILG)  !< Depth of ponded water on soil surface [m]
  real, intent(in) :: DT    (ILG)  !< Time period over which water movement takes place [s]
  !
  !     * WORK FIELDS (FOR ALL CALLED ROUTINES AS WELL).
  !
  real, intent(in) :: ZMAT  (ILG,IGP2,IGP1)
  !
  real, intent(inout) :: WMOVE (ILG,IGP2), TMOVE(ILG,IGP2)
  real :: GRKSATF(ILG,IG), THPORF(ILG,IG)
  !
  real, intent(inout) :: THLIQX(ILG,IGP1), THICEX(ILG,IGP1), TBARWX(ILG,IGP1), &
                         DELZX (ILG,IGP1), ZBOTX (ILG,IGP1), FDT   (ILG,IGP1), &
                         TFDT  (ILG,IGP1), PSIF  (ILG,IGP1), THLINF(ILG,IGP1), &
                         GRKINF(ILG,IGP1), THLMAX(ILG,IG),   THTEST(ILG,IG), &
                         ZRMDR (ILG,IGP1), FDUMMY(ILG,IGP1), TDUMMY(ILG,IGP1), &
                         THLDUM(ILG,IG),   THIDUM(ILG,IG),   TDUMW (ILG,IG), &
                         RDUMMY(ILG),      TRMDR(ILG),       ZF(ILG)
  !
  real, intent(in) :: FMAX  (ILG), TUSED (ILG), &
                      ZERO  (ILG), WEXCES(ILG), FDTBND(ILG), &
                      WADD  (ILG), TADD  (ILG), WADJ  (ILG), TIMPND(ILG), &
                      DZF   (ILG), DTFLOW(ILG), THLNLZ(ILG), THLQLZ(ILG), &
                      DZDISP(ILG), WDISP (ILG), WABS  (ILG)
  !
  !     * SOIL INFORMATION ARRAYS.
  !
  real, intent(in) :: THPOR (ILG,IG)   !< Pore volume in soil layer \f$[m^3 m^{-3}] (\theta_p)\f$
  real, intent(in) :: THLRET(ILG,IG)   !< Liquid water retention capacity for organic
  !< soil \f$[m^3 m^{-3}] (\theta_{l,ret})\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !< remaining after freezing or evaporation
  !< \f$[m^3 m^{-3}] (theta_{l,min})\f$
  real, intent(in) :: BI    (ILG,IG)   !< Clapp and Hornberger empirical "b" parameter [ ] (b)
  real, intent(in) :: PSISAT(ILG,IG)   !< Soil moisture suction at saturation \f$[m] (\Psi_{sat})\f$
  real, intent(in) :: GRKSAT(ILG,IG)   !< Hydraulic conductivity of soil at
  !< saturation \f$[m s^{-1}] (K_{sat})\f$
  real, intent(in) :: THLRAT(ILG,IG)   !< Fractional saturation of soil behind the
  !< wetting front \f$[ ] (f_{inf})\f$
  real, intent(in) :: THFC  (ILG,IG)   !< Field capacity \f$[m^3 m^{-3}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable depth of soil layer \f$[m] (\Delta z_{g,w})\f$
  real, intent(in) :: ZBOTW (ILG,IG)   !< Depth to permeable bottom of soil layer [m]
  real, intent(in) :: XDRAIN(ILG)      !< Drainage index for water flow at bottom of soil profile [ ]
  real, intent(in) :: DELZ(IG)         !< Overall thickness of soil layer [m]
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: PSIINF, GRK, PSI
  !
  !     * VARIOUS INTEGER ARRAYS.
  !
  integer, intent(in)    :: ISAND(ILG,IG) !< Sand content flag
  integer, intent(inout) :: IGRN (ILG)    !< Flag to indicate whether calculations in this subroutine are to be done
  integer, intent(in)    :: IGRD (ILG)    !< Flag to indicate whether calculations in subroutine waterFlowNonInfiltrate are to be done
  integer, intent(in)    :: IZERO(ILG), IFIND(ILG), &
                            ITER (ILG), NEND (ILG), ISIMP (ILG)

  integer, intent(inout) :: IFILL(ILG), LZF(ILG), NINF(ILG)
  integer, intent(in)    :: IGDR (ILG)    !< Index of soil layer in which bedrock is encountered
  !
  !-----------------------------------------------------------------------

  !     * DETERMINE POINTS WHICH SATISFY CONDITIONS FOR THESE CALCULATIONS
  !     * AND STORE THEM AS HAVING NON-ZERO VALUES FOR WORK ARRAY "IGRN".
  !
  !>
  !! In loop 50, the flag IGRN is first set to 1 for all grid cells
  !! where the calculations in this subroutine are to be performed.
  !! The necessary conditions are: that the surface being modelled is
  !! not a glacier or ice sheet (ISAND > -4); that the time period DT
  !! is greater than zero; and that either the rainfall rate is
  !! greater than zero or that ponded water exists on the surface. If
  !! any of these conditions is not met, IGRN is set to zero.
  !!
  do I = IL1,IL2 ! loop 75
    if (FI(I) > 0. .and. &
        ISAND(I,1) > - 4 .and. DT(I) > 0. .and. &
        (R(I) > 0. .or. ZPOND(I) > 0.)) then
      IGRN(I) = 1
      RDUMMY(I) = 0.
    else
      IGRN(I) = 0
    end if
  end do ! loop 75
  !
  !     * ADJUST GRKSAT FOR VISCOSITY OF WATER AND PRESENCE OF ICE
  !     * ADJUST THPOR FOR PRESENCE OF ICE.
  !     * INITIALIZATION; DETERMINATION OF SOIL HYDRAULIC CONDUCTIVITIES
  !     * AND SOIL MOISTURE SUCTION ACROSS WETTING FRONT.
  !
  !>
  !! In loop 100, a series of local arrays THLIQX, THICEX, TBARWX,
  !! DELZX and ZBOTX are defined for use in the subsequent
  !! infiltration calculations, with a “y” dimension of IG+1, where IG
  !! is the number of soil layers. The entries from 1 to IG are set
  !! equal to the corresponding values in THLIQ, THICE, TBARW, DELZW
  !! and ZBOTW. An effective saturated hydraulic conductivity GRKSATF
  !! is calculated for each soil layer, applying an empirical
  !! correction for the presence of ice. This ice content factor,
  !! \f$f_{ice}\f$, is calculated from Zhao and Gray (1997) as:
  !!
  !! \f$f_{ice} = [1.0 – min(1.0, \theta_i / \theta_p)]^2\f$
  !!
  !! To further account for the presence of ice, a modified pore
  !! volume THPORF for the soil layer is calculated as the maximum of
  !! the available pore volume \f$\theta_p\f$ – \f$\theta_i\f$ (where \f$\theta_p\f$ is the total
  !! pore volume and \f$\theta_i\f$ is the ice content of the layer), the
  !! actual liquid water content \f$\theta_l\f$, and the minimum residual
  !! liquid water content \f$\theta_{l, min}\f$. (The last two conditions are
  !! required because in the case of a saturated soil undergoing
  !! freezing, since water expands when frozen, the sum of the liquid
  !! and frozen volumetric water contents may be greater than the pore
  !! volume, and thus \f$\theta_l\f$ or \f$\theta_{l, min}\f$ may be greater than
  !! \f$\theta_p\f$ – \f$\theta_i\f$.) Finally, the water content THLINF and the
  !! hydraulic conductivity GRKINF behind the wetting front are
  !! evaluated, following the analysis of Green and Ampt (1911) \cite Green1911-gy and
  !! treating the change in soil moisture due to infiltration as a
  !! downward-propagating square wave. THLINF is calculated as the
  !! maximum of \f$f_{inf} ( \theta_p – \theta_i)\f$, \f$\theta_l\f$, and \f$\theta_{l, min}\f$, where \f$f_{inf}\f$
  !! represents the fractional saturation of the soil behind the
  !! wetting front, corresponding to a hydraulic conductivity of half
  !! the saturated value (this correction is applied in order to
  !! account for the fact that as infiltration occurs, a small amount
  !! of air is usually trapped in the soil). GRKINF is calculated from
  !! GRKSATF, THLINF and THPORF, using the classic Clapp and
  !! Hornberger (1978) \cite Clapp1978-898 equation
  !!
  !! \f$K(z) = K_{sat} ( \theta_l / \theta_p)^{(2b + 3)}\f$
  !!
  !! where K(z) is the hydraulic conductivity at a depth z, \f$K_{sat}\f$ is
  !! the hydraulic conductivity at saturation, \f$\theta_p\f$ is the pore
  !! volume and b is an empirical coefficient.
  !!
  do J = 1,IG ! loop 100
    do I = IL1,IL2
      if (IGRN(I) > 0) then
        THLIQX(I,J) = THLIQ(I,J)
        THICEX(I,J) = THICE(I,J)
        TBARWX(I,J) = TBARW(I,J)
        DELZX(I,J) = DELZW(I,J)
        ZBOTX(I,J) = ZBOTW(I,J)
        FDT (I,J) = 0.0
        TFDT(I,J) = 0.0
        if (ISAND(I,J) > - 3) then
          GRKSATF(I,J) = GRKSAT(I,J) * (1.0 - MAX(0.0,MIN(1.0, &
                         THICE(I,J) / THPOR(I,J)))) ** 2
          THPORF(I,J) = MAX((THPOR(I,J) - THICE(I,J) - 0.00001), &
                        THLIQ(I,J),THLMIN(I,J))
          THLINF(I,J) = MAX(THLIQ(I,J),THLMIN(I,J), &
                        THLRAT(I,J) * (THPOR(I,J) - &
                        THICE(I,J) - 0.00001))
          GRKINF(I,J) = GRKSATF(I,J) * (THLINF(I,J) / THPORF(I,J)) &
                        ** (2. * BI(I,J) + 3.)
        else
          GRKSATF(I,J) = 0.0
          THPORF(I,J) = 0.0
          THLINF(I,J) = 0.0
          GRKINF(I,J) = 0.0
        end if
      end if
    end do
  end do ! loop 100
  !
  !>
  !! In loop 150, values for the IG+1 entries in each of the above
  !! arrays are assigned. If the permeable thickness DELZW of the IG
  !! layer is less than its overall thickness DELZ (indicating the
  !! presence of bedrock) the IG+1 entries are set to 0; otherwise
  !! they are set to the values for the IG layer in the case of
  !! THLIQX, THICEX, TBARWX and THLINF, and for DELZX and ZBOTX they
  !! are set to large positive numbers. GRKINF is set to the value for
  !! the IG layer multiplied by the drainage parameter XDRAIN, which
  !! takes a value of 0 if the soil is underlain by an impermeable
  !! layer (as in a bog), and a value of 1 otherwise.
  !!
  do I = IL1,IL2 ! loop 150
    if (IGRN(I) > 0) then
      if (DELZW(I,IG) < DELZ(IG)) then
        THLIQX(I,IG + 1) = 0.0
        THICEX(I,IG + 1) = 0.0
        TBARWX(I,IG + 1) = 0.0
        DELZX(I,IG + 1) = 0.0
        THLINF(I,IG + 1) = 0.0
        GRKINF(I,IG + 1) = 0.0
      else
        THLIQX(I,IG + 1) = THLIQX(I,IG)
        THICEX(I,IG + 1) = THICEX(I,IG)
        TBARWX(I,IG + 1) = TBARWX(I,IG)
        DELZX(I,IG + 1) = 999999.
        THLINF(I,IG + 1) = THLINF(I,IG)
        GRKINF(I,IG + 1) = GRKINF(I,IG) * XDRAIN(I)
      end if
      ZBOTX (I,IG + 1) = ZBOTX(I,IG) + DELZX(I,IG + 1)
      FDT   (I,IG + 1) = 0.0
      TFDT  (I,IG + 1) = 0.0
    end if
  end do ! loop 150
  !
  !>
  !! In loop 200, the soil water suction across the wetting front
  !! psi_f is calculated for each soil layer, using equation 25 in
  !! Verseghy (1991):
  !!
  !! \f$\Psi_f = b[\Psi_{inf} K_{inf} - \Psi(z) K(z)]/[K_{inf} (b+3)]\f$
  !!
  !! where \f$\Psi_{inf}\f$ and \f$K_{inf}\f$ are the soil moisture suction and the
  !! hydraulic conductivity behind the wetting front, and \f$\Psi(z)\f$ and
  !! \f$K(z)\f$ are the soil moisture suction and the hydraulic conductivity
  !! ahead of the wetting front. The soil moisture suction values are
  !! obtained using the classic Clapp and Hornberger (1978) equation:
  !!
  !! \f$\Psi(z) = \Psi_{sat} (\theta_l / \theta_p)^{(-b)}\f$
  !!
  !! where \f$\Psi_{sat}\f$ is the soil moisture suction at saturation.
  !!
  do J = 1,IG ! loop 200
    do I = IL1,IL2
      if (IGRN(I) > 0) then
        if (THPOR(I,J) > 0.0001) then
          PSIINF = MAX(PSISAT(I,J) * (THLINF(I,J) / THPORF(I,J)) ** &
                   ( - BI(I,J)),PSISAT(I,J))
          GRK = MIN(GRKSATF(I,J) * (THLIQ(I,J) / THPORF(I,J)) ** &
                (2. * BI(I,J) + 3.),GRKSATF(I,J))
          PSI = MAX(PSISAT(I,J) * (THLIQ(I,J) / THPORF(I,J)) ** &
                ( - BI(I,J)),PSISAT(I,J))
        else
          PSIINF = PSISAT(I,J)
          GRK = GRKSATF(I,J)
          PSI = PSISAT(I,J)
        end if
        if (THLINF(I,J) > THLIQ(I,J)) then
          PSIF (I,J) = MAX(BI(I,J) * (GRKINF(I,J) * PSIINF - GRK * PSI) / &
                       (GRKINF(I,J) * (BI(I,J) + 3.)),0.0)
        else
          PSIF (I,J) = 0.0
        end if
      end if
    end do
  end do ! loop 200

  do I = IL1,IL2 ! loop 250
    if (IGRN(I) > 0) then
      PSIF (I,IG + 1) = PSIF (I,IG)
      TRMDR(I) = DELT
    else
      TRMDR(I) = 0.
    end if
  end do ! loop 250

  do J = 1,IGP2 ! loop 300
    do I = IL1,IL2
      if (IGRN(I) > 0) then
        WMOVE(I,J) = 0.0
        TMOVE(I,J) = 0.0
      end if
    end do
  end do ! loop 300
  !
  !     * DETERMINE STARTING POSITION OF WETTING FRONT; INITIALIZATION
  !     * FOR SATURATED INFILTRATION.
  !
  !>
  !! In loop 400, a test is carried out to determine whether saturated
  !! conditions are already present in the soil. Generally it is
  !! assumed that a period of unsaturated flow occurs initially, so
  !! the flag IFILL is set to 1, the depth of the wetting front ZF is
  !! set to zero, and the flag LZF indicating the index of the soil
  !! layer in which the wetting front occurs is set to 1. If a pond is
  !! present on the soil surface or if the hydraulic conductivity if
  !! the first soil layer is very small, it is deemed that saturated
  !! flow begins immediately, so IFILL is set to zero. Then each soil
  !! layer is checked in turn, to ascertain whether the liquid water
  !! content is greater than THLINF. If so, it is concluded that
  !! saturated flow is occurring in this layer; ZF is set to ZBOTW,
  !! the depth of the bottom of the soil layer, LZF is set to the
  !! index of the next layer, and IFILL is set to zero. For each layer
  !! found to be undergoing saturated flow, its water content is
  !! stored in the J+1 level of the water movement matrix WMOVE, and
  !! its water temperature is stored in the J+1 level of the matrix
  !! TMOVE. The counter NINF is set to the number of soil layers
  !! behind and including the one with the wetting front, plus 2.
  !!
  do I = IL1,IL2 ! loop 400
    if (IGRN(I) > 0) then
      IFILL(I) = 1
      ZF(I) = 0.0
      LZF(I) = 1
      if (ZPOND(I) > 0. .or. GRKINF(I,1) < 1.0E-12) then
        NINF(I) = 2
        TMOVE(I,2) = TBARWX(I,1)
        IFILL(I) = 0
      end if
      do J = 1,IG ! loop 350
        if (THLIQ(I,J) >= (THLINF(I,J) - 1.0E-6) .and. &
            THLIQ(I,J) > 0.0001 .and. LZF(I) == J) then
          ZF(I) = ZBOTW(I,J)
          LZF(I) = J + 1
          NINF(I) = J + 2
          WMOVE(I,J + 1) = THLIQ(I,J) * DELZW(I,J)
          TMOVE(I,J + 1) = TBARWX(I,J)
          TMOVE(I,J + 2) = TBARWX(I,J + 1)
          IFILL(I) = 0
        end if
      end do ! loop 350
    else
      IFILL(I) = 0
      LZF(I) = 0
      NINF(I) = 0
    end if

  end do ! loop 400
  !>
  !! A series of subroutines is called to complete the infiltration
  !! calculations. Subroutine waterInfiltrateUnsat treats the period of unsaturated
  !! flow up the point when saturated flow begins. Subroutine waterInfiltrateSat
  !! treats the period of saturated flow. Subroutine waterBaseflow reassigns
  !! liquid water contents and temperatures at the conclusion of the
  !! infiltration period. THLIQ, THICE and TBARW are then set to the
  !! corresponding updated values of THLIQX, THICEX and TBARWX. The
  !! average depth of the wetting front in the soil layer in which it
  !! occurs, the index of the soil layer, and the moisture content
  !! behind the wetting front are stored in output variables ZFAV,
  !! LZFAV and THLINV respectively. Finally, if there is any time
  !! remaining in the current time step following the infiltration
  !! period, subroutine waterFlowNonInfiltrate is called to treat the movement of soil
  !! water over the remainder of the time step.
  !!
  !
  !     * IF SATURATED INFILTRATION CONDITIONS ARE NOT PRESENT AT ONCE
  !     * (IFILL=1), CALL "waterInfiltrateUnsat" TO DO PROCESSING FOR PERIOD OF
  !     * UNSATURATED INFILTRATION.
  !
  call waterInfiltrateUnsat(WMOVE, TMOVE, LZF, NINF, ZF, TRMDR, R, TR, & ! Formerly WFILL
                            PSIF, GRKINF, THLINF, THLIQX, TBARWX, &
                            DELZX, ZBOTX, DZF, TIMPND, WADJ, WADD, &
                            IFILL, IFIND, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !
  !     * CALL "waterInfiltrateSat" TO DO PROCESSING FOR PERIOD OF SATURATED
  !     * INFILTRATION.
  !
  call waterInfiltrateSat(WMOVE, TMOVE, LZF, NINF, TRMDR, TPOND, ZPOND, & ! Formerly WFLOW
                          R, TR, EVAP, PSIF, GRKINF, THLINF, THLIQX, TBARWX, &
                          DELZX, ZBOTX, FMAX, ZF, DZF, DTFLOW, THLNLZ, &
                          THLQLZ, DZDISP, WDISP, WABS, ITER, NEND, ISIMP, &
                          IGRN, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !
  !     * RECALCULATE TEMPERATURES AND LIQUID MOISTURE CONTENTS OF
  !     * SOIL LAYERS FOLLOWING INFILTRATION.
  !
  call waterBaseflow(THLIQX, THICEX, TBARWX, ZPOND, TPOND, & ! Formerly WEND
                     BASFLW, TBASFL, RUNOFF, TRUNOF, FI, &
                     WMOVE, TMOVE, LZF, NINF, TRMDR, THLINF, DELZX, &
                     ZMAT, ZRMDR, FDTBND, WADD, TADD, FDT, TFDT, &
                     THLMAX, THTEST, THLDUM, THIDUM, TDUMW, &
                     TUSED, RDUMMY, ZERO, WEXCES, XDRAIN, &
                     THPOR, THLRET, THLMIN, BI, PSISAT, GRKSAT, &
                     THFC, DELZW, ISAND, IGRN, IGRD, IGDR, IZERO, &
                     IVEG, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !
  do J = 1,IG ! loop 800
    do I = IL1,IL2
      if (IGRN(I) > 0) then
        THLIQ(I,J) = THLIQX(I,J)
        THICE(I,J) = THICEX(I,J)
        TBARW(I,J) = TBARWX(I,J)
      end if
    end do
  end do ! loop 800
  !
  do I = IL1,IL2 ! loop 850
    if (IGRN(I) > 0 .and. LZF(I) < IG + 1) then
      ZFAV(I) = (ZF(I) + MAX(ZBOTW(I,LZF(I)) - DELZW(I,LZF(I)),0.0)) / &
                2.0
      LZFAV(I) = LZF(I)
      THLINV(I) = THLINF(I,LZF(I))
    else
      ZFAV(I) = 0.0
      LZFAV(I) = 0
      THLINV(I) = 0.0
    end if
  end do ! loop 850
  !
  !     * IF TIME REMAINS IN THE CURRENT MODEL STEP AFTER INFILTRATION
  !     * HAS CEASED (TRMDR>0), CALL "waterFlowNonInfiltrate" TO CALCULATE WATER FLOWS
  !     * BETWEEN LAYERS FOR THE REMAINDER OF THE TIME STEP.
  !
  call waterFlowNonInfiltrate(IVEG, THLIQ, THICE, TBARW, FDUMMY, TDUMMY, BASFLW, & ! Formerly GRDRAN
                              TBASFL, RUNOFF, TRUNOF, QFG, WLOST, FI, EVAP, ZERO, ZERO, &
                              TRMDR, WEXCES, THLMAX, THTEST, THPOR, THLRET, THLMIN, &
                              BI, PSISAT, GRKSAT, THFC, DELZW, XDRAIN, ISAND, IZERO, &
                              IZERO, IGRD, IGDR, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !
  return
end subroutine waterFlowInfiltrate
