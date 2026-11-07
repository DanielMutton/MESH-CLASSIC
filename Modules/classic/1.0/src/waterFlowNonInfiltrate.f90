!> \file
!! Quantifies movement of liquid water between soil layers
!! under non-infiltrating conditions, in response to gravity and
!! tension forces.
!! @author D. Verseghy, M. Lazare, P. Bartlett, R. Soulis, F. Seglenieks, V. Fortin, Y. Delage, L. Spacek
!
subroutine waterFlowNonInfiltrate (IVEG, THLIQ, THICE, TBARW, FDT, TFDT, BASFLW, TBASFL, & ! Formerly GRDRAN
                                   RUNOFF, TRUNOF, QFG, WLOST, FI, EVAP, R, ZPOND, DT, &
                                   WEXCES, THLMAX, THTEST, THPOR, THLRET, THLMIN, &
                                   BI, PSISAT, GRKSAT, THFC, DELZW, XDRAIN, ISAND, LZF, &
                                   IGRN, IGRD, IGDR, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !

  !     * OCT 18/11 - M.LAZARE.   PASS IN "IGDR" AS AN INPUT FIELD
  !     *                         (ORIGINATING IN soilProperties) RATHER
  !     *                         THAN REPEATING THE CALCULATION HERE
  !     *                         AS AN INTERNAL WORK FIELD.
  !     * SEP 11/11 - D.VERSEGHY. CHANGE IF CONDITION ON "J<IG"
  !     *                         TO "J<IGDR(I)" IN LOOPS 400
  !     *                         AND 500, TO BE CONSISTENT WITH
  !     *                         OTHERS.
  !     * DEC 10/10 - D.VERSEGHY. ALLOW DRAINAGE AT BEDROCK SURFACE
  !     *                         ANYWHERE IN SOIL PROFILE.
  !     * DEC 23/09 - V.FORTIN.   NEW CALCULATION OF BASEFLOW.
  !     * MAR 31/09 - D.VERSEGHY. PASS IN LZF, AND ZERO OUT FLOWS AT
  !     *                         TOP AND BOTTOM OF SOIL LAYERS DOWN
  !     *                         TO LAYER CONTAINING WETTING FRONT
  !     *                         IN CASES WHERE INFILTRATION IS
  !     *                         OCCURRING.
  !     * JAN 06/09 - D.VERSEGHY. MODIFIED CALCULATION OF GRKSATF
  !     *                         ADJUSTMENTS TO WATER FLUX
  !     *                         CORRECTIONS IN 500 LOOP.
  !     * MAR 27/08 - D.VERSEGHY. MOVE VISCOSITY ADJUSTMENT TO waterCalcPrep.
  !     * OCT 31/06 - R.SOULIS.   ADJUST GRKSAT FOR VISCOSITY OF
  !     *                         WATER AND PRESENCE OF ICE; ADJUST
  !     *                         THPOR FOR PRESENCE OF ICE.
  !     * JUN 06/06 - F.SEGLENIEKS. CHANGE CALCULATION OF GRSBND
  !     *                           TO USE HARMONIC MEAN INSTEAD OF
  !     *                           GEOMETRIC MEAN.
  !     * MAY 17/06 - D.VERSEGHY. MODIFY CALCULATION OF THLMAX TO
  !     *                         ALLOW FOR OVERSATURATED CONDITIONS.
  !     * MAR 21/06 - D.VERSEGHY. PROTECT CALCULATIONS OF TBASFL AND
  !     *                         TRUNOF AGAINST DIVISION BY ZERO.
  !     * MAR 23/05 - R.SOULIS/D.VERSEGHY. CALCULATE GRKSAT AND
  !     *                         PSISAT AT LAYER BOUNDARIES USING
  !     *                         GEOMETRIC MEAN; SET BASEFLOW TO
  !     *                         GRKSAT IF THLIQ > THFC; ADD
  !     *                         CALCULATION OF RUNOFF TEMPERATURE.
  !     * MAR 16/05 - D.VERSEGHY. TREAT FROZEN SOIL WATER AS ICE
  !     *                         VOLUME RATHER THAN AS EQUIVALENT
  !     *                         LIQUID WATER VOLUME.
  !     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 27/04 - D.VERSEGHY, Y.DELAGE. PROTECT SENSITIVE
  !     *                         CALCULATIONS AGAINST ROUNDOFF ERRORS.
  !     * DEC 03/03 - D.VERSEGHY. IMPROVE HANDLING OF CAPILLARY RISE
  !     *                         VS. GRAVITY DRAINAGE (ESPECIALLY
  !     *                         FOR ORGANIC SOILS).
  !     * JUL 31/03 - D.VERSEGHY. ALWAYS CALCULATE THLMAX IN 100 LOOP.
  !     * OCT 23/02 - D.VERSEGHY. REFINEMENT OF TEST IN 400 LOOP.
  !     * JUN 21/02 - D.VERSEGHY. BUGFIX IN CALCULATION OF FDT'S IN
  !     *                         400 LOOP; UPDATE SUBROUTINE CALL
  !     *                         SHORTENED CLASS4 COMMON BLOCK.
  !     * MAY 21/02 - D.VERSEGHY. STREAMLINE CALCULATIONS FOR ORGANIC
  !     *                         SOILS AND MODIFY CHECK ON EVAPORATION
  !     *                         RATE.
  !     * DEC 12/01 - D.VERSEGHY. ADD SEPARATE CALCULATION OF BASEFLOW
  !     *                         AT BOTTOM OF SOIL COLUMN.
  !     * SEP 28/00 - P.BARTLETT/D.VERSEGHY. BUG FIX IN CALCULATION
  !     *                                    OF PSI IN LOOP 200.
  !     * FEB 08/00 - D.VERSEGHY/L.SPACEK. MINOR BUG FIX IN LOOP 600
  !     *                                  RE. ADDRESSING OF THLIQ.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * DEC 30/96 - D.VERSEGHY. CLASS - VERSION 2.6.
  !     *                         BUGFIX IN CALCULATION OF QFG.
  !     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         ADDITIONAL DIAGNOSTIC CALCULATIONS.
  !     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. UPDATE SOIL LAYER TEMPERATURES AND
  !     *                         LIQUID MOISTURE CONTENTS FOR
  !     *                         NON-INFILTRATING CONDITIONS (I.E.
  !     *                         NO PONDED WATER AND NO RAINFALL
  !     *                         OCCURRING WITHIN CURRENT TIMESTEP).
  !
  use classicParams, only : DELT, TFREZ, HCPW, RHOW, RHOICE, CLHMLT

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IVEG  !< Subarea type flag
  integer, intent(in) :: IG, IGP1, IGP2, ILG, IL1, IL2, JL, N
  integer :: I, J, K, IPTBAD
  !
  !     * INPUT/OUTPUT FIELDS.
  !
  real, intent(inout) :: THLIQ (ILG,IG)   !< Volumetric liquid water content of soil
  !< layer \f$[m^3 m^{-3}] (\theta_l)\f$
  real, intent(inout) :: THICE (ILG,IG)   !< Volumetric frozen water content of soil
  !< layer \f$[m^3 m^{-3}] (\theta_i)\f$
  real, intent(inout) :: TBARW (ILG,IG)   !< Temperature of water in soil layer [C]
  real, intent(inout) :: FDT  (ILG,IGP1)  !< Water flow at soil layer interfaces during
  !< current time step [m]
  real, intent(inout) :: TFDT  (ILG,IGP1) !< Temperature of water flowing between soil layers [C]
  !
  real, intent(inout) :: BASFLW(ILG)  !< Base flow from bottom of soil column [m]
  real, intent(inout) :: TBASFL (ILG) !< Temperature of base flow from bottom of soil column [K]
  real, intent(inout) :: RUNOFF(ILG)  !< Total runoff from soil column [m]
  real, intent(inout) :: TRUNOF (ILG) !< Temperature of total runoff from soil column [K]
  real, intent(inout) :: QFG    (ILG) !< Evaporation from soil surface (diagnostic)
  !< \f$[kg m^{-2} s^{-1}]\f$
  real, intent(inout) :: WLOST (ILG)  !< Residual amount of water that cannot be
  !< supplied by surface stores \f$[kg m^{-2}]\f$
  !
  !     * INPUT FIELDS.
  !
  real, intent(in) :: FI    (ILG)  !< Fractional coverage of subarea in question on
  !< modelled area [ ] \f$(X_i)\f$
  real, intent(in) :: EVAP  (ILG)  !< Evaporation rate from ground surface \f$[m s^{-1}]\f$
  real, intent(in) :: R     (ILG)  !< Rainfall rate at ground surface \f$[m s^{-1}]\f$
  real, intent(in) :: ZPOND (ILG)  !< Depth of ponded water on soil surface [m]
  real, intent(in) :: DT    (ILG)  !< Time period over which water movement takes place [s]
  !
  integer, intent(in) :: IGRN  (ILG)  !< Flag to indicate whether calculations in subroutine waterFlowInfiltrate are done
  integer, intent(in) :: LZF   (ILG)  !< Index of soil layer in which wetting front is located
  integer, intent(inout) :: IGRD  (ILG)  !< Flag to indicate whether calculations in this subroutine are to be done
  integer, intent(in) :: IGDR  (ILG)  !< Index of soil layer in which bedrock is encountered
  !
  !     * WORK FIELDS.
  !
  real, intent(inout) :: WEXCES(ILG), THLMAX(ILG,IG), THTEST(ILG,IG)
  real                :: GRKSATF(ILG,IG), THPORF(ILG,IG)
  !
  !     * SOIL INFORMATION ARRAYS.
  !
  real, intent(in) :: THPOR (ILG,IG)   !< Pore volume in soil layer \f$[m^3 m^{-3}] (\theta_p)\f$
  real, intent(in) :: THLRET(ILG,IG)   !< Liquid water retention capacity for organic
  !< soil \f$[m^3 m^{-3} ] (\theta_{l,ret})\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content
  !< remaining after freezing or evaporation
  !< \f$[m^3 m^{-3}] (\theta_{l,min})\f$
  real, intent(in) :: BI    (ILG,IG)   !< Clapp and Hornberger empirical "b" parameter [ ] (b)
  real, intent(in) :: PSISAT(ILG,IG)   !< Soil moisture suction at saturation \f$[m] (\Psi_{sat})\f$
  real, intent(in) :: GRKSAT(ILG,IG)   !< Hydraulic conductivity of soil at
  !< saturation \f$[m s^{-1}] (K_{sat})\f$
  real, intent(in) :: THFC  (ILG,IG)   !< Field capacity \f$[m^3 m^{-3}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable depth of soil layer \f$[m] (\Delta_{zg,w})\f$
  real, intent(in) :: XDRAIN(ILG)      !< Drainage index for water flow at bottom of soil profile [ ]
  !
  integer, intent(in) :: ISAND (ILG,IG)!< Sand content flag
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: THPBND, THLBND, DTHLDZ, DTHPDZ, BBND, GRSBND, PSSBND, GRK, PSI, &
          WLIMIT, THSUBL, THLTHR, CCH, ASAT, ASATC, SATB
  !
  !-----------------------------------------------------------------------

  !     * DETERMINE POINTS WHICH SATISFY CONDITIONS FOR THESE CALCULATIONS
  !     * AND STORE THEM AS HAVING NON-ZERO VALUES FOR WORK ARRAY "IGRD".
  !     * NOTE THAT POINTS WHICH GO THROUGH THE ROUTINE "waterFlowInfiltrate" SHOULD
  !     * NOT GO THROUGH THIS ROUTINE WHEN IT IS CALLED FROM waterBudgetDriver.
  !     * THE INPUT ARRAY "IGRN" HANDLES THIS CONDITION (PASSED AS
  !     * "IZERO" ARRAY WHEN CALLED FROM "waterBaseflow" OR THE END OF "waterFlowInfiltrate").
  !
  !>
  !! In loop 50, the flag IGRD is first set to 1 for all grid cells
  !! where the calculations in this subroutine are to be performed.
  !! The necessary conditions are: that the surface being modelled is
  !! not a glacier or ice sheet (ISAND > -4); that the time period DT
  !! is greater than zero; that the infiltration calculations in
  !! subroutine waterFlowInfiltrate are not simultaneously being performed
  !! (IGRN = 0); and that the rainfall rate and the depth of ponded
  !! water are both vanishingly small. If any of these conditions is
  !! not met, IGRD is set to zero.
  !!
  do I = IL1,IL2
    if (FI (I) > 0. .and. &
        ISAND(I,1) > - 4 .and. DT(I) > 0. .and. IGRN(I) == 0 .and. &
        (R(I) < 1.0E-12 .and. ZPOND(I) < 1.0E-12)) then
      IGRD(I) = 1
    else
      IGRD(I) = 0
    end if
  end do ! loop 50
  !
  !     * CALCULATE MAXIMUM LIQUID WATER CONTENT OF EACH SOIL LAYER
  !     * ADJUST GRKSAT FOR VISCOSITY OF WATER AND PRESENCE OF ICE
  !     * ADJUST THPOR FOR PRESENCE OF ICE.
  !
  !>
  !! In loop 100, if the surface being modelled is not an ice sheet
  !! (ISAND = -4) and if the soil layer in question is not completely
  !! rock (ISAND = -3), the maximum possible water content of each
  !! soil layer, THLMAX, is calculated as the maximum of the available
  !! pore volume THPOR – THICE (where THPOR is the total pore volume
  !! and THICE is the ice content of the layer), the actual liquid
  !! water content THLIQ, and the minimum residual liquid water
  !! content THLMIN. The last two conditions are required because in
  !! the case of a saturated soil undergoing freezing, since water
  !! expands when frozen, the sum of the liquid and frozen volumetric
  !! water contents may be greater than the pore volume, and thus
  !! THLIQ or THLMIN may be greater than THPOR – THICE. An effective
  !! saturated hydraulic conductivity GRKSATF of the soil layer is
  !! also defined, applying an empirical correction for the presence
  !! of ice. This ice content factor, fice, is calculated from Zhao
  !! and Gray (1997) as:
  !!
  !! \f$f_{ice} = [1.0 – min((\theta_p - \theta_{l, min} )/\theta_p , \theta_i / \theta_p)]^2\f$
  !!
  !! The pore volume of the soil layer is corrected for the presence
  !! of ice by defining the effective porevolume THPORF as equivalent
  !! to THLMAX.
  !!
  do J = 1,IG ! loop 100
    do I = IL1,IL2
      if (IGRD(I) > 0) then
        if (ISAND(I,J) > - 3) then
          THLMAX(I,J) = MAX((THPOR(I,J) - THICE(I,J) - 0.00001), &
                        THLIQ(I,J),THLMIN(I,J))
          GRKSATF(I,J) = GRKSAT(I,J) * (1.0 - MAX(0.0,MIN((THPOR(I,J) - &
                         THLMIN(I,J)) / THPOR(I,J),THICE(I,J) / THPOR(I,J)))) ** 2
          THPORF(I,J) = THLMAX(I,J)
        else
          THLMAX(I,J) = 0.0
          GRKSATF(I,J) = 0.0
          THPORF(I,J) = 0.0
        end if
      end if
    end do
  end do ! loop 100
  !
  !     * CALCULATE THEORETICAL FLOW RATES AT BOTTOM OF PERMEABLE SOIL
  !     * DEPTH AND BETWEEN SOIL LAYERS.
  !
  !>
  !! In loops 150 and 200, the theoretical amounts of water FDT
  !! flowing across the soil layer boundaries are calculated. FDT is
  !! evaluated as the product of the flow rate F(z) at the given depth
  !! z, and the period of time (DT) over which the flow is occurring.
  !! At z=0, the water flux is simply equal to the soil surface
  !! evaporation rate. Within the soil, the flow rate F at a given
  !! depth z is obtained using equation 21 from Verseghy (1991):
  !!
  !! \f$F(z) = K(z) [-b \Psi(z)/\theta_l(z) \times d\theta_l / dz + 1]\f$
  !!
  !! where K(z) is the hydraulic conductivity and psi(z) is the soil
  !! moisture suction at depth z, and b is an empirical parameter
  !! developed by Clapp and Hornberger (1978) \cite Clapp1978-898. K(z) and \f$\Psi(z)\f$ are
  !! calculated following Clapp and Hornberger as
  !!
  !! \f$K(z) = K_{sat} (\theta_l/\theta_p)^{(2b + 3)}\f$
  !!
  !! \f$\Psi(z) = \Psi_{sat} (\theta_l/\theta_p)^{(-b)}\f$
  !!
  !! where \f$K_{sat}\f$ and \f$\Psi_{sat}\f$ are the values of \f$K\f$ and \f$\Psi\f$ respectively
  !! at saturation.
  !!
  !! At the bottom of the permeable soil depth, in layer IGDR, if the
  !! liquid water content of the soil is greater than the field
  !! capacity, the vertical flow out of the bottom of the soil profile
  !! is calculated using a relation derived from Soulis et al. (2010):
  !!
  !! \f$F(z_b) = K_{sat} \times min{1, (\theta_l /\theta_p)/[1 – 1/(2b + 3)]^{(2b + 3)}}\f$
  !!
  !! This flow rate is multiplied by a drainage parameter XDRAIN,
  !! which is set to 0 if the soil is underlain by an impermeable
  !! layer (as in a bog), and to 1 otherwise.
  !!

  do I = IL1,IL2 ! loop 150
    if (IGRD(I) > 0) then
      FDT(I,1) = - EVAP(I) * DT(I)
      if (DELZW(I,IGDR(I)) > 0.0001) then
        if (THLIQ(I,IGDR(I)) > THFC(I,IGDR(I))) then
          CCH = 2.0 * BI(I,IGDR(I)) + 3.0
          ASATC = 1.0 - (1.0 / CCH)
          ASAT = THLIQ(I,IGDR(I)) / THPORF(I,IGDR(I))
          SATB = MIN(1.0,ASAT / ASATC)
          FDT(I,IGDR(I) + 1) = GRKSATF(I,IGDR(I)) * DT(I) * &
                               XDRAIN(I) * SATB ** CCH
        else
          FDT(I,IGDR(I) + 1) = 0.0
        end if
      else
        FDT(I,IGDR(I) + 1) = 0.0
      end if
    end if
  end do ! loop 150
  !
  !>
  !! Between soil layers, values of \f$K(z)\f$ and \f$\Psi(z)\f$ must be
  !! determined. This requires the estimation of soil properties at
  !! the layer interfaces. For \f$\theta_l\f$, \f$\theta_p\f$ and \f$d(\theta_l)/dz\f$, slightly
  !! different approaches are followed if the permeable soil layer
  !! thickness DELZW increases with depth (the normal case), or if it
  !! decreases between one soil layer and the next (indicating the
  !! presence of an impermeable barrier at some depth in the second
  !! layer). In the first case, the pore volume THPBND, liquid water
  !! content THLBND, and liquid water gradient DTHLDZ are simply
  !! calculated as arithmetic averages of the values above and below
  !! the interface. This has the effect of weighting the interface
  !! values towards the upper layer values, roughly emulating a
  !! surfaceward exponential decay curve of liquid water content. In
  !! the second case, in order to avoid a spurious weighting towards
  !! the lower layer, the gradients of liquid water content and pore
  !! volume between the upper and lower layers are calculated as
  !! linear relations, and then solved for the values at the
  !! interface.
  !!
  !! The Clapp and Hornberger b parameter at the interface is
  !! calculated as a simple average of the values in the upper and
  !! lower soil layers. The saturated hydraulic conductivity \f$K_{sat, bnd}\f$ is
  !! calculated as a harmonic mean, and the saturated soil moisture
  !! suction \f$\Psi_{sat, bnd}\f$ as a geometric mean, of the top and bottom
  !! layer values:
  !!
  !! \f$
  !! K_{sat, bnd} = K_{sat, t} K_{sat, b}
  !! ( \Delta z_{g, w, , t} + \Delta z_{g, w, , b} ) / ( K_{sat, t} \Delta z_{g, w, b} + K_{sat, b} \Delta z_{g, w, t} )
  !! \f$
  !!
  !! \f$
  !! \Psi_{sat, bnd} = \Psi_{sat, t}^{\Delta zg, w, t/(\Delta zg, w, t + \Delta zg, w, b)}
  !! \Psi_{sat, b}^{\Delta zg, w, b/(\Delta zg, w, t + \Delta zg, w, b)}
  !! \f$
  !!
  !! Finally, K(z), \f$\Psi(z)\f$ and F(z) at the interface are calculated
  !! from the equations given above. At the end of the loop, if
  !! set to zero behind the soil layer LZF containing the wetting
  !! front, and the flow rate at the bottom of LZF is constrained to
  !! be \f$\geq\f$ 0.
  !!
  do J = 1,IG - 1 ! loop 200
    do I = IL1,IL2
      if (IGRD(I) > 0) then
        if (J < IGDR(I)) then
          if (THPOR(I,J) > 0.0 .and. THPOR(I,J + 1) > 0.0 .and. &
              ISAND(I,J + 1) > - 3) then
            if (DELZW(I,J + 1) > DELZW(I,J)) then
              THPBND = (THPORF(I,J) + THPORF(I,J + 1)) / 2.0
              THLBND = (THLIQ(I,J) + THLIQ(I,J + 1)) / 2.0
              DTHLDZ = (THLIQ(I,J + 1) - THLBND) / DELZW(I,J + 1) + &
                       (THLBND - THLIQ(I,J)) / DELZW(I,J)
            else
              DTHLDZ = 2.0 * (THLIQ(I,J + 1) - THLIQ(I,J)) / &
                       (DELZW(I,J + 1) + DELZW(I,J))
              THLBND = THLIQ(I,J) + 0.5 * DTHLDZ * DELZW(I,J)
              DTHPDZ = 2.0 * (THPORF(I,J + 1) - THPORF(I,J)) / &
                       (DELZW(I,J + 1) + DELZW(I,J))
              THPBND = THPORF(I,J) + 0.5 * DTHPDZ * DELZW(I,J)
            end if
            BBND = (BI(I,J) + BI(I,J + 1)) / 2.0
            !                  GRSBND=GRKSAT(I,J)**(DELZW(I,J)/(DELZW(I,J)+
            !     1                DELZW(I,J+1)))*GRKSAT(I,J+1)**(DELZW(I,J+1)/
            !     2                (DELZW(I,J)+DELZW(I,J+1)))
            GRSBND = GRKSATF(I,J) * GRKSATF(I,J + 1) * (DELZW(I,J) + &
                     DELZW(I,J + 1)) / (GRKSATF(I,J) * DELZW(I,J + 1) + &
                     GRKSATF(I,J + 1) * DELZW(I,J))
            PSSBND = PSISAT(I,J) ** (DELZW(I,J) / (DELZW(I,J) + &
                     DELZW(I,J + 1))) * PSISAT(I,J + 1) ** (DELZW(I,J + 1) / &
                     (DELZW(I,J) + DELZW(I,J + 1)))
            GRK = MIN(GRSBND * (THLBND / THPBND) ** (2. * BBND + 3.), &
                  GRSBND)
            PSI = MAX(PSSBND * (THLBND / THPBND) ** ( - BBND),PSSBND)
            FDT(I,J + 1) = GRK * DT(I) * (( - BBND * PSI * DTHLDZ / THLBND) + 1.)
          else
            FDT(I,J + 1) = 0.0
          end if
          if (ABS(THLIQ(I,J) - THLIQ(I,J + 1)) < 0.05 .and. &
              FDT(I,J) < 0.0) FDT(I,J + 1) = 0.0
          if (LZF(I) > 0 .and. J < LZF(I)) FDT(I,J + 1) = 0.0
          if (LZF(I) > 0 .and. J == LZF(I) .and. FDT(I,J + 1) &
              < 0.0) FDT(I,J + 1) = 0.0
        end if
      end if
    end do
  end do ! loop 200
  !
  !     * CHECK FOR SUSTAINABLE EVAPORATION RATE FROM TOP SOIL LAYER; IF
  !     * LIQUID WATER SUPPLY IS INSUFFICIENT, TRY TO REMOVE WATER FROM
  !     * FROZEN SOIL MOISTURE.
  !
  !>
  !! In the next several loops, checks are carried out to ascertain
  !! whether the theoretically determined flow rates at the soil layer
  !! interfaces can be supported by the water contents in the layers.
  !! At the beginning of loop 250, the soil surface evaporation rate
  !! is addressed. A trial liquid water content THTEST is calculated
  !! for the first soil layer, reflecting the removal of the
  !! evaporated water. If THTEST is less than the minimum soil water
  !! content, F(0) is set to zero and \f$\theta_l\f$ is set to \f$\theta_{l, min}\f$. The
  !! excess surface flux that was not met by the removed liquid water
  !! is converted to a frozen water content, and an attempt is made to
  !! remove it from the frozen water in the layer. (The energy
  !! involved in converting the required frozen water mass to liquid
  !! water is obtained by decreasing the water temperature of the
  !! layer.) The frozen water content of the layer is adjusted to
  !! reflect the removal of the required amount. If the demand exceeds
  !! the supply of frozen water available, the frozen water content is
  !! set to zero, the diagnostic evaporative flux QFG at the soil
  !! surface is adjusted, and the remaining water flux not able to be
  !! met by the first soil layer is assigned to the variable WLOST.
  !!
  !! In the remainder of the 250 loop, checks are carried out for soil
  !! layers with liquid water contents already effectively equal to
  !! \f$\theta_{l, min}\f$. If the flux at the top of the layer is upward and that at
  !! the bottom of the layer is downward, both are set to zero. If
  !! both are downward, only the flux at the bottom is set to zero; if
  !! both are upward, only the flux at the top is set to zero. If the
  !! soil layer is an organic one and the liquid water content is less
  !! than the layer retention capacity THLRET and the flux at the
  !! bottom of the layer is downward, it is set to zero.
  !!
  IPTBAD = 0
  do J = 1,IG ! loop 250
    do I = IL1,IL2
      if (IGRD(I) > 0 .and. J == 1 .and. FDT(I,J) < 0. .and. &
          DELZW(I,J) > 0.0) then
        THTEST(I,J) = THLIQ(I,J) + FDT(I,J) / DELZW(I,J)
        if (THTEST(I,J) < THLMIN(I,J)) then
          FDT(I,J) = FDT(I,J) + (THLIQ(I,J) - THLMIN(I,J)) * &
                     DELZW(I,J)
          THLIQ(I,J) = THLMIN(I,J)
          WEXCES(I) = - FDT(I,J)
          FDT(I,J) = 0.0
          THSUBL = WEXCES(I) * RHOW / (RHOICE * DELZW(I,J))
          if (THLIQ(I,J) > 0.0) then
            TBARW(I,J) = TBARW(I,J) - (CLHMLT * RHOICE * THSUBL) / &
                         (HCPW * THLIQ(I,J))
          end if
          if (THSUBL <= THICE(I,J)) then
            THICE(I,J) = THICE(I,J) - THSUBL
          else
            THSUBL = THSUBL - THICE(I,J)
            THICE(I,J) = 0.0
            QFG(I) = QFG(I) - FI(I) * THSUBL * RHOICE * DELZW(I,J) / &
                     DELT
            WLOST(I) = WLOST(I) + THSUBL * RHOICE * DELZW(I,J)
          end if
        end if
        if (THICE(I,J) < 0.) IPTBAD = I
      end if
      !
      !     * ENSURE THAT CALCULATED WATER FLOWS BETWEEN SOIL LAYERS DO NOT
      !     * CAUSE LIQUID MOISTURE CONTENT OF ANY LAYER TO FALL BELOW THE
      !     * RESIDUAL VALUE OR TO EXCEED THE CALCULATED MAXIMUM.
      !
      if (IGRD(I) > 0) then
        if (THLIQ(I,J) <= (THLMIN(I,J) + 0.001) &
            .and. J <= IGDR(I)) then
          if (FDT(I,J) <= 0. .and. FDT(I,J + 1) >= 0.) then
            FDT(I,J) = 0.0
            FDT(I,J + 1) = 0.0
          else if (FDT(I,J) >= 0. .and. FDT(I,J + 1) > 0.) then
            FDT(I,J + 1) = 0.0
          else if (FDT(I,J) < 0. .and. FDT(I,J + 1) <= 0.) then
            FDT(I,J) = 0.0
          end if
        end if
      end if
      if (IGRD(I) > 0 .and. ISAND(I,J) == - 2 .and. &
          THLIQ(I,J) <= THLRET(I,J)) then
        if (FDT(I,J + 1) > 0.0) FDT(I,J + 1) = 0.0
      end if
    end do
  end do ! loop 250
  !
  if (IPTBAD /= 0) then
    write(6,6500) IPTBAD,JL,IVEG,THICE(IPTBAD,1)
6500 format('0AT (I,J) = (',I3,',',I3,'),IVEG = ',I2,' THICE(1) = ', &
             E13.5)
    call errorHandler('waterFlowNonInfiltrate', - 1)
  end if
  !
  !>
  !! In the 300 loop, checks are carried out for soil layers with
  !! liquid water contents already effectively equal to THLMAX. If the
  !! flux at the top of the layer is downward and that at the bottom
  !! of the layer is upward, both are set to zero. If both are
  !! downward, then if the top flux is greater than the bottom flux,
  !! it is set equal to the bottom flux. If both are upward, then if
  !! magnitude of the bottom flux is greater than that of the top
  !! flux, it is set equal to the top flux.
  !!
  do J = IG,1, - 1 ! loop 300
    do I = IL1,IL2
      if (IGRD(I) > 0) then
        if (THLIQ(I,J) >= (THLMAX(I,J) - 0.001) &
            .and. J <= IGDR(I)) then
          if (FDT(I,J) >= 0. .and. FDT(I,J + 1) <= 0.) then
            FDT(I,J) = 0.0
            FDT(I,J + 1) = 0.0
          else if (FDT(I,J) > 0. .and. FDT(I,J + 1) >= 0.) then
            if (FDT(I,J) > FDT(I,J + 1)) FDT(I,J) = FDT(I,J + 1)
          else if (FDT(I,J) <= 0. .and. FDT(I,J + 1) < 0.) then
            if (FDT(I,J + 1) < FDT(I,J)) FDT(I,J + 1) = FDT(I,J)
          end if
        end if
      end if
    end do
  end do ! loop 300
  !
  !>
  !! In the 400 loop, for each soil layer THTEST is recalculated as
  !! the liquid water content resulting from the updated boundary
  !! fluxes, and is compared with a residual value THLTHR. For organic
  !! soil layers deeper than the first layer, THLTHR is set to the
  !! minimum of \f$\theta_{l, ret}\f$ and \f$\theta_l\f$, since only evapotranspiration can
  !! cause the soil moisture to fall below the retention capacity.
  !! (The first layer is excepted because surface evaporation can
  !! drive its moisture content below \f$\theta_{l, ret}\f$.) For mineral soils,
  !! \f$\theta_{l, min}\f$ is set to THLMIN. If THTEST < THLTHR, then if the flow at
  !! the bottom of the soil layer is downward, it is recalculated as
  !! the sum of the flow at the top of the layer plus the amount of
  !! water that must be removed to make the liquid water content of
  !! the layer equal to THLTHR. If the flow at the bottom of the layer
  !! is upward, the flow at the top of the layer is recalculated as
  !! the sum of the flow at the bottom of the layer minus the amount
  !! of water that must be removed to make the liquid water content of
  !! the layer equal to THLTHR. THTEST of the current layer is then
  !! reset to THLTHR, and THTEST of the overlying and underlying
  !! layers (if any) are recalculated using the new interface fluxes.
  !!
  do J = 1,IG ! loop 400
    do I = IL1,IL2
      if (IGRD(I) > 0) then
        if (J <= IGDR(I) .and. ISAND(I,J) /= - 3) then
          THTEST(I,J) = THLIQ(I,J) + (FDT(I,J) - FDT(I,J + 1)) / DELZW(I,J)
          if (ISAND(I,J) == - 2 .and. J /= 1) then
            THLTHR = MIN(THLRET(I,J),THLIQ(I,J))
          else
            THLTHR = THLMIN(I,J)
          end if
          if (THTEST(I,J) < THLTHR) then
            if (FDT(I,J + 1) > 0.) then
              FDT(I,J + 1) = FDT(I,J) + (THLIQ(I,J) - THLTHR) * &
                             DELZW(I,J)
            else
              FDT(I,J) = FDT(I,J + 1) - (THLIQ(I,J) - THLTHR) * &
                         DELZW(I,J)
            end if
            THTEST(I,J) = THLTHR
            if (J < IGDR(I)) then
              if (DELZW(I,J + 1) > 0.0) THTEST(I,J + 1) = THLIQ(I,J + 1) &
                  + (FDT(I,J + 1) - FDT(I,J + 2)) / DELZW(I,J + 1)
            end if
            if (J > 1) then
              if (DELZW(I,J - 1) > 0.0) THTEST(I,J - 1) = THLIQ(I,J - 1) &
                  + (FDT(I,J - 1) - FDT(I,J)) / DELZW(I,J - 1)
            end if
          end if
        else
          THTEST(I,J) = 0.0
        end if
      end if
    end do
  end do ! loop 400
  !
  !>
  !! In the 500 loop, for each soil layer THTEST is compared with
  !! THLMAX. If THTEST > THLMAX, two temporary variables are defined:
  !! WLIMIT as the amount of water required to raise the liquid water
  !! content of the soil layer to THLMAX, and WEXCES as the amount of
  !! water representing the excess of THTEST over THLMAX. If the flux
  !! at the top of the layer is downward and the flux at the bottom of
  !! the layer is upward, then if the negative of the flux at the
  !! bottom of the layer is greater than WLIMIT, it is set equal to
  !! -WLIMIT and the flux at the top is set to zero; otherwise WEXCES
  !! is subtracted from the flux at the top. If both the flux at the
  !! top and the flux at the bottom are downward, then WEXCES is
  !! subtracted from the flux at the top. If both are upward, then
  !! WEXCES is added to the flux at the bottom, and a correction is
  !! applied to the flux at the bottom of the underlying layer.
  !! Finally, THTEST is recalculated for all the soil layers using the
  !! new interface fluxes.
  !!
  do J = IG,1, - 1 ! loop 500
    do I = IL1,IL2
      if (IGRD(I) > 0) then
        if (THTEST(I,J) > THLMAX(I,J) .and. J <= IGDR(I)) then
          WLIMIT = MAX((THLMAX(I,J) - THLIQ(I,J)),0.0) * DELZW(I,J)
          WEXCES(I) = (THTEST(I,J) - THLMAX(I,J)) * DELZW(I,J)
          if (FDT(I,J) > 0. .and. FDT(I,J + 1) <= 0.) then
            if ( - FDT(I,J + 1) > WLIMIT) then
              FDT(I,J + 1) = - WLIMIT
              FDT(I,J) = 0.0
            else
              FDT(I,J) = FDT(I,J) - WEXCES(I)
            end if
            !                IF (FDT(I,J)>=WLIMIT)             THEN
            !                   FDT(I,J)=WLIMIT
            !                   FDT(I,J+1)=0.0
            !                ELSE
            !                   FDT(I,J+1)=FDT(I,J)-WLIMIT
            !                END IF
          else if (FDT(I,J) > 0. .and. FDT(I,J + 1) >= 0.) then
            FDT(I,J) = FDT(I,J) - WEXCES(I)
          else if (FDT(I,J) <= 0. .and. FDT(I,J + 1) < 0.) then
            FDT(I,J + 1) = FDT(I,J + 1) + WEXCES(I)
            if (J < IGDR(I)) then
              if (FDT(I,J + 2) < 0.) FDT(I,J + 2) = 0.0
            end if
          end if
          do K = 1,IG ! loop 450
            if (DELZW(I,K) > 0.0) then
              THTEST(I,K) = THLIQ(I,K) + (FDT(I,K) - FDT(I,K + 1)) / &
                            DELZW(I,K)
            end if
          end do ! loop 450
        end if
      end if
    end do
  end do ! loop 500
  !
  !>
  !! In the first part of the 600 loop, a correction is performed in
  !! the unlikely event that as a result of the above adjustments, the
  !! flux at the bottom of the last permeable soil layer has become
  !! upward. If this occurs, all of the fluxes above are corrected by
  !! the inverse of this same amount. However, this will result in a
  !! small spurious upward water flux at the surface. Since the liquid
  !! water flows are now accounted for, an attempt is made, as in loop
  !! 250, to remove the required water from the frozen water store in
  !! the first layer. The excess that cannot be supplied from this
  !! source is assigned to WLOST.
  !!
  !! In the second part of the 600 loop, the temperatures TFDT of the
  !! water fluxes at the top and bottom of the layers in the permeable
  !! soil profile are set equal to the temperatures of the water in
  !! the first and last layers respectively. The flow at the bottom of
  !! the profile and the temperature of the flow are used to update
  !! BASFLW and TBASFL, the baseflow from the current subarea and its
  !! temperature, and RUNOFF and TRUNOF, the total runoff from the
  !! grid cell in question and its temperature.
  !!
  IPTBAD = 0
  do I = IL1,IL2
    if (IGRD(I) > 0) then
      if (FDT(I,IGDR(I) + 1) < 0.) then
        WEXCES(I) = - FDT(I,IGDR(I) + 1)
        do J = 1,IGDR(I) + 1
          FDT(I,J) = FDT(I,J) + WEXCES(I)
        end do ! loop 550
        THSUBL = WEXCES(I) * RHOW / (RHOICE * DELZW(I,1))
        if (THLIQ(I,1) > 0.0) then
          TBARW(I,1) = TBARW(I,1) - (CLHMLT * RHOICE * THSUBL) / &
                       (HCPW * THLIQ(I,1))
        end if
        if (THSUBL <= THICE(I,1)) then
          THICE(I,1) = THICE(I,1) - THSUBL
        else
          THSUBL = THSUBL - THICE(I,1)
          THICE(I,1) = 0.0
          QFG(I) = QFG(I) - FI(I) * THSUBL * RHOICE * DELZW(I,1) / &
                   DELT
          WLOST(I) = WLOST(I) + THSUBL * RHOICE * DELZW(I,1)
        end if
        if (THICE(I,1) < 0.0) IPTBAD = I
      end if
      !
      !     * CALCULATE DRAINAGE FROM BOTTOM OF SOIL COLUMN AND RE-EVALUATE
      !     * SOIL LAYER TEMPERATURES AND LIQUID MOISTURE CONTENTS AFTER
      !     * WATER MOVEMENT.
      !
      TFDT(I,1) = TBARW(I,1)
      TFDT(I,IGDR(I) + 1) = TBARW(I,IGDR(I))
      if (FDT(I,IGDR(I) + 1) > 0.0) then
        if ((BASFLW(I) + FI(I) * FDT(I,IGDR(I) + 1)) > 0.0) &
            TBASFL(I) = (TBASFL(I) * BASFLW(I) + (TFDT(I,IGDR(I) + 1) + &
            TFREZ) * FI(I) * FDT(I,IGDR(I) + 1)) / (BASFLW(I) + FI(I) * &
            FDT(I,IGDR(I) + 1))
        BASFLW(I) = BASFLW(I) + FI(I) * FDT(I,IGDR(I) + 1)
        if ((RUNOFF(I) + FDT(I,IGDR(I) + 1)) > 0.0) &
            TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + (TFDT(I,IGDR(I) + 1) + &
            TFREZ) * FDT(I,IGDR(I) + 1)) / (RUNOFF(I) + &
            FDT(I,IGDR(I) + 1))
        RUNOFF(I) = RUNOFF(I) + FDT(I,IGDR(I) + 1)
      end if
    end if
  end do ! loop 600
  !
  if (IPTBAD /= 0) then
    write(6,6500) IPTBAD,JL,IVEG,THICE(IPTBAD,1)
    call errorHandler('waterFlowNonInfiltrate', - 3)
  end if
  !
  !>
  !! In the 700 loop, for each successive soil layer the temperature
  !! of the water flux at the bottom of the layer is first assigned.
  !! If the flux is downward, TFDT is set to the temperature of the
  !! water in the current layer; if it is upward, TFDT is set to the
  !! temperature of the water in the layer below. The temperature of
  !! the water in the current layer is then updated using the FDT and
  !! TFDT values at the top and bottom of the layer, and the liquid
  !! water content of the layer is set to THTEST.
  !!
  do J = 1,IG ! loop 700
    do I = IL1,IL2
      if (IGRD(I) > 0) then
        if (J <= IGDR(I)) then
          if (J < IGDR(I)) then
            if (FDT(I,J + 1) > 0.) then
              TFDT(I,J + 1) = TBARW(I,J)
            else
              TFDT(I,J + 1) = TBARW(I,J + 1)
            end if
          end if
          if (THTEST(I,J) > 0.0 .and. DELZW(I,J) > 0.0) then
            TBARW(I,J) = (THLIQ(I,J) * TBARW(I,J) + (FDT(I,J) * TFDT(I,J) - &
                         FDT(I,J + 1) * TFDT(I,J + 1)) / DELZW(I,J)) / &
                         THTEST(I,J)
          end if
          THLIQ(I,J) = THTEST(I,J)
        end if
      end if
    end do
  end do ! loop 700

  return
end subroutine waterFlowNonInfiltrate
