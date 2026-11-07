!> \file
!! Recalculates liquid water content of soil layers after
!! infiltration, and evaluates baseflow.
!! @author D. Verseghy, M. Lazare
!
subroutine waterBaseflow (THLIQX, THICEX, TBARWX, ZPOND, TPOND, & ! Formerly WEND
                          BASFLW, TBASFL, RUNOFF, TRUNOF, FI, &
                          WMOVE, TMOVE, LZF, NINF, TRMDR, THLINF, DELZX, &
                          ZMAT, ZRMDR, FDTBND, WADD, TADD, FDT, TFDT, &
                          THLMAX, THTEST, THLDUM, THIDUM, TDUMW, &
                          TUSED, RDUMMY, ZERO, WEXCES, XDRAIN, &
                          THPOR, THLRET, THLMIN, BI, PSISAT, GRKSAT, &
                          THFC, DELZW, ISAND, IGRN, IGRD, IGDR, IZERO, &
                          IVEG, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !
  !     * OCT 18/11 - M.LAZARE.   PASS IN "IGDR" AS AN INPUT FIELD
  !     *                         (ORIGINATING IN soilProperties) RATHER
  !     *                         THAN REPEATING THE CALCULATION HERE
  !     *                         AS AN INTERNAL WORK FIELD.
  !     * DEC 15/10 - D.VERSEGHY. ALLOW FOR BASEFLOW WHEN BEDROCK
  !     *                         LIES WITHIN SOIL PROFILE.
  !     * JAN 06/09 - D.VERSEGHY. ADD ZPOND AND TPOND TO SUBROUTINE
  !     *                         CALL; ASSIGN RESIDUAL OF WMOVE TO
  !     *                         PONDED WATER; REVISE LOOP 550
  !     *                         DELETE CALCULATION OF FDTBND.
  !     * MAY 17/06 - D.VERSEGHY. PROTECT AGAINST DIVISIONS BY ZERO.
  !     * OCT 21/05 - D.VERSEGHY. FIX MINOR BUGS IN CLEANUP AND
  !     *                         RUNOFF TEMPERATURE CALCULATION.
  !     * MAR 23/05 - D.VERSEGHY. ADD VARIABLES TO waterFlowNonInfiltrate CALL
  !     *                         ADD CALCULATION OF RUNOFF
  !     *                         TEMPERATURE.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * APR 24/03 - D.VERSEGHY. ADD CHECK FOR OVERFLOW IN SOIL
  !     *                         LAYER CONTAINING WETTING FRONT.
  !     * OCT 15/02 - D.VERSEGHY. BUGFIX IN CALCULATION OF FDTBND
  !     *                         (PRESENT ONLY IN PROTOTYPE
  !     *                         VERSIONS OF CLASS VERSION 3.0).
  !     * JUN 21/02 - D.VERSEGHY. UPDATE SUBROUTINE CALL; SHORTENED
  !     *                         CLASS4 COMMON BLOCK.
  !     * DEC 12/01 - D.VERSEGHY. ADD SEPARATE CALCULATION OF BASEFLOW
  !     *                         AT BOTTOM OF SOIL COLUMN.
  !     * OCT 20/97 - D.VERSEGHY. APPLY ACCURACY LIMIT ON FLOWS IN AND
  !     *                         OUT OF LAYER CONTAINING WETTING FRONT,
  !     *                         IN ORDER TO ENSURE MOISTURE CONSERVATION.
  !     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
  !     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
  !     *                         SOIL PERMEABLE DEPTH.
  !     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
  !     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
  !     *                         BETWEEN SOIL LAYERS.
  !     * APR 24/92 - D.VERSEGHY, M.LAZARE. CLASS - VERSION 2.1.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. RECALCULATE LIQUID MOISTURE CONTENT
  !     *                         OF SOIL LAYERS AFTER INFILTRATION
  !     *                         AND EVALUATE FLOW ("RUNOFF") FROM
  !     *                         BOTTOM OF SOIL COLUMN.
  !
  use classicParams, only : DELT, TFREZ

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IVEG, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N
  integer :: I, J, K
  !
  !     * OUTPUT FIELDS.
  !
  real, intent(inout) :: THLIQX(ILG,IGP1) !< Volumetric liquid water content of soil
  !< layer \f$[m^3 m^{-3}]\f$
  real, intent(in)    :: THICEX(ILG,IGP1) !< Volumetric frozen water content of soil
  !< layer \f$[m^3 m^{-3}]\f$
  real, intent(inout) :: TBARWX(ILG,IGP1) !< Temperature of water in soil layer [C]
  real, intent(inout) :: ZPOND (ILG)      !< Depth of ponded water [m]
  real, intent(inout) :: TPOND (ILG)      !< Temperature of ponded water [C]
  real, intent(inout) :: BASFLW(ILG)      !< Base flow from bottom of soil column \f$[kg m^{-2}]\f$
  real, intent(inout) :: TBASFL(ILG)      !< Temperature of base flow from bottom of soil column [K]
  real, intent(inout) :: RUNOFF(ILG)      !< Total runoff from soil column [m]
  real, intent(inout) :: TRUNOF(ILG)      !< Temperature of total runoff from soil column [K]
  !
  !     * INPUT FIELDS.
  !
  real, intent(inout) :: WMOVE (ILG,IGP2) !< Water movement matrix \f$[m^3 m^{-2}]\f$
  real, intent(in)    :: TMOVE (ILG,IGP2) !< Temperature matrix associated with ground water movement [C]
  real, intent(in)    :: THLINF(ILG,IGP1) !< Volumetric liquid water content behind the
  !< wetting front \f$[m^3 m^{-3}]\f$
  real, intent(in)    :: FI    (ILG)      !< Fractional coverage of subarea in question on modelled area [ ]
  real, intent(in)    :: TRMDR (ILG)      !< Time remaining in current time step [s]
  real, intent(in)    :: DELZX (ILG,IGP1) !< Permeable depth of soil layer [m]
  !
  integer, intent(in) :: LZF(ILG), IGRN(ILG)
  integer, intent(inout) :: NINF(ILG)
  !
  !     * INTERNAL WORK ARRAYS.
  !
  real, intent(inout) :: ZMAT(ILG,IGP2,IGP1), ZRMDR(ILG,IGP1), &
                         FDTBND(ILG), WADD(ILG), TADD(ILG)
  !
  !     * INTERNAL ARRAYS USED IN CALLING waterFlowNonInfiltrate.
  !
  real, intent(in)    :: TFDT(ILG,IGP1)
  real, intent(inout) :: FDT (ILG,IGP1)
  !
  real, intent(in)    :: THLMAX(ILG,IG), THTEST(ILG,IG)
  real, intent(inout) :: THLDUM(ILG,IG), THIDUM(ILG,IG), TDUMW(ILG,IG)
  !
  real, intent(in)    :: RDUMMY(ILG), ZERO(ILG), WEXCES(ILG)
  real, intent(inout) :: TUSED (ILG)
  !
  integer, intent(in) :: IGRD(ILG), IZERO(ILG), IGDR(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: WREM, TREM, THDRAN, THINFL, WDRA, TDRA
  !
  !     * SOIL INFORMATION ARRAYS.
  !
  real, intent(in) :: THPOR (ILG,IG)   !< Pore volume in soil layer \f$[m^3 m^{-3}] (\theta_p)\f$
  real, intent(in) :: THLRET(ILG,IG)   !< Liquid water retention capacity for organic
  !< soil \f$[m^3 m^{-3} ] (\theta_{l,ret})\f$
  real, intent(in) :: THLMIN(ILG,IG)   !< Residual soil liquid water content (variable, accounting for unfrozen
  !< water in frozen soils) remaining after freezing or evaporation \f$[m^3 m^{-3}]\f$
  real, intent(in) :: BI    (ILG,IG)   !< Clapp and Hornberger empirical "b" parameter [ ] (b)
  real, intent(in) :: PSISAT(ILG,IG)   !< Soil moisture suction at saturation \f$[m] (\Psi_{sat})\f$
  real, intent(in) :: GRKSAT(ILG,IG)   !< Hydraulic conductivity of soil at
  !< saturation \f$[m s^{-1}] (K_{sat})\f$
  real, intent(in) :: THFC  (ILG,IG)   !< Field capacity \f$[m^3 m^{-3}]\f$
  real, intent(in) :: DELZW (ILG,IG)   !< Permeable depth of soil layer \f$[m] (\Delta_{zg,w})\f$
  real, intent(in) :: XDRAIN(ILG)      !< Drainage index for water flow at bottom of soil profile [ ]
  !
  integer, intent(in) :: ISAND (ILG,IG)
  !
  !-----------------------------------------------------------------------
  !
  !>
  !! At levels in the soil profile lower than the bottom of the
  !! wetting front, redistribution of soil liquid water proceeds in
  !! response to normal gravity and suction forces. Subroutine waterFlowNonInfiltrate
  !! is called to calculate these redistributions. The time period
  !! TUSED that is passed to waterFlowNonInfiltrate is set in the 100 loop to the time
  !! period over which infiltration was occurring during the current
  !! time step, except if the wetting front has passed the bottom of
  !! the lowest soil layer, in which case TUSED is set to zero (since
  !! the flows calculated by waterFlowNonInfiltrate are not required). waterFlowNonInfiltrate is
  !! called using dummy variables THLDUM, THIDUM and TDUMW, which are
  !! set in loop 125 to the liquid water content, the frozen water
  !! content and the water temperature of the soil layers
  !! respectively.
  !!
  !
  !     * DETERMINE AMOUNT OF TIME OUT OF CURRENT MODEL STEP DURING WHICH
  !     * INFILTRATION WAS OCCURRING.
  !     * SET WORK ARRAY "TUSED" TO ZERO FOR POINTS WHERE WETTING FRONT
  !     * IS BELOW BOTTOM OF LOWEST SOIL LAYER TO SUPPRESS CALCULATIONS
  !     * DONE IN "waterFlowNonInfiltrate".
  !
  do I = IL1,IL2 ! loop 100
    if (IGRN(I) > 0 .and. LZF(I) <= IG) then
      TUSED(I) = DELT - TRMDR(I)
    else
      TUSED(I) = 0.
    end if
  end do ! loop 100
  !
  !     * INITIALIZATION.
  !
  do J = 1,IG ! loop 125
    do I = IL1,IL2
      if (IGRN(I) > 0) then
        THLDUM(I,J) = THLIQX(I,J)
        THIDUM(I,J) = THICEX(I,J)
        TDUMW (I,J) = TBARWX(I,J)
      end if
    end do
  end do ! loop 125
  !
  !     * CALL "waterFlowNonInfiltrate" WITH COPIES OF CURRENT LIQUID AND FROZEN SOIL
  !     * MOISTURE CONTENTS AND LAYER TEMPERATURES TO DETERMINE MOISTURE
  !     * FLOW BETWEEN LAYERS BELOW THE WETTING FRONT.
  !
  call waterFlowNonInfiltrate(IVEG, THLDUM, THIDUM, TDUMW, FDT, TFDT, RDUMMY, RDUMMY, & ! Formerly GRDRAN
                              RDUMMY, RDUMMY, RDUMMY, RDUMMY, FI, ZERO, ZERO, ZERO, &
                              TUSED, WEXCES, THLMAX, THTEST, THPOR, THLRET, THLMIN, &
                              BI, PSISAT, GRKSAT, THFC, DELZW, XDRAIN, ISAND, LZF, &
                              IZERO, IGRD, IGDR, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !>
  !! After waterFlowNonInfiltrate has been called, the maximum value of the water
  !! movement index NINF is set to the number of soil layers plus 1.
  !! The values in the matrix ZRMDR, representing for each soil layer
  !! the depth that has not been affected by infiltration, are
  !! initialized to the soil permeable layer thicknesses DELZX. The
  !! water flows FDT coming out of waterFlowNonInfiltrate are set to zero at the soil
  !! layer interfaces above the wetting front. For the layer
  !! containing the wetting front, if the flow at the bottom of the
  !! layer is upward, it is set to zero (to avoid possible overflows
  !! in liquid water content).
  !!
  !
  !     * INITIALIZATION OF ARRAYS IN PREPARATION FOR RE-ALLOCATION OF
  !     * MOISTURE STORES WITHIN SOIL LAYERS; SUPPRESS WATER FLOWS
  !     * CALCULATED IN waterFlowNonInfiltrate ABOVE WETTING FRONT; CONSISTENCY CHECK
  !     * FOR WATER FLOWS INTO LAYER CONTAINING WETTING FRONT.
  !
  do I = IL1,IL2 ! loop 150
    if (IGRN(I) > 0) then
      NINF(I) = MIN(NINF(I),IGP1)
    end if
  end do ! loop 150
  !
  do J = IGP1,1, - 1 ! loop 200
    do I = IL1,IL2
      if (IGRN(I) > 0) then
        ZRMDR(I,J) = DELZX(I,J)
        if (J <= LZF(I)) then
          FDT(I,J) = 0.0
        end if
        if (J < IGP1) then
          if (J == LZF(I) .and. FDT(I,J + 1) < 0.) then
            FDT(I,J + 1) = 0.0
          end if
        end if
      end if
    end do
  end do ! loop 200
  !
  !>
  !! The values in the three-dimensional matrix ZMAT are initialized
  !! to zero. This matrix contains the depth of each soil layer J
  !! (including the dummy soil layer below the lowest layer) that is
  !! filled by water from level K in the water movement matrix WMOVE.
  !! In WMOVE, the first level contains the amount of water that has
  !! infiltrated at the surface during the time period in question,
  !! and each successive level K contains the amount of water in soil
  !! layer K-1 that has been displaced during the infiltration, down
  !! to the soil layer containing the wetting front. Thus, the number
  !! of levels in WMOVE that are used in the current infiltration
  !! calculations, NINF, is equal to LZF+1, where LZF is the index of
  !! the soil layer containing the wetting front; or to IGP1, the
  !! number of soil layers IG plus 1, if LZF is greater than IG, i.e.
  !! if the wetting front has penetrated below the bottom of the
  !! lowest soil layer into the underlying dummy layer. In the 400
  !! loop, starting at the top of the soil profile, an attempt is made
  !! to assign each layer of WMOVE, converted into a depth by using
  !! the volumetric water content THLINF behind the wetting front for
  !! the soil layer, in turn to the K, J level of ZMAT. If the
  !! calculated value of ZMAT is greater than the available depth
  !! ZRMDR of the layer, ZMAT is set to ZRMDR, WMOVE is decremented by
  !! ZRMDR converted back to a water amount, and ZRMDR is set to zero.
  !! Otherwise the calculated value of ZMAT is accepted, ZRMDR is
  !! decremented by ZMAT, and WMOVE is set to zero. At the end of
  !! these calculations, any remaining residual amounts in the WMOVE
  !! matrix are assigned to ponded water, and the ponded water
  !! temperature is updated accordingly.
  !!
  do J = 1,IGP1 ! loop 300                                                              <<< unterminated do loop
    do K = 1,IGP1
      do I = IL1,IL2
        if (IGRN(I) > 0 .and. K <= NINF(I)) ZMAT(I,K,J) = 0.0
      end do
    end do
  end do ! loop 300
  !
  !     * ASSIGN VALUES IN MATRIX "ZMAT": DETERMINE DEPTH OUT OF EACH
  !     * SOIL LAYER J WHICH IS FILLED BY WATER FROM RESERVOIR K
  !     * IN "WMOVE"; FIND THE DEPTH "ZRMDR" LEFT OVER WITHIN EACH
  !     * SOIL LAYER.
  !
  do K = 1,IGP1 ! loop 400
    do J = 1,IGP1
      do I = IL1,IL2
        if (IGRN(I) > 0 .and. K <= NINF(I)) then
          if (ZRMDR(I,J) > 1.0E-5 .and. WMOVE(I,K) > 0.) then
            ZMAT(I,K,J) = WMOVE(I,K) / THLINF(I,J)
            if (ZMAT(I,K,J) > ZRMDR(I,J)) then
              ZMAT(I,K,J) = ZRMDR(I,J)
              WMOVE(I,K) = WMOVE(I,K) - ZRMDR(I,J) * THLINF(I,J)
              ZRMDR(I,J) = 0.0
            else
              ZRMDR(I,J) = ZRMDR(I,J) - ZMAT(I,K,J)
              WMOVE(I,K) = 0.0
            end if
          end if
        end if
      end do
    end do
  end do ! loop 400
  !
  do J = 1,IGP1 ! loop 450
    do I = IL1,IL2
      if (IGRN(I) > 0 .and. WMOVE(I,J) > 0.0) then
        TPOND(I) = (TPOND(I) * ZPOND(I) + TMOVE(I,J) * WMOVE(I,J)) / &
                   (ZPOND(I) + WMOVE(I,J))
        ZPOND(I) = ZPOND(I) + WMOVE(I,J)
      end if
    end do
  end do ! loop 450
  !>
  !! As a result of the above processes, the liquid water content and
  !! temperature of each soil layer (excluding the bottom, dummy
  !! layer) will be a combined result of infiltration processes (WADD,
  !! TADD), redistribution processes (WDRA, TDRA), and water in the
  !! layer that has remained unaffected (WREM, TREM). For each layer,
  !! WADD is calculated by summing over the respective ZMAT values
  !! corresponding to that layer multiplied by THLINF, and TADD by
  !! summing over the ZMAT and THLINF values multiplied by the
  !! respective TMOVE values. WREM is obtained as the product of the
  !! original water content THLIQX multiplied by ZRMDR, and TREM as
  !! the product of the water temperature TBARWX, THLIQX and ZRMDR.
  !! For the soil layer containing the wetting front, a check is
  !! carried out to determine whether the liquid moisture content
  !! resulting from the infiltration and drainage processes, THINFL,
  !! is less than the residual liquid moisture content THLMIN. If so,
  !! the flow FDT at the bottom of the layer is recalculated as the
  !! value required to keep THLIQX at THLMIN. WDRA is obtained from
  !! the difference between the water fluxes FDT at the top and bottom
  !! of the layer, supplied by waterFlowNonInfiltrate, and TDRA is obtained from the
  !! water fluxes FDT and their corresponding temperatures TFDT.
  !! Finally, THLIQX is calculated as the sum of WADD, WREM and WDRA
  !! normalized by DELZX, and TBARWX as the sum of TADD, TREM and TDRA
  !! normalized by the product of THLIQX and DELZX.
  !!
  !
  !     * ADD WATER CONTENT AND TEMPERATURE CHANGES DUE TO INFILTRATION
  !     * (WADD, TADD) AND DRAINAGE (WDRA, TDRA) TO WATER REMAINING IN
  !     * EACH SOIL LAYER AFTER THESE PROCESSES (WREM, TREM).
  !
  do J = IG,1, - 1 ! loop 600
    do I = IL1,IL2 ! loop 500
      if (IGRN(I) > 0) then
        WADD(I) = 0.
        TADD(I) = 0.
      end if
    end do ! loop 500
    !
    do K = 1,IGP1 ! loop 525
      do I = IL1,IL2
        if (IGRN(I) > 0 .and. K <= NINF(I)) then
          WADD(I) = WADD(I) + THLINF(I,J) * ZMAT(I,K,J)
          TADD(I) = TADD(I) + TMOVE(I,K) * THLINF(I,J) * ZMAT(I,K,J)
        end if
      end do
    end do ! loop 525
    !
    do I = IL1,IL2 ! loop 550
      if (IGRN(I) > 0 .and. DELZW(I,J) > 1.0E-4) then
        if (ZRMDR(I,J) > 1.0E-5) then
          WREM = THLIQX(I,J) * ZRMDR(I,J)
          TREM = TBARWX(I,J) * THLIQX(I,J) * ZRMDR(I,J)
        else
          WREM = 0.0
          TREM = 0.0
        end if
        if (J == LZF(I)) then
          THINFL = (WADD(I) + WREM + FDT(I,J) - FDT(I,J + 1)) / DELZW(I,J)
          if (THINFL < THLMIN(I,J)) then
            FDT(I,J + 1) = WADD(I) + WREM + FDT(I,J) - THLMIN(I,J) * &
                           DELZW(I,J)
          end if
        end if
        WDRA = FDT(I,J) - FDT(I,J + 1)
        TDRA = FDT(I,J) * TFDT(I,J) - FDT(I,J + 1) * TFDT(I,J + 1)
        THLIQX(I,J) = (WADD(I) + WREM + WDRA) / DELZW(I,J)
        THLIQX(I,J) = MAX(THLIQX(I,J),THLMIN(I,J))
        TBARWX(I,J) = (TADD(I) + TREM + TDRA) / (THLIQX(I,J) * &
                      DELZW(I,J))
      end if
    end do ! loop 550
  end do ! loop 600

  !>
  !! Lastly, the base flow BASFLW at the bottom of the soil profile
  !! and its temperature TBASFL are calculated. If the wetting front
  !! is located in the dummy soil layer below the soil profile, BASFLW
  !! is obtained by summing over the ZMAT values for the IGP1 level,
  !! multiplied by the dummy layer THLINF value and the fractional
  !! coverage FI of the modelled subarea. TBASFL is similarly obtained
  !! as the weighted average of the original TBASFL and the values of
  !! TMOVE corresponding to the ZMAT values. The overall subarea
  !! runoff RUNOFF and its temperature TRUNOF are calculated in the
  !! same manner without the FI weightings. Otherwise, the baseflow
  !! and total runoff are obtained from the value of FDT at the bottom
  !! of the IGDR layer, and their temperatures from the values of TFDT
  !! and FDT.
  !!
  !
  !     * CALCULATE FLOW OUT OF BOTTOM OF SOIL COLUMN DUE TO INFILTRATION
  !     * AND GRAVITY DRAINAGE AND ADD TO TOTAL RUNOFF AND BASEFLOW.
  !
  do K = 1,IGP1 ! loop 700
    do I = IL1,IL2
      if (IGRN(I) > 0) then
        if (LZF(I) == IGP1 .and. K <= NINF(I) .and. &
            THLINF(I,IGP1) * ZMAT(I,K,IGP1) > 0.0) then
          TBASFL(I) = (TBASFL(I) * BASFLW(I) + FI(I) * (TMOVE(I,K) + &
                      TFREZ) * THLINF(I,IGP1) * ZMAT(I,K,IGP1)) / (BASFLW(I) + &
                      FI(I) * THLINF(I,IGP1) * ZMAT(I,K,IGP1))
          BASFLW(I) = BASFLW(I) + FI(I) * THLINF(I,IGP1) * &
                      ZMAT(I,K,IGP1)
          TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + (TMOVE(I,K) + TFREZ) * &
                      THLINF(I,IGP1) * ZMAT(I,K,IGP1)) / (RUNOFF(I) + &
                      THLINF(I,IGP1) * ZMAT(I,K,IGP1))
          RUNOFF(I) = RUNOFF(I) + THLINF(I,IGP1) * ZMAT(I,K,IGP1)
        else if (K == (IGDR(I) + 1) .and. FDT(I,K) > 1.0E-8) then
          TBASFL(I) = (TBASFL(I) * BASFLW(I) + FI(I) * (TFDT(I,K) + &
                      TFREZ) * FDT(I,K)) / (BASFLW(I) + FI(I) * FDT(I,K))
          BASFLW(I) = BASFLW(I) + FI(I) * FDT(I,K)
          TRUNOF(I) = (TRUNOF(I) * RUNOFF(I) + (TFDT(I,K) + TFREZ) * &
                      FDT(I,K)) / (RUNOFF(I) + FDT(I,K))
          RUNOFF(I) = RUNOFF(I) + FDT(I,K)
        end if
      end if
    end do
  end do ! loop 700
  !
  return
end subroutine waterBaseflow
