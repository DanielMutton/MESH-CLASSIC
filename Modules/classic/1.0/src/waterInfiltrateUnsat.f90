!> \file
!! Evaluates infiltration of water into soil under
!! unsaturated conditions.
!! @author D. Verseghy, M. Lazare
!
subroutine waterInfiltrateUnsat (WMOVE, TMOVE, LZF, NINF, ZF, TRMDR, R, TR, & ! Formerly WFILL
                                 PSIF, GRKINF, THLINF, THLIQX, TBARWX, &
                                 DELZX, ZBOTX, DZF, TIMPND, WADJ, WADD, &
                                 IFILL, IFIND, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)
  !
  !     * JAN 06/09 - D.VERSEGHY. CORRECT LZF AND ZF ASSIGNMENTS IN LOOP
  !     *                         100; ADDITIONAL DZF CHECK IN LOOP 400.
  !     * MAR 22/06 - D.VERSEGHY. MOVE IFILL TEST OUTSIDE ALL IF BLOCKS.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 29/04 - D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS
  !     *                         AGAINST ROUNDOFF ERRORS.
  !     * JUN 21/02 - D.VERSEGHY. UPDATE SUBROUTINE CALL.
  !     * MAY 17/99 - D.VERSEGHY. PUT LIMIT ON CONDITION BASED ON "GRKINF"
  !     *                         SO THAT "LZF" IS ALWAYS INITIALIZED.
  !     * NOV 30/94 - M.LAZARE.   BRACKET TERMS IN "WADJ" CALCULATION IN
  !     *                         LOOP 200 TO AVOID OPTIMIZATION LEADING
  !     *                         TO RE-ORDERING OF CALCULATION AND
  !     *                         RARE VERY LARGE ITERATION LIMITS.
  !     * AUG 16/93 - D.VERSEGHY/M.LAZARE. ADD MISSING OUTER LOOP ON "J"
  !     *                                  IN 200 LOOP.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. UNSATURATED FLOW OF WATER INTO SOIL.
  !
  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IG, IGP1, IGP2, ILG, IL1, IL2, JL, N
  integer :: I, J
  !
  !     * OUTPUT FIELDS.
  !
  real, intent(inout) :: WMOVE  (ILG,IGP2) !< Water movement matrix \f$[m^3 m^{-2}]\f$
  real, intent(out)   :: TMOVE  (ILG,IGP2) !< Temperature matrix associated with ground water movement [C}
  real, intent(inout) :: ZF     (ILG)      !< Depth of the wetting front [m]
  real, intent(inout) :: TRMDR  (ILG)      !< Remainder of time step after unsaturated infiltration ceases [s]
  integer, intent(inout) :: LZF (ILG)      !< Index of soil layer in which wetting front is located
  integer, intent(out)   :: NINF(ILG)      !< Number of levels involved in water movement
  !
  !     * INPUT FIELDS.
  !
  real, intent(in) :: R     (ILG)      !< Rainfall rate at ground surface \f$[m s^{-1}]\f$
  real, intent(in) :: TR    (ILG)      !< Temperature of rainfall [C]
  real, intent(in) :: PSIF  (ILG,IGP1) !< Soil water suction across the wetting front [m]
  real, intent(in) :: GRKINF(ILG,IGP1) !< Hydraulic conductivity of soil behind the wetting front \f$[m s^{-1}]\f$
  real, intent(in) :: THLINF(ILG,IGP1) !< Volumetric liquid water content behind the wetting front \f$[m^3 m^{-3}]\f$
  real, intent(in) :: THLIQX(ILG,IGP1) !< Volumetric liquid water content of soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in) :: TBARWX(ILG,IGP1) !< Temperature of water in soil layer [C]
  real, intent(in) :: DELZX (ILG,IGP1) !< Permeable depth of soil layer \f$[m] (\Delta z_{g,w})\f$
  real, intent(in) :: ZBOTX (ILG,IGP1) !< Depth of bottom of soil layer [m]
  integer, intent(in) :: IFILL (ILG)      !< Flag indicating whether unsaturated infiltration is occurring
  !
  !
  !     * INTERNAL WORK FIELDS.
  !
  real, intent(inout)    :: DZF  (ILG), TIMPND(ILG), WADJ(ILG), WADD(ILG)
  integer, intent(inout) :: IFIND(ILG)
  !
  !     * TEMPORARY VARIABLE.
  !
  real :: THLADD
  !-----------------------------------------------------------------------
  !>
  !! The infiltration rate \f$F_{inf}\f$ under conditions of a constant water
  !! supply can be expressed, e.g. in Mein and Larson (1973), as
  !!
  !! \f$F_{inf} = \f$K_{inf}\f$ [(\Psi_f + z_f)/ z_f ]\f$
  !!
  !! where \f$K_{inf}\f$ is the hydraulic conductivity of the soil behind the
  !! wetting front, \f$\Psi_f\f$ is the soil moisture suction across the
  !! wetting front, and \f$z_f\f$ is the depth of the wetting front. It can
  !! be seen that \f$F_{inf}\f$ decreases with increasing \f$z_f\f$ to an asymptotic
  !! value of \f$K_{inf}\f$. Thus, if the rainfall rate r is less than
  !! \f$K_{inf}\f$, the actual infiltration rate is limited by r, i.e.
  !! \f$F_{inf} = r\f$. Otherwise, \f$F_{inf}\f$ will be equal to R until the right-hand
  !! side of the above equation becomes less than r, after which point
  !! the above equation applies and ponding of excess water begins on
  !! the surface. The depth of the wetting front at this time \f$t_p\f$ can
  !! be calculated by setting \f$F_{inf}\f$ equal to r in the above equation
  !! and solving for \f$z_f\f$. This results in:
  !!
  !! \f$z_f = \Psi_f /[r/K_{inf} – 1]\f$
  !!
  !! The amount of water added to the soil up to the time of ponding
  !! is \f$t_p r\f$, or \f$z_f (\theta_{inf} – \theta_l)\f$, where \f$\theta_l\f$ and \f$\theta_{inf}\f$ are
  !! respectively the liquid water content of the soil before and
  !! after the wetting front has passed. Setting these two equal and
  !! solving for \f$t_p\f$ results in
  !!
  !! \f$t_p = z_f (\theta_{inf} – \theta_l)/r\f$
  !!
  !     * INITIALIZATION.
  !
  do I = IL1,IL2
    IFIND(I) = 0
    WADJ(I) = 0.
  end do ! loop 50
  !
  !     * TEST SUCCESSIVE SOIL LAYERS TO FIND DEPTH OF WETTING FRONT
  !     * AT THE TIME PONDING BEGINS, I.E. AT THE TIME THE DECREASING
  !     * INFILTRATION RATE EQUALS THE RAINFALL RATE.
  !
  !>
  !! In the 100 loop, a check is done for each successive soil layer
  !! to compare the infiltration rate in the layer with the rainfall
  !! rate. If \f$K_{inf} < R\f$, a test calculation is performed to determine
  !! where the depth of the wetting front would theoretically occur at
  !! the ponding time \f$t_p\f$. If the calculated value of \f$z_f\f$ is less than
  !! the depth of the top of the soil layer, \f$z_f\f$ is set to the depth of
  !! the top of the layer; if \f$z_f\f$ falls within the soil layer, that
  !! value of \f$z_f\f$ is accepted. In both cases, the index LZF is set to
  !! the index of the layer, and the flag IFIND, indicating that \f$z_f\f$,
  !! has been successfully located, is set to 1. If the infiltration
  !! rate in the soil layer is greater than the rainfall rate, \f$z_f\f$ is
  !! provisionally set to the bottom of the current layer, and LZF to
  !! the index of the next layer. IFIND remains zero. If the
  !! infiltration rate in the layer is vanishingly small, \f$z_f\f$ is set to
  !! the depth of the top of the current layer, LZF to the index of
  !! the overlying layer, and IFIND to 1.
  !!
  do J = 1,IGP1 ! loop 100
    do I = IL1,IL2
      if (ifILL(I) > 0 .and. ifIND(I) == 0) then
        if (GRKINF(I,J) > 1.0E-12 .and. &
            GRKINF(I,J) < (R(I) - 1.0E-8)) then
          ZF(I) = PSIF (I,J) / (R(I) / GRKINF(I,J) - 1.0)
          if (ZF(I) < (ZBOTX(I,J) - DELZX(I,J))) then
            ZF(I) = MAX(ZBOTX(I,J) - DELZX(I,J),0.0)
            LZF(I) = J
            IFIND(I) = 1
          else if (ZF(I) < ZBOTX(I,J)) then
            LZF(I) = J
            IFIND(I) = 1
          end if
        else if (GRKINF(I,J) > 1.0E-12) then
          ZF(I) = ZBOTX(I,J)
          LZF(I) = MIN(J + 1,IGP1)
        else if (GRKINF(I,J) <= 1.0E-12) then
          if (J == 1) then
            ZF(I) = 0.0
            LZF(I) = 1
          else
            ZF(I) = ZBOTX(I,J - 1)
            LZF(I) = J - 1
          end if
          IFIND(I) = 1
        end if
      end if
    end do
  end do ! loop 100
  !>
  !! If LZF is greater than 1, some adjustment to the equation for \f$t_p\f$
  !! above is required to account for the fact that the values of
  !! \f$\theta_{inf}\f$ and \f$\theta_l\f$ in the layer containing the wetting front may
  !! differ from those in the overlying layers. The equation for \f$t_p\f$
  !! above can be rewritten as
  !!
  !! \f$t_p = [z_f [\theta_{inf}(z_f) – \theta_l(z_f)] + w_{adj}]/r\f$
  !!
  !! where \f$w_{adj}\f$ is calculated as
  !!
  !! \f$LZF-1\f$
  !!
  !! \f$w_{adj} = \Sigma[(\theta_{inf. i} – \theta_{l, i} ) – (\theta_{inf} (z_f) – \theta_l (z_f))] \Delta z_{g, w}\f$
  !!
  !! \f$i=1\f$
  !!
  !! The adjusting volume WADJ is calculated in loop 200, and the time
  !! to ponding TIMPND in loop 250. If TIMPND is greater than the
  !! amount of time remaining in the current time step TRMDR, then
  !! unsaturated infiltration is deemed to be occurring over the
  !! entire time step. In this case, the amount of water infiltrating
  !! over the time step is assigned to the first level of the water
  !! movement matrix WMOVE and to the accounting variable WADD, and
  !! the temperature of the infiltrating water is assigned to the
  !! first level of the matrix TMOVE.
  !!
  !
  !     * FIND THE VOLUME OF WATER NEEDED TO CORRECT FOR THE DIFFERENCE
  !     * (IF ANY) BETWEEN THE LIQUID MOISTURE CONTENTS OF THE LAYERS
  !     * OVERLYING THE WETTING FRONT AND THAT OF THE LAYER CONTAINING
  !     * THE WETTING FRONT.
  !
  do J = 1,IGP1 ! loop 200
    do I = IL1,IL2
      if (ifILL(I) > 0 .and. LZF(I) > 1 .and. J < LZF(I)) then
        WADJ(I) = WADJ(I) + DELZX(I,J) * ( (THLINF(I,J) - THLIQX(I,J)) - &
                  (THLINF(I,LZF(I)) - THLIQX(I,LZF(I))) )
      end if
    end do
  end do ! loop 200
  !
  !     * CALCULATE THE TIME TO PONDING, GIVEN THE DEPTH REACHED BY THE
  !     * WETTING FRONT AT THAT TIME.

  do I = IL1,IL2 ! loop 250
    if (ifILL(I) > 0) then
      TIMPND(I) = (ZF(I) * (THLINF(I,LZF(I)) - THLIQX(I,LZF(I))) + &
                  WADJ(I)) / R(I)
      TIMPND(I) = MAX(TIMPND(I),0.0)
      if (ZF(I) > 10.0) TIMPND(I) = 1.0E+8
      !
      !     * IN THE CASE WHERE THE TIME TO PONDING EXCEEDS OR EQUALS THE
      !     * TIME REMAINING IN THE CURRENT MODEL STEP, RECALCULATE THE
      !     * ACTUAL DEPTH ATTAINED BY THE WETTING FRONT OVER THE CURRENT
      !     * MODEL STEP; ASSIGN VALUES IN THE WATER MOVEMENT MATRIX.
      !
      if (TIMPND(I) >= TRMDR(I)) then
        TMOVE(I,1) = TR(I)
        WMOVE(I,1) = R(I) * TRMDR(I)
        WADD(I) = WMOVE(I,1)
      end if
    end if
  end do ! loop 250
  !
  !>
  !! In loop 300 WADD is partitioned over the soil profile by
  !! comparing in turn the liquid water content of each soil layer
  !! with the calculated liquid water content behind the wetting front
  !! THLINF, and decrementing WADD layer by layer until a layer is
  !! reached in which the remainder of WADD is insufficient to raise
  !! the liquid water content to THLINF. If this condition is reached,
  !! LZF is set to the index of the soil layer; the depth of the
  !! wetting front DZF within the layer, obtained as
  !! WADD/(THLINF-THLIQX), is added to the depth of the bottom of the
  !! overlying layer to obtain ZF.
  !!
  do J = 1,IGP1 ! loop 300
    do I = IL1,IL2
      if (ifILL(I) > 0) then
        if (TIMPND(I) >= TRMDR(I) .and. WADD(I) > 0.) then
          THLADD = MAX(THLINF(I,J) - THLIQX(I,J),0.0)
          if (THLADD > 0.) then
            DZF(I) = WADD(I) / THLADD
          else
            DZF(I) = 1.0E+8
          end if
          if (DZF(I) > (DELZX(I,J) + 1.0E-5)) then
            WADD(I) = WADD(I) - THLADD * DELZX(I,J)
          else
            DZF(I) = MIN(DZF(I),DELZX(I,J))
            LZF(I) = J
            if (J == 1) then
              ZF(I) = DZF(I)
            else
              ZF(I) = ZBOTX(I,J - 1) + DZF(I)
            end if
            WADD(I) = 0.
          end if
        end if
      end if
    end do
  end do ! loop 300
  !
  !>
  !! In loop 400, the water content in each soil layer J existing
  !! above ZF is assigned to the J+1 level of the water movement
  !! matrix WMOVE, and the respective water temperatures are assigned
  !! to TMOVE.
  !!
  do J = 1,IGP1 ! loop 400
    do I = IL1,IL2
      if (ifILL(I) > 0) then
        if (TIMPND(I) >= TRMDR(I) .and. J <= LZF(I)) then
          TMOVE(I,J + 1) = TBARWX(I,J)
          if (J == LZF(I) .and. DZF(I) < DELZX(I,J)) then
            WMOVE(I,J + 1) = THLIQX(I,J) * DZF(I)
          else
            WMOVE(I,J + 1) = THLIQX(I,J) * DELZX(I,J)
          end if
        end if
      end if
    end do
  end do ! loop 400
  !>
  !! If TIMPND < TRMDR, the amount of water infiltrating between the
  !! start of the time step and TIMPND is again assigned to the first
  !! level of the water movement matrix WMOVE, and the temperature of
  !! the infiltrating water is assigned to the first level of the
  !! matrix TMOVE. The depth DZF of the wetting front within the layer
  !! containing it is calculated by subtracting the depth of the
  !! bottom of the overlying layer from ZF.
  !!
  !
  !     * IN THE CASE WHERE THE TIME TO PONDING IS LESS THAN THE TIME
  !     * REMAINING IN THE CURRENT MODEL STEP, ACCEPT THE DEPTH OF THE
  !     * WETTING FRONT FROM LOOP 100; ASSIGN VALUES IN THE WATER
  !     * MOVEMENT MATRIX.
  !
  do I = IL1,IL2 ! loop 450
    if (ifILL(I) > 0) then
      if (TIMPND(I) < TRMDR(I)) then
        TMOVE(I,1) = TR(I)
        WMOVE(I,1) = R(I) * TIMPND(I)
        if (LZF(I) == 1) then
          DZF(I) = ZF(I)
        else
          DZF(I) = ZF(I) - ZBOTX(I,LZF(I) - 1)
        end if
      end if
    end if
  end do ! loop 450
  !
  !>
  !! In loop 500, the water content in each soil layer J existing
  !! above ZF is assigned to the J+1 level of the water movement
  !! matrix WMOVE, and the respective water temperatures are assigned
  !! to TMOVE.
  !!
  do J = 1,IGP1 ! loop 500
    do I = IL1,IL2
      if (ifILL(I) > 0) then
        if (TIMPND(I) < TRMDR(I) .and. J <= LZF(I)) then
          TMOVE(I,J + 1) = TBARWX(I,J)
          if (J == LZF(I)) then
            WMOVE(I,J + 1) = THLIQX(I,J) * DZF(I)
          else
            WMOVE(I,J + 1) = THLIQX(I,J) * DELZX(I,J)
          end if
        end if
      end if
    end do
  end do ! loop 500
  !>
  !! Finally, the time remaining in the current time step after the
  !! period of unsaturated infiltration is recalculated, and the
  !! counter NINF is set to LZF+1.
  !!
  !
  !     * CALCULATE TIME REMAINING IN CURRENT MODEL STEP AFTER
  !     * UNSATURATED FLOW.
  !
  do I = IL1,IL2 ! loop 600
    if (ifILL(I) > 0) then
      if (TIMPND(I) >= TRMDR(I)) then
        TRMDR(I) = 0.
      else
        TRMDR(I) = TRMDR(I) - TIMPND(I)
      end if
      NINF(I) = LZF(I) + 1
    end if
  end do ! loop 600
  !
  return
end subroutine waterInfiltrateUnsat
