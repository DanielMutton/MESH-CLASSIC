!> \file
!! Evaluates infiltration of water into soil under saturated
!! conditions.
!! @author D. Verseghy, M. Lazare, Y. Delage, J. P. Blanchette
!
subroutine waterInfiltrateSat (WMOVE, TMOVE, LZF, NINF, TRMDR, TPOND, ZPOND, & ! Formerly WFLOW
                               R, TR, EVAP, PSIF, GRKINF, THLINF, THLIQX, TBARWX, &
                               DELZX, ZBOTX, FMAX, ZF, DZF, DTFLOW, THLNLZ, &
                               THLQLZ, DZDISP, WDISP, WABS, ITER, NEND, ISIMP, &
                               IGRN, IG, IGP1, IGP2, ILG, IL1, IL2, JL, N)

  !     * FEB 09/09 - J.P.BLANCHETTE. INCREASE LIMITING VALUE OF NEND.
  !     * JAN 06/09 - D.VERSEGHY. CHECKS ON WETTING FRONT LOCATION
  !     *                         IN 500 LOOP.
  !     * AUG 07/07 - D.VERSEGHY. INCREASE ITERATION COUNTER NEND
  !     *                         MOVE CALCULATION OF FMAX FROM
  !     *                         waterFlowInfiltrate TO THIS ROUTINE.
  !     * MAY 17/06 - D.VERSEGHY. PROTECT AGAINST DIVISIONS BY ZERO.
  !     * SEP 13/05 - D.VERSEGHY. REPLACE HARD-CODED 4 WITH IGP1.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * JUL 30/04 - D.VERSEGHY/Y.DELAGE. PROTECT SENSITIVE
  !     *                         CALCULATIONS AGAINST ROUNDOFF ERRORS.
  !     * JUN 21/02 - D.VERSEGHY. UPDATE SUBROUTINE CALL.
  !     * MAR 04/02 - D.VERSEGHY. DEFINE "NEND" FOR ALL CASES.
  !     * DEC 16/94 - D.VERSEGHY. BUG FIX - SPECIFY TMOVE BEHIND
  !     *                         WETTING FRONT AFTER ANY FULL-LAYER
  !     *                         JUMP DOWNWARD.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
  !     *                                  REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
  !     *                         CLASS VERSION 2.0 (WITH CANOPY).
  !     * APR 11/89 - D.VERSEGHY. SATURATED FLOW OF WATER THROUGH SOIL.
  !
  use classicParams, only : DELT

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: IG, IGP1, IGP2, ILG, IL1, IL2, JL, N
  integer             :: I, J, NPNTS, NIT
  !
  !     * INPUT/OUTPUT FIELDS.
  !
  real, intent(inout)    :: WMOVE (ILG,IGP2) !< Water movement matrix \f$[m^3 m^{-2}]\f$
  real, intent(inout)    :: TMOVE (ILG,IGP2) !< Temperature matrix associated with ground water movement [C]
  integer, intent(inout) :: LZF   (ILG)      !< Index of soil layer in which wetting front is located
  integer, intent(inout) :: NINF  (ILG)      !< Number of levels involved in water movement
  real, intent(inout)    :: TRMDR (ILG)      !< Time remaining in current time step [s]
  real, intent(inout)    :: TPOND (ILG)      !< Temperature of ponded water [C]
  real, intent(inout)    :: ZPOND (ILG)      !< Depth of ponded water \f$[m] (z_p)\f$
  !
  !     * INPUT FIELDS.
  !
  real, intent(in)    :: R     (ILG)      !< Rainfall rate at ground surface \f$[m s^{-1}]\f$
  real, intent(in)    :: TR    (ILG)      !< Temperature of rainfall [C]
  real, intent(in)    :: EVAP  (ILG)      !< Surface evaporation rate \f$[m s^{-1}]\f$
  real, intent(in)    :: PSIF  (ILG,IGP1) !< Soil water suction across the wetting front \f$[m] (\Psi_f)\f$
  real, intent(in)    :: GRKINF(ILG,IGP1) !< Hydraulic conductivity of soil behind the wetting front \f$[m s^{-1}] (K)\f$
  real, intent(in)    :: THLINF(ILG,IGP1) !< Volumetric liquid water content behind the wetting front \f$[m^3 m^{-3}]\f$
  real, intent(in)    :: THLIQX(ILG,IGP1) !< Volumetric liquid water content of soil layer \f$[m^3 m^{-3}]\f$
  real, intent(in)    :: TBARWX(ILG,IGP1) !< Temperature of water in soil layer [C]
  real, intent(in)    :: DELZX (ILG,IGP1) !< Permeable depth of soil layer \f$[m] (\Delta z_{z,w})\f$
  real, intent(in)    :: ZBOTX (ILG,IGP1) !< Depth of bottom of soil layer [m]
  real, intent(inout) :: FMAX  (ILG)      !< Maximum infiltration rate, defined as minimum value of GRKINF
  real, intent(inout) :: ZF    (ILG)      !< Depth of the wetting front \f$[m] (z_f)\f$
  integer, intent(in) :: IGRN  (ILG)      !< Flag to indicate whether infiltration is occurring
  !
  !     * INTERNAL WORK FIELDS.
  !
  real, intent(inout)    :: DZF(ILG), DTFLOW(ILG), THLNLZ(ILG), THLQLZ(ILG), &
                            DZDISP(ILG), WDISP(ILG), WABS(ILG)
  integer, intent(inout) :: ITER(ILG), NEND(ILG), ISIMP(ILG)
  !
  !     * TEMPORARY VARIABLES.
  !
  real :: RESID, FINF, ZPTEST, WINF
  !
  !-----------------------------------------------------------------------
  !>
  !! General calculations performed: The infiltration rate \f$F_{inf}\f$ under saturated conditions is calculated as
  !!
  !! \f$F_{inf} = K_{inf} [(\Psi_f + z_f + z_p)/z_f ]\f$
  !!
  !! where \f$K_{inf}\f$ is the hydraulic conductivity of the soil behind the wetting front, \f$\Psi_f\f$ is the soil moisture
  !! suction across the wetting front, \f$z_f\f$ is the depth of the wetting front, and \f$z_p\f$ is the depth of water ponded
  !! on the surface. Since \f$\Psi_f\f$ will vary by soil layer, and \f$z_f\f$ and \f$z_p\f$ will vary with time, the period of
  !! infiltration is divided into two-minute segments. A maximum iteration flag NEND is defined as the number of
  !! iteration segments plus 100, as a safeguard against runaway iterations. Since the surface infiltration rate
  !! will be limited by the lowest infiltration rate encountered, a maximum infiltration rate FMAX is defined as
  !! the minimum value of GRKINF, the hydraulic conductivity behind the wetting front, in all of the soil layers
  !! above and including that containing the wetting front.
  !!
  !
  !     * CALCULATE ITERATION ENDPOINT "NEND", AND SET SWITCH "ITER" TO 1
  !     * FOR POINTS OVER WHICH THIS SUBROUTINE IS TO BE PERFORMED.
  !
  do I = IL1,IL2 ! loop 50
    if (IGRN(I) > 0 .and. TRMDR(I) > 0.0) then
      RESID = MOD(TRMDR(I),120.)
      if (RESID > 0.) then
        NEND(I) = NINT(TRMDR(I) / 120. + 0.5) + 100
      else
        NEND(I) = NINT(TRMDR(I) / 120.) + 100
      end if
      ITER(I) = 1
      FMAX(I) = 999999.
    else
      NEND(I) = 0
      ITER(I) = 0
    end if
  end do ! loop 50
  NIT = 1
  !
  do J = 1,IGP1 ! loop 75
    do I = IL1,IL2
      if (ITER(I) > 0 .and. LZF(I) > 1 .and. J < LZF(I)) then
        FMAX(I) = MIN(GRKINF(I,J),FMAX(I))
      end if
    end do
  end do ! loop 75
  !
  !ignoreLint(1)
100 continue ! anchor for a GO TO statement that needs to be removed

  !     * BEGINNING OF ITERATION SEQUENCE.
  !     * SET OR RESET NUMBER OF POINTS TO BE PROCESSED ON THE CURRENT
  !     * LATITUDE CIRCLE(S).
  !
  NPNTS = 0
  !>
  !! In loop 200, a check is performed to determine whether the
  !! liquid water content of the current soil layer, THLIQX, equals
  !! or exceeds that behind the wetting front, THLINF. If so, the
  !! depth of the wetting front is reset to the depth of the soil
  !! layer. The water content of the soil layer is assigned to the
  !! level of the water movement matrix WMOVE corresponding to the
  !! soil layer index plus 1, and the temperature of the water in the
  !! layer is assigned to the same level of the matrix TMOVE. FMAX is
  !! recalculated, and the flags LZF and NINF are each incremented
  !! by 1. The flag ISIMP is set to 1, indicating that the following
  !! calculations are to be bypassed for this iteration.
  !!
  !! If ISIMP is not 1, the time period DTFLOW applying to the
  !! current iteration loop is set to two minutes or to the
  !! remainder of the time step, whichever is less. If GRKINF for
  !! the current soil layer is not vanishingly small, the flag ISIMP
  !! is set to -2. Otherwise, it is deemed that infiltration is
  !! suppressed, and the temperature and depth of water ponded on
  !! the surface are simply updated. The temperature of the ponded
  !! water is calculated as the weighted average of the current pond
  !! temperature and the rainfall added to it. The ponded water
  !! remaining after rainfall and evaporation have taken place is
  !! calculated as ZPTEST. If ZPTEST is less than zero, it is
  !! deduced that evaporation must have removed the ponded water
  !! before the end of the current iteration loop. If this is the
  !! case, the time period of the current iteration loop is
  !! recalculated as the amount of time required for the difference
  !! between the evaporation
  !! and the rainfall rates to consume the ponded water, and the
  !! pond depth ZPOND is set to zero. Otherwise, ZPOND is set to
  !! ZPTEST. Finally, ISIMP is set to -1.
  !!

  !
  !     * IF THE WATER CONTENT OF THE CURRENT SOIL LAYER EQUALS OR EXCEEDS
  !     * THE WATER CONTENT BEHIND THE WETTING FRONT, INSTANTANEOUSLY
  !     * RELOCATE WETTING FRONT TO BOTTOM OF CURRENT SOIL LAYER
  !     * RE-EVALUATE INFILTRATION PARAMETERS, UPDATE WATER MOVEMENT
  !     * MATRIX, SET SWITCH "ISIMP" TO 1 AND DROP TO END OF ITERATION
  !     * LOOP.
  !     * (SOME ARRAYS ARE GATHERED (ON LZF) TO AVOID MULTIPLE INDIRECT-
  !     * ADDRESSING REFERENCES IN THE ENSUING LOOPS.)
  !
  do I = IL1,IL2 ! loop 200
    if (ITER(I) == 1) then
      THLNLZ(I) = THLINF(I,LZF(I))
      THLQLZ(I) = THLIQX(I,LZF(I))
      if (THLQLZ(I) > (THLNLZ(I) - 1.0E-6) .and. LZF(I) < IGP1 &
          .and. THLQLZ(I) > 0.0001) then
        ZF(I) = ZBOTX(I,LZF(I))
        WMOVE(I,NINF(I)) = THLQLZ(I) * DELZX(I,LZF(I))
        TMOVE(I,NINF(I)) = TBARWX(I,LZF(I))
        FINF = GRKINF(I,LZF(I)) * (ZF(I) + ZPOND(I) + PSIF (I,LZF(I))) / &
               ZF(I)
        FMAX(I) = MIN(FMAX(I),FINF)
        LZF(I) = LZF(I) + 1
        NINF(I) = NINF(I) + 1
        ISIMP(I) = 1
      else
        ISIMP(I) = 0
      end if
    else
      ISIMP(I) = 0
    end if
  end do ! loop 200
  !
  !     * INFILTRATION CALCULATIONS TAKING FINITE TIME. SET TIMESTEP OF
  !     * CURRENT ITERATION PASS AND CHECK HYDRAULIC CONDUCTIVITY OF
  !     * CURRENT SOIL LAYER. IF ZERO, RECALCULATE POND DEPTH AND POND
  !     * TEMPERATURE AND SET "ISIMP" TO -1; ELSE SET "ISIMP" TO -2.
  !
  do I = IL1,IL2 ! loop 300
    if (ITER(I) == 1 .and. ISIMP(I) /= 1) then
      DTFLOW(I) = MIN(TRMDR(I),120.)
      if (GRKINF(I,LZF(I)) > 1.0E-12) then
        ISIMP(I) = - 2
      else
        TPOND(I) = (ZPOND(I) * TPOND(I) + R(I) * DTFLOW(I) * TR(I)) / &
                   (ZPOND(I) + R(I) * DTFLOW(I))
        if (ABS(R(I) - EVAP(I)) > 1.0E-11) then
          ZPTEST = ZPOND(I) + (R(I) - EVAP(I)) * DTFLOW(I)
        else
          ZPTEST = ZPOND(I)
        end if
        if (ZPTEST < 0.) then
          DTFLOW(I) = ZPOND(I) / (EVAP(I) - R(I))
          ZPOND(I) = 0.0
        else
          ZPOND(I) = ZPTEST
        end if
        ISIMP(I) = - 1
      end if
    end if
  end do ! loop 300
  !>
  !! The 400 loop addresses the infiltration process under saturated conditions. Such infiltration is modelled
  !! as “piston flow”. First the current infiltration rate FINF is calculated using the equation given above.
  !! If the wetting front has passed the bottom of the lowest soil layer, PSIF Cis neglected. If FINF is greater
  !! than the rainfall rate R, FINF is set equal to R; if FINF is greater than CFMAX, FINF is set equal to FMAX.
  !! If the wetting front has not passed the bottom of the lowest soil layer, the change in depth of the wetting
  !! front DZF over the current time interval is calculated from the amount of infiltrating water WINF and the
  !! volumetric water content behind the wetting front, THLINF. The amount of soil water WDISP displaced by this
  !! movement of the wetting front is calculated from DZF and THLIQX, and the depth to which this water penetrates,
  !! DZDISP, is calculated as WDISP/(THLINF – THLIQX). The amount of soil water WABS entrained by the movement of
  !! WDISP itself is calculated as the product of DZDISP and THLIQX. The change in depth of the wetting front,
  !! behind which the liquid water content of the   soil layer is THLINF, is now the sum of the depth represented
  !! by infiltration of water at the surface, DZF, and the depth represented by displacement of the pre-existing
  !! soil water, DZDISP. If this change in depth causes the wetting front to move beyond the bottom of the current
  !! soil layer, DTFLOW is recalculated as the amount of time required for the composite wetting front to reach the
  !! bottom of the soil layer, and DZF, WDISP, DZDISP and WABS are likewise recalculated. As in the case for ISIMP = -1,
  !! the temperature of the ponded water on the surface is calculated as the weighted average of the current pond temperature
  !! and the rainfall added to it. The ponded water remaining after rainfall, infiltration and evaporation have taken place
  !! is calculated as ZPTEST. If ZPTEST is less than zero, it is deduced that infiltration and evaporation must have removed
  !! the ponded water before the end of the current iteration loop. If this is the case, the time period of the current iteration
  !! loop is recalculated as the amount of time required for the infiltration and evaporation minus the rainfall to
  !! consume the ponded water. The pond depth ZPOND is set to zero; if the wetting front has not passed the bottom of the
  !! lowest soil layer, DZF, WDISP, DZDISP and WABS are recalculated. Otherwise, ZPOND is set to ZPTEST. Finally, the first
  !! layer of the water movement matrix WMOVE is incremented with the amount of the infiltrated water, and the first layer
  !! of the matrix TMOVE with the temperature of the infiltrated water. The last layer of WMOVE is incremented by WDISP+WABS,
  !! and the depth of the wetting front by DZF+DZDISP.
  !!
  !     * "ISIMP"=-2: NORMAL SATURATED INFILTRATION UNDER PISTON-FLOW
  !     * CONDITIONS. CALCULATE CURRENT INFILTRATION RATE (FINF); WATER
  !     * INFILTRATING DURING CURRENT ITERATION PASS (WINF); SOIL WATER
  !     * DISPLACED INTO EMPTY PORES AHEAD OF WETTING FRONT (WDISP)
  !     * AND SOIL WATER OVERTAKEN BY DISPLACED SOIL WATER (WABS).
  !     * RE-EVALUATE POND TEMPERATURE AND POND DEPTH; UPDATE WATER
  !     * MOVEMENT MATRIX; ADJUST CURRENT POSITION OF WETTING FRONT.
  !
  do I = IL1,IL2 ! loop 400
    if (ITER(I) == 1 .and. ISIMP(I) == - 2) then
      if (LZF(I) < IGP1) then
        if (ZF(I) > 1.0E-7) then
          FINF = GRKINF(I,LZF(I)) * (ZF(I) + ZPOND(I) + &
                 PSIF (I,LZF(I))) / ZF(I)
        else
          FINF = GRKINF(I,1)
        end if
      else
        FINF = GRKINF(I,LZF(I)) * (ZF(I) + ZPOND(I)) / ZF(I)
      end if
      if (ZPOND(I) < 1.0E-8 .and. FINF > R(I))      FINF = R(I)
      if (FINF > FMAX(I)) FINF = FMAX(I)
      WINF = FINF * DTFLOW(I)
      if (LZF(I) < IGP1) then
        DZF(I) = WINF / THLNLZ(I)
        WDISP(I) = DZF(I) * THLQLZ(I)
        DZDISP(I) = WDISP(I) / (THLNLZ(I) - THLQLZ(I))
        WABS(I) = DZDISP(I) * THLQLZ(I)
        if ((ZF(I) + DZF(I) + DZDISP(I)) > ZBOTX(I,LZF(I))) then
          DTFLOW(I) = (ZBOTX(I,LZF(I)) - ZF(I)) / &
                      (FINF / THLNLZ(I) + &
                      (FINF * THLQLZ(I)) / &
                      (THLNLZ(I) * &
                      (THLNLZ(I) - THLQLZ(I))))
          WINF = FINF * DTFLOW(I)
          DZF(I) = WINF / THLNLZ(I)
          WDISP(I) = DZF(I) * THLQLZ(I)
          DZDISP(I) = WDISP(I) / (THLNLZ(I) - THLQLZ(I))
          WABS(I) = DZDISP(I) * THLQLZ(I)
        end if
      end if
      if (ZPOND(I) + R(I) * DTFLOW(I) > 1.0E-8) then
        TPOND(I) = (ZPOND(I) * TPOND(I) + R(I) * DTFLOW(I) * TR(I)) / &
                   (ZPOND(I) + R(I) * DTFLOW(I))
      end if
      if (ABS(R(I) - FINF - EVAP(I)) > 1.0E-11) then
        ZPTEST = ZPOND(I) + R(I) * DTFLOW(I) - WINF - EVAP(I) * &
                 DTFLOW(I)
      else
        ZPTEST = ZPOND(I)
      end if
      if (ZPTEST < 0.) then
        DTFLOW(I) = ZPOND(I) / (FINF + EVAP(I) - R(I))
        WINF = FINF * DTFLOW(I)
        ZPOND(I) = 0.0
        if (LZF(I) < IGP1) then
          DZF(I) = WINF / THLNLZ(I)
          WDISP(I) = DZF(I) * THLQLZ(I)
          DZDISP(I) = WDISP(I) / (THLNLZ(I) - THLQLZ(I))
          WABS(I) = DZDISP(I) * THLQLZ(I)
        end if
      else
        ZPOND(I) = ZPTEST
      end if
      if ((WMOVE(I,1) + WINF) > 0.) then
        TMOVE(I,1) = (WMOVE(I,1) * TMOVE(I,1) + WINF * TPOND(I)) / &
                     (WMOVE(I,1) + WINF)
      end if
      WMOVE(I,1) = WMOVE(I,1) + WINF
    end if
  end do ! loop 400
  !
  !     * (THIS PORTION OF THE ABOVE DO-LOOP WAS SPLIT OFF ON THE CRAY
  !     * BECAUSE IT WOULD NOT VECTORIZE. ONE MIGHT TRY AND RE-COMBINE
  !     * IT ON THE SX-3 (GOES IN FIRST PART OF IF BLOCK)).
  !
  do I = IL1,IL2 ! loop 450
    if (ITER(I) == 1 .and. ISIMP(I) == - 2 .and. &
        LZF(I) < IGP1) then
      WMOVE(I,NINF(I)) = WMOVE(I,NINF(I)) + WDISP(I) + WABS(I)
      ZF(I) = ZF(I) + DZF(I) + DZDISP(I)
    end if
  end do ! loop 450
  !>
  !! In the 500 loop, if ISIMP < 0 (i.e. water movement has occurred), the time remaining in the
  !! current time step is recalculated. If the wetting front is at a soil layer boundary, if the
  !! time remaining is non-zero, if the wetting front is still within the modelled soil column, and
  !! if the calculated water content behind the wetting front is greater than zero, FMAX is updated,
  !! and LZF and NINF are incremented by 1. The NINF level of the TMOVE matrix is set to the water
  !! temperature of the new soil layer.
  !!
  !     * CALCULATE REMAINING ITERATION TIME; RE-EVALUATE INFILTRATION
  !     * PARAMETERS.
  !
  do I = IL1,IL2 ! loop 500
    if (ITER(I) == 1 .and. ISIMP(I) /= 1) then
      TRMDR(I) = TRMDR(I) - DTFLOW(I)
      if (ABS(ZF(I) - ZBOTX(I,LZF(I))) < 1.0E-6 .and. &
          TRMDR(I) > 0. .and. LZF(I) < IGP1 .and. &
          THLQLZ(I) > 0.0001) then
        FINF = GRKINF(I,LZF(I)) * (ZBOTX(I,LZF(I)) + ZPOND(I) + &
               PSIF (I,LZF(I))) / ZBOTX(I,LZF(I))
        FMAX(I) = MIN(FMAX(I),FINF)
        LZF(I) = LZF(I) + 1
        NINF(I) = NINF(I) + 1
        TMOVE(I,NINF(I)) = TBARWX(I,LZF(I))
      end if
    end if
  end do ! loop 500
  !
  !     * INCREMENT ITERATION COUNTER ("NIT") AND SEE IF ANY POINTS STILL
  !     * REMAIN TO BE DONE (USING "NPNTS"). IF SO, RETURN TO BEGINNING
  !     * TO COMPLETE THESE REMAINING POINTS.
  !
  NIT = NIT + 1
  !>
  !! At the end of the iteration pass, checks are done in the 600 loop to ascertain whether the number
  !! of iteration passes is still less than NEND, whether either ponded water still exists or rain is
  !! still falling, and whether there is still time remaining in the current time step. If these conditions
  !! are all fulfilled, the counter NPNTS representing the number of points in the current vector of mosaic
  !! tiles for which infiltration is still occurring is incremented by 1. Otherwise, the iteration flag for
  !! the current tile is changed from 1 to 0, signaling the end of the saturated infiltration calculations
  !! for that tile.
  !!
  do I = IL1,IL2 ! loop 600
    if (IGRN(I) > 0) then
      if (NIT <= NEND(I) .and. ITER(I) == 1 .and. &
          (ZPOND(I) > 1.0E-8 .or. R(I) > 0.)  .and. &
          TRMDR(I) > 0.) then
        NPNTS = NPNTS + 1
      else
        ITER(I) = 0
      end if
    end if
  end do ! loop 600
  !
  if (NPNTS > 0) GO TO 100
  !
  return
end subroutine waterInfiltrateSat
