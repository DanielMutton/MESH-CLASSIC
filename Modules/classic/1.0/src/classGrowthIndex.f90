!> \file
!! Evaluates growth index used in calculating vegetation
!! parameters for forests.
!! @author D. Verseghy, M. Lazare, V. Arora
!>
!! The growth index that is calculated here varies from a value of 1
!! for periods when the trees are mature and fully leaved, to 0 for
!! dormant and leafless periods, with a linear transition between
!! the two. The transition periods are assumed to last for sixty
!! days; therefore during these periods the growth index is
!! incremented by \f$ \Delta t /5.184x10^6 \f$ where \f$\Delta t\f$ is the time step in
!! seconds.
!!
!! The transition period from dormant to fully leafed is triggered
!! when both the air temperature and the temperature of the first
!! soil layer are above 2 C. If one of these conditions is not met
!! afterwards, the growth index is reset back to 0. Increments are
!! added continuously thereafter until the index reaches 1.
!!
!! The transition from fully leafed to dormant is triggered when
!! either the air temperature or the temperature of the first soil
!! layer falls below 2 C. When this first happens at the end of the
!! fully-leafed period, the growth index is set instantaneously to
!! -1 and increments are continuously added from that point until
!! the index reaches 0.
!!
!! The absolute value of this growth index is utilized for
!! performing calculations of various forest vegetation parameters
!! in subroutine calcLandSurfParams; thus its shape as used there is that of a
!! symmetrical trapezoidal function.
!!
subroutine classGrowthIndex (GROWTH, TBAR, TA, FC, FCS, ILG, IG, IL1, IL2, JL) ! Formerly CGROW


  !     * MAR 09/07 - D.VERSEGHY. CHANGE SENESCENCE THRESHOLD FROM
  !     *                         0.10 TO 0.90.
  !     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * DEC 03/03 - V.ARORA. CHANGE THRESHOLD OF LEAF-OUT FROM 0 C TO
  !     *                      2 C.
  !     * SEP 25/97 - M.LAZARE. CLASS - VERSION 2.7.
  !     *                       INSERT "IF" CONDITION TO PERFORM THESE
  !     *                       CALCULATIONS ONLY IF CANOPY IS PRESENT.
  !     * APR 24/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
  !     *                                  FOR MODEL VERSION GCM7.
  !     * APR 11/89 - D.VERSEGHY. INCREMENT/DECREMENT GROWTH INDEX FOR
  !     *                         VEGETATION TYPES 1 AND 2 (NEEDLELEAF
  !     *                         AND BROADLEAF TREES).
  !
  use classicParams,       only : DELT, TFREZ

  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(in) :: ILG, IG, IL1, IL2, JL
  integer :: I
  !
  !     * OUTPUT ARRAY.
  !
  real, intent(inout) :: GROWTH(ILG)  !< Tree growth index [ ]
  !
  !     * INPUT ARRAYS.
  !
  real, intent(in) :: TBAR  (ILG,IG)   !< Temperature of soil layers [K]
  !
  real, intent(in) :: TA(ILG)  !< Air temperature [K]
  real, intent(in) :: FC(ILG)  !< Fractional coverage of vegetation without snow on
  !! modelled area [ ]
  real, intent(in) :: FCS(ILG) !< Fractional coverage of vegetation with underlying
  !! snow pack on modelled area [ ]
  !-----------------------------------------------------------------------

  do I = IL1,IL2 ! loop 100
    if ((FC(I) + FCS(I)) > 0.0) then
      if (GROWTH(I) < 0.0) then
        GROWTH(I) = MIN(0.0,(GROWTH(I) + DELT / 5.184E6))
      else
        if (TA(I) > (TFREZ + 2.0) .and. TBAR(I,1) > (TFREZ + 2.0)) &
            then
          GROWTH(I) = MIN(1.0,(GROWTH(I) + DELT / 5.184E6))
        else
          if (GROWTH(I) > 0.90) then
            GROWTH(I) = - GROWTH(I)
          else
            GROWTH(I) = 0.0
          end if
        end if
      end if
    end if
  end do ! loop 100
  !
  return
end subroutine classGrowthIndex
